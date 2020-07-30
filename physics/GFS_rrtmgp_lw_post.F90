module GFS_rrtmgp_lw_post 
  use machine,                    only: kind_phys
  use GFS_typedefs,               only: GFS_coupling_type,  &
                                        GFS_control_type,   &
                                        GFS_grid_type,      &
                                        GFS_radtend_type,   &
                                        GFS_statein_type,   &
                                        GFS_diag_type
  use module_radiation_aerosols, only: NSPC1
  use module_radlw_parameters,   only: topflw_type, sfcflw_type, proflw_type
  ! RRTMGP DDT's
  use mo_gas_optics_rrtmgp,      only: ty_gas_optics_rrtmgp
  use mo_fluxes_byband,          only: ty_fluxes_byband
  use mo_heating_rates,          only: compute_heating_rate
  use rrtmgp_aux,                only: check_error_msg
  implicit none
  
  public GFS_rrtmgp_lw_post_init,GFS_rrtmgp_lw_post_run,GFS_rrtmgp_lw_post_finalize

contains
  ! #########################################################################################
  ! SUBROUTINE GFS_rrtmgp_lw_post_init
  ! #########################################################################################
  subroutine GFS_rrtmgp_lw_post_init()
  end subroutine GFS_rrtmgp_lw_post_init

 ! #########################################################################################
  ! SUBROUTINE GFS_rrtmgp_lw_post_run
  ! #########################################################################################
!> \section arg_table_GFS_rrtmgp_lw_post_run
!! \htmlinclude GFS_rrtmgp_lw_post.html
!!
  subroutine GFS_rrtmgp_lw_post_run (Model, Grid, Radtend, Coupling, Diag,  Statein, im,   &
       p_lev, tsfa, fluxlwUP_allsky, fluxlwDOWN_allsky, fluxlwUP_clrsky, fluxlwDOWN_clrsky,&
       raddt, aerodp, cldsa, mtopa, mbota, cld_frac, cldtaulw,                             &
       flxprf_lw, errmsg, errflg)

    ! Inputs
    type(GFS_control_type), intent(in) :: &
         Model                ! Fortran DDT: FV3-GFS model control parameters
    type(GFS_grid_type), intent(in) :: &
         Grid                 ! Fortran DDT: FV3-GFS grid and interpolation related data  
    type(GFS_statein_type), intent(in) :: &
         Statein              ! Fortran DDT: FV3-GFS prognostic state data in from dycore  
    integer, intent(in) :: &
         im                   ! Horizontal loop extent 
    real(kind_phys), dimension(size(Grid%xlon,1)), intent(in) ::  &
         tsfa                 ! Lowest model layer air temperature for radiation (K)
    real(kind_phys), dimension(size(Grid%xlon,1), Model%levs+1), intent(in) :: &
         p_lev                ! Pressure @ model layer-interfaces    (hPa)
    real(kind_phys), dimension(size(Grid%xlon,1), Model%levs+1), intent(in) :: &
         fluxlwUP_allsky,   & ! RRTMGP longwave all-sky flux (W/m2)
         fluxlwDOWN_allsky, & ! RRTMGP longwave all-sky flux (W/m2)
         fluxlwUP_clrsky,   & ! RRTMGP longwave clear-sky flux (W/m2)
         fluxlwDOWN_clrsky    ! RRTMGP longwave clear-sky flux (W/m2)
    real(kind_phys), intent(in) :: &
         raddt                ! Radiation time step
    real(kind_phys), dimension(im,NSPC1), intent(in) :: &
         aerodp               ! Vertical integrated optical depth for various aerosol species  
    real(kind_phys), dimension(im,5), intent(in) :: &
         cldsa                ! Fraction of clouds for low, middle, high, total and BL 
    integer,         dimension(im,3), intent(in) ::&
         mbota,             & ! vertical indices for low, middle and high cloud tops 
         mtopa                ! vertical indices for low, middle and high cloud bases
    real(kind_phys), dimension(im,Model%levs), intent(in) :: &
         cld_frac,          & ! Total cloud fraction in each layer
         cldtaulw             ! approx 10.mu band layer cloud optical depth  
    real(kind_phys),dimension(size(Grid%xlon,1), Model%levs)  :: &
         hlwc,              & ! Longwave all-sky heating-rate (K/sec)  
         hlw0                 ! Longwave clear-sky heating-rate (K/sec)

    ! Outputs (mandatory)
    character(len=*), intent(out) :: &
         errmsg
    integer, intent(out) :: &
         errflg
    type(GFS_coupling_type), intent(inout) :: &
         Coupling             ! Fortran DDT: FV3-GFS fields to/from coupling with other components 
    type(GFS_radtend_type), intent(inout) :: &
         Radtend              ! Fortran DDT: FV3-GFS radiation tendencies 
    type(GFS_diag_type), intent(inout) :: &
         Diag                 ! Fortran DDT: FV3-GFS diagnotics data 

    ! Outputs (optional)
    type(proflw_type), dimension(size(Grid%xlon,1), Model%levs+1), optional, intent(inout) :: &
         flxprf_lw            ! 2D radiative fluxes, components:
                              ! upfxc - total sky upward flux (W/m2)
                              ! dnfxc - total sky dnward flux (W/m2)
                              ! upfx0 - clear sky upward flux (W/m2)
                              ! dnfx0 - clear sky dnward flux (W/m2)

    ! Local variables
    integer :: i, j, k, iSFC, iTOA, itop, ibtc
    logical :: l_fluxeslw2d, top_at_1
    real(kind_phys) :: tem0d, tem1, tem2

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    if (.not. Model%lslwr) return

    ! Are any optional outputs requested?
    l_fluxeslw2d    = present(flxprf_lw)

    ! #######################################################################################
    ! What is vertical ordering?
    ! #######################################################################################
    top_at_1 = (p_lev(1,1) .lt. p_lev(1, Model%levs))
    if (top_at_1) then 
       iSFC = Model%levs+1
       iTOA = 1
    else
       iSFC = 1
       iTOA = Model%levs+1
    endif

    ! #######################################################################################
    ! Compute LW heating-rates. 
    ! #######################################################################################
    if (Model%lslwr) then
       ! Clear-sky heating-rate (optional)
       if (Model%lwhtr) then
          call check_error_msg('GFS_rrtmgp_post',compute_heating_rate(  &
               fluxlwUP_clrsky,   & ! IN  - RRTMGP upward longwave clear-sky flux profiles (W/m2)
               fluxlwDOWN_clrsky, & ! IN  - RRTMGP downward longwave clear-sky flux profiles (W/m2)
               p_lev,             & ! IN  - Pressure @ layer-interfaces (Pa)
               hlw0))               ! OUT - Longwave clear-sky heating rate (K/sec)
       endif
       ! All-sky heating-rate (mandatory)
       call check_error_msg('GFS_rrtmgp_post',compute_heating_rate(     &
            fluxlwUP_allsky,      & ! IN  - RRTMGP upward longwave all-sky flux profiles (W/m2)
            fluxlwDOWN_allsky,    & ! IN  - RRTMGP downward longwave all-sky flux profiles (W/m2)
            p_lev,                & ! IN  - Pressure @ layer-interfaces (Pa)
            hlwc))                  ! OUT - Longwave all-sky heating rate (K/sec)
       
       ! Copy fluxes from RRTGMP types into model radiation types.
       ! Mandatory outputs
       Diag%topflw(:)%upfxc    = fluxlwUP_allsky(:,iTOA)
       Diag%topflw(:)%upfx0    = fluxlwUP_clrsky(:,iTOA)
       Radtend%sfcflw(:)%upfxc = fluxlwUP_allsky(:,iSFC)
       Radtend%sfcflw(:)%upfx0 = fluxlwUP_clrsky(:,iSFC)
       Radtend%sfcflw(:)%dnfxc = fluxlwDOWN_allsky(:,iSFC)
       Radtend%sfcflw(:)%dnfx0 = fluxlwDOWN_clrsky(:,iSFC)
       
       ! Optional outputs
       if(l_fluxeslw2d) then
          flxprf_lw%upfxc = fluxlwUP_allsky
          flxprf_lw%dnfxc = fluxlwDOWN_allsky
          flxprf_lw%upfx0 = fluxlwUP_clrsky
          flxprf_lw%dnfx0 = fluxlwDOWN_clrsky
       endif
    endif
    
    ! #######################################################################################
    !  Save LW outputs.
    ! #######################################################################################
    if (Model%lslwr) then
       ! Save surface air temp for diurnal adjustment at model t-steps
       Radtend%tsflw (:) = tsfa(:)
       
       ! All-sky heating rate profile
       do k = 1, model%levs
          Radtend%htrlw(1:im,k) = hlwc(1:im,k)
       enddo
       if (Model%lwhtr) then
          do k = 1, model%levs
             Radtend%lwhc(1:im,k) = hlw0(1:im,k)
          enddo
       endif
       
       ! Radiation fluxes for other physics processes
       Coupling%sfcdlw(:) = Radtend%sfcflw(:)%dnfxc
    endif 

    ! #######################################################################################
    ! Save LW diagnostics
    ! - For time averaged output quantities (including total-sky and clear-sky SW and LW 
    !   fluxes at TOA and surface; conventional 3-domain cloud amount, cloud top and base 
    !   pressure, and cloud top temperature; aerosols AOD, etc.), store computed results in
    !   corresponding slots of array fluxr with appropriate time weights.
    ! - Collect the fluxr data for wrtsfc
    ! #######################################################################################
    if (Model%lssav) then
       if (Model%lslwr) then
          do i=1,im
             ! LW all-sky fluxes
             Diag%fluxr(i,1 ) = Diag%fluxr(i,1 ) + Model%fhlwr * fluxlwUP_allsky(  i,iTOA)   ! total sky top lw up
             Diag%fluxr(i,19) = Diag%fluxr(i,19) + Model%fhlwr * fluxlwDOWN_allsky(i,iSFC)   ! total sky sfc lw dn
             Diag%fluxr(i,20) = Diag%fluxr(i,20) + Model%fhlwr * fluxlwUP_allsky(  i,iSFC)   ! total sky sfc lw up
             ! LW clear-sky fluxes
             Diag%fluxr(i,28) = Diag%fluxr(i,28) + Model%fhlwr * fluxlwUP_clrsky(  i,iTOA)   ! clear sky top lw up
             Diag%fluxr(i,30) = Diag%fluxr(i,30) + Model%fhlwr * fluxlwDOWN_clrsky(i,iSFC)   ! clear sky sfc lw dn
             Diag%fluxr(i,33) = Diag%fluxr(i,33) + Model%fhlwr * fluxlwUP_clrsky(  i,iSFC)   ! clear sky sfc lw up
          enddo
          
          do i=1,im
             Diag%fluxr(i,17) = Diag%fluxr(i,17) + raddt * cldsa(i,4)
             Diag%fluxr(i,18) = Diag%fluxr(i,18) + raddt * cldsa(i,5)
          enddo
          
          ! Save cld frac,toplyr,botlyr and top temp, note that the order of h,m,l cloud is reversed for 
          ! the fluxr output. save interface pressure (pa) of top/bot
          do j = 1, 3
             do i = 1, IM
                tem0d = raddt * cldsa(i,j)
                itop  = mtopa(i,j)
                ibtc  = mbota(i,j)
                Diag%fluxr(i, 8-j) = Diag%fluxr(i, 8-j) + tem0d
                Diag%fluxr(i,11-j) = Diag%fluxr(i,11-j) + tem0d * Statein%prsi(i,itop)
                Diag%fluxr(i,14-j) = Diag%fluxr(i,14-j) + tem0d * Statein%prsi(i,ibtc)
                Diag%fluxr(i,17-j) = Diag%fluxr(i,17-j) + tem0d * Statein%tgrs(i,itop)
                
                ! Add optical depth and emissivity output
                tem2 = 0.
                do k=ibtc,itop
                   tem2 = tem2 + cldtaulw(i,k)      ! approx 10. mu channel
                enddo
                Diag%fluxr(i,46-j) = Diag%fluxr(i,46-j) + tem0d * (1.0-exp(-tem2))
             enddo
          enddo
       endif       
    endif

  end subroutine GFS_rrtmgp_lw_post_run

  ! #########################################################################################
  ! SUBROUTINE GFS_rrtmgp_lw_post_finalize
  ! #########################################################################################
  subroutine GFS_rrtmgp_lw_post_finalize ()
  end subroutine GFS_rrtmgp_lw_post_finalize

end module GFS_rrtmgp_lw_post

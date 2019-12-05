!>\file GFS_rrtmgp_sw_post
!!This file contains
module GFS_rrtmgp_sw_post 
  use machine,                    only: kind_phys
  use GFS_typedefs,               only: GFS_coupling_type,  &
                                        GFS_control_type,   &
                                        GFS_grid_type,      &
                                        GFS_radtend_type,   &
                                        GFS_diag_type,      &
                                        GFS_statein_type,   &
                                        GFS_interstitial_type
  use module_radiation_aerosols, only: NSPC1
  use module_radsw_parameters,   only: topfsw_type, sfcfsw_type, profsw_type, cmpfsw_type
  ! RRTMGP DDT's
  use mo_gas_optics_rrtmgp,      only: ty_gas_optics_rrtmgp
  use mo_fluxes_byband,          only: ty_fluxes_byband
  use mo_heating_rates,          only: compute_heating_rate
  use rrtmgp_aux,                only: check_error_msg
  implicit none
  
  public GFS_rrtmgp_sw_post_init,GFS_rrtmgp_sw_post_run,GFS_rrtmgp_sw_post_finalize

contains

  ! #########################################################################################
  ! SUBROUTINE GFS_rrtmgp_sw_post_init
  ! #########################################################################################
  subroutine GFS_rrtmgp_sw_post_init()
  end subroutine GFS_rrtmgp_sw_post_init

  ! #########################################################################################
  ! SUBROUTINE GFS_rrtmgp_sw_post_run
  ! #########################################################################################
!> \section arg_table_GFS_rrtmgp_sw_post_run
!! \htmlinclude GFS_rrtmgp_sw_post.html
!!
  subroutine GFS_rrtmgp_sw_post_run (Model, Interstitial, Grid, Diag, Radtend, Coupling,    &
       Statein, scmpsw, im, p_lev, sw_gas_props, nday, idxday, fluxswUP_allsky,             &
       fluxswDOWN_allsky, fluxswUP_clrsky, fluxswDOWN_clrsky, raddt, aerodp, cldsa, mbota,  &
       mtopa, cld_frac, cldtausw, flxprf_sw, hsw0, errmsg, errflg)

    ! Inputs
    type(GFS_control_type), intent(in) :: &
         Model                ! Fortran DDT: FV3-GFS model control parameters
    type(GFS_Interstitial_type), intent(in) :: &
         Interstitial         ! Fortran DDT: FV3-GFS interstitial arrays
    type(GFS_grid_type), intent(in) :: &
         Grid                 ! Fortran DDT: FV3-GFS grid and interpolation related data 
    type(GFS_coupling_type), intent(inout) :: &
         Coupling             ! Fortran DDT: FV3-GFS fields to/from coupling with other components 
    type(GFS_radtend_type), intent(inout) :: &
         Radtend              ! Fortran DDT: FV3-GFS radiation tendencies 
    type(GFS_diag_type), intent(inout) :: &
         Diag                 ! Fortran DDT: FV3-GFS diagnotics data  
    type(GFS_statein_type), intent(in) :: &
         Statein              ! Fortran DDT: FV3-GFS prognostic state data in from dycore  
    integer, intent(in) :: &
         im,                & ! Horizontal loop extent 
         nDay                 ! Number of daylit columns
    integer, intent(in), dimension(nday) :: &
         idxday               ! Index array for daytime points
    type(ty_gas_optics_rrtmgp),intent(in) :: &
         sw_gas_props         ! DDT containing SW spectral information
    real(kind_phys), dimension(size(Grid%xlon,1), Model%levs+1), intent(in) :: &
         p_lev                ! Pressure @ model layer-interfaces    (hPa)
    real(kind_phys), dimension(size(Grid%xlon,1), Model%levs+1), intent(in) :: &
         fluxswUP_allsky,   & ! SW All-sky flux                    (W/m2)
         fluxswDOWN_allsky, & ! SW All-sky flux                    (W/m2)
         fluxswUP_clrsky,   & ! SW Clear-sky flux                  (W/m2)
         fluxswDOWN_clrsky    ! SW All-sky flux                    (W/m2)
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
         cldtausw             ! approx .55mu band layer cloud optical depth  
    real(kind_phys),dimension(size(Grid%xlon,1), Model%levs) :: & 
         hswc                 ! All-sky heating rates (K/s)

    ! Outputs (mandatory)
    character(len=*), intent(out) :: &
         errmsg
    integer, intent(out) :: &
         errflg

    ! Outputs (optional)
    real(kind_phys), dimension(size(Grid%xlon,1), Model%levs), optional, intent(inout) :: &
         hsw0             ! Shortwave clear-sky heating-rate         (K/sec)
    type(profsw_type), dimension(size(Grid%xlon,1), Model%levs+1), intent(inout), optional :: &
         flxprf_sw        ! 2D radiative fluxes, components:
                          ! upfxc - total sky upward flux            (W/m2)
                          ! dnfxc - total sky dnward flux            (W/m2)
                          ! upfx0 - clear sky upward flux            (W/m2)
                          ! dnfx0 - clear sky dnward flux            (W/m2)
    type(cmpfsw_type), dimension(size(Grid%xlon,1)), intent(inout), optional :: &
         scmpsw           ! 2D surface fluxes, components:
                          ! uvbfc - total sky downward uv-b flux at  (W/m2)
                          ! uvbf0 - clear sky downward uv-b flux at  (W/m2)
                          ! nirbm - downward nir direct beam flux    (W/m2)
                          ! nirdf - downward nir diffused flux       (W/m2)
                          ! visbm - downward uv+vis direct beam flux (W/m2)
                          ! visdf - downward uv+vis diffused flux    (W/m2)
    ! Local variables
    integer :: i, j, k, iSFC, iTOA, itop, ibtc
    real(kind_phys) :: tem0d, tem1, tem2
    real(kind_phys), dimension(nDay, Model%levs) :: thetaTendClrSky, thetaTendAllSky
    logical :: l_clrskysw_hr, l_fluxessw2d, top_at_1, l_sfcFluxessw1D

   ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    ! Are any optional outputs requested?
    l_clrskysw_hr   = present(hsw0)
    l_fluxessw2d    = present(flxprf_sw)
    l_sfcfluxessw1D = present(scmpsw)

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
    ! Compute SW heating-rates
    ! #######################################################################################
    ! Initialize
    hswc = 0
    Diag%topfsw    = topfsw_type ( 0., 0., 0. )
    Radtend%sfcfsw = sfcfsw_type ( 0., 0., 0., 0. )
    if (l_clrskysw_hr) then
       hsw0(:,:) = 0.
    endif
    if (l_fluxessw2D) then
       flxprf_sw = profsw_type ( 0., 0., 0., 0. )
    endif
    if (l_sfcfluxessw1D) then
       scmpsw = cmpfsw_type (0.,0.,0.,0.,0.,0.)
    endif

    if (Model%lsswr .and. nDay .gt. 0) then
       ! Clear-sky heating-rate (optional)
       if (l_clrskysw_HR) then
          call check_error_msg('GFS_rrtmgp_post',compute_heating_rate( &
               fluxswUP_clrsky(idxday,:),                &
               fluxswDOWN_clrsky(idxday,:),                &
               p_lev(idxday,:),     &
               thetaTendClrSky))
          hsw0(idxday,:)=thetaTendClrSky
       endif
       ! All-sky heating-rate (mandatory)
       call check_error_msg('GFS_rrtmgp_post',compute_heating_rate(    &
            fluxswUP_allsky(idxday,:),                   &
            fluxswDOWN_allsky(idxday,:),                   &
            p_lev(idxday,:),        &
            thetaTendAllSky))
       hswc(idxday,:) = thetaTendAllSky

       ! Copy fluxes from RRTGMP types into model radiation types.
       ! Mandatory outputs
       Diag%topfsw(:)%upfxc    = fluxswUP_allsky(:,iTOA)
       Diag%topfsw(:)%upfx0    = fluxswUP_clrsky(:,iTOA)
       Diag%topfsw(:)%dnfxc    = fluxswDOWN_allsky(:,iTOA)
       Radtend%sfcfsw(:)%upfxc = fluxswUP_allsky(:,iSFC)
       Radtend%sfcfsw(:)%upfx0 = fluxswUP_clrsky(:,iSFC)
       Radtend%sfcfsw(:)%dnfxc = fluxswDOWN_allsky(:,iSFC)
       Radtend%sfcfsw(:)%dnfx0 = fluxswDOWN_clrsky(:,iSFC)

       ! Optional output
       if(l_fluxessw2D) then
          flxprf_sw(:,:)%upfxc = fluxswUP_allsky(:,:)
          flxprf_sw(:,:)%dnfxc = fluxswDOWN_allsky(:,:)
          flxprf_sw(:,:)%upfx0 = fluxswUP_clrsky(:,:)
          flxprf_sw(:,:)%dnfx0 = fluxswDOWN_clrsky(:,:)
       endif
    endif
    ! #######################################################################################
    ! Save SW outputs
    ! #######################################################################################
    if (Model%lsswr) then
       if (nday > 0) then
          ! All-sky heating rate
          do k = 1, Model%levs
             Radtend%htrsw(1:im,k) = hswc(1:im,k)
          enddo
          ! Clear-sky heating rate
          if (Model%swhtr) then
             do k = 1, Model%levs
                Radtend%swhc(1:im,k) = hsw0(1:im,k)
             enddo
          endif
          
          ! Surface down and up spectral component fluxes
          ! - Save two spectral bands' surface downward and upward fluxes for output.
          do i=1,im
             Coupling%nirbmdi(i) = scmpsw(i)%nirbm
             Coupling%nirdfdi(i) = scmpsw(i)%nirdf
             Coupling%visbmdi(i) = scmpsw(i)%visbm
             Coupling%visdfdi(i) = scmpsw(i)%visdf
             
             Coupling%nirbmui(i) = scmpsw(i)%nirbm * Interstitial%sfc_alb_nir_dir(1,i)
             Coupling%nirdfui(i) = scmpsw(i)%nirdf * Interstitial%sfc_alb_nir_dif(1,i)
             Coupling%visbmui(i) = scmpsw(i)%visbm * Interstitial%sfc_alb_uvvis_dir(1,i)
             Coupling%visdfui(i) = scmpsw(i)%visdf * Interstitial%sfc_alb_uvvis_dif(1,i)
          enddo
       else                   ! if_nday_block
          Radtend%htrsw(:,:) = 0.0
          Radtend%sfcfsw     = sfcfsw_type( 0.0, 0.0, 0.0, 0.0 )
          Diag%topfsw        = topfsw_type( 0.0, 0.0, 0.0 )
          scmpsw             = cmpfsw_type( 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 )
          
          do i=1,im
             Coupling%nirbmdi(i) = 0.0
             Coupling%nirdfdi(i) = 0.0
             Coupling%visbmdi(i) = 0.0
             Coupling%visdfdi(i) = 0.0
             
             Coupling%nirbmui(i) = 0.0
             Coupling%nirdfui(i) = 0.0
             Coupling%visbmui(i) = 0.0
             Coupling%visdfui(i) = 0.0
          enddo
          
          if (Model%swhtr) then
             Radtend%swhc(:,:) = 0
          endif
       endif                  ! end_if_nday
       
       ! Radiation fluxes for other physics processes
       do i=1,im
          Coupling%sfcnsw(i) = Radtend%sfcfsw(i)%dnfxc - Radtend%sfcfsw(i)%upfxc
          Coupling%sfcdsw(i) = Radtend%sfcfsw(i)%dnfxc
       enddo  
    endif                                ! end_if_lsswr

    ! #######################################################################################
    ! Save SW diagnostics
    ! - For time averaged output quantities (including total-sky and clear-sky SW and LW 
    !   fluxes at TOA and surface; conventional 3-domain cloud amount, cloud top and base 
    !   pressure, and cloud top temperature; aerosols AOD, etc.), store computed results in
    !   corresponding slots of array fluxr with appropriate time weights.
    ! - Collect the fluxr data for wrtsfc
    ! #######################################################################################
    if (Model%lssav) then
       if (Model%lsswr) then
          do i=1,im
             Diag%fluxr(i,34) = Diag%fluxr(i,34) + Model%fhswr*aerodp(i,1)  ! total aod at 550nm
             Diag%fluxr(i,35) = Diag%fluxr(i,35) + Model%fhswr*aerodp(i,2)  ! DU aod at 550nm
             Diag%fluxr(i,36) = Diag%fluxr(i,36) + Model%fhswr*aerodp(i,3)  ! BC aod at 550nm
             Diag%fluxr(i,37) = Diag%fluxr(i,37) + Model%fhswr*aerodp(i,4)  ! OC aod at 550nm
             Diag%fluxr(i,38) = Diag%fluxr(i,38) + Model%fhswr*aerodp(i,5)  ! SU aod at 550nm
             Diag%fluxr(i,39) = Diag%fluxr(i,39) + Model%fhswr*aerodp(i,6)  ! SS aod at 550nm
             if (Radtend%coszen(i) > 0.) then
                ! SW all-sky fluxes
                tem0d = Model%fhswr * Radtend%coszdg(i) / Radtend%coszen(i)
                Diag%fluxr(i,2 ) = Diag%fluxr(i,2)  +    Diag%topfsw(i)%upfxc * tem0d  ! total sky top sw up
                Diag%fluxr(i,3 ) = Diag%fluxr(i,3)  + Radtend%sfcfsw(i)%upfxc * tem0d  
                Diag%fluxr(i,4 ) = Diag%fluxr(i,4)  + Radtend%sfcfsw(i)%dnfxc * tem0d  ! total sky sfc sw dn
                ! SW uv-b fluxes
                Diag%fluxr(i,21) = Diag%fluxr(i,21) + scmpsw(i)%uvbfc * tem0d          ! total sky uv-b sw dn
                Diag%fluxr(i,22) = Diag%fluxr(i,22) + scmpsw(i)%uvbf0 * tem0d          ! clear sky uv-b sw dn
                ! SW TOA incoming fluxes
                Diag%fluxr(i,23) = Diag%fluxr(i,23) + Diag%topfsw(i)%dnfxc * tem0d     ! top sw dn 
                ! SW SFC flux components
                Diag%fluxr(i,24) = Diag%fluxr(i,24) + scmpsw(i)%visbm * tem0d          ! uv/vis beam sw dn
                Diag%fluxr(i,25) = Diag%fluxr(i,25) + scmpsw(i)%visdf * tem0d          ! uv/vis diff sw dn
                Diag%fluxr(i,26) = Diag%fluxr(i,26) + scmpsw(i)%nirbm * tem0d          ! nir beam sw dn
                Diag%fluxr(i,27) = Diag%fluxr(i,27) + scmpsw(i)%nirdf * tem0d          ! nir diff sw dn
                ! SW clear-sky fluxes
                Diag%fluxr(i,29) = Diag%fluxr(i,29) + Diag%topfsw(i)%upfx0 * tem0d
                Diag%fluxr(i,31) = Diag%fluxr(i,31) + Radtend%sfcfsw(i)%upfx0 * tem0d 
                Diag%fluxr(i,32) = Diag%fluxr(i,32) + Radtend%sfcfsw(i)%dnfx0 * tem0d

             endif
          enddo

          ! Save total and boundary-layer clouds
          do i=1,im
             Diag%fluxr(i,17) = Diag%fluxr(i,17) + raddt * cldsa(i,4)
             Diag%fluxr(i,18) = Diag%fluxr(i,18) + raddt * cldsa(i,5)
          enddo

          ! Save cld frac,toplyr,botlyr and top temp, note that the order of h,m,l cloud 
          ! is reversed for the fluxr output. save interface pressure (pa) of top/bot
          do j = 1, 3
             do i = 1, im
                tem0d = raddt * cldsa(i,j)
                itop  = mtopa(i,j)
                ibtc  = mbota(i,j)
                Diag%fluxr(i, 8-j) = Diag%fluxr(i, 8-j) + tem0d
                Diag%fluxr(i,11-j) = Diag%fluxr(i,11-j) + tem0d * Statein%prsi(i,itop)
                Diag%fluxr(i,14-j) = Diag%fluxr(i,14-j) + tem0d * Statein%prsi(i,ibtc)
                Diag%fluxr(i,17-j) = Diag%fluxr(i,17-j) + tem0d * Statein%tgrs(i,itop)
                
                ! Add optical depth and emissivity output
                tem1 = 0.
                do k=ibtc,itop
                   tem1 = tem1 + cldtausw(i,k)      ! approx .55 mu channel
                enddo
                Diag%fluxr(i,43-j) = Diag%fluxr(i,43-j) + tem0d * tem1
             enddo
          enddo
       endif
    endif

 
  end subroutine GFS_rrtmgp_sw_post_run

  ! #########################################################################################
  ! SUBROUTINE GFS_rrtmgp_sw_post_finalize
  ! #########################################################################################
  subroutine GFS_rrtmgp_sw_post_finalize ()
  end subroutine GFS_rrtmgp_sw_post_finalize

end module GFS_rrtmgp_sw_post

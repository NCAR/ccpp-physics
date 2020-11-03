module GFS_rrtmgp_lw_post 
  use machine,                   only: kind_phys
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
  ! ########################################################################################
!> \section arg_table_GFS_rrtmgp_lw_post_run
!! \htmlinclude GFS_rrtmgp_lw_post.html
!!
  subroutine GFS_rrtmgp_lw_post_run (nCol, nLev, lslwr, do_lw_clrsky_hr, save_diag, fhlwr, &
       p_lev, t_lay, tsfa, fluxlwUP_allsky, fluxlwDOWN_allsky, fluxlwUP_clrsky,            &
       fluxlwDOWN_clrsky, raddt, aerodp, cldsa, mtopa, mbota, cld_frac, cldtaulw, fluxr,   &
       sfcdlw, sfcflw, tsflw, htrlw, topflw, flxprf_lw, htrlwc, errmsg, errflg)

    ! Inputs                    
    integer, intent(in) :: &
         nCol,              & ! Horizontal loop extent 
         nLev                 ! Number of vertical layers
    logical, intent(in) :: & 
         lslwr,             & ! Logical flags for lw radiation calls
         do_lw_clrsky_hr,   & ! Output clear-sky SW heating-rate?         
    	 save_diag            ! Output radiation diagnostics?
    real(kind_phys), intent(in) :: &
         fhlwr                ! Frequency for SW radiation    	 
    real(kind_phys), dimension(nCol), intent(in) ::  &
         tsfa                 ! Lowest model layer air temperature for radiation (K)
    real(kind_phys), dimension(nCol, nLev), intent(in) :: &
         t_lay                ! Temperature @ model layer centers (K)
    real(kind_phys), dimension(nCol, nLev+1), intent(in) :: &
         p_lev,             & ! Pressure @ model layer-interfaces    (Pa)
         fluxlwUP_allsky,   & ! RRTMGP longwave all-sky flux (W/m2)
         fluxlwDOWN_allsky, & ! RRTMGP longwave all-sky flux (W/m2)
         fluxlwUP_clrsky,   & ! RRTMGP longwave clear-sky flux (W/m2)
         fluxlwDOWN_clrsky    ! RRTMGP longwave clear-sky flux (W/m2)
    real(kind_phys), intent(in) :: &
         raddt                ! Radiation time step
    real(kind_phys), dimension(nCol,NSPC1), intent(in) :: &
         aerodp               ! Vertical integrated optical depth for various aerosol species  
    real(kind_phys), dimension(nCol,5), intent(in) :: &
         cldsa                ! Fraction of clouds for low, middle, high, total and BL 
    integer,         dimension(nCol,3), intent(in) ::&
         mbota,             & ! vertical indices for low, middle and high cloud tops 
         mtopa                ! vertical indices for low, middle and high cloud bases
    real(kind_phys), dimension(nCol,nLev), intent(in) :: &
         cld_frac,          & ! Total cloud fraction in each layer
         cldtaulw             ! approx 10.mu band layer cloud optical depth  
    
    real(kind=kind_phys), dimension(:,:), intent(inout) :: fluxr
    
    ! Outputs (mandatory)
	real(kind_phys), dimension(nCol), intent(out) :: &
	     sfcdlw,            & ! Total sky sfc downward lw flux (W/m2)
	     tsflw                ! surface air temp during lw calculation (K)
	type(sfcflw_type), dimension(nCol), intent(out) :: &
	     sfcflw               ! LW radiation fluxes at sfc    
    real(kind_phys), dimension(nCol,nLev), intent(out) :: &
         htrlw                ! LW all-sky heating rate
    type(topflw_type), dimension(nCol), intent(out) :: &
         topflw               ! lw_fluxes_top_atmosphere
    character(len=*), intent(out) :: &
         errmsg
    integer, intent(out) :: &
         errflg
                  
    ! Outputs (optional)
    type(proflw_type), dimension(nCol, nLev+1), optional, intent(out) :: &
         flxprf_lw            ! 2D radiative fluxes, components:
                              ! upfxc - total sky upward flux (W/m2)
                              ! dnfxc - total sky dnward flux (W/m2)
                              ! upfx0 - clear sky upward flux (W/m2)
                              ! dnfx0 - clear sky dnward flux (W/m2)
    real(kind_phys),dimension(nCol, nLev),intent(out),optional  :: &
         htrlwc               ! Longwave clear-sky heating-rate (K/sec)
         
    ! Local variables
    integer :: i, j, k, iSFC, iTOA, itop, ibtc
    logical :: l_fluxeslw2d, top_at_1
    real(kind_phys) :: tem0d, tem1, tem2
    real(kind_phys),dimension(nCol,nLev) :: hlwc

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    if (.not. lslwr) return

    ! Are any optional outputs requested?
    l_fluxeslw2d    = present(flxprf_lw)

    ! #######################################################################################
    ! What is vertical ordering?
    ! #######################################################################################
    top_at_1 = (p_lev(1,1) .lt. p_lev(1, nLev))
    if (top_at_1) then 
       iSFC = nLev+1
       iTOA = 1
    else
       iSFC = 1
       iTOA = nLev+1
    endif

    ! #######################################################################################
    ! Compute LW heating-rates. 
    ! #######################################################################################
    ! Clear-sky heating-rate (optional)
    if (do_lw_clrsky_hr) then
        call check_error_msg('GFS_rrtmgp_post',compute_heating_rate(  &
            fluxlwUP_clrsky,   & ! IN  - RRTMGP upward longwave clear-sky flux profiles (W/m2)
            fluxlwDOWN_clrsky, & ! IN  - RRTMGP downward longwave clear-sky flux profiles (W/m2)
            p_lev,             & ! IN  - Pressure @ layer-interfaces (Pa)
            htrlwc))               ! OUT - Longwave clear-sky heating rate (K/sec)
    endif
    
    ! All-sky heating-rate (mandatory)
    call check_error_msg('GFS_rrtmgp_post',compute_heating_rate(     &
        fluxlwUP_allsky,      & ! IN  - RRTMGP upward longwave all-sky flux profiles (W/m2)
        fluxlwDOWN_allsky,    & ! IN  - RRTMGP downward longwave all-sky flux profiles (W/m2)
        p_lev,                & ! IN  - Pressure @ layer-interfaces (Pa)
        htrlw))                 ! OUT - Longwave all-sky heating rate (K/sec)

    ! #######################################################################################
    ! Save LW outputs.
    ! #######################################################################################
    ! Copy fluxes from RRTGMP types into model radiation types.
    ! Mandatory outputs
    topflw(:)%upfxc = fluxlwUP_allsky(:,iTOA)
    topflw(:)%upfx0 = fluxlwUP_clrsky(:,iTOA)
    sfcflw(:)%upfxc = fluxlwUP_allsky(:,iSFC)
    sfcflw(:)%upfx0 = fluxlwUP_clrsky(:,iSFC)
    sfcflw(:)%dnfxc = fluxlwDOWN_allsky(:,iSFC)
    sfcflw(:)%dnfx0 = fluxlwDOWN_clrsky(:,iSFC)
       
    ! Optional outputs
    if(l_fluxeslw2d) then
        flxprf_lw%upfxc = fluxlwUP_allsky
        flxprf_lw%dnfxc = fluxlwDOWN_allsky
        flxprf_lw%upfx0 = fluxlwUP_clrsky
        flxprf_lw%dnfx0 = fluxlwDOWN_clrsky
    endif
    
    ! Save surface air temp for diurnal adjustment at model t-steps
    tsflw (:) = tsfa(:)

    ! Radiation fluxes for other physics processes
    sfcdlw(:) = sfcflw(:)%dnfxc

    ! #######################################################################################
    ! Save LW diagnostics
    ! - For time averaged output quantities (including total-sky and clear-sky SW and LW 
    !   fluxes at TOA and surface; conventional 3-domain cloud amount, cloud top and base 
    !   pressure, and cloud top temperature; aerosols AOD, etc.), store computed results in
    !   corresponding slots of array fluxr with appropriate time weights.
    ! - Collect the fluxr data for wrtsfc
    ! #######################################################################################
    if (save_diag) then
          do i=1,nCol
             ! LW all-sky fluxes
             fluxr(i,1 ) = fluxr(i,1 ) + fhlwr * fluxlwUP_allsky(  i,iTOA)   ! total sky top lw up
             fluxr(i,19) = fluxr(i,19) + fhlwr * fluxlwDOWN_allsky(i,iSFC)   ! total sky sfc lw dn
             fluxr(i,20) = fluxr(i,20) + fhlwr * fluxlwUP_allsky(  i,iSFC)   ! total sky sfc lw up
             ! LW clear-sky fluxes
             fluxr(i,28) = fluxr(i,28) + fhlwr * fluxlwUP_clrsky(  i,iTOA)   ! clear sky top lw up
             fluxr(i,30) = fluxr(i,30) + fhlwr * fluxlwDOWN_clrsky(i,iSFC)   ! clear sky sfc lw dn
             fluxr(i,33) = fluxr(i,33) + fhlwr * fluxlwUP_clrsky(  i,iSFC)   ! clear sky sfc lw up
          enddo
          
          do i=1,nCol
             fluxr(i,17) = fluxr(i,17) + raddt * cldsa(i,4)
             fluxr(i,18) = fluxr(i,18) + raddt * cldsa(i,5)
          enddo

          ! Save cld frac,toplyr,botlyr and top temp, note that the order of h,m,l cloud is reversed for 
          ! the fluxr output. save interface pressure (pa) of top/bot
          do j = 1, 3
             do i = 1, nCol
                tem0d = raddt * cldsa(i,j)
                itop  = mtopa(i,j)
                ibtc  = mbota(i,j)
                fluxr(i, 8-j) = fluxr(i, 8-j) + tem0d
                fluxr(i,11-j) = fluxr(i,11-j) + tem0d * p_lev(i,itop)
                fluxr(i,14-j) = fluxr(i,14-j) + tem0d * p_lev(i,ibtc)
                fluxr(i,17-j) = fluxr(i,17-j) + tem0d * t_lay(i,itop)
                
                ! Add optical depth and emissivity output
                tem2 = 0.
                do k=ibtc,itop
                   tem2 = tem2 + cldtaulw(i,k)      ! approx 10. mu channel
                enddo
                fluxr(i,46-j) = fluxr(i,46-j) + tem0d * (1.0-exp(-tem2))
             enddo
          enddo
    endif

  end subroutine GFS_rrtmgp_lw_post_run

  ! #########################################################################################
  ! SUBROUTINE GFS_rrtmgp_lw_post_finalize
  ! #########################################################################################
  subroutine GFS_rrtmgp_lw_post_finalize ()
  end subroutine GFS_rrtmgp_lw_post_finalize

end module GFS_rrtmgp_lw_post

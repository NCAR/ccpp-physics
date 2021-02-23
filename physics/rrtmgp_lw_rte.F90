! ###########################################################################################
! ###########################################################################################
module rrtmgp_lw_rte
  use machine,                only: kind_phys
  use mo_rte_kind,            only: wl
  use mo_gas_optics_rrtmgp,   only: ty_gas_optics_rrtmgp
  use mo_cloud_optics,        only: ty_cloud_optics
  use mo_optical_props,       only: ty_optical_props_1scl, ty_optical_props_2str
  use mo_rte_lw,              only: rte_lw
  use mo_fluxes_byband,       only: ty_fluxes_byband
  use mo_source_functions,    only: ty_source_func_lw
  use rrtmgp_aux,             only: check_error_msg
  use rrtmgp_lw_gas_optics,   only: lw_gas_props
  implicit none

  public rrtmgp_lw_rte_init, rrtmgp_lw_rte_run, rrtmgp_lw_rte_finalize
contains

  ! #########################################################################################
  ! SUBROUTINE rrtmgp_lw_rte_init
  ! #########################################################################################
  subroutine rrtmgp_lw_rte_init()
  end subroutine rrtmgp_lw_rte_init

  ! #########################################################################################
  ! SUBROUTINE rrtmgp_lw_rte_run
  ! #########################################################################################
!! \section arg_table_rrtmgp_lw_rte_run
!! \htmlinclude rrtmgp_lw_rte_run.html
!!
  subroutine rrtmgp_lw_rte_run(doLWrad, doLWclrsky, use_LW_jacobian, doGP_lwscat, nCol,    &
       nLev, p_lev, sfc_emiss_byband, sources, lw_optical_props_clrsky,                    &
       lw_optical_props_clouds, lw_optical_props_aerosol, nGauss_angles, fluxlwUP_allsky,  &
       fluxlwDOWN_allsky, fluxlwUP_clrsky, fluxlwDOWN_clrsky, sfculw_jac, errmsg, errflg)

    ! Inputs
    logical, intent(in) :: &
         doLWrad,                 & ! Logical flag for longwave radiation call
         doLWclrsky,              & ! Compute clear-sky fluxes for clear-sky heating-rate?
         use_LW_jacobian,         & ! Compute Jacobian of LW to update radiative fluxes between radiation calls?
         doGP_lwscat                ! Include scattering in LW cloud-optics?
    integer, intent(in) :: &
         nCol,                    & ! Number of horizontal gridpoints
         nLev,                    & ! Number of vertical levels
         nGauss_angles              ! Number of angles used in Gaussian quadrature
    real(kind_phys), dimension(ncol,nLev+1), intent(in) :: &
         p_lev                      ! Pressure @ model layer-interfaces (hPa)
    real(kind_phys), dimension(lw_gas_props%get_nband(),ncol), intent(in) :: &
         sfc_emiss_byband           ! Surface emissivity in each band
    type(ty_source_func_lw),intent(in) :: &
         sources                    ! RRTMGP DDT: longwave source functions
    type(ty_optical_props_1scl),intent(inout) :: &
         lw_optical_props_aerosol, &! RRTMGP DDT: longwave aerosol radiative properties
         lw_optical_props_clrsky    ! RRTMGP DDT: longwave clear-sky radiative properties 
    type(ty_optical_props_2str),intent(inout) :: &
         lw_optical_props_clouds    ! RRTMGP DDT: longwave cloud radiative properties          
         
    ! Outputs
    real(kind_phys), dimension(ncol), intent(inout) :: &
         sfculw_jac                  ! Jacobian of upwelling LW surface radiation (W/m2/K)
    real(kind_phys), dimension(ncol,nLev+1), intent(inout) :: &
         fluxlwUP_allsky,          & ! All-sky flux (W/m2)
         fluxlwDOWN_allsky,        & ! All-sky flux (W/m2)
         fluxlwUP_clrsky,          & ! Clear-sky flux (W/m2)
         fluxlwDOWN_clrsky           ! All-sky flux (W/m2)
    character(len=*), intent(out) :: & 
         errmsg                      ! CCPP error message
    integer, intent(out) :: & 
         errflg                      ! CCPP error flag         

    ! Local variables
    type(ty_fluxes_byband) :: &
         flux_allsky, flux_clrsky
    real(kind_phys), dimension(ncol,nLev+1,lw_gas_props%get_nband()),target :: &
         fluxLW_up_allsky, fluxLW_up_clrsky, fluxLW_dn_allsky, fluxLW_dn_clrsky
    real(kind_phys), dimension(nCol,nLev+1) :: fluxlwUP_jac,fluxlwDOWN_jac
    logical :: &
         top_at_1
    integer :: iSFC, iTOA

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    if (.not. doLWrad) return

    ! Vertical ordering?
    top_at_1 = (p_lev(1,1) .lt. p_lev(1, nLev))
    if (top_at_1) then
       iSFC = nLev+1
       iTOA = 1
    else
       iSFC = 1
       iTOA = nLev+1
    endif

    ! Initialize RRTMGP DDT containing 2D(3D) fluxes
    flux_allsky%bnd_flux_up => fluxLW_up_allsky
    flux_allsky%bnd_flux_dn => fluxLW_dn_allsky
    flux_clrsky%bnd_flux_up => fluxLW_up_clrsky
    flux_clrsky%bnd_flux_dn => fluxLW_dn_clrsky

    !
    ! Compute clear-sky fluxes (if requested)
    !
    ! Add aerosol optics to gas optics
    call check_error_msg('rrtmgp_lw_rte_run',lw_optical_props_aerosol%increment(lw_optical_props_clrsky))
    call lw_optical_props_aerosol%finalize()

    ! Call RTE solver
    if (doLWclrsky) then
       call check_error_msg('rrtmgp_lw_rte_run',rte_lw(           &
            lw_optical_props_clrsky,         & ! IN  - optical-properties
            top_at_1,                        & ! IN  - veritcal ordering flag
            sources,                         & ! IN  - source function
            sfc_emiss_byband,                & ! IN  - surface emissivity in each LW band
            flux_clrsky,                     & ! OUT - Fluxes
            n_gauss_angles = nGauss_angles))   ! IN  - Number of angles in Gaussian quadrature

       ! Store fluxes
       fluxlwUP_clrsky   = sum(flux_clrsky%bnd_flux_up,dim=3)
       fluxlwDOWN_clrsky = sum(flux_clrsky%bnd_flux_dn,dim=3)
    else
       fluxlwUP_clrsky   = 0.0
       fluxlwDOWN_clrsky = 0.0   
    endif
    
    !
    ! All-sky fluxes
    !

    ! Include LW cloud-scattering?
    if (doGP_lwscat) then 
       ! Add clear-sky optics to cloud-optics (2-stream)
       call check_error_msg('rrtmgp_lw_rte_run',lw_optical_props_clrsky%increment(lw_optical_props_clouds))
       call lw_optical_props_clrsky%finalize() 
       
       if (use_LW_jacobian) then
          ! Compute LW Jacobians
          call check_error_msg('rrtmgp_lw_rte_run',rte_lw(           &
               lw_optical_props_clouds,         & ! IN  - optical-properties
               top_at_1,                        & ! IN  - veritcal ordering flag
               sources,                         & ! IN  - source function
               sfc_emiss_byband,                & ! IN  - surface emissivity in each LW band
               flux_allsky,                     & ! OUT - Flxues 
               n_gauss_angles = nGauss_angles,  & ! IN  - Number of angles in Gaussian quadrature
               flux_up_Jac    = fluxlwUP_jac,   & ! OUT - surface temperature flux (upward) Jacobian (W/m2/K)
               flux_dn_Jac    = fluxlwDOWN_jac))  ! OUT - surface temperature flux (downward) Jacobian (W/m2/K)
          sfculw_jac = fluxlwUP_jac(:,iSFC)
       else
          call check_error_msg('rrtmgp_lw_rte_run',rte_lw(           &
               lw_optical_props_clouds,         & ! IN  - optical-properties
               top_at_1,                        & ! IN  - veritcal ordering flag
               sources,                         & ! IN  - source function
               sfc_emiss_byband,                & ! IN  - surface emissivity in each LW band
               flux_allsky,                     & ! OUT - Flxues 
               n_gauss_angles = nGauss_angles))   ! IN  - Number of angles in Gaussian quadrature    
       end if
    ! No scattering in LW clouds.   
    else
       ! Add cloud optics to clear-sky optics (scalar)
       call check_error_msg('rrtmgp_lw_rte_run',lw_optical_props_clouds%increment(lw_optical_props_clrsky))
       call lw_optical_props_clouds%finalize() 
    
       if (use_LW_jacobian) then
          ! Compute LW Jacobians
          call check_error_msg('rrtmgp_lw_rte_run',rte_lw(           &
               lw_optical_props_clrsky,         & ! IN  - optical-properties
               top_at_1,                        & ! IN  - veritcal ordering flag
               sources,                         & ! IN  - source function
               sfc_emiss_byband,                & ! IN  - surface emissivity in each LW band
               flux_allsky,                     & ! OUT - Flxues 
               n_gauss_angles = nGauss_angles,  & ! IN  - Number of angles in Gaussian quadrature
               flux_up_Jac    = fluxlwUP_jac,   & ! OUT - surface temperature flux (upward) Jacobian (W/m2/K)
               flux_dn_Jac    = fluxlwDOWN_jac))  ! OUT - surface temperature flux (downward) Jacobian (W/m2/K)
          sfculw_jac = fluxlwUP_jac(:,iSFC)
       else
          call check_error_msg('rrtmgp_lw_rte_run',rte_lw(           &
               lw_optical_props_clrsky,         & ! IN  - optical-properties
               top_at_1,                        & ! IN  - veritcal ordering flag
               sources,                         & ! IN  - source function
               sfc_emiss_byband,                & ! IN  - surface emissivity in each LW band
               flux_allsky,                     & ! OUT - Flxues 
               n_gauss_angles = nGauss_angles))   ! IN  - Number of angles in Gaussian quadrature    
       end if    
    endif
    
    ! Store fluxes
    fluxlwUP_allsky   = sum(flux_allsky%bnd_flux_up,dim=3)
    fluxlwDOWN_allsky = sum(flux_allsky%bnd_flux_dn,dim=3) 

  end subroutine rrtmgp_lw_rte_run
  
  ! #########################################################################################
  ! SUBROUTINE rrtmgp_lw_rte_finalize
  ! #########################################################################################
  subroutine rrtmgp_lw_rte_finalize()
  end subroutine rrtmgp_lw_rte_finalize


end module rrtmgp_lw_rte

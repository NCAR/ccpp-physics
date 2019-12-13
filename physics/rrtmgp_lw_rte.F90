! ###########################################################################################
! ###########################################################################################
module rrtmgp_lw_rte
  use machine,                only: kind_phys
  use mo_rte_kind,            only: wl
  use mo_gas_optics_rrtmgp,   only: ty_gas_optics_rrtmgp
  use mo_cloud_optics,        only: ty_cloud_optics
  use mo_optical_props,       only: ty_optical_props_1scl
  use mo_rte_lw,              only: rte_lw
  use mo_fluxes_byband,       only: ty_fluxes_byband
  use mo_source_functions,    only: ty_source_func_lw
  use rrtmgp_aux,             only: check_error_msg

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
!! \htmlinclude rrtmgp_lw_rte.html
!!
  subroutine rrtmgp_lw_rte_run(doLWrad, nCol, nLev, p_lay, t_lay, p_lev, skt, lw_gas_props, &
       sfc_emiss_byband, sources, lw_optical_props_clrsky, lw_optical_props_clouds,         &
       lw_optical_props_aerosol, fluxlwUP_allsky, fluxlwDOWN_allsky, fluxlwUP_clrsky,       &
       fluxlwDOWN_clrsky, hlw0, hlwb, errmsg, errflg)

    ! Inputs
    logical, intent(in) :: &
         doLWrad                    ! Logical flag for longwave radiation call
    integer, intent(in) :: &
         nCol,                    & ! Number of horizontal gridpoints
         nLev                       ! Number of vertical levels
    real(kind_phys), dimension(ncol,nLev), intent(in) :: &
         p_lay,                   & ! Pressure @ model layer-centers (hPa)
         t_lay                      ! Temperature (K)
    real(kind_phys), dimension(ncol,nLev+1), intent(in) :: &
         p_lev                      ! Pressure @ model layer-interfaces (hPa)
    real(kind_phys), dimension(ncol), intent(in) :: &
         skt                        ! Surface(skin) temperature (K)
    type(ty_gas_optics_rrtmgp),intent(in) :: &
         lw_gas_props               ! DDT containing LW spectral information
    real(kind_phys), dimension(lw_gas_props%get_nband(),ncol), intent(in) :: &
         sfc_emiss_byband           ! Surface emissivity in each band
    type(ty_source_func_lw),intent(in) :: &
         sources
    type(ty_optical_props_1scl),intent(inout) :: &
         lw_optical_props_clrsky    ! RRTMGP DDT: longwave clear-sky radiative properties 
    type(ty_optical_props_1scl),intent(in) :: &
         lw_optical_props_clouds, & ! RRTMGP DDT: longwave cloud radiative properties 
         lw_optical_props_aerosol   ! RRTMGP DDT: longwave aerosol radiative properties

    ! Outputs
    real(kind_phys), dimension(ncol,nLev+1), intent(out) :: &
         fluxlwUP_allsky,          & ! All-sky flux (W/m2)
         fluxlwDOWN_allsky,        & ! All-sky flux (W/m2)
         fluxlwUP_clrsky,          & ! Clear-sky flux (W/m2)
         fluxlwDOWN_clrsky           ! All-sky flux (W/m2)
    character(len=*), intent(out) :: & 
         errmsg                      ! CCPP error message
    integer, intent(out) :: & 
         errflg                      ! CCPP error flag

    ! Outputs (optional)
    real(kind_phys), dimension(ncol,nLev,lw_gas_props%get_nband()), optional, intent(inout) :: &
         hlwb                        ! All-sky heating rate, by band (K/sec)
    real(kind_phys), dimension(ncol,nLev), optional, intent(inout) :: &
         hlw0                        ! Clear-sky heating rate (K/sec)

    ! Local variables
    type(ty_fluxes_byband) :: &
         flux_allsky, flux_clrsky
    real(kind_phys), dimension(ncol,nLev+1,lw_gas_props%get_nband()),target :: &
         fluxLW_up_allsky, fluxLW_up_clrsky, fluxLW_dn_allsky, fluxLW_dn_clrsky
    logical :: &
         l_ClrSky_HR, l_AllSky_HR_byband, top_at_1

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

    if (.not. doLWrad) return

    ! Vertical ordering?
    top_at_1 = (p_lev(1,1) .lt. p_lev(1, nLev))

    ! Are any optional outputs requested? Need to know now to compute correct fluxes.
    l_ClrSky_HR        = present(hlw0)
    l_AllSky_HR_byband = present(hlwb)

    ! Initialize RRTMGP DDT containing 2D(3D) fluxes
    flux_allsky%bnd_flux_up => fluxLW_up_allsky
    flux_allsky%bnd_flux_dn => fluxLW_dn_allsky
    flux_clrsky%bnd_flux_up => fluxLW_up_clrsky
    flux_clrsky%bnd_flux_dn => fluxLW_dn_clrsky

    ! Compute clear-sky fluxes (if requested)
    ! Clear-sky fluxes are gas+aerosol
    call check_error_msg('rrtmgp_lw_rte_run',lw_optical_props_aerosol%increment(lw_optical_props_clrsky))
    if (l_ClrSky_HR) then
       call check_error_msg('rrtmgp_lw_rte_run',rte_lw(           &
            lw_optical_props_clrsky, & ! IN  - optical-properties
            top_at_1,                & ! IN  - veritcal ordering flag
            sources,                 & ! IN  - source function
            sfc_emiss_byband,        & ! IN  - surface emissivity in each LW band
            flux_clrsky))
       ! Store fluxes
       fluxlwUP_clrsky   = sum(flux_clrsky%bnd_flux_up,dim=3)
       fluxlwDOWN_clrsky = sum(flux_clrsky%bnd_flux_dn,dim=3)
    endif

    ! All-sky fluxes
    ! Clear-sky fluxes are (gas+aerosol)+clouds
    call check_error_msg('rrtmgp_lw_rte_run',lw_optical_props_clouds%increment(lw_optical_props_clrsky))
    call check_error_msg('rrtmgp_lw_rte_run',rte_lw(           &
         lw_optical_props_clrsky, & ! IN  - optical-properties
         top_at_1,                & ! IN  - veritcal ordering flag
         sources,                 & ! IN  - source function
         sfc_emiss_byband,        & ! IN  - surface emissivity in each LW band
         flux_allsky))
    ! Store fluxes
    fluxlwUP_allsky   = sum(flux_allsky%bnd_flux_up,dim=3)
    fluxlwDOWN_allsky = sum(flux_allsky%bnd_flux_dn,dim=3) 

    ! Only output fluxes by-band when heating-rate profiles by band are requested.
    !if (l_AllSky_HR_byband) then
    !endif

  end subroutine rrtmgp_lw_rte_run
  
  ! #########################################################################################
  ! SUBROUTINE rrtmgp_lw_rte_finalize
  ! #########################################################################################
  subroutine rrtmgp_lw_rte_finalize()
  end subroutine rrtmgp_lw_rte_finalize


end module rrtmgp_lw_rte

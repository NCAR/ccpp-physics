! ###########################################################################################
! ###########################################################################################
module rrtmgp_lw_clrallsky_driver
  use machine,                only: kind_phys
  use GFS_typedefs,           only: GFS_control_type, GFS_radtend_type
  use mo_rte_kind,            only: wl
  use mo_gas_optics_rrtmgp,   only: ty_gas_optics_rrtmgp
  use mo_cloud_optics,        only: ty_cloud_optics
  use mo_optical_props,       only: ty_optical_props_1scl
  use mo_rrtmgp_clr_all_sky,  only: rte_lw
  use mo_gas_concentrations,  only: ty_gas_concs
  use mo_fluxes_byband,       only: ty_fluxes_byband
  use rrtmgp_aux,             only: check_error_msg

  public rrtmgp_lw_clrallsky_driver_init, rrtmgp_lw_clrallsky_driver_run, rrtmgp_lw_clrallsky_driver_finalize
contains

  ! #########################################################################################
  ! SUBROUTINE rrtmgp_lw_clrallsky_driver_init
  ! #########################################################################################
  subroutine rrtmgp_lw_clrallsky_driver_init()
  end subroutine rrtmgp_lw_clrallsky_driver_init

  ! #########################################################################################
  ! SUBROUTINE rrtmgp_lw_clrallsky_driver_run
  ! #########################################################################################
!! \section arg_table_rrtmgp_lw_clrallsky_driver_run
!! \htmlinclude rrtmgp_lw_clrallsky_driver.html
!!
  subroutine rrtmgp_lw_clrallsky_driver_run(Model, Radtend, ncol, lw_gas_props, p_lay, t_lay,&
       p_lev, skt, gas_concentrations, lw_optical_props_clouds, lw_optical_props_aerosol,    &
       lslwr, fluxlwUP_allsky, fluxlwDOWN_allsky, fluxlwUP_clrsky, fluxlwDOWN_clrsky, hlw0,  &
       hlwb, errmsg, errflg)

    ! Inputs
    type(GFS_control_type),   intent(in) :: &
         Model
    type(GFS_radtend_type), intent(in) :: &
         Radtend                 ! Fortran DDT containing FV3-GFS radiation tendencies 
    integer, intent(in) :: &
         ncol                    ! Number of horizontal gridpoints
    real(kind_phys), dimension(ncol,model%levs), intent(in) :: &
         p_lay,                & ! Pressure @ model layer-centers         (hPa)
         t_lay                   ! Temperature                            (K)
    real(kind_phys), dimension(ncol,model%levs+1), intent(in) :: &
         p_lev                   ! Pressure @ model layer-interfaces      (hPa)
    real(kind_phys), dimension(ncol), intent(in) :: &
         skt                     ! Surface(skin) temperature              (K)
    type(ty_gas_optics_rrtmgp),intent(in) :: &
         lw_gas_props                ! DDT containing LW spectral information
    type(ty_optical_props_1scl),intent(in) :: &
         lw_optical_props_clouds, & ! RRTMGP DDT: longwave cloud radiative properties 
         lw_optical_props_aerosol ! RRTMGP DDT: longwave aerosol radiative properties
    type(ty_gas_concs),intent(in) :: &
         gas_concentrations      ! RRTMGP DDT: trace gas concentrations   (vmr)
    logical, intent(in) :: &
         lslwr                   ! Flag to calculate LW irradiances
 
    ! Outputs
    character(len=*), intent(out) :: errmsg
    integer, intent(out) :: errflg
    real(kind_phys), dimension(ncol,model%levs), intent(out) :: &
         fluxlwUP_allsky,   & ! All-sky flux                    (W/m2)
         fluxlwDOWN_allsky, & ! All-sky flux                    (W/m2)
         fluxlwUP_clrsky,   & ! Clear-sky flux                  (W/m2)
         fluxlwDOWN_clrsky    ! All-sky flux                    (W/m2)

    ! Outputs (optional)
    real(kind_phys), dimension(ncol,model%levs,lw_gas_props%get_nband()), optional, intent(inout) :: &
         hlwb             ! All-sky heating rate, by band     (K/sec)
    real(kind_phys), dimension(ncol,model%levs), optional, intent(inout) :: &
         hlw0             ! Clear-sky heating rate            (K/sec)

    ! Local variables
    type(ty_fluxes_byband) :: &
         flux_allsky, & ! All-sky flux                        (W/m2)
         flux_clrsky    ! Clear-sky flux                      (W/m2)
    real(kind_phys), dimension(ncol,model%levs+1),target :: &
         fluxLW_up_allsky, fluxLW_up_clrsky, fluxLW_dn_allsky, fluxLW_dn_clrsky
    real(kind_phys), dimension(ncol,model%levs+1,lw_gas_props%get_nband()),target :: &
         fluxLWBB_up_allsky, fluxLWBB_dn_allsky
    logical :: l_ClrSky_HR, l_AllSky_HR_byband

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0
    if (.not. lslwr) return
 
    ! Are any optional outputs requested? Need to know now to compute correct fluxes.
    l_ClrSky_HR        = present(hlw0)
    l_AllSky_HR_byband = present(hlwb)

    ! Initialize RRTMGP DDT containing 2D(3D) fluxes
    flux_allsky%flux_up => fluxLW_up_allsky
    flux_allsky%flux_dn => fluxLW_dn_allsky
    flux_clrsky%flux_up => fluxLW_up_clrsky
    flux_clrsky%flux_dn => fluxLW_dn_clrsky
    ! Only calculate fluxes by-band, only when heating-rate profiles by band are requested.
    if (l_AllSky_HR_byband) then
       flux_allsky%bnd_flux_up => fluxLWBB_up_allsky
       flux_allsky%bnd_flux_dn => fluxLWBB_dn_allsky
    endif

    ! Call RRTMGP LW scheme
    call check_error_msg('rrtmgp_lw_clrallsky_driver_run',rte_lw(           &
         lw_gas_props,                       & ! IN  - spectral information 
         gas_concentrations,                 & ! IN  - gas concentrations (vmr)
         p_lay,                              & ! IN  - pressure at layer interfaces (Pa)
         t_lay,                              & ! IN  - temperature at layer interfaes (K)
         p_lev,                              & ! IN  - pressure at layer centers (Pa)
         skt,                                & ! IN  - skin temperature (K)
         Radtend%sfc_emiss_byband,           & ! IN  - surface emissivity in each LW band
         lw_optical_props_clouds,            & ! IN  - DDT containing cloud optical information 
         flux_allsky,                        & ! OUT - Fluxes, all-sky, 3D (nCol,model%levs,nBand) 
         flux_clrsky,                        & ! OUT - Fluxes, clear-sky, 3D (nCol,model%levs,nBand) 
         aer_props = lw_optical_props_aerosol)) ! IN(optional) - DDT containing aerosol optical information
    fluxlwUP_allsky   = flux_allsky%flux_up
    fluxlwDOWN_allsky = flux_allsky%flux_dn 
    fluxlwUP_clrsky   = flux_clrsky%flux_up
    fluxlwDOWN_clrsky = flux_clrsky%flux_dn

  end subroutine rrtmgp_lw_clrallsky_driver_run
  
  ! #########################################################################################
  ! SUBROUTINE rrtmgp_lw_clrallsky_driver_finalize
  ! #########################################################################################
  subroutine rrtmgp_lw_clrallsky_driver_finalize()
  end subroutine rrtmgp_lw_clrallsky_driver_finalize


end module rrtmgp_lw_clrallsky_driver

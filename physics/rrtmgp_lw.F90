! ###########################################################################################
! ###########################################################################################
module rrtmgp_lw
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

  public rrtmgp_lw_init, rrtmgp_lw_run, rrtmgp_lw_finalize
contains

  ! #########################################################################################
  ! SUBROUTINE rrtmgp_lw_init
  ! #########################################################################################
  subroutine rrtmgp_lw_init()
  end subroutine rrtmgp_lw_init

  ! #########################################################################################
  ! SUBROUTINE rrtmgp_lw_run
  ! #########################################################################################
!! \section arg_table_rrtmgp_lw_run Argument Table
!! | local_name              | standard_name                                                                                 | long_name                                                          | units | rank | type                  |    kind   | intent | optional |
!! |-------------------------|-----------------------------------------------------------------------------------------------|--------------------------------------------------------------------|-------|------|-----------------------|-----------|--------|----------|
!! | Model                   | GFS_control_type_instance                                                                     | Fortran DDT containing FV3-GFS model control parameters            | DDT   |    0 | GFS_control_type      |           | in     | F        |
!! | Radtend                 | GFS_radtend_type_instance                                                                     | Fortran DDT containing FV3-GFS radiation tendencies                | DDT   |    0 | GFS_radtend_type      |           | in     | F        |
!! | ncol                    | horizontal_loop_extent                                                                        | horizontal dimension                                               | count |    0 | integer               |           | in     | F        |
!! | p_lay                   | air_pressure_at_layer_for_RRTMGP_in_hPa                                                       | air pressure layer                                                 | hPa   |    2 | real                  | kind_phys | in     | F        |
!! | p_lev                   | air_pressure_at_interface_for_RRTMGP_in_hPa                                                   | air pressure level                                                 | hPa   |    2 | real                  | kind_phys | in     | F        |
!! | t_lay                   | air_temperature_at_layer_for_RRTMGP                                                           | air temperature layer                                              | K     |    2 | real                  | kind_phys | in     | F        |
!! | skt                     | surface_ground_temperature_for_radiation                                                      | surface ground temperature for radiation                           | K     |    1 | real                  | kind_phys | in     | F        |
!! | lw_gas_props            | coefficients_for_lw_gas_optics                                                                | DDT containing spectral information for RRTMGP LW radiation scheme | DDT   |    0 | ty_gas_optics_rrtmgp  |           | in     | F        |
!! | optical_propsLW_clds    | longwave_optical_properties_for_cloudy_atmosphere                                             | Fortran DDT containing RRTMGP optical properties                   | DDT   |    0 | ty_optical_props_1scl |           | in     | F        |
!! | optical_propsLW_aerosol | longwave_optical_properties_for_aerosols                                                      | Fortran DDT containing RRTMGP optical properties                   | DDT   |    0 | ty_optical_props_1scl |           | in     | F        |
!! | gas_concentrations      | Gas_concentrations_for_RRTMGP_suite                                                           | DDT containing gas concentrations for RRTMGP radiation scheme      | DDT   |    0 | ty_gas_concs          |           | in     | F        |
!! | lslwr                   | flag_to_calc_lw                                                                               | flag to calculate LW irradiances                                   | flag  |    0 | logical               |           | in     | F        |
!! | hlw0                    | tendency_of_air_temperature_due_to_longwave_heating_assuming_clear_sky_on_radiation_time_step | longwave clear sky heating rate                                    | K s-1 |    2 | real                  | kind_phys | in     | T        |
!! | hlwb                    | lw_heating_rate_spectral                                                                      | longwave total sky heating rate (spectral)                         | K s-1 |    3 | real                  | kind_phys | in     | T        |
!! | fluxUP_allsky           | lw_flux_profile_upward_allsky                                                                 | RRTMGP upward longwave all-sky flux profile                        | W m-2 |    2 | real                  | kind_phys | out    | F        |
!! | fluxDOWN_allsky         | lw_flux_profile_downward_allsky                                                               | RRTMGP downward longwave all-sky flux profile                      | W m-2 |    2 | real                  | kind_phys | out    | F        |
!! | fluxUP_clrsky           | lw_flux_profile_upward_clrsky                                                                 | RRTMGP upward longwave clr-sky flux profile                        | W m-2 |    2 | real                  | kind_phys | out    | F        |
!! | fluxDOWN_clrsky         | lw_flux_profile_downward_clrsky                                                               | RRTMGP downward longwave clr-sky flux profile                      | W m-2 |    2 | real                  | kind_phys | out    | F        |
!! | errmsg                  | ccpp_error_message                                                                            | error message for error handling in CCPP                           | none  |    0 | character             | len=*     | out    | F        |
!! | errflg                  | ccpp_error_flag                                                                               | error flag for error handling in CCPP                              | flag  |    0 | integer               |           | out    | F        |
!!
  subroutine rrtmgp_lw_run(Model, Radtend, ncol, lw_gas_props, p_lay, t_lay, p_lev, skt, &
       gas_concentrations, optical_propsLW_clds, optical_propsLW_aerosol,&
       lslwr, fluxUP_allsky, fluxDOWN_allsky, fluxUP_clrsky, fluxDOWN_clrsky, hlw0, hlwb, errmsg, errflg)

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
         optical_propsLW_clds, & ! RRTMGP DDT: longwave cloud radiative properties 
         optical_propsLW_aerosol ! RRTMGP DDT: longwave aerosol radiative properties
    type(ty_gas_concs),intent(in) :: &
         gas_concentrations      ! RRTMGP DDT: trace gas concentrations   (vmr)
    logical, intent(in) :: &
         lslwr                   ! Flag to calculate LW irradiances
 
    ! Outputs
    character(len=*), intent(out) :: errmsg
    integer, intent(out) :: errflg
    real(kind_phys), dimension(ncol,model%levs), intent(out) :: &
         fluxUP_allsky,   & ! All-sky flux                    (W/m2)
         fluxDOWN_allsky, & ! All-sky flux                    (W/m2)
         fluxUP_clrsky,   & ! Clear-sky flux                  (W/m2)
         fluxDOWN_clrsky    ! All-sky flux                    (W/m2)

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
    call check_error_msg('rrtmgp_lw_run',rte_lw(           &
         lw_gas_props,                       & ! IN  - spectral information 
         gas_concentrations,                 & ! IN  - gas concentrations (vmr)
         p_lay,                              & ! IN  - pressure at layer interfaces (Pa)
         t_lay,                              & ! IN  - temperature at layer interfaes (K)
         p_lev,                              & ! IN  - pressure at layer centers (Pa)
         skt,                                & ! IN  - skin temperature (K)
         Radtend%sfc_emiss_byband,           & ! IN  - surface emissivity in each LW band
         optical_propsLW_clds,               & ! IN  - DDT containing cloud optical information 
         flux_allsky,                        & ! OUT - Fluxes, all-sky, 3D (nCol,model%levs,nBand) 
         flux_clrsky,                        & ! OUT - Fluxes, clear-sky, 3D (nCol,model%levs,nBand) 
         aer_props = optical_propsLW_aerosol)) ! IN(optional) - DDT containing aerosol optical information
    fluxUP_allsky   = flux_allsky%flux_up
    fluxDOWN_allsky = flux_allsky%flux_dn 
    fluxUP_clrsky   = flux_clrsky%flux_up
    fluxDOWN_clrsky = flux_clrsky%flux_dn

  end subroutine rrtmgp_lw_run
  
  ! #########################################################################################
  ! SUBROUTINE rrtmgp_lw_finalize
  ! #########################################################################################
  subroutine rrtmgp_lw_finalize()
  end subroutine rrtmgp_lw_finalize


end module rrtmgp_lw

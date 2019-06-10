! ###########################################################################################
! ###########################################################################################
module rrtmgp_lw
  use machine,                only: kind_phys
  use GFS_typedefs,           only: GFS_control_type, GFS_radtend_type, GFS_statein_type
  use mo_rte_kind,            only: wl
  use mo_gas_optics_rrtmgp,   only: ty_gas_optics_rrtmgp
  use mo_cloud_optics,        only: ty_cloud_optics
  use mo_optical_props,       only: ty_optical_props_1scl
  use mo_rte_lw,              only: rte_lw
  use mo_fluxes_byband,       only: ty_fluxes_byband
  use mo_source_functions,    only: ty_source_func_lw
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
!! | local_name            | standard_name                                                                                 | long_name                                                           | units | rank | type                  |    kind   | intent | optional |
!! |-----------------------|-----------------------------------------------------------------------------------------------|---------------------------------------------------------------------|-------|------|-----------------------|-----------|--------|----------|
!! | Model                 | GFS_control_type_instance                                                                     | Fortran DDT containing FV3-GFS model control parameters             | DDT   |    0 | GFS_control_type      |           | in     | F        |
!! | Radtend               | GFS_radtend_type_instance                                                                     | Fortran DDT containing FV3-GFS radiation tendencies                 | DDT   |    0 | GFS_radtend_type      |           | in     | F        |
!! | Statein               | GFS_statein_type_instance                                                                     | Fortran DDT containing FV3-GFS prognostic state data in from dycore | DDT   |    0 | GFS_statein_type      |           | in     | F        |
!! | ncol                  | horizontal_loop_extent                                                                        | horizontal dimension                                                | count |    0 | integer               |           | in     | F        |
!! | p_lay                 | air_pressure_at_layer_for_RRTMGP_in_hPa                                                       | air pressure layer                                                  | hPa   |    2 | real                  | kind_phys | in     | F        |
!! | p_lev                 | air_pressure_at_interface_for_RRTMGP_in_hPa                                                   | air pressure level                                                  | hPa   |    2 | real                  | kind_phys | in     | F        |
!! | t_lay                 | air_temperature_at_layer_for_RRTMGP                                                           | air temperature layer                                               | K     |    2 | real                  | kind_phys | in     | F        |
!! | skt                   | surface_ground_temperature_for_radiation                                                      | surface ground temperature for radiation                            | K     |    1 | real                  | kind_phys | in     | F        |
!! | lw_gas_props          | coefficients_for_lw_gas_optics                                                                | DDT containing spectral information for RRTMGP LW radiation scheme  | DDT   |    0 | ty_gas_optics_rrtmgp  |           | in     | F        |
!! | optical_props_clrsky  | longwave_optical_properties_for_clear_sky                                                     | Fortran DDT containing RRTMGP optical properties                    | DDT   |    0 | ty_optical_props_1scl |           | inout  | F        |
!! | optical_props_cloud   | longwave_optical_properties_for_cloudy_atmosphere                                             | Fortran DDT containing RRTMGP optical properties                    | DDT   |    0 | ty_optical_props_1scl |           | in     | F        |
!! | optical_props_aerosol | longwave_optical_properties_for_aerosols                                                      | Fortran DDT containing RRTMGP optical properties                    | DDT   |    0 | ty_optical_props_1scl |           | in     | F        |
!! | sources               | longwave_source_function                                                                      | Fortran DDT containing RRTMGP source functions                      | DDT   |    0 | ty_source_func_lw     |           | in     | F        |
!! | lslwr                 | flag_to_calc_lw                                                                               | flag to calculate LW irradiances                                    | flag  |    0 | logical               |           | in     | F        |
!! | hlw0                  | tendency_of_air_temperature_due_to_longwave_heating_assuming_clear_sky_on_radiation_time_step | longwave clear sky heating rate                                     | K s-1 |    2 | real                  | kind_phys | in     | T        |
!! | hlwb                  | lw_heating_rate_spectral                                                                      | longwave total sky heating rate (spectral)                          | K s-1 |    3 | real                  | kind_phys | in     | T        |
!! | fluxUP_allsky         | lw_flux_profile_upward_allsky                                                                 | RRTMGP upward longwave all-sky flux profile                         | W m-2 |    2 | real                  | kind_phys | out    | F        |
!! | fluxDOWN_allsky       | lw_flux_profile_downward_allsky                                                               | RRTMGP downward longwave all-sky flux profile                       | W m-2 |    2 | real                  | kind_phys | out    | F        |
!! | fluxUP_clrsky         | lw_flux_profile_upward_clrsky                                                                 | RRTMGP upward longwave clr-sky flux profile                         | W m-2 |    2 | real                  | kind_phys | out    | F        |
!! | fluxDOWN_clrsky       | lw_flux_profile_downward_clrsky                                                               | RRTMGP downward longwave clr-sky flux profile                       | W m-2 |    2 | real                  | kind_phys | out    | F        |
!! | errmsg                | ccpp_error_message                                                                            | error message for error handling in CCPP                            | none  |    0 | character             | len=*     | out    | F        |
!! | errflg                | ccpp_error_flag                                                                               | error flag for error handling in CCPP                               | flag  |    0 | integer               |           | out    | F        |
!!
  subroutine rrtmgp_lw_run(Model, Statein, Radtend, ncol, lw_gas_props, p_lay, t_lay, p_lev, &
       skt, sources, optical_props_clrsky, optical_props_cloud, optical_props_aerosol, lslwr,&
       fluxUP_allsky, fluxDOWN_allsky, fluxUP_clrsky, fluxDOWN_clrsky, hlw0, hlwb, errmsg, errflg)

    ! Inputs
    type(GFS_control_type), intent(in) :: &
         Model                   ! Fortran DDT containing FV3-GFS model control parameters 
    type(GFS_radtend_type), intent(in) :: &
         Radtend                 ! Fortran DDT containing FV3-GFS radiation tendencies 
    type(GFS_statein_type), intent(in) :: &
         Statein                 ! Fortran DDT containing FV3-GFS prognostic state data in from dycore 
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
         lw_gas_props            ! DDT containing LW spectral information
    type(ty_optical_props_1scl),intent(inout) :: &
         optical_props_clrsky    ! RRTMGP DDT: longwave clear-sky radiative properties 
    type(ty_optical_props_1scl),intent(in) :: &
         optical_props_cloud,  & ! RRTMGP DDT: longwave cloud radiative properties 
         optical_props_aerosol   ! RRTMGP DDT: longwave aerosol radiative properties
    type(ty_source_func_lw),intent(in) :: &
         sources
    logical, intent(in) :: &
         lslwr                   ! Flag to calculate LW irradiances
 
    ! Outputs
    character(len=*), intent(out) :: & 
         errmsg                  ! CCPP error message
    integer, intent(out) :: & 
         errflg                  ! CCPP error flag
    real(kind_phys), dimension(ncol,model%levs), intent(out) :: &
         fluxUP_allsky,        & ! All-sky flux                    (W/m2)
         fluxDOWN_allsky,      & ! All-sky flux                    (W/m2)
         fluxUP_clrsky,        & ! Clear-sky flux                  (W/m2)
         fluxDOWN_clrsky         ! All-sky flux                    (W/m2)

    ! Outputs (optional)
    real(kind_phys), dimension(ncol,model%levs,lw_gas_props%get_nband()), optional, intent(inout) :: &
         hlwb                    ! All-sky heating rate, by band   (K/sec)
    real(kind_phys), dimension(ncol,model%levs), optional, intent(inout) :: &
         hlw0                    ! Clear-sky heating rate          (K/sec)

    ! Local variables
    type(ty_fluxes_byband) :: &
         flux_allsky, flux_clrsky
    real(kind_phys), dimension(ncol,model%levs+1),target :: &
         fluxLW_up_allsky, fluxLW_up_clrsky, fluxLW_dn_allsky, fluxLW_dn_clrsky
    real(kind_phys), dimension(ncol,model%levs+1,lw_gas_props%get_nband()),target :: &
         fluxLWBB_up_allsky, fluxLWBB_dn_allsky
    logical :: &
         l_ClrSky_HR, l_AllSky_HR_byband, top_at_1

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0
    if (.not. lslwr) return

    ! Vertical ordering?
    top_at_1 = (Statein%prsi(1,1) .lt.  Statein%prsi(1, Model%levs))

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

    ! Compute clear-sky fluxes (if requested)
    ! Clear-sky fluxes are gas+aerosol
    call check_error_msg('rrtmgp_lw_run',optical_props_aerosol%increment(optical_props_clrsky))
    if (l_ClrSky_HR) then
       call check_error_msg('rrtmgp_lw_run',rte_lw(           &
            optical_props_clrsky,               & ! IN  - optical-properties
            top_at_1,                           & ! IN  - veritcal ordering flag
            sources,                            & ! IN  - source function
            Radtend%sfc_emiss_byband,           & ! IN  - surface emissivity in each LW band
            flux_clrsky))
       ! Store fluxes
       fluxUP_clrsky   = flux_clrsky%flux_up
       fluxDOWN_clrsky = flux_clrsky%flux_dn
    endif

    ! All-sky fluxes
    ! Clear-sky fluxes are (gas+aerosol)+clouds
    call check_error_msg('rrtmgp_lw_run',optical_props_cloud%increment(optical_props_clrsky))
    call check_error_msg('rrtmgp_lw_run',rte_lw(           &
         optical_props_clrsky,               & ! IN  - optical-properties
         top_at_1,                           & ! IN  - veritcal ordering flag
         sources,                            & ! IN  - source function
         Radtend%sfc_emiss_byband,           & ! IN  - surface emissivity in each LW band
         flux_allsky))
    ! Store fluxes
    fluxUP_allsky   = flux_allsky%flux_up
    fluxDOWN_allsky = flux_allsky%flux_dn 

  end subroutine rrtmgp_lw_run
  
  ! #########################################################################################
  ! SUBROUTINE rrtmgp_lw_finalize
  ! #########################################################################################
  subroutine rrtmgp_lw_finalize()
  end subroutine rrtmgp_lw_finalize


end module rrtmgp_lw

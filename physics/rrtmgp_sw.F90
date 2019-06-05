! ###########################################################################################
! ###########################################################################################
module rrtmgp_sw
  use machine,                 only: kind_phys
  use GFS_typedefs,            only: GFS_control_type, GFS_radtend_type
  use mo_rte_kind,             only: wl
  use mo_gas_optics_rrtmgp,    only: ty_gas_optics_rrtmgp
  use mo_cloud_optics,         only: ty_cloud_optics
  use mo_optical_props,        only: ty_optical_props_2str
  use mo_rrtmgp_clr_all_sky,   only: rte_sw
  use mo_gas_concentrations,   only: ty_gas_concs
  use mo_fluxes_byband,        only: ty_fluxes_byband
  use module_radsw_parameters, only: cmpfsw_type
  use rrtmgp_sw_cloud_optics,  only: rrtmgp_sw_cloud_optics_init
  use rrtmgp_sw_gas_optics,    only: rrtmgp_sw_gas_optics_init, check_error_msg

  public rrtmgp_sw_init, rrtmgp_sw_run, rrtmgp_sw_finalize

contains

  ! #########################################################################################
  ! SUBROUTINE rrtmgp_sw_init
  ! #########################################################################################
  subroutine rrtmgp_sw_init()
  end subroutine rrtmgp_sw_init

  ! #########################################################################################
  ! SUBROUTINE rrtmgp_sw_run
  ! #########################################################################################
!! \section arg_table_rrtmgp_sw_run Argument Table
!! | local_name              | standard_name                                                                                  | long_name                                                                | units | rank | type                  |    kind   | intent | optional |
!! |-------------------------|------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------|-------|------|-----------------------|-----------|--------|----------|
!! | Model                   | GFS_control_type_instance                                                                      | Fortran DDT containing FV3-GFS model control parameters                  | DDT   |    0 | GFS_control_type      |           | in     | F        |
!! | Radtend                 | GFS_radtend_type_instance                                                                      | Fortran DDT containing FV3-GFS radiation tendencies                      | DDT   |    0 | GFS_radtend_type      |           | in     | F        |
!! | ncol                    | horizontal_loop_extent                                                                         | horizontal dimension                                                     | count |    0 | integer               |           | in     | F        |
!! | p_lay                   | air_pressure_at_layer_for_RRTMGP_in_hPa                                                        | air pressure layer                                                       | hPa   |    2 | real                  | kind_phys | in     | F        |
!! | p_lev                   | air_pressure_at_interface_for_RRTMGP_in_hPa                                                    | air pressure level                                                       | hPa   |    2 | real                  | kind_phys | in     | F        |
!! | t_lay                   | air_temperature_at_layer_for_RRTMGP                                                            | air temperature layer                                                    | K     |    2 | real                  | kind_phys | in     | F        |
!! | sw_gas_props            | coefficients_for_sw_gas_optics                                                                 | DDT containing spectral information for RRTMGP SW radiation scheme       | DDT   |    0 | ty_gas_optics_rrtmgp  |           | in     | F        |
!! | optical_props_clds      | shortwave_optical_properties_for_cloudy_atmosphere                                             | Fortran DDT containing RRTMGP optical properties                         | DDT   |    0 | ty_optical_props_2str |           | in     | F        |
!! | optical_props_aerosol   | shortwave_optical_properties_for_aerosols                                                      | Fortran DDT containing RRTMGP optical properties                         | DDT   |    0 | ty_optical_props_2str |           | in     | F        |
!! | gas_concentrations      | Gas_concentrations_for_RRTMGP_suite                                                            | DDT containing gas concentrations for RRTMGP radiation scheme            | DDT   |    0 | ty_gas_concs          |           | in     | F        |
!! | lsswr                   | flag_to_calc_sw                                                                                | flag to calculate SW irradiances                                         | flag  |    0 | logical               |           | in     | F        |
!! | sfcalb_nir_dir          | surface_shortwave_albedo_near_infrared_direct_in_each_band                                     | surface sw near-infrared direct albedo in each SW band                   | frac  |    2 | real                  | kind_phys | in     | F        |
!! | sfcalb_nir_dif          | surface_shortwave_albedo_near_infrared_diffuse_in_each_band                                    | surface sw near-infrared diffuse albedo in each SW band                  | frac  |    2 | real                  | kind_phys | in     | F        |
!! | cossza                  | cosine_of_zenith_angle                                                                         | cosine of the solar zenit angle                                          | none  |    1 | real                  | kind_phys | in     | F        |
!! | nday                    | daytime_points_dimension                                                                       | daytime points dimension                                                 | count |    0 | integer               |           | in     | F        |
!! | idxday                  | daytime_points                                                                                 | daytime points                                                           | index |    1 | integer               |           | in     | F        |
!! | hsw0                    | tendency_of_air_temperature_due_to_shortwave_heating_assuming_clear_sky_on_radiation_time_step | shortwave clear sky heating rate                                         | K s-1 |    2 | real                  | kind_phys | inout  | T        |
!! | hswb                    | sw_heating_rate_spectral                                                                       | shortwave total sky heating rate (spectral)                              | K s-1 |    3 | real                  | kind_phys | inout  | T        |
!! | scmpsw                  | components_of_surface_downward_shortwave_fluxes                                                | derived type for special components of surface downward shortwave fluxes | W m-2 |    1 | cmpfsw_type           |           | inout  | T        |
!! | fluxUP_allsky           | sw_flux_profile_upward_allsky                                                                  | RRTMGP upward shortwave all-sky flux profile                             | W m-2 |    2 | real                  | kind_phys | out    | F        |
!! | fluxDOWN_allsky         | sw_flux_profile_downward_allsky                                                                | RRTMGP downward shortwave all-sky flux profile                           | W m-2 |    2 | real                  | kind_phys | out    | F        |
!! | fluxUP_clrsky           | sw_flux_profile_upward_clrsky                                                                  | RRTMGP upward shortwave clr-sky flux profile                             | W m-2 |    2 | real                  | kind_phys | out    | F        |
!! | fluxDOWN_clrsky         | sw_flux_profile_downward_clrsky                                                                | RRTMGP downward shortwave clr-sky flux profile                           | W m-2 |    2 | real                  | kind_phys | out    | F        |
!! | errmsg                  | ccpp_error_message                                                                             | error message for error handling in CCPP                                 | none  |    0 | character             | len=*     | out    | F        |
!! | errflg                  | ccpp_error_flag                                                                                | error flag for error handling in CCPP                                    | flag  |    0 | integer               |           | out    | F        |
!!
  subroutine rrtmgp_sw_run(Model, Radtend, ncol, sw_gas_props, p_lay, t_lay, p_lev, gas_concentrations, &
       optical_props_clds, optical_props_aerosol, &
       lsswr, sfcalb_nir_dir, sfcalb_nir_dif, cossza,  nday, idxday, hsw0, hswb, scmpsw, &
       fluxUP_allsky, fluxDOWN_allsky, fluxUP_clrsky, fluxDOWN_clrsky, errmsg, errflg)

    ! Inputs
    type(GFS_control_type),   intent(in)    :: Model
    type(GFS_radtend_type),   intent(in)    :: Radtend
    integer, intent(in) :: &
         ncol,                 & ! Number of horizontal gridpoints
         nday                    ! Number of daytime points
    integer, intent(in), dimension(nday) :: &
         idxday                  ! Index array for daytime points
    real(kind_phys), dimension(ncol,Model%levs), intent(in) :: &
         p_lay,                & ! Pressure @ model layer-centers         (hPa)
         t_lay                   ! Temperature                            (K)
    real(kind_phys), dimension(ncol,Model%levs+1), intent(in) :: &
         p_lev                   ! Pressure @ model layer-interfaces      (hPa)
    type(ty_gas_optics_rrtmgp),intent(in) :: &
         sw_gas_props                ! DDT containing SW spectral information
    real(kind_phys), dimension(sw_gas_props%get_nband(),ncol), intent(in) :: &
         sfcalb_nir_dir,  & ! Surface albedo direct (near-IR)         (1)
         sfcalb_nir_dif     ! Surface albedo diffuse (near-IR)        (1)
    real(kind_phys), dimension(ncol), intent(in) :: &
         cossza             ! Cosine of solar zenith angle            (1)
    type(ty_optical_props_2str),intent(in) :: &
         optical_props_clds, & ! RRTMGP DDT: longwave cloud radiative properties 
         optical_props_aerosol ! RRTMGP DDT: longwave aerosol radiative properties

    type(ty_gas_concs),intent(in) :: &
         gas_concentrations      ! RRTMGP DDT: trace gas concentrations   (vmr)
    logical, intent(in) :: &
         lsswr                   ! Flag to calculate SW irradiances

    ! Outputs
    character(len=*), intent(out) :: errmsg
    integer, intent(out) :: errflg
    real(kind_phys), dimension(ncol,Model%levs), intent(out) :: &
         fluxUP_allsky,   & ! All-sky flux                    (W/m2)
         fluxDOWN_allsky, & ! All-sky flux                    (W/m2)
         fluxUP_clrsky,   & ! Clear-sky flux                  (W/m2)
         fluxDOWN_clrsky    ! All-sky flux                    (W/m2)

    ! Inputs (optional) (NOTE. We only need the optional arguments to know what fluxes to output, HR's are computed later)
    real(kind_phys), dimension(ncol,Model%levs), optional, intent(inout) :: &
         hsw0             ! Clear-sky heating rate            (K/sec)
    real(kind_phys), dimension(ncol,Model%levs,sw_gas_props%get_nband()), intent(inout), optional :: &
         hswb             ! All-sky heating rate, by band     (K/sec)
    ! Outputs (optional)
    type(cmpfsw_type), dimension(ncol), intent(inout),optional :: &
         scmpsw           ! 2D surface fluxes, components:
                          ! uvbfc - total sky downward uv-b flux at  (W/m2)
                          ! uvbf0 - clear sky downward uv-b flux at  (W/m2)
                          ! nirbm - downward nir direct beam flux    (W/m2)
                          ! nirdf - downward nir diffused flux       (W/m2)
                          ! visbm - downward uv+vis direct beam flux (W/m2)
                          ! visdf - downward uv+vis diffused flux    (W/m2)

    ! Local variables
    type(ty_fluxes_byband) :: &
         flux_allsky, & ! All-sky flux                      (W/m2)
         flux_clrsky    ! Clear-sky flux                    (W/m2)
    real(kind_phys), dimension(nday,Model%levs+1),target :: &
         fluxSW_up_allsky, fluxSW_up_clrsky, fluxSW_dn_allsky, fluxSW_dn_clrsky
    real(kind_phys), dimension(nday,Model%levs+1,sw_gas_props%get_nband()),target :: &
         fluxSWBB_up_allsky, fluxSWBB_dn_allsky
    real(kind_phys), dimension(ncol,Model%levs) :: vmrTemp
    logical :: l_ClrSky_HR=.false., l_AllSky_HR_byband=.false., l_scmpsw=.false.
    integer :: iGas
    type(ty_optical_props_2str)  :: &
         optical_props_clds_daylit, & ! RRTMGP DDT: longwave cloud radiative properties 
         optical_props_aerosol_daylit ! RRTMGP DDT: longwave aerosol radiative properties
    type(ty_gas_concs) :: &
         gas_concentrations_daylit    ! RRTMGP DDT: trace gas concentrations   (vmr)

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg  = 0

    if (.not. lsswr) return

    ! Are any optional outputs requested? Need to know now to compute correct fluxes.
    l_ClrSky_HR        = present(hsw0)
    l_AllSky_HR_byband = present(hswb)
    l_scmpsw           = present(scmpsw)
    if ( l_scmpsw ) then
       scmpsw = cmpfsw_type (0., 0., 0., 0., 0., 0.)
    endif
    fluxUP_allsky(:,:)   = 0._kind_phys
    fluxDOWN_allsky(:,:) = 0._kind_phys
    fluxUP_clrsky(:,:)   = 0._kind_phys
    fluxDOWN_clrsky(:,:) = 0._kind_phys

    if (nDay .gt. 0) then

       ! Subset the cloud and aerosol radiative properties over daylit points.
       ! Cloud optics [nDay,Model%levs,nBands]
       call check_error_msg('rrtmgp_sw_run',optical_props_clds_daylit%alloc_2str(nday, Model%levs, sw_gas_props))
       optical_props_clds_daylit%tau    = optical_props_clds%tau(idxday,:,:)
       optical_props_clds_daylit%ssa    = optical_props_clds%ssa(idxday,:,:)
       optical_props_clds_daylit%g      = optical_props_clds%g(idxday,:,:)
       ! Aerosol optics [nDay,Model%levs,nBands]
       call check_error_msg('rrtmgp_sw_run',optical_props_aerosol_daylit%alloc_2str(nday, Model%levs, sw_gas_props%get_band_lims_wavenumber()))
       optical_props_aerosol_daylit%tau = optical_props_aerosol%tau(idxday,:,:)
       optical_props_aerosol_daylit%ssa = optical_props_aerosol%ssa(idxday,:,:)
       optical_props_aerosol_daylit%g   = optical_props_aerosol%g(idxday,:,:)
      
       ! Similarly, subset the gas concentrations.
       do iGas=1,Model%nGases
          call check_error_msg('rrtmgp_sw_run',gas_concentrations%get_vmr(trim(Radtend%active_gases(iGas)),vmrTemp))
          call check_error_msg('rrtmgp_sw_run',gas_concentrations_daylit%set_vmr(trim(Radtend%active_gases(iGas)),vmrTemp(idxday,:)))
       enddo

       ! Initialize RRTMGP DDT containing 2D(3D) fluxes
       flux_allsky%flux_up => fluxSW_up_allsky
       flux_allsky%flux_dn => fluxSW_dn_allsky
       flux_clrsky%flux_up => fluxSW_up_clrsky
       flux_clrsky%flux_dn => fluxSW_dn_clrsky
       ! Only calculate fluxes by-band, only when heating-rate profiles by band are requested.
       if (l_AllSky_HR_byband) then
          flux_allsky%bnd_flux_up => fluxSWBB_up_allsky
          flux_allsky%bnd_flux_dn => fluxSWBB_dn_allsky
       endif
       
       ! Call RRTMGP SW scheme
       call check_error_msg('rrtmgp_sw_run',rte_sw(               &
            sw_gas_props,                             & ! IN  - spectral information 
            gas_concentrations_daylit,                & ! IN  - gas concentrations (vmr)
            p_lay(idxday,1:Model%levs),               & ! IN  - pressure at layer interfaces (Pa)  
            t_lay(idxday,1:Model%levs),               & ! IN  - temperature at layer interfaes (K)
            p_lev(idxday,1:Model%levs+1),             & ! IN  - pressure at layer centers (Pa)
            cossza(idxday),                           & ! IN  - Cosine of solar zenith angle
            sfcalb_nir_dir(:,idxday),                 & ! IN  - Shortwave surface albedo (direct)
            sfcalb_nir_dif(:,idxday),                 & ! IN  - Shortwave surface albedo (diffuse)
            optical_props_clds_daylit,                & ! IN  - DDT containing cloud optical information 
            flux_allsky,                              & ! OUT - Fluxes, all-sky, 3D (nCol,Model%levs,nBand) 
            flux_clrsky,                              & ! OUT - Fluxes, clear-sky, 3D (nCol,Model%levs,nBand) 
            aer_props = optical_props_aerosol_daylit))  ! IN(optional) - DDT containing aerosol optical information
       fluxUP_allsky(idxday,:)   = flux_allsky%flux_up
       fluxDOWN_allsky(idxday,:) = flux_allsky%flux_dn 
       fluxUP_clrsky(idxday,:)   = flux_clrsky%flux_up
       fluxDOWN_clrsky(idxday,:) = flux_clrsky%flux_dn
    endif
  end subroutine rrtmgp_sw_run
  
  ! #########################################################################################
  ! SUBROUTINE rrtmgp_sw_finalize
  ! #########################################################################################
  subroutine rrtmgp_sw_finalize()
  end subroutine rrtmgp_sw_finalize

end module rrtmgp_sw

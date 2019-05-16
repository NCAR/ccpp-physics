! ###########################################################################################
! ###########################################################################################
module rrtmgp_lw
  use machine,               only: kind_phys
  use mo_gas_optics_rrtmgp,  only: ty_gas_optics_rrtmgp
  use mo_cloud_optics,       only: ty_cloud_optics
  use mo_optical_props,      only: ty_optical_props_1scl
  use mo_rrtmgp_clr_all_sky, only: rte_lw
  use mo_gas_concentrations, only: ty_gas_concs
  use mo_fluxes_byband,      only: ty_fluxes_byband

  public rrtmgp_lw_init, rrtmgp_lw_run, rrtmgp_lw_finalize
contains

  subroutine rrtmgp_lw_init()
  end subroutine rrtmgp_lw_init
  
  ! #########################################################################################
  ! #########################################################################################
!! \section arg_table_rrtmgp_lw_run Argument Table
!! | local_name            | standard_name                                   | long_name                                                          | units | rank | type                  |    kind   | intent | optional |
!! |-----------------------|-------------------------------------------------|--------------------------------------------------------------------|-------|------|-----------------------|-----------|--------|----------|
!! | ncol                  | horizontal_loop_extent                          | horizontal dimension                                               | count |    0 | integer               |           | in     | F        |
!! | nlay                  | adjusted_vertical_layer_dimension_for_radiation | number of vertical layers for radiation                            | count |    0 | integer               |           | in     | F        |
!! | p_lay                 | air_pressure_at_layer_for_radiation_in_hPa      | air pressure layer                                                 | hPa   |    2 | real                  | kind_phys | in     | F        |
!! | p_lev                 | air_pressure_at_interface_for_radiation_in_hPa  | air pressure level                                                 | hPa   |    2 | real                  | kind_phys | in     | F        |
!! | t_lay                 | air_temperature_at_layer_for_radiation          | air temperature layer                                              | K     |    2 | real                  | kind_phys | in     | F        |
!! | skt                   | surface_ground_temperature_for_radiation        | surface ground temperature for radiation                           | K     |    1 | real                  | kind_phys | in     | F        |
!! | sfc_emiss             | surface_longwave_emissivity_in_each_band        | surface lw emissivity in fraction in each LW band                  | frac  |    2 | real                  | kind_phys | in     | F        |
!! | kdist_lw              | K_distribution_file_for_RRTMGP_LW_scheme        | DDT containing spectral information for RRTMGP LW radiation scheme | DDT   |    0 | ty_gas_optics_rrtmgp  |           | in     | F        |
!! | optical_props_clds    | optical_properties_for_cloudy_atmosphere        | Fortran DDT containing RRTMGP optical properties                   | DDT   |    0 | ty_optical_props_1scl |           | in     | F        |
!! | optical_props_aerosol | optical_properties_for_aerosols                 | Fortran DDT containing RRTMGP optical properties                   | DDT   |    0 | ty_optical_props_1scl |           | in     | F        |
!! | gas_concentrations    | Gas_concentrations_for_RRTMGP_suite             | DDT containing gas concentrations for RRTMGP radiation scheme      | DDT   |    0 | ty_gas_concs          |           | in     | F        |
!! | fluxLW_allsky         | lw_flux_profiles_byband_allsky                  | Fortran DDT containing RRTMGP 3D fluxes                            | DDT   |    0 | ty_fluxes_byband      |           | out    | F        |
!! | fluxLW_clrsky         | lw_flux_profiles_byband_clrsky                  | Fortran DDT containing RRTMGP 3D fluxes                            | DDT   |    0 | ty_fluxes_byband      |           | out    | F        |
!! | errmsg                | ccpp_error_message                              | error message for error handling in CCPP                           | none  |    0 | character             | len=*     | out    | F        |
!! | errflg                | ccpp_error_flag                                 | error flag for error handling in CCPP                              | flag  |    0 | integer               |           | out    | F        |
!!
  subroutine rrtmgp_lw_run(ncol, nlay, kdist_lw, p_lay, t_lay, p_lev, skt, &
       sfc_emiss, gas_concentrations, optical_props_clds, optical_props_aerosol,&
       fluxLW_allsky, fluxLW_clrsky, errmsg, errflg)

    ! Inputs
    integer, intent(in) :: &
         ncol,         & ! Number of horizontal gridpoints
         nlay            ! Number of vertical layers
    real(kind_phys), dimension(ncol,nlay), intent(in) :: &
         p_lay,        & ! Pressure @ model layer-centers         (hPa)
         t_lay           ! Temperature                            (K)
    real(kind_phys), dimension(ncol,nlay+1), intent(in) :: &
         p_lev           ! Pressure @ model layer-interfaces      (hPa)
    real(kind_phys), dimension(ncol), intent(in) :: &
         skt             ! Surface(skin) temperature              (K)
    type(ty_gas_optics_rrtmgp),intent(in) :: &
         kdist_lw        ! DDT containing LW spectral information
    real(kind_phys), dimension(kdist_lw%get_nband(),ncol) :: &
         sfc_emiss       ! Surface emissivity                     (1)
    type(ty_optical_props_1scl),intent(in) :: &
         optical_props_clds, &
         optical_props_aerosol
    type(ty_gas_concs),intent(in) :: &
         gas_concentrations
    type(ty_fluxes_byband),intent(out) :: &
         fluxLW_allsky, & ! All-sky flux                      (W/m2)
         fluxLW_clrsky    ! Clear-sky flux                    (W/m2)

    ! Outputs
    character(len=*), intent(out) :: errmsg
    integer, intent(out) :: errflg
    
    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0
    
    ! Call RRTMGP LW scheme
    call check_error_msg(rte_lw(           &
         kdist_lw,                         & ! IN  - spectral information 
         gas_concentrations,               & ! IN  - gas concentrations (vmr)
         p_lay,                            & ! IN  - pressure at layer interfaces (Pa)
         t_lay,                            & ! IN  - temperature at layer interfaes (K)
         p_lev,                            & ! IN  - pressure at layer centers (Pa)
         skt,                              & ! IN  - skin temperature (K)
         sfc_emiss,                        & ! IN  - surface emissivity in each LW band
         optical_props_clds,               & ! IN  - DDT containing cloud optical information 
         fluxLW_allsky,                    & ! OUT - Fluxes, all-sky, 3D (nCol,nLay,nBand) 
         fluxLW_clrsky,                    & ! OUT - Fluxes, clear-sky, 3D (nCol,nLay,nBand) 
         aer_props = optical_props_aerosol)) ! IN(optional) - DDT containing aerosol optical information

  end subroutine rrtmgp_lw_run
  
  subroutine rrtmgp_lw_finalize()
  end subroutine rrtmgp_lw_finalize
  subroutine check_error_msg(error_msg)
    character(len=*), intent(in) :: error_msg
    
    if(error_msg /= "") then
       print*,"ERROR(rrtmgp_sw_main.F90): "
       print*,trim(error_msg)
       return
    end if
  end subroutine check_error_msg    


end module rrtmgp_lw

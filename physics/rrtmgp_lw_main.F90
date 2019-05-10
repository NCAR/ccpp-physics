! ###########################################################################################
! ###########################################################################################
module rrtmgp_lw_main
  use machine,              only: kind_phys
  use mo_gas_optics_rrtmgp, only: ty_gas_optics_rrtmgp
  use mo_cloud_optics,      only: ty_cloud_optics
  use mo_optical_props,     only: ty_optical_props_1scl

  public rrtmgp_lw_main_init, rrtmgp_lw_main_run, rrtmgp_lw_main_finalize
contains

  subroutine rrtmgp_lw_main_init()
  end subroutine rrtmgp_lw_main_init
  
  ! #########################################################################################
  ! #########################################################################################
!! \section arg_table_rrtmgp_lw_main_run Argument Table
!! | local_name         | standard_name                                   | long_name                                                          | units | rank | type                  |    kind   | intent | optional |
!! |--------------------|-------------------------------------------------|--------------------------------------------------------------------|-------|------|-----------------------|-----------|--------|----------|
!! | ncol               | horizontal_loop_extent                          | horizontal dimension                                               | count |    0 | integer               |           | in     | F        |
!! | nlay               | adjusted_vertical_layer_dimension_for_radiation | number of vertical layers for radiation                            | count |    0 | integer               |           | in     | F        |
!! | p_lay              | air_pressure_at_layer_for_radiation_in_hPa      | air pressure layer                                                 | hPa   |    2 | real                  | kind_phys | in     | F        |
!! | p_lev              | air_pressure_at_interface_for_radiation_in_hPa  | air pressure level                                                 | hPa   |    2 | real                  | kind_phys | in     | F        |
!! | t_lay              | air_temperature_at_layer_for_radiation          | air temperature layer                                              | K     |    2 | real                  | kind_phys | in     | F        |
!! | skt                | surface_ground_temperature_for_radiation        | surface ground temperature for radiation                           | K     |    1 | real                  | kind_phys | in     | F        |
!! | sfc_emiss          | surface_longwave_emissivity_in_each_band        | surface lw emissivity in fraction in each LW band                  | frac  |    2 | real                  | kind_phys | in     | F        |
!! | kdist_lw           | K_distribution_file_for_RRTMGP_LW_scheme        | DDT containing spectral information for RRTMGP LW radiation scheme | DDT   |    0 | ty_gas_optics_rrtmgp  |           | in     | F        |
!! | optical_props_clds | optical_properties_for_cloudy_atmosphere        | Fortran DDT containing RRTMGP optical properties                   | DDT   |    0 | ty_optical_props_1scl |           | in     | F        |
!! | errmsg             | ccpp_error_message                              | error message for error handling in CCPP                           | none  |    0 | character             | len=*     | out    | F        |
!! | errflg             | ccpp_error_flag                                 | error flag for error handling in CCPP                              | flag  |    0 | integer               |           | out    | F        |
!!
  subroutine rrtmgp_lw_main_run(ncol, nlay, kdist_lw, p_lay, t_lay,     &
       p_lev, skt, sfc_emiss, optical_props_clds, errmsg, errflg)

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
         optical_props_clds
 
    ! Outputs
    character(len=*), intent(out) :: errmsg
    integer, intent(out) :: errflg

    ! Initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

  end subroutine rrtmgp_lw_main_run
  
  subroutine rrtmgp_lw_main_finalize()
  end subroutine rrtmgp_lw_main_finalize
  


end module rrtmgp_lw_main

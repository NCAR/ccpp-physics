!> \file GFS_surface_diag.F90
!!  Contains code related to the surface diagnostic scheme.

      module sfc_diag_post

      contains

      subroutine sfc_diag_post_init ()
      end subroutine sfc_diag_post_init

      subroutine sfc_diag_post_finalize()
      end subroutine sfc_diag_post_finalize
#if 0
!> \section arg_table_sfc_diag_post_run Argument Table
!! | local_name     | standard_name                                                                                                       | long_name                                                                           | units       | rank | type       |    kind   | intent | optional |
!! |----------------|---------------------------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------|-------------|------|------------|-----------|--------|----------|
!! | im             | horizontal_loop_extent                                                                                              | horizontal loop extent                                                              | count       |    0 | integer    |           | in     | F        |
!! | lssav          | flag_diagnostics                                                                                                    | logical flag for storing diagnostics                                                | flag        |    0 | logical    |           | in     | F        |
!! | dtf            | time_step_for_dynamics                                                                                              | dynamics timestep                                                                   | s           |    0 | real       | kind_phys | in     | F        |
!! | con_eps        | ratio_of_dry_air_to_water_vapor_gas_constants                                                                       | rd/rv                                                                               | none        |    0 | real       | kind_phys | in     | F        |
!! | con_epsm1      | ratio_of_dry_air_to_water_vapor_gas_constants_minus_one                                                             | (rd/rv) - 1                                                                         | none        |    0 | real       | kind_phys | in     | F        |
!! | pgr            | surface_air_pressure                                                                                                | surface pressure                                                                    | Pa          |    1 | real       | kind_phys | in     | F        |
!! | t2m            | temperature_at_2m                                                                                                   | 2 meter temperature                                                                 | K           |    1 | real       | kind_phys | in     | F        |
!! | q2m            | specific_humidity_at_2m                                                                                             | 2 meter specific humidity                                                           | kg kg-1     |    1 | real       | kind_phys | in     | F        |
!! | u10m           | x_wind_at_10m                                                                                                       | 10 meter u wind speed                                                               | m s-1       |    1 | real       | kind_phys | in     | F        |
!! | v10m           | y_wind_at_10m                                                                                                       | 10 meter v wind speed                                                               | m s-1       |    1 | real       | kind_phys | in     | F        |
!! | tmpmin         | minimum_temperature_at_2m                                                                                           | min temperature at 2m height                                                        | K           |    1 | real       | kind_phys | inout  | F        |
!! | tmpmax         | maximum_temperature_at_2m                                                                                           | max temperature at 2m height                                                        | K           |    1 | real       | kind_phys | inout  | F        |
!! | spfhmin        | minimum_specific_humidity_at_2m                                                                                     | minimum specific humidity at 2m height                                              | kg kg-1     |    1 | real       | kind_phys | inout  | F        |
!! | spfhmax        | maximum_specific_humidity_at_2m                                                                                     | maximum specific humidity at 2m height                                              | kg kg-1     |    1 | real       | kind_phys | inout  | F        |
!! | wind10mmax     | maximum_wind_at_10m                                                                                                 | maximum wind speed at 10 m                                                          | m s-1       |    1 | real       | kind_phys | inout  | F        |
!! | u10mmax        | maximum_x_wind_at_10m                                                                                               | maximum x wind at 10 m                                                              | m s-1       |    1 | real       | kind_phys | inout  | F        |
!! | v10mmax        | maximum_y_wind_at_10m                                                                                               | maximum y wind at 10 m                                                              | m s-1       |    1 | real       | kind_phys | inout  | F        |
!! | dpt2m          | dewpoint_temperature_at_2m                                                                                          | 2 meter dewpoint temperature                                                        | K           |    1 | real       | kind_phys | inout  | F        |
!! | errmsg         | ccpp_error_message                                                                                                  | error message for error handling in CCPP                                            | none        |    0 | character  | len=*     | out    | F        |
!! | errflg         | ccpp_error_flag                                                                                                     | error flag for error handling in CCPP                                               | flag        |    0 | integer    |           | out    | F        |
!!
#endif
      subroutine sfc_diag_post_run (im, lssav, dtf, con_eps, con_epsm1, pgr,    &
                         t2m, q2m, u10m, v10m, tmpmin, tmpmax, spfhmin, spfhmax,&
                         wind10mmax, u10mmax, v10mmax, dpt2m, errmsg, errflg)

        use machine,               only: kind_phys

        implicit none

        integer,                              intent(in) :: im
        logical,                              intent(in) :: lssav
        real(kind=kind_phys),                 intent(in) :: dtf, con_eps, con_epsm1
        real(kind=kind_phys), dimension(im),  intent(in) :: pgr, t2m, q2m, u10m, v10m
        real(kind=kind_phys), dimension(im),  intent(inout) :: tmpmin, tmpmax, spfhmin, spfhmax
        real(kind=kind_phys), dimension(im),  intent(inout) :: wind10mmax, u10mmax, v10mmax, dpt2m

        character(len=*),                     intent(out) :: errmsg
        integer,                              intent(out) :: errflg

        integer :: i
        real(kind=kind_phys) :: tem

        ! Initialize CCPP error handling variables
        errmsg = ''
        errflg = 0

        if (lssav) then
          do i=1,im
            tmpmax(i)  = max(tmpmax(i),t2m(i))
            tmpmin(i)  = min(tmpmin(i),t2m(i))
            spfhmax(i) = max(spfhmax(i),q2m(i))
            spfhmin(i) = min(spfhmin(i),q2m(i))
          enddo
          ! Find max wind speed then decompose
          do i=1, im
             tem = sqrt(u10m(i)*u10m(i) + v10m(i)*v10m(i))
             if (tem > wind10mmax(i)) then
                wind10mmax(i) = tem
                u10mmax(i)    = u10m(i)
                v10mmax(i)    = v10m(i)
             endif
             ! Compute dew point, first using vapor pressure
             tem = max(pgr(i) * q2m(i) / ( con_eps - con_epsm1 *q2m(i)), 1.e-8)
             dpt2m(i) = 243.5 / ( ( 17.67 / log(tem/611.2) ) - 1.) + 273.14
          enddo
        endif

      end subroutine sfc_diag_post_run

      end module sfc_diag_post

!>  \file cires_ugwp_post.F90
!! This file contains
module cires_ugwp_post

contains

!>\defgroup cires_ugwp_post CIRES UGWP Scheme Post
!! @{
!> \section arg_table_cires_ugwp_post_init Argument Table
!!
    subroutine cires_ugwp_post_init ()
    end subroutine cires_ugwp_post_init

!>@brief The subroutine initializes the CIRES UGWP
#if 0
!> \section arg_table_cires_ugwp_post_run Argument Table
!! | local_name       | standard_name                                                                  | long_name                                                              | units     | rank |  type     |   kind    | intent | optional |
!! |------------------|--------------------------------------------------------------------------------|------------------------------------------------------------------------|-----------|------|-----------|-----------|--------|----------|
!! | ldiag_ugwp       | diag_ugwp_flag                                                                 | flag for CIRES UGWP Diagnostics                                        | flag      | 0    | logical   |           | in     | F        |
!! | dtf              | time_step_for_dynamics                                                         | dynamics timestep                                                      | s         | 0    | real      | kind_phys | none   | F        |
!! | im               | horizontal_loop_extent                                                         | horizontal loop extent                                                 | count     | 0    | integer   |           | in     | F        |
!! | levs             | vertical_dimension                                                             | number of vertical levels                                              | count     | 0    | integer   |           | in     | F        |
!! | gw_dudt          | tendency_of_x_wind_due_to_ugwp                                                 | zonal wind tendency due to UGWP                                        | m s-2     | 2    | real      | kind_phys | none   | F        |
!! | tau_tofd         | instantaneous_momentum_flux_due_to_turbulent_orographic_form_drag              | momentum flux or stress due to TOFD                                    | Pa        | 2    | real      | kind_phys | none   | F        |
!! | tau_mtb          | instantaneous_momentum_flux_due_to_mountain_blocking_drag                      | momentum flux or stress due to mountain blocking drag                  | Pa        | 2    | real      | kind_phys | none   | F        |
!! | tau_ogw          | instantaneous_momentum_flux_due_to_orographic_gravity_wave_drag                | momentum flux or stress due to orographic gravity wave drag            | Pa        | 2    | real      | kind_phys | none   | F        |
!! | tau_ngw          | instantaneous_momentum_flux_due_to_nonstationary_gravity_wave                  | momentum flux or stress due to nonstationary gravity waves             | Pa        | 2    | real      | kind_phys | none   | F        |
!! | zmtb             | height_of_mountain_blocking                                                    | height of mountain blocking drag                                       | m         | 1    | real      | kind_phys | none   | F        |
!! | zlwb             | height_of_low_level_wave_breaking                                              | height of low level wave breaking                                      | m         | 1    | real      | kind_phys | none   | F        |
!! | zogw             | height_of_launch_level_of_orographic_gravity_wave                              | height of launch level of orographic gravity wave                      | m         | 1    | real      | kind_phys | none   | F        |
!! | dudt_mtb         | instantaneous_change_in_x_wind_due_to_mountain_blocking_drag                   | instantaneous change in x wind due to mountain blocking drag           | m s-2     | 2    | real      | kind_phys | none   | F        |
!! | dudt_ogw         | instantaneous_change_in_x_wind_due_to_orographic_gravity_wave_drag             | instantaneous change in x wind due to orographic gw drag               | m s-2     | 2    | real      | kind_phys | none   | F        |
!! | dudt_tms         | instantaneous_change_in_x_wind_due_to_turbulent_orographic_form_drag           | instantaneous change in x wind due to TOFD                             | m s-2     | 2    | real      | kind_phys | none   | F        |
!! | cnvgwd           | flag_convective_gravity_wave_drag                                              | flag for conv gravity wave drag                                        | flag      | 0    | logical   |           | in     | F        |
!! | tot_zmtb         | time_integral_of_height_of_mountain_blocking                                   | time integral of height of mountain blocking drag                      | m         | 1    | real      | kind_phys | none   | F        |
!! | tot_zlwb         | time_integral_of_height_of_low_level_wave_breaking                             | time integral of height of drag due to low level wave breaking         | m         | 1    | real      | kind_phys | none   | F        |
!! | tot_zogw         | time_integral_of_height_of_launch_level_of_orographic_gravity_wave             | time integral of height of launch level of orographic gravity wave     | m         | 1    | real      | kind_phys | none   | F        |
!! | tot_tofd         | time_integral_of_momentum_flux_due_to_turbulent_orographic_form_drag           | time integral of momentum flux due to TOFD                             | Pa        | 2    | real      | kind_phys | none   | F        |
!! | tot_mtb          | time_integral_of_momentum_flux_due_to_mountain_blocking_drag                   | time integral of momentum flux due to mountain blocking drag           | Pa        | 2    | real      | kind_phys | none   | F        |
!! | tot_ogw          | time_integral_of_momentum_flux_due_to_orographic_gravity_wave_drag             | time integral of momentum flux due to orographic gravity wave drag     | Pa        | 2    | real      | kind_phys | none   | F        |
!! | tot_ngw          | time_integral_of_momentum_flux_due_to_nonstationary_gravity_wave               | time integral of momentum flux due to nonstationary gravity waves      | Pa        | 2    | real      | kind_phys | none   | F        |
!! | du3dt_mtb        | time_integral_of_change_in_x_wind_due_to_mountain_blocking_drag                | time integral of change in x wind due to mountain blocking drag        | m s-2     | 2    | real      | kind_phys | none   | F        |
!! | du3dt_ogw        | time_integral_of_change_in_x_wind_due_to_orographic_gravity_wave_drag          | time integral of change in x wind due to orographic gw drag            | m s-2     | 2    | real      | kind_phys | none   | F        |
!! | du3dt_tms        | time_integral_of_change_in_x_wind_due_to_turbulent_orographic_form_drag        | time integral of change in x wind due to TOFD                          | m s-2     | 2    | real      | kind_phys | none   | F        |
!! | du3dt_ngw        | time_integral_of_change_in_x_wind_due_to_nonstationary_gravity_wave            | time integral of change in x wind due to NGW                           | m s-2     | 2    | real      | kind_phys | none   | F        |
!! | errmsg           | ccpp_error_message                                                             | error message for error handling in CCPP                               | none      | 0    | character | len=*     | out    | F        |
!! | errflg           | ccpp_error_flag                                                                | error flag for error handling in CCPP                                  | flag      | 0    | integer   |           | out    | F        |
!!
#endif


     subroutine cires_ugwp_post_run (ldiag_ugwp, dtf, im, levs,     &
         gw_dudt, tau_tofd, tau_mtb, tau_ogw, tau_ngw,              &
         zmtb, zlwb, zogw, dudt_mtb, dudt_ogw, dudt_tms,            &
         tot_zmtb, tot_zlwb, tot_zogw,                              &
         tot_tofd, tot_mtb, tot_ogw, tot_ngw,                       &
         du3dt_mtb,du3dt_ogw, du3dt_tms, du3dt_ngw,                 &
         cnvgwd, errmsg, errflg)

        use machine,                only: kind_phys

        implicit none

        ! Interface variables
        integer,              intent(in) :: im, levs
        real(kind=kind_phys), intent(in) :: dtf
        logical,              intent(in) :: ldiag_ugwp      !< flag for CIRES UGWP Diagnostics
        logical,              intent(inout) :: cnvgwd       !< flag to turn on/off convective gwd

        real(kind=kind_phys), intent(in),  dimension(im)       :: zmtb, zlwb, zogw
        real(kind=kind_phys), intent(in),  dimension(im)       :: tau_mtb, tau_ogw, tau_tofd, tau_ngw
        real(kind=kind_phys), intent(out), dimension(im)       :: tot_mtb, tot_ogw, tot_tofd, tot_ngw
        real(kind=kind_phys), intent(out), dimension(im)       :: tot_zmtb, tot_zlwb, tot_zogw
        real(kind=kind_phys), intent(in),  dimension(im, levs) :: gw_dudt, dudt_mtb, dudt_ogw, dudt_tms 
        real(kind=kind_phys), intent(out), dimension(im, levs) :: du3dt_mtb, du3dt_ogw, du3dt_tms, du3dt_ngw

        character(len=*),        intent(out) :: errmsg
        integer,                 intent(out) :: errflg


        ! Initialize CCPP error handling variables
        errmsg = ''
        errflg = 0

        if (.not. (ldiag_ugwp)) return


        if (ldiag_ugwp) then
            tot_zmtb =  tot_zmtb + dtf *zmtb
            tot_zlwb =  tot_zlwb + dtf *zlwb
            tot_zogw =  tot_zogw + dtf *zogw
    
            tot_tofd  = tot_tofd + dtf *tau_tofd
            tot_mtb   = tot_mtb +  dtf *tau_mtb
            tot_ogw   = tot_ogw +  dtf *tau_ogw
            tot_ngw   = tot_ngw +  dtf *tau_ngw
    
            du3dt_mtb = du3dt_mtb + dtf *dudt_mtb
            du3dt_tms = du3dt_tms + dtf *dudt_tms
            du3dt_ogw = du3dt_ogw + dtf *dudt_ogw
            du3dt_ngw = du3dt_ngw + dtf *gw_dudt 
         endif          


         cnvgwd = .false.

      end subroutine cires_ugwp_post_run

!> \section arg_table_cires_ugwp_post_finalize Argument Table
!!
      subroutine cires_ugwp_post_finalize ()
      end subroutine cires_ugwp_post_finalize

!! @}
end module cires_ugwp_post

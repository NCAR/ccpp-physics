!>  \file GFS_surface_loop_control.F90
!!  This file contains the GFS_surface_loop_control scheme.

!> \defgroup GFS_surface_loop_control GFS_surface_loop_control scheme
!! @{

      module GFS_surface_loop_control_part1
      contains

      subroutine GFS_surface_loop_control_part1_init
      end subroutine GFS_surface_loop_control_part1_init

      subroutine GFS_surface_loop_control_part1_finalize
      end subroutine GFS_surface_loop_control_part1_finalize

!> \brief Brief description of the subroutine
!!
#if 0
!! \section arg_table_GFS_surface_loop_control_part1_run Arguments
!! | local_name     | standard_name                                          | long_name                                  | units      | rank | type      |    kind   | intent | optional |
!! |----------------|--------------------------------------------------------|--------------------------------------------|------------|------|-----------|-----------|--------|----------|
!! | im             | horizontal_loop_extent                                 | horizontal loop extent                     | count      |    0 | integer   |           | in     | F        |
!! | iter           | ccpp_loop_counter                                      | loop counter for subcycling loops in CCPP  | index      |    0 | integer   |           | in     | F        |
!! | wind           | wind_speed_at_lowest_model_layer                       | wind speed at lowest model level           | m s-1      |    1 | real      | kind_phys | in     | F        |
!! | flag_guess     | flag_for_guess_run                                     | flag for guess run                         | flag       |    1 | logical   |           | inout  | F        |
!! | errmsg         | ccpp_error_message                                     | error message for error handling in CCPP   | none       |    0 | character | len=*     | out    | F        |
!! | errflg         | ccpp_error_flag                                        | error flag for error handling in CCPP      | flag       |    0 | integer   |           | out    | F        |
!!
#endif
!!  \section general General Algorithm
!!  \section detailed Detailed Algorithm
!!  @{

      subroutine GFS_surface_loop_control_part1_run (im, iter, wind, flag_guess, errmsg, errflg)

      use machine,           only: kind_phys

      implicit none

      ! Interface variables
      integer, intent(in)                                :: im
      integer, intent(in)                                :: iter
      real(kind=kind_phys), dimension(im), intent(in)    :: wind
      logical,              dimension(im), intent(inout) :: flag_guess

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      ! Local variables
      integer :: i

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      do i=1,im
        if (iter == 1 .and. wind(i) < 2.0) then
          flag_guess(i) = .true.
        endif
      enddo

      end subroutine GFS_surface_loop_control_part1_run
!> @}
      end module  GFS_surface_loop_control_part1
!> @}

!> \defgroup GFS_surface_loop_control GFS_surface_loop_control scheme
!! @{

      module GFS_surface_loop_control_part2
      contains

      subroutine GFS_surface_loop_control_part2_init
      end subroutine GFS_surface_loop_control_part2_init

      subroutine GFS_surface_loop_control_part2_finalize
      end subroutine GFS_surface_loop_control_part2_finalize

!> \brief Brief description of the subroutine
!!
#if 0
!! \section arg_table_GFS_surface_loop_control_part2_run Arguments
!! | local_name     | standard_name                                          | long_name                                                            | units      | rank | type      |    kind   | intent | optional |
!! |----------------|--------------------------------------------------------|----------------------------------------------------------------------|------------|------|-----------|-----------|--------|----------|
!! | im             | horizontal_loop_extent                                 | horizontal loop extent                                               | count      |    0 | integer   |           | in     | F        |
!! | iter           | ccpp_loop_counter                                      | loop counter for subcycling loops in CCPP                            | index      |    0 | integer   |           | in     | F        |
!! | wind           | wind_speed_at_lowest_model_layer                       | wind speed at lowest model level                                     | m s-1      |    1 | real      | kind_phys | in     | F        |
!! | flag_guess     | flag_for_guess_run                                     | flag for guess run                                                   | flag       |    1 | logical   |           | inout  | F        |
!! | flag_iter      | flag_for_iteration                                     | flag for iteration                                                   | flag       |    1 | logical   |           | inout  | F        |
!! | dry            | flag_nonzero_land_surface_fraction                     | flag indicating presence of some land surface area fraction          | flag       |    1 | logical   |           | in     | F        |
!! | wet            | flag_nonzero_wet_surface_fraction                      | flag indicating presence of some ocean or lake surface area fraction | flag       |    1 | logical   |           | in     | F        |
!! | icy            | flag_nonzero_sea_ice_surface_fraction                  | flag indicating presence of some sea ice surface area fraction       | flag       |    1 | logical   |           | in     | F        |
!! | nstf_name1     | flag_for_nsstm_run                                     | NSSTM flag: off/uncoupled/coupled=0/1/2                              | flag       |    0 | integer   |           | in     | F        |
!! | errmsg         | ccpp_error_message                                     | error message for error handling in CCPP                             | none       |    0 | character | len=*     | out    | F        |
!! | errflg         | ccpp_error_flag                                        | error flag for error handling in CCPP                                | flag       |    0 | integer   |           | out    | F        |
!!
#endif
!!  \section general General Algorithm
!!  \section detailed Detailed Algorithm
!!  @{

      subroutine GFS_surface_loop_control_part2_run (im, iter,  wind, &
             flag_guess, flag_iter, dry, wet, icy, nstf_name1, errmsg, errflg)

      use machine,           only: kind_phys

      implicit none

      ! Interface variables
      integer,                             intent(in)    :: im
      integer,                             intent(in)    :: iter
      real(kind=kind_phys), dimension(im), intent(in)    :: wind
      logical,              dimension(im), intent(inout) :: flag_guess
      logical,              dimension(im), intent(inout) :: flag_iter
      logical,              dimension(im), intent(in)    :: dry, wet, icy
      integer,                             intent(in)    :: nstf_name1

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      ! Local variables
      integer :: i

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      do i = 1, im
        flag_iter(i)  = .false.
        flag_guess(i) = .false.

        if (iter == 1 .and. wind(i) < 2.0) then
          if (dry(i) .or. (wet(i) .and. .not.icy(i) .and. nstf_name1 > 0)) then
            flag_iter(i) = .true.
          endif
        endif

      enddo

      end subroutine GFS_surface_loop_control_part2_run
!> @}

      end module GFS_surface_loop_control_part2
!> @}

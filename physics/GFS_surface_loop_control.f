!>  \file GFS_surface_loop_control.f
!!  This file contains the GFS_surface_loop_control scheme.

!> \defgroup GFS_surface_loop_control GFS_surface_loop_control scheme
!! @{
!!  \brief Brief description of the parameterization
!!  \section diagram Calling Hierarchy Diagram
!!  \section intraphysics Intraphysics Communication

      module GFS_surface_loop_control_part1
      contains

      subroutine GFS_surface_loop_control_part1_init
      end subroutine GFS_surface_loop_control_part1_init

      subroutine GFS_surface_loop_control_part1_finalize
      end subroutine GFS_surface_loop_control_part1_finalize

!> \brief Brief description of the subroutine
!!
!! \section arg_table_lsmnoah_pre_run Arguments
!!| local var name | longname                                                    | description                                | units      | rank | type    |    kind   | intent | optional |
!!|----------------|-------------------------------------------------------------|--------------------------------------------|------------|------|---------|-----------|--------|----------|
!!| im             | horizontal_loop_extent                                      | horizontal loop extent, start at 1         | index      |    0 | integer |           | in     | F        |
!!| iter           | iteration_number                                            | number of iteration                        | index      |    0 | real    | kind_phys | in     | F        |
!!| wind           | surface_wind_speed                                          | surface wind speed                         | m s-1      | 1    | real    | kind_phys | in     | F        |
!!| flag_guess     | flag_for_guess_run                                          | flag for guess run                         | flag       | 1    | logical |           | inout  | F        |
!!
!!  \section general General Algorithm
!!  \section detailed Detailed Algorithm
!!  @{

      subroutine GFS_surface_loop_control_part1_run                     &
     & ( im,iter,wind,flag_guess
     & ) 

      use machine,           only: kind_phys

!  ---  interface variables
      integer, intent(in) :: im, iter
      real(kind=kind_phys), dimension(im), intent(in)  ::               &
     &   wind 
      logical, dimension(im), intent(inout)  ::                         &
     &   flag_guess 

      do i = 1, im
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
!!  \brief Brief description of the parameterization
!!  \section diagram Calling Hierarchy Diagram
!!  \section intraphysics Intraphysics Communication

      module GFS_surface_loop_control_part2
      contains

      subroutine GFS_surface_loop_control_part2_init
      end subroutine GFS_surface_loop_control_part2_init

      subroutine GFS_surface_loop_control_part2_finalize
      end subroutine GFS_surface_loop_control_part2_finalize

!> \brief Brief description of the subroutine
!!
!! \section arg_table_lsmnoah_pre_run Arguments
!!| local var name | longname                                                    | description                                | units      | rank | type    |    kind   | intent | optional |
!!|----------------|-------------------------------------------------------------|--------------------------------------------|------------|------|---------|-----------|--------|----------|
!!| im             | horizontal_loop_extent                                      | horizontal loop extent, start at 1         | index      |    0 | integer |           | in     | F        |
!!| iter           | iteration_number                                            | number of iteration                        | index      |    0 | real    | kind_phys | in     | F        |
!!| wind           | surface_wind_speed                                          | surface wind speed                         | m s-1      | 1    | real    | kind_phys | in     | F        |
!!| flag_guess     | flag_for_guess_run                                          | flag for guess run                         | flag       | 1    | logical |           | inout  | F        |
!!| flag_iter      | flag_for_iteration                                          | flag for iteration                         | flag       | 1    | logical |           | in     | F        |
!!| islmsk         | sea_land_ice_mask                                           | landmask: sea/land/ice=0/1/2               | flag       | 1    | integer |           | in     | F        |
!!| nstf_name1     | flag_for_nsstm_run                                          | NSSTM flag: off/uncoupled/coupled=0/1/2    | flag       |    0 | integer |           | in     | F        |
!!
!!  \section general General Algorithm
!!  \section detailed Detailed Algorithm
!!  @{

      subroutine GFS_surface_loop_control_part2_run                     $
     $ (im,iter,wind,flag_guess,flag_iter,islmsk,nstf_name1
     $ ) 

      use machine,           only: kind_phys

!  ---  interface variables
      integer, intent(in) :: im, iter, nstf_name1
      integer, dimension(im), intent(in) :: islmsk
!      integer, dimension(im) :: islmsk
      real(kind=kind_phys), dimension(im), intent(in)  ::               &
     &   wind
      logical, dimension(im), intent(inout)  ::                         &
     &   flag_guess,flag_iter 

      do i = 1, im
        flag_iter(i)  = .false.
        flag_guess(i) = .false.

        if (iter == 1 .and. wind(i) < 2.0) then
          if ((islmsk(i) == 1) .or. ((islmsk(i) == 0) .and.             &
     &                               (nstf_name1 > 0))) then
            flag_iter(i) = .true.
          endif
        endif

      enddo

      end subroutine GFS_surface_loop_control_part2_run
!> @}

      end module GFS_surface_loop_control_part2
!> @}

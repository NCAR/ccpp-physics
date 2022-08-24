!>  \file GFS_surface_loop_control_part1.F90
!!  This file contains the GFS_surface_loop_control_part1 scheme.

!> \defgroup GFS_surface_loop_control GFS_surface_loop_control_part1 scheme
!! This module contains the GFS_surface_loop_control_part1 scheme.
!! @{
      module GFS_surface_loop_control_part1
      contains

!> \brief Brief description of the subroutine
!!
!! \section arg_table_GFS_surface_loop_control_part1_run Arguments
!! \htmlinclude GFS_surface_loop_control_part1_run.html
!!
!!  \section gen_loop1 General Algorithm
!!  \section detailed_loop1 Detailed Algorithm
      subroutine GFS_surface_loop_control_part1_run (im, iter,       &
                                   wind, flag_guess, errmsg, errflg)

      use machine,           only: kind_phys

      implicit none

      ! Interface variables
      integer, intent(in)                               :: im
      integer, intent(in)                               :: iter
      real(kind=kind_phys), dimension(:), intent(in)    :: wind
      logical,              dimension(:), intent(inout) :: flag_guess

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      ! Local variables
      integer :: i

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      do i=1,im
        if (iter == 1 .and. wind(i) < 2.0d0) then
          flag_guess(i) = .true.
        endif
      enddo

      end subroutine GFS_surface_loop_control_part1_run
      end module  GFS_surface_loop_control_part1
!> @}

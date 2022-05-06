!>  \file GFS_surface_loop_control_part2.F90
!!  This file contains the GFS_surface_loop_control_part2 scheme.

!> \defgroup GFS_surface_loop_control2 GFS_surface_loop_control_part2 Module
!! This module contains the GFS_surface_loop_control_part2 scheme.
!> @{
      module GFS_surface_loop_control_part2
      contains

#if 0
!> \section arg_table_GFS_surface_loop_control_part2_run Arguments
!! \htmlinclude GFS_surface_loop_control_part2_run.html
!!
#endif
!>  \section looptwo_general General Algorithm
      subroutine GFS_surface_loop_control_part2_run (im, lsm, lsm_noahmp, iter,&
       wind, flag_guess, flag_iter, dry, wet, icy, nstf_name1, errmsg, errflg)

      use machine,           only: kind_phys

      implicit none

      ! Interface variables
      integer,                             intent(in)    :: im
      integer,                             intent(in)    :: iter
      integer,                             intent(in)    :: lsm
      integer,                             intent(in)    :: lsm_noahmp
      real(kind=kind_phys), dimension(:),  intent(in)    :: wind
      logical,              dimension(:),  intent(inout) :: flag_guess
      logical,              dimension(:),  intent(inout) :: flag_iter
      logical,              dimension(:),  intent(in)    :: dry, wet, icy
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

        if (iter == 1 .and. wind(i) < 2.0d0) then
          !if (dry(i) .or. (wet(i) .and. .not.icy(i) .and. nstf_name1 > 0)) then
          if((dry(i) .and. lsm /= lsm_noahmp) .or. (wet(i) .and. nstf_name1 > 0)) then
            flag_iter(i) = .true.
          endif
        endif

      enddo

      end subroutine GFS_surface_loop_control_part2_run
      end module GFS_surface_loop_control_part2
!> @}

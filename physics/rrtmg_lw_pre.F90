!>\file rrtmg_lw_pre.F90
!!
      module rrtmg_lw_pre
      contains

!>\defgroup rrtmg_lw_pre GFS RRTMG-LW scheme pre
!! This module contains RRTMG-LW pre module.
!> @{
!> \section arg_table_rrtmg_lw_pre_run Argument Table
!! \htmlinclude rrtmg_lw_pre_run.html
!!
      subroutine rrtmg_lw_pre_run (errmsg, errflg)

      implicit none

      character(len=*), intent(  out) :: errmsg
      integer,          intent(  out) :: errflg

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      end subroutine rrtmg_lw_pre_run

!> @}
      end module rrtmg_lw_pre

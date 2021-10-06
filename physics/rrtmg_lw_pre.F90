!>\file rrtmg_lw_pre.f90
      module rrtmg_lw_pre
      contains

!>\defgroup rrtmg_lw_pre GFS RRTMG scheme pre
!! @{
      subroutine rrtmg_lw_pre_init ()
      end subroutine rrtmg_lw_pre_init

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

      subroutine rrtmg_lw_pre_finalize ()
      end subroutine rrtmg_lw_pre_finalize
!! @}
      end module rrtmg_lw_pre

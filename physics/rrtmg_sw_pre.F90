!>\file rrtmg_sw_pre.f90
!! This file contains a subroutine to module_radiation_surface::setalb() to
!! setup surface albedo for SW radiation.
      module rrtmg_sw_pre
      contains

!>\defgroup rrtmg_sw_pre GFS RRTMG scheme Pre
!! @{
      subroutine rrtmg_sw_pre_init ()
      end subroutine rrtmg_sw_pre_init

!> \section arg_table_rrtmg_sw_pre_run Argument Table
!! \htmlinclude rrtmg_sw_pre_run.html
!!
      subroutine rrtmg_sw_pre_run (im, lsswr, coszen, nday, idxday, errmsg, errflg)

      use machine,                   only: kind_phys

      implicit none

      integer,                              intent(in)    :: im
      logical,                              intent(in)    :: lsswr
      real(kind=kind_phys), dimension(im),  intent(in)    :: coszen
      integer,                              intent(out)   :: nday
      integer, dimension(:),                intent(out)   :: idxday
      character(len=*),                     intent(out)   :: errmsg
      integer,                              intent(out)   :: errflg

      ! Local variables
      integer :: i

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

!  --- ...  start radiation calculations
!           remember to set heating rate unit to k/sec!

!> -# Start SW radiation calculations
      if (lsswr) then
!>  - Check for daytime points for SW radiation.
        nday = 0
        idxday = 0
        do i = 1, IM
          if (coszen(i) >= 0.0001) then
            nday = nday + 1
            idxday(nday) = i
          endif
        enddo
      else
        nday   = 0
        idxday = 0
      endif

      end subroutine rrtmg_sw_pre_run

      subroutine rrtmg_sw_pre_finalize ()
      end subroutine rrtmg_sw_pre_finalize

!! @}
      end module rrtmg_sw_pre

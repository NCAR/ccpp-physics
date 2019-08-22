!> \file noahmp_pre.F90
!! This file contains subroutines that prepare data for the NoahMP land surface model scheme
!! as part of the GFS physics suite.
      module noahmp_pre

      implicit none

      contains

      subroutine noahmp_pre_init()
      end subroutine noahmp_pre_init

      subroutine noahmp_pre_finalize()
      end subroutine noahmp_pre_finalize

!> \section arg_table_noahmp_pre_run Argument Table
!! \htmlinclude noahmp_pre_run.html
!!      
      subroutine noahmp_pre_run (jdat, julian, yearlen, errmsg, errflg)

      use machine, only : kind_phys
      implicit none
      
      integer,          intent(in)  :: jdat(1:8)
      
      real(kind=kind_phys), intent(out) :: julian
      integer             , intent(out) :: yearlen
      
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      
      integer ::  iw3jdn      
      integer :: jd0, jd1
      real    :: fjd
      
      
      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0
      
      ! Julian day calculation (fcst day of the year)
      ! we need yearln and julian to
      ! pass to noah mp sflx, idate is init, jdat is fcst;idate = jdat when kdt=1
      ! jdat is changing
      !

      jd1    = iw3jdn(jdat(1),jdat(2),jdat(3))
      jd0    = iw3jdn(jdat(1),1,1)
      fjd    = float(jdat(5))/24.0 + float(jdat(6))/1440.0

      julian = float(jd1-jd0) + fjd

      !
      ! Year length
      !
      ! what if the integration goes from one year to another?
      ! iyr or jyr ? from 365 to 366 or from 366 to 365
      !
      ! is this against model's noleap yr assumption?
      if (mod(jdat(1),4) == 0) then
        yearlen = 366
        if (mod(jdat(1),100) == 0) then
          yearlen = 365
          if (mod(jdat(1),400) == 0) then
            yearlen = 366
          endif
        endif
      endif

    end subroutine noahmp_pre_run

    end module noahmp_pre
!> \file GFS_time_vary_pre.scm.F90
!!  Contains code related to GFS physics suite setup (generic part of time_vary_step)

   module GFS_time_vary_pre

      use funcphys, only: gfuncphys

      implicit none

      private

      public GFS_time_vary_pre_init, GFS_time_vary_pre_timestep_init, GFS_time_vary_pre_finalize

      logical :: is_initialized = .false.

      contains

!> \section arg_table_GFS_time_vary_pre_init Argument Table
!! \htmlinclude GFS_time_vary_pre_init.html
!!
      subroutine GFS_time_vary_pre_init (errmsg, errflg)

         implicit none

         character(len=*),                 intent(out)   :: errmsg
         integer,                          intent(out)   :: errflg

         ! Initialize CCPP error handling variables
         errmsg = ''
         errflg = 0

         if (is_initialized) return
 
         !--- Call gfuncphys (funcphys.f) to compute all physics function tables.
         call gfuncphys ()

         is_initialized = .true.

      end subroutine GFS_time_vary_pre_init


!> \section arg_table_GFS_time_vary_pre_finalize Argument Table
!! \htmlinclude GFS_time_vary_pre_finalize.html
!!
      subroutine GFS_time_vary_pre_finalize(errmsg, errflg)

         implicit none

         character(len=*),                 intent(out)   :: errmsg
         integer,                          intent(out)   :: errflg

         ! Initialize CCPP error handling variables
         errmsg = ''
         errflg = 0

         if (.not. is_initialized) return

         ! DH* this is the place to deallocate whatever is allocated by gfuncphys() in GFS_time_vary_pre_init

         is_initialized = .false.

      end subroutine GFS_time_vary_pre_finalize


!> \section arg_table_GFS_time_vary_pre_timestep_init Argument Table
!! \htmlinclude GFS_time_vary_pre_timestep_init.html
!!
      subroutine GFS_time_vary_pre_timestep_init (jdat, idat, dtp, lsm, lsm_noahmp, nsswr, &
        nslwr, idate, debug, me, master, nscyc, sec, phour, zhour, fhour, kdt,   &
        julian, yearlen, ipt, lprnt, lssav, lsswr, lslwr, solhr, errmsg, errflg)

        use machine,               only: kind_phys

        implicit none
        
        integer,                          intent(in)    :: idate(4)
        integer,                          intent(in)    :: jdat(1:8), idat(1:8)
        integer,                          intent(in)    :: lsm, lsm_noahmp,      &
                                                           nsswr, nslwr, me,     &
                                                           master, nscyc
        logical,                          intent(in)    :: debug
        real(kind=kind_phys),             intent(in)    :: dtp
        
        integer,                          intent(out)   :: kdt, yearlen, ipt
        logical,                          intent(out)   :: lprnt, lssav, lsswr,  &
                                                           lslwr
        real(kind=kind_phys),             intent(out)   :: sec, phour, zhour,    &
                                                           fhour, julian, solhr
        
        character(len=*),                 intent(out)   :: errmsg
        integer,                          intent(out)   :: errflg

        real(kind=kind_phys), parameter :: con_24  =   24.0_kind_phys
        real(kind=kind_phys), parameter :: con_hr  = 3600.0_kind_phys
        real(kind=kind_phys) :: rinc(5)
        
        integer ::  iw3jdn      
        integer :: jd0, jd1
        real    :: fjd

        ! Initialize CCPP error handling variables
        errmsg = ''
        errflg = 0

        ! Check initialization status
        if (.not.is_initialized) then
           write(errmsg,'(*(a))') "Logic error: GFS_time_vary_pre_timestep_init called &
                                  &before GFS_time_vary_pre_init"
           errflg = 1
           return
        end if

        !--- jdat is being updated directly inside of the time integration
        !--- loop of gmtb_scm.F90
        !--- update calendars and triggers
        rinc(1:5)   = 0
        call w3difdat(jdat,idat,4,rinc)
        sec = rinc(4)
        phour = sec/con_hr
        !--- set current bucket hour
        zhour = phour
        fhour = (sec + dtp)/con_hr
        kdt   = nint((sec + dtp)/dtp)
        
        !GJF* These calculations were originally in GFS_physics_driver.F90 for 
        !     NoahMP. They were moved to this routine since they only depends 
        !     on time (not space). Note that this code is included as-is from 
        !     GFS_physics_driver.F90, but it may be simplified by using more 
        !     NCEP W3 library calls (e.g., see W3DOXDAT, W3FS13 for Julian day 
        !     of year and W3DIFDAT to determine the integer number of days in 
        !     a given year). *GJF
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

        ipt    = 1
        lprnt  = .false.
        lssav  = .true.

        !--- radiation triggers
        lsswr  = (mod(kdt, nsswr) == 1)
        lslwr  = (mod(kdt, nslwr) == 1)
        !--- allow for radiation to be called on every physics time step, if needed
        if (nsswr == 1)  lsswr = .true.
        if (nslwr == 1)  lslwr = .true.

        !--- set the solar hour based on a combination of phour and time initial hour
        solhr  = mod(phour+idate(1),con_24)

        if ((debug) .and. (me == master)) then
          print *,'   sec ', sec
          print *,'   kdt ', kdt
          print *,' nsswr ', nsswr
          print *,' nslwr ', nslwr
          print *,' nscyc ', nscyc
          print *,' lsswr ', lsswr
          print *,' lslwr ', lslwr
          print *,' fhour ', fhour
          print *,' phour ', phour
          print *,' solhr ', solhr
        endif

      end subroutine GFS_time_vary_pre_timestep_init

    end module GFS_time_vary_pre

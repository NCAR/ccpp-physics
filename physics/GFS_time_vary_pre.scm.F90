!> \file GFS_time_vary_pre.F90
!!  Contains code related to GFS physics suite setup (generic part of time_vary_step)

   module GFS_time_vary_pre

      use funcphys, only: gfuncphys

      implicit none

      private

      public GFS_time_vary_pre_init, GFS_time_vary_pre_run, GFS_time_vary_pre_finalize

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


!> \section arg_table_GFS_time_vary_pre_run Argument Table
!! \htmlinclude GFS_time_vary_pre_run.html
!!
      subroutine GFS_time_vary_pre_run (Model, errmsg, errflg)

        use machine,               only: kind_phys
        use GFS_typedefs,          only: GFS_control_type

        implicit none

        type(GFS_control_type),           intent(inout) :: Model
        character(len=*),                 intent(out)   :: errmsg
        integer,                          intent(out)   :: errflg

        real(kind=kind_phys), parameter :: con_24  =   24.0_kind_phys
        real(kind=kind_phys), parameter :: con_hr  = 3600.0_kind_phys
        real(kind=kind_phys) :: rinc(5)

        integer              :: fjd, iyear, imon, iday, ihr, imin, jd0, jd1
        integer              :: iw3jdn

        ! Initialize CCPP error handling variables
        errmsg = ''
        errflg = 0

        ! Check initialization status
        if (.not.is_initialized) then
           write(errmsg,'(*(a))') "Logic error: GFS_time_vary_pre_run called before GFS_time_vary_pre_init"
           errflg = 1
           return
        end if

        !--- Model%jdat is being updated directly inside of FV3GFS_cap.F90
        !--- update calendars and triggers
        rinc(1:5)   = 0
        call w3difdat(Model%jdat,Model%idat,4,rinc)
        Model%sec = rinc(4)
        Model%phour = Model%sec/con_hr
        !--- set current bucket hour
        Model%zhour = Model%phour
        Model%fhour = (Model%sec + Model%dtp)/con_hr
        Model%kdt   = nint((Model%sec + Model%dtp)/Model%dtp)

        Model%ipt    = 1
        Model%lprnt  = .false.
        Model%lssav  = .true.

        !--- radiation triggers
        Model%lsswr  = (mod(Model%kdt, Model%nsswr) == 1)
        Model%lslwr  = (mod(Model%kdt, Model%nslwr) == 1)
        !--- allow for radiation to be called on every physics time step, if needed
        if (Model%nsswr == 1)  Model%lsswr = .true.
        if (Model%nslwr == 1)  Model%lslwr = .true.

        !--- set the solar hour based on a combination of phour and time initial hour
        Model%solhr  = mod(Model%phour+Model%idate(1),con_24)

        if (Model%lsm == Model%lsm_noahmp) then
!
! Julian day calculation (fcst day of the year)
! we need imn to init lai and sai and yearln and julian to
! pass to noah mp sflx, idate is init, jdat is fcst;idate = jdat when kdt=1
! jdat is changing
!
          Model%imn = Model%idate(2)

          iyear = Model%jdat(1)
          imon  = Model%jdat(2)
          iday  = Model%jdat(3)
          ihr   = Model%jdat(5)
          imin  = Model%jdat(6)

          jd1   = iw3jdn(iyear,imon,iday)
          jd0   = iw3jdn(iyear,1,1)
          fjd   = float(ihr)/24.0 + float(imin)/1440.0

          Model%julian = float(jd1-jd0) + fjd

!
! Year length
!
! what if the integration goes from one year to another?
! iyr or jyr ? from 365 to 366 or from 366 to 365
!
! is this against model's noleap yr assumption?

          if (mod(iyear,400) == 0) then
            Model%yearlen = 366
          elseif (mod(iyear,100) == 0) then
            Model%yearlen = 365
          elseif (mod(iyear,4) == 0) then
            Model%yearlen = 366
          else
            Model%yearlen = 365
          endif
        endif !  if (Model%lsm == Model%lsm_noahmp)
!

        if ((Model%debug) .and. (Model%me == Model%master)) then
          print *,'   sec ', Model%sec
          print *,'   kdt ', Model%kdt
          print *,' nsswr ', Model%nsswr
          print *,' nslwr ', Model%nslwr
          print *,' nscyc ', Model%nscyc
          print *,' lsswr ', Model%lsswr
          print *,' lslwr ', Model%lslwr
          print *,' fhour ', Model%fhour
          print *,' phour ', Model%phour
          print *,' solhr ', Model%solhr
        endif

      end subroutine GFS_time_vary_pre_run

    end module GFS_time_vary_pre

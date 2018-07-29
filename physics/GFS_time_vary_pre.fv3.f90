!> \file GFS_time_vary_pre.f90
!!  Contains code related to GFS physics suite setup (generic part of time_vary_step)

   module GFS_time_vary_pre

      implicit none

      private

      public GFS_time_vary_pre_init, GFS_time_vary_pre_run, GFS_time_vary_pre_finalize

      contains

      subroutine GFS_time_vary_pre_init ()
      end subroutine GFS_time_vary_pre_init

      subroutine GFS_time_vary_pre_finalize()
      end subroutine GFS_time_vary_pre_finalize

!> \section arg_table_GFS_time_vary_pre_run Argument Table
!! | local_name     | standard_name                                          | long_name                                                               | units    | rank |  type                 |   kind    | intent | optional |
!! |----------------|--------------------------------------------------------|-------------------------------------------------------------------------|----------|------|-----------------------|-----------|--------|----------|
!! | Model          | FV3-GFS_Control_type                                   | Fortran DDT containing FV3-GFS model control parameters                 | DDT      |    0 | GFS_control_type      |           | inout  | F        |
!! | errmsg         | ccpp_error_message                                     | error message for error handling in CCPP                                | none     |    0 | character             | len=*     | out    | F        |
!! | errflg         | ccpp_error_flag                                        | error flag for error handling in CCPP                                   | flag     |    0 | integer               |           | out    | F        |
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

        ! Initialize CCPP error handling variables
        errmsg = ''
        errflg = 0

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

        !--- set the solar hour based on a combination of phour and time initial hour
        Model%solhr  = mod(Model%phour+Model%idate(1),con_24)

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

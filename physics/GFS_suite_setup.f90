!> \file GFS_suite_setup.f90
!!  Contains code related to GFS physics suite setup.

    module GFS_suite_setup_1

      contains

      subroutine GFS_suite_setup_1_init ()
      end subroutine GFS_suite_setup_1_init

      subroutine GFS_suite_setup_1_finalize()
      end subroutine GFS_suite_setup_1_finalize

!> \section arg_table_GFS_suite_setup_1_run Argument Table
!! | local var name | longname                                               | description                                                           | units         | rank | type                          |    kind   | intent | optional |
!! |----------------|--------------------------------------------------------|-----------------------------------------------------------------------|---------------|------|-------------------------------|-----------|--------|----------|
!! | Model          | FV3-GFS_Control_type                                   | Fortran DDT containing FV3-GFS model control parameters               | DDT           |    0 | GFS_control_type              |           | inout  | F        |
!! | sec            | seconds_elapsed_since_model_initialization             | seconds elapsed since model initialization                            | s             |    0 | real                          | kind_phys | inout  | F        |
!! | blkno          | block_number                                           | for explicit data blocking: block number of this block                | index         |    0 | integer                       |           | in     | F        |
!!
      subroutine GFS_suite_setup_1_run (Model, sec, blkno)

        use machine,               only: kind_phys
        use GFS_typedefs,          only: GFS_control_type

        implicit none

        type(GFS_control_type),           intent(inout) :: Model
        real(kind=kind_phys),             intent(inout) :: sec
        integer,                          intent(in)    :: blkno

        real(kind=kind_phys), parameter :: con_24  =   24.0_kind_phys
        real(kind=kind_phys), parameter :: con_hr  = 3600.0_kind_phys
        real(kind=kind_phys) :: rinc(5)

        ! DH* check if it is sufficient to do this only for blkno=1 - sec will be part of IPD_Data%Interstitial in the full CCPP code
        !--- Model%jdat is being updated directly inside of FV3GFS_cap.F90
        !--- update calendars and triggers
        rinc(1:5)   = 0
        call w3difdat(Model%jdat,Model%idat,4,rinc)
        sec = rinc(4)

        if (blkno==1) then
          Model%phour = sec/con_hr
          !--- set current bucket hour
          Model%zhour = Model%phour
          Model%fhour = (sec + Model%dtp)/con_hr
          Model%kdt   = nint((sec + Model%dtp)/Model%dtp)

          Model%ipt    = 1
          Model%lprnt  = .false.
          Model%lssav  = .true.

          !--- radiation triggers
          Model%lsswr  = (mod(Model%kdt, Model%nsswr) == 1)
          Model%lslwr  = (mod(Model%kdt, Model%nslwr) == 1)

          !--- set the solar hour based on a combination of phour and time initial hour
          Model%solhr  = mod(Model%phour+Model%idate(1),con_24)

          if ((Model%debug) .and. (Model%me == Model%master)) then
            print *,'   sec ', sec
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
        endif

      end subroutine GFS_suite_setup_1_run

    end module GFS_suite_setup_1

    module GFS_suite_setup_2

      contains

      subroutine GFS_suite_setup_2_init ()
      end subroutine GFS_suite_setup_2_init

      subroutine GFS_suite_setup_2_finalize()
      end subroutine GFS_suite_setup_2_finalize

!> \section arg_table_GFS_suite_setup_2_run Argument Table
!! | local var name | longname                                                                | description                                                             | units         | rank | type                          |    kind   | intent | optional |
!! |----------------|-------------------------------------------------------------------------|-------------------------------------------------------------------------|---------------|------|-------------------------------|-----------|--------|----------|
!! | Grid           | FV3-GFS_Grid_type                                                       | Fortran DDT containing FV3-GFS grid and interpolation related data      | DDT           |    0 | GFS_grid_type                 |           | in     | F        |
!! | Model          | FV3-GFS_Control_type                                                    | Fortran DDT containing FV3-GFS model control parameters                 | DDT           |    0 | GFS_control_type              |           | inout  | F        |
!! | Tbd            | FV3-GFS_Tbd_type                                                        | Fortran DDT containing FV3-GFS miscellaneous data                       | DDT           |    0 | GFS_tbd_type                  |           | inout  | F        |
!! | Sfcprop        | FV3-GFS_Sfcprop_type                                                    | Fortran DDT containing FV3-GFS surface fields                           | DDT           |    0 | GFS_sfcprop_type              |           | inout  | F        |
!! | Cldprop        | FV3-GFS_Cldprop_type                                                    | Fortran DDT containing FV3-GFS cloud fields                             | DDT           |    0 | GFS_cldprop_type              |           | inout  | F        |
!! | Diag           | FV3-GFS_Diag_type                                                       | Fortran DDT containing FV3-GFS fields targeted for diagnostic output    | DDT           |    0 | GFS_diag_type                 |           | inout  | F        |
!!
      subroutine GFS_suite_setup_2_run (Grid, Model, Tbd, Sfcprop, Cldprop, Diag, Sfccycle)

        use mersenne_twister, only: random_setseed, random_number
        use machine,               only: kind_phys
        use physcons,              only: dxmin, dxinv
        use GFS_typedefs,          only: GFS_control_type, GFS_grid_type, &
                                         GFS_Tbd_type, GFS_sfcprop_type,  &
                                         GFS_cldprop_type, GFS_diag_type, &
                                         GFS_sfccycle_type

        implicit none

        type(GFS_grid_type),              intent(in)    :: Grid
        type(GFS_control_type),           intent(inout) :: Model
        type(GFS_tbd_type),               intent(inout) :: Tbd
        type(GFS_sfcprop_type),           intent(inout) :: Sfcprop
        type(GFS_cldprop_type),           intent(inout) :: Cldprop
        type(GFS_diag_type),              intent(inout) :: Diag
        type(GFS_sfccycle_type),          intent(inout) :: Sfccycle

        real(kind=kind_phys), parameter :: con_hr  = 3600.0_kind_phys
        real(kind=kind_phys), parameter :: con_99  =   99.0_kind_phys
        real(kind=kind_phys), parameter :: con_100 =  100.0_kind_phys

        integer :: i, j, k, iseed, iskip, ix, nb
        real(kind=kind_phys) :: wrk(1)
        real(kind=kind_phys) :: rannie(Model%cny)
        real(kind=kind_phys) :: rndval(Model%cnx*Model%cny*Model%nrcm)

        if (Tbd%blkno==1) then
          !--- switch for saving convective clouds - cnvc90.f
          !--- aka Ken Campana/Yu-Tai Hou legacy
          if ((mod(Model%kdt,Model%nsswr) == 0) .and. (Model%lsswr)) then
            !--- initialize,accumulate,convert
            Model%clstp = 1100 + min(Model%fhswr/con_hr,Model%fhour,con_99)
          elseif (mod(Model%kdt,Model%nsswr) == 0) then
            !--- accumulate,convert
            Model%clstp = 0100 + min(Model%fhswr/con_hr,Model%fhour,con_99)
          elseif (Model%lsswr) then
            !--- initialize,accumulate
            Model%clstp = 1100
          else
            !--- accumulate
            Model%clstp = 0100
          endif
        endif

        !--- random number needed for RAS and old SAS and when cal_pre=.true.
        if ( ((Model%imfdeepcnv <= 0) .or. (Model%cal_pre)) .and. (Model%random_clds) ) then
          iseed = mod(con_100*sqrt(Model%fhour*con_hr),1.0d9) + Model%seed0
          call random_setseed(iseed)
          call random_number(wrk)
          do i = 1,Model%cnx*Model%nrcm
            iseed = iseed + nint(wrk(1)) * i
            call random_setseed(iseed)
            call random_number(rannie)
            rndval(1+(i-1)*Model%cny:i*Model%cny) = rannie(1:Model%cny)
          enddo

          ! DH* TODO - this could be sped up by saving jsc, jec, isc, iec in Tbd (for example)
          ! and looping just over them; ix would then run from 1 to blksz(nb); one could also
          ! use OpenMP to speed up this loop or the inside loops *DH
          do k = 1,Model%nrcm
            iskip = (k-1)*Model%cnx*Model%cny
            ix = 0
            nb = 1
            do j = 1,Model%ny
              do i = 1,Model%nx
                ix = ix + 1
                if (ix .gt. Tbd%blksz(nb)) then
                  ix = 1
                  nb = nb + 1
                endif
                if (nb == Tbd%blkno) then
                  Tbd%rann(ix,k) = rndval(i+Model%isc-1 + (j+Model%jsc-2)*Model%cnx + iskip)
                endif
              enddo
            enddo
          enddo
        endif  ! imfdeepcnv, cal_re, random_clds

        !--- o3 interpolation
        if (Model%ntoz > 0) then
          call ozinterpol (Model%me, Tbd%blksz(Tbd%blkno), Model%idate, Model%fhour, &
                           Grid%jindx1_o3, Grid%jindx2_o3, Tbd%ozpl, Grid%ddy_o3)
        endif

        !--- h2o interpolation
        if (Model%h2o_phys) then
          call h2ointerpol (Model%me, Tbd%blksz(Tbd%blkno), Model%idate, Model%fhour, &
                            Grid%jindx1_h, Grid%jindx2_h, Tbd%h2opl, Grid%ddy_h)
        endif

        !--- repopulate specific time-varying sfc properties for AMIP/forecast runs
        if (Model%nscyc >  0) then
          if (mod(Model%kdt,Model%nscyc) == 1) THEN
            call gcycle (Tbd%blksz(Tbd%blkno), Model, Grid, Sfcprop, Cldprop, Sfccycle)
          endif
        endif

        !--- determine if diagnostics buckets need to be cleared
        if (mod(Model%kdt,Model%nszero) == 1) then
          call Diag%rad_zero  (Model)
          call Diag%phys_zero (Model)
        !!!!  THIS IS THE POINT AT WHICH DIAG%ZHOUR NEEDS TO BE UPDATED
        endif

      end subroutine GFS_suite_setup_2_run

    end module GFS_suite_setup_2

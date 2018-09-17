!> \file GFS_phys_time_vary.F90
!!  Contains code related to GFS physics suite setup (physics part of time_vary_step)

   module GFS_phys_time_vary

      use ozne_def, only : levozp, oz_coeff, oz_lat, oz_pres, oz_time, ozplin
      use ozinterp, only : read_o3data, setindxoz, ozinterpol

      use h2o_def,   only : levh2o, h2o_coeff, h2o_lat, h2o_pres, h2o_time, h2oplin
      use h2ointerp, only : read_h2odata, setindxh2o, h2ointerpol

      implicit none

      private

      public GFS_phys_time_vary_init, GFS_phys_time_vary_run, GFS_phys_time_vary_finalize

      logical :: is_initialized = .false.

      contains

!> \section arg_table_GFS_phys_time_vary_init Argument Table
!! | local_name     | standard_name                                          | long_name                                                               | units    | rank |  type                 |   kind    | intent | optional |
!! |----------------|--------------------------------------------------------|-------------------------------------------------------------------------|----------|------|-----------------------|-----------|--------|----------|
!! | Data           | FV3-GFS_Data_type_all_blocks                           | Fortran DDT containing FV3-GFS data                                     | DDT      |    1 | GFS_data_type         |           | inout  | F        |
!! | Model          | FV3-GFS_Control_type                                   | Fortran DDT containing FV3-GFS model control parameters                 | DDT      |    0 | GFS_control_type      |           | inout  | F        |
!! | Interstitial   | FV3-GFS_Interstitial_type_all_threads                  | Fortran DDT containing FV3-GFS interstitial data                        | DDT      |    1 | GFS_interstitial_type |           | inout  | F        |
!! | errmsg         | ccpp_error_message                                     | error message for error handling in CCPP                                | none     |    0 | character             | len=*     | out    | F        |
!! | errflg         | ccpp_error_flag                                        | error flag for error handling in CCPP                                   | flag     |    0 | integer               |           | out    | F        |
!!
      subroutine GFS_phys_time_vary_init (Data, Model, Interstitial, errmsg, errflg)

         use GFS_typedefs,          only: GFS_control_type, GFS_data_type, GFS_interstitial_type

         implicit none

         ! Interface variables
         type(GFS_data_type),              intent(inout) :: Data(:)
         type(GFS_control_type),           intent(inout) :: Model
         type(GFS_interstitial_type),      intent(inout) :: Interstitial(:)
         character(len=*),                 intent(out)   :: errmsg
         integer,                          intent(out)   :: errflg

         ! Local variables
         integer :: nb, nblks, nt, nthrds

         ! Initialize CCPP error handling variables
         errmsg = ''
         errflg = 0

         if (is_initialized) return

         nblks = size(Model%blksz)
         nthrds = size(Interstitial)

         call read_o3data  (Model%ntoz, Model%me, Model%master)

         ! Consistency check that the hardcoded values for levozp and
         ! oz_coeff in GFS_typedefs.F90 match what is set by read_o3data
         ! in GFS_typedefs.F90: allocate (Tbd%ozpl  (IM,levozp,oz_coeff))
         if (size(Data(1)%Tbd%ozpl, dim=2).ne.levozp) then
            write(errmsg,'(2a,i0,a,i0)') "Value error in GFS_phys_time_vary_init: ",    &
                  "levozp from read_o3data does not match value in GFS_typedefs.F90: ", &
                  levozp, " /= ", size(Data(1)%Tbd%ozpl, dim=2)
            errflg = 1
            return
         end if
         if (size(Data(1)%Tbd%ozpl, dim=3).ne.oz_coeff) then
            write(errmsg,'(2a,i0,a,i0)') "Value error in GFS_phys_time_vary_init: ",      &
                  "oz_coeff from read_o3data does not match value in GFS_typedefs.F90: ", &
                  oz_coeff, " /= ", size(Data(1)%Tbd%ozpl, dim=3)
            errflg = 1
            return
         end if

         call read_h2odata (Model%h2o_phys, Model%me, Model%master)

         ! Consistency check that the hardcoded values for levh2o and
         ! h2o_coeff in GFS_typedefs.F90 match what is set by read_o3data
         ! in GFS_typedefs.F90: allocate (Tbd%h2opl (IM,levh2o,h2o_coeff))
         if (size(Data(1)%Tbd%h2opl, dim=2).ne.levh2o) then
            write(errmsg,'(2a,i0,a,i0)') "Value error in GFS_phys_time_vary_init: ",     &
                  "levh2o from read_h2odata does not match value in GFS_typedefs.F90: ", &
                  levh2o, " /= ", size(Data(1)%Tbd%h2opl, dim=2)
            errflg = 1
            return
         end if
         if (size(Data(1)%Tbd%h2opl, dim=3).ne.h2o_coeff) then
            write(errmsg,'(2a,i0,a,i0)') "Value error in GFS_phys_time_vary_init: ",       &
                  "h2o_coeff from read_h2odata does not match value in GFS_typedefs.F90: ", &
                  h2o_coeff, " /= ", size(Data(1)%Tbd%h2opl, dim=3)
            errflg = 1
            return
         end if

         ! DH* OpenMP parallel region with OpenMP do loops for the three do loops?
         ! there is no threading on the outside, i.e. can use CCPP's omp_threads here

         ! Update values of oz_pres in Interstitial data type for all threads
         if (Model%ntoz > 0) then
            do nt=1,nthrds
               Interstitial(nt)%oz_pres = oz_pres
            end do
         end if

         ! Update values of h2o_pres in Interstitial data type for all threads
         if (Model%h2o_phys) then
            do nt=1,nthrds
               Interstitial(nt)%h2o_pres = h2o_pres
            end do
         end if

         !--- read in and initialize ozone and water
         if (Model%ntoz > 0) then
           do nb = 1, nblks
             call setindxoz (Model%blksz(nb), Data(nb)%Grid%xlat_d, Data(nb)%Grid%jindx1_o3, &
                             Data(nb)%Grid%jindx2_o3, Data(nb)%Grid%ddy_o3)
           enddo
         endif

         if (Model%h2o_phys) then
           do nb = 1, nblks
             call setindxh2o (Model%blksz(nb), Data(nb)%Grid%xlat_d, Data(nb)%Grid%jindx1_h, &
                              Data(nb)%Grid%jindx2_h, Data(nb)%Grid%ddy_h)
           enddo
         endif

         ! *DH

         is_initialized = .true.

      end subroutine GFS_phys_time_vary_init


!> \section arg_table_GFS_phys_time_vary_finalize Argument Table
!! | local_name     | standard_name                                          | long_name                                                               | units    | rank |  type                 |   kind    | intent | optional |
!! |----------------|--------------------------------------------------------|-------------------------------------------------------------------------|----------|------|-----------------------|-----------|--------|----------|
!! | errmsg         | ccpp_error_message                                     | error message for error handling in CCPP                                | none     |    0 | character             | len=*     | out    | F        |
!! | errflg         | ccpp_error_flag                                        | error flag for error handling in CCPP                                   | flag     |    0 | integer               |           | out    | F        |
!!
      subroutine GFS_phys_time_vary_finalize(errmsg, errflg)

         implicit none

         ! Interface variables
         character(len=*),                 intent(out)   :: errmsg
         integer,                          intent(out)   :: errflg

         ! Initialize CCPP error handling variables
         errmsg = ''
         errflg = 0

         if (.not.is_initialized) return

         ! Deallocate ozone arrays
         if (allocated(oz_lat)  ) deallocate(oz_lat)
         if (allocated(oz_pres) ) deallocate(oz_pres)
         if (allocated(oz_time) ) deallocate(oz_time)
         if (allocated(ozplin)  ) deallocate(ozplin)

         ! Deallocate h2o arrays
         if (allocated(h2o_lat) ) deallocate(h2o_lat)
         if (allocated(h2o_pres)) deallocate(h2o_pres)
         if (allocated(h2o_time)) deallocate(h2o_time)
         if (allocated(h2oplin) ) deallocate(h2oplin)

         is_initialized = .false.

      end subroutine GFS_phys_time_vary_finalize


!> \section arg_table_GFS_phys_time_vary_run Argument Table
!! | local_name     | standard_name                                          | long_name                                                               | units    | rank |  type                 |   kind    | intent | optional |
!! |----------------|--------------------------------------------------------|-------------------------------------------------------------------------|----------|------|-----------------------|-----------|--------|----------|
!! | Data           | FV3-GFS_Data_type_all_blocks                           | Fortran DDT containing FV3-GFS data                                     | DDT      |    1 | GFS_data_type         |           | inout  | F        |
!! | Model          | FV3-GFS_Control_type                                   | Fortran DDT containing FV3-GFS model control parameters                 | DDT      |    0 | GFS_control_type      |           | inout  | F        |
!! | errmsg         | ccpp_error_message                                     | error message for error handling in CCPP                                | none     |    0 | character             | len=*     | out    | F        |
!! | errflg         | ccpp_error_flag                                        | error flag for error handling in CCPP                                   | flag     |    0 | integer               |           | out    | F        |
!!
      subroutine GFS_phys_time_vary_run (Data, Model, errmsg, errflg)

        use mersenne_twister,      only: random_setseed, random_number
        use machine,               only: kind_phys
        use GFS_typedefs,          only: GFS_control_type, GFS_data_type

        implicit none

        ! Interface variables
        type(GFS_data_type),              intent(in)    :: Data(:)
        type(GFS_control_type),           intent(inout) :: Model
        character(len=*),                 intent(out)   :: errmsg
        integer,                          intent(out)   :: errflg

        ! Local variables
        real(kind=kind_phys), parameter :: con_hr  = 3600.0_kind_phys
        real(kind=kind_phys), parameter :: con_99  =   99.0_kind_phys
        real(kind=kind_phys), parameter :: con_100 =  100.0_kind_phys

        integer :: i, j, k, iseed, iskip, ix, nb, nblks
        real(kind=kind_phys) :: wrk(1)
        real(kind=kind_phys) :: rannie(Model%cny)
        real(kind=kind_phys) :: rndval(Model%cnx*Model%cny*Model%nrcm)

        ! Initialize CCPP error handling variables
        errmsg = ''
        errflg = 0

        ! Check initialization status
        if (.not.is_initialized) then
           write(errmsg,'(*(a))') "Logic error: GFS_phys_time_vary_run called before GFS_phys_time_vary_init"
           errflg = 1
           return
        end if

        nblks = size(Model%blksz)

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
                if (ix .gt. Model%blksz(nb)) then
                  ix = 1
                  nb = nb + 1
                endif
                Data(nb)%Tbd%rann(ix,k) = rndval(i+Model%isc-1 + (j+Model%jsc-2)*Model%cnx + iskip)
              enddo
            enddo
          enddo
        endif  ! imfdeepcnv, cal_re, random_clds

        !--- o3 interpolation
        if (Model%ntoz > 0) then
! DH* OpenMP?
          do nb = 1, nblks
            call ozinterpol (Model%me, Model%blksz(nb), Model%idate, Model%fhour, &
                             Data(nb)%Grid%jindx1_o3, Data(nb)%Grid%jindx2_o3,              &
                             Data(nb)%Tbd%ozpl, Data(nb)%Grid%ddy_o3)
          enddo
        endif

        !--- h2o interpolation
        if (Model%h2o_phys) then
! DH* OpenMP?
          do nb = 1, nblks
            call h2ointerpol (Model%me, Model%blksz(nb), Model%idate, Model%fhour, &
                              Data(nb)%Grid%jindx1_h, Data(nb)%Grid%jindx2_h,                &
                              Data(nb)%Tbd%h2opl, Data(nb)%Grid%ddy_h)
          enddo
        endif

        !--- repopulate specific time-varying sfc properties for AMIP/forecast runs
        if (Model%nscyc >  0) then
          if (mod(Model%kdt,Model%nscyc) == 1) THEN
            call gcycle (nblks, Model, Data(:)%Grid, Data(:)%Sfcprop, Data(:)%Cldprop)
          endif
        endif

        !--- determine if diagnostics buckets need to be cleared
        if (mod(Model%kdt,Model%nszero) == 1) then
! DH* OpenMP?
          do nb = 1,nblks
            call Data(nb)%Intdiag%rad_zero  (Model)
            call Data(nb)%Intdiag%phys_zero (Model)
        !!!!  THIS IS THE POINT AT WHICH DIAG%ZHOUR NEEDS TO BE UPDATED
          enddo
        endif

      end subroutine GFS_phys_time_vary_run

   end module GFS_phys_time_vary

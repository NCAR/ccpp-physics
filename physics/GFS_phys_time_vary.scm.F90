!> \file GFS_phys_time_vary.F90
!!  Contains code related to GFS physics suite setup (physics part of time_vary_step)

   module GFS_phys_time_vary

     use ozne_def, only : levozp, oz_coeff, oz_lat, oz_pres, oz_time, ozplin
     use ozinterp, only : read_o3data, setindxoz, ozinterpol

     use h2o_def,   only : levh2o, h2o_coeff, h2o_lat, h2o_pres, h2o_time, h2oplin
     use h2ointerp, only : read_h2odata, setindxh2o, h2ointerpol

     use aerclm_def, only : aerin, aer_pres, ntrcaer, ntrcaerm
     use aerinterp,  only : read_aerdata, setindxaer, aerinterpol

     use iccn_def,   only : ciplin, ccnin, ci_pres
     use iccninterp, only : read_cidata, setindxci, ciinterpol

#if 0
      !--- variables needed for calculating 'sncovr'
      use namelist_soilveg, only: salp_data, snupx
#endif

      implicit none

      private

      public GFS_phys_time_vary_init, GFS_phys_time_vary_run, GFS_phys_time_vary_finalize

      logical :: is_initialized = .false.

      contains

!> \section arg_table_GFS_phys_time_vary_init Argument Table
!! \htmlinclude GFS_phys_time_vary_init.html
!!
      subroutine GFS_phys_time_vary_init (Grid, Model, Interstitial, Tbd, errmsg, errflg)

         use GFS_typedefs,          only: GFS_control_type, GFS_grid_type, &
                                          GFS_Tbd_type, GFS_interstitial_type

         implicit none

         ! Interface variables
         type(GFS_grid_type),              intent(inout) :: Grid
         type(GFS_control_type),           intent(in)    :: Model
         type(GFS_interstitial_type),      intent(inout) :: Interstitial
         type(GFS_tbd_type),               intent(in)    :: Tbd
         character(len=*),                 intent(out)   :: errmsg
         integer,                          intent(out)   :: errflg

         ! Local variables
         integer :: i, j, ix, nb, nt

         ! Initialize CCPP error handling variables
         errmsg = ''
         errflg = 0

         if (is_initialized) return

         nb = 1
         nt = 1
         
         call read_o3data  (Model%ntoz, Model%me, Model%master)

         ! Consistency check that the hardcoded values for levozp and
         ! oz_coeff in GFS_typedefs.F90 match what is set by read_o3data
         ! in GFS_typedefs.F90: allocate (Tbd%ozpl  (IM,levozp,oz_coeff))
         if (size(Tbd%ozpl, dim=2).ne.levozp) then
            write(errmsg,'(2a,i0,a,i0)') "Value error in GFS_phys_time_vary_init: ",    &
                  "levozp from read_o3data does not match value in GFS_typedefs.F90: ", &
                  levozp, " /= ", size(Tbd%ozpl, dim=2)
            errflg = 1
         end if
         if (size(Tbd%ozpl, dim=3).ne.oz_coeff) then
            write(errmsg,'(2a,i0,a,i0)') "Value error in GFS_phys_time_vary_init: ",      &
                  "oz_coeff from read_o3data does not match value in GFS_typedefs.F90: ", &
                  oz_coeff, " /= ", size(Tbd%ozpl, dim=3)
            errflg = 1
         end if
           
         call read_h2odata (Model%h2o_phys, Model%me, Model%master)
         
         ! Consistency check that the hardcoded values for levh2o and
         ! h2o_coeff in GFS_typedefs.F90 match what is set by read_o3data
         ! in GFS_typedefs.F90: allocate (Tbd%h2opl (IM,levh2o,h2o_coeff))
         if (size(Tbd%h2opl, dim=2).ne.levh2o) then
            write(errmsg,'(2a,i0,a,i0)') "Value error in GFS_phys_time_vary_init: ",     &
                  "levh2o from read_h2odata does not match value in GFS_typedefs.F90: ", &
                  levh2o, " /= ", size(Tbd%h2opl, dim=2)
            errflg = 1
         end if
         if (size(Tbd%h2opl, dim=3).ne.h2o_coeff) then
            write(errmsg,'(2a,i0,a,i0)') "Value error in GFS_phys_time_vary_init: ",       &
                  "h2o_coeff from read_h2odata does not match value in GFS_typedefs.F90: ", &
                  h2o_coeff, " /= ", size(Tbd%h2opl, dim=3)
            errflg = 1
         end if 
                       
         if (Model%iaerclm) then
           ! Consistency check that the value for ntrcaerm set in GFS_typedefs.F90
           ! and used to allocate Tbd%aer_nm matches the value defined in aerclm_def
           if (size(Tbd%aer_nm, dim=3).ne.ntrcaerm) then
              write(errmsg,'(2a,i0,a,i0)') "Value error in GFS_phys_time_vary_init: ",     &
                    "ntrcaerm from aerclm_def does not match value in GFS_typedefs.F90: ", &
                    ntrcaerm, " /= ", size(Tbd%aer_nm, dim=3)
              errflg = 1
           else
              ! Update the value of ntrcaer in aerclm_def with the value defined
              ! in GFS_typedefs.F90 that is used to allocate the Tbd DDT.
              ! If Model%iaerclm is .true., then ntrcaer == ntrcaerm
              ntrcaer = size(Tbd%aer_nm, dim=3)
              ! Read aerosol climatology
              call read_aerdata (Model%me,Model%master,Model%iflip,Model%idate,errmsg,errflg)
              if (errflg/=0) return
           endif
         else
            ! Update the value of ntrcaer in aerclm_def with the value defined
            ! in GFS_typedefs.F90 that is used to allocate the Tbd DDT.
            ! If Model%iaerclm is .false., then ntrcaer == 1
            ntrcaer = size(Tbd%aer_nm, dim=3)
         endif
         
         if (Model%iccn == 1) then
            call read_cidata  ( Model%me, Model%master)
            ! No consistency check needed for in/ccn data, all values are
            ! hardcoded in module iccn_def.F and GFS_typedefs.F90
         endif
         
         ! Update values of oz_pres in Interstitial data type for all threads
         if (Model%ntoz > 0) then
            Interstitial%oz_pres = oz_pres
         end if

         ! Update values of h2o_pres in Interstitial data type for all threads
         if (Model%h2o_phys) then
            Interstitial%h2o_pres = h2o_pres
         end if
         
         
         !--- read in and initialize ozone
         if (Model%ntoz > 0) then
            call setindxoz (Model%blksz(nb), Grid%xlat_d, Grid%jindx1_o3, &
                            Grid%jindx2_o3, Grid%ddy_o3)
         endif

         !--- read in and initialize stratospheric water
         if (Model%h2o_phys) then
            call setindxh2o (Model%blksz(nb), Grid%xlat_d, Grid%jindx1_h, &
                             Grid%jindx2_h, Grid%ddy_h)
         endif

         !--- read in and initialize aerosols
         if (Model%iaerclm) then
           call setindxaer (Model%blksz(nb), Grid%xlat_d, Grid%jindx1_aer,           &
                              Grid%jindx2_aer, Grid%ddy_aer, Grid%xlon_d,     &
                              Grid%iindx1_aer, Grid%iindx2_aer, Grid%ddx_aer, &
                              Model%me, Model%master)
         endif
          !--- read in and initialize IN and CCN
         if (Model%iccn == 1) then
             call setindxci (Model%blksz(nb), Grid%xlat_d, Grid%jindx1_ci,       &
                             Grid%jindx2_ci, Grid%ddy_ci, Grid%xlon_d,  &
                             Grid%iindx1_ci, Grid%iindx2_ci, Grid%ddx_ci)
         endif
         
        !--- initial calculation of maps local ix -> global i and j, store in Tbd
         ix = 0
         nb = 1
         do j = 1,Model%ny
           do i = 1,Model%nx
             ix = ix + 1
             if (ix .gt. Model%blksz(nb)) then
               ix = 1
               nb = nb + 1
             endif
             Tbd%jmap(ix) = j
             Tbd%imap(ix) = i
           enddo
         enddo

         is_initialized = .true.
         
      end subroutine GFS_phys_time_vary_init

!> \section arg_table_GFS_phys_time_vary_finalize Argument Table
!! \htmlinclude GFS_phys_time_vary_finalize.html
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

        ! Deallocate aerosol arrays
        if (allocated(aerin)   ) deallocate(aerin)
        if (allocated(aer_pres)) deallocate(aer_pres)

        ! Deallocate IN and CCN arrays
        if (allocated(ciplin)  ) deallocate(ciplin)
        if (allocated(ccnin)   ) deallocate(ccnin)
        if (allocated(ci_pres) ) deallocate(ci_pres)

        is_initialized = .false.

      end subroutine GFS_phys_time_vary_finalize

!> \section arg_table_GFS_phys_time_vary_run Argument Table
!! \htmlinclude GFS_phys_time_vary_run.html
!!
      subroutine GFS_phys_time_vary_run (Grid, Statein, Model, Tbd, Sfcprop, Cldprop, Diag, first_time_step, errmsg, errflg)

        use mersenne_twister,      only: random_setseed, random_number
        use machine,               only: kind_phys
        use GFS_typedefs,          only: GFS_control_type, GFS_grid_type, &
                                         GFS_Tbd_type, GFS_sfcprop_type,  &
                                         GFS_cldprop_type, GFS_diag_type, &
                                         GFS_statein_type

        implicit none

        type(GFS_grid_type),              intent(in)    :: Grid
        type(GFS_statein_type),           intent(in)    :: Statein
        type(GFS_control_type),           intent(inout) :: Model
        type(GFS_tbd_type),               intent(inout) :: Tbd
        type(GFS_sfcprop_type),           intent(inout) :: Sfcprop
        type(GFS_cldprop_type),           intent(inout) :: Cldprop
        type(GFS_diag_type),              intent(inout) :: Diag
        logical,                          intent(in)    :: first_time_step
        character(len=*),                 intent(out)   :: errmsg
        integer,                          intent(out)   :: errflg

        real(kind=kind_phys), parameter :: con_hr  = 3600.0_kind_phys
        real(kind=kind_phys), parameter :: con_99  =   99.0_kind_phys
        real(kind=kind_phys), parameter :: con_100 =  100.0_kind_phys

        integer :: i, j, k, iseed, iskip, ix, nb, kdt_rad, vegtyp
        real(kind=kind_phys) :: sec_zero, rsnow
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

        nb = 1

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
        if ( (Model%imfdeepcnv <= 0 .or. Model%cal_pre) .and. Model%random_clds ) then
          iseed = mod(con_100*sqrt(Model%fhour*con_hr),1.0d9) + Model%seed0
          call random_setseed(iseed)
          call random_number(wrk)
          do i = 1,Model%cnx*Model%nrcm
            iseed = iseed + nint(wrk(1)*1000.0) * i
            call random_setseed(iseed)
            call random_number(rannie)
            rndval(1+(i-1)*Model%cny:i*Model%cny) = rannie(1:Model%cny)
          enddo

          do k = 1,Model%nrcm
            iskip = (k-1)*Model%cnx*Model%cny
            do ix=1,Model%blksz(nb)
                j = Tbd%jmap(ix)
                i = Tbd%imap(ix)
                Tbd%rann(ix,k) = rndval(i+Model%isc-1 + (j+Model%jsc-2)*Model%cnx + iskip)
              enddo
          enddo
        endif  ! imfdeepcnv, cal_re, random_clds

        !--- o3 interpolation
        if (Model%ntoz > 0) then
          call ozinterpol (Model%me, Model%blksz(nb), Model%idate, Model%fhour, &
                           Grid%jindx1_o3, Grid%jindx2_o3, Tbd%ozpl, Grid%ddy_o3)
        endif

        !--- h2o interpolation
        if (Model%h2o_phys) then
          call h2ointerpol (Model%me, Model%blksz(nb), Model%idate, Model%fhour, &
                            Grid%jindx1_h, Grid%jindx2_h, Tbd%h2opl, Grid%ddy_h)
        endif

        !--- aerosol interpolation
        if (Model%iaerclm) then
          call aerinterpol (Model%me, Model%master, Model%blksz(nb),             &
                             Model%idate, Model%fhour,                            &
                             Grid%jindx1_aer, Grid%jindx2_aer,  &
                             Grid%ddy_aer,Grid%iindx1_aer,      &
                             Grid%iindx2_aer,Grid%ddx_aer,      &
                             Model%levs,Statein%prsl,                    &
                             Tbd%aer_nm)
        endif
         !--- ICCN interpolation
        if (Model%iccn == 1) then
            call ciinterpol (Model%me, Model%blksz(nb), Model%idate, Model%fhour, &
                             Grid%jindx1_ci, Grid%jindx2_ci,    &
                             Grid%ddy_ci,Grid%iindx1_ci,        &
                             Grid%iindx2_ci,Grid%ddx_ci,        &
                             Model%levs,Statein%prsl,                    &
                             Tbd%in_nm, Tbd%ccn_nm)
        endif

        !--- original FV3 code, not needed for SCM; also not compatible with the way
        !    the time vary steps are run (over each block) --> cannot use
        !--- repopulate specific time-varying sfc properties for AMIP/forecast runs
        !if (Model%nscyc >  0) then
        !  if (mod(Model%kdt,Model%nscyc) == 1) THEN
        !    call gcycle (nblks, Model, Grid(:), Sfcprop(:), Cldprop(:))
        !  endif
        !endif

        !--- determine if diagnostics buckets need to be cleared
        sec_zero = nint(Model%fhzero*con_hr)
        if (sec_zero >= nint(max(Model%fhswr,Model%fhlwr))) then
          if (mod(Model%kdt,Model%nszero) == 1) then
              call Diag%rad_zero  (Model)
              call Diag%phys_zero (Model)
        !!!!  THIS IS THE POINT AT WHICH DIAG%ZHOUR NEEDS TO BE UPDATED
          endif
        else
          if (mod(Model%kdt,Model%nszero) == 1) then
              call Diag%phys_zero (Model)
        !!!!  THIS IS THE POINT AT WHICH DIAG%ZHOUR NEEDS TO BE UPDATED
          endif
          kdt_rad = nint(min(Model%fhswr,Model%fhlwr)/Model%dtp)
          if (mod(Model%kdt, kdt_rad) == 1) then
              call Diag%rad_zero  (Model)
        !!!!  THIS IS THE POINT AT WHICH DIAG%ZHOUR NEEDS TO BE UPDATED
          endif
        endif

#if 0
        !Calculate sncovr if it was read in but empty (from FV3/io/FV3GFS_io.F90/sfc_prop_restart_read)
        if (first_time_step) then
          if (nint(Sfcprop%sncovr(1)) == -9999) then
            !--- compute sncovr from existing variables
            !--- code taken directly from read_fix.f
              do ix = 1, Model%blksz(nb)
                Sfcprop%sncovr(ix) = 0.0
                if (Sfcprop%slmsk(ix) > 0.001) then
                  vegtyp = Sfcprop%vtype(ix)
                  if (vegtyp == 0) vegtyp = 7
                  rsnow  = 0.001*Sfcprop%weasd(ix)/snupx(vegtyp)
                  if (0.001*Sfcprop%weasd(ix) < snupx(vegtyp)) then
                    Sfcprop%sncovr(ix) = 1.0 - (exp(-salp_data*rsnow) - rsnow*exp(-salp_data))
                  else
                    Sfcprop%sncovr(ix) = 1.0
                  endif
                endif
              enddo
              ! DH* 20201104: don't forget snocvr_ice for RUC LSM (see FV3GFS_io.F90)
          endif
        endif
#endif

      end subroutine GFS_phys_time_vary_run

    end module GFS_phys_time_vary

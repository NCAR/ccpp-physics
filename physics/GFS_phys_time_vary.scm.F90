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
      subroutine GFS_phys_time_vary_init (im, nx, ny, me, master, nblks, ntoz, iflip, &
        iccn, levh2o_int, levozp_int, idate, blksz, h2o_phys, iaerclm, xlat_d, xlon_d,&
        ozpl, h2opl, aer_nm, imap, jmap, jindx1_o3, jindx2_o3, jindx1_h, jindx2_h,    &
        jindx1_aer, jindx2_aer, iindx1_aer, iindx2_aer, jindx1_ci, jindx2_ci,         &
        iindx1_ci, iindx2_ci, ddy_o3, ddy_h, ddy_aer, ddx_aer, ddy_ci, ddx_ci,        &
        oz_pres_int, h2o_pres_int, errmsg, errflg)

         use machine,                   only: kind_phys

         implicit none

         ! Interface variables
         integer,                                intent(in) :: im, nx, ny, me, master,  &
                                                               nblks, ntoz, iflip, iccn,&
                                                               levh2o_int, levozp_int 
         integer, dimension(4),                  intent(in) :: idate
         integer, dimension(nblks),              intent(in) :: blksz
         logical,                                intent(in) :: h2o_phys, iaerclm
         real(kind=kind_phys), dimension(im),    intent(in) :: xlat_d, xlon_d
         real(kind=kind_phys), dimension(:,:,:), intent(in) :: ozpl, h2opl, aer_nm
                  
         integer, dimension(im),              intent(inout) :: imap, jmap
         integer, dimension(:),               intent(inout) :: jindx1_o3, jindx2_o3,    &
                                                               jindx1_h, jindx2_h,      &
                                                               jindx1_aer, jindx2_aer,  &
                                                               iindx1_aer, iindx2_aer,  &
                                                               jindx1_ci, jindx2_ci,    &
                                                               iindx1_ci, iindx2_ci
         real(kind=kind_phys), dimension(:),  intent(inout) :: ddy_o3, ddy_h, ddy_aer,  &
                                                               ddx_aer, ddy_ci, ddx_ci
         real(kind=kind_phys), dimension(levozp_int), intent(inout) :: oz_pres_int
         real(kind=kind_phys), dimension(levh2o_int), intent(inout) :: h2o_pres_int
         
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
         
         call read_o3data  (ntoz, me, master)

         ! Consistency check that the hardcoded values for levozp and
         ! oz_coeff in GFS_typedefs.F90 match what is set by read_o3data
         ! in GFS_typedefs.F90: allocate (Tbd%ozpl  (IM,levozp,oz_coeff))
         if (size(ozpl, dim=2).ne.levozp) then
            write(errmsg,'(2a,i0,a,i0)') "Value error in GFS_phys_time_vary_init: ",    &
                  "levozp from read_o3data does not match value in GFS_typedefs.F90: ", &
                  levozp, " /= ", size(ozpl, dim=2)
            errflg = 1
         end if
         if (size(ozpl, dim=3).ne.oz_coeff) then
            write(errmsg,'(2a,i0,a,i0)') "Value error in GFS_phys_time_vary_init: ",      &
                  "oz_coeff from read_o3data does not match value in GFS_typedefs.F90: ", &
                  oz_coeff, " /= ", size(ozpl, dim=3)
            errflg = 1
         end if
           
         call read_h2odata (h2o_phys, me, master)
         
         ! Consistency check that the hardcoded values for levh2o and
         ! h2o_coeff in GFS_typedefs.F90 match what is set by read_o3data
         ! in GFS_typedefs.F90: allocate (Tbd%h2opl (IM,levh2o,h2o_coeff))
         if (size(h2opl, dim=2).ne.levh2o) then
            write(errmsg,'(2a,i0,a,i0)') "Value error in GFS_phys_time_vary_init: ",     &
                  "levh2o from read_h2odata does not match value in GFS_typedefs.F90: ", &
                  levh2o, " /= ", size(h2opl, dim=2)
            errflg = 1
         end if
         if (size(h2opl, dim=3).ne.h2o_coeff) then
            write(errmsg,'(2a,i0,a,i0)') "Value error in GFS_phys_time_vary_init: ",       &
                  "h2o_coeff from read_h2odata does not match value in GFS_typedefs.F90: ", &
                  h2o_coeff, " /= ", size(h2opl, dim=3)
            errflg = 1
         end if 
                       
         if (iaerclm) then
           ! Consistency check that the value for ntrcaerm set in GFS_typedefs.F90
           ! and used to allocate Tbd%aer_nm matches the value defined in aerclm_def
           if (size(aer_nm, dim=3).ne.ntrcaerm) then
              write(errmsg,'(2a,i0,a,i0)') "Value error in GFS_phys_time_vary_init: ",     &
                    "ntrcaerm from aerclm_def does not match value in GFS_typedefs.F90: ", &
                    ntrcaerm, " /= ", size(aer_nm, dim=3)
              errflg = 1
           else
              ! Update the value of ntrcaer in aerclm_def with the value defined
              ! in GFS_typedefs.F90 that is used to allocate the Tbd DDT.
              ! If iaerclm is .true., then ntrcaer == ntrcaerm
              ntrcaer = size(aer_nm, dim=3)
              ! Read aerosol climatology
              call read_aerdata (me, master, iflip, idate, errmsg, errflg)
              if (errflg/=0) return
           endif
         else
            ! Update the value of ntrcaer in aerclm_def with the value defined
            ! in GFS_typedefs.F90 that is used to allocate the Tbd DDT.
            ! If iaerclm is .false., then ntrcaer == 1
            ntrcaer = size(aer_nm, dim=3)
         endif
         
         if (iccn == 1) then
            call read_cidata  (me, master)
            ! No consistency check needed for in/ccn data, all values are
            ! hardcoded in module iccn_def.F and GFS_typedefs.F90
         endif
         
         ! Update values of oz_pres in Interstitial data type for all threads
         if (ntoz > 0) then
            oz_pres_int = oz_pres
         end if

         ! Update values of h2o_pres in Interstitial data type for all threads
         if (h2o_phys) then
            h2o_pres_int = h2o_pres
         end if
         
         
         !--- read in and initialize ozone
         if (ntoz > 0) then
            call setindxoz (blksz(nb), xlat_d, jindx1_o3, &
                            jindx2_o3, ddy_o3)
         endif

         !--- read in and initialize stratospheric water
         if (h2o_phys) then
            call setindxh2o (blksz(nb), xlat_d, jindx1_h, &
                             jindx2_h, ddy_h)
         endif

         !--- read in and initialize aerosols
         if (iaerclm) then
           call setindxaer (blksz(nb), xlat_d, jindx1_aer,           &
                              jindx2_aer, ddy_aer, xlon_d,     &
                              iindx1_aer, iindx2_aer, ddx_aer, &
                              me, master)
         endif
          !--- read in and initialize IN and CCN
         if (iccn == 1) then
             call setindxci (blksz(nb), xlat_d, jindx1_ci,       &
                             jindx2_ci, ddy_ci, xlon_d,  &
                             iindx1_ci, iindx2_ci, ddx_ci)
         endif
         
        !--- initial calculation of maps local ix -> global i and j, store in Tbd
         ix = 0
         nb = 1
         do j = 1, ny
           do i = 1, nx
             ix = ix + 1
             if (ix .gt. blksz(nb)) then
               ix = 1
               nb = nb + 1
             endif
             jmap(ix) = j
             imap(ix) = i
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
      subroutine GFS_phys_time_vary_run (levs, cnx, cny, isc, jsc, me, master,   &
        ntoz, iccn, nrcm, nsswr, nszero, kdt, imfdeepcnv, seed0, first_time_step,&
        lsswr, cal_pre, random_clds, h2o_phys, iaerclm, fhswr, fhlwr, fhour,     &
        fhzero, dtp, idate, jindx1_o3, jindx2_o3, jindx1_h, jindx2_h, jindx1_aer,&
        jindx2_aer, iindx1_aer, iindx2_aer, jindx1_ci, jindx2_ci, iindx1_ci,     &
        iindx2_ci, blksz, imap, jmap, ddy_o3, ddy_h, ddy_aer, ddx_aer, ddy_ci,   &
        ddx_ci, slmsk, vtype, weasd, prsl, Model, clstp, sncovr, rann, in_nm,    &
        ccn_nm, ozpl, h2opl, aer_nm, Diag, errmsg, errflg)

        use mersenne_twister,      only: random_setseed, random_number
        use machine,               only: kind_phys
        use GFS_typedefs,          only: GFS_control_type, GFS_diag_type
        
        implicit none

        integer,                              intent(in) :: levs, cnx, cny, isc, jsc, &
                                                            me, master, ntoz, iccn,   &
                                                            nrcm, nsswr, nszero, kdt, &
                                                            imfdeepcnv, seed0 
        logical,                              intent(in) :: first_time_step, lsswr,   &
                                                            cal_pre, random_clds,     &
                                                            h2o_phys, iaerclm
        real(kind=kind_phys),                 intent(in) :: fhswr, fhlwr, fhour,      &
                                                            fhzero, dtp
        
        integer, dimension(4),                intent(in) :: idate
        integer, dimension(:),                intent(in) :: jindx1_o3, jindx2_o3,     &
                                                            jindx1_h, jindx2_h,       &
                                                            jindx1_aer, jindx2_aer,   &
                                                            iindx1_aer, iindx2_aer,   &
                                                            jindx1_ci, jindx2_ci,     &
                                                            iindx1_ci, iindx2_ci,     &
                                                            blksz, imap, jmap
        real(kind=kind_phys), dimension(:),   intent(in) :: ddy_o3, ddy_h, ddy_aer,   &
                                                            ddx_aer, ddy_ci, ddx_ci,  &
                                                            slmsk, vtype, weasd
        real(kind=kind_phys), dimension(:,:), intent(in) :: prsl
        
        type(GFS_control_type),               intent(in) :: Model
        
        real(kind=kind_phys),                   intent(inout) :: clstp
        real(kind=kind_phys), dimension(:),     intent(inout) :: sncovr
        real(kind=kind_phys), dimension(:,:),   intent(inout) :: rann, in_nm, ccn_nm
        real(kind=kind_phys), dimension(:,:,:), intent(inout) :: ozpl, h2opl, aer_nm
        
        type(GFS_diag_type),              intent(inout) :: Diag

        character(len=*),                 intent(out)   :: errmsg
        integer,                          intent(out)   :: errflg

        real(kind=kind_phys), parameter :: con_hr  = 3600.0_kind_phys
        real(kind=kind_phys), parameter :: con_99  =   99.0_kind_phys
        real(kind=kind_phys), parameter :: con_100 =  100.0_kind_phys

        integer :: i, j, k, iseed, iskip, ix, nb, kdt_rad, vegtyp
        real(kind=kind_phys) :: sec_zero, rsnow
        real(kind=kind_phys) :: wrk(1)
        real(kind=kind_phys) :: rannie(cny)
        real(kind=kind_phys) :: rndval(cnx*cny*nrcm)

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
        if ((mod(kdt,nsswr) == 0) .and. (lsswr)) then
          !--- initialize,accumulate,convert
          clstp = 1100 + min(fhswr/con_hr,fhour,con_99)
        elseif (mod(kdt,nsswr) == 0) then
          !--- accumulate,convert
          clstp = 0100 + min(fhswr/con_hr,fhour,con_99)
        elseif (lsswr) then
          !--- initialize,accumulate
          clstp = 1100
        else
          !--- accumulate
          clstp = 0100
        endif

        !--- random number needed for RAS and old SAS and when cal_pre=.true.
        if ( (imfdeepcnv <= 0 .or. cal_pre) .and. random_clds ) then
          iseed = mod(con_100*sqrt(fhour*con_hr),1.0d9) + seed0
          call random_setseed(iseed)
          call random_number(wrk)
          do i = 1,cnx*nrcm
            iseed = iseed + nint(wrk(1)*1000.0) * i
            call random_setseed(iseed)
            call random_number(rannie)
            rndval(1+(i-1)*cny:i*cny) = rannie(1:cny)
          enddo

          do k = 1,nrcm
            iskip = (k-1)*cnx*cny
            do ix=1,blksz(nb)
                j = jmap(ix)
                i = imap(ix)
                rann(ix,k) = rndval(i+isc-1 + (j+jsc-2)*cnx + iskip)
              enddo
          enddo
        endif  ! imfdeepcnv, cal_re, random_clds

        !--- o3 interpolation
        if (ntoz > 0) then
          call ozinterpol (me, blksz(nb), idate, fhour, &
                           jindx1_o3, jindx2_o3, ozpl, ddy_o3)
        endif

        !--- h2o interpolation
        if (h2o_phys) then
          call h2ointerpol (me, blksz(nb), idate, fhour, &
                            jindx1_h, jindx2_h, h2opl, ddy_h)
        endif

        !--- aerosol interpolation
        if (iaerclm) then
          call aerinterpol (me, master, blksz(nb),             &
                             idate, fhour,                            &
                             jindx1_aer, jindx2_aer,  &
                             ddy_aer,iindx1_aer,      &
                             iindx2_aer,ddx_aer,      &
                             levs,prsl,                    &
                             aer_nm)
        endif
         !--- ICCN interpolation
        if (iccn == 1) then
            call ciinterpol (me, blksz(nb), idate, fhour, &
                             jindx1_ci, jindx2_ci,    &
                             ddy_ci,iindx1_ci,        &
                             iindx2_ci,ddx_ci,        &
                             levs,prsl,                    &
                             in_nm, ccn_nm)
        endif

        !--- original FV3 code, not needed for SCM; also not compatible with the way
        !    the time vary steps are run (over each block) --> cannot use
        !--- repopulate specific time-varying sfc properties for AMIP/forecast runs
        !if (Model%nscyc >  0) then
        !  if (mod(kdt,Model%nscyc) == 1) THEN
        !    call gcycle (nblks, Model, Grid(:), Sfcprop(:), Cldprop(:))
        !  endif
        !endif

        !--- determine if diagnostics buckets need to be cleared
        sec_zero = nint(fhzero*con_hr)
        if (sec_zero >= nint(max(fhswr,fhlwr))) then
          if (mod(kdt,nszero) == 1) then
              call Diag%rad_zero  (Model)
              call Diag%phys_zero (Model)
        !!!!  THIS IS THE POINT AT WHICH DIAG%ZHOUR NEEDS TO BE UPDATED
          endif
        else
          if (mod(kdt,nszero) == 1) then
              call Diag%phys_zero (Model)
        !!!!  THIS IS THE POINT AT WHICH DIAG%ZHOUR NEEDS TO BE UPDATED
          endif
          kdt_rad = nint(min(fhswr,fhlwr)/dtp)
          if (mod(kdt, kdt_rad) == 1) then
              call Diag%rad_zero  (Model)
        !!!!  THIS IS THE POINT AT WHICH DIAG%ZHOUR NEEDS TO BE UPDATED
          endif
        endif

#if 0
        !Calculate sncovr if it was read in but empty (from FV3/io/FV3GFS_io.F90/sfc_prop_restart_read)
        if (first_time_step) then
          if (nint(sncovr(1)) == -9999) then
            !--- compute sncovr from existing variables
            !--- code taken directly from read_fix.f
              do ix = 1, blksz(nb)
                sncovr(ix) = 0.0
                if (slmsk(ix) > 0.001) then
                  vegtyp = vtype(ix)
                  if (vegtyp == 0) vegtyp = 7
                  rsnow  = 0.001*weasd(ix)/snupx(vegtyp)
                  if (0.001*weasd(ix) < snupx(vegtyp)) then
                    sncovr(ix) = 1.0 - (exp(-salp_data*rsnow) - rsnow*exp(-salp_data))
                  else
                    sncovr(ix) = 1.0
                  endif
                endif
              enddo
          endif
        endif
#endif

      end subroutine GFS_phys_time_vary_run

    end module GFS_phys_time_vary

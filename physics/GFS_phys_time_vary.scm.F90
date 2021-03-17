!> \file GFS_phys_time_vary.scm.F90
!!  Contains code related to GFS physics suite setup (physics part of time_vary_step)

!>\defgroup mod_GFS_phys_time_vary GFS Physics Time Update
!! This module contains GFS physics time vary subroutines including ozone, stratospheric water vapor, 
!! aerosol, IN&CCN and surface properties updates. 
!> @{
   module GFS_phys_time_vary
     
      use machine, only : kind_phys

      use mersenne_twister, only: random_setseed, random_number

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

      public GFS_phys_time_vary_init, GFS_phys_time_vary_timestep_init, GFS_phys_time_vary_timestep_finalize, GFS_phys_time_vary_finalize

      logical :: is_initialized = .false.

      real(kind=kind_phys), parameter :: con_hr  = 3600.0_kind_phys
      real(kind=kind_phys), parameter :: con_99  =   99.0_kind_phys
      real(kind=kind_phys), parameter :: con_100 =  100.0_kind_phys

      contains

!> \section arg_table_GFS_phys_time_vary_init Argument Table
!! \htmlinclude GFS_phys_time_vary_init.html
!!
!>\section gen_GFS_phys_time_vary_init GFS_phys_time_vary_init General Algorithm
!! @{
      subroutine GFS_phys_time_vary_init (                                                         &
              me, master, ntoz, h2o_phys, iaerclm, iccn, iflip, im, nx, ny, idate, xlat_d, xlon_d, &
              jindx1_o3, jindx2_o3, ddy_o3, ozpl, jindx1_h, jindx2_h, ddy_h, h2opl,                &
              jindx1_aer, jindx2_aer, ddy_aer, iindx1_aer, iindx2_aer, ddx_aer, aer_nm,            &
              jindx1_ci, jindx2_ci, ddy_ci, iindx1_ci, iindx2_ci, ddx_ci, imap, jmap,              &
              errmsg, errflg)

         implicit none

         ! Interface variables
         integer,              intent(in)    :: me, master, ntoz, iccn, iflip, im, nx, ny
         logical,              intent(in)    :: h2o_phys, iaerclm
         integer,              intent(in)    :: idate(:)
         real(kind_phys),      intent(in)    :: xlat_d(:), xlon_d(:)

         integer,              intent(inout) :: jindx1_o3(:), jindx2_o3(:), jindx1_h(:), jindx2_h(:)
         real(kind_phys),      intent(inout) :: ddy_o3(:),  ddy_h(:)
         real(kind_phys),      intent(in)    :: ozpl(:,:,:), h2opl(:,:,:)
         integer,              intent(inout) :: jindx1_aer(:), jindx2_aer(:), iindx1_aer(:), iindx2_aer(:)
         real(kind_phys),      intent(inout) :: ddy_aer(:), ddx_aer(:)
         real(kind_phys),      intent(in)    :: aer_nm(:,:,:)
         integer,              intent(inout) :: jindx1_ci(:), jindx2_ci(:), iindx1_ci(:), iindx2_ci(:)
         real(kind_phys),      intent(inout) :: ddy_ci(:), ddx_ci(:)
         integer,              intent(inout) :: imap(:), jmap(:)

         character(len=*),     intent(out)   :: errmsg
         integer,              intent(out)   :: errflg

         ! Local variables
         integer :: i, j, ix

         ! Initialize CCPP error handling variables
         errmsg = ''
         errflg = 0

         if (is_initialized) return
         
!> - Call read_o3data() to read ozone data 
         call read_o3data (ntoz, me, master)

         ! Consistency check that the hardcoded values for levozp and
         ! oz_coeff in GFS_typedefs.F90 match what is set by read_o3data
         ! in GFS_typedefs.F90: allocate (Tbd%ozpl (IM,levozp,oz_coeff))
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

!> - Call read_h2odata() to read stratospheric water vapor data
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

!> - Call read_aerdata() to read aerosol climatology
         if (iaerclm) then
            ! Consistency check that the value for ntrcaerm set in GFS_typedefs.F90
            ! and used to allocate aer_nm matches the value defined in aerclm_def
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
               call read_aerdata (me,master,iflip,idate,errmsg,errflg)
            endif
         else
            ! Update the value of ntrcaer in aerclm_def with the value defined
            ! in GFS_typedefs.F90 that is used to allocate the Tbd DDT.
            ! If iaerclm is .false., then ntrcaer == 1
            ntrcaer = size(aer_nm, dim=3)
         endif

!> - Call read_cidata() to read IN and CCN data
         if (iccn == 1) then
           call read_cidata (me,master)
           ! No consistency check needed for in/ccn data, all values are
           ! hardcoded in module iccn_def.F and GFS_typedefs.F90
         endif

!> - Call setindxoz() to initialize ozone data
         if (ntoz > 0) then
           call setindxoz (im, xlat_d, jindx1_o3, jindx2_o3, ddy_o3)
         endif

!> - Call setindxh2o() to initialize stratospheric water vapor data
         if (h2o_phys) then
           call setindxh2o (im, xlat_d, jindx1_h, jindx2_h, ddy_h)
         endif

!> - Call setindxaer() to initialize aerosols data
         if (iaerclm) then
           call setindxaer (im, xlat_d, jindx1_aer,          &
                            jindx2_aer, ddy_aer, xlon_d,     &
                            iindx1_aer, iindx2_aer, ddx_aer, &
                            me, master)
         endif

!> - Call setindxci() to initialize IN and CCN data
         if (iccn == 1) then
           call setindxci (im, xlat_d, jindx1_ci,      &
                           jindx2_ci, ddy_ci, xlon_d,  &
                           iindx1_ci, iindx2_ci, ddx_ci)
         endif

         !--- initial calculation of maps local ix -> global i and j
         ix = 0
         do j = 1,ny
           do i = 1,nx
             ix = ix + 1
             jmap(ix) = j
             imap(ix) = i
           enddo
         enddo

#if 0
        !Calculate sncovr if it was read in but empty (from FV3/io/FV3GFS_io.F90/sfc_prop_restart_read)
        ! if (first_time_step) then
        !   if (nint(Sfcprop%sncovr(1)) == -9999) then
        !     !--- compute sncovr from existing variables
        !     !--- code taken directly from read_fix.f
        !       do ix = 1, im
        !         Sfcprop%sncovr(ix) = 0.0
        !         if (Sfcprop%slmsk(ix) > 0.001) then
        !           vegtyp = Sfcprop%vtype(ix)
        !           if (vegtyp == 0) vegtyp = 7
        !           rsnow  = 0.001*Sfcprop%weasd(ix)/snupx(vegtyp)
        !           if (0.001*Sfcprop%weasd(ix) < snupx(vegtyp)) then
        !             Sfcprop%sncovr(ix) = 1.0 - (exp(-salp_data*rsnow) - rsnow*exp(-salp_data))
        !           else
        !             Sfcprop%sncovr(ix) = 1.0
        !           endif
        !         endif
        !       enddo
        !       ! DH* 20201104: don't forget snocvr_ice for RUC LSM (see FV3GFS_io.F90)
        !   endif
        ! endif
#endif

         is_initialized = .true.

      end subroutine GFS_phys_time_vary_init
!! @}

!> \section arg_table_GFS_phys_time_vary_timestep_init Argument Table
!! \htmlinclude GFS_phys_time_vary_timestep_init.html
!!
!>\section gen_GFS_phys_time_vary_timestep_init GFS_phys_time_vary_timestep_init General Algorithm
!! @{
      subroutine GFS_phys_time_vary_timestep_init (                                                 &
            me, master, cnx, cny, isc, jsc, nrcm, im, levs, kdt, idate, nsswr, fhswr, lsswr, fhour, &
            imfdeepcnv, cal_pre, random_clds,        ntoz, h2o_phys, iaerclm, iccn, clstp,          &
            jindx1_o3, jindx2_o3, ddy_o3, ozpl, jindx1_h, jindx2_h, ddy_h, h2opl,                   &
            jindx1_aer, jindx2_aer, ddy_aer, iindx1_aer, iindx2_aer, ddx_aer, aer_nm,               &
            jindx1_ci, jindx2_ci, ddy_ci, iindx1_ci, iindx2_ci, ddx_ci, in_nm, ccn_nm,              &
            imap, jmap, prsl, seed0, rann, errmsg, errflg)

         implicit none

         ! Interface variables
         integer,              intent(in)    :: me, master, cnx, cny, isc, jsc, nrcm, im, levs, kdt, &
                                                nsswr, imfdeepcnv, iccn, ntoz
         integer,              intent(in)    :: idate(:)
         real(kind_phys),      intent(in)    :: fhswr, fhour
         logical,              intent(in)    :: lsswr, cal_pre, random_clds, h2o_phys, iaerclm
         real(kind_phys),      intent(out)   :: clstp
         integer,              intent(in)    :: jindx1_o3(:), jindx2_o3(:), jindx1_h(:), jindx2_h(:)
         real(kind_phys),      intent(in)    :: ddy_o3(:),  ddy_h(:)
         real(kind_phys),      intent(inout) :: ozpl(:,:,:), h2opl(:,:,:)
         integer,              intent(in)    :: jindx1_aer(:), jindx2_aer(:), iindx1_aer(:), iindx2_aer(:)
         real(kind_phys),      intent(in)    :: ddy_aer(:), ddx_aer(:)
         real(kind_phys),      intent(inout) :: aer_nm(:,:,:)
         integer,              intent(in)    :: jindx1_ci(:), jindx2_ci(:), iindx1_ci(:), iindx2_ci(:)
         real(kind_phys),      intent(in)    :: ddy_ci(:), ddx_ci(:)
         real(kind_phys),      intent(inout) :: in_nm(:,:), ccn_nm(:,:)
         integer,              intent(in)    :: imap(:), jmap(:)
         real(kind_phys),      intent(in)    :: prsl(:,:)
         integer,              intent(in)    :: seed0
         real(kind_phys),      intent(inout) :: rann(:,:)
         !
         character(len=*),     intent(out)   :: errmsg
         integer,              intent(out)   :: errflg

         ! Local variables
         integer :: i, j, k, iseed, iskip, ix, kdt_rad
         real(kind=kind_phys) :: sec_zero, rsnow
         real(kind=kind_phys) :: wrk(1)
         real(kind=kind_phys) :: rannie(cny)
         real(kind=kind_phys) :: rndval(cnx*cny*nrcm)

         ! Initialize CCPP error handling variables
         errmsg = ''
         errflg = 0

         ! Check initialization status
         if (.not.is_initialized) then
            write(errmsg,'(*(a))') "Logic error: GFS_phys_time_vary_timestep_init called before GFS_phys_time_vary_init"
            errflg = 1
            return
         end if

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
         !    imfdeepcnv < 0 when ras = .true.
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
            do ix=1,im
              j = jmap(ix)
              i = imap(ix)
              rann(ix,k) = rndval(i+isc-1 + (j+jsc-2)*cnx + iskip)
            enddo
          enddo

         endif  ! imfdeepcnv, cal_re, random_clds

!> - Call ozinterpol() to make ozone interpolation
         if (ntoz > 0) then
           call ozinterpol (me, im, idate, fhour, &
                            jindx1_o3, jindx2_o3, &
                            ozpl, ddy_o3)
         endif

!> - Call h2ointerpol() to make stratospheric water vapor data interpolation
         if (h2o_phys) then
           call h2ointerpol (me, im, idate, fhour, &
                             jindx1_h, jindx2_h,   &
                             h2opl, ddy_h)
         endif

!> - Call aerinterpol() to make aerosol interpolation
         if (iaerclm) then
           call aerinterpol (me, master, im, idate, fhour, &
                             jindx1_aer, jindx2_aer,       &
                             ddy_aer, iindx1_aer,          &
                             iindx2_aer, ddx_aer,          &
                             levs, prsl, aer_nm)
         endif

!> - Call ciinterpol() to make IN and CCN data interpolation
         if (iccn == 1) then
           call ciinterpol (me, im, idate, fhour,    &
                            jindx1_ci, jindx2_ci,    &
                            ddy_ci, iindx1_ci,       &
                            iindx2_ci, ddx_ci,       &
                            levs, prsl, in_nm, ccn_nm)
         endif

!       Not needed for SCM:
!> - Call gcycle() to repopulate specific time-varying surface properties for AMIP/forecast runs
        !if (nscyc >  0) then
        !   if (mod(kdt,nscyc) == 1) THEN
        !     call gcycle (me, nthrds, nx, ny, isc, jsc, nsst, tile_num, nlunit,       &
        !         input_nml_file, lsoil, lsoil_lsm, kice, idate, ialb, isot, ivegsrc,  &
        !         use_ufo, nst_anl, fhcyc, phour, lakefrac, min_seaice, min_lakeice,   &
        !         frac_grid, smc, slc, stc, smois, sh2o, tslb, tiice, tg3, tref, tsfc, &
        !         tsfco, tisfc, hice, fice, facsf, facwf, alvsf, alvwf, alnsf, alnwf,  &
        !         zorli, zorll, zorlo, weasd, slope, snoalb, canopy, vfrac, vtype,     &
        !         stype, shdmin, shdmax, snowd, cv, cvb, cvt, oro, oro_uf,             &
        !         xlat_d, xlon_d, slmsk, imap, jmap)
        !   endif
        !endif

      end subroutine GFS_phys_time_vary_timestep_init
!! @}

!> \section arg_table_GFS_phys_time_vary_timestep_finalize Argument Table
!! \htmlinclude GFS_phys_time_vary_timestep_finalize.html
!!
!>\section gen_GFS_phys_time_vary_timestep_finalize GFS_phys_time_vary_timestep_finalize General Algorithm
!! @{
      subroutine GFS_phys_time_vary_timestep_finalize (errmsg, errflg)

         implicit none

         ! Interface variables
         character(len=*),                 intent(out)   :: errmsg
         integer,                          intent(out)   :: errflg

         ! Initialize CCPP error handling variables
         errmsg = ''
         errflg = 0

      end subroutine GFS_phys_time_vary_timestep_finalize
!! @}

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

   end module GFS_phys_time_vary
!> @}

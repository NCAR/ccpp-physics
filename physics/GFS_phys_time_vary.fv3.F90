!> \file GFS_phys_time_vary.fv3.F90
!!  Contains code related to GFS physics suite setup (physics part of time_vary_step)

!>\defgroup mod_GFS_phys_time_vary GFS Physics Time Update
!! This module contains GFS physics time vary subroutines including ozone, stratospheric water vapor,
!! aerosol, IN&CCN and surface properties updates.
!> @{
   module GFS_phys_time_vary

#ifdef OPENMP
      use omp_lib
#endif

      use machine, only : kind_phys

      use mersenne_twister, only: random_setseed, random_number

      use ozne_def, only : levozp, oz_coeff, oz_lat, oz_pres, oz_time, ozplin
      use ozinterp, only : read_o3data, setindxoz, ozinterpol

      use h2o_def,   only : levh2o, h2o_coeff, h2o_lat, h2o_pres, h2o_time, h2oplin
      use h2ointerp, only : read_h2odata, setindxh2o, h2ointerpol

      use aerclm_def, only : aerin, aer_pres, ntrcaer, ntrcaerm
      use aerinterp,  only : read_aerdata, setindxaer, aerinterpol, read_aerdataf

      use iccn_def,   only : ciplin, ccnin, ci_pres
      use iccninterp, only : read_cidata, setindxci, ciinterpol

      use gcycle_mod, only : gcycle

      use cires_tauamf_data,   only:  cires_indx_ugwp,  read_tau_amf, tau_amf_interp
      use cires_tauamf_data,   only:  tau_limb,  days_limb, ugwp_taulat

      !--- variables needed for calculating 'sncovr'
      use namelist_soilveg, only: salp_data, snupx
      use set_soilveg_mod, only: set_soilveg

      ! --- needed for Noah MP init
      use noahmp_tables, only: laim_table,saim_table,sla_table,      &
                               bexp_table,smcmax_table,smcwlt_table, &
                               dwsat_table,dksat_table,psisat_table, &
                               isurban_table,isbarren_table,         &
                               isice_table,iswater_table

      implicit none

      private

      public GFS_phys_time_vary_init, GFS_phys_time_vary_timestep_init, GFS_phys_time_vary_timestep_finalize, GFS_phys_time_vary_finalize

      logical :: is_initialized = .false.

      real(kind=kind_phys), parameter :: con_hr        =  3600.0_kind_phys
      real(kind=kind_phys), parameter :: con_99        =    99.0_kind_phys
      real(kind=kind_phys), parameter :: con_100       =   100.0_kind_phys
      real(kind=kind_phys), parameter :: missing_value = 9.99e20_kind_phys
      real(kind=kind_phys), parameter :: drythresh     =   1.e-4_kind_phys
      real(kind=kind_phys), parameter :: zero          =     0.0_kind_phys
      real(kind=kind_phys), parameter :: one           =     1.0_kind_phys

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
              do_ugwp_v1, jindx1_tau, jindx2_tau, ddy_j1tau, ddy_j2tau,                            &
              isot, ivegsrc, nlunit, sncovr, sncovr_ice, lsm, lsm_noahmp, lsm_ruc, min_seaice,     &
              fice, landfrac, vtype, weasd, lsoil, zs, dzs, lsnow_lsm_lbound, lsnow_lsm_ubound,    &
              tvxy, tgxy, tahxy, canicexy, canliqxy, eahxy, cmxy, chxy, fwetxy, sneqvoxy, alboldxy,&
              qsnowxy, wslakexy, albdvis, albdnir, albivis, albinir, emiss, taussxy, waxy, wtxy,   &
              zwtxy, xlaixy, xsaixy, lfmassxy, stmassxy, rtmassxy, woodxy, stblcpxy, fastcpxy,     &
              smcwtdxy, deeprechxy, rechxy, snowxy, snicexy, snliqxy, tsnoxy , smoiseq, zsnsoxy,   &
              slc, smc, stc, tsfcl, snowd, canopy, tg3, stype, con_t0c, nthrds, errmsg, errflg)

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
         logical,              intent(in)    :: do_ugwp_v1
         real(kind_phys),      intent(inout) :: ddy_j1tau(:), ddy_j2tau(:)
         integer,              intent(inout) :: jindx1_tau(:), jindx2_tau(:)

         integer,              intent(in)    :: isot, ivegsrc, nlunit
         real(kind_phys),      intent(inout) :: sncovr(:), sncovr_ice(:)
         integer,              intent(in)    :: lsm, lsm_noahmp, lsm_ruc
         real(kind_phys),      intent(in)    :: min_seaice, fice(:)
         real(kind_phys),      intent(in)    :: landfrac(:), vtype(:)
         real(kind_phys),      intent(inout) :: weasd(:)

         ! NoahMP - only allocated when NoahMP is used
         integer, intent(in) :: lsoil, lsnow_lsm_lbound, lsnow_lsm_ubound
         real(kind_phys),      intent(in)    :: zs(:)
         real(kind_phys),      intent(in)    :: dzs(:)
         real(kind_phys),      intent(inout) :: tvxy(:)
         real(kind_phys),      intent(inout) :: tgxy(:)
         real(kind_phys),      intent(inout) :: tahxy(:)
         real(kind_phys),      intent(inout) :: canicexy(:)
         real(kind_phys),      intent(inout) :: canliqxy(:)
         real(kind_phys),      intent(inout) :: eahxy(:)
         real(kind_phys),      intent(inout) :: cmxy(:)
         real(kind_phys),      intent(inout) :: chxy(:)
         real(kind_phys),      intent(inout) :: fwetxy(:)
         real(kind_phys),      intent(inout) :: sneqvoxy(:)
         real(kind_phys),      intent(inout) :: alboldxy(:)
         real(kind_phys),      intent(inout) :: qsnowxy(:)
         real(kind_phys),      intent(inout) :: wslakexy(:)
         real(kind_phys),      intent(inout) :: albdvis(:)
         real(kind_phys),      intent(inout) :: albdnir(:)
         real(kind_phys),      intent(inout) :: albivis(:)
         real(kind_phys),      intent(inout) :: albinir(:)
         real(kind_phys),      intent(inout) :: emiss(:)
         real(kind_phys),      intent(inout) :: taussxy(:)
         real(kind_phys),      intent(inout) :: waxy(:)
         real(kind_phys),      intent(inout) :: wtxy(:)
         real(kind_phys),      intent(inout) :: zwtxy(:)
         real(kind_phys),      intent(inout) :: xlaixy(:)
         real(kind_phys),      intent(inout) :: xsaixy(:)
         real(kind_phys),      intent(inout) :: lfmassxy(:)
         real(kind_phys),      intent(inout) :: stmassxy(:)
         real(kind_phys),      intent(inout) :: rtmassxy(:)
         real(kind_phys),      intent(inout) :: woodxy(:)
         real(kind_phys),      intent(inout) :: stblcpxy(:)
         real(kind_phys),      intent(inout) :: fastcpxy(:)
         real(kind_phys),      intent(inout) :: smcwtdxy(:)
         real(kind_phys),      intent(inout) :: deeprechxy(:)
         real(kind_phys),      intent(inout) :: rechxy(:)
         real(kind_phys),      intent(inout) :: snowxy(:)
         real(kind_phys),      intent(inout) :: snicexy(:,lsnow_lsm_lbound:)
         real(kind_phys),      intent(inout) :: snliqxy(:,lsnow_lsm_lbound:)
         real(kind_phys),      intent(inout) :: tsnoxy (:,lsnow_lsm_lbound:)
         real(kind_phys),      intent(inout) :: smoiseq(:,:)
         real(kind_phys),      intent(inout) :: zsnsoxy(:,lsnow_lsm_lbound:)
         real(kind_phys),      intent(inout) :: slc(:,:)
         real(kind_phys),      intent(inout) :: smc(:,:)
         real(kind_phys),      intent(inout) :: stc(:,:)
         real(kind_phys),      intent(in)    :: tsfcl(:)
         real(kind_phys),      intent(in)    :: snowd(:)
         real(kind_phys),      intent(in)    :: canopy(:)
         real(kind_phys),      intent(in)    :: tg3(:)
         real(kind_phys),      intent(in)    :: stype(:)
         real(kind_phys),      intent(in)    :: con_t0c

         integer,              intent(in)    :: nthrds
         character(len=*),     intent(out)   :: errmsg
         integer,              intent(out)   :: errflg

         ! Local variables
         integer :: i, j, ix, vegtyp, iamin, iamax, jamin, jamax
         real(kind_phys) :: rsnow

         !--- Noah MP
         integer              :: soiltyp, isnow, is, imn
         real(kind=kind_phys) :: masslai, masssai, snd
         real(kind=kind_phys) :: bexp, ddz, smcmax, smcwlt, dwsat, dksat, psisat

         real(kind=kind_phys), dimension(:), allocatable :: dzsno
         real(kind=kind_phys), dimension(:), allocatable :: dzsnso

         ! Initialize CCPP error handling variables
         errmsg = ''
         errflg = 0

         if (is_initialized) return
         iamin=999
         iamax=-999
         jamin=999
         jamax=-999

!$OMP parallel num_threads(nthrds) default(none)                                    &
!$OMP          shared (me,master,ntoz,h2o_phys,im,nx,ny,idate)                      &
!$OMP          shared (xlat_d,xlon_d,imap,jmap,errmsg,errflg)                       &
!$OMP          shared (levozp,oz_coeff,oz_pres,ozpl)                                &
!$OMP          shared (levh2o,h2o_coeff,h2o_pres,h2opl)                             &
!$OMP          shared (iamin, iamax, jamin, jamax)                                  &
!$OMP          shared (iaerclm,ntrcaer,aer_nm,iflip,iccn)                           &
!$OMP          shared (jindx1_o3,jindx2_o3,ddy_o3,jindx1_h,jindx2_h,ddy_h)          &
!$OMP          shared (jindx1_aer,jindx2_aer,ddy_aer,iindx1_aer,iindx2_aer,ddx_aer) &
!$OMP          shared (jindx1_ci,jindx2_ci,ddy_ci,iindx1_ci,iindx2_ci,ddx_ci)       &
!$OMP          shared (do_ugwp_v1,jindx1_tau,jindx2_tau,ddy_j1tau,ddy_j2tau)        &
!$OMP          shared (isot,ivegsrc,nlunit,sncovr,sncovr_ice,lsm,lsm_ruc)           &
!$OMP          shared (min_seaice,fice,landfrac,vtype,weasd,snupx,salp_data)        &
!$OMP          private (ix,i,j,rsnow,vegtyp)

!$OMP sections

!$OMP section
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

!$OMP section
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

!$OMP section
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

!$OMP section
!> - Call read_cidata() to read IN and CCN data
         if (iccn == 1) then
           call read_cidata (me,master)
           ! No consistency check needed for in/ccn data, all values are
           ! hardcoded in module iccn_def.F and GFS_typedefs.F90
         endif

!$OMP section
!> - Call tau_amf dats for  ugwp_v1
         if (do_ugwp_v1) then
            call read_tau_amf(me, master, errmsg, errflg)
         endif

!$OMP section
!> - Initialize soil vegetation (needed for sncovr calculation further down)
         call set_soilveg(me, isot, ivegsrc, nlunit)

!$OMP end sections

! Need an OpenMP barrier here (implicit in "end sections")

!$OMP sections

!$OMP section
!> - Call setindxoz() to initialize ozone data
         if (ntoz > 0) then
           call setindxoz (im, xlat_d, jindx1_o3, jindx2_o3, ddy_o3)
         endif

!$OMP section
!> - Call setindxh2o() to initialize stratospheric water vapor data
         if (h2o_phys) then
           call setindxh2o (im, xlat_d, jindx1_h, jindx2_h, ddy_h)
         endif

!> - Call setindxaer() to initialize aerosols data
!$OMP section
         if (iaerclm) then
           call setindxaer (im, xlat_d, jindx1_aer,          &
                            jindx2_aer, ddy_aer, xlon_d,     &
                            iindx1_aer, iindx2_aer, ddx_aer, &
                            me, master)
           iamin=min(minval(iindx1_aer), iamin)
           iamax=max(maxval(iindx2_aer), iamax)
           jamin=min(minval(jindx1_aer), jamin)
           jamax=max(maxval(jindx2_aer), jamax)
         endif
!$OMP section

!> - Call setindxci() to initialize IN and CCN data
         if (iccn == 1) then
           call setindxci (im, xlat_d, jindx1_ci,      &
                           jindx2_ci, ddy_ci, xlon_d,  &
                           iindx1_ci, iindx2_ci, ddx_ci)
         endif

!$OMP section
!> - Call  cires_indx_ugwp to read monthly-mean GW-tau diagnosed from FV3GFS-runs that can resolve GWs
         if (do_ugwp_v1) then
            call cires_indx_ugwp (im, me, master, xlat_d, jindx1_tau, jindx2_tau,  &
                                  ddy_j1tau, ddy_j2tau)
         endif

!$OMP section
         !--- initial calculation of maps local ix -> global i and j
         ix = 0
         do j = 1,ny
           do i = 1,nx
             ix = ix + 1
             jmap(ix) = j
             imap(ix) = i
           enddo
         enddo

!$OMP section
         !--- if sncovr does not exist in the restart, need to create it
         if (all(sncovr < zero)) then
           if (me == master ) write(0,'(a)') 'GFS_phys_time_vary_init: compute sncovr from weasd and soil vegetation parameters'
           !--- compute sncovr from existing variables
           !--- code taken directly from read_fix.f
           sncovr(:) = zero
           do ix=1,im
             if (landfrac(ix) >= drythresh .or. fice(ix) >= min_seaice) then
               vegtyp = vtype(ix)
               if (vegtyp == 0) vegtyp = 7
               rsnow  = 0.001_kind_phys*weasd(ix)/snupx(vegtyp)
               if (0.001_kind_phys*weasd(ix) < snupx(vegtyp)) then
                 sncovr(ix) = one - (exp(-salp_data*rsnow) - rsnow*exp(-salp_data))
               else
                 sncovr(ix) = one
               endif
             endif
           enddo
         endif

         !--- For RUC LSM: create sncovr_ice from sncovr
         if (lsm == lsm_ruc) then
           if (all(sncovr_ice < zero)) then
             if (me == master ) write(0,'(a)') 'GFS_phys_time_vary_init: fill sncovr_ice with sncovr for RUC LSM'
             sncovr_ice(:) = sncovr(:)
           endif
         endif

!$OMP end sections

!$OMP end parallel
         if (iaerclm) then
           call read_aerdataf (iamin, iamax, jamin, jamax, me,master,iflip,            &
                              idate,errmsg,errflg)
         endif

         if (lsm == lsm_noahmp) then
           if (all(tvxy < zero)) then

             allocate(dzsno (lsnow_lsm_lbound:lsnow_lsm_ubound))
             allocate(dzsnso(lsnow_lsm_lbound:lsoil)           )
             dzsno(:)    = missing_value
             dzsnso(:)   = missing_value

             tvxy(:)     = missing_value
             tgxy(:)     = missing_value
             tahxy(:)    = missing_value
             canicexy(:) = missing_value
             canliqxy(:) = missing_value
             eahxy(:)    = missing_value
             cmxy(:)     = missing_value
             chxy(:)     = missing_value
             fwetxy(:)   = missing_value
             sneqvoxy(:) = missing_value
             alboldxy(:) = missing_value
             qsnowxy(:)  = missing_value
             wslakexy(:) = missing_value
             albdvis(:)  = missing_value
             albdnir(:)  = missing_value
             albivis(:)  = missing_value
             albinir(:)  = missing_value
             emiss(:)    = missing_value
             taussxy(:)  = missing_value
             waxy(:)     = missing_value
             wtxy(:)     = missing_value
             zwtxy(:)    = missing_value
             xlaixy(:)   = missing_value
             xsaixy(:)   = missing_value

             lfmassxy(:)   = missing_value
             stmassxy(:)   = missing_value
             rtmassxy(:)   = missing_value
             woodxy(:)     = missing_value
             stblcpxy(:)   = missing_value
             fastcpxy(:)   = missing_value
             smcwtdxy(:)   = missing_value
             deeprechxy(:) = missing_value
             rechxy(:)     = missing_value

             snowxy (:)    = missing_value
             snicexy(:,:)  = missing_value
             snliqxy(:,:)  = missing_value
             tsnoxy (:,:)  = missing_value
             smoiseq(:,:)  = missing_value
             zsnsoxy(:,:)  = missing_value

             do ix=1,im
               if (landfrac(ix) >= drythresh) then
                 tvxy(ix)     = tsfcl(ix)
                 tgxy(ix)     = tsfcl(ix)
                 tahxy(ix)    = tsfcl(ix)

                 if (snowd(ix) > 0.01_kind_phys .and. tsfcl(ix) > con_t0c ) tvxy(ix)  = con_t0c
                 if (snowd(ix) > 0.01_kind_phys .and. tsfcl(ix) > con_t0c ) tgxy(ix)  = con_t0c
                 if (snowd(ix) > 0.01_kind_phys .and. tsfcl(ix) > con_t0c ) tahxy(ix) = con_t0c

                 canicexy(ix) = 0.0_kind_phys
                 canliqxy(ix) = canopy(ix)

                 eahxy(ix)    = 2000.0_kind_phys

                 cmxy(ix)     = zero
                 chxy(ix)     = zero
                 fwetxy(ix)   = zero
                 sneqvoxy(ix) = weasd(ix)     ! mm
                 alboldxy(ix) = 0.65_kind_phys
                 qsnowxy(ix)  = zero

!                 if (srflag(ix) > 0.001) qsnowxy(ix) = tprcp(ix)/dtp
                 ! already set to 0.0
                 wslakexy(ix) = zero
                 taussxy(ix)  = zero
                 albdvis(ix)  = 0.2_kind_phys
                 albdnir(ix)  = 0.2_kind_phys
                 albivis(ix)  = 0.2_kind_phys
                 albinir(ix)  = 0.2_kind_phys
                 emiss(ix)    = 0.95_kind_phys


                 waxy(ix)     = 4900.0_kind_phys
                 wtxy(ix)     = waxy(ix)
                 zwtxy(ix)    = (25.0_kind_phys + 2.0_kind_phys) - waxy(ix) / 1000.0_kind_phys / 0.2_kind_phys

                 vegtyp       = vtype(ix)
                 if (vegtyp == 0) vegtyp = 7
                 imn          = idate(2)

                 if ((vegtyp == isbarren_table) .or. (vegtyp == isice_table) .or. (vegtyp == isurban_table) .or. (vegtyp == iswater_table)) then

                   xlaixy(ix)   = zero
                   xsaixy(ix)   = zero

                   lfmassxy(ix) = zero
                   stmassxy(ix) = zero
                   rtmassxy(ix) = zero

                   woodxy   (ix) = zero
                   stblcpxy (ix) = zero
                   fastcpxy (ix) = zero

                 else

                   xlaixy(ix)   = max(laim_table(vegtyp, imn),0.05_kind_phys)
!                   xsaixy(ix)   = max(saim_table(vegtyp, imn),0.05)
                   xsaixy(ix)   = max(xlaixy(ix)*0.1_kind_phys,0.05_kind_phys)

                   masslai      = 1000.0_kind_phys / max(sla_table(vegtyp),one)
                   lfmassxy(ix) = xlaixy(ix)*masslai
                   masssai      = 1000.0_kind_phys / 3.0_kind_phys
                   stmassxy(ix) = xsaixy(ix)* masssai

                   rtmassxy(ix) = 500.0_kind_phys

                   woodxy(ix)   = 500.0_kind_phys
                   stblcpxy(ix) = 1000.0_kind_phys
                   fastcpxy(ix) = 1000.0_kind_phys

                 endif  ! non urban ...

                 if (vegtyp == isice_table) then
                   do is = 1,lsoil
                     stc(ix,is) = min(stc(ix,is),min(tg3(ix),263.15_kind_phys))
                     smc(ix,is) = one
                     slc(ix,is) = zero
                   enddo
                 endif

                 snd = snowd(ix)/1000.0_kind_phys  ! go to m from snwdph

                 if (weasd(ix) /= zero .and. snd == zero ) then
                   snd = weasd(ix)/1000.0
                 endif

                 if (vegtyp == 15) then                      ! land ice in MODIS/IGBP
                   if (weasd(ix) < 0.1_kind_phys) then
                     weasd(ix) = 0.1_kind_phys
                     snd       = 0.01_kind_phys
                   endif
                 endif

                 if (snd < 0.025_kind_phys ) then
                   snowxy(ix)   = zero
                   dzsno(-2:0)  = zero
                 elseif (snd >= 0.025_kind_phys .and. snd <= 0.05_kind_phys ) then
                   snowxy(ix)   = -1.0_kind_phys
                   dzsno(0)     = snd
                 elseif (snd > 0.05_kind_phys .and. snd <= 0.10_kind_phys ) then
                   snowxy(ix)   = -2.0_kind_phys
                   dzsno(-1)    = 0.5_kind_phys*snd
                   dzsno(0)     = 0.5_kind_phys*snd
                 elseif (snd > 0.10_kind_phys .and. snd <= 0.25_kind_phys ) then
                   snowxy(ix)   = -2.0_kind_phys
                   dzsno(-1)    = 0.05_kind_phys
                   dzsno(0)     = snd - 0.05_kind_phys
                 elseif (snd > 0.25_kind_phys .and. snd <= 0.45_kind_phys ) then
                   snowxy(ix)   = -3.0_kind_phys
                   dzsno(-2)    = 0.05_kind_phys
                   dzsno(-1)    = 0.5_kind_phys*(snd-0.05_kind_phys)
                   dzsno(0)     = 0.5_kind_phys*(snd-0.05_kind_phys)
                 elseif (snd > 0.45_kind_phys) then
                   snowxy(ix)   = -3.0_kind_phys
                   dzsno(-2)    = 0.05_kind_phys
                   dzsno(-1)    = 0.20_kind_phys
                   dzsno(0)     = snd - 0.05_kind_phys - 0.20_kind_phys
                 else
                   errmsg = 'Error in GFS_phys_time_vary.fv3.F90: Problem with the logic assigning snow layers in Noah MP initialization'
                   errflg = 1
                   return
                 endif

! Now we have the snowxy field
! snice + snliq + tsno allocation and compute them from what we have

                 tsnoxy(ix,:)  = zero
                 snicexy(ix,:) = zero
                 snliqxy(ix,:) = zero
                 zsnsoxy(ix,:) = zero

                 isnow = nint(snowxy(ix))+1 ! snowxy <=0.0, dzsno >= 0.0

                 do is = isnow,0
                   tsnoxy(ix,is)  = tgxy(ix)
                   snliqxy(ix,is) = zero
                   snicexy(ix,is) = one * dzsno(is) * weasd(ix)/snd
                 enddo
!
!zsnsoxy, all negative ?
!
                 do is = isnow,0
                   dzsnso(is) = -dzsno(is)
                 enddo

                 do is = 1,4
                   dzsnso(is) = -dzs(is)
                 enddo
!
! Assign to zsnsoxy
!
                 zsnsoxy(ix,isnow) = dzsnso(isnow)
                 do is = isnow+1,4
                   zsnsoxy(ix,is) = zsnsoxy(ix,is-1) + dzsnso(is)
                 enddo
!
! smoiseq
! Init water table related quantities here
!
                 soiltyp  = stype(ix)
                 if (soiltyp /= 0) then
                   bexp   = bexp_table(soiltyp)
                   smcmax = smcmax_table(soiltyp)
                   smcwlt = smcwlt_table(soiltyp)
                   dwsat  = dwsat_table(soiltyp)
                   dksat  = dksat_table(soiltyp)
                   psisat = -psisat_table(soiltyp)
                 endif

                 if (vegtyp == isurban_table) then
                   smcmax = 0.45_kind_phys
                   smcwlt = 0.40_kind_phys
                 endif

                 if ((bexp > zero) .and. (smcmax > zero) .and. (-psisat > zero)) then
                   do is = 1, lsoil
                     if ( is == 1 )then
                       ddz = -zs(is+1) * 0.5_kind_phys
                     elseif ( is < lsoil ) then
                       ddz = ( zs(is-1) - zs(is+1) ) * 0.5_kind_phys
                     else
                       ddz = zs(is-1) - zs(is)
                     endif
                     smoiseq(ix,is) = min(max(find_eq_smc(bexp, dwsat, dksat, ddz, smcmax),1.e-4_kind_phys),smcmax*0.99_kind_phys)
                   enddo
                 else                                    ! bexp <= 0.0
                   smoiseq(ix,1:4) = smcmax
                 endif                                   ! end the bexp condition

                 smcwtdxy(ix)   = smcmax
                 deeprechxy(ix) = zero
                 rechxy(ix)     = zero

               endif

             enddo ! ix

             deallocate(dzsno)
             deallocate(dzsnso)

           endif
         endif   !if Noah MP cold start ends

         is_initialized = .true.

      contains

!
! Use newton-raphson method to find eq soil moisture
!
         function find_eq_smc(bexp, dwsat, dksat, ddz, smcmax) result(smc)
            implicit none
            real(kind=kind_phys), intent(in) :: bexp, dwsat, dksat, ddz, smcmax
            real(kind=kind_phys) :: smc
            real(kind=kind_phys) :: expon, aa, bb, func, dfunc, dx
            integer :: iter
            !
            expon = bexp + 1.
            aa    = dwsat / ddz
            bb    = dksat / smcmax ** expon
            smc = 0.5 * smcmax
            !
            do iter = 1,100
              func  = (smc - smcmax) * aa +  bb * smc ** expon
              dfunc = aa + bb * expon * smc ** bexp
              dx    = func / dfunc
              smc   = smc - dx
              if ( abs (dx) < 1.e-6_kind_phys) return
            enddo
         end function find_eq_smc

      end subroutine GFS_phys_time_vary_init
!! @}

!> \section arg_table_GFS_phys_time_vary_timestep_init Argument Table
!! \htmlinclude GFS_phys_time_vary_timestep_init.html
!!
!>\section gen_GFS_phys_time_vary_timestep_init GFS_phys_time_vary_timestep_init General Algorithm
!! @{
      subroutine GFS_phys_time_vary_timestep_init (                                                 &
            me, master, cnx, cny, isc, jsc, nrcm, im, levs, kdt, idate, nsswr, fhswr, lsswr, fhour, &
            imfdeepcnv, cal_pre, random_clds, nscyc, ntoz, h2o_phys, iaerclm, iccn, clstp,          &
            jindx1_o3, jindx2_o3, ddy_o3, ozpl, jindx1_h, jindx2_h, ddy_h, h2opl,                   &
            jindx1_aer, jindx2_aer, ddy_aer, iindx1_aer, iindx2_aer, ddx_aer, aer_nm,               &
            jindx1_ci, jindx2_ci, ddy_ci, iindx1_ci, iindx2_ci, ddx_ci, in_nm, ccn_nm,              &
            imap, jmap, prsl, seed0, rann, nthrds, nx, ny, nsst, tile_num, nlunit, lsoil, lsoil_lsm,&
            kice, ialb, isot, ivegsrc, input_nml_file, use_ufo, nst_anl, frac_grid, fhcyc, phour,   &
            lakefrac, min_seaice, min_lakeice, smc, slc, stc, smois, sh2o, tslb, tiice, tg3, tref,  &
            tsfc, tsfco, tisfc, hice, fice, facsf, facwf, alvsf, alvwf, alnsf, alnwf, zorli, zorll, &
            zorlo, weasd, slope, snoalb, canopy, vfrac, vtype, stype, shdmin, shdmax, snowd,        &
            cv, cvb, cvt, oro, oro_uf, xlat_d, xlon_d, slmsk,                                       &
            do_ugwp_v1, jindx1_tau, jindx2_tau, ddy_j1tau, ddy_j2tau, tau_amf, errmsg, errflg)

         implicit none

         ! Interface variables
         integer,              intent(in)    :: me, master, cnx, cny, isc, jsc, nrcm, im, levs, kdt, &
                                                nsswr, imfdeepcnv, iccn, nscyc, ntoz
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

         logical,              intent(in)    :: do_ugwp_v1
         integer,              intent(in)    :: jindx1_tau(:), jindx2_tau(:)
         real(kind_phys),      intent(in)    :: ddy_j1tau(:), ddy_j2tau(:)
         real(kind_phys),      intent(inout) :: tau_amf(:)

         ! For gcycle only
         integer,              intent(in)    :: nthrds, nx, ny, nsst, tile_num, nlunit, lsoil
         integer,              intent(in)    :: lsoil_lsm, kice, ialb, isot, ivegsrc
         character(len=*),     intent(in)    :: input_nml_file(:)
         logical,              intent(in)    :: use_ufo, nst_anl, frac_grid
         real(kind_phys),      intent(in)    :: fhcyc, phour, lakefrac(:), min_seaice, min_lakeice,  &
                                                xlat_d(:), xlon_d(:)
         real(kind_phys),      intent(inout) :: smc(:,:), slc(:,:), stc(:,:), smois(:,:), sh2o(:,:), &
                                      tslb(:,:), tiice(:,:), tg3(:), tref(:),                        &
                                      tsfc(:), tsfco(:), tisfc(:), hice(:), fice(:),                 &
                                      facsf(:), facwf(:), alvsf(:), alvwf(:), alnsf(:), alnwf(:),    &
                                      zorli(:), zorll(:), zorlo(:), weasd(:), slope(:), snoalb(:),   &
                                      canopy(:), vfrac(:), vtype(:), stype(:), shdmin(:), shdmax(:), &
                                      snowd(:), cv(:), cvb(:), cvt(:), oro(:), oro_uf(:), slmsk(:)

         character(len=*),     intent(out)   :: errmsg
         integer,              intent(out)   :: errflg

         ! Local variables
         integer :: i, j, k, iseed, iskip, ix
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
           call ciinterpol (me, im, idate, fhour,     &
                            jindx1_ci, jindx2_ci,     &
                            ddy_ci, iindx1_ci,        &
                            iindx2_ci, ddx_ci,        &
                            levs, prsl, in_nm, ccn_nm)
         endif

!> - Call  cires_indx_ugwp to read monthly-mean GW-tau diagnosed from FV3GFS-runs that resolve GW-activ
         if (do_ugwp_v1) then
           call tau_amf_interp(me, master, im, idate, fhour, &
                               jindx1_tau, jindx2_tau,       &
                               ddy_j1tau, ddy_j2tau, tau_amf)
         endif

!> - Call gcycle() to repopulate specific time-varying surface properties for AMIP/forecast runs
         if (nscyc >  0) then
           if (mod(kdt,nscyc) == 1) THEN
             call gcycle (me, nthrds, nx, ny, isc, jsc, nsst, tile_num, nlunit,       &
                 input_nml_file, lsoil, lsoil_lsm, kice, idate, ialb, isot, ivegsrc,  &
                 use_ufo, nst_anl, fhcyc, phour, lakefrac, min_seaice, min_lakeice,   &
                 frac_grid, smc, slc, stc, smois, sh2o, tslb, tiice, tg3, tref, tsfc, &
                 tsfco, tisfc, hice, fice, facsf, facwf, alvsf, alvwf, alnsf, alnwf,  &
                 zorli, zorll, zorlo, weasd, slope, snoalb, canopy, vfrac, vtype,     &
                 stype, shdmin, shdmax, snowd, cv, cvb, cvt, oro, oro_uf,             &
                 xlat_d, xlon_d, slmsk, imap, jmap)
           endif
         endif

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

         ! Deallocate UGWP-input arrays
         if (allocated(ugwp_taulat)) deallocate(ugwp_taulat)
         if (allocated(tau_limb   )) deallocate(tau_limb)
         if (allocated(days_limb  )) deallocate(days_limb)

         is_initialized = .false.

      end subroutine GFS_phys_time_vary_finalize

   end module GFS_phys_time_vary
!> @}

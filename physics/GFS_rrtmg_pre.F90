!> \file GFS_rrtmg_pre.f90
!! This file contains
      module GFS_rrtmg_pre

      public GFS_rrtmg_pre_run

      contains

!> \defgroup GFS_rrtmg_pre GFS RRTMG Scheme Pre
!! @{
      subroutine GFS_rrtmg_pre_init ()
      end subroutine GFS_rrtmg_pre_init

!> \section arg_table_GFS_rrtmg_pre_run Argument Table
!! \htmlinclude GFS_rrtmg_pre_run.html
!!
      ! Attention - the output arguments lm, im, lmk, lmp must not be set
      ! in the CCPP version - they are defined in the interstitial_create routine
      subroutine GFS_rrtmg_pre_run (im, levs, lm, lmk, lmp, n_var_lndp,        &
        imfdeepcnv, imfdeepcnv_gf, me, ncnd, ntrac, num_p3d, npdf3d, ncnvcld3d,&
        ntqv, ntcw,ntiw, ntlnc, ntinc, ncld, ntrw, ntsw, ntgl, ntwa, ntoz,     &
        ntclamt, nleffr, nieffr, nseffr, lndp_type, kdt, imp_physics,          &
        imp_physics_thompson, imp_physics_gfdl, imp_physics_zhao_carr,         &
        imp_physics_zhao_carr_pdf, imp_physics_mg, imp_physics_wsm6,           &
        imp_physics_fer_hires, julian, yearlen, lndp_var_list, lsswr, lslwr,   &
        ltaerosol, lgfdlmprad, uni_cld, effr_in, do_mynnedmf, lmfshal,         &
        lmfdeep2, fhswr, fhlwr, solhr, sup, con_eps, epsm1, fvirt,             &
        rog, rocp, con_rd, xlat_d, xlat, xlon, coslat, sinlat, tsfc, slmsk,    &
        prsi, prsl, prslk, tgrs, sfc_wts, mg_cld, effrr_in, pert_clds,sppt_wts,&
        sppt_amp, cnvw_in, cnvc_in, qgrs, aer_nm, dx, icloud,                  & !inputs from here and above
        coszen, coszdg, effrl_inout, effri_inout, effrs_inout,                 &
        clouds1, clouds2, clouds3, clouds4, clouds5,                           & !in/out from here and above
        kd, kt, kb, mtopa, mbota, raddt, tsfg, tsfa, de_lgth, alb1d, delp, dz, & !output from here and below
        plvl, plyr, tlvl, tlyr, qlyr, olyr, gasvmr_co2, gasvmr_n2o, gasvmr_ch4,&
        gasvmr_o2, gasvmr_co, gasvmr_cfc11, gasvmr_cfc12, gasvmr_cfc22,        &
        gasvmr_ccl4,  gasvmr_cfc113, aerodp, clouds6, clouds7, clouds8,        &
        clouds9, cldsa, cldfra, faersw1, faersw2, faersw3, faerlw1, faerlw2,   &
        faerlw3, alpha, errmsg, errflg)

      use machine,                   only: kind_phys

      use physparam

      use radcons,                   only: itsfc,ltp, lextop, qmin,  &
                                           qme5, qme6, epsq, prsmin
      use funcphys,                  only: fpvs

      use module_radiation_astronomy,only: coszmn                      ! sol_init, sol_update
      use module_radiation_gases,    only: NF_VGAS, getgases, getozn   ! gas_init, gas_update,
      use module_radiation_aerosols, only: NF_AESW, NF_AELW, setaer, & ! aer_init, aer_update,
     &                                     NSPC1
      use module_radiation_clouds,   only: NF_CLDS,                  & ! cld_init
     &                                     progcld1, progcld3,       &
     &                                     progcld2,                 &
     &                                     progcld4, progcld5,       &
     &                                     progcld6,                 &
     &                                     progclduni,               &
     &                                     cal_cldfra3,              &
     &                                     find_cloudLayers,         &
     &                                     adjust_cloudIce,          &
     &                                     adjust_cloudH2O,          &
     &                                     adjust_cloudFinal

      use module_radsw_parameters,   only: topfsw_type, sfcfsw_type, &
     &                                     profsw_type, NBDSW
      use module_radlw_parameters,   only: topflw_type, sfcflw_type, &
     &                                     proflw_type, NBDLW
      use surface_perturbation,      only: cdfnor,ppfbet

      ! For Thompson MP
      use module_mp_thompson,        only: calc_effectRad, Nt_c,     &
                                           re_qc_min, re_qc_max,     &
                                           re_qi_min, re_qi_max,     &
                                           re_qs_min, re_qs_max
      use module_mp_thompson_make_number_concentrations, only:       &
                                           make_IceNumber,           &
                                           make_DropletNumber,       &
                                           make_RainNumber

      implicit none

      integer,              intent(in)  :: im, levs, lm, lmk, lmp, n_var_lndp, &
                                           imfdeepcnv,                         &
                                           imfdeepcnv_gf, me, ncnd, ntrac,     &
                                           num_p3d, npdf3d, ncnvcld3d, ntqv,   &
                                           ntcw, ntiw, ntlnc, ntinc, ncld,     &
                                           ntrw, ntsw, ntgl, ntwa, ntoz,       &
                                           ntclamt, nleffr, nieffr, nseffr,    &
                                           lndp_type,                          &
                                           kdt, imp_physics,                   &
                                           imp_physics_thompson,               &
                                           imp_physics_gfdl,                   &
                                           imp_physics_zhao_carr,              &
                                           imp_physics_zhao_carr_pdf,          &
                                           imp_physics_mg, imp_physics_wsm6,   &
                                           imp_physics_fer_hires,              &
                                           yearlen, icloud

      character(len=3), dimension(:), intent(in) :: lndp_var_list

      logical,              intent(in) :: lsswr, lslwr, ltaerosol, lgfdlmprad, &
                                          uni_cld, effr_in, do_mynnedmf,       &
                                          lmfshal, lmfdeep2, pert_clds

      real(kind=kind_phys), intent(in) :: fhswr, fhlwr, solhr, sup, julian, sppt_amp
      real(kind=kind_phys), intent(in) :: con_eps, epsm1, fvirt, rog, rocp, con_rd

      real(kind=kind_phys), dimension(:), intent(in) :: xlat_d, xlat, xlon,    &
                                                        coslat, sinlat, tsfc,  &
                                                        slmsk, dx

      real(kind=kind_phys), dimension(:,:), intent(in) :: prsi, prsl, prslk,   &
                                                          tgrs, sfc_wts,       &
                                                          mg_cld, effrr_in,    &
                                                          cnvw_in, cnvc_in,    &
                                                          sppt_wts

      real(kind=kind_phys), dimension(:,:,:), intent(in) :: qgrs, aer_nm

      real(kind=kind_phys), dimension(:),   intent(inout) :: coszen, coszdg

      real(kind=kind_phys), dimension(:,:), intent(inout) :: effrl_inout,      &
                                                             effri_inout,      &
                                                             effrs_inout
      real(kind=kind_phys), dimension(im,lm+LTP), intent(inout) :: clouds1,    &
                                                             clouds2, clouds3, &
                                                             clouds4, clouds5

      integer,                                      intent(out) :: kd, kt, kb

      integer, dimension(im,3),                     intent(out) :: mbota, mtopa

      real(kind=kind_phys),                         intent(out) :: raddt

      real(kind=kind_phys), dimension(im),          intent(out) :: tsfg, tsfa
      real(kind=kind_phys), dimension(im),          intent(out) :: de_lgth,    &
                                                                   alb1d

      real(kind=kind_phys), dimension(im,lm+LTP),   intent(out) :: delp, dz,   &
                                                                   plyr, tlyr, &
                                                                   qlyr, olyr

      real(kind=kind_phys), dimension(im,lm+1+LTP), intent(out) :: plvl, tlvl



      real(kind=kind_phys), dimension(im,lm+LTP),   intent(out) :: gasvmr_co2, &
                                                                   gasvmr_n2o, &
                                                                   gasvmr_ch4, &
                                                                   gasvmr_o2,  &
                                                                   gasvmr_co,  &
                                                                   gasvmr_cfc11,&
                                                                   gasvmr_cfc12,&
                                                                   gasvmr_cfc22,&
                                                                   gasvmr_ccl4,&
                                                                   gasvmr_cfc113
      real(kind=kind_phys), dimension(im,NSPC1),     intent(out) :: aerodp
      real(kind=kind_phys), dimension(im,lm+LTP),    intent(out) :: clouds6,   &
                                                                    clouds7,   &
                                                                    clouds8,   &
                                                                    clouds9,   &
                                                                    cldfra
      real(kind=kind_phys), dimension(im,5),            intent(out) :: cldsa

      real(kind=kind_phys), dimension(im,lm+LTP,NBDSW), intent(out) :: faersw1,&
                                                                       faersw2,&
                                                                       faersw3

      real(kind=kind_phys), dimension(im,lm+LTP,NBDLW), intent(out) :: faerlw1,&
                                                                       faerlw2,&
                                                                       faerlw3
      real(kind=kind_phys), dimension(im,lm+LTP),       intent(out) :: alpha

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      ! Local variables
      integer :: ncndl

      integer :: i, j, k, k1, k2, lsk, lv, n, itop, ibtc, LP1, lla, llb, lya,lyb

      real(kind=kind_phys) :: es, qs, delt, tem0d, gridkm, pfac

      real(kind=kind_phys), dimension(im) :: cvt1, cvb1, tem1d, tskn, xland

      real(kind=kind_phys), dimension(im,lm+LTP) ::         &
                          htswc, htlwc, gcice, grain, grime, htsw0, htlw0, &
                          rhly, tvly,qstl, vvel, clw, ciw, prslk1, tem2da, &
                          dzb, hzb, cldcov, deltaq, cnvc, cnvw,            &
                          effrl, effri, effrr, effrs, rho, orho, plyrpa

      ! for Thompson MP
      real(kind=kind_phys), dimension(im,lm+LTP) ::         &
                                  re_cloud, re_ice, re_snow, qv_mp, qc_mp, &
                                  qi_mp, qs_mp, nc_mp, ni_mp, nwfa

      ! for F-A MP
      real(kind=kind_phys), dimension(im,lm+LTP)   :: qc_save, qi_save, qs_save
      real(kind=kind_phys), dimension(im,lm+LTP+1) :: tem2db, hz

      real(kind=kind_phys), dimension(im,lm+LTP,min(4,ncnd))   :: ccnd
      real(kind=kind_phys), dimension(im,lm+LTP,2:ntrac)       :: tracer1
      real(kind=kind_phys), dimension(im,lm+LTP,NF_CLDS)       :: clouds
      real(kind=kind_phys), dimension(im,lm+LTP,NF_VGAS)       :: gasvmr
      real(kind=kind_phys), dimension(im,lm+LTP,NBDSW,NF_AESW) :: faersw
      real(kind=kind_phys), dimension(im,lm+LTP,NBDLW,NF_AELW) :: faerlw
      
      ! for stochastic cloud perturbations
      real(kind=kind_phys), dimension(im) :: cldp1d
      real (kind=kind_phys) :: alpha0,beta0,m,s,cldtmp,tmp_wt,cdfz
      integer  :: iflag

      integer :: ids, ide, jds, jde, kds, kde, &
                 ims, ime, jms, jme, kms, kme, &
                 its, ite, jts, jte, kts, kte

      real(kind=kind_phys) :: qvs
!
!===> ...  begin here
!
      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      if (.not. (lsswr .or. lslwr)) return

      !--- set commonly used integers
      ncndl = min(ncnd,4)

      LP1 = LM + 1               ! num of in/out levels


!  --- ...  set local /level/layer indexes corresponding to in/out
!  variables

      if ( lextop ) then
        if ( ivflip == 1 ) then    ! vertical from sfc upward
          kd = 0                   ! index diff between in/out and local
          kt = 1                   ! index diff between lyr and upper bound
          kb = 0                   ! index diff between lyr and lower bound
          lla = LMK                ! local index at the 2nd level from top
          llb = LMP                ! local index at toa level
          lya = LM                 ! local index for the 2nd layer from top
          lyb = LP1                ! local index for the top layer
        else                       ! vertical from toa downward
          kd = 1                   ! index diff between in/out and local
          kt = 0                   ! index diff between lyr and upper bound
          kb = 1                   ! index diff between lyr and lower bound
          lla = 2                  ! local index at the 2nd level from top
          llb = 1                  ! local index at toa level
          lya = 2                  ! local index for the 2nd layer from top
          lyb = 1                  ! local index for the top layer
        endif                      ! end if_ivflip_block
      else
        kd = 0
        if ( ivflip == 1 ) then    ! vertical from sfc upward
          kt = 1                   ! index diff between lyr and upper bound
          kb = 0                   ! index diff between lyr and lower bound
        else                       ! vertical from toa downward
          kt = 0                   ! index diff between lyr and upper bound
          kb = 1                   ! index diff between lyr and lower bound
        endif                      ! end if_ivflip_block
      endif   ! end if_lextop_block

      raddt = min(fhswr, fhlwr)
!     print *,' in grrad : raddt=',raddt


!> -# Setup surface ground temperature and ground/air skin temperature
!! if required.

      if ( itsfc == 0 ) then            ! use same sfc skin-air/ground temp
        do i = 1, IM
          tskn(i) = tsfc(i)
          tsfg(i) = tsfc(i)
        enddo
      else                              ! use diff sfc skin-air/ground temp
        do i = 1, IM
          tskn(i) = tsfc(i)
          tsfg(i) = tsfc(i)
        enddo
      endif


!> -# Prepare atmospheric profiles for radiation input.
!

      lsk = 0
      if (ivflip == 0 .and. lm < levs) lsk = levs - lm

!     convert pressure unit from pa to mb
      do k = 1, LM
        k1 = k + kd
        k2 = k + lsk
        do i = 1, IM
          plvl(i,k1+kb) = prsi(i,k2+kb) * 0.01   ! pa to mb (hpa)
          plyr(i,k1)    = prsl(i,k2)    * 0.01   ! pa to mb (hpa)
          tlyr(i,k1)    = tgrs(i,k2)
          prslk1(i,k1)  = prslk(i,k2)

!>  - Compute relative humidity.
          es  = min( prsl(i,k2),  fpvs( tgrs(i,k2) ) )  ! fpvs and prsl in pa
          qs  = max( QMIN, con_eps * es / (prsl(i,k2) + epsm1*es) )
          rhly(i,k1) = max( 0.0, min( 1.0, max(QMIN, qgrs(i,k2,ntqv))/qs ) )
          qstl(i,k1) = qs
        enddo
      enddo

      !--- recast remaining all tracers (except sphum) forcing them all to be positive
      do j = 2, ntrac
        do k = 1, LM
          k1 = k + kd
          k2 = k + lsk
          tracer1(:,k1,j) = max(0.0, qgrs(:,k2,j))
        enddo
      enddo
!
      if (ivflip == 0) then                                ! input data from toa to sfc
        if (lsk > 0) then
          k1 = 1 + kd
          k2 = k1 + kb
          do i = 1, IM
            plvl(i,k2)   = 0.01 * prsi(i,1+kb)          ! pa to mb (hpa)
            plyr(i,k1)   = 0.5 * (plvl(i,k2+1) + plvl(i,k2))
            prslk1(i,k1) = (plyr(i,k1)*0.001) ** rocp
          enddo
        else
          k1 = 1 + kd
          do i = 1, IM
            plvl(i,k1) = prsi(i,1) * 0.01   ! pa to mb (hpa)
          enddo
        endif
      else                                                 ! input data from sfc to top
        if (levs > lm) then
          k1 = lm + kd
          do i = 1, IM
            plvl(i,k1+1) = 0.01 * prsi(i,levs+1)  ! pa to mb (hpa)
            plyr(i,k1)   = 0.5 * (plvl(i,k1+1) + plvl(i,k1))
            prslk1(i,k1) = (plyr(i,k1)*0.001) ** rocp
          enddo
        else
          k1 = lp1 + kd
          do i = 1, IM
            plvl(i,k1) = prsi(i,lp1) * 0.01   ! pa to mb (hpa)
          enddo
        endif
      endif
!
      if ( lextop ) then                 ! values for extra top layer
        do i = 1, IM
          plvl(i,llb) = prsmin
          if ( plvl(i,lla) <= prsmin ) plvl(i,lla) = 2.0*prsmin
          plyr(i,lyb)   = 0.5 * plvl(i,lla)
          tlyr(i,lyb)   = tlyr(i,lya)
          prslk1(i,lyb) = (plyr(i,lyb)*0.001) ** rocp ! plyr in Pa
          rhly(i,lyb)   = rhly(i,lya)
          qstl(i,lyb)   = qstl(i,lya)
        enddo

!  ---  note: may need to take care the top layer amount
        tracer1(:,lyb,:) = tracer1(:,lya,:)
      endif


!>  - Get layer ozone mass mixing ratio (if use ozone climatology data,
!!    call getozn()).

      if (ntoz > 0) then            ! interactive ozone generation
        do k=1,lmk
          do i=1,im
            olyr(i,k) = max( QMIN, tracer1(i,k,ntoz) )
          enddo
        enddo
      else                                ! climatological ozone
        call getozn (prslk1, xlat, im, lmk,    &     !  ---  inputs
                     olyr)                           !  ---  outputs
      endif                               ! end_if_ntoz

!>  - Call coszmn(), to compute cosine of zenith angle (only when SW is called)
      if (lsswr) then
        call coszmn (xlon,sinlat,coslat,solhr,im,me, &     !  ---  inputs
                     coszen, coszdg)                       !  ---  outputs
      endif

!>  - Call getgases(), to set up non-prognostic gas volume mixing
!!    ratioes (gasvmr).
!  - gasvmr(:,:,1)  -  co2 volume mixing ratio
!  - gasvmr(:,:,2)  -  n2o volume mixing ratio
!  - gasvmr(:,:,3)  -  ch4 volume mixing ratio
!  - gasvmr(:,:,4)  -  o2  volume mixing ratio
!  - gasvmr(:,:,5)  -  co  volume mixing ratio
!  - gasvmr(:,:,6)  -  cf11 volume mixing ratio
!  - gasvmr(:,:,7)  -  cf12 volume mixing ratio
!  - gasvmr(:,:,8)  -  cf22 volume mixing ratio
!  - gasvmr(:,:,9)  -  ccl4 volume mixing ratio
!  - gasvmr(:,:,10) -  cfc113 volumne mixing ratio

!  --- ...  set up non-prognostic gas volume mixing ratioes

      call getgases (plvl, xlon, xlat, IM, LMK, & !  --- inputs
                     gasvmr)                      !  --- outputs

!CCPP: re-assign gasvmr(:,:,NF_VGAS) to gasvmr_X(:,:)
      do k = 1, LMK
        do i = 1, IM
           gasvmr_co2    (i,k)  = gasvmr(i,k,1)
           gasvmr_n2o    (i,k)  = gasvmr(i,k,2)
           gasvmr_ch4    (i,k)  = gasvmr(i,k,3)
           gasvmr_o2     (i,k)  = gasvmr(i,k,4)
           gasvmr_co     (i,k)  = gasvmr(i,k,5)
           gasvmr_cfc11  (i,k)  = gasvmr(i,k,6)
           gasvmr_cfc12  (i,k)  = gasvmr(i,k,7)
           gasvmr_cfc22  (i,k)  = gasvmr(i,k,8)
           gasvmr_ccl4   (i,k)  = gasvmr(i,k,9)
           gasvmr_cfc113 (i,k)  = gasvmr(i,k,10)
         enddo
      enddo

!>  - Get temperature at layer interface, and layer moisture.
      do k = 2, LMK
        do i = 1, IM
          tem2da(i,k) = log( plyr(i,k) )
          tem2db(i,k) = log( plvl(i,k) )
        enddo
      enddo

      if (ivflip == 0) then              ! input data from toa to sfc
        do i = 1, IM
          tem1d (i)   = QME6
          tem2da(i,1) = log( plyr(i,1) )
          tem2db(i,1) = log( max(prsmin, plvl(i,1)) )
          tem2db(i,LMP) = log( plvl(i,LMP) )
          tsfa  (i)   = tlyr(i,LMK)                  ! sfc layer air temp
          tlvl(i,1)   = tlyr(i,1)
          tlvl(i,LMP) = tskn(i)
        enddo

        do k = 1, LM
          k1 = k + kd
          do i = 1, IM
            qlyr(i,k1) = max( tem1d(i), qgrs(i,k,ntqv) )
            tem1d(i)   = min( QME5, qlyr(i,k1) )
            tvly(i,k1) = tgrs(i,k) * (1.0 + fvirt*qlyr(i,k1)) ! virtual T (K)
            delp(i,k1) = plvl(i,k1+1) - plvl(i,k1)
          enddo
        enddo

        if ( lextop ) then
          do i = 1, IM
            qlyr(i,lyb) = qlyr(i,lya)
            tvly(i,lyb) = tvly(i,lya)
            delp(i,lyb) = plvl(i,lla) - plvl(i,llb)
          enddo
        endif

        do k = 2, LMK
          do i = 1, IM
            tlvl(i,k) = tlyr(i,k) + (tlyr(i,k-1) - tlyr(i,k))           &
     &                * (tem2db(i,k)   - tem2da(i,k))                   &
     &                / (tem2da(i,k-1) - tem2da(i,k))
          enddo
        enddo

!  ---  ...  level height and layer thickness (km)
! dz:  Layer thickness between layer boundaries
! dzb: Layer thickness between layer centers (lowest is from surface to lowest layer center)
! hz:  Height of each level (i.e. layer boundary)
! hzb: Height of each layer center

        tem0d = 0.001 * rog
        do i = 1, IM
          do k = 1, LMK
            dz(i,k) = tem0d * (tem2db(i,k+1) - tem2db(i,k)) * tvly(i,k)
          enddo

          hz(i,LMP) = 0.0
          do k = LMK, 1, -1
            hz(i,k) = hz(i,k+1) + dz(i,k)
          enddo

          do k = LMK, 1, -1
            pfac = (tem2db(i,k+1) - tem2da(i,k)) / (tem2db(i,k+1) - tem2db(i,k))
            hzb(i,k) = hz(i,k+1) + pfac * (hz(i,k) - hz(i,k+1))
          enddo

          do k = LMK-1, 1, -1
            dzb(i,k) = hzb(i,k) - hzb(i,k+1)
          enddo
          dzb(i,LMK) = hzb(i,LMK) - hz(i,LMP)
        enddo

      else                               ! input data from sfc to toa

        do i = 1, IM
          tem1d (i)   = QME6
          tem2da(i,1) = log( plyr(i,1) )
          tem2db(i,1) = log( plvl(i,1) )
          tem2db(i,LMP) = log( max(prsmin, plvl(i,LMP)) )
          tsfa  (i)   = tlyr(i,1)                    ! sfc layer air temp
          tlvl(i,1)   = tskn(i)
          tlvl(i,LMP) = tlyr(i,LMK)
        enddo

        do k = LM, 1, -1
          do i = 1, IM
            qlyr(i,k) = max( tem1d(i), qgrs(i,k,ntqv) )
            tem1d(i)  = min( QME5, qlyr(i,k) )
            tvly(i,k) = tgrs(i,k) * (1.0 + fvirt*qlyr(i,k)) ! virtual T (K)
            delp(i,k) = plvl(i,k) - plvl(i,k+1)
          enddo
        enddo

        if ( lextop ) then
          do i = 1, IM
            qlyr(i,lyb) = qlyr(i,lya)
            tvly(i,lyb) = tvly(i,lya)
            delp(i,lyb) = plvl(i,lla) - plvl(i,llb)
          enddo
        endif

        do k = 1, LMK-1
          do i = 1, IM
            tlvl(i,k+1) = tlyr(i,k) + (tlyr(i,k+1) - tlyr(i,k))         &
     &                  * (tem2db(i,k+1) - tem2da(i,k))                 &
     &                  / (tem2da(i,k+1) - tem2da(i,k))
          enddo
        enddo

!  ---  ...  level height and layer thickness (km)
! dz:  Layer thickness between layer boundaries
! dzb: Layer thickness between layer centers (lowest is from surface to lowest layer center)
! hz:  Height of each level (i.e. layer boundary)
! hzb: Height of each layer center

        tem0d = 0.001 * rog
        do i = 1, IM
          do k = LMK, 1, -1
            dz(i,k) = tem0d * (tem2db(i,k) - tem2db(i,k+1)) * tvly(i,k)
          enddo

          hz(i,1) = 0.0
          do k = 1, LMK
            hz(i,k+1) = hz(i,k) + dz(i,k)
          enddo

          do k = 1, LMK
            pfac = (tem2db(i,k) - tem2da(i,k)) / (tem2db(i,k) - tem2db(i,k+1))
            hzb(i,k) = hz(i,k) + pfac * (hz(i,k+1) - hz(i,k))
          enddo

          do k = 2, LMK
            dzb(i,k) = hzb(i,k) - hzb(i,k-1)
          enddo
          dzb(i,1) = hzb(i,1) - hz(i,1)
        enddo

      endif                              ! end_if_ivflip

!>  - Call module_radiation_aerosols::setaer(),to setup aerosols
!! property profile for radiation.

!check  print *,' in grrad : calling setaer '

      call setaer (plvl, plyr, prslk1, tvly, rhly, slmsk,    & !  ---  inputs
                   tracer1, aer_nm, xlon, xlat, IM, LMK, LMP,&
                   lsswr,lslwr,                              &
                   faersw,faerlw,aerodp)                       !  ---  outputs

! CCPP
      do j = 1,NBDSW
        do k = 1, LMK
          do i = 1, IM
            ! NF_AESW = 3
            faersw1(i,k,j) = faersw(i,k,j,1)
            faersw2(i,k,j) = faersw(i,k,j,2)
            faersw3(i,k,j) = faersw(i,k,j,3)
          enddo
        enddo
       enddo

      do j = 1,NBDLW
        do k = 1, LMK
          do i = 1, IM
            ! NF_AELW = 3
            faerlw1(i,k,j) = faerlw(i,k,j,1)
            faerlw2(i,k,j) = faerlw(i,k,j,2)
            faerlw3(i,k,j) = faerlw(i,k,j,3)
          enddo
        enddo
       enddo

!>  - Obtain cloud information for radiation calculations
!!    (clouds,cldsa,mtopa,mbota)
!!\n   for  prognostic cloud:
!!    - For Zhao/Moorthi's prognostic cloud scheme,
!!      call module_radiation_clouds::progcld1()
!!    - For Zhao/Moorthi's prognostic cloud+pdfcld,
!!      call module_radiation_clouds::progcld3()
!!      call module_radiation_clouds::progclduni() for unified cloud and ncld=2

!  --- ...  obtain cloud information for radiation calculations

!      if (ntcw > 0) then                            ! prognostic cloud schemes
        ccnd = 0.0_kind_phys
        if (ncnd == 1) then                          ! Zhao_Carr_Sundqvist
          do k=1,LMK
            do i=1,IM
              ccnd(i,k,1) = tracer1(i,k,ntcw)        ! liquid water/ice
            enddo
          enddo
        elseif (ncnd == 2) then                      ! MG
          do k=1,LMK
            do i=1,IM
              ccnd(i,k,1) = tracer1(i,k,ntcw)        ! liquid water
              ccnd(i,k,2) = tracer1(i,k,ntiw)        ! ice water
            enddo
          enddo
        elseif (ncnd == 4) then                      ! MG2
          do k=1,LMK
            do i=1,IM
              ccnd(i,k,1) = tracer1(i,k,ntcw)                     ! liquid water
              ccnd(i,k,2) = tracer1(i,k,ntiw)                     ! ice water
              ccnd(i,k,3) = tracer1(i,k,ntrw)                     ! rain water
              ccnd(i,k,4) = tracer1(i,k,ntsw)                     ! snow water
            enddo
          enddo
        elseif (ncnd == 5) then                         ! GFDL MP, Thompson, MG3, FA
          do k=1,LMK
            do i=1,IM
              ccnd(i,k,1) = tracer1(i,k,ntcw)                     ! liquid water
              ccnd(i,k,2) = tracer1(i,k,ntiw)                     ! ice water
              ccnd(i,k,3) = tracer1(i,k,ntrw)                     ! rain water
              if (imp_physics == imp_physics_fer_hires ) then
                  ccnd(i,k,4) = 0.0
              else
                  ccnd(i,k,4) = tracer1(i,k,ntsw) + tracer1(i,k,ntgl) ! snow + graupel
              endif
            enddo
          enddo
          ! for Thompson MP - prepare variables for calc_effr
          if (imp_physics == imp_physics_thompson .and. ltaerosol) then
            do k=1,LMK
              do i=1,IM
                qvs = qgrs(i,k,ntqv)
                qv_mp (i,k) = qvs/(1.-qvs)
                rho   (i,k) = con_eps*prsl(i,k)/(con_rd*tgrs(i,k)*(qv_mp(i,k)+con_eps))
                orho  (i,k) = 1.0/rho(i,k)
                qc_mp (i,k) = tracer1(i,k,ntcw)/(1.-qvs)
                qi_mp (i,k) = tracer1(i,k,ntiw)/(1.-qvs)
                qs_mp (i,k) = tracer1(i,k,ntsw)/(1.-qvs)
                nc_mp (i,k) = tracer1(i,k,ntlnc)/(1.-qvs)
                ni_mp (i,k) = tracer1(i,k,ntinc)/(1.-qvs)
                nwfa  (i,k) = tracer1(i,k,ntwa)
              enddo
            enddo
          elseif (imp_physics == imp_physics_thompson) then
            do k=1,LMK
              do i=1,IM
                qvs = qgrs(i,k,ntqv)
                qv_mp (i,k) = qvs/(1.-qvs)
                rho   (i,k) = con_eps*prsl(i,k)/(con_rd*tgrs(i,k)*(qv_mp(i,k)+con_eps))
                orho  (i,k) = 1.0/rho(i,k)
                qc_mp (i,k) = tracer1(i,k,ntcw)/(1.-qvs)
                qi_mp (i,k) = tracer1(i,k,ntiw)/(1.-qvs)
                qs_mp (i,k) = tracer1(i,k,ntsw)/(1.-qvs)
                nc_mp (i,k) = nt_c*orho(i,k)
                ni_mp (i,k) = tracer1(i,k,ntinc)/(1.-qvs)
              enddo
            enddo
          endif
        endif
        do n=1,ncndl
          do k=1,LMK
            do i=1,IM
              if (ccnd(i,k,n) < epsq) ccnd(i,k,n) = 0.0
            enddo
          enddo
        enddo
        if (imp_physics == imp_physics_gfdl ) then
          if (.not. lgfdlmprad) then


! rsun the  summation methods and order make the difference in calculation

!            clw(:,:) = clw(:,:) + tracer1(:,1:LMK,ntcw)   &
!                                + tracer1(:,1:LMK,ntiw)   &
!                                + tracer1(:,1:LMK,ntrw)   &
!                                + tracer1(:,1:LMK,ntsw)   &
!                                + tracer1(:,1:LMK,ntgl)
            ccnd(:,:,1) =               tracer1(:,1:LMK,ntcw)
            ccnd(:,:,1) = ccnd(:,:,1) + tracer1(:,1:LMK,ntrw)
            ccnd(:,:,1) = ccnd(:,:,1) + tracer1(:,1:LMK,ntiw)
            ccnd(:,:,1) = ccnd(:,:,1) + tracer1(:,1:LMK,ntsw)
            ccnd(:,:,1) = ccnd(:,:,1) + tracer1(:,1:LMK,ntgl)

!          else
!            do j=1,ncld
!              ccnd(:,:,1) = ccnd(:,:,1) + tracer1(:,1:LMK,ntcw+j-1) ! cloud condensate amount
!            enddo
          endif
          do k=1,LMK
            do i=1,IM
              if (ccnd(i,k,1) < EPSQ ) ccnd(i,k,1) = 0.0
            enddo
          enddo
        endif
!
        if (uni_cld) then
          if (effr_in) then
            do k=1,lm
              k1 = k + kd
              do i=1,im
                cldcov(i,k1) = mg_cld(i,k)
                effrl(i,k1)  = effrl_inout(i,k)
                effri(i,k1)  = effri_inout(i,k)
                effrr(i,k1)  = effrr_in(i,k)
                effrs(i,k1)  = effrs_inout(i,k)
              enddo
            enddo
          else
            do k=1,lm
              k1 = k + kd
              do i=1,im
                cldcov(i,k1) = mg_cld(i,k)
              enddo
            enddo
          endif
        elseif (imp_physics == imp_physics_gfdl) then            ! GFDL MP
          if (do_mynnedmf .and. kdt>1) THEN
            do k=1,lm
              k1 = k + kd
              do i=1,im
                if (tracer1(i,k1,ntrw)>1.0e-7 .OR. tracer1(i,k1,ntsw)>1.0e-7) then
                ! GFDL cloud fraction
                  cldcov(i,k1) = tracer1(I,k1,ntclamt)
                else
                ! MYNN sub-grid cloud fraction
                  cldcov(i,k1) = clouds1(i,k1)
                endif
              enddo
            enddo
          else
            ! GFDL cloud fraction
            cldcov(1:IM,1+kd:LM+kd) = tracer1(1:IM,1:LM,ntclamt)
          endif
          if(effr_in) then
            do k=1,lm
              k1 = k + kd
              do i=1,im
                effrl(i,k1) = effrl_inout(i,k)
                effri(i,k1) = effri_inout(i,k)
                effrr(i,k1) = effrr_in(i,k)
                effrs(i,k1) = effrs_inout(i,k)
!                if(me==0) then
!                  if(effrl(i,k1)> 5.0) then
!                    write(6,*) 'rad driver:cloud radii:',kdt, i,k1,       &
!                    effrl(i,k1)
!                  endif
!                  if(effrs(i,k1)==0.0) then
!                    write(6,*) 'rad driver:snow mixing ratio:',kdt, i,k1, &
!                    tracer1(i,k,ntsw)
!                  endif
!                endif
              enddo
            enddo
          endif
        elseif (imp_physics == imp_physics_thompson) then       !  Thompson MP
          !
          ! Compute effective radii for QC, QI, QS with (GF, MYNN) or without (all others) sub-grid clouds
          !
          ! Update number concentration, consistent with sub-grid clouds (GF, MYNN) or without (all others)
          do k=1,lm
            do i=1,im
              if (ltaerosol .and. qc_mp(i,k)>1.e-12 .and. nc_mp(i,k)<100.) then
                nc_mp(i,k) = make_DropletNumber(qc_mp(i,k)*rho(i,k), nwfa(i,k)*rho(i,k)) * orho(i,k)
              endif
              if (qi_mp(i,k)>1.e-12 .and. ni_mp(i,k)<100.) then
                ni_mp(i,k) = make_IceNumber(qi_mp(i,k)*rho(i,k), tlyr(i,k)) * orho(i,k)
              endif
            end do
          end do
          ! Call Thompson's subroutine to compute effective radii
          do i=1,im
            ! Effective radii [m] are now intent(out), bounds applied in calc_effectRad
            !tgs: progclduni has different limits for ice radii (10.0-150.0) than
            !     calc_effectRad (4.99-125.0 for WRFv3.8.1; 2.49-125.0 for WRFv4+)
            !     it will raise the low limit from 5 to 10, but the high limit will remain 125.
            call calc_effectRad (tlyr(i,:), plyr(i,:)*100., qv_mp(i,:), qc_mp(i,:),   &
                                 nc_mp(i,:), qi_mp(i,:), ni_mp(i,:), qs_mp(i,:), &
                                 re_cloud(i,:), re_ice(i,:), re_snow(i,:), 1, lm )
            do k=1,lm
              re_cloud(i,k) = MAX(re_qc_min, MIN(re_cloud(i,k), re_qc_max))
              re_ice(i,k)   = MAX(re_qi_min, MIN(re_ice(i,k),   re_qi_max))
              re_snow(i,k)  = MAX(re_qs_min, MIN(re_snow(i,k),  re_qs_max))
            end do
          end do
          ! Scale Thompson's effective radii from meter to micron
          do k=1,lm
            do i=1,im
              re_cloud(i,k) = re_cloud(i,k)*1.e6
              re_ice(i,k)   = re_ice(i,k)*1.e6
              re_snow(i,k)  = re_snow(i,k)*1.e6
            end do
          end do
          do k=1,lm
            k1 = k + kd
            do i=1,im
              effrl(i,k1) = re_cloud (i,k)
              effri(i,k1) = re_ice (i,k)
              effrr(i,k1) = 1000. ! rrain_def=1000.
              effrs(i,k1) = re_snow(i,k)
            enddo
          enddo
          ! Update global arrays
          do k=1,lm
            k1 = k + kd
            do i=1,im
              effrl_inout(i,k) = effrl(i,k1)
              effri_inout(i,k) = effri(i,k1)
              effrs_inout(i,k) = effrs(i,k1)
            enddo
          enddo
        else                                                           ! all other cases
          cldcov = 0.0
        endif

!
!  --- add suspended convective cloud water to grid-scale cloud water
!      only for cloud fraction & radiation computation
!      it is to enhance cloudiness due to suspended convec cloud water
!      for zhao/moorthi's (imp_phys=99) &
!          ferrier's (imp_phys=5) microphysics schemes

        if ((num_p3d == 4) .and. (npdf3d == 3)) then       ! same as imp_physics = 98
          do k=1,lm
            k1 = k + kd
            do i=1,im
              !GJF: this is not consistent with GFS_typedefs,
              !     but it looks like the Zhao-Carr-PDF scheme is not in the CCPP
              deltaq(i,k1) = 0.0!Tbd%phy_f3d(i,k,5)      !GJF: this variable is not in phy_f3d anymore
              cnvw  (i,k1) = cnvw_in(i,k)
              cnvc  (i,k1) = cnvc_in(i,k)
            enddo
          enddo
        elseif ((npdf3d == 0) .and. (ncnvcld3d == 1)) then ! same as imp_physics=99
          do k=1,lm
            k1 = k + kd
            do i=1,im
              deltaq(i,k1) = 0.0
              cnvw  (i,k1) = cnvw_in(i,k)
              cnvc  (i,k1) = 0.0
            enddo
          enddo
        else                                                      ! all the rest
          do k=1,lmk
            do i=1,im
              deltaq(i,k) = 0.0
              cnvw  (i,k) = 0.0
              cnvc  (i,k) = 0.0
            enddo
          enddo
        endif

        !mz HWRF physics: icloud=3
        if(icloud == 3) then

          ! Set internal dimensions
          ids = 1
          ims = 1
          its = 1
          ide = size(xlon,1)
          ime = size(xlon,1)
          ite = size(xlon,1)
          jds = 1
          jms = 1
          jts = 1
          jde = 1
          jme = 1
          jte = 1
          kds = 1
          kms = 1
          kts = 1
          kde = lm+LTP ! should this be lmk instead of lm? no, or?
          kme = lm+LTP
          kte = lm+LTP

          do k = 1, LMK
            do i = 1, IM
              rho(i,k)=plyr(i,k)*100./(con_rd*tlyr(i,k))
              plyrpa(i,k)=plyr(i,k)*100.    !hPa->Pa
            end do
          end do

          do i=1,im
            if (slmsk(i)==1. .or. slmsk(i)==2.) then ! sea/land/ice mask (=0/1/2) in FV3
               xland(i)=1.0                          ! but land/water = (1/2) in HWRF
            else
               xland(i)=2.0
            endif
          enddo

          gridkm = sqrt(2.0)*sqrt(dx(1)*0.001*dx(1)*0.001)

          do i =1, im
            do k =1, lmk
               qc_save(i,k) = ccnd(i,k,1)  
               qi_save(i,k) = ccnd(i,k,2) 
               qs_save(i,k) = ccnd(i,k,4)
            enddo
          enddo


          call cal_cldfra3(cldcov,qlyr,ccnd(:,:,1),ccnd(:,:,2),      &
                           ccnd(:,:,4),plyrpa,tlyr,rho,xland,gridkm, &
                           ids,ide,jds,jde,kds,kde,                  &
                           ims,ime,jms,jme,kms,kme,                  &
                           its,ite,jts,jte,kts,kte)

          !mz* back to micro-only qc  qi,qs
          do i =1, im
            do k =1, lmk
              ccnd(i,k,1) = qc_save(i,k)
              ccnd(i,k,2) = qi_save(i,k)
              ccnd(i,k,4) = qs_save(i,k)
            enddo
          enddo

        endif ! icloud == 3

        if (lextop) then
          do i=1,im
            cldcov(i,lyb) = cldcov(i,lya)
            deltaq(i,lyb) = deltaq(i,lya)
            cnvw  (i,lyb) = cnvw  (i,lya)
            cnvc  (i,lyb) = cnvc  (i,lya)
          enddo
          if (effr_in) then
            do i=1,im
              effrl(i,lyb) = effrl(i,lya)
              effri(i,lyb) = effri(i,lya)
              effrr(i,lyb) = effrr(i,lya)
              effrs(i,lyb) = effrs(i,lya)
            enddo
          endif
        endif

        if (imp_physics == imp_physics_zhao_carr) then
          ccnd(1:IM,1:LMK,1) = ccnd(1:IM,1:LMK,1) + cnvw(1:IM,1:LMK)
        endif

        if (imp_physics == imp_physics_zhao_carr .or. imp_physics == imp_physics_mg) then ! zhao/moorthi's prognostic cloud scheme
                                         ! or unified cloud and/or with MG microphysics

          if (uni_cld .and. ncld >= 2) then
            call progclduni (plyr, plvl, tlyr, tvly, ccnd, ncndl,         & !  ---  inputs
                             xlat, xlon, slmsk, dz, delp,                 &
                             IM, LMK, LMP, cldcov,                        &
                             effrl, effri, effrr, effrs, effr_in,         &
                             dzb, xlat_d, julian, yearlen,                &
                             clouds, cldsa, mtopa, mbota, de_lgth, alpha)   !  ---  outputs
          else
            call progcld1 (plyr ,plvl, tlyr, tvly, qlyr, qstl, rhly,      & !  ---  inputs
                           ccnd(1:IM,1:LMK,1), xlat, xlon, slmsk, dz,     &
                           delp, IM, LMK, LMP, uni_cld, lmfshal, lmfdeep2,&
                           cldcov, effrl, effri, effrr, effrs, effr_in,   &
                           dzb, xlat_d, julian, yearlen,                  &
                           clouds, cldsa, mtopa, mbota, de_lgth, alpha)     !  ---  outputs
          endif

        elseif(imp_physics == imp_physics_zhao_carr_pdf) then      ! zhao/moorthi's prognostic cloud+pdfcld

          call progcld3 (plyr, plvl, tlyr, tvly, qlyr, qstl, rhly,        &  !  ---  inputs
                         ccnd(1:IM,1:LMK,1), cnvw, cnvc, xlat, xlon,      &
                         slmsk, dz, delp, im, lmk, lmp, deltaq, sup, kdt, &
                         me, dzb, xlat_d, julian, yearlen,                &
                         clouds, cldsa, mtopa, mbota, de_lgth, alpha)        !  ---  outputs

        elseif (imp_physics == imp_physics_gfdl) then           ! GFDL cloud scheme

          if (.not. lgfdlmprad) then
            call progcld4 (plyr, plvl, tlyr, tvly, qlyr, qstl, rhly,      &    !  ---  inputs
                           ccnd(1:IM,1:LMK,1), cnvw, cnvc, xlat, xlon,    &
                           slmsk, cldcov, dz, delp, im, lmk, lmp,         &
                           dzb, xlat_d, julian, yearlen,                  &
                           clouds, cldsa, mtopa, mbota, de_lgth, alpha)        !  ---  outputs
          else

            call progclduni (plyr, plvl, tlyr, tvly, ccnd, ncndl, xlat,   &    !  ---  inputs
                            xlon, slmsk, dz,delp, IM, LMK, LMP, cldcov,   &
                            effrl, effri, effrr, effrs, effr_in,          &
                            dzb, xlat_d, julian, yearlen,                 &
                            clouds, cldsa, mtopa, mbota, de_lgth, alpha)       !  ---  outputs
!           call progcld4o (plyr, plvl, tlyr, tvly, qlyr, qstl, rhly,       &   !  ---  inputs
!                           tracer1, xlat, xlon, slmsk, dz, delp,           &
!                           ntrac-1, ntcw-1,ntiw-1,ntrw-1,                  &
!                           ntsw-1,ntgl-1,ntclamt-1,                        &
!                           im, lmk, lmp,                                   &
!                           dzb, xlat_d, julian, yearlen,                   &
!                           clouds, cldsa, mtopa, mbota, de_lgth, alpha)        !  ---  outputs
          endif

        elseif(imp_physics == imp_physics_fer_hires) then
          if (kdt == 1) then
            effrl_inout(:,:) = 10.
            effri_inout(:,:) = 50.
            effrs_inout(:,:) = 250.
          endif

          call progcld5 (plyr,plvl,tlyr,tvly,qlyr,qstl,rhly,tracer1,       &  !  --- inputs
                         xlat,xlon,slmsk,dz,delp,                          &
                         ntrac-1, ntcw-1,ntiw-1,ntrw-1,                    &
!mz                       ntsw-1,ntgl-1,                                   &
                         im, lmk, lmp, icloud, uni_cld, lmfshal, lmfdeep2, &
                         cldcov(:,1:LMK),effrl_inout(:,:),                 &
                         effri_inout(:,:), effrs_inout(:,:),               &
                         dzb, xlat_d, julian, yearlen,                     &
                         clouds,cldsa,mtopa,mbota, de_lgth, alpha)            !  --- outputs

        elseif(imp_physics == imp_physics_thompson) then                              ! Thompson MP

          if(do_mynnedmf .or. imfdeepcnv == imfdeepcnv_gf ) then ! MYNN PBL or GF conv
              !-- MYNN PBL or convective GF
              !-- use cloud fractions with SGS clouds
              do k=1,lmk
                do i=1,im
                  clouds(i,k,1)  = clouds1(i,k)
                enddo
              enddo

                ! --- use clduni as with the GFDL microphysics.
                ! --- make sure that effr_in=.true. in the input.nml!
                call progclduni (plyr, plvl, tlyr, tvly, ccnd, ncndl,   & !  ---  inputs
                         xlat, xlon, slmsk, dz, delp, IM, LMK, LMP,     &
                         clouds(:,1:LMK,1),                             &
                         effrl, effri, effrr, effrs, effr_in ,          &
                         dzb, xlat_d, julian, yearlen,                  &
                         clouds, cldsa, mtopa, mbota, de_lgth, alpha)     !  ---  outputs

          else
            ! MYNN PBL or GF convective are not used
            call progcld6 (plyr,plvl,tlyr,qlyr,qstl,rhly,tracer1,   & !  --- inputs
                         xlat,xlon,slmsk,dz,delp,                   &
                         ntrac-1, ntcw-1,ntiw-1,ntrw-1,             &
                         ntsw-1,ntgl-1,                             &
                         im, lmk, lmp, uni_cld, lmfshal, lmfdeep2,  &
                         cldcov(:,1:LMK), effrl_inout(:,:),         &
                         effri_inout(:,:), effrs_inout(:,:),        &
                         dzb, xlat_d, julian, yearlen,              &
                         clouds, cldsa, mtopa ,mbota, de_lgth, alpha) !  --- outputs
          endif ! MYNN PBL or GF

        endif                            ! end if_imp_physics

!      endif                             ! end_if_ntcw

! perturb cld cover
       if (pert_clds) then
          do i=1,im
             tmp_wt= -1*log( ( 2.0 / ( sppt_wts(i,38) ) ) - 1 )
              call cdfnor(tmp_wt,cdfz)
              cldp1d(i) = cdfz
          enddo
          do k = 1, LMK
             do i = 1, IM
                ! compute beta distribution parameters
                m = clouds(i,k,1)
                if (m<0.99 .AND. m > 0.01) then
                   s = sppt_amp*m*(1.-m)
                   alpha0 = m*m*(1.-m)/(s*s)-m
                   beta0  = alpha0*(1.-m)/m
           ! compute beta distribution value corresponding
           ! to the given percentile albPpert to use as new albedo
                   call ppfbet(cldp1d(i),alpha0,beta0,iflag,cldtmp)
                   clouds(i,k,1) = cldtmp
                else
                   clouds(i,k,1) = m
                endif
             enddo     ! end_do_i_loop
          enddo     ! end_do_k_loop
       endif
       do k = 1, LMK
         do i = 1, IM
            clouds1(i,k)  = clouds(i,k,1)
            clouds2(i,k)  = clouds(i,k,2)
            clouds3(i,k)  = clouds(i,k,3)
            clouds4(i,k)  = clouds(i,k,4)
            clouds5(i,k)  = clouds(i,k,5)
            clouds6(i,k)  = clouds(i,k,6)
            clouds7(i,k)  = clouds(i,k,7)
            clouds8(i,k)  = clouds(i,k,8)
            clouds9(i,k)  = clouds(i,k,9)
            cldfra(i,k)   = clouds(i,k,1)
         enddo
       enddo

! mg, sfc-perts
!  ---  scale random patterns for surface perturbations with
!  perturbation size
!  ---  turn vegetation fraction pattern into percentile pattern
      alb1d(:) = 0.
      if (lndp_type==1) then
          do k =1,n_var_lndp
            if (lndp_var_list(k) == 'alb') then
              do i=1,im
                call cdfnor(sfc_wts(i,k),alb1d(i))
                !lndp_alb = lndp_prt_list(k)
              enddo
            endif
          enddo
      endif
! mg, sfc-perts

      end subroutine GFS_rrtmg_pre_run

      subroutine GFS_rrtmg_pre_finalize ()
      end subroutine GFS_rrtmg_pre_finalize

!! @}
      end module GFS_rrtmg_pre

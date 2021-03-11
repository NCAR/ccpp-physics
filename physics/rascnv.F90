!>  \file rascnv.F90
!!  This file contains the entire Relaxed Arakawa-Schubert convection
!!  parameteriztion

      module rascnv

      USE machine , ONLY : kind_phys
      implicit none
      public :: rascnv_init, rascnv_run, rascnv_finalize
      private
      logical :: is_initialized = .False.
!
      integer,               parameter :: kp = kind_phys
      integer,               parameter :: nrcmax=32 ! Maximum # of random clouds per 1200s

      integer,               parameter :: idnmax=999
      real (kind=kind_phys), parameter :: delt_c=1800.0_kp/3600.0_kp          &
!     Adjustment time scales in hrs for deep and shallow clouds
!    &,                                   adjts_d=3.0, adjts_s=0.5
!    &,                                   adjts_d=2.5, adjts_s=0.5
     &,                                   adjts_d=2.0_kp, adjts_s=0.5_kp
!
      logical,               parameter :: fix_ncld_hr=.true.

!
      real (kind=kind_phys), parameter :: ZERO=0.0_kp,      HALF=0.5_kp       &
     &,                                   pt25=0.25_kp,     ONE=1.0_kp        &
     &,                                   TWO=2.0_kp,       FOUR=4.0_kp       &
     &,                                   twoo3=two/3.0_kp                    &
     &,                                   FOUR_P2=4.0e2_kp, ONE_M10=1.0e-10_kp&
     &,                                   ONE_M6=1.0e-6_kp, ONE_M5=1.0e-5_kp  &
     &,                                   ONE_M2=1.0e-2_kp, ONE_M1=1.0e-1_kp  &
     &,                                   oneolog10=one/log(10.0_kp)          &
     &,                                   rain_min=1.0e-13_kp                 &
     &,                                   facmb=0.01_kp                       & ! conversion factor from Pa to hPa (or mb)
     &,                                   cmb2pa=100.0_kp                       ! Conversion from hPa to Pa
!
!     real (kind=kind_phys), parameter :: frac=0.5_kp,    crtmsf=0.0_kp     &
      real (kind=kind_phys), parameter :: frac=0.1_kp,    crtmsf=0.0_kp     &
     &,                                   tfrac_max=0.15_kp                 &
     &,                                   rhfacs=0.75_kp, rhfacl=0.75_kp    &
     &,                                   face=5.0_kp,    delx=10000.0_kp   &
     &,                                   ddfac=face*delx*0.001_kp          &
     &,                                   max_neg_bouy=0.15_kp              &
!    &,                                   max_neg_bouy=pt25_kp              &
     &,                                   testmb=0.1_kp, testmbi=one/testmb &
     &,                                   dpd=0.5_kp, rknob=1.0_kp, eknob=1.0_kp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     logical, parameter :: aw_scal=.false., cumfrc=.true.              &
      logical, parameter :: aw_scal=.true., cumfrc=.true.               &
     &,                     updret=.false., vsmooth=.false.             &
     &,                     wrkfun=.false., crtfun=.true.               &
     &,                     calkbl=.true.,  botop=.true.,  revap=.true. &
     &,                     advcld=.true.,  advups=.false.,advtvd=.true.
!    &,                     advcld=.true.,  advups=.true., advtvd=.false.
!    &,                     advcld=.true.,  advups=.false.,advtvd=.false.


      real(kind=kind_phys), parameter :: TF=233.16_kp, TCR=273.16_kp    &
     &,                                  TCRF=one/(TCR-TF), TCL=2.0_kp

!
!    For pressure gradient force in momentum mixing
!     real (kind=kind_phys), parameter :: pgftop=0.80, pgfbot=0.30       &
!    No pressure gradient force in momentum mixing
      real (kind=kind_phys), parameter :: pgftop=0.0_kp, pgfbot=0.0_kp   &
!     real (kind=kind_phys), parameter :: pgftop=0.55, pgfbot=0.55       &
     &,                                  pgfgrad=(pgfbot-pgftop)*0.001_kp&
     &,                                  cfmax=0.1_kp
!
!     For Tilting Angle Specification
!
      real(kind=kind_phys)       :: REFP(6), REFR(6), TLAC(8), PLAC(8), &
                                    TLBPL(7), drdp(5)
!
      DATA PLAC/100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0/
      DATA TLAC/ 35.0,  25.0,  20.0,  17.5,  15.0,  12.5,  10.0,  7.5/
      DATA REFP/500.0, 300.0, 250.0, 200.0, 150.0, 100.0/
      DATA REFR/  1.0,   2.0,   3.0,   4.0,   6.0,   8.0/
!
      real(kind=kind_phys)       :: AC(16), AD(16)
!
      integer, parameter :: nqrp=500001
      real(kind=kind_phys)       ::  C1XQRP, C2XQRP, TBQRP(NQRP),       &
                                     TBQRA(NQRP), TBQRB(NQRP)
!
      integer, parameter :: nvtp=10001
      real(kind=kind_phys)       ::  C1XVTP, C2XVTP, TBVTP(NVTP)
!
      real(kind=kind_phys)       ::  afc, facdt,                        &
                            grav, cp, alhl, alhf, rgas, rkap, nu, pi,   &
                            t0c,  rv, cvap, cliq, csol, ttp, eps, epsm1,&
!
                            ONEBG,   GRAVCON,  onebcp, GRAVFAC, ELOCP,  &
                            ELFOCP,  oneoalhl, CMPOR,  picon,   zfac,   &
                            deg2rad, PIINV,    testmboalhl,             &
                            rvi,     facw,     faci, hsub, tmix, DEN

      contains
 
! -----------------------------------------------------------------------
! CCPP entry points for gfdl cloud microphysics
! -----------------------------------------------------------------------

!>\brief The subroutine initializes rascnv
!!
!> \section arg_table_rascnv_init Argument Table
!! \htmlinclude rascnv_init.html
!!
      subroutine rascnv_init(me, dt,   con_g,    con_cp,   con_rd,      &
                             con_rv,   con_hvap, con_hfus, con_fvirt,   &
                             con_t0c,  con_ttp,  con_cvap, con_cliq,    &
                             con_csol, con_eps,  con_epsm1,             &
                             errmsg,   errflg)
!
      Implicit none
!
      integer,              intent(in)  :: me
      real(kind=kind_phys), intent(in)  :: dt,                          &
                      con_g,    con_cp,    con_rd,  con_rv,   con_hvap, &
                      con_hfus, con_fvirt, con_t0c, con_cvap, con_cliq, &
                      con_csol, con_ttp,   con_eps, con_epsm1

      character(len=*),     intent(out) :: errmsg
      integer,              intent(out) :: errflg
!
      real(kind=kind_phys), parameter ::  actp=1.7_kp, facm=1.00_kp
!
      real(kind=kind_phys) :: PH(15), A(15)
!
      DATA PH/150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0    &
     &,       550.0, 600.0, 650.0, 700.0, 750.0, 800.0, 850.0/
!
       DATA A/ 1.6851, 1.1686, 0.7663, 0.5255, 0.4100, 0.3677           &
     &,        0.3151, 0.2216, 0.1521, 0.1082, 0.0750, 0.0664           &
     &,        0.0553, 0.0445, 0.0633/
!
      real(kind=kind_phys) tem,  actop, tem1, tem2
      integer              i, l
!
! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0
      if (is_initialized) return
!                           set critical workfunction arrays
      ACTOP = ACTP*FACM
      DO L=1,15
        A(L) = A(L)*FACM
      ENDDO
      DO L=2,15
        TEM   = one / (PH(L) - PH(L-1))
        AC(L) = (PH(L)*A(L-1) - PH(L-1)*A(L)) * TEM
        AD(L) = (A(L) - A(L-1)) * TEM
      ENDDO
      AC(1)  = ACTOP
      AC(16) = A(15)
      AD(1)  = zero
      AD(16) = zero
!
      CALL SETQRP
      CALL SETVTP
!
      do i=1,7
        tlbpl(i) = (tlac(i)-tlac(i+1)) / (plac(i)-plac(i+1))
      enddo
      do i=1,5
        drdp(i)  = (REFR(i+1)-REFR(i)) / (REFP(i+1)-REFP(i))
      enddo
!
!     VTP = 36.34*SQRT(1.2)* (0.001)**0.1364
!
      AFC = -(1.01097e-4_kp*DT)*(3600.0_kp/DT)**0.57777778_kp
!
      if (fix_ncld_hr) then
        facdt = delt_c / dt
      else
        facdt = one / 3600.0_kp
      endif
!
      grav = con_g     ; cp    = con_cp   ; alhl = con_hvap
      alhf = con_hfus  ; rgas  = con_rd
      nu   = con_FVirt ; t0c   = con_t0c
      rv   = con_rv    ; cvap  = con_cvap
      cliq = con_cliq  ; csol  = con_csol ; ttp  = con_ttp
      eps  = con_eps   ; epsm1 = con_epsm1
!
      pi       = four*atan(one) ; PIINV   = one/PI
      ONEBG    = ONE / GRAV     ; GRAVCON = cmb2pa * ONEBG
      onebcp   = one / cp       ; GRAVFAC = GRAV / CMB2PA
      rkap     = rgas * onebcp  ; deg2rad = pi/180.0_kp
      ELOCP    = ALHL * onebcp  ; ELFOCP  = (ALHL+ALHF) * onebcp
      oneoalhl = one/alhl       ; CMPOR   = CMB2PA / RGAS
      picon    = half*pi*onebg  ; zfac    = 0.28888889e-4_kp * ONEBG
      testmboalhl = testmb/alhl
!
      rvi  = one / rv       ; facw = CVAP - CLIQ 
      faci = CVAP - CSOL    ; hsub = alhl + alhf
      tmix = TTP  - 20.0_kp ; DEN  = one / (TTP-TMIX)
!

      if (me == 0) write(0,*) ' NO DOWNDRAFT FOR CLOUD TYPES'           &
     &,                ' DETRAINING AT NORMALIZED PRESSURE ABOVE ',DPD
!
      is_initialized = .true.

!
      end subroutine rascnv_init
!
!! \section arg_table_rascnv_finalize Argument Table
!! \htmlinclude rascnv_finalize.html
!!
      subroutine rascnv_finalize (errmsg, errflg)

      implicit none

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      end subroutine rascnv_finalize
!!
!!
!!===================================================================== !
!! rascnv_run:                                                          !
!!                                                                      !
!! program history log:                                                 !
!!   Oct  2019  -- shrinivas moorthi                                    !
!!                                                                      !
!!                                                                      !
!! ====================  defination of variables  ====================
!! !
!!                                                                      !
!! inputs:                                                       size
!! !
!!    im       - integer, horiz dimension and num of used pts      1    !
!!    k        - integer, vertical dimension                       1    !
!!    dt       - real, time step in seconds                        1    !
!!    dtf      - real, dynamics time step in seconds               1    !
!!    rannum   - real, array holding random numbers between 0 an 1 (im,nrcm)  !
!!    tin      - real, input temperature (K)
!!    qin      - real, input specific humidity (kg/kg)
!!    uin      - real, input zonal wind component
!!    vin      - real, input meridional wind component
!!    ccin     - real, input condensates+tracers 
!!    fscav    - real 
!!    prsi     - real, layer interface pressure
!!    prsl     - real, layer mid pressure
!!    prsik    - real, layer interface Exner function
!!    prslk    - real, layer mid Exner function
!!    phil     - real, layer mid geopotential height
!!    phii     - real, layer interface geopotential height
!!    kpbl     - integer pbl top index
!!    cdrag    - real, drag coefficient
!!    rainc    - real, convectinve rain (m/sec)
!!    kbot     - integer, cloud bottom index
!!    ktop     - integer, cloud top index
!!    knv      - integer, 0 - no convvection; 1 - convection
!!    ddvel    - downdraft induced surface wind
!!    flipv    - logical, true if input data from bottom to top
!!    me       - integer, current pe number
!!    area     - real, grid area
!!    ccwf     - real, multiplication factor for critical workfunction
!!    nrcm     - integer, number of random numbers at each grid point
!!    rhc      - real, critical relative humidity
!!    ud_mf    - real, updraft mass flux
!!    dd_mf    - real, downdraft mass flux
!!    dt_mf    - real, detrained mass flux
!!    qw0      - real, min cloud water before autoconversion
!!    qi0      - real, min cloud ice before autoconversion
!!    dlqfac   - real,fraction of condensated detrained in layers
!!    kdt      - integer, current teime step
!!    revap    - logial,  when true reevaporate falling rain/snow
!!    qlcn     - real
!!    qicn     - real
!!    w_upi    - real
!!    cf_upi   - real
!!    cnv_mfd  - real
!!    cnv_dqldt- real
!!    clcn     - real
!!    cnv_fice - real
!!    cnv_ndrop- real
!!    cnv_nice - real
!!    mp_phys  - integer, microphysics option
!!    mp_phys_mg - integer, flag for MG microphysics option
!!    trcmin   - real, floor value for tracers
!!    ntk      - integer, index representing TKE in the tracer array
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! \section arg_table_rascnv_run Argument Table
!! \htmlinclude rascnv_run.html
!!
      subroutine rascnv_run(IM,     k,      ntr,   dt,   dtf            &
     &,                     ccwf,   area,   dxmin, dxinv                &
     &,                     psauras, prauras, wminras, dlqf, flipv      &
     &,                     me,     rannum, nrcm,  mp_phys, mp_phys_mg  &
     &,                     ntk,    kdt,                            rhc &
     &,                     tin,    qin,    uin,   vin,     ccin, fscav &
     &,                     prsi,   prsl,   prsik, prslk,   phil, phii  &
     &,                     KPBL,   CDRAG,  RAINC, kbot,    ktop, kcnv  &
     &,                     DDVEL,  ud_mf,  dd_mf, dt_mf                &
     &,                     QLCN,   QICN,   w_upi,  cf_upi, CNV_MFD     &
     &,                     CNV_DQLDT,CLCN,CNV_FICE,CNV_NDROP,CNV_NICE  &
     &,                     errmsg,  errflg)
!
!*********************************************************************
!*********************************************************************
!************         Relaxed Arakawa-Schubert      ******************
!************             Parameterization          ******************
!************          Plug Compatible Driver       ******************
!************               23 May 2002             ******************
!************                                       ******************
!************               Developed By            ******************
!************                                       ******************
!************             Shrinivas Moorthi         ******************
!************                                       ******************
!************                  EMC/NCEP             ******************
!*********************************************************************
!*********************************************************************
!
!
      Implicit none
!
      LOGICAL FLIPV
!
!      input
!
      integer, intent(in) :: im, k, ntr, me, nrcm, ntk, kdt             &
     &,                      mp_phys, mp_phys_mg
      integer, dimension(im) :: kbot, ktop, kcnv, kpbl
!
      real(kind=kind_phys), intent(in)        :: dxmin, dxinv, ccwf(2)  &
     &,                                          psauras(2), prauras(2) &
     &,                                          wminras(2), dlqf(2)
!
      real(kind=kind_phys), dimension(im,k+1) :: prsi, prsik, phii

      real(kind=kind_phys), dimension(im,k)   :: tin, qin,  uin, vin    &
     &,                                          prsl, prslk, phil      &

     &,                                          ud_mf, dd_mf, dt_mf    &
     &,                                          rhc, qlcn, qicn, w_upi &
     &,                                          cnv_mfd                &
     &,                                          cnv_dqldt, clcn        &
     &,                                          cnv_fice, cnv_ndrop    &
     &,                                          cnv_nice, cf_upi
      real(kind=kind_phys), dimension(im)     :: area,  cdrag           &
     &,                                          rainc, ddvel
      real(kind=kind_phys), dimension(im,nrcm):: rannum
      real(kind=kind_phys)                       ccin(im,k,ntr+2)
      real(kind=kind_phys)                       trcmin(ntr+2)

      real(kind=kind_phys) DT, dtf
!
!     Added for aerosol scavenging for GOCART
!
      real(kind=kind_phys), intent(in)  :: fscav(ntr)

!    &,                                   ctei_r(im), ctei_rm
      character(len=*),     intent(out) :: errmsg
      integer,              intent(out) :: errflg
!
!     locals
!
      real(kind=kind_phys), dimension(k)   :: toi,    qoi, tcu, qcu     &
     &,                                       pcu,    clw, cli, qii, qli&
     &,                                       phi_l,  prsm,psjm         &
     &,                                       alfinq, alfind, rhc_l     &
     &,                                       qoi_l, qli_l, qii_l
      real(kind=kind_phys), dimension(k+1) :: prs, psj, phi_h, flx, flxd


      integer, dimension(100)         :: ic
      real(kind=kind_phys), parameter :: clwmin=1.0e-10_kp
!
      real(kind=kind_phys), allocatable ::  ALFINT(:,:), uvi(:,:)       &
     &,                                     trcfac(:,:), rcu(:,:)
      real(kind=kind_phys)                  dtvd(2,4)
!    &,                                     DPI(K)
      real(kind=kind_phys) CFAC, TEM,  sgc, ccwfac, tem1, tem2, rain    &
     &,                    wfnc,tla,pl,qiid,qlid, c0, c0i, dlq_fac, sumq&
     &,                    rainp
!     integer                          :: nrcmax    ! Maximum # of random clouds per 1200s
!
      Integer              KCR,  KFX, NCMX, NC,  KTEM, I,   ii, Lm1, l  &
     &,                    ntrc,      ll,   km1, kp1,  ipt, lv, KBL, n  &
     &,                    KRMIN, KRMAX, KFMAX, kblmx, irnd,ib          &
     &,                    kblmn, ksfc, ncrnd
      real(kind=kind_phys) sgcs(k)
!
!  Scavenging related parameters
!
      real                fscav_(ntr+2)  ! Fraction scavenged per km
!
      fscav_ = -999.0_kp                 ! By default no scavenging
      if (ntr > 0 .and. fscav(1) > zero) then
        do i=1,ntr
          fscav_(i) = fscav(i)
        enddo
      endif
      trcmin = -99999.0_kp
      if (ntk-2 > 0) trcmin(ntk-2) = 1.0e-4_kp

!> - Initialize CCPP error handling variables

      errmsg = ''
      errflg = 0
!
      km1 = k - 1
      kp1 = k + 1
      if (flipv) then
        ksfc = 1
      else
        ksfc = kp1
      endif
!
      ntrc = ntr
      IF (CUMFRC) THEN
        ntrc = ntrc + 2
      ENDIF
      if (ntrc > 0) then
        if (.not. allocated(trcfac)) allocate (trcfac(k,ntrc))
        if (.not. allocated(uvi))    allocate (uvi(k,ntrc))
        if (.not. allocated(rcu))    allocate (rcu(k,ntrc))
        do n=1, ntrc
          do l=1,k
            trcfac(l,n) = one         !  For other tracers
            rcu(l,n)    = zero
          enddo
        enddo
      endif
!
!!!!! initialization for microphysics ACheng
      if(mp_phys == mp_phys_mg) then
        do l=1,K
          do i=1,im
            QLCN(i,l)      = zero
            QICN(i,l)      = zero
            w_upi(i,l)     = zero
            cf_upi(i,l)    = zero
            CNV_MFD(i,l)   = zero
!           CNV_PRC3(i,l)  = zero
            CNV_DQLDT(i,l) = zero
            CLCN(i,l)      = zero
            CNV_FICE(i,l)  = zero
            CNV_NDROP(i,l) = zero
            CNV_NICE(i,l)  = zero
          enddo
        enddo
      endif
!
      if (.not. allocated(alfint)) allocate(alfint(k,ntrc+4))
!
!     call set_ras_afc(dt)
!     AFC = -(1.04E-4*DT)*(3600./DT)**0.578
!     AFC = -(1.01097E-4*DT)*(3600./DT)**0.57777778
!
      do l=1,k
        do i=1,im
          ud_mf(i,l) = zero
          dd_mf(i,l) = zero
          dt_mf(i,l) = zero
        enddo
      enddo
      DO IPT=1,IM

        tem1    = max(zero, min(one, (log(area(ipt)) - dxmin) * dxinv))
        tem2    = one - tem1
        ccwfac  = ccwf(1)*tem1 + ccwf(2)*tem2
        dlq_fac = dlqf(1)*tem1 + dlqf(2)*tem2
        tem     = one  + dlq_fac
        c0i     = (psauras(1)*tem1 + psauras(2)*tem2) * tem
        c0      = (prauras(1)*tem1 + prauras(2)*tem2) * tem
        if (ccwfac == zero) ccwfac = half
!
!       ctei = .false.
!       if (ctei_r(ipt) > ctei_rm) ctei = .true.
!
!     Compute NCRND  : 
!                      if flipv is true, then input variables are from bottom
!                      to top while RAS goes top to bottom
!
        tem = one / prsi(ipt,ksfc)

        KRMIN = 1
        KRMAX = km1
        KFMAX = KRMAX
        kblmx = 1
        kblmn = 1
        sgcs(k) = one
        DO L=1,KM1
          ll = l
          if (flipv) ll = kp1 -l ! Input variables are bottom to top!
          SGC = prsl(ipt,ll) * tem
          sgcs(l) = sgc
          IF (SGC <= 0.050_kp) KRMIN = L
!         IF (SGC <= 0.700_kp) KRMAX = L
!         IF (SGC <= 0.800_kp) KRMAX = L
          IF (SGC <= 0.760_kp) KRMAX = L
!         IF (SGC <= 0.930_kp) KFMAX = L
          IF (SGC <= 0.970_kp) KFMAX = L    ! Commented on 20060202
!         IF (SGC <= 0.700_kp) kblmx = L    ! Commented on 20101015
          IF (SGC <= 0.600_kp) kblmx = L    ! 
!         IF (SGC <= 0.650_kp) kblmx = L    ! Commented on 20060202
          IF (SGC <= 0.980_kp) kblmn = L    ! 
        ENDDO
        krmin = max(krmin,2)

!
        if (fix_ncld_hr) then
!!!       NCRND = min(nrcmax, (KRMAX-KRMIN+1)) * (DTF/1200) + 0.50001
          NCRND = min(nrcmax, (KRMAX-KRMIN+1)) * (DTF/1800) + 0.10001_kp
!         NCRND = min(nrcmax, (KRMAX-KRMIN+1)) * (DTF/1200) + 0.10001
!         NCRND = min(nrcmax, (KRMAX-KRMIN+1)) * (DTF/900) + 0.50001
!         NCRND = min(nrcmax, (KRMAX-KRMIN+1)) * (DTF/600) + 0.50001
!         NCRND = min(nrcmax, (KRMAX-KRMIN+1)) * (DTF/360) + 0.50001
!    &                                         + 0.50001
!         NCRND = min(nrcmax, (KRMAX-KRMIN+1)) * min(1.0,DTF/360) + 0.1
        else
          NCRND = min(nrcmax, (KRMAX-KRMIN+1))
        endif
        NCRND   = min(nrcm,max(NCRND, 1))
!
        KCR     = MIN(K,KRMAX)
        KTEM    = MIN(K,KFMAX)
        KFX     = KTEM - KCR

        IF (KFX > 0) THEN
          IF (BOTOP) THEN
            DO NC=1,KFX
              IC(NC) = KTEM + 1 - NC
            ENDDO
          ELSE
            DO NC=KFX,1,-1
              IC(NC) = KTEM + 1 - NC
            ENDDO
          ENDIF
        ENDIF
!
        NCMX  = KFX + NCRND
        IF (NCRND > 0) THEN
          DO I=1,NCRND
            II = mod(i-1,nrcm) + 1
            IRND = (RANNUM(ipt,II)-0.0005_kp)*(KCR-KRMIN+1)
            IC(KFX+I) = IRND + KRMIN
          ENDDO
        ENDIF
!
        do l=1,k
          CLW(l)     = zero
          CLI(l)     = zero
                                 ! to be zero i.e. no environmental condensate!!!
          QII(l)     = zero
          QLI(l)     = zero
!                          Initialize heating, drying, cloudiness etc.
          tcu(l)     = zero
          qcu(l)     = zero
          pcu(l)     = zero
          flx(l)     = zero
          flxd(l)    = zero
          do n=1,ntrc
            rcu(l,n) = zero
          enddo
        enddo
        flx(kp1)  = zero
        flxd(kp1) = zero
        rain      = zero
!
        if (flipv) then   ! Input variables are bottom to top!
          do l=1,k
            ll = kp1 - l
                          ! Transfer input prognostic data into local variable
            toi(l)   = tin(ipt,ll)
            qoi(l)   = qin(ipt,ll)

            PRSM(L)  = prsl(ipt,ll) * facmb
            PSJM(L)  = prslk(ipt,ll)
            phi_l(L) = phil(ipt,ll)
            rhc_l(L) = rhc(ipt,ll)
!
            if (ntrc > ntr) then               ! CUMFRC is true
              uvi(l,ntr+1) = uin(ipt,ll)
              uvi(l,ntr+2) = vin(ipt,ll)
            endif
!
            if (ntr > 0) then                  ! tracers such as O3, dust etc
              do n=1,ntr
                uvi(l,n) = ccin(ipt,ll,n+2)
                if (abs(uvi(l,n)) < 1.0e-20_kp) uvi(l,n) = zero
              enddo
            endif
          enddo
          do l=1,kp1
            ll = kp1 + 1 - l      ! Input variables are bottom to top!
            PRS(LL)   = prsi(ipt,L) * facmb
            PSJ(LL)   = prsik(ipt,L)
            phi_h(LL) = phii(ipt,L)
          enddo
!
          if (ccin(ipt,1,2) <= -998.0_kp) then  ! input ice/water are together 
            do l=1,k
              ll = kp1 -l
              tem = ccin(ipt,ll,1)                                      &
     &            * MAX(ZERO, MIN(ONE, (TCR-toi(L))*TCRF))
              ccin(ipt,ll,2) = ccin(ipt,ll,1) - tem
              ccin(ipt,ll,1) = tem
            enddo
          endif
          if (advcld) then
            do l=1,k
              ll = kp1 -l ! Input variables are bottom to top!
              QII(L) = ccin(ipt,ll,1)
              QLI(L) = ccin(ipt,ll,2)
            enddo
          endif
          KBL  = MAX(MIN(k, kp1-KPBL(ipt)), k/2)
!
        else              ! Input variables are top to bottom!

          do l=1,k
                          ! Transfer input prognostic data into local variable
            toi(l)   = tin(ipt,l)
            qoi(l)   = qin(ipt,l)

            PRSM(L)  = prsl(ipt, L) * facmb
            PSJM(L)  = prslk(ipt,L)
            phi_l(L) = phil(ipt,L)
            rhc_l(L) = rhc(ipt,L)
!
            if (ntrc > ntr) then               ! CUMFRC is true
              uvi(l,ntr+1) = uin(ipt,l)
              uvi(l,ntr+2) = vin(ipt,l)
            endif
!
            if (ntr > 0) then                  ! tracers such as O3, dust etc
              do n=1,ntr
                uvi(l,n) = ccin(ipt,l,n+2)
                if (abs(uvi(l,n)) < 1.0e-20_kp) uvi(l,n) = zero
              enddo
            endif
          enddo
          DO L=1,kp1
            PRS(L)   = prsi(ipt,L) * facmb
            PSJ(L)   = prsik(ipt,L)
            phi_h(L) = phii(ipt,L)
          ENDDO
!
          if (ccin(ipt,1,2) <= -998.0_kp) then  ! input ice/water are together
            do l=1,k
              tem = ccin(ipt,l,1)                                       &
     &            * MAX(ZERO, MIN(ONE, (TCR-toi(L))*TCRF))
              ccin(ipt,l,2) = ccin(ipt,l,1) - tem
              ccin(ipt,l,1) = tem
            enddo
          endif
          if (advcld) then
            do l=1,k
              QII(L) = ccin(ipt,l,1)
              QLI(L) = ccin(ipt,l,2)
            enddo
          endif
!
          KBL  = KPBL(ipt)
!
        endif      ! end of if (flipv) then
!
!       do l=k,kctop(1),-1
!!        DPI(L)  = 1.0 / (PRS(L+1) - PRS(L))
!       enddo
!
!     print *,' ipt=',ipt

        if (advups) then               ! For first order upstream for updraft
          alfint(:,:) = one
        elseif (advtvd) then           ! TVD flux limiter scheme for updraft
          alfint(:,:) = one
          l   = krmin
          lm1 = l - 1
          dtvd(1,1) = cp*(toi(l)-toi(lm1)) + phi_l(l)-phi_l(lm1)        &
     &              + alhl*(qoi(l)-qoi(lm1))
          dtvd(1,2) = qoi(l) - qoi(lm1)
          dtvd(1,3) = qli(l) - qli(lm1)
          dtvd(1,4) = qii(l) - qii(lm1)
          do l=krmin+1,k
            lm1 = l - 1

!     write(0,*)' toi=',toi(l),toi(lm1),' phi_l=',phi_l(l),phi_l(lm1)
!    &,' qoi=',qoi(l),qoi(lm1),' cp=',cp,' alhl=',alhl

            dtvd(2,1)   = cp*(toi(l)-toi(lm1)) + phi_l(l)-phi_l(lm1)    &
     &                  + alhl*(qoi(l)-qoi(lm1))

!     write(0,*)' l=',l,' dtvd=',dtvd(:,1)

            if (abs(dtvd(2,1)) > 1.0e-10_kp) then
              tem1        = dtvd(1,1) / dtvd(2,1)
              tem2        = abs(tem1)
              alfint(l,1) = one - half*(tem1 + tem2)/(one + tem2)   ! for h
            endif

!     write(0,*)' alfint=',alfint(l,1),' l=',l,' ipt=',ipt

            dtvd(1,1) = dtvd(2,1)
!
            dtvd(2,2) = qoi(l) - qoi(lm1)

!     write(0,*)' l=',l,' dtvd2=',dtvd(:,2)

            if (abs(dtvd(2,2)) > 1.0e-10_kp) then
              tem1        = dtvd(1,2) / dtvd(2,2)
              tem2        = abs(tem1)
              alfint(l,2) = one - half*(tem1 + tem2)/(one + tem2)   ! for q
            endif
            dtvd(1,2) = dtvd(2,2)
!
            dtvd(2,3) = qli(l) - qli(lm1)

!     write(0,*)' l=',l,' dtvd3=',dtvd(:,3)

            if (abs(dtvd(2,3)) > 1.0e-10_kp) then
              tem1        = dtvd(1,3) / dtvd(2,3)
              tem2        = abs(tem1)
              alfint(l,3) = one - half*(tem1 + tem2)/(one + tem2)   ! for ql
            endif
            dtvd(1,3) = dtvd(2,3)
!
            dtvd(2,4) = qii(l) - qii(lm1)

!     write(0,*)' l=',l,' dtvd4=',dtvd(:,4)

            if (abs(dtvd(2,4)) > 1.0e-10_kp) then
              tem1        = dtvd(1,4) / dtvd(2,4)
              tem2        = abs(tem1)
              alfint(l,4) = one - half*(tem1 + tem2)/(one + tem2)   ! for qi
            endif
            dtvd(1,4) = dtvd(2,4)
          enddo
!
          if (ntrc > 0) then
            do n=1,ntrc
              l = krmin
              dtvd(1,1) = uvi(l,n) - uvi(l-1,n)
              do l=krmin+1,k
                dtvd(2,1) = uvi(l,n) - uvi(l-1,n)

!     write(0,*)' l=',l,' dtvdn=',dtvd(:,1),' n=',n,' l=',l

                if (abs(dtvd(2,1)) > 1.0e-10_kp) then
                  tem1          = dtvd(1,1) / dtvd(2,1)
                  tem2          = abs(tem1)
                  alfint(l,n+4) = one - half*(tem1 + tem2)/(one + tem2) ! for tracers
                endif
                dtvd(1,1) = dtvd(2,1)
              enddo
            enddo
          endif
        else
          alfint(:,:) = half              ! For second order scheme
        endif
        alfind(:)   = half
!
!     write(0,*)' after alfint for ipt=',ipt

! Resolution dependent press grad correction momentum mixing

        if (CUMFRC) then
          do l=krmin,k
            tem = one - max(pgfbot, min(pgftop, pgftop+pgfgrad*prsm(l)))
            trcfac(l,ntr+1) = tem
            trcfac(l,ntr+2) = tem
          enddo
        endif
!
!       if (calkbl) kbl = k

        if (calkbl) then
          kbl = kblmn
        else
          kbl = min(kbl, kblmn)
        endif
!
        DO NC=1,NCMX     ! multi cloud loop
!
          IB = IC(NC)    ! cloud top level index
          if (ib > kbl-1) cycle
!
!****************************************************************************
!         if (advtvd) then           ! TVD flux limiter scheme for updraft
!           l   = ib
!           lm1 = l - 1
!           dtvd(1,1) = cp*(toi(l)-toi(lm1)) + phi_l(l)-phi_l(lm1)
!    &                + alhl*(qoi(l)-qoi(lm1))
!           dtvd(1,2) = qoi(l) - qoi(lm1)
!           dtvd(1,3) = qli(l) - qli(lm1)
!           dtvd(1,4) = qii(l) - qii(lm1)
!           do l=ib+1,k
!             lm1 = l - 1
!             dtvd(2,1)   = cp*(toi(l)-toi(lm1)) + phi_l(l)-phi_l(lm1)
!    &                    + alhl*(qoi(l)-qoi(lm1))
!             if (abs(dtvd(2,1)) > 1.0e-10) then
!               tem1        = dtvd(1,1) / dtvd(2,1)
!               tem2        = abs(tem1)
!               alfint(l,1) = 1.0 - 0.5*(tem1 + tem2)/(1.0 + tem2)   ! for h
!             endif
!             dtvd(1,1)   = dtvd(2,1)
!
!             dtvd(2,2)   = qoi(l) - qoi(lm1)
!             if (abs(dtvd(2,2)) > 1.0e-10) then
!               tem1        = dtvd(1,2) / dtvd(2,2)
!               tem2        = abs(tem1)
!               alfint(l,2) = 1.0 - 0.5*(tem1 + tem2)/(1.0 + tem2)   ! for q
!             endif
!             dtvd(1,2)   = dtvd(2,2)
!
!             dtvd(2,3)   = qli(l) - qli(lm1)
!             if (abs(dtvd(2,3)) > 1.0e-10) then
!               tem1        = dtvd(1,3) / dtvd(2,3)
!               tem2        = abs(tem1)
!               alfint(l,3) = 1.0 - 0.5*(tem1 + tem2)/(1.0 + tem2)   ! for ql
!             endif
!             dtvd(1,3)   = dtvd(2,3)
!
!             dtvd(2,4)   = qii(l) - qii(lm1)
!             if (abs(dtvd(2,4)) > 1.0e-10) then
!               tem1        = dtvd(1,4) / dtvd(2,4)
!               tem2        = abs(tem1)
!               alfint(l,4) = 1.0 - 0.5*(tem1 + tem2)/(1.0 + tem2)   ! for qi
!             endif
!             dtvd(1,4)   = dtvd(2,4)
!           enddo
!
!           if (ntrc > 0) then
!             do n=1,ntrc
!               l = ib
!               dtvd(1,1)   = uvi(l,n) - uvi(l-1,n)
!               do l=ib+1,k
!                 dtvd(2,1)     = uvi(l,n) - uvi(l-1,n)
!                 if (abs(dtvd(2,1)) > 1.0e-10) then
!                   tem1        = dtvd(1,1) / dtvd(2,1)
!                   tem2          = abs(tem1)
!                   alfint(l,n+4) = 1.0 - 0.5*(tem1 + tem2)/(1.0 + tem2) ! for tracers
!                 endif
!                 dtvd(1,1)     = dtvd(2,1)
!               enddo
!             enddo
!           endif
!         endif
!****************************************************************************
!
          WFNC = zero
          do L=IB,KP1
            FLX(L)  = zero
            FLXD(L) = zero
          enddo
!
          TLA  = -10.0_kp
!
          qiid = qii(ib)         ! cloud top level ice before convection
          qlid = qli(ib)         ! cloud top level water before convection
!
          rainp = rain

          CALL CLOUD(K,  KP1, IB, ntrc, kblmx, kblmn                    &
     &,              FRAC,  MAX_NEG_BOUY, vsmooth, aw_scal              &
     &,              REVAP, WRKFUN, CALKBL, CRTFUN                      &
     &,              DT, KDT, TLA, DPD                                  &
     &,              ALFINT, rhfacl, rhfacs, area(ipt)                  &
     &,              ccwfac, CDRAG(ipt), trcfac                         &
     &,              alfind, rhc_l, phi_l, phi_h, PRS, PRSM,sgcs        &
     &,              TOI, QOI, UVI, QLI, QII, KBL, DDVEL(ipt)           &
     &,              TCU, QCU, RCU, PCU, FLX, FLXD, RAIN, WFNC, fscav_  &
     &,              trcmin, ntk-2, c0, wminras(1), c0i, wminras(2)     &
     &,              dlq_fac)
!    &,              ctei)

!
          if (flipv) then
            do L=IB,K
              ll = kp1 -l    ! Input variables are bottom to top!
              ud_mf(ipt,ll) = ud_mf(ipt,ll) + flx(l+1)
              dd_mf(ipt,ll) = dd_mf(ipt,ll) + flxd(l+1)
            enddo
            ll  = kp1 - ib
            dt_mf(ipt,ll) = dt_mf(ipt,ll) + flx(ib)

            if (mp_phys == mp_phys_mg) then !         Anning Cheng for microphysics 11/14/2015

              CNV_MFD(ipt,ll)   = CNV_MFD(ipt,ll)   + flx(ib)/dt

!!            CNV_DQLDT(ipt,ll) = CNV_DQLDT(ipt,ll)
!!   &                          + max(0.,(QLI(ib)+QII(ib)-qiid-qlid))/dt
              CNV_DQLDT(ipt,ll) = CNV_DQLDT(ipt,ll) + flx(ib)*          &
     &                            max(0.,(QLI(ib)+QII(ib)-qiid-qlid))/dt
!!   &                                max(0.,(QLI(ib)+QII(ib)))/dt/3.
              if(flx(ib)<0) write(*,*)"AAA666", flx(ib),QLI(ib),QII(ib) &
     &                                       ,ipt,ll
            endif

          else

            do L=IB,K
              ud_mf(ipt,l) = ud_mf(ipt,l) + flx(l+1)
              dd_mf(ipt,l) = dd_mf(ipt,l) + flxd(l+1)
            enddo
            dt_mf(ipt,ib) = dt_mf(ipt,ib) + flx(ib)

            if (mp_phys == mp_phys_mg) then !         Anning Cheng for microphysics 11/14/2015
              CNV_MFD(ipt,ib)   = CNV_MFD(ipt,ib)   + flx(ib)/dt
!!            CNV_DQLDT(ipt,ib) = CNV_DQLDT(ipt,ib)
!!   &                          + max(0.,(QLI(ib)+QII(ib)-qiid-qlid))/dt
              CNV_DQLDT(ipt,ib) = CNV_DQLDT(ipt,ib) + flx(ib)*          &
     &                            max(zero,(QLI(ib)+QII(ib)-qiid-qlid))/dt
!!   &                                max(0.,(QLI(ib)+QII(ib)))/dt/3.
              if(flx(ib)<0) write(*,*)"AAA666", flx(ib),QLI(ib),QII(ib) &
     &                                       ,ipt,ib
            endif

          endif
! 
!
!   Warning!!!!
!   ------------
!   By doing the following, CLOUD does not contain environmental
!   condensate!
!
          if (.not. advcld) then
            do l=1,K
              clw(l) = clw(l) + QLI(L)
              cli(l) = cli(l) + QII(L)
              QLI(L) = zero
              QII(L) = zero
            enddo
          endif
!
        ENDDO                        ! End of the NC loop!
!
        RAINC(ipt) = rain * 0.001_kp    ! Output rain is in meters
        if (rainc(ipt) < rain_min) rainc(ipt) = zero

        ktop(ipt) = kp1
        kbot(ipt) = 0

        kcnv(ipt) = 0


        do l=k,1,-1
!         qli(l) = max(qli(l), zero)
!         qii(l) = max(qii(l), zero)
!         clw(i) = max(clw(i), zero)
!         cli(i) = max(cli(i), zero)

          if (sgcs(l) < 0.93_kp .and. abs(tcu(l)) > one_m10) then
!         if (sgcs(l) < 0.90_kp .and. tcu(l) .ne. zero) then
!         if (sgcs(l) < 0.85_kp .and. tcu(l) .ne. zero) then
             kcnv(ipt) = 1
          endif
!  New test for convective clouds ! added in 08/21/96
          if (clw(l)+cli(l) > zero .OR.                                 &
     &        qli(l)+qii(l) > clwmin) ktop(ipt) = l
        enddo
        do l=1,km1
          if (clw(l)+cli(l) > zero .OR.                                 &
     &        qli(l)+qii(l) > clwmin) kbot(ipt) = l
        enddo
!
        if (flipv) then
          do l=1,k
            ll = kp1 - l
            tin(ipt,ll) = toi(l)                  ! Temperature
            qin(ipt,ll) = qoi(l)                  ! Specific humidity
            uin(ipt,ll) = uvi(l,ntr+1)            ! U momentum
            vin(ipt,ll) = uvi(l,ntr+2)            ! V momentum

!!        for 2M microphysics, always output these variables
            if (mp_phys == mp_phys_mg) then
              if (advcld) then
                QLCN(ipt,ll)     = max(qli(l)-ccin(ipt,ll,2), zero)
                QICN(ipt,ll)     = max(qii(l)-ccin(ipt,ll,1), zero)
                CNV_FICE(ipt,ll) = QICN(ipt,ll)                         &
     &                           / max(1.0e-10_kp,QLCN(ipt,ll)+QICN(ipt,ll))
              else
                QLCN(ipt,ll)     = qli(l)
                QICN(ipt,ll)     = qii(l)
                CNV_FICE(ipt,ll) = qii(l)/max(1.0e-10_kp,qii(l)+qli(l))
              endif
              cf_upi(ipt,ll)   = max(zero,min(0.02_kp*log(one+           &
     &                             500.0_kp*ud_mf(ipt,ll)/dt), cfmax))
!    &                             500*ud_mf(ipt,ll)/dt), 0.60))
              CLCN(ipt,ll)     = cf_upi(ipt,ll)  !downdraft is below updraft
              w_upi(ipt,ll)    = ud_mf(ipt,ll)*toi(l)*rgas /            &
     &                      (dt*max(cf_upi(ipt,ll),1.0e-12_kp)*prsl(ipt,ll))
            endif

            if (ntr > 0) then
              do n=1,ntr
                ccin(ipt,ll,n+2) = uvi(l,n)           ! Tracers
              enddo
            endif
          enddo
          if (advcld) then
            do l=1,k
              ll  = kp1 - l
              ccin(ipt,ll,1) = qii(l)          ! Cloud ice
              ccin(ipt,ll,2) = qli(l)          ! Cloud water
            enddo
          else
            do l=1,k
              ll  = kp1 - l
              ccin(ipt,ll,1) = ccin(ipt,ll,1) + cli(l)
              ccin(ipt,ll,2) = ccin(ipt,ll,2) + clw(l)
            enddo
          endif
!
          ktop(ipt) = kp1 - ktop(ipt)
          kbot(ipt) = kp1 - kbot(ipt)
!
        else

          do l=1,k
            tin(ipt,l) = toi(l)                   ! Temperature
            qin(ipt,l) = qoi(l)                   ! Specific humidity
            uin(ipt,l) = uvi(l,ntr+1)             ! U momentum
            vin(ipt,l) = uvi(l,ntr+2)             ! V momentum

!!        for 2M microphysics, always output these variables
            if (mp_phys == mp_phys_mg) then
              if (advcld) then
                QLCN(ipt,l)     = max(qli(l)-ccin(ipt,l,2), zero)
                QICN(ipt,l)     = max(qii(l)-ccin(ipt,l,1), zero)
                CNV_FICE(ipt,l) = QICN(ipt,l)                           &
     &                          / max(1.0e-10_kp,QLCN(ipt,l)+QICN(ipt,l))
              else
                QLCN(ipt,l)     = qli(l)
                QICN(ipt,l)     = qii(l)
                CNV_FICE(ipt,l) = qii(l)/max(1.0e-10_kp,qii(l)+qli(l))
              endif
!!            CNV_PRC3(ipt,l) = PCU(l)/dt
!             CNV_PRC3(ipt,l) = zero
!             if(PCU(l) < zero) write(*,*)"AAA777",PCU(l),ipt,l
              cf_upi(ipt,l)   = max(zero,min(0.02_kp*log(one+           &
     &                             500.0_kp*ud_mf(ipt,l)/dt), cfmax))
!    &                             500*ud_mf(ipt,l)/dt), 0.60))
              CLCN(ipt,l)     = cf_upi(ipt,l)  !downdraft is below updraft
              w_upi(ipt,l)    = ud_mf(ipt,l)*toi(l)*rgas /              &
     &                        (dt*max(cf_upi(ipt,l),1.0e-12_kp)*prsl(ipt,l))
            endif

            if (ntr > 0) then
              do n=1,ntr
                ccin(ipt,l,n+2) = uvi(l,n)           ! Tracers
              enddo
            endif
          enddo
          if (advcld) then
            do l=1,k
              ccin(ipt,l,1) = qii(l)          ! Cloud ice
              ccin(ipt,l,2) = qli(l)          ! Cloud water
            enddo
          else
            do l=1,k
              ccin(ipt,l,1) = ccin(ipt,l,1) + cli(l)
              ccin(ipt,l,2) = ccin(ipt,l,2) + clw(l)
            enddo
          endif
        endif
!
!     Velocity scale from the downdraft!
!
        DDVEL(ipt) = DDVEL(ipt) * DDFAC * GRAV / (prs(KP1)-prs(K))
!
      ENDDO                            ! End of the IPT Loop!

      deallocate (alfint, uvi, trcfac, rcu)
!
      RETURN
      end subroutine rascnv_run
      SUBROUTINE CLOUD(                                                 &
     &                  K, KP1, KD, NTRC, KBLMX, kblmn                  &
     &,                 FRACBL, MAX_NEG_BOUY, vsmooth, aw_scal          &
     &,                 REVAP, WRKFUN, CALKBL, CRTFUN                   &
     &,                 DT, KDT, TLA, DPD                               &
     &,                 ALFINT, RHFACL, RHFACS, area, ccwf, cd, trcfac  &
     &,                 alfind, rhc_ls, phil, phih, prs, prsm, sgcs     &
     &,                 TOI, QOI, ROI,  QLI, QII, KPBL, DSFC            &
     &,                 TCU, QCU, RCU, PCU, FLX, FLXD, CUP, WFNC,fscav_ &
     &,                 trcmin, ntk, c0, qw0, c0i, qi0, dlq_fac)
!    &,                 ctei)

!
!***********************************************************************
!******************** Relaxed  Arakawa-Schubert ************************
!****************** Plug Compatible Scalar Version *********************
!************************ SUBROUTINE CLOUD  ****************************
!************************  October 2004     ****************************
!********************  VERSION 2.0  (modified) *************************
!************* Shrinivas.Moorthi@noaa.gov (301) 683-3718  ***** ********
!***********************************************************************
!*References:
!-----------
!     NOAA Technical Report NWS/NCEP 99-01:
!     Documentation of Version 2 of Relaxed-Arakawa-Schubert
!     Cumulus Parameterization with Convective Downdrafts, June 1999.
!     by S. Moorthi and M. J. Suarez.
!
!     Relaxed Arakawa-Schubert Cumulus Parameterization (Version 2)
!     with Convective Downdrafts - Unpublished Manuscript (2002)
!     by Shrinivas Moorthi and Max J. Suarez.
!
!***********************************************************************
!
!===>    UPDATES CLOUD TENDENCIES DUE TO A SINGLE CLOUD
!===>    DETRAINING AT LEVEL KD.
!
!***********************************************************************
!
!===>  TOI(K)     INOUT   TEMPERATURE             KELVIN
!===>  QOI(K)     INOUT   SPECIFIC HUMIDITY       NON-DIMENSIONAL
!===>  ROI(K,NTRC)INOUT   TRACER                  ARBITRARY
!===>  QLI(K)     INOUT   LIQUID WATER            NON-DIMENSIONAL
!===>  QII(K)     INOUT   ICE                     NON-DIMENSIONAL

!===>  PRS(KP1)   INPUT   PRESSURE @ EDGES        MB
!===>  PRSM(K)    INPUT   PRESSURE @ LAYERS       MB
!===>  SGCS(K)    INPUT   Local sigma
!===>  PHIH(KP1)  INPUT   GEOPOTENTIAL @ EDGES  IN MKS units
!===>  PHIL(K)    INPUT   GEOPOTENTIAL @ LAYERS IN MKS units
!===>  PRJ(KP1)   INPUT   (P/P0)^KAPPA  @ EDGES   NON-DIMENSIONAL
!===>  PRJM(K)    INPUT   (P/P0)^KAPPA  @ LAYERS  NON-DIMENSIONAL

!===>  K          INPUT   THE RISE & THE INDEX OF THE SUBCLOUD LAYER
!===>  KD         INPUT   DETRAINMENT LEVEL ( 1<= KD < K )          
!===>  NTRC       INPUT   NUMBER OF TRACERS. MAY BE ZERO.
!===>  kblmx      INPUT   highest level the pbl can take
!===>  kblmn      INPUT   lowest  level the pbl can take
!===>  DPD        INPUT   Critical normalized pressure (i.e. sigma) at the cloud top
!                         No downdraft calculation if the cloud top pressure is higher
!                         than DPD*PRS(KP1)
!
!===>  TCU(K  )   UPDATE  TEMPERATURE TENDENCY       DEG
!===>  QCU(K  )   UPDATE  WATER VAPOR TENDENCY       (G/G)
!===>  RCU(K,NTRC)UPDATE  TRACER TENDENCIES          ND
!===>  PCU(K)     UPDATE  PRECIP @ BASE OF LAYER     KG/M^2
!===>  FLX(K  )   UPDATE  MASS FLUX @ TOP OF LAYER   KG/M^2
!===>  CUP        UPDATE  PRECIPITATION AT THE SURFACE KG/M^2
!
      IMPLICIT NONE
!
      real (kind=kind_phys), parameter :: RHMAX=1.0_kp                   & ! MAX RELATIVE HUMIDITY
     &,                                   QUAD_LAM=1.0_kp                & ! MASK FOR QUADRATIC LAMBDA
     &,                                   RHRAM=0.05_kp                  & ! PBL RELATIVE HUMIDITY RAMP
!    &,                                   RHRAM=0.15_kp                 !& ! PBL RELATIVE HUMIDITY RAMP
     &,                                   HCRITD=4000.0_kp               & ! Critical Moist Static Energy for Deep clouds
!    &,                                   HCRITS=2000.0_kp               & ! Critical Moist Static Energy for Shallow clouds
     &,                                   HCRITS=2500.0_kp               & ! Critical Moist Static Energy for Shallow clouds
     &,                                   pcrit_lcl=250.0_kp             & ! Critical pressure difference between boundary layer top
                                                                          ! layer top and lifting condensation level (hPa)
!    &,                                   hpert_fac=1.01_kp             !& ! Perturbation on hbl when ctei=.true.
!    &,                                   hpert_fac=1.005_kp            !& ! Perturbation on hbl when ctei=.true.
     &,                                   qudfac=quad_lam*half           &
     &,                                   shalfac=3.0_kp                 &
!    &,                                   qudfac=quad_lam*pt25, shalfac=3.0_kp       !& !  Yogesh's
!    &,                                   c0ifac=0.07_kp                 & ! following Han et al, 2016 MWR
!    &,                                   c0ifac=0.001_kp                & ! following Han et al, 2017 Weather and Forecasting
     &,                                   c0ifac=0.0_kp                  &
     &,                                   dpnegcr  = 150.0_kp
!    &,                                   dpnegcr  = 100.0_kp
!    &,                                   dpnegcr  = 200.0_kp
!
      real(kind=kind_phys), parameter :: ERRMIN=0.0001_kp                &
     &,                                  ERRMI2=0.1_kp*ERRMIN            &
!    &,                                  rainmin=1.0e-9_kp              !&
     &,                                  rainmin=1.0e-8_kp               &
     &,                                  oneopt9=one/0.09_kp             &
     &,                                  oneopt4=one/0.04_kp
      real(kind=kind_phys), parameter :: almax=1.0e-2_kp                 &
     &,                                  almin1=0.0_kp,   almin2=0.0_kp
      real(kind=kind_phys), parameter :: bldmax=300.0_kp, bldmin=25.0_kp
!
!  INPUT ARGUMENTS

!     LOGICAL REVAP, WRKFUN, CALKBL, CRTFUN, CALCUP, ctei
      LOGICAL REVAP, WRKFUN, CALKBL, CRTFUN, CALCUP
      logical vsmooth, aw_scal
      INTEGER K, KP1, KD, NTRC, kblmx, kblmn, ntk


      real(kind=kind_phys), dimension(K)   ::  TOI,  QOI, PRSM, QLI, QII&
     &,                                        PHIL, SGCS, rhc_ls       &
     &,                                        alfind
      real(kind=kind_phys), dimension(KP1) ::  PRS,  PHIH
      real(kind=kind_phys), dimension(K,NTRC) :: ROI, trcfac
      real(kind=kind_phys), dimension(ntrc)   :: trcmin
      real(kind=kind_phys)                    :: CD, DSFC
      INTEGER                                 :: KPBL, KBL, KB1, kdt

      real(kind=kind_phys) ALFINT(K,NTRC+4)
      real(kind=kind_phys) FRACBL, MAX_NEG_BOUY, DPD                    &
     &,                    RHFACL, RHFACS, area, ccwf                   &
     &,                    c0, qw0, c0i, qi0, dlq_fac
 
!  UPDATE ARGUMENTS

      real(kind=kind_phys), dimension(K)      :: TCU, QCU, TCD, QCD, PCU
      real(kind=kind_phys), dimension(KP1)    :: FLX, FLXD
      real(kind=kind_phys), dimension(K,NTRC) :: RCU
      real(kind=kind_phys)                    :: CUP
!
!  TEMPORARY WORK SPACE

      real(kind=kind_phys), dimension(KD:K)   :: HOL, QOL, HST, QST     &
     &,                           TOL, GMH, AKT, AKC, BKC, LTL, RNN     &
     &,                           FCO, PRI, QIL, QLL, ZET, XI, RNS      &
     &,                           Q0U, Q0D, vtf, CIL, CLL, ETAI, dlq    &
     &,                           wrk1, wrk2, dhdp, qrb, qrt, evp       &
     &,                           ghd, gsd, etz, cldfr,  sigf, rho

      real(kind=kind_phys), dimension(KD:KP1) :: GAF, GMS, GAM, DLB     &
     &,                           DLT, ETA, PRL, BUY, ETD, HOD, QOD, wvl
      real(kind=kind_phys), dimension(KD:K-1) :: etzi

      real(kind=kind_phys) fscav_(ntrc)

      LOGICAL ep_wfn, cnvflg, LOWEST, DDFT, UPDRET

      real(kind=kind_phys) ALM,   DET,    HCC,  CLP                     &
     &,                    HSU,   HSD,    QTL,  QTV                     &
     &,                    AKM,   WFN,    HOS,  QOS                     &
     &,                    AMB,   TX1,    TX2,  TX3                     &
     &,                    TX4,   TX5,    QIS,  QLS                     &
     &,                    HBL,   QBL,    RBL(NTRC), wcbase             &
     &,                    QLB,   QIB,    PRIS                          &
     &,                    WFNC,  TX6,    ACR                           &
     &,                    TX7,   TX8,    TX9,  RHC                     &
     &,                    hstkd, qstkd,  ltlkd, q0ukd, q0dkd, dlbkd    &
     &,                    qtp, qw00, qi00, qrbkd                       &
     &,                    hstold, rel_fac, prism                       &
     &,                    TL, PL, QL, QS, DQS, ST1, SGN, TAU,          &
     &                     QTVP, HB, QB, TB, QQQ,                       &
     &                     HCCP, DS, DH, AMBMAX, X00, EPP, QTLP,        &
     &                     DPI, DPHIB, DPHIT, DEL_ETA, DETP,            &
     &                     TEM, TEM1, TEM2, TEM3, TEM4,                 &
     &                     ST2, ST3, ST4, ST5,                          &
     &                     ERRH, ERRW, ERRE, TEM5,                      &
     &                     TEM6, HBD, QBD, st1s, shal_fac, hmax, hmin,  &
     &                     dhdpmn, avt, avq, avr, avh                   &
     &,                    TRAIN, DOF, CLDFRD, tla, gmf                 &
     &,                    FAC, RSUM1, RSUM2, RSUM3, dpneg, hcrit       &
     &,                    ACTEVAP,AREARAT,DELTAQ,MASS,MASSINV,POTEVAP  &
     &,                    TEQ,QSTEQ,DQDT,QEQ                           &
     &,                    CLFRAC, DT,      clvfr, delzkm, fnoscav, delp
!    &,                    CLFRAC, DT, clf, clvfr, delzkm, fnoscav, delp
!    &,                    almin1, almin2

      INTEGER I, L,  N,  KD1, II, iwk, idh, lcon                        &
     &,       IT, KM1, KTEM, KK, KK1, LM1, LL, LP1, kbls, kmxh          &
     &,       kblh, kblm, kblpmn, kmax, kmaxm1, kmaxp1, klcl, kmin, kmxb
!
!***********************************************************************
!
!     almin2 = 0.2 * sqrt(pi/area)
!     almin1 = almin2

      KM1 = K  - 1
      KD1 = KD + 1

      do l=1,K
        tcd(L) = zero
        qcd(L) = zero
      enddo
!
      CLDFRD   = zero
      DOF      = zero
      PRL(KP1) = PRS(KP1)
!
      DO L=KD,K
        RNN(L)   = zero
        ZET(L)   = zero
        XI(L)    = zero
!
        TOL(L)   = TOI(L)
        QOL(L)   = QOI(L)
        PRL(L)   = PRS(L)
        CLL(L)   = QLI(L)
        CIL(L)   = QII(L)
        BUY(L)   = zero

        wvl(l)   = zero
      ENDDO
      wvl(kp1) = zero
!
      if (vsmooth) then
        do l=kd,k
          wrk1(l) = tol(l)
          wrk2(l) = qol(l)
        enddo
        do l=kd1,km1
          tol(l) = pt25*wrk1(l-1) + half*wrk1(l) + pt25*wrk1(l+1)
          qol(l) = pt25*wrk2(l-1) + half*wrk2(l) + pt25*wrk2(l+1)
        enddo
      endif
!
      DO L=KD, K
        DPI    = ONE / (PRL(L+1) - PRL(L))
        PRI(L) = GRAVFAC * DPI
!
        PL     = PRSM(L)
        TL     = TOL(L)

        rho(l) = cmb2pa * pl / (rgas*tl*(one+nu*qol(l)))

        AKT(L) = (PRL(L+1) - PL) * DPI
!
        CALL QSATCN(TL, PL, QS, DQS)
!
        QST(L) = QS
        GAM(L) = DQS * ELOCP
        ST1    = ONE + GAM(L)
        GAF(L) = ONEOALHL * GAM(L) / ST1
 
        QL     = MAX(MIN(QS*RHMAX,QOL(L)), ONE_M10)
        QOL(L) = QL
 
        TEM    = CP * TL
        LTL(L) = TEM * ST1 / (ONE+NU*(QST(L)+TL*DQS))
        vtf(L) = one + NU * QL
        ETA(L) = ONE / (LTL(L) * VTF(L))

        HOL(L) = TEM + QL * ALHL
        HST(L) = TEM + QS * ALHL
!
      ENDDO
!
      ETA(KP1) = ZERO
      GMS(K)   = ZERO
!
      AKT(KD)  = HALF
      GMS(KD)  = ZERO
!
      CLP      = ZERO
!
      GAM(KP1) = GAM(K)
      GAF(KP1) = GAF(K)
!
      DO L=K,KD1,-1
        DPHIB  = PHIL(L) - PHIH(L+1)
        DPHIT  = PHIH(L) - PHIL(L)
!
        DLB(L) = DPHIB * ETA(L) ! here eta contains 1/(L*(1+nu*q))
        DLT(L) = DPHIT * ETA(L)
!
        QRB(L) = DPHIB
        QRT(L) = DPHIT
!
        ETA(L) = ETA(L+1) + DPHIB

        HOL(L) = HOL(L) + ETA(L)
        hstold = hst(l)
        HST(L) = HST(L) + ETA(L)
!
        ETA(L) = ETA(L) + DPHIT
      ENDDO
!
!     For the cloud top layer
!
      L = KD

      DPHIB  = PHIL(L) - PHIH(L+1)
!
      DLB(L) = DPHIB * ETA(L)
!
      QRB(L) = DPHIB
      QRT(L) = DPHIB
!
      ETA(L) = ETA(L+1) + DPHIB

      HOL(L) = HOL(L) + ETA(L)
      HST(L) = HST(L) + ETA(L)
!
!     To determine KBL internally -- If KBL is defined externally
!     the following two loop should be skipped
!
      if (sgcs(kd) < 0.5_kp) then
         hcrit = hcritd
      elseif (sgcs(kd) > 0.65_kp) then
         hcrit = hcrits
      else
         hcrit = (hcrits*(sgcs(kd)-0.5_kp) + hcritd*(0.65_kp-sgcs(kd)))&
     &         * (one/0.15_kp)
      endif
      IF (CALKBL) THEN
         KTEM = MAX(KD+1, KBLMX)
         hmin = hol(k)
         kmin = k
         do l=km1,kd,-1
           if (hmin > hol(l)) then
             hmin = hol(l)
             kmin = l
           endif
         enddo
         if (kmin == k) return
         hmax = hol(k)
         kmax = k
         do l=km1,ktem,-1
           if (hmax < hol(l)) then
             hmax = hol(l)
             kmax = l
           endif
         enddo
         kmxb = kmax
         if (kmax < kmin) then
           kmax = k
           kmxb = k
           hmax = hol(kmax)
         elseif (kmax < k) then
           do l=kmax+1,k
!            if (abs(hol(kmax)-hol(l)) > half * hcrit) then
             if (abs(hol(kmax)-hol(l)) >        hcrit) then
               kmxb = l - 1
               exit
             endif
           enddo
         endif
         kmaxm1 = kmax - 1
         kmaxp1 = kmax + 1
         kblpmn = kmax
!
         dhdp(kmax:k) = zero
         dhdpmn = dhdp(kmax)
         do l=kmaxm1,ktem,-1
           dhdp(l) = (HOL(L)-HOL(L+1)) / (PRL(L+2)-PRL(L))
           if (dhdp(l) < dhdpmn) then
             dhdpmn = dhdp(l)
             kblpmn = l + 1
           elseif (dhdp(l) > zero .and. l <= kmin) then
             exit
           endif
         enddo
         kbl = kmax
         if (kblpmn < kmax) then
           do l=kblpmn,kmaxm1
             if (hmax-hol(l) < half*hcrit) then
               kbl = l
               exit
             endif
           enddo
         endif
!
         klcl = kd1
         if (kmax > kd1) then
           do l=kmaxm1,kd1,-1
             if (hmax > hst(l)) then
               klcl = l+1
               exit
             endif
           enddo
         endif

!        if (klcl == kd .or. klcl < ktem) return

!        This is to handle mid-level convection from quasi-uniform h

         if (kmax < kmxb) then
           kmax   = max(kd1, min(kmxb,k))
           kmaxm1 = kmax - 1
           kmaxp1 = kmax + 1
         endif


!        if (prl(Kmaxp1) - prl(klcl) > 250.0 ) return

         ii  = max(kbl,kd1)
         kbl = max(klcl,kd1)
         tem = min(50.0_kp,max(10.0_kp,(prl(kmaxp1)-prl(kd))*0.10_kp))
         if (prl(kmaxp1) - prl(ii) > tem .and. ii > kbl) kbl = ii

         if (kbl .ne. ii) then
           if (PRL(kmaxp1)-PRL(KBL) > bldmax) kbl = max(kbl,ii)
         endif
         if (kbl < ii) then
           if (hol(ii)-hol(ii-1) > half*hcrit) kbl = ii
         endif

         if (prl(kbl) - prl(klcl) > pcrit_lcl) return
!
!        KBL  = min(kmax, MAX(KBL,KBLMX))
         KBL  = min(kblmn, MAX(KBL,KBLMX))
!        kbl  = min(kblh,kbl)
!!!
!        tem1 = max(prl(kP1)-prl(k),                                    &
!    &                     min((prl(kbl) - prl(kd))*0.05, 10.0))
!!   &                     min((prl(kbl) - prl(kd))*0.05, 20.0))
!!   &                     min((prl(kbl) - prl(kd))*0.05, 30.0))
!        if (prl(kp1)-prl(kbl) < tem1) then
!          KTEM = MAX(KD+1, KBLMX)
!          do l=k,KTEM,-1
!            tem = prl(kp1) - prl(l)
!            if (tem > tem1) then
!              kbl = min(kbl,l)
!              exit
!            endif
!          enddo
!        endif
!        if (kbl == kblmx .and. kmax >= km1) kbl = k - 1
!!!

         KPBL = KBL

      ELSE
         KBL = KPBL
      ENDIF
!
      KBL = min(kmax,MAX(KBL,KD+2))
      KB1 = KBL - 1
!
      if (PRL(Kmaxp1)-PRL(KBL) > bldmax .or. kb1 <= kd ) then
!    &          .or. PRL(Kmaxp1)-PRL(KBL) < bldmin) then
        return
      endif
!
      PRIS     = ONE / (PRL(KP1)-PRL(KBL))
      PRISM    = ONE / (PRL(Kmaxp1)-PRL(KBL))
      TX1      = ETA(KBL)         ! geopotential height at KBL
!
      GMS(KBL) = zero
      XI(KBL)  = zero
      ZET(KBL) = zero
!
      shal_fac = one
!     if (prl(kbl)-prl(kd) < 300.0 .and. kmax == k) shal_fac = shalfac
      if (prl(kbl)-prl(kd) < 350.0_kp .and. kmax == k) shal_fac = shalfac
      DO L=Kmax,KD,-1
        IF (L >= KBL) THEN
          ETA(L) = (PRL(Kmaxp1)-PRL(L)) * PRISM
        ELSE
          ZET(L) = (ETA(L) - TX1) * ONEBG
          XI(L)  =  ZET(L) * ZET(L) * (QUDFAC*shal_fac)
          ETA(L) =  ZET(L) - ZET(L+1)
          GMS(L) =  XI(L)  - XI(L+1)
        ENDIF
      ENDDO
      if (kmax < k) then
        do l=kmaxp1,kp1
          eta(l) = zero
        enddo
      endif
!
      HBL = HOL(Kmax) * ETA(Kmax)
      QBL = QOL(Kmax) * ETA(Kmax)
      QLB = CLL(Kmax) * ETA(Kmax)
      QIB = CIL(Kmax) * ETA(Kmax)
      TX1 = QST(Kmax) * ETA(Kmax)
!
      DO L=Kmaxm1,KBL,-1
         TEM = ETA(L) - ETA(L+1)
         HBL = HBL + HOL(L) * TEM
         QBL = QBL + QOL(L) * TEM
         QLB = QLB + CLL(L) * TEM
         QIB = QIB + CIL(L) * TEM
         TX1 = TX1 + QST(L) * TEM
      ENDDO

!     if (ctei .and. sgcs(kd) > 0.65) then
!        hbl = hbl * hpert_fac
!        qbl = qbl * hpert_fac
!     endif

!                                   Find Min value of HOL in TX2
      TX2 = HOL(KD)
      IDH = KD1
      DO L=KD1,KB1
        IF (HOL(L) < TX2) THEN
           TX2 = HOL(L)
           IDH = L             ! Level of minimum moist static energy!
        ENDIF
      ENDDO
      IDH = 1
!     IDH = MAX(KD1, IDH)
      IDH = MAX(KD, IDH)       ! Moorthi May, 31, 2019
!
      TEM1 = HBL - HOL(KD)
      TEM  = HBL - HST(KD1) - LTL(KD1) * NU *(QOL(KD1)-QST(KD1))
      LOWEST = KD == KB1

      lcon = kd
      do l=kb1,kd1,-1
        if (hbl >= hst(l)) then
          lcon = l
          exit
        endif
      enddo
!
      if (lcon == kd .or. kbl <= kd .or. prl(kbl)-prsm(lcon) > 150.0_kp) &
     &                                    return
!
      TX1    = RHFACS - QBL / TX1       !     Average RH

      cnvflg = (TEM > ZERO .OR. (LOWEST .AND. TEM1 >= ZERO))            &
     &         .AND. TX1 < RHRAM

      IF (.NOT. cnvflg) RETURN
!
      RHC = MAX(ZERO, MIN(ONE, EXP(-20.0_kp*TX1) ))
!
      wcbase = 0.1_kp
      if (ntrc > 0) then
        DO N=1,NTRC
          RBL(N) = ROI(Kmax,N) * ETA(Kmax)
        ENDDO
        DO N=1,NTRC
          DO L=KmaxM1,KBL,-1
            RBL(N) = RBL(N) + ROI(L,N)*(ETA(L)-ETA(L+1))
          ENDDO
        ENDDO
!
!       if (ntk > 0 .and. aw_scal) then
        if (ntk > 0) then
          if (rbl(ntk) > zero) then
            wcbase = min(two, max(wcbase, sqrt(twoo3*rbl(ntk))))
!           wcbase = min(one, max(wcbase, sqrt(twoo3*rbl(ntk))))
          endif
        endif

      endif
!
      TX4 = zero
      TX5 = zero
!
      TX3 = QST(KBL) - GAF(KBL) * HST(KBL)
      DO L=KBL,K
        QIL(L) = MAX(ZERO, MIN(ONE, (TCR-TCL-TOL(L))*TCRF))
      ENDDO
!
      DO L=KB1,KD1,-1
        lp1      = l + 1
        TEM      = QST(L) - GAF(L) * HST(L)
        TEM1     = (TX3 + TEM) * half
        ST2      = (GAF(L)+GAF(LP1)) * half
!
        FCO(LP1) =            TEM1 + ST2 * HBL

        RNN(LP1) = ZET(LP1) * TEM1 + ST2 * TX4
        GMH(LP1) = XI(LP1)  * TEM1 + ST2 * TX5
!
        TX3      = TEM
        TX4      = TX4 + ETA(L) * HOL(L)
        TX5      = TX5 + GMS(L) * HOL(L)
!
        QIL(L)   = MAX(ZERO, MIN(ONE, (TCR-TCL-TOL(L))*TCRF))
        QLL(LP1) = (half*ALHF) * ST2 * (QIL(L)+QIL(LP1)) + ONE
      ENDDO
!
!     FOR THE CLOUD TOP -- L=KD
!
      L = KD
!
      lp1      = l + 1
      TEM      = QST(L) - GAF(L) * HST(L)
      TEM1     = (TX3 + TEM) * half
      ST2      = (GAF(L)+GAF(LP1)) * half
!
      FCO(LP1) =            TEM1 + ST2 * HBL
      RNN(LP1) = ZET(LP1) * TEM1 + ST2 * TX4
      GMH(LP1) = XI(LP1)  * TEM1 + ST2 * TX5
!
      FCO(L)   = TEM + GAF(L) * HBL
      RNN(L)   = TEM * ZET(L) + (TX4 + ETA(L)*HOL(L)) * GAF(L)
      GMH(L)   = TEM * XI(L)  + (TX5 + GMS(L)*HOL(L)) * GAF(L)
!
!   Replace FCO for the Bottom
!
      FCO(KBL) = QBL
      RNN(KBL) = zero
      GMH(KBL) = zero
!
      QIL(KD)  =  MAX(ZERO, MIN(ONE, (TCR-TCL-TOL(KD))*TCRF))
      QLL(KD1) = (half*ALHF) * ST2 * (QIL(KD) + QIL(KD1)) + ONE
      QLL(KD ) = ALHF * GAF(KD) * QIL(KD) + ONE
!
      st1  = qil(kd)
      st2  = c0i * st1
      if (c0ifac > 1.0e-6_kp) st2  = st2 * exp(c0ifac*min(tol(kd)-t0c,zero))
      tem  = c0  * (one-st1)
      tem2 = st2*qi0 + tem*qw0
!
      DO L=KD,KB1
         lp1    = l + 1
         tx2    = akt(l) * eta(l)
         tx1    = tx2 * tem2
         q0u(l) = tx1
         FCO(L) = FCO(LP1) - FCO(L) + tx1
         RNN(L) = RNN(LP1) - RNN(L)                                     &
     &          + ETA(L)*(QOL(L)+CLL(L)+CIL(L)) + tx1*zet(l)
         GMH(L) = GMH(LP1) - GMH(L)                                     &
     &          + GMS(L)*(QOL(L)+CLL(L)+CIL(L)) + tx1*xi(l)
!
         tem1   = (one-akt(l)) * eta(l)

         AKT(L) = QLL(L)   + (st2 + tem) * tx2

         AKC(L) = one / AKT(L)
!
         st1    = half * (qil(l)+qil(lp1))
         st2    = c0i * st1
         if (c0ifac > 1.0e-6_kp) st2 = st2 * exp(c0ifac*min(tol(lp1)-t0c,zero))
         tem    = c0  * (one-st1)
         tem2   = st2*qi0 + tem*qw0
!
         BKC(L) = QLL(LP1) - (st2 + tem) * tem1
!
         tx1    = tem1*tem2
         q0d(l) = tx1
         FCO(L) = FCO(L) + tx1
         RNN(L) = RNN(L) + tx1*zet(lp1)
         GMH(L) = GMH(L) + tx1*xi(lp1)
      ENDDO

      qw00 = qw0
      qi00 = qi0
      ii = 0
  777 continue
!
      ep_wfn = .false.
      RNN(KBL) = zero
      TX3      = bkc(kb1) * (QIB + QLB)
      TX4      = zero
      TX5      = zero
      DO L=KB1,KD1,-1
        TEM    = BKC(L-1)       * AKC(L)
        TX3    = (TX3 + FCO(L)) * TEM
        TX4    = (TX4 + RNN(L)) * TEM
        TX5    = (TX5 + GMH(L)) * TEM
      ENDDO
      IF (KD < KB1) THEN
         HSD   = HST(KD1) + LTL(KD1) *  NU *(QOL(KD1)-QST(KD1))
      ELSE
         HSD   = HBL
      ENDIF
!
      TX3 = (TX3 + FCO(KD)) * AKC(KD)
      TX4 = (TX4 + RNN(KD)) * AKC(KD)
      TX5 = (TX5 + GMH(KD)) * AKC(KD)
      ALM = ALHF*QIL(KD) - LTL(KD) * VTF(KD)
!
      HSU = HST(KD) + LTL(KD) * NU * (QOL(KD)-QST(KD))
!
!===> VERTICAL INTEGRALS NEEDED TO COMPUTE THE ENTRAINMENT PARAMETER
!
      TX1 = ALM * TX4
      TX2 = ALM * TX5

      DO L=KD,KB1
        TAU = HOL(L) - HSU
        TX1 = TX1 + TAU * ETA(L)
        TX2 = TX2 + TAU * GMS(L)
      ENDDO
!
!     MODIFY HSU TO INCLUDE CLOUD LIQUID WATER AND ICE TERMS
!
      HSU    = HSU - ALM * TX3
!
      CLP    = ZERO
      ALM    = -100.0_kp
      HOS    = HOL(KD)
      QOS    = QOL(KD)
      QIS    = CIL(KD)
      QLS    = CLL(KD)

      cnvflg = HBL > HSU .and. abs(tx1) > 1.0e-4_kp

!***********************************************************************

      ST1 = HALF*(HSU + HSD)

      IF (cnvflg) THEN
!
!  STANDARD CASE:
!   CLOUD CAN BE NEUTRALLY BOUYANT AT MIDDLE OF LEVEL KD W/ +VE LAMBDA.
!   EPP < .25 IS REQUIRED TO HAVE REAL ROOTS.
!
        clp = one
        st2 = hbl - hsu
!
        if (tx2 == zero) then
          alm = - st2 / tx1
          if (alm > almax) alm = -100.0_kp
        else
          x00 = tx2 + tx2
          epp = tx1 * tx1 - (x00+x00)*st2
          if (epp > zero) then
            x00  = one / x00
            tem  = sqrt(epp)
            tem1 = (-tx1-tem)*x00
            tem2 = (-tx1+tem)*x00
            if (tem1 > almax) tem1 = -100.0_kp
            if (tem2 > almax) tem2 = -100.0_kp
            alm  = max(tem1,tem2)

          endif
        endif
!
!  CLIP CASE:
!   NON-ENTRAINIG CLOUD DETRAINS IN LOWER HALF OF TOP LAYER.
!   NO CLOUDS ARE ALLOWED TO DETRAIN BELOW THE TOP LAYER.
!
      ELSEIF (HBL <= HSU .AND. HBL > ST1) THEN
        ALM = ZERO
!       CLP = (HBL-ST1) / (HSU-ST1)    ! commented on Jan 16, 2010
      ENDIF
!
      cnvflg = .TRUE.
      IF (ALMIN1 > zero) THEN
        IF (ALM >= ALMIN1) cnvflg = .FALSE.
      ELSE
        LOWEST = KD == KB1
        IF ( (ALM > ZERO) .OR.                                          &
     &     (.NOT. LOWEST .AND. ALM == ZERO) ) cnvflg = .FALSE.
      ENDIF
!
!===>  IF NO SOUNDING MEETS SECOND CONDITION, RETURN
!
      IF (cnvflg) THEN
         IF (ii > 0 .or. (qw00 == zero .and. qi00 == zero)) RETURN
         CLP = one
         ep_wfn = .true.
         GO TO 888
      ENDIF
!
      st1s = ONE
      IF(CLP > ZERO .AND. CLP < ONE) THEN
        ST1     = HALF*(ONE+CLP)
        ST2     = ONE - ST1
        st1s    = st1
        hstkd   = hst(kd)
        qstkd   = qst(kd)
        ltlkd   = ltl(kd)
        q0ukd   = q0u(kd)
        q0dkd   = q0d(kd)
        dlbkd   = dlb(kd)
        qrbkd   = qrb(kd)
!
        HST(KD) = HST(KD)*ST1 + HST(KD1)*ST2
        HOS     = HOL(KD)*ST1 + HOL(KD1)*ST2
        QST(KD) = QST(KD)*ST1 + QST(KD1)*ST2
        QOS     = QOL(KD)*ST1 + QOL(KD1)*ST2
        QLS     = CLL(KD)*ST1 + CLL(KD1)*ST2
        QIS     = CIL(KD)*ST1 + CIL(KD1)*ST2
        LTL(KD) = LTL(KD)*ST1 + LTL(KD1)*ST2
!
        DLB(KD) = DLB(KD)*CLP
        qrb(KD) = qrb(KD)*CLP
        ETA(KD) = ETA(KD)*CLP
        GMS(KD) = GMS(KD)*CLP
        Q0U(KD) = Q0U(KD)*CLP
        Q0D(KD) = Q0D(KD)*CLP
      ENDIF
!
!
!***********************************************************************
!
!    Critical workfunction is included in this version
!
      ACR  = zero
      TEM  = PRL(KD1) - (PRL(KD1)-PRL(KD)) * CLP * HALF
      tx1  = PRL(KBL) - TEM
      tx2  = min(900.0_kp, max(tx1,100.0_kp))
      tem1 = log(tx2*0.01_kp) * oneolog10
      tem2 = one - tem1
      if ( kdt == 1 ) then
!       rel_fac = (dt * facdt)  / (tem1*12.0_kp + tem2*3.0)
        rel_fac = (dt * facdt)  / (tem1*6.0_kp + tem2*adjts_s)
      else
        rel_fac = (dt * facdt) / (tem1*adjts_d + tem2*adjts_s)
      endif
!
!     rel_fac = max(zero, min(one,rel_fac))
      rel_fac = max(zero, min(half,rel_fac))
      
      IF (CRTFUN) THEN
        iwk = tem*0.02_kp - 0.999999999_kp
        iwk = MAX(1, MIN(iwk, 16))
        ACR = tx1 * (AC(iwk) + tem * AD(iwk)) * CCWF
      ENDIF
!
!===>  NORMALIZED MASSFLUX
!
!  ETA IS THE THICKNESS COMING IN AND normalized  MASS FLUX GOING OUT.
!  GMS IS THE THICKNESS SQUARE ; IT IS LATER REUSED FOR GAMMA_S
!
!     ETA(K) = ONE

      DO L=KB1,KD,-1
        ETA(L)  = ETA(L+1) + ALM * (ETA(L) + ALM * GMS(L))
        ETAI(L) = one / ETA(L)
      ENDDO
      ETAI(KBL) = one
!
!===>  CLOUD WORKFUNCTION
!
      WFN    = ZERO
      AKM    = ZERO
      DET    = ZERO
      HCC    = HBL
      cnvflg = .FALSE.
      QTL    = QST(KB1) - GAF(KB1)*HST(KB1)
      TX1    = HBL
!
      qtv    = qbl
      det    = qlb + qib
!
      tx2    = zero
      dpneg  = zero
!
      DO L=KB1,KD1,-1
         lm1 = l - 1
         lp1 = l + 1
         DEL_ETA = ETA(L) - ETA(LP1)
         HCCP = HCC + DEL_ETA*HOL(L)
!
         QTLP = QST(LM1) - GAF(LM1)*HST(LM1)
         QTVP = half * ((QTLP+QTL)*ETA(L)                               &
     &               + (GAF(L)+GAF(LM1))*HCCP)
         ST1  = ETA(L)*Q0U(L) + ETA(LP1)*Q0D(L)
         DETP = (BKC(L)*DET - (QTVP-QTV)                                &
     &        + DEL_ETA*(QOL(L)+CLL(L)+CIL(L)) + ST1)  * AKC(L)

         TEM1   = AKT(L)   - QLL(L)
         TEM2   = QLL(LP1) - BKC(L)
         RNS(L) = TEM1*DETP  + TEM2*DET - ST1

         qtp    = half * (qil(L)+qil(LM1))
         tem2   = min(qtp*(detp-eta(l)*qw00),                           &
     &               (one-qtp)*(detp-eta(l)*qi00))
         st1    = min(tx2,tem2)
         tx2    = tem2
!
         IF (rns(l) < zero .or. st1 < zero) ep_wfn = .TRUE.
         IF (DETP <= ZERO) cnvflg = .TRUE.

         ST1  = HST(L) - LTL(L)*NU*(QST(L)-QOL(L))


         TEM2 = HCCP   + DETP   * QTP * ALHF
!
         ST2  = LTL(L) * VTF(L)
         TEM5 = CLL(L) + CIL(L)
         TEM3 = (TX1  - ETA(LP1)*ST1 - ST2*(DET-TEM5*eta(lp1))) * DLB(L)
         TEM4 = (TEM2 - ETA(L  )*ST1 - ST2*(DETP-TEM5*eta(l)))  * DLT(L)
!
         ST1  = TEM3 + TEM4

         WFN = WFN + ST1       
         AKM = AKM - min(ST1,ZERO)

         if (st1 < zero .and. wfn < zero) then
           dpneg = dpneg + prl(lp1) - prl(l)
         endif

         BUY(L) = half * (tem3/(eta(lp1)*qrb(l)) + tem4/(eta(l)*qrt(l)))
!
         HCC = HCCP
         DET = DETP
         QTL = QTLP
         QTV = QTVP
         TX1 = TEM2

      ENDDO

      DEL_ETA = ETA(KD) - ETA(KD1)
      HCCP    = HCC + DEL_ETA*HOS
!
      QTLP    = QST(KD) - GAF(KD)*HST(KD)
      QTVP    = QTLP*ETA(KD) + GAF(KD)*HCCP
      ST1     = ETA(KD)*Q0U(KD) + ETA(KD1)*Q0D(KD)
      DETP    = (BKC(KD)*DET - (QTVP-QTV)                               &
     &        + DEL_ETA*(QOS+QLS+QIS) + ST1) * AKC(KD)
!
      TEM1    = AKT(KD)  - QLL(KD)
      TEM2    = QLL(KD1) - BKC(KD)
      RNS(KD) = TEM1*DETP  + TEM2*DET - ST1
!
      IF (rns(kd) < zero) ep_wfn = .TRUE.
      IF (DETP   <= ZERO) cnvflg = .TRUE.
!
  888 continue

      if (ep_wfn) then
        IF ((qw00 == zero .and. qi00 == zero)) RETURN
        if (ii == 0) then
          ii  = 1
          if (clp > zero .and. clp < one) then
            hst(kd) = hstkd
            qst(kd) = qstkd
            ltl(kd) = ltlkd
            q0u(kd) = q0ukd
            q0d(kd) = q0dkd
            dlb(kd) = dlbkd
            qrb(kd) = qrbkd
          endif
          do l=kd,kb1
            lp1    = l + 1
            FCO(L) = FCO(L) - q0u(l) - q0d(l)
            RNN(L) = RNN(L) - q0u(l)*zet(l) - q0d(l)*zet(lp1)
            GMH(L) = GMH(L) - q0u(l)*xi(l)  - q0d(l)*zet(lp1)
            ETA(L) = ZET(L) - ZET(LP1)
            GMS(L) = XI(L)  - XI(LP1)
            Q0U(L) = zero
            Q0D(L) = zero
          ENDDO
          qw00 = zero
          qi00 = zero

          go to 777
        else
          cnvflg = .true.
        endif
      endif
!
!
!     ST1 = 0.5 * (HST(KD)  - LTL(KD)*NU*(QST(KD)-QOS)                  &
!    &          +  HST(KD1) - LTL(KD1)*NU*(QST(KD1)-QOL(KD1)))
!
      ST1  = HST(KD) - LTL(KD)*NU*(QST(KD)-QOS)
      ST2  = LTL(KD) * VTF(KD)
      TEM5 = (QLS + QIS) * eta(kd1)
      ST1  = HALF * (TX1-ETA(KD1)*ST1-ST2*(DET-TEM5))*DLB(KD)
!
      WFN  = WFN + ST1
      AKM  = AKM - min(ST1,ZERO)   ! Commented on 08/26/02 - does not include top
!

      BUY(KD) = ST1 / (ETA(KD1)*qrb(kd))
!
      DET = DETP
      HCC = HCCP
      AKM = AKM / WFN


!***********************************************************************
!
      IF (WRKFUN) THEN ! If only to calculate workfunction save it and return
        IF (WFN >= zero) WFNC = WFN
        RETURN
      ELSEIF (.NOT. CRTFUN) THEN
        ACR = WFNC
      ENDIF
!
!===>  THIRD CHECK BASED ON CLOUD WORKFUNCTION
!
      CALCUP = .FALSE.

      TEM = max(0.05_kp, MIN(CD*200.0_kp, MAX_NEG_BOUY))
      IF (.not. cnvflg    .and. WFN > ACR .and.                         &
     &    dpneg < dpnegcr .and. AKM <= TEM)     CALCUP = .TRUE.
!
!===>  IF NO SOUNDING MEETS THIRD CONDITION, RETURN
!
      IF (.NOT. CALCUP) RETURN
!
! This is for not LL - 20050601
!     IF (ALMIN2 .NE. zero) THEN
!       IF (ALMIN1 .NE. ALMIN2) ST1 = one / max(ONE_M10,(ALMIN2-ALMIN1))
!       IF (ALM < ALMIN2) THEN
!          CLP = CLP * max(zero, min(one,(0.3 + 0.7*(ALM-ALMIN1)*ST1)))
!!         CLP = CLP * max(0.0, min(1.0,(0.2 + 0.8*(ALM-ALMIN1)*ST1)))
!!         CLP = CLP * max(0.0, min(1.0,(0.1 + 0.9*(ALM-ALMIN1)*ST1)))
!       ENDIF
!     ENDIF
!
      CLP = CLP * RHC
      dlq = zero
      tem = one / (one + dlq_fac)
      do l=kd,kb1
        rnn(l) = rns(l) * tem
        dlq(l) = rns(l) * tem * dlq_fac
      enddo
      DO L=KBL,K 
        RNN(L) = zero
      ENDDO
!
!     If downdraft is to be invoked, do preliminary check to see
!     if enough rain is available and then call DDRFT.
!
      DDFT = .FALSE.
      IF (dpd > zero) THEN
        TRAIN = zero
        IF (CLP > zero) THEN
          DO L=KD,KB1
            TRAIN = TRAIN + RNN(L)
          ENDDO
        ENDIF

        PL = (PRL(KD1) + PRL(KD))*HALF
        IF (TRAIN > 1.0e-4_kp .AND. PL <= dpd*prl(kp1)) DDFT  = .TRUE.
      ENDIF
!
      IF (DDFT) THEN ! Downdraft scheme based on (Cheng and Arakawa, 1997)
        CALL DDRFT(                                                     &
     &              K,   KP1,    KD                                     &
     &,             TLA, ALFIND, wcbase                                 &
     &,             TOL, QOL, HOL,   PRL, QST, HST, GAM, GAF            &
!    &,             TOL, QOL, HOL,   PRL, QST, HST, GAM, GAF, HBL, QBL  &
     &,             QRB, QRT, BUY,   KBL, IDH, ETA, RNN, ETAI           &
     &,             ALM, WFN, TRAIN, DDFT                               &
     &,             ETD, HOD, QOD,   EVP, DOF, CLDFR, ETZ               &
     &,             GMS, GSD, GHD,   wvl)               

      ENDIF
!
!  No Downdraft case (including case with no downdraft solution)
!  ---------------------------------------------------------
!
      IF (.NOT. DDFT) THEN
        DO L=KD,KP1
          ETD(L) = zero
          HOD(L) = zero
          QOD(L) = zero
          wvl(l) = zero
        ENDDO
        DO L=KD,K
          EVP(L) = zero
          ETZ(L) = zero
        ENDDO

      ENDIF

!
!===> CALCULATE GAMMAS  i.e. TENDENCIES PER UNIT CLOUD BASE MASSFLUX
!           Includes downdraft terms!

      avh = zero

!
!     Fraction of detrained condensate evaporated
!
!     tem1 = max(ZERO, min(HALF, (prl(kd)-FOUR_P2)*ONE_M2))
!     tem1 = max(ZERO, min(HALF, (prl(kd)-300.0)*0.005))
      tem1 = zero
!     tem1 = 1.0
!     if (kd1 == kbl) tem1 = 0.0
!
      tem2 = one - tem1
      TEM  = DET * QIL(KD)


      st1 = (HCC+ALHF*TEM-ETA(KD)*HST(KD)) / (one+gam(KD))
      DS  = ETA(KD1) * (HOS- HOL(KD)) - ALHL*(QOS - QOL(KD))
      DH  = ETA(KD1) * (HOS- HOL(KD))


      GMS(KD) = (DS + st1 - tem1*det*alhl-tem*alhf) * PRI(KD)
      GMH(KD) = PRI(KD) * (HCC-ETA(KD)*HOS + DH)
!
!      TENDENCY FOR SUSPENDED ENVIRONMENTAL ICE AND/OR LIQUID WATER
!
      QLL(KD) = (tem2*(DET-TEM) + ETA(KD1)*(QLS-CLL(KD))                &
     &        + (one-QIL(KD))*dlq(kd) - ETA(KD)*QLS ) * PRI(KD)

      QIL(KD) =     (tem2*TEM + ETA(KD1)*(QIS-CIL(KD))                  &
     &        + QIL(KD)*dlq(kd) - ETA(KD)*QIS ) * PRI(KD)
!
      GHD(KD) = zero
      GSD(KD) = zero
!
      DO L=KD1,K
         lm1 = l - 1
         ST1 = ONE - ALFINT(L,1)
         ST2 = ONE - ALFINT(L,2)
         ST3 = ONE - ALFINT(L,3)
         ST4 = ONE - ALFINT(L,4)
         ST5 = ONE - ALFIND(L)
         HB       = ALFINT(L,1)*HOL(LM1) + ST1*HOL(L)
         QB       = ALFINT(L,2)*QOL(LM1) + ST2*QOL(L)

         TEM      = ALFINT(L,4)*CIL(LM1) + ST4*CIL(L)
         TEM2     = ALFINT(L,3)*CLL(LM1) + ST3*CLL(L)
 
         TEM1     = ETA(L) * (TEM - CIL(L))
         TEM3     = ETA(L) * (TEM2 - CLL(L))

         HBD      = ALFIND(L)*HOL(LM1) + ST5*HOL(L)
         QBD      = ALFIND(L)*QOL(LM1) + ST5*QOL(L)

         TEM5     = ETD(L) * (HOD(L) - HBD)
         TEM6     = ETD(L) * (QOD(L) - QBD)
!
         DH       = ETA(L) * (HB - HOL(L)) + TEM5
         DS       = DH - ALHL * (ETA(L) * (QB - QOL(L)) + TEM6)

         GMH(L)   = DH * PRI(L)
         GMS(L)   = DS * PRI(L)
!
         GHD(L)   = TEM5 * PRI(L)
         GSD(L)   = (TEM5 - ALHL * TEM6) * PRI(L)
!
         QLL(L)   = (TEM3 + (one-QIL(L))*dlq(l)) * PRI(L)
         QIL(L)   = (TEM1 + QIL(L)*dlq(l)) * PRI(L)

         TEM1     = ETA(L) * (CIL(LM1) - TEM)
         TEM3     = ETA(L) * (CLL(LM1) - TEM2)

         DH       = ETA(L) * (HOL(LM1) - HB) - TEM5
         DS       = DH - ALHL * ETA(L) * (QOL(LM1) - QB)                &
     &                 + ALHL * (TEM6 - EVP(LM1))

         GMH(LM1) = GMH(LM1) + DH * PRI(LM1)
         GMS(LM1) = GMS(LM1) + DS * PRI(LM1)
!
         GHD(LM1) = GHD(LM1) - TEM5 * PRI(LM1)
         GSD(LM1) = GSD(LM1) - (TEM5-ALHL*(TEM6-EVP(LM1))) * PRI(LM1)

         QIL(LM1) = QIL(LM1) + TEM1 * PRI(LM1)
         QLL(LM1) = QLL(LM1) + TEM3 * PRI(LM1)
!
        avh = avh + gmh(lm1)*(prs(l)-prs(lm1))

      ENDDO
!
      HBD    = HOL(K)
      QBD    = QOL(K)
      TEM5   =  ETD(KP1) * (HOD(KP1) - HBD)
      TEM6   =  ETD(KP1) * (QOD(KP1) - QBD)
      DH     = - TEM5
      DS     = DH + ALHL * TEM6
      TEM1   = DH * PRI(K)
      TEM2   = (DS - ALHL * EVP(K)) * PRI(K)
      GMH(K) = GMH(K) + TEM1
      GMS(K) = GMS(K) + TEM2
      GHD(K) = GHD(K) + TEM1
      GSD(K) = GSD(K) + TEM2

!
      avh    = avh + gmh(K)*(prs(KP1)-prs(K))
!
      tem4   = - GRAVFAC * pris
      TX1    = DH * tem4
      TX2    = DS * tem4
!
      DO L=KBL,K
        GMH(L) = GMH(L) + TX1
        GMS(L) = GMS(L) + TX2
        GHD(L) = GHD(L) + TX1
        GSD(L) = GSD(L) + TX2
!
        avh = avh + tx1*(prs(l+1)-prs(l))
      ENDDO
!
!***********************************************************************
!***********************************************************************

!===>  KERNEL (AKM) CALCULATION BEGINS

!===>  MODIFY SOUNDING WITH UNIT MASS FLUX
!
      DO L=KD,K

         TEM1   = GMH(L)
         TEM2   = GMS(L)
         HOL(L) = HOL(L) +  TEM1*TESTMB
         QOL(L) = QOL(L) + (TEM1-TEM2)  * TESTMBOALHL
         HST(L) = HST(L) +  TEM2*(ONE+GAM(L))*TESTMB
         QST(L) = QST(L) +  TEM2*GAM(L) * TESTMBOALHL
         CLL(L) = CLL(L) + QLL(L) * TESTMB
         CIL(L) = CIL(L) + QIL(L) * TESTMB
      ENDDO
!
      if (alm > zero) then
        HOS     = HOS + GMH(KD)  * TESTMB
        QOS     = QOS + (GMH(KD)-GMS(KD)) * TESTMBOALHL
        QLS     = QLS + QLL(KD) * TESTMB
        QIS     = QIS + QIL(KD) * TESTMB
      else
        st2 = one - st1s
        HOS     = HOS + (st1s*GMH(KD)+st2*GMH(KD1))  * TESTMB
        QOS     = QOS + (st1s * (GMH(KD)-GMS(KD))                       &
     &                +  st2  * (GMH(KD1)-GMS(KD1))) * TESTMBOALHL
        HST(kd) = HST(kd) + (st1s*GMS(kd)*(ONE+GAM(kd))                 &
     &                    +  st2*gms(kd1)*(ONE+GAM(kd1))) * TESTMB
        QST(kd) = QST(kd) + (st1s*GMS(kd)*GAM(kd)                       &
     &                    +  st2*gms(kd1)*gam(kd1)) * TESTMBOALHL

        QLS     = QLS + (st1s*QLL(KD)+st2*QLL(KD1)) * TESTMB
        QIS     = QIS + (st1s*QIL(KD)+st2*QIL(KD1)) * TESTMB
      endif

!
      TEM = PRL(Kmaxp1) - PRL(Kmax)
      HBL = HOL(Kmax) * TEM
      QBL = QOL(Kmax) * TEM
      QLB = CLL(Kmax) * TEM
      QIB = CIL(Kmax) * TEM
      DO L=KmaxM1,KBL,-1
        TEM = PRL(L+1) - PRL(L)
        HBL = HBL + HOL(L) * TEM
        QBL = QBL + QOL(L) * TEM
        QLB = QLB + CLL(L) * TEM
        QIB = QIB + CIL(L) * TEM
      ENDDO
      HBL = HBL * PRISM
      QBL = QBL * PRISM
      QLB = QLB * PRISM
      QIB = QIB * PRISM

!     if (ctei .and. sgcs(kd) > 0.65) then
!        hbl = hbl * hpert_fac
!        qbl = qbl * hpert_fac
!     endif
 
!***********************************************************************

!===>  CLOUD WORKFUNCTION FOR MODIFIED SOUNDING, THEN KERNEL (AKM)
!
      AKM = ZERO
      TX1 = ZERO
      QTL = QST(KB1) - GAF(KB1)*HST(KB1)
      QTV = QBL
      HCC = HBL
      TX2 = HCC
      TX4 = (ALHF*half)*MAX(ZERO,MIN(ONE,(TCR-TCL-TOL(KB1))*TCRF))
!
      qtv = qbl
      tx1 = qib + qlb
!

      DO L=KB1,KD1,-1
         lm1 = l - 1
         lp1 = l + 1
         DEL_ETA = ETA(L) - ETA(LP1)
         HCCP = HCC + DEL_ETA*HOL(L)
!
         QTLP = QST(LM1) - GAF(LM1)*HST(LM1)
         QTVP = half * ((QTLP+QTL)*ETA(L) + (GAF(L)+GAF(LM1))*HCCP)

         DETP = (BKC(L)*TX1 - (QTVP-QTV)                                &
     &        +  DEL_ETA*(QOL(L)+CLL(L)+CIL(L))                         &
     &        +  ETA(L)*Q0U(L) + ETA(LP1)*Q0D(L)) * AKC(L)
         IF (DETP <= ZERO) cnvflg = .TRUE.

         ST1  = HST(L) - LTL(L)*NU*(QST(L)-QOL(L))

         TEM2 = (ALHF*half)*MAX(ZERO,MIN(ONE,(TCR-TCL-TOL(LM1))*TCRF))
         TEM1 = HCCP + DETP * (TEM2+TX4)

         ST2  = LTL(L) * VTF(L)
         TEM5 = CLL(L) + CIL(L)
         AKM  = AKM +                                                   &
     &     (  (TX2  -ETA(LP1)*ST1-ST2*(TX1-TEM5*eta(lp1))) * DLB(L)     &
     &      + (TEM1 -ETA(L  )*ST1-ST2*(DETP-TEM5*eta(l)))  * DLT(L) )
!
         HCC  = HCCP
         TX1  = DETP
         TX2  = TEM1
         QTL  = QTLP
         QTV  = QTVP
         TX4  = TEM2
      ENDDO
!
      if (cnvflg) return
!
!  Eventhough we ignore the change in lambda, we still assume
!  that the cLoud-top contribution is zero; as though we still
!  had non-bouyancy there.
!
!
      ST1  = HST(KD)  - LTL(KD)*NU*(QST(KD)-QOS)
      ST2  = LTL(KD)  * VTF(KD)
      TEM5 = (QLS + QIS) * eta(kd1)
      AKM  = AKM + HALF * (TX2-ETA(KD1)*ST1-ST2*(TX1-TEM5)) * DLB(KD)
!
      AKM  = (AKM - WFN) * TESTMBI


!***********************************************************************

!===>   MASS FLUX
!
      AMB  = - (WFN-ACR) / AKM
!
!===>   RELAXATION AND CLIPPING FACTORS
!
      AMB  = AMB * CLP * rel_fac

!!!   if (DDFT) AMB = MIN(AMB, ONE/CLDFRD)
       
!===>   SUB-CLOUD LAYER DEPTH LIMIT ON MASS FLUX

      AMBMAX = (PRL(KMAXP1)-PRL(KBL))*(FRACBL*GRAVCON)
      AMB    = MAX(MIN(AMB, AMBMAX),ZERO)

!***********************************************************************
!*************************RESULTS***************************************
!***********************************************************************

!===>  PRECIPITATION AND CLW DETRAINMENT
!
      if (amb > zero) then

!
!       if (wvl(kd) > zero) then
!         tx1 = one - amb * eta(kd) / (rho(kd)*wvl(kd))
!         sigf(kd) =  max(zero, min(one, tx1 * tx1))
!       endif
        if (aw_scal) then
          tx1 = (0.2_kp / max(alm, 1.0e-5_kp))
          tx2 = one - min(one, pi * tx1 * tx1 / area)

          tx2 = tx2 * tx2 

! comnet out the following for now - 07/23/18
!         do l=kd1,kbl
!           lp1 = min(K, l+1)
!           if (wvl(l) > zero .and. wvl(lp1) > zero) then
!             tx1     = one - amb *  (eta(l)+eta(lp1))
!    &                      / ((wvl(l)+wvl(lp1))*rho(l)*grav)
!             sigf(l) = max(zero, min(one, tx1 * tx1))
!           else
!             sigf(l) = min(one,tx2)
!           endif
!           sigf(l) = max(sigf(l), tx2)
!         enddo
!         sigf(kd) = sigf(kd1)
!         if (kbl < k) then
!           sigf(kbl+1:k) = sigf(kbl)
!         endif
          sigf(kd:k) = tx2
        else
          sigf(kd:k) = one
        endif

        tx1 = max(1.0e-6_kp, abs(gms(kd) * onebcp * sigf(kd)))
        amb = min(tx1*amb, tfrac_max*toi(kd)) / tx1

!
        avt = zero
        avq = zero
        avr = dof * sigf(kbl)
!
        DSFC = DSFC + AMB * ETD(K) * (one/DT) * sigf(kbl)
!
        DO L=K,KD,-1
          PCU(L) = PCU(L) + AMB*RNN(L)*sigf(l)      !  (A40)
          avr    = avr + rnn(l) * sigf(l)
        ENDDO
        pcu(k) = pcu(k) + amb * dof * sigf(kbl)
!
!===> TEMPARATURE AND Q CHANGE AND CLOUD MASS FLUX DUE TO CLOUD TYPE KD
!
        TX1 = AMB * ONEBCP
        TX2 = AMB * ONEOALHL
        DO L=KD,K
          delp   = prs(l+1) - prs(l)
          tx3    = amb * sigf(l)
          ST1    = GMS(L) * TX1 * sigf(l)
          TOI(L) = TOI(L) + ST1
          TCU(L) = TCU(L) + ST1
          TCD(L) = TCD(L) + GSD(L) * TX1 * sigf(l)
!
          st1 = st1 - ELOCP * (QIL(L) + QLL(L)) * tx3

          avt = avt + st1 * delp

          FLX(L)  = FLX(L)  + ETA(L) * tx3
          FLXD(L) = FLXD(L) + ETD(L) * tx3
!
          QII(L)  = QII(L) + QIL(L) * tx3
          TEM     = zero

          QLI(L)  = QLI(L) + QLL(L) * tx3 + TEM

          ST1     = (GMH(L)-GMS(L)) * TX2 * sigf(l)

          QOI(L)  = QOI(L) + ST1
          QCU(L)  = QCU(L) + ST1
          QCD(L)  = QCD(L) + (GHD(L)-GSD(L)) * TX2 * sigf(l)
!
          avq = avq + (st1 + (QLL(L)+QIL(L))*tx3) * delp
!         avq = avq + st1 * (prs(l+1)-prs(l))
!         avr = avr + (QLL(L) + QIL(L)*(1+alhf/alhl))
          avr = avr + (QLL(L) + QIL(L)) * delp * sigf(l) * gravcon

!      Correction for negative condensate!
          if (qii(l) < zero) then
            tem    = qii(l) * elfocp
            QOI(L) = QOI(L) + qii(l)
            qcu(l) = qcu(l) + qii(l)
            toi(l) = toi(l) - tem
            tcu(l) = tcu(l) - tem
            qii(l) = zero
          endif
          if (qli(l) < zero) then
            tem    = qli(l) * elocp
            QOI(L) = QOI(L) + qli(l)
            qcu(l) = qcu(l) + qli(l)
            toi(l) = toi(l) - tem
            tcu(l) = tcu(l) - tem
            qli(l) = zero
          endif

        ENDDO
        avr = avr * amb
!
!      Correction for negative condensate!
!     if (advcld) then
!       do l=kd,k
!         if (qli(l) < zero) then
!           qoi(l) = qoi(l) + qli(l)
!           toi(l) = toi(l) - (alhl/cp) * qli(l)
!           qli(l) = zero
!         endif
!         if (qii(l) < zero) then
!           qoi(l) = qoi(l) + qii(l)
!           toi(l) = toi(l) - ((alhl+alhf)/cp) * qii(l)
!           qii(l) = zero
!         endif
!       enddo
!     endif

!
!
        TX1 = zero
        TX2 = zero
!
        IF (REVAP) THEN !     REEVAPORATION OF FALLING CONVECTIVE RAIN
!
         tem = zero
         do l=kd,kbl
           IF (L < IDH .or. (.not. DDFT)) THEN
             tem = tem + amb * rnn(l) * sigf(l)
           endif
         enddo
         tem = tem + amb * dof * sigf(kbl)
         tem = tem * (3600.0_kp/dt)
         tem1 = sqrt(max(one, min(100.0_kp,(6.25e10_kp/max(area,one))))) ! 20110530

         clfrac = max(ZERO, min(half, rknob*clf(tem)*tem1))
         cldfrd = clfrac
!
         DO L=KD,KBL         ! Testing on 20070926
!                                                 for L=KD,K
           IF (L >= IDH .AND. DDFT) THEN
             tem    = amb * sigf(l)
             TX2    = TX2 + tem * RNN(L)
             CLDFRD = MIN(tem*CLDFR(L), clfrac)
           ELSE
             TX1 = TX1 + AMB * RNN(L) * sigf(l)
           ENDIF
           tx4 = zfac * phil(l)
           tx4 = (one - tx4 * (one - half*tx4)) * afc
!
           IF (TX1 > zero .OR. TX2 > zero) THEN
             TEQ     = TOI(L)
             QEQ     = QOI(L)
             PL      = half * (PRL(L+1)+PRL(L))

             ST1     = MAX(ZERO, MIN(ONE, (TCR-TEQ)*TCRF))
             ST2     = ST1*ELFOCP + (one-ST1)*ELOCP

             CALL QSATCN ( TEQ,PL,QSTEQ,DQDT)
!
             DELTAQ = half * (QSTEQ*rhc_ls(l)-QEQ) / (one+ST2*DQDT)
!
             QEQ    = QEQ + DELTAQ
             TEQ    = TEQ - DELTAQ*ST2
!
             TEM1   = MAX(ZERO, MIN(ONE, (TCR-TEQ)*TCRF))
             TEM2   = TEM1*ELFOCP + (one-TEM1)*ELOCP

             CALL QSATCN ( TEQ,PL,QSTEQ,DQDT)
!
             DELTAQ = (QSTEQ*rhc_ls(l)-QEQ) / (one+TEM2*DQDT)
!
             QEQ    = QEQ + DELTAQ
             TEQ    = TEQ - DELTAQ*TEM2

             IF (QEQ > QOI(L)) THEN
               POTEVAP = (QEQ-QOI(L))*(PRL(L+1)-PRL(L))*GRAVCON

               tem4    = zero
               if (tx1 > zero)                                          &
     &         TEM4    = POTEVAP * (one - EXP( tx4*TX1**0.57777778_kp ))
               ACTEVAP = MIN(TX1, TEM4*CLFRAC)

               if (tx1 < rainmin*dt) actevap = min(tx1, potevap)
!
               tem4    = zero
               if (tx2 > zero)                                          &
     &         TEM4    = POTEVAP * (one - EXP( tx4*TX2**0.57777778_kp ))
               TEM4    = min(MIN(TX2, TEM4*CLDFRD), potevap-actevap)
               if (tx2 < rainmin*dt) tem4 = min(tx2, potevap-actevap)
!
               TX1     = TX1 - ACTEVAP
               TX2     = TX2 - TEM4
               ST1     = (ACTEVAP+TEM4) * PRI(L)
               QOI(L)  = QOI(L) + ST1
               QCU(L)  = QCU(L) + ST1
!

               ST1     = ST1 * ELOCP
               TOI(L)  = TOI(L) - ST1 
               TCU(L)  = TCU(L) - ST1
             ENDIF
           ENDIF
         ENDDO
!
          CUP = CUP + TX1 + TX2 + DOF * AMB * sigf(kbl)
        ELSE
          DO L=KD,K
            TX1 = TX1 + AMB * RNN(L) * sigf(l)
          ENDDO
          CUP = CUP + TX1 + DOF * AMB * sigf(kbl)
        ENDIF
!
!    Convective transport (mixing) of passive tracers
!
        if (NTRC > 0) then
          do l=kd,km1
            if (etz(l) /= zero) etzi(l) = one / etz(l)
          enddo
          DO N=1,NTRC        ! Tracer loop ; first two are u and v

            DO L=KD,K
              HOL(L) = ROI(L,N)
            ENDDO
!
            HCC     = RBL(N)
            HOD(KD) = HOL(KD)
!      Compute downdraft properties for the tracer
            DO L=KD1,K
              lm1 = l - 1
              ST1 = ONE - ALFIND(L)
              HB  = ALFIND(L)  * HOL(LM1) + ST1 * HOL(L)
              IF (ETZ(LM1) /= ZERO) THEN
                TEM = ETZI(LM1)
                IF (ETD(L)  > ETD(LM1)) THEN
                 HOD(L) = (ETD(LM1)*(HOD(LM1)-HOL(LM1))                 &
     &                   +  ETD(L)  *(HOL(LM1)-HB) +  ETZ(LM1)*HB) * TEM
                ELSE
                 HOD(L) = (ETD(LM1)*(HOD(LM1)-HB) + ETZ(LM1)*HB) * TEM
                ENDIF
               ELSE
                HOD(L) = HB
              ENDIF
            ENDDO
             
            DO L=KB1,KD,-1
              HCC = HCC + (ETA(L)-ETA(L+1))*HOL(L)
            ENDDO
!
!         Scavenging -- fscav   - fraction scavenged [km-1]
!                       delz    - distance from the entrainment to detrainment layer [km]
!                       fnoscav - the fraction not scavenged
!                                 following Liu et al. [JGR,2001] Eq 1

            if (FSCAV_(N) > zero) then
              DELZKM = ( PHIL(KD) - PHIH(KD1) ) *(onebg*0.001_kp)
              FNOSCAV = exp(- FSCAV_(N) * DELZKM)
            else
              FNOSCAV = one
            endif

            GMH(KD) = PRI(KD) * (HCC-ETA(KD)*HOL(KD)) * trcfac(kd,n)    &
     &                                                * FNOSCAV
            DO L=KD1,K
              if (FSCAV_(N) > zero) then
                DELZKM = ( PHIL(KD) - PHIH(L+1) ) *(onebg*0.001_kp)
                FNOSCAV = exp(- FSCAV_(N) * DELZKM)
              endif
              lm1      = l - 1
              ST1      = ONE - ALFINT(L,N+4)
              ST2      = ONE - ALFIND(L)
              HB       = ALFINT(L,N+4) * HOL(LM1) + ST1 * HOL(L)
              HBD      = ALFIND(L) * HOL(LM1) + ST2 * HOL(L)
              TEM5     = ETD(L)    * (HOD(L) - HBD)
              DH       = ETA(L)    * (HB - HOL(L))   * FNOSCAV + TEM5
              GMH(L  ) = DH * PRI(L) * trcfac(l,n)
              DH       = ETA(L)    * (HOL(LM1) - HB) * FNOSCAV - TEM5
              GMH(LM1) = GMH(LM1)  + DH * PRI(LM1) * trcfac(l,n)
            ENDDO
!
            st2 = zero
            DO L=KD,K
              ST1 = GMH(L)*AMB*sigf(l) + st2
              st3 = HOL(L) + st1
              st2 = st3 - trcmin(n) ! if trcmin is defined limit change
              if (st2 < zero) then
                ROI(L,N) = trcmin(n)
                RCU(L,N) = RCU(L,N) + ST1
                if (l < k)                                              &
     &          st2 = st2 * (prl(l+1)-prl(l))*pri(l+1) * (cmb2pa/grav)
              else
                ROI(L,N) = ST3
                RCU(L,N) = RCU(L,N) + ST1
                st2 = zero
              endif
              
            ENDDO
          ENDDO                             ! Tracer loop NTRC
        endif
      endif             ! amb > zero

      RETURN
      end subroutine cloud

      SUBROUTINE DDRFT(                                                 &
     &                  K,   KP1, KD                                    &
     &,                 TLA, ALFIND, wcbase                             &
     &,                 TOL, QOL, HOL, PRL, QST, HST, GAM, GAF          &
!    &,                 TOL, QOL, HOL, PRL, QST, HST, GAM, GAF, HBL, QBL&
     &,                 QRB, QRT, BUY, KBL, IDH, ETA, RNN, ETAI         &
     &,                 ALM, WFN, TRAIN, DDFT                           &
     &,                 ETD, HOD, QOD, EVP, DOF, CLDFRD, WCB            &
     &,                 GMS, GSD, GHD, wvlu)                   

!
!***********************************************************************
!******************** Cumulus Downdraft Subroutine *********************
!****************** Based on Cheng and Arakawa (1997)  ****** **********
!************************ SUBROUTINE DDRFT  ****************************
!*************************  October 2004  ******************************
!***********************************************************************
!***********************************************************************
!************* Shrinivas.Moorthi@noaa.gov (301) 683-3718 ***************
!***********************************************************************
!***********************************************************************
!23456789012345678901234567890123456789012345678901234567890123456789012
!
!===>  TOL(K)     INPUT   TEMPERATURE            KELVIN
!===>  QOL(K)     INPUT   SPECIFIC HUMIDITY      NON-DIMENSIONAL

!===>  PRL(KP1)   INPUT   PRESSURE @ EDGES       MB

!===>  K     INPUT   THE RISE & THE INDEX OF THE SUBCLOUD LAYER
!===>  KD    INPUT   DETRAINMENT LEVEL ( 1<= KD < K )          
!     
      IMPLICIT NONE
!
!  INPUT ARGUMENTS
!
      INTEGER K, KP1, KD, KBL
      real(kind=kind_phys) ALFIND(K), wcbase

      real(kind=kind_phys), dimension(kd:k)   ::  HOL, QOL, HST, QST    &
     &,                                           TOL, QRB, QRT, RNN    &
     &,                                           RNS, ETAI
      real(kind=kind_phys), dimension(kd:kp1) ::  GAF, BUY, GAM, ETA    &
     &,                                           PRL
!
!     real(kind=kind_phys)    HBL,     QBL,        PRIS                 &
!    &,                       TRAIN,   WFN,        ALM
!
!     TEMPORARY WORK SPACE
!
      real(kind=kind_phys), dimension(KD:K)   :: RNF, WCB,  EVP, STLT   &
     &,                                          GHD, GSD,  CLDFRD      &
     &,                                          GQW, QRPI, QRPS, BUD

      real(kind=kind_phys), dimension(KD:KP1) :: QRP, WVL, WVLU, ETD    &
     &,                                          HOD, QOD, ROR,  GMS

      real(kind=kind_phys) TL,     PL,     QL,      QS,   DQS,  ST1     &
     &,                    QQQ,            DEL_ETA, HB,   QB,   TB      &
     &,                    TEM,    TEM1,   TEM2,    TEM3, TEM4, ST2     &
     &,                    ERRMIN, ERRMI2, ERRH,    ERRW, ERRE, TEM5    &
     &,                    TEM6,   HBD,    QBD,     TX1,  TX2,  TX3     &
     &,                    TX4,    TX5,    TX6,     TX7,  TX8,  TX9     &
     &,                    WFN,    ALM,             AL2                 &
     &,                    TRAIN,  GMF,    ONPG,    CTLA, VTRM          &
     &,                    RPART,  QRMIN,  AA1,     BB1,  CC1,   DD1    &
!    &,                    WC2MIN, WCMIN,  WCBASE,  F2,   F3,    F5     &
     &,                    WC2MIN, WCMIN,           F2,   F3,    F5     &
     &,                    GMF1,   GMF5,   QRAF,    QRBF, del_tla       &
     &,                    TLA,    STLA,   CTL2,    CTL3                &
!    &,                    TLA,    STLA,   CTL2,    CTL3, ASIN          &
!    &,                    RNT,    RNB,    ERRQ,    RNTP, QRPF,  VTPF   &
     &,                    RNT,    RNB,    ERRQ,    RNTP                &
     &,                    EDZ,    DDZ,    CE,      QHS,  FAC,   FACG   &
     &,                    RSUM1,  RSUM2,  RSUM3,   CEE,  DOF,   DOFW
!    &,                    sialf

      INTEGER              I, L,  N, IX, KD1, II, kb1, IP1, JJ, ntla    &
     &,                    IT, KM1, KTEM, KK, KK1, LM1, LL, LP1         &
     &,                    IDW, IDH, IDN(K), idnm, itr
!
      parameter (ERRMIN=0.0001_kp, ERRMI2=0.1_kp*ERRMIN)
!     parameter (ERRMIN=0.00001, ERRMI2=0.1*ERRMIN)
!
!     real (kind=kind_phys), parameter :: PIINV=one/PI, pio2=half*pi
!
      parameter (ONPG=one+half, GMF=one/ONPG, RPART=zero)
!     parameter (ONPG=1.0+0.5, GMF=1.0/ONPG, RPART=1.0)
!     parameter (ONPG=1.0+0.5, GMF=1.0/ONPG, RPART=0.5)
!     PARAMETER (AA1=1.0, BB1=1.5, CC1=1.1, DD1=0.85, F3=CC1, F5=2.5)
!     PARAMETER (AA1=2.0, BB1=1.5, CC1=1.1, DD1=0.85, F3=CC1, F5=2.5)
      PARAMETER (AA1=1.0_kp, BB1=1.0_kp, CC1=1.0_kp, DD1=1.0_kp,        &
     &           F3=CC1,  F5=1.0_kp)
      parameter (QRMIN=1.0e-6_kp, WC2MIN=0.01_kp, GMF1=GMF/AA1, GMF5=GMF/F5)
!     parameter (QRMIN=1.0E-6, WC2MIN=1.00, GMF1=GMF/AA1, GMF5=GMF/F5)
      parameter (WCMIN=sqrt(wc2min))
!     parameter (sialf=0.5)
!
      integer, parameter :: itrmu=25,  itrmd=25                         &
     &,                     itrmin=15, itrmnd=12, numtla=2

!     uncentering for vvel in dd
      real(kind=kind_phys), parameter :: ddunc1=0.25_kp                 &
     &,                                  ddunc2=one-ddunc1              &
!    &,                                  ddunc1=0.4,  ddunc2=one-ddunc1 &
!    &,                                  ddunc1=0.3,  ddunc2=one-ddunc1 &
     &,                      VTPEXP=-0.3636_kp                          &
     &,                      VTP=36.34_kp*SQRT(1.2_kp)*(0.001_kp)**0.1364_kp
!
!     real(kind=kind_phys) EM(K*K), ELM(K)
      real(kind=kind_phys) ELM(K), AA(KD:K,KD:KP1), QW(KD:K,KD:K)       &
     &,                    VT(2),  VRW(2), TRW(2), QA(3), WA(3)

      LOGICAL SKPUP, cnvflg, DDFT, UPDRET, DDLGK

!***********************************************************************


      KD1    = KD + 1
      KM1    = K  - 1
      KB1    = KBL - 1
!
!     VTP    = 36.34*SQRT(1.2)* (0.001)**0.1364
!     VTPEXP = -0.3636
!     PIINV  = 1.0 / PI
!     PICON  = PIO2 * ONEBG
!
!     Compute Rain Water Budget of the Updraft (Cheng and Arakawa, 1997)
!
      CLDFRD = zero
      RNTP   = zero
      DOF    = zero
      ERRQ   = 10.0_kp
      RNB    = zero
      RNT    = zero
      TX2    = PRL(KBL)
!
      TX1      = (PRL(KD) + PRL(KD1)) * half
      ROR(KD)  = CMPOR*TX1 / (TOL(KD)*(one+NU*QOL(KD)))
!     GMS(KD)  = VTP * ROR(KD) ** VTPEXP
      GMS(KD)  = VTP * VTPF(ROR(KD))
!
      QRP(KD)  = QRMIN
!
      TEM      = TOL(K) * (one + NU * QOL(K))
      ROR(KP1) = half * CMPOR * (PRL(KP1)+PRL(K)) / TEM
      GMS(KP1) = VTP * VTPF(ROR(KP1))
      QRP(KP1) = QRMIN
!
      kk = kbl
      DO L=KD1,K
        TEM = half * (TOL(L)+TOL(L-1))                                  &
     &      * (one + (half*NU) * (QOL(L)+QOL(L-1)))
        ROR(L) = CMPOR * PRL(L) / TEM
!       GMS(L) = VTP * ROR(L) ** VTPEXP
        GMS(L) = VTP * VTPF(ROR(L))
        QRP(L) = QRMIN
        if (buy(l) <= zero .and. kk == KBL) then
          kk = l
        endif
      ENDDO
      if (kk /= kbl) then
        do l=kk,kbl
          buy(l) = 0.9_kp * buy(l-1)
        enddo
      endif
!
      do l=kd,k
        qrpi(l) = buy(l)
      enddo
      do l=kd1,kb1
        buy(l) = 0.25_kp * (qrpi(l-1)+qrpi(l)+qrpi(l)+qrpi(l+1))
      enddo
      
!
!     CALL ANGRAD(TX1, ALM, STLA, CTL2, AL2, PI, TLA, TX2, WFN, TX3)
      tx1 = 1000.0_kp + tx1 - prl(kp1)
!     CALL ANGRAD(TX1, ALM,  AL2, TLA, TX2, WFN, TX3)
      CALL ANGRAD(TX1, ALM,  AL2, TLA)
!
!    Following Ucla approach for rain profile
!
      F2      = (BB1+BB1)*ONEBG/(PI*0.2_kp)
!     WCMIN   = SQRT(WC2MIN)
!     WCBASE  = WCMIN
!
!     del_tla = TLA * 0.2
!     del_tla = TLA * 0.25
      del_tla = TLA * 0.3_kp
      TLA     = TLA - DEL_TLA
!
      DO L=KD,K
        RNF(L)   = zero
        RNS(L)   = zero
        STLT(L)  = zero
        GQW(L)   = zero
        QRP(L)   = QRMIN
        DO N=KD,K
          QW(N,L) = zero
        ENDDO
      ENDDO
!     DO L=KD,KP1
!       WVL(L)   = zero
!     ENDDO
!
!-----QW(N,L) = D(W(N)*W(N))/DQR(L)
!
      KK = KBL
      QW(KD,KD) = -QRB(KD)  * GMF1
      GHD(KD)   = ETA(KD)   * ETA(KD)
      GQW(KD)   = QW(KD,KD) * GHD(KD)
      GSD(KD)   = ETAI(KD)  * ETAI(KD)
!
      GQW(KK)   = - QRB(KK-1) * (GMF1+GMF1)
!
      WCB(KK)   = WCBASE * WCBASE

      TX1       = WCB(KK)
      GSD(KK)   = one
      GHD(KK)   = one
!
      TEM       = GMF1 + GMF1
      DO L=KB1,KD1,-1
        GHD(L)  = ETA(L)  * ETA(L)
        GSD(L)  = ETAI(L) * ETAI(L)
        GQW(L)  = - GHD(L) * (QRB(L-1)+QRT(L)) * TEM
        QW(L,L) = - QRT(L) * TEM
!
        st1     = half * (eta(l) + eta(l+1))
        TX1     = TX1 + BUY(L) * TEM * (qrb(l)+qrt(l)) * st1 * st1
        WCB(L)  = TX1 * GSD(L)
      ENDDO
!
      TEM1        = (QRB(KD) + QRT(KD1) + QRT(KD1)) * GMF1
      GQW(KD1)    = - GHD(KD1) * TEM1
      QW(KD1,KD1) = - QRT(KD1) * TEM
      st1         = half * (eta(kd) + eta(kd1))
      WCB(KD)     = (TX1 + BUY(KD)*TEM*qrb(kd)*st1*st1) * GSD(KD)
!
      DO L=KD1,KBL
        DO N=KD,L-1
          QW(N,L) = GQW(L) * GSD(N)
        ENDDO
      ENDDO
      QW(KBL,KBL) = zero
!
      do ntla=1,numtla    ! numtla is the the maximimu number of tilting angle tries
                          ! ------
!       if (errq < 1.0 .or. tla > 45.0) cycle
        if (errq < 0.1_kp .or. tla > 45.0_kp) cycle
!
        tla  = tla + del_tla
        STLA = SIN(TLA*deg2rad)     ! sine of tilting angle
        CTL2 = one - STLA * STLA    ! cosine square of tilting angle
!
        STLA = F2        * STLA * AL2
        CTL2 = DD1       * CTL2
        CTL3 = 0.1364_kp * CTL2
!
        DO L=KD,K
          RNF(L)  = zero
          STLT(L) = zero
          QRP(L)  = QRMIN
        ENDDO
        DO L=KD,KP1
          WVL(L)  = zero
        ENDDO
        WVL(KBL)  = WCBASE
        STLT(KBL) = one / WCBASE
!
        DO L=KD,KP1
          DO N=KD,K
            AA(N,L) = zero
          ENDDO
        ENDDO
!
        SKPUP = .FALSE.
!
        DO ITR=1,ITRMU               ! Rain Profile Iteration starts!
          IF (.NOT. SKPUP) THEN
!            wvlu = wvl
!
!-----CALCULATING THE VERTICAL VELOCITY
!
            TX1      = zero
            QRPI(KBL) = one / QRP(KBL)
            DO L=KB1,KD,-1
              TX1 = TX1    + QRP(L+1)*GQW(L+1)
              ST1 = WCB(L) + QW(L,L)*QRP(L) + TX1*GSD(L)
!             if (st1 > wc2min) then
              if (st1 > zero) then
                WVL(L) = max(ddunc1*SQRT(ST1) + ddunc2*WVL(L), wcmin)
!               WVL(L) = SQRT(ST1)
!               WVL(L) = max(half * (SQRT(ST1) + WVL(L)), wcmin)
!               qrp(l) = half*((wvl(l)*wvl(l)-wcb(l)-tx1*gsd(l))/qw(l,l)&
!    &                      + qrp(l))
              else

!               wvl(l) = 0.5*(wcmin+wvl(l))
!               wvl(l) = max(half*(wvl(l) + wvl(l+1)), wcmin)
                wvl(l) = max(wvl(l),wcmin)
                qrp(l) = (wvl(l)*wvl(l) - wcb(l) - tx1*gsd(l))/qw(l,l)
!               qrp(l) = half*((wvl(l)*wvl(l)-wcb(l)-tx1*gsd(l))/qw(l,l)&
!    &                      + qrp(l))
              endif
              qrp(l) = max(qrp(l), qrmin)

              STLT(L) = one / WVL(L)
              QRPI(L) = one / QRP(L)
            ENDDO
!
!-----CALCULATING TRW, VRW AND OF
!
!           VT(1)   = GMS(KD) * QRP(KD)**0.1364
            VT(1)   = GMS(KD) * QRPF(QRP(KD))
            TRW(1)  = ETA(KD) * QRP(KD) * STLT(KD)
            TX6     = TRW(1) * VT(1)
            VRW(1)  = F3*WVL(KD) - CTL2*VT(1)
            BUD(KD) = STLA * TX6 * QRB(KD) * half
            RNF(KD) = BUD(KD)
            DOF     = 1.1364_kp * BUD(KD) * QRPI(KD)
            DOFW    = -BUD(KD)  * STLT(KD)
!
            RNT     = TRW(1) * VRW(1)
            TX2     = zero
            TX4     = zero
            RNB     = RNT
            TX1     = half
            TX8     = zero
!
            IF (RNT >= zero) THEN
              TX3 = (RNT-CTL3*TX6) * QRPI(KD)
              TX5 = CTL2 * TX6 * STLT(KD)
            ELSE
              TX3 = zero
              TX5 = zero
              RNT = zero
              RNB = zero
            ENDIF
!
            DO L=KD1,KB1
              KTEM    = MAX(L-2, KD)
              LL      = L - 1
! 
!             VT(2)   = GMS(L) * QRP(L)**0.1364
              VT(2)   = GMS(L) * QRPF(QRP(L))
              TRW(2)  = ETA(L) * QRP(L) * STLT(L)
              VRW(2)  = F3*WVL(L) - CTL2*VT(2)
              QQQ     = STLA * TRW(2) * VT(2)
              ST1     = TX1  * QRB(LL)
              BUD(L)  = QQQ * (ST1 + QRT(L))
!
              QA(2)   = DOF
              WA(2)   = DOFW
              DOF     = 1.1364_kp * BUD(L) * QRPI(L)
              DOFW    = -BUD(L)   * STLT(L)
!
              RNF(LL) = RNF(LL) + QQQ * ST1
              RNF(L)  =           QQQ * QRT(L)
!
              TEM3    = VRW(1) + VRW(2)
              TEM4    = TRW(1) + TRW(2)
!
              TX6     = pt25 * TEM3 * TEM4
              TEM4    = TEM4 * CTL3
!
!-----BY QR ABOVE
!
!             TEM1    = pt25*(TRW(1)*TEM3 - TEM4*VT(1))*TX7
              TEM1    = pt25*(TRW(1)*TEM3 - TEM4*VT(1))*QRPI(LL)
              ST1     = pt25*(TRW(1)*(CTL2*VT(1)-VRW(2))                &
     &                      * STLT(LL) + F3*TRW(2))
!-----BY QR BELOW
              TEM2    = pt25*(TRW(2)*TEM3 - TEM4*VT(2))*QRPI(L)
              ST2     = pt25*(TRW(2)*(CTL2*VT(2)-VRW(1))                &
     &                      * STLT(L)  + F3*TRW(1))
!
!      From top to  the KBL-2 layer
!
              QA(1)   = TX2
              QA(2)   = QA(2) + TX3 - TEM1
              QA(3)   = -TEM2
!
              WA(1)   = TX4
              WA(2)   = WA(2) + TX5 - ST1
              WA(3)   = -ST2
!
              TX2     = TEM1
              TX3     = TEM2
              TX4     = ST1
              TX5     = ST2
!
              VT(1)   = VT(2)
              TRW(1)  = TRW(2)
              VRW(1)  = VRW(2)
!
              IF (WVL(KTEM) == WCMIN) WA(1) = zero
              IF (WVL(LL)   == WCMIN) WA(2) = zero
              IF (WVL(L)    == WCMIN) WA(3) = zero
              DO N=KTEM,KBL
                AA(LL,N) = (WA(1)*QW(KTEM,N) * STLT(KTEM)               &
     &                   +  WA(2)*QW(LL,N)   * STLT(LL)                 &
     &                   +  WA(3)*QW(L,N)    * STLT(L) ) * half
              ENDDO
              AA(LL,KTEM) = AA(LL,KTEM) + QA(1)
              AA(LL,LL)   = AA(LL,LL)   + QA(2)
              AA(LL,L)    = AA(LL,L)    + QA(3)
              BUD(LL)     = (TX8 + RNN(LL)) * half                      &
     &                      - RNB + TX6 - BUD(LL)
              AA(LL,KBL+1) = BUD(LL)
              RNB = TX6
              TX1 = one
              TX8 = RNN(LL)
            ENDDO
            L  = KBL
            LL = L - 1
!           VT(2)   = GMS(L) * QRP(L)**0.1364
            VT(2)   = GMS(L) * QRPF(QRP(L))
            TRW(2)  = ETA(L) * QRP(L) * STLT(L)
            VRW(2)  = F3*WVL(L) - CTL2*VT(2)
            ST1     = STLA * TRW(2) * VT(2) * QRB(LL)
            BUD(L)  = ST1

            QA(2)   = DOF
            WA(2)   = DOFW
            DOF     = 1.1364_kp * BUD(L) * QRPI(L)
            DOFW    = -BUD(L)   * STLT(L)
!
            RNF(LL) = RNF(LL) + ST1
!
            TEM3    = VRW(1) + VRW(2)
            TEM4    = TRW(1) + TRW(2)
!
            TX6     = pt25 * TEM3 * TEM4
            TEM4    = TEM4 * CTL3
!
!-----BY QR ABOVE
!
            TEM1    = pt25*(TRW(1)*TEM3 - TEM4*VT(1))*QRPI(LL)
            ST1     = pt25*(TRW(1)*(CTL2*VT(1)-VRW(2))                  &
     &                    * STLT(LL) + F3*TRW(2))
!-----BY QR BELOW
            TEM2    = pt25*(TRW(2)*TEM3 - TEM4*VT(2))*QRPI(L)
            ST2     = pt25*(TRW(2)*(CTL2*VT(2)-VRW(1))                  &
     &                    * STLT(L)  + F3*TRW(1))
!
!      For the layer next to the top of the boundary layer
!
            QA(1)   = TX2
            QA(2)   = QA(2) + TX3 - TEM1
            QA(3)   = -TEM2
!
            WA(1)   = TX4
            WA(2)   = WA(2) + TX5 - ST1
            WA(3)   = -ST2
!
            TX2     = TEM1
            TX3     = TEM2
            TX4     = ST1
            TX5     = ST2
!
            IDW     = MAX(L-2, KD)
!
            IF (WVL(IDW) == WCMIN) WA(1) = zero
            IF (WVL(LL)  == WCMIN) WA(2) = zero
            IF (WVL(L)   == WCMIN) WA(3) = zero
!
            KK = IDW
            DO N=KK,L
              AA(LL,N) = (WA(1)*QW(KK,N) * STLT(KK)                     &
     &                 +  WA(2)*QW(LL,N) * STLT(LL)                     &
     &                 +  WA(3)*QW(L,N)  * STLT(L) ) * half

            ENDDO
!
            AA(LL,IDW) = AA(LL,IDW) + QA(1)
            AA(LL,LL)  = AA(LL,LL)  + QA(2)
            AA(LL,L)   = AA(LL,L)   + QA(3)
            BUD(LL)    = (TX8+RNN(LL)) * half - RNB + TX6 - BUD(LL)
!
            AA(LL,L+1) = BUD(LL)
!
            RNB        = TRW(2) * VRW(2)
!
!      For the top of the boundary layer
!
            IF (RNB < zero) THEN
               KK    = KBL
               TEM   = VT(2) * TRW(2)
               QA(2) = (RNB - CTL3*TEM) * QRPI(KK)
               WA(2) = CTL2 * TEM * STLT(KK)
            ELSE
               RNB   = zero
               QA(2) = zero
               WA(2) = zero
            ENDIF
!
            QA(1) = TX2
            QA(2) = DOF + TX3 - QA(2)
            QA(3) = zero
!
            WA(1) = TX4
            WA(2) = DOFW + TX5 - WA(2)
            WA(3) = zero
!
            KK = KBL
            IF (WVL(KK-1) == WCMIN) WA(1) = zero
            IF (WVL(KK)   == WCMIN) WA(2) = zero
!
            DO II=1,2
               N = KK + II - 2
               AA(KK,N) = (WA(1)*QW(KK-1,N) * STLT(KK-1)                &
     &                  +  WA(2)*QW(KK,N)   * STLT(KK)) * half
            ENDDO
            FAC = half
            LL  = KBL
            L   = LL + 1
            LM1 = LL - 1
            AA(LL,LM1)  = AA(LL,LM1) + QA(1)
            AA(LL,LL)   = AA(LL,LL)  + QA(2)
            BUD(LL)     = half*RNN(LM1) - TX6 + RNB - BUD(LL)
            AA(LL,LL+1) = BUD(LL)
!
!-----SOLVING THE BUDGET EQUATIONS FOR DQR
!
            DO L=KD1,KBL
              LM1  = L - 1
              cnvflg = ABS(AA(LM1,LM1)) < ABS(AA(L,LM1))
              DO  N=LM1,KBL+1
                 IF (cnvflg) THEN
                    TX1       = AA(LM1,N)
                    AA(LM1,N) = AA(L,N)
                    AA(L,N)   = TX1
                 ENDIF
              ENDDO
              TX1 = AA(L,LM1) / AA(LM1,LM1)
              DO  N=L,KBL+1
                AA(L,N) = AA(L,N) - TX1 * AA(LM1,N)
              ENDDO
            ENDDO     
!
!-----BACK SUBSTITUTION AND CHECK IF THE SOLUTION CONVERGES
!
            KK = KBL
            KK1 = KK + 1
            AA(KK,KK1) = AA(KK,KK1) / AA(KK,KK)      !   Qr correction !
            TX2        = ABS(AA(KK,KK1)) * QRPI(KK)  !   Error Measure !
!
            KK = KBL + 1
            DO L=KB1,KD,-1
               LP1   = L + 1
               TX1  = zero
               DO N=LP1,KBL
                 TX1  = TX1 + AA(L,N) * AA(N,KK)
               ENDDO
               AA(L,KK) = (AA(L,KK) - TX1) / AA(L,L)       ! Qr correction !
               TX2      = MAX(TX2, ABS(AA(L,KK))*QRPI(L))  ! Error Measure !
            ENDDO
!
!           tem = 0.5
            if (tx2 > one .and. abs(errq-tx2) > 0.1_kp) then
              tem = half
!!          elseif (tx2 < 0.1) then
!!            tem = 1.2
            else
              tem = one
            endif
!
            DO L=KD,KBL
!              QRP(L) = MAX(QRP(L)+AA(L,KBL+1), QRMIN)
               QRP(L) = MAX(QRP(L)+AA(L,KBL+1)*tem, QRMIN)
            ENDDO
!
            IF (ITR < ITRMIN) THEN
               TEM = ABS(ERRQ-TX2) 
               IF (TEM >= ERRMI2 .AND. TX2 >= ERRMIN) THEN 
                 ERRQ  = TX2                              ! Further iteration !
               ELSE 
                 SKPUP = .TRUE.                           ! Converges      !
                 ERRQ  = zero                             ! Rain profile exists!
               ENDIF 
            ELSE
               TEM = ERRQ - TX2
!              IF (TEM < ZERO .AND. ERRQ > 0.1_kp) THEN
               IF (TEM < ZERO .AND. ERRQ > 0.5_kp) THEN
!              IF (TEM < ZERO .and.                                    &
!    &            (ntla < numtla .or. ERRQ > 0.5_kp)) THEN
                 SKPUP = .TRUE.                           ! No convergence !
                 ERRQ  = 10.0_kp                          ! No rain profile!
!!!!           ELSEIF (ABS(TEM) < ERRMI2 .OR. TX2 < ERRMIN) THEN
               ELSEIF (TX2 < ERRMIN) THEN
                 SKPUP = .TRUE.                           ! Converges      !
                 ERRQ = zero                              ! Rain profile exists!
               elseif (tem < zero .and. errq < 0.1_kp) then
                 skpup = .true.
!                if (ntla == numtla .or. tem > -0.003) then
                   errq  = zero
!                else
!                  errq = 10.0
!                endif
               ELSE
                 ERRQ = TX2                               ! Further iteration !

!              if (itr == itrmu .and. ERRQ > ERRMIN*10                  &
!    &            .and. ntla == 1) ERRQ = 10.0 
               ENDIF
            ENDIF
!
          ENDIF                                           ! SKPUP  ENDIF!
!
        ENDDO                                             ! End of the ITR Loop!!
!
        IF (ERRQ < 0.1_kp) THEN
          DDFT = .TRUE.
          RNB  = - RNB
!         do l=kd1,kb1-1
!           if (wvl(l)-wcbase < 1.0E-9) ddft = .false.
!         enddo
        ELSE
          DDFT = .FALSE.
        ENDIF

      enddo                                          ! End of ntla loop
!
!     Caution !! Below is an adjustment to rain flux to maintain
!                conservation of precip!
!
      IF (DDFT) THEN
        TX1 = zero
        DO L=KD,KB1
          TX1 = TX1 + RNF(L)
        ENDDO
        TX1 = TRAIN / (TX1+RNT+RNB)
        IF (ABS(TX1-one) < 0.2_kp) THEN
          RNT = MAX(RNT*TX1,ZERO)
          RNB = RNB * TX1
          DO L=KD,KB1
            RNF(L) = RNF(L) * TX1
          ENDDO
!      rain flux adjustment is over

        ELSE
          DDFT = .FALSE.
          ERRQ = 10.0_kp
        ENDIF
      ENDIF
!
      DOF = zero
      IF (.NOT. DDFT) then
        wvlu(kd:kp1) = zero
        RETURN     ! Rain profile did not converge!
                   ! No down draft for this case - rerurn
                   ! ------------------------------------
!
      else         ! rain profile converged - do downdraft calculation
                   ! ------------------------------------------------

        wvlu(kd:kp1) = wvl(kd:kp1) ! save updraft vertical velocity for output

!
!     Downdraft calculation begins
!     ----------------------------
!
        DO L=KD,K
          WCB(L) = zero
        ENDDO
!
        ERRQ  = 10.0_kp
! At this point stlt contains inverse of updraft vertical velocity 1/Wu.

        KK = MAX(KB1,KD1)
        DO L=KK,K
          STLT(L) = STLT(L-1)
        ENDDO
        TEM = stla / BB1       ! this is 2/(pi*radius*grav)
!
        DO L=KD,K
          IF (L <= KBL) THEN
            STLT(L) = ETA(L) * STLT(L) * TEM / ROR(L)
          ELSE
            STLT(L) = zero
          ENDIF
        ENDDO

        rsum1 = zero
        rsum2 = zero
!
        IDN(:) = idnmax
        DO L=KD,KP1
          ETD(L)  = zero
          WVL(L)  = zero
!         QRP(L)  = zero
        ENDDO
        DO L=KD,K
          EVP(L)   = zero
          BUY(L)   = zero
          QRP(L+1) = zero
        ENDDO
        HOD(KD)  = HOL(KD)
        QOD(KD)  = QOL(KD)
        TX1      = zero
!!!     TX1      = STLT(KD)*QRB(KD)*ONE              ! sigma at the top
!       TX1      = MIN(STLT(KD)*QRB(KD)*ONE, ONE)    ! sigma at the top
!       TX1      = MIN(STLT(KD)*QRB(KD)*0.5, ONE)    ! sigma at the top
        RNTP     = zero
        TX5      = TX1
        QA(1)    = zero
!
!       Here we assume RPART of detrained rain RNT goes to Pd
!
        IF (RNT > zero) THEN
          if (TX1 > zero) THEN
            QRP(KD) = (RPART*RNT / (ROR(KD)*TX1*GMS(KD)))               &
     &                                          ** (one/1.1364_kp)
          else
            tx1 = RPART*RNT / (ROR(KD)*GMS(KD)*QRP(KD)**1.1364_kp)
          endif
          RNTP    = (one - RPART) * RNT
          BUY(KD) = - ROR(KD) * TX1 * QRP(KD)
        ELSE
          QRP(KD) = zero
        ENDIF
!
!     L-loop for the downdraft iteration from KD1 to KP1 (bottom surface)
!
!       BUD(KD) = ROR(KD)
        idnm = 1
        DO L=KD1,KP1

          QA(1) = zero
          ddlgk = idn(idnm) == idnmax
          if (.not. ddlgk) cycle
          IF (L <= K) THEN
            ST1   = one - ALFIND(L)
            WA(1) = ALFIND(L)*HOL(L-1) + ST1*HOL(L)
            WA(2) = ALFIND(L)*QOL(L-1) + ST1*QOL(L)
            WA(3) = ALFIND(L)*TOL(L-1) + ST1*TOL(L)
            QA(2) = ALFIND(L)*HST(L-1) + ST1*HST(L)
            QA(3) = ALFIND(L)*QST(L-1) + ST1*QST(L)
          ELSE
            WA(1) = HOL(K)
            WA(2) = QOL(K)
            WA(3) = TOL(K)
            QA(2) = HST(K)
            QA(3) = QST(K)
          ENDIF
!
          FAC = two
          IF (L == KD1) FAC = one

          FACG    = FAC * half * GMF5     !  12/17/97
!
!         DDLGK   =  IDN(idnm) == 99

          BUD(KD) = ROR(L)

          TX1    = TX5
          WVL(L) = MAX(WVL(L-1),ONE_M1)

          QRP(L) = MAX(QRP(L-1),QRP(L))
!
!         VT(1)  = GMS(L-1) * QRP(L-1) ** 0.1364
          VT(1)  = GMS(L-1) * QRPF(QRP(L-1))
          RNT    = ROR(L-1) * (WVL(L-1)+VT(1))*QRP(L-1)
!
!         TEM    = MAX(ALM, 2.5E-4) * MAX(ETA(L), 1.0)
          TEM    = MAX(ALM,ONE_M6) * MAX(ETA(L), ONE)
!         TEM    = MAX(ALM, 1.0E-5) * MAX(ETA(L), 1.0)
          TRW(1) = PICON*TEM*(QRB(L-1)+QRT(L-1))
          TRW(2) = one / TRW(1)
!
          VRW(1) = half * (GAM(L-1) + GAM(L))
          VRW(2) = one / (VRW(1) + VRW(1))
!
          TX4    =  (QRT(L-1)+QRB(L-1))*(ONEBG*FAC*500.0_kp*EKNOB)
!
          DOFW   = one / (WA(3) * (one + NU*WA(2)))      !  1.0 / TVbar!
!
          ETD(L) = ETD(L-1)
          HOD(L) = HOD(L-1)
          QOD(L) = QOD(L-1)
!
          ERRQ   = 10.0_kp

!
          IF (L <= KBL) THEN
            TX3 = STLT(L-1) * QRT(L-1) * (half*FAC)
            TX8 = STLT(L)   * QRB(L-1) * (half*FAC)
            TX9 = TX8 + TX3
          ELSE
            TX3 = zero
            TX8 = zero
            TX9 = zero
          ENDIF
!
          TEM  = WVL(L-1) + VT(1)
          IF (TEM > zero) THEN
            TEM1 = one / (TEM*ROR(L-1))
            TX3 = VT(1) * TEM1 * ROR(L-1) * TX3
            TX6 = TX1 * TEM1
          ELSE
            TX6 = one
          ENDIF
!
          IF (L == KD1) THEN
            IF (RNT > zero) THEN
              TEM    = MAX(QRP(L-1),QRP(L))
              WVL(L) = TX1 * TEM * QRB(L-1)*(FACG*5.0_kp)
            ENDIF
            WVL(L) = MAX(ONE_M2, WVL(L))
            TRW(1) = TRW(1) * half
            TRW(2) = TRW(2) + TRW(2)
          ELSE
            IF (DDLGK) EVP(L-1) = EVP(L-2)
          ENDIF
!
!       No downdraft above level IDH
!

          IF (L < IDH) THEN

            ETD(L)   = zero
            HOD(L)   = WA(1)
            QOD(L)   = WA(2)
            EVP(L-1) = zero
            WVL(L)   = zero
            QRP(L)   = zero
            BUY(L)   = zero
            TX5      = TX9
            ERRQ     = zero
            RNTP     = RNTP + RNT * TX1
            RNT      = zero
            WCB(L-1) = zero

!         ENDIF
!         BUD(KD) = ROR(L)
!
!       Iteration loop for a given level L begins
!
          else
            DO ITR=1,ITRMD
!
!             cnvflg =  DDLGK .AND. (ERRQ > ERRMIN)
              cnvflg =  ERRQ > ERRMIN
              IF (cnvflg) THEN
!
!               VT(1) = GMS(L) * QRP(L) ** 0.1364
                VT(1) = GMS(L) * QRPF(QRP(L))
                TEM   =  WVL(L) + VT(1)
!
                IF (TEM > zero) THEN
                  ST1  = ROR(L) * TEM * QRP(L) + RNT
                  IF (ST1 /= zero) ST1 = two * EVP(L-1) / ST1
                  TEM1 = one / (TEM*ROR(L))
                  TEM2 = VT(1) * TEM1 * ROR(L) * TX8
                ELSE
                  TEM1 = zero
                  TEM2 = TX8
                  ST1  = zero
                ENDIF
!
                st2 = tx5
                TEM = ROR(L)*WVL(L) - ROR(L-1)*WVL(L-1)
                if (tem > zero) then
                  TX5 = (TX1 - ST1 + TEM2 + TX3)/(one+tem*tem1)
                else
                  TX5 = TX1 - tem*tx6 - ST1 + TEM2 + TX3
                endif
                TX5   = MAX(TX5,ZERO)
                tx5   = half * (tx5 + st2)
!
!               qqq = 1.0 + tem * tem1 * (1.0 - sialf)
!
!               if (qqq > 0.0) then
!                 TX5   = (TX1 - sialf*tem*tx6 - ST1 + TEM2 + TX3) / qqq
!               else
!                 TX5   = (TX1 - tem*tx6 - ST1 + TEM2 + TX3)
!               endif
!
                TEM1   = ETD(L)
                ETD(L) = ROR(L) * TX5 * MAX(WVL(L),ZERO)
!
                if (etd(l) > zero) etd(l) = half * (etd(l) + tem1)
!

                DEL_ETA = ETD(L) - ETD(L-1)

!                 TEM       = DEL_ETA * TRW(2)
!                 TEM2      = MAX(MIN(TEM, 1.0), -1.0)
!                 IF (ABS(TEM) > 1.0 .AND. ETD(L) > 0.0 ) THEN
!                   DEL_ETA = TEM2 * TRW(1)
!                   ETD(L)  = ETD(L-1) + DEL_ETA
!                 ENDIF
!                 IF (WVL(L) > 0.0) TX5 = ETD(L) / (ROR(L)*WVL(L))
!
                  ERRE  = ETD(L) - TEM1
!
                  tem  = max(abs(del_eta), trw(1))
                  tem2 = del_eta / tem
                  TEM1 = SQRT(MAX((tem+DEL_ETA)*(tem-DEL_ETA),ZERO))
!                 TEM1 = SQRT(MAX((TRW(1)+DEL_ETA)*(TRW(1)-DEL_ETA),0.0))

                  EDZ  = (half + ASIN(TEM2)*PIINV)*DEL_ETA + TEM1*PIINV

                DDZ   = EDZ - DEL_ETA
                WCB(L-1) = ETD(L) + DDZ
!
                TEM1  = HOD(L)
                IF (DEL_ETA > zero) THEN
                  QQQ    = one / (ETD(L) + DDZ)
                  HOD(L) = (ETD(L-1)*HOD(L-1) + DEL_ETA*HOL(L-1)        &
     &                                              + DDZ*WA(1)) * QQQ
                  QOD(L) = (ETD(L-1)*QOD(L-1) + DEL_ETA*QOL(L-1)        &
     &                                              + DDZ*WA(2)) * QQQ
                ELSEif((ETD(L-1) + EDZ) > zero) then
                  QQQ    = one / (ETD(L-1) + EDZ)
                  HOD(L) = (ETD(L-1)*HOD(L-1) + EDZ*WA(1)) * QQQ
                  QOD(L) = (ETD(L-1)*QOD(L-1) + EDZ*WA(2)) * QQQ
                ENDIF
                ERRH  = HOD(L) - TEM1
                ERRQ  = ABS(ERRH/HOD(L))  + ABS(ERRE/MAX(ETD(L),ONE_M5))

                DOF   = DDZ
                VT(2) = QQQ
!
                DDZ  = DOF
                TEM4 = QOD(L)
                TEM1 = VRW(1)
!
                QHS  = QA(3) + half * (GAF(L-1)+GAF(L)) * (HOD(L)-QA(2))
!
!                                           First iteration       !
!
                ST2  = PRL(L) * (QHS + TEM1 * (QHS-QOD(L)))
                TEM2 = ROR(L) * QRP(L)
                CALL QRABF(TEM2,QRAF,QRBF)
                TEM6 = TX5 * (1.6_kp + 124.9_kp * QRAF) * QRBF * TX4
!
                CE   = TEM6 * ST2 / ((5.4e5_kp*ST2 + 2.55e6_kp)*(ETD(L)+DDZ))
!
                TEM2   = - ((one+TEM1)*(QHS+CE) + TEM1*QOD(L))
                TEM3   = (one + TEM1) * QHS * (QOD(L)+CE)
                TEM    = MAX(TEM2*TEM2 - four*TEM1*TEM3,ZERO)
                QOD(L) = MAX(TEM4, (- TEM2 - SQRT(TEM)) * VRW(2))
!
!
!                                            second iteration   !
!
                ST2  = PRL(L) * (QHS + TEM1 * (QHS-QOD(L)))
                CE   = TEM6 * ST2 / ((5.4e5_kp*ST2 + 2.55e6_kp)*(ETD(L)+DDZ))
!               CEE  = CE * (ETD(L)+DDZ)
!


                TEM2   = - ((one+TEM1)*(QHS+CE) + TEM1*tem4)
                TEM3   = (one + TEM1) * QHS * (tem4+CE)
                TEM    = MAX(TEM2*TEM2 - four*TEM1*TEM3,ZERO)
                QOD(L) = MAX(TEM4, (- TEM2 - SQRT(TEM)) * VRW(2))
!                                              Evaporation in Layer L-1
!
                EVP(L-1) = (QOD(L)-TEM4) * (ETD(L)+DDZ)
!                                              Calculate Pd (L+1/2)
                QA(1)    = TX1*RNT + RNF(L-1) - EVP(L-1)
!
                if (qa(1) > zero) then
                  IF (ETD(L) > zero) THEN
                    TEM    = QA(1) / (ETD(L)+ROR(L)*TX5*VT(1))
                    QRP(L) = MAX(TEM,ZERO)
                  ELSEIF (TX5 > zero) THEN
                    QRP(L) = (MAX(ZERO,QA(1)/(ROR(L)*TX5*GMS(L))))      &
     &                                            ** (one/1.1364_kp)
                  ELSE
                    QRP(L) = zero
                  ENDIF
                else
                  qrp(l) = half * qrp(l)
                endif
!                                              Compute Buoyancy
                TEM1   = WA(3) + (HOD(L)-WA(1)-ALHL*(QOD(L)-WA(2)))     &
     &                         * onebcp

                TEM1   = TEM1 * (one + NU*QOD(L))
                ROR(L) = CMPOR * PRL(L) / TEM1
                TEM1   = TEM1 * DOFW
!!!             TEM1   = TEM1 * (1.0 + NU*QOD(L)) * DOFW

                BUY(L) = (TEM1 - one - QRP(L)) * ROR(L) * TX5
!                                              Compute W (L+1/2)

                TEM1   = WVL(L)
                WVL(L) = VT(2) * (ETD(L-1)*WVL(L-1) - FACG              &
     &                   * (BUY(L-1)*QRT(L-1)+BUY(L)*QRB(L-1)))

!
                if (wvl(l) < zero) then
!                 WVL(L) = max(wvl(l), 0.1*tem1)
!                 WVL(L) = 0.5*tem1
!                 WVL(L) = 0.1*tem1
!                 WVL(L) = 0.0
                  WVL(L) = 1.0e-10_kp
                else
                  WVL(L) = half*(WVL(L)+TEM1)
                endif

!
!               WVL(L) = max(0.5*(WVL(L)+TEM1), 0.0)

                ERRW   = WVL(L) - TEM1
!
                ERRQ   = ERRQ + ABS(ERRW/MAX(WVL(L),ONE_M5))


!               IF (ITR >= MIN(ITRMIN,ITRMD/2)) THEN
                IF (ITR >= MIN(ITRMND,ITRMD/2)) THEN

                  IF (ETD(L-1) == zero .AND. ERRQ > 0.2_kp) THEN

                    ROR(L) = BUD(KD)
                    ETD(L) = zero
                    WVL(L) = zero
                    ERRQ   = zero
                    HOD(L) = WA(1)
                    QOD(L) = WA(2)
!                   TX5      = TX1 + TX9
                    if (L <= KBL) then
                      TX5 = TX9
                    else
                      TX5 = (STLT(KB1) * QRT(KB1)                       &
     &                    +  STLT(KBL) * QRB(KB1)) * (0.5_kp*FAC)
                    endif

                    EVP(L-1) = zero
                    TEM      = MAX(TX1*RNT+RNF(L-1),ZERO)
                    QA(1)    = TEM - EVP(L-1)
!                   IF (QA(1) > 0.0) THEN

                    QRP(L)   = (QA(1) / (ROR(L)*TX5*GMS(L)))            &
     &                                              ** (one/1.1364_kp)
!                   endif
                    BUY(L)   = - ROR(L) * TX5 * QRP(L)
                    WCB(L-1) = zero
                  ENDIF
!
                  DEL_ETA = ETD(L) - ETD(L-1)
                  IF(DEL_ETA < zero .AND. ERRQ > 0.1_kp) THEN
                    ROR(L)   = BUD(KD)
                    ETD(L)   = zero
                    WVL(L)   = zero
!!!!!               TX5      = TX1 + TX9
                    CLDFRD(L-1) = TX5
!
                    DEL_ETA  = - ETD(L-1)
                    EDZ      = zero
                    DDZ      = -DEL_ETA
                    WCB(L-1) = DDZ
!
                    HOD(L)   = HOD(L-1)
                    QOD(L)   = QOD(L-1)
!
                    TEM4     = QOD(L)
                    TEM1     = VRW(1)
!
                    QHS      = QA(3) + half * (GAF(L-1)+GAF(L))         &
     &                                      * (HOD(L)-QA(2))

!
!                                           First iteration       !
!
                    ST2  = PRL(L) * (QHS + TEM1 * (QHS-QOD(L)))
                    TEM2 = ROR(L) * QRP(L-1)
                    CALL QRABF(TEM2,QRAF,QRBF)
                    TEM6 = TX5 * (1.6_kp + 124.9_kp * QRAF) * QRBF * TX4
!
                    CE   = TEM6*ST2/((5.4e5_kp*ST2 + 2.55e6_kp)*(ETD(L)+DDZ))
!

                    TEM2   = - ((one+TEM1)*(QHS+CE) + TEM1*QOD(L))
                    TEM3   = (one + TEM1) * QHS * (QOD(L)+CE)
                    TEM    = MAX(TEM2*TEM2 -FOUR*TEM1*TEM3,ZERO)
                    QOD(L) = MAX(TEM4, (- TEM2 - SQRT(TEM)) * VRW(2))
!
!                                            second iteration   !
!
                    ST2  = PRL(L) * (QHS + TEM1 * (QHS-QOD(L)))
                    CE   = TEM6*ST2/((5.4e5_kp*ST2 + 2.55e6_kp)*(ETD(L)+DDZ))
!                   CEE  = CE * (ETD(L)+DDZ)
!


                    TEM2   = - ((one+TEM1)*(QHS+CE) + TEM1*tem4)
                    TEM3   = (one + TEM1) * QHS * (tem4+CE)
                    TEM    = MAX(TEM2*TEM2 -FOUR*TEM1*TEM3,ZERO)
                    QOD(L) = MAX(TEM4, (- TEM2 - SQRT(TEM)) * VRW(2))

!                                              Evaporation in Layer L-1
!
                    EVP(L-1) = (QOD(L)-TEM4) * (ETD(L)+DDZ)

!                                               Calculate Pd (L+1/2)
!                   RNN(L-1) = TX1*RNT + RNF(L-1) - EVP(L-1)

                    QA(1)    = TX1*RNT + RNF(L-1)
                    EVP(L-1) = min(EVP(L-1), QA(1))
                    QA(1)    = QA(1) - EVP(L-1)
                    qrp(l)   = zero

!
!                   IF (QA(1) > 0.0) THEN
!!                    RNS(L-1) = QA(1)
!!!                   tx5      = tx9
!                     QRP(L) = (QA(1) / (ROR(L)*TX5*GMS(L)))              &
!    &                                           ** (1.0/1.1364)
!                   endif
!                   ERRQ   = 0.0
!                                              Compute Buoyancy
!                   TEM1   = WA(3)+(HOD(L)-WA(1)-ALHL*(QOD(L)-WA(2)))     &
!    &                                                    * (1.0/CP)
!                   TEM1   = TEM1 * (1.0 + NU*QOD(L)) * DOFW
!                   BUY(L) = (TEM1 - 1.0 - QRP(L)) * ROR(L) * TX5
!
!                   IF (QA(1) > 0.0) RNS(L) = QA(1)

                    IF (L .LE. K) THEN
                       RNS(L) = QA(1)
                       QA(1)  = zero
                    ENDIF
                    tx5    = tx9
                    ERRQ   = zero
                    QRP(L) = zero
                    BUY(L) = zero
!
                  ENDIF
                ENDIF
              ENDIF
!
            ENDDO                ! End of the iteration loop  for a given L!
            IF (L <= K) THEN
              IF (ETD(L-1) == zero .AND. ERRQ > 0.1_kp .and. l <= kbl) THEN
!!!  &           .AND. ERRQ > ERRMIN*10.0 .and. l <= kbl) THEN
!    &           .AND. ERRQ > ERRMIN*10.0) THEN
                 ROR(L)   = BUD(KD)
                 HOD(L)   = WA(1)
                 QOD(L)   = WA(2)
                 TX5      =       TX9     ! Does not make too much difference!
!                TX5      = TX1 + TX9
                 EVP(L-1) = zero
!                EVP(L-1) = CEE * (1.0 - qod(l)/qa(3))
                 QA(1)    = TX1*RNT + RNF(L-1)
                 EVP(L-1) = min(EVP(L-1), QA(1))
                 QA(1)    = QA(1) - EVP(L-1)

!                QRP(L)   = 0.0
!              if (tx5 == 0.0 .or. gms(l) == 0.0) then
!                write(0,*)' Ctx5=',tx5,' gms=',gms(l),' ror=',ror(l)   &
!    &,          ' L=',L,' QA=',QA(1),' tx1=',tx1,' tx9=',tx9           &
!    &,          ' kbl=',kbl,' etd1=',etd(l-1),' DEL_ETA=',DEL_ETA
!              endif
!                IF (QA(1) > 0.0) THEN

                   QRP(L) = (QA(1) / (ROR(L)*TX5*GMS(L)))               &
     &                                           ** (one/1.1364_kp)
!                ENDIF
                 ETD(L)   = zero
                 WVL(L)   = zero
                 ST1      = one - ALFIND(L)

                 ERRQ     = zero
                 BUY(L)   = - ROR(L) * TX5 * QRP(L)
                 WCB(L-1) = zero
              ENDIF
            ENDIF
!
            LL = MIN(IDN(idnm), KP1)
            IF (ERRQ < one .AND. L <= LL) THEN
              IF (ETD(L-1) > zero .AND. ETD(L) == zero) THEN
               IDN(idnm) = L
               wvl(l)    = zero
               if (L < KBL .or. tx5 > zero) idnm  = idnm + 1
               errq = zero
              ENDIF
              if (etd(l) == zero .and. l > kbl) then
                idn(idnm) = l
                if (tx5 > zero) idnm  = idnm + 1
              endif
            ENDIF
! 
!     If downdraft properties are not obtainable, (i.e.solution does
!      not converge) , no downdraft is assumed
!
!           IF (ERRQ > ERRMIN*100.0 .AND. IDN(idnm) == 99)                &
            IF (ERRQ > 0.1_kp .AND. IDN(idnm) == idnmax) DDFT = .FALSE.
!
            DOF = zero
            IF (.NOT. DDFT) RETURN
!
!         if (ddlgk .or. l .le. idn(idnm)) then
!           rsum2 = rsum2 + evp(l-1)
!           write(0,*)' rsum1=',rsum1,' rsum2=',rsum2,' L=',L,' qa=',qa(1)&
!    &,   ' evp=',evp(l-1)
!         else
!           rsum1 = rsum1 + rnf(l-1)
!           write(0,*)' rsum1=',rsum1,' rsum2=',rsum2,' L=',L,' rnf=',    &
!     &     rnf(l-1)
!         endif

          endif         ! if (l < idh)
        ENDDO                      ! End of the L Loop of downdraft !

        TX1 = zero

        DOF = QA(1)
!
!       write(0,*)' dof=',dof,' rntp=',rntp,' rnb=',rnb
!       write(0,*)' total=',(rsum1+dof+rntp+rnb)
!
        dof     = max(dof, zero)
        RNN(KD) = RNTP
        TX1     = EVP(KD)
        TX2     = RNTP + RNB + DOF

        II = IDH
        IF (II >= KD1+1) THEN
           RNN(KD)   = RNN(KD) + RNF(KD)
           TX2       = TX2 + RNF(KD)
           RNN(II-1) = zero
           TX1       = EVP(II-1)
        ENDIF
        DO L=KD,K
          II = IDH

          IF (L > KD1 .AND. L < II) THEN
            RNN(L-1) = RNF(L-1)
            TX2      = TX2 + RNN(L-1)
          ELSEIF (L >= II .AND. L < IDN(idnm)) THEN
            rnn(l)   = rns(l)
            tx2      = tx2 + rnn(l)
            TX1      = TX1 + EVP(L)
          ELSEIF (L >= IDN(idnm)) THEN
            ETD(L+1) = zero
            HOD(L+1) = zero
            QOD(L+1) = zero
            EVP(L)   = zero
            RNN(L)   = RNF(L) + RNS(L)
            TX2      = TX2    + RNN(L)
          ENDIF
        ENDDO
!
!      For Downdraft case the rain is that falls thru the bottom

        L = KBL

        RNN(L)    = RNN(L) + RNB
        CLDFRD(L) = TX5

!
!     Caution !! Below is an adjustment to rain flux to maintain
!                conservation of precip!

!
        IF (TX1 > zero) THEN
          TX1 = (TRAIN - TX2) / TX1
        ELSE
          TX1 = zero
        ENDIF

        DO L=KD,K
          EVP(L) = EVP(L) * TX1
        ENDDO

      ENDIF                       ! if (.not. DDFT) loop   endif
!
!***********************************************************************
!***********************************************************************

      RETURN
      end subroutine ddrft

      SUBROUTINE QSATCN(TT,P,Q,DQDT)
!
      USE FUNCPHYS , ONLY : fpvs

      implicit none
!
      real(kind=kind_phys) TT, P, Q, DQDT
!
!     real(kind=kind_phys), parameter :: ZERO=0.0, ONE=1.0              &
!    &,                                  rvi=one/rv,     facw=CVAP-CLIQ &
!    &,                                  faci=CVAP-CSOL, hsub=alhl+alhf &
!    &,                                  tmix=TTP-20.0                  &
!    &,                                  DEN=one/(TTP-TMIX)
!
      real(kind=kind_phys) es, d, hlorv, W
!
!     es = 10.0 * fpvs(tt)                ! fpvs is in centibars!
      es = min(p, 0.01_kp * fpvs(tt))     ! fpvs is in Pascals!
!     D  = one / max(p+epsm1*es,ONE_M10)
      D  = one / (p+epsm1*es)
!
      q  = MIN(eps*es*D, ONE)
!
      W     = max(ZERO, min(ONE, (TT - TMIX)*DEN))
      hlorv = ( W      * (alhl + FACW * (tt-ttp))                       &
     &       + (one-W) * (hsub + FACI * (tt-ttp)) ) * RVI
      dqdt  = p * q * hlorv *  D / (tt*tt)
!
      return
      end subroutine qsatcn

      SUBROUTINE ANGRAD(PRES, ALM,  AL2, TLA)
      implicit none

      real(kind=kind_phys) PRES, ALM,  AL2,  TLA,  TEM
!
      integer i
!
      IF (TLA < 0.0_kp) THEN
          IF (PRES <= PLAC(1)) THEN
            TLA = TLAC(1)
          ELSEIF (PRES <= PLAC(2)) THEN
            TLA = TLAC(2) + (PRES-PLAC(2))*tlbpl(1)
          ELSEIF (PRES <= PLAC(3)) THEN
            TLA = TLAC(3) + (PRES-PLAC(3))*tlbpl(2)
          ELSEIF (PRES <= PLAC(4)) THEN
            TLA = TLAC(4) + (PRES-PLAC(4))*tlbpl(3)
          ELSEIF (PRES <= PLAC(5)) THEN
            TLA = TLAC(5) + (PRES-PLAC(5))*tlbpl(4)
          ELSEIF (PRES <= PLAC(6)) THEN
            TLA = TLAC(6) + (PRES-PLAC(6))*tlbpl(5)
          ELSEIF (PRES <= PLAC(7)) THEN
            TLA = TLAC(7) + (PRES-PLAC(7))*tlbpl(6)
          ELSEIF (PRES <= PLAC(8)) THEN
            TLA = TLAC(8) + (PRES-PLAC(8))*tlbpl(7)
          ELSE
            TLA = TLAC(8)
          ENDIF
      ENDIF
        IF (PRES >= REFP(1)) THEN
          TEM = REFR(1)
        ELSEIF (PRES >= REFP(2)) THEN
          TEM = REFR(1) + (PRES-REFP(1)) * drdp(1)
        ELSEIF (PRES >= REFP(3)) THEN
          TEM = REFR(2) + (PRES-REFP(2)) * drdp(2)
        ELSEIF (PRES >= REFP(4)) THEN
          TEM = REFR(3) + (PRES-REFP(3)) * drdp(3)
        ELSEIF (PRES >= REFP(5)) THEN
          TEM = REFR(4) + (PRES-REFP(4)) * drdp(4)
        ELSEIF (PRES >= REFP(6)) THEN
          TEM = REFR(5) + (PRES-REFP(5)) * drdp(5)
        ELSE
          TEM = REFR(6)
        ENDIF
!
        tem = 2.0e-4_kp / tem
        al2 = min(4.0_kp*tem, max(alm, tem))
!
      RETURN
      end subroutine angrad

      SUBROUTINE SETQRP
      implicit none

      real(kind=kind_phys) tem2,tem1,x,xinc,xmax,xmin
      integer jx
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!     XMIN   = 1.0E-6
      XMIN   = 0.0_kp
      XMAX   = 5.0_kp
      XINC   = (XMAX-XMIN)/(NQRP-1)
      C2XQRP = one / XINC
      C1XQRP = one - XMIN*C2XQRP
      TEM1   = 0.001_kp ** 0.2046_kp
      TEM2   = 0.001_kp ** 0.525_kp
      DO JX=1,NQRP
        X         = XMIN + (JX-1)*XINC
        TBQRP(JX) =        X ** 0.1364_kp
        TBQRA(JX) = TEM1 * X ** 0.2046_kp
        TBQRB(JX) = TEM2 * X ** 0.525_kp
      ENDDO    
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      RETURN
      end subroutine setqrp

      SUBROUTINE QRABF(QRP,QRAF,QRBF)
      implicit none
!
      real(kind=kind_phys) QRP, QRAF, QRBF, XJ, REAL_NQRP
      INTEGER JX
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      REAL_NQRP = REAL(NQRP)
      XJ   = MIN(MAX(C1XQRP+C2XQRP*QRP,ONE),REAL_NQRP)
      JX   = MIN(XJ,NQRP-ONE)
      XJ   = XJ - JX
      QRAF = TBQRA(JX)  + XJ * (TBQRA(JX+1)-TBQRA(JX))
      QRBF = TBQRB(JX)  + XJ * (TBQRB(JX+1)-TBQRB(JX))
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      RETURN
      end subroutine qrabf

      SUBROUTINE SETVTP
      implicit none

      real(kind=kind_phys), parameter :: vtpexp=-0.3636_kp
      real(kind=kind_phys) xinc,x,xmax,xmin
      integer jx
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      XMIN   = 0.05_kp
      XMAX   = 1.5_kp
      XINC   = (XMAX-XMIN)/(NVTP-1)
      C2XVTP = one / XINC
      C1XVTP = one - XMIN*C2XVTP
      DO JX=1,NVTP
        X         = XMIN + (JX-1)*XINC
        TBVTP(JX) =        X ** VTPEXP
      ENDDO
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      RETURN
      end subroutine setvtp
!
      real(kind=kind_phys) FUNCTION QRPF(QRP)
!
      implicit none

      real(kind=kind_phys) QRP, XJ, REAL_NQRP
      INTEGER JX
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      REAL_NQRP = REAL(NQRP)
      XJ   = MIN(MAX(C1XQRP+C2XQRP*QRP,ONE),REAL_NQRP)
!     XJ   = MIN(MAX(C1XQRP+C2XQRP*QRP,ONE),FLOAT(NQRP))
      JX   = MIN(XJ,NQRP-ONE)
      QRPF = TBQRP(JX)  + (XJ-JX) * (TBQRP(JX+1)-TBQRP(JX))
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      RETURN
      end function qrpf

      real(kind=kind_phys) FUNCTION VTPF(ROR)
!
      implicit none
      real(kind=kind_phys) ROR, XJ, REAL_NVTP
      INTEGER JX
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      REAL_NVTP = REAL(NVTP)
      XJ   = MIN(MAX(C1XVTP+C2XVTP*ROR,ONE),REAL_NVTP)
      JX   = MIN(XJ,NVTP-ONE)
      VTPF = TBVTP(JX)  + (XJ-JX) * (TBVTP(JX+1)-TBVTP(JX))
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      RETURN
      end function vtpf

      real(kind=kind_phys) FUNCTION CLF(PRATE)
!
      implicit none
      real(kind=kind_phys) PRATE
!
      real (kind=kind_phys), parameter :: ccf1=0.30_kp, ccf2=0.09_kp    &
     &,                                   ccf3=0.04_kp, ccf4=0.01_kp    &
     &,                                   pr1=1.0_kp,   pr2=5.0_kp      &
     &,                                   pr3=20.0_kp
!
      if (prate < pr1) then
        clf = ccf1
      elseif (prate < pr2) then
        clf = ccf2
      elseif (prate < pr3) then
        clf = ccf3
      else
        clf = ccf4
      endif
!
      RETURN
      end function clf
      end module rascnv

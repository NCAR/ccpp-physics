!> \file canopy_driver.F

!> This file contains ...
      module canopy_driver
      use mfpbltq_mod
      use tridi_mod
      use mfscuq_mod
      !PCC_CANOPY_utilities
      use canopy_utils_mod
!IVAI
      use satmedmfvdifq_can_mod
      use canopy_mask_mod

! =====================
!     contains: canopy_mask_init, canopy_mask_run
! =====================

      use canopy_levs_mod

! =====================
!     contains: canopy_levs_init,  canopy_levs_run
!
! !Layers in reverse order!
!   1   is top    resolved layer
!   km  is bottom model hybrid layer
!   nkt is bottom canopy layer
!      ZMID_CAN3 (:,:,NLAYT)
!      ZMOM_CAN3 (:,:,NLAYT+1)
! =====================

      use canopy_transfer_mod

! =====================
!     contains: canopy_transfer_run
!  In:
!   Q1    (:,:, NLAYS, NSPCSD) :  Chemical tracers conc. ppmv on model levels
!
!  Output:
!   Q1_CAN3(:,:, NLAYT, NSPCSD) :  Chemical tracers conc. ppmv on combined canopy+resolved layers
!                                    ! CANOPY COLUMNS ONLY !
! ================================
!IVAI

      contains

!> \defgroup module_canopy_driver GFS TKE-EDMF PBL Module
!! This file contains the CCPP-compliant SATMEDMF scheme (updated version) which
!! computes subgrid vertical turbulence mixing using scale-aware TKE-based moist
!! eddy-diffusion mass-flux (TKE-EDMF) parameterization (by Jongil Han).
!> @{
!! \brief This subroutine contains all of the logic for the
!! scale-aware TKE-based moist eddy-diffusion mass-flux (TKE-EDMF, updated version) scheme.
!! For local turbulence mixing, a TKE closure model is used.
!! Updated version of satmedmfvdif.f to have better low level
!! inversion, to reduce the cold bias in lower troposphere,
!! and to reduce the negative wind speed bias in upper troposphere
!!
!! Incorporate the LES-based changes for TC simulation
!! (Chen et al.,2022 \cite Chen_2022)
!! with additional improvements on MF working with Cu schemes.
!!
!> \section arg_table_canopy_driver_init Argument Table
!! \htmlinclude canopy_driver_init.html
!!
      subroutine canopy_driver_init (satmedmf,                          &
     &                               isatmedmf,isatmedmf_vdifq,         &
     &                               errmsg,errflg)

      logical, intent(in   ) :: satmedmf
      integer, intent(in) :: isatmedmf,isatmedmf_vdifq

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

! Consistency checks
      if (.not. satmedmf) then
        write(errmsg,fmt='(*(a))') 'Logic error: satmedmf = .false.'
        errflg = 1
        return
      end if

      if (.not. isatmedmf==isatmedmf_vdifq) then
        write(errmsg,fmt='(*(a))') 'Logic error: satmedmfvdif is ',
     &                   'called, but isatmedmf/=isatmedmf_vdifq.'
        errflg = 1
        return
      end if

      end subroutine canopy_driver_init

!> \section arg_table_canopy_driver_run Argument Table
!! \htmlinclude canopy_driver_run.html
!!
!!\section gen_canopy_driver GFS canopy_driver General Algorithm
!! canopy_driver_run() computes subgrid vertical turbulence mixing
!! using the scale-aware TKE-based moist eddy-diffusion mass-flux (EDMF) parameterization of
!! Han and Bretherton (2019) \cite Han_2019 .
!! -# The local turbulent mixing is represented by an eddy-diffusivity scheme which 
!! is a function of a prognostic TKE.
!! -# For the convective boundary layer, nonlocal transport by large eddies 
!! (mfpbltq.f), is represented using a mass flux approach (Siebesma et al.(2007) \cite Siebesma_2007 ).
!! -# A mass-flux approach is also used to represent the stratocumulus-top-induced turbulence
!! (mfscuq.f).
!! \section detail_canopy GFS canopy_driver Detailed Algorithm 
      subroutine canopy_driver_run(im,km,ntrac,ntcw,ntrw,               &
     &     ntiw,ntke,ntchm,ntchs,ntche, ntoz,nto3,                      &
     &     ndtend,                                                      & !add ndtend
     &     con_rocp,                                                    & 
     &     grav,pi,rd,cp,rv,hvap,hfus,fv,eps,epsm1,                     &
!The following three variables are for SA-3D-TKE
     &     def_1,def_2,def_3,sa3dtke,dku3d_h,dku3d_e,                   &
!  dv=GFS_Interstitial(cdata%thrd_no)%dvdt
!  du=GFS_Interstitial(cdata%thrd_no)%dudt
! tdt=GFS_Interstitial(cdata%thrd_no)%dtdt
! rtg=GFS_Interstitial(cdata%thrd_no)%dvdftra
!  q1=GFS_Interstitial(cdata%thrd_no)%vdftra
!  u1=GFS_Statein%ugrs
!  v1=GFS_Statein%vgrs
!  t1=GFS_Statein%tgrs
!  q1=GFS_Interstitial(cdata%thrd_no)%vdftra
     &     dv,du,tdt,rtg,                                               & ! inout: dv,du,tdt,rtg
!  dv=GFS_Interstitial(cdata%thrd_no)%dvdt_can
!  du=GFS_Interstitial(cdata%thrd_no)%dudt_can
! tdt=GFS_Interstitial(cdata%thrd_no)%dtdt_can
! rtg=GFS_Interstitial(cdata%thrd_no)%dvdftra_can
     &     u1,v1,t1,q1,usfco,vsfco,use_oceanuv,                         & ! in
! swh=GFS_Radtend%htrsw
! hlw=GFS_Radtend%htrlw
     &     swh,hlw,xmu,garea,zvfun,sigmaf,                              & ! in
! u10m=GFS_Intdiag%u10m
! v10m=GFS_Intdiag%v10m
! t2m=GFS_Sfcprop%t2m
! q2m=GFS_Sfcprop%q2m
! fm=GFS_Sfcprop%ffmm
! fh=GFS_Sfcprop%ffhh
     &     psk,rbsoil,zorl,u10m,v10m,t2m, q2m,fm,fh,                    &
!  kpbl=GFS_Interstitial(cdata%thrd_no)%kpbl 2D
     &     tsea,heat,evap,stress,spd1,
     &     kpbl,                                                        & ! inout: kpbl
!  pgr=GFS_Statein%pgr ???
     &     pgr, 
!  del=GFS_Interstitial(cdata%thrd_no)%del
! prsi=GFS_Statein%prsi
! prsl=GFS_Statein%prsl
! prslk=GFS_Statein%prslk
! phii=GFS_Statein%phii
! phil=GFS_Statein%phil
     &     prsi,del,prsl,prslk,phii,phil,delt,tte_edmf,                 &
!dusfc=GFS_Interstitial(cdata%thrd_no)%dusfc1
!dvsfc=GFS_Interstitial(cdata%thrd_no)%dvsfc1
!dtsfc=GFS_Interstitial(cdata%thrd_no)%dtsfc1
!dqsfc=GFS_Interstitial(cdata%thrd_no)%dqsfc1
! hpbl=GFS_Tbd%hpbl
!  dkt=GFS_Intdiag%dkt
!  dku=GFS_Intdiag%dku
! tkeh=GFS_Interstitial(cdata%thrd_no)%tkeh
     &     dspheat,                                                     & ! in
     &     dusfc,dvsfc,dtsfc,dqsfc,hpbl,                                & ! in: dusfc,dvsfc,dtsfc,dqsfc,hpbl
     &     dkt,dku,tkeh,                                                & ! inout: dkt,dku, tkeh
     &     dkt_can,dku_can,                                             & ! out
!kinver=GFS_Interstitial 2d 
     &     kinver,xkzm_m,xkzm_h,xkzm_s,dspfac,bl_upfr,bl_dnfr,          & ! in 
     &     rlmx,elmx,sfc_rlm,tc_pbl,use_lpt,                            & ! in
!IVAI: canopy inputs from AQM                      
     &     do_canopy, cplaqm, claie, cfch, cfrt, cclu, cpopu,           & ! in
     &     ntqv,dtend,dtidx,index_of_temperature,index_of_x_wind,       & ! inout: dtend (.ldiag3d.)
     &     index_of_y_wind,index_of_process_pbl,gen_tend,ldiag3d,       &
     &     errmsg,errflg,                                               &
!IVAI: aux arrays
     &     naux2d,naux3d,aux2d,aux3d)

!
      use machine  , only : kind_phys
      use funcphys , only : fpvs
!
      implicit none
!
!----------------------------------------------------------------------
      integer, intent(in)  :: im, km,                                   &
     &                        ntrac, ntcw, ntrw, ntiw, ntke, ntqv,      &
     &                        ntchm,ntchs,ntche, ntoz,nto3, ndtend      & !IVAI 
      integer, intent(in)  :: sfc_rlm
      integer, intent(in)  :: tc_pbl
      integer, intent(in)  :: use_lpt
      integer, intent(in)  :: kinver(:)
      integer, intent(inout)  :: kpbl(:)
      logical, intent(in)  :: gen_tend,ldiag3d
!
      real(kind=kind_phys), intent(in) :: grav,pi,rd,cp,rv,hvap,hfus,fv,&
     &                                    eps,epsm1,                    &
     &                                    con_rocp                       !IVAI
      real(kind=kind_phys), intent(in) :: delt, xkzm_m, xkzm_h, xkzm_s
      real(kind=kind_phys), intent(in) :: dspfac, bl_upfr, bl_dnfr
      real(kind=kind_phys), intent(in) :: rlmx, elmx
!PCC CANOPY------------------------------------
      logical, intent(in)  :: do_canopy, cplaqm

!IVAI: canopy inputs
      real(kind=kind_phys), optional, intent(in)  ::
! 2D
     &       claie(:) , cfch(:),                                        &
     &       cfrt(:)  , cclu(:), cpopu(:)

      !----------------------------------------------
      real(kind=kind_phys), intent(inout) ::                            & 
     &                      dv(:,:),     du(:,:),                       &
     &                     tdt(:,:),   tkeh(:,:),
     &                     rtg(:,:,:)
      real(kind=kind_phys), intent(in) ::                               &
     &                     u1(:,:),     v1(:,:),                        &
     &                     usfco(:),    vsfco(:),                       &
     &                     t1(:,:),     q1(:,:,:),                      &
!The following two variables are for SA-3D-TKE
     &                     def_1(:,:),  def_2(:,:), def_3(:,:),         &
     &                     swh(:,:),    hlw(:,:),                       &
     &                     xmu(:),       garea(:),                      &
     &                     zvfun(:),     sigmaf(:),                     &
     &                     psk(:),       rbsoil(:),                     &
     &                     zorl(:),      tsea(:),                       &
     &                     u10m(:),      v10m(:),                       &
     &                     t2m(:),       q2m(:),                        & !IVAI
     &                     fm(:),        fh(:),                         &
     &                     evap(:),      heat(:),                       &
     &                     stress(:),    spd1(:),                       &
     &                     pgr(:),                                      & !IVAI: pgr=surface air pressure
     &                     prsi(:,:), del(:,:),                         &
     &                     prsl(:,:),   prslk(:,:),                     &
     &                     phii(:,:), phil(:,:)

      real(kind=kind_phys), intent(inout), dimension(:,:,:), optional ::&
     &     dtend

      integer, intent(in) :: dtidx(:,:), index_of_temperature,          &
     &       index_of_x_wind, index_of_y_wind, index_of_process_pbl
      logical, intent(in) :: use_oceanuv
      real(kind=kind_phys), intent(in) ::                               &
     &                     dusfc(:),     dvsfc(:),                      &
     &                     dtsfc(:),     dqsfc(:)                       &
      real(kind=kind_phys), intent(inout) ::
     &                     hpbl(:) ! use resolved hpbl in non-canopy columns
      real(kind=kind_phys), intent(inout) ::                            &
     &                     dkt(:,:), dku(:,:)
      real(kind=kind_phys), intent(out) ::                              &
     &                     dkt_can(:,:), dku_can(:,:)

      logical, intent(in)  :: sa3dtke !flag for SA-3D-TKE scheme
!
!
!  flag for tke dissipative heating
      logical, intent(in)  :: dspheat
!  flag for TTE-EDMF scheme
      logical, intent(in)  :: tte_edmf

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

!For passing dku to the dyn_core (SA-3D-TKE scheme)
       real(kind=kind_phys), intent(in) ::
     &            dku3d_h(:,:),dku3d_e(:,:)

!IVAI
      integer, intent(in) :: naux2d,naux3d
      real(kind_phys), intent(inout) :: aux2d(:,:)
      real(kind_phys), intent(inout) :: aux3d(:,:,:)
!IVAI

!          flag for tke dissipative heating
!
!----------------------------------------------------------------------
!***
!***  local variables
      real(kind=kind_phys) spd1_m
!***
      integer i,is,k,n,ndt,km1,kmpbl,kmscu,ntrac1,
     &                     idtend
!      integer kps,kbx,kmx
!      integer lcld(im),kcld(im),krad(im),mrad(im)
!
!      real(kind=kind_phys) tke(im,km), tkei(im,km-1), e2(im,0:km)
!
!      real(kind=kind_phys) theta(im,km),thvx(im,km),  thlvx(im,km),
!     &                     qlx(im,km),  thetae(im,km),thlx(im,km),
!     &                     slx(im,km),  svx(im,km),   qtx(im,km),
!     &                     tvx(im,km),  pix(im,km),   radx(im,km-1),
!     &                     dkq(im,km-1),cku(im,km-1), ckt(im,km-1)
!
!      real(kind=kind_phys) plyr(im,km), rhly(im,km),  cfly(im,km),
!     &                     qstl(im,km)
!
      real(kind=kind_phys) dtdz1(im), gdx(im),
     &                     phih(im),  phim(im),
     &                     phims(im), prn(im,km-1),
     &                     rbdn(im),  rbup(im),    thermal(im),
     &                     ustar(im), wstar(im),   
     &                     ust3(im),  wst3(im),     ! rho_a(im),
     &                     z0(im),    crb(im),     tkemean(im),
     &                     hgamt(im), hgamq(im),
     &                     wscale(im),vpert(im),
     &                     zol(im),   sflux(im),
     &                     sumx(im),  tx1(im),     tx2(im)
!
!      real(kind=kind_phys) radmin(im)
!
      real(kind=kind_phys) zi(im,km+1),  zl(im,km),   zm(im,km),
     &                     xkzo(im,km-1),xkzmo(im,km-1),
     &                     xkzm_hx(im),  xkzm_mx(im), tkmnz(im,km-1),
     &                     rdzt(im,km-1),rlmnz(im,km),
     &                     al(im,km-1),  ad(im,km),   au(im,km-1),
     &                     f1(im,km),    f2(im,km*(ntrac-1))
!
      real(kind=kind_phys) elm(im,km),   ele(im,km),
     &                     ckz(im,km),   chz(im,km),
     &                     diss(im,km-1),prod(im,km-1), 
     &                     bf(im,km-1),  shr2(im,km-1), 
     &                     xlamue(im,km-1), xlamde(im,km-1),
     &                     gotvx(im,km), rlam(im,km-1)
!
!   variables for updrafts (thermals)
!
!      real(kind=kind_phys) tcko(im,km),  qcko(im,km,ntrac),
!     &                     ucko(im,km),  vcko(im,km)
!
!   variables for stratocumulus-top induced downdrafts
!
!      real(kind=kind_phys) tcdo(im,km),  qcdo(im,km,ntrac),
!     &                     ucdo(im,km),  vcdo(im,km)
!
!  variables for Total Variation Diminishing (TVD) flux-limiter scheme
!
!      real(kind=kind_phys) e_half(im,km-1), e_diff(im,0:km-1),
!     &                     q_half(im,km-1,ntrac-1),
!     &                     qh(im,km-1,ntrac-1),
!     &                     q_diff(im,0:km-1,ntrac-1)
!      real(kind=kind_phys) rrkp, phkp
!      real(kind=kind_phys) tsumn(im), tsump(im), rtnp(im)
!      real(kind=kind_phys) sfcpbl(im), vez0fun(im)
!
!      logical  pblflg(im), sfcflg(im), flg(im)
!      logical  scuflg(im), pcnvflg(im)
!      logical  mlenflg
!
!  pcnvflg: true for unstable pbl
!
      real(kind=kind_phys) aphi16,  aphi5,
     &                     wfac,    cfac,
     &                     gamcrt,  gamcrq, sfcfrac,
!    &                     conq,    cont,   conw,
     &                     dsdz2,   dsdzt,  dkmax,
     &                     dsig,    dt2,    dtodsd,
     &                     dtodsu,  g,      factor, dz,
     &                     gocp,    gravi,  zol1,   zolcru,
     &                     concmin,                                !IVAI
     &                     buop,    shrp,   dtn,
     &                     prnum,   prmax,  prmin,  prtke,
     &                     prscu,   pr0,    ri,
     &                     dw2,     dw2min, zk,
     &                     elmfac,  elefac, dspmax,
     &                     alp,     clwt,   cql,
     &                     f0,      robn,   crbmin, crbmax,
     &                     es,      qs,     value,  onemrh,
     &                     cfh,     gamma,  elocp,  el2orc,
     &                     epsi,    beta,   chx,    cqx,
     &                     rdt,     rdz,    qmin,   qlmin,
     &                     rimin,   rbcr,   rbint,  tdzmin,
     &                     rlmn,    rlmn0,  rlmn1,  rlmn2,
     &                     ttend,   utend,  vtend,  qtend,
     &                     zfac,    zfmin,  vk,     spdk2,
     &                     tkmin,   tkbmx,  xkgdx,
     &                     xkinv1,  xkinv2,
     &                     zlup,    zldn,   cs0,    csmf,
     &                     tem,     tem1,   tem2,   tem3,
     &                     ptem,    ptem0,  ptem1,  ptem2
!
      real(kind=kind_phys) slfac
!
      real(kind=kind_phys) vegflo, vegfup, z0lo, z0up, vc0, zc0
!
      real(kind=kind_phys) ck0, ck1, ch0, ch1, ce0, rchck
!
      real(kind=kind_phys) qlcr, zstblmax, hcrinv
!
      real(kind=kind_phys) h1

      real(kind=kind_phys) bfac, mffac

      real(kind=kind_phys) q1_new(im,km,ntrac-1)

!IVAI
      integer kount 

!PCC_CANOPY------------------------------------
      integer COUNTCAN,KCAN
      ! Declare local maximum canopy layers
      integer, parameter :: MAXCAN = 1000

      real(kind=kind_phys) :: FCH,
     &                     ZI05, ZFL, DZFL, BOTCAN, FZI05, TTCORR,
     &               UUCAN, VVCAN, TTCAN, TKECAN, TKEHCAN,
     &               UU1, VV1, TT1, TKE1, TKEH1,
     &               UU_INT,VV_INT,TT_INT,TKE_INT, TKEH_INT,
     &               TT_SUM, ZZ_INT

      real(kind=kind_phys) :: ZCANX (MAXCAN), ZOOOX(MAXCAN),
     &                        UUX(MAXCAN),   VVX(MAXCAN), TTX(MAXCAN),
     &                       TKEX(MAXCAN), TKEHX(MAXCAN)
!PCC_CANOPY------------------------------------

      real(kind=kind_phys) ::
     &                      dv_can(im,km),  du_can(im,km),              & ! size (km)
     &                     duv_can(im,km),                              & ! size (km)
     &                     tdt_can(im,km), rtg_can(im,km,ntrac)           ! size (km)

      real(kind=kind_phys) :: 
     &                        dum3d_h (im, km)  , dum3d_e(im, km), 
     &                       tkeh_mod (im, km)
 
      real(kind=kind_phys) ::  dtend_can(im, km , ndtend), 
     &        aux3d_02(im,km), aux3d_04(im,km), aux3d_06(im,km)

! Out: list sat_canopy call
      real(kind=kind_phys) :: dusfc_can(im),  dvsfc_can(im),            &
     &                        dtsfc_can(im),  dqsfc_can(im),            &
     &                         hpbl_can(im), 
     &       aux2d_02(im),     aux2d_04(im), aux2d_06(im)

      integer              ::  kpbl_can(im)

      real(kind=kind_phys) :: qv_vmr(im,km), 
     &                   rho1(im,km), 
     &                   t2  (im,km), 
     &                   u2  (im,km), v2(im,km), ws2(im,km),
     &                   rho2(im,km), duv(im,km)
      real(kind=kind_phys) :: wind_dir_to_degrees,  
     &                   wind_dir_from_degrees, wind_dir_from_rad, 
     &                   wind_dir_cardinal, wind_dir_cardinal_rad

      integer
     &         kcan1, kc_can,
     &         kc, nkt1                ,
     &         kmod   (im, km)         ,
     &         kcan3  (im, nkc)        ,
     &         nfrct  (km+nkc,    im)     , !nkt
     &         ifrct  (km+nkc, 2, im)       !nkt

      real(kind=kind_phys) ::
! 2D arrays
     &       FRT_MASK  (im)            ,
! met2d arrays
     &         u2_mod (im, km)  ,  v2_mod(im, km),
     &        ws2_mod (im, km)  ,  t2_mod(im, km),
     &         du_mod (im, km)  ,  dv_mod(im, km),
     &        duv_mod (im, km)  , tdt_mod(im, km),
     &       U10M_CAN (im)      ,
     &       V10M_CAN (im)      ,
! all gas-phase species array
! NB. mfpbltq_mod: q1(ix,km,ntrac1) with  ntrac1 = ntrac - 1
     &       RTG1_MOD(im, km,     ntrac), !     before diffusion
     &       Q1_MOD  (im, km,     ntrac), !     before diffusion
     &       Q2_MOD  (im, km,     ntrac), !     after  diffusion
     &       Q2      (im, km,     ntrac), !     after  diffusion
     &       RTG_MOD (im, km,     ntrac), !     after  diffusion
! sat_can inputs
     &       Q1_CAN  (im, km    , ntrac), ! size (km)  before diffusion
     &       Q2_CAN  (im, km    , ntrac), ! size (km)  after  diffusion
     &       U1_CAN  (im, km )          , ! size (km)
     &       V1_CAN  (im, km )          , ! size (km)
     &      WS1_CAN  (im, km )          , ! size (km)
     &       T1_CAN  (im, km )          , ! size (km)
     &     TKEH_CAN  (im, km )          ,
!
     &        swh_can (im, km)        ,
     &        hlw_can (im, km)        ,
! sat_can inputs
     &       phii_can (im, km+1)      ,
     &       prsi_can (im, km+1)      ,
     &       prsl_can (im, km)        ,
     &        del_can (im, km)        ,
     &      prslk_can (im, km)        ,
     &       phil_can (im, km)        ,
! Canopy layers
     &      DKT_CAN3 (im, km+nkc)       ,
     &      DKU_CAN3 (im, km+nkc)       ,
     &       QV_CAN3 (im, km+nkc)       , ! nkt before diffusion
     &       Q1_CAN3 (im, km+nkc, ntrac), ! nkt before diffusion
     &       Q2_CAN3 (im, km+nkc, ntrac), ! nkt after  diffsion
! before diffusion
     &      WS1_CAN3 (im, km+nkc)        , ! using km for now only
     &       U1_CAN3 (im, km+nkc)        , ! using km for now only
     &       V1_CAN3 (im, km+nkc)        , ! using km for now only
     &       T1_CAN3 (im, km+nkc)        , ! using km for now only
     &     RHO1_CAN3 (im, km+nkc)        ,
     &     TKEH_CAN3 (im, km+nkc)        ,
! after diffusion
     &      WS2_CAN3 (im, km+nkc)        ,
     &       U2_CAN3 (im, km+nkc)        ,
     &       V2_CAN3 (im, km+nkc)        ,
     &       T2_CAN3 (im, km+nkc)        ,
     &     RHO2_CAN3 (im, km+nkc)        ,
! 2D
     &       Q1_2M   (im,      ntrac),    !     before diffusion
     &       Q2_2M   (im,      ntrac),    !     after  diffusion
! met3d arrays
     &       phii_can3 (im, km+nkc+1)      , !nkt
     &       phil_can3 (im, km+nkc)        ,
     &       prsi_can3 (im, km+nkc+1)      ,
     &       prsl_can3 (im, km+nkc)        ,
     &      prslk_can3 (im, km+nkc)        ,
     &        del_can3 (im, km+nkc)        ,
!
     &       zi_can3   (im, km+nkc+1)      , !nkt
     &       ZL_CAN3   (im, km+nkc)        , ! zl_can is ZH_CAN
     &       ZM_CAN3   (im, km+nkc)        , ! zm_can is ZF_CAN
     &       dz_can3   (im, km+nkc)        ,
!
! model layers
     &  wind_dir_to_rad(im, km)         ,
     &       ws1       (im, km)         ,
     &       wdir      (im, km)         ,

! layer height arrays !layers are in reverse order!
!     1   is top    resolved layer
!     km  is bottom model hybrid layer
!     km+nkc=nkt is bottom canopy layer
     &       zmom_can3 (im, km+nkc)        ,  ! dim(im, nkt+1)
     &       zmid_can3 (im, km+nkc)        ,
     &     sigmom_can3 (im, km+nkc)        ,  ! ~zm (nkt) or ~zi (nkt+1)
     &     sigmid_can3 (im, km+nkc)        ,  ! ~zl
     &    massair_can3 (im, km+nkc)        ,
     &    massair      (im, km)         ,
     &     mmr_o3_can3 (im, km+nkc)     ,
     &       frctr2c   (km+nkc, 2, im)     ,
     &       frctc2r   (km+nkc, 2, im)
!IVAI

!!
      parameter(bfac=100.)
      parameter(wfac=7.0,cfac=4.5)
      parameter(gamcrt=3.,gamcrq=0.,sfcfrac=0.1)
      parameter(vk=0.4,rimin=-100.,slfac=0.1)
      parameter(rbcr=0.25,zolcru=-0.02,tdzmin=1.e-3)
      parameter(rlmn=30.,rlmn0=5.,rlmn1=5.,rlmn2=10.)
      parameter(prmin=0.25,prmax=4.0)
      parameter(pr0=1.0,prtke=1.0,prscu=0.67)
      parameter(f0=1.e-4,crbmin=0.15,crbmax=0.35)
      parameter(tkmin=1.e-9,tkbmx=0.2,dspmax=10.0)
      parameter(qmin=1.e-8,qlmin=1.e-12,zfmin=1.e-8)
!IVAI:rbdriver.F CONCMIN   = 1.0E-30
      PARAMETER(concmin = 1.0E-30)        ! Minimum conc
!IVAI
      parameter(aphi5=5.,aphi16=16.)
      parameter(elmfac=1.0,elefac=1.0,cql=100.)
      parameter(dw2min=1.e-4,dkmax=1000.,xkgdx=1000.)
      parameter(qlcr=3.5e-5,zstblmax=2500.)
      parameter(xkinv1=0.15,xkinv2=0.3)
      parameter(h1=0.33333333,hcrinv=250.)
      parameter(vegflo=0.1,vegfup=1.0,z0lo=0.1,z0up=1.0)
      parameter(vc0=1.0,zc0=1.0)
      parameter(ck1=0.15,ch1=0.15)
      parameter(cs0=0.4,csmf=0.5)
      parameter(rchck=1.5,ndt=20)

      if (tc_pbl == 0) then
        ck0 = 0.4
        ch0 = 0.4
        ce0 = 0.4
      else if (tc_pbl == 1) then
        ck0 = 0.55
        ch0 = 0.55
        ce0 = 0.12
      endif
      gravi = 1.0 / grav
      g = grav
      gocp = g / cp
!     cont=cp/g
!     conq=hvap/g
!     conw=1.0/g  ! for del in pa
!!    parameter(cont=1000.*cp/g,conq=1000.*hvap/g,conw=1000./g) !kpa
      elocp  = hvap / cp
      el2orc = hvap * hvap / (rv * cp)

!
!************************************************************************
! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

!> ## Compute preliminary variables from input arguments
      dt2  = delt
      rdt = 1. / dt2
!
! the code is written assuming ntke=ntrac
! if ntrac > ntke, the code needs to be modified
!
      ntrac1 = ntrac - 1
      km1 = km - 1
      kmpbl = km / 2
      kmscu = km / 2
!>  - Compute physical height of the layer centers and interfaces from 
!! the geopotential height (\p zi and \p zl)
      do k=1,km
        do i=1,im
          zi(i,k) = phii(i,k) * gravi
          zl(i,k) = phil(i,k) * gravi
          ckz(i,k) = ck1
          chz(i,k) = ch1
          rlmnz(i,k) = rlmn0
        enddo
      enddo
      do i=1,im
        zi(i,km+1) = phii(i,km+1) * gravi
      enddo
      do k=1,km
        do i=1,im
          zm(i,k) = zi(i,k+1)
        enddo
      enddo
!>  - Compute horizontal grid size (\p gdx)
      do i=1,im
        gdx(i) = sqrt(garea(i))
      enddo

!>  - Some output variables and logical flags are initialized
!      do i = 1,im
!         z0(i)    = 0.01 * zorl(i)
!         rho_a(i) = prsl(i,1)/(rd*t1(i,1)*(1.+fv*max(q1(i,1,1),qmin)))

!         pblflg(i)= .true.
!         sfcflg(i)= .true.
!         if(rbsoil(i) > 0.) sfcflg(i) = .false.
!         pcnvflg(i)= .false.
!         scuflg(i)= .true.
!         if(scuflg(i)) then
!           radmin(i)= 0.
!           mrad(i)  = km1
!           krad(i)  = 1
!           lcld(i)  = km1
!           kcld(i)  = km1
!         endif
!      enddo

!IVAI
! Initialize canopy layers concentrtions with values before diffusion
      if (do_canopy .and. cplaqm) then

! TKE
        aux3d(:,:, 7) = q1(:,:, ntke) ! ntke=198 "tke"

        do k = 1,km
          do i = 1,im
            rho1(i,k) = prsl(i,k)/
     &                (rd*t1(i,k)*
     &         (1.+fv*max(q1(i,k, 1),qmin)))

          enddo
        enddo

! convert mass m.r. to volume m.r.
!        qv_vmr(:,:) = q1(:,:,1)   ! ntqv=1 water vapor 
! q2m in ppmv ???
!       q2m_vmr(:) = q2m(:)


        CALL canopy_mask_init( im, im, km, ! nkt, 
     &     claie, cfch, cfrt, cclu, cpopu, !in:
     &     FRT_MASK,                      !out: nkt defined here!
     &     errmsg, errflg)

        if (errflg /= 0) return

!      nkt=km+nkc defined in "canopy_mask_init"
      nkt1 = nkt - 1

        CALL canopy_levs_init( im, im, km, nkc, nkt,
     &       ntrac-1, ntqv, ntke,                 ! in
     &       errmsg, errflg,                      &
     &       zi, zl, zm, !in: 3D
     &       prsl, prsi, !in: 3D
     &       dv, du, tdt, rtg,                    ! in: 3D
     &       u1, v1, t1, q1,                      ! in 3D / 4D
     &       rho1, dkt, dku,                      ! in 3D
     &       dtend, 
     &     zmom_can3,   zmid_can3,                !out 3D
     &   sigmom_can3, sigmid_can3,                !out 3D
     &       ZL_CAN3,     ZM_CAN3,                !out 3D ZH_CAN, ZF_CAN rename half- and full-layer height
     &     PRSL_CAN3,   PRSI_CAN3,                !out 3D set to zero
     &       dv_can, du_can, tdt_can, rtg_can,    !out: 3D size (km)
     &       T1_CAN3,     QV_CAN3,rho1_CAN3,      !out 3D set to zero
     &      WS1_CAN3,    DKT_CAN3, DKU_CAN3,      !out 3D set to WS DKT
     &       Q1_CAN3, Q1_2M,                      !out 4D set to Q1
     &       DTEND_CAN )

        if (errflg /= 0) return

! ================
! In;  Q1 (im,km,NTRAC-1)
! Out:
!      Q1_CAN3(:,1:3,NTRAC-1) <= Q1(:,1,NTRAC-1) ! ALL CANOPY & NON-CANOPY COLUMNS !
! ===============

        CALL canopy_transfer_init(im, im, km,        !in
     &        massair_can3, massair,                  !out
     &        mmr_o3_can3,                            !inout
     &        nfrct, ifrct,                          !out
     &        frctr2c, frctc2r,                      !out
     &        errmsg, errflg )

        if (errflg /= 0) return


      endif ! (do_canopy .and. cplaqm)
!IVAI


      !PCC_CANOPY------------------------------------
      kount=0 !IVAI
      if (do_canopy .and. cplaqm) then 

! NB. Call canopy routines after eddy diffusivities are calculated!!!

        CALL canopy_mask_run( im, im, km,           !in
     &     claie, cfch, cfrt, cclu, cpopu,           !in
     &     FRT_MASK,                                !out
     &     errmsg, errflg)

        if (errflg /= 0) return

! Temporary 2D output
!         aux2d(:, 1) = FRT_mask(:)
        
! Wind direction, degrees
! ATAN2(Y, X) computes the principal value of the argument function of the complex number X + i Y.
!  This function can be used to transform from Cartesian into polar coordinates and allows to determine the angle in the correct quadrant.

         do k=1,km
           do i=1,im

             ws1(i,k) = sqrt(u1(i,k)**2+v1(i,k)**2)

             wind_dir_to_rad(i,k) =
!    &                       atan2(v1(i,k)/ws1(i,k),u1(i,k)/ws1(i,k))   ! to radians
     &                       atan2(u1(i,k)/ws1(i,k),v1(i,k)/ws1(i,k))   ! to radians

             wind_dir_to_degrees = wind_dir_to_rad(i,k) * 180./pi       ! to degrees
             wind_dir_from_degrees = wind_dir_to_degrees + 180.         ! from degrees
             wind_dir_cardinal = 90. - wind_dir_from_degrees            ! cardinal degrees

             wind_dir_from_rad = wind_dir_from_degrees * pi/180.        ! from radians
             wind_dir_cardinal_rad = wind_dir_cardinal * pi/180.        ! cardinal radians
           enddo
         enddo

! ===============
! In canopy_levs_run, the vertical layers are going down (reversed to phot)
!
! k = 1, NLAYT from  top to bottom of augmented canopy+resolved model layers
!      1 = top model layer
!  NLAYS = bottom resolved model layer
!  NLAYC = 3 canopy layers
!  NLAYT = NLAYS + NLAYC, augmented canopy plus resolved model layers
! NLAYS+1= top    canopy layer at 1.0*CH
! NLAYS+2= middle canopy layer at 0.5*CH
! NLAYS+3= bottom canopy layer at 0.2*CH
! =================

         CALL canopy_levs_run(im, im, km, nkc, nkt, ! in
     &       ntrac1, ntqv, ntke,                   ! in
     &       errmsg, errflg,                       ! out
     &       RD, PI,                               ! in gry gas constant
     &       zi, zl, zm,                           ! in
     &       prsl, prsi, pgr,                      ! in  (Pa)
     &       cfch,                                 ! in: canopy data input
     &       garea, u10m,v10m, fm, fh,             ! in: 2D
     &       rbsoil,                               ! in: 2D
     &       t2m, q2m,                             ! in 2D
     &       stress, spd1,                         ! in: 2D
     &       dv, du, tdt, rtg,                     ! in: 3D
     &       u1, v1, t1, q1,                       ! in: 3D " 4D
     &       rho1, dkt, dku,                       ! in 3D
     &       FRT_MASK,                             ! in 2D canopy_mask
     &       kmod, kcan3,                          ! out
     &      zmom_can3,   zmid_can3,                ! out
     &    sigmom_can3, sigmid_can3,                ! out
     &        ZL_CAN3,     ZM_CAN3,                ! out: zl=ZH_CAN, zm=ZF_CAN rename half- and full-layer height
     &      PRSL_CAN3,   PRSI_CAN3,                ! out: mean layer pressure; air pressure at model layer interfaces
     &       dv_can, du_can, tdt_can, rtg_can,     ! out: 3D size (km)
     &        T1_CAN3,     QV_CAN3, rho1_CAN3,     ! out 3D: 2-m interpolated T1 QV rho1
     &       WS1_CAN3,    DKT_CAN3,  DKU_CAN3,     ! out 3D: 10-m interpolated WS1
     &        Q1_CAN3,  Q1_2M)                     ! inout kg kg-1

         if (errflg /= 0) return

! ================
! Out: 
!      T1_CAN3   (:,:,NLAYT)
!      QV_CAN3   (:,:,NLAYT) Q2m interpolated
!      PRSI_CAN3 (   NLAYT+1)
!      PRSL_CAN3 (:, NLAYT)
!      rho1_CAN3 (:, NLAYT)
!      Q1_CAN3(1,2,3) <= Q1(1)   ! ALL CANOPY & NON-CANOPY COLUMNS !

! !Layers in reverse order!
!     1   is top    resolved layer
!     km  is bottom model hybrid layer
!     nkt is bottom canopy layer
!  zmid_can3   (:,:, NLAYT)   layers are in reverse order!
!  zmom_can3   (:,:, NLAYT+1) layers are in reverse order!
!  massair_can3(:,:, NLAYT)    :  mass of air in canopy layers (kg)
!  massair     (:,:, NLAYS)    :  mass of air in model layers (kg)
!                                        (gathered canopy + resolved scale columns)
!   nfrct  (NLAYT,  :,:)        :  Number of original model levels contributing to canopy level k
!   ifrct  (NLAYT,2,:,:)        :  Index of the original model level contributing to canopy level k
!   frctr2c(NLAYT,2,:,:)        :  Fractional contribution of the original model level to canopy level k
!   frctc2r(NLAYT,2,:,:)        :  Fractional contribution of the canopy level to the original model level
! ================

         do i = 1,im
!           rho1_can_sfc(i) = prsl_can3(i,1)/
!     &               (rd*t1_can3(i,1)*(1.+fv*max(q1_can3(i,1, ntqv),qmin))) ! ntqv=1
         enddo 

        do k = 1,km ! ntk to do
          do i = 1,im

            q1_can3(i,k, ntqv) = qv_can3(i,k) ! 2-m interpolated humidity

            rho1_can3(i,k) = prsl_can3(i,k)/
     &                     (rd*t1_can3(i,k)*
     &              (1.+fv*max(q1_can3(i,k, ntqv),qmin))) ! ntqv=1 2-m interpolated

          enddo
        enddo

! Humidity on canopy layers, interpolated from q2m qv_can3
! Humidity overwritten in "resolved_to_canopy" mass transfer
!         q1_can3(:,1:km, ntqv ) = qv_can3(:,1:km)

! Output 2D pbl diags
!        aux2d(:, 1) =  rho_a (:)      ! GOOD
!        aux2d(:, 1) =  q1_can3(:,1, 1) ! GOOD
!        aux2d(:, 2) =  rho_a_can(:)   ! GOOD

!        aux2d(:, 5) =  rho1    (:,1) GOOD
!        aux2d(:, 6) =  rho1_can3(:,1) GOOD

! Output 3D pbl diags
!        aux3d(:,:,5) = rho1    (:,1:km) 
!        aux3d(:,:,6) = rho1_can3(:,1:km) !  "2-m interpolated" air density

!        aux3d(:,:,3) = t1    (:,1:km) ! save in sat routine
!        aux3d(:,:,4) = t1_can3(:,1:km) ! "2-m interpolated" temperature

!        aux3d(:,:,2) =   qv_can3(:,1:km) ! "2-m interpolated" humidity 

!        aux3d(:,:,2) = prsl_can3(:,1:km) 

!        aux3d(:,:,1) = rho1     (:,1:km)
!        aux3d(:,:,2) = rho1_can3(:,1:km)

!        aux3d(:,:,6) =   zh_can(:,1:km)
!        aux3d(:,:,6) =   zf_can(:,1:km)

!        aux3d(:,:,5) = ws1    (:,1:km) ! save in sat routine
!        aux3d(:,:,6) = WS1_CAN3(:,1:km) ! 10-m interpolated 

!        aux3d(:,:, 4) = dkt_can3 (:,1:km)

! del    (:,k) = prsi    (:,k) - prsi    (:,k+1)
! del_can3(:,k) = prsi_can3(:,k) - prsi_can3(:,k+1)
! Above canopy layers
         do k=1,km 
           do i=1,im
             ! kc = 4,5,6.. 67
             kc = nkc + k
             del_can3  (i,kc) =  prsi_can3(i,kc) - prsi_can3(i,kc+1) 

! Exner function canopy layers
! !< exner function = (p/p0)**rocp
             prslk_can3(i,kc) = (prsl_can3(i,kc) /pgr(i)) ** con_rocp

           enddo
         enddo

! Canopy layers: kc = 1, 2, 3
         do kc = 1, nkc
           do i = 1,im

! Canopy columns
             IF (FRT_MASK(i) >  0.) THEN
               del_can3(i,kc) = prsi_can3(i, kc) - prsi_can3(i, kc+1)

! Non-canopy columns set to del
            ELSE IF (FRT_MASK(i) <=  0.) THEN
              del_can3(i,kc) = del(i,1) 

            ENDIF

! Exner function canopy layers
! !< exner function = (p/p0)**rocp
             prslk_can3(i,kc) = (prsl_can3(i,kc) /pgr(i)) ** con_rocp

           end do
         end do

!--------
! zmid(1)  = zl(km) is top    height model  layer centers
! zmid(km) = zl(1)  is bottom height model  layer centers
! NB.
! ZH_CAN  (i, nkt + 1 - k) = zmid_can(i, k) (k = 1, nkt)  combined layer centers => rename zl_can
! zi height model layer interfaces
!  dz_can = zl_can3(i,k+1) - zl_can3(i,k)
! so zm(i,k) = zi(i,k+1) = zl_can3(i,k) + dz_can3(i,k)/2
!      zm (:,k) = zi(:,k+1), so zm_can (i, k) = zi_can3(i,k+1)  (k=1,km) 

! Above canopy layers
         do k = 1,km
           do i = 1,im
             ! kc = 4,5,6.. 67
             kc = nkc + k
! dim zi_can3(im, nkt+1)
             zi_can3(i,kc+1) = zm_can3(i, kc)                    ! upper  interface
             dz_can3(i,kc)   = zi_can3(i, kc+1) - zi_can3(i, kc)
           end do
         end do

! Canopy layers: kc = 1, 2, 3
         zi_can3(:, 1) = 0.
         do kc = 1, nkc
           do i = 1,im
             zi_can3(i,kc+1) = zm_can3(i, kc)
             dz_can3(i,kc)   = zi_can3(i, kc+1) - zi_can3(i, kc)
           end do
         end do

!>  - Compute geopotential physical height of the layer centers and interfaces from
!! the physical height (\p zi and \p zl)
         do k=1,nkt
           do i=1,im
             phil_can3(i,k) = zl_can3(i,k) * grav
             phii_can3(i,k) = zi_can3(i,k) * grav
           enddo
         enddo
         do i=1,im
           phii_can3(i,nkt+1) = zi_can3(i,nkt+1) * grav
         enddo

! Print

! First, ...
         do k = 1,km

            ! nkt is top
            ! nkc+1 is bottom
            kc=nkc+k

            do i=1,im

              u1_can3(i,kc) = ws1_can3(i,kc) * sin(wind_dir_to_rad(i,k))      ! m/s
              v1_can3(i,kc) = ws1_can3(i,kc) * cos(wind_dir_to_rad(i,k))      ! m/s

            end do
         end do

! Canopy Layers: use 1hy resolved model layer wind direction
         do kc = 1,nkc
            do i=1,im

             u1_can3(i,kc) = ws1_can3(i,kc) * sin(wind_dir_to_rad(i,1))    ! m/s
             v1_can3(i,kc) = ws1_can3(i,kc) * cos(wind_dir_to_rad(i,1))    ! m/s

            end do
         end do

! All columns & all layers
         q1_mod (:,1:km, :) = q1   (:,1:km, :)         ! before "resolved_to_canopy"
!        Q1_2M  ... move out of "canopy_levs"

!===============================================================================
! Distribute tracer concentration from model resolved layers into canopy layers
! flag = 0 "resolved_to_canopy"
!===============================================================================

         CALL canopy_transfer_run(im, im, km,        !in
     &        ntrac1, ntoz,                          !in
     &        garea,                                 !in
     &        zi, zl, zm,                            !in
     &        q1, rho1,                              !in kg kg-1
     &        0,                                     !in 0 = "resolved_to_canopy"
     &        FRT_MASK,                              !in
     &        kmod, kcan3,                           !in
     &        zmom_can3, zmid_can3,                  ! in
     &        PRSL_CAN3, rho1_CAN3,                  !in:    before diffusion
     &        Q1_MOD, Q1_CAN3, Q1_2M,                !inout: kg kg-1 before diffusion
     &     massair_can3, massair,                  !inout
     &      mmr_o3_can3,                            !inout
     &        nfrct, ifrct,                          !inout
     &        frctr2c, frctc2r,                      !inout
     &        errmsg, errflg )

         if (errflg /= 0) return

! ==============
!  Input:
!   Q1    (:,:, NLAYS, ntrac1) :  Chemical tracers conc. ppmv on model levels
!
!  Output:
!   Q1_CAN3(:,:, NLAYT, ntrac1) :  Chemical tracers conc. ppmv on combined canopy+resolved layers
!                                    ! CANOPY COLUMNS ONLY !
!   Q1_2M (:,:       , ntrac1) :  2M Chemical tracers conc. ppmv Diagnostics
!
! ================================

! Comment out to use 2-m interpolated value from "canopy_levs" qv_can3(:,:)
!          q1_can3 (:,1:km, ntqv) = qv_can3(:,1:km) ! ntqv=1

! Subset the canopy tracers/arrays for input to "sat_can", since  routine is coded on dim(km) 
         prsi_can  (:,1:km+1) =  prsi_can3 (:,1:km+1)
         prsl_can  (:,1:km)   =  prsl_can3 (:,1:km)
         prslk_can (:,1:km)   = prslk_can3 (:,1:km) 
         del_can   (:,1:km)   =   del_can3 (:,1:km)

         phii_can  (:,1:km+1) = phii_can3(:,1:km+1)
         phil_can  (:,1:km)   = phil_can3(:,1:km)

! NB. Using 10-m interpolated values creates shear and gives very high TKE tendencies

         t1_can  (:,1:km)     = t1_can3  (:,1:km)

         dku_can (:,1:km)     = dku_can3 (:,1:km) ! "canopy_levs"
         dkt_can (:,1:km)     = dkt_can3 (:,1:km) ! "canopy_levs"

! Prepare arrays for "Resolved_to_Canopy" transfer
! Mass tracers (ntrac1) except TKE
         q1_can (:,1:km, 1:ntrac1) = q1_can3(:,1:km, 1:ntrac1) ! ntrac1 "resolved_to_canopy"
! TKE tracer
         q1_can (:,1:km,   ntke  ) = q1_can3 (:,1:km,   ntke  ) ! ntke   "resolved_to_canopy"

      endif !do_canopy .and. cplaqm

      if (do_canopy .and. cplaqm) then

! Save a copy of dtend for the canopy call , before adding vdiff tendecies on model layers
! dtend(im,km,ndtend)

! 3D array on combined canopy plus resolved model layers

! Test with km combined canopy plus resolved layers, so skip the top combined 3 layers
! This should be nkt layers...       dv_can(:,nkc+1:nkt) = dv(:,1:km) ! nkt combined canopy plus resolved layers

! Sub-Canopy
         swh_can (:, nkc+1:km ) = swh(:,1:km)
         swh_can (:, 3 )        = swh(:,1   )
         swh_can (:, 2 )        = swh(:,1   )
         swh_can (:, 1 )        = swh(:,1   )

         hlw_can (:, nkc+1:km ) = hlw(:,1:km)
         hlw_can (:, 3 )        = hlw(:,1   )
         hlw_can (:, 2 )        = hlw(:,1   )
         hlw_can (:, 1 )        = hlw(:,1   )

! Output pbl diags
!        aux3d(:,:, 2) = q1_can3(:,1:km, ntke)

! Subset combined layers (minus top nkc layers)
         do k = km-nkc, 1, -1 ! top to 1hy model layer
           ! km     is top combined subset
           ! nkc+1  is bot combined
           kc= nkc+k     ! 4th from top (nkt) to nkc+1 combined canopy plus resolved model layer

! Above-canopy TKE tracer set to hybrid model layers ("resoved_to_canopy" only does mass tranfer of mass tracers)
           q1_can (:,kc, ntke    ) = q1(:,k, ntke    ) ! ntke always on

! Above-canopy wind components set to hybrid model layers
           u1_can (:,kc) =  u1 (:,k)
           v1_can (:,kc) =  v1 (:,k)

           ws1_can(:,kc) = sqrt(u1_can(:,kc)**2+v1_can(:,kc)**2)

         end do

         do kc = 1, nkc ! 3-nkc canopy layers

! Sub-canopy values of TKE  set to 1hy model layer
! ("canopy_transfer" only does mass transfer to mass conc. tracers)
            q1_can (:,kc, ntke    ) = q1(:,1, ntke    ) ! ntke always on

!Sub-canopy values of wind components  set to 1hy model layer
            u1_can (:,kc) =  u1 (:,1)
            v1_can (:,kc) =  v1 (:,1)

            ws1_can(:,kc) = sqrt(u1_can(:,kc)**2+v1_can(:,kc)**2)
         end do

!!! BEFORE SAT CANOPY CALL!!

! Output canopy pbl tendency of QV
         if(ldiag3d) then

! Output pbl diffusivities
!          aux3d(:,:, 5) = dku (:,1:km) ! InOut  GOOD
!          aux3d(:,:, 3) = dkt (:,1:km) ! InOut  GOOD

! GFSv17_p8:
!  rtg_no2_index   = 10 ! "mp_thompson"
!  rtg_no_index    = 11 ! "mp_thompson"
!  rtg_o3_index    = 12 ! "mp_thompson"
!  rtg_no3_index   = 13 ! "mp_thompson"
!          aux3d(:,:, 1) =  dv (:,:)       ! before "canopy_to_resolved"
!          aux3d(:,:, 1) =  du (:,:)       ! before "canopy_to_resolved"
! duv is below after u2 & v2

! Output pbl diags 2D
           aux2d(:, 3) = float(kpbl(:)) ! before canopy
           aux2d(:, 1) = hpbl(:)        ! before canopy

         endif

!
!>  - Call satmedmfvdifq_can(), which is ...
!!    to take into account ...
         CALL satmedmfvdifq_can(im,km,ntrac,ntcw,ntrw,                  &
     &     ntiw,ntke,grav,pi,rd,cp,rv,hvap,hfus,fv,eps,epsm1,           &
!The following three variables are for SA-3D-TKE
     &     def_1,def_2,def_3,sa3dtke,dum3d_h,dum3d_e,                   &
     &     dv_can,du_can,tdt_can,rtg_can,                               & ! InOut
     &     u1_can,v1_can,t1_can, q1_can,                                & ! In: canopy inputs
     &     usfco,vsfco,use_oceanuv,                                     &
     &     swh_can,hlw_can,                                             & ! In: canopy inputs
     &     xmu,garea,zvfun,sigmaf,                                      &
     &     psk,rbsoil,zorl,u10m,v10m,fm,fh,                             &
     &     tsea,heat,evap,stress,spd1,
     &     kpbl_can,                                                    & ! Out
     &     prsi_can,del_can,prsl_can,prslk_can,phii_can,phil_can,       & ! In: canopy inputs
     &     delt,tte_edmf,                                               &
     &     dspheat,
     &     dusfc_can,dvsfc_can,dtsfc_can,dqsfc_can,hpbl_can,            & ! Out
     &     dkt,     dku,     tkeh_can,                                  & ! Out/Out:tkeh_can
     &     dkt_can,dku_can,                                             & ! In: canopy inputs
     &     kinver,xkzm_m,xkzm_h,xkzm_s,dspfac,bl_upfr,bl_dnfr,          &
     &     rlmx,elmx,sfc_rlm,tc_pbl,use_lpt,                            &
!IVAI: canopy inputs from AQM
     &     do_canopy, cplaqm, claie, cfch, cfrt, cclu, cpopu,           &
!IVAI
     &     ntqv,
     &     dtend_can,                                                   & !inout: dtend (.ldiag3d.)
     &     dtidx,index_of_temperature,index_of_x_wind,                  &
     &     index_of_y_wind,index_of_process_pbl,gen_tend,ldiag3d,       &
     &     errmsg,errflg, 
!IVAI: aux arrays
     &     naux2d,naux3d,aux2d,aux3d)

! Set non-canopy columns to resolved values
! NB. Only vars not ALREADY defined in non-canopy columns
         do i = 1,im
         IF (FRT_MASK(i) <=  0.) THEN

! Non-canopy columns 
!           kpbl_can(i) = kpbl(i) ! kpbl zero before original sat call 
            hpbl_can(i) = hpbl(i)

! Non-canopy columns
            dusfc_can(i) = dusfc(i) ! dusfc is zero before the main sat call
            dvsfc_can(i) = dvsfc(i) ! dvsfc is zero before the main sat call
            dtsfc_can(i) = dtsfc(i) ! dtsfc is zero before the main sat call
            dqsfc_can(i) = dqsfc(i) ! dqsfc is zero before the main sat call

         END IF !(FRT_MASK)
         end do

        aux2d(:, 4) = float(kpbl_can(:))  ! after canopy
        aux2d(:, 2) = hpbl_can (:)        ! after canopy

! Output 3D pbl diags
!       aux3d(:,:, 6) = dku_can  (:,1:km) !
!       aux3d(:,:, 4) = dkt_can  (:,1:km) !

c
!> - Apply the tendencies of heat and moisture on canopy layers
! NB. before doing "canopy_to_resolved" mass transfer
c

!U-wind/V-wind on original model layers after diffusion
          u2    (:,1:km)    = u1    (:,1:km) +
     &                        du    (:,1:km) * dt2      ! before "canopy_to_resolved"
          v2    (:,1:km)    = v1    (:,1:km) +
     &                        dv    (:,1:km) * dt2      ! before "canopy_to_resolved"

          ws2   (:,1:km)    = sqrt(u2(:,1:km)**2+v2(:,1:km)**2)
          duv   (:,1:km)    =    (ws2(:,1:km) - ws1(:,1:km) )*rdt  ! before "canopy_to_resolved"

! Air temperature on original model layers after diffusion
          t2    (:,1:km)    = t1    (:,1:km) +
     &                        tdt   (:,1:km) * dt2      ! before "canopy_to_resolved"

! All mass tracers (excepr TKE) on original model layers after diffusion for use in "canopy_to_resolved"
          q2    (:,1:km, 1:ntrac1) = q1 (:,1:km, 1:ntrac1) +
     &                              rtg (:,1:km, 1:ntrac1) * dt2   ! before "canopy_to_resolved"

! TKE tracers on original model layers after diffusion for use in "canopy_to_resolved"
          q2    (:,1:km,   ntke ) = q1 (:,1:km,   ntke ) +
     &                             rtg (:,1:km,   ntke ) * dt2   ! before "canopy_to_resolved"
          q2 (:,:, ntke) = max(q2 (:,:, ntke), tkmin)            ! before "canopy_to_resolved"

! Temperature & wind components
          t2_mod  (:,1:km)   = t2  (:,1:km)
          u2_mod  (:,1:km)   = u2  (:,1:km)
          v2_mod  (:,1:km)   = v2  (:,1:km)

! Wind speed
          ws2_mod (:,1:km)   = sqrt(u2_mod(:,1:km)**2+v2_mod(:,1:km)**2) ! after diffusion

! Tracers
          q2_mod  (:,1:km,:) = q2  (:,1:km, :)          ! before "canopy_to_resolved"
          tkeh_mod(:,1:km)   = tkeh(:,1:km)             ! before "canopy_to_resolved"

! Tendencies
          rtg_mod(:,1:km, :) = rtg(:,1:km, :)           ! before "canopy_to_resolved" (km)

          tdt_mod(:,1:km)    = tdt(:,1:km)              ! before "canopy_to_resolved" (km)

          du_mod (:,1:km)    = du (:,1:km)              ! before "canopy_to_resolved" (km)
          dv_mod (:,1:km)    = dv (:,1:km)              ! before "canopy_to_resolved" (km)
          duv_mod(:,1:km)    = duv(:,1:km)              ! before "canopy_to_resolved" (km)

! Before "canopy_to_resolved"
!          aux3d(:,:,6) = Q2_MOD(:,:, ntke) ! ntke=198 "tke  before "canopy_to_resolved"

!          aux3d(:,:,4) = Q2_MOD(:,:, 11  ) ! n=11   "no"    before "canopy_to_resolved"
!          aux3d(:,:,4) = Q2_MOD(:,:, 12  ) ! n=12   "o3"    before "canopy_to_resolved"
!          aux3d(:,:,4) = Q2_MOD(:,:, 10  ) ! n=10   "no2"   before "canopy_to_resolved"

!          aux3d(:,:,6) = Q2_MOD(:,:, ntoz) ! ntoz=7 "o3mr"  before "canopy_to_resolved"
!          aux3d(:,:,2) = Q2_MOD(:,:, ntqv) ! ntqv=1 humidity before "canopy_to_resolved"

! Air Density after diffusion model layers
          rho2 (:,1:km) = prsl (:,1:km)/
     &                  (rd*t2 (:,1:km)*
     &           (1.+fv*max(q2 (:,1:km, ntqv),qmin))) ! ntqv=1 before "canopy_to_resolved"

! Set non-canopy columns to resolved values
! NB. Only vars not ALREADY defined in non-canopy columns

! Above-canopy layers (canopy layers below) non-canopy columns
         do k = 1, km-nkc
           ! km     is top combined subset
           ! nkc+1  is bot combined
           kc= nkc+k     ! 4th from top (nkt) to nkc+1 combined canopy plus resolved model layer

           do i = 1,im
           IF (FRT_MASK(i) <=  0.) THEN

! Above-canopy U-Wind/V-Wind/Temp     non-canopy columns
            u2_can3 (i,kc) = u2 (i,k)
            v2_can3 (i,kc) = v2 (i,k)
            t2_can3 (i,kc) = t2 (i,k)

! Above-canopy TKE tracer half layers non-canopy columns
            TKEH_CAN(i,kc) = TKEH(i,k)       ! J s-1 after diffusion (km)

! Tendencies
            TDT_CAN (i,kc) = TDT (i,k)       ! K s-1 after diffusion (km)
            DU_CAN  (i,kc) = DU  (i,k)       ! m s-2 after diffusion (km)
            DV_CAN  (i,kc) = DV  (i,k)       ! m s-2 after diffusion (km)
            DUV_CAN (i,kc) = DUV (i,k)       ! m s-2 after diffusion (km)

! All tendencies except TKE non-canopy columns
             RTG_CAN (i,kc, 1:ntrac1) = RTG (i,k, 1:ntrac1) ! kg kg-1 s-1
! TKE Tendency non-canopy columns
             RTG_CAN (i,kc, ntke) = RTG (i,k, ntke ) !         s-1

           END IF ! (FRT_MASK)
           end do ! i=1,im
         end do ! k = 1, km-nkc

! Canopy layers non-canopy columns
         do kc = 1, nkc ! 3-nkc canopy layers
           do i = 1,im
           IF (FRT_MASK(i) <=  0.) THEN

! Canopy layers U-wind/V-wind/Temp non-canopy columns
            u2_can3 (i,kc) = u2 (i,1)
            v2_can3 (i,kc) = v2 (i,1)
            t2_can3 (i,kc) = t2 (i,1)

! Canopy layers TKE half layers non-canopy columns
            TKEH_CAN(i,kc) = TKEH(i,1)       ! J S-1 (km)

! Canopy layers Tendencies momentum and heat non-canopy columns
            TDT_CAN (i,kc) = TDT (i,1)       ! K s-1 (km)
            DU_CAN  (i,kc) = DU  (i,1)       ! m s-2 (km)
            DV_CAN  (i,kc) = DV  (i,1)       ! m s-2 (km)
            DUV_CAN (i,kc) = DUV (i,1)       ! m s-2 (km)

! Canopy layers Tendencies mass tracers non-canopy columns
            RTG_CAN (i,kc, 1:ntrac1) = RTG (i,1, 1:ntrac1) ! kg kg-1 s-1 (km)
! Canopy layers Tendency TKE tracer     non-canopy columns
            RTG_CAN (i,kc,   ntke  ) = RTG (i,1,   ntke)   ! J s-1 s-1 (km)

           ENDIF ! (FRT_MASK)
           end do  ! do i=1,im
         end do ! kc=1,nkc

! Combined layers (1,km) Humidity after diffusion
! Apply minimum value on humidity qmin  before "canopy_to_resolved"
         q2_can (:,:, ntqv) = q1_can (:,:, ntqv) +
     &                       rtg_can (:,:, ntqv) * dt2                ! after diffusion (km)
         q2_can (:,:, ntqv) = max(q2_can (:,:, ntqv),qmin)            ! after diffusion (km)

! Output pbl diags
! Reproducing v16 but values are too large, check q2_can (:,:, ntke)
        aux3d(:,:, 4) = tdt_can (:,:)       ! after diffusion
!       aux3d(:,:, 2) =  dv_can (:,:)       ! after diffusion 
!       aux3d(:,:, 2) =  du_can (:,:)       ! after diffusion 
!
! duv_can calculated below

! Winds & temperature on combined layers after diffusion

! Top 3 combined layers NOT calculated in "sat_can": fill with top 3 model layers
          u2_can3 (:,km+nkc) = u2 (:,km  )  ! after diffusion (nkt)
          u2_can3 (:,km+2  ) = u2 (:,km-1)  ! after diffusion (nkt)
          u2_can3 (:,km+1  ) = u2 (:,km-2)  ! after diffusion (nkt)

! Top 3 combined layers NOT calculated in "sat_can": fill with top 3 model layers
          v2_can3 (:,km+nkc) = v2 (:,km  )  ! after diffusion (nkt)
          v2_can3 (:,km+2  ) = v2 (:,km-1)  ! after diffusion (nkt)
          v2_can3 (:,km+1  ) = v2 (:,km-2)  ! after diffusion (nkt)

! Top 3 combined layers NOT calculated in "sat_can": fill with top 3 model layers
          t2_can3 (:,km+nkc) = t2 (:,km  )  ! after diffusion (nkt)
          t2_can3 (:,km+2  ) = t2 (:,km-1)  ! after diffusion (nkt)
          t2_can3 (:,km+1  ) = t2 (:,km-2)  ! after diffusion (nkt)

! All tracers on combined layers after diffusion, for use in "canopy_to_resolved"

! Top 3 combined layers NOT calculated in "sat_can": fill with top 3 model layers
          q2_can3 (:,km+nkc ,   ntke  ) = q2 (:,km   ,   ntke  )         ! after diffusion (nkt)
          q2_can3 (:,km+2   ,   ntke  ) = q2 (:,km-1 ,   ntke  )         ! after diffusion (nkt)
          q2_can3 (:,km+1   ,   ntke  ) = q2 (:,km-2 ,   ntke  )         ! after diffusion (nkt)

! Top 3 combined layers NOT calculated in "sat_can": fill with top 3 model layers
          q2_can3 (:,km+nkc , 1:ntrac1) = q2 (:,km   , 1:ntrac1)         ! nkt after diffusion (nkt)
          q2_can3 (:,km+2   , 1:ntrac1) = q2 (:,km-1 , 1:ntrac1)         ! nkt after diffusion (nkt)
          q2_can3 (:,km+1   , 1:ntrac1) = q2 (:,km-2 , 1:ntrac1)         ! nkt after diffusion (nkt)

! Wind and temperature after diffusion

! Combined layers (1,km) U-Wind/V-Wind/Temp after diffusion
          u2_can3 (:,1:km) = u1_can (:,1:km) +  du_can(:,1:km) * dt2    ! after diffusion (nkt)
          v2_can3 (:,1:km) = v1_can (:,1:km) +  dv_can(:,1:km) * dt2    ! after diffusion (nkt)

! Combined layers (1,km) Temperature after diffusion
          t2_can3 (:,1:km) = t1_can (:,1:km) + tdt_can(:,1:km) * dt2    ! after diffusion (nkt)
!          t2_can (:,1:km) = t1_can (:,1:km) + tdt_can(:,1:km) * dt2    ! after diffusion (km)

! Wind Speed after diffusion on canopy layers
          ws2_can3(:,1:km) = sqrt(u2_can3(:,1:km)**2+v2_can3(:,1:km)**2)

! 10m-interpolated wind ws1_can3
!         duv_can (:,1:km) = (ws2_can3(:,1:km) - ws1_can3(:,1:km)) * rdt ! 10m-interpolated ws1_can3
! 1hy model layer wind  ws1_can (1-4cy)
          duv_can (:,1:km) = (ws2_can3(:,1:km) - ws1_can (:,1:km)) * rdt ! ws1_can is using 1hy model layer u1&v1(:,1)

! Output pbl diags
          aux3d(:,:, 2) =  duv_can (:,1:km)   ! after diffusion

! Tracers after diffusion
! Subset (1:km) combined layers out of total ntk layers (NB. dim(:,nkt,:) <= dim(:,km,:)
          q2_can3 (:,1:km, 1:ntrac1) = q1_can (:,1:km, 1:ntrac1) +
     &                                rtg_can (:,1:km, 1:ntrac1) * dt2   ! after diffusion (nkt)
! TKE after diffusion
          q2_can3 (:,1:km,   ntke  ) = q1_can (:,1:km,   ntke  ) +
     &                                rtg_can (:,1:km,   ntke  ) * dt2   ! after diffusion (nkt)
          q2_can  (:,1:km,   ntke  ) = q1_can (:,1:km,   ntke  ) +
     &                                rtg_can (:,1:km,   ntke  ) * dt2   ! after diffusion (km)

! Apply minimum value on TKE tracer before "canopy_to_resolved"
          q2_can3 (:,:, ntke) = max(q2_can3 (:,:, ntke), tkmin)           ! after diffusion (nkt)
          q2_can  (:,:, ntke) = max(q2_can  (:,:, ntke), tkmin)           ! after diffusion (km)

! Apply minimum value on humidity qmin before "canopy_to_resolved" and the tendency update
          q2_can3 (:,:, ntqv) = max(q2_can3 (:,:, ntqv), qmin   )  ! ntqv=1

! Apply minimum value on chemical conc  before "canopy_to_resolved"
          q2_can3 (:,:, ntoz:ntrac1 ) =
     &                         max(q2_can3 (:,:, ntoz:ntrac1), concmin)

! Output pbl diags
!         aux3d(:,:, 6) = q2_can3 (:,1:km, ntoz)  ! ntoz=7   "o3mr" after diffusion
!         aux3d(:,:, 4) = q2_can3 (:,1:km, 10  )  ! n=10     "o3"   after diffusion
!         aux3d(:,:, 2) = q2_can3 (:,1:km, ntqv)  ! ntqv=1   humidity

! Top 3 combined layers set to resolved
!    NB. Q2_can3 tracers array & t2_can3  after diffusion only updated 1:km
          rho2_can3 (:,km+nkc) = prsl (:,km  )/          ! after diffusion
     &                        (rd*t2 (:,km  )*
     &                 (1.+fv*max(q2 (:,km  , ntqv),qmin))) ! ntqv=1
          rho2_can3 (:,km+2  ) = prsl (:,km-1)/          ! after diffusion
     &                        (rd*t2 (:,km-1)*
     &                 (1.+fv*max(q2 (:,km-1, ntqv),qmin))) ! ntqv=1

          rho2_can3 (:,km+1  ) = prsl (:,km-2)/          ! after diffusion
     &                        (rd*t2 (:,km-2)*
     &                 (1.+fv*max(q2 (:,km-2, ntqv),qmin))) ! ntqv=1

! Air density after diffusion on canopy layers
          rho2_can3 (:,1:km   ) = prsl_can3 (:,1:km)/           ! after diffusion
     &                          (rd*t2_can3 (:,1:km)*
     &                   (1.+fv*max(q2_can3 (:,1:km, ntqv),qmin))) ! ntqv=1

! Output pbl diags
!         aux3d(:,:, 2) = rho2_can3 (:,1:km   )


      endif !do_canopy .and. cplaqm
!IVAI

!IVAI
      if (do_canopy) then

!===============================================================================
! Gather tracer concentration from canopy layers into model resolved layers (flag = 1)
!===============================================================================

!!!??????split "canopy_to_resolved" and "resolved_to_canopy" as separate routines !!!!!!!!!

        CALL canopy_transfer_run(im, im, km,         !in
     &        ntrac1, ntoz,                          !in
     &        garea,                                 !in
     &        zi, zl, zm,                            !in
     &        q2, rho2,                              !in kg kg-1
     &        1,                                     !in 1 = "canopy_to_resolved"
     &        FRT_MASK,                              !in
     &        kmod, kcan3,                           !in
     &        zmom_can3, zmid_can3,                  !in
     &        PRSL_CAN3, rho2_CAN3,                  !in
     &        Q2_MOD, Q2_CAN3, Q2_2M,                !inout kg kg-1 after diffusion
     &        massair_can3, massair,                 !inout
     &        mmr_o3_can3,                           !inout
     &        nfrct, ifrct,                          !inout
     &        frctr2c, frctc2r,                      !inout
     &        errmsg, errflg )

        if (errflg /= 0) return

! ==============
!  Input:
!   Q2_CAN3(:,:, NLAYT, NSPCSD) :  Chemical tracers mass conc. kg kg-1 on combined canopy+resolved layers after diffusion
!
!  InOutput
!   Q2_MOD(:,:, NLAYS, NSPCSD) :  Chemical tracers mass conc. kg kg-1 on model levels after diffusion
!   Q2                            Canopy columns only!!!
!
! ================================


! 
! Canopy columns calculated in "canopy_to_resolved"
!
! Non-Canopy columns filled with resolved values
!
        do k = 1, km
          do i = 1,im

! Non-canopy columns set to resolved values
          IF (FRT_MASK(i) <=  0.) THEN

! Non-canopy columns after "canopy_to_resolved"  Line 1843: forrtl: error (78): process killed (SIGTERM)
            q2_mod (i,k, 1:ntrac1) = q2(i,k, 1:ntrac1)       ! after diffusion & after "canopy_to_resolved"
            q2_mod (i,k,   ntke  ) = q2(i,k,   ntke)         ! after diffusion

            t2_mod (i,k) = t2 (i,k)   ! after diffusion
            u2_mod (i,k) = u2 (i,k)   ! after diffusion
            v2_mod (i,k) = v2 (i,k)   ! after diffusion

            ws2_mod(i,k) = sqrt(u2_mod(i,k)**2+v2_mod(i,k)**2) ! after diffusion

          END IF ! (FRT_MASK)
          end do ! i=1,im
        end do ! k = 1, km-nkc

! Apply minimum value on chemical conc  after "canopy_to_resolved"
        q2_mod  (:,:, ntoz:ntrac1 ) =
     &                         max(q2_mod  (:,:, ntoz:ntrac1), concmin)  ! after "canopy_to_resolved"

! Apply minimum value on humidity qmin  before doing canopy_transfer & update tendency
        q2_mod(:,:, ntqv) = max(q2_mod(:,:, ntqv),qmin)

! Output 2D pbl diags
!        aux2d(:, 6) = rtg_2m(:, 11  )         ! after diffusion n=11 "no"
!        aux2d(:, 4) = rtg_2m(:, 12  )         ! after diffusion n=12 "o3"
!        aux2d(:, 2) = rtg_2m(:, 10  )         ! after diffusion n=10 "no2"

!        aux2d(:, 6) = Q2_2m (:, 11  )         ! after diffusion n=11 "no"
!        aux2d(:, 4) = Q2_2m (:, 12  )         ! after diffusion n=12 "o3"
!        aux2d(:, 2) = Q2_2m (:, 10  )         ! after diffusion n=10 "no2"

! 2-m diag is always 1cy layer
! GFSv16 - GFDL misrophysics

! GFSv17 - MP microphysics

! Output 3D pbl diags

!       aux3d(:,:,6) = ws2_mod(:,:)   ! ws2
!       aux3d(:,:,4) = t2_mod(:,:)       ! t2
!       aux3d(:,:,2) = q2_mod(:,:, ntke) ...? ! ntke=198 TKE

!       aux3d(:,:,4) = Q2_MOD(:,:, 11  ) ! n=11 no after "canopy_to_resolved"
!       aux3d(:,:,4) = Q2_MOD(:,:, 12  ) ! n=12 o3 after "canopy_to_resolved"
!       aux3d(:,:,4) = Q2_MOD(:,:, 10  ) ! n=10 no2 after "canopy_to_resolved"

!       aux3d(:,:,6) = Q2_MOD(:,:, ntoz) ! ntoz=7  after "canopy_to_resolved" 
!       aux3d(:,:,2) = Q2_MOD(:,:, ntqv) ! ntqv=1  after "canopy_to_resolved"
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do k = 1, km1-1     ! from bottom to top resolved model levels

          kount=0 !IVAI
          ZOOOX(:) = 1.
          do i = 1, im

! Canopy columns/grid cells
! There is a contiguous forest canopy, calculated integrated canopy correction over model layers
          IF ( FRT_MASK(i) >  0. ) THEN

            FCH  = cfch(i)

! Determine if canopy inside the model layer (kcan=1) or not (kcan=0)
            IF (k .EQ. 1) THEN !use model layer interfaces
              KCAN = 1
            ELSE
              IF ( cfch(i) > zi(i,k) .AND. cfch(i) <= zi(i,k+1) ) THEN
                KCAN = 1
              ELSE
                KCAN= 0
              END IF
            END IF

            IF (KCAN == 1 ) THEN !canopy inside model layer k

              ZI05 = zi(i,k+1) ! Initialize each model layer top that contains canopy (m)
                               ! Integrate across total model interface
              ZFL =  zi(i,k+1) ! Set ZFL = ZI05
              COUNTCAN = 0   ! Initialize canopy layers

              IF (k .EQ. 1) THEN  !Find bottom in each model layer
                  BOTCAN = 0.0
              ELSE
                  BOTCAN = zi(i,k)
              END IF

              DZFL =  zi(i,k+1) - BOTCAN

! Canopy inside model layer
              DO WHILE (ZI05 .GE. BOTCAN)

!!!!!!!!!!!!!!!!!!!!!!!!
! Steping down in-canopy zi_can3(:,1)= 0.

                 do kc = nkt-1, 1, -1         ! from top to bottom combined canopy plus resolved  model layers

!!!!!!!!!!!!!!!!!!!!!!!!

! Between two canopy layers
                 IF ( ZI05 >  zi_can3(i,kc  ) .and.
     &                ZI05 <= zi_can3(i,kc+1) ) THEN

                   kc_can   = kc
                   FZI05    = 1./(zi_can3 (i,kc+1) - zi_can3 (i,kc)) *
     &                           (ZI05             - zi_can3 (i,kc))

                   TTCORR   = (t2_can3 (i,kc+1) - t2_can3 (i,kc))/
     &                        (zi_can3 (i,kc+1) - zi_can3 (i,kc)) *
     &                        (ZI05             - zi_can3 (i,kc))

! U-Wind/V-Wind after diffusion on canopy layers
                   UUCAN  =  u2_can3 (i,kc  ) +
     &                      (u2_can3 (i,kc+1) - u2_can3 (i,kc))/
     &                      (zi_can3 (i,kc+1) - zi_can3 (i,kc)) *
     &                      (ZI05             - zi_can3 (i,kc))

                   VVCAN  =  v2_can3 (i,kc  ) +
     &                      (v2_can3 (i,kc+1) - v2_can3 (i,kc))/
     &                      (zi_can3 (i,kc+1) - zi_can3 (i,kc)) *
     &                      (ZI05             - zi_can3 (i,kc))
! Temperature after diffusion on canopy layers
                   TTCAN  =  t2_can3 (i,kc  ) +
     &                      (t2_can3 (i,kc+1) - t2_can3 (i,kc))/
     &                      (zi_can3 (i,kc+1) - zi_can3 (i,kc)) *
     &                      (ZI05             - zi_can3 (i,kc))


! TKE after diffusion on canopy layers
                   TKECAN =  q2_can3 (i,kc, ntke) +
     &                      (q2_can3 (i,kc+1, ntke) -
     &                       q2_can3 (i,kc  , ntke))/
     &                      (zi_can3 (i,kc+1) - zi_can3 (i,kc)) *
     &                      (ZI05             - zi_can3 (i,kc))

! TKE half layers
                   TKEHCAN =  tkeh_can (i,kc  ) +
     &                      (tkeh_can (i,kc+1) - tkeh_can (i,kc))/
     &                      (  zi_can3(i,kc+1) -   zi_can3(i,kc)) *
     &                      (ZI05              -   zi_can3(i,kc))

                 END IF !  "zi_can3(kc) < ZI05 <= zi_can3(kc+1)"

!!!!!!!!!!!!!! 
! End steping down in-canopy
!
                 end do  ! kc = nkt-1, 1, -1         ! from top to bottom combined canopy plus resolved  model layers
!
!!!!!!!!!!!!!!

                 IF ( ZI05 .LE. ZFL )  THEN  ! Model layers

                   COUNTCAN = COUNTCAN + 1
                   ZCANX  (COUNTCAN) = ZI05

! U-Wind/V-Wind after diffusion on model layer
                   UUX  (COUNTCAN) = UUCAN
                   VVX  (COUNTCAN) = VVCAN
! Temperature after diffusion on model layers
                   TTX  (COUNTCAN) = TTCAN
! TKE on model layers
                   TKEX (COUNTCAN) = TKECAN
! TKE half layers on model layers
                   TKEHX(COUNTCAN) = TKEHCAN

                 END IF  ! ( ZI05 .LE. ZFL )

                 ZI05 = ZI05-0.5  !step down in-canopy resolution of 0.5m

              END DO ! DO WHILE (ZI05.GE.BOTCAN)

! IntegrateTrapezoid = sum((y(1+1:n-0) + y(1+0:n-1))*
!     &                            (x(1+1:n-0) - x(1+0:n-1)))/2
              ZZ_INT= IntegrateTrapezoid(
     &                  ZCANX(COUNTCAN:1:-1), ZOOOX(COUNTCAN:1:-1))

! U-wind
              UU_INT = IntegrateTrapezoid(
     &                   ZCANX(COUNTCAN:1:-1) ,
     &                     UUX(COUNTCAN:1:-1) ) /
     &                     ZZ_INT

! V-wind
              VV_INT = IntegrateTrapezoid(
     &                   ZCANX(COUNTCAN:1:-1) ,
     &                     VVX(COUNTCAN:1:-1) ) /
     &                     ZZ_INT

! Temp
              TT_INT = IntegrateTrapezoid(
     &                   ZCANX(COUNTCAN:1:-1) ,
     &                    TTX(COUNTCAN:1:-1)   ) /
     &                    ZZ_INT

! TKE
              TKE_INT = IntegrateTrapezoid(
     &                    ZCANX(COUNTCAN:1:-1) ,
     &                     TKEX(COUNTCAN:1:-1)  ) /
     &                     ZZ_INT

! TKEH
              TKEH_INT= IntegrateTrapezoid(
     &                    ZCANX(COUNTCAN:1:-1) ,
     &                    TKEHX(COUNTCAN:1:-1) ) /
     &                    ZZ_INT

! Sum
              TT_SUM = sum(  TTX(COUNTCAN:1:-1))/COUNTCAN

! U-wind/V-wind Canopy Columns
              u2_mod (i,k)   = UU_INT     ! after "canopy-to-resolved"
              v2_mod (i,k)   = VV_INT     ! after "canopy-to-resolved"

              ws2_mod(i,k) = sqrt(u2_mod(i,k)**2+v2_mod(i,k)**2) ! after "canopy-to-resolved"

! Temperature Canopy Columns
              t2_mod (i,k)   = TT_INT     ! after "canopy-to-resolved"

! TKE Canopy Columns
              q2_mod (i,k, ntke) = TKE_INT    ! after "canopy-to-resolved"
              q2_mod (i,k, ntke) = max(q2_mod (i,k, ntke), tkmin)           ! after "canopy-to-resolved"

! TKEH Canopy Columns
              tkeh_mod(i,k) = TKEH_INT   ! after "canopy-to-resolved"
              tkeh_mod(i,k) = max(tkeh_mod(i,k), tkmin)
! Apply minimum value on TKE tracer before "canopy_to_resolved"

            END IF  ! (KCAN .EQ. 1) model layer(s) containing canopy

          END IF  ! contigous canopy conditions
          enddo   ! i = 1, im

          kount = kount + 1

        enddo ! k = 1, km1-1           ! from bottom to top resolved model levels



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Update ALL tracers with in-canopy tendencies (average sub-canopy values )
! Here just wind components, temperature TKE, and interstitial tracers 
! The chemical tracers are updated below with values from "canopy_to_resolved" transfer
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Update ALL model layer
       do k = 1,km

          do i = 1,im
          ! Non-Canopy columns set to resolved
          IF (FRT_MASK(i) <= 0.) THEN

! Non-canopy columns after "canopy_to_resolved"
! already assigned before "canopy_to_resolved"
!            du_mod (i,k) = du (i,k)                      ! after "canopy_to_resolved"
!            dv_mod (i,k) = dv (i,k)                      ! after "canopy_to_resolved"
!            tdt_mod(i,k) = tdt(i,k)                      ! after "canopy_to_resolved"

! Non-canopy columns
! already assigned before "canopy_to_resolved"
!            rtg_mod(i,k, ntke) = rtg(i,k, ntke)          ! after "canopy_to_resolved"

! Non-canopy columns
! already assigned before "canopy_to_resolved"
!            rtg_mod(i,k, 1:ntrac1) = rtg(i,k, 1:ntrac1)  ! after "canopy_to_resolved"

          ! Update Canopy columns only
          ELSE IF (FRT_MASK(i) >  0.) THEN

! Canopy Columns
! U-Wind/V-Wind after sub-canopy diffusion
            du_mod (i,k)  = (u2_mod (i,k) - u1 (i,k) )*rdt  ! after "canopy_to_resolved" dim(km)
            dv_mod (i,k)  = (v2_mod (i,k) - v1 (i,k) )*rdt  ! after "canopy_to_resolved" dim(km)

            duv_mod(i,k)  = (ws2_mod(i,k) - ws1(i,k) )*rdt  ! after "canopy_to_resolved" dim(km)

! Temperature after sub-canopy diffusion
            tdt_mod(i,k) = (t2_mod(i,k) - t1(i,k) )*rdt  ! after "canopy_to_resolved" dim(km)

! Tendency mass tracers after sub-canopy diffusion
            rtg_mod(i,k, 1:ntrac1) = (q2_mod(i,k, 1:ntrac1) -
     &                                q1    (i,k, 1:ntrac1))*rdt  ! after "canopy_to_resolved"
! Tendency TKE  tracer  after sub-canopy diffusion
            rtg_mod(i,k, ntke)     = (q2_mod(i,k,   ntke)   -
     &                                q1    (i,k,   ntke)  )*rdt  ! after "canopy_to_resolved"

!!!!!!!!!!!!!!!!!!!!!

          ENDIF ! (FRT_MASK)
          enddo ! do i=1,im
        enddo ! k=1,km


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Update Met & TKE & MP (microphysics) cloud fields
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Update all model layers
        do k = 1,km
          do i = 1,im

          IF (FRT_MASK(i) >  0.) THEN
! Tendency Wind components
!             du (i,k) = du_mod (i,k)            ! <<<<<<<========== UPDATE UU  TEND =========>>>>>>>
!             dv (i,k) = dv_mod (i,k)            ! <<<<<<<========== UPDATE VV  TEND =========>>>>>>>
! Tendency Temperature
!             tdt(i,k) = tdt_mod(i,k)            ! <<<<<<<========== UPDATE TT  TEND =========>>>>>>>
! Tendency TKE (ntke=198)
!             rtg(i,k, ntke) = rtg_mod(i,k, ntke) ! <<<<<<<========== UPDATE TKE TEND =========>>>>>>>

! TKE half layers
!             tkeh(i,k) = tkeh_mod(i,k)           ! <<<<<<<========== UPDATE TKEH     =========>>>>>>>

          ENDIF ! Contiguous canopy
          enddo ! i
        enddo !k



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Update ONLY chemical tracers (n=8, ntrac1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! All chemical tracers (n=9, ntrac1)
! ntqv=1-8 are cloud/rain and "sgs_tke"
! n=8   "no2" (ntchs=9) GFDL
! n=9   "no"            GFDL
! n=10  "o3"            GFDL
! ...
! NTRAC1 = 196
! ---NO PBL TEND --------
! n=197 ntche
! n=198 ntke
! -----------------------
        DO n =  ntchs-1, ntche-1 ! 10, ntche-1 (same as NTRAC1)

! Update all model layers
        do k = 1,km

          do i = 1,im
          IF (FRT_MASK(i) >  0.) THEN

!!!!!!!!!!!!!!!!!!!!!!!!!
! Select all chemical tracers (selected tracers above)
!!!!!!!!!!!!!!!!!!!!!!!!!
            rtg(i,k, n) = rtg_mod(i,k, n)   ! <<<<<<<========== UPDATE TEND =========>>>>>>>
          END IF
          end do

        end do ! k = 1,km
        end do ! n = 1, NTRAC1

      endif !if(do_canopy)

!IVAI

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!> ## Save PBL height for diagnostic purpose
!
      if (do_canopy) then

      do i = 1, im
         hpbl(i) = hpbl_can(i)
         kpbl(i) = kpbl_can(i)
      enddo

! Output pbl diags 2D
!       aux2d(:, 2) = hpbl_can(:) ! after canopy

      endif  !do_canopy
!IVAI
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      return
      end subroutine canopy_driver_run
!> @}
      end module canopy_driver

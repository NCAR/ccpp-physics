!> \file module_mynnsfc_wrapper.F90
!!  Contains all of the code related to running the MYNN surface layer scheme 

      MODULE mynnsfc_wrapper

          USE module_sf_mynn

          !Global variables:
          INTEGER, PARAMETER :: psi_opt = 0   !0: MYNN
                                              !1: GFS

      contains

!! \section arg_table_mynnsfc_wrapper_init Argument Table
!! \htmlinclude mynnsfc_wrapper_init.html

!!
      subroutine mynnsfc_wrapper_init(errmsg, errflg)

         character(len=*), intent(out) :: errmsg
         integer, intent(out) :: errflg

         ! Initialize CCPP error handling variables
         errmsg = ''
         errflg = 0

         ! initialize tables for psih and psim (stable and unstable)
         CALL PSI_INIT(psi_opt,errmsg,errflg)

         IF (debug_code >= 1) THEN
           print*,"CHECK INITIALIZATION OF PSI:"
           print*,"psim_stab(0-1):",psim_stab(0),psim_stab(1)
           print*,"psih_stab(0-1):",psih_stab(0),psih_stab(1)
           print*,"psim_unstab(0-1):",psim_unstab(0),psim_unstab(1)
           print*,"psih_unstab(0-1):",psih_unstab(0),psih_unstab(1)
         ENDIF

      end subroutine mynnsfc_wrapper_init

      subroutine mynnsfc_wrapper_finalize ()
      end subroutine mynnsfc_wrapper_finalize

!>\defgroup gsd_mynn_sfc GSD MYNN Surface Layer Scheme Module
!> \brief This scheme (1) performs pre-mynnsfc work, (2) runs the mynn sfc layer scheme, and (3) performs post-mynnsfc work
!! \section arg_table_mynnsfc_wrapper_run Argument Table
!! \htmlinclude mynnsfc_wrapper_run.html
!!
!###===================================================================
SUBROUTINE mynnsfc_wrapper_run(            &
     &  im,levs,                           &
     &  itimestep,iter,                    &
     &  flag_init,flag_restart,lsm,lsm_ruc,&
     &  sigmaf,vegtype,shdmax,ivegsrc,     &  !intent(in)
     &  z0pert,ztpert,                     &  !intent(in)
     &  redrag,sfc_z0_type,                &  !intent(in)
     &  delt,dx,                           &
     &  u, v, t3d, qvsh, qc, prsl, phii,   &
     &  exner, ps, PBLH, slmsk,            &
     &         wet,       dry,       icy,  &  !intent(in)
     &   tskin_ocn, tskin_lnd, tskin_ice,  &  !intent(in)
     &   tsurf_ocn, tsurf_lnd, tsurf_ice,  &  !intent(in)
     &    qsfc_ocn,  qsfc_lnd,  qsfc_ice,  &  !intent(in)
     &   snowh_ocn, snowh_lnd, snowh_ice,  &  !intent(in)
     &     znt_ocn,   znt_lnd,   znt_ice,  &  !intent(inout)
     &     ust_ocn,   ust_lnd,   ust_ice,  &  !intent(inout)
     &      cm_ocn,    cm_lnd,    cm_ice,  &  !intent(inout)
     &      ch_ocn,    ch_lnd,    ch_ice,  &  !intent(inout)
     &      rb_ocn,    rb_lnd,    rb_ice,  &  !intent(inout)
     &  stress_ocn,stress_lnd,stress_ice,  &  !intent(inout)
     &      fm_ocn,    fm_lnd,    fm_ice,  &  !intent(inout)
     &      fh_ocn,    fh_lnd,    fh_ice,  &  !intent(inout)
     &    fm10_ocn,  fm10_lnd,  fm10_ice,  &  !intent(inout)
     &     fh2_ocn,   fh2_lnd,   fh2_ice,  &  !intent(inout)
     &    hflx_ocn,  hflx_lnd,  hflx_ice,  &
     &    qflx_ocn,  qflx_lnd,  qflx_ice,  &
     &  QSFC, qsfc_lnd_ruc, qsfc_ice_ruc,  &
     &  USTM, ZOL, MOL,                    &
     &  RMOL, WSPD, ch, HFLX, QFLX, LH,    &
     &  FLHC, FLQC,                        &
     &  U10, V10, TH2, T2, Q2,             &
     &  wstar, CHS2, CQS2,                 &
!     &  CP, G, ROVCP, R, XLV,           &
!     &  SVP1, SVP2, SVP3, SVPT0,        &
!     &  EP1,EP2,KARMAN,                 &
     &  lprnt, errmsg, errflg           )


! should be moved to inside the mynn:
      use machine , only : kind_phys

!      use physcons, only : cp     => con_cp,              &
!     &                     g      => con_g,               &
!     &                     r_d    => con_rd,              &
!     &                     r_v    => con_rv,              &
!     &                     cpv    => con_cvap,            &
!     &                     cliq   => con_cliq,            &
!     &                     Cice   => con_csol,            &
!     &                     rcp    => con_rocp,            &
!     &                     XLV    => con_hvap,            &
!     &                     XLF    => con_hfus,            &
!     &                     EP_1   => con_fvirt,           &
!     &                     EP_2   => con_eps

!      USE module_sf_mynn, only : SFCLAY_mynn 

!------------------------------------------------------------------- 
      implicit none
!------------------------------------------------------------------- 
!  ---  constant parameters:
!      real(kind=kind_phys), parameter :: rvovrd  = r_v/r_d
      real(kind=kind_phys), parameter :: karman  = 0.4
!      real(kind=kind_phys), parameter :: XLS     = 2.85E6
!      real(kind=kind_phys), parameter :: p1000mb = 100000.
      real(kind=kind_phys), parameter :: SVP1    = 0.6112
      real(kind=kind_phys), parameter :: SVP2    = 17.67
      real(kind=kind_phys), parameter :: SVP3    = 29.65
      real(kind=kind_phys), parameter :: SVPT0   = 273.15

  REAL, PARAMETER :: xlvcp=xlv/cp, xlscp=(xlv+xlf)/cp, ev=xlv, rd=r_d, &
       &rk=cp/rd, svp11=svp1*1.e3, p608=ep_1, ep_3=1.-ep_2, g_inv=1./g


  character(len=*), intent(out) :: errmsg
  integer, intent(out) :: errflg

!MISC CONFIGURATION OPTIONS
      INTEGER, PARAMETER ::       &
     &       spp_pbl  = 0,        &
     &       isftcflx = 0,        & !control: 0
     &       iz0tlnd  = 0,        & !control: 0
     &       isfflx   = 1

      integer, intent(in) :: ivegsrc
      integer, intent(in) :: sfc_z0_type ! option for calculating surface roughness length over ocean
      logical, intent(in) :: redrag ! reduced drag coeff. flag for high wind over sea (j.han)

!Input data
      integer, dimension(im), intent(in) :: vegtype
      real(kind=kind_phys), dimension(im), intent(in)    ::       &
     &                    sigmaf,shdmax,z0pert,ztpert

!MYNN-1D
      REAL    :: delt
      INTEGER :: im, levs
      INTEGER :: iter, k, i, itimestep, lsm, lsm_ruc
      LOGICAL :: flag_init,flag_restart,lprnt
      INTEGER :: IDS,IDE,JDS,JDE,KDS,KDE,                   &
     &            IMS,IME,JMS,JME,KMS,KME,                  &
     &            ITS,ITE,JTS,JTE,KTS,KTE

      real(kind=kind_phys), dimension(im,levs+1),           &
     &      intent(in)  ::                  phii
      real(kind=kind_phys), dimension(im,levs),             &
     &      intent(in)  ::         exner, PRSL,             &
     &                     u, v, t3d, qvsh, qc

      real(kind=kind_phys), dimension(im,levs) ::           &
     &        pattern_spp_pbl, dz, th, qv

      logical, dimension(im), intent(in) :: wet, dry, icy

      real(kind=kind_phys), dimension(im), intent(in)    :: &
     &                    tskin_ocn, tskin_lnd, tskin_ice,  &
     &                    tsurf_ocn, tsurf_lnd, tsurf_ice,  &
     &                    snowh_ocn, snowh_lnd, snowh_ice

      real(kind=kind_phys), dimension(im), intent(inout) :: &
     &                      znt_ocn,   znt_lnd,   znt_ice,  &
     &                      ust_ocn,   ust_lnd,   ust_ice,  &
     &                       cm_ocn,    cm_lnd,    cm_ice,  &
     &                       ch_ocn,    ch_lnd,    ch_ice,  &
     &                       rb_ocn,    rb_lnd,    rb_ice,  &
     &                   stress_ocn,stress_lnd,stress_ice,  &
     &                       fm_ocn,    fm_lnd,    fm_ice,  &
     &                       fh_ocn,    fh_lnd,    fh_ice,  &
     &                     fm10_ocn,  fm10_lnd,  fm10_ice,  &
     &                      fh2_ocn,   fh2_lnd,   fh2_ice,  &
     &                     hflx_ocn,  hflx_lnd,  hflx_ice,  &
     &                     qflx_ocn,  qflx_lnd,  qflx_ice,  &
     &                     qsfc_ocn,  qsfc_lnd,  qsfc_ice

!MYNN-2D
      real(kind=kind_phys), dimension(:), intent(in)    ::  &
     &        dx, pblh, slmsk, ps,                          &
     &        qsfc_lnd_ruc, qsfc_ice_ruc

      real(kind=kind_phys), dimension(im), intent(inout) :: &
     &        ustm, hflx, qflx, wspd, qsfc,                 &
     &        FLHC, FLQC, U10, V10, TH2, T2, Q2,            &
     &        CHS2, CQS2, rmol, zol, mol, ch,               &
     &        lh, wstar
     !LOCAL
      real, dimension(im) ::                                &
     &        hfx, znt, psim, psih,                         &
     &        chs, ck, cd, mavail, xland, GZ1OZ0,           &
     &        cpm, qgh, qfx, qsfc_ruc

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

!      if (lprnt) then
!         write(0,*)"=============================================="
!         write(0,*)"in mynn surface layer wrapper..."
!         write(0,*)"flag_init=",flag_init
!         write(0,*)"flag_restart=",flag_restart
!         write(0,*)"iter=",iter
!      endif

      ! prep MYNN-only variables
      do k=1,2 !levs
        do i=1,im
           dz(i,k)=(phii(i,k+1) - phii(i,k))*g_inv
           th(i,k)=t3d(i,k)/exner(i,k)
           !qc(i,k)=MAX(qgrs(i,k,ntcw),0.0)
           qv(i,k)=qvsh(i,k)/(1.0 - qvsh(i,k))
           pattern_spp_pbl(i,k)=0.0
        enddo
      enddo
      do i=1,im
          if (slmsk(i)==1. .or. slmsk(i)==2.)then !sea/land/ice mask (=0/1/2) in FV3
            xland(i)=1.0                          !but land/water = (1/2) in SFCLAY_mynn
          else
            xland(i)=2.0
          endif
          qgh(i)=0.0
          mavail(i)=1.0
          !snowh(i)=snowd(i)*800. !mm -> m
          !znt_lnd(i)=znt_lnd(i)*0.01  !cm -> m
          !znt_ocn(i)=znt_ocn(i)*0.01  !cm -> m
          !znt_ice(i)=znt_ice(i)*0.01  !cm -> m
          cpm(i)=cp
      enddo

      ! cm -> m
      where (dry) znt_lnd=znt_lnd*0.01
      where (wet) znt_ocn=znt_ocn*0.01
      where (icy) znt_ice=znt_ice*0.01

      ! qsfc ruc
      qsfc_ruc = 0.0
      if (lsm==lsm_ruc) then
        where (dry) qsfc_ruc = qsfc_lnd_ruc
        where (icy) qsfc_ruc = qsfc_ice_ruc
      end if

!      if (lprnt) then
!          write(0,*)"CALLING SFCLAY_mynn; input:"
!          write(0,*)"T:",t3d(1,1),t3d(1,2),t3d(1,3)
!          write(0,*)"TH:",th(1,1),th(1,2),th(1,3)
!          write(0,*)"u:",u(1,1:3)
!          write(0,*)"v:",v(1,1:3)
!          !write(0,*)"qv:",qv(1,1:3,1)
!          write(0,*)"p:",prsl(1,1)
!          write(0,*)"dz:",dz(1,1)," qsfc=",qsfc(1)," rmol:",rmol(1)
!          write(0,*)"         land      water    ice"
!          write(0,*)dry(1),wet(1),icy(1)
!          write(0,*)"ust:",ust_lnd(1),ust_ocn(1),ust_ice(1)
!          write(0,*)"Tsk:",tskin_lnd(1),tskin_ocn(1),tskin_ice(1)
!          write(0,*)"Tsurf:",tsurf_lnd(1),tsurf_ocn(1),tsurf_ice(1)
!          write(0,*)"Qsfc:",qsfc_lnd(1),qsfc_ocn(1),qsfc_ice(1)
!          write(0,*)"sno:",snowh_lnd(1),snowh_ocn(1),snowh_ice(1)
!          write(0,*)"znt:",znt_lnd(1),znt_ocn(1),znt_ice(1)
!          !write(0,*)"HFX:",hfx(1)," qfx",qfx(1)
!          write(0,*)"qsfc:",qsfc(1)," ps:",ps(1)
!          write(0,*)"wspd:",wspd(1),"rb=",rb_ocn(1)
!          write(0,*)"delt=",delt," im=",im," levs=",levs
!          write(0,*)"flag_init=",flag_init 
!          write(0,*)"flag_restart=",flag_restart 
!          write(0,*)"iter=",iter
!          write(0,*)"zlvl(1)=",dz(1,1)*0.5
!          write(0,*)"PBLH=",pblh(1)," xland=",xland(1)
!       endif


        CALL SFCLAY_mynn(                                                     &
             u3d=u,v3d=v,t3d=t3d,qv3d=qv,p3d=prsl,dz8w=dz,                    &
             th3d=th,pi3d=exner,qc3d=qc,                                      &
             PSFCPA=ps,PBLH=pblh,MAVAIL=mavail,XLAND=xland,DX=dx,             &
             CP=cp,G=g,ROVCP=rcp,R=r_d,XLV=xlv,                               &
             SVP1=svp1,SVP2=svp2,SVP3=svp3,SVPT0=svpt0,                       &
             EP1=ep_1,EP2=ep_2,KARMAN=karman,                                 &
             ISFFLX=isfflx,isftcflx=isftcflx,LSM=lsm,                         &
             iz0tlnd=iz0tlnd,psi_opt=psi_opt,                                 &
    &        sigmaf=sigmaf,vegtype=vegtype,shdmax=shdmax,ivegsrc=ivegsrc,     & !intent(in)
    &        z0pert=z0pert,ztpert=ztpert,                                     & !intent(in)
    &        redrag=redrag,sfc_z0_type=sfc_z0_type,                           & !intent(in)
             itimestep=itimestep,iter=iter,                                   &
                         wet=wet,              dry=dry,              icy=icy, &  !intent(in)
             tskin_ocn=tskin_ocn,  tskin_lnd=tskin_lnd,  tskin_ice=tskin_ice, &  !intent(in)
             tsurf_ocn=tsurf_ocn,  tsurf_lnd=tsurf_lnd,  tsurf_ice=tsurf_ice, &  !intent(in)
               qsfc_ocn=qsfc_ocn,    qsfc_lnd=qsfc_lnd,    qsfc_ice=qsfc_ice, &  !intent(in)
             snowh_ocn=snowh_ocn,  snowh_lnd=snowh_lnd,  snowh_ice=snowh_ice, &  !intent(in)
                 znt_ocn=znt_ocn,      znt_lnd=znt_lnd,      znt_ice=znt_ice, &  !intent(inout)
                 ust_ocn=ust_ocn,      ust_lnd=ust_lnd,      ust_ice=ust_ice, &  !intent(inout)
                   cm_ocn=cm_ocn,        cm_lnd=cm_lnd,        cm_ice=cm_ice, &  !intent(inout)
                   ch_ocn=ch_ocn,        ch_lnd=ch_lnd,        ch_ice=ch_ice, &  !intent(inout)
                   rb_ocn=rb_ocn,        rb_lnd=rb_lnd,        rb_ice=rb_ice, &  !intent(inout)
           stress_ocn=stress_ocn,stress_lnd=stress_lnd,stress_ice=stress_ice, &  !intent(inout)
                   fm_ocn=fm_ocn,        fm_lnd=fm_lnd,        fm_ice=fm_ice, &  !intent(inout)
                   fh_ocn=fh_ocn,        fh_lnd=fh_lnd,        fh_ice=fh_ice, &  !intent(inout)
               fm10_ocn=fm10_ocn,    fm10_lnd=fm10_lnd,    fm10_ice=fm10_ice, &  !intent(inout)
                 fh2_ocn=fh2_ocn,      fh2_lnd=fh2_lnd,      fh2_ice=fh2_ice, &  !intent(inout)
               hflx_ocn=hflx_ocn,    hflx_lnd=hflx_lnd,    hflx_ice=hflx_ice, &
               qflx_ocn=qflx_ocn,    qflx_lnd=qflx_lnd,    qflx_ice=qflx_ice, &
             ch=ch,CHS=chs,CHS2=chs2,CQS2=cqs2,CPM=cpm,                       &
             ZNT=znt,USTM=ustm,ZOL=zol,MOL=mol,RMOL=rmol,                     &
             psim=psim,psih=psih,                                             &
             HFLX=hflx,HFX=hfx,QFLX=qflx,QFX=qfx,LH=lh,FLHC=flhc,FLQC=flqc,   &
             QGH=qgh,QSFC=qsfc,QSFC_RUC=qsfc_ruc,                             &
             U10=u10,V10=v10,TH2=th2,T2=t2,Q2=q2,                             &
             GZ1OZ0=GZ1OZ0,WSPD=wspd,wstar=wstar,                             &
             spp_pbl=spp_pbl,pattern_spp_pbl=pattern_spp_pbl,                 &
             ids=1,ide=im, jds=1,jde=1, kds=1,kde=levs,                       &
             ims=1,ime=im, jms=1,jme=1, kms=1,kme=levs,                       &
             its=1,ite=im, jts=1,jte=1, kts=1,kte=levs                        )


        !! POST MYNN SURFACE LAYER (INTERSTITIAL) WORK:
        !do i = 1, im
        !   !* Taken from sfc_nst.f
        !   !* ch         = surface exchange coeff heat & moisture(m/s) im
        !   !* rch(i)     = rho_a(i) * cp * ch(i) * wind(i)
        !   !* hflx(i)    = rch(i) * (tsurf(i) - theta1(i))  !K m s-1
        !   !* hflx(i)=hfx(i)/(rho(i,1)*cp) - now calculated inside module_sf_mynn.F90
        !   !* Taken from sfc_nst.f
        !   !* evap(i)    = elocp * rch(i) * (qss(i) - q0(i)) !kg kg-1 m s-1
        !   !NOTE: evap & qflx will be solved for later
        !   !qflx(i)=QFX(i)/
        !   !evap(i)=QFX(i)   !or /rho ??
        !   ! DH* note - this could be automated (CCPP knows how to convert m to cm)
        !   znt_lnd(i)=znt_lnd(i)*100.   !m -> cm
        !   znt_ocn(i)=znt_ocn(i)*100.
        !   znt_ice(i)=znt_ice(i)*100.
        !enddo

        ! m -> cm
        where (dry) znt_lnd=znt_lnd*100.
        where (wet) znt_ocn=znt_ocn*100.
        where (icy) znt_ice=znt_ice*100.

!      if (lprnt) then
!         write(0,*)
!         write(0,*)"finished with mynn_surface layer; output:"
!         write(0,*)"         land      water    ice"
!         write(0,*)dry(1),wet(1),icy(1)
!         write(0,*)"ust:",ust_lnd(1),ust_ocn(1),ust_ice(1)
!         write(0,*)"Tsk:",tskin_lnd(1),tskin_ocn(1),tskin_ice(1)
!         write(0,*)"Tsurf:",tsurf_lnd(1),tsurf_ocn(1),tsurf_ice(1)
!         write(0,*)"Qsfc:",qsfc_lnd(1),qsfc_ocn(1),qsfc_ice(1)
!         write(0,*)"sno:",snowh_lnd(1),snowh_ocn(1),snowh_ice(1)
!         write(0,*)"znt (cm):",znt_lnd(1),znt_ocn(1),znt_ice(1)
!         write(0,*)"cm:",cm_lnd(1),cm_ocn(1),cm_ice(1)
!         write(0,*)"ch:",ch_lnd(1),ch_ocn(1),ch_ice(1)
!         write(0,*)"fm:",fm_lnd(1),fm_ocn(1),fm_ice(1)
!         write(0,*)"fh:",fh_lnd(1),fh_ocn(1),fh_ice(1) 
!         write(0,*)"rb:",rb_lnd(1),rb_ocn(1),rb_ice(1)
!         write(0,*)"xland=",xland(1)," wstar:",wstar(1)
!         write(0,*)"HFX:",hfx(1)," qfx:",qfx(1)
!         write(0,*)"HFLX:",hflx(1)," evap:",evap(1)
!         write(0,*)"qsfc:",qsfc(1)," ps:",ps(1)," wspd:",wspd(1)
!         write(0,*)"ZOL:",ZOL(1)," rmol=",rmol(1)
!         write(0,*)"psim:",psim(1)," psih=",psih(1)," pblh:",pblh(1)
!         write(0,*)"FLHC=",FLHC(1)," CHS=",CHS(1)
!         write(0,*)
!      endif


  END SUBROUTINE mynnsfc_wrapper_run

!###=================================================================

END MODULE mynnsfc_wrapper

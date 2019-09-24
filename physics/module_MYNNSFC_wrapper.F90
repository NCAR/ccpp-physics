!> \file module_mynnsfc_wrapper.F90
!!  Contains all of the code related to running the MYNN surface layer scheme 

      MODULE mynnsfc_wrapper

      contains

      subroutine mynnsfc_wrapper_init ()
      end subroutine mynnsfc_wrapper_init

      subroutine mynnsfc_wrapper_finalize ()
      end subroutine mynnsfc_wrapper_finalize

!>\defgroup gsd_mynn_sfc GSD MYNN Surface Layer Scheme Module
!> \brief This scheme (1) performs pre-mynnsfc work, (2) runs the mynn sfc layer scheme, and (3) performs post-mynnsfc work
#if 0
!! \section arg_table_mynnsfc_wrapper_run Argument Table
!! \htmlinclude mynnsfc_wrapper_run.html
!!
#endif
!###===================================================================
SUBROUTINE mynnsfc_wrapper_run(         &
     &  ix,im,levs,                     &
     &  iter,flag_init,flag_restart,    &
     &  delt,dx,                        &
     &  u, v, t3d, qvsh, qc, prsl, phii,&
     &  exner, tsq, qsq, cov, sh3d,     &
     &  el_pbl, qc_bl, cldfra_bl,       &
     &  ps, PBLH, slmsk, TSK,           &
     &  QSFC, snowd,                    &
     &  zorl,UST,USTM, ZOL,MOL,RMOL,    &
     &  fm, fh, fm10, fh2, WSPD, br, ch,&
     &  HFLX, QFX, LH, FLHC, FLQC,      &
     &  U10, V10, TH2, T2, Q2,          &
     &  wstar, CHS2, CQS2,              &
     &  cda, cka, stress,               &
!     &  CP, G, ROVCP, R, XLV,           &
!     &  SVP1, SVP2, SVP3, SVPT0,        &
!     &  EP1,EP2,KARMAN,                 &
     &  icloud_bl, bl_mynn_cloudpdf,    &
     &  lprnt, errmsg, errflg           )


! should be moved to inside the mynn:
      use machine , only : kind_phys
!      use funcphys, only : fpvs

      use physcons, only : cp     => con_cp,              &
     &                     g      => con_g,               &
     &                     r_d    => con_rd,              &
     &                     r_v    => con_rv,              &
     &                     cpv    => con_cvap,            &
     &                     cliq   => con_cliq,            &
     &                     Cice   => con_csol,            &
     &                     rcp    => con_rocp,            &
     &                     XLV    => con_hvap,            &
     &                     XLF    => con_hfus,            &
     &                     EP_1   => con_fvirt,           &
     &                     EP_2   => con_eps

      USE module_sf_mynn, only : SFCLAY_mynn 

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

!-------------------------------------------------------------------
!For WRF:
!-------------------------------------------------------------------
!  USE module_model_constants, only: &
!       &karman, g, p1000mb, &
!       &cp, r_d, r_v, rcp, xlv, xlf, xls, &
!      &svp1, svp2, svp3, svpt0, ep_1, ep_2, rvovrd, &
!       &cpv, cliq, cice

!-------------------------------------------------------------------
!For reference
!   REAL    , PARAMETER :: karman       = 0.4
!   REAL    , PARAMETER :: g            = 9.81
!   REAL    , PARAMETER :: r_d          = 287.
!   REAL    , PARAMETER :: cp           = 7.*r_d/2.
!   REAL    , PARAMETER :: r_v          = 461.6
!   REAL    , PARAMETER :: cpv          = 4.*r_v
!   REAL    , PARAMETER :: cliq         = 4190.
!   REAL    , PARAMETER :: Cice         = 2106.
!   REAL    , PARAMETER :: rcp          = r_d/cp
!   REAL    , PARAMETER :: XLS          = 2.85E6
!   REAL    , PARAMETER :: XLV          = 2.5E6
!   REAL    , PARAMETER :: XLF          = 3.50E5
!   REAL    , PARAMETER :: p1000mb      = 100000.
!   REAL    , PARAMETER :: rvovrd       = r_v/r_d
!   REAL    , PARAMETER :: SVP1         = 0.6112
!   REAL    , PARAMETER :: SVP2         = 17.67
!   REAL    , PARAMETER :: SVP3         = 29.65
!   REAL    , PARAMETER :: SVPT0        = 273.15
!   REAL    , PARAMETER :: EP_1         = R_v/R_d-1.
!   REAL    , PARAMETER :: EP_2         = R_d/R_v

  REAL, PARAMETER :: xlvcp=xlv/cp, xlscp=(xlv+xlf)/cp, ev=xlv, rd=r_d, &
       &rk=cp/rd, svp11=svp1*1.e3, p608=ep_1, ep_3=1.-ep_2, g_inv=1/g


  character(len=*), intent(out) :: errmsg
  integer, intent(out) :: errflg

! NAMELIST OPTIONS (INPUT):
      INTEGER, INTENT(IN) ::                                &
     &       bl_mynn_cloudpdf,                              &
     &       icloud_bl

!MISC CONFIGURATION OPTIONS
      INTEGER, PARAMETER ::                                 &
     &       spp_pbl  = 0,                                  &
     &       isftcflx = 0,                                  &
     &       iz0tlnd  = 0,                                  &
     &       isfflx   = 1

!MYNN-1D
      REAL    :: delt
      INTEGER :: im, ix, levs
      INTEGER :: iter, k, i, itimestep
      LOGICAL :: flag_init,flag_restart,lprnt
      INTEGER :: IDS,IDE,JDS,JDE,KDS,KDE,                   &
     &            IMS,IME,JMS,JME,KMS,KME,                  &
     &            ITS,ITE,JTS,JTE,KTS,KTE

!MYNN-3D
      real(kind=kind_phys), dimension(im,levs+1) :: phii
      real(kind=kind_phys), dimension(im,levs) ::           &
     &        exner, PRSL,                                  &
     &        u, v, t3d, qvsh, qc,                          &
     &        Sh3D, EL_PBL, EXCH_H,                         &
     &        qc_bl, cldfra_bl,                             &
     &        Tsq, Qsq, Cov
     !LOCAL
      real(kind=kind_phys), dimension(im,levs) ::           &
     &        dz, rho, th, qv,                              &
     &        pattern_spp_pbl

!MYNN-2D
      real(kind=kind_phys), dimension(im) ::                &
     &        dx, pblh, slmsk, tsk, qsfc, ps,               &
     &        zorl, ust, ustm, hflx, qfx, br, wspd, snowd,  &
     &        FLHC, FLQC, U10, V10, TH2, T2, Q2,            &
     &        CHS2, CQS2, rmol, zol, mol, ch,               &
     &        fm, fh, fm10, fh2,                            &
     &        lh, cda, cka, stress, wstar
     !LOCAL
      real, dimension(im) ::                                &
     &        qcg, hfx, znt, ts, snowh, psim, psih,         &
     &        chs, ck, cd, mavail, regime, xland, GZ1OZ0       

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      if (lprnt) then
         write(0,*)"=============================================="
         write(0,*)"in mynn surface layer wrapper..."
         write(0,*)"flag_init=",flag_init
         write(0,*)"flag_restart=",flag_restart
         write(0,*)"iter=",iter
      endif

      ! If initialization is needed and mynnsfc_wrapper is called
      ! in a subcycling loop, then test for (flag_init==.T. .and. iter==1);
      ! initialization in sfclay_mynn is triggered by itimestep == 1
      ! DH* TODO: Use flag_restart to distinguish which fields need
      ! to be initialized and which are read from restart files
      if (flag_init.and.iter==1) then
          itimestep = 1
      else
          itimestep = 2
      endif

      !prep MYNN-only variables
            do k=1,levs
              do i=1,im
                 dz(i,k)=(phii(i,k+1) - phii(i,k))*g_inv
                 th(i,k)=t3d(i,k)/exner(i,k)
                 !qc(i,k)=MAX(qgrs(i,k,ntcw),0.0)
                 qv(i,k)=qvsh(i,k)/(1.0 - qvsh(i,k))
                 rho(i,k)=prsl(i,k)/(r_d*t3d(i,k)) !gt0(i,k))
                 pattern_spp_pbl(i,k)=0.0
              enddo
            enddo
            do i=1,im
                if (slmsk(i)==1. .or. slmsk(i)==2.)then !sea/land/ice mask (=0/1/2) in FV3
                  xland(i)=1.0                          !but land/water = (1/2) in SFCLAY_mynn
                else
                  xland(i)=2.0
                endif
!                ust(i) = sqrt(stress(i))
                !ch(i)=0.0
                HFX(i)=hflx(i)*rho(i,1)*cp
                !QFX(i)=evap(i)
                !wstar(i)=0.0
                qcg(i)=0.0
                snowh(i)=snowd(i)*800. !mm -> m
                znt(i)=zorl(i)*0.01    !cm -> m?
                ts(i)=tsk(i)/exner(i,1)  !theta
!                qsfc(i)=qss(i)
!                ps(i)=pgr(i)
!                wspd(i)=wind(i)
                mavail(i)=1.0  !????
            enddo

      if (lprnt) then
          write(0,*)"CALLING SFCLAY_mynn; input:"
          print*,"T:",t3d(1,1),t3d(1,2),t3d(1,3)
          print*,"TH:",th(1,1),th(1,2),th(1,3)
          print*,"rho:",rho(1,1),rho(1,2),rho(1,3)
          print*,"u:",u(1,1:3)
          !print*,"qv:",qv(1,1:3,1)
          print*,"p:",prsl(1,1)," snowh=",snowh(1)
          print*,"dz:",dz(1,1)," qsfc=",qsfc(1)
          print*,"rmol:",rmol(1)," ust:",ust(1)
          print*,"Tsk:",tsk(1)," Thetasurf:",ts(1)
          print*,"HFX:",hfx(1)," qfx",qfx(1)
          print*,"qsfc:",qsfc(1)," ps:",ps(1)
          print*,"wspd:",wspd(1),"br=",br(1)
          print*,"znt:",znt(1)," delt=",delt
          print*,"im=",im," levs=",levs
          print*,"flag_init=",flag_init !," ntcw=",ntcw!," ntk=",ntk
          print*,"flag_restart=",flag_restart !," ntcw=",ntcw!," ntk=",ntk
          print*,"iter=",iter
          !print*,"ncld=",ncld," ntrac(gq0)=",ntrac
          print*,"zlvl(1)=",dz(1,1)*0.5
          print*,"PBLH=",pblh(1)," xland=",xland(1)
       endif


        CALL SFCLAY_mynn(                                                 &
                     u3d=u,v3d=v,t3d=t3d,qv3d=qv,p3d=prsl,dz8w=dz,        &
                     CP=cp,G=g,ROVCP=rcp,R=r_d,XLV=xlv,                   &
                     PSFCPA=ps,CHS=chs,CHS2=chs2,CQS2=cqs2,               &
                     ZNT=znt,UST=ust,PBLH=pblh,MAVAIL=mavail,             &
                     ZOL=zol,MOL=mol,REGIME=regime,psim=psim,psih=psih,   &
                     psix=fm,psit=fh,psix10=fm10,psit2=fh2,               &
!                     fm=psix,fh=psit,fm10=psix10,fh2=psit2,               &
                     XLAND=xland,HFX=hfx,QFX=qfx,LH=lh,TSK=tsk,           &
                     FLHC=flhc,FLQC=flqc,QSFC=qsfc,RMOL=rmol,             &
                     U10=u10,V10=v10,TH2=th2,T2=t2,Q2=q2,SNOWH=snowh,     &
                     GZ1OZ0=GZ1OZ0,WSPD=wspd,BR=br,ISFFLX=isfflx,DX=dx,   &
                     SVP1=svp1,SVP2=svp2,SVP3=svp3,SVPT0=svpt0,           &
                     EP1=ep_1,EP2=ep_2,KARMAN=karman,                     &
                     itimestep=itimestep,ch=ch,                           &
                     th3d=th,pi3d=exner,qc3d=qc,rho3d=rho,                &
                     tsq=tsq,qsq=qsq,cov=cov,sh3d=sh3d,el_pbl=el_pbl,     &
                     qcg=qcg,wstar=wstar,                                 &
                     icloud_bl=icloud_bl,qc_bl=qc_bl,cldfra_bl=cldfra_bl, &
                     spp_pbl=spp_pbl,pattern_spp_pbl=pattern_spp_pbl,     &
                     ids=1,ide=im, jds=1,jde=1, kds=1,kde=levs,           &
                     ims=1,ime=im, jms=1,jme=1, kms=1,kme=levs,           &
                     its=1,ite=im, jts=1,jte=1, kts=1,kte=levs,           &
                     ustm=ustm, ck=ck, cka=cka, cd=cd, cda=cda,           &
                     isftcflx=isftcflx, iz0tlnd=iz0tlnd,                  &
                     bl_mynn_cloudpdf=bl_mynn_cloudpdf      )


     ! POST MYNN SURFACE LAYER (INTERSTITIAL) WORK:
        do i = 1, im
           hflx(i)=hfx(i)/(rho(i,1)*cp)
           !QFX(i)=evap(i)                                                                                                                                                                                                      
           zorl(i)=znt(i)*100.             !m -> cm
           stress(i) = ust(i)**2
        enddo


      if (lprnt) then
         print*
         print*,"finished with mynn_surface layer; output:"
         print*,"xland=",xland(1)," cda=",cda(1)
         print*,"rmol:",rmol(1)," ust:",ust(1)
         print*,"Tsk:",tsk(1)," Thetasurf:",ts(1)
         print*,"HFX:",hfx(1)," qfx",qfx(1)
         print*,"qsfc:",qsfc(1)," ps:",ps(1)
         print*,"wspd:",wspd(1)," br=",br(1)
         print*,"znt:",znt(1),"pblh:",pblh(1)
         print*,"FLHC=",FLHC(1)," CHS=",CHS(1)
         print*
      endif


  END SUBROUTINE mynnsfc_wrapper_run

!###=================================================================

END MODULE mynnsfc_wrapper

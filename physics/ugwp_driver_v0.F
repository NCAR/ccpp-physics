!!23456	
       module sso_coorde
!
! specific to COORDE-2019 project OGW switches/sensitivity
! to diagnose SSO effects pgwd=1 (OGW is on) =0 (off)
!                         pgd4=4 (4 timse taub, control pgwd=1)
!     
       use machine,      only: kind_phys
       real(kind=kind_phys),parameter :: pgwd  = 1._kind_phys
       real(kind=kind_phys),parameter :: pgwd4 = 1._kind_phys
       end module sso_coorde 
!
!
! Routine cires_ugwp_driver_v0 is replaced with cires_ugwp.F90/cires_ugwp_run in CCPP
#if 0
       subroutine cires_ugwp_driver_v0(me,  master,
     &    im,  levs, nmtvr, dtp, kdt, imx, do_ugwp, do_tofd,
     &    cdmbgwd,  xlat, xlatd, sinlat, coslat, spgrid,
     &    ugrs, vgrs, tgrs, qgrs, prsi, prsl, prslk,
     &    phii, phil, del, hprime, oc, oa4, clx, theta,
     &    gamm, sigma, elvmax, sgh30, kpbl, 
     &    dusfcg, dvsfcg, gw_dudt, gw_dvdt, gw_dtdt, gw_kdis,
     &    tau_tofd, tau_mtb, tau_ogw, tau_ngw,
     &    zmtb, zlwb, zogw, du3dt_mtb,du3dt_ogw, du3dt_tms,rdxzb,
     &    rain, ntke, tke, lprnt, ipr)
!-----------------------------------------------------------
! Part 1 "old-revised" gfs-gwdps_v0 or "old" gwdps (if do_ugwp=.false.)
! Part 2  non-stationary multi-wave GWs FV3GFS-v0
! Part 3  Dissipative version of UGWP-tendency application
!         (similar to WAM-2017)
!----------------------------------------------------------- 
       use machine,          only : kind_phys
       use physcons,         only : con_cp, con_g, con_rd, con_rv
 
       use ugwp_wmsdis_init, only : tamp_mpa, ilaunch
       use sso_coorde,       only : pgwd, pgwd4
       implicit none
!input
 
       integer, intent(in) :: me,  master
       integer, intent(in) :: im, levs, kdt, imx, nmtvr, ntke, ipr

       real(kind=kind_phys), intent(in) :: dtp, cdmbgwd(4)
       logical             :: do_ugwp, do_tofd, lprnt
       integer, intent(in) :: kpbl(im)
       real(kind=kind_phys), intent(in), dimension(im) :: xlat, xlatd
     &,                            sgh30, sinlat, coslat, spgrid    ! spgrid = tile-area
     &,                            rain

       real(kind=kind_phys), intent(in), dimension(im,levs) :: ugrs
     &,        vgrs, tgrs, qgrs, prsl, prslk, phil, del
       real(kind=kind_phys), intent(in), dimension(im,levs+1) :: 
     &         phii, prsi

!      real(kind=kind_phys), intent(in) :: oro_stat(im,nmtvr)
       real(kind=kind_phys), intent(in), dimension(im) :: hprime, oc
     &,                        theta, gamm, sigma, elvmax
       real(kind=kind_phys), intent(in), dimension(im,4) :: oa4, clx
       real(kind=kind_phys), intent(in)                  :: tke(im,levs)
!out
       real(kind=kind_phys), dimension(im,levs) :: gw_dudt, gw_dvdt
     &,                                            gw_dTdt, gw_kdis

!-----locals	+ diagnostics output

       real(kind=kind_phys), dimension(im,levs)   :: Pdvdt, Pdudt
     &,                      Pdtdt, Pkdis, ed_dudt, ed_dvdt, ed_dTdt

       real(kind=kind_phys), dimension(im)      :: dusfcg, dvsfcg

       real(kind=kind_phys), dimension(im)      :: rdxzb, zmtb,
     &    zlwb, zogw, tau_mtb, tau_ogw, tau_tofd,  tau_ngw, turb_fac
       real(kind=kind_phys), dimension(im,levs) :: du3dt_mtb, du3dt_ogw
     &,                                            du3dt_tms
       real(kind=kind_phys), dimension(im)      :: tem

! locals
       real(kind=kind_phys) :: rfac, tx1
       integer              :: i, j, k, ix
!
! define hprime, oc, oa4, clx, theta, sigma, gamm, elvmax
!
!      real(kind=kind_phys), dimension(im)     :: hprime,
!    &       oc, theta, sigma, gamm, elvmax
!      real(kind=kind_phys), dimension(im, 4)  :: clx, oa4
!	
! switches that activate impact of OGWs and NGWs along with eddy diffusion
!
       real(kind=kind_phys), parameter :: pogw=1.0, pngw=1.0, pked=1.0
     &,                                   ompked=1.0-pked
!
! switches for GW-effects: pogw=1 (OGWs)  pngw=1 (NGWs) pked=1 (eddy mixing)
!
       if (me == master .and. kdt < 2) then
         print *
         write(6,*) 'FV3GFS execute ugwp_driver_v0 '
!        write(6,*) 'FV3GFS execute ugwp_driver_v0 nmtvr=', nmtvr
         write(6,*) ' COORDE EXPER pogw = ' , pogw
         write(6,*) ' COORDE EXPER pgwd = ' , pgwd
         write(6,*) ' COORDE EXPER pgwd4 = ', pgwd4
         print *
       endif
 
       do i=1,im
         zlwb(i) = 0.
       enddo
!
! 1) ORO stationary GWs
!    ------------------
       
       if (do_ugwp .and. nmtvr == 14) then   ! calling revised old GFS gravity wave drag
         CALL GWDPS_V0(IM,    levs,  imx,   do_tofd,
     &                 Pdvdt, Pdudt, Pdtdt, Pkdis,
     &                 ugrs , vgrs, tgrs, qgrs,KPBL, prsi,del,prsl,
     &                 prslk, phii, phil, DTP,KDT,
     &                 sgh30, HPRIME, OC, OA4, CLX, THETA,
     &                 SIGMA, GAMM,   ELVMAX,
     &                 DUSFCg, DVSFCg, xlatd, sinlat, coslat, spgrid,
     &                 cdmbgwd(1:2), me, master, rdxzb,
     &                 zmtb, zogw, tau_mtb, tau_ogw, tau_tofd,
     &                 du3dt_mtb, du3dt_ogw, du3dt_tms)
!
         if (me == master .and. kdt < 2) then
           print *
           write(6,*) 'FV3GFS finished gwdps_v0 in ugwp_driver_v0 '
           print *
         endif  
       else                                  ! calling old GFS gravity wave drag as is
         do k=1,levs
           do i=1,im
             pdvdt(i,k) = 0.0
             pdudt(i,k) = 0.0
             pdtdt(i,k) = 0.0
             pkdis(i,k) = 0.0
           enddo
         enddo
         if (cdmbgwd(1) > 0.0 .or. cdmbgwd(2) > 0.0) then
           call gwdps(im, im, im, levs, Pdvdt, Pdudt, Pdtdt             &
     &,               ugrs, vgrs, tgrs, qgrs                            &
     &,               kpbl, prsi, del, prsl, prslk, phii, phil, dtp, kdt&
     &,               hprime, oc, oa4, clx, theta, sigma, gamm          &
     &,               elvmax, dusfcg, dvsfcg                            &
     &,               con_g,  con_cp, con_rd, con_rv, imx               &
     &,               nmtvr,  cdmbgwd(1:2), me, lprnt, ipr, rdxzb)
         endif

         tau_mtb   = 0.0  ; tau_ogw   = 0.0 ;  tau_tofd = 0.0
         du3dt_mtb = 0.0  ; du3dt_ogw = 0.0 ;  du3dt_tms= 0.0
       endif
!
       if (cdmbgwd(3) > 0.0) then
! 2) non-stationary GWs with GEOS-5/MERRA GW-forcing
!    ----------------------------------------------
!--------	
! GMAO GEOS-5/MERRA GW-forcing	lat-dep
!--------
         call slat_geos5_tamp(im, tamp_mpa, xlatd, tau_ngw)
 
!        call slat_geos5(im, xlatd, tau_ngw)
!
         if (abs(1.0-cdmbgwd(3)) > 1.0e-6) then
           if (cdmbgwd(4) > 0.0) then
             do i=1,im
               turb_fac(i) = 0.0
             enddo
             if (ntke > 0) then
               do k=1,(levs+levs)/3
                 do i=1,im
                   turb_fac(i) = turb_fac(i) + del(i,k) * tke(i,k)
                   tem(i)      = tem(i)      + del(i,k)
                 enddo
               enddo
               do i=1,im
                 turb_fac(i) = turb_fac(i) / tem(i)
               enddo
             endif
             rfac = 86400000 / dtp
             do i=1,im
               tx1 = cdmbgwd(4)*min(10.0, max(turb_fac(i),rain(i)*rfac))
               tau_ngw(i) = tau_ngw(i) * max(0.1, min(5.0, tx1))
             enddo
           endif
           do i=1,im
             tau_ngw(i) = tau_ngw(i) * cdmbgwd(3)
           enddo
         endif
!
         call fv3_ugwp_solv2_v0(im,   levs,  dtp,
     &                          tgrs, ugrs,  vgrs,   qgrs, prsl, prsi,
     &                          phil, xlatd, sinlat, coslat, 
     &                          gw_dudt, gw_dvdt, gw_dTdt, gw_kdis,
     &                          tau_ngw, me, master, kdt)

         if (me == master .and. kdt < 2) then
           print *
           write(6,*)'FV3GFS finished fv3_ugwp_v0 in ugwp_driver_v0 '
           write(6,*) ' non-stationary GWs with GMAO/MERRA GW-forcing '
           print *
         endif
         do k=1,levs
           do i=1,im
             gw_dtdt(i,k) = pngw*gw_dtdt(i,k) + pogw*Pdtdt(i,k)
             gw_dudt(i,k) = pngw*gw_dudt(i,k) + pogw*Pdudt(i,k)
             gw_dvdt(i,k) = pngw*gw_dvdt(i,k) + pogw*Pdvdt(i,k)
             gw_kdis(i,k) = pngw*gw_kdis(i,k) + pogw*Pkdis(i,k)
           enddo
         enddo
       else
         do k=1,levs
           do i=1,im
             gw_dtdt(i,k) = Pdtdt(i,k)
             gw_dudt(i,k) = Pdudt(i,k)
             gw_dvdt(i,k) = Pdvdt(i,k)
             gw_kdis(i,k) = Pkdis(i,k)
           enddo
         enddo
       endif
 
       if (pogw == 0.0) then
!        zmtb = 0.;  zogw =0.
         tau_mtb   = 0.0  ; tau_ogw   = 0.0 ;  tau_tofd = 0.0
         du3dt_mtb = 0.0  ; du3dt_ogw = 0.0 ;  du3dt_tms= 0.0
       endif

       return
 
!=============================================================================
! make "ugwp eddy-diffusion" update for gw_dtdt/gw_dudt/gw_dvdt by solving
! vert diffusion equations & update "Statein%tgrs, Statein%ugrs, Statein%vgrs"
!=============================================================================
!
! 3) application of "eddy"-diffusion to "smooth" UGWP-related tendencies
!------------------------------------------------------------------------------
       do k=1,levs
         do i=1,im
           ed_dudt(i,k) = 0.0 ; ed_dvdt(i,k) = 0.0 ; ed_dtdt(i,k) = 0.0
         enddo
       enddo

       call edmix_ugwp_v0(im,   levs, dtp,
     &                    tgrs, ugrs, vgrs, qgrs, del,
     &                    prsl, prsi, phil, prslk,
     &                    gw_dudt, gw_dvdt, gw_dTdt, gw_kdis,
     &                    ed_dudt, ed_dvdt, ed_dTdt,
     &                    me, master, kdt )

      do k=1,levs
        do i=1,im
          gw_dtdt(i,k) = gw_dtdt(i,k)*ompked + ed_dtdt(i,k)*pked
          gw_dvdt(i,k) = gw_dvdt(i,k)*ompked + ed_dvdt(i,k)*pked
          gw_dudt(i,k) = gw_dudt(i,k)*ompked + ed_dudt(i,k)*pked
        enddo
      enddo

      end subroutine cires_ugwp_driver_v0
#endif
!	
!=====================================================================	
!
!ugwp-v0 subroutines: GWDPS_V0 and fv3_ugwp_solv2_v0	
!  
!=====================================================================
!>\ingroup cires_ugwp_run
!> @{
!!<b>Note for the sub-grid scale orography scheme in UGWP-v0:</b> Due to degraded forecast scores of simulations with revised schemes for subgrid-scale orography effects in FV3GFS, EMC reinstalled the original gwdps-code with updated efficiency factors for the mountain blocking and OGW drag. The GFS OGW is described in the separate section (\ref GFS_GWDPS) and its “call” moved into UGWP-driver subroutine. This combination of NGW and OGW schemes was tested in the FV3GFS-L127 medium-range forecasts (15-30 days) for C96, C192, C384 and C768 resolutions and work in progress to introduce the optimal choice for the scale-aware representations of the efficiency factors that will reflect the better simulations of GW activity by FV3 dynamical core at higher horizontal resolutions. With the MERRA-2 VMF function for NGWs (\ref slat_geos5_tamp) and operational OGW drag scheme (\ref GFS_GWDPS), FV3GFS simulations can successfully forecast the recent major mid-winter sudden stratospheric warming (SSW) events of 2018-02-12 and 2018-12-31 (10-14 days before the SSW onset; Yudin et al. 2019 \cite yudin_et_al_2019). The first multi-year (2015-2018) FV3GFS simulations with UGWP-v0 also produce the equatorial QBO-like oscillations in the zonal wind and temperature anomalies.
!!
      SUBROUTINE GWDPS_V0(IM,  km,    imx, do_tofd,
     &    Pdvdt, Pdudt, Pdtdt, Pkdis, U1,V1,T1,Q1,KPBL,
     &    PRSI,DEL,PRSL,PRSLK,PHII, PHIL,DTP,KDT,
     &    sgh30, HPRIME,OC,OA4,CLX4,THETA,vSIGMA,vGAMMA,ELVMAXD,
     &    DUSFC, DVSFC,  xlatd, sinlat, coslat, sparea,
     $    cdmbgwd, me, master, rdxzb,
     &    zmtb, zogw, tau_mtb, tau_ogw, tau_tofd,
     &    dudt_mtb, dudt_ogw, dudt_tms)
!----------------------------------------
! ugwp_v0
!
! modified/revised version of gwdps.f (with bug fixes, tofd, appropriate
!   computation of kref for OGW + COORDE diagnostics
!   all constants/parameters inside cires_ugwp_initialize.F90 
!----------------------------------------    

      USE MACHINE ,      ONLY : kind_phys
      use ugwp_common ,  only : rgrav,   grav,  cpd, rd, rv, rcpd, rcpd2
     &,                         pi,      rad_to_deg, deg_to_rad, pi2
     &,                         rdi,     gor,    grcp, gocp,  fv, gr2
     &,                         bnv2min, dw2min, velmin, arad
 
      use ugwp_oro_init, only : rimin,  ric,     efmin,     efmax
     &,                         hpmax,  hpmin,   sigfaci => sigfac
     &,                         dpmin,  minwnd,  hminmt,    hncrit
     &,                         RLOLEV, GMAX,    VELEPS,    FACTOP
     &,                         FRC,    CE,      CEOFRC,    frmax, CG
     &,                         FDIR,   MDIR,    NWDIR
     &,                         cdmb,   cleff,   fcrit_gfs, fcrit_mtb
     &,                         n_tofd, ze_tofd, ztop_tofd

      use cires_ugwp_module, only : kxw,  max_kdis, max_axyz
      use sso_coorde,        only : pgwd, pgwd4
!----------------------------------------
      implicit none
      character(len=8)    :: strsolver='PSS-1986'  ! current operational solver or  'WAM-2017'
      integer, intent(in) :: im, km, imx, kdt
      integer, intent(in) :: me, master
      logical, intent(in) :: do_tofd
      real(kind=kind_phys), parameter  :: sigfac = 3, sigfacS = 0.5
      real(kind=kind_phys)             :: ztopH,zlowH,ph_blk, dz_blk
      integer, intent(in)              :: KPBL(IM)    ! Index for the PBL top layer!
      real(kind=kind_phys), intent(in) :: dtp         !  time step
      real(kind=kind_phys), intent(in) :: cdmbgwd(2)
          
      real(kind=kind_phys), intent(in), dimension(im,km) ::
     &                                   u1,  v1,   t1,    q1,
     &                                   del, prsl, prslk, phil
      real(kind=kind_phys), intent(in),dimension(im,km+1):: prsi, phii
      real(kind=kind_phys), intent(in) :: xlatd(im),sinlat(im),
     &                                    coslat(im)
      real(kind=kind_phys), intent(in) :: sparea(im)

      real(kind=kind_phys), intent(in) :: OC(IM), OA4(im,4), CLX4(im,4)
      real(kind=kind_phys), intent(in) :: HPRIME(IM),  sgh30(IM)
      real(kind=kind_phys), intent(in) :: ELVMAXD(IM), THETA(IM)
      real(kind=kind_phys), intent(in) :: vSIGMA(IM),  vGAMMA(IM)
      real(kind=kind_phys)             :: SIGMA(IM),   GAMMA(IM)
      
!output -phys-tend
      real(kind=kind_phys),dimension(im,km),intent(out) ::
     &                      Pdvdt,    Pdudt,    Pkdis, Pdtdt
! output - diag-coorde
     &,                     dudt_mtb, dudt_ogw, dudt_tms
!                     
      real(kind=kind_phys),dimension(im) :: RDXZB,   zmtb,    zogw
     &,                                     tau_ogw, tau_mtb, tau_tofd
     &,                                     dusfc,   dvsfc
!
!---------------------------------------------------------------------
! # of permissible sub-grid orography hills for "any" resolution  < 25
!    correction for "elliptical" hills based on shilmin-area =sgrid/25 
!     4.*gamma*b_ell*b_ell  >=  shilmin
!     give us limits on [b_ell & gamma *b_ell] > 5 km =sso_min
!     gamma_min = 1/4*shilmin/sso_min/sso_min
!23.01.2019:  cdmb = 4.*192/768_c192=1 x 0.5
!     192: cdmbgwd        = 0.5, 2.5
!     cleff = 2.5*0.5e-5 * sqrt(192./768.) => Lh_eff = 1004. km
!      6*dx = 240 km 8*dx = 320. ~ 3-5 more effective
!---------------------------------------------------------------------
      real(kind=kind_phys)            :: gammin = 0.00999999
      real(kind=kind_phys), parameter :: nhilmax = 25.
      real(kind=kind_phys), parameter :: sso_min = 3000.
      logical, parameter              :: do_adjoro = .true.
!
      real(kind=kind_phys)            :: shilmin, sgrmax, sgrmin
      real(kind=kind_phys)            :: belpmin, dsmin,  dsmax
!     real(kind=kind_phys)            :: arhills(im)              ! not used why do we need?
      real(kind=kind_phys)            :: xlingfs
         
! 
! locals       
! mean flow
      real(kind=kind_phys), dimension(im,km) :: RI_N, BNV2, RO
     &,                                         VTK, VTJ, VELCO
!mtb     
      real(kind=kind_phys), dimension(im)    :: OA,  CLX , elvmax, wk
     &,                                         PE, EK, UP
      
      real(kind=kind_phys), dimension(im,km) :: DB, ANG, UDS

      real(kind=kind_phys) :: ZLEN, DBTMP, R, PHIANG, DBIM, ZR
      real(kind=kind_phys) :: ENG0, ENG1, COSANG2, SINANG2
      real(kind=kind_phys) :: bgam, cgam, gam2, rnom, rdem       
!
! TOFD
!     Some constants now in "use ugwp_oro_init" +   "use ugwp_common"
!
!==================
      real(kind=kind_phys)   :: unew, vnew,  zpbl,  sigflt, zsurf
      real(kind=kind_phys), dimension(km)    :: utofd1, vtofd1
     &,                                         epstofd1, krf_tofd1
     &,                                         up1, vp1, zpm
      real(kind=kind_phys),dimension(im, km) :: axtms, aytms
! 
! OGW
!
      LOGICAL ICRILV(IM)
!
      real(kind=kind_phys), dimension(im) :: XN, YN, UBAR, VBAR, ULOW,
     &               ROLL,  bnv2bar, SCOR, DTFAC, XLINV, DELKS, DELKS1
!
      real(kind=kind_phys) :: TAUP(IM,km+1), TAUD(IM,km)
      real(kind=kind_phys) :: taub(im), taulin(im), heff, hsat, hdis

      integer, dimension(im) :: kref, idxzb, ipt, kreflm, 
     &                          iwklm, iwk, izlow
!      
!check what we need
!
      real(kind=kind_phys) :: bnv,  fr, ri_gw
     &,                       brvf,   tem,   tem1,  tem2, temc, temv
     &,                       ti,     rdz,   dw2,   shr2, bvf2
     &,                       rdelks, efact, coefm, gfobnv
     &,                       scork,  rscor, hd,    fro,  sira
     &,                       dtaux,  dtauy, pkp1log, pklog
     &,                       grav2, rcpdt, windik, wdir
     &,                       sigmin, dxres,sigres,hdxres
     &,                       cdmb4, mtbridge
     &,                       kxridge, inv_b2eff, zw1, zw2
     &,                       belps, aelps, nhills, selps
     
      integer ::          kmm1, kmm2, lcap, lcapp1
     &,            npt,   kbps, kbpsp1,kbpsm1
     &,            kmps,  idir, nwd,  klcap, kp1, kmpbl, kmll
     &,            k_mtb, k_zlow, ktrial, klevm1, i, j, k
! 
      rcpdt = 1.0 / (cpd*dtp)
      grav2 = grav + grav
!       
! mtb-blocking  sigma_min and dxres => cires_initialize
!  
      sgrmax = maxval(sparea) ; sgrmin = minval(sparea)
      dsmax  = sqrt(sgrmax)   ; dsmin  = sqrt(sgrmin)

      dxres   = pi2*arad/float(IMX)
      hdxres  = 0.5*dxres
!     shilmin = sgrmin/nhilmax            ! not used - Moorthi

!     gammin = min(sso_min/dsmax, 1.)     ! Moorthi - with this results are not reproducible
      gammin = min(sso_min/dxres, 1.)     ! Moorthi

!     sigmin = 2.*hpmin/dsmax      !dxres ! Moorthi - this will not reproduce
      sigmin = 2.*hpmin/dxres      !dxres

!     if (kdt == 1) then
!       print *, sgrmax, sgrmin , ' min-max sparea '
!       print *, 'sigmin-hpmin-dsmax', sigmin, hpmin, dsmax
!       print *, 'dxres/dsmax ', dxres, dsmax
!       print *, ' shilmin gammin ', shilmin, gammin
!     endif

      kxridge = float(IMX)/arad * cdmbgwd(2)

      if (me == master .and. kdt == 1) then
        print *, ' gwdps_v0 kxridge ', kxridge
        print *, ' gwdps_v0 scale2 ', cdmbgwd(2)
        print *, ' gwdps_v0 IMX ', imx
        print *, ' gwdps_v0 GAM_MIN ', gammin
        print *, ' gwdps_v0 SSO_MIN ', sso_min
      endif

      do i=1,im
        idxzb(i)    = 0
        zmtb(i)     = 0.0
        zogw(i)     = 0.0
        rdxzb(i)    = 0.0      
        tau_ogw(i)  = 0.0
        tau_mtb(i)  = 0.0 
        dusfc(i)    = 0.0
        dvsfc(i)    = 0.0
        tau_tofd(i) = 0.0
!
        ipt(i) = 0
        sigma(i) = max(vsigma(i), sigmin)
        gamma(i) = max(vgamma(i), gammin)
      enddo

      do k=1,km
        do i=1,im
          pdvdt(i,k)    = 0.0
          pdudt(i,k)    = 0.0
          pdtdt(i,k)    = 0.0
          pkdis(i,k)    = 0.0
          dudt_mtb(i,k) = 0.0
          dudt_ogw(i,k) = 0.0
          dudt_tms(i,k) = 0.0
        enddo
      enddo
 
! ----  for lm and gwd calculation points
      
      npt = 0
      do i = 1,im
        if ( elvmaxd(i) >= hminmt .and. hprime(i)  >= hpmin ) then
          
          npt      = npt + 1
          ipt(npt) = i
!         arhills(i) = 1.0
!
          sigres = max(sigmin, sigma(i))
!         if (sigma(i) < sigmin) sigma(i)=  sigmin
          dxres = sqrt(sparea(i))
          if (2.*hprime(i)/sigres > dxres) sigres=2.*hprime(i)/dxres
          aelps = min(2.*hprime(i)/sigres, 0.5*dxres)
          if (gamma(i) > 0.0 ) belps = min(aelps/gamma(i),.5*dxres)
!
! small-scale "turbulent" oro-scales < sso_min
!
          if( aelps < sso_min .and. do_adjoro) then
 
! a, b > sso_min upscale ellipse  a/b > 0.1 a>sso_min & h/b=>new_sigm
!
            aelps = sso_min 
            if (belps < sso_min ) then
              gamma(i) = 1.0
              belps = aelps*gamma(i)
            else
              gamma(i) = min(aelps/belps, 1.0)
            endif
            sigma(i) = 2.*hprime(i)/aelps
            gamma(i) = min(aelps/belps, 1.0)
          endif
 
          selps      = belps*belps*gamma(i)*4.    ! ellipse area of the el-c hill 
          nhills     = min(nhilmax, sparea(i)/selps)
!         arhills(i) = max(nhills, 1.0)

!333   format( ' nhil: ', I6, 4(2x, F9.3), 2(2x, E9.3))	    
!      if (kdt==1 )
!     & write(6,333) nint(nhills)+1,xlatd(i), hprime(i),aelps*1.e-3,
!     &   belps*1.e-3, sigma(i),gamma(i)

        endif
      enddo
        
      IF (npt == 0) then
!         print *,  'oro-npt = 0 elvmax ', maxval(elvmaxd), hminmt
!         print *,  'oro-npt = 0 hprime ', maxval(hprime), hpmin	     	    
        RETURN      ! No gwd/mb calculation done
      endif


      do i=1,npt
        iwklm(i)  = 2
        IDXZB(i)  = 0 
        kreflm(i) = 0
      enddo
 
      do k=1,km
        do i=1,im
          db(i,k)  = 0.0
          ang(i,k) = 0.0
          uds(i,k) = 0.0 
        enddo
      enddo

      KMM1 = km - 1 ;  KMM2   = km - 2 ; KMLL   = kmm1
      LCAP = km     ;  LCAPP1 = LCAP + 1 
 
      DO I = 1, npt
        j = ipt(i)
        ELVMAX(J) = min (ELVMAXd(J)*0. + sigfac * hprime(j), hncrit)
        izlow(i)  = 1          ! surface-level
      ENDDO
!
      DO K = 1, kmm1
        DO I = 1, npt
          j = ipt(i)
          ztopH   = sigfac * hprime(j)
          zlowH   = sigfacs* hprime(j) 
          pkp1log =  phil(j,k+1) * rgrav
          pklog   =  phil(j,k)   * rgrav
!         if (( ELVMAX(j) <= pkp1log) .and. (ELVMAX(j).ge.pklog) )
!     &      iwklm(I)  =  MAX(iwklm(I), k+1 ) 
          if (( ztopH <= pkp1log) .and. (zTOPH >= pklog) )
     &        iwklm(I)  =  MAX(iwklm(I), k+1 )
!
          if (zlowH <= pkp1log .and. zlowH >= pklog)
     &        izlow(I)  =  MAX(izlow(I),k)
        ENDDO
      ENDDO
!
      DO K = 1,km
        DO I =1,npt
          J         = ipt(i)
          VTJ(I,K)  = T1(J,K)  * (1.+FV*Q1(J,K))
          VTK(I,K)  = VTJ(I,K) / PRSLK(J,K)
          RO(I,K)   = RDI * PRSL(J,K) / VTJ(I,K)       ! DENSITY mid-levels
          TAUP(I,K) = 0.0
        ENDDO
      ENDDO
!
! check RI_N or RI_MF computation
!
      DO K = 1,kmm1
        DO I =1,npt
          J         = ipt(i)
          RDZ       = grav   / (phil(j,k+1) - phil(j,k))
          TEM1      = U1(J,K) - U1(J,K+1)
          TEM2      = V1(J,K) - V1(J,K+1)
          DW2       = TEM1*TEM1 + TEM2*TEM2
          SHR2      = MAX(DW2,DW2MIN) * RDZ * RDZ
!          TI        = 2.0 / (T1(J,K)+T1(J,K+1))
!          BVF2      = Grav*(GOCP+RDZ*(VTJ(I,K+1)-VTJ(I,K)))* TI
!          RI_N(I,K) = MAX(BVF2/SHR2,RIMIN)   ! Richardson number
!
          BVF2 = grav2 * RDZ * (VTK(I,K+1)-VTK(I,K))
     &                       / (VTK(I,K+1)+VTK(I,K))
          bnv2(i,k+1) = max( BVF2, bnv2min )
          RI_N(I,K+1) = Bnv2(i,k)/SHR2        ! Richardson number consistent with BNV2	
!
! add here computation for Ktur and OGW-dissipation fro VE-GFS
!	    
        ENDDO
      ENDDO
      K = 1
      DO I = 1, npt
        bnv2(i,k) = bnv2(i,k+1)
      ENDDO
!		
! level iwklm =>phil(j,k)/g < sigfac * hprime(j) < phil(j,k+1)/g 
!
      DO I = 1, npt
        J   = ipt(i)
        k_zlow = izlow(I)
        if (k_zlow == iwklm(i)) k_zlow = 1
        DELKS(I)   = 1.0 / (PRSI(J,k_zlow) - PRSI(J,iwklm(i)))
!       DELKS1(I)  = 1.0 /(PRSL(J,k_zlow) - PRSL(J,iwklm(i)))
        UBAR (I)   = 0.0
        VBAR (I)   = 0.0
        ROLL (I)   = 0.0
        PE   (I)   = 0.0
        EK   (I)   = 0.0
        BNV2bar(I) = 0.0   
      ENDDO
!
      DO I = 1, npt
        k_zlow = izlow(I)
        if (k_zlow == iwklm(i)) k_zlow = 1
        DO K = k_zlow, iwklm(I)-1                        ! Kreflm(I)= iwklm(I)-1 
          J       = ipt(i)                               ! laye-aver Rho, U, V
          RDELKS  = DEL(J,K) * DELKS(I)
          UBAR(I) = UBAR(I)  + RDELKS * U1(J,K)          ! trial Mean U below 
          VBAR(I) = VBAR(I)  + RDELKS * V1(J,K)          ! trial Mean V below 
          ROLL(I) = ROLL(I)  + RDELKS * RO(I,K)          ! trial Mean RO below 
!   
          BNV2bar(I) = BNV2bar(I) + .5*(BNV2(I,K)+BNV2(I,K+1))* RDELKS
        ENDDO
      ENDDO
!
      DO I = 1, npt
        J = ipt(i)
!
! integrate from Ztoph = sigfac*hprime  down to Zblk if exists
! find ph_blk, dz_blk like in LM-97 and IFS
!	
        ph_blk =0.  
        DO K = iwklm(I), 1, -1
          PHIANG   =  atan2(V1(J,K),U1(J,K))*RAD_TO_DEG
          ANG(I,K) = ( THETA(J) - PHIANG )
          if ( ANG(I,K) >  90. ) ANG(I,K) = ANG(I,K) - 180.
          if ( ANG(I,K) < -90. ) ANG(I,K) = ANG(I,K) + 180.
          ANG(I,K) = ANG(I,K) * DEG_TO_RAD
          UDS(I,K) = 
     &        MAX(SQRT(U1(J,K)*U1(J,K) + V1(J,K)*V1(J,K)), velmin)
!
          IF (IDXZB(I) == 0 ) then
            dz_blk = ( PHII(J,K+1) - PHII(J,K) ) *rgrav
            PE(I)  =  PE(I) + BNV2(I,K) * 
     &         ( ELVMAX(J) - phil(J,K)*rgrav ) * dz_blk

            UP(I)  =  max(UDS(I,K) * cos(ANG(I,K)), velmin)  
            EK(I)  = 0.5 *  UP(I) * UP(I) 

            ph_blk = ph_blk + dz_blk*sqrt(BNV2(I,K))/UP(I)

! --- Dividing Stream lime  is found when PE =exceeds EK. oper-l GFS
!           IF ( PE(I) >=  EK(I) ) THEN
            IF ( ph_blk >=  fcrit_gfs ) THEN
               IDXZB(I) = K
               zmtb (J) = PHIL(J, K)*rgrav
               RDXZB(J) = real(k, kind=kind_phys)
            ENDIF

          ENDIF
        ENDDO
!
! Alternative expression: ZMTB = max(Heff*(1. -Fcrit_gfs/Fr), 0)
! fcrit_gfs/fr	 
!
        goto 788

        BNV     = SQRT( BNV2bar(I) )
        heff    = 2.*min(HPRIME(J),hpmax)
        zw2     = UBAR(I)*UBAR(I)+VBAR(I)*VBAR(I)
        Ulow(i) = sqrt(max(zw2,dw2min))
        Fr      = heff*bnv/Ulow(i)
        ZW1     = max(Heff*(1. -fcrit_gfs/fr), 0.0)
        zw2     = phil(j,2)*rgrav
        if (Fr > fcrit_gfs .and. zw1 > zw2 ) then 
          do k=2, kmm1
            pkp1log =  phil(j,k+1) * rgrav
            pklog   =  phil(j,k)   * rgrav
            if (zw1 <= pkp1log .and. zw1 >= pklog)  exit
          enddo
            IDXZB(I) = K
            zmtb (J) = PHIL(J, K)*rgrav
        else
           zmtb (J) = 0.
           IDXZB(I) = 0
        endif
788     continue
      ENDDO

!
! --- The drag for mtn blocked flow
! 
      cdmb4 = 0.25*cdmb 
      DO I = 1, npt
        J = ipt(i)
!
        IF ( IDXZB(I) > 0 ) then
! (4.16)-IFS	  
          gam2 = gamma(j)*gamma(j)
          BGAM = 1.0 - 0.18*gamma(j) - 0.04*gam2
          CGAM =       0.48*gamma(j) + 0.30*gam2
          DO K = IDXZB(I)-1, 1, -1

            ZLEN = SQRT( ( PHIL(J,IDXZB(I)) - PHIL(J,K) ) / 
     &                   ( PHIL(J,K ) + Grav * hprime(J) ) )

            tem     = cos(ANG(I,K))
            COSANG2 = tem * tem
            SINANG2 = 1.0 - COSANG2 
!	      
!  cos =1 sin =0 =>   1/R= gam     ZR = 2.-gam 
!  cos =0 sin =1 =>   1/R= 1/gam   ZR = 2.- 1/gam
!
            rdem = COSANG2      +  GAM2 * SINANG2
            rnom = COSANG2*GAM2 +         SINANG2
!	       
! metOffice Dec 2010
! correction of H. Wells & A. Zadra for the
! aspect ratio  of the hill seen by MF
! (1/R , R-inverse below: 2-R)

            rdem = max(rdem, 1.e-6)       
            R    = sqrt(rnom/rdem)
            ZR   =  MAX( 2. - R, 0. )

            sigres = max(sigmin, sigma(J))
            if (hprime(J)/sigres > dxres) sigres = hprime(J)/dxres
            mtbridge = ZR * sigres*ZLEN / hprime(J)
! (4.15)-IFS 	   
!           DBTMP = CDmb4 * mtbridge *
!     &             MAX(cos(ANG(I,K)), gamma(J)*sin(ANG(I,K)))
! (4.16)-IFS
            DBTMP  = CDmb4*mtbridge*(bgam* COSANG2 +cgam* SINANG2)
            DB(I,K)= DBTMP * UDS(I,K)
          ENDDO
!                  
        endif
      ENDDO
! 
!.............................
!.............................
! end  mtn blocking section
!.............................
!.............................
!
!--- Orographic Gravity Wave Drag Section
!     
!  Scale cleff between IM=384*2 and 192*2 for T126/T170 and T62
!  inside "cires_ugwp_initialize.F90" now
!
      KMPBL  = km / 2 
      iwk(1:npt) = 2
!
! METO-scheme: 
! k_mtb = max(k_zmtb, k_n*hprime/2] to reduce diurnal variations taub_ogw 
!     
      DO K=3,KMPBL
        DO I=1,npt
          j   = ipt(i)
          tem = (prsi(j,1) - prsi(j,k))
          if (tem < dpmin) iwk(i) = k    ! dpmin=50 mb
  
!===============================================================	  
! lev=111      t=311.749     hkm=0.430522     Ps-P(iwk)=52.8958 
!           below "Hprime" - source of OGWs  and below Zblk !!!
!           27           2  kpbl ~ 1-2 km   < Hprime
!===============================================================	  
        enddo
      enddo
!
! iwk - adhoc GFS-parameter to select OGW-launch level between
!      LEVEL ~0.4-0.5 KM from surface or/and  PBL-top
! in UGWP-V1: options to modify as  Htop ~ (2-3)*Hprime > Zmtb
! in UGWP-V0 we ensured that : Zogw > Zmtb
!

      KBPS  = 1
      KMPS  = km
      K_mtb = 1
      DO I=1,npt
        J         = ipt(i)
        K_mtb     = max(1, idxzb(i))

        kref(I)   = MAX(IWK(I), KPBL(J)+1 )             ! reference level PBL or smt-else ????
        kref(I)   = MAX(kref(i), iwklm(i) )             ! iwklm => sigfac*hprime

        if (kref(i) <= idxzb(i)) kref(i) = idxzb(i) + 1 ! layer above zmtb
        KBPS      = MAX(KBPS,  kref(I))
        KMPS      = MIN(KMPS,  kref(I))
!
        DELKS(I)  = 1.0 / (PRSI(J,k_mtb) - PRSI(J,kref(I)))
        UBAR (I)  = 0.0
        VBAR (I)  = 0.0
        ROLL (I)  = 0.0
        BNV2bar(I)= 0.0
      ENDDO
!
      KBPSP1 = KBPS + 1
      KBPSM1 = KBPS - 1
      K_mtb  = 1
!
      DO I = 1,npt
        K_mtb = max(1, idxzb(i))
        DO K = k_mtb,KBPS            !KBPS = MAX(kref) ;KMPS= MIN(kref)
          IF (K < kref(I)) THEN
            J          = ipt(i)
            RDELKS     = DEL(J,K) * DELKS(I)
            UBAR(I)    = UBAR(I)  + RDELKS * U1(J,K)   ! Mean U below kref
            VBAR(I)    = VBAR(I)  + RDELKS * V1(J,K)   ! Mean V below kref
            ROLL(I)    = ROLL(I)  + RDELKS * RO(I,K)   ! Mean RO below kref
            BNV2bar(I) = BNV2bar(I) + .5*(BNV2(I,K)+BNV2(I,K+1))* RDELKS
          ENDIF
        ENDDO
      ENDDO
!
! orographic asymmetry parameter (OA), and (CLX)
      DO I = 1,npt
        J      = ipt(i)
        wdir   = atan2(UBAR(I),VBAR(I)) + pi
        idir   = mod(nint(fdir*wdir),mdir) + 1
        nwd    = nwdir(idir)
        OA(I)  = (1-2*INT( (NWD-1)/4 )) * OA4(J,MOD(NWD-1,4)+1)
        CLX(I) = CLX4(J,MOD(NWD-1,4)+1)
      ENDDO
!
      DO I = 1,npt
       DTFAC(I)  = 1.0
       ICRILV(I) = .FALSE.                      ! INITIALIZE CRITICAL LEVEL CONTROL VECTOR
       ULOW(I) = MAX(SQRT(UBAR(I)*UBAR(I)+VBAR(I)*VBAR(I)),velmin)
       XN(I)  = UBAR(I) / ULOW(I)
       YN(I)  = VBAR(I) / ULOW(I)
      ENDDO
!
      DO  K = 1, kmm1
        DO  I = 1,npt
          J            = ipt(i)
          VELCO(I,K)   = 0.5 * ((U1(J,K)+U1(J,K+1))*XN(I)
     &                       +  (V1(J,K)+V1(J,K+1))*YN(I))
        ENDDO
      ENDDO
!
!------------------
! v0: incorporates latest modifications for kxridge and heff/hsat
!             and taulin for Fr <=fcrit_gfs 
!             and concept of "clipped" hill if zmtb > 0. to make
! the integrated "tau_sso = tau_ogw +tau_mtb" close to reanalysis data
!      it is still used the "single-OGWave"-approach along ULOW-upwind
!
!      in contrast to the 2-orthogonal wave (2OTW) schemes of IFS/METO/E-CANADA  
! 2OTW scheme requires "aver angle" and wind projections on 2 axes of ellipse a-b
!     with 2-stresses:  taub_a & taub_b from AS of Phillips et al. (1984)
!------------------
      taub(:)  = 0. ; taulin(:)= 0.
      DO I = 1,npt
        J    = ipt(i)
        BNV  = SQRT( BNV2bar(I) )
        heff = min(HPRIME(J),hpmax)

        if( zmtb(j) > 0.) heff = max(sigfac*heff-zmtb(j), 0.)/sigfac  
        if (heff <= 0) cycle

        hsat = fcrit_gfs*ULOW(I)/bnv
        heff = min(heff, hsat)

        FR   = min(BNV  * heff /ULOW(I), FRMAX)
!
        EFACT    = (OA(I) + 2.) ** (CEOFRC*FR)
        EFACT    = MIN( MAX(EFACT,EFMIN), EFMAX )
!
        COEFM    = (1. + CLX(I)) ** (OA(I)+1.)
!
        XLINV(I) = COEFM * CLEFF           ! effective kxw for Lin-wave
        XLINGFS  = COEFM * CLEFF 
!
        TEM      = FR    * FR * OC(J)
        GFOBNV   = GMAX  * TEM / ((TEM + CG)*BNV) 
!
!new specification of XLINV(I) & taulin(i)

        sigres = max(sigmin, sigma(J))
        if (heff/sigres > hdxres) sigres = heff/hdxres
        inv_b2eff =  0.5*sigres/heff
        kxridge   =  1.0 / sqrt(sparea(J))
        XLINV(I)  = XLINGFS    !or max(kxridge, inv_b2eff)  ! 6.28/Lx ..0.5*sigma(j)/heff = 1./Lridge  
        taulin(i) = 0.5*ROLL(I)*XLINV(I)*BNV*ULOW(I)*
     &                 heff*heff*pgwd4

        if ( FR > fcrit_gfs ) then
          TAUB(I)  = XLINV(I) * ROLL(I) * ULOW(I) * ULOW(I)
     &             * ULOW(I)  * GFOBNV  * EFACT *pgwd4      ! nonlinear FLUX Tau0...XLINV(I)
!
        else
!	 
          TAUB(I)  = XLINV(I) * ROLL(I) * ULOW(I) * ULOW(I)
     &             * ULOW(I)  * GFOBNV  * EFACT *pgwd4  
!     
!         TAUB(I)  = taulin(i)                             !  linear flux for FR <= fcrit_gfs
!   
        endif
!		 
!
        K       = MAX(1, kref(I)-1)
        TEM     = MAX(VELCO(I,K)*VELCO(I,K), dw2min)
        SCOR(I) = BNV2(I,K) / TEM  ! Scorer parameter below ref level
!
! diagnostics for zogw > zmtb
!
        zogw(J) = PHII(j, kref(I)) *rgrav
      ENDDO
!                                                                       
!----SET UP BOTTOM VALUES OF STRESS
!
      DO K = 1, KBPS
        DO I = 1,npt
          IF (K <= kref(I)) TAUP(I,K) = TAUB(I)
        ENDDO
      ENDDO
      
      if (strsolver == 'PSS-1986') then     
       
!======================================================
!   V0-GFS OROGW-solver of Palmer et al 1986 -"PSS-1986"
!   in V1-OROGW  LINSATDIS of                 "WAM-2017"
!     with LLWB-mechanism for
!     rotational/non-hydrostat OGWs important for 
!     HighRES-FV3GFS with dx < 10 km
!======================================================
      
        DO K = KMPS, KMM1                   ! Vertical Level Loop
          KP1 = K + 1
          DO I = 1, npt
!
            IF (K >= kref(I)) THEN
              ICRILV(I) = ICRILV(I) .OR. ( RI_N(I,K) < RIC)
     &                              .OR. (VELCO(I,K) <= 0.0)
            ENDIF
          ENDDO
!
          DO I = 1,npt
            IF (K >= kref(I))   THEN
              IF (.NOT.ICRILV(I) .AND. TAUP(I,K) > 0.0 ) THEN
                TEMV = 1.0 / max(VELCO(I,K), velmin)
!
                IF (OA(I) > 0. .AND. kp1 < kref(i)) THEN
                  SCORK   = BNV2(I,K) * TEMV * TEMV
                  RSCOR   = MIN(1.0, SCORK / SCOR(I))
                  SCOR(I) = SCORK
                ELSE 
                  RSCOR   = 1.
                ENDIF
!
                BRVF = SQRT(BNV2(I,K))        ! Brent-Vaisala Frequency interface
!               TEM1 = XLINV(I)*(RO(I,KP1)+RO(I,K))*BRVF*VELCO(I,K)*0.5

                TEM1 = XLINV(I)*(RO(I,KP1)+RO(I,K))*BRVF*0.5
     &                         * max(VELCO(I,K), velmin)
                HD   = SQRT(TAUP(I,K) / TEM1)
                FRO  = BRVF * HD * TEMV
!
!    RIM is the  "WAVE"-RICHARDSON NUMBER BY PALMER,Shutts, Swinbank 1986
!

                TEM2   = SQRT(ri_n(I,K))
                TEM    = 1. + TEM2 * FRO
                RI_GW  = ri_n(I,K) * (1.0-FRO) / (TEM * TEM)
!
!    CHECK STABILITY TO EMPLOY THE 'dynamical SATURATION HYPOTHESIS'
!    OF PALMER,Shutts, Swinbank 1986
!                                       ----------------------
                IF (RI_GW <= RIC .AND.
     &             (OA(I) <= 0. .OR.  kp1 >= kref(i) )) THEN
                   TEMC = 2.0 + 1.0 / TEM2
                   HD   = VELCO(I,K) * (2.*SQRT(TEMC)-TEMC) / BRVF
                   TAUP(I,KP1) = TEM1 * HD * HD
                ELSE 
                   TAUP(I,KP1) = TAUP(I,K) * RSCOR
                ENDIF
                taup(i,kp1) = min(taup(i,kp1), taup(i,k))
              ENDIF
            ENDIF
          ENDDO
        ENDDO
!     
!  zero momentum deposition at the top model layer
!      
        taup(1:npt,km+1) = taup(1:npt,km)      
!
!     Calculate wave acc-n: - (grav)*d(tau)/d(p) = taud
!
        DO K = 1,KM
          DO I = 1,npt
            TAUD(I,K) = GRAV*(TAUP(I,K+1) - TAUP(I,K))/DEL(ipt(I),K)
          ENDDO
        ENDDO
!
!------scale MOMENTUM DEPOSITION  AT TOP TO 1/2 VALUE
! it is zero now      
!       DO I = 1,npt
!         TAUD(I, km) = TAUD(I,km) * FACTOP
!       ENDDO
    
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!------IF THE GRAVITY WAVE DRAG WOULD FORCE A CRITICAL LINE IN THE
!------LAYERS BELOW SIGMA=RLOLEV DURING THE NEXT DELTIM TIMESTEP,
!------THEN ONLY APPLY DRAG UNTIL THAT CRITICAL LINE IS REACHED.
! Empirical implementation of the LLWB-mechanism: Lower Level Wave Breaking
! by limiting "Ax = Dtfac*Ax" due to possible LLWB around Kref and 500 mb
! critical line [V - Ax*dtp = 0.] is smt like "LLWB" for stationary OGWs
!2019:  this option limits sensitivity of taux/tauy to  increase/decreaseof TAUB
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        DO K = 1,KMM1
          DO I = 1,npt
           IF (K  >= kref(I) .and. PRSI(ipt(i),K) >= RLOLEV) THEN

               IF(TAUD(I,K) /= 0.) THEN
                 TEM = DTP * TAUD(I,K)
                 DTFAC(I) = MIN(DTFAC(I),ABS(VELCO(I,K)/TEM))
!	         DTFAC(I) = 1.0
               ENDIF
            ENDIF
          ENDDO
        ENDDO
!
!--------------------------- OROGW-solver of GFS PSS-1986
!  
      else 
!
!--------------------------- OROGW-solver of WAM2017
!
!       sigres = max(sigmin, sigma(J))
!	if (heff/sigres.gt.dxres) sigres=heff/dxres
!         inv_b2eff =  0.5*sigres/heff 	
!       XLINV(I) = max(kxridge, inv_b2eff)           ! 0.5*sigma(j)/heff = 1./Lridge      
        dtfac(:) =  1.0
       
        call oro_wam_2017(im, km, npt, ipt, kref, kdt, me, master,
     &       dtp, dxres, taub, u1, v1, t1, xn, yn, bnv2, ro, prsi,prsL, 
     &       del, sigma, hprime, gamma, theta,
     &       sinlat, xlatd, taup, taud, pkdis)
     
      endif            !  oro_wam_2017 - LINSATDIS-solver of WAM-2017
!      
!--------------------------- OROGW-solver of WAM2017
! 
! TOFD as in BELJAARS-2004
!
! ---------------------------    
      IF( do_tofd ) then
        axtms(:,:) = 0.0 ; aytms(:,:) = 0.0 
        if ( kdt == 1 .and. me == 0) then
          print *, 'VAY do_tofd  from surface to ', ztop_tofd 
        endif
        DO I = 1,npt          
          J = ipt(i)
          zpbl  =rgrav*phil( j, kpbl(j) )
 
          sigflt = min(sgh30(j), 0.3*hprime(j)) ! cannot exceed 30% of LS-SSO
 
          zsurf = phii(j,1)*rgrav
          do k=1,km
            zpm(k) = phiL(j,k)*rgrav
            up1(k) = u1(j,k)
            vp1(k) = v1(j,k)
          enddo
 
          call ugwp_tofd1d(km, sigflt, elvmaxd(j), zsurf, zpbl, 
     &             up1, vp1, zpm,  utofd1, vtofd1, epstofd1, krf_tofd1)
     
          do k=1,km
            axtms(j,k) = utofd1(k)
            aytms(j,k) = vtofd1(k)
!	 
! add TOFD to GW-tendencies
!	 
            pdvdt(J,k)  = pdvdt(J,k) + aytms(j,k)
            pdudt(J,k)  = pdudt(J,k) + axtms(j,k)
          enddo
!2018-diag
          tau_tofd(J) = sum( utofd1(1:km)* del(j,1:km))
        enddo
      ENDIF             ! do_tofd 

!---------------------------
! combine oro-drag effects
!---------------------------  
! +  diag-3d

      dudt_tms = axtms 
      tau_ogw  = 0.
      tau_mtb  = 0.
 
      DO K = 1,KM
        DO I = 1,npt
          J    = ipt(i)
!
          ENG0 = 0.5*(U1(j,K)*U1(j,K)+V1(J,K)*V1(J,K))
!
          if ( K < IDXZB(I) .AND. IDXZB(I) /= 0 ) then
!
! if blocking layers -- no OGWs
!	  
            DBIM       = DB(I,K) / (1.+DB(I,K)*DTP)
            Pdvdt(j,k) = - DBIM * V1(J,K) +Pdvdt(j,k)
            Pdudt(j,k) = - DBIM * U1(J,K) +Pdudt(j,k)
            ENG1       = ENG0*(1.0-DBIM*DTP)*(1.-DBIM*DTP)
            
            DUSFC(J)   = DUSFC(J) - DBIM * U1(J,K) * DEL(J,K)
            DVSFC(J)   = DVSFC(J) - DBIM * V1(J,K) * DEL(J,K)
!2018-diag 
            dudt_mtb(j,k) = -DBIM * U1(J,K)
            tau_mtb(j)    = tau_mtb(j) + dudt_mtb(j,k)* DEL(J,K)

          else
!
! OGW-s above blocking height
!	  
            TAUD(I,K)  = TAUD(I,K) * DTFAC(I)
            DTAUX      = TAUD(I,K) * XN(I) * pgwd
            DTAUY      = TAUD(I,K) * YN(I) * pgwd
    
            Pdvdt(j,k)   = DTAUY  +Pdvdt(j,k)
            Pdudt(j,k)   = DTAUX  +Pdudt(j,k)
    
            unew   = U1(J,K) + DTAUX*dtp       !   Pdudt(J,K)*DTP
            vnew   = V1(J,K) + DTAUY*dtp       !   Pdvdt(J,K)*DTP
            ENG1   = 0.5*(unew*unew + vnew*vnew)
!
            DUSFC(J) = DUSFC(J)  + DTAUX * DEL(J,K)
            DVSFC(J) = DVSFC(J)  + DTAUY * DEL(J,K)
!2018-diag
            dudt_ogw(j,k) = DTAUX
            tau_ogw(j)    = tau_ogw(j) +DTAUX*DEL(j,k)
          endif
!	  
! local energy deposition SSO-heat
!	
          Pdtdt(j,k) = max(ENG0-ENG1,0.)*rcpdt  
        ENDDO
      ENDDO
! dusfc w/o tofd  sign as in the ERA-I, MERRA  and CFSR
      DO I = 1,npt
         J           = ipt(i)
         DUSFC(J)    = -rgrav * DUSFC(J)
         DVSFC(J)    = -rgrav * DVSFC(J)
         tau_mtb(j)  = -rgrav * tau_mtb(j)
         tau_ogw(j)  = -rgrav * tau_ogw(j)
         tau_tofd(J) = -rgrav * tau_tofd(j)
       ENDDO

       RETURN


!============ debug ------------------------------------------------
       if (kdt <= 2 .and. me == 0) then
        print *, 'vgw-oro done gwdps_v0 in ugwp-v0 step-proc ', kdt, me
!
        print *, maxval(pdudt)*86400.,  minval(pdudt)*86400, 'vgw_axoro'
        print *, maxval(pdvdt)*86400.,  minval(pdvdt)*86400, 'vgw_ayoro'
!       print *, maxval(kdis),  minval(kdis),  'vgw_kdispro m2/sec'
        print *, maxval(pdTdt)*86400.,  minval(pdTdt)*86400,'vgw_epsoro'
        print *, maxval(zmtb), ' z_mtb ',  maxval(tau_mtb), ' tau_mtb '
        print *, maxval(zogw), ' z_ogw ',  maxval(tau_ogw), ' tau_ogw '
!       print *, maxval(tau_tofd),  ' tau_tofd '
!       print *, maxval(axtms)*86400.,  minval(axtms)*86400, 'vgw_axtms'
!       print *,maxval(dudt_mtb)*86400.,minval(dudt_mtb)*86400,'vgw_axmtb'
        if (maxval(abs(pdudt))*86400. > 100.) then

          print *, maxval(u1),  minval(u1),  ' u1 gwdps-v0 '
          print *, maxval(v1),  minval(v1),  ' v1 gwdps-v0 '
          print *, maxval(t1),  minval(t1),  ' t1 gwdps-v0 '
          print *, maxval(q1),  minval(q1),  ' q1 gwdps-v0 '
          print *, maxval(del), minval(del), ' del gwdps-v0 '
          print *, maxval(phil)*rgrav,minval(phil)*rgrav, 'zmet'
          print *, maxval(phii)*rgrav,minval(phii)*rgrav, 'zmeti'
          print *, maxval(prsi), minval(prsi), ' prsi '
          print *, maxval(prsL), minval(prsL), ' prsL '
          print *, maxval(RO), minval(RO), ' RO-dens '
          print *, maxval(bnv2(1:npt,:)), minval(bnv2(1:npt,:)),' BNV2 '
          print *, maxval(kpbl), minval(kpbl), ' kpbl '
          print *, maxval(sgh30), maxval(hprime), maxval(elvmax),'oro-d'
          print *
          do i =1, npt
            j= ipt(i)
            print *,zogw(J)/hprime(j), zmtb(j)/hprime(j),
     &              phil(j,1)/9.81, nint(hprime(j)/sigma(j))
!
!....................................................................
!
!   zogw/hp=5.9   zblk/hp=10.7    zm=11.1m   ridge/2=2,489m/9,000m
!   from 5 to 20 km , we need to count for "ridges" > dx/4 ~ 15 km
!   we must exclude blocking by small ridges
!   VAY-kref < iblk           zogw-lev 15        block-level:  39
!
! velmin => 1.0, 0.01, 0.1 etc.....unification of wind limiters
! MAX(SQRT(U1(J,K)*U1(J,K) + V1(J,K)*V1(J,K)), minwnd)
! MAX(DW2,DW2MIN) * RDZ * RDZ
! ULOW(I) = MAX(SQRT(UBAR(I)*UBAR(I) + VBAR(I)*VBAR(I)), 1.0)
! TEM      = MAX(VELCO(I,K)*VELCO(I,K), 0.1)
!              TEMV = 1.0 / max(VELCO(I,K), 0.01)
!     &                   * max(VELCO(I,K),0.01) 
!....................................................................     
          enddo
          print *
          stop
        endif
       endif
      
!
      RETURN
!---------------------------------------------------------------
! review of OLD-GFS code 2017/18 most substantial changes
!  a) kref > idxzb if idxzb > KPBL "OK" clipped-hill for OGW
!  b) tofd -sgh30                   "OK"
!
!  c) FR < Frc linear theory for taub-specification
!
!  d) solver of Palmer et al. (1987) => Linsat of McFarlane
!
!---------------------------------------------------------------   
      end subroutine gwdps_v0  
      
      
      
!===============================================================================
!    use fv3gfs-v0
!    first beta version of ugwp for fv3gfs-128
!    cires/swpc - jan 2018
!    non-tested wam ugwp-solvers in fv3gfs: "lsatdis", "dspdis", "ado99dis"
!    they reqiure extra-work to put them in with intializtion and namelists
!          next will be lsatdis for both fv3wam & fv3gfs-128l implementations
!    with (a) stochastic-deterministic propagation solvers for wave packets/spectra
!         (b) gw-sources: oro/convection/dyn-instability (fronts/jets/pv-anomalies)
!         (c) guidance from high-res runs for GW sources and res-aware tune-ups  
!23456 
!
!      call gwdrag_wam(1,  im,   ix,  km,   ksrc, dtp,
!     & xlat, gw_dudt, gw_dvdt,  taux,  tauy)
!        call fv3_ugwp_wms17(kid1, im,  ix,  km,  ksrc_ifs, dtp,
!     &  adt,adu,adv,prsl,prsi,phil,xlat, gw_dudt, gw_dvdt, gw_dtdt, gw_ked,
!     &  taux,tauy,grav, amol_i, me, lstep_first )
!
!
!23456==============================================================================

!>\ingroup cires_ugwp_run
!> @{
!! 
      subroutine fv3_ugwp_solv2_v0(klon, klev, dtime,
     &           tm1 , um1, vm1, qm1,
     &           prsl, prsi,   philg, xlatd, sinlat, coslat,
     &           pdudt, pdvdt, pdtdt, dked, tau_ngw,
     &           mpi_id, master, kdt)
!


!=======================================================
!
!      nov 2015 alternative gw-solver for nggps-wam
!      nov 2017 nh/rotational gw-modes for nh-fv3gfs
! ---------------------------------------------------------------------------------
!  
  
      use ugwp_common ,     only : rgrav,      grav,  cpd,    rd,  rv
     &,                            omega2,     rcpd2,  pi,    pi2, fv
     &,                            rad_to_deg, deg_to_rad
     &,                            rdi,        gor,    grcp,   gocp
     &,                            bnv2min,    dw2min, velmin, gr2
!
      use ugwp_wmsdis_init, only : hpscale, rhp2,    bv2min,   gssec
     &,                            v_kxw,   v_kxw2,  tamp_mpa, zfluxglob
     &,                            maxdudt, gw_eff,  dked_min
     &,                            nslope,  ilaunch, zmsi
     &,                            zci,     zdci,    zci4, zci3, zci2
     &,                            zaz_fct, zcosang, zsinang
     &,                            nwav,    nazd,    zcimin, zcimax
! 
      implicit none
!23456 
       
      integer, intent(in) :: klev             ! vertical level
      integer, intent(in) :: klon             ! horiz tiles

      real,    intent(in) :: dtime            ! model time step       
      real,    intent(in) :: vm1(klon,klev)   ! meridional wind 
      real,    intent(in) :: um1(klon,klev)   ! zonal wind
      real,    intent(in) :: qm1(klon,klev)   ! spec. humidity
      real,    intent(in) :: tm1(klon,klev)   ! kin temperature 

      real,    intent(in) :: prsl(klon,klev)  ! mid-layer pressure
      real,    intent(in) :: philg(klon,klev) ! m2/s2-phil => meters !!!!!       phil =philg/grav
      real,    intent(in) :: prsi(klon,klev+1)! prsi interface pressure
      real,    intent(in) :: xlatd(klon)      ! lat was in radians, now with xlat_d in degrees
      real,    intent(in) :: sinlat(klon)
      real,    intent(in) :: coslat(klon)
      real,    intent(in) :: tau_ngw(klon)

      integer, intent(in) :: mpi_id, master, kdt
!  
!
! out-gw effects
!
      real, intent(out) :: pdudt(klon,klev)     ! zonal momentum tendency
      real, intent(out) :: pdvdt(klon,klev)     ! meridional momentum tendency
      real, intent(out) :: pdtdt(klon,klev)     ! gw-heating (u*ax+v*ay)/cp
      real, intent(out) :: dked(klon,klev)      ! gw-eddy diffusion	
      real, parameter   :: minvel = 0.5         !      
      real, parameter   :: epsln  = 1.0d-12     !      
                 
!vay-2018
   
      real              :: taux(klon,klev+1)    ! EW component of vertical momentum flux (pa)
      real              :: tauy(klon,klev+1)    ! NS component of vertical momentum flux (pa)
      real              :: phil(klon,klev)      ! gphil/grav	
!
! local ===============================================================================================
!
 
!      real  :: zthm1(klon,klev)                       ! temperature interface levels   	 
       real  :: zthm1                                  ! 1.0 / temperature interface levels   	 
       real  :: zbvfhm1(klon,ilaunch:klev)             ! interface BV-frequency 
       real  :: zbn2(klon,ilaunch:klev)                ! interface BV-frequency 
       real  :: zrhohm1(klon,ilaunch:klev)             ! interface density 
       real  :: zuhm1(klon,ilaunch:klev)               ! interface zonal wind
       real  :: zvhm1(klon,ilaunch:klev)               ! meridional wind
       real  :: v_zmet(klon,ilaunch:klev)
       real  :: vueff(klon,ilaunch:klev)
       real  :: zbvfl(klon)                            ! BV at launch level	 
       real  :: c2f2(klon)

!23456
       real  :: zul(klon,nazd)                ! velocity in azimuthal direction at launch level
       real  :: zci_min(klon,nazd)
!      real  :: zcrt(klon,klev,nazd)          ! not used - do we need it?  Moorthi
       real  :: zact(klon, nwav, nazd)        ! if =1 then critical level encountered => c-u
!      real  :: zacc(klon, nwav, nazd)        ! not used!
!
       real  :: zpu(klon,klev,  nazd)         ! momentum flux
!      real  :: zdfl(klon,klev, nazd)
       real  :: zfct(klon,klev)
       real  :: zfnorm(klon)                  ! normalisation factor

       real  ::  zfluxlaun(klon)
       real  ::  zui(klon, klev,nazd)
!
       real  :: zdfdz_v(klon,klev, nazd)   ! axj = -df*rho/dz       directional momentum depositiom
       real  :: zflux(klon, nwav, nazd)   ! momentum flux at each level   stored as ( ix, mode, iazdim)

       real  :: zflux_z (klon, nwav,klev)    !momentum flux at each azimuth stored as ( ix, mode, klev)
!
       real  :: vm_zflx_mode, vc_zflx_mode
       real  :: kzw2, kzw3, kdsat, cdf2, cdf1, wdop2

!      real  :: zang, znorm, zang1, ztx
       real  :: zu, zcin, zcpeak, zcin4, zbvfl4
       real  :: zcin2, zbvfl2, zcin3, zbvfl3, zcinc
       real  :: zatmp, zfluxs, zdep, zfluxsq, zulm, zdft, ze1, ze2

!  
       real  :: zdelp,zrgpts
       real  :: zthstd,zrhostd,zbvfstd
       real  :: tvc1,  tvm1, tem1, tem2, tem3
       real  :: zhook_handle
       real  :: delpi(klon,ilaunch:klev)


!      real  :: rcpd, grav2cpd
       real, parameter ::  rcpdl    = cpd/grav     ! 1/[g/cp]  == cp/g
     &,                    grav2cpd = grav/rcpdl   ! g*(g/cp)= g^2/cp
     &,                    cpdi     = 1.0d0/cpd
       
       real :: expdis, fdis
!      real :: fmode, expdis, fdis
       real :: v_kzi, v_kzw, v_cdp, v_wdp, sc, tx1

       integer :: j, k, inc, jk, jl, iazi
!       
!--------------------------------------------------------------------------
!
        do k=1,klev
          do j=1,klon
            pdvdt(j,k) = 0.0
            pdudt(j,k) = 0.0
            pdtdt(j,k) = 0.0
            dked(j,k)  = 0.0
            phil(j,k)  = philg(j,k) * rgrav
          enddo
        enddo
!-----------------------------------------------------------	
! also other options to alter tropical values
! tamp = 100.e-3*1.e3 = 100 mpa
! vay-2017   zfluxglob=> lat-dep here from geos-5/merra-2 
!-----------------------------------------------------------
!        call slat_geos5_tamp(klon, tamp_mpa, xlatd, tau_ngw)	

     
!        phil = philg*rgrav
 
!        rcpd     = 1.0/(grav/cpd)     ! 1/[g/cp]
!        grav2cpd = grav*grav/cpd      ! g*(g/cp)= g^2/cp

        if (kdt ==1 .and. mpi_id == master) then
          print *,  maxval(tm1), minval(tm1), 'vgw: temp-res '
          print *,  'ugwp-v0: zcimin=' , zcimin
          print *,  'ugwp-v0: zcimax=' , zcimax 
          print *
        endif
!
!=================================================
       do iazi=1, nazd
         do jk=1,klev
           do jl=1,klon
             zpu(jl,jk,iazi)  = 0.0
!            zcrt(jl,jk,iazi) = 0.0
!            zdfl(jl,jk,iazi) = 0.0
           enddo
         enddo
       enddo

!   
! set initial min Cxi for critical level absorption
       do iazi=1,nazd
         do jl=1,klon
           zci_min(jl,iazi) = zcimin
         enddo
       enddo
!       define half model level winds and temperature
!       ---------------------------------------------
       do jk=max(ilaunch,2),klev
         do jl=1,klon
           tvc1 = tm1(jl,jk)   * (1. +fv*qm1(jl,jk))
           tvm1 = tm1(jl,jk-1) * (1. +fv*qm1(jl,jk-1))
!          zthm1(jl,jk)   = 0.5 *(tvc1+tvm1)
           zthm1          = 2.0 / (tvc1+tvm1)
           zuhm1(jl,jk)   = 0.5 *(um1(jl,jk-1)+um1(jl,jk))
           zvhm1(jl,jk)   = 0.5 *(vm1(jl,jk-1)+vm1(jl,jk))
!          zrhohm1(jl,jk) = prsi(jl,jk)*rdi/zthm1(jl,jk)   !  rho = p/(RTv) 
           zrhohm1(jl,jk) = prsi(jl,jk)*rdi*zthm1          !  rho = p/(RTv) 
           zdelp          = phil(jl,jk)-phil(jl,jk-1)      !>0 ...... dz-meters
           v_zmet(jl,jk)  = zdelp + zdelp
           delpi(jl,jk)  = grav / (prsi(jl,jk-1) - prsi(jl,jk))
           vueff(jl,jk)   =
     &     2.e-5*exp( (phil(jl,jk)+phil(jl,jk-1))*rhp2)+dked_min
!
!          zbn2(jl,jk)    =  grav2cpd/zthm1(jl,jk)
           zbn2(jl,jk)    =  grav2cpd*zthm1
     &                    * (1.0+rcpdl*(tm1(jl,jk)-tm1(jl,jk-1))/zdelp)
           zbn2(jl,jk)    = max(min(zbn2(jl,jk), gssec), bv2min)
           zbvfhm1(jl,jk) = sqrt(zbn2(jl,jk))       ! bn = sqrt(bn2)          			
         enddo
       enddo

       if (ilaunch == 1) then
         jk = 1
         do jl=1,klon
!          zthm1(jl,jk)  = tm1(jl,jk) * (1. +fv*qm1(jl,jk))   ! not used
           zuhm1(jl,jk)  = um1(jl,jk)
           zvhm1(jl,jk)  = vm1(jl,jk)
           ZBVFHM1(JL,1) = ZBVFHM1(JL,2)
           V_ZMET(JL,1)  = V_ZMET(JL,2)
           VUEFF(JL,1)   = DKED_MIN
           ZBN2(JL,1)    = ZBN2(JL,2)
         enddo
       endif
       do jl=1,klon
         tx1       = OMEGA2 * SINLAT(JL) / V_KXW
         C2F2(JL)  = tx1 * tx1
         zbvfl(jl) = zbvfhm1(jl,ilaunch)
       enddo
!	
!        define intrinsic velocity (relative to launch level velocity) u(z)-u(zo), and coefficinets
!       ------------------------------------------------------------------------------------------        
        do iazi=1, nazd
          do jl=1,klon
            zul(jl,iazi) = zcosang(iazi) * zuhm1(jl,ilaunch)
     &                   + zsinang(iazi) * zvhm1(jl,ilaunch)
          enddo
        enddo
!
         do jk=ilaunch, klev-1     ! from z-launch up   model level from which gw spectrum is launched
           do iazi=1, nazd
             do jl=1,klon
               zu = zcosang(iazi)*zuhm1(jl,jk)
     &            + zsinang(iazi)*zvhm1(jl,jk)
               zui(jl,jk,iazi) =  zu - zul(jl,iazi)
             enddo
          enddo

        enddo
!                                         define rho(zo)/n(zo)
!       ------------------- 
      do jk=ilaunch, klev-1               
        do jl=1,klon
           zfct(jl,jk) = zrhohm1(jl,jk) / zbvfhm1(jl,jk)
        enddo
      enddo

!      ----------------------------------------- 
!       set launch momentum flux spectral density
!       ----------------------------------------- 

      if(nslope == 1) then                   ! s=1 case
                                             ! --------
        do inc=1,nwav
          zcin  = zci(inc)
          zcin4 = zci4(inc)
          do jl=1,klon
!n4
            zbvfl4 = zbvfl(jl) * zbvfl(jl)
            zbvfl4 = zbvfl4 * zbvfl4
            zflux(jl,inc,1) = zfct(jl,ilaunch)*zbvfl4*zcin
     &                      / (zbvfl4+zcin4)
          enddo
         enddo
      elseif(nslope == 2) then               ! s=2 case
                                             ! --------
        do inc=1, nwav
          zcin  = zci(inc)
          zcin4 = zci4(inc)
          do jl=1,klon
            zbvfl4 = zbvfl(jl) * zbvfl(jl)
            zbvfl4 = zbvfl4    * zbvfl4
            zcpeak = zbvfl(jl) * zmsi
            zflux(jl,inc,1) = zfct(jl,ilaunch)*
     &                     zbvfl4*zcin*zcpeak/(zbvfl4*zcpeak+zcin4*zcin)
          enddo
        enddo
      elseif(nslope == -1) then              ! s=-1 case
                                             ! --------
        do inc=1,nwav
          zcin  = zci(inc)
          zcin2 = zci2(inc)
          do jl=1,klon
            zbvfl2 = zbvfl(jl)*zbvfl(jl)
            zflux(jl,inc,1) = zfct(jl,ilaunch)*zbvfl2*zcin
     &                      / (zbvfl2+zcin2)
          enddo
        enddo
      elseif(nslope == 0) then               ! s=0 case
                                             ! --------

        do inc=1, nwav
           zcin  = zci(inc)
           zcin3 = zci3(inc)
          do jl=1,klon
            zbvfl3 = zbvfl(jl)**3
            zflux(jl,inc,1) = zfct(jl,ilaunch)*zbvfl3*zcin
     &                      / (zbvfl3+zcin3)
          enddo
        enddo

      endif  ! for slopes
!
! normalize momentum flux at the src-level
!       ------------------------------
! integrate (zflux x dx)
      do inc=1, nwav
        zcinc = zdci(inc)
        do jl=1,klon
          zpu(jl,ilaunch,1) = zpu(jl,ilaunch,1) + zflux(jl,inc,1)*zcinc
        enddo
      enddo
!      
!       normalize and include lat-dep  (precip or merra-2)
!       -----------------------------------------------------------
! also other options to alter tropical values
!
      do jl=1,klon
        zfluxlaun(jl) = tau_ngw(jl)     !*(.5+.75*coslat(JL))      !zfluxglob/2  on poles
        zfnorm(jl)    = zfluxlaun(jl) / zpu(jl,ilaunch,1)
      enddo
!
      do iazi=1,nazd
        do jl=1,klon
          zpu(jl,ilaunch,iazi) = zfluxlaun(jl)
        enddo
      enddo

!       adjust constant zfct

      do jk=ilaunch, klev-1
        do jl=1,klon
          zfct(jl,jk) = zfnorm(jl)*zfct(jl,jk)
        enddo
      enddo
!       renormalize each spectral mode

      do inc=1, nwav
        do jl=1,klon
          zflux(jl,inc,1) = zfnorm(jl)*zflux(jl,inc,1)
        enddo
      enddo

!     copy zflux into all other azimuths
!     --------------------------------
!     zact(:,:,:) = 1.0 ; zacc(:,:,:) = 1.0
      zact(:,:,:) = 1.0
      do iazi=2, nazd
        do inc=1,nwav
          do jl=1,klon
            zflux(jl,inc,iazi) = zflux(jl,inc,1)
          enddo
        enddo
      enddo

! -------------------------------------------------------------      
!                                        azimuth do-loop
! --------------------
      do iazi=1, nazd

!     write(0,*)' iazi=',iazi,' ilaunch=',ilaunch
!                                       vertical do-loop
! ----------------
        do jk=ilaunch, klev-1
! first check for critical levels
! ------------------------
          do jl=1,klon
            zci_min(jl,iazi) = max(zci_min(jl,iazi),zui(jl,jk,iazi))
          enddo
! set zact to zero if critical level encountered
! ----------------------------------------------
          do inc=1, nwav
!           zcin = zci(inc)
            do jl=1,klon
!             zatmp = minvel + sign(minvel,zcin-zci_min(jl,iazi))
!             zacc(jl,inc,iazi) = zact(jl,inc,iazi)-zatmp
!             zact(jl,inc,iazi) = zatmp
              zact(jl,inc,iazi) =  minvel
     &                          + sign(minvel,zci(inc)-zci_min(jl,iazi))
            enddo
          enddo
!
!   zdfl not used! - do we need it? Moorthi
! integrate to get critical-level contribution to mom deposition
! ---------------------------------------------------------------
!         do inc=1, nwav
!           zcinc = zdci(inc)
!           do jl=1,klon
!             zdfl(jl,jk,iazi) = zdfl(jl,jk,iazi) +
!    &                zacc(jl,inc,iazi)*zflux(jl,inc,iazi)*zcinc
!           enddo
!         enddo
! --------------------------------------------
! get weighted average of phase speed in layer  zcrt is not used - do we need it? Moorthi
! --------------------------------------------
!         do jl=1,klon
!     write(0,*)' jk=',jk,' jl=',jl,' iazi=',iazi, zdfl(jl,jk,iazi)
!           if(zdfl(jl,jk,iazi) > epsln ) then
!             zatmp = zcrt(jl,jk,iazi)
!             do inc=1, nwav
!                 zatmp = zatmp + zci(inc) *
!    &                    zacc(jl,inc,iazi)*zflux(jl,inc,iazi)*zdci(inc)
!             enddo
!
!             zcrt(jl,jk,iazi)  = zatmp / zdfl(jl,jk,iazi)
!           else
!              zcrt(jl,jk,iazi) = zcrt(jl,jk-1,iazi)
!           endif
!         enddo

!
          do inc=1, nwav
            zcin  = zci(inc)
            if (abs(zcin) > epsln) then
              zcinc = 1.0 / zcin
            else
              zcinc = 1.0
            endif
            do jl=1,klon
!=======================================================================
! saturated limit    wfit = kzw*kzw*kt; wfdt = wfit/(kxw*cx)*betat
! & dissipative      kzi = 2.*kzw*(wfdm+wfdt)*dzpi(k)
!           define   kxw = 
!=======================================================================  
              v_cdp =  abs(zcin-zui(jL,jk,iazi))
              v_wdp = v_kxw*v_cdp
              wdop2 = v_wdp* v_wdp
              cdf2  = v_cdp*v_cdp - c2f2(jL)
              if (cdf2 > 0) then
                kzw2 = (zBn2(jL,jk)-wdop2)/Cdf2  - v_kxw2
              else
                kzw2 = 0.0
              endif
              if ( kzw2 > 0 ) then
                v_kzw = sqrt(kzw2)
!	       
!linsatdis:  kzw2, kzw3, kdsat, c2f2,  cdf2, cdf1
!	             
!kzw2 = (zBn2(k)-wdop2)/Cdf2  - rhp4 - v_kx2w  ! full lin DS-NiGW (N2-wd2)*k2=(m2+k2+[1/2H]^2)*(wd2-f2) 
!              Kds = kxw*Cdf1*rhp2/kzw3
!
                v_cdp  = sqrt( cdf2 )
                v_wdp  = v_kxw *  v_cdp
                v_kzi  = abs(v_kzw*v_kzw*vueff(jl,jk)/v_wdp*v_kzw)
                expdis = exp(-v_kzi*v_zmet(jl,jk))
              else
                v_kzi  = 0.
                expdis = 1.0
                v_kzw  = 0.
                v_cdp  = 0.   ! no effects of reflected waves
              endif
       
!             fmode =  zflux(jl,inc,iazi)
!             fdis  =  fmode*expdis
              fdis  =  expdis * zflux(jl,inc,iazi)
!
! saturated flux + wave dissipation - Keddy_gwsat in UGWP-V1
!  linsatdis = 1.0 , here:   u'^2 ~ linsatdis* [v_cdp*v_cdp]
!
              zfluxs = zfct(jl,jk)*v_cdp*v_cdp*zcinc
!                                     
!             zfluxs= zfct(jl,jk)*(zcin-zui(jl,jk,iazi))**2/zcin
! flux_tot - sat.flux
! 
              zdep = zact(jl,inc,iazi)* (fdis-zfluxs)
              if(zdep > 0.0 ) then
! subs on sat-limit
                zflux(jl,inc,iazi) = zfluxs
                zflux_z(jl,inc,jk) = zfluxs
              else
! assign dis-ve flux
                zflux(jl,inc,iazi) = fdis
                zflux_z(jl,inc,jk) = fdis
              endif
            enddo
          enddo
!
! integrate over spectral modes  zpu(y, z, azimuth)    zact(jl,inc,iazi)*zflux(jl,inc,iazi)*[d("zcinc")]
!
          zdfdz_v(:,jk,iazi) = 0.0
 
          do inc=1, nwav
            zcinc = zdci(inc)                    ! dc-integration
            do jl=1,klon
              vc_zflx_mode    = zact(jl,inc,iazi)*zflux(jl,inc,iazi)
              zpu(jl,jk,iazi) = zpu(jl,jk,iazi) + vc_zflx_mode*zcinc
              
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! check monotonic decrease
!     (heat deposition integration over spectral mode for each azimuth
!      later sum over selected azimuths as "non-negative" scalars)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              if (jk > ilaunch)then
!               zdelp = grav/(prsi(jl,jk-1)-prsi(jl,jk))*
!    &                 abs(zcin-zui(jl,jk,iazi)) *zcinc
                zdelp = delpi(jl,jk) * abs(zcin-zui(jl,jk,iazi)) *zcinc
                vm_zflx_mode = zact(jl,inc,iazi)* zflux_z(jl,inc,jk-1)

                if (vc_zflx_mode > vm_zflx_mode)
     &              vc_zflx_mode = vm_zflx_mode ! no-flux increase
                zdfdz_v( jl,jk,iazi) = zdfdz_v( jl,jk,iazi) +
     &                             (vm_zflx_mode-vc_zflx_mode)*zdelp    ! heating >0
!
!
               endif
            enddo                          !jl=1,klon             
          enddo                            !waves inc=1,nwav

! --------------
        enddo                              ! end jk do-loop vertical loop
! ---------------
      enddo                               ! end nazd do-loop
! ----------------------------------------------------------------------------    
!       sum contribution for total zonal and meridional flux +
!           energy dissipation
!       ---------------------------------------------------
!      
       do jk=1,klev+1
         do jl=1,klon
           taux(jl,jk) = 0.0 
           tauy(jl,jk) = 0.0 
         enddo
       enddo     
    
       tem3 = zaz_fct*cpdi
       do iazi=1,nazd
         tem1 = zaz_fct*zcosang(iazi)
         tem2 = zaz_fct*zsinang(iazi)
         do jk=ilaunch,  klev-1
           do jl=1,klon
             taux(jl,jk)  = taux(jl,jk)  + tem1 * zpu(jl,jk,iazi)       ! zaz_fct - "azimuth"-norm-n
             tauy(jl,jk)  = tauy(jl,jk)  + tem2 * zpu(jl,jk,iazi)
             pdtdt(jl,jk) = pdtdt(jl,jk) + tem3 * zdfdz_v(jl,jk,iazi)   ! eps_dis =sum( +d(flux_e)/dz) > 0.
           enddo
         enddo

       enddo
!
!    update du/dt and dv/dt tendencies   ..... no contribution to heating => keddy/tracer-mom-heat
!    ----------------------------   
!

       do jk=ilaunch,klev
         do jl=1, klon
!          zdelp = grav / (prsi(jl,jk-1)-prsi(jl,jk))
           zdelp = delpi(jl,jk)
           ze1   = (taux(jl,jk)-taux(jl,jk-1))*zdelp
           ze2   = (tauy(jl,jk)-tauy(jl,jk-1))*zdelp
           if (abs(ze1) >= maxdudt ) then
             ze1 = sign(maxdudt, ze1)
           endif 
           if (abs(ze2) >= maxdudt ) then
             ze2 = sign(maxdudt, ze2)
           endif
           pdudt(jl,jk) = -ze1
           pdvdt(jl,jk) = -ze2
!
! Cx =0 based Cx=/= 0. above
!
           pdtdt(jl,jk) = (ze1*um1(jl,jk) + ze2*vm1(jl,jk)) * cpdi
!
           dked(jl,jk)  = max(dked_min, pdtdt(jl,jk)/zbn2(jl,jk))
!          if (dked(jl,jk) < 0)  dked(jl,jk) = dked_min
         enddo
       enddo
!	
! add limiters/efficiency for "unbalanced ics" if it is needed
!       
       do jk=ilaunch,klev
         do jl=1, klon
           pdudt(jl,jk) = gw_eff * pdudt(jl,jk)
           pdvdt(jl,jk) = gw_eff * pdvdt(jl,jk)
           pdtdt(jl,jk) = gw_eff * pdtdt(jl,jk)
           dked(jl,jk)  = gw_eff * dked(jl,jk)
         enddo
       enddo
!  
!--------------------------------------------------------------------------- 
!
       if (kdt == 1 .and. mpi_id == master) then
         print *, 'vgw done  '
!
         print *, maxval(pdudt)*86400.,  minval(pdudt)*86400, 'vgw ax'
         print *, maxval(pdvdt)*86400.,  minval(pdvdt)*86400, 'vgw ay'
         print *, maxval(dked)*1.,  minval(dked)*1,  'vgw keddy m2/sec'
         print *, maxval(pdtdt)*86400.,  minval(pdtdt)*86400,'vgw eps'
!
!        print *, ' ugwp -heating rates '
       endif

       return
       end subroutine fv3_ugwp_solv2_v0
!-------------------------------------------------------------------------------
!
! Part-3 of UGWP-V01 Dissipative (eddy) effects of UGWP it will be activated
!           after tests of OGW (new revision) and NGW with MERRA-2 forcing.
!
!-------------------------------------------------------------------------------
      subroutine edmix_ugwp_v0(im, levs,   dtp,
     &                         t1, u1, v1, q1, del,
     &                         prsl,  prsi, phil, prslk,
     &                         pdudt, pdvdt, pdTdt, pkdis,
     &                         ed_dudt, ed_dvdt, ed_dTdt,
     &                         me, master, kdt )
!
      use machine,           only : kind_phys
      use ugwp_common ,      only : rgrav, grav, cpd, rd, rdi, fv
!    &,                             pi, rad_to_deg, deg_to_rad, pi2
     &,                             bnv2min, velmin, arad

      implicit none

      integer, intent(in) :: me,  master, kdt
      integer, intent(in) :: im, levs
      real(kind=kind_phys), intent(in) :: dtp 
      real(kind=kind_phys), intent(in), dimension(im,levs) :: 
     &       u1,  v1, t1,   q1, del,  prsl, prslk, phil  
!       
      real(kind=kind_phys), intent(in),dimension(im,levs+1):: prsi
      real(kind=kind_phys),dimension(im,levs) :: pdudt, pdvdt, pdTdt
      real(kind=kind_phys),dimension(im,levs) :: pkdis
!
! out
!
      real(kind=kind_phys),dimension(im,levs) :: ed_dudt, ed_dvdt
      real(kind=kind_phys),dimension(im,levs) :: ed_dTdt 
!         
! locals
!
       integer :: i, j, k
!------------------------------------------------------------------------       
!      solving 1D-vertical eddy diffusion to "smooth" 
!      GW-related tendencies:   du/dt, dv/dt, d(PT)/dt
!      we need to use sum of molecular + eddy terms including turb-part
!       of PBL extended to the model top, because "phys-tend" dx/dt 
!       should be smoothed as "entire" fields therefore one should
!       first estimate and collect "effective" diffusion and applied
!       it to each part of tendency or "sum of tendencies + Xdyn"
!       this "diffusive-way" is tested with UGWP-tendencies
!       forced by various       wave sources. X' =dx/dt *dt
!       d(X + X')/dt = K*diff(X + X') =>
!   
!   wave1 dX'/dt = Kw * diff(X')... eddy part  "Kwave" on wave-part
!   turb2 dX/dt  = Kturb * diff(X) ... resolved scale mixing "Kturb"  like PBL
!       we may assume "zero-GW"-tendency at the top lid and "zero" flux
!       or "vertical gradient" near the surface
!
!  1-st trial w/o PBL interactions: add dU, dV dT tendencies
!  compute BV, SHR2, Ri => Kturb, Kturb + Kwave => Apply it to "X_Tend +X "
!  ed_X = X_ed - X => final eddy tendencies
!--------------------------------------------------------------------------- 
!    rzs=30m              dk       = rzs*rzs*sqrt(shr2(i,k))
!                  Ktemp     = dk/(1+5.*ri)**2  Kmom = Pr*Ktemp
!
      real(kind=kind_phys) :: Sw(levs), Sw1(levs), Fw(levs), Fw1(levs)
      real(kind=kind_phys) :: Km(levs), Kpt(levs), Pt(levs), Ptmap(levs)
      real(kind=kind_phys) :: rho(levs), rdp(levs), rdpm(levs-1)
      real(kind=kind_phys),dimension(levs) :: ktur, vumol, up, vp, tp
      real(kind=kind_phys),dimension(levs) :: bn2,  shr2, ksum
      real(kind=kind_phys) ::  eps_shr, eps_bn2, eps_dis
      real(kind=kind_phys) ::  rdz , uz, vz, ptz
! -------------------------------------------------------------------------     
! Prw*Lsat2 =1, for GW-eddy diffusion Pr_wave = Kv/Kt
!  Pr_wave ~1/Lsat2 = 1/Frcit2 = 2. => Lsat2 = 1./2 (Frc ~0.7)
!  m*u'/N = u'/{c-U) = h'N/(c-U) = Lsat = Fcrit 
! > PBL:  0.25 < prnum = 1.0 + 2.1*ri < 4
!     monin-edmf parameter(rlam=30.0,vk=0.4,vk2=vk*vk) rlamun=150.0
!
      real(kind=kind_phys), parameter :: iPr_pt = 0.5, dw2min = 1.e-4
      real(kind=kind_phys), parameter :: lturb = 30., sc2 = lturb*lturb 
      real(kind=kind_phys), parameter :: ulturb=150.,sc2u=ulturb* ulturb
      real(kind=kind_phys), parameter :: ric =0.25
      real(kind=kind_phys), parameter :: rimin = -10., prmin = 0.25
      real(kind=kind_phys), parameter :: prmax = 4.0
      real(kind=kind_phys), parameter :: hps = 7000., h4 = 0.25/hps
      real(kind=kind_phys), parameter :: kedmin  = 0.01, kedmax = 250.
      
            
      real(kind=kind_phys) :: rdtp, rineg, kamp, zmet, zgrow
      real(kind=kind_phys) ::  stab, stab_dt, dtstab, ritur
      integer              ::  nstab   
      real(kind=kind_phys) ::  w1, w2, w3
      rdtp  = 1./dtp
       nstab = 1
       stab_dt = 0.9999

      do i =1, im

        rdp(1:levs) = grav/del(i, 1:levs)
 
        up(1:levs) = u1(i,1:levs) +pdudt(i,1:levs)*dtp
        vp(1:levs) = v1(i,1:levs) +pdvdt(i,1:levs)*dtp
        tp(1:levs) = t1(i,1:levs) +pdTdt(i,1:levs)*dtp
        Ptmap(1:levs) = (1.+fv*q1(i,1:levs))/prslk(i,1:levs)
        rho(1:levs) = rdi*prsl(i, 1:levs)/tp(1:levs)
        Pt(1:levs) = tp(1:levs)*Ptmap(1:levs)

        do k=1, levs-1
          rdpm(k) = grav/(prsl(i,k)-prsl(i,k+1))
          rdz = .5*rdpm(k)*(rho(k)+rho(k+1))
          uz = up(k+1)-up(k)
          vz = vp(k+1)-vp(k)
          ptz =2.*(pt(k+1)-pt(k))/(pt(k+1)+pt(k))
          shr2(k) = rdz*rdz*(max(uz*uz+vz*vz, dw2min)) 
          bn2(k) = grav*rdz*ptz  
          zmet = phil(j,k)*rgrav
          zgrow = exp(zmet*h4)
          if ( bn2(k) < 0. ) then
!	    
! adjust PT-profile to bn2(k) = bnv2min	-- neutral atmosphere
!  adapt "pdtdt = (Ptadj-Ptdyn)/Ptmap" 
!
             print *,' UGWP-V0 unstab PT(z) via gwdTdt ', bn2(k), k
    
             rineg = bn2(k)/shr2(k)
             bn2(k) = max(bn2(k), bnv2min)
             kamp = sqrt(shr2(k))*sc2u *zgrow
             ktur(k) =kamp* (1+8.*(-rineg)/(1+1.746*sqrt(-rineg)))
          endif
          ritur = max(bn2(k)/shr2(k), rimin)
          if (ritur > 0. ) then
            kamp = sqrt(shr2(k))*sc2 *zgrow
            w1 = 1./(1. + 5*ritur)
            ktur(k)= kamp * w1 * w1
          endif
          vumol(k) = 2.e-5 *exp(zmet/hps)
          ksum(k) =ktur(k)+Pkdis(i,k)+vumol(k)
          ksum(k) = max(ksum(k), kedmin)
          ksum(k) = min(ksum(k), kedmax)
          stab = 2.*ksum(k)*rdz*rdz*dtp
          if ( stab >= 1.0 ) then
               stab_dt = max(stab_dt, stab)
          endif 
        enddo
        nstab = max(1, nint(stab_dt)+1)
        dtstab = dtp / float(nstab)
        ksum(levs) = ksum(levs-1)
        Fw(1:levs) =  pdudt(i, 1:levs)
        Fw1(1:levs) = pdvdt(i, 1:levs)
        Km(1:levs) = ksum(1:levs) * rho(1:levs)* rho(1:levs)
 
        do j=1, nstab
          call diff_1d_wtend(levs, dtstab, Fw, Fw1, Km,
     &                       rdp, rdpm, Sw, Sw1)
          Fw  = Sw
          Fw1 = Sw1
        enddo

        ed_dudt(i,:) = Sw
        ed_dvdt(i,:) = Sw1
 
        Pt(1:levs) = t1(i,1:levs)*Ptmap(1:levs)
        Kpt = Km*iPr_pt
        Fw(1:levs) =  pdTdt(i, 1:levs)*Ptmap(1:levs)
        do j=1, nstab
          call diff_1d_ptend(levs, dtstab, Fw, Kpt, rdp, rdpm, Sw)
          Fw  = Sw
        enddo
         ed_dtdt(i,1:levs) = Sw(1:levs)/Ptmap(1:levs)

      enddo 
  
      end subroutine edmix_ugwp_v0

      subroutine diff_1d_wtend(levs, dt, F, F1, Km, rdp, rdpm, S, S1)
      use machine,      only: kind_phys
      implicit none
      integer :: levs
      real(kind=kind_phys) :: dt
      real(kind=kind_phys) :: S(levs), S1(levs), F(levs), F1(levs)
      real(kind=kind_phys) :: Km(levs),  rdp(levs), rdpm(levs-1)
      integer              :: i, k
      real(kind=kind_phys) ::      Kp1, ad, cd, bd 
!     real(kind=kind_phys) :: km1, Kp1, ad, cd, bd 
!      S(:) = 0.0 ; S1(:) = 0.0
!         
! explicit diffusion solver
!
      k = 1
!     km1 = 0. ; ad =0.
                 ad =0.
      kp1 = .5*(Km(k)+Km(k+1))
      cd = rdp(1)*rdpm(1)*kp1*dt
      bd = 1. - cd - ad
!     S(k) = cd*F(k+1) + ad *F(k-1) + bd *F(k)
      S(K)  = F(k)
      S1(K) = F1(k)
      do k=2, levs-1
       ad = cd
       kp1 = .5*(Km(k)+Km(k+1))
       cd = rdp(k)*rdpm(k)*kp1*dt
       bd = 1.-(ad +cd)
       S(k)  = cd*F(k+1)  + ad *F(k-1)  + bd *F(k)
       S1(k) = cd*F1(k+1) + ad *F1(k-1) + bd *F1(k)
      enddo
       k     = levs
       S(k)  = F(k)
       S1(k) = F1(k)
      end subroutine diff_1d_wtend

      subroutine diff_1d_ptend(levs, dt, F, Km, rdp, rdpm, S)
       use machine,      only: kind_phys
       implicit none
       integer :: levs
      real(kind=kind_phys) :: dt
      real(kind=kind_phys) :: S(levs),   S1(levs),  F(levs), F1(levs)
      real(kind=kind_phys) :: Km(levs),  rdp(levs), rdpm(levs-1)
      integer              :: i, k
      real(kind=kind_phys) ::      Kp1, ad, cd, bd
!     real(kind=kind_phys) :: km1, Kp1, ad, cd, bd
!        
! explicit "eddy" smoother for tendencies
!           

      k = 1
!     km1 = 0. ; ad =0.
                 ad =0.
      kp1 = .5*(Km(k)+Km(k+1))
      cd = rdp(1)*rdpm(1)*kp1*dt
      bd = 1. -(cd +ad)
!      S(k) = cd*F(k+1) + ad *F(k-1) + bd *F(k)
      S(K) = F(k)
      do k=2, levs-1
       ad = cd
       kp1 = .5*(Km(k)+Km(k+1))
       cd = rdp(k)*rdpm(k)*kp1*dt
       bd = 1.-(ad +cd)
       S(k) = cd*F(k+1) + ad *F(k-1) + bd *F(k)
      enddo
       k = levs
       S(k) = F(k)
      end subroutine diff_1d_ptend

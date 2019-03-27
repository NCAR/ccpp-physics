!23456	 

       subroutine cires_ugwp_driver_v0(me,  master,
     &    im,  levs, ntrac, nvoro, dtp, kdt, imx,
     &    nmtvr,  do_tofd,  cdmbgwd,  xlat, xlatd, sinlat, coslat,  
     &    Statein,  del, oro_stat, sgh30, kpbl,
     &    dusfcg, dvsfcg, gw_dudt, gw_dvdt, gw_dtdt, gw_kdis,
     &    tau_tofd, tau_mtb, tau_ogw, tau_ngw,
     &    zmtb, zlwb, zogw, du3dt_mtb,du3dt_ogw, du3dt_tms,rdxzb )
!-----------------------------------------------------------
! Part 1 "old-revised" gfs-gwdps_v0
! Part 2  non-stationary multi-wave GWs FV3GFS-v0 
! Part 3  Dissipative version of UGWP-tendency application
!         (similar to WAM-2017)
!----------------------------------------------------------- 
        use machine,      only: kind_phys
!        use physcons,     only: con_cp, con_fvirt, con_g, con_rd,    
!     &                    con_rv, con_rerth, con_pi
     
        use GFS_typedefs,     only : GFS_statein_type
        use ugwp_wmsdis_init, only : tamp_mpa	 
		
	implicit none
!input
	type(GFS_statein_type),  intent(in) :: Statein 
	 
        integer, intent(in) :: me,  master
        integer, intent(in) :: im, levs, ntrac, nvoro
        integer, intent(in) :: kdt, imx,  nmtvr
	
	real(kind=kind_phys), intent(in) :: dtp, cdmbgwd(2)
	logical             :: do_tofd 
        integer, intent(in) :: kpbl(im)
	real(kind=kind_phys), intent(in), dimension(im) :: xlat, xlatd
	real(kind=kind_phys), intent(in), dimension(im) :: sinlat, coslat
	
	real(kind=kind_phys), intent(in) :: del(im, levs)
	
	real(kind=kind_phys), intent(in) :: oro_stat(im, nvoro), sgh30(im)
!out
        real(kind=kind_phys), dimension(im, levs) :: gw_dudt, gw_dvdt
        real(kind=kind_phys), dimension(im, levs) :: gw_dTdt, gw_kdis	
	
!-----locals	+ diagnostics output

        real(kind=kind_phys), dimension(im, levs) :: Pdvdt,Pdudt
        real(kind=kind_phys), dimension(im, levs) :: Pdtdt,Pkdis				
!	
        real(kind=kind_phys), dimension(im, levs) :: ed_dudt, ed_dvdt
        real(kind=kind_phys), dimension(im, levs) :: ed_dTdt
	
        real(kind=kind_phys), dimension(im)       :: dusfcg, dvsfcg	
	
        real(kind=kind_phys), dimension(im)       ::  rdxzb
        real(kind=kind_phys), dimension(im)       ::  zmtb,  
     &     zlwb, zogw, tau_mtb, tau_ogw, tau_tofd, tau_ngw
        real(kind=kind_phys), dimension(im, levs) :: du3dt_mtb,du3dt_ogw
        real(kind=kind_phys), dimension(im, levs) :: du3dt_tms     
! locals     
         integer :: i, j, k, ix
!
! define hprime, oc, oa4, clx, theta, sigma, gamm, elvmax
!
        real(kind=kind_phys), dimension(im)     ::  hprime, 
     &       oc, theta, sigma, gamm, elvmax  
        real(kind=kind_phys), dimension(im, 4)  :: clx, oa4 
!	
! switches that activate impact of OGWs and NGWs
!	
        real(kind=kind_phys),parameter :: pogw=1.,pngw=1.,pked=1. 
!
!	   
        if (me == master .and. kdt < 2) then
	 print *
	 write(6,*) 'FV3GFS execute ugwp_driver_v0 '
	 print *
	endif 
	
        ix = im
!         print *, ' NMTVR in driver ', nmtvr
       if(nmtvr == 14) then         ! current operational - as of 2014
         do i=1,im
	  hprime(i) = oro_stat(i,1)
          oc(i)     = oro_stat(i,2)
          oa4(i,1)  = oro_stat(i,3)
          oa4(i,2)  = oro_stat(i,4)
          oa4(i,3)  = oro_stat(i,5)
          oa4(i,4)  = oro_stat(i,6)
          clx(i,1)  = oro_stat(i,7)
          clx(i,2)  = oro_stat(i,8)
          clx(i,3)  = oro_stat(i,9)
          clx(i,4)  = oro_stat(i,10)
          theta(i)  = oro_stat(i,11)
          gamm(i)   = oro_stat(i,12)
          sigma(i)  = oro_stat(i,13)
          elvmax(i) = oro_stat(i,14) 
	 enddo
	endif   
!
! 1) ORO stationary GWs
!
!       pdvdt(:,:) = 0. ; pdudt(:,:) = 0. 
!       pkdis(:,:) = 0. ; pdtdt(:,:) = 0.
        zlwb(:)   = 0.
       
       CALL GWDPS_V0(IM, levs, imx, do_tofd,                    
     &               Pdvdt, Pdudt, Pdtdt, Pkdis,
     &  Statein%ugrs, Statein%vgrs, Statein%tgrs,                        
     &  Statein%qgrs(:,:,1),KPBL, Statein%prsi,del,Statein%prsl,                          
     &  Statein%prslk, Statein%phii, Statein%phil, DTP,KDT,         
     &  sgh30, HPRIME,OC,OA4, CLX, THETA,SIGMA,GAMM,ELVMAX,       
     &  DUSFCg, DVSFCg,                   
     &  nmtvr, cdmbgwd, me, master, rdxzb, 
     &  zmtb, zogw, tau_mtb, tau_ogw, tau_tofd,
     &  du3dt_mtb, du3dt_ogw, du3dt_tms) 
!       
!      
! non-stationary GW-scheme  with GMAO/MERRA GW-forcing
!     
        if (me == master .and. kdt < 2) then
	 print *
	 write(6,*) 'FV3GFS finished gwdps_v0 in ugwp_driver_v0 '
	 print *
	endif  
!--------	
! GMAO GEOS-5/MERRA GW-forcing	lat-dep
!--------
       call slat_geos5_tamp(im, tamp_mpa, xlatd, tau_ngw)
       
!       call slat_geos5(im, xlatd, tau_ngw)
!
! 2) non-stationary GWs with GEOS-5/MERRA GW-forcing
!       	               
       call fv3_ugwp_solv2_v0(im, levs, dtp,      
     &   Statein%tgrs, Statein%ugrs, Statein%vgrs,
     &   Statein%qgrs(:,:,1),      
     &   Statein%prsl, Statein%prsi, Statein%phil, xlatd,
     &   sinlat, coslat, gw_dudt, gw_dvdt, gw_dTdt, gw_kdis,
     &   tau_ngw, me, master, kdt )
	 
        if (me == master .and. kdt < 2) then
	 print *
	 write(6,*)'FV3GFS finished fv3_ugwp_v0 in ugwp_driver_v0 '
	 write(6,*) ' non-stationary GWs with GMAO/MERRA GW-forcing '
	 print *
	endif 
	 	 
        do k=1,levs
         do i=1,im      
           gw_dtdt(i,k) = pngw*gw_dtdt(i,k)+ pogw*Pdtdt(i,k)
           gw_dudt(i,k) = pngw*gw_dudt(i,k)+ pogw*Pdudt(i,k)
           gw_dvdt(i,k) = pngw*gw_dvdt(i,k)+ pogw*Pdvdt(i,k)
	   gw_kdis(i,k) = pngw*gw_kdis(i,k)+ pogw*Pkdis(i,k)
          enddo
         enddo	
	 
	 return
	 
!=================================================================================
! make "ugwp eddy-diffusion" update for gw_dtdt/gw_dudt/gw_dvdt by solving
! vertical diffusion equations & update "Statein%tgrs, Statein%ugrs, Statein%vgrs"
!=================================================================================
!
! 3) application of "eddy"-diffusion to "smooth" UGWP-related tendencies
!
        ed_dudt(:,:) =0.; ed_dvdt(:,:) = 0. ; ed_dtdt(:,:) = 0. 
	
!       call edmix_ugwp_betav0(im, levs, dtp,      
!     &   Statein%tgrs, Statein%ugrs, Statein%vgrs,       
!     &   Statein%prsl, Statein%prsi, Statein%phil,
!     &   gw_dudt, gw_dvdt, gw_dTdt,
!     &   ed_dudt, ed_dvdt, ed_dTdt, gw_kdis,
!     &   me, master, kdt )
     	
	gw_dtdt = gw_dtdt +  ed_dtdt*pked
	gw_dvdt = gw_dvdt +  ed_dvdt*pked
	gw_dudt = gw_dudt +  ed_dudt*pked
			   			
        end subroutine cires_ugwp_driver_v0
!=====================================================================	
!
!ugwp-v0 subroutines: GWDPS_V0 and fv3_ugwp_solv2_v0	
!  
!=====================================================================
      SUBROUTINE GWDPS_V0(IM, levs, imx, do_tofd,                   
     &    Pdvdt, Pdudt, Pdtdt, Pkdis, U1,V1,T1,Q1,KPBL,              
     &    PRSI,DEL,PRSL,PRSLK,PHII, PHIL,DTP,KDT,         
     &    sgh30, HPRIME,OC,OA4,CLX4,THETA,SIGMA,GAMMA,ELVMAXD,       
     &    DUSFC, DVSFC, nmtvr, cdmbgwd, me, master, rdxzb,                  
     &    zmtb, zogw, tau_mtb, tau_ogw, tau_tofd,
     &    dudt_mtb, dudt_ogw, dudt_tms)
!----------------------------------------
! ugwp_v0
!
! modified/revised version of gwdps.f (with bug fixes, tofd, appropriate
!   computation of kref for OGW + COORDE diagnostics
!   all constants/parameters inside cires_ugwp_initialize.F90 
!----------------------------------------    

      USE MACHINE , ONLY : kind_phys
      use ugwp_common ,      only : rgrav, grav, cpd, rd, rv, rcpd, rcpd2      
      use ugwp_common ,      only : pi, rad_to_deg, deg_to_rad, pi2       
      use ugwp_common ,      only : rdi, gor,  grcp,  gocp,  fv,    gr2
      use ugwp_common ,      only : bnv2min, dw2min, velmin, arad 
      
      use ugwp_oro_init, only: rimin, ric, efmin, efmax,hpmax,hpmin
      use ugwp_oro_init, only: dpmin, minwnd, hminmt,hncrit, sigfac    
      use ugwp_oro_init, only: RLOLEV,GMAX, VELEPS, FACTOP 
      use ugwp_oro_init, only: FRC, CE, CEOFRC, frmax, CG                  
      use ugwp_oro_init, only: FDIR, MDIR, NWDIR
      use ugwp_oro_init, only: cdmb, cleff, fcrit_gfs,fcrit_mtb 
      use ugwp_oro_init, only: n_tofd, ze_tofd, ztop_tofd    
                    
      use cires_ugwp_module, only  : kxw,  max_kdis, max_axyz 
      
!----------------------------------------                  
      implicit none
      
      integer, intent(in) :: im, levs, imx, kdt
      integer, intent(in) :: me, master, nmtvr
      integer             :: km
      logical, intent(in) :: do_tofd
           
      integer, intent(in)              :: KPBL(IM)    ! Index for the PBL top layer!
      real(kind=kind_phys), intent(in) :: dtp         !  time step
      real(kind=kind_phys), intent(in) :: cdmbgwd(2) 
          
      real(kind=kind_phys), intent(in), dimension(im,levs) :: 
     &                                   u1,  v1, t1,   q1,
     &                                   del,  prsl, prslk, phil   
      real(kind=kind_phys), intent(in),dimension(im,levs+1):: prsi, phii
      
      real(kind=kind_phys), intent(in) :: OC(IM), OA4(im,4), CLX4(im,4) 
      real(kind=kind_phys), intent(in) :: HPRIME(IM), sgh30(IM)                     
      real(kind=kind_phys), intent(in) :: ELVMAXD(IM), THETA(IM)
      real(kind=kind_phys), intent(in) :: SIGMA(IM),GAMMA(IM)  
      
!output -phys-tend
      real(kind=kind_phys),dimension(im,levs),intent(out):: Pdvdt,Pdudt
      real(kind=kind_phys),dimension(im,levs),intent(out):: Pkdis,Pdtdt
!     
!
! output          
! diag-coorde
      real(kind=kind_phys), dimension(im,levs), intent(out) ::
     &                      dudt_mtb, dudt_ogw, dudt_tms
!                     
      real(kind=kind_phys),dimension(im) :: RDXZB,  zmtb,   zogw
      real(kind=kind_phys),dimension(im) :: tau_ogw, tau_mtb, tau_tofd
      real(kind=kind_phys),dimension(im) :: dusfc, dvsfc                    
         
! 
! locals       
! mean flow
      real(kind=kind_phys) :: RI_N(IM,levs), BNV2(IM,levs), RO(IM,levs)
      real(kind=kind_phys) :: VTK(IM,levs),VTJ(IM,levs),VELCO(IM,levs)  
!mtb     
      real(kind=kind_phys) ::  OA(IM),  CLX(IM) , elvmax(im)       
      real(kind=kind_phys) ::  wk(IM)
      real(kind=kind_phys), dimension(im) :: PE, EK, UP      
      
      real(kind=kind_phys) :: DB(IM,levs),ANG(IM,levs),UDS(IM, levs)
      real(kind=kind_phys) :: ZLEN, DBTMP, R, PHIANG, DBIM, ZR
      real(kind=kind_phys) :: ENG0, ENG1, COSANG2, SINANG2
      real(kind=kind_phys) :: bgam, cgam, gam2       
!
! TOFD
!     Some constants now in "use ugwp_oro_init" +   "use ugwp_common"
!
!==================
      real(kind=kind_phys)   :: unew, vnew,  zpbl,  sigflt  
      real(kind=kind_phys), dimension(levs)    ::  utofd1, vtofd1
      real(kind=kind_phys), dimension(levs)    :: epstofd1, krf_tofd1 
      real(kind=kind_phys), dimension(levs)    ::  up1, vp1, zpm
      real(kind=kind_phys)                     ::  zsurf
      real(kind=kind_phys),dimension(im, levs) :: axtms, aytms 
! 
! OGW
!                     
      LOGICAL ICRILV(IM)
!
      real(kind=kind_phys) :: XN(IM), YN(IM), UBAR(IM),     
     &               VBAR(IM),  ULOW(IM),     
     &               ROLL(IM),  ULOI(IM),  bnv2bar(im), SCOR(IM), 
     &               DTFAC(IM), XLINV(IM), DELKS(IM), DELKS1(IM)
!
      real(kind=kind_phys) :: TAUP(IM,levs+1), TAUD(IM,levs)         
      real(kind=kind_phys) :: taub(im), taulin(im), heff, hsat, hdis       
         
      integer ::  kref(IM), idxzb(im), ipt(im), k_mtb 
      integer ::  kreflm(IM), iwklm(im), iwk(im)
      integer ::  ktrial, klevm1
!      
!check what we need
!
      real(kind=kind_phys) :: bnv,  fr, ri_gw ,          
     &                    brvf,   tem,   tem1,  tem2, temc, temv,  
     &                    ti,    rdz,   dw2,   shr2, bvf2,        
     &                    rdelks, efact, coefm, gfobnv,                   
     &                    scork,  rscor, hd,  fro,  sira,        
     &                    dtaux,  dtauy, pkp1log, pklog
     
      integer ::   kmm1, kmm2, lcap, lcapp1
      integer ::   npt       
      integer ::   kbps, kbpsp1,kbpsm1,               
     &       kmps, idir, nwd,  klcap, kp1, kmpbl, kmll                                
! 
       integer  ::  i, j, k
       real(kind=kind_phys)     :: grav2, rcpdt, windik, wdir
       real(kind=kind_phys)     :: sigmin, dxres, sigres
       real(kind=kind_phys)     :: cdmb4
       real(kind=kind_phys)     :: kxridge, inv_b2eff      
       rcpdt =1./cpd/dtp       
       grav2 = 2.*grav
!       
! mtb-blocking  sigma_min and dxres => cires_initialize
!       
       dxres = pi2*arad/float(IMX)
       sigmin = hpmin/dxres       
       kxridge = float(IMX)/arad * cdmbgwd(2)
       
       if (me == master .and. kdt==1) then
        print *, ' gwdps_v0 kxridge ', kxridge
        print *, ' gwdps_v0 scale2 ', cdmbgwd(2)
        print *, ' gwdps_v0 IMX ', imx	 	
       endif
              
       idxzb(:) = 0 	
       zmtb(:) = 0.    ; zogw(:) =0.
       rdxzb(:) =0.       
       tau_ogw(:) =0. ; tau_mtb(:) =0.; dusfc(:)=0. ; dvsfc(:)=0. 
       tau_tofd(:) =0.
       
       pdvdt  =0.   ; pdudt =0.    ; pdtdt =0. ; pkdis =0.
       dudt_mtb =0. ; dudt_ogw =0. ; dudt_tms =0. 
          
! ----  for lm and gwd calculation points
      
        ipt(:) = 0
        npt = 0
        do i = 1,im
          if ( (elvmaxd(i) .gt. hminmt) 
     &       .and. (hprime(i) .gt. hpmin) )  then
             npt      = npt + 1
             ipt(npt) = i
          endif
        enddo         

        IF (npt .eq. 0) then 	    
!        print *,  'oro-npt = 0 elvmax ', maxval(elvmaxd), hminmt
!        print *,  'oro-npt = 0 hprime ', maxval(hprime), hpmin	     	    
	RETURN     ! No gwd/mb calculation done
	endif
	
	
	  iwklm(1:npt)  = 2
          IDXZB(1:npt)  = 0 
          kreflm(1:npt) = 0
	  	
	db =0. ; ang = 0. ; uds =0. 
        km = levs
	KMM1   = levs- 1 ; KMM2   = levs - 2 ; KMLL   = kmm1
        LCAP   = levs   ;  LCAPP1 = LCAP + 1       
        
          DO I = 1, npt
            j = ipt(i)
            ELVMAX(J) = min (ELVMAXd(J)*0. + sigfac * hprime(j), hncrit)
          ENDDO
!
        DO K = 1, levs-1
          DO I = 1, npt
            j = ipt(i)
            pkp1log =  phil(j,k+1) *rgrav
            pklog =    phil(j,k)   *rgrav
        if (( ELVMAX(j).le.pkp1log) .and. (ELVMAX(j).ge.pklog) )
     &      iwklm(I)  =  MAX(iwklm(I), k+1 ) 

          ENDDO
        ENDDO	    
!
       RO = rdi*PRSL/T1

      DO K = 1,levs
        DO I =1,npt
          J         = ipt(i)
          VTJ(I,K)  = T1(J,K)  * (1.+FV*Q1(J,K))
          VTK(I,K)  = VTJ(I,K) / PRSLK(J,K)
          RO(I,K)   = RDI * PRSL(J,K) / VTJ(I,K)       ! DENSITY 
          TAUP(I,K) = 0.0
        ENDDO
      ENDDO
      bnv2(:,:) = 4.e-4
!
! check RI_N or RI_MF computation
!      
      DO K = 1,levs-1
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
          bnv2(i,k) = max( BVF2, bnv2min )
          RI_N(I,K) = Bnv2(i,k)/SHR2        ! Richardson number consistent with BNV2	  
        ENDDO
      ENDDO         
        K = levs 
        DO I = 1, npt
	   bnv2(i,k) =bnv2(i,k-1) 
	ENDDO	
!
        DO I = 1, npt
          J   = ipt(i)
          DELKS(I)  = 1.0 / (PRSI(J,1) - PRSI(J,iwklm(i)))
          DELKS1(I) = 1.0 / (PRSL(J,1) - PRSL(J,iwklm(i)))
          UBAR (I)  = 0.0
          VBAR (I)  = 0.0
          ROLL (I)  = 0.0
          PE   (I)  = 0.0
          EK   (I)  = 0.0
          BNV2bar(I) = (PRSL(J,1)-PRSL(J,2)) * DELKS1(I) * BNV2(I,1)
        ENDDO

! 
        DO Ktrial = KMLL, 1, -1
          DO I = 1, npt
             IF ( Ktrial .LT. iwklm(I) .and. kreflm(I) .eq. 0 ) then
                kreflm(I) = Ktrial
             ENDIF
          ENDDO
        ENDDO
! 
        DO I = 1, npt
          DO K = 1, Kreflm(I)
            J        = ipt(i)
            RDELKS     = DEL(J,K) * DELKS(I)
            UBAR(I)    = UBAR(I)  + RDELKS * U1(J,K) ! trial Mean U below 
            VBAR(I)    = VBAR(I)  + RDELKS * V1(J,K) ! trial Mean V below 
            ROLL(I)    = ROLL(I)  + RDELKS * RO(I,K) ! trial Mean RO below 
            RDELKS     = (PRSL(J,K)-PRSL(J,K+1)) * DELKS1(I)
            BNV2bar(I) = BNV2bar(I) + BNV2(I,K) * RDELKS
          ENDDO
        ENDDO
!
        DO I = 1, npt
          J = ipt(i)
          DO K = iwklm(I), 1, -1
            PHIANG   =  atan2(V1(J,K),U1(J,K))*RAD_TO_DEG
            ANG(I,K) = ( THETA(J) - PHIANG )
            if ( ANG(I,K) .gt.  90. ) ANG(I,K) = ANG(I,K) - 180.
            if ( ANG(I,K) .lt. -90. ) ANG(I,K) = ANG(I,K) + 180.
            ANG(I,K) = ANG(I,K) * DEG_TO_RAD
            UDS(I,K) = 
     &          MAX(SQRT(U1(J,K)*U1(J,K) + V1(J,K)*V1(J,K)), velmin)
!
            IF (IDXZB(I) .eq. 0 ) then
              PE(I) = PE(I) + BNV2(I,K) * 
     &           ( ELVMAX(J) - phil(J,K)*rgrav ) * 
     &           ( PHII(J,K+1) - PHII(J,K) ) *rgrav

              UP(I)  =  UDS(I,K) * cos(ANG(I,K))
              EK(I)  = 0.5 *  UP(I) * UP(I) 

! --- Dividing Stream lime  is found when PE =exceeds EK.
              IF ( PE(I) .ge.  EK(I) ) THEN
                 IDXZB(I) = K
		 zmtb (J) =PHII(J, K)*rgrav
!	  if (zmtb (J). gt. 3*hprime(j)) print *, 'ZBLK > 3*HP'
                 RDXZB(J) =real(k, kind=kind_phys)
		 
              ENDIF

            ENDIF
          ENDDO
        ENDDO

!
! --- The drag for mtn blocked flow
! 
        DO I = 1, npt
          J = ipt(i)
!
          IF ( IDXZB(I) .gt. 0 ) then
	  gam2 = gamma(j)*gamma(j)
	  BGAM = 1.0-0.18*gamma(j)-0.04*gam2
	  CGAM =    0.48*gamma(j) +0.30*gam2	  
            DO K = IDXZB(I)-1, 1, -1

                ZLEN = SQRT( ( PHIL(J,IDXZB(I)) - PHIL(J,K) ) / 
     &                       ( PHIL(J,K ) + Grav * hprime(J) ) )

               COSANG2 = cos(ANG(I,K))*cos(ANG(I,K))
               SINANG2 =1.0 -COSANG2  
               if ( abs(GAMMA(J) * COSANG2 + SINANG2) 
     &              .lt. 1.e-06 ) then
                 ZR = 2.0
               else
                 R = (COSANG2 + GAMMA(J) * SINANG2) /
     &              (GAMMA(J) * COSANG2 + SINANG2)
                 ZR =  MAX( 2. - 1. / R, 0. )
               endif
	   
           sigres = max(sigmin, sigma(J))
	   if (hprime(J)/sigres.gt.dxres) sigres=hprime(J)/dxres
	   
! (4.15)-IFS 	   
!           DBTMP = 0.25 *  CDmb * ZR * sigres*ZLEN /hprime(J) *
!     &             MAX(cos(ANG(I,K)), gamma(J)*sin(ANG(I,K))) 
! (4.16)-IFS     
           DBTMP = 0.25 *  CDmb * ZR * sigres*ZLEN / hprime(J) *
     &                  (bgam* COSANG2 +cgam* SINANG2)  
          
           DB(I,K) =  DBTMP * UDS(I,K)  
           ENDDO	        
!                  
         endif	       
        ENDDO	
! 
!.............................
!.............................
! end  mtn blocking section
!
!
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
     
      DO K=3,KMPBL
        DO I=1,npt
          j   = ipt(i)
          tem = (prsi(j,1) - prsi(j,k))
          if (tem .lt. dpmin) iwk(i) = k    ! dpmin=50 mb
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
      KBPS = 1
      KMPS = levs
      K_mtb =    1 
      DO I=1,npt
        J         = ipt(i)
	K_mtb     =    max(1, idxzb(i))
        kref(I)   = MAX(IWK(I), KPBL(J)+1 )             ! reference level 
      if (kref(i) .lt. idxzb(i)) kref(i) = idxzb(i) + 2 ! ~2-layers above zmtb
        KBPS      = MAX(KBPS,  kref(I))
        KMPS      = MIN(KMPS,  kref(I))      	
!        print *, 'VAY-kref < iblk ', kref(i), idxzb(i)
!	endif
        DELKS(I)  = 1.0 / (PRSI(J,k_mtb) - PRSI(J,kref(I)))
        DELKS1(I) = 1.0 / (PRSL(J,k_mtb) - PRSL(J,kref(I)))
        UBAR (I)  = 0.0
        VBAR (I)  = 0.0
        ROLL (I)  = 0.0
        BNV2bar(I)=(PRSL(J,k_mtb)-PRSL(J,k_mtb+1))* 
     &	           DELKS1(I)* BNV2(I,k_mtb)
      ENDDO
!  
      KBPSP1 = KBPS + 1
      KBPSM1 = KBPS - 1
        K_mtb  = 1    
      DO I = 1,npt
	K_mtb = max(1, idxzb(i))
       DO K = k_mtb,KBPS            !KBPS = MAX(kref) ;KMPS= MIN(kref)
          IF (K .LT. kref(I)) THEN
            J          = ipt(i)
            RDELKS     = DEL(J,K) * DELKS(I)
            UBAR(I)    = UBAR(I)  + RDELKS * U1(J,K)   ! Mean U below kref
            VBAR(I)    = VBAR(I)  + RDELKS * V1(J,K)   ! Mean V below kref
            ROLL(I)    = ROLL(I)  + RDELKS * RO(I,K)   ! Mean RO below kref
            RDELKS     = (PRSL(J,K)-PRSL(J,K+1)) * DELKS1(I)
            BNV2bar(I) = BNV2bar(I) + BNV2(I,K) * RDELKS
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
       ICRILV(I) = .FALSE. ! INITIALIZE CRITICAL LEVEL CONTROL VECTOR 
       ULOW(I) = MAX(SQRT(UBAR(I)*UBAR(I)+VBAR(I)*VBAR(I)),velmin)
       ULOI(I) = 1.0 / ULOW(I)
       XN(I)  = UBAR(I) * ULOI(I)
       YN(I)  = VBAR(I) * ULOI(I)       
      ENDDO
!
      DO  K = 1, levs-1
        DO  I = 1,npt
          J            = ipt(i)
          VELCO(I,K)   = 0.5 * ((U1(J,K)+U1(J,K+1))*XN(I)
     &                       + (V1(J,K)+V1(J,K+1))*YN(I))
        ENDDO
      ENDDO

!
!------------------
! v0: incorporates latest modifications for kxridge and heff/hsat
!             and taulin for Fr <=fcrit_gfs to make
!             and concept of "clipped" hill if zntb > 0.
! the integrated "tau_sso = tau_ogw +tau_mtb" close to reanalysis data
!------------------
      taub(:)  = 0. ; taulin(:)= 0.
      DO I = 1,npt
        J      = ipt(i)
        BNV    = SQRT( BNV2bar(I) )
        heff  = min(HPRIME(J),hpmax)
	
        if( zmtb(j) > 0.)heff=max(sigfac*heff-zmtb(j), 0.)/sigfac  
	if (heff .le. 0) cycle	 
	
	hsat     = fcrit_gfs*ULOI(I)/bnv
	heff     = min(heff, hsat)
	taulin(i)= CLEFF*BNV*0.5*ULOI(I)*heff*heff	
	
        FR     = BNV     * ULOI(I) * heff
        FR     = MIN(FR, FRMAX)
!
        EFACT    = (OA(I) + 2.) ** (CEOFRC*FR)
        EFACT    = MIN( MAX(EFACT,EFMIN), EFMAX )
!
        COEFM    = (1. + CLX(I)) ** (OA(I)+1.)
!
        XLINV(I) = COEFM * CLEFF           ! effective kxw for Lin-wave
!
        TEM      = FR    * FR * OC(J)
				
        GFOBNV   = GMAX  * TEM / ((TEM + CG)*BNV) 

!
!new specification of XLINV(I) & taulin(i)
!          sigres = max(sigmin, sigma(J))
!	   if (heff/sigres.gt.dxres) sigres=heff/dxres
!          inv_b2eff =  0.5*sigres/heff 	
!        XLINV(I) = max(kxridge, inv_b2eff)           ! 0.5*sigma(j)/heff = 1./Lridge  
!	taulin(i)= XLINV(I)*BNV*0.5*ULOI(I)*heff*heff

       if ( FR > fcrit_gfs ) then
        TAUB(I)  = XLINV(I) * ROLL(I) * ULOW(I) * ULOW(I)
     &           * ULOW(I)  * GFOBNV  * EFACT         ! nonlinear FLUX Tau0
!
         else
!
!        TAUB(I)  = XLINV(I) * ROLL(I) * ULOW(I) * BNV * heff*heff	 
!        TAUB(I)  = XLINV(I) * ROLL(I) * ULOW(I) * ULOW(I)
!     &           * ULOW(I)  * GFOBNV  * EFACT 
        TAUB(I)  = taulin(i)            
        endif	 
!
        K        = MAX(1, kref(I)-1)
        TEM      = MAX(VELCO(I,K)*VELCO(I,K), dw2min)
        SCOR(I)  = BNV2(I,K) / TEM  ! Scorer parameter below ref level
!
! diagnostics for zogw > zmtb
!
	zogw(J)    = PHII(j, kref(I)) *rgrav
      ENDDO
!    
!                                                                       
!----SET UP BOTTOM VALUES OF STRESS
!
      DO K = 1, KBPS
        DO I = 1,npt
          IF (K .LE. kref(I)) TAUP(I,K) = TAUB(I)
        ENDDO
      ENDDO
!================================================
!   V0-GFS OROGW-solver of Palmer et al 1986
!   V1-alternative LINSATDIS of WAM
!     with LLWB-mechanism for
!     rotational/non-hydrostat OGWs for 
!     HRES-FV3GFS with dx < 10 km
!=================================================

      DO K = KMPS, KMM1                   ! Vertical Level Loop
        KP1 = K + 1
        DO I = 1, npt
!
          IF (K .GE. kref(I)) THEN
            ICRILV(I) = ICRILV(I) .OR. ( RI_N(I,K) .LT. RIC)
     &                            .OR. (VELCO(I,K) .LE. 0.0)
          ENDIF
        ENDDO
!
        DO I = 1,npt
          IF (K .GE. kref(I))   THEN
            IF (.NOT.ICRILV(I) .AND. TAUP(I,K) .GT. 0.0 ) THEN
              TEMV = 1.0 / max(VELCO(I,K), velmin)
!
              IF (OA(I).GT.0. .AND. kp1 .lt. kref(i)) THEN
                SCORK   = BNV2(I,K) * TEMV * TEMV
                RSCOR   = MIN(1.0, SCORK / SCOR(I))
                SCOR(I) = SCORK
              ELSE 
                RSCOR   = 1.
              ENDIF
!
              BRVF = SQRT(BNV2(I,K))        ! Brunt-Vaisala Frequency
!             TEM1 = XLINV(I)*(RO(I,KP1)+RO(I,K))*BRVF*VELCO(I,K)*0.5

              TEM1 = XLINV(I)*(RO(I,KP1)+RO(I,K))*BRVF*0.5
     &                       * max(VELCO(I,K), velmin)
              HD   = SQRT(TAUP(I,K) / TEM1)
              FRO  = BRVF * HD * TEMV
!
!    RIM is the  "WAVE"-RICHARDSON NUMBER BY PALMER,Shutts, Swinbank 1986
!

              TEM2   = SQRT(ri_n(I,K))
              TEM    = 1. + TEM2 * FRO
              RI_GW    = ri_n(I,K) * (1.-FRO) / (TEM * TEM)
!
!    CHECK STABILITY TO EMPLOY THE 'dynamical SATURATION HYPOTHESIS'
!    OF PALMER,Shutts, Swinbank 1986
!                                       ----------------------
              IF (RI_GW .LE. RIC .AND.
     &           (OA(I) .LE. 0. .OR.  kp1 .ge. kref(i) )) THEN
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
         taup(1:npt,levs+1) = taup(1:npt,levs)      
!
!     Calculate - (grav)*d(tau)/d(p) and Decel terms DTAUX, DTAUY
!
      DO K = 1,KM
        DO I = 1,npt
          TAUD(I,K)=GRAV*(TAUP(I,K+1) - TAUP(I,K))/DEL(ipt(I),K)
        ENDDO
      ENDDO
!
!------scale MOMENTUM DEPOSITION  AT TOP TO 1/2 VALUE
!
      DO KLCAP = LCAP, KM
         DO I = 1,npt
            TAUD(I,KLCAP) = TAUD(I,KLCAP) * FACTOP
         ENDDO
      ENDDO
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!------IF THE GRAVITY WAVE DRAG WOULD FORCE A CRITICAL LINE IN THE
!------LAYERS BELOW SIGMA=RLOLEV DURING THE NEXT DELTIM TIMESTEP,
!------THEN ONLY APPLY DRAG UNTIL THAT CRITICAL LINE IS REACHED.
! Empirical implementation of the LLWB-mechanism: Lower Level Wave Breaking
! by limiting "Ax = Dtfac*Ax" due to possible LLWB around Kref and 500 mb
! critical line [V - Ax*dtp = 0.] is smt like "LLWB" for stationary OGWs
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      DO K = 1,KMM1
        DO I = 1,npt
         IF (K .GE. kref(I) .and. PRSI(ipt(i),K) .GE. RLOLEV) THEN
         
             IF(TAUD(I,K).NE.0.) THEN
               TEM = DTP * TAUD(I,K)
               DTFAC(I) = MIN(DTFAC(I),ABS(VELCO(I,K)/TEM))
!	       DTFAC(I) = 1.0
             ENDIF
          ENDIF
        ENDDO
      ENDDO
!
!---------------------------
!  
      IF( do_tofd ) then
      axtms(:,:) = 0.0 ; aytms(:,:) = 0.0 
       if ( kdt == 1) then
           print *, 'VAY do_tofd  from surface to ', ztop_tofd 
       endif
       DO I = 1,npt          
	 J = ipt(i)
	 zpbl  =rgrav*phil( j, kpbl(j) )
	 
         sigflt = min(sgh30(j), 0.2*hprime(j)) ! cannot exceed 20% of LS-SSO
	 
	 zsurf = phii(j,1)*rgrav
	 zpm(1:levs) =phiL(j,1:levs)*rgrav
	 up1(1:levs) = u1(j,1:levs)
	 vp1(1:levs) = v1(j,1:levs)
	 
         call ugwp_tofd1d(levs, sigflt, elvmaxd(j), zsurf, zpbl, 
     &	    up1, vp1, zpm,  utofd1, vtofd1, epstofd1, krf_tofd1)
     
         axtms(j,1:levs) = utofd1
         aytms(j,1:levs) = vtofd1
!	 
! add TOFD to GW-tendencies
!	 
	 pdvdt(J,:)  = pdvdt(J,:) + aytms(j,:)
	 pdudt(J,:)  = pdudt(J,:) + axtms(j,:)	     
!2018-diag
          tau_tofd(J) = sum( utofd1(1:levs)* del(j,1:levs))
        enddo
       ENDIF	! do_tofd 

!---------------------------
! combine oro-drag effects
!---------------------------  
! +  diag-3d

	 dudt_tms = axtms 
	 tau_ogw =0.
	 tau_mtb =0.
	 
      DO K = 1,KM
        DO I = 1,npt
          J          = ipt(i)
!

          ENG0       = 0.5*(U1(j,K)*U1(j,K)+V1(J,K)*V1(J,K))
!
          if ( K .lt. IDXZB(I) .AND. IDXZB(I) .ne. 0 ) then
!
! if blocking layers -- no OGWs
!	  
            DBIM = DB(I,K) / (1.+DB(I,K)*DTP)	    	    
            Pdvdt(j,k) = - DBIM * V1(J,K) +Pdvdt(j,k)
	    Pdudt(j,k) = - DBIM * U1(J,K) +Pdudt(j,k)	    
            ENG1    = ENG0*(1.0-DBIM*DTP)*(1.-DBIM*DTP)
            
            DUSFC(J)   = DUSFC(J) - DBIM * U1(J,K) * DEL(J,K)
            DVSFC(J)   = DVSFC(J) - DBIM * V1(J,K) * DEL(J,K)
!2018-diag 
            dudt_mtb(j,k) = -DBIM * U1(J,K)	    
	    tau_mtb(j) = tau_mtb(j) + dudt_mtb(j,k)* DEL(J,K)
	    
          else
!
! GW-s above blocking height
!	  
            TAUD(I,K)  = TAUD(I,K) * DTFAC(I)
            DTAUX      = TAUD(I,K) * XN(I)
            DTAUY      = TAUD(I,K) * YN(I)	  
	    
            Pdvdt(j,k)   = DTAUY  +Pdvdt(j,k)
	    Pdudt(j,k)   = DTAUX  +Pdudt(j,k)
	    
	    unew = U1(J,K) + DTAUX*dtp       !   Pdudt(J,K)*DTP
	    vnew = V1(J,K) + DTAUY*dtp       !   Pdvdt(J,K)*DTP
            ENG1   = 0.5*(unew*unew + vnew*vnew)	    	    	    
!
            DUSFC(J)   = DUSFC(J)  + DTAUX * DEL(J,K)
            DVSFC(J)   = DVSFC(J)  + DTAUY * DEL(J,K)
!2018-diag
	    dudt_ogw(j,k) =DTAUX 	    
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
         J          = ipt(i)
         DUSFC(J)   = -rgrav * DUSFC(J)
         DVSFC(J)   = -rgrav * DVSFC(J)
	 tau_mtb(j) = -rgrav * tau_mtb(j)
	 tau_ogw(j) = -rgrav * tau_ogw(j)	 
	 tau_tofd(J)= -rgrav * tau_tofd(j)
       ENDDO
  
       RETURN


!============ debug ------------------------------------------------      
       if (kdt .le. 2 ) then
        print *, 'vgw-oro done gwdps_v0 in ugwp-v0 step-proc ', kdt, me
!
        print *, maxval(pdudt)*86400.,  minval(pdudt)*86400, 'vgw_axoro'
        print *, maxval(pdvdt)*86400.,  minval(pdvdt)*86400, 'vgw_ayoro'
!        print *, maxval(kdis),  minval(kdis),  'vgw_kdispro m2/sec'
        print *, maxval(pdTdt)*86400.,  minval(pdTdt)*86400,'vgw_epsoro'
        print *, maxval(zmtb), ' z_mtb ',  maxval(tau_mtb), ' tau_mtb '
        print *, maxval(zogw), ' z_ogw ',  maxval(tau_ogw), ' tau_ogw '
!        print *, maxval(tau_tofd),  ' tau_tofd '
!        print *, maxval(axtms)*86400.,  minval(axtms)*86400, 'vgw_axtms'
!      print *,maxval(dudt_mtb)*86400.,minval(dudt_mtb)*86400,'vgw_axmtb'
        if (maxval(abs(pdudt))*86400. > 100.) then
	
	    print *,  maxval(u1), minval(u1), ' u1 gwdps-v0 '
	    print *,  maxval(v1), minval(v1), ' v1 gwdps-v0 '
	    print *,  maxval(t1), minval(t1), ' t1 gwdps-v0 '
	    print *,  maxval(q1), minval(q1), ' q1 gwdps-v0 '
	    print *,  maxval(del), minval(del), ' del gwdps-v0 '
	    print *,  maxval(phil)*rgrav,minval(phil)*rgrav, 'zmet'
	    print *,  maxval(phii)*rgrav,minval(phii)*rgrav, 'zmeti'
	    print *,  maxval(prsi), minval(prsi), ' prsi '
	    print *,  maxval(prsL), minval(prsL), ' prsL '
	    print *,  maxval(RO), minval(RO), ' RO-dens '
        print*,maxval(bnv2(1:npt,:)), minval(bnv2(1:npt,:)),' BNV2 '
	    print *,  maxval(kpbl), minval(kpbl), ' kpbl ' 
	print *, maxval(sgh30), maxval(hprime), maxval(elvmax),'oro-d'
	print *	 
	do i =1, npt
	   j= ipt(i)
      print *,zogw(J)/hprime(j), zmtb(j)/hprime(j),
     &	 phil(j,1)/9.81, nint(hprime(j)/sigma(j))
!
!....................................................................
!
!   zogw/hp=5.9   zblk/hp=10.7    zm=11.1m   ridge/2=2,489m/9,000m
!   from 5 to 20 km , we need to count for "ridges" > dx/4 ~ 15 km
!   we must exclude blocking by small ridges
!   VAY-kref < iblk           zogw-lev 15        block-level:  39
!
! velmin => 1.0, 0.01, 0.1 etc.....
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
!    with (a) stochastic-deterministic propagation solvers for wave packet/spectra
!         (b) gw-sources: oro/convection/dyn-instability (fronts/jets/pv-anomalies)
!         (c) guidance from high-res runs for gw sources and res-aware tune-ups  
!23456 
!
!      call gwdrag_wam(1,  im,   ix,  levs,   ksrc, dtp,
!     & xlat, gw_dudt, gw_dvdt,  taux,  tauy)
!        call fv3_ugwp_wms17(kid1, im,  ix,  levs,  ksrc_ifs, dtp,
!     &  adt,adu,adv,prsl,prsi,phil,xlat, gw_dudt, gw_dvdt, gw_dtdt, gw_ked,
!     &  taux,tauy,grav, amol_i, me, lstep_first )
!
!
!23456==============================================================================


      subroutine fv3_ugwp_solv2_v0(klon, klev, dtime,
     & tm1 , um1, vm1, qm1, 
     & prsl, prsi,   philg, xlatd, sinlat, coslat,
     & pdudt, pdvdt, pdtdt, dked, tau_ngw, mpi_id, master, kdt)
!


!=======================================================
!
!      nov 2015 alternative gw-solver for nggps-wam
!      nov 2017 nh/rotational gw-modes for nh-fv3gfs
! ---------------------------------------------------------------------------------
!  
  
      use ugwp_common ,      only : rgrav, grav, cpd, rd, rv
      use ugwp_common ,      only : omega2, rcpd2      
      use ugwp_common ,      only : pi, rad_to_deg, deg_to_rad, pi2       
      use ugwp_common ,      only : rdi, gor,  grcp,  gocp,  fv,  gr2
      use ugwp_common ,      only : bnv2min, dw2min, velmin                
!
      use ugwp_wmsdis_init, only : hpscale, rhp2, bv2min, gssec 
      use ugwp_wmsdis_init, only : v_kxw, v_kxw2, tamp_mpa, zfluxglob    
      use ugwp_wmsdis_init, only : maxdudt, gw_eff, dked_min
      use ugwp_wmsdis_init, only : nslope, ilaunch  
      use ugwp_wmsdis_init, only : zms, zci,  zdci, zci4, zci3, zci2
      use ugwp_wmsdis_init, only : zaz_fct, zcosang, zsinang 
      use ugwp_wmsdis_init, only : nwav, nazd  , zcimin, zcimax 
! 
       implicit none
!23456 
       
        integer  ,intent(in)  :: klev              ! vertical level
        integer  ,intent(in)  :: klon              ! horiz tiles
	
        real ,intent(in)   :: dtime               ! model time step       
        real ,intent(in)   :: vm1(klon,klev)      ! meridional wind 
        real ,intent(in)   :: um1(klon,klev)      ! zonal wind
        real ,intent(in)   :: qm1(klon,klev)      ! spec. humidity
        real ,intent(in)   :: tm1(klon,klev)      ! kin temperature 

        real ,intent(in) :: prsl(klon,klev)       ! mid-layer pressure
        real ,intent(in) :: philg(klon,klev)      ! m2/s2-phil => meters !!!!!       phil =philg/grav
        real ,intent(in) :: prsi(klon,klev+1)     !   prsi interface pressure
        real ,intent(in) :: xlatd(klon)           ! lat was in radians, now with xlat_d in degrees
        real ,intent(in) :: sinlat(klon)
        real ,intent(in) :: coslat(klon)          		
        real ,intent(in) :: tau_ngw(klon)
		 	
        integer, intent(in)   :: mpi_id, master, kdt
!  
!
! out-gw effects
!
        real ,intent(out):: pdudt(klon,klev)     ! zonal momentum tendency
        real ,intent(out):: pdvdt(klon,klev)     ! meridional momentum tendency
        real ,intent(out):: pdtdt(klon,klev)     ! gw-heating (u*ax+v*ay)/cp
        real ,intent(out):: dked(klon,klev)      ! gw-eddy diffusion	
        real, parameter  :: minvel = 0.5    !      
	                 
!vay-2018
   
        real             :: taux(klon,klev+1)  ! EW component of vertical momentum flux (pa)
        real             :: tauy(klon,klev+1)  ! NS component of vertical momentum flux (pa)
        real             :: phil(klon,klev)    !   gphil/grav	
!
! local ===============================================================================================
!
	 
         real  :: zbvfhm1(klon,klev), zbn2(klon,klev)    ! interface BV-frequency 
         real  :: zbvfl(klon)                            ! BV at launch level	 
         real  :: zrhohm1(klon,klev)                     ! interface density 
         real  :: zuhm1(klon,klev)                       ! interface zonal wind
         real  :: zvhm1(klon,klev)                       ! meridional wind
         real  :: zthm1(klon,klev)                       !temperature interface levels   	 
         real  :: v_zmet(klon, klev)	 	 
         real  :: vueff(klon, klev)	 
	 real  :: c2f2(klon)	 
	
!23456
         real  :: zul(klon,nazd)                !velocity in azimuthal direction at launch level
         real  :: zci_min(klon,nazd)
         real  :: zcrt(klon,klev,nazd)
         real  :: zact(klon, nwav, nazd)        !if =1 then critical level encountered => c-u
         real  :: zacc(klon, nwav, nazd)
!
         real  :: zpu(klon,klev,  nazd)         !momentum flux
         real  :: zdfl(klon,klev, nazd)	 
         real  :: zfct(klon,klev)
         real  :: zfnorm(klon)                 !normalisation factor
	 
       real  ::  zfluxlaun(klon)
       real  ::  zui(klon, klev,nazd)
!
       real  :: zdfdz_v(klon,klev, nazd)   ! axj = -df*rho/dz       directional momentum depositiom
       real  :: zflux(klon, nwav, nazd)   ! momentum flux at each level   stored as ( ix, mode, iazdim)

       real  :: zflux_z (klon, nwav,klev)    !momentum flux at each azimuth stored as ( ix, mode, klev)
!
       real  :: vm_zflx_mode, vc_zflx_mode
       real  :: kzw2, kzw3, kdsat, cdf2, cdf1, wdop2

       real  :: zang, znorm, zang1, ztx
       real  :: zu, zcin, zcpeak, zcin4, zbvfl4
       real  :: zcin2, zbvfl2, zcin3, zbvfl3, zcinc
       real  :: zatmp, zfluxs, zdep, zfluxsq, zulm, zdft, ze1, ze2

!  
       real  :: zdelp,zrgpts
       real  :: zthstd,zrhostd,zbvfstd
       real  :: tvc1,  tvm1       
       real  :: zhook_handle


       real  :: rcpd, grav2cpd
       
       real :: fmode, expdis, fdis
       real :: v_kzi, v_kzw, v_cdp, v_wdp, sc
                  
       integer :: j, k, inc, jk, jl, iazi             
!       
!--------------------------------------------------------------------------
!
        pdvdt  =0.   ; pdudt =0.    ; pdtdt =0. ; dked =0. 
!-----------------------------------------------------------	
! also other options to alter tropical values
! tamp = 100.e-3*1.e3 = 100 mpa
! vay-2017   zfluxglob=> lat-dep here from geos-5/merra-2 
!-----------------------------------------------------------
!        call slat_geos5_tamp(klon, tamp_mpa, xlatd, tau_ngw)	
	
     
         phil =philg*rgrav
                              	
         rcpd = 1./(grav/cpd)        ! 1/[g/cp]
         grav2cpd=grav*grav/cpd      ! g*(g/cp)= g^2/cp

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
         zpu(jl,jk,iazi) =0.0 
         zcrt(jl,jk,iazi)=0.0 
         zdfl(jl,jk,iazi)=0.0 
         enddo
         enddo      
       enddo  

!   
! set initial min Cxi for critical level absorption
       do iazi=1,nazd
        do jl=1,klon
            zci_min(jl,iazi)=zcimin
        enddo
       enddo
!       define half model level winds and temperature
!       -----------------------------------
       do jk=2,klev
        do jl=1,klon  
	tvc1 = tm1(jl,jk)*(1. +fv*qm1(jl,jk)) 
	tvm1 = tm1(jl,jk-1)*(1. +fv*qm1(jl,jk-1))                                     
        zthm1(jl,jk) =0.5 *(tvc1+tvm1) 
        zuhm1(jl,jk) =0.5 *(um1(jl,jk-1)+um1(jl,jk))
        zvhm1(jl,jk) =0.5 *(vm1(jl,jk-1)+vm1(jl,jk))
        zrhohm1(jl,jk)=prsi(jl,jk)*rdi/zthm1(jl,jk)                !  rho = p/(RTv) 
        zdelp=phil(jl,jk)-phil(jl,jk-1)      !>0 ...... dz-meters
        v_zmet(jl,jk) = 2.*zdelp
        vueff(jl,jk) = 
     &  2.e-5*exp( (phil(jl,jk)+phil(jl,jk-1))*rhp2)+dked_min
!     
        zbn2(jl,jk)= grav2cpd/zthm1(jl,jk)*
     &    (1.0 + rcpd*(tm1(jl,jk)-tm1(jl,jk-1))/zdelp)
         zbn2(jl,jk)=min(zbn2(jl,jk),gssec)
         zbn2(jl,jk)=max(zbn2(jl,jk),bv2min)
         zbvfhm1(jl,jk)=sqrt(zbn2(jl,jk))       ! bn = sqrt(bn2)          			
        enddo
       enddo

          jk=1
       do jl=1,klon                                        
         zthm1(jl,jk)=tm1(jl,jk)*(1. +fv*qm1(jl,jk))
         zuhm1(jl,jk)=um1(jl,jk)
         zvhm1(jl,jk)=vm1(jl,jk)
         ZBVFHM1(JL,1)  = ZBVFHM1(JL,2)
         V_ZMET(JL,1)  = V_ZMET(JL,2)
         VUEFF(JL,1)   = DKED_MIN
	 ZBN2(JL,1)   = ZBN2(JL,2)
	 C2F2(JL) = (OMEGA2*SINLAT(JL))**2/V_KXW2
         zbvfl(jl)=zbvfhm1(jl,ilaunch)	 	
        enddo
!	
!        define intrinsic velocity (relative to launch level velocity) u(z)-u(zo), and coefficinets
!       ------------------------------------------------------------------------------------------        
        do iazi=1, nazd
         do jl=1,klon
         zul(jl,iazi)=zcosang(iazi)*zuhm1(jl,ilaunch)+
     &                zsinang(iazi)*zvhm1(jl,ilaunch)
         enddo
        enddo
!
         do jk=ilaunch, klev-1     ! from z-launch up   model level from which gw spectrum is launched
          do iazi=1, nazd
           do jl=1,klon
          zu=zcosang(iazi)*zuhm1(jl,jk)+zsinang(iazi)*zvhm1(jl,jk)
          zui(jl,jk,iazi) =  zu - zul(jl,iazi)
          enddo
        enddo

        enddo
!                                         define rho(zo)/n(zo)
!       ------------------- 
      do jk=ilaunch, klev-1               
        do jl=1,klon
        zfct(jl,jk)=zrhohm1(jl,jk)/zbvfhm1(jl,jk)
        enddo
      enddo

!      ----------------------------------------- 
!       set launch momentum flux spectral density
!       ----------------------------------------- 

      if(nslope==1) then
! s=1 case
       do inc=1,nwav
             zcin=zci(inc)
             zcin4= zci4(inc)
       do jl=1,klon
!n4
         zbvfl4=zbvfl(jl)**4
         zflux(jl,inc,1)=zfct(jl,ilaunch)*zbvfl4*zcin/(zbvfl4+zcin4)
      enddo
      enddo
           elseif(nslope==2) then
! s=2 case
      do inc=1, nwav
      zcin=zci(inc)
      zcin4= zci4(inc)
        do jl=1,klon
         zbvfl4=zbvfl(jl)**4
         zcpeak=zbvfl(jl)/zms
        zflux(jl,inc,1)=zfct(jl,ilaunch)*
     &     zbvfl4*zcin*zcpeak/(zbvfl4*zcpeak+zcin4*zcin)
         enddo
       enddo
          elseif(nslope==-1) then
! s=-1 case
       do inc=1,nwav
         zcin=zci(inc)
         zcin2= zci2(inc)
       do jl=1,klon
        zbvfl2=zbvfl(jl)**2
        zflux(jl,inc,1)=zfct(jl,ilaunch)*zbvfl2*zcin/(zbvfl2+zcin2) 
       enddo
       enddo
!s=0 case
           elseif(nslope==0) then

       do inc=1, nwav
           zcin=zci(inc)
           zcin3= zci3(inc)
       do jl=1,klon
        zbvfl3=zbvfl(jl)**3
        zflux(jl,inc,1)=zfct(jl,ilaunch)*zbvfl3*zcin/(zbvfl3+zcin3)
       enddo
       enddo

       endif  ! for slopes
!
! normalize momentum flux at the src-level
!       ------------------------------
! integrate (zflux x dx)
      do inc=1, nwav
         zcinc=zdci(inc)
      do jl=1,klon
      zpu(jl,ilaunch,1)=zpu(jl,ilaunch,1)+zflux(jl,inc,1)*zcinc
      enddo    
      enddo
!      
!       normalize and include lat-dep  (precip or merra-2)
!       -----------------------------------------------------------
! also other options to alter tropical values
!
       do jl=1,klon       
        zfluxlaun(jl)=tau_ngw(jl)     !*(.5+.75*coslat(JL))      !zfluxglob/2  on poles
        zfnorm(jl)= zfluxlaun(jl)/zpu(jl,ilaunch,1)
       enddo
!
      do iazi=1,nazd
      do jl=1,klon
      zpu(jl,ilaunch,iazi)=zfluxlaun(jl)
      enddo
      enddo

!       adjust constant zfct

         do jk=ilaunch, klev-1
        do jl=1,klon
            zfct(jl,jk)=zfnorm(jl)*zfct(jl,jk)
        enddo
        enddo
!               renormalize each spectral mode

        do inc=1, nwav
         do jl=1,klon
          zflux(jl,inc,1)=zfnorm(jl)*zflux(jl,inc,1)
         enddo
         enddo

!       copy zflux into all other azimuths
!       --------------------------------
          zact(:,:,:) = 1.0 ; zacc(:,:,:) = 1.0
          do iazi=2, nazd
            do inc=1,nwav
            do jl=1,klon
            zflux(jl,inc,iazi)=zflux(jl,inc,1)
           enddo
          enddo
          enddo

! -------------------------------------------------------------      
!                                        azimuth do-loop
! --------------------
        do iazi=1, nazd
!                                       vertical do-loop
! ----------------
          do jk=ilaunch, klev-1                     
! first check for critical levels
! ------------------------
           do jl=1,klon
            zci_min(jl,iazi)=max(zci_min(jl,iazi),zui(jl,jk,iazi))               
           enddo
! set zact to zero if critical level encountered
! ----------------------------------------------
           do inc=1, nwav
            zcin=zci(inc)
            do jl=1,klon
               zatmp= minvel+sign(minvel,zcin-zci_min(jl,iazi))
               zacc(jl,inc,iazi)=zact(jl,inc,iazi)-zatmp
               zact(jl,inc,iazi)=zatmp
            enddo
           enddo
!
! integrate to get critical-level contribution to mom deposition
! ---------------------------------------------------------------
           do inc=1, nwav
            zcinc=zdci(inc)
            do jl=1,klon
               zdfl(jl,jk,iazi)=zdfl(jl,jk,iazi)+
     &                zacc(jl,inc,iazi)*zflux(jl,inc,iazi)*zcinc
            enddo
           enddo
! --------------------------------------------
! get weighted average of phase speed in layer
! --------------------------------------------
          do jl=1,klon
            if(zdfl(jl,jk,iazi)>0.0 ) then
               zatmp=zcrt(jl,jk,iazi)
             do inc=1, nwav
                  zatmp=zatmp+zci(inc)*
     &                   zacc(jl,inc,iazi)*zflux(jl,inc,iazi)*zdci(inc)
             enddo
!
               zcrt(jl,jk,iazi)=zatmp/zdfl(jl,jk,iazi)
            else
               zcrt(jl,jk,iazi)=zcrt(jl,jk-1,iazi)
            endif
         enddo

!
             do inc=1, nwav
                zcin=zci(inc)
                zcinc=1.0 /zcin
                do jl=1,klon
!=======================================================================
! saturated limit    wfit = kzw*kzw*kt; wfdt = wfit/(kxw*cx)*betat
! & dissipative      kzi = 2.*kzw*(wfdm+wfdt)*dzpi(k)
!           define   kxw = 
!=======================================================================  
               v_cdp =  abs(zcin-zui(jL,jk,iazi))
               v_wdp = v_kxw*v_cdp
	       wdop2 = v_wdp* v_wdp
	       cdf2 = v_cdp*v_cdp - c2f2(jL) 
	       if (cdf2 .gt. 0) then
	           kzw2 = (zBn2(jL,jk)-wdop2)/Cdf2  - v_kxw2  
		else
		 kzw2 = 0.0
	       endif	        	       
               if ( kzw2 .gt. 0 ) then 
	       v_kzw = sqrt(kzw2)
!	       
!linsatdis:  kzw2, kzw3, kdsat, c2f2,  cdf2, cdf1
!	             
!kzw2 = (zBn2(k)-wdop2)/Cdf2  - rhp4 - v_kx2w  ! full lin DS-NiGW (N2-wd2)*k2=(m2+k2+[1/2H]^2)*(wd2-f2) 
!              Kds = kxw*Cdf1*rhp2/kzw3
!
                v_cdp = sqrt( cdf2 ) 
	        v_wdp = v_kxw *  v_cdp    
                v_kzi  = abs(v_kzw*v_kzw*vueff(jl,jk)/v_wdp*v_kzw)
                  expdis = exp(-v_kzi*v_zmet(jl,jk))
               else
                  v_kzi = 0.
                  expdis = 1.0  
		  v_kzw = 0.
		  v_cdp = 0.   ! no effects of reflected waves
               endif
	       
               fmode =  zflux(jl,inc,iazi) 
               fdis  =  fmode*expdis
!
! saturated flux + wave dissipation - Keddy_gwsat in UGWP-V1
!  linsatdis = 1.0 , here:   u'^2 ~ linsatdis* [v_cdp*v_cdp]
!
               zfluxs= zfct(jl,jk)*v_cdp*v_cdp*zcinc
!                                     
!               zfluxs= zfct(jl,jk)*(zcin-zui(jl,jk,iazi))**2/zcin
! flux_tot - sat.flux
! 

               zdep=zact(jl,inc,iazi)* (fdis-zfluxs)		   
                  if(zdep>0.0 ) then
! subs on sat-limit
                      zflux(jl,inc,iazi)=zfluxs
                      zflux_z(jl,inc,jk)=zfluxs
                   else
! assign dis-ve flux
                      zflux(jl,inc,iazi)=fdis 
                      zflux_z(jl,inc,jk)=fdis 
                   endif  
               enddo                 
             enddo
!
! integrate over spectral modes  zpu(y, z, azimuth)    zact(jl,inc,iazi)*zflux(jl,inc,iazi)*[d("zcinc")]
!
           zdfdz_v(:,jk, iazi) =0.
        
          do inc=1, nwav
                 zcinc=zdci(inc)                    ! dc-integration
            do jl=1,klon

               vc_zflx_mode = zact(jl,inc,iazi)*zflux(jl,inc,iazi)
               zpu(jl,jk,iazi)=zpu(jl,jk,iazi) + vc_zflx_mode*zcinc
              
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
! check monotonic decrease
!     (heat deposition integration over spectral mode for each azimuth
!      later sum over selected azimuths as "non-negative" scalars)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if (jk.gt.ilaunch)then
       zdelp= grav/(prsi(jl,jk-1)-prsi(jl,jk))*
     &        abs(zcin-zui(jl,jk,iazi)) *zcinc
       vm_zflx_mode=zact(jl,inc,iazi)* zflux_z(jl,inc,jk-1)

      if (vc_zflx_mode.gt.vm_zflx_mode) vc_zflx_mode =vm_zflx_mode ! no-flux increase
      zdfdz_v( jl,jk,iazi) =zdfdz_v( jl,jk,iazi) +
     &           (vm_zflx_mode-vc_zflx_mode)*zdelp                 ! heating >0            
!
!              
       endif              
            enddo                          !jl=1,klon             
         enddo                             !waves inc=1,nwav

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
         taux(jl,jk)=0.0 
         tauy(jl,jk)=0.0 
       enddo
       enddo     
    
      do iazi=1,nazd
        do jk=ilaunch,  klev-1      
        do jl=1,klon
       taux(jl,jk)=taux(jl,jk)+zpu(jl,jk,iazi)*zaz_fct*zcosang(iazi)       ! zaz_fct - "azimuth"-norm-n
       tauy(jl,jk)=tauy(jl,jk)+zpu(jl,jk,iazi)*zaz_fct*zsinang(iazi)
      pdtdt(jl,jk) =pdtdt(jl,jk)+ zdfdz_v(jl,jk,iazi)*zaz_fct/cpd      ! eps_dis =sum( +d(flux_e)/dz) > 0.
        enddo
        enddo

      enddo
!
!    update du/dt and dv/dt tendencies   ..... no contribution to heating => keddy/tracer-mom-heat
!    ----------------------------   
! 
 

        do jk=ilaunch,klev
        do jl=1, klon
        zdelp= grav/(prsi(jl,jk-1)-prsi(jl,jk))            
        ze1=(taux(jl,jk)-taux(jl,jk-1))*zdelp
        ze2=(tauy(jl,jk)-tauy(jl,jk-1))*zdelp  
	if (abs(ze1) .ge. maxdudt ) then
	    ze1 = sign(maxdudt, ze1)
	endif 
	if (abs(ze2) .ge. maxdudt ) then
	    ze2 = sign(maxdudt, ze2)
	endif 	
        pdudt(jl,jk)=-ze1
        pdvdt(jl,jk)=-ze2
!	
! Cx =0 based Cx=/= 0. above
!
        pdtdt(jl,jk) = (ze1*um1(jl,jk) + ze2*vm1(jl,jk))/cpd
!
        dked(jl,jk)=pdtdt(jl,jk)/zbn2(jl,jk)
        if (dked(jl,jk).lt.0)  dked(jl,jk) = dked_min
        enddo
        enddo
!	
! add limiters/efficiency for "unbalanced ics"
!       
	pdudt = gw_eff*pdudt
	pdvdt = gw_eff*pdvdt
	pdtdt = gw_eff*pdtdt		
	dked =  gw_eff* dked	        	
!  
!--------------------------------------------------------------------------- 
!
       if (kdt == 1 .and. mpi_id == master) then
             print *, 'vgw done gwdrag_wms in ifs '
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

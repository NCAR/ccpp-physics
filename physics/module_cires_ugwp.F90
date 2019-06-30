!
module  cires_ugwp_module

!
!   driver is called after pbl & before chem-parameterizations
!
!....................................................................................
!  order = dry-adj=>conv=mp-aero=>radiation -sfc/land- chem -> vertdiff-> [rf-gws]=> ion-re
!...................................................................................
!
!
    implicit none 
    logical ::    module_is_initialized
!    logical            :: do_ugwp   = .false.               ! control => true - ugwp false old gws + rayeleigh friction

    logical            :: do_physb_gwsrcs   = .false.      ! control for physics-based GW-sources
    logical            :: do_rfdamp         = .false.      ! control for Rayleigth friction inside ugwp_driver  
  
    real, parameter    :: arad=6370.e3
    real, parameter    :: pi = atan(1.0)
    real, parameter    :: pi2 = 2.*pi  
    real, parameter    :: hps   = 7000.
    real, parameter    :: hpskm = hps/1000. 
!             
    real               ::  kxw = 6.28e-3/100.               ! single horizontal wavenumber of ugwp schemes  
    real, parameter    ::  ricrit = 0.25
    real, parameter    ::  frcrit = 0.50
    real, parameter    ::  linsat = 1.00
    real, parameter    ::  linsat2 = linsat*linsat     
!  

    integer               :: knob_ugwp_solver=1             ! 1, 2, 3, 4 - (linsat, ifs_2010, ad_gfdl, dsp_dis)
    integer, dimension(4) :: knob_ugwp_source               ! [1,1,1,0]  - (oro, fronts, conv, imbf-owp]
    integer, dimension(4) :: knob_ugwp_wvspec               !  number of waves for- (oro, fronts, conv, imbf-owp]
    integer, dimension(4) :: knob_ugwp_azdir                !   number of wave azimuths for- (oro, fronts, conv, imbf-owp]
    integer, dimension(4) :: knob_ugwp_stoch                !  1 - deterministic ; 0 - stochastic
    real,    dimension(4) :: knob_ugwp_effac                !  efficiency factors for- (oro, fronts, conv, imbf-owp]
    
    integer               :: knob_ugwp_doaxyz=1             ! 1 -gwdrag
    integer               :: knob_ugwp_doheat=1             ! 1 -gwheat
    integer               :: knob_ugwp_dokdis=0             ! 1 -gwmixing
    integer               :: knob_ugwp_ndx4lh = 2           ! n-number  of  "unresolved" "n*dx" for lh_gw 
!    
    integer  :: ugwp_azdir
    integer  :: ugwp_stoch

    integer  :: ugwp_src    
    integer  :: ugwp_nws 
    real     :: ugwp_effac    
     
!                   
    data knob_ugwp_source / 1,0, 1, 0 /                     !  oro-conv-fjet-okw-taub_lat:      1-active 0-off
    data knob_ugwp_wvspec /1,32,32,32/                      !  number of waves for- (oro, fronts, conv, imbf-owp, taulat]
    data knob_ugwp_azdir  /2, 4, 4,4/                       !  number of wave azimuths for- (oro, fronts, conv, imbf-okwp]
    data knob_ugwp_stoch  /0, 0, 0,0/                       !  0 - deterministic ; 1 - stochastic, non-activated option
    data knob_ugwp_effac  /1.,1.,1.,1./                     !  efficiency factors for- (oro, fronts, conv, imbf-owp] 
    integer  :: knob_ugwp_version = 0    
    namelist /cires_ugwp_nml/ knob_ugwp_solver, knob_ugwp_source,knob_ugwp_wvspec, knob_ugwp_azdir,  &
            knob_ugwp_stoch,  knob_ugwp_effac,knob_ugwp_doaxyz,  knob_ugwp_doheat, knob_ugwp_dokdis, &
            knob_ugwp_ndx4lh, knob_ugwp_version
!&cires_ugwp_nml 
! knob_ugwp_solver=2 
! knob_ugwp_source=1,1,1,0
! knob_ugwp_wvspec=1,32,32,32 
! knob_ugwp_azdir =2, 4, 4,4   
! knob_ugwp_stoch =0, 0, 0,0
! knob_ugwp_effac=1, 1, 1,1   
! knob_ugwp_doaxyz=1 
! knob_ugwp_doheat=1 
! knob_ugwp_dokdis=0 
! knob_ugwp_ndx4lh=4
!/ 
!
! allocatable arrays, initilized during "cires_ugwp_init" &
!                     released   during "cires_ugwp_finalize" 
!   	  
   real, allocatable :: kvg(:), ktg(:), krad(:), kion(:)
   real, allocatable :: zkm(:), pmb(:)   
   real, allocatable :: rfdis(:), rfdist(:) 
   integer           :: levs_rf
   real              :: pa_rf, tau_rf
!
! limiters
!   
   real, parameter ::  max_kdis = 400.              ! 400 m2/s
   real, parameter ::  max_axyz = 400.e-5           ! 400 m/s/day
   real, parameter ::  max_eps =  max_kdis*4.e-7    ! ~16   K/day    
!
!======================================================================   
   real, parameter :: F_coriol=1                    ! Coriolis effects
   real, parameter :: F_nonhyd=1                    ! Nonhydrostatic waves
   real, parameter :: F_kds   =0                    ! Eddy mixing due to GW-unstable below      
   real, parameter :: iPr_ktgw =1./3., iPr_spgw=iPr_ktgw 
   real, parameter :: iPr_turb =1./3., iPr_mol =1.95
   real, parameter ::  rhp1=1./hps, rhp2=0.5*rhp1, rhp4 = rhp2*rhp2
   real, parameter ::  khp =  0.287*rhp1            ! R/Cp/Hp 
   real, parameter ::  cd_ulim = 1.0                ! critical level precision or Lz ~ 0 ~dz of model   
      
   contains
!   
! -----------------------------------------------------------------------
!
! init  of cires_ugwp   (_init)  called from GFS_driver.F90
!
! -----------------------------------------------------------------------
subroutine cires_ugwp_mod_init (me, master, nlunit, logunit, fn_nml2, &
              lonr, latr, levs, ak, bk, pref, dtp, cdmvgwd, cgwf)
!
!  input_nml_file ='input.nml'=fn_nml 
!
    use  ugwp_oro_init,     only :  init_oro_gws
    use  ugwp_conv_init,    only :  init_conv_gws  
    use  ugwp_fjet_init,    only :  init_fjet_gws
    use  ugwp_okw_init,     only :  init_okw_gws
    use  ugwp_wmsdis_init,  only :  initsolv_wmsdis
    use  ugwp_lsatdis_init, only :  initsolv_lsatdis             
    implicit none
    
    integer, intent (in) :: me
    integer, intent (in) :: master
    integer, intent (in) :: nlunit
    integer, intent (in) :: logunit
    integer, intent (in) :: lonr
    integer, intent (in) :: levs
    integer, intent (in) :: latr   
    real,    intent (in) :: ak(levs+1), bk(levs+1), pref
    real,    intent (in) :: dtp   
    real,    intent (in) :: cdmvgwd(2), cgwf(2)             ! "scaling" controls for "old" GFS-GW schemes     
         
    character(len=64), intent (in) :: fn_nml2
    character(len=64), parameter   :: fn_nml='input.nml' 
       
!    character,  intent (in) :: input_nml_file
!    integer, parameter :: logunit =  6   
    integer :: ios
    logical :: exists
    real    :: dxsg 
    integer :: k
!       
    if (me == master) print *, trim (fn_nml), ' GW-namelist file '
    inquire (file =trim (fn_nml) , exist = exists)
!    
    if (.not. exists) then
       if (me == master) &
        write (6, *) 'separate ugwp :: namelist file: ', trim (fn_nml), ' does not exist'
    else
        open (unit = nlunit, file = trim(fn_nml), readonly, status = 'old', iostat = ios)
    endif
    rewind (nlunit)
    read   (nlunit, nml = cires_ugwp_nml)
    close  (nlunit)
!
    ! write version number and namelist to log file
    if (me == master) then
        write (logunit, *) " ================================================================== "
        write (logunit, *) "cires_ugwp_cires"
        write (logunit, nml = cires_ugwp_nml)
    endif
!
! effective kxw - resolution-aware
!    
    dxsg =  pi2*arad/float(lonr) * knob_ugwp_ndx4lh
!
!    kxw  =  pi2/dxsg
!
! init global background dissipation for ugwp -> 4d-variable for fv3wam linked with pbl-vert_diff 
!    

!    allocate(fcor(latr), fcor2(latr)  )
!
    allocate( kvg(levs+1),   ktg(levs+1)  ) 
    allocate( krad(levs+1),  kion(levs+1) )        
    allocate( zkm(levs),   pmb(levs) )  
    allocate( rfdis(levs), rfdist(levs) )      
!    
! ak -pa  bk-dimensionless  from surf => tol_lid_pressure =0
!    
    do k=1, levs 
       pmb(k) = 1.e0*(ak(k) + pref*bk(k))    ! Pa -unit  Pref = 1.e5
       zkm(k) = -hpskm*alog(pmb(k)/pref)
    enddo   
! 
! Part-1 :init_global_gwdis
!    
    call init_global_gwdis(levs, zkm, pmb, kvg, ktg, krad, kion)
    call rf_damp_init     (levs, pa_rf, tau_rf, dtp, pmb, rfdis, rfdist, levs_rf)   
! 
! Part-2 :init_SOURCES_gws
!    
    
!    
! call init-solver for "stationary" multi-wave spectra and sub-grid oro
!  
    call init_oro_gws( knob_ugwp_wvspec(1), knob_ugwp_azdir(1), &
         knob_ugwp_stoch(1), knob_ugwp_effac(1), lonr, kxw, cdmvgwd )
!    
! call init-sources for "non-sationary" multi-wave spectra 
!
    do_physb_gwsrcs=.true.
    
    IF (do_physb_gwsrcs) THEN

     if (me == master) print *, ' do_physb_gwsrcs ',  do_physb_gwsrcs, ' in cires_ugwp_init '
    if (knob_ugwp_wvspec(4) > 0) then 	
! okw     
       call init_okw_gws(knob_ugwp_wvspec(4), knob_ugwp_azdir(4), &
         knob_ugwp_stoch(4), knob_ugwp_effac(4), lonr, kxw )
       if (me == master) print *, ' init_okw_gws '
    endif	
    
    if (knob_ugwp_wvspec(3) > 0) then 
! fronts    	
       call init_fjet_gws(knob_ugwp_wvspec(3), knob_ugwp_azdir(3), &
         knob_ugwp_stoch(3), knob_ugwp_effac(3), lonr, kxw )
       if (me == master) print *, ' init_fjet_gws '	 
    endif 
    
     if (knob_ugwp_wvspec(2) > 0) then 
! conv    	
       call init_conv_gws(knob_ugwp_wvspec(2), knob_ugwp_azdir(2), &
         knob_ugwp_stoch(2), knob_ugwp_effac(2), lonr, kxw, cgwf )
       if (me == master)   &
           print *, ' init_convective GWs cgwf', knob_ugwp_wvspec(2), knob_ugwp_azdir(2)
       
    endif 
    
    ENDIF   !IF (do_physb_gwsrcs)

!======================    
! Part-3 :init_SOLVERS
! =====================       	      		      	 
	 
!
! call init-solvers for "broad" non-stationary multi-wave spectra
!    
    if   (knob_ugwp_solver==1) then
!   
      call initsolv_lsatdis(me, master, knob_ugwp_wvspec(2), knob_ugwp_azdir(2), &
           knob_ugwp_stoch(2), knob_ugwp_effac(2), do_physb_gwsrcs, kxw )
	  
    endif
    	   
     if   (knob_ugwp_solver==2) then 

       call initsolv_wmsdis(me, master, knob_ugwp_wvspec(2), knob_ugwp_azdir(2), &
         knob_ugwp_stoch(2), knob_ugwp_effac(2), do_physb_gwsrcs, kxw)
	 
     endif 
!
! other solvers not yet tested for fv3gfs
!     
!<    if   (knob_ugwp_solver==3) call init_dspdis
!<    if   (knob_ugwp_solver==4) call init_adodis
! 

!======================                
    module_is_initialized = .true.
    if (me == master)   print *, ' VAY-ugwp is initialized ', module_is_initialized

end subroutine cires_ugwp_mod_init

!   
!=====================================================================  
!
!ugwp-v0 subroutines: GWDPS_V0 and fv3_ugwp_solv2_v0    
!  
!=====================================================================
      SUBROUTINE GWDPS_V0(IM, levs, imx, do_tofd,   &
         Pdvdt, Pdudt, Pdtdt, Pkdis, U1,V1,T1,Q1,KPBL, &
         PRSI,DEL,PRSL,PRSLK,PHII, PHIL,DTP,KDT, &
         sgh30, HPRIME,OC,OA4,CLX4,THETA,vSIGMA,vGAMMA,ELVMAXD, &
         DUSFC, DVSFC,  xlatd, sinlat, coslat, sparea, &
         nmtvr, cdmbgwd, me, master, rdxzb, &
         zmtb, zogw, tau_mtb, tau_ogw, tau_tofd, &
         dudt_mtb, dudt_ogw, dudt_tms)
!       
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
      use ugwp_oro_init, only: dpmin, minwnd, hminmt,hncrit
      use ugwp_oro_init, only: sigfaci => sigfac    
      use ugwp_oro_init, only: RLOLEV,GMAX, VELEPS, FACTOP 
      use ugwp_oro_init, only: FRC, CE, CEOFRC, frmax, CG                  
      use ugwp_oro_init, only: FDIR, MDIR, NWDIR
      use ugwp_oro_init, only: cdmb, cleff, fcrit_gfs,fcrit_mtb 
      use ugwp_oro_init, only: n_tofd, ze_tofd, ztop_tofd    
                    
!      use cires_ugwp_module, only  : kxw,  max_kdis, max_axyz 
!      use sso_coorde,       only   : pgwd, pgwd4      
!----------------------------------------                  
      implicit none
      character(len=8)    :: strsolver='PSS-1986'  ! current operational solver or  'WAM-2017'
      integer, intent(in) :: im, levs, imx, kdt
      integer, intent(in) :: me, master, nmtvr
      integer             :: km
      logical, intent(in) :: do_tofd
      real(kind=kind_phys), parameter  :: sigfac =3, sigfacS = 0.5    
      real(kind=kind_phys)             :: ztopH,zlowH,ph_blk, dz_blk            
      integer, intent(in)              :: KPBL(IM)    ! Index for the PBL top layer!
      real(kind=kind_phys), intent(in) :: dtp         !  time step
      real(kind=kind_phys), intent(in) :: cdmbgwd(2) 
          
      real(kind=kind_phys), intent(in), dimension(im,levs) :: u1,  v1, t1,   q1, del,  prsl, prslk, phil   
      real(kind=kind_phys), intent(in),dimension(im,levs+1):: prsi, phii
      real(kind=kind_phys), intent(in) :: xlatd(im),sinlat(im),coslat(im)
      real(kind=kind_phys), intent(in) :: sparea(im)      
            
      real(kind=kind_phys), intent(in) :: OC(IM), OA4(im,4), CLX4(im,4) 
      real(kind=kind_phys), intent(in) :: HPRIME(IM), sgh30(IM)                     
      real(kind=kind_phys), intent(in) :: ELVMAXD(IM), THETA(IM)
      real(kind=kind_phys), intent(in) :: vSIGMA(IM),vGAMMA(IM)  
      real(kind=kind_phys)             :: SIGMA(IM),GAMMA(IM)        
      
!output -phys-tend
      real(kind=kind_phys),dimension(im,levs),intent(out):: Pdvdt,Pdudt
      real(kind=kind_phys),dimension(im,levs),intent(out):: Pkdis,Pdtdt
!     
!
! output          
! diag-coorde
      real(kind=kind_phys), dimension(im,levs), intent(out) :: dudt_mtb, dudt_ogw, dudt_tms
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
      real(kind=kind_phys) :: bgam, cgam, gam2, rnom, rdem       
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
      real(kind=kind_phys) :: XN(IM), YN(IM), UBAR(IM),    & 
                     VBAR(IM),  ULOW(IM),                  &
                     ROLL(IM),  bnv2bar(im), SCOR(IM),     &
                    DTFAC(IM), XLINV(IM), DELKS(IM), DELKS1(IM)
!
      real(kind=kind_phys) :: TAUP(IM,levs+1), TAUD(IM,levs)         
      real(kind=kind_phys) :: taub(im), taulin(im), heff, hsat, hdis       
         
      integer ::  kref(IM), idxzb(im), ipt(im), k_mtb,k_zlow 
      integer ::  kreflm(IM), iwklm(im), iwk(im), izlow(im)
      integer ::  ktrial, klevm1
!      
!check what we need
!
      real(kind=kind_phys) :: bnv,  fr, ri_gw , brvf,   tem,   tem1,  tem2, temc, temv, & 
                         ti,    rdz,   dw2,   shr2, bvf2, rdelks, efact, coefm, gfobnv,& 
                         scork,  rscor, hd,  fro,  sira, dtaux,  dtauy, pkp1log, pklog
     
      integer ::   kmm1, kmm2, lcap, lcapp1
      integer ::   npt       
      integer ::   kbps, kbpsp1,kbpsm1, kmps, idir, nwd,  klcap, kp1, kmpbl, kmll 
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
      real(kind=kind_phys)            :: shilmin, sgrmax, sgrmin
      real(kind=kind_phys)            :: belpmin, dsmin,  dsmax        
      real(kind=kind_phys), parameter :: nhilmax =  25.
      real(kind=kind_phys), parameter :: sso_min = 3000. 
      real(kind=kind_phys)            :: xlingfs
      real(kind=kind_phys)            :: arhills(im)       
      logical, parameter :: do_adjoro = .true.                     
! 
       integer  ::  i, j, k
       real(kind=kind_phys)     :: grav2, rcpdt, windik, wdir
       real(kind=kind_phys)     :: sigmin, dxres,sigres,hdxres
       real(kind=kind_phys)     :: cdmb4, mtbridge
       real(kind=kind_phys)     :: kxridge, inv_b2eff, zw1, zw2
       real(kind=kind_phys)     :: belps, aelps, nhills, selps             
       rcpdt = 1./cpd/dtp       
       grav2 = 2.*grav
!       
! mtb-blocking  sigma_min and dxres => cires_initialize
!  
       sgrmax = maxval(sparea) ; sgrmin = minval(sparea)
        dsmax = sqrt(sgrmax)   ; dsmin = sqrt(sgrmin)      

        dxres = pi2*arad/float(IMX)
        hdxres = .5*dxres   
        shilmin = sgrmin/nhilmax

    gammin = min(sso_min/dsmax, 1.)
    
        sigmin = 2.*hpmin/dsmax      !dxres 
       if (kdt==1) then 
    print *, sgrmax, sgrmin , ' min-max sparea '
    print *, 'sigmin-hpmin-dsmax',  sigmin, hpmin, dsmax
    print *, 'dxres/dsmax ', dxres, dsmax
    print *, ' shilmin gammin ', shilmin, gammin    
       endif
                  
        kxridge = float(IMX)/arad * cdmbgwd(2)
       
       if (me == master .and. kdt==1) then
        print *, ' gwdps_v0 kxridge ', kxridge
        print *, ' gwdps_v0 scale2 ', cdmbgwd(2)
        print *, ' gwdps_v0 IMX ', imx
        print *, ' gwdps_v0 GAM_MIN ',  gammin
    print *, ' gwdps_v0 SSO_MIN ', sso_min      
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
    sigma = vsigma
    gamma = vgamma
        do i = 1,im
    
          if ( (elvmaxd(i) .ge. hminmt)   .and. (gamma(i)   .ge. gammin)   .and. (hprime(i)  .ge. hpmin) ) then
          
             npt      = npt + 1
             ipt(npt) = i
         arhills(i) = 1.0
!
          if (gamma(i) .lt. gammin) gamma(i)= gammin         
          sigres = max(sigmin, sigma(i))
          if (sigma(i) .lt. sigmin) sigma(i)=  sigmin         
      dxres = sqrt(sparea(i))
       if (2.*hprime(i)/sigres.gt.dxres) sigres=2.*hprime(i)/dxres       
      aelps = min(2.*hprime(i)/sigres, 0.5*dxres)
      if (gamma(i) > 0.0 ) belps=min(aelps/gamma(i),.5*dxres) 
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
       
      selps = belps*belps*gamma(i)*4.    ! ellipse area of the el-c hill 
      nhills = sparea(i)/selps
      if (nhills > nhilmax) nhills=nhilmax
      arhills(i)= max(nhills, 1.0)
!333   format( ' nhil: ', I6, 4(2x, F9.3), 2(2x, E9.3))     
!      if (kdt==1 )
!     & write(6,333) nint(nhills)+1,xlatd(i), hprime(i),aelps*1.e-3,
!     &   belps*1.e-3, sigma(i),gamma(i)
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
         izlow(:) =1                       ! surface-level
        DO K = 1, levs-1
          DO I = 1, npt
            j = ipt(i)
        ztopH  = sigfac * hprime(j)
        zlowH  = sigfacs* hprime(j) 
            pkp1log =  phil(j,k+1) *rgrav
            pklog =    phil(j,k)   *rgrav
!        if (( ELVMAX(j).le.pkp1log) .and. (ELVMAX(j).ge.pklog) )
!     &      iwklm(I)  =  MAX(iwklm(I), k+1 ) 
        if (( ztopH.le.pkp1log) .and. (zTOPH.ge.pklog) )      iwklm(I)  =  MAX(iwklm(I), k+1 )
!     
        if (( zlowH.le.pkp1log) .and. (zlowH.ge.pklog) )      izlow(I)  =  MAX(izlow(I), k )      
          ENDDO
        ENDDO       
!
      DO K = 1,levs
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
          BVF2 = grav2 * RDZ * (VTK(I,K+1)-VTK(I,K))  / (VTK(I,K+1)+VTK(I,K))
          bnv2(i,k+1) = max( BVF2, bnv2min )
          RI_N(I,K+1) = Bnv2(i,k)/SHR2        ! Richardson number consistent with BNV2  
!
! add here computation for Ktur and OGW-dissipation fro VE-GFS
!       
        ENDDO
      ENDDO         
        K = 1 
        DO I = 1, npt
       bnv2(i,k) =bnv2(i,k+1) 
    ENDDO
!       
! level iwklm =>phil(j,k)/g < sigfac * hprime(j) < phil(j,k+1)/g 
!
        DO I = 1, npt
          J   = ipt(i)
      k_zlow = izlow(I)
      if (k_zlow == iwklm(i)) k_zlow = 1
          DELKS(I)  = 1.0 / (PRSI(J,k_zlow) - PRSI(J,iwklm(i)))
!          DELKS1(I) = 1.0 /(PRSL(J,k_zlow) - PRSL(J,iwklm(i)))
          UBAR (I)  = 0.0
          VBAR (I)  = 0.0
          ROLL (I)  = 0.0
          PE   (I)  = 0.0
          EK   (I)  = 0.0
          BNV2bar(I) = 0.0   
        ENDDO
!
        DO I = 1, npt     
      k_zlow = izlow(I)
      if (k_zlow == iwklm(i)) k_zlow = 1
          DO K = k_zlow, iwklm(I)-1                                !Kreflm(I)= iwklm(I)-1 
            J        = ipt(i)                                 ! laye-aver Rho, U, V
            RDELKS     = DEL(J,K) * DELKS(I)
            UBAR(I)    = UBAR(I)  + RDELKS * U1(J,K)          ! trial Mean U below 
            VBAR(I)    = VBAR(I)  + RDELKS * V1(J,K)          ! trial Mean V below 
            ROLL(I)    = ROLL(I)  + RDELKS * RO(I,K)          !   trial Mean RO below 
!       
            BNV2bar(I)=BNV2bar(I)+.5*(BNV2(I,K)+BNV2(I,K+1))* RDELKS
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
            if ( ANG(I,K) .gt.  90. ) ANG(I,K) = ANG(I,K) - 180.
            if ( ANG(I,K) .lt. -90. ) ANG(I,K) = ANG(I,K) + 180.
            ANG(I,K) = ANG(I,K) * DEG_TO_RAD
            UDS(I,K) =    MAX(SQRT(U1(J,K)*U1(J,K) + V1(J,K)*V1(J,K)), velmin)
!
            IF (IDXZB(I) .eq. 0 ) then
          dz_blk=( PHII(J,K+1) - PHII(J,K) ) *rgrav
              PE(I) = PE(I) + BNV2(I,K) *   ( ELVMAX(J) - phil(J,K)*rgrav ) * dz_blk

              UP(I)  =  max(UDS(I,K) * cos(ANG(I,K)), velmin)  
              EK(I)  = 0.5 *  UP(I) * UP(I) 
          
          ph_blk = ph_blk + dz_blk*sqrt(BNV2(I,K))/UP(I)

! --- Dividing Stream lime  is found when PE =exceeds EK. oper-l GFS
!              IF ( PE(I) .ge.  EK(I) ) THEN
              IF ( ph_blk .ge.  fcrit_gfs ) THEN
                 IDXZB(I) = K
         zmtb (J) =PHIL(J, K)*rgrav
                 RDXZB(J) =real(k, kind=kind_phys)       
              ENDIF

            ENDIF
          ENDDO
!
! Alternative expression: ZMTB = max(Heff*(1. -Fcrit_gfs/Fr), 0)
! fcrit_gfs/fr   
!
         goto 788 
      
     BNV    = SQRT( BNV2bar(I) )
         heff  = 2.*min(HPRIME(J),hpmax)
     zw2 = UBAR(I)*UBAR(I)+VBAR(I)*VBAR(I)
     Ulow(i)=sqrt(max(zw2,dw2min))
     Fr = heff*bnv/Ulow(i)
     ZW1 =max(Heff*(1. -fcrit_gfs/fr), 0.)  
     zw2 =  phil(j,2)*rgrav
     if (Fr > fcrit_gfs .and. zw1 > zw2 ) then 
       do k=2, levs-1
            pkp1log =  phil(j,k+1) *rgrav
            pklog =    phil(j,k)   *rgrav       
             if (zw1.le.pkp1log .and. zw1.ge.pklog)  exit         
       enddo
        IDXZB(I) =K
        zmtb (J) =PHIL(J, K)*rgrav
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
          IF ( IDXZB(I) .gt. 0 ) then
! (4.16)-IFS      
      gam2 = gamma(j)*gamma(j)
      BGAM = 1.0-0.18*gamma(j)-0.04*gam2
      CGAM =    0.48*gamma(j) +0.30*gam2      
            DO K = IDXZB(I)-1, 1, -1

                ZLEN = SQRT( ( PHIL(J,IDXZB(I)) - PHIL(J,K) ) /   ( PHIL(J,K ) + Grav * hprime(J) ) )

               COSANG2 = cos(ANG(I,K))*cos(ANG(I,K))
               SINANG2 =1.0 -COSANG2 
!         
!  cos =1 sin =0 =>   1/R= gam     ZR = 2.-gam 
!  cos =0 sin =1 =>   1/R= 1/gam   ZR = 2.- 1/gam
!
           rdem = COSANG2 +  GAM2 * SINANG2
           rnom = COSANG2*GAM2 +    SINANG2
!          
! metOffice Dec 2010
! correction of H. Wells & A. Zadra for the
! aspect ratio  of the hill seen by MF
! (1/R , R-inverse below: 2-R)

                 rdem = max(rdem, 1.e-6)       
                 R = sqrt(rnom/rdem)
                 ZR =  MAX( 2. - R, 0. )
                   
           sigres = max(sigmin, sigma(J))
       if (hprime(J)/sigres.gt.dxres) sigres=hprime(J)/dxres
           mtbridge = ZR * sigres*ZLEN / hprime(J)     
! (4.15)-IFS       
!           DBTMP = CDmb4 * mtbridge *
!     &           MAX(cos(ANG(I,K)), gamma(J)*sin(ANG(I,K))) 
! (4.16)-IFS     
           DBTMP = CDmb4*mtbridge*(bgam* COSANG2 +cgam* SINANG2)          
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
        
    kref(I)   = MAX(IWK(I), KPBL(J)+1 )             ! reference level PBL or smt-else ????
    kref(I)   = MAX(kref(i), iwklm(i) )             ! iwklm => sigfac*hprime
    
      if (kref(i) .lt. idxzb(i)) kref(i) = idxzb(i) + 1 ! layer above zmtb
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
          IF (K .LT. kref(I)) THEN
            J          = ipt(i)
            RDELKS     = DEL(J,K) * DELKS(I)
            UBAR(I)    = UBAR(I)  + RDELKS * U1(J,K)   ! Mean U below kref
            VBAR(I)    = VBAR(I)  + RDELKS * V1(J,K)   ! Mean V below kref
            ROLL(I)    = ROLL(I)  + RDELKS * RO(I,K)   ! Mean RO below kref
          BNV2bar(I)=BNV2bar(I)+.5*(BNV2(I,K)+BNV2(I,K+1))* RDELKS  
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
      DO  K = 1, levs-1
        DO  I = 1,npt
          J            = ipt(i)
          VELCO(I,K)   = 0.5 * ((U1(J,K)+U1(J,K+1))*XN(I) + (V1(J,K)+V1(J,K+1))*YN(I))
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
        J      = ipt(i)
        BNV    = SQRT( BNV2bar(I) )
        heff  = min(HPRIME(J),hpmax)
    
        if( zmtb(j) > 0.)heff=max(sigfac*heff-zmtb(j), 0.)/sigfac  
    if (heff .le. 0) cycle   
    
    hsat     = fcrit_gfs*ULOW(I)/bnv
    heff     = min(heff, hsat)
        
        FR = min(BNV  * heff /ULOW(I), FRMAX)
!
        EFACT    = (OA(I) + 2.) ** (CEOFRC*FR)
        EFACT    = MIN( MAX(EFACT,EFMIN), EFMAX )
!
        COEFM    = (1. + CLX(I)) ** (OA(I)+1.)
!
        XLINV(I) = COEFM * CLEFF           ! effective kxw for Lin-wave
    XLINGFS =  COEFM * CLEFF 
!
        TEM      = FR    * FR * OC(J)               
        GFOBNV   = GMAX  * TEM / ((TEM + CG)*BNV) 
!
!new specification of XLINV(I) & taulin(i)

        sigres = max(sigmin, sigma(J))
    if (heff/sigres.gt.hdxres) sigres=heff/hdxres
        inv_b2eff =  0.5*sigres/heff    
    kxridge =    1./sqrt(sparea(J))
        XLINV(I) = XLINGFS    !or max(kxridge, inv_b2eff)         ! 6.28/Lx ..0.5*sigma(j)/heff = 1./Lridge  
      taulin(i)=.5*ROLL(I)*XLINV(I)*BNV*ULOW(I)* heff*heff

       if ( FR > fcrit_gfs ) then
        TAUB(I)  = XLINV(I) * ROLL(I) * ULOW(I) * ULOW(I) * ULOW(I)  * GFOBNV  * EFACT       ! nonlinear FLUX Tau0...XLINV(I)
!
         else
!    
        TAUB(I)  = XLINV(I) * ROLL(I) * ULOW(I) * ULOW(I)  * ULOW(I)  * GFOBNV  * EFACT  
!     
!        TAUB(I)  = taulin(i)                             !  linear flux for FR <= fcrit_gfs
!   
        endif
!        
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
      
      if (strsolver=='PSS-1986') then     
       
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
          IF (K .GE. kref(I)) THEN
            ICRILV(I) = ICRILV(I) .OR. ( RI_N(I,K) .LT. RIC)   .OR. (VELCO(I,K) .LE. 0.0)
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
              BRVF = SQRT(BNV2(I,K))        ! Brent-Vaisala Frequency interface
!             TEM1 = XLINV(I)*(RO(I,KP1)+RO(I,K))*BRVF*VELCO(I,K)*0.5

              TEM1 = XLINV(I)*(RO(I,KP1)+RO(I,K))*BRVF*0.5  * max(VELCO(I,K), velmin)
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
              IF (RI_GW .LE. RIC .AND.  (OA(I) .LE. 0. .OR.  kp1 .ge. kref(i) )) THEN
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
!     Calculate wave acc-n: - (grav)*d(tau)/d(p) = taud
!
      DO K = 1,KM
        DO I = 1,npt
          TAUD(I,K)=GRAV*(TAUP(I,K+1) - TAUP(I,K))/DEL(ipt(I),K)
        ENDDO
      ENDDO
!
!------scale MOMENTUM DEPOSITION  AT TOP TO 1/2 VALUE
! it is zero now      
!         DO I = 1,npt
!            TAUD(I, levs) = TAUD(I,levs) * FACTOP
!         ENDDO
    
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
         IF (K .GE. kref(I) .and. PRSI(ipt(i),K) .GE. RLOLEV) THEN
         
             IF(TAUD(I,K).NE.0.) THEN
               TEM = DTP * TAUD(I,K)
               DTFAC(I) = MIN(DTFAC(I),ABS(VELCO(I,K)/TEM))
!          DTFAC(I) = 1.0
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
!          sigres = max(sigmin, sigma(J))
!      if (heff/sigres.gt.dxres) sigres=heff/dxres
!          inv_b2eff =  0.5*sigres/heff     
!        XLINV(I) = max(kxridge, inv_b2eff)           ! 0.5*sigma(j)/heff = 1./Lridge      
       dtfac(:) =  1.0
       
      call oro_wam_2017(im, levs, npt, ipt, kref, kdt, me, master,&
      dtp, dxres, taub, u1, v1, t1, xn, yn, bnv2, ro, prsi,prsL, &
        del, sigma, hprime, gamma, theta,&
        sinlat, xlatd, taup, taud, pkdis)
     
      endif            !  oro_wam_2017 - LINSATDIS-solver of WAM-2017
!      
!--------------------------- OROGW-solver of WAM2017
! 
! TOFD as in BELJAARS-2004
!
! ---------------------------    
      IF( do_tofd ) then
      axtms(:,:) = 0.0 ; aytms(:,:) = 0.0 
       if ( kdt == 1) then
           print *, 'VAY do_tofd  from surface to ', ztop_tofd 
       endif
       DO I = 1,npt          
     J = ipt(i)
     zpbl  =rgrav*phil( j, kpbl(j) )
     
         sigflt = min(sgh30(j), 0.3*hprime(j)) ! cannot exceed 30% of LS-SSO
     
     zsurf = phii(j,1)*rgrav
     zpm(1:levs) =phiL(j,1:levs)*rgrav
     up1(1:levs) = u1(j,1:levs)
     vp1(1:levs) = v1(j,1:levs)
     
         call ugwp_tofd1d(levs, sigflt, elvmaxd(j), zsurf, zpbl,& 
           up1, vp1, zpm,  utofd1, vtofd1, epstofd1, krf_tofd1)
     
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
       ENDIF    ! do_tofd 

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
! OGW-s above blocking height
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
      print *,zogw(J)/hprime(j), zmtb(j)/hprime(j), phil(j,1)/9.81, nint(hprime(j)/sigma(j))
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
!      call gwdrag_wam(1,  im,   ix,  levs,   ksrc, dtp,
!     & xlat, gw_dudt, gw_dvdt,  taux,  tauy)
!        call fv3_ugwp_wms17(kid1, im,  ix,  levs,  ksrc_ifs, dtp,
!     &  adt,adu,adv,prsl,prsi,phil,xlat, gw_dudt, gw_dvdt, gw_dtdt, gw_ked,
!     &  taux,tauy,grav, amol_i, me, lstep_first )
!
!
!23456==============================================================================

      subroutine slat_geos5_tamp(im, tau_amp, xlatdeg, tau_gw)
!=================
! GEOS-5 & MERRA-2 lat-dependent GW-source function  tau(z=Zlaunch) =rho*<u'w'>
!=================
      implicit none
      integer :: im
      real  :: xlatdeg(im),  tau_amp
      real  :: tau_gw(im)
      real  :: latdeg
!      real, parameter  :: tau_amp = 100.e-3
      real             :: trop_gw,  flat_gw
      integer :: i

!
! if-lat
      trop_gw = 0.75
      do i=1, im
      latdeg = xlatdeg(i)
      if (-15.3 < latdeg .and. latdeg < 15.3) then
          flat_gw = trop_gw*exp(-( (abs(latdeg)-3.)/8.0)**2)
          if (flat_gw < 1.2 .and. abs(latdeg) <= 3.) flat_gw = trop_gw
      else if (latdeg > -31. .and. latdeg <= -15.3) then
          flat_gw =  0.10
      else if (latdeg <  31. .and. latdeg >=  15.3) then
           flat_gw =  0.10
      else if (latdeg > -60. .and. latdeg <= -31.) then
        flat_gw =  0.50*exp(-((abs(latdeg)-60.)/23.)**2)
      else if (latdeg <  60. .and. latdeg >=  31.) then
        flat_gw =  0.50*exp(-((abs(latdeg)-60.)/23.)**2)
      else if (latdeg <= -60.) then
         flat_gw =  0.50*exp(-((abs(latdeg)-60.)/70.)**2)
      else if (latdeg >=  60.) then
         flat_gw =  0.50*exp(-((abs(latdeg)-60.)/70.)**2)
      end if
      tau_gw(i) = tau_amp*flat_gw
      enddo
!      
      end subroutine slat_geos5_tamp



      subroutine fv3_ugwp_solv2_v0(klon, klev, dtime, &
      tm1 , um1, vm1, qm1, &
      prsl, prsi,   philg, xlatd, sinlat, coslat, &
      pdudt, pdvdt, pdtdt, dked, tau_ngw, mpi_id, master, kdt)
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
        vueff(jl,jk) =  2.e-5*exp( (phil(jl,jk)+phil(jl,jk-1))*rhp2)+dked_min
!     
        zbn2(jl,jk)= grav2cpd/zthm1(jl,jk)*    (1.0 + rcpd*(tm1(jl,jk)-tm1(jl,jk-1))/zdelp)
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
         zul(jl,iazi)=zcosang(iazi)*zuhm1(jl,ilaunch)+  zsinang(iazi)*zvhm1(jl,ilaunch)
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
        zflux(jl,inc,1)=zfct(jl,ilaunch)* zbvfl4*zcin*zcpeak/(zbvfl4*zcpeak+zcin4*zcin)
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
               zdfl(jl,jk,iazi)=zdfl(jl,jk,iazi)+ zacc(jl,inc,iazi)*zflux(jl,inc,iazi)*zcinc
            enddo
           enddo
! --------------------------------------------
! get weighted average of phase speed in layer
! --------------------------------------------
          do jl=1,klon
            if(zdfl(jl,jk,iazi)>0.0 ) then
               zatmp=zcrt(jl,jk,iazi)
             do inc=1, nwav
                  zatmp=zatmp+zci(inc)*  zacc(jl,inc,iazi)*zflux(jl,inc,iazi)*zdci(inc)
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
       zdelp= grav/(prsi(jl,jk-1)-prsi(jl,jk))*  abs(zcin-zui(jl,jk,iazi)) *zcinc
       vm_zflx_mode=zact(jl,inc,iazi)* zflux_z(jl,inc,jk-1)

      if (vc_zflx_mode.gt.vm_zflx_mode) vc_zflx_mode =vm_zflx_mode ! no-flux increase
      zdfdz_v( jl,jk,iazi) =zdfdz_v( jl,jk,iazi) + (vm_zflx_mode-vc_zflx_mode)*zdelp                 ! heating >0            
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
! add limiters/efficiency for "unbalanced ics" if it is needed
!       
    pdudt = gw_eff*pdudt
    pdvdt = gw_eff*pdvdt
    pdtdt = gw_eff*pdtdt        
    dked =  gw_eff* dked                
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
      subroutine edmix_ugwp_v0(im, levs, dtp, & 
        t1, u1, v1, q1, del,       &
        prsl,  prsi, phil, prslk,&
        pdudt, pdvdt, pdTdt, pkdis,&
        ed_dudt, ed_dvdt, ed_dTdt,&
        me, master, kdt )
!
      use machine,           only : kind_phys
      use ugwp_common ,      only : rgrav, grav, cpd, rd, rdi, fv      
!      use ugwp_common ,      only : pi, rad_to_deg, deg_to_rad, pi2       
      use ugwp_common ,      only : bnv2min, velmin, arad
      
      implicit none
         
        integer, intent(in) :: me,  master, kdt
        integer, intent(in) :: im, levs 
    real(kind=kind_phys), intent(in) :: dtp 
        real(kind=kind_phys), intent(in), dimension(im,levs) ::  u1,  v1, t1,   q1, del,  prsl, prslk, phil  
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
      real(kind=kind_phys),parameter :: ulturb=150.,sc2u=ulturb* ulturb    
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
! adjust PT-profile to bn2(k) = bnv2min -- neutral atmosphere
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
        if ( stab .ge. 1.0 ) then
            stab_dt = max(stab_dt, stab)
        endif    
       enddo   
       nstab = max(1, nint(stab_dt)+1)
       dtstab = dtp /float(nstab)     
         ksum(levs) = ksum(levs-1)
         Fw(1:levs) =  pdudt(i, 1:levs)
         Fw1(1:levs) = pdvdt(i, 1:levs)      
     Km(1:levs) = ksum(1:levs) * rho(1:levs)* rho(1:levs)
     
     do j=1, nstab   
     !call diff_1d_wtend(levs, dtstab, Fw, Fw1, Km, del(i,:), Sw, Sw1)
     !WL*
      call diff_1d_wtend(levs, dtstab, Fw, Fw1, Km, rdp, rdpm, Sw, Sw1)
     !*WL
     Fw  = Sw
     Fw1 = Sw1   
     enddo
         
      ed_dudt(i,:) = Sw
      ed_dvdt(i,:) = Sw1
     
      Pt(1:levs) = t1(i,1:levs)*Ptmap(1:levs)      
      Kpt = Km*iPr_pt
          Fw(1:levs) =  pdTdt(i, 1:levs)*Ptmap(1:levs)
     do j=1, nstab                 
     ! call diff_1d_ptend(levs, dtstab, Fw, Kpt, del(i,:), Sw)
     !WL*
         call diff_1d_ptend(levs, dtstab, Fw, Kpt, rdp, rdpm, Sw)
     !*WL
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
      real(kind=kind_phys) :: Km1, Kp1, ad, cd, bd 
!      S(:) = 0.0 ; S1(:) = 0.0
!         
! explicit diffusion solver
!
      k = 1 
      km1 = 0. ; ad =0.
      kp1 = .5*(Km(k)+Km(k+1)) 
      cd = rdp(1)*rdpm(1)*kp1*dt
      bd = 1. - cd -ad
!      S(k) = cd*F(k+1) + ad *F(k-1) + bd *F(k)
      S(K)  = F(k)
      S1(K) = F1(k)      
      do k=2, levs-1
       ad = cd
       kp1 = .5*(Km(k)+Km(k+1))    
       cd = rdp(k)*rdpm(k)*kp1*dt
       bd = 1.-(ad +cd) 
       S(k) = cd*F(k+1) +   ad *F(k-1)  + bd *F(k)
       S1(k)= cd*F1(k+1)  + ad *F1(k-1) + bd *F1(k)     
      enddo 
       k =levs
       S(k) = F(k)
       S1(k) = F1(k)            
      end subroutine diff_1d_wtend
      
      subroutine diff_1d_ptend(levs, dt, F, Km, rdp, rdpm, S)
       use machine,      only: kind_phys
       implicit none
       integer :: levs
      real(kind=kind_phys) :: dt      
      real(kind=kind_phys) :: S(levs), S1(levs), F(levs), F1(levs) 
      real(kind=kind_phys) :: Km(levs),  rdp(levs), rdpm(levs-1)
      integer              :: i, k     
      real(kind=kind_phys) :: Km1, Kp1, ad, cd, bd 
!        
! explicit "eddy" smoother for tendencies
!           

      k = 1 
      km1 = 0. ; ad =0.
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
       k =levs
       S(k) = F(k)  
      end subroutine diff_1d_ptend  
        

     subroutine ugwp_tofd1d(levs, sigflt, elvmax, zsurf, zpbl,  u, v, zmid, utofd, vtofd, epstofd, krf_tofd)
       use machine ,          only : kind_phys
       use ugwp_common      , only :  rcpd2
       use ugwp_oro_init, only : n_tofd, const_tofd, ze_tofd, a12_tofd, ztop_tofd
!       
     implicit none
      integer ::  levs
      real(kind_phys), dimension(levs)        ::   u, v, zmid
      real(kind_phys)                         ::  sigflt, elvmax, zpbl, zsurf
      real(kind_phys), dimension(levs)        ::  utofd, vtofd, epstofd, krf_tofd     
!
! locals
!
     integer :: i, k
     real    :: sghmax = 5.
     real    :: sgh2, ekin
     real    :: zdec, rzdec
     real    :: umag, zmet, zarg, ztexp, krf
!     
     utofd =0.0 ; vtofd = 0.0 ;  epstofd =0.0 ; krf_tofd =0.0
!         
       zdec = max(n_tofd*sigflt, zpbl)             ! ntimes*sgh_turb or Zpbl
       zdec = min(ze_tofd, zdec)                ! cannot exceed 18 km
       rzdec =1./zdec
       sgh2 = max(sigflt*sigflt, sghmax*sghmax) ! 25 meters dz-of the first layer

      do k=1, levs
         zmet = zmid(k)-zsurf
      if (zmet > ztop_tofd) cycle
         ekin = u(k)*u(k) + v(k)*v(k)
         umag = sqrt(ekin)
         zarg = zmet*rzdec
         ztexp = exp(-zarg*sqrt(zarg))
         krf         = const_tofd* a12_tofd *sgh2* zmet ** (-1.2) *ztexp
         utofd(k)  = -krf*u(k)
         vtofd(k)  = -krf*v(k)
         epstofd(k)= rcpd2*krf*ekin  ! more accurate heat/mom form using "implicit tend-solver"
                                     ! to update momentum and temp-re; epstofd(k) can be skipped
         krf_tofd(k) = krf
     enddo
!                
     end subroutine ugwp_tofd1d

!===================================================  
!   
!part-1 init =>   wave dissipation + rfriction
!
!===================================================
     subroutine init_global_gwdis(levs, zkm, pmb, kvg, ktg, krad, kion)
     implicit none

     integer                              :: levs
     real, intent(in)                     :: zkm(levs), pmb(levs)
     real, intent(out), dimension(levs+1) :: kvg, ktg, krad, kion
!
!locals + data 
!  
     integer :: k
     real, parameter :: vusurf = 2.e-5
     real, parameter :: musurf = vusurf/1.95
     real, parameter :: hpmol  = 8.5
!      
     real, parameter :: kzmin  =  0.1
     real, parameter :: kturbo = 100.
     real, parameter :: zturbo = 130.
     real, parameter :: zturw  = 30.
     real, parameter :: inv_pra = 3.  !kt/kv =inv_pr     
!     
     real, parameter :: alpha  = 1./86400./15.
!     
     real, parameter :: kdrag  = 1./86400./10.
     real, parameter :: zdrag  = 100.
     real, parameter :: zgrow  = 50.
!
     real ::    vumol, mumol, keddy, ion_drag
!     
     do k=1, levs
      vumol =vusurf*exp(-zkm(k)/hpmol)
      mumol =musurf*exp(-zkm(k)/hpmol)

      keddy = kturbo*exp(-((zkm(k)-zturbo) /zturw)**2)

      kvg(k) = vumol +keddy
      ktg(k) = mumol +keddy*inv_pra

      krad(k) = alpha
!
      ion_drag = kdrag
!       
      kion(k) = ion_drag
     enddo

      k= levs+1
      kion(k) = kion(k-1)
      krad(k) = krad(k-1)
      kvg(k) =  kvg(k-1)
      ktg(k) =  ktg(k-1)
!                                                          
     end subroutine init_global_gwdis

     subroutine  rf_damp_init(levs, pa_rf, tau_rf, dtp, pmb, rfdis, rfdist, levs_rf)
     implicit none

     integer         ::   levs
     real            ::   pa_rf, tau_rf
     real            ::   dtp

     real            ::   pmb(levs)
     real            ::   rfdis(levs),  rfdist(levs)
     integer         ::   levs_rf

     real    ::            krf, krfz
     integer ::            k
!
     rfdis(1:levs) =  1.0
     rfdist(1:levs) = 0.0
     levs_rf = levs
     if (tau_rf <= 0.0 .or. pa_rf == 0.0) return

      krf = 1.0/(tau_rf*86400.0)

      do k=levs, 1, -1
        if(pmb(k) < pa_rf ) then                ! applied only on constant pressure surfaces fixed pmb in "Pa"
         krfz = krf*alog(pa_rf/pmb(k))
         rfdis(k)  =  1.0/(1.+krfz*dtp)
     rfdist(k)  =  (rfdis(k) -1.0)/dtp  ! du/dtp  
         levs_rf = k
        endif
      enddo
        
     end subroutine  rf_damp_init

      subroutine oro_wam_2017(im, levs,npt,ipt, kref,kdt,me,master, &
        dtp,dxres, taub, u1, v1, t1, xn, yn, bn2, rho, prsi, prsL, &
        del, sigma, hprime, gamma, theta, sinlat, xlatd, taup, taud, pkdis) 
!     
      USE MACHINE , ONLY          : kind_phys
      use ugwp_common ,      only : rgrav, grav, cpd, rd, rv, rcpd, rcpd2 
      use ugwp_common ,      only : rdi   
      use ugwp_common ,      only : pi, rad_to_deg, deg_to_rad, pi2 
      use ugwp_common ,      only : omega2, velmin, dw2min    
!      
      implicit none    
      
      integer :: im, levs
      integer :: npt
      integer :: kdt, me, master     
      integer :: kref(im), ipt(im)
      real(kind=kind_phys), intent(in) :: dtp, dxres 
      real(kind=kind_phys), intent(in) :: taub(im)
      
      real(kind=kind_phys), intent(in) :: sinlat(im), xlatd(im)
      real(kind=kind_phys), intent(in), dimension(im) :: sigma, hprime, gamma, theta
     
      real(kind=kind_phys), intent(in), dimension(im) :: xn, yn
      
      real(kind=kind_phys), intent(in), dimension(im, levs) :: u1, v1, t1,  bn2,  rho,   prsl, del
     
      real(kind=kind_phys), intent(in), dimension(im, levs+1) :: prsi
!
! out : taup,  taud, pkdis
!
      real(kind=kind_phys), intent(inout), dimension(im, levs+1) :: taup
      real(kind=kind_phys), intent(inout), dimension(im, levs)  :: taud
      real(kind=kind_phys), intent(inout), dimension(im, levs)  :: pkdis
      real(kind=kind_phys)     :: belps, aelps, nhills, selps
!
! multiwave oro-spectra
! locals
!
      integer :: i, j, k, isp, iw
      
      integer, parameter               :: nworo = 30     
      real(kind=kind_phys), parameter  :: fc_flag = 0.
      real(kind=kind_phys), parameter  :: mkzmin = 6.28e-3/50.
      real(kind=kind_phys), parameter  :: mkz2min = mkzmin* mkzmin    
      real(kind=kind_phys), parameter  :: kedmin = 1.e-3 
      real(kind=kind_phys), parameter  :: kedmax = 350.,axmax=250.e-5
      real(kind=kind_phys), parameter  :: rtau   = 0.01   ! nonlin-OGW scale 1/10sec  
      real(kind=kind_phys), parameter  :: Linsat2 =0.5 
      real(kind=kind_phys), parameter  :: kxmin = 6.28e-3/100.
      real(kind=kind_phys), parameter  :: kxmax = 6.28e-3/5.0
      real(kind=kind_phys), parameter  :: dkx = (kxmax -kxmin)/(nworo-1)
      real(kind=kind_phys), parameter  :: kx_slope= -5./3. 
      real(kind=kind_phys), parameter  :: hps =7000., rhp2 = .5/hps  
      real(kind=kind_phys), parameter  :: cxmin=0.5, cxmin2=cxmin*cxmin
                          
      real  :: akx(nworo),  cxoro(nworo), akx2(nworo)
      real  :: aspkx(nworo), c2f2(nworo) , cdf2(nworo)        
      real  :: tau_sp(nworo,levs+1), wkdis(nworo, levs+1)
      real  :: tau_kx(nworo),taub_kx(nworo)      
      real, dimension(nworo, levs+1)  :: wrms, akzw  
          
      real  :: tauz(levs+1), rms_wind(levs+1)    
      real  :: wave_act(nworo,levs+1) 
      
      real  :: kxw, kzw, kzw2, kzw3, kzi, dzmet, rhoint 
      real  ::  rayf, kturb                
      real  ::  uz, bv, bv2,kxsp, fcor2, cf2 
      
      real  ::  fdis       
      real  ::  wfdm, wfdt, wfim, wfit   
      real  ::  betadis, betam, betat, kds, cx, rhofac  
      real  ::  etwk, etws, tauk, cx2sat
      real  ::  cdf1, tau_norm     
!
! mean flow
!      
      real, dimension(levs+1) :: uzi,rhoi,ktur, kalp, dzi
      
      integer :: nw, nzi, ksrc
      taud (:, :) = 0.0 ; pkdis(:,:) = 0.0 ; taup (:,:) = 0.0
      tau_sp (:,:) = 0.0 ; wrms(:,:) = 0.0
      nw  =  nworo 
      nzi = levs+1  
       
      do iw = 1, nw
!                                 !kxw =  0.25/(dxres)*iw
     kxw = kxmin+(iw-1)*dkx 
     akx(iw)  = kxw
     akx2(iw) = kxw*kxw
     aspkx(iw) = kxw ** (kx_slope)
     tau_kx(iw) = aspkx(iw)*dkx  
      enddo
      
         tau_norm = sum(tau_kx) 
     tau_kx(:) =  tau_kx(:)/tau_norm 
      
    if (kdt == 1) then
771     format( 'vay-oro19 ', 3(2x,F8.3))   
        write(6,771)    maxval(tau_kx)*maxval(taub)*1.e3,   minval(tau_kx), maxval(tau_kx)
        endif         
!
! main loop over oro-points
!
      do i =1, npt      
         j = ipt(i)
    
!
! estimate "nhills" => stochastic choices for OGWs
!   
    if (taub(i) > 0.) then
!       
! max_kxridge =min( .5*sigma(j)/hprime(j), kmax)
! ridge-dependent dkx = (max_kxridge -kxmin)/(nw-1)
! option to make grid-box variable kx-spectra kxw = kxmin+(iw-1)*dkx
!
      wave_act(1:nw, 1:levs+1) = 1.0    
      ksrc = kref(i)    
      tauz(1:ksrc) = taub(i)    
      taub_kx(1:nw) =tau_kx(1:nw) * taub(i)
          wkdis(:,:) = kedmin   
     
    call oro_meanflow(levs, nzi, u1(j,:), v1(j,:), t1(j,:), prsi(j,:), &
       prsL(j,:),del(j,:),rho(i,:),bn2(i,:), uzi, rhoi,ktur, kalp,dzi,&
       xn(i), yn(i))    
     
      fcor2 = (omega2*sinlat(j))*(omega2*sinlat(j))*fc_flag

      k = ksrc
      
      bv2 =bn2(i,k)
      uz = uzi(k)           !u1(j,ksrc)*xn(i)+v1(j,ksrc)*yn(i)!
          kturb =ktur(k) 
      rayf  =kalp(k)
      rhoint = rhoi(k)
      dzmet = dzi(k)
          kzw =max(sqrt(bv2)/max(cxmin, uz), mkzmin) 
!             
! specify oro-kx spectra and related variables k=ksrc
!     
      do iw = 1, nw
        kxw =  akx(iw)      
        cxoro(iw) = 0.0 - uz 
        c2f2(iw) = fcor2/akx2(iw)       
            wrms(iw,k)= taub_kx(iw)/rhoint*kzw/kxw    
        tau_sp(iw, k) = taub_kx(iw)    
!           
!               
           if ( cxoro(iw) .gt. cxmin) then
            wave_act(iw,k:levs+1) = 0.                 ! crit-level     
        else
        cdf2(iw) =  cxoro(iw)*cxoro(iw) -c2f2(iw) 
         if ( cdf2(iw) .lt. cxmin2)  then 
              wave_act(iw,k:levs+1) = 0.               ! coriolis cut-off        
          else     
          kzw2 = max(Bv2/Cdf2(iw) - akx2(iw), mkz2min)
          kzw = sqrt(kzw2)
          akzw(iw,k)= kzw
          wrms(iw,k)= taub_kx(iw)/rhoint * kzw/kxw
         endif
        endif 
       enddo             !  nw-spectral loop 
!
! defined abobe, k = ksrc: akx(nworo),  cxoro(nworo), tau_sp(ksrc, nworo)
! propagate upward multiwave-spectra are filtered by dissipation & instability
!
!          tau_sp(:,ksrc+1:levs+1) = tau_sp(:, ksrc)                  
      do k= ksrc+1, levs
          uz = uzi(k)
          bv2 =bn2(i,k)
          bv = sqrt(bv2)
             rayf = kalp(k)
             rhoint= rhoi(k) 
             dzmet = dzi(k)
         rhofac = rhoi(k-1)/rhoi(k)       
         
       do iw = 1, nworo
!       
        if (wave_act(iw, k-1) .le. 0.) cycle
           cxoro(iw)= 0.-uz
             if ( cxoro(iw) .gt. cxmin) then
               wave_act(iw,k:levs+1) = 0.     ! crit-level
          else                  
             cdf2(iw) =  cxoro(iw)*cxoro(iw) -c2f2(iw)
         if ( cdf2(iw) .lt. cxmin2)  wave_act(iw,k:levs+1) = 0.
             endif       
             if ( wave_act(iw,k) .le. 0.) cycle 
!          
! upward propagation
!          
             kzw2 = Bv2/Cdf2(iw) - akx2(iw) 
               
             if (kzw2 .lt. mkz2min) then 
                 wave_act(iw,k:levs+1) =0.
             else  
!          
! upward propagation w/o reflection
!                    
           kxw = akx(iw)         
               kzw =sqrt(kzw2)
           akzw(iw,k)= kzw
               kzw3 =kzw2*kzw
         
         cx = cxoro(iw)
         betadis = cdf2(iw)/(Cx*Cx+c2f2(iw)) 
         betaM =1./(1.+betadis)
         betaT =1.- BetaM        
         kds  = wkdis(iw,k-1)
            
         etws =  wrms(iw,k-1)*rhofac * kzw/akzw(iw,k-1)

         kturb = ktur(k)+pkdis(j,k-1)
         wfiM = kturb*kzw2 +rayf
             wfiT = wfiM                   ! do updates with Pr-numbers Kv/Kt
         cdf1 = sqrt(Cdf2(iw))
         wfdM = wfiM/(kxw*Cdf1)*BetaM
             wfdT = wfiT/(kxw*Cdf1)*BetaT
         kzi = 2.*kzw*(wfdM+wfdT)*dzmet
         Fdis  = exp(-kzi)  
         
         etwk = etws*Fdis   
         Cx2sat = Linsat2*Cdf2(iw)
                          
       if (etwk .gt. cx2sat) then       
        Kds =  kxw*Cdf1*rhp2/kzw3
        etwk =cx2sat
        wfiM = kds*kzw2
        wfdM = wfiM/(kxw*Cdf1)
        kzi = 2.*kzw*(wfdm + wfdm)*dzmet
        etwk = cx2sat*exp(-kzi)     
       endif 
!                                if( lat(j) eq 40.5 ) then stop
        wkdis(iw,k) = kds
        wrms(iw,k) =  etwk 
            tauk = etwk*kxw/kzw 
        tau_sp(iw,k) = tauk *rhoint 
    if ( tau_sp(iw,k)>tau_sp(iw,k-1)) tau_sp(iw,k)=tau_sp(iw,k-1)  
    
        ENDIF    ! upward 
        ENDDO    ! spectral
    
!......... do spectral sum of rms, wkdis, tau    

        tauz(k)    = sum( tau_sp(:,k)*wave_act(:,k) )   
    rms_wind(k)= sum(   wrms(:,k)*wave_act(:,k) )
    
    pkdis(j,k)=sum(wkdis(:,k)*wave_act(:,k))+rms_wind(k)*rtau
            
    if (pkdis(j,k) .gt. kedmax) pkdis(j,k)= kedmax
    
        ENDDO   ! k=ksrc+1, levs
    
        k = ksrc
        tauz(k)      = sum(tau_sp(:,k)*wave_act(:,k))
    tauz(k)      = tauz(k+1)   ! zero momentum dep-n at k=ksrc 
    
    pkdis(j,k) =  sum(wkdis(:,k)*wave_act(:,k)) 
    rms_wind(k)  = sum(wrms(:,k)*wave_act(:,k)) 
          tauz(levs+1) = tauz(levs)
      taup(i, 1:levs+1) = tauz(1:levs+1)     
        do  k=ksrc, levs  
         taud(i,k) = ( tauz(k+1) - tauz(k))*grav/del(j,k)  
!    if (taud(i,k) .gt. 0)taud(i,k)=taud(i,k)*.01
!    if (abs(taud(i,k)).ge.axmax)taud(i,k)=sign(taud(i,k),axmax)    
    enddo       
      endif                    ! taub > 0
      enddo                    ! oro-points (i, j, ipt)
!23456      
      end subroutine oro_wam_2017

!-------------------------------------------------------------
!
! define mean flow  and dissipation for OGW-kx spectrum
!
!-------------------------------------------------------------      
      subroutine oro_meanflow(nz, nzi, u1, v1, t1, pint, pmid, &
       delp, rho, bn2, uzi, rhoi, ktur, kalp, dzi, xn, yn)
        
      use ugwp_common ,  only : grav, rgrav, rdi,  velmin, dw2min   
      implicit none
      
      integer :: nz, nzi
      real, dimension(nz  ) ::  u1,  v1, t1, delp, rho, pmid
      real, dimension(nz  ) ::  bn2  ! define at the interfaces
      real, dimension(nz+1) ::  pint
      real                  ::  xn, yn      
! output     
 
      real, dimension(nz+1) ::  dzi,  uzi, rhoi, ktur, kalp  
      
! locals      
      integer :: i, j, k
      real :: ui, vi, ti, uz, vz, shr2, rdz, kamp
      real :: zgrow, zmet, rdpm, ritur, kmol, w1
! paremeters
      real, parameter :: hps = 7000., rpspa = 1.e-5
      real, parameter :: rhps=1./hps
      real, parameter :: h4= 0.25/hps
      real, parameter :: rimin = 1./8., kedmin = 0.01
      real, parameter :: lturb = 30. ,  uturb = 150.    
      real, parameter :: lsc2 = lturb*lturb,usc2 = uturb*uturb         
      kalp(1:nzi) = 2.e-7        ! radiative damping
      
      do k=2, nz
        rdpm = grav/(pmid(k-1)-pmid(k)) 
        ui = .5*(u1(k-1)+u1(k))
        vi = .5*(v1(k-1)+v1(k))
    uzi(k) = Ui*xn + Vi*yn
        ti = .5*(t1(k-1)+t1(k))
    rhoi(k) = rdi*pint(k)/ti
    rdz = rdpm *rhoi(k)
    dzi(k) = 1./rdz
    uz = u1(k)-u1(k-1)
    vz = v1(k)-v1(k-1)
    shr2 = rdz*rdz*(max(uz*uz+vz*vz, dw2min))
    zmet = -hps*alog(pint(k)*rpspa)
    zgrow = exp(zmet*h4)    
    kmol =  2.e-5*exp(zmet*rhps)+kedmin
    ritur = max(bn2(k)/shr2, rimin)
        kamp = sqrt(shr2)*lsc2 *zgrow
     w1 = 1./(1. + 5*ritur)
     ktur(k)= kamp * w1 * w1 +kmol
      enddo 
       
      k=1 
      uzi(k) = uzi(k+1)
      ktur(k) = ktur(k+1)
      rhoi(k) = rdi*pint(k)/t1(k+1)
      dzi(k) = rgrav*delp(k)/rhoi(k)
      
      k = nzi
      uzi(k) = uzi(k-1)
      ktur(k) = ktur(k-1)
      rhoi(k) = rhoi(k-1)*.5
      dzi(k) =  dzi(k-1)   
         
      end subroutine oro_meanflow  

subroutine cires_ugwp_advance
!-----------------------------------------------------------------------
!
!   options for the day-to-day variable sources/spectra + diagnostics
!           for stochastic "triggers"
!   diagnose GW-source functions * FGF + OKWP + SGO/CONV from IAU-fields
!     or use for stochastic GWP-sources "memory"
!-----------------------------------------------------------------------
      implicit none
!
! update sources
!  a) physics-based triggers for multi-wave
!  b) stochastic-based spectra and amplitudes
!  c) use "memory" on GW-spectra from previous time-step
!  d) update "background" GW  dissipation as needed 
!      
end subroutine cires_ugwp_advance   
 
!      
! -----------------------------------------------------------------------
! finalize  of cires_ugwp   (_finalize) 
! -----------------------------------------------------------------------
subroutine cires_ugwp_mod_finalize
!
! deallocate sources/spectra & some diagnostics need to find where "deaalocate them"
! before "end" of the FV3GFS
!
    implicit none
!    
!   deallocate arrays employed in:
!     cires_ugwp_advance / cires_ugwp_driver / cires_ugwp_init
!
    deallocate( kvg,   ktg  ) 
    deallocate( krad,  kion )        
    deallocate( zkm,   pmb  )  
    deallocate( rfdis,   rfdist)   
            
end subroutine cires_ugwp_mod_finalize    
!



end module cires_ugwp_module


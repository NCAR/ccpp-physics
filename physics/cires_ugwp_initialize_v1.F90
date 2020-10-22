!===============================
! cu-cires ugwp-scheme
!    initialization of selected 
!  init gw-solvers (1,2,3,4)
!  init gw-source specifications
!  init gw-background dissipation
!==============================
!
! Part-0   specifications of common constants, limiters and "criiical" values
! 
!
    
    module ugwp_common_v1
!
!     use machine,  only : kind_phys
!     use physcons, only : pi => con_pi, grav => con_g, rd => con_rd,   &
!                          rv => con_rv, cpd => con_cp, fv => con_fvirt,&
!                          arad => con_rerth
     implicit none

      real, parameter ::  grav =9.81, cpd = 1004.
      real, parameter ::  rd = 287.0 , rv =461.5      
      real, parameter ::  grav2 = grav + grav
      real, parameter ::  rgrav = 1.0/grav, rgrav2= rgrav*rgrav
           
      real, parameter ::  fv   = rv/rd - 1.0  
      real, parameter ::  rdi  = 1.0 / rd, rcpd = 1./cpd, rcpd2 = 0.5/cpd  
      real, parameter ::  gor  = grav/rd
      real, parameter ::  gr2  = grav*gor
      real, parameter ::  grcp = grav*rcpd, gocp = grcp
      real, parameter ::  rcpdl    = cpd*rgrav                         ! 1/[g/cp]  == cp/g
      real, parameter ::  grav2cpd = grav*grcp                         !  g*(g/cp)= g^2/cp

      real, parameter ::  pi = 4.*atan(1.0), pi2 = 2.*pi, pih = .5*pi  
      real, parameter ::  rad_to_deg=180.0/pi, deg_to_rad=pi/180.0

      real, parameter ::  arad = 6370.e3
!               
      real, parameter ::  bnv2min = (pi2/1800.)*(pi2/1800.)
      real, parameter ::  bnv2max = (pi2/30.)*(pi2/30.)
       
      real, parameter :: dw2min=1.0,  velmin=sqrt(dw2min), minvel = 0.5     
      real, parameter :: omega1 = pi2/86400.
      real, parameter :: omega2 = 2.*omega1, omega3 = 3.*omega1     
      real, parameter :: hpscale= 7000., rhp=1./hpscale, rhp2=.5*rhp, rh4 = 0.25*rhp 
      real, parameter :: mkzmin = pi2/80.0e3, mkz2min = mkzmin*mkzmin
      real, parameter :: mkzmax = pi2/500.,  mkz2max = mkzmax*mkzmax 
      real, parameter :: cdmin  = 2.e-2/mkzmax   
     end module ugwp_common_v1
!
!
!===================================================
!
!Part-1 init =>   wave dissipation + RFriction
!
!===================================================
     subroutine init_global_gwdis_v1(levs, zkm, pmb, kvg, ktg, krad, kion, con_pi,  &
                                     pa_rf, tau_rf, me, master)

      
     implicit none
     integer , intent(in)                 :: me, master
     integer , intent(in)                 :: levs
     real, intent(in)                     :: con_pi, pa_rf, tau_rf    
     real, intent(in)                     :: zkm(levs), pmb(levs)    ! in km-Pa
     real, intent(out), dimension(levs+1) :: kvg, ktg, krad, kion
!
!locals + data
!
     integer :: k
     real, parameter :: vusurf  = 2.e-5
     real, parameter :: musurf  = vusurf/1.95
     real, parameter :: hpmol   = 8.5
!
     real, parameter :: kzmin   =  0.1
     real, parameter :: kturbo  = 100.
     real, parameter :: zturbo  = 130.
     real, parameter :: zturw   = 30.
     real, parameter :: inv_pra = 3.  !kt/kv =inv_pr
!
     real, parameter :: alpha  = 1./86400./15.   ! height variable see Zhu-1993 from 60-days => 6 days
     real  ::  pa_alp  = 750.                     ! super-RF parameters
     real  :: tau_alp  = 10.                     ! days   (750 Pa /10days)
!     
     real, parameter :: kdrag  = 1./86400./30.   !parametrization for WAM for FV3GFS SuperRF
     real, parameter :: zdrag  = 100.
     real, parameter :: zgrow  = 50.
!
     real ::    vumol, mumol, keddy, ion_drag
     real ::   rf_fv3, rtau_fv3, ptop, pih_dlog       
!
     real ::  ae1 ,ae2
     real :: pih

     pih = 0.5*con_pi

     pa_alp = pa_rf
     tau_alp = tau_rf
      
     ptop = pmb(levs)  
     rtau_fv3 = 1./86400./tau_alp
     pih_dlog = pih/log(pa_alp/ptop)

     do k=1, levs
      ae1 = -zkm(k)/hpmol
      vumol    = vusurf*exp(ae1)
      mumol    = musurf*exp(ae1)
      ae2 = -((zkm(k)-zturbo) /zturw)**2
      keddy    = kturbo*exp(ae2)

      kvg(k)   = vumol + keddy
      ktg(k)   = mumol + keddy*inv_pra

      krad(k)  = alpha
!
      ion_drag = kdrag
!
      kion(k)  = ion_drag!      
! add  Rayleigh_Super of FV3  for pmb < pa_alp
!
      if (pmb(k) .le. pa_alp) then
       rf_fv3=rtau_fv3*sin(pih_dlog*log(pa_alp/pmb(k)))**2
       krad(k) = krad(k) + rf_fv3 
       kion(k) = kion(k) + rf_fv3   

      endif

!      write(6,132) zkm(k), kvg(k), kvg(k)*(6.28/5000.)**2,  kion(k) 
     enddo

      k= levs+1
      kion(k) = kion(k-1)
      krad(k) = krad(k-1)
      kvg(k)  =  kvg(k-1)
      ktg(k)  =  ktg(k-1)
      if (me == master) then
	  write(6, * ) '  zkm(k), kvg(k), kvg(k)*(6.28/5000.)**2,  kion(k) '     
        do k=1, levs, 1
	  write(6,132) zkm(k), kvg(k), kvg(k)*(6.28/5000.)**2,  kion(k), pmb(k) 
	enddo      
      endif
!
 132  format( 2x, F8.3,' dis-scales:', 4(2x, E10.3))    
                                                          
     end subroutine init_global_gwdis_v1
!
!
     subroutine  rf_damp_init_v1(levs, pa_rf, tau_rf, dtp, pmb, rfdis, rfdist, levs_rf)
     implicit none

     integer         ::   levs
     real            ::   pa_rf, tau_rf
     real            ::   dtp

     real            ::   pmb(levs)
     real            ::   rfdis(levs),  rfdist(levs)
     integer         ::   levs_rf
 
     real            ::   krf, krfz
     integer         ::   k
!
     rfdis(1:levs)  = 1.0
     rfdist(1:levs) = 0.0
     levs_rf = levs
     if (tau_rf <= 0.0 .or. pa_rf == 0.0) return
 
      krf = 1.0/(tau_rf*86400.0)

      do k=levs, 1, -1
        if(pmb(k) < pa_rf ) then               ! applied only on constant pressure surfaces fixed pmb in "Pa"
         krfz = krf*log(pa_rf/pmb(k))
         rfdis(k)  =  1.0/(1.+krfz*dtp)
         rfdist(k) =  (rfdis(k) -1.0)/dtp      ! du/dtp
         levs_rf = k
        endif
      enddo

     end subroutine  rf_damp_init_v1
! ========================================================================
! Part 2 - sources
!      wave  sources
! ========================================================================
!
!    ugwp_oro_init_v1
!
!=========================================================================
     module ugwp_oro_init_v1

     use ugwp_common_v1, only : bnv2min, grav, grcp, fv, grav, cpd, grcp, pi
     use ugwp_common_v1, only : mkzmin, mkz2min
     implicit none
!  
! constants and "crirtical" values to run oro-mtb_gw physics
!
! choice of oro-scheme: strver = 'vay_2018' , 'gfs_2018', 'kdn_2005', 'smc_2000'
!
!
      real,      parameter :: hncrit=9000.   ! max value in meters for elvmax
      real,      parameter :: hminmt=50.     ! min mtn height (*j*)
      real,      parameter :: sigfac=4.0     ! mb3a expt test for elvmax factor
!
!
      real,      parameter :: minwnd=1.0     ! min wind component (*j*)
      real,      parameter :: dpmin=5000.0   ! minimum thickness of the reference layer in pa
      real,      parameter :: hpmax=2400.0, hpmin=25.0

      character(len=8)  ::  strver  = 'gfs_2018'
      character(len=8)  ::  strbase = 'gfs_2018'
      real, parameter   ::  rimin=-10., ric=0.25

!
      real, parameter :: efmin=0.5,    efmax=10.0
      

      real, parameter :: sigma_std=1./100., gamm_std=1.0

      real, parameter :: frmax=10., frc =1.0, frmin =0.01
!

       real, parameter :: ce=0.8,   ceofrc=ce/frc, cg=0.5
       real, parameter :: gmax=1.0, veleps=1.0, factop=0.5
!
       real, parameter :: rlolev=50000.0
!
 
 
!  hncrit set to 8000m and sigfac added to enhance elvmax mtn hgt


 
       real,      parameter ::   kxoro=6.28e-3/200.    !
       real,      parameter ::   coro = 0.0
       integer,   parameter ::   nridge=2
 
      real    ::  cdmb                      ! scale factors for mtb
      real    ::  cleff                     ! scale factors for orogw
      integer ::  nworo                     ! number of waves
      integer ::  nazoro                    ! number of azimuths
      integer ::  nstoro                    ! flag for stochastic launch above SG-peak

      integer, parameter ::  mdir = 8
      real,    parameter ::  fdir=.5*mdir/pi

      integer nwdir(mdir)
      data nwdir/6,7,5,8,2,3,1,4/
      save nwdir

      real,    parameter   ::   odmin = 0.1, odmax = 10.0
!------------------------------------------------------------------------------
!    small-scale orography parameters  for TOFD of Beljaars et al., 2004, QJRMS
!------------------------------------------------------------------------------

     integer, parameter  :: n_tofd = 2                      ! depth of SSO for TOFD compared with Zpbl
     real, parameter     :: const_tofd = 0.0759             ! alpha*beta*Cmd*Ccorr*2.109 = 12.*1.*0.005*0.6*2.109 = 0.0759
     real, parameter     :: ze_tofd    = 1500.0             ! BJ's z-decay in meters
     real, parameter     :: a12_tofd   = 0.0002662*0.005363 ! BJ's k-spect const for sigf2 * a1*a2*exp(-[z/zdec]**1.5]
     real, parameter     :: ztop_tofd  = 10.*ze_tofd        ! no TOFD > this height too higher 15 km
!------------------------------------------------------------------------------
!
      real, parameter :: fcrit_sm  = 0.7, fcrit_sm2  = fcrit_sm * fcrit_sm
      real, parameter :: fcrit_gfs = 0.7
      real, parameter :: fcrit_mtb = 0.7

      real,  parameter :: zbr_pi  = (1.0/2.0)*pi
      real,  parameter :: zbr_ifs = 0.5*pi

      contains
!
      subroutine init_oro_gws(nwaves, nazdir, nstoch, effac, &
                              lonr, kxw, cdmbgwd )
!
!
      integer :: nwaves, nazdir, nstoch
      integer :: lonr
      real    :: cdmbgwd(2)   ! scaling factors for MTb (1) & (2) for cleff = cleff  * cdmbgwd(2)
                              ! high res-n "larger" MTB and "less-active" cleff in GFS-2018   
      real    :: cdmbX
      real    :: kxw
      real    :: effac        ! it is analog of cdmbgwd(2) for GWs, off for now
!-----------------------------! GFS-setup for cdmb & cleff
! cdmb =  4.0 * (192.0/IMX)  
! cleff = 0.5E-5 / SQRT(IMX/192.0)  = 0.5E-5*SQRT(192./IMX)
!
      real, parameter :: lonr_refmb =  4.0 * 192.0
      real, parameter :: lonr_refgw =  192.0

! copy  to "ugwp_oro_init_v1"  =>  nwaves, nazdir, nstoch
 
      nworo  =  nwaves
      nazoro =  nazdir
      nstoro =  nstoch

      cdmbX = lonr_refmb/float(lonr)
      cdmb  = cdmbX
      if (cdmbgwd(1) >= 0.0) cdmb = cdmb * cdmbgwd(1)
 
       cleff = 0.5e-5 * sqrt(lonr_refgw/float(lonr))  !* effac

!!!    cleff = kxw    * sqrt(lonr_refgw/float(lonr))  !* effac

      if (cdmbgwd(2) >= 0.0) cleff = cleff  * cdmbgwd(2)
! 
!....................................................................
! higher res => smaller h' ..&.. higher kx
! flux_gwd ~ 'u'^2*kx/kz ~kxu/n ~1/dx *u/n    tau ~ h'*h'*kx*kx = const (h'-less kx-grow)
!....................................................................
!
!      print *, ' init_oro_gws 2-1cdmb',  cdmbgwd(2), cdmbgwd(1)
      end subroutine init_oro_gws
!

    end module ugwp_oro_init_v1
! =========================================================================
!
!    ugwp_conv_init_v1
!
!=========================================================================
    module ugwp_conv_init_v1

     implicit none
      real    ::  eff_con                   ! scale factors for conv GWs
      integer ::  nwcon                     ! number of waves
      integer ::  nazcon                    ! number of azimuths
      integer ::  nstcon                    ! flag for stochastic choice of launch level above Conv-cloud
      real    ::  con_dlength
      real    ::  con_cldf

      real, parameter    :: cmin  =  5  !2.5
      real, parameter    :: cmax  = 95. !82.5
      real, parameter    :: cmid  = 22.5
      real, parameter    :: cwid  = cmid
      real, parameter    :: bns   = 2.e-2, bns2 = bns*bns, bns4=bns2*bns2
      real, parameter    :: mstar = 6.28e-3/2.  !  2km
      real               :: dc

      real, allocatable  :: ch_conv(:),  spf_conv(:)
      real, allocatable  :: xaz_conv(:), yaz_conv(:)
     contains
!
     subroutine init_conv_gws(nwaves, nazdir, nstoch, effac, &
                              con_pi, arad, lonr, kxw, cgwf)

     implicit none
 
      integer :: nwaves, nazdir, nstoch
      integer :: lonr
      real    :: con_pi, arad
      real    :: cgwf(2)
      real    :: kxw,  effac
      real    :: work1 = 0.5
      real    :: chk, tn4, snorm
      integer :: k

      nwcon    = nwaves
      nazcon   = nazdir
      nstcon   = nstoch
      eff_con  = effac

      con_dlength = 2.0*con_pi*arad/float(lonr)
      con_cldf    = cgwf(1) * work1 + cgwf(2) *(1.-work1)
!
! allocate & define spectra in "selected direction": "dc" "ch(nwaves)"
!
       if (.not. allocated(ch_conv))  allocate (ch_conv(nwaves))
       if (.not. allocated(spf_conv)) allocate (spf_conv(nwaves))
       if (.not. allocated(xaz_conv)) allocate (xaz_conv(nazdir))
       if (.not. allocated(yaz_conv)) allocate (yaz_conv(nazdir))

      dc = (cmax-cmin)/float(nwaves-1)
!
! we may use different spectral "shapes"
! for example FVS-93 "Desabeius"
!  E(s=1, t=3,m, w, k) ~ m^s/(m*^4 + m^4) ~ m^-3 saturated tail
!
      do k = 1,nwaves
         chk         = cmin + (k-1)*dc
         tn4         = (mstar*chk)**4
         ch_conv(k)  =  chk
         spf_conv(k) =  bns4*chk/(bns4+tn4)
      enddo

      snorm = sum(spf_conv)
      spf_conv = spf_conv/snorm*1.5
 
      call init_nazdir(con_pi, nazdir,  xaz_conv,  yaz_conv)
     end subroutine init_conv_gws


    end module ugwp_conv_init_v1
!=========================================================================
!
!    ugwp_fjet_init_v1
!
!=========================================================================

   module ugwp_fjet_init_v1
      implicit none
      real    ::  eff_fj                     ! scale factors for conv GWs
      integer ::  nwfj                       ! number of waves
      integer ::  nazfj                      ! number of azimuths
      integer ::  nstfj                      ! flag for stochastic choice of launch level above Conv-cloud
!
      real, parameter    ::  fjet_trig=0.    ! if ( abs(frgf) > fjet_trig ) launch GW-packet

      
      real, parameter    :: cmin =  2.5
      real, parameter    :: cmax = 67.5
      real               :: dc
      real, allocatable  :: ch_fjet(:) , spf_fjet(:)
      real, allocatable  :: xaz_fjet(:), yaz_fjet(:)
     contains
     subroutine init_fjet_gws(nwaves, nazdir, nstoch, effac, &
                              con_pi, lonr, kxw)
     implicit none

      integer :: nwaves, nazdir, nstoch
      integer :: lonr
      real    :: con_pi
      real    :: kxw,  effac , chk
 
      integer :: k

      nwfj   =  nwaves
      nazfj  =  nazdir
      nstfj  =  nstoch
      eff_fj =  effac

       if (.not. allocated(ch_fjet))  allocate (ch_fjet(nwaves))
       if (.not. allocated(spf_fjet)) allocate (spf_fjet(nwaves))
       if (.not. allocated(xaz_fjet)) allocate (xaz_fjet(nazdir))
       if (.not. allocated(yaz_fjet)) allocate (yaz_fjet(nazdir))
 
      dc = (cmax-cmin)/float(nwaves-1)
      do k = 1,nwaves
         chk         = cmin + (k-1)*dc
         ch_fjet(k)  =  chk
         spf_fjet(k) =  1.0
      enddo
      call init_nazdir(con_pi, nazdir,  xaz_fjet,  yaz_fjet)

     end subroutine init_fjet_gws

    end module ugwp_fjet_init_v1
!
!=========================================================================
!
!
     module ugwp_okw_init_v1
!=========================================================================
     implicit none

      real    ::  eff_okw                     ! scale factors for conv GWs
      integer ::  nwokw                       ! number of waves
      integer ::  nazokw                      ! number of azimuths
      integer ::  nstokw                      ! flag for stochastic choice of launch level above Conv-cloud
!
      real, parameter    ::  okw_trig=0.      ! if ( abs(okwp) > okw_trig ) launch GW-packet

      real, parameter    :: cmin =  2.5
      real, parameter    :: cmax = 67.5
      real               :: dc
      real, allocatable  :: ch_okwp(:),   spf_okwp(:)
      real, allocatable  :: xaz_okwp(:),  yaz_okwp(:)

     contains
!
     subroutine init_okw_gws(nwaves, nazdir, nstoch, effac, &
                             con_pi, lonr, kxw)

     implicit none

      integer :: nwaves, nazdir, nstoch
      integer :: lonr
      real    :: con_pi
      real    :: kxw,  effac , chk
 
      integer :: k

      nwokw   =  nwaves
      nazokw  =  nazdir
      nstokw  =  nstoch
      eff_okw =  effac

       if (.not. allocated(ch_okwp))  allocate (ch_okwp(nwaves))
       if (.not. allocated(spf_okwp)) allocate (spf_okwp(nwaves))
       if (.not. allocated(xaz_okwp)) allocate (xaz_okwp(nazdir))
       if (.not. allocated(yaz_okwp)) allocate (yaz_okwp(nazdir))
      dc = (cmax-cmin)/float(nwaves-1)
      do k = 1,nwaves
         chk =  cmin + (k-1)*dc
         ch_okwp(k) = chk
         spf_okwp(k) = 1.
      enddo

      call init_nazdir(con_pi, nazdir,  xaz_okwp,  yaz_okwp)

     end subroutine init_okw_gws
 
     end module ugwp_okw_init_v1

!=============================== end of GW  sources
!
!  init specific  gw-solvers (1,2,3,4)
!

!===============================
!  Part -3  init  wave solvers
!===============================

  module ugwp_lsatdis_init_v1
     implicit none

      integer  :: nwav, nazd
      integer  :: nst
      real     :: eff
      integer, parameter  :: incdim = 4, iazdim = 4
!
     contains

     subroutine initsolv_lsatdis(me, master,  nwaves, nazdir, nstoch, effac, do_physb, kxw)

     implicit none
!
     integer  :: me, master
     integer  :: nwaves, nazdir
     integer  :: nstoch
     real     :: effac
     logical  :: do_physb
     real     :: kxw
!
!locals: define azimuths and Ch(nwaves) - domain when physics-based soureces
!                                          are not actibve
!
     integer :: inc, jk, jl, iazi, i, j, k

     if( nwaves == 0 .or. nstoch == 1 ) then
!                                redefine from the default
       nwav = incdim
       nazd = iazdim
       nst  = 0
       eff  = 1.0
     else
!                                from input_nml multi-wave spectra
       nwav = nwaves
       nazd = nazdir
       nst  = nstoch
       eff  = effac
     endif
!
     end subroutine initsolv_lsatdis
!
  end module ugwp_lsatdis_init_v1
!
!
  module ugwp_wmsdis_init_v1

     use ugwp_common_v1, only : arad, pi, pi2, hpscale, rhp, rhp2, rh4, omega2 
     use ugwp_common_v1, only : bnv2max,  bnv2min, minvel 
     use ugwp_common_v1, only : mkzmin,  mkz2min, mkzmax, mkz2max, cdmin
    implicit none

      real,     parameter   :: maxdudt = 250.e-5, maxdtdt=15.e-2
      real,     parameter   :: dked_min =0.01,    dked_max=250.0
   
      real,     parameter   :: gptwo=2.0
 
      real ,     parameter  :: bnfix  = pi2/300., bnfix2= bnfix * bnfix
      real ,     parameter  :: bnfix4 =  bnfix2 * bnfix2  
      real ,     parameter  :: bnfix3 =  bnfix2 * bnfix        
!
! make parameter list that will be passed to SOLVER
!
!     integer, parameter    :: klaunch=55      ! 32 - ~ 1km ;55 - 5.5 km ; 52 4.7km ; 60-7km index for selecting launch level
!     integer, parameter    :: ilaunch=klaunch
 
      integer  , parameter  :: iazidim=4       ! number of azimuths
      integer  , parameter  :: incdim=25       ! number of discrete cx - spectral elements in launch spectrum
      real     , parameter  :: ucrit=cdmin
 
      real ,     parameter  :: zcimin = 2.5
      real ,     parameter  :: zcimax = 125.0
      real ,     parameter  :: zgam   = 0.25
!
! Verical spectra
!
      real ,     parameter  :: pind_wd = 5./3.
      real ,     parameter  :: sind_kz = 1.
      real ,     parameter  :: tind_kz = 3.    
      real ,     parameter  :: stind_kz = sind_kz + tind_kz  
!
! from kmob_ugwp namelist
!      
      real                  :: nslope              ! the GW sprctral slope at small-m             
      real                  :: lzstar
      real                  :: lzmin
      real                  :: lzmax      
      real                  :: lhmet
      real                  :: tamp_mpa            !amplitude for GEOS-5/MERRA-2
      real                  :: tau_min             ! min of GW MF 0.25 mPa                        
      integer               :: ilaunch
      real                  :: gw_eff

      real                  :: v_kxw, rv_kxw, v_kxw2 
     
     
 
!===========================================================================
      integer  :: nwav, nazd, nst
      real     :: eff
 
      real                :: zaz_fct, zms
      real, allocatable   :: zci(:), zci4(:), zci3(:),zci2(:), zdci(:)
      real, allocatable   :: zcosang(:), zsinang(:)
      real, allocatable   :: lzmet(:), czmet(:), mkzmet(:), dczmet(:), dmkz(:) 

!
! GW-eddy constants for wave-mode dissipation by background and stability of
!         "final" flow after application of GW-effects
! 
      real, parameter :: iPr_pt = 0.5
      real, parameter :: lturb = 30., sc2  = lturb*lturb     ! stable on 80-km TL lmix ~ 500 met.
      real, parameter :: ulturb=150., sc2u = ulturb* ulturb  ! unstable
      real, parameter :: ric =0.25
      real, parameter :: rimin = -10., prmin = 0.25
      real, parameter :: prmax = 4.0
!
      contains
!============================================================================
     subroutine initsolv_wmsdis(me, master,  nwaves, nazdir, nstoch, effac, do_physb, kxw)
 
!        call initsolv_wmsdis(me, master, knob_ugwp_wvspec(2), knob_ugwp_azdir(2), &
!         knob_ugwp_stoch(2), knob_ugwp_effac(2), do_physb_gwsrcs, kxw)
!
     implicit none
!
!input -control for solvers:
!      nwaves, nazdir, nstoch, effac, do_physb, kxw
!
!
     integer  :: me, master, nwaves, nazdir, nstoch
     real     :: effac, kxw
     logical  :: do_physb
     real     :: dlzmet   
!
!locals
!
      integer :: inc, jk, jl, iazi
!
      real :: zang, zang1, znorm
      real :: zx1, zx2, ztx, zdx, zxran, zxmin, zxmax, zx, zpexp
      real :: fpc, fpc_dc
      real :: ae1,ae2
     if( nwaves == 0) then
!
!     redefine from the deafault
!
       nwav   = incdim
       nazd   = iazidim
       nst    = 0
       eff    = 1.0      
       gw_eff = eff
     else
!
! from input.nml
!
       nwav   = nwaves
       nazd   = nazdir
       nst    = nstoch
       gw_eff = effac
     endif  

      
      v_kxw  = pi2/lhmet ; v_kxw2 = v_kxw*v_kxw 
      rv_kxw = 1./v_kxw

      allocate ( zci(nwav),  zci4(nwav), zci3(nwav),zci2(nwav), zdci(nwav)  )
      allocate ( zcosang(nazd), zsinang(nazd) )
      allocate (lzmet(nwav), czmet(nwav), mkzmet(nwav), dczmet(nwav), dmkz(nwav) )   

      if (me == master) then
         print *, 'ugwp_v1: init_gw_wmsdis_control '
!  
         print *, 'ugwp_v1: WMS_DIS launch layer ',    ilaunch
         print *, 'ugwp_v1: WMS_DIS tot_mflux in mpa', tamp_mpa*1000.
         print *, 'ugwp_v1: WMS_DIS lhmet in km  ' , lhmet*1.e-3      	 
       endif

       zpexp = gptwo * 0.5                    ! gptwo=2 , zpexp = 1.

!
!      set up azimuth directions and some trig factors
!
!
       zang = pi2 / float(nazd)

! get normalization factor to ensure that the same amount of momentum
! flux is directed (n,s,e,w) no mater how many azimuths are selected.
!
       znorm = 0.0
       do iazi=1, nazd
         zang1         = (iazi-1)*zang
         zcosang(iazi) = cos(zang1)
         zsinang(iazi) = sin(zang1)
         znorm         = znorm + abs(zcosang(iazi))
       enddo
!      zaz_fct = 1.0
       zaz_fct = 2.0 / znorm            ! correction factor for azimuthal sums

!       define coordinate transform for "Ch"   ....x = 1/c stretching transform
!       ----------------------------------------------- 
! 
! x=1/Cphase transform
! see eq. 28-30 Scinocca 2003.   x = 1/c stretching transform
!
          zxmax = 1.0 / zcimin
          zxmin = 1.0 / zcimax
          zxran = zxmax - zxmin
          zdx   = zxran / real(nwav-1)                            ! dkz
!
          ae1=zxran/zgam
          zx1   = zxran/(exp(ae1)-1.0 )                    ! zgam =1./4.
          zx2   = zxmin - zx1

!
! computations for zci =1/zx, stretching "accuracy" is not "accurate" spectra transform
!  it represents additional "empirical" redistribution of "spectral" mode in C-space
!                               
         zms = pi2 / lzstar

        do inc=1, nwav
          ztx = real(inc-1)*zdx+zxmin
	  ae1 = (ztx-zxmin)/zgam
          zx  = zx1*exp(ae1)+zx2                            !eq.(29-30),Scinocca-2003
          zci(inc)  = 1.0 /zx                                            !
          zdci(inc) = zci(inc)**2*(zx1/zgam)*exp(ae1)*zdx   !
          zci4(inc) = (zms*zci(inc))**4
          zci2(inc) = (zms*zci(inc))**2
          zci3(inc) = (zms*zci(inc))**3
        enddo
!
!
!  alternatuve lzmax-lzmin
!
!
	dlzmet = (lzmax-lzmin)/ real(nwav-1)    
         do inc=1, nwav
          lzmet(inc) = lzmin + (inc-1)*dlzmet
	  mkzmet(inc) = pi2/lzmet(inc)
          zci(inc) =lzmet(inc)/(pi2/bnfix)
          zci4(inc) = (zms*zci(inc))**4
          zci2(inc) = (zms*zci(inc))**2
          zci3(inc) = (zms*zci(inc))**3	
	  
	 enddo

	  zdx = (zci(nwav)-zci(1))/ real(nwav-1)  
    

         if (me == master) then
           print *
           print *,  'ugwp_v0: zcimin=' , zcimin
           print *,  'ugwp_v0: zcimax=' , zcimax
           print *,  'ugwp_v0: zgam=   ', zgam                
           print *
           
!          print *, ' ugwp_v1  nslope=', nslope
           print *
           print *,  'ugwp_v1: zcimin/zci=' , maxval(zci)
           print *,  'ugwp_v1: zcimax/zci=' , minval(zci)
           print *,  'ugwp_v1: cd_crit=', ucrit            
           print *,  'ugwp_v1: launch_level',  ilaunch
           print *, ' ugwp_v1  lzstar=', lzstar
           print *, ' ugwp_v1  nslope=', nslope

           print *
        do inc=1, nwav	    	  
          zdci(inc) = zdx 
          if (nslope == 1) fpc = bnfix4*zci(inc)/ (bnfix4+zci4(inc))
          if (nslope == 0) fpc = bnfix3*zci(inc)/ (bnfix3+zci3(inc))	  
          fpc_dc = fpc * zdci(inc)
	  write(6,111)  inc, zci(inc), zdci(inc),ucrit, fpc, fpc_dc, 6.28e-3/bnfix*zci(inc)
        enddo
         endif
 111     format( 'wms-zci', i4, 7 (3x, F8.3))

     end subroutine initsolv_wmsdis
!
! make a list of  all-initilized parameters needed for "gw_solver_wmsdis"
!

  end module ugwp_wmsdis_init_v1
!=========================================================================
!
! work TODO for 2-extra WAM-solvers:
!           DSPDIS (Hines)+ADODIS (Alexander-Dunkerton-Ortland)
!
!========================================================================= 
     subroutine init_dspdis_v1
     implicit none
     end subroutine init_dspdis_v1

     subroutine init_adodis_v1
     implicit none
     end subroutine init_adodis_v1
     

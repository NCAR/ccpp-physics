!===============================
! cu-cires ugwp-scheme
!  initialization of selected 
!  init gw-solvers (1,2,3,4)
!  init gw-source specifications
!  init gw-background dissipation
!==============================
!
! Part-0   specifications of common constants, limiters and "criiical" values
! 
!
    
    module ugwp_common
!
     use machine,  only : kind_phys
			 
     implicit none
     
      real(kind=kind_phys)   ::  pi, pi2, pih, rad_to_deg, deg_to_rad
      real(kind=kind_phys)   ::  arad, p0s     
      real(kind=kind_phys)   ::  grav, grav2, rgrav, rgrav2
      real(kind=kind_phys)   ::  cpd,  rd, rv, fv            
      real(kind=kind_phys)   ::  rdi,  rcpd, rcpd2  
      
      real(kind=kind_phys)   ::  gor,  gr2,  grcp, gocp, rcpdl, grav2cpd
      real(kind=kind_phys)   ::  bnv2min, bnv2max
      real(kind=kind_phys)   ::  dw2min,  velmin, minvel        
      real(kind=kind_phys)   ::  omega1, omega2,   omega3   
      real(kind=kind_phys)   ::  hpscale, rhp, rhp2, rh4, rhp4, khp, hpskm
      real(kind=kind_phys)   ::  mkzmin, mkz2min,  mkzmax, mkz2max, cdmin    
      real(kind=kind_phys)   ::  rcpdt      
           
!      real(kind=kind_phys), parameter ::  grav2 = grav + grav
!      real(kind=kind_phys), parameter ::  rgrav = 1.0/grav, rgrav2= rgrav*rgrav           
!      real(kind=kind_phys), parameter ::  rdi  = 1.0 / rd, rcpd = 1./cpd, rcpd2 = 0.5/cpd  
!      real(kind=kind_phys), parameter ::  gor  = grav/rd, rcpdt = 1./(cp*dtp)
      
!      real(kind=kind_phys), parameter ::  gr2  = grav*gor
!      real(kind=kind_phys), parameter ::  grcp = grav*rcpd, gocp = grcp
!      real(kind=kind_phys), parameter ::  rcpdl    = cpd*rgrav                         ! 1/[g/cp]  == cp/g
!      real(kind=kind_phys), parameter ::  grav2cpd = grav*grcp                         !  g*(g/cp)= g^2/cp
!      real(kind=kind_phys), parameter ::  pi2 = 2.*pi, pih = .5*pi  
!      real(kind=kind_phys), parameter ::  rad_to_deg=180.0/pi, deg_to_rad=pi/180.0
!               
!      real(kind=kind_phys), parameter ::  bnv2min = (pi2/1800.)*(pi2/1800.)
!      real(kind=kind_phys), parameter ::  bnv2max = (pi2/30.)*(pi2/30.)      
!      real(kind=kind_phys), parameter :: dw2min=1.0,  velmin=sqrt(dw2min), minvel = 0.5     
!      real(kind=kind_phys), parameter :: omega1 = pi2/86400., omega2 = 2.*omega1, omega3 = 3.*omega1     
!   
!      real(kind=kind_phys), parameter :: hpscale= 7000., rhp=1./hpscale, rhp2=.5*rhp, rh4 = 0.25*rhp 
!      real(kind=kind_phys), parameter :: mkzmin = pi2/80.0e3, mkz2min = mkzmin*mkzmin
!      real(kind=kind_phys), parameter :: mkzmax = pi2/500.,   mkz2max = mkzmax*mkzmax 
!      real(kind=kind_phys), parameter :: cdmin  = 2.e-2/mkzmax 
!      real(kind=kind_phys), parameter ::  pi = 4.*atan(1.0),
!      real(kind=kind_phys), parameter ::  grav =9.81, cpd = 1004.
!      real(kind=kind_phys), parameter ::  rd = 287.0 , rv =461.5  
!      real(kind=kind_phys), parameter ::  fv   = rv/rd - 1.0    
!      real(kind=kind_phys), parameter ::  arad = 6370.e3
        
     end module ugwp_common
     
      subroutine init_nazdir(naz,  xaz,  yaz)
      
      use machine,     only : kind_phys
      use ugwp_common, only :  pi2
      
      implicit none
      
      integer :: naz
      real(kind=kind_phys), dimension(naz) :: xaz,  yaz
      integer :: idir
      real(kind=kind_phys)    :: phic, drad
      
      drad  = pi2/float(naz)
      if (naz.ne.4) then     
        do idir =1, naz
         Phic = drad*(float(idir)-1.0)
         xaz(idir) = cos(Phic)
         yaz(idir) = sin(Phic)
        enddo
      else 
!                                  if (naz.eq.4) then
          xaz(1) = 1.0     !E
          yaz(1) = 0.0
          xaz(2) = 0.0     
          yaz(2) = 1.0     !N
          xaz(3) =-1.0     !W
          yaz(3) = 0.0
          xaz(4) = 0.0
          yaz(4) =-1.0     !S
      endif      
      end  subroutine init_nazdir     
!
!
!===================================================
!
!Part-1 init =>   wave dissipation + RFriction
!
!===================================================
     subroutine init_global_gwdis(levs, zkm, pmb, kvg, ktg, krad, kion, me, master)
!				      
     use machine ,      only : kind_phys 
     use ugwp_common,   only : pih, pi  
     
     implicit none
     integer , intent(in)                 :: me, master
     integer , intent(in)                 :: levs
     real(kind=kind_phys), intent(in)                     :: zkm(levs), pmb(levs)    ! in km-Pa
     real(kind=kind_phys), intent(out), dimension(levs+1) :: kvg, ktg, krad, kion
!
!locals + data
!
     integer :: k
     real(kind=kind_phys), parameter :: vusurf  = 2.e-5
     real(kind=kind_phys), parameter :: musurf  = vusurf/1.95
     real(kind=kind_phys), parameter :: hpmol   = 7.0
!
     real(kind=kind_phys), parameter :: kzmin   =  0.1
     real(kind=kind_phys), parameter :: kturbo  = 100.
     real(kind=kind_phys), parameter :: zturbo  = 130.
     real(kind=kind_phys), parameter :: zturw   = 30.
     real(kind=kind_phys), parameter :: inv_pra = 3.  !kt/kv =inv_pr
!
     real(kind=kind_phys), parameter :: alpha  = 1./86400./15.   ! height variable see Zhu-1993 from 60-days => 6 days
     real(kind=kind_phys)  ::  pa_alp  = 750.                    ! super-RF parameters from FV3-dycore GFSv17/16 sett
     real(kind=kind_phys)  :: tau_alp  = 10.                     ! days   (750 Pa /10days)
!     
     real(kind=kind_phys), parameter :: kdrag  = 1./86400./30.   !parametrization for WAM ion drag as e-density function
     real(kind=kind_phys), parameter :: zdrag  = 100.
     real(kind=kind_phys), parameter :: zgrow  = 50.
!
     real(kind=kind_phys) ::   vumol, mumol, keddy, ion_drag
     real(kind=kind_phys) ::   rf_fv3, rtau_fv3, ptop, pih_dlog       
!
     real(kind=kind_phys) ::  ae1 ,ae2
!     
      
     ptop = pmb(levs)  
     rtau_fv3 = 1./86400./tau_alp
     pih_dlog = pih/log(pa_alp/ptop)

     do k=1, levs
      ae1 =    zkm(k)/hpmol
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
      
!      if (me == master) then
!	  write(6, * ) '  zkm(k), kvg(k), kvg(k)*(6.28/5000.)**2,  kion(k) ... init_global_gwdis'     
!        do k=1, levs, 1
!	  write(6,132) zkm(k), kvg(k), kvg(k)*(6.28/5000.)**2,  kion(k), pmb(k) 
!	enddo      
!      endif
!
! 132  format( 2x, F8.3,' dis-scales:', 4(2x, E10.3))    
                                                          
     end subroutine init_global_gwdis
!
! ========================================================================
! Part 2 - sources
!      wave  sources
! ========================================================================
!
!    ugwp_oro_init
!
!=========================================================================
     module ugwp_oro_init
     use machine ,    only : kind_phys
     use ugwp_common, only : bnv2min, grav, grcp, fv, grav, cpd, grcp, pi
     use ugwp_common, only : mkzmin, mkz2min
     implicit none
!  
! constants and "crirtical" values to run oro-mtb_gw physics
!
!
      real(kind=kind_phys),      parameter :: hncrit=9000.   ! max value in meters for elvmax
      real(kind=kind_phys),      parameter :: hminmt=50.     ! min mtn height (*j*)
      real(kind=kind_phys),      parameter :: sigfac=4.0     ! mb3a expt test for elvmax factor
      real(kind=kind_phys),      parameter :: hpmax=2500.0
      real(kind=kind_phys),      parameter :: hpmin=25.0      
!
!
      real(kind=kind_phys),      parameter :: minwnd=1.0     ! min wind component (*j*)
      real(kind=kind_phys),      parameter :: dpmin=5000.0   ! minimum thickness of the reference layer in pa

      
      character(len=8)  ::  strver  = 'gfs_2018'
      character(len=8)  ::  strbase = 'gfs_2018'
      
      real(kind=kind_phys), parameter   ::  rimin=-10., ric=0.25
      
      real(kind=kind_phys), parameter   ::  frmax=10., frc =1.0, frmin =0.01
      real(kind=kind_phys), parameter   ::  ce=0.8,   ceofrc=ce/frc, cg=0.5
      real(kind=kind_phys), parameter   ::  gmax=1.0, veleps=1.0, factop=0.5!
      real(kind=kind_phys), parameter   ::  efmin=0.5,    efmax=10.0
      
      real(kind=kind_phys), parameter   ::  rlolev=50000.0
      integer,              parameter   ::  mdir = 8
      real(kind=kind_phys), parameter   ::  fdir=mdir/(8.*atan(1.0))
      real(kind=kind_phys), parameter   ::  zpgeo=2.*atan(1.0)
      
      integer nwdir(mdir)
      data nwdir/6,7,5,8,2,3,1,4/
      save nwdir
      
      real(kind=kind_phys), parameter   :: odmin = 0.1, odmax = 10.0
      real(kind=kind_phys), parameter   :: fcrit_sm  = 0.7, fcrit_sm2  = fcrit_sm * fcrit_sm
      real(kind=kind_phys), parameter   :: fcrit_gfs = 0.7, fcrit_v1 = 0.7
      real(kind=kind_phys), parameter   :: fcrit_mtb = 0.7

      real(kind=kind_phys),  parameter  :: zbr_pi  = zpgeo
      real(kind=kind_phys),  parameter  :: zbr_ifs = zpgeo
 
!  
 
      real(kind=kind_phys),   parameter ::   kxoro=6.28e-3/200.    !
      real(kind=kind_phys),   parameter ::   coro = 0.0
      integer,parameter ::    nridge=2
      real(kind=kind_phys),   parameter ::   sigma_std=1./100., gamm_std=1.0       
 
      real(kind=kind_phys)    ::  cdmb                      ! scale factors for mtb
      real(kind=kind_phys)    ::  cleff                     ! scale factors for orogw
      
      integer ::  nworo                     ! number of waves
      integer ::  nazoro                    ! number of azimuths
      integer ::  nstoro                    ! flag for stochastic launch above SG-peak

      
!------------------------------------------------------------------------------
!    small-scale orography parameters  for TOFD of Beljaars et al., 2004, QJRMS
!    SA-option can be controlled by  Integral limits of fluxes
!    in B2004: klow = 0.003 1/m ~ 2km and kinf ~ 6.28/10/(Z1)~< 1 km => meters
!    these limits can change strength of TOFD... choice of k0tr ~1/10 km (10km ~dx of C768)
!    kmax = kdis_pbl
!------------------------------------------------------------------------------
     real(kind=kind_phys), parameter     :: kmax = 6.28/(10.*25.)           ! max k-tofd
     real(kind=kind_phys), parameter     :: k1tr = 6.28/(2100)              ! max k-transition from -1.9/slope to -2.8/slope
     real(kind=kind_phys), parameter     :: kflt = 6.28/(18.e3)             !
     real(kind=kind_phys), parameter     :: k0tr = 6.28/(10.e3)             ! min k-tofd     
     real(kind=kind_phys), parameter     :: nk1tr = 2.8
     real(kind=kind_phys), parameter     :: nk0tr = 1.9
     real(kind=kind_phys), parameter     :: a1_tofd =  kflt ** nk1tr *1.e3
     real(kind=kind_phys), parameter     :: a2_tofd =  k1tr ** (nk0tr-nk1tr)    
     real(kind=kind_phys), parameter     :: fix_tofd = 2.* 0.005 * 12 *0.6     !value= 0.072
!     
! B2004 scheme is based on the empirical vertical profile of the tofd divergence: 
!       Ax_tofd(Z)=exp(-[Z/ze_tofd]^3/2) / Z^1.2.....
! TOFD-flux/TMS-flux must dissipate due to PBL-diffusion with spectral damping
! Here we can enhance TOFD-impact by selecting k0tr and kmax limits
!      as functions of resolution and PBL-dissipation
!                            
     integer, parameter  :: n_tofd = 2                                      ! depth of SSO for TOFD compared with Zpbl
     real(kind=kind_phys), parameter     :: const_tofd = 0.0759             ! alpha*beta*Cmd*Ccorr*2.109 = 12.*1.*0.005*0.6*2.109 = 0.0759
     real(kind=kind_phys), parameter     :: ze_tofd    = 1500.0             ! BJ's z-decay in meters, 1.5 km
     real(kind=kind_phys), parameter     :: a12_tofd   = 0.0002662*0.005363 ! BJ's k-spect const for sigf2 * a1*a2*exp(-[z/zdec]**1.5]
     real(kind=kind_phys), parameter     :: ztop_tofd  = 3.*ze_tofd         ! no TOFD > this height 4.5 km
!------------------------------------------------------------------------------
!

      contains
!
      subroutine init_oro_gws(nwaves, nazdir, nstoch, effac, &
                              lonr, kxw )
!
!
      integer :: nwaves, nazdir, nstoch
      integer :: lonr			       
      real(kind=kind_phys)    :: cdmbX
      real(kind=kind_phys)    :: kxw
      real(kind=kind_phys)    :: effac        ! it is analog of cdmbgwd(2) for GWs, off for now
!-----------------------------! GFS-setup for cdmb & cleff
! cdmb =  4.0 * (192.0/IMX)  
! cleff = 0.5E-5 / SQRT(IMX/192.0)  = 0.5E-5*SQRT(192./IMX)
!
      real(kind=kind_phys), parameter :: lonr_refmb =  4.0 * 192.0
      real(kind=kind_phys), parameter :: lonr_refgw =  192.0
      real(kind=kind_phys), parameter :: cleff_ref  = 0.5e-5           !       1256 km = 10 * 125 km ???
      
! copy  to "ugwp_oro_init"  =>  nwaves, nazdir, nstoch
 
      nworo  =  nwaves
      nazoro =  nazdir
      nstoro =  nstoch

      cdmbX = lonr_refmb/float(lonr)
      
      cdmb  = cdmbX
      cleff = cleff_ref * sqrt(lonr_refgw/float(lonr))   !* effac
!
      end subroutine init_oro_gws
!

    end module ugwp_oro_init
! =========================================================================
!
!    ugwp_conv_init
!
!=========================================================================
    module ugwp_conv_init
    
     use machine ,              only : kind_phys

     
     implicit none
      real(kind=kind_phys)    ::  eff_con                   ! scale factors for conv GWs
      integer ::  nwcon                     ! number of waves
      integer ::  nazcon                    ! number of azimuths
      integer ::  nstcon                    ! flag for stochastic choice of launch level above Conv-cloud
      real(kind=kind_phys)    ::  con_dlength
      real(kind=kind_phys)    ::  con_cldf

      real(kind=kind_phys), parameter    :: cmin  =  5  !2.5
      real(kind=kind_phys), parameter    :: cmax  = 95. !82.5
      real(kind=kind_phys), parameter    :: cmid  = 22.5
      real(kind=kind_phys), parameter    :: cwid  = cmid
      real(kind=kind_phys), parameter    :: bns   = 2.e-2, bns2 = bns*bns, bns4=bns2*bns2
      real(kind=kind_phys), parameter    :: mstar = 6.28e-3/2.  !  2km
      real(kind=kind_phys)               :: dc

      real(kind=kind_phys), allocatable  :: ch_conv(:),  spf_conv(:)
      real(kind=kind_phys), allocatable  :: xaz_conv(:), yaz_conv(:)
     contains
!
     subroutine init_conv_gws(nwaves, nazdir, nstoch, effac, lonr, kxw)
!
     use ugwp_common,  only : pi2, arad

     implicit none
     
 
      integer :: nwaves, nazdir, nstoch
      integer :: lonr
!      
! ccpp
!      
      
      real(kind=kind_phys)    :: kxw,  effac
      real(kind=kind_phys)    :: work1 = 0.5
      real(kind=kind_phys)    :: chk, tn4, snorm
      integer :: k

      nwcon    = nwaves
      nazcon   = nazdir
      nstcon   = nstoch
      eff_con  = effac

      con_dlength = pi2*arad/float(lonr)      
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
 
      call init_nazdir(nazdir,  xaz_conv,  yaz_conv)
     end subroutine init_conv_gws


    end module ugwp_conv_init
!=========================================================================
!
!    ugwp_fjet_init
!
!=========================================================================

   module ugwp_fjet_init
   use machine ,              only : kind_phys   

   
   
      implicit none
      real(kind=kind_phys)    ::  eff_fj                     ! scale factors for conv GWs
      integer ::  nwfj                       ! number of waves
      integer ::  nazfj                      ! number of azimuths
      integer ::  nstfj                      ! flag for stochastic choice of launch level above Conv-cloud
!
      real(kind=kind_phys), parameter    ::  fjet_trig=0.    ! if ( abs(frgf) > fjet_trig ) launch GW-packet

      
      real(kind=kind_phys), parameter    :: cmin =  2.5
      real(kind=kind_phys), parameter    :: cmax = 67.5
      real(kind=kind_phys)               :: dc
      real(kind=kind_phys), allocatable  :: ch_fjet(:) , spf_fjet(:)
      real(kind=kind_phys), allocatable  :: xaz_fjet(:), yaz_fjet(:)
     contains
     
     subroutine init_fjet_gws(nwaves, nazdir, nstoch, effac,lonr, kxw) 

     use ugwp_common,  only : pi2, arad			      
			      
     implicit none

      integer :: nwaves, nazdir, nstoch
      integer :: lonr
      real(kind=kind_phys)    :: kxw,  effac , chk
 
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
      call init_nazdir(nazdir,  xaz_fjet,  yaz_fjet)

     end subroutine init_fjet_gws

    end module ugwp_fjet_init
!
!=========================================================================
!
!
     module ugwp_okw_init
!=========================================================================
       use machine ,              only : kind_phys      
     
       implicit none

      real(kind=kind_phys)    ::  eff_okw                     ! scale factors for conv GWs
      integer ::  nwokw                       ! number of waves
      integer ::  nazokw                      ! number of azimuths
      integer ::  nstokw                      ! flag for stochastic choice of launch level above Conv-cloud
!
      real(kind=kind_phys), parameter    ::  okw_trig=0.      ! if ( abs(okwp) > okw_trig ) launch GW-packet

      real(kind=kind_phys), parameter    :: cmin =  2.5
      real(kind=kind_phys), parameter    :: cmax = 67.5
      real(kind=kind_phys)               :: dc
      real(kind=kind_phys), allocatable  :: ch_okwp(:),   spf_okwp(:)
      real(kind=kind_phys), allocatable  :: xaz_okwp(:),  yaz_okwp(:)

     contains
!

			     
     subroutine init_okw_gws(nwaves, nazdir, nstoch, effac, lonr, kxw)
     use ugwp_common,  only : pi2, arad

     implicit none

      integer :: nwaves, nazdir, nstoch
      integer :: lonr
      real(kind=kind_phys)    :: kxw,  effac , chk
 
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
  
      call init_nazdir(nazdir,  xaz_okwp,  yaz_okwp)
!
     end subroutine init_okw_gws
 
     end module ugwp_okw_init

!=============================== end of GW  sources
!
!  init specific  gw-solvers (1,2,3,4)
!
!===============================
!  Part -3  init  wave solvers
!===============================

  module ugwp_lsatdis_init
     use machine ,              only : kind_phys    
     implicit none

      integer  :: nwav, nazd
      integer  :: nst
      real(kind=kind_phys)     :: eff
      integer, parameter  :: incdim = 4, iazdim = 4
!
     contains

     subroutine initsolv_lsatdis(me, master,  nwaves, nazdir, nstoch, effac, do_physb, kxw)

     implicit none
!
     integer  :: me, master
     integer  :: nwaves, nazdir
     integer  :: nstoch
     real(kind=kind_phys)     :: effac
     logical  :: do_physb
     real(kind=kind_phys)     :: kxw
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
  end module ugwp_lsatdis_init
!
!
  module ugwp_wmsdis_init
  
     use machine ,    only : kind_phys    
     use ugwp_common, only : arad, pi, pi2, hpscale, rhp, rhp2, rh4, omega2 
     use ugwp_common, only : bnv2max,  bnv2min, minvel 
     use ugwp_common, only : mkzmin,  mkz2min, mkzmax, mkz2max, ucrit => cdmin
     
    implicit none

      real(kind=kind_phys),     parameter   :: maxdudt = 250.e-5, maxdtdt=15.e-2
      real(kind=kind_phys),     parameter   :: dked_min =0.01,    dked_max=250.0
   
      real(kind=kind_phys),     parameter   :: gptwo=2.0
 
      real(kind=kind_phys) ,     parameter  :: bnfix  =  6.28/300., bnfix2= bnfix * bnfix
      real(kind=kind_phys) ,     parameter  :: bnfix4 =  bnfix2 * bnfix2  
      real(kind=kind_phys) ,     parameter  :: bnfix3 =  bnfix2 * bnfix        
!
! make parameter list that will be passed to SOLVER
! 
      integer  , parameter  :: iazidim=4       ! number of azimuths
      integer  , parameter  :: incdim=25       ! number of discrete cx - spectral elements in launch spectrum
 
      real(kind=kind_phys) ,     parameter  :: zcimin = 2.5
      real(kind=kind_phys) ,     parameter  :: zcimax = 125.0
      real(kind=kind_phys) ,     parameter  :: zgam   = 0.25
!
! Verical spectra
!
      real(kind=kind_phys) ,     parameter  :: pind_wd = 5./3.
      real(kind=kind_phys) ,     parameter  :: sind_kz = 1.
      real(kind=kind_phys) ,     parameter  :: tind_kz = 3.    
      real(kind=kind_phys) ,     parameter  :: stind_kz = sind_kz + tind_kz  
!
! copies from kmob_ugwp namelist
!      
      real(kind=kind_phys)                  :: nslope              ! the GW sprctral slope at small-m             
      real(kind=kind_phys)                  :: lzstar
      real(kind=kind_phys)                  :: lzmin
      real(kind=kind_phys)                  :: lzmax      
      real(kind=kind_phys)                  :: lhmet
      real(kind=kind_phys)                  :: tamp_mpa            !amplitude for GEOS-5/MERRA-2
      real(kind=kind_phys)                  :: tau_min             ! min of GW MF 0.25 mPa                        
      integer               :: ilaunch
      real(kind=kind_phys)                  :: gw_eff

      real(kind=kind_phys)                  :: v_kxw, rv_kxw, v_kxw2 
     
     
 
!===========================================================================
      integer  :: nwav, nazd, nst
      real(kind=kind_phys)     :: eff
 
      real(kind=kind_phys)                :: zaz_fct, zms
      real(kind=kind_phys), allocatable   :: zci(:), zci4(:), zci3(:),zci2(:), zdci(:)
      real(kind=kind_phys), allocatable   :: zcosang(:), zsinang(:)
      real(kind=kind_phys), allocatable   :: lzmet(:), czmet(:), mkzmet(:), dczmet(:), dmkz(:) 

!
! GW-eddy constants for wave-mode dissipation by background and stability of
!         "final" flow after application of GW-effects
! 
      real(kind=kind_phys), parameter :: iPr_pt = 0.5
      real(kind=kind_phys), parameter :: lturb = 30., sc2  = lturb*lturb     ! stable on 80-km TL lmix ~ 500 met.
      real(kind=kind_phys), parameter :: ulturb=150., sc2u = ulturb* ulturb  ! unstable
      real(kind=kind_phys), parameter :: ric =0.25
      real(kind=kind_phys), parameter :: rimin = -10., prmin = 0.25
      real(kind=kind_phys), parameter :: prmax = 4.0
!
      contains
!============================================================================
     subroutine initsolv_wmsdis(me, master,  nwaves, nazdir, nstoch, effac, do_physb, kxw, version)
 
!        call initsolv_wmsdis(me, master, knob_ugwp_wvspec(2), knob_ugwp_azdir(2), &
!         knob_ugwp_stoch(2), knob_ugwp_effac(2), do_physb_gwsrcs, kxw,version)
!
     implicit none
!
!input-control for solvers:
!      nwaves, nazdir, nstoch, effac, do_physb, kxw
!
!
     integer, intent(in)  :: me, master, nwaves, nazdir, nstoch
     integer, intent(in)  :: version 
        
     real(kind=kind_phys), intent(in)     :: effac, kxw
     logical,  intent(in) :: do_physb
      
!
!locals
!
      real(kind=kind_phys)     :: dlzmet   
      real(kind=kind_phys)     :: cstar,rcstar, nslope3, fnorm, zcin     
            
      integer  :: inc, jk, jl, iazi
!
      real(kind=kind_phys) :: zang, zang1, znorm
      real(kind=kind_phys) :: zx1, zx2, ztx, zdx, zxran, zxmin, zxmax, zx, zpexp
      real(kind=kind_phys) :: fpc, fpc_dc
      real(kind=kind_phys) :: ae1,ae2
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

      
      v_kxw  = kxw ; v_kxw2 = v_kxw*v_kxw 
      rv_kxw = 1./v_kxw

      allocate ( zci(nwav),  zci4(nwav), zci3(nwav),zci2(nwav), zdci(nwav)  )
      allocate ( zcosang(nazd), zsinang(nazd) )
      allocate (lzmet(nwav), czmet(nwav), mkzmet(nwav), dczmet(nwav), dmkz(nwav) )   

!      if (me == master) then
!         print *, 'ugwp_v1/v0: init_gw_wmsdis_control '
!  
!         print *, 'ugwp_v1/v0: WMS_DIS launch layer ',    ilaunch
!         print *, 'ugwp_v1/v0: WMS_DIS tot_mflux in mpa', tamp_mpa*1000.
!         print *, 'ugwp_v1/v0: WMS_DIS lhmet in km  ' ,   lhmet*1.e-3      	 
!       endif

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
! Scinocca 2003.   x = 1/c stretching transform
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
!  alternatuve lzmax-lzmin without  x=1/c transform
!
!
        if (version == 1) then 
	
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
          do inc=1, nwav	    	  
          zdci(inc) = zdx     
	  enddo   
	 
	  cstar = bnfix/zms
	  rcstar = 1./cstar	
	ENDIF                               !   if (version == 1) then 	 
	 
	RETURN	  	 
!===================  Diag prints after return  ====================
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
	   nslope3=nslope+3.0
         do inc=1, nwav	    	  
           zcin =zci(inc)*rcstar
	   fpc =  rcstar*(zcin*zcin)/(1.+ zcin**nslope3)	   	    
           fpc_dc = fpc * zdci(inc)
	   write(6,111)  inc, zci(inc), zdci(inc),ucrit, fpc, fpc_dc, 6.28e-3/bnfix*zci(inc)
          enddo
         endif
	 

	 
 111     format( 'wms-zci', i4, 7 (3x, F8.3))

     end subroutine initsolv_wmsdis
!
!
  end module ugwp_wmsdis_init

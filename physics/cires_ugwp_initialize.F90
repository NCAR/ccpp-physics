!===============================
! cu-cires ugwp-scheme
!    initialization of selected 
!  init gw-solvers (1,2,3,4)
!  init gw-source specifications
!  init gw-background dissipation
!==============================
!
! Part-0   specifications of common constants, limiters and "criiical" values
 

!   module oro_state
    
!   integer, parameter :: kind_phys=8 
!   integer, parameter :: nvaroro=14   
!   real (kind=kind_phys), allocatable :: oro_stat(:, :)
!   contains

!   subroutine fill_oro_stat(nx, oc, oa4, clx4, theta, gamm, sigma, elvmax, hprime)
     
!    real  (kind=kind_phys),dimension(nx) :: oc, theta, gamm, sigma, elvmax, hprime  
!    real(kind=kind_phys),dimension(nx,4) :: oa4, clx4
!    integer :: i
!    do i=1, nx
!      oro_stat(i,1)    = hprime(i)
!      oro_stat(i,2)    =  oc(i)
!      oro_stat(i,3:6)  = oa4(i,1:4)
!      oro_stat(i,7:10) = clx4(i,1:4)
!      oro_stat(i,11)   = theta(i)
!      oro_stat(i,12)   = gamm(i)
!      oro_stat(i,13)   = sigma(i)
!      oro_stat(i,14)   = elvmax(i)
!    enddo   
!   end   subroutine fill_oro_stat

!   end module oro_state
    
    module ugwp_common
!
     implicit none

      real, parameter ::  grav =9.80665, cpd = 1004.6, grcp = grav/cpd
      real, parameter ::  rd = 287.05 , rv =461.5
      real, parameter ::  rgrav = 1.0/grav

      real, parameter ::  fv   = rv/rd - 1.0
      real, parameter ::  rdi  = 1.0 / rd
      real, parameter ::  gor  = grav/rd
      real, parameter ::  gr2  = grav*gor
      real, parameter ::  gocp = grav/cpd
      real, parameter ::  pi = 4.*atan(1.0), pi2 = 2.*pi
!
      real, parameter ::  rad_to_deg=180.0/pi, deg_to_rad=pi/180.0

      real, parameter ::  arad = 6370.e3
      real, parameter ::  rcpd2 = 0.5/cpd,  rcpd = 1./cpd
      real, parameter ::  dw2min=1.0
      real, parameter ::  bnv2min=1.e-6
      real, parameter ::  velmin=sqrt(dw2min)
      real, parameter ::  omega1 = pi2/86400.
      real, parameter ::  omega2 = 2.*omega1
     end module ugwp_common
!
!
!===================================================
!
!Part-1 init =>   wave dissipation + RFriction
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
     real, parameter :: alpha  = 1./86400./15.
!     
     real, parameter :: kdrag  = 1./86400./10.
     real, parameter :: zdrag  = 100.
     real, parameter :: zgrow  = 50.
!
     real ::    vumol, mumol, keddy, ion_drag
!
     do k=1, levs
      vumol    = vusurf*exp(-zkm(k)/hpmol)
      mumol    = musurf*exp(-zkm(k)/hpmol)

      keddy    = kturbo*exp(-((zkm(k)-zturbo) /zturw)**2)

      kvg(k)   = vumol + keddy
      ktg(k)   = mumol + keddy*inv_pra

      krad(k)  = alpha
!
      ion_drag = kdrag
!
      kion(k)  = ion_drag
     enddo

      k= levs+1
      kion(k) = kion(k-1)
      krad(k) = krad(k-1)
      kvg(k)  =  kvg(k-1)
      ktg(k)  =  ktg(k-1)
!                                                          
     end subroutine init_global_gwdis
!
!
     subroutine  rf_damp_init(levs, pa_rf, tau_rf, dtp, pmb, rfdis, rfdist, levs_rf)
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

     end subroutine  rf_damp_init
! ========================================================================
! Part 2 - sources
!      wave  sources
! ========================================================================
!
!    ugwp_oro_init
!
!=========================================================================
     module ugwp_oro_init

     use ugwp_common, only : bnv2min, grav, grcp, fv, grav, cpd, grcp, pi

     implicit none
!  
! constants and "crirtical" values to run oro-mtb_gw physics
!
! choice of oro-scheme: strver = 'vay_2018' , 'gfs_2018', 'kdn_2005', 'smc_2000'
!
      character(len=8)  ::  strver  = 'gfs_2018'
      character(len=8)  ::  strbase = 'gfs_2018'
      real, parameter   ::  rimin=-10., ric=0.25

!
      real, parameter :: efmin=0.5,    efmax=10.0
      real, parameter :: hpmax=2400.0, hpmin=25.0
      real, parameter :: sigma_std=1./100., gamm_std=1.0

      real, parameter :: frmax=10., frc =1.0, frmin =0.01
!

       real, parameter :: ce=0.8, ceofrc=ce/frc, cg=0.5
       real, parameter :: gmax=1.0, veleps=1.0, factop=0.5
!
       real, parameter :: rlolev=50000.0
!
       real,      parameter :: hncrit=9000.   ! max value in meters for elvmax
 
!  hncrit set to 8000m and sigfac added to enhance elvmax mtn hgt

       real,      parameter :: sigfac=4.0     ! mb3a expt test for elvmax factor
       real,      parameter :: hminmt=50.     ! min mtn height (*j*)
       real,      parameter :: minwnd=1.0     ! min wind component (*j*)
       real,      parameter :: dpmin=5000.0   ! minimum thickness of the reference layer in pa
 
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

     real,    parameter   ::   odmin  = 0.1, odmax = 10.0
!------------------------------------------------------------------------------
!    small-scale orography parameters  for TOFD of Beljaars et al., 2004, QJRMS
!------------------------------------------------------------------------------

     integer, parameter  ::   n_tofd=2                  ! depth of SSO for TOFD compared with Zpbl
     real, parameter     ::  const_tofd =  0.0759       ! alpha*beta*Cmd*Ccorr*2.109 = 12.*1.*0.005*0.6*2.109 = 0.0759
     real, parameter     ::     ze_tofd =1500.0         ! BJ's z-decay in meters
     real, parameter     ::  a12_tofd =0.0002662*0.005363   ! BJ's k-spect const for sigf2 * a1*a2*exp(-[z/zdec]**1.5]
     real, parameter     ::      ztop_tofd =10.*ze_tofd     ! no TOFD > this height too higher 15 km
!------------------------------------------------------------------------------
!
      real, parameter :: fcrit_sm  = 0.7, fcrit_sm2  = fcrit_sm * fcrit_sm
      real, parameter :: fcrit_gfs = 0.7
      real, parameter :: fcrit_mtb = 0.7

      real,  parameter :: lzmax  = 18.e3                      ! 18 km
      real,  parameter :: mkzmin = 6.28/lzmax
      real,  parameter :: mkz2min = mkzmin*mkzmin
      real,  parameter :: zbr_pi  = 3./2.*4.*atan(1.0)        ! 3pi/2
      real,  parameter :: zbr_ifs = 2.*atan(1.0)              ! pi/2

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

! copy  to "ugwp_oro_init"  =>  nwaves, nazdir, nstoch
 
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

    end module ugwp_oro_init
! =========================================================================
!
!    ugwp_conv_init
!
!=========================================================================
    module ugwp_conv_init

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
                              lonr, kxw, cgwf)
     use ugwp_common,  only : pi2, arad
     implicit none
 
      integer :: nwaves, nazdir, nstoch
      integer :: lonr
      real    :: cgwf(2)
      real    :: kxw,  effac
      real    :: work1 = 0.5
      real    :: chk, tn4, snorm
      integer :: k

      nwcon    = nwaves
      nazcon   = nazdir
      nstcon   = nstoch
      eff_con  = effac

      con_dlength = pi2*arad/float(lonr)
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
 
      call init_nazdir(nazdir,  xaz_conv,  yaz_conv)
     end subroutine init_conv_gws


    end module ugwp_conv_init
!=========================================================================
!
!    ugwp_fjet_init
!
!=========================================================================

   module ugwp_fjet_init
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
     subroutine init_fjet_gws(nwaves, nazdir, nstoch, effac, lonr, kxw)
     use ugwp_common,  only : pi2, arad
     implicit none

      integer :: nwaves, nazdir, nstoch
      integer :: lonr
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
      call init_nazdir(nazdir,  xaz_fjet,  yaz_fjet)

     end subroutine init_fjet_gws

    end module ugwp_fjet_init
!
!=========================================================================
!
!
     module ugwp_okw_init
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
     subroutine init_okw_gws(nwaves, nazdir, nstoch, effac, lonr, kxw)

     use ugwp_common,  only : pi2, arad
     implicit none

      integer :: nwaves, nazdir, nstoch
      integer :: lonr
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

      call init_nazdir(nazdir,  xaz_okwp,  yaz_okwp)

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
  end module ugwp_lsatdis_init
!
!
  module ugwp_wmsdis_init
 
    implicit none

      real,     parameter   :: maxdudt = 250.e-5

      real,     parameter   :: hpscale= 7000., rhp2   =  0.5/hpscale
      real,     parameter   :: omega2 = 2.*6.28/86400
      real,     parameter   :: gptwo=2.0

      real,     parameter   :: dked_min =0.01
      real,     parameter   :: gssec = (6.28/30.)**2        ! max-value for bn2
      real,     parameter   :: bv2min = (6.28/60./120.)**2  ! min-value for bn2  7.6(-7)  2 hrs
      real,     parameter   :: minvel = 0.5
 
!
! make parameter list that will be passed to SOLVER
!

      real, parameter       :: v_kxw  = 6.28e-3/200.
      real, parameter       :: v_kxw2 = v_kxw*v_kxw
      real, parameter       :: tamp_mpa = 30.e-3
      real, parameter       :: zfluxglob= 3.75e-3

      real ,     parameter  :: nslope=1        ! the GW sprctral slope at small-m
!     integer, parameter    :: klaunch=55      ! 32 - ~ 1km ;55 - 5.5 km ; 52 4.7km ; 60-7km index for selecting launch level
!     integer, parameter    :: ilaunch=klaunch
 
      integer  , parameter  :: iazidim=4       ! number of azimuths
      integer  , parameter  :: incdim=25       ! number of discrete cx - spectral elements in launch spectrum
      real ,     parameter  :: ucrit2=0.5
 
      real ,     parameter  :: zcimin = ucrit2
      real ,     parameter  :: zcimax = 125.0
      real ,     parameter  :: zgam   =   0.25
      real ,     parameter  :: zms_l  = 2000.0

      integer               :: ilaunch
      real                  :: gw_eff
 
!===========================================================================
      integer  :: nwav, nazd, nst
      real     :: eff
 
      real                :: zaz_fct , zms
      real, allocatable   :: zci(:), zci4(:), zci3(:),zci2(:), zdci(:)
      real, allocatable   :: zcosang(:), zsinang(:)
      contains
!============================================================================
     subroutine initsolv_wmsdis(me, master,  nwaves, nazdir, nstoch, effac, do_physb, kxw)
 
!        call initsolv_wmsdis(me, master, knob_ugwp_wvspec(2), knob_ugwp_azdir(2), &
!         knob_ugwp_stoch(2), knob_ugwp_effac(2), do_physb_gwsrcs, kxw)
!
     use ugwp_common, only :   pi, pi2
     implicit none
!
!input -control for solvers:
!      nwaves, nazdir, nstoch, effac, do_physb, kxw
!
!
     integer  :: me, master, nwaves, nazdir, nstoch
     real     :: effac, kxw
     logical  :: do_physb
!
!locals
!
      integer :: inc, jk, jl, iazi
!
      real :: zang, zang1, znorm
      real :: zx1, zx2, ztx, zdx, zxran, zxmin, zxmax, zx, zpexp

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

      allocate ( zci(nwav),  zci4(nwav), zci3(nwav),zci2(nwav), zdci(nwav)  )
      allocate ( zcosang(nazd), zsinang(nazd) )

      if (me == master) then
         print *, 'ugwp_v0: init_gw_wmsdis_control '
!        print *, 'ugwp_v0: WMSDIS launch layer ',  klaunch
         print *, 'ugwp_v0: WMSDIS launch layer ',  ilaunch
         print *, 'ugwp_v0: WMSDID tot_mflux in mpa', tamp_mpa*1000.
       endif

       zpexp = gptwo * 0.5                    ! gptwo=2 , zpexp = 1.

!
!      set up azimuth directions and some trig factors
!
!
       zang=pi2/float(nazd)

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
       zaz_fct = 1.0
       zaz_fct = 2.0 / znorm            ! correction factot for azimuthal sums

!       define coordinate transform for "Ch"   ....x = 1/c stretching transform
!       ----------------------------------------------- 
! note that this is expresed in terms of the intrinsic phase speed
! at launch ci=c-u_o so that the transformation is identical
! see eq. 28-30 of scinocca 2003.   x = 1/c stretching transform
!
          zxmax = 1.0 / zcimin
          zxmin = 1.0 / zcimax
          zxran = zxmax - zxmin
          zdx   = zxran / real(nwav-1)                            ! dkz
!
          zx1   = zxran/(exp(zxran/zgam)-1.0 )                    ! zgam =1./4.
          zx2   = zxmin - zx1

!
! computations for zci =1/zx
!                                   if(lgacalc) zgam=(zxmax-zxmin)/log(zxmax/zxmin)
!                                   zx1=zxran/(exp(zxran/zgam)-1.0_jprb)
!                                   zx2=zxmin-zx1
        zms = 2.*pi/zms_l
        do inc=1, nwav
          ztx = real(inc-1)*zdx+zxmin
          zx  = zx1*exp((ztx-zxmin)/zgam)+zx2                            !eq. 29 of scinocca 2003
          zci(inc)  = 1.0 /zx                                            !eq. 28 of scinocca 2003
          zdci(inc) = zci(inc)**2*(zx1/zgam)*exp((ztx-zxmin)/zgam)*zdx   !eq. 30 of scinocca 2003
          zci4(inc) = (zms*zci(inc))**4
          zci2(inc) = (zms*zci(inc))**2
          zci3(inc) = (zms*zci(inc))**3
        enddo
!
!
!  all done and print-out
!
!
         if (me == master) then
           print *
           print *,  'ugwp_v0: zcimin=' , zcimin
           print *,  'ugwp_v0: zcimax=' , zcimax
           print *,  'ugwp_v0: cd_crit=', zgam                        ! m/s precision for crit-level
           print *,  'ugwp_v0: launch_level',  ilaunch
           print *, ' ugwp_v0 zms_l=', zms_l
           print *, ' ugwp_vgw  nslope=', nslope

           print *
         endif
 

     end subroutine initsolv_wmsdis
!
! make a list of  all-initilized parameters needed for "gw_solver_wmsdis"
!

  end module ugwp_wmsdis_init
!=========================================================================
!
! work TODO for 2-extra WAM-solvers:
!           DSPDIS (Hines)+ADODIS (Alexander-Dunkerton-Ortland)
!
!========================================================================= 
     subroutine init_dspdis
     implicit none
     end subroutine init_dspdis

     subroutine init_adodis
     implicit none
     end subroutine init_adodis
     

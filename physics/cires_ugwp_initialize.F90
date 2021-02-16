!===============================
! cu-cires ugwp-scheme
!  initialization of ugwp_common_v0
!  init gw-solvers (1,2) .. no UFS-funds for (3,4) tests
!  init gw-source specifications
!  init gw-background dissipation
!===============================    
    module ugwp_common_v0
!
     use machine,  only: kind_phys
     use physcons, only : pi => con_pi, grav => con_g, rd => con_rd,   &
                          rv => con_rv, cpd => con_cp, fv => con_fvirt,&
                          arad => con_rerth
     implicit none

      real(kind=kind_phys), parameter ::  grcp = grav/cpd, rgrav = 1.0d0/grav, &
                          rdi  = 1.0d0/rd,                                     &
                          gor  = grav/rd,  gr2   = grav*gor, gocp = grav/cpd,  &
                          rcpd = 1./cpd,   rcpd2 = 0.5*rcpd,                   &
                          pi2  = pi + pi,  omega1 = pi2/86400.0,               &
                          omega2 = omega1+omega1,                              &
                          rad_to_deg=180.0/pi, deg_to_rad=pi/180.0,            &
                          dw2min=1.0, bnv2min=1.e-6, velmin=sqrt(dw2min)


     end module ugwp_common_v0
!
!
!===================================================
!
!Part-1 init =>   wave dissipation + RFriction
!
!===================================================
     subroutine init_global_gwdis_v0(levs, zkm, pmb, kvg, ktg, krad, kion)
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
     end subroutine init_global_gwdis_v0

     
! ========================================================================
! Part 2 - sources
!      wave  sources
! ========================================================================
!
!    ugwpv0_oro_init
!
!=========================================================================
     module ugwpv0_oro_init

     use ugwp_common_v0, only : bnv2min, grav, grcp, fv, grav, cpd, grcp, pi

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

       real, parameter :: ce=0.8,   ceofrc=ce/frc, cg=0.5
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

      real,  parameter :: lzmax   = 18.e3                      ! 18 km
      real,  parameter :: mkzmin  = 6.28/lzmax
      real,  parameter :: mkz2min = mkzmin*mkzmin
      real,  parameter :: zbr_pi  = (3.0/2.0)*pi
      real,  parameter :: zbr_ifs = 0.5*pi

      contains
!
      subroutine init_oro_gws_v0(nwaves, nazdir, nstoch, effac, &
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
      end subroutine init_oro_gws_v0
!

    end module ugwpv0_oro_init
!=============================== end of GW  sources
!
!  init specific  gw-solvers (1,2,3,4)
!

!===============================
!  Part -3  init  wave solvers
!===============================

  module ugwpv0_lsatdis_init
     implicit none

      integer  :: nwav, nazd
      integer  :: nst
      real     :: eff
      integer, parameter  :: incdim = 4, iazdim = 4
!
     contains

     subroutine initsolv_lsatdis_v0(me, master,  nwaves, nazdir, nstoch, effac, do_physb, kxw)

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
     end subroutine initsolv_lsatdis_v0
!
  end module ugwpv0_lsatdis_init
!
!
  module ugwpv0_wmsdis_init
 
    use ugwp_common_v0, only :   pi, pi2
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
 
      integer  , parameter  :: iazidim=4       ! number of azimuths
      integer  , parameter  :: incdim=25       ! number of discrete cx - spectral elements in launch spectrum
      real ,     parameter  :: ucrit2=0.5
 
      real ,     parameter  :: zcimin = ucrit2
      real ,     parameter  :: zcimax = 125.0
      real ,     parameter  :: zgam   =   0.25
      real ,     parameter  :: zms_l  = 2000.0, zms = pi2 / zms_l, zmsi = 1.0 / zms

      integer               :: ilaunch
      real                  :: gw_eff
 
!===========================================================================
      integer  :: nwav, nazd, nst
      real     :: eff
 
      real                :: zaz_fct
      real, allocatable   :: zci(:), zci4(:), zci3(:),zci2(:), zdci(:)
      real, allocatable   :: zcosang(:), zsinang(:)
      contains
!============================================================================
     subroutine initsolv_wmsdis_v0(me, master,  nwaves, nazdir, nstoch, effac, do_physb, kxw)
 
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
!       zms = pi2 / zms_l
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

     end subroutine initsolv_wmsdis_v0
!
  end module ugwpv0_wmsdis_init

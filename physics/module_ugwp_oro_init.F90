! ========================================================================= 
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
           
     real,    parameter   ::   odmin  =0.1, odmax =10.0      
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

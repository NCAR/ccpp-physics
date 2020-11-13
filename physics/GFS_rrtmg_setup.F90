!> \file GFS_rrtmg_setup.f90
!! This file contains
module GFS_rrtmg_setup

   use physparam, only : isolar , ictmflg, ico2flg, ioznflg, iaerflg,&
!  &             iaermdl, laswflg, lalwflg, lavoflg, icldflg,         &
   &             iaermdl,                            icldflg,         &
   &             iovrRad=>iovr, lcrick , lcnorm , lnoprec,            &
   &             ialbflg, iemsflg, isubcsw, isubclw, ivflip , ipsd0,  &
   &             iswcliq,                                             &
   &             kind_phys

   use radcons, only: ltp, lextop

   implicit none

   public GFS_rrtmg_setup_init, GFS_rrtmg_setup_run, GFS_rrtmg_setup_finalize

   private

   logical :: is_initialized = .false.

   !  ---  version tag and last revision date
   character(40), parameter ::                                       &
        &   VTAGRAD='NCEP-Radiation_driver    v5.2  Jan 2013 '
   !    &   VTAGRAD='NCEP-Radiation_driver    v5.1  Nov 2012 '
   !    &   VTAGRAD='NCEP-Radiation_driver    v5.0  Aug 2012 '

   !> new data input control variables (set/reset in subroutines radinit/radupdate):
   integer :: month0 = 0
   integer :: iyear0 = 0
   integer :: monthd = 0

   !> control flag for the first time of reading climatological ozone data
   !! (set/reset in subroutines radinit/radupdate, it is used only if the
   !! control parameter ioznflg=0)
   logical :: loz1st = .true.

   contains

!> \defgroup GFS_rrtmg_setup GFS RRTMG Scheme Setup
!! @{
!! \section arg_table_GFS_rrtmg_setup_init Argument Table
!! \htmlinclude GFS_rrtmg_setup_init.html
!!
   subroutine GFS_rrtmg_setup_init (                                    &
          si, levr, ictm, isol, ico2, iaer, ialb, iems, ntcw,  num_p2d, &
          num_p3d, npdf3d, ntoz, iovr, isubc_sw, isubc_lw,              &
          icliq_sw, crick_proof, ccnorm,                                &
          imp_physics,                                                  &
          norad_precip, idate, iflip,                                   &
          im, faerlw, faersw, aerodp,                                   & ! for consistency checks
          me, errmsg, errflg)
! =================   subprogram documentation block   ================ !
!                                                                       !
! subprogram:   GFS_rrtmg_setup_init - a subprogram to initialize radiation !
!                                                                       !
! usage:        call GFS_rrtmg_setup_init                               !
!                                                                       !
! attributes:                                                           !
!   language:  fortran 90                                               !
!                                                                       !
! program history:                                                      !
!   mar 2012  - yu-tai hou   create the program to initialize fixed     !
!                 control variables for radiaion processes.  this       !
!                 subroutine is called at the start of model run.       !
!   nov 2012  - yu-tai hou   modified control parameter through         !
!                 module 'physparam'.                                   !
!   mar 2014  - sarah lu  iaermdl is determined from iaer               !
!   jul 2014  - s moorthi add npdf3d for pdf clouds                     !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
! input parameters:                                                     !
!   si               : model vertical sigma interface or equivalence    !
!   levr             : number of model vertical layers                  !
!   ictm             :=yyyy#, external data time/date control flag      !
!                     =   -2: same as 0, but superimpose seasonal cycle !
!                             from climatology data set.                !
!                     =   -1: use user provided external data for the   !
!                             forecast time, no extrapolation.          !
!                     =    0: use data at initial cond time, if not     !
!                             available, use latest, no extrapolation.  !
!                     =    1: use data at the forecast time, if not     !
!                             available, use latest and extrapolation.  !
!                     =yyyy0: use yyyy data for the forecast time,      !
!                             no further data extrapolation.            !
!                     =yyyy1: use yyyy data for the fcst. if needed, do !
!                             extrapolation to match the fcst time.     !
!   isol             := 0: use the old fixed solar constant in "physcon"!
!                     =10: use the new fixed solar constant in "physcon"!
!                     = 1: use noaa ann-mean tsi tbl abs-scale data tabl!
!                     = 2: use noaa ann-mean tsi tbl tim-scale data tabl!
!                     = 3: use cmip5 ann-mean tsi tbl tim-scale data tbl!
!                     = 4: use cmip5 mon-mean tsi tbl tim-scale data tbl!
!   ico2             :=0: use prescribed global mean co2 (old  oper)    !
!                     =1: use observed co2 annual mean value only       !
!                     =2: use obs co2 monthly data with 2-d variation   !
!   iaer             : 4-digit aerosol flag (dabc for aermdl,volc,lw,sw)!
!                     d: =0 or none, opac-climatology aerosol scheme    !
!                        =1 use gocart climatology aerosol scheme       !
!                        =2 use gocart progostic aerosol scheme         !
!                     a: =0 use background stratospheric aerosol        !
!                        =1 incl stratospheric vocanic aeros            !
!                     b: =0 no topospheric aerosol in lw radiation      !
!                        =1 include tropspheric aerosols for lw         !
!                     c: =0 no topospheric aerosol in sw radiation      !
!                        =1 include tropspheric aerosols for sw         !
!   ialb             : control flag for surface albedo schemes          !
!                     =0: climatology, based on surface veg types       !
!                     =1: modis retrieval based surface albedo scheme   !
!   iems             : ab 2-digit control flag                          !
!                     a: =0 set sfc air/ground t same for lw radiation  !
!                        =1 set sfc air/ground t diff for lw radiation  !
!                     b: =0 use fixed sfc emissivity=1.0 (black-body)   !
!                        =1 use varying climtology sfc emiss (veg based)!
!                        =2 future development (not yet)                !
!   ntcw             :=0 no cloud condensate calculated                 !
!                     >0 array index location for cloud condensate      !
!   num_p3d          :=3: ferrier's microphysics cloud scheme           !
!                     =4: zhao/carr/sundqvist microphysics cloud        !
!   npdf3d            =0 no pdf clouds                                  !
!                     =3 (when num_p3d=4) pdf clouds with zhao/carr/    !
!                        sundqvist scheme                               !
!   ntoz             : ozone data control flag                          !
!                     =0: use climatological ozone profile              !
!                     >0: use interactive ozone profile                 !
!   icliq_sw         : sw optical property for liquid clouds            !
!                     =0:input cld opt depth, ignoring iswcice setting  !
!                     =1:cloud optical property scheme based on Hu and  !
!                        Stamnes(1993) \cite hu_and_stamnes_1993 method !
!                     =2:cloud optical property scheme based on Hu and  !
!                        Stamnes(1993) -updated                         !
!   iovr             : control flag for cloud overlap (sw/lw rad)       !
!                     =0: random overlapping clouds                     !
!                     =1: max/ran overlapping clouds                    !
!                     =2: maximum overlap clouds       (mcica only)     !
!                     =3: decorrelation-length overlap (mcica only)     !
!                     =4: exponential overlap clouds
!   isubc_sw/isubc_lw: sub-column cloud approx control flag (sw/lw rad) !
!                     =0: with out sub-column cloud approximation       !
!                     =1: mcica sub-col approx. prescribed random seed  !
!                     =2: mcica sub-col approx. provided random seed    !
!   crick_proof      : control flag for eliminating CRICK               !
!   ccnorm           : control flag for in-cloud condensate mixing ratio!
!   norad_precip     : control flag for not using precip in radiation   !
!   idate(4)         : ncep absolute date and time of initial condition !
!                      (hour, month, day, year)                         !
!   iflip            : control flag for direction of vertical index     !
!                     =0: index from toa to surface                     !
!                     =1: index from surface to toa                     !
!   me               : print control flag                               !
!                                                                       !
!  subroutines called: radinit                                          !
!                                                                       !
!  ===================================================================  !
!
      use module_radsw_parameters,  only: NBDSW
      use module_radlw_parameters,  only: NBDLW
      use module_radiation_aerosols,only: NF_AELW, NF_AESW, NSPC1
      use module_radiation_clouds,  only: NF_CLDS
      use module_radiation_gases,   only: NF_VGAS
      use module_radiation_surface, only: NF_ALBD

      implicit none

      ! interface variables
      real (kind=kind_phys), intent(in) :: si(levr+1)
      integer, intent(in) :: levr
      integer, intent(in) :: ictm
      integer, intent(in) :: isol
      integer, intent(in) :: ico2
      integer, intent(in) :: iaer
      integer, intent(in) :: ialb
      integer, intent(in) :: iems
      integer, intent(in) :: ntcw
      integer, intent(in) :: num_p2d
      integer, intent(in) :: num_p3d
      integer, intent(in) :: npdf3d
      integer, intent(in) :: ntoz
      integer, intent(in) :: iovr
      integer, intent(in) :: isubc_sw
      integer, intent(in) :: isubc_lw
      integer, intent(in) :: icliq_sw
      logical, intent(in) :: crick_proof
      logical, intent(in) :: ccnorm
      integer, intent(in) :: imp_physics
      logical, intent(in) :: norad_precip
      integer, intent(in) :: idate(4)
      integer, intent(in) :: iflip
      ! For consistency checks
      integer, intent(in)         :: im
      real(kind_phys), intent(in) :: faerlw(:,:,:,:)
      real(kind_phys), intent(in) :: faersw(:,:,:,:)
      real(kind_phys), intent(in) :: aerodp(:,:)
      ! End for consistency checks
      integer, intent(in)           :: me
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      ! For consistency checks
      real(kind_phys), dimension(im,levr+ltp,NBDLW,NF_AELW) :: faerlw_check
      real(kind_phys), dimension(im,levr+ltp,NBDSW,NF_AESW) :: faersw_check
      real(kind_phys), dimension(im,NSPC1)                  :: aerodp_check
      ! End for consistency checks

      ! Initialize the CCPP error handling variables
      errmsg = ''
      errflg = 0

      if (is_initialized) return

      ! Consistency checks for dimensions of arrays, this is required
      ! to detect differences in FV3's parameters that are used to
      ! dimension certain arrays and the values in ccpp-physics
      if (size(faerlw(1,:,:,:)).ne.size(faerlw_check(1,:,:,:))) then
         write(errmsg,"(3a,4i4,a,4i4)") &
               "Runtime error: dimension mismatch for faerlw,",        &
               " check definitions of levr, ltp, nbdlw, nf_aelw:",     &
               " expected shape ", shape(faerlw_check(:,:,:,:)),       &
               " but got ", shape(faerlw(:,:,:,:))
         errflg = 1
         return
      end if
      if (size(faersw(1,:,:,:)).ne.size(faersw_check(1,:,:,:))) then
         write(errmsg,"(3a,4i4,a,4i4)") &
               "Runtime error: dimension mismatch for faersw,",        &
               " check definitions of levr, ltp, nbdsw, nf_aesw:",     &
               " expected shape ", shape(faersw_check(:,:,:,:)),       &
               " but got ", shape(faersw(:,:,:,:))
         errflg = 1
         return
      end if
      if (size(aerodp(1,:)).ne.size(aerodp_check(1,:))) then
         write(errmsg,"(3a,2i4,a,2i4)") &
               "Runtime error: dimension mismatch for aerodp,",        &
               " check definitions of nspc1:",                         &
               " expected shape ", shape(aerodp_check(:,:)),           &
               " but got ", shape(aerodp(:,:))
         errflg = 1
         return
      end if
      
      ! End of consistency checks

      isolar = isol                     ! solar constant control flag

      ictmflg= ictm                     ! data ic time/date control flag
      ico2flg= ico2                     ! co2 data source control flag
      ioznflg= ntoz                     ! ozone data source control flag

      if ( ictm==0 .or. ictm==-2 ) then
        iaerflg = mod(iaer, 100)        ! no volcanic aerosols for clim hindcast
      else
        iaerflg = mod(iaer, 1000)   
      endif
      iaermdl = iaer/1000               ! control flag for aerosol scheme selection
      if ( iaermdl < 0 .or.  (iaermdl>2 .and. iaermdl/=5) ) then
         print *, ' Error -- IAER flag is incorrect, Abort'
         stop 7777
      endif

!     if ( ntcw > 0 ) then
        icldflg = 1                     ! prognostic cloud optical prop scheme
!     else
!       icldflg = 0                     ! no support for diag cloud opt prop scheme
!     endif

      iswcliq = icliq_sw                ! optical property for liquid clouds for sw

      ! iovr comes from the model. In the RRTMG implementation this is stored in phyrparam.f,
      ! it comes in from the host-model and is set here. 
      ! In GP, iovr is passed directly into the routines.
      iovrRAD = iovr
      lcrick  = crick_proof             ! control flag for eliminating CRICK 
      lcnorm  = ccnorm                  ! control flag for in-cld condensate 
      lnoprec = norad_precip            ! precip effect on radiation flag (ferrier microphysics)
      isubcsw = isubc_sw                ! sub-column cloud approx flag in sw radiation
      isubclw = isubc_lw                ! sub-column cloud approx flag in lw radiation

      ialbflg= ialb                     ! surface albedo control flag
      iemsflg= iems                     ! surface emissivity control flag

      ivflip = iflip                    ! vertical index direction control flag

!  ---  assign initial permutation seed for mcica cloud-radiation
      if ( isubc_sw>0 .or. isubc_lw>0 ) then
!       ipsd0 = 17*idate(1)+43*idate(2)+37*idate(3)+23*idate(4) + ipsd0
        ipsd0 = 17*idate(1)+43*idate(2)+37*idate(3)+23*idate(4)
      endif

      if ( me == 0 ) then
        print *,'  In rad_initialize (GFS_rrtmg_setup_init), before calling radinit'
        print *,' si =',si
        print *,' levr=',levr,' ictm=',ictm,' isol=',isol,' ico2=',ico2,&
     &          ' iaer=',iaer,' ialb=',ialb,' iems=',iems,' ntcw=',ntcw
        print *,' np3d=',num_p3d,' ntoz=',ntoz,                         &
     &          ' iovr=',iovr,' isubc_sw=',isubc_sw,                    &
     &          ' isubc_lw=',isubc_lw,' icliq_sw=',icliq_sw,            &
     &          ' iflip=',iflip,'  me=',me
        print *,' crick_proof=',crick_proof,                            &
     &          ' ccnorm=',ccnorm,' norad_precip=',norad_precip
      endif

      call radinit                                                      &
!  ---  inputs:
     &     ( si, levr, imp_physics, me )
!  ---  outputs:
!          ( none )

      if ( me == 0 ) then
        print *,'  Radiation sub-cloud initial seed =',ipsd0,           &
     &          ' IC-idate =',idate
        print *,' return from rad_initialize (GFS_rrtmg_setup_init) - after calling radinit'
      endif
!
      is_initialized = .true.
!
      return

   end subroutine GFS_rrtmg_setup_init

!> \section arg_table_GFS_rrtmg_setup_run Argument Table
!! \htmlinclude GFS_rrtmg_setup_run.html
!!
   subroutine GFS_rrtmg_setup_run (                &
          idate, jdate, deltsw, deltim, lsswr, me, &
          slag, sdec, cdec, solcon, errmsg, errflg)

      implicit none

      ! interface variables
      integer,              intent(in)  :: idate(:)
      integer,              intent(in)  :: jdate(:)
      real(kind=kind_phys), intent(in)  :: deltsw
      real(kind=kind_phys), intent(in)  :: deltim
      logical,              intent(in)  :: lsswr
      integer,              intent(in)  :: me
      real(kind=kind_phys), intent(out) :: slag
      real(kind=kind_phys), intent(out) :: sdec
      real(kind=kind_phys), intent(out) :: cdec
      real(kind=kind_phys), intent(out) :: solcon
      character(len=*),     intent(out) :: errmsg
      integer,              intent(out) :: errflg

      ! Check initialization state
      if (.not.is_initialized) then
         write(errmsg, fmt='((a))') 'GFS_rrtmg_setup_run called before GFS_rrtmg_setup_init'
         errflg = 1
         return
      end if

      ! Initialize the CCPP error handling variables
      errmsg = ''
      errflg = 0

      call radupdate(idate,jdate,deltsw,deltim,lsswr,me, &
                     slag,sdec,cdec,solcon)

   end subroutine GFS_rrtmg_setup_run

!> \section arg_table_GFS_rrtmg_setup_finalize Argument Table
!! \htmlinclude GFS_rrtmg_setup_finalize.html
!!
   subroutine GFS_rrtmg_setup_finalize (errmsg, errflg)

      implicit none

      character(len=*),          intent(  out) :: errmsg
      integer,                   intent(  out) :: errflg

      ! Initialize the CCPP error handling variables
      errmsg = ''
      errflg = 0

      if (.not.is_initialized) return

      ! do finalization stuff if needed

      is_initialized = .false.

   end subroutine GFS_rrtmg_setup_finalize


! Private functions


   subroutine radinit( si, NLAY, imp_physics, me )
!...................................

!  ---  inputs:
!     &     ( si, NLAY, imp_physics, me )
!  ---  outputs:
!          ( none )

! =================   subprogram documentation block   ================ !
!                                                                       !
! subprogram:   radinit     initialization of radiation calculations    !
!                                                                       !
! usage:        call radinit                                            !
!                                                                       !
! attributes:                                                           !
!   language:  fortran 90                                               !
!   machine:   wcoss                                                    !
!                                                                       !
!  ====================  definition of variables  ====================  !
!                                                                       !
! input parameters:                                                     !
!   si               : model vertical sigma interface                   !
!   NLAY             : number of model vertical layers                  !
!   imp_physics      : MP identifier                                    !
!   me               : print control flag                               !
!                                                                       !
!  outputs: (none)                                                      !
!                                                                       !
!  external module variables:  (in module physparam)                     !
!   isolar   : solar constant cntrol flag                               !
!              = 0: use the old fixed solar constant in "physcon"       !
!              =10: use the new fixed solar constant in "physcon"       !
!              = 1: use noaa ann-mean tsi tbl abs-scale with cycle apprx!
!              = 2: use noaa ann-mean tsi tbl tim-scale with cycle apprx!
!              = 3: use cmip5 ann-mean tsi tbl tim-scale with cycl apprx!
!              = 4: use cmip5 mon-mean tsi tbl tim-scale with cycl apprx!
!   iaerflg  : 3-digit aerosol flag (abc for volc, lw, sw)              !
!              a:=0 use background stratospheric aerosol                !
!                =1 include stratospheric vocanic aeros                 !
!              b:=0 no topospheric aerosol in lw radiation              !
!                =1 compute tropspheric aero in 1 broad band for lw     !
!                =2 compute tropspheric aero in multi bands for lw      !
!              c:=0 no topospheric aerosol in sw radiation              !
!                =1 include tropspheric aerosols for sw                 !
!   ico2flg  : co2 data source control flag                             !
!              =0: use prescribed global mean co2 (old  oper)           !
!              =1: use observed co2 annual mean value only              !
!              =2: use obs co2 monthly data with 2-d variation          !
!   ictmflg  : =yyyy#, external data ic time/date control flag          !
!              =   -2: same as 0, but superimpose seasonal cycle        !
!                      from climatology data set.                       !
!              =   -1: use user provided external data for the          !
!                      forecast time, no extrapolation.                 !
!              =    0: use data at initial cond time, if not            !
!                      available, use latest, no extrapolation.         !
!              =    1: use data at the forecast time, if not            !
!                      available, use latest and extrapolation.         !
!              =yyyy0: use yyyy data for the forecast time,             !
!                      no further data extrapolation.                   !
!              =yyyy1: use yyyy data for the fcst. if needed, do        !
!                      extrapolation to match the fcst time.            !
!   ioznflg  : ozone data source control flag                           !
!              =0: use climatological ozone profile                     !
!              =1: use interactive ozone profile                        !
!   ialbflg  : albedo scheme control flag                               !
!              =0: climatology, based on surface veg types              !
!              =1: modis retrieval based surface albedo scheme          !
!   iemsflg  : emissivity scheme cntrl flag (ab 2-digit integer)        !
!              a:=0 set sfc air/ground t same for lw radiation          !
!                =1 set sfc air/ground t diff for lw radiation          !
!              b:=0 use fixed sfc emissivity=1.0 (black-body)           !
!                =1 use varying climtology sfc emiss (veg based)        !
!                =2 future development (not yet)                        !
!   icldflg  : cloud optical property scheme control flag               !
!              =0: use diagnostic cloud scheme                          !
!              =1: use prognostic cloud scheme (default)                !
!   imp_physics  : cloud microphysics scheme control flag               !
!              =99 zhao/carr/sundqvist microphysics scheme              !
!              =98 zhao/carr/sundqvist microphysics+pdf cloud&cnvc,cnvw !
!              =11 GFDL cloud microphysics                              !
!              =8 Thompson microphysics scheme                          !
!              =6 WSM6 microphysics scheme                              !
!              =10 MG microphysics scheme                               !
!   iovr     : control flag for cloud overlap in radiation              !
!              =0: random overlapping clouds                            !
!              =1: max/ran overlapping clouds                           !
!   isubcsw  : sub-column cloud approx control flag in sw radiation     !
!   isubclw  : sub-column cloud approx control flag in lw radiation     !
!              =0: with out sub-column cloud approximation              !
!              =1: mcica sub-col approx. prescribed random seed         !
!              =2: mcica sub-col approx. provided random seed           !
!   lcrick   : control flag for eliminating CRICK                       !
!              =t: apply layer smoothing to eliminate CRICK             !
!              =f: do not apply layer smoothing                         !
!   lcnorm   : control flag for in-cld condensate                       !
!              =t: normalize cloud condensate                           !
!              =f: not normalize cloud condensate                       !
!   lnoprec  : precip effect in radiation flag (ferrier microphysics)   !
!              =t: snow/rain has no impact on radiation                 !
!              =f: snow/rain has impact on radiation                    !
!   ivflip   : vertical index direction control flag                    !
!              =0: index from toa to surface                            !
!              =1: index from surface to toa                            !
!                                                                       !
!  subroutines called: sol_init, aer_init, gas_init, cld_init,          !
!                      sfc_init, rlwinit, rswinit                       !
!                                                                       !
!  usage:       call radinit                                            !
!                                                                       !
!  ===================================================================  !
!

      use module_radiation_astronomy, only : sol_init
      use module_radiation_aerosols,  only : aer_init
      use module_radiation_gases,     only : gas_init
      use module_radiation_surface,   only : sfc_init
      use module_radiation_clouds,    only : cld_init
      ! DH* these should be called by rrtmg_lw_init and rrtmg_sw_init!
      use rrtmg_lw,                   only : rlwinit
      use rrtmg_sw,                   only : rswinit

      implicit none

!  ---  inputs:
      integer, intent(in) :: NLAY, me, imp_physics 

      real (kind=kind_phys), intent(in) :: si(:)

!  ---  outputs: (none, to module variables)

!  ---  locals:

!
!===> ...  begin here
!
!> -# Set up control variables and external module variables in
!!    module physparam
#if 0
      ! GFS_radiation_driver.F90 may in the future initialize air/ground
      ! temperature differently; however, this is not used at the moment
      ! and as such we avoid the difficulty of dealing with exchanging
      ! itsfc between GFS_rrtmg_setup and a yet-to-be-created/-used
      ! interstitial routine (or GFS_radiation_driver.F90)
      itsfc  = iemsflg / 10             ! sfc air/ground temp control
#endif
      loz1st = (ioznflg == 0)           ! first-time clim ozone data read flag
      month0 = 0
      iyear0 = 0
      monthd = 0

      if (me == 0) then
!       print *,' NEW RADIATION PROGRAM STRUCTURES -- SEP 01 2004'
        print *,' NEW RADIATION PROGRAM STRUCTURES BECAME OPER. ',      &
     &          '  May 01 2007'
        print *, VTAGRAD                !print out version tag
        print *,' - Selected Control Flag settings: ICTMflg=',ictmflg,  &
     &    ' ISOLar =',isolar, ' ICO2flg=',ico2flg,' IAERflg=',iaerflg,  &
     &    ' IALBflg=',ialbflg,' IEMSflg=',iemsflg,' ICLDflg=',icldflg,  &
     &    ' IMP_PHYSICS=',imp_physics,' IOZNflg=',ioznflg
        print *,' IVFLIP=',ivflip,' IOVR=',iovrRad,                     &
     &    ' ISUBCSW=',isubcsw,' ISUBCLW=',isubclw
        print *,' LCRICK=',lcrick,' LCNORM=',lcnorm,' LNOPREC=',lnoprec
        print *,' LTP =',ltp,', add extra top layer =',lextop

        if ( ictmflg==0 .or. ictmflg==-2 ) then
          print *,'   Data usage is limited by initial condition!'
          print *,'   No volcanic aerosols'
        endif

        if ( isubclw == 0 ) then
          print *,' - ISUBCLW=',isubclw,' No McICA, use grid ',         &
     &            'averaged cloud in LW radiation'
        elseif ( isubclw == 1 ) then
          print *,' - ISUBCLW=',isubclw,' Use McICA with fixed ',       &
     &            'permutation seeds for LW random number generator'
        elseif ( isubclw == 2 ) then
          print *,' - ISUBCLW=',isubclw,' Use McICA with random ',      &
     &            'permutation seeds for LW random number generator'
        else
          print *,' - ERROR!!! ISUBCLW=',isubclw,' is not a ',          &
     &            'valid option '
          stop
        endif

        if ( isubcsw == 0 ) then
          print *,' - ISUBCSW=',isubcsw,' No McICA, use grid ',         &
     &            'averaged cloud in SW radiation'
        elseif ( isubcsw == 1 ) then
          print *,' - ISUBCSW=',isubcsw,' Use McICA with fixed ',       &
     &            'permutation seeds for SW random number generator'
        elseif ( isubcsw == 2 ) then
          print *,' - ISUBCSW=',isubcsw,' Use McICA with random ',      &
     &            'permutation seeds for SW random number generator'
        else
          print *,' - ERROR!!! ISUBCSW=',isubcsw,' is not a ',          &
     &            'valid option '
          stop
        endif

        if ( isubcsw /= isubclw ) then
          print *,' - *** Notice *** ISUBCSW /= ISUBCLW !!!',           &
     &            isubcsw, isubclw
        endif
      endif

!> -# Initialization
!! - astronomy initialization routine:
!! call module_radiation_astronomy::sol_init()
!! - aerosols initialization routine:
!! call module_radiation_aerosols::aer_init()
!! - CO2 and other gases intialization routine:
!! call module_radiation_gases::gas_init()
!! - surface intialization routine:
!! call module_radiation_surface::sfc_init()
!! - cloud initialization routine:
!! call module_radiation_clouds::cld_init()
!! - LW radiation initialization routine:
!! call module_radlw_main::rlwinit()
!! - SW radiation initialization routine:
!! call module_radsw_main::rswinit()
!     Initialization

      call sol_init ( me )          !  --- ...  astronomy initialization routine

      call aer_init ( NLAY, me )    !  --- ...  aerosols initialization routine

      call gas_init ( me )          !  --- ...  co2 and other gases initialization routine

      call sfc_init ( me )          !  --- ...  surface initialization routine

      call cld_init ( si, NLAY, imp_physics, me) !  --- ...  cloud initialization routine

      call rlwinit ( me )           !  --- ...  lw radiation initialization routine

      call rswinit ( me )           !  --- ...  sw radiation initialization routine
!
      return
!...................................
      end subroutine radinit
      !-----------------------------------

!> This subroutine checks and updates time sensitive data used by
!! radiation computations. This subroutine needs to be placed inside
!! the time advancement loop but outside of the horizontal grid loop.
!! It is invoked at radiation calling frequncy but before any actual
!! radiative transfer computations.
!! \param idate          NCEP absolute date and time of intial condition
!!                       (year,month,day,time-zone,hour,minute,second,
!!                        mil-second)
!! \param jdate          NCEP absolute date and time at forecast time
!!                       (year,month,day,time-zone,hour,minute,second,
!!                        mil-second)
!! \param deltsw         SW radiation calling time interval in seconds
!! \param deltim         model advancing time-step duration in seconds
!! \param lsswr          logical control flag for SW radiation calculations
!! \param me             print control flag
!! \param slag           equation of time in radians
!! \param sdec,cdec      sine and cosine of the solar declination angle
!! \param solcon         solar constant adjusted by sun-earth distance \f$(W/m^2)\f$
!> \section gen_radupdate General Algorithm
!> @{
!-----------------------------------
      subroutine radupdate( idate,jdate,deltsw,deltim,lsswr, me,        &
     &                      slag,sdec,cdec,solcon)
!...................................

! =================   subprogram documentation block   ================ !
!                                                                       !
! subprogram:   radupdate   calls many update subroutines to check and  !
!   update radiation required but time varying data sets and module     !
!   variables.                                                          !
!                                                                       !
! usage:        call radupdate                                          !
!                                                                       !
! attributes:                                                           !
!   language:  fortran 90                                               !
!   machine:   ibm sp                                                   !
!                                                                       !
!  ====================  definition of variables  ====================  !
!                                                                       !
! input parameters:                                                     !
!   idate(8)       : ncep absolute date and time of initial condition   !
!                    (yr, mon, day, t-zone, hr, min, sec, mil-sec)      !
!   jdate(8)       : ncep absolute date and time at fcst time           !
!                    (yr, mon, day, t-zone, hr, min, sec, mil-sec)      !
!   deltsw         : sw radiation calling frequency in seconds          !
!   deltim         : model timestep in seconds                          !
!   lsswr          : logical flags for sw radiation calculations        !
!   me             : print control flag                                 !
!                                                                       !
!  outputs:                                                             !
!   slag           : equation of time in radians                        !
!   sdec, cdec     : sin and cos of the solar declination angle         !
!   solcon         : sun-earth distance adjusted solar constant (w/m2)  !
!                                                                       !
!  external module variables:                                           !
!   isolar   : solar constant cntrl  (in module physparam)              !
!              = 0: use the old fixed solar constant in "physcon"       !
!              =10: use the new fixed solar constant in "physcon"       !
!              = 1: use noaa ann-mean tsi tbl abs-scale with cycle apprx!
!              = 2: use noaa ann-mean tsi tbl tim-scale with cycle apprx!
!              = 3: use cmip5 ann-mean tsi tbl tim-scale with cycl apprx!
!              = 4: use cmip5 mon-mean tsi tbl tim-scale with cycl apprx!
!   ictmflg  : =yyyy#, external data ic time/date control flag          !
!              =   -2: same as 0, but superimpose seasonal cycle        !
!                      from climatology data set.                       !
!              =   -1: use user provided external data for the          !
!                      forecast time, no extrapolation.                 !
!              =    0: use data at initial cond time, if not            !
!                      available, use latest, no extrapolation.         !
!              =    1: use data at the forecast time, if not            !
!                      available, use latest and extrapolation.         !
!              =yyyy0: use yyyy data for the forecast time,             !
!                      no further data extrapolation.                   !
!              =yyyy1: use yyyy data for the fcst. if needed, do        !
!                      extrapolation to match the fcst time.            !
!                                                                       !
!  module variables:                                                    !
!   loz1st   : first-time clim ozone data read flag                     !
!                                                                       !
!  subroutines called: sol_update, aer_update, gas_update               !
!                                                                       !
!  ===================================================================  !
!
      use module_radiation_astronomy, only : sol_update
      use module_radiation_aerosols,  only : aer_update
      use module_radiation_gases,     only : gas_update

      implicit none

!  ---  inputs:
      integer, intent(in) :: idate(:), jdate(:), me
      logical, intent(in) :: lsswr

      real (kind=kind_phys), intent(in) :: deltsw, deltim

!  ---  outputs:
      real (kind=kind_phys), intent(out) :: slag, sdec, cdec, solcon

!  ---  locals:
      integer :: iyear, imon, iday, ihour
      integer :: kyear, kmon, kday, khour

      logical :: lmon_chg       ! month change flag
      logical :: lco2_chg       ! cntrl flag for updating co2 data
      logical :: lsol_chg       ! cntrl flag for updating solar constant
!
!===> ...  begin here
!
!> -# Set up time stamp at fcst time and that for green house gases
!! (currently co2 only)
!  --- ...  time stamp at fcst time

      iyear = jdate(1)
      imon  = jdate(2)
      iday  = jdate(3)
      ihour = jdate(5)

!  --- ...  set up time stamp used for green house gases (** currently co2 only)

      if ( ictmflg==0 .or. ictmflg==-2 ) then  ! get external data at initial condition time
        kyear = idate(1)
        kmon  = idate(2)
        kday  = idate(3)
        khour = idate(5)
      else                           ! get external data at fcst or specified time
        kyear = iyear
        kmon  = imon
        kday  = iday
        khour = ihour
      endif   ! end if_ictmflg_block

      if ( month0 /= imon ) then
        lmon_chg = .true.
        month0   = imon
      else
        lmon_chg = .false.
      endif

!> -# Call module_radiation_astronomy::sol_update(), yearly update, no
!! time interpolation.
      if (lsswr) then

        if ( isolar == 0 .or. isolar == 10 ) then
          lsol_chg = .false.
        elseif ( iyear0 /= iyear ) then
          lsol_chg = .true.
        else
          lsol_chg = ( isolar==4 .and. lmon_chg )
        endif
        iyear0 = iyear

        call sol_update                                                 &
!  ---  inputs:
     &     ( jdate,kyear,deltsw,deltim,lsol_chg, me,                    &
!  ---  outputs:
     &       slag,sdec,cdec,solcon                                      &
     &     )

      endif  ! end_if_lsswr_block

!> -# Call module_radiation_aerosols::aer_update(), monthly update, no
!! time interpolation
      if ( lmon_chg ) then
        call aer_update ( iyear, imon, me )
      endif

!> -# Call co2 and other gases update routine:
!! module_radiation_gases::gas_update()
      if ( monthd /= kmon ) then
        monthd = kmon
        lco2_chg = .true.
      else
        lco2_chg = .false.
      endif

      call gas_update ( kyear,kmon,kday,khour,loz1st,lco2_chg, me )

      if ( loz1st ) loz1st = .false.

!> -# Call surface update routine (currently not needed)
!     call sfc_update ( iyear, imon, me )

!> -# Call clouds update routine (currently not needed)
!     call cld_update ( iyear, imon, me )
!
      return
!...................................
      end subroutine radupdate
!-----------------------------------

!! @}
end module GFS_rrtmg_setup

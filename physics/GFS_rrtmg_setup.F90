!> \file GFS_rrtmg_setup.F90
!! This file contains

!> \defgroup GFS_rrtmg_setup_mod GFS RRTMG Scheme Setup
!! This subroutine initializes RRTMG. 
!> @{
module GFS_rrtmg_setup

   use machine, only:  kind_phys

   implicit none

   public GFS_rrtmg_setup_init, GFS_rrtmg_setup_timestep_init, GFS_rrtmg_setup_finalize

   private

   logical :: is_initialized = .false.

   !  ---  version tag and last revision date
   character(40), parameter ::                                       &
        &   VTAGRAD='NCEP-Radiation_driver    v5.2  Jan 2013 '
   !    &   VTAGRAD='NCEP-Radiation_driver    v5.1  Nov 2012 '
   !    &   VTAGRAD='NCEP-Radiation_driver    v5.0  Aug 2012 '

   !> new data input control variables (set/reset in subroutine radupdate):
   integer :: month0 = 0
   integer :: iyear0 = 0
   integer :: monthd = 0

   !> control flag for the first time of reading climatological ozone data
   !! (set/reset in subroutines radinit/radupdate, it is used only if the
   !! control parameter ntoz=0)
   logical :: loz1st = .true.

   contains

!> \section arg_table_GFS_rrtmg_setup_init Argument Table
!! \htmlinclude GFS_rrtmg_setup_init.html
!!
   subroutine GFS_rrtmg_setup_init ( si, levr, ictm, isol, solar_file, ico2, &
        iaer, ntcw, num_p3d, npdf3d, ntoz, iovr, iovr_rand, iovr_maxrand,    &
        iovr_max, iovr_dcorr, iovr_exp, iovr_exprand, icliq_sw, lcrick,      &
        lcnorm, imp_physics, lnoprec, idate, iflip, do_RRTMGP, me, lalw1bd,  &
        iaermdl, iaerflg, aeros_file, con_pi, con_t0c, con_c, con_boltz,     &
        con_plnk, con_solr_2008, con_solr_2002, con_g, con_rd, co2usr_file,  &
        co2cyc_file, rad_hr_units, inc_minor_gas, icliq_lw, isubcsw, isubclw,&
        iswmode, ipsd0, ltp, lextop, errmsg, errflg)
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
!                     =4: exponential overlap clouds                    !
!   isubcsw/isubclw: sub-column cloud approx control flag (sw/lw rad)   !
!                     =0: with out sub-column cloud approximation       !
!                     =1: mcica sub-col approx. prescribed random seed  !
!                     =2: mcica sub-col approx. provided random seed    !
!   lcrick           : control flag for eliminating CRICK               !
!   lcnorm           : control flag for in-cloud condensate mixing ratio!
!   lnoprec          : control flag for not using precip in radiation   !
!   idate(4)         : ncep absolute date and time of initial condition !
!                      (hour, month, day, year)                         !
!   iflip            : control flag for direction of vertical index     !
!                     =0: index from toa to surface                     !
!                     =1: index from surface to toa                     !
!   me               : print control flag                               !
!   ltp              : number of radiation extra top layers             !
!   lextop           : control flag to denote extra top layers are used !
!                                                                       !
!  subroutines called: radinit                                          !
!                                                                       !
!  ===================================================================  !
!
      use module_radiation_astronomy, only : sol_init
      use module_radiation_aerosols,  only : aer_init
      use module_radiation_gases,     only : gas_init
      use module_radiation_clouds,    only : cld_init
      use rrtmg_lw,                   only : rlwinit
      use rrtmg_sw,                   only : rswinit
      implicit none

      ! interface variables
      real (kind=kind_phys), intent(in) :: si(:)
      integer, intent(in) :: levr, ictm, isol, ico2, iaer, ntcw, num_p3d, &
           ltp, npdf3d, ntoz, iovr, iovr_rand, iovr_maxrand, iovr_max,    &
           iovr_dcorr, iovr_exp, iovr_exprand, icliq_sw, imp_physics,     &
           iflip, me, rad_hr_units, icliq_lw, isubcsw, isubclw, iswmode
      integer, intent(in) :: idate(:)
      logical, intent(in) :: lcrick, lcnorm, lnoprec, do_RRTMGP, lalw1bd, &
           inc_minor_gas, lextop
      character(len=26),intent(in)  :: aeros_file, solar_file, co2usr_file,&
           co2cyc_file
      real(kind_phys),  intent(in)  :: con_pi, con_t0c, con_c, con_boltz,  &
           con_plnk, con_solr_2008, con_solr_2002, con_g, con_rd
      integer,          intent(inout) :: ipsd0
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
      integer,          intent(out) :: iaermdl, iaerflg

      ! Initialize the CCPP error handling variables
      errmsg = ''
      errflg = 0

      if (is_initialized) return
      
      if (do_RRTMGP) then
        write(errmsg,'(*(a))') "Logic error: do_RRTMGP must be set to .false."
        errflg = 1
        return
      end if

      if ( ictm==0 .or. ictm==-2 ) then
        iaerflg = mod(iaer, 100)        ! no volcanic aerosols for clim hindcast
      else
        iaerflg = mod(iaer, 1000)   
      endif
      iaermdl = iaer/1000               ! control flag for aerosol scheme selection
      if ( iaermdl < 0 .or.  (iaermdl>2 .and. iaermdl/=5) ) then
         print *, ' Error -- IAER flag is incorrect, Abort'
         errflg = 1
         errmsg = 'ERROR(GFS_rrtmg_setup): IAER flag is incorrect'
         return
      endif

!  ---  assign initial permutation seed for mcica cloud-radiation
      if ( isubcsw>0 .or. isubclw>0 ) then
!       ipsd0 = 17*idate(1)+43*idate(2)+37*idate(3)+23*idate(4) + ipsd0
        ipsd0 = 17*idate(1)+43*idate(2)+37*idate(3)+23*idate(4)
      endif

      if ( me == 0 ) then
         print *,' In rad_initialize (GFS_rrtmg_setup_init), before calling RRTMG initialization'
         print *,' si =',si
         print *,' levr=',levr,' ictm=',ictm,' isol=',isol,' ico2=',ico2,&
                 ' iaermdl=',iaermdl,' iaerflg=',iaerflg
         print *,' np3d=',num_p3d,' ntoz=',ntoz,                         &
                 ' iovr=',iovr,' isubcsw=',isubcsw,                      &
                 ' isubclw=',isubclw,' icliq_sw=',icliq_sw,              &
                 ' iflip=',iflip,'  me=',me
         print *,' lcrick=',lcrick,                                      &
                 ' lcnorm=',lcnorm,' lnoprec=',lnoprec
         print *, 'lextop=',lextop, ' ltp=',ltp
      endif

      ! Call initialization routines
      call sol_init ( me, isol, solar_file, con_solr_2008,con_solr_2002,&
           con_pi )
      call aer_init ( levr, me, iaermdl, iaerflg, lalw1bd, aeros_file,  &
           con_pi, con_t0c, con_c, con_boltz, con_plnk, errflg, errmsg)
      call gas_init ( me, co2usr_file, co2cyc_file, ico2, ictm, ntoz,   &
           con_pi, errflg, errmsg)
      call cld_init ( si, levr, imp_physics, me, con_g, con_rd, errflg, errmsg)
      call rlwinit ( me, rad_hr_units, inc_minor_gas, icliq_lw, isubcsw, &
           iovr, iovr_rand, iovr_maxrand, iovr_max, iovr_dcorr,         &
           iovr_exp, iovr_exprand, errflg, errmsg )
      call rswinit ( me, rad_hr_units, inc_minor_gas, icliq_sw, isubclw, &
           iovr, iovr_rand, iovr_maxrand, iovr_max, iovr_dcorr,         &
           iovr_exp, iovr_exprand,iswmode, errflg, errmsg )

      if ( me == 0 ) then
        print *,'  Radiation sub-cloud initial seed =',ipsd0,           &
     &          ' IC-idate =',idate
        print *,' return from rad_initialize (GFS_rrtmg_setup_init) - after calling RRTMG initialization'
      endif
!
      is_initialized = .true.
!
      return

   end subroutine GFS_rrtmg_setup_init

!> \section arg_table_GFS_rrtmg_setup_timestep_init Argument Table
!! \htmlinclude GFS_rrtmg_setup_timestep_init.html
!!
   subroutine GFS_rrtmg_setup_timestep_init (idate, jdate, deltsw, deltim, &
        lsswr, me, iaermdl, iaerflg, isol, aeros_file, slag, sdec, cdec,   &
        solcon, con_pi, co2dat_file, co2gbl_file, ictm, ico2, ntoz, errmsg, errflg)

      implicit none

      ! interface variables
      integer,              intent(in)  :: idate(:)
      integer,              intent(in)  :: jdate(:)
      real(kind=kind_phys), intent(in)  :: deltsw
      real(kind=kind_phys), intent(in)  :: deltim
      real(kind=kind_phys), intent(in)  :: con_pi
      logical,              intent(in)  :: lsswr
      integer,              intent(in)  :: me
      integer,              intent(in)  :: iaermdl, iaerflg, isol, ictm, ico2, ntoz
      character(len=26),    intent(in)  :: aeros_file, co2dat_file, co2gbl_file
      real(kind=kind_phys), intent(out) :: slag
      real(kind=kind_phys), intent(out) :: sdec
      real(kind=kind_phys), intent(out) :: cdec
      real(kind=kind_phys), intent(out) :: solcon
      character(len=*),     intent(out) :: errmsg
      integer,              intent(out) :: errflg

      ! Check initialization state
      if (.not.is_initialized) then
         write(errmsg, fmt='((a))') 'GFS_rrtmg_setup_timestep_init called before GFS_rrtmg_setup_init'
         errflg = 1
         return
      end if

      ! Initialize the CCPP error handling variables
      errmsg = ''
      errflg = 0

      call radupdate(idate,jdate,deltsw,deltim,lsswr,me,iaermdl, iaerflg,isol,aeros_file,&
           slag,sdec,cdec,solcon,con_pi,co2dat_file,co2gbl_file,ictm,ico2,ntoz,errflg,errmsg)

   end subroutine GFS_rrtmg_setup_timestep_init

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
!-----------------------------------
      subroutine radupdate( idate,jdate,deltsw,deltim,lsswr,me, iaermdl,&
           iaerflg, isol, aeros_file, slag,sdec,cdec,solcon, con_pi,    &
           co2dat_file,co2gbl_file, ictm, ico2, ntoz, errflg, errmsg)
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
!  subroutines called: sol_update, aer_update, gas_update               !
!                                                                       !
!  ===================================================================  !
!
      use module_radiation_astronomy, only : sol_update
      use module_radiation_aerosols,  only : aer_update
      use module_radiation_gases,     only : gas_update

      implicit none

!  ---  inputs:
      integer, intent(in) :: idate(:), jdate(:), me, iaermdl, iaerflg, isol, ictm, ntoz, ico2
      logical, intent(in) :: lsswr
      character(len=26),intent(in) :: aeros_file,co2dat_file,co2gbl_file

      real (kind=kind_phys), intent(in) :: deltsw, deltim, con_pi

!  ---  outputs:
      real (kind=kind_phys), intent(out) :: slag, sdec, cdec, solcon
      character(len=*),     intent(out) :: errmsg
      integer,              intent(out) :: errflg

!  ---  locals:
      integer :: iyear, imon, iday, ihour
      integer :: kyear, kmon, kday, khour

      logical :: lmon_chg       ! month change flag
      logical :: lco2_chg       ! cntrl flag for updating co2 data
      logical :: lsol_chg       ! cntrl flag for updating solar constant
!
!===> ...  begin here
!

      ! Initialize the CCPP error handling variables
      errmsg = ''
      errflg = 0

!> -# Set up time stamp at fcst time and that for green house gases
!! (currently co2 only)
!  --- ...  time stamp at fcst time

      iyear = jdate(1)
      imon  = jdate(2)
      iday  = jdate(3)
      ihour = jdate(5)

!  --- ...  set up time stamp used for green house gases (** currently co2 only)

      if ( ictm==0 .or. ictm==-2 ) then  ! get external data at initial condition time
        kyear = idate(1)
        kmon  = idate(2)
        kday  = idate(3)
        khour = idate(5)
      else                           ! get external data at fcst or specified time
        kyear = iyear
        kmon  = imon
        kday  = iday
        khour = ihour
      endif   ! end if_ictm_block

      if ( month0 /= imon ) then
        lmon_chg = .true.
        month0   = imon
      else
        lmon_chg = .false.
      endif

!> -# Call module_radiation_astronomy::sol_update(), yearly update, no
!! time interpolation.
      if (lsswr) then

        if ( isol == 0 .or. isol == 10 ) then
          lsol_chg = .false.
        elseif ( iyear0 /= iyear ) then
          lsol_chg = .true.
        else
          lsol_chg = ( isol==4 .and. lmon_chg )
        endif
        iyear0 = iyear

        call sol_update                                                 &
!  ---  inputs:
     &     ( jdate,kyear,deltsw,deltim,lsol_chg, me,                    &
!  ---  outputs:
     &       slag,sdec,cdec,solcon,con_pi,errmsg,errflg                 &
     &     )

      endif  ! end_if_lsswr_block

!> -# Call module_radiation_aerosols::aer_update(), monthly update, no
!! time interpolation
      if ( lmon_chg ) then
        call aer_update ( iyear, imon, me, iaermdl, aeros_file, errflg, errmsg )
      endif

!> -# Call co2 and other gases update routine:
!! module_radiation_gases::gas_update()
      if ( monthd /= kmon ) then
        monthd = kmon
        lco2_chg = .true.
      else
        lco2_chg = .false.
      endif

      call gas_update ( kyear,kmon,kday,khour,loz1st,lco2_chg, me, co2dat_file, &
           co2gbl_file, ictm, ico2, ntoz, errflg, errmsg )

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
!> @}

end module GFS_rrtmg_setup

!> \file GFS_rrtmgp_setup.f90
!! This file contains
module GFS_rrtmgp_setup
  use machine,                    only : kind_phys
  use module_radiation_astronomy, only : sol_init, sol_update
  use module_radiation_aerosols,  only : aer_init, aer_update
  use module_radiation_gases,     only : gas_init, gas_update
  use module_radiation_surface,   only : sfc_init
  use GFS_cloud_diagnostics,      only : hml_cloud_diagnostics_initialize
  ! *NOTE* These parameters below are required radiation_****** modules. They are not
  !        directly used by the RRTMGP routines.
  use physparam,                  only : isolar,  ictmflg, ico2flg, ioznflg, iaerflg,    &
                                         iaermdl, ialbflg, iemsflg, ivflip
  implicit none
  
  public GFS_rrtmgp_setup_init, GFS_rrtmgp_setup_run, GFS_rrtmgp_setup_finalize
  
  ! Version tag and last revision date
  character(40), parameter ::                                       &
       VTAGRAD='NCEP-RRTMGP_driver       v1.0  Sep 2019 '
  
  ! Module paramaters
  integer ::  &
       month0 = 0, &
       iyear0 = 0, &
       monthd = 0
  logical ::  &
       is_initialized = .false.
  ! Control flag for the first time of reading climatological ozone data
  ! (set/reset in subroutines GFS_rrtmgp_setup_init/GFS_rrtmgp_setuup_run, it is used only if 
  ! the control parameter ioznflg=0)
  logical :: loz1st = .true.
  
contains
  ! #########################################################################################  
  ! SUBROUTINE GFS_rrtmgp_setup_init
  ! #########################################################################################  
!> \defgroup GFS_rrtmgp_setup GFS RRTMGP Scheme Setup
!! @{
!! \section arg_table_GFS_rrtmgp_setup_init
!! \htmlinclude GFS_rrtmgp_setup_init.html
!!
  subroutine GFS_rrtmgp_setup_init(imp_physics, imp_physics_fer_hires, imp_physics_gfdl,    &
       imp_physics_thompson, imp_physics_wsm6, imp_physics_zhao_carr,                       &
       imp_physics_zhao_carr_pdf, imp_physics_mg,  si, levr, ictm, isol, ico2, iaer, ialb,  &
       iems, ntcw, num_p3d,  ntoz, iovr, isubc_sw, isubc_lw, icliq_sw, crick_proof, ccnorm, &
       norad_precip, idate, iflip, me, errmsg, errflg)

    ! Inputs
    integer, intent(in) :: &
         imp_physics,               & ! Flag for MP scheme
         imp_physics_fer_hires,     & ! Flag for fer-hires scheme
         imp_physics_gfdl,          & ! Flag for gfdl scheme
         imp_physics_thompson,      & ! Flag for thompsonscheme
         imp_physics_wsm6,          & ! Flag for wsm6 scheme
         imp_physics_zhao_carr,     & ! Flag for zhao-carr scheme
         imp_physics_zhao_carr_pdf, & ! Flag for zhao-carr+PDF scheme
         imp_physics_mg               ! Flag for MG scheme
    real(kind_phys), dimension(levr+1), intent(in) :: &
         si
    integer, intent(in) :: levr, ictm, isol, ico2, iaer, ialb, iems,   & 
         ntcw, num_p3d, ntoz, iovr, isubc_sw, isubc_lw,                &
         icliq_sw, iflip, me 
    logical, intent(in) :: &
         crick_proof, ccnorm, norad_precip
    integer, intent(in), dimension(4) :: &
         idate

    ! Outputs
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg
    
    ! Initialize the CCPP error handling variables
    errmsg = ''
    errflg = 0
    
    if (is_initialized) return
    
    ! Set radiation parameters
    isolar  = isol                     ! solar constant control flag
    ictmflg = ictm                     ! data ic time/date control flag
    ico2flg = ico2                     ! co2 data source control flag
    ioznflg = ntoz                     ! ozone data source control flag
    ialbflg = ialb                     ! surface albedo control flag
    iemsflg = iems                     ! surface emissivity control flag
    ivflip  = iflip                    ! vertical index direction control flag
    
    if ( ictm==0 .or. ictm==-2 ) then
       iaerflg = mod(iaer, 100)        ! no volcanic aerosols for clim hindcast
    else
       iaerflg = mod(iaer, 1000)   
    endif
    iaermdl = iaer/1000               ! control flag for aerosol scheme selection
    if ( iaermdl < 0 .or.  (iaermdl>2 .and. iaermdl/=5) ) then
       errmsg = trim(errmsg) // ' Error -- IAER flag is incorrect, Abort'
       errflg = 1
       return
    endif
    
    if ( me == 0 ) then
       print *,'  In rad_initialize (GFS_rrtmgp_setup_init), before calling radinit'
       print *,' si       = ',si
       print *,' levr     = ',levr,      &
               ' ictm     = ',ictm,      &
               ' isol     = ',isol,      &
               ' ico2     = ',ico2,      &
               ' iaer     = ',iaer,      &
               ' ialb     = ',ialb,      &
               ' iems     = ',iems,      &
               ' ntcw     = ',ntcw
       print *,' np3d     = ',num_p3d,   &
               ' ntoz     = ',ntoz,      &
               ' iovr     = ',iovr,      &
               ' isubc_sw = ',isubc_sw,  &
               ' isubc_lw = ',isubc_lw,  &
               ' icliq_sw = ',icliq_sw,  &
               ' iflip    = ',iflip,     &
               ' me       = ',me
    endif
    
#if 0
    ! GFS_radiation_driver.F90 may in the future initialize air/ground
    ! temperature differently; however, this is not used at the moment
    ! and as such we avoid the difficulty of dealing with exchanging
    ! itsfc between GFS_rrtmgp_setup and a yet-to-be-created/-used
    ! interstitial routine (or GFS_radiation_driver.F90)
    itsfc  = iemsflg / 10             ! sfc air/ground temp control
#endif
    loz1st = (ioznflg == 0)           ! first-time clim ozone data read flag
    month0 = 0
    iyear0 = 0
    monthd = 0

    ! Call initialization routines..
    call sol_init ( me )
    call aer_init ( levr, me )
    call gas_init ( me )
    call sfc_init ( me )
    call hml_cloud_diagnostics_initialize(imp_physics, imp_physics_fer_hires,           &
         imp_physics_gfdl, imp_physics_thompson, imp_physics_wsm6,                      &
         imp_physics_zhao_carr, imp_physics_zhao_carr_pdf, imp_physics_mg, levr, me, si,&
         errflg)

    if ( me == 0 ) then
       print *,' return from rad_initialize (GFS_rrtmgp_setup_init) - after calling radinit'
    endif
    
    is_initialized = .true.

    return
  end subroutine GFS_rrtmgp_setup_init

  ! #########################################################################################
  ! SUBROUTINE GFS_rrtmgp_setup_run
  ! #########################################################################################
!> \section arg_table_GFS_rrtmgp_setup_run
!! \htmlinclude GFS_rrtmgp_setup_run.html
!!
  subroutine GFS_rrtmgp_setup_run (idate, jdate, deltsw, deltim, lsswr, me, &
       slag, sdec, cdec, solcon, errmsg, errflg)
     
    ! Inputs
    integer,         intent(in)  :: idate(:)
    integer,         intent(in)  :: jdate(:)
    real(kind_phys), intent(in)  :: deltsw
    real(kind_phys), intent(in)  :: deltim
    logical,         intent(in)  :: lsswr
    integer,         intent(in)  :: me

    ! Outputs
    real(kind_phys), intent(out) :: slag
    real(kind_phys), intent(out) :: sdec
    real(kind_phys), intent(out) :: cdec
    real(kind_phys), intent(out) :: solcon
    character(len=*),intent(out) :: errmsg
    integer,         intent(out) :: errflg

    ! Locals
    integer :: iyear, imon, iday, ihour
    integer :: kyear, kmon, kday, khour
    logical :: lmon_chg       ! month change flag
    logical :: lco2_chg       ! cntrl flag for updating co2 data
    logical :: lsol_chg       ! cntrl flag for updating solar constant
    
    ! Initialize the CCPP error handling variables
    errmsg = ''
    errflg = 0

    ! Check initialization state
    if (.not.is_initialized) then
       write(errmsg, fmt='((a))') 'GFS_rrtmgp_setup_run called before GFS_rrtmgp_setup_init'
       errflg = 1
       return
    end if

    ! Set up time stamp at fcst time and that for green house gases
    iyear = jdate(1)
    imon  = jdate(2)
    iday  = jdate(3)
    ihour = jdate(5)

    ! Set up time stamp used for green house gases (** currently co2 only)
    ! get external data at initial condition time 
    if ( ictmflg==0 .or. ictmflg==-2 ) then 
       kyear = idate(1)
       kmon  = idate(2)
       kday  = idate(3)
       khour = idate(5)
    ! get external data at fcst or specified time 
    else
       kyear = iyear
       kmon  = imon
       kday  = iday
       khour = ihour
    endif
    
    if ( month0 /= imon ) then
       lmon_chg = .true.
       month0   = imon
    else
       lmon_chg = .false.
    endif

    ! Update solar forcing...
    if (lsswr) then
       if ( isolar == 0 .or. isolar == 10 ) then
          lsol_chg = .false.
       elseif ( iyear0 /= iyear ) then
          lsol_chg = .true.
       else
          lsol_chg = ( isolar==4 .and. lmon_chg )
       endif
       iyear0 = iyear
       call sol_update(jdate, kyear, deltsw, deltim, lsol_chg, me, slag, sdec, cdec, solcon)
    endif

    ! Update aerosols...
    if ( lmon_chg ) then
       call aer_update ( iyear, imon, me )
    endif

    ! Update trace gases (co2 only)...
    if ( monthd /= kmon ) then
       monthd = kmon
       lco2_chg = .true.
    else
       lco2_chg = .false.
    endif
    call gas_update (kyear, kmon, kday, khour, loz1st, lco2_chg, me )
    
    if ( loz1st ) loz1st = .false.

    return
  end subroutine GFS_rrtmgp_setup_run

  ! #########################################################################################
  ! SUBROUTINE GFS_rrtmgp_setup_finalize
  ! ######################################################################################### 
!> \section arg_table_GFS_rrtmgp_setup_finalize
!! \htmlinclude GFS_rrtmgp_setup_finalize.html
!!
  subroutine GFS_rrtmgp_setup_finalize (errmsg, errflg)
    character(len=*),          intent(  out) :: errmsg
    integer,                   intent(  out) :: errflg
    
    ! Initialize the CCPP error handling variables
    errmsg = ''
    errflg = 0
     
    if (.not.is_initialized) return
    
    ! do finalization stuff if needed
    is_initialized = .false.
    
  end subroutine GFS_rrtmgp_setup_finalize
end module GFS_rrtmgp_setup

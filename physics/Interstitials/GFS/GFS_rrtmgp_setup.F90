!> \file GFS_rrtmgp_setup.F90
!! This file initializes the RRTMGP radiation scheme

module GFS_rrtmgp_setup
  use machine,                    only : kind_phys
  use module_radiation_astronomy, only : sol_init, sol_update
  use module_radiation_aerosols,  only : aer_init, aer_update
  use module_radiation_gases,     only : gas_init, gas_update
  implicit none
  
  public GFS_rrtmgp_setup_init, GFS_rrtmgp_setup_timestep_init, GFS_rrtmgp_setup_finalize

  private

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
  ! (set/reset in subroutines GFS_rrtmgp_setup_init/GFS_rrtmgp_setup_timestep_init, it is used only if 
  ! the control parameter ntoz=0)
  logical :: loz1st = .true.
  
contains

!> \defgroup GFS_rrtmgp_setup_mod GFS RRTMGP Scheme Setup Module
!! \section arg_table_GFS_rrtmgp_setup_init
!! \htmlinclude GFS_rrtmgp_setup_init.html
!!
  subroutine GFS_rrtmgp_setup_init(do_RRTMGP, imp_physics, imp_physics_fer_hires,        &
       imp_physics_gfdl, imp_physics_thompson, imp_physics_wsm6, imp_physics_zhao_carr,  &
       imp_physics_zhao_carr_pdf, imp_physics_mg,  si, levr, ictm, isol, ico2, iaer,     &
       ntcw, ntoz, iovr, isubc_sw, isubc_lw, lalw1bd, idate, me, aeros_file,             &
       iaermdl, iaerflg, con_pi, con_t0c, con_c, con_boltz, con_plnk, solar_file,        &
       con_solr_2008, con_solr_2002, co2usr_file, co2cyc_file, ipsd0, errmsg, errflg)

    ! Inputs
    logical, intent(in) :: do_RRTMGP
    integer, intent(in) :: &
         imp_physics,               & ! Flag for MP scheme
         imp_physics_fer_hires,     & ! Flag for fer-hires scheme
         imp_physics_gfdl,          & ! Flag for gfdl scheme
         imp_physics_thompson,      & ! Flag for thompsonscheme
         imp_physics_wsm6,          & ! Flag for wsm6 scheme
         imp_physics_zhao_carr,     & ! Flag for zhao-carr scheme
         imp_physics_zhao_carr_pdf, & ! Flag for zhao-carr+PDF scheme
         imp_physics_mg               ! Flag for MG scheme
    real(kind_phys), intent(in) :: &
         con_pi, con_t0c, con_c, con_boltz, con_plnk, con_solr_2008, con_solr_2002
    real(kind_phys), dimension(:), intent(in) :: &
         si
    integer, intent(in) :: levr, ictm, isol, ico2, iaer, & 
         ntcw, ntoz, iovr, isubc_sw, isubc_lw,  &
         me
    logical, intent(in) :: &
         lalw1bd
    integer, intent(in), dimension(:) :: &
         idate
    character(len=26),intent(in) :: aeros_file, solar_file, co2usr_file, co2cyc_file

    ! Outputs
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errflg
    integer,          intent(inout) :: ipsd0
    integer,          intent(out)   :: iaermdl, iaerflg
    
    ! Initialize the CCPP error handling variables
    errmsg = ''
    errflg = 0
    
    if (is_initialized) return

    ! Consistency checks
    if (.not. do_RRTMGP) then
      write(errmsg,'(*(a))') "Logic error: do_RRTMGP must be set to .true."
      errflg = 1
      return
    end if

    ! Set radiation parameters
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

    ! Assign initial permutation seed for mcica cloud-radiation
    if ( isubc_sw>0 .or. isubc_lw>0 ) then
       ipsd0 = 17*idate(1)+43*idate(2)+37*idate(3)+23*idate(4)
    endif
    
    if ( me == 0 ) then
       print *,'  In rad_initialize (GFS_rrtmgp_setup_init), before calling radinit'
       print *,' si       = ',si
       print *,' levr     = ',levr,      &
               ' ictm     = ',ictm,      &
               ' isol     = ',isol,      &
               ' ico2     = ',ico2,      &
               ' iaermdl  = ',iaermdl,   &
               ' iaerflg  = ',iaerflg,   &
               ' ntcw     = ',ntcw,      &
               ' ntoz     = ',ntoz,      &
               ' iovr     = ',iovr,      &
               ' isubc_sw = ',isubc_sw,  &
               ' isubc_lw = ',isubc_lw,  &
               ' ipsd0    = ',ipsd0,     &
               ' me       = ',me
    endif

    loz1st = (ntoz == 0)           ! first-time clim ozone data read flag
    month0 = 0
    iyear0 = 0
    monthd = 0

    ! Call initialization routines..
    call sol_init ( me, isol, solar_file, con_solr_2008, con_solr_2002, con_pi )
    call aer_init ( levr, me, iaermdl, iaerflg, lalw1bd, aeros_file, con_pi, con_t0c,    &
         con_c, con_boltz, con_plnk, errflg, errmsg)
    call gas_init ( me, co2usr_file, co2cyc_file, ico2, ictm, ntoz, con_pi, errflg, errmsg )

    if ( me == 0 ) then
       print *,' return from rad_initialize (GFS_rrtmgp_setup_init) - after calling radinit'
    endif
    
    is_initialized = .true.

    return
  end subroutine GFS_rrtmgp_setup_init

  ! #########################################################################################
  ! SUBROUTINE GFS_rrtmgp_setup_timestep_init
  ! #########################################################################################
!> \section arg_table_GFS_rrtmgp_setup_timestep_init
!! \htmlinclude GFS_rrtmgp_setup_timestep_init.html
!!
  subroutine GFS_rrtmgp_setup_timestep_init (idate, jdate, deltsw, deltim, doSWrad, me,     &
       iaermdl, aeros_file, isol, slag, sdec, cdec, solcon, con_pi, co2dat_file,            &
       co2gbl_file, ictm, ico2, ntoz, errmsg, errflg)
     
    ! Inputs
    integer,         intent(in)  :: idate(:)
    integer,         intent(in)  :: jdate(:)
    real(kind_phys), intent(in)  :: deltsw
    real(kind_phys), intent(in)  :: deltim
    logical,         intent(in)  :: doSWrad
    real(kind_phys), intent(in)  :: con_pi
    integer,         intent(in)  :: me
    integer,         intent(in)  :: iaermdl,isol,ictm,ico2,ntoz
    character(len=26), intent(in) :: aeros_file,co2dat_file,co2gbl_file

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
       write(errmsg, fmt='((a))') 'GFS_rrtmgp_setup_timestep_init called before GFS_rrtmgp_setup_init'
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
    if ( ictm==0 .or. ictm==-2 ) then 
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
    if (doSWrad) then
       if ( isol == 0 .or. isol == 10 ) then
          lsol_chg = .false.
       elseif ( iyear0 /= iyear ) then
          lsol_chg = .true.
       else
          lsol_chg = ( isol==4 .and. lmon_chg )
       endif
       iyear0 = iyear
       call sol_update(jdate, kyear, deltsw, deltim, lsol_chg, me, slag, sdec, cdec, solcon, con_pi, errmsg, errflg)
    endif

    ! Update aerosols...
    if ( lmon_chg ) then
       call aer_update ( iyear, imon, me, iaermdl, aeros_file, errflg, errmsg)
    endif

    ! Update trace gases (co2 only)...
    if ( monthd /= kmon ) then
       monthd = kmon
       lco2_chg = .true.
    else
       lco2_chg = .false.
    endif
    call gas_update (kyear, kmon, kday, khour, loz1st, lco2_chg, me, co2dat_file,           &
         co2gbl_file, ictm, ico2, ntoz, errflg, errmsg )
    
    if ( loz1st ) loz1st = .false.

    return
  end subroutine GFS_rrtmgp_setup_timestep_init

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

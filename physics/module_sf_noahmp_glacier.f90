!>  \file module_sf_noahmp_glacier.f90
!!  This file contains the NoahMP Glacier scheme.

!>\ingroup NoahMP_LSM
module noahmp_glacier_globals

  implicit none

! ==================================================================================================
!------------------------------------------------------------------------------------------!
! physical constants:                                                                      !
!------------------------------------------------------------------------------------------!

  real, parameter :: grav   = 9.80616   !acceleration due to gravity (m/s2)
  real, parameter :: sb     = 5.67e-08  !stefan-boltzmann constant (w/m2/k4)
  real, parameter :: vkc    = 0.40      !von karman constant
  real, parameter :: tfrz   = 273.16    !freezing/melting point (k)
  real, parameter :: hsub   = 2.8440e06 !latent heat of sublimation (j/kg)
  real, parameter :: hvap   = 2.5104e06 !latent heat of vaporization (j/kg)
  real, parameter :: hfus   = 0.3336e06 !latent heat of fusion (j/kg)
  real, parameter :: cwat   = 4.188e06  !specific heat capacity of water (j/m3/k)
  real, parameter :: cice   = 2.094e06  !specific heat capacity of ice (j/m3/k)
  real, parameter :: cpair  = 1004.64   !heat capacity dry air at const pres (j/kg/k)
  real, parameter :: tkwat  = 0.6       !thermal conductivity of water (w/m/k)
  real, parameter :: tkice  = 2.2       !thermal conductivity of ice (w/m/k)
  real, parameter :: tkair  = 0.023     !thermal conductivity of air (w/m/k)
  real, parameter :: rair   = 287.04    !gas constant for dry air (j/kg/k)
  real, parameter :: rw     = 461.269   !gas constant for  water vapor (j/kg/k)
  real, parameter :: denh2o = 1000.     !density of water (kg/m3)
  real, parameter :: denice = 917.      !density of ice (kg/m3)

! =====================================options for different schemes================================
! options for dynamic vegetation: 
! 1 -> off (use table lai; use fveg = shdfac from input)
! 2 -> on (together with opt_crs = 1)
! 3 -> off (use table lai; calculate fveg)
! 4 -> off (use table lai; use maximum vegetation fraction)

  integer :: dveg    != 2   !

! options for canopy stomatal resistance
! 1-> ball-berry; 2->jarvis

  integer :: opt_crs != 1    !(must 1 when dveg = 2)

! options for soil moisture factor for stomatal resistance
! 1-> noah (soil moisture) 
! 2-> clm  (matric potential)
! 3-> ssib (matric potential)

  integer :: opt_btr != 1    !(suggested 1)

! options for runoff and groundwater
! 1 -> topmodel with groundwater (niu et al. 2007 jgr) ;
! 2 -> topmodel with an equilibrium water table (niu et al. 2005 jgr) ;
! 3 -> original surface and subsurface runoff (free drainage)
! 4 -> bats surface and subsurface runoff (free drainage)

  integer :: opt_run != 1    !(suggested 1)

! options for surface layer drag coeff (ch & cm)
! 1->m-o ; 2->original noah (chen97); 3->myj consistent; 4->ysu consistent. 

  integer :: opt_sfc != 1    !(1 or 2 or 3 or 4)

! options for supercooled liquid water (or ice fraction)
! 1-> no iteration (niu and yang, 2006 jhm); 2: koren's iteration 

  integer :: opt_frz != 1    !(1 or 2)

! options for frozen soil permeability
! 1 -> linear effects, more permeable (niu and yang, 2006, jhm)
! 2 -> nonlinear effects, less permeable (old)

  integer :: opt_inf != 1    !(suggested 1)

! options for radiation transfer
! 1 -> modified two-stream (gap = f(solar angle, 3d structure ...)<1-fveg)
! 2 -> two-stream applied to grid-cell (gap = 0)
! 3 -> two-stream applied to vegetated fraction (gap=1-fveg)

  integer :: opt_rad != 1    !(suggested 1)

! options for ground snow surface albedo
! 1-> bats; 2 -> class

  integer :: opt_alb != 2    !(suggested 2)

! options for partitioning  precipitation into rainfall & snowfall
! 1 -> jordan (1991); 2 -> bats: when sfctmp<tfrz+2.2 ; 3-> sfctmp<tfrz

  integer :: opt_snf != 1    !(suggested 1)

! options for lower boundary condition of soil temperature
! 1 -> zero heat flux from bottom (zbot and tbot not used)
! 2 -> tbot at zbot (8m) read from a file (original noah)

  integer :: opt_tbot != 2   !(suggested 2)

! options for snow/soil temperature time scheme (only layer 1)
! 1 -> semi-implicit; 2 -> full implicit (original noah)

  integer :: opt_stc != 1    !(suggested 1)

! adjustable parameters for snow processes

  real, parameter :: z0sno  = 0.002  !snow surface roughness length (m) (0.002)
  real, parameter :: ssi    = 0.03   !liquid water holding capacity for snowpack (m3/m3) (0.03)
  real, parameter :: swemx  = 1.00   !new snow mass to fully cover old snow (mm)
                                     !equivalent to 10mm depth (density = 100 kg/m3)

!------------------------------------------------------------------------------------------!
end module noahmp_glacier_globals
!------------------------------------------------------------------------------------------!

!>\ingroup NoahMP_LSM
module noahmp_glacier_routines
  use noahmp_glacier_globals
#ifndef CCPP
  use  module_wrf_utl
#endif
  implicit none

  public  :: noahmp_options_glacier
  public  :: noahmp_glacier

  private :: atm_glacier
  private :: energy_glacier
  private ::       thermoprop_glacier
  private ::               csnow_glacier
  private ::       radiation_glacier
  private ::               snow_age_glacier
  private ::               snowalb_bats_glacier  
  private ::               snowalb_class_glacier
  private ::       glacier_flux
  private ::               sfcdif1_glacier                  
  private ::       tsnosoi_glacier
  private ::               hrt_glacier
  private ::               hstep_glacier   
  private ::                         rosr12_glacier
  private ::       phasechange_glacier

  private :: water_glacier
  private ::       snowwater_glacier
  private ::               snowfall_glacier
  private ::               combine_glacier
  private ::               divide_glacier
  private ::                         combo_glacier
  private ::               compact_glacier
  private ::               snowh2o_glacier

  private :: error_glacier

contains
!
! ==================================================================================================

!>\ingroup NoahMP_LSM
  subroutine noahmp_glacier (&
                   iloc    ,jloc    ,cosz    ,nsnow   ,nsoil   ,dt      , & ! in : time/space/model-related
                   sfctmp  ,sfcprs  ,uu      ,vv      ,q2      ,soldn   , & ! in : forcing
                   prcp    ,lwdn    ,tbot    ,zlvl    ,ficeold ,zsoil   , & ! in : forcing
                   qsnow   ,sneqvo  ,albold  ,cm      ,ch      ,isnow   , & ! in/out : 
                   sneqv   ,smc     ,zsnso   ,snowh   ,snice   ,snliq   , & ! in/out :
                   tg      ,stc     ,sh2o    ,tauss   ,qsfc    ,          & ! in/out : 
                   fsa     ,fsr     ,fira    ,fsh     ,fgev    ,ssoil   , & ! out : 
                   trad    ,edir    ,runsrf  ,runsub  ,sag     ,albedo  , & ! out :
                   qsnbot  ,ponding ,ponding1,ponding2,t2m     ,q2e     , & ! out :
#ifdef CCPP
                   emissi,  fpice   ,ch2b    , esnow, errmsg, errflg) 
#else
                   emissi,  fpice   ,ch2b    , esnow) 
#endif
                   

! --------------------------------------------------------------------------------------------------
! initial code: guo-yue niu, oct. 2007
! modified to glacier: michael barlage, june 2012
! --------------------------------------------------------------------------------------------------
  implicit none
! --------------------------------------------------------------------------------------------------
! input
  integer                        , intent(in)    :: iloc   !grid index
  integer                        , intent(in)    :: jloc   !grid index
  real                           , intent(in)    :: cosz   !cosine solar zenith angle [0-1]
  integer                        , intent(in)    :: nsnow  !maximum no. of snow layers        
  integer                        , intent(in)    :: nsoil  !no. of soil layers        
  real                           , intent(in)    :: dt     !time step [sec]
  real                           , intent(in)    :: sfctmp !surface air temperature [k]
  real                           , intent(in)    :: sfcprs !pressure (pa)
  real                           , intent(in)    :: uu     !wind speed in eastward dir (m/s)
  real                           , intent(in)    :: vv     !wind speed in northward dir (m/s)
  real                           , intent(in)    :: q2     !mixing ratio (kg/kg) lowest model layer
  real                           , intent(in)    :: soldn  !downward shortwave radiation (w/m2)
  real                           , intent(in)    :: prcp   !precipitation rate (kg m-2 s-1)
  real                           , intent(in)    :: lwdn   !downward longwave radiation (w/m2)
  real                           , intent(in)    :: tbot   !bottom condition for soil temp. [k]
  real                           , intent(in)    :: zlvl   !reference height (m)
  real, dimension(-nsnow+1:    0), intent(in)    :: ficeold!ice fraction at last timestep
  real, dimension(       1:nsoil), intent(in)    :: zsoil  !layer-bottom depth from soil surf (m)


! input/output : need arbitary intial values
  real                           , intent(inout) :: qsnow  !snowfall [mm/s]
  real                           , intent(inout) :: sneqvo !snow mass at last time step (mm)
  real                           , intent(inout) :: albold !snow albedo at last time step (class type)
  real                           , intent(inout) :: cm     !momentum drag coefficient
  real                           , intent(inout) :: ch     !sensible heat exchange coefficient

! prognostic variables
  integer                        , intent(inout) :: isnow  !actual no. of snow layers [-]
  real                           , intent(inout) :: sneqv  !snow water eqv. [mm]
  real, dimension(       1:nsoil), intent(inout) :: smc    !soil moisture (ice + liq.) [m3/m3]
  real, dimension(-nsnow+1:nsoil), intent(inout) :: zsnso  !layer-bottom depth from snow surf [m]
  real                           , intent(inout) :: snowh  !snow height [m]
  real, dimension(-nsnow+1:    0), intent(inout) :: snice  !snow layer ice [mm]
  real, dimension(-nsnow+1:    0), intent(inout) :: snliq  !snow layer liquid water [mm]
  real                           , intent(inout) :: tg     !ground temperature (k)
  real, dimension(-nsnow+1:nsoil), intent(inout) :: stc    !snow/soil temperature [k]
  real, dimension(       1:nsoil), intent(inout) :: sh2o   !liquid soil moisture [m3/m3]
  real                           , intent(inout) :: tauss  !non-dimensional snow age
  real                           , intent(inout) :: qsfc   !mixing ratio at lowest model layer

! output
  real                           , intent(out)   :: fsa    !total absorbed solar radiation (w/m2)
  real                           , intent(out)   :: fsr    !total reflected solar radiation (w/m2)
  real                           , intent(out)   :: fira   !total net lw rad (w/m2)  [+ to atm]
  real                           , intent(out)   :: fsh    !total sensible heat (w/m2) [+ to atm]
  real                           , intent(out)   :: fgev   !ground evap heat (w/m2) [+ to atm]
  real                           , intent(out)   :: ssoil  !ground heat flux (w/m2)   [+ to soil]
  real                           , intent(out)   :: trad   !surface radiative temperature (k)
  real                           , intent(out)   :: edir   !soil surface evaporation rate (mm/s]
  real                           , intent(out)   :: runsrf !surface runoff [mm/s] 
  real                           , intent(out)   :: runsub !baseflow (saturation excess) [mm/s]
  real                           , intent(out)   :: sag    !solar rad absorbed by ground (w/m2)
  real                           , intent(out)   :: albedo !surface albedo [-]
  real                           , intent(out)   :: qsnbot !snowmelt [mm/s]
  real                           , intent(out)   :: ponding!surface ponding [mm]
  real                           , intent(out)   :: ponding1!surface ponding [mm]
  real                           , intent(out)   :: ponding2!surface ponding [mm]
  real                           , intent(out)   :: t2m     !2-m air temperature over bare ground part [k]
  real                           , intent(out)   :: q2e
  real                           , intent(out)   :: emissi
  real                           , intent(out)   :: fpice
  real                           , intent(out)   :: ch2b
  real                           , intent(out)   :: esnow

#ifdef CCPP  
  character(len=*), intent(inout)    :: errmsg
  integer,          intent(inout)    :: errflg
#endif
  
! local
  integer                                        :: iz     !do-loop index
  integer, dimension(-nsnow+1:nsoil)             :: imelt  !phase change index [1-melt; 2-freeze]
  real                                           :: rhoair !density air (kg/m3)
  real, dimension(-nsnow+1:nsoil)                :: dzsnso !snow/soil layer thickness [m]
  real                                           :: thair  !potential temperature (k)
  real                                           :: qair   !specific humidity (kg/kg) (q2/(1+q2))
  real                                           :: eair   !vapor pressure air (pa)
  real, dimension(       1:    2)                :: solad  !incoming direct solar rad (w/m2)
  real, dimension(       1:    2)                :: solai  !incoming diffuse solar rad (w/m2)
  real, dimension(       1:nsoil)                :: sice   !soil ice content (m3/m3)
  real, dimension(-nsnow+1:    0)                :: snicev !partial volume ice of snow [m3/m3]
  real, dimension(-nsnow+1:    0)                :: snliqv !partial volume liq of snow [m3/m3]
  real, dimension(-nsnow+1:    0)                :: epore  !effective porosity [m3/m3]
  real                                           :: qdew   !ground surface dew rate [mm/s]
  real                                           :: qvap   !ground surface evap. rate [mm/s]
  real                                           :: lathea !latent heat [j/kg]
  real                                           :: qmelt  !internal pack melt
  real                                           :: swdown !downward solar [w/m2]
  real                                           :: beg_wb !beginning water for error check
  real                                           :: zbot = -8.0 

  character*256 message

! --------------------------------------------------------------------------------------------------
! re-process atmospheric forcing

   call atm_glacier (sfcprs ,sfctmp ,q2     ,soldn  ,cosz   ,thair  , & 
                     qair   ,eair   ,rhoair ,solad  ,solai  ,swdown )

   beg_wb = sneqv

! snow/soil layer thickness (m); interface depth: zsnso < 0; layer thickness dzsnso > 0

     do iz = isnow+1, nsoil
         if(iz == isnow+1) then
           dzsnso(iz) = - zsnso(iz)
         else
           dzsnso(iz) = zsnso(iz-1) - zsnso(iz)
         end if
     end do

! compute energy budget (momentum & energy fluxes and phase changes) 

    call energy_glacier (nsnow  ,nsoil  ,isnow  ,dt     ,qsnow  ,rhoair , & !in
                         eair   ,sfcprs ,qair   ,sfctmp ,lwdn   ,uu     , & !in
                         vv     ,solad  ,solai  ,cosz   ,zlvl   ,         & !in
                         tbot   ,zbot   ,zsnso  ,dzsnso ,                 & !in
                         tg     ,stc    ,snowh  ,sneqv  ,sneqvo ,sh2o   , & !inout
                         smc    ,snice  ,snliq  ,albold ,cm     ,ch     , & !inout
#ifdef CCPP
                         tauss  ,qsfc   ,errmsg ,errflg ,                 & !inout
#else
                         tauss  ,qsfc   ,                                 & !inout
#endif
                         imelt  ,snicev ,snliqv ,epore  ,qmelt  ,ponding, & !out
                         sag    ,fsa    ,fsr    ,fira   ,fsh    ,fgev   , & !out
                         trad   ,t2m    ,ssoil  ,lathea ,q2e    ,emissi, ch2b )   !out

#ifdef CCPP
    if (errflg /= 0) return
#endif

    sice = max(0.0, smc - sh2o)   
    sneqvo  = sneqv

    qvap = max( fgev/lathea, 0.)       ! positive part of fgev [mm/s] > 0
    qdew = abs( min(fgev/lathea, 0.))  ! negative part of fgev [mm/s] > 0
    edir = qvap - qdew

! compute water budgets (water storages, et components, and runoff)

     call water_glacier (nsnow  ,nsoil  ,imelt  ,dt     ,prcp   ,sfctmp , & !in
                         qvap   ,qdew   ,ficeold,zsoil  ,                 & !in
                         isnow  ,snowh  ,sneqv  ,snice  ,snliq  ,stc    , & !inout
                         dzsnso ,sh2o   ,sice   ,ponding,zsnso  ,         & !inout
                         runsrf ,runsub ,qsnow  ,ponding1       ,ponding2,qsnbot,fpice,esnow &  !out
                        )

!    if(maxval(sice) < 0.0001) then
!      write(message,*) "glacier has melted at:",iloc,jloc," are you sure this should be a glacier point?"
!      call wrf_debug(10,trim(message))
!    end if
     
! water and energy balance check

     call error_glacier (iloc   ,jloc   ,swdown ,fsa    ,fsr    ,fira   , &
                         fsh    ,fgev   ,ssoil  ,sag    ,prcp   ,edir   , &
#ifdef CCPP
                         runsrf ,runsub ,sneqv  ,dt     ,beg_wb, errmsg, errflg )
#else
                         runsrf ,runsub ,sneqv  ,dt     ,beg_wb )
#endif

#ifdef CCPP
     if (errflg /= 0) return
#endif

    if(snowh <= 1.e-6 .or. sneqv <= 1.e-3) then
     snowh = 0.0
     sneqv = 0.0
    end if

    if(swdown.ne.0.) then
      albedo = fsr / swdown
    else
      albedo = -999.9
    end if
    

  end subroutine noahmp_glacier
! ==================================================================================================
!>\ingroup NoahMP_LSM
  subroutine atm_glacier (sfcprs ,sfctmp ,q2     ,soldn  ,cosz   ,thair  , &
                          qair   ,eair   ,rhoair ,solad  ,solai  , &
                          swdown )     
! --------------------------------------------------------------------------------------------------
! re-process atmospheric forcing
! --------------------------------------------------------------------------------------------------
  implicit none
! --------------------------------------------------------------------------------------------------
! inputs

  real                          , intent(in)  :: sfcprs !pressure (pa)
  real                          , intent(in)  :: sfctmp !surface air temperature [k]
  real                          , intent(in)  :: q2     !mixing ratio (kg/kg)
  real                          , intent(in)  :: soldn  !downward shortwave radiation (w/m2)
  real                          , intent(in)  :: cosz   !cosine solar zenith angle [0-1]

! outputs

  real                          , intent(out) :: thair  !potential temperature (k)
  real                          , intent(out) :: qair   !specific humidity (kg/kg) (q2/(1+q2))
  real                          , intent(out) :: eair   !vapor pressure air (pa)
  real, dimension(       1:   2), intent(out) :: solad  !incoming direct solar radiation (w/m2)
  real, dimension(       1:   2), intent(out) :: solai  !incoming diffuse solar radiation (w/m2)
  real                          , intent(out) :: rhoair !density air (kg/m3)
  real                          , intent(out) :: swdown !downward solar filtered by sun angle [w/m2]

!locals

  real                                        :: pair   !atm bottom level pressure (pa)
! --------------------------------------------------------------------------------------------------

       pair   = sfcprs                   ! atm bottom level pressure (pa)
       thair  = sfctmp * (sfcprs/pair)**(rair/cpair) 
!       qair   = q2 / (1.0+q2)           ! mixing ratio to specific humidity [kg/kg]
       qair   = q2                       ! in wrf, driver converts to specific humidity

       eair   = qair*sfcprs / (0.622+0.378*qair)
       rhoair = (sfcprs-0.378*eair) / (rair*sfctmp)

       if(cosz <= 0.) then 
          swdown = 0.
       else
          swdown = soldn
       end if 

       solad(1) = swdown*0.7*0.5     ! direct  vis
       solad(2) = swdown*0.7*0.5     ! direct  nir
       solai(1) = swdown*0.3*0.5     ! diffuse vis
       solai(2) = swdown*0.3*0.5     ! diffuse nir

  end subroutine atm_glacier
! ==================================================================================================
! --------------------------------------------------------------------------------------------------
!>\ingroup NoahMP_LSM
  subroutine energy_glacier (nsnow  ,nsoil  ,isnow  ,dt     ,qsnow  ,rhoair , & !in
                             eair   ,sfcprs ,qair   ,sfctmp ,lwdn   ,uu     , & !in
                             vv     ,solad  ,solai  ,cosz   ,zref   ,         & !in
                             tbot   ,zbot   ,zsnso  ,dzsnso ,                 & !in
                             tg     ,stc    ,snowh  ,sneqv  ,sneqvo ,sh2o   , & !inout
                             smc    ,snice  ,snliq  ,albold ,cm     ,ch     , & !inout
#ifdef CCPP
                             tauss  ,qsfc   ,errmsg, errflg,                  & !inout
#else
                             tauss  ,qsfc   ,                                 & !inout
#endif
                             imelt  ,snicev ,snliqv ,epore  ,qmelt  ,ponding, & !out
                             sag    ,fsa    ,fsr    ,fira   ,fsh    ,fgev   , & !out
                             trad   ,t2m    ,ssoil  ,lathea ,q2e    ,emissi, ch2b )   !out

! --------------------------------------------------------------------------------------------------
! --------------------------------------------------------------------------------------------------
!  use noahmp_veg_parameters
!  use noahmp_rad_parameters
! --------------------------------------------------------------------------------------------------
  implicit none
! --------------------------------------------------------------------------------------------------
! inputs
  integer                           , intent(in)    :: nsnow  !maximum no. of snow layers        
  integer                           , intent(in)    :: nsoil  !number of soil layers
  integer                           , intent(in)    :: isnow  !actual no. of snow layers
  real                              , intent(in)    :: dt     !time step [sec]
  real                              , intent(in)    :: qsnow  !snowfall on the ground (mm/s)
  real                              , intent(in)    :: rhoair !density air (kg/m3)
  real                              , intent(in)    :: eair   !vapor pressure air (pa)
  real                              , intent(in)    :: sfcprs !pressure (pa)
  real                              , intent(in)    :: qair   !specific humidity (kg/kg)
  real                              , intent(in)    :: sfctmp !air temperature (k)
  real                              , intent(in)    :: lwdn   !downward longwave radiation (w/m2)
  real                              , intent(in)    :: uu     !wind speed in e-w dir (m/s)
  real                              , intent(in)    :: vv     !wind speed in n-s dir (m/s)
  real   , dimension(       1:    2), intent(in)    :: solad  !incoming direct solar rad. (w/m2)
  real   , dimension(       1:    2), intent(in)    :: solai  !incoming diffuse solar rad. (w/m2)
  real                              , intent(in)    :: cosz   !cosine solar zenith angle (0-1)
  real                              , intent(in)    :: zref   !reference height (m)
  real                              , intent(in)    :: tbot   !bottom condition for soil temp. (k) 
  real                              , intent(in)    :: zbot   !depth for tbot [m]
  real   , dimension(-nsnow+1:nsoil), intent(in)    :: zsnso  !layer-bottom depth from snow surf [m]
  real   , dimension(-nsnow+1:nsoil), intent(in)    :: dzsnso !depth of snow & soil layer-bottom [m]

! input & output
  real                              , intent(inout) :: tg     !ground temperature (k)
  real   , dimension(-nsnow+1:nsoil), intent(inout) :: stc    !snow/soil temperature [k]
  real                              , intent(inout) :: snowh  !snow height [m]
  real                              , intent(inout) :: sneqv  !snow mass (mm)
  real                              , intent(inout) :: sneqvo !snow mass at last time step (mm)
  real   , dimension(       1:nsoil), intent(inout) :: sh2o   !liquid soil moisture [m3/m3]
  real   , dimension(       1:nsoil), intent(inout) :: smc    !soil moisture (ice + liq.) [m3/m3]
  real   , dimension(-nsnow+1:    0), intent(inout) :: snice  !snow ice mass (kg/m2)
  real   , dimension(-nsnow+1:    0), intent(inout) :: snliq  !snow liq mass (kg/m2)
  real                              , intent(inout) :: albold !snow albedo at last time step(class type)
  real                              , intent(inout) :: cm     !momentum drag coefficient
  real                              , intent(inout) :: ch     !sensible heat exchange coefficient
  real                              , intent(inout) :: tauss  !snow aging factor
  real                              , intent(inout) :: qsfc   !mixing ratio at lowest model layer
  
#ifdef CCPP  
  character(len=*)                  , intent(inout) :: errmsg
  integer                           , intent(inout) :: errflg
#endif

! outputs
  integer, dimension(-nsnow+1:nsoil), intent(out)   :: imelt  !phase change index [1-melt; 2-freeze]
  real   , dimension(-nsnow+1:    0), intent(out)   :: snicev !partial volume ice [m3/m3]
  real   , dimension(-nsnow+1:    0), intent(out)   :: snliqv !partial volume liq. water [m3/m3]
  real   , dimension(-nsnow+1:    0), intent(out)   :: epore  !effective porosity [m3/m3]
  real                              , intent(out)   :: qmelt  !snowmelt [mm/s]
  real                              , intent(out)   :: ponding!pounding at ground [mm]
  real                              , intent(out)   :: sag    !solar rad. absorbed by ground (w/m2)
  real                              , intent(out)   :: fsa    !tot. absorbed solar radiation (w/m2)
  real                              , intent(out)   :: fsr    !tot. reflected solar radiation (w/m2)
  real                              , intent(out)   :: fira   !total net lw. rad (w/m2)   [+ to atm]
  real                              , intent(out)   :: fsh    !total sensible heat (w/m2) [+ to atm]
  real                              , intent(out)   :: fgev   !ground evaporation (w/m2)  [+ to atm]
  real                              , intent(out)   :: trad   !radiative temperature (k)
  real                              , intent(out)   :: t2m    !2 m height air temperature (k)
  real                              , intent(out)   :: ssoil  !ground heat flux (w/m2)   [+ to soil]
  real                              , intent(out)   :: lathea !latent heat vap./sublimation (j/kg)
  real                              , intent(out)   :: q2e
  real                              , intent(out)   :: emissi
  real                              , intent(out)   :: ch2b   !sensible heat conductance, canopy air to zlvl air (m/s)


! local
  real                                              :: ur     !wind speed at height zlvl (m/s)
  real                                              :: zlvl   !reference height (m)
  real                                              :: rsurf  !ground surface resistance (s/m)
  real                                              :: zpd    !zero plane displacement (m)
  real                                              :: z0mg   !z0 momentum, ground (m)
  real                                              :: emg    !ground emissivity
  real                                              :: fire   !emitted ir (w/m2)
  real, dimension(-nsnow+1:nsoil)                   :: fact   !temporary used in phase change
  real, dimension(-nsnow+1:nsoil)                   :: df     !thermal conductivity [w/m/k]
  real, dimension(-nsnow+1:nsoil)                   :: hcpct  !heat capacity [j/m3/k]
  real                                              :: gamma  !psychrometric constant (pa/k)
  real                                              :: rhsur  !raltive humidity in surface soil/snow air space (-)

! ---------------------------------------------------------------------------------------------------

! wind speed at reference height: ur >= 1

    ur = max( sqrt(uu**2.+vv**2.), 1. )

! roughness length and displacement height

     z0mg = z0sno
     zpd  = snowh

     zlvl = zpd + zref

! thermal properties of soil, snow, lake, and frozen soil

  call thermoprop_glacier (nsoil   ,nsnow   ,isnow   ,dzsnso  ,          & !in
                           dt      ,snowh   ,snice   ,snliq   ,          & !in
                           df      ,hcpct   ,snicev  ,snliqv  ,epore   , & !out
                           fact    )                                       !out

! solar radiation: absorbed & reflected by the ground

  call  radiation_glacier (dt      ,tg      ,sneqvo  ,sneqv   ,cosz    , & !in
                           qsnow   ,solad   ,solai   ,                   & !in
                           albold  ,tauss   ,                            & !inout
                           sag     ,fsr     ,fsa)                          !out

! vegetation and ground emissivity

     emg = 0.98

! soil surface resistance for ground evap.

     rhsur = 1.0
     rsurf = 1.0

! set psychrometric constant

     lathea = hsub
     gamma = cpair*sfcprs/(0.622*lathea)

! surface temperatures of the ground and energy fluxes

    call glacier_flux (nsoil   ,nsnow   ,emg     ,isnow   ,df      ,dzsnso  ,z0mg    , & !in
                       zlvl    ,zpd     ,qair    ,sfctmp  ,rhoair  ,sfcprs  , & !in
                       ur      ,gamma   ,rsurf   ,lwdn    ,rhsur   ,smc     , & !in
                       eair    ,stc     ,sag     ,snowh   ,lathea  ,sh2o    , & !in
#ifdef CCPP
                       cm      ,ch      ,tg      ,qsfc    ,errmsg  ,errflg  , & !inout
#else
                       cm      ,ch      ,tg      ,qsfc    ,          & !inout
#endif
                       fira    ,fsh     ,fgev    ,ssoil   ,          & !out
                       t2m     ,q2e     ,ch2b)                         !out 

!energy balance at surface: sag=(irb+shb+evb+ghb)

    fire = lwdn + fira

    if(fire <=0.) then
#ifdef CCPP
      errflg = 1
      errmsg = "stop in noah-mp: emitted longwave <0"
      return 
#else
      call wrf_error_fatal("stop in noah-mp: emitted longwave <0")
#endif
    end if

    ! compute a net emissivity
    emissi = emg

    ! when we're computing a trad, subtract from the emitted ir the
    ! reflected portion of the incoming lwdn, so we're just
    ! considering the ir originating in the canopy/ground system.
    
    trad = ( ( fire - (1-emissi)*lwdn ) / (emissi*sb) ) ** 0.25

! 3l snow & 4l soil temperatures

    call tsnosoi_glacier (nsoil   ,nsnow   ,isnow   ,dt      ,tbot    , & !in
                          ssoil   ,snowh   ,zbot    ,zsnso   ,df      , & !in
		          hcpct   ,                                     & !in
                          stc     )                                       !inout

! adjusting snow surface temperature
     if(opt_stc == 2) then
      if (snowh > 0.05 .and. tg > tfrz) tg = tfrz
     end if

! energy released or consumed by snow & frozen soil

 call phasechange_glacier (nsnow   ,nsoil   ,isnow   ,dt      ,fact    , & !in
                           dzsnso  ,                                     & !in
                           stc     ,snice   ,snliq   ,sneqv   ,snowh   , & !inout
                           smc     ,sh2o    ,                            & !inout
                           qmelt   ,imelt   ,ponding )                     !out


  end subroutine energy_glacier
! ==================================================================================================
!>\ingroup NoahMP_LSM
  subroutine thermoprop_glacier (nsoil   ,nsnow   ,isnow   ,dzsnso  , & !in
                                 dt      ,snowh   ,snice   ,snliq   , & !in
                                 df      ,hcpct   ,snicev  ,snliqv  ,epore   , & !out
                                 fact    )                                       !out
! ------------------------------------------------------------------------------------------------- 
! -------------------------------------------------------------------------------------------------
  implicit none
! --------------------------------------------------------------------------------------------------
! inputs
  integer                        , intent(in)  :: nsoil   !number of soil layers
  integer                        , intent(in)  :: nsnow   !maximum no. of snow layers        
  integer                        , intent(in)  :: isnow   !actual no. of snow layers
  real                           , intent(in)  :: dt      !time step [s]
  real, dimension(-nsnow+1:    0), intent(in)  :: snice   !snow ice mass (kg/m2)
  real, dimension(-nsnow+1:    0), intent(in)  :: snliq   !snow liq mass (kg/m2)
  real, dimension(-nsnow+1:nsoil), intent(in)  :: dzsnso  !thickness of snow/soil layers [m]
  real                           , intent(in)  :: snowh   !snow height [m]

! outputs
  real, dimension(-nsnow+1:nsoil), intent(out) :: df      !thermal conductivity [w/m/k]
  real, dimension(-nsnow+1:nsoil), intent(out) :: hcpct   !heat capacity [j/m3/k]
  real, dimension(-nsnow+1:    0), intent(out) :: snicev  !partial volume of ice [m3/m3]
  real, dimension(-nsnow+1:    0), intent(out) :: snliqv  !partial volume of liquid water [m3/m3]
  real, dimension(-nsnow+1:    0), intent(out) :: epore   !effective porosity [m3/m3]
  real, dimension(-nsnow+1:nsoil), intent(out) :: fact    !computing energy for phase change
! --------------------------------------------------------------------------------------------------
! locals

  integer :: iz, iz2
  real, dimension(-nsnow+1:    0)              :: cvsno   !volumetric specific heat (j/m3/k)
  real, dimension(-nsnow+1:    0)              :: tksno   !snow thermal conductivity (j/m3/k)
  real                                         :: zmid    !mid-point soil depth
! --------------------------------------------------------------------------------------------------

! compute snow thermal conductivity and heat capacity

    call csnow_glacier (isnow   ,nsnow   ,nsoil   ,snice   ,snliq   ,dzsnso  , & !in
                        tksno   ,cvsno   ,snicev  ,snliqv  ,epore   )   !out

    do iz = isnow+1, 0
      df   (iz) = tksno(iz)
      hcpct(iz) = cvsno(iz)
    end do

! compute soil thermal properties (using noah glacial ice approximations)

    do  iz = 1, nsoil
       zmid      = 0.5 * (dzsnso(iz))
       do iz2 = 1, iz-1
         zmid = zmid + dzsnso(iz2)
       end do
       hcpct(iz) = 1.e6 * ( 0.8194 + 0.1309*zmid )
       df(iz)    = 0.32333 + ( 0.10073 * zmid )
    end do
       
! combine a temporary variable used for melting/freezing of snow and frozen soil

    do iz = isnow+1,nsoil
     fact(iz) = dt/(hcpct(iz)*dzsnso(iz))
    end do

! snow/soil interface

    if(isnow == 0) then
       df(1) = (df(1)*dzsnso(1)+0.35*snowh)      / (snowh    +dzsnso(1)) 
    else
       df(1) = (df(1)*dzsnso(1)+df(0)*dzsnso(0)) / (dzsnso(0)+dzsnso(1))
    end if


  end subroutine thermoprop_glacier
! ==================================================================================================
! --------------------------------------------------------------------------------------------------
!>\ingroup NoahMP_LSM  
  subroutine csnow_glacier (isnow   ,nsnow   ,nsoil   ,snice   ,snliq   ,dzsnso  , & !in
                            tksno   ,cvsno   ,snicev  ,snliqv  ,epore   )   !out
! --------------------------------------------------------------------------------------------------
! snow bulk density,volumetric capacity, and thermal conductivity
!---------------------------------------------------------------------------------------------------
  implicit none
!---------------------------------------------------------------------------------------------------
! inputs

  integer,                          intent(in) :: isnow  !number of snow layers (-)            
  integer                        ,  intent(in) :: nsnow  !maximum no. of snow layers        
  integer                        ,  intent(in) :: nsoil  !number of soil layers
  real, dimension(-nsnow+1:    0),  intent(in) :: snice  !snow ice mass (kg/m2)
  real, dimension(-nsnow+1:    0),  intent(in) :: snliq  !snow liq mass (kg/m2) 
  real, dimension(-nsnow+1:nsoil),  intent(in) :: dzsnso !snow/soil layer thickness [m]

! outputs

  real, dimension(-nsnow+1:    0), intent(out) :: cvsno  !volumetric specific heat (j/m3/k)
  real, dimension(-nsnow+1:    0), intent(out) :: tksno  !thermal conductivity (w/m/k)
  real, dimension(-nsnow+1:    0), intent(out) :: snicev !partial volume of ice [m3/m3]
  real, dimension(-nsnow+1:    0), intent(out) :: snliqv !partial volume of liquid water [m3/m3]
  real, dimension(-nsnow+1:    0), intent(out) :: epore  !effective porosity [m3/m3]

! locals

  integer :: iz
  real, dimension(-nsnow+1:    0) :: bdsnoi  !bulk density of snow(kg/m3)

!---------------------------------------------------------------------------------------------------
! thermal capacity of snow

  do iz = isnow+1, 0
      snicev(iz)   = min(1., snice(iz)/(dzsnso(iz)*denice) )
      epore(iz)    = 1. - snicev(iz)
      snliqv(iz)   = min(epore(iz),snliq(iz)/(dzsnso(iz)*denh2o))
  enddo

  do iz = isnow+1, 0
      bdsnoi(iz) = (snice(iz)+snliq(iz))/dzsnso(iz)
      cvsno(iz) = cice*snicev(iz)+cwat*snliqv(iz)
!      cvsno(iz) = 0.525e06                          ! constant
  enddo

! thermal conductivity of snow

  do iz = isnow+1, 0
     tksno(iz) = 3.2217e-6*bdsnoi(iz)**2.           ! stieglitz(yen,1965)
!    tksno(iz) = 2e-2+2.5e-6*bdsnoi(iz)*bdsnoi(iz)   ! anderson, 1976
!    tksno(iz) = 0.35                                ! constant
!    tksno(iz) = 2.576e-6*bdsnoi(iz)**2. + 0.074    ! verseghy (1991)
!    tksno(iz) = 2.22*(bdsnoi(iz)/1000.)**1.88      ! douvill(yen, 1981)
  enddo

  end subroutine csnow_glacier
!===================================================================================================
!>\ingroup NoahMP_LSM
  subroutine radiation_glacier (dt      ,tg      ,sneqvo  ,sneqv   ,cosz    , & !in
                                qsnow   ,solad   ,solai   ,                   & !in
                                albold  ,tauss   ,                            & !inout
                                sag     ,fsr     ,fsa)                          !out
! --------------------------------------------------------------------------------------------------
  implicit none
! --------------------------------------------------------------------------------------------------
! input
  real, intent(in)                     :: dt     !time step [s]
  real, intent(in)                     :: tg     !ground temperature (k)
  real, intent(in)                     :: sneqvo !snow mass at last time step(mm)
  real, intent(in)                     :: sneqv  !snow mass (mm)
  real, intent(in)                     :: cosz   !cosine solar zenith angle (0-1)
  real, intent(in)                     :: qsnow  !snowfall (mm/s)
  real, dimension(1:2)    , intent(in) :: solad  !incoming direct solar radiation (w/m2)
  real, dimension(1:2)    , intent(in) :: solai  !incoming diffuse solar radiation (w/m2)

! inout
  real,                  intent(inout) :: albold !snow albedo at last time step (class type)
  real,                  intent(inout) :: tauss  !non-dimensional snow age

! output
  real, intent(out)                    :: sag    !solar radiation absorbed by ground (w/m2)
  real, intent(out)                    :: fsr    !total reflected solar radiation (w/m2)
  real, intent(out)                    :: fsa    !total absorbed solar radiation (w/m2)

! local
  integer                              :: ib     !number of radiation bands
  integer                              :: nband  !number of radiation bands
  real                                 :: fage   !snow age function (0 - new snow)
  real, dimension(1:2)                 :: albsnd !snow albedo (direct)
  real, dimension(1:2)                 :: albsni !snow albedo (diffuse)
  real                                 :: alb    !current class albedo
  real                                 :: abs    !temporary absorbed rad
  real                                 :: ref    !temporary reflected rad
  real                                 :: fsno   !snow-cover fraction, = 1 if any snow
  real, dimension(1:2)                 :: albice !albedo land ice: 1=vis, 2=nir

  real,parameter :: mpe = 1.e-6

! --------------------------------------------------------------------------------------------------

  nband = 2
  albsnd = 0.0
  albsni = 0.0
  albice(1) = 0.80    !albedo land ice: 1=vis, 2=nir
  albice(2) = 0.55

! snow age

  call snow_age_glacier (dt,tg,sneqvo,sneqv,tauss,fage)

! snow albedos: age even when sun is not present

  if(opt_alb == 1) &
     call snowalb_bats_glacier (nband,cosz,fage,albsnd,albsni)
  if(opt_alb == 2) then
     call snowalb_class_glacier(nband,qsnow,dt,alb,albold,albsnd,albsni)
     albold = alb
  end if

! zero summed solar fluxes

   sag = 0.
   fsa = 0.
   fsr = 0.
   
   fsno = 0.0
   if(sneqv > 0.0) fsno = 1.0

! loop over nband wavebands

  do ib = 1, nband

    albsnd(ib) = albice(ib)*(1.-fsno) + albsnd(ib)*fsno
    albsni(ib) = albice(ib)*(1.-fsno) + albsni(ib)*fsno

! solar radiation absorbed by ground surface

    abs = solad(ib)*(1.-albsnd(ib)) + solai(ib)*(1.-albsni(ib))
    sag = sag + abs
    fsa = fsa + abs
    
    ref = solad(ib)*albsnd(ib) + solai(ib)*albsni(ib)
    fsr = fsr + ref
    
  end do

  end subroutine radiation_glacier
! ==================================================================================================
!>\ingroup NoahMP_LSM
  subroutine snow_age_glacier (dt,tg,sneqvo,sneqv,tauss,fage)
! --------------------------------------------------------------------------------------------------
  implicit none
! ------------------------ code history ------------------------------------------------------------
! from bats
! ------------------------ input/output variables --------------------------------------------------
!input
   real, intent(in) :: dt        !main time step (s)
   real, intent(in) :: tg        !ground temperature (k)
   real, intent(in) :: sneqvo    !snow mass at last time step(mm)
   real, intent(in) :: sneqv     !snow water per unit ground area (mm)

! inout
  real,  intent(inout) :: tauss  !non-dimensional snow age

!output
   real, intent(out) :: fage     !snow age

!local
   real            :: tage       !total aging effects
   real            :: age1       !effects of grain growth due to vapor diffusion
   real            :: age2       !effects of grain growth at freezing of melt water
   real            :: age3       !effects of soot
   real            :: dela       !temporary variable
   real            :: sge        !temporary variable
   real            :: dels       !temporary variable
   real            :: dela0      !temporary variable
   real            :: arg        !temporary variable
! see yang et al. (1997) j.of climate for detail.
!---------------------------------------------------------------------------------------------------

   if(sneqv.le.0.0) then
          tauss = 0.
   else if (sneqv.gt.800.) then
          tauss = 0.
   else
!          tauss = 0.
          dela0 = 1.e-6*dt
          arg   = 5.e3*(1./tfrz-1./tg)
          age1  = exp(arg)
          age2  = exp(amin1(0.,10.*arg))
          age3  = 0.3
          tage  = age1+age2+age3
          dela  = dela0*tage
          dels  = amax1(0.0,sneqv-sneqvo) / swemx
          sge   = (tauss+dela)*(1.0-dels)
          tauss = amax1(0.,sge)
   endif

   fage= tauss/(tauss+1.)

  end subroutine snow_age_glacier
! ==================================================================================================
! --------------------------------------------------------------------------------------------------
!>\ingroup NoahMP_LSM
  subroutine snowalb_bats_glacier (nband,cosz,fage,albsnd,albsni)
! --------------------------------------------------------------------------------------------------
  implicit none
! --------------------------------------------------------------------------------------------------
! input

  integer,intent(in) :: nband  !number of waveband classes

  real,intent(in) :: cosz    !cosine solar zenith angle
  real,intent(in) :: fage    !snow age correction

! output

  real, dimension(1:2),intent(out) :: albsnd !snow albedo for direct(1=vis, 2=nir)
  real, dimension(1:2),intent(out) :: albsni !snow albedo for diffuse
! ---------------------------------------------------------------------------------------------

  real :: fzen                 !zenith angle correction
  real :: cf1                  !temperary variable
  real :: sl2                  !2.*sl
  real :: sl1                  !1/sl
  real :: sl                   !adjustable parameter
  real, parameter :: c1 = 0.2  !default in bats 
  real, parameter :: c2 = 0.5  !default in bats
!  real, parameter :: c1 = 0.2 * 2. ! double the default to match sleepers river's
!  real, parameter :: c2 = 0.5 * 2. ! snow surface albedo (double aging effects)
! ---------------------------------------------------------------------------------------------
! zero albedos for all points

        albsnd(1: nband) = 0.
        albsni(1: nband) = 0.

! when cosz > 0

        sl=2.0
        sl1=1./sl
        sl2=2.*sl
        cf1=((1.+sl1)/(1.+sl2*cosz)-sl1)
        fzen=amax1(cf1,0.)

        albsni(1)=0.95*(1.-c1*fage)         
        albsni(2)=0.65*(1.-c2*fage)        

        albsnd(1)=albsni(1)+0.4*fzen*(1.-albsni(1))    !  vis direct
        albsnd(2)=albsni(2)+0.4*fzen*(1.-albsni(2))    !  nir direct

  end subroutine snowalb_bats_glacier
! ==================================================================================================
! --------------------------------------------------------------------------------------------------
!>\ingroup NoahMP_LSM
  subroutine snowalb_class_glacier (nband,qsnow,dt,alb,albold,albsnd,albsni)
! --------------------------------------------------------------------------------------------------
  implicit none
! --------------------------------------------------------------------------------------------------
! input

  integer,intent(in) :: nband  !number of waveband classes

  real,intent(in) :: qsnow     !snowfall (mm/s)
  real,intent(in) :: dt        !time step (sec)
  real,intent(in) :: albold    !snow albedo at last time step

! in & out

  real,                intent(inout) :: alb        ! 
! output

  real, dimension(1:2),intent(out) :: albsnd !snow albedo for direct(1=vis, 2=nir)
  real, dimension(1:2),intent(out) :: albsni !snow albedo for diffuse
! ---------------------------------------------------------------------------------------------

! ---------------------------------------------------------------------------------------------
! zero albedos for all points

        albsnd(1: nband) = 0.
        albsni(1: nband) = 0.

! when cosz > 0

         alb = 0.55 + (albold-0.55) * exp(-0.01*dt/3600.)

! 1 mm fresh snow(swe) -- 10mm snow depth, assumed the fresh snow density 100kg/m3
! here assume 1cm snow depth will fully cover the old snow

         if (qsnow > 0.) then
           alb = alb + min(qsnow*dt,swemx) * (0.84-alb)/(swemx)
         endif

         albsni(1)= alb         ! vis diffuse
         albsni(2)= alb         ! nir diffuse
         albsnd(1)= alb         ! vis direct
         albsnd(2)= alb         ! nir direct

  end subroutine snowalb_class_glacier
! ==================================================================================================
!>\ingroup NoahMP_LSM
  subroutine glacier_flux (nsoil   ,nsnow   ,emg     ,isnow   ,df      ,dzsnso  ,z0m     , & !in
                           zlvl    ,zpd     ,qair    ,sfctmp  ,rhoair  ,sfcprs  , & !in
                           ur      ,gamma   ,rsurf   ,lwdn    ,rhsur   ,smc     , & !in
                           eair    ,stc     ,sag     ,snowh   ,lathea  ,sh2o    , & !in
#ifdef CCPP
                           cm      ,ch      ,tgb     ,qsfc    ,errmsg  ,errflg  , & !inout
#else
                           cm      ,ch      ,tgb     ,qsfc    ,          & !inout
#endif
                           irb     ,shb     ,evb     ,ghb     ,          & !out
                           t2mb    ,q2b     ,ehb2)                         !out 

! --------------------------------------------------------------------------------------------------
! use newton-raphson iteration to solve ground (tg) temperature
! that balances the surface energy budgets for glacier.

! bare soil:
! -sab + irb[tg] + shb[tg] + evb[tg] + ghb[tg] = 0
! ----------------------------------------------------------------------
!  use module_model_constants
! ----------------------------------------------------------------------
  implicit none
! ----------------------------------------------------------------------
! input
  integer, intent(in)                         :: nsnow  !maximum no. of snow layers        
  integer, intent(in)                         :: nsoil  !number of soil layers
  real,                            intent(in) :: emg    !ground emissivity
  integer,                         intent(in) :: isnow  !actual no. of snow layers
  real, dimension(-nsnow+1:nsoil), intent(in) :: df     !thermal conductivity of snow/soil (w/m/k)
  real, dimension(-nsnow+1:nsoil), intent(in) :: dzsnso !thickness of snow/soil layers (m)
  real,                            intent(in) :: z0m    !roughness length, momentum, ground (m)
  real,                            intent(in) :: zlvl   !reference height (m)
  real,                            intent(in) :: zpd    !zero plane displacement (m)
  real,                            intent(in) :: qair   !specific humidity at height zlvl (kg/kg)
  real,                            intent(in) :: sfctmp !air temperature at reference height (k)
  real,                            intent(in) :: rhoair !density air (kg/m3)
  real,                            intent(in) :: sfcprs !density air (kg/m3)
  real,                            intent(in) :: ur     !wind speed at height zlvl (m/s)
  real,                            intent(in) :: gamma  !psychrometric constant (pa/k)
  real,                            intent(in) :: rsurf  !ground surface resistance (s/m)
  real,                            intent(in) :: lwdn   !atmospheric longwave radiation (w/m2)
  real,                            intent(in) :: rhsur  !raltive humidity in surface soil/snow air space (-)
  real,                            intent(in) :: eair   !vapor pressure air at height (pa)
  real, dimension(-nsnow+1:nsoil), intent(in) :: stc    !soil/snow temperature (k)
  real, dimension(       1:nsoil), intent(in) :: smc    !soil moisture
  real, dimension(       1:nsoil), intent(in) :: sh2o   !soil liquid water
  real,                            intent(in) :: sag    !solar radiation absorbed by ground (w/m2)
  real,                            intent(in) :: snowh  !actual snow depth [m]
  real,                            intent(in) :: lathea !latent heat of vaporization/subli (j/kg)

! input/output
  real,                         intent(inout) :: cm     !momentum drag coefficient
  real,                         intent(inout) :: ch     !sensible heat exchange coefficient
  real,                         intent(inout) :: tgb    !ground temperature (k)
  real,                         intent(inout) :: qsfc   !mixing ratio at lowest model layer
  
#ifdef CCPP  
  character(len=*),             intent(inout) :: errmsg
  integer,                      intent(inout) :: errflg
#endif
  
! output
! -sab + irb[tg] + shb[tg] + evb[tg] + ghb[tg] = 0
  real,                           intent(out) :: irb    !net longwave rad (w/m2)   [+ to atm]
  real,                           intent(out) :: shb    !sensible heat flux (w/m2) [+ to atm]
  real,                           intent(out) :: evb    !latent heat flux (w/m2)   [+ to atm]
  real,                           intent(out) :: ghb    !ground heat flux (w/m2)  [+ to soil]
  real,                           intent(out) :: t2mb   !2 m height air temperature (k)
  real,                           intent(out) :: q2b    !bare ground heat conductance
  real,                           intent(out) :: ehb2   !sensible heat conductance for diagnostics


! local variables 
  integer :: niterb  !number of iterations for surface temperature
  real    :: mpe     !prevents overflow error if division by zero
  real    :: dtg        !change in tg, last iteration (k)
  integer :: mozsgn  !number of times moz changes sign
  real    :: mozold     !monin-obukhov stability parameter from prior iteration
  real    :: fm2          !monin-obukhov momentum adjustment at 2m
  real    :: fh2          !monin-obukhov heat adjustment at 2m
  real    :: ch2          !surface exchange at 2m
  real    :: h          !temporary sensible heat flux (w/m2)
  real    :: fv         !friction velocity (m/s)
  real    :: cir        !coefficients for ir as function of ts**4
  real    :: cgh        !coefficients for st as function of ts
  real    :: csh        !coefficients for sh as function of ts
  real    :: cev        !coefficients for ev as function of esat[ts]
  real    :: cq2b       !
  integer :: iter    !iteration index
  real    :: z0h        !roughness length, sensible heat, ground (m)
  real    :: moz        !monin-obukhov stability parameter
  real    :: fm         !momentum stability correction, weighted by prior iters
  real    :: fh         !sen heat stability correction, weighted by prior iters
  real    :: ramb       !aerodynamic resistance for momentum (s/m)
  real    :: rahb       !aerodynamic resistance for sensible heat (s/m)
  real    :: rawb       !aerodynamic resistance for water vapor (s/m)
  real    :: estg       !saturation vapor pressure at tg (pa)
  real    :: destg      !d(es)/dt at tg (pa/k)
  real    :: esatw      !es for water
  real    :: esati      !es for ice
  real    :: dsatw      !d(es)/dt at tg (pa/k) for water
  real    :: dsati      !d(es)/dt at tg (pa/k) for ice
  real    :: a          !temporary calculation
  real    :: b          !temporary calculation
  real    :: t, tdc     !kelvin to degree celsius with limit -50 to +50
  real, dimension(       1:nsoil) :: sice   !soil ice

  tdc(t)   = min( 50., max(-50.,(t-tfrz)) )

! -----------------------------------------------------------------
! initialization variables that do not depend on stability iteration
! -----------------------------------------------------------------
        niterb = 5
        mpe    = 1e-6
        dtg    = 0.
        mozsgn = 0
        mozold = 0.
        h      = 0.
        fv     = 0.1

        cir = emg*sb
        cgh = 2.*df(isnow+1)/dzsnso(isnow+1)

! -----------------------------------------------------------------
      loop3: do iter = 1, niterb  ! begin stability iteration

        z0h = z0m 

!       for now, only allow sfcdif1 until others can be fixed

        call sfcdif1_glacier(iter   ,zlvl   ,zpd    ,z0h    ,z0m    , & !in
                     qair   ,sfctmp ,h      ,rhoair ,mpe    ,ur     , & !in
#ifdef CCPP
       &             moz ,mozsgn ,fm ,fh ,fm2 ,fh2   ,errmsg, errflg, & !inout
#else 
       &             moz ,mozsgn ,fm ,fh ,fm2 ,fh2                  , & !inout
#endif
       &             fv     ,cm     ,ch     ,ch2)                       !out

#ifdef CCPP
        if (errflg /= 0) return
#endif
        ramb = max(1.,1./(cm*ur))
        rahb = max(1.,1./(ch*ur))
        rawb = rahb

! es and d(es)/dt evaluated at tg

        t = tdc(tgb)
        call esat(t, esatw, esati, dsatw, dsati)
        if (t .gt. 0.) then
            estg  = esatw
            destg = dsatw
        else
            estg  = esati
            destg = dsati
        end if

        csh = rhoair*cpair/rahb
        cev = rhoair*cpair/gamma/(rsurf+rawb)

! surface fluxes and dtg

        irb   = cir * tgb**4 - emg*lwdn
        shb   = csh * (tgb        - sfctmp      )
        evb   = cev * (estg*rhsur - eair        )
        ghb   = cgh * (tgb        - stc(isnow+1))

        b     = sag-irb-shb-evb-ghb
        a     = 4.*cir*tgb**3 + csh + cev*destg + cgh
        dtg   = b/a

        irb = irb + 4.*cir*tgb**3*dtg
        shb = shb + csh*dtg
        evb = evb + cev*destg*dtg
        ghb = ghb + cgh*dtg

! update ground surface temperature
        tgb = tgb + dtg

! for m-o length
        h = csh * (tgb - sfctmp)

        t = tdc(tgb)
        call esat(t, esatw, esati, dsatw, dsati)
        if (t .gt. 0.) then
            estg  = esatw
        else
            estg  = esati
        end if
        qsfc = 0.622*(estg*rhsur)/(sfcprs-0.378*(estg*rhsur))

     end do loop3 ! end stability iteration
! -----------------------------------------------------------------

! if snow on ground and tg > tfrz: reset tg = tfrz. reevaluate ground fluxes.

     sice = smc - sh2o
     if(opt_stc == 1) then
     if ((maxval(sice) > 0.0 .or. snowh > 0.0) .and. tgb > tfrz) then
          tgb = tfrz
          irb = cir * tgb**4 - emg*lwdn
          shb = csh * (tgb        - sfctmp)
          evb = cev * (estg*rhsur - eair )          !estg reevaluate ?
          ghb = sag - (irb+shb+evb)
     end if
     end if

! 2m air temperature
     ehb2  = fv*vkc/(log((2.+z0h)/z0h)-fh2)
     cq2b  = ehb2
     if (ehb2.lt.1.e-5 ) then
       t2mb  = tgb
       q2b   = qsfc
     else
       t2mb  = tgb - shb/(rhoair*cpair) * 1./ehb2
       q2b   = qsfc - evb/(lathea*rhoair)*(1./cq2b + rsurf)
     endif

! update ch 
     ch = 1./rahb

  end subroutine glacier_flux
!  ==================================================================================================
!>\ingroup NoahMP_LSM
  subroutine esat(t, esw, esi, desw, desi)
!---------------------------------------------------------------------------------------------------
! use polynomials to calculate saturation vapor pressure and derivative with
! respect to temperature: over water when t > 0 c and over ice when t <= 0 c
  implicit none
!---------------------------------------------------------------------------------------------------
! in

  real, intent(in)  :: t              !temperature

!out

  real, intent(out) :: esw            !saturation vapor pressure over water (pa)
  real, intent(out) :: esi            !saturation vapor pressure over ice (pa)
  real, intent(out) :: desw           !d(esat)/dt over water (pa/k)
  real, intent(out) :: desi           !d(esat)/dt over ice (pa/k)

! local

  real :: a0,a1,a2,a3,a4,a5,a6  !coefficients for esat over water
  real :: b0,b1,b2,b3,b4,b5,b6  !coefficients for esat over ice
  real :: c0,c1,c2,c3,c4,c5,c6  !coefficients for dsat over water
  real :: d0,d1,d2,d3,d4,d5,d6  !coefficients for dsat over ice

  parameter (a0=6.107799961    , a1=4.436518521e-01,  &
             a2=1.428945805e-02, a3=2.650648471e-04,  &
             a4=3.031240396e-06, a5=2.034080948e-08,  &
             a6=6.136820929e-11)

  parameter (b0=6.109177956    , b1=5.034698970e-01,  &
             b2=1.886013408e-02, b3=4.176223716e-04,  &
             b4=5.824720280e-06, b5=4.838803174e-08,  &
             b6=1.838826904e-10)

  parameter (c0= 4.438099984e-01, c1=2.857002636e-02,  &
             c2= 7.938054040e-04, c3=1.215215065e-05,  &
             c4= 1.036561403e-07, c5=3.532421810e-10,  &
             c6=-7.090244804e-13)

  parameter (d0=5.030305237e-01, d1=3.773255020e-02,  &
             d2=1.267995369e-03, d3=2.477563108e-05,  &
             d4=3.005693132e-07, d5=2.158542548e-09,  &
             d6=7.131097725e-12)

  esw  = 100.*(a0+t*(a1+t*(a2+t*(a3+t*(a4+t*(a5+t*a6))))))
  esi  = 100.*(b0+t*(b1+t*(b2+t*(b3+t*(b4+t*(b5+t*b6))))))
  desw = 100.*(c0+t*(c1+t*(c2+t*(c3+t*(c4+t*(c5+t*c6))))))
  desi = 100.*(d0+t*(d1+t*(d2+t*(d3+t*(d4+t*(d5+t*d6))))))

  end subroutine esat
! ==================================================================================================
!>\ingroup NoahMP_LSM
  subroutine sfcdif1_glacier(iter   ,zlvl   ,zpd    ,z0h    ,z0m    , & !in
                     qair   ,sfctmp ,h      ,rhoair ,mpe    ,ur     , & !in
#ifdef CCPP
       &             moz ,mozsgn ,fm ,fh ,fm2 ,fh2 ,errmsg ,errflg  , & !inout
#else
       &             moz ,mozsgn ,fm ,fh ,fm2 ,fh2                  , & !inout
#endif
       &             fv     ,cm     ,ch     ,ch2     )                  !out
! -------------------------------------------------------------------------------------------------
! computing surface drag coefficient cm for momentum and ch for heat
! -------------------------------------------------------------------------------------------------
    implicit none
! -------------------------------------------------------------------------------------------------
! inputs
    integer,              intent(in) :: iter   !iteration index
    real,                 intent(in) :: zlvl   !reference height  (m)
    real,                 intent(in) :: zpd    !zero plane displacement (m)
    real,                 intent(in) :: z0h    !roughness length, sensible heat, ground (m)
    real,                 intent(in) :: z0m    !roughness length, momentum, ground (m)
    real,                 intent(in) :: qair   !specific humidity at reference height (kg/kg)
    real,                 intent(in) :: sfctmp !temperature at reference height (k)
    real,                 intent(in) :: h      !sensible heat flux (w/m2) [+ to atm]
    real,                 intent(in) :: rhoair !density air (kg/m**3)
    real,                 intent(in) :: mpe    !prevents overflow error if division by zero
    real,                 intent(in) :: ur     !wind speed (m/s)

! in & out
    real,              intent(inout) :: moz    !monin-obukhov stability (z/l)
    integer,           intent(inout) :: mozsgn !number of times moz changes sign
    real,              intent(inout) :: fm     !momentum stability correction, weighted by prior iters
    real,              intent(inout) :: fh     !sen heat stability correction, weighted by prior iters
    real,              intent(inout) :: fm2    !sen heat stability correction, weighted by prior iters
    real,              intent(inout) :: fh2    !sen heat stability correction, weighted by prior iters

#ifdef CCPP  
    character(len=*),  intent(inout) :: errmsg
    integer,           intent(inout) :: errflg
#endif

! outputs
    real,                intent(out) :: fv     !friction velocity (m/s)
    real,                intent(out) :: cm     !drag coefficient for momentum
    real,                intent(out) :: ch     !drag coefficient for heat
    real,                intent(out) :: ch2    !drag coefficient for heat

! locals
    real    :: mozold                   !monin-obukhov stability parameter from prior iteration
    real    :: tmpcm                    !temporary calculation for cm
    real    :: tmpch                    !temporary calculation for ch
    real    :: mol                      !monin-obukhov length (m)
    real    :: tvir                     !temporary virtual temperature (k)
    real    :: tmp1,tmp2,tmp3           !temporary calculation
    real    :: fmnew                    !stability correction factor, momentum, for current moz
    real    :: fhnew                    !stability correction factor, sen heat, for current moz
    real    :: moz2                     !2/l
    real    :: tmpcm2                   !temporary calculation for cm2
    real    :: tmpch2                   !temporary calculation for ch2
    real    :: fm2new                   !stability correction factor, momentum, for current moz
    real    :: fh2new                   !stability correction factor, sen heat, for current moz
    real    :: tmp12,tmp22,tmp32        !temporary calculation

    real    :: cmfm, chfh, cm2fm2, ch2fh2


! -------------------------------------------------------------------------------------------------
! monin-obukhov stability parameter moz for next iteration

    mozold = moz
  
    if(zlvl <= zpd) then
       write(*,*) 'critical glacier problem: zlvl <= zpd; model stops', zlvl, zpd
#ifdef CCPP
       errflg = 1
       errmsg = "stop in noah-mp glacier"
       return 
#else
       call wrf_error_fatal("stop in noah-mp glacier")
#endif
    endif

    tmpcm = log((zlvl-zpd) / z0m)
    tmpch = log((zlvl-zpd) / z0h)
    tmpcm2 = log((2.0 + z0m) / z0m)
    tmpch2 = log((2.0 + z0h) / z0h)

    if(iter == 1) then
       fv   = 0.0
       moz  = 0.0
       mol  = 0.0
       moz2 = 0.0
    else
       tvir = (1. + 0.61*qair) * sfctmp
       tmp1 = vkc * (grav/tvir) * h/(rhoair*cpair)
       if (abs(tmp1) .le. mpe) tmp1 = mpe
       mol  = -1. * fv**3 / tmp1
       moz  = min( (zlvl-zpd)/mol, 1.)
       moz2  = min( (2.0 + z0h)/mol, 1.)
    endif

! accumulate number of times moz changes sign.

    if (mozold*moz .lt. 0.) mozsgn = mozsgn+1
    if (mozsgn .ge. 2) then
       moz = 0.
       fm = 0.
       fh = 0.
       moz2 = 0.
       fm2 = 0.
       fh2 = 0.
    endif

! evaluate stability-dependent variables using moz from prior iteration
    if (moz .lt. 0.) then
       tmp1 = (1. - 16.*moz)**0.25
       tmp2 = log((1.+tmp1*tmp1)/2.)
       tmp3 = log((1.+tmp1)/2.)
       fmnew = 2.*tmp3 + tmp2 - 2.*atan(tmp1) + 1.5707963
       fhnew = 2*tmp2

! 2-meter
       tmp12 = (1. - 16.*moz2)**0.25
       tmp22 = log((1.+tmp12*tmp12)/2.)
       tmp32 = log((1.+tmp12)/2.)
       fm2new = 2.*tmp32 + tmp22 - 2.*atan(tmp12) + 1.5707963
       fh2new = 2*tmp22
    else
       fmnew = -5.*moz
       fhnew = fmnew
       fm2new = -5.*moz2
       fh2new = fm2new
    endif

! except for first iteration, weight stability factors for previous
! iteration to help avoid flip-flops from one iteration to the next

    if (iter == 1) then
       fm = fmnew
       fh = fhnew
       fm2 = fm2new
       fh2 = fh2new
    else
       fm = 0.5 * (fm+fmnew)
       fh = 0.5 * (fh+fhnew)
       fm2 = 0.5 * (fm2+fm2new)
       fh2 = 0.5 * (fh2+fh2new)
    endif

! exchange coefficients

    fh = min(fh,0.9*tmpch)
    fm = min(fm,0.9*tmpcm)
    fh2 = min(fh2,0.9*tmpch2)
    fm2 = min(fm2,0.9*tmpcm2)

    cmfm = tmpcm-fm
    chfh = tmpch-fh
    cm2fm2 = tmpcm2-fm2
    ch2fh2 = tmpch2-fh2
    if(abs(cmfm) <= mpe) cmfm = mpe
    if(abs(chfh) <= mpe) chfh = mpe
    if(abs(cm2fm2) <= mpe) cm2fm2 = mpe
    if(abs(ch2fh2) <= mpe) ch2fh2 = mpe
    cm  = vkc*vkc/(cmfm*cmfm)
    ch  = vkc*vkc/(cmfm*chfh)
    ch2  = vkc*vkc/(cm2fm2*ch2fh2)
        
! friction velocity

    fv = ur * sqrt(cm)
    ch2  = vkc*fv/ch2fh2

  end subroutine sfcdif1_glacier
! ==================================================================================================
!>\ingroup NoahMP_LSM
  subroutine tsnosoi_glacier (nsoil   ,nsnow   ,isnow   ,dt      ,tbot    , & !in
                              ssoil   ,snowh   ,zbot    ,zsnso   ,df      , & !in
			      hcpct   ,                                     & !in
                              stc     )                                       !inout
! --------------------------------------------------------------------------------------------------
! compute snow (up to 3l) and soil (4l) temperature. note that snow temperatures
! during melting season may exceed melting point (tfrz) but later in phasechange
! subroutine the snow temperatures are reset to tfrz for melting snow.
! --------------------------------------------------------------------------------------------------
  implicit none
! --------------------------------------------------------------------------------------------------
!input

    integer,                         intent(in)  :: nsoil  !no of soil layers (4)
    integer,                         intent(in)  :: nsnow  !maximum no of snow layers (3)
    integer,                         intent(in)  :: isnow  !actual no of snow layers

    real,                            intent(in)  :: dt     !time step (s)
    real,                            intent(in)  :: tbot   !
    real,                            intent(in)  :: ssoil  !ground heat flux (w/m2)
    real,                            intent(in)  :: snowh  !snow depth (m)
    real,                            intent(in)  :: zbot   !from soil surface (m)
    real, dimension(-nsnow+1:nsoil), intent(in)  :: zsnso  !layer-bot. depth from snow surf.(m)
    real, dimension(-nsnow+1:nsoil), intent(in)  :: df     !thermal conductivity
    real, dimension(-nsnow+1:nsoil), intent(in)  :: hcpct  !heat capacity (j/m3/k)

!input and output

    real, dimension(-nsnow+1:nsoil), intent(inout) :: stc

!local

    integer                                      :: iz
    real                                         :: zbotsno   !zbot from snow surface
    real, dimension(-nsnow+1:nsoil)              :: ai, bi, ci, rhsts
    real                                         :: eflxb !energy influx from soil bottom (w/m2)
    real, dimension(-nsnow+1:nsoil)              :: phi   !light through water (w/m2)

! ----------------------------------------------------------------------

! prescribe solar penetration into ice/snow

    phi(isnow+1:nsoil) = 0.

! adjust zbot from soil surface to zbotsno from snow surface

    zbotsno = zbot - snowh    !from snow surface

! compute ice temperatures

      call hrt_glacier   (nsnow     ,nsoil     ,isnow     ,zsnso     , &
                          stc       ,tbot      ,zbotsno   ,df        , &
                          hcpct     ,ssoil     ,phi       ,            &
                          ai        ,bi        ,ci        ,rhsts     , &
                          eflxb     )

      call hstep_glacier (nsnow     ,nsoil     ,isnow     ,dt        , &
                          ai        ,bi        ,ci        ,rhsts     , &
                          stc       ) 

  end subroutine tsnosoi_glacier
! ==================================================================================================
! ----------------------------------------------------------------------
!>\ingroup NoahMP_LSM
  subroutine hrt_glacier (nsnow     ,nsoil     ,isnow     ,zsnso     , & !in
                          stc       ,tbot      ,zbot      ,df        , & !in
                          hcpct     ,ssoil     ,phi       ,            & !in
                          ai        ,bi        ,ci        ,rhsts     , & !out
                          botflx    )                                    !out
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! calculate the right hand side of the time tendency term of the soil
! thermal diffusion equation.  also to compute ( prepare ) the matrix
! coefficients for the tri-diagonal matrix of the implicit time scheme.
! ----------------------------------------------------------------------
    implicit none
! ----------------------------------------------------------------------
! input

    integer,                         intent(in)  :: nsoil  !no of soil layers (4)
    integer,                         intent(in)  :: nsnow  !maximum no of snow layers (3)
    integer,                         intent(in)  :: isnow  !actual no of snow layers
    real,                            intent(in)  :: tbot   !bottom soil temp. at zbot (k)
    real,                            intent(in)  :: zbot   !depth of lower boundary condition (m)
                                                           !from soil surface not snow surface
    real,                            intent(in)  :: ssoil  !ground heat flux (w/m2)
    real, dimension(-nsnow+1:nsoil), intent(in)  :: zsnso  !depth of layer-bottom of snow/soil (m)
    real, dimension(-nsnow+1:nsoil), intent(in)  :: stc    !snow/soil temperature (k)
    real, dimension(-nsnow+1:nsoil), intent(in)  :: df     !thermal conductivity [w/m/k]
    real, dimension(-nsnow+1:nsoil), intent(in)  :: hcpct  !heat capacity [j/m3/k]
    real, dimension(-nsnow+1:nsoil), intent(in)  :: phi    !light through water (w/m2)

! output

    real, dimension(-nsnow+1:nsoil), intent(out) :: rhsts  !right-hand side of the matrix
    real, dimension(-nsnow+1:nsoil), intent(out) :: ai     !left-hand side coefficient
    real, dimension(-nsnow+1:nsoil), intent(out) :: bi     !left-hand side coefficient
    real, dimension(-nsnow+1:nsoil), intent(out) :: ci     !left-hand side coefficient
    real,                            intent(out) :: botflx !energy influx from soil bottom (w/m2)

! local

    integer                                      :: k
    real, dimension(-nsnow+1:nsoil)              :: ddz
    real, dimension(-nsnow+1:nsoil)              :: denom
    real, dimension(-nsnow+1:nsoil)              :: dtsdz
    real, dimension(-nsnow+1:nsoil)              :: eflux
    real                                         :: temp1
! ----------------------------------------------------------------------

    do k = isnow+1, nsoil
        if (k == isnow+1) then
           denom(k)  = - zsnso(k) * hcpct(k)
           temp1     = - zsnso(k+1)
           ddz(k)    = 2.0 / temp1
           dtsdz(k)  = 2.0 * (stc(k) - stc(k+1)) / temp1
           eflux(k)  = df(k) * dtsdz(k) - ssoil - phi(k)
        else if (k < nsoil) then
           denom(k)  = (zsnso(k-1) - zsnso(k)) * hcpct(k)
           temp1     = zsnso(k-1) - zsnso(k+1)
           ddz(k)    = 2.0 / temp1
           dtsdz(k)  = 2.0 * (stc(k) - stc(k+1)) / temp1
           eflux(k)  = (df(k)*dtsdz(k) - df(k-1)*dtsdz(k-1)) - phi(k)
        else if (k == nsoil) then
           denom(k)  = (zsnso(k-1) - zsnso(k)) * hcpct(k)
           temp1     =  zsnso(k-1) - zsnso(k)
           if(opt_tbot == 1) then
               botflx     = 0. 
           end if
           if(opt_tbot == 2) then
               dtsdz(k)  = (stc(k) - tbot) / ( 0.5*(zsnso(k-1)+zsnso(k)) - zbot)
               botflx    = -df(k) * dtsdz(k)
           end if
           eflux(k)  = (-botflx - df(k-1)*dtsdz(k-1) ) - phi(k)
        end if
    end do

    do k = isnow+1, nsoil
        if (k == isnow+1) then
           ai(k)    =   0.0
           ci(k)    = - df(k)   * ddz(k) / denom(k)
           if (opt_stc == 1) then
              bi(k) = - ci(k)
           end if                                        
           if (opt_stc == 2) then
              bi(k) = - ci(k) + df(k)/(0.5*zsnso(k)*zsnso(k)*hcpct(k))
           end if
        else if (k < nsoil) then
           ai(k)    = - df(k-1) * ddz(k-1) / denom(k) 
           ci(k)    = - df(k  ) * ddz(k  ) / denom(k) 
           bi(k)    = - (ai(k) + ci (k))
        else if (k == nsoil) then
           ai(k)    = - df(k-1) * ddz(k-1) / denom(k) 
           ci(k)    = 0.0
           bi(k)    = - (ai(k) + ci(k))
        end if
           rhsts(k)  = eflux(k)/ (-denom(k))
    end do

  end subroutine hrt_glacier
! ==================================================================================================
! ----------------------------------------------------------------------
!>\ingroup NoahMP_LSM
  subroutine hstep_glacier (nsnow     ,nsoil     ,isnow     ,dt        ,  & !in
                            ai        ,bi        ,ci        ,rhsts     ,  & !inout
                            stc       )                                     !inout
! ----------------------------------------------------------------------
! calculate/update the soil temperature field.
! ----------------------------------------------------------------------
    implicit none
! ----------------------------------------------------------------------
! input

    integer,                         intent(in)    :: nsoil
    integer,                         intent(in)    :: nsnow
    integer,                         intent(in)    :: isnow
    real,                            intent(in)    :: dt

! output & input
    real, dimension(-nsnow+1:nsoil), intent(inout) :: ai
    real, dimension(-nsnow+1:nsoil), intent(inout) :: bi
    real, dimension(-nsnow+1:nsoil), intent(inout) :: ci
    real, dimension(-nsnow+1:nsoil), intent(inout) :: stc
    real, dimension(-nsnow+1:nsoil), intent(inout) :: rhsts

! local
    integer                                        :: k
    real, dimension(-nsnow+1:nsoil)                :: rhstsin
    real, dimension(-nsnow+1:nsoil)                :: ciin
! ----------------------------------------------------------------------

    do k = isnow+1,nsoil
       rhsts(k) =   rhsts(k) * dt
       ai(k)    =      ai(k) * dt
       bi(k)    = 1. + bi(k) * dt
       ci(k)    =      ci(k) * dt
    end do

! copy values for input variables before call to rosr12

    do k = isnow+1,nsoil
       rhstsin(k) = rhsts(k)
       ciin(k)    = ci(k)
    end do

! solve the tri-diagonal matrix equation

    call rosr12_glacier (ci,ai,bi,ciin,rhstsin,rhsts,isnow+1,nsoil,nsnow)

! update snow & soil temperature

    do k = isnow+1,nsoil
       stc (k) = stc (k) + ci (k)
    end do

  end subroutine hstep_glacier
! ==================================================================================================
!>\ingroup NoahMP_LSM
  subroutine rosr12_glacier (p,a,b,c,d,delta,ntop,nsoil,nsnow)
! ----------------------------------------------------------------------
! subroutine rosr12
! ----------------------------------------------------------------------
! invert (solve) the tri-diagonal matrix problem shown below:
! ###                                            ### ###  ###   ###  ###
! #b(1), c(1),  0  ,  0  ,  0  ,   . . .  ,    0   # #      #   #      #
! #a(2), b(2), c(2),  0  ,  0  ,   . . .  ,    0   # #      #   #      #
! # 0  , a(3), b(3), c(3),  0  ,   . . .  ,    0   # #      #   # d(3) #
! # 0  ,  0  , a(4), b(4), c(4),   . . .  ,    0   # # p(4) #   # d(4) #
! # 0  ,  0  ,  0  , a(5), b(5),   . . .  ,    0   # # p(5) #   # d(5) #
! # .                                          .   # #  .   # = #   .  #
! # .                                          .   # #  .   #   #   .  #
! # .                                          .   # #  .   #   #   .  #
! # 0  , . . . , 0 , a(m-2), b(m-2), c(m-2),   0   # #p(m-2)#   #d(m-2)#
! # 0  , . . . , 0 ,   0   , a(m-1), b(m-1), c(m-1)# #p(m-1)#   #d(m-1)#
! # 0  , . . . , 0 ,   0   ,   0   ,  a(m) ,  b(m) # # p(m) #   # d(m) #
! ###                                            ### ###  ###   ###  ###
! ----------------------------------------------------------------------
    implicit none

    integer, intent(in)   :: ntop           
    integer, intent(in)   :: nsoil,nsnow
    integer               :: k, kk

    real, dimension(-nsnow+1:nsoil),intent(in):: a, b, d
    real, dimension(-nsnow+1:nsoil),intent(inout):: c,p,delta

! ----------------------------------------------------------------------
! initialize eqn coef c for the lowest soil layer
! ----------------------------------------------------------------------
    c (nsoil) = 0.0
    p (ntop) = - c (ntop) / b (ntop)
! ----------------------------------------------------------------------
! solve the coefs for the 1st soil layer
! ----------------------------------------------------------------------
    delta (ntop) = d (ntop) / b (ntop)
! ----------------------------------------------------------------------
! solve the coefs for soil layers 2 thru nsoil
! ----------------------------------------------------------------------
    do k = ntop+1,nsoil
       p (k) = - c (k) * ( 1.0 / (b (k) + a (k) * p (k -1)) )
       delta (k) = (d (k) - a (k)* delta (k -1))* (1.0/ (b (k) + a (k)&
            * p (k -1)))
    end do
! ----------------------------------------------------------------------
! set p to delta for lowest soil layer
! ----------------------------------------------------------------------
    p (nsoil) = delta (nsoil)
! ----------------------------------------------------------------------
! adjust p for soil layers 2 thru nsoil
! ----------------------------------------------------------------------
    do k = ntop+1,nsoil
       kk = nsoil - k + (ntop-1) + 1
       p (kk) = p (kk) * p (kk +1) + delta (kk)
    end do
! ----------------------------------------------------------------------
  end subroutine rosr12_glacier
! ----------------------------------------------------------------------
! ==================================================================================================
!>\ingroup NoahMP_LSM
  subroutine phasechange_glacier (nsnow   ,nsoil   ,isnow   ,dt      ,fact    , & !in
                                  dzsnso  ,                                     & !in
                                  stc     ,snice   ,snliq   ,sneqv   ,snowh   , & !inout
                                  smc     ,sh2o    ,                            & !inout
                                  qmelt   ,imelt   ,ponding )                     !out
! ----------------------------------------------------------------------
! melting/freezing of snow water and soil water
! ----------------------------------------------------------------------
  implicit none
! ----------------------------------------------------------------------
! inputs

  integer, intent(in)                             :: nsnow  !maximum no. of snow layers [=3]
  integer, intent(in)                             :: nsoil  !no. of soil layers [=4]
  integer, intent(in)                             :: isnow  !actual no. of snow layers [<=3]
  real, intent(in)                                :: dt     !land model time step (sec)
  real, dimension(-nsnow+1:nsoil), intent(in)     :: fact   !temporary
  real, dimension(-nsnow+1:nsoil), intent(in)     :: dzsnso !snow/soil layer thickness [m]

! inputs/outputs

  real, dimension(-nsnow+1:nsoil), intent(inout)  :: stc    !snow/soil layer temperature [k]
  real, dimension(-nsnow+1:0)    , intent(inout)  :: snice  !snow layer ice [mm]
  real, dimension(-nsnow+1:0)    , intent(inout)  :: snliq  !snow layer liquid water [mm]
  real, intent(inout)                             :: sneqv
  real, intent(inout)                             :: snowh
  real, dimension(       1:nsoil), intent(inout)  :: sh2o   !soil liquid water [m3/m3]
  real, dimension(       1:nsoil), intent(inout)  :: smc    !total soil water [m3/m3]

! outputs
  real,                               intent(out) :: qmelt  !snowmelt rate [mm/s]
  integer, dimension(-nsnow+1:nsoil), intent(out) :: imelt  !phase change index
  real,                               intent(out) :: ponding!snowmelt when snow has no layer [mm]

! local

  integer                         :: j,k         !do loop index
  real, dimension(-nsnow+1:nsoil) :: hm        !energy residual [w/m2]
  real, dimension(-nsnow+1:nsoil) :: xm        !melting or freezing water [kg/m2]
  real, dimension(-nsnow+1:nsoil) :: wmass0
  real, dimension(-nsnow+1:nsoil) :: wice0 
  real, dimension(-nsnow+1:nsoil) :: wliq0 
  real, dimension(-nsnow+1:nsoil) :: mice      !soil/snow ice mass [mm]
  real, dimension(-nsnow+1:nsoil) :: mliq      !soil/snow liquid water mass [mm]
  real, dimension(-nsnow+1:nsoil) :: heatr     !energy residual or loss after melting/freezing
  real                            :: temp1     !temporary variables [kg/m2]
  real                            :: propor
  real                            :: xmf       !total latent heat of phase change

! ----------------------------------------------------------------------
! initialization

    qmelt   = 0.
    ponding = 0.
    xmf     = 0.

    do j = isnow+1,0           ! all snow layers
         mice(j) = snice(j)
         mliq(j) = snliq(j)
    end do

    do j = 1, nsoil            ! all soil layers
         mliq(j) =  sh2o(j)            * dzsnso(j) * 1000.
         mice(j) = (smc(j) - sh2o(j))  * dzsnso(j) * 1000.
    end do

    do j = isnow+1,nsoil       ! all layers
         imelt(j)    = 0
         hm(j)       = 0.
         xm(j)       = 0.
         wice0(j)    = mice(j)
         wliq0(j)    = mliq(j)
         wmass0(j)   = mice(j) + mliq(j)
    enddo
    
    do j = isnow+1,nsoil
         if (mice(j) > 0. .and. stc(j) >= tfrz) then  ! melting 
             imelt(j) = 1
         endif
         if (mliq(j) > 0. .and. stc(j)  < tfrz) then  ! freezing 
             imelt(j) = 2
         endif

         ! if snow exists, but its thickness is not enough to create a layer
         if (isnow == 0 .and. sneqv > 0. .and. j == 1) then
             if (stc(j) >= tfrz) then
                imelt(j) = 1
             endif
         endif
    enddo

! calculate the energy surplus and loss for melting and freezing

    do j = isnow+1,nsoil
         if (imelt(j) > 0) then
             hm(j) = (stc(j)-tfrz)/fact(j)
             stc(j) = tfrz
         endif

         if (imelt(j) == 1 .and. hm(j) < 0.) then
            hm(j) = 0.
            imelt(j) = 0
         endif
         if (imelt(j) == 2 .and. hm(j) > 0.) then
            hm(j) = 0.
            imelt(j) = 0
         endif
         xm(j) = hm(j)*dt/hfus                           
    enddo

! the rate of melting and freezing for snow without a layer, needs more work.

    if (isnow == 0 .and. sneqv > 0. .and. xm(1) > 0.) then  
        temp1  = sneqv
        sneqv  = max(0.,temp1-xm(1))  
        propor = sneqv/temp1
        snowh  = max(0.,propor * snowh)
        heatr(1)  = hm(1) - hfus*(temp1-sneqv)/dt  
        if (heatr(1) > 0.) then
              xm(1) = heatr(1)*dt/hfus             
              hm(1) = heatr(1) 
	      imelt(1) = 1                   
        else
              xm(1) = 0.
              hm(1) = 0.
	      imelt(1) = 0                   
        endif
        qmelt   = max(0.,(temp1-sneqv))/dt
        xmf     = hfus*qmelt
        ponding = temp1-sneqv
    endif

! the rate of melting and freezing for snow and soil

    do j = isnow+1,nsoil
      if (imelt(j) > 0 .and. abs(hm(j)) > 0.) then

         heatr(j) = 0.
         if (xm(j) > 0.) then                            
            mice(j) = max(0., wice0(j)-xm(j))
            heatr(j) = hm(j) - hfus*(wice0(j)-mice(j))/dt
         else if (xm(j) < 0.) then                      
            mice(j) = min(wmass0(j), wice0(j)-xm(j))  
            heatr(j) = hm(j) - hfus*(wice0(j)-mice(j))/dt
         endif

         mliq(j) = max(0.,wmass0(j)-mice(j))

         if (abs(heatr(j)) > 0.) then
            stc(j) = stc(j) + fact(j)*heatr(j)
            if (j <= 0) then                             ! snow
               if (mliq(j)*mice(j)>0.) stc(j) = tfrz
            end if
         endif

         if (j > 0) xmf = xmf + hfus * (wice0(j)-mice(j))/dt

         if (j < 1) then
            qmelt = qmelt + max(0.,(wice0(j)-mice(j)))/dt
         endif
      endif
    enddo
    heatr = 0.0
    xm = 0.0

! deal with residuals in ice/soil

! first remove excess heat by reducing temperature of layers

    if (any(stc(1:4) > tfrz) .and. any(stc(1:4) < tfrz)) then
      do j = 1,nsoil
        if ( stc(j) > tfrz ) then                                       
	  heatr(j) = (stc(j)-tfrz)/fact(j)
          do k = 1,nsoil
	    if (j .ne. k .and. stc(k) < tfrz .and. heatr(j) > 0.1) then
	      heatr(k) = (stc(k)-tfrz)/fact(k)
	      if (abs(heatr(k)) > heatr(j)) then  ! layer absorbs all
	        heatr(k) = heatr(k) + heatr(j)
		stc(k) = tfrz + heatr(k)*fact(k)
		heatr(j) = 0.0
              else
	        heatr(j) = heatr(j) + heatr(k)
		heatr(k) = 0.0
		stc(k) = tfrz
              end if
	    end if
	  end do
          stc(j) = tfrz + heatr(j)*fact(j)
        end if
      end do
    end if

! now remove excess cold by increasing temperature of layers (may not be necessary with above loop)

    if (any(stc(1:4) > tfrz) .and. any(stc(1:4) < tfrz)) then
      do j = 1,nsoil
        if ( stc(j) < tfrz ) then                                       
	  heatr(j) = (stc(j)-tfrz)/fact(j)
          do k = 1,nsoil
	    if (j .ne. k .and. stc(k) > tfrz .and. heatr(j) < -0.1) then
	      heatr(k) = (stc(k)-tfrz)/fact(k)
	      if (heatr(k) > abs(heatr(j))) then  ! layer absorbs all
	        heatr(k) = heatr(k) + heatr(j)
		stc(k) = tfrz + heatr(k)*fact(k)
		heatr(j) = 0.0
              else
	        heatr(j) = heatr(j) + heatr(k)
		heatr(k) = 0.0
		stc(k) = tfrz
              end if
	    end if
	  end do
          stc(j) = tfrz + heatr(j)*fact(j)
        end if
      end do
    end if

! now remove excess heat by melting ice

    if (any(stc(1:4) > tfrz) .and. any(mice(1:4) > 0.)) then
      do j = 1,nsoil
        if ( stc(j) > tfrz ) then                                       
	  heatr(j) = (stc(j)-tfrz)/fact(j)
          xm(j) = heatr(j)*dt/hfus                           
          do k = 1,nsoil
	    if (j .ne. k .and. mice(k) > 0. .and. xm(j) > 0.1) then
	      if (mice(k) > xm(j)) then  ! layer absorbs all
	        mice(k) = mice(k) - xm(j)
		xmf = xmf + hfus * xm(j)/dt
		stc(k) = tfrz
		xm(j) = 0.0
              else
	        xm(j) = xm(j) - mice(k)
		xmf = xmf + hfus * mice(k)/dt
		mice(k) = 0.0
		stc(k) = tfrz
              end if
              mliq(k) = max(0.,wmass0(k)-mice(k))
	    end if
	  end do
	  heatr(j) = xm(j)*hfus/dt
          stc(j) = tfrz + heatr(j)*fact(j)
        end if
      end do
    end if

! now remove excess cold by freezing liquid of layers (may not be necessary with above loop)

    if (any(stc(1:4) < tfrz) .and. any(mliq(1:4) > 0.)) then
      do j = 1,nsoil
        if ( stc(j) < tfrz ) then                                       
	  heatr(j) = (stc(j)-tfrz)/fact(j)
          xm(j) = heatr(j)*dt/hfus                           
          do k = 1,nsoil
	    if (j .ne. k .and. mliq(k) > 0. .and. xm(j) < -0.1) then
	      if (mliq(k) > abs(xm(j))) then  ! layer absorbs all
	        mice(k) = mice(k) - xm(j)
		xmf = xmf + hfus * xm(j)/dt
		stc(k) = tfrz
		xm(j) = 0.0
              else
	        xm(j) = xm(j) + mliq(k)
		xmf = xmf - hfus * mliq(k)/dt
		mice(k) = wmass0(k)
		stc(k) = tfrz
              end if
              mliq(k) = max(0.,wmass0(k)-mice(k))
	    end if
	  end do
	  heatr(j) = xm(j)*hfus/dt
          stc(j) = tfrz + heatr(j)*fact(j)
        end if
      end do
    end if

    do j = isnow+1,0             ! snow
       snliq(j) = mliq(j)
       snice(j) = mice(j)
    end do

    do j = 1, nsoil              ! soil
       sh2o(j) =  mliq(j)            / (1000. * dzsnso(j))
       sh2o(j) =  max(0.0,min(1.0,sh2o(j)))
!       smc(j)  = (mliq(j) + mice(j)) / (1000. * dzsnso(j))
       smc(j)  = 1.0 
    end do
   
  end subroutine phasechange_glacier
! ==================================================================================================
!>\ingroup NoahMP_LSM
  subroutine water_glacier (nsnow  ,nsoil  ,imelt  ,dt     ,prcp   ,sfctmp , & !in
                            qvap   ,qdew   ,ficeold,zsoil  ,                 & !in
                            isnow  ,snowh  ,sneqv  ,snice  ,snliq  ,stc    , & !inout
                            dzsnso ,sh2o   ,sice   ,ponding,zsnso  ,         & !inout
                            runsrf ,runsub ,qsnow  ,ponding1 ,ponding2,qsnbot,fpice,esnow     &   !out
                            )  !out
! ----------------------------------------------------------------------  
! code history:
! initial code: guo-yue niu, oct. 2007
! ----------------------------------------------------------------------
  implicit none
! ----------------------------------------------------------------------
! input
  integer,                         intent(in)    :: nsnow   !maximum no. of snow layers
  integer,                         intent(in)    :: nsoil   !no. of soil layers
  integer, dimension(-nsnow+1:0) , intent(in)    :: imelt   !melting state index [1-melt; 2-freeze]
  real,                            intent(in)    :: dt      !main time step (s)
  real,                            intent(in)    :: prcp    !precipitation (mm/s)
  real,                            intent(in)    :: sfctmp  !surface air temperature [k]
  real,                            intent(in)    :: qvap    !soil surface evaporation rate[mm/s]
  real,                            intent(in)    :: qdew    !soil surface dew rate[mm/s]
  real, dimension(-nsnow+1:    0), intent(in)    :: ficeold !ice fraction at last timestep
  real, dimension(       1:nsoil), intent(in)    :: zsoil  !layer-bottom depth from soil surf (m)

! input/output
  integer,                         intent(inout) :: isnow   !actual no. of snow layers
  real,                            intent(inout) :: snowh   !snow height [m]
  real,                            intent(inout) :: sneqv   !snow water eqv. [mm]
  real, dimension(-nsnow+1:    0), intent(inout) :: snice   !snow layer ice [mm]
  real, dimension(-nsnow+1:    0), intent(inout) :: snliq   !snow layer liquid water [mm]
  real, dimension(-nsnow+1:nsoil), intent(inout) :: stc     !snow/soil layer temperature [k]
  real, dimension(-nsnow+1:nsoil), intent(inout) :: dzsnso  !snow/soil layer thickness [m]
  real, dimension(       1:nsoil), intent(inout) :: sh2o    !soil liquid water content [m3/m3]
  real, dimension(       1:nsoil), intent(inout) :: sice    !soil ice content [m3/m3]
  real                           , intent(inout) :: ponding ![mm]
  real, dimension(-nsnow+1:nsoil), intent(inout) :: zsnso   !layer-bottom depth from snow surf [m]

! output
  real,                            intent(out)   :: runsrf  !surface runoff [mm/s] 
  real,                            intent(out)   :: runsub  !baseflow (sturation excess) [mm/s]
  real,                            intent(out)   :: qsnow   !snow at ground srf (mm/s) [+]
  real,                            intent(out)   :: ponding1
  real,                            intent(out)   :: ponding2
  real,                            intent(out)   :: qsnbot  !melting water out of snow bottom [mm/s]
  real,                            intent(out)   :: fpice   !precipitation frozen fraction
  real,                            intent(out)   :: esnow   !

! local
  real                                           :: qrain   !rain at ground srf (mm) [+]
  real                                           :: qseva   !soil surface evap rate [mm/s]
  real                                           :: qsdew   !soil surface dew rate [mm/s]
  real                                           :: qsnfro  !snow surface frost rate[mm/s]
  real                                           :: qsnsub  !snow surface sublimation rate [mm/s]
  real                                           :: snowhin !snow depth increasing rate (m/s)
  real                                           :: snoflow !glacier flow [mm/s]
  real                                           :: bdfall  !density of new snow (mm water/m snow)
  real                                           :: replace !replacement water due to sublimation of glacier
  real, dimension(       1:nsoil)                :: sice_save  !soil ice content [m3/m3]
  real, dimension(       1:nsoil)                :: sh2o_save  !soil liquid water content [m3/m3]
  integer :: ilev


! ----------------------------------------------------------------------
! initialize

   snoflow         = 0.
   runsub          = 0.
   runsrf          = 0.
   sice_save       = sice
   sh2o_save       = sh2o

! --------------------------------------------------------------------
! partition precipitation into rain and snow (from canwater)

! jordan (1991)

     if(opt_snf == 1 .or. opt_snf == 4) then
       if(sfctmp > tfrz+2.5)then
           fpice = 0.
       else
         if(sfctmp <= tfrz+0.5)then
           fpice = 1.0
         else if(sfctmp <= tfrz+2.)then
           fpice = 1.-(-54.632 + 0.2*sfctmp)
         else
           fpice = 0.6
         endif
       endif
     endif

     if(opt_snf == 2) then
       if(sfctmp >= tfrz+2.2) then
           fpice = 0.
       else
           fpice = 1.0
       endif
     endif

     if(opt_snf == 3) then
       if(sfctmp >= tfrz) then
           fpice = 0.
       else
           fpice = 1.0
       endif
     endif
!     print*, 'fpice: ',fpice

! hedstrom nr and jw pomeroy (1998), hydrol. processes, 12, 1611-1625
! fresh snow density

     bdfall = min(120.,67.92+51.25*exp((sfctmp-tfrz)/2.59)) !mb: change to min v3.7

     qrain   = prcp * (1.-fpice)
     qsnow   = prcp * fpice
     snowhin = qsnow/bdfall
!     print *, 'qrain, qsnow',qrain,qsnow,qrain*dt,qsnow*dt

! sublimation, frost, evaporation, and dew

!     qsnsub = 0.
!     if (sneqv > 0.) then
!       qsnsub = min(qvap, sneqv/dt)
!     endif
!     qseva = qvap-qsnsub

!     qsnfro = 0.
!     if (sneqv > 0.) then
!        qsnfro = qdew
!     endif
!     qsdew = qdew - qsnfro

     qsnsub = qvap  ! send total sublimation/frost to snowwater and deal with it there
     qsnfro = qdew
     esnow = qsnsub*2.83e+6


!     print *, 'qvap',qvap,qvap*dt
!     print *, 'qsnsub',qsnsub,qsnsub*dt
!     print *, 'qseva',qseva,qseva*dt
!     print *, 'qsnfro',qsnfro,qsnfro*dt
!     print *, 'qdew',qdew,qdew*dt
!     print *, 'qsdew',qsdew,qsdew*dt
!print *, 'before snowwater', sneqv,snowh,snice,snliq,sh2o,sice
     call snowwater_glacier (nsnow  ,nsoil  ,imelt  ,dt     ,sfctmp , & !in
                             snowhin,qsnow  ,qsnfro ,qsnsub ,qrain  , & !in
                             ficeold,zsoil  ,                         & !in
                             isnow  ,snowh  ,sneqv  ,snice  ,snliq  , & !inout
                             sh2o   ,sice   ,stc    ,dzsnso ,zsnso  , & !inout
                             qsnbot ,snoflow,ponding1       ,ponding2)  !out
!print *, 'after snowwater', sneqv,snowh,snice,snliq,sh2o,sice
!print *, 'ponding', ponding,ponding1,ponding2

    !ponding: melting water from snow when there is no layer
    
    runsrf = (ponding+ponding1+ponding2)/dt

    if(isnow == 0) then
      runsrf = runsrf + qsnbot + qrain
    else
      runsrf = runsrf + qsnbot
    endif

    
    replace = 0.0
    do ilev = 1,nsoil
       replace = replace + dzsnso(ilev)*(sice(ilev) - sice_save(ilev) + sh2o(ilev) - sh2o_save(ilev))
    end do
    replace = replace * 1000.0 / dt     ! convert to [mm/s]
    
    sice = min(1.0,sice_save)
    sh2o = 1.0 - sice
!print *, 'replace', replace
    
    ! use runsub as a water balancer, snoflow is snow that disappears, replace is
    !   water from below that replaces glacier loss

    runsub       = snoflow + replace

  end subroutine water_glacier
! ==================================================================================================
! ----------------------------------------------------------------------
!>\ingroup NoahMP_LSM
  subroutine snowwater_glacier (nsnow  ,nsoil  ,imelt  ,dt     ,sfctmp , & !in
                                snowhin,qsnow  ,qsnfro ,qsnsub ,qrain  , & !in
                                ficeold,zsoil  ,                         & !in
                                isnow  ,snowh  ,sneqv  ,snice  ,snliq  , & !inout
                                sh2o   ,sice   ,stc    ,dzsnso ,zsnso  , & !inout
                                qsnbot ,snoflow,ponding1       ,ponding2)  !out
! ----------------------------------------------------------------------
  implicit none
! ----------------------------------------------------------------------
! input
  integer,                         intent(in)    :: nsnow  !maximum no. of snow layers
  integer,                         intent(in)    :: nsoil  !no. of soil layers
  integer, dimension(-nsnow+1:0) , intent(in)    :: imelt  !melting state index [0-no melt;1-melt]
  real,                            intent(in)    :: dt     !time step (s)
  real,                            intent(in)    :: sfctmp !surface air temperature [k]
  real,                            intent(in)    :: snowhin!snow depth increasing rate (m/s)
  real,                            intent(in)    :: qsnow  !snow at ground srf (mm/s) [+]
  real,                            intent(in)    :: qsnfro !snow surface frost rate[mm/s]
  real,                            intent(in)    :: qsnsub !snow surface sublimation rate[mm/s]
  real,                            intent(in)    :: qrain  !snow surface rain rate[mm/s]
  real, dimension(-nsnow+1:0)    , intent(in)    :: ficeold!ice fraction at last timestep
  real, dimension(       1:nsoil), intent(in)    :: zsoil  !layer-bottom depth from soil surf (m)

! input & output
  integer,                         intent(inout) :: isnow  !actual no. of snow layers
  real,                            intent(inout) :: snowh  !snow height [m]
  real,                            intent(inout) :: sneqv  !snow water eqv. [mm]
  real, dimension(-nsnow+1:    0), intent(inout) :: snice  !snow layer ice [mm]
  real, dimension(-nsnow+1:    0), intent(inout) :: snliq  !snow layer liquid water [mm]
  real, dimension(       1:nsoil), intent(inout) :: sh2o   !soil liquid moisture (m3/m3)
  real, dimension(       1:nsoil), intent(inout) :: sice   !soil ice moisture (m3/m3)
  real, dimension(-nsnow+1:nsoil), intent(inout) :: stc    !snow layer temperature [k]
  real, dimension(-nsnow+1:nsoil), intent(inout) :: dzsnso !snow/soil layer thickness [m]
  real, dimension(-nsnow+1:nsoil), intent(inout) :: zsnso  !layer-bottom depth from snow surf [m]

! output
  real,                              intent(out) :: qsnbot !melting water out of snow bottom [mm/s]
  real,                              intent(out) :: snoflow!glacier flow [mm]
  real,                              intent(out) :: ponding1
  real,                              intent(out) :: ponding2

! local
  integer :: iz
  real    :: bdsnow  !bulk density of snow (kg/m3)
! ----------------------------------------------------------------------
   snoflow = 0.0
   ponding1 = 0.0
   ponding2 = 0.0

   call snowfall_glacier (nsoil  ,nsnow  ,dt     ,qsnow  ,snowhin, & !in
                          sfctmp ,                                 & !in
                          isnow  ,snowh  ,dzsnso ,stc    ,snice  , & !inout
                          snliq  ,sneqv  )                           !inout

   if(isnow < 0) then        !when more than one layer
     call  compact_glacier (nsnow  ,nsoil  ,dt     ,stc    ,snice  , & !in
                            snliq  ,imelt  ,ficeold,                 & !in
                            isnow  ,dzsnso )                           !inout

     call  combine_glacier (nsnow  ,nsoil  ,                         & !in
                            isnow  ,sh2o   ,stc    ,snice  ,snliq  , & !inout
                            dzsnso ,sice   ,snowh  ,sneqv  ,         & !inout
                            ponding1       ,ponding2)                  !out

     call   divide_glacier (nsnow  ,nsoil  ,                         & !in
                            isnow  ,stc    ,snice  ,snliq  ,dzsnso )   !inout
   end if

!set empty snow layers to zero

   do iz = -nsnow+1, isnow
        snice(iz) = 0.
        snliq(iz) = 0.
        stc(iz)   = 0.
        dzsnso(iz)= 0.
        zsnso(iz) = 0.
   enddo

   call  snowh2o_glacier (nsnow  ,nsoil  ,dt     ,qsnfro ,qsnsub , & !in 
                          qrain  ,                                 & !in
                          isnow  ,dzsnso ,snowh  ,sneqv  ,snice  , & !inout
                          snliq  ,sh2o   ,sice   ,stc    ,         & !inout
			  ponding1       ,ponding2       ,         & !inout
                          qsnbot )                                   !out

!to obtain equilibrium state of snow in glacier region
       
   if(sneqv > 2000.) then   ! 2000 mm -> maximum water depth
      bdsnow      = snice(0) / dzsnso(0)
      snoflow     = (sneqv - 2000.)
      snice(0)    = snice(0)  - snoflow 
      dzsnso(0)   = dzsnso(0) - snoflow/bdsnow
      snoflow     = snoflow / dt
   end if

! sum up snow mass for layered snow

   if(isnow /= 0) then
       sneqv = 0.
       do iz = isnow+1,0
             sneqv = sneqv + snice(iz) + snliq(iz)
       enddo
   end if

! reset zsnso and layer thinkness dzsnso

   do iz = isnow+1, 0
        dzsnso(iz) = -dzsnso(iz)
   end do

   dzsnso(1) = zsoil(1)
   do iz = 2,nsoil
        dzsnso(iz) = (zsoil(iz) - zsoil(iz-1))
   end do

   zsnso(isnow+1) = dzsnso(isnow+1)
   do iz = isnow+2 ,nsoil
       zsnso(iz) = zsnso(iz-1) + dzsnso(iz)
   enddo

   do iz = isnow+1 ,nsoil
       dzsnso(iz) = -dzsnso(iz)
   end do

  end subroutine snowwater_glacier
! ==================================================================================================
!>\ingroup NoahMP_LSM
  subroutine snowfall_glacier (nsoil  ,nsnow  ,dt     ,qsnow  ,snowhin , & !in
                               sfctmp ,                                  & !in
                               isnow  ,snowh  ,dzsnso ,stc    ,snice   , & !inout
                               snliq  ,sneqv  )                            !inout
! ----------------------------------------------------------------------
! snow depth and density to account for the new snowfall.
! new values of snow depth & density returned.
! ----------------------------------------------------------------------
    implicit none
! ----------------------------------------------------------------------
! input

  integer,                            intent(in) :: nsoil  !no. of soil layers
  integer,                            intent(in) :: nsnow  !maximum no. of snow layers
  real,                               intent(in) :: dt     !main time step (s)
  real,                               intent(in) :: qsnow  !snow at ground srf (mm/s) [+]
  real,                               intent(in) :: snowhin!snow depth increasing rate (m/s)
  real,                               intent(in) :: sfctmp !surface air temperature [k]

! input and output

  integer,                         intent(inout) :: isnow  !actual no. of snow layers
  real,                            intent(inout) :: snowh  !snow depth [m]
  real,                            intent(inout) :: sneqv  !swow water equivalent [m]
  real, dimension(-nsnow+1:nsoil), intent(inout) :: dzsnso !thickness of snow/soil layers (m)
  real, dimension(-nsnow+1:nsoil), intent(inout) :: stc    !snow layer temperature [k]
  real, dimension(-nsnow+1:    0), intent(inout) :: snice  !snow layer ice [mm]
  real, dimension(-nsnow+1:    0), intent(inout) :: snliq  !snow layer liquid water [mm]

! local

  integer :: newnode            ! 0-no new layers, 1-creating new layers
! ----------------------------------------------------------------------
    newnode  = 0

! shallow snow / no layer

    if(isnow == 0 .and. qsnow > 0.)  then
      snowh = snowh + snowhin * dt
      sneqv = sneqv + qsnow * dt
    end if

! creating a new layer
 
    if(isnow == 0  .and. qsnow>0. .and. snowh >= 0.05) then
      isnow    = -1
      newnode  =  1
      dzsnso(0)= snowh
      snowh    = 0.
      stc(0)   = min(273.16, sfctmp)   ! temporary setup
      snice(0) = sneqv
      snliq(0) = 0.
    end if

! snow with layers

    if(isnow <  0 .and. newnode == 0 .and. qsnow > 0.) then
         snice(isnow+1)  = snice(isnow+1)   + qsnow   * dt
         dzsnso(isnow+1) = dzsnso(isnow+1)  + snowhin * dt
    endif

! ----------------------------------------------------------------------
  end subroutine snowfall_glacier
! ==================================================================================================
! ----------------------------------------------------------------------
!>\ingroup NoahMP_LSM
  subroutine compact_glacier (nsnow  ,nsoil  ,dt     ,stc    ,snice , & !in
                              snliq  ,imelt  ,ficeold,                & !in
                              isnow  ,dzsnso )                          !inout
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
  implicit none
! ----------------------------------------------------------------------
! input
   integer,                         intent(in)    :: nsoil  !no. of soil layers [ =4]
   integer,                         intent(in)    :: nsnow  !maximum no. of snow layers [ =3]
   integer, dimension(-nsnow+1:0) , intent(in)    :: imelt  !melting state index [0-no melt;1-melt]
   real,                            intent(in)    :: dt     !time step (sec)
   real, dimension(-nsnow+1:nsoil), intent(in)    :: stc    !snow layer temperature [k]
   real, dimension(-nsnow+1:    0), intent(in)    :: snice  !snow layer ice [mm]
   real, dimension(-nsnow+1:    0), intent(in)    :: snliq  !snow layer liquid water [mm]
   real, dimension(-nsnow+1:    0), intent(in)    :: ficeold!ice fraction at last timestep

! input and output
   integer,                         intent(inout) :: isnow  ! actual no. of snow layers
   real, dimension(-nsnow+1:nsoil), intent(inout) :: dzsnso ! snow layer thickness [m]

! local
   real, parameter     :: c2 = 21.e-3   ![m3/kg] ! default 21.e-3
   real, parameter     :: c3 = 2.5e-6   ![1/s]  
   real, parameter     :: c4 = 0.04     ![1/k]
   real, parameter     :: c5 = 2.0      !
   real, parameter     :: dm = 100.0    !upper limit on destructive metamorphism compaction [kg/m3]
   real, parameter     :: eta0 = 0.8e+6 !viscosity coefficient [kg-s/m2] 
                                        !according to anderson, it is between 0.52e6~1.38e6
   real :: burden !pressure of overlying snow [kg/m2]
   real :: ddz1   !rate of settling of snow pack due to destructive metamorphism.
   real :: ddz2   !rate of compaction of snow pack due to overburden.
   real :: ddz3   !rate of compaction of snow pack due to melt [1/s]
   real :: dexpf  !expf=exp(-c4*(273.15-stc)).
   real :: td     !stc - tfrz [k]
   real :: pdzdtc !nodal rate of change in fractional-thickness due to compaction [fraction/s]
   real :: void   !void (1 - snice - snliq)
   real :: wx     !water mass (ice + liquid) [kg/m2]
   real :: bi     !partial density of ice [kg/m3]
   real, dimension(-nsnow+1:0) :: fice   !fraction of ice at current time step

   integer  :: j

! ----------------------------------------------------------------------
    burden = 0.0

    do j = isnow+1, 0

        wx      = snice(j) + snliq(j)
        fice(j) = snice(j) / wx
        void    = 1. - (snice(j)/denice + snliq(j)/denh2o) / dzsnso(j)

        ! allow compaction only for non-saturated node and higher ice lens node.
        if (void > 0.001 .and. snice(j) > 0.1) then
           bi = snice(j) / dzsnso(j)
           td = max(0.,tfrz-stc(j))
           dexpf = exp(-c4*td)

           ! settling as a result of destructive metamorphism

           ddz1 = -c3*dexpf

           if (bi > dm) ddz1 = ddz1*exp(-46.0e-3*(bi-dm))

           ! liquid water term

           if (snliq(j) > 0.01*dzsnso(j)) ddz1=ddz1*c5

           ! compaction due to overburden

           ddz2 = -(burden+0.5*wx)*exp(-0.08*td-c2*bi)/eta0 ! 0.5*wx -> self-burden

           ! compaction occurring during melt

           if (imelt(j) == 1) then
              ddz3 = max(0.,(ficeold(j) - fice(j))/max(1.e-6,ficeold(j)))
              ddz3 = - ddz3/dt           ! sometimes too large
           else
              ddz3 = 0.
           end if

           ! time rate of fractional change in dz (units of s-1)

           pdzdtc = (ddz1 + ddz2 + ddz3)*dt
           pdzdtc = max(-0.5,pdzdtc)

           ! the change in dz due to compaction

           dzsnso(j) = dzsnso(j)*(1.+pdzdtc)
        end if

        ! pressure of overlying snow

        burden = burden + wx

    end do

  end subroutine compact_glacier
! ==================================================================================================
!>\ingroup NoahMP_LSM
  subroutine combine_glacier (nsnow  ,nsoil  ,                         & !in
                              isnow  ,sh2o   ,stc    ,snice  ,snliq  , & !inout
                              dzsnso ,sice   ,snowh  ,sneqv  ,         & !inout
                              ponding1       ,ponding2)                  !inout
! ----------------------------------------------------------------------
    implicit none
! ----------------------------------------------------------------------
! input

    integer, intent(in)     :: nsnow                        !maximum no. of snow layers
    integer, intent(in)     :: nsoil                        !no. of soil layers

! input and output

    integer,                         intent(inout) :: isnow !actual no. of snow layers
    real, dimension(       1:nsoil), intent(inout) :: sh2o  !soil liquid moisture (m3/m3)
    real, dimension(       1:nsoil), intent(inout) :: sice  !soil ice moisture (m3/m3)
    real, dimension(-nsnow+1:nsoil), intent(inout) :: stc   !snow layer temperature [k]
    real, dimension(-nsnow+1:    0), intent(inout) :: snice !snow layer ice [mm]
    real, dimension(-nsnow+1:    0), intent(inout) :: snliq !snow layer liquid water [mm]
    real, dimension(-nsnow+1:nsoil), intent(inout) :: dzsnso!snow layer depth [m]
    real,                            intent(inout) :: sneqv !snow water equivalent [m]
    real,                            intent(inout) :: snowh !snow depth [m]
    real,                            intent(inout) :: ponding1
    real,                            intent(inout) :: ponding2

! local variables:

    integer :: i,j,k,l               ! node indices
    integer :: isnow_old             ! number of top snow layer
    integer :: mssi                  ! node index
    integer :: neibor                ! adjacent node selected for combination
    real    :: zwice                 ! total ice mass in snow
    real    :: zwliq                 ! total liquid water in snow
    real    :: dzmin(3)              ! minimum of top snow layer
    data dzmin /0.045, 0.05, 0.2/
!    data dzmin /0.025, 0.025, 0.1/  ! mb: change limit
!-----------------------------------------------------------------------

       isnow_old = isnow

       do j = isnow_old+1,0
          if (snice(j) <= .1) then
             if(j /= 0) then
                snliq(j+1) = snliq(j+1) + snliq(j)
                snice(j+1) = snice(j+1) + snice(j)
             else
               if (isnow_old < -1) then
                snliq(j-1) = snliq(j-1) + snliq(j)
                snice(j-1) = snice(j-1) + snice(j)
               else
                ponding1 = ponding1 + snliq(j)       ! isnow will get set to zero below
                sneqv = snice(j)                     ! ponding will get added to ponding from
                snowh = dzsnso(j)                    ! phasechange which should be zero here
                snliq(j) = 0.0                       ! because there it was only calculated
                snice(j) = 0.0                       ! for thin snow
                dzsnso(j) = 0.0
               endif
!                sh2o(1) = sh2o(1)+snliq(j)/(dzsnso(1)*1000.)
!                sice(1) = sice(1)+snice(j)/(dzsnso(1)*1000.)
             endif

             ! shift all elements above this down by one.
             if (j > isnow+1 .and. isnow < -1) then
                do i = j, isnow+2, -1
                   stc(i)   = stc(i-1)
                   snliq(i) = snliq(i-1)
                   snice(i) = snice(i-1)
                   dzsnso(i)= dzsnso(i-1)
                end do
             end if
             isnow = isnow + 1
          end if
       end do

! to conserve water in case of too large surface sublimation

       if(sice(1) < 0.) then
          sh2o(1) = sh2o(1) + sice(1)
          sice(1) = 0.
       end if

       if(isnow ==0) return   ! mb: get out if no longer multi-layer

       sneqv  = 0.
       snowh  = 0.
       zwice  = 0.
       zwliq  = 0.

       do j = isnow+1,0
             sneqv = sneqv + snice(j) + snliq(j)
             snowh = snowh + dzsnso(j)
             zwice = zwice + snice(j)
             zwliq = zwliq + snliq(j)
       end do

! check the snow depth - all snow gone
! the liquid water assumes ponding on soil surface.

!       if (snowh < 0.025 .and. isnow < 0 ) then ! mb: change limit
       if (snowh < 0.05 .and. isnow < 0 ) then
          isnow  = 0
          sneqv = zwice
          ponding2 = ponding2 + zwliq           ! limit of isnow < 0 means input ponding
          if(sneqv <= 0.) snowh = 0.            ! should be zero; see above
       end if

!       if (snowh < 0.05 ) then
!          isnow  = 0
!          sneqv = zwice
!          sh2o(1) = sh2o(1) + zwliq / (dzsnso(1) * 1000.)
!          if(sneqv <= 0.) snowh = 0.
!       end if

! check the snow depth - snow layers combined

       if (isnow < -1) then

          isnow_old = isnow
          mssi     = 1

          do i = isnow_old+1,0
             if (dzsnso(i) < dzmin(mssi)) then

                if (i == isnow+1) then
                   neibor = i + 1
                else if (i == 0) then
                   neibor = i - 1
                else
                   neibor = i + 1
                   if ((dzsnso(i-1)+dzsnso(i)) < (dzsnso(i+1)+dzsnso(i))) neibor = i-1
                end if

                ! node l and j are combined and stored as node j.
                if (neibor > i) then
                   j = neibor
                   l = i
                else
                   j = i
                   l = neibor
                end if

                call combo_glacier (dzsnso(j), snliq(j), snice(j), &
                   stc(j), dzsnso(l), snliq(l), snice(l), stc(l) )

                ! now shift all elements above this down one.
                if (j-1 > isnow+1) then
                   do k = j-1, isnow+2, -1
                      stc(k)   = stc(k-1)
                      snice(k) = snice(k-1)
                      snliq(k) = snliq(k-1)
                      dzsnso(k) = dzsnso(k-1)
                   end do
                end if

                ! decrease the number of snow layers
                isnow = isnow + 1
                if (isnow >= -1) exit
             else

                ! the layer thickness is greater than the prescribed minimum value
                mssi = mssi + 1

             end if
          end do

       end if

  end subroutine combine_glacier
! ==================================================================================================

! ----------------------------------------------------------------------
!>\ingroup NoahMP_LSM
  subroutine combo_glacier(dz,  wliq,  wice, t, dz2, wliq2, wice2, t2)
! ----------------------------------------------------------------------
    implicit none
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------s
! input

    real, intent(in)    :: dz2   !nodal thickness of 2 elements being combined [m]
    real, intent(in)    :: wliq2 !liquid water of element 2 [kg/m2]
    real, intent(in)    :: wice2 !ice of element 2 [kg/m2]
    real, intent(in)    :: t2    !nodal temperature of element 2 [k]
    real, intent(inout) :: dz    !nodal thickness of 1 elements being combined [m]
    real, intent(inout) :: wliq  !liquid water of element 1
    real, intent(inout) :: wice  !ice of element 1 [kg/m2]
    real, intent(inout) :: t     !node temperature of element 1 [k]

! local 

    real                :: dzc   !total thickness of nodes 1 and 2 (dzc=dz+dz2).
    real                :: wliqc !combined liquid water [kg/m2]
    real                :: wicec !combined ice [kg/m2]
    real                :: tc    !combined node temperature [k]
    real                :: h     !enthalpy of element 1 [j/m2]
    real                :: h2    !enthalpy of element 2 [j/m2]
    real                :: hc    !temporary

!-----------------------------------------------------------------------

    dzc = dz+dz2
    wicec = (wice+wice2)
    wliqc = (wliq+wliq2)
    h = (cice*wice+cwat*wliq) * (t-tfrz)+hfus*wliq
    h2= (cice*wice2+cwat*wliq2) * (t2-tfrz)+hfus*wliq2

    hc = h + h2
    if(hc < 0.)then
       tc = tfrz + hc/(cice*wicec + cwat*wliqc)
    else if (hc.le.hfus*wliqc) then
       tc = tfrz
    else
       tc = tfrz + (hc - hfus*wliqc) / (cice*wicec + cwat*wliqc)
    end if

    dz = dzc
    wice = wicec
    wliq = wliqc
    t = tc

  end subroutine combo_glacier
! ==================================================================================================
!>\ingroup NoahMP_LSM
  subroutine divide_glacier (nsnow  ,nsoil  ,                         & !in
                             isnow  ,stc    ,snice  ,snliq  ,dzsnso  )  !inout
! ----------------------------------------------------------------------
    implicit none
! ----------------------------------------------------------------------
! input

    integer, intent(in)                            :: nsnow !maximum no. of snow layers [ =3]
    integer, intent(in)                            :: nsoil !no. of soil layers [ =4]

! input and output

    integer                        , intent(inout) :: isnow !actual no. of snow layers 
    real, dimension(-nsnow+1:nsoil), intent(inout) :: stc   !snow layer temperature [k]
    real, dimension(-nsnow+1:    0), intent(inout) :: snice !snow layer ice [mm]
    real, dimension(-nsnow+1:    0), intent(inout) :: snliq !snow layer liquid water [mm]
    real, dimension(-nsnow+1:nsoil), intent(inout) :: dzsnso!snow layer depth [m]

! local variables:

    integer                                        :: j     !indices
    integer                                        :: msno  !number of layer (top) to msno (bot)
    real                                           :: drr   !thickness of the combined [m]
    real, dimension(       1:nsnow)                :: dz    !snow layer thickness [m]
    real, dimension(       1:nsnow)                :: swice !partial volume of ice [m3/m3]
    real, dimension(       1:nsnow)                :: swliq !partial volume of liquid water [m3/m3]
    real, dimension(       1:nsnow)                :: tsno  !node temperature [k]
    real                                           :: zwice !temporary
    real                                           :: zwliq !temporary
    real                                           :: propor!temporary
    real                                           :: dtdz  !temporary
! ----------------------------------------------------------------------

    do j = 1,nsnow
          if (j <= abs(isnow)) then
             dz(j)    = dzsnso(j+isnow)
             swice(j) = snice(j+isnow)
             swliq(j) = snliq(j+isnow)
             tsno(j)  = stc(j+isnow)
          end if
    end do

       msno = abs(isnow)

       if (msno == 1) then
          ! specify a new snow layer
          if (dz(1) > 0.05) then
             msno = 2
             dz(1)    = dz(1)/2.
             swice(1) = swice(1)/2.
             swliq(1) = swliq(1)/2.
             dz(2)    = dz(1)
             swice(2) = swice(1)
             swliq(2) = swliq(1)
             tsno(2)  = tsno(1)
          end if
       end if

       if (msno > 1) then
          if (dz(1) > 0.05) then
             drr      = dz(1) - 0.05
             propor   = drr/dz(1)
             zwice    = propor*swice(1)
             zwliq    = propor*swliq(1)
             propor   = 0.05/dz(1)
             swice(1) = propor*swice(1)
             swliq(1) = propor*swliq(1)
             dz(1)    = 0.05

             call combo_glacier (dz(2), swliq(2), swice(2), tsno(2), drr, &
                  zwliq, zwice, tsno(1))

             ! subdivide a new layer
!             if (msno <= 2 .and. dz(2) > 0.20) then  ! mb: change limit
             if (msno <= 2 .and. dz(2) > 0.10) then
                msno = 3
                dtdz = (tsno(1) - tsno(2))/((dz(1)+dz(2))/2.)
                dz(2)    = dz(2)/2.
                swice(2) = swice(2)/2.
                swliq(2) = swliq(2)/2.
                dz(3)    = dz(2)
                swice(3) = swice(2)
                swliq(3) = swliq(2)
                tsno(3) = tsno(2) - dtdz*dz(2)/2.
                if (tsno(3) >= tfrz) then
                   tsno(3)  = tsno(2)
                else
                   tsno(2) = tsno(2) + dtdz*dz(2)/2.
                endif

             end if
          end if
       end if

       if (msno > 2) then
          if (dz(2) > 0.2) then
             drr = dz(2) - 0.2
             propor   = drr/dz(2)
             zwice    = propor*swice(2)
             zwliq    = propor*swliq(2)
             propor   = 0.2/dz(2)
             swice(2) = propor*swice(2)
             swliq(2) = propor*swliq(2)
             dz(2)    = 0.2
             call combo_glacier (dz(3), swliq(3), swice(3), tsno(3), drr, &
                  zwliq, zwice, tsno(2))
          end if
       end if

       isnow = -msno

    do j = isnow+1,0
             dzsnso(j) = dz(j-isnow)
             snice(j) = swice(j-isnow)
             snliq(j) = swliq(j-isnow)
             stc(j)   = tsno(j-isnow)
    end do


!    do j = isnow+1,nsoil
!    write(*,'(i5,7f10.3)') j, dzsnso(j), snice(j), snliq(j),stc(j)
!    end do

  end subroutine divide_glacier
! ==================================================================================================
!>\ingroup NoahMP_LSM
  subroutine snowh2o_glacier (nsnow  ,nsoil  ,dt     ,qsnfro ,qsnsub , & !in 
                              qrain  ,                                 & !in
                              isnow  ,dzsnso ,snowh  ,sneqv  ,snice  , & !inout
                              snliq  ,sh2o   ,sice   ,stc    ,         & !inout
                              ponding1       ,ponding2       ,         & !inout
                              qsnbot )                                   !out
! ----------------------------------------------------------------------
! renew the mass of ice lens (snice) and liquid (snliq) of the
! surface snow layer resulting from sublimation (frost) / evaporation (dew)
! ----------------------------------------------------------------------
   implicit none
! ----------------------------------------------------------------------
! input

   integer,                         intent(in)    :: nsnow  !maximum no. of snow layers[=3]
   integer,                         intent(in)    :: nsoil  !no. of soil layers[=4]
   real,                            intent(in)    :: dt     !time step
   real,                            intent(in)    :: qsnfro !snow surface frost rate[mm/s]
   real,                            intent(in)    :: qsnsub !snow surface sublimation rate[mm/s]
   real,                            intent(in)    :: qrain  !snow surface rain rate[mm/s]

! output

   real,                            intent(out)   :: qsnbot !melting water out of snow bottom [mm/s]

! input and output

   integer,                         intent(inout) :: isnow  !actual no. of snow layers
   real, dimension(-nsnow+1:nsoil), intent(inout) :: dzsnso ! snow layer depth [m]
   real,                            intent(inout) :: snowh  !snow height [m]
   real,                            intent(inout) :: sneqv  !snow water eqv. [mm]
   real, dimension(-nsnow+1:0),     intent(inout) :: snice  !snow layer ice [mm]
   real, dimension(-nsnow+1:0),     intent(inout) :: snliq  !snow layer liquid water [mm]
   real, dimension(       1:nsoil), intent(inout) :: sh2o   !soil liquid moisture (m3/m3)
   real, dimension(       1:nsoil), intent(inout) :: sice   !soil ice moisture (m3/m3)
   real, dimension(-nsnow+1:nsoil), intent(inout) :: stc    !snow layer temperature [k]
   real,                            intent(inout) :: ponding1
   real,                            intent(inout) :: ponding2

! local variables:

   integer                     :: j         !do loop/array indices
   real                        :: qin       !water flow into the element (mm/s)
   real                        :: qout      !water flow out of the element (mm/s)
   real                        :: wgdif     !ice mass after minus sublimation
   real, dimension(-nsnow+1:0) :: vol_liq   !partial volume of liquid water in layer
   real, dimension(-nsnow+1:0) :: vol_ice   !partial volume of ice lens in layer
   real, dimension(-nsnow+1:0) :: epore     !effective porosity = porosity - vol_ice
   real :: propor, temp
! ----------------------------------------------------------------------

!for the case when sneqv becomes '0' after 'combine'

   if(sneqv == 0.) then
      sice(1) =  sice(1) + (qsnfro-qsnsub)*dt/(dzsnso(1)*1000.)
   end if

! for shallow snow without a layer
! snow surface sublimation may be larger than existing snow mass. to conserve water,
! excessive sublimation is used to reduce soil water. smaller time steps would tend 
! to aviod this problem.

   if(isnow == 0 .and. sneqv > 0.) then
      temp   = sneqv
      sneqv  = sneqv - qsnsub*dt + qsnfro*dt
      propor = sneqv/temp
      snowh  = max(0.,propor * snowh)

      if(sneqv < 0.) then
         sice(1) = sice(1) + sneqv/(dzsnso(1)*1000.)
         sneqv   = 0.
         snowh   = 0.
      end if
      if(sice(1) < 0.) then
         sh2o(1) = sh2o(1) + sice(1)
         sice(1) = 0.
      end if
   end if

   if(snowh <= 1.e-8 .or. sneqv <= 1.e-6) then
     snowh = 0.0
     sneqv = 0.0
   end if

! for deep snow

   if ( isnow < 0 ) then !kwm added this if statement to prevent out-of-bounds array references

      wgdif = snice(isnow+1) - qsnsub*dt + qsnfro*dt
      snice(isnow+1) = wgdif
      if (wgdif < 1.e-6 .and. isnow <0) then
         call  combine_glacier (nsnow  ,nsoil  ,                         & !in
                                isnow  ,sh2o   ,stc    ,snice  ,snliq  , & !inout
                                dzsnso ,sice   ,snowh  ,sneqv  ,         & !inout
                               ponding1, ponding2 )                        !inout
      endif
      !kwm:  subroutine combine can change isnow to make it 0 again?
      if ( isnow < 0 ) then !kwm added this if statement to prevent out-of-bounds array references
         snliq(isnow+1) = snliq(isnow+1) + qrain * dt
         snliq(isnow+1) = max(0., snliq(isnow+1))
      endif
      
   endif !kwm  -- can the endif be moved toward the end of the subroutine (just set qsnbot=0)?

! porosity and partial volume

   !kwm looks to me like loop index / if test can be simplified.

   do j = -nsnow+1, 0
      if (j >= isnow+1) then
         vol_ice(j)      = min(1., snice(j)/(dzsnso(j)*denice))
         epore(j)        = 1. - vol_ice(j)
         vol_liq(j)      = min(epore(j),snliq(j)/(dzsnso(j)*denh2o))
      end if
   end do

   qin = 0.
   qout = 0.

   !kwm looks to me like loop index / if test can be simplified.

   do j = -nsnow+1, 0
      if (j >= isnow+1) then
         snliq(j) = snliq(j) + qin
         if (j <= -1) then
            if (epore(j) < 0.05 .or. epore(j+1) < 0.05) then
               qout = 0.
            else
               qout = max(0.,(vol_liq(j)-ssi*epore(j))*dzsnso(j))
               qout = min(qout,(1.-vol_ice(j+1)-vol_liq(j+1))*dzsnso(j+1))
            end if
         else
            qout = max(0.,(vol_liq(j) - ssi*epore(j))*dzsnso(j))
         end if
         qout = qout*1000.
         snliq(j) = snliq(j) - qout
         qin = qout
      end if
   end do

! liquid water from snow bottom to soil

   qsnbot = qout / dt           ! mm/s

  end subroutine snowh2o_glacier
! ********************* end of water subroutines ******************************************
! ==================================================================================================
!>\ingroup NoahMP_LSM
  subroutine error_glacier (iloc   ,jloc   ,swdown ,fsa    ,fsr    ,fira   , &
                            fsh    ,fgev   ,ssoil  ,sag    ,prcp   ,edir   , &
#ifdef CCPP
                           runsrf ,runsub ,sneqv  ,dt     ,beg_wb, errmsg, errflg )
#else 
                           runsrf ,runsub ,sneqv  ,dt     ,beg_wb )
#endif
! --------------------------------------------------------------------------------------------------
! check surface energy balance and water balance
! --------------------------------------------------------------------------------------------------
  implicit none
! --------------------------------------------------------------------------------------------------
! inputs
  integer                        , intent(in) :: iloc   !grid index
  integer                        , intent(in) :: jloc   !grid index
  real                           , intent(in) :: swdown !downward solar filtered by sun angle [w/m2]
  real                           , intent(in) :: fsa    !total absorbed solar radiation (w/m2)
  real                           , intent(in) :: fsr    !total reflected solar radiation (w/m2)
  real                           , intent(in) :: fira   !total net longwave rad (w/m2)  [+ to atm]
  real                           , intent(in) :: fsh    !total sensible heat (w/m2)     [+ to atm]
  real                           , intent(in) :: fgev   !ground evaporation heat (w/m2) [+ to atm]
  real                           , intent(in) :: ssoil  !ground heat flux (w/m2)        [+ to soil]
  real                           , intent(in) :: sag

  real                           , intent(in) :: prcp   !precipitation rate (kg m-2 s-1)
  real                           , intent(in) :: edir   !soil surface evaporation rate[mm/s]
  real                           , intent(in) :: runsrf !surface runoff [mm/s] 
  real                           , intent(in) :: runsub !baseflow (saturation excess) [mm/s]
  real                           , intent(in) :: sneqv  !snow water eqv. [mm]
  real                           , intent(in) :: dt     !time step [sec]
  real                           , intent(in) :: beg_wb !water storage at begin of a timesetp [mm]

#ifdef CCPP  
  character(len=*)               , intent(inout) :: errmsg
  integer                        , intent(inout) :: errflg
#endif

  real                                        :: end_wb !water storage at end of a timestep [mm]
  real                                        :: errwat !error in water balance [mm/timestep]
  real                                        :: erreng !error in surface energy balance [w/m2]
  real                                        :: errsw  !error in shortwave radiation balance [w/m2]
  character(len=256)                          :: message
! --------------------------------------------------------------------------------------------------
   errsw   = swdown - (fsa + fsr)
   if (errsw > 0.01) then            ! w/m2
     write(*,*) "sag    =",sag
     write(*,*) "fsa    =",fsa
     write(*,*) "fsr    =",fsr
     write(message,*) 'errsw =',errsw
#ifdef CCPP
     errflg = 1
     errmsg = trim(message)//NEW_LINE('A')//"radiation budget problem in noahmp glacier"
     return 
#else
     call wrf_message(trim(message))
     call wrf_error_fatal("radiation budget problem in noahmp glacier")
#endif
   end if

   erreng = sag-(fira+fsh+fgev+ssoil)
   if(erreng > 0.01) then
      write(message,*) 'erreng =',erreng
#ifdef CCPP
      errmsg = trim(message)    
#else
      call wrf_message(trim(message))
#endif  
      write(message,'(i6,1x,i6,1x,5f10.4)')iloc,jloc,sag,fira,fsh,fgev,ssoil
#ifdef CCPP
     errflg = 1
     errmsg = trim(errmsg)//NEW_LINE('A')//"energy budget problem in noahmp glacier"
     return 
#else
     call wrf_message(trim(message))
     call wrf_error_fatal("energy budget problem in noahmp glacier")
#endif
   end if

   end_wb = sneqv
   errwat = end_wb-beg_wb-(prcp-edir-runsrf-runsub)*dt


 end subroutine error_glacier
! ==================================================================================================

!>\ingroup NoahMP_LSM
  subroutine noahmp_options_glacier(idveg     ,iopt_crs  ,iopt_btr  ,iopt_run  ,iopt_sfc  ,iopt_frz , & 
                             iopt_inf  ,iopt_rad  ,iopt_alb  ,iopt_snf  ,iopt_tbot, iopt_stc )

  implicit none

  integer,  intent(in) :: idveg     !dynamic vegetation (1 -> off ; 2 -> on) with opt_crs = 1
  integer,  intent(in) :: iopt_crs  !canopy stomatal resistance (1-> ball-berry; 2->jarvis)
  integer,  intent(in) :: iopt_btr  !soil moisture factor for stomatal resistance (1-> noah; 2-> clm; 3-> ssib)
  integer,  intent(in) :: iopt_run  !runoff and groundwater (1->simgm; 2->simtop; 3->schaake96; 4->bats)
  integer,  intent(in) :: iopt_sfc  !surface layer drag coeff (ch & cm) (1->m-o; 2->chen97)
  integer,  intent(in) :: iopt_frz  !supercooled liquid water (1-> ny06; 2->koren99)
  integer,  intent(in) :: iopt_inf  !frozen soil permeability (1-> ny06; 2->koren99)
  integer,  intent(in) :: iopt_rad  !radiation transfer (1->gap=f(3d,cosz); 2->gap=0; 3->gap=1-fveg)
  integer,  intent(in) :: iopt_alb  !snow surface albedo (1->bats; 2->class)
  integer,  intent(in) :: iopt_snf  !rainfall & snowfall (1-jordan91; 2->bats; 3->noah)
  integer,  intent(in) :: iopt_tbot !lower boundary of soil temperature (1->zero-flux; 2->noah)

  integer,  intent(in) :: iopt_stc  !snow/soil temperature time scheme (only layer 1)
                                    ! 1 -> semi-implicit; 2 -> full implicit (original noah)

! -------------------------------------------------------------------------------------------------

  dveg = idveg
  
  opt_crs  = iopt_crs  
  opt_btr  = iopt_btr  
  opt_run  = iopt_run  
  opt_sfc  = iopt_sfc  
  opt_frz  = iopt_frz  
  opt_inf  = iopt_inf  
  opt_rad  = iopt_rad  
  opt_alb  = iopt_alb  
  opt_snf  = iopt_snf  
  opt_tbot = iopt_tbot 
  opt_stc  = iopt_stc
  
  end subroutine noahmp_options_glacier
 
end module noahmp_glacier_routines
! ==================================================================================================

module module_sf_noahmp_glacier

  use noahmp_glacier_routines
  use noahmp_glacier_globals

end module module_sf_noahmp_glacier


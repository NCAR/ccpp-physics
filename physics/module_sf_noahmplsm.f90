!>  \file module_sf_noahmplsm.f90
!!  This file contains the NoahMP land surface model.

!>\ingroup NoahMP_LSM
module module_sf_noahmplsm
#ifndef CCPP  
  use  module_wrf_utl
#endif

  implicit none

  public  :: noahmp_options
  public  :: noahmp_sflx

  private :: atm
  private :: phenology
  private :: precip_heat
  private :: energy
  private ::       thermoprop
  private ::               csnow
  private ::               tdfcnd
  private ::       radiation
  private ::               albedo
  private ::                         snow_age
  private ::                         snowalb_bats  
  private ::                         snowalb_class
  private ::                         groundalb
  private ::                         twostream
  private ::               surrad
  private ::       vege_flux
  private ::               sfcdif1                  
  private ::               sfcdif2                
  private ::               stomata                  
  private ::               canres                  
  private ::               esat
  private ::               ragrb
  private ::       bare_flux
  private ::       tsnosoi
  private ::               hrt
  private ::               hstep   
  private ::                         rosr12
  private ::       phasechange
  private ::               frh2o           

  private :: water
  private ::       canwater
  private ::       snowwater
  private ::               snowfall
  private ::               combine
  private ::               divide
  private ::                         combo
  private ::               compact
  private ::               snowh2o
  private ::       soilwater
  private ::               zwteq
  private ::               infil
  private ::               srt
  private ::                         wdfcnd1        
  private ::                         wdfcnd2       
  private ::               sstep
  private ::       groundwater
  private ::       shallowwatertable

  private :: carbon
  private ::       co2flux
!  private ::       bvocflux
!  private ::       ch4flux

  private :: error

! =====================================options for different schemes================================
! **recommended

  integer :: dveg     ! options for dynamic vegetation: 
                      !   1 -> off (use table lai; use fveg = shdfac from input)
                      !   2 -> on  (together with opt_crs = 1)
                      !   3 -> off (use table lai; calculate fveg)
                      ! **4 -> off (use table lai; use maximum vegetation fraction)
                      ! **5 -> on  (use maximum vegetation fraction)

  integer :: opt_crs  ! options for canopy stomatal resistance
                      ! **1 -> ball-berry
		      !   2 -> jarvis

  integer :: opt_btr  ! options for soil moisture factor for stomatal resistance
                      ! **1 -> noah (soil moisture) 
                      !   2 -> clm  (matric potential)
                      !   3 -> ssib (matric potential)

  integer :: opt_run  ! options for runoff and groundwater
                      ! **1 -> topmodel with groundwater (niu et al. 2007 jgr) ;
                      !   2 -> topmodel with an equilibrium water table (niu et al. 2005 jgr) ;
                      !   3 -> original surface and subsurface runoff (free drainage)
                      !   4 -> bats surface and subsurface runoff (free drainage)
                      !   5 -> miguez-macho&fan groundwater scheme (miguez-macho et al. 2007 jgr; fan et al. 2007 jgr)
		      !          (needs further testing for public use)

  integer :: opt_sfc  ! options for surface layer drag coeff (ch & cm)
                      ! **1 -> m-o
		      ! **2 -> original noah (chen97)
		      ! **3 -> myj consistent; 4->ysu consistent. mb: removed in v3.7 for further testing

  integer :: opt_frz  ! options for supercooled liquid water (or ice fraction)
                      ! **1 -> no iteration (niu and yang, 2006 jhm)
		      !   2 -> koren's iteration 

  integer :: opt_inf  ! options for frozen soil permeability
                      ! **1 -> linear effects, more permeable (niu and yang, 2006, jhm)
                      !   2 -> nonlinear effects, less permeable (old)

  integer :: opt_rad  ! options for radiation transfer
                      !   1 -> modified two-stream (gap = f(solar angle, 3d structure ...)<1-fveg)
                      !   2 -> two-stream applied to grid-cell (gap = 0)
                      ! **3 -> two-stream applied to vegetated fraction (gap=1-fveg)

  integer :: opt_alb  ! options for ground snow surface albedo
                      !   1 -> bats
		      ! **2 -> class

  integer :: opt_snf  ! options for partitioning  precipitation into rainfall & snowfall
                      ! **1 -> jordan (1991)
		      !   2 -> bats: when sfctmp<tfrz+2.2 
		      !   3 -> sfctmp < tfrz
		      !   4 -> use wrf microphysics output

  integer :: opt_tbot ! options for lower boundary condition of soil temperature
                      !   1 -> zero heat flux from bottom (zbot and tbot not used)
                      ! **2 -> tbot at zbot (8m) read from a file (original noah)

  integer :: opt_stc  ! options for snow/soil temperature time scheme (only layer 1)
                      ! **1 -> semi-implicit; flux top boundary condition
		      !   2 -> full implicit (original noah); temperature top boundary condition
                      !   3 -> same as 1, but fsno for ts calculation (generally improves snow; v3.7)

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
  real, parameter :: tkair  = 0.023     !thermal conductivity of air (w/m/k) (not used mb: 20140718)
  real, parameter :: rair   = 287.04    !gas constant for dry air (j/kg/k)
  real, parameter :: rw     = 461.269   !gas constant for  water vapor (j/kg/k)
  real, parameter :: denh2o = 1000.     !density of water (kg/m3)
  real, parameter :: denice = 917.      !density of ice (kg/m3)

  integer, private, parameter :: mband = 2

  type noahmp_parameters ! define a noahmp parameters type

!------------------------------------------------------------------------------------------!
! from the veg section of mptable.tbl
!------------------------------------------------------------------------------------------!

    logical :: urban_flag
    integer :: iswater
    integer :: isbarren
    integer :: isice
    integer :: eblforest

    real :: ch2op              !maximum intercepted h2o per unit lai+sai (mm)
    real :: dleaf              !characteristic leaf dimension (m)
    real :: z0mvt              !momentum roughness length (m)
    real :: hvt                !top of canopy (m)
    real :: hvb                !bottom of canopy (m)
    real :: den                !tree density (no. of trunks per m2)
    real :: rc                 !tree crown radius (m)
    real :: mfsno              !snowmelt m parameter ()
    real :: saim(12)           !monthly stem area index, one-sided
    real :: laim(12)           !monthly leaf area index, one-sided
    real :: sla                !single-side leaf area per kg [m2/kg]
    real :: dilefc             !coeficient for leaf stress death [1/s]
    real :: dilefw             !coeficient for leaf stress death [1/s]
    real :: fragr              !fraction of growth respiration  !original was 0.3 
    real :: ltovrc             !leaf turnover [1/s]

    real :: c3psn              !photosynthetic pathway: 0. = c4, 1. = c3
    real :: kc25               !co2 michaelis-menten constant at 25c (pa)
    real :: akc                !q10 for kc25
    real :: ko25               !o2 michaelis-menten constant at 25c (pa)
    real :: ako                !q10 for ko25
    real :: vcmx25             !maximum rate of carboxylation at 25c (umol co2/m**2/s)
    real :: avcmx              !q10 for vcmx25
    real :: bp                 !minimum leaf conductance (umol/m**2/s)
    real :: mp                 !slope of conductance-to-photosynthesis relationship
    real :: qe25               !quantum efficiency at 25c (umol co2 / umol photon)
    real :: aqe                !q10 for qe25
    real :: rmf25              !leaf maintenance respiration at 25c (umol co2/m**2/s)
    real :: rms25              !stem maintenance respiration at 25c (umol co2/kg bio/s)
    real :: rmr25              !root maintenance respiration at 25c (umol co2/kg bio/s)
    real :: arm                !q10 for maintenance respiration
    real :: folnmx             !foliage nitrogen concentration when f(n)=1 (%)
    real :: tmin               !minimum temperature for photosynthesis (k)
       
    real :: xl                 !leaf/stem orientation index
    real :: rhol(mband)        !leaf reflectance: 1=vis, 2=nir
    real :: rhos(mband)        !stem reflectance: 1=vis, 2=nir
    real :: taul(mband)        !leaf transmittance: 1=vis, 2=nir
    real :: taus(mband)        !stem transmittance: 1=vis, 2=nir

    real :: mrp                !microbial respiration parameter (umol co2 /kg c/ s)
    real :: cwpvt              !empirical canopy wind parameter

    real :: wrrat              !wood to non-wood ratio
    real :: wdpool             !wood pool (switch 1 or 0) depending on woody or not [-]
    real :: tdlef              !characteristic t for leaf freezing [k]

  integer :: nroot              !number of soil layers with root present
     real :: rgl                !parameter used in radiation stress function
     real :: rsmin              !minimum stomatal resistance [s m-1]
     real :: hs                 !parameter used in vapor pressure deficit function
     real :: topt               !optimum transpiration air temperature [k]
     real :: rsmax              !maximal stomatal resistance [s m-1]

     real :: slarea
     real :: eps(5)

!------------------------------------------------------------------------------------------!
! from the rad section of mptable.tbl
!------------------------------------------------------------------------------------------!

     real :: albsat(mband)       !saturated soil albedos: 1=vis, 2=nir
     real :: albdry(mband)       !dry soil albedos: 1=vis, 2=nir
     real :: albice(mband)       !albedo land ice: 1=vis, 2=nir
     real :: alblak(mband)       !albedo frozen lakes: 1=vis, 2=nir
     real :: omegas(mband)       !two-stream parameter omega for snow
     real :: betads              !two-stream parameter betad for snow
     real :: betais              !two-stream parameter betad for snow
     real :: eg(2)               !emissivity

!------------------------------------------------------------------------------------------!
! from the globals section of mptable.tbl
!------------------------------------------------------------------------------------------!
 
     real :: co2          !co2 partial pressure
     real :: o2           !o2 partial pressure
     real :: timean       !gridcell mean topgraphic index (global mean)
     real :: fsatmx       !maximum surface saturated fraction (global mean)
     real :: z0sno        !snow surface roughness length (m) (0.002)
     real :: ssi          !liquid water holding capacity for snowpack (m3/m3)
     real :: swemx        !new snow mass to fully cover old snow (mm)

!------------------------------------------------------------------------------------------!
! from the soilparm.tbl tables, as functions of soil category.
!------------------------------------------------------------------------------------------!
     real :: bexp         !b parameter
     real :: smcdry       !dry soil moisture threshold where direct evap from top
                          !layer ends (volumetric) (not used mb: 20140718)
     real :: smcwlt       !wilting point soil moisture (volumetric)
     real :: smcref       !reference soil moisture (field capacity) (volumetric)
     real :: smcmax       !porosity, saturated value of soil moisture (volumetric)
     real :: f1           !soil thermal diffusivity/conductivity coef (not used mb: 20140718)
     real :: psisat       !saturated soil matric potential
     real :: dksat        !saturated soil hydraulic conductivity
     real :: dwsat        !saturated soil hydraulic diffusivity
     real :: quartz       !soil quartz content
!------------------------------------------------------------------------------------------!
! from the genparm.tbl file
!------------------------------------------------------------------------------------------!
     real :: slope       !slope index (0 - 1)
     real :: csoil       !vol. soil heat capacity [j/m3/k]
     real :: zbot        !depth (m) of lower boundary soil temperature
     real :: czil        !calculate roughness length of heat

     real :: kdt         !used in compute maximum infiltration rate (in infil)
     real :: frzx        !used in compute maximum infiltration rate (in infil)

  end type noahmp_parameters

contains
!
!== begin noahmp_sflx ==============================================================================

!>\ingroup NoahMP_LSM
  subroutine noahmp_sflx (parameters, &
                   iloc    , jloc    , lat     , yearlen , julian  , cosz    , & ! in : time/space-related
                   dt      , dx      , dz8w    , nsoil   , zsoil   , nsnow   , & ! in : model configuration 
                   shdfac  , shdmax  , vegtyp  , ice     , ist     ,           & ! in : vegetation/soil characteristics
                   smceq   ,                                                   & ! in : vegetation/soil characteristics
                   sfctmp  , sfcprs  , psfc    , uu      , vv      , q2      , & ! in : forcing
                   qc      , soldn   , lwdn    ,                               & ! in : forcing
	           prcpconv, prcpnonc, prcpshcv, prcpsnow, prcpgrpl, prcphail, & ! in : forcing
                   tbot    , co2air  , o2air   , foln    , ficeold , zlvl    , & ! in : forcing
                   lheatstrg                                                 , & ! in : canopy heat storage
                   albold  , sneqvo  ,                                         & ! in/out : 
                   stc     , sh2o    , smc     , tah     , eah     , fwet    , & ! in/out : 
                   canliq  , canice  , tv      , tg      , qsfc    , qsnow   , & ! in/out : 
                   isnow   , zsnso   , snowh   , sneqv   , snice   , snliq   , & ! in/out : 
                   zwt     , wa      , wt      , wslake  , lfmass  , rtmass  , & ! in/out : 
                   stmass  , wood    , stblcp  , fastcp  , lai     , sai     , & ! in/out : 
                   cm      , ch      , tauss   ,                               & ! in/out : 
                   smcwtd  ,deeprech , rech    , cpfac                       , & ! in/out :
		   z0wrf   , &
                   fsa     , fsr     , fira    , fshx    , ssoil   , fcev    , & ! out : 
                   fgev    , fctr    , ecan    , etran   , edir    , trad    , & ! out :
                   tgb     , tgv     , t2mv    , t2mb    , q2v     , q2b     , & ! out :
                   runsrf  , runsub  , apar    , psn     , sav     , sag     , & ! out :
                   fsno    , nee     , gpp     , npp     , fveg    , albedo  , & ! out :
                   qsnbot  , ponding , ponding1, ponding2, rssun   , rssha   , & ! out :
                   bgap    , wgap    , chv     , chb     , emissi  ,           & ! out :
		   shg     , shc     , shb     , evg     , evb     , ghv     , & ! out :
		   ghb     , irg     , irc     , irb     , tr      , evc     , & ! out :
		   chleaf  , chuc    , chv2    , chb2    , fpice   , pahv    , &
#ifdef CCPP
                   pahg    , pahb    , pah     , esnow, errmsg, errflg)
#else
                   pahg    , pahb    , pah     , esnow)
#endif

! --------------------------------------------------------------------------------------------------
! initial code: guo-yue niu, oct. 2007
! --------------------------------------------------------------------------------------------------
  implicit none
! --------------------------------------------------------------------------------------------------
! input
  type (noahmp_parameters), intent(in) :: parameters

  integer                        , intent(in)    :: ice    !ice (ice = 1)
  integer                        , intent(in)    :: ist    !surface type 1->soil; 2->lake
  integer                        , intent(in)    :: vegtyp !vegetation type 
  integer                        , intent(in)    :: nsnow  !maximum no. of snow layers        
  integer                        , intent(in)    :: nsoil  !no. of soil layers        
  integer                        , intent(in)    :: iloc   !grid index
  integer                        , intent(in)    :: jloc   !grid index
  real                           , intent(in)    :: dt     !time step [sec]
  real, dimension(       1:nsoil), intent(in)    :: zsoil  !layer-bottom depth from soil surf (m)
  real                           , intent(in)    :: q2     !mixing ratio (kg/kg) lowest model layer
  real                           , intent(in)    :: sfctmp !surface air temperature [k]
  real                           , intent(in)    :: uu     !wind speed in eastward dir (m/s)
  real                           , intent(in)    :: vv     !wind speed in northward dir (m/s)
  real                           , intent(in)    :: soldn  !downward shortwave radiation (w/m2)
  real                           , intent(in)    :: lwdn   !downward longwave radiation (w/m2)
  real                           , intent(in)    :: sfcprs !pressure (pa)
  real                           , intent(inout) :: zlvl   !reference height (m)
  logical                        , intent(in)    :: lheatstrg ! flag for canopy heat storage parameterization       
  real                           , intent(in)    :: cosz   !cosine solar zenith angle [0-1]
  real                           , intent(in)    :: tbot   !bottom condition for soil temp. [k]
  real                           , intent(in)    :: foln   !foliage nitrogen (%) [1-saturated]
  real                           , intent(in)    :: shdfac !green vegetation fraction [0.0-1.0]
  integer                        , intent(in)    :: yearlen!number of days in the particular year.
  real                           , intent(in)    :: julian !julian day of year (floating point)
  real                           , intent(in)    :: lat    !latitude (radians)
  real, dimension(-nsnow+1:    0), intent(in)    :: ficeold!ice fraction at last timestep
  real, dimension(       1:nsoil), intent(in)    :: smceq  !equilibrium soil water  content [m3/m3]
  real                           , intent(in)    :: prcpconv ! convective precipitation entering  [mm/s]    ! mb/an : v3.7
  real                           , intent(in)    :: prcpnonc ! non-convective precipitation entering [mm/s] ! mb/an : v3.7
  real                           , intent(in)    :: prcpshcv ! shallow convective precip entering  [mm/s]   ! mb/an : v3.7
  real                           , intent(in)    :: prcpsnow ! snow entering land model [mm/s]              ! mb/an : v3.7
  real                           , intent(in)    :: prcpgrpl ! graupel entering land model [mm/s]           ! mb/an : v3.7
  real                           , intent(in)    :: prcphail ! hail entering land model [mm/s]              ! mb/an : v3.7

!jref:start; in 
  real                           , intent(in)    :: qc     !cloud water mixing ratio
  real                           , intent(inout)    :: qsfc   !mixing ratio at lowest model layer
  real                           , intent(in)    :: psfc   !pressure at lowest model layer
  real                           , intent(in)    :: dz8w   !thickness of lowest layer
  real                           , intent(in)    :: dx
  real                           , intent(in)    :: shdmax  !yearly max vegetation fraction
!jref:end


! input/output : need arbitary intial values
  real                           , intent(inout) :: qsnow  !snowfall [mm/s]
  real                           , intent(inout) :: fwet   !wetted or snowed fraction of canopy (-)
  real                           , intent(inout) :: sneqvo !snow mass at last time step (mm)
  real                           , intent(inout) :: eah    !canopy air vapor pressure (pa)
  real                           , intent(inout) :: tah    !canopy air tmeperature (k)
  real                           , intent(inout) :: albold !snow albedo at last time step (class type)
  real                           , intent(inout) :: cm     !momentum drag coefficient
  real                           , intent(inout) :: ch     !sensible heat exchange coefficient
  real                           , intent(inout) :: tauss  !non-dimensional snow age

! prognostic variables
  integer                        , intent(inout) :: isnow  !actual no. of snow layers [-]
  real                           , intent(inout) :: canliq !intercepted liquid water (mm)
  real                           , intent(inout) :: canice !intercepted ice mass (mm)
  real                           , intent(inout) :: sneqv  !snow water eqv. [mm]
  real, dimension(       1:nsoil), intent(inout) :: smc    !soil moisture (ice + liq.) [m3/m3]
  real, dimension(-nsnow+1:nsoil), intent(inout) :: zsnso  !layer-bottom depth from snow surf [m]
  real                           , intent(inout) :: snowh  !snow height [m]
  real, dimension(-nsnow+1:    0), intent(inout) :: snice  !snow layer ice [mm]
  real, dimension(-nsnow+1:    0), intent(inout) :: snliq  !snow layer liquid water [mm]
  real                           , intent(inout) :: tv     !vegetation temperature (k)
  real                           , intent(inout) :: tg     !ground temperature (k)
  real, dimension(-nsnow+1:nsoil), intent(inout) :: stc    !snow/soil temperature [k]
  real, dimension(       1:nsoil), intent(inout) :: sh2o   !liquid soil moisture [m3/m3]
  real                           , intent(inout) :: zwt    !depth to water table [m]
  real                           , intent(inout) :: wa     !water storage in aquifer [mm]
  real                           , intent(inout) :: wt     !water in aquifer&saturated soil [mm]
  real                           , intent(inout) :: wslake !lake water storage (can be neg.) (mm)
  real,                            intent(inout) :: smcwtd !soil water content between bottom of the soil and water table [m3/m3]
  real,                            intent(inout) :: deeprech !recharge to or from the water table when deep [m]
  real,                            intent(inout) :: rech !recharge to or from the water table when shallow [m] (diagnostic)
  real,                            intent(inout) :: cpfac  ! heat capacity enhancement factor due to heat storage

! output
  real                           , intent(out)   :: z0wrf  !combined z0 sent to coupled model
  real                           , intent(out)   :: fsa    !total absorbed solar radiation (w/m2)
  real                           , intent(out)   :: fsr    !total reflected solar radiation (w/m2)
  real                           , intent(out)   :: fira   !total net lw rad (w/m2)  [+ to atm]
  real                           , intent(out)   :: fshx   !total sensible heat (w/m2) [+ to atm]
  real                           , intent(out)   :: fcev   !canopy evap heat (w/m2) [+ to atm]
  real                           , intent(out)   :: fgev   !ground evap heat (w/m2) [+ to atm]
  real                           , intent(out)   :: fctr   !transpiration heat (w/m2) [+ to atm]
  real                           , intent(out)   :: ssoil  !ground heat flux (w/m2)   [+ to soil]
  real                           , intent(out)   :: trad   !surface radiative temperature (k)
  real                                           :: ts     !surface temperature (k)
  real                           , intent(out)   :: ecan   !evaporation of intercepted water (mm/s)
  real                           , intent(out)   :: etran  !transpiration rate (mm/s)
  real                           , intent(out)   :: edir   !soil surface evaporation rate (mm/s]
  real                           , intent(out)   :: runsrf !surface runoff [mm/s] 
  real                           , intent(out)   :: runsub !baseflow (saturation excess) [mm/s]
  real                           , intent(out)   :: psn    !total photosynthesis (umol co2/m2/s) [+]
  real                           , intent(out)   :: apar   !photosyn active energy by canopy (w/m2)
  real                           , intent(out)   :: sav    !solar rad absorbed by veg. (w/m2)
  real                           , intent(out)   :: sag    !solar rad absorbed by ground (w/m2)
  real                           , intent(out)   :: fsno   !snow cover fraction on the ground (-)
  real                           , intent(out)   :: fveg   !green vegetation fraction [0.0-1.0]
  real                           , intent(out)   :: albedo !surface albedo [-]
  real                                           :: errwat !water error [kg m{-2}]
  real                           , intent(out)   :: qsnbot !snowmelt out bottom of pack [mm/s]
  real                           , intent(out)   :: ponding!surface ponding [mm]
  real                           , intent(out)   :: ponding1!surface ponding [mm]
  real                           , intent(out)   :: ponding2!surface ponding [mm]
  real                           , intent(out)   :: esnow

!jref:start; output
  real                           , intent(out)     :: t2mv   !2-m air temperature over vegetated part [k]
  real                           , intent(out)     :: t2mb   !2-m air temperature over bare ground part [k]
  real, intent(out) :: rssun        !sunlit leaf stomatal resistance (s/m)
  real, intent(out) :: rssha        !shaded leaf stomatal resistance (s/m)
  real, intent(out) :: bgap
  real, intent(out) :: wgap
  real, intent(out) :: tgv
  real, intent(out) :: tgb
  real              :: q1
  real, intent(out) :: emissi
!jref:end
#ifdef CCPP
  character(len=*), intent(inout)    :: errmsg
  integer,          intent(inout)    :: errflg
#endif

! local
  integer                                        :: iz     !do-loop index
  integer, dimension(-nsnow+1:nsoil)             :: imelt  !phase change index [1-melt; 2-freeze]
  real                                           :: cmc    !intercepted water (canice+canliq) (mm)
  real                                           :: taux   !wind stress: e-w (n/m2)
  real                                           :: tauy   !wind stress: n-s (n/m2)
  real                                           :: rhoair !density air (kg/m3)
  real                                           :: fsh    !total sensible heat (w/m2) [+ to atm]
!  real, dimension(       1:    5)                :: vocflx !voc fluxes [ug c m-2 h-1]
  real, dimension(-nsnow+1:nsoil)                :: dzsnso !snow/soil layer thickness [m]
  real                                           :: thair  !potential temperature (k)
  real                                           :: qair   !specific humidity (kg/kg) (q2/(1+q2))
  real                                           :: eair   !vapor pressure air (pa)
  real, dimension(       1:    2)                :: solad  !incoming direct solar rad (w/m2)
  real, dimension(       1:    2)                :: solai  !incoming diffuse solar rad (w/m2)
  real                                           :: qprecc !convective precipitation (mm/s)
  real                                           :: qprecl !large-scale precipitation (mm/s)
  real                                           :: igs    !growing season index (0=off, 1=on)
  real                                           :: elai   !leaf area index, after burying by snow
  real                                           :: esai   !stem area index, after burying by snow
  real                                           :: bevap  !soil water evaporation factor (0 - 1)
  real, dimension(       1:nsoil)                :: btrani !soil water transpiration factor (0 - 1)
  real                                           :: btran  !soil water transpiration factor (0 - 1)
  real                                           :: qin    !groundwater recharge [mm/s]
  real                                           :: qdis   !groundwater discharge [mm/s]
  real, dimension(       1:nsoil)                :: sice   !soil ice content (m3/m3)
  real, dimension(-nsnow+1:    0)                :: snicev !partial volume ice of snow [m3/m3]
  real, dimension(-nsnow+1:    0)                :: snliqv !partial volume liq of snow [m3/m3]
  real, dimension(-nsnow+1:    0)                :: epore  !effective porosity [m3/m3]
  real                                           :: totsc  !total soil carbon (g/m2)
  real                                           :: totlb  !total living carbon (g/m2)
  real                                           :: t2m    !2-meter air temperature (k)
  real                                           :: qdew   !ground surface dew rate [mm/s]
  real                                           :: qvap   !ground surface evap. rate [mm/s]
  real                                           :: lathea !latent heat [j/kg]
  real                                           :: swdown !downward solar [w/m2]
  real                                           :: qmelt  !snowmelt [mm/s]
  real                                           :: beg_wb !water storage at begin of a step [mm]
  real,intent(out)                                              :: irc    !canopy net lw rad. [w/m2] [+ to atm]
  real,intent(out)                                              :: irg    !ground net lw rad. [w/m2] [+ to atm]
  real,intent(out)                                              :: shc    !canopy sen. heat [w/m2]   [+ to atm]
  real,intent(out)                                              :: shg    !ground sen. heat [w/m2]   [+ to atm]
  real,intent(out)                                              :: evg    !ground evap. heat [w/m2]  [+ to atm]
  real,intent(out)                                              :: ghv    !ground heat flux [w/m2]  [+ to soil]
  real,intent(out)                                              :: irb    !net longwave rad. [w/m2] [+ to atm]
  real,intent(out)                                              :: shb    !sensible heat [w/m2]     [+ to atm]
  real,intent(out)                                              :: evb    !evaporation heat [w/m2]  [+ to atm]
  real,intent(out)                                              :: ghb    !ground heat flux [w/m2] [+ to soil]
  real,intent(out)                                              :: evc    !canopy evap. heat [w/m2]  [+ to atm]
  real,intent(out)                                              :: tr     !transpiration heat [w/m2] [+ to atm]
  real, intent(out)   :: fpice   !snow fraction in precipitation
  real, intent(out)   :: pahv    !precipitation advected heat - vegetation net (w/m2)
  real, intent(out)   :: pahg    !precipitation advected heat - under canopy net (w/m2)
  real, intent(out)   :: pahb    !precipitation advected heat - bare ground net (w/m2)
  real, intent(out)                                           :: pah     !precipitation advected heat - total (w/m2)

!jref:start 
  real                                           :: fsrv
  real                                           :: fsrg
  real,intent(out)                               :: q2v
  real,intent(out)                               :: q2b
  real :: q2e
  real :: qfx
  real,intent(out)                               :: chv    !sensible heat exchange coefficient over vegetated fraction
  real,intent(out)                               :: chb    !sensible heat exchange coefficient over bare-ground
  real,intent(out)                               :: chleaf !leaf exchange coefficient
  real,intent(out)                               :: chuc   !under canopy exchange coefficient
  real,intent(out)                               :: chv2    !sensible heat exchange coefficient over vegetated fraction
  real,intent(out)                               :: chb2    !sensible heat exchange coefficient over bare-ground
!jref:end  

! carbon
! inputs
  real                           , intent(in)    :: co2air !atmospheric co2 concentration (pa)
  real                           , intent(in)    :: o2air  !atmospheric o2 concentration (pa)

! inputs and outputs : prognostic variables
  real                        , intent(inout)    :: lfmass !leaf mass [g/m2]
  real                        , intent(inout)    :: rtmass !mass of fine roots [g/m2]
  real                        , intent(inout)    :: stmass !stem mass [g/m2]
  real                        , intent(inout)    :: wood   !mass of wood (incl. woody roots) [g/m2]
  real                        , intent(inout)    :: stblcp !stable carbon in deep soil [g/m2]
  real                        , intent(inout)    :: fastcp !short-lived carbon, shallow soil [g/m2]
  real                        , intent(inout)    :: lai    !leaf area index [-]
  real                        , intent(inout)    :: sai    !stem area index [-]

! outputs
  real                          , intent(out)    :: nee    !net ecosys exchange (g/m2/s co2)
  real                          , intent(out)    :: gpp    !net instantaneous assimilation [g/m2/s c]
  real                          , intent(out)    :: npp    !net primary productivity [g/m2/s c]
  real                                           :: autors !net ecosystem respiration (g/m2/s c)
  real                                           :: heters !organic respiration (g/m2/s c)
  real                                           :: troot  !root-zone averaged temperature (k)
  real                                           :: bdfall   !bulk density of new snow (kg/m3)    ! mb/an: v3.7
  real                                           :: rain     !rain rate                   (mm/s)  ! mb/an: v3.7
  real                                           :: snow     !liquid equivalent snow rate (mm/s)  ! mb/an: v3.7
  real                                           :: fp                                            ! mb/an: v3.7
  real                                           :: prcp                                          ! mb/an: v3.7
!more local variables for precip heat mb
  real                                           :: qintr   !interception rate for rain (mm/s)
  real                                           :: qdripr  !drip rate for rain (mm/s)
  real                                           :: qthror  !throughfall for rain (mm/s)
  real                                           :: qints   !interception (loading) rate for snowfall (mm/s)
  real                                           :: qdrips  !drip (unloading) rate for intercepted snow (mm/s)
  real                                           :: qthros  !throughfall of snowfall (mm/s)
  real                                           :: qrain   !rain at ground srf (mm/s) [+]
  real                                           :: snowhin !snow depth increasing rate (m/s)
  real                                 :: latheav !latent heat vap./sublimation (j/kg)
  real                                 :: latheag !latent heat vap./sublimation (j/kg)
  logical                             :: frozen_ground ! used to define latent heat pathway
  logical                             :: frozen_canopy ! used to define latent heat pathway
  
  ! intent (out) variables need to be assigned a value.  these normally get assigned values
  ! only if dveg == 2.
  nee = 0.0
  npp = 0.0
  gpp = 0.0
      pahv  = 0.
      pahg  = 0.
      pahb  = 0.
      pah  = 0.

! --------------------------------------------------------------------------------------------------
! re-process atmospheric forcing

   call atm (parameters,sfcprs  ,sfctmp   ,q2      ,                            &
             prcpconv, prcpnonc,prcpshcv,prcpsnow,prcpgrpl,prcphail, &
             soldn   ,cosz     ,thair   ,qair    ,                   & 
             eair    ,rhoair   ,qprecc  ,qprecl  ,solad   ,solai   , &
             swdown  ,bdfall   ,rain    ,snow    ,fp      ,fpice   , prcp )     

! snow/soil layer thickness (m)

     do iz = isnow+1, nsoil
         if(iz == isnow+1) then
           dzsnso(iz) = - zsnso(iz)
         else
           dzsnso(iz) = zsnso(iz-1) - zsnso(iz)
         end if
     end do

! root-zone temperature

     troot  = 0.
     do iz=1,parameters%nroot
        troot = troot + stc(iz)*dzsnso(iz)/(-zsoil(parameters%nroot))
     enddo

! total water storage for water balance check
    
     if(ist == 1) then
     beg_wb = canliq + canice + sneqv + wa
     do iz = 1,nsoil
        beg_wb = beg_wb + smc(iz) * dzsnso(iz) * 1000.
     end do
     end if

! vegetation phenology

     call phenology (parameters,vegtyp , snowh  , tv     , lat   , yearlen , julian , & !in
                     lai    , sai    , troot  , elai    , esai   ,igs)

!input gvf should be consistent with lai
     if(dveg == 1) then
        fveg = shdfac
        if(fveg <= 0.05) fveg = 0.05
     else if (dveg == 2 .or. dveg == 3) then
        fveg = 1.-exp(-0.52*(lai+sai))
        if(fveg <= 0.05) fveg = 0.05
     else if (dveg == 4 .or. dveg == 5) then
        fveg = shdmax
        if(fveg <= 0.05) fveg = 0.05
     else
        write(*,*) "-------- fatal called in sflx -----------"
#ifdef CCPP
        errflg = 1
        errmsg = "namelist parameter dveg unknown"
        return 
#else
        call wrf_error_fatal("namelist parameter dveg unknown") 
#endif
     endif
     if(parameters%urban_flag .or. vegtyp == parameters%isbarren) fveg = 0.0
     if(elai+esai == 0.0) fveg = 0.0

    call precip_heat(parameters,iloc   ,jloc   ,vegtyp ,dt     ,uu     ,vv     , & !in
                     elai   ,esai   ,fveg   ,ist    ,                 & !in
                     bdfall ,rain   ,snow   ,fp     ,                 & !in
                     canliq ,canice ,tv     ,sfctmp ,tg     ,         & !in
                     qintr  ,qdripr ,qthror ,qints  ,qdrips ,qthros , & !out
                     pahv   ,pahg   ,pahb   ,qrain  ,qsnow  ,snowhin, & !out
	             fwet   ,cmc                                    )   !out

! compute energy budget (momentum & energy fluxes and phase changes) 

    call energy (parameters,ice    ,vegtyp ,ist    ,nsnow  ,nsoil  , & !in
                 isnow  ,dt     ,rhoair ,sfcprs ,qair   , & !in
                 sfctmp ,thair  ,lwdn   ,uu     ,vv     ,zlvl   , & !in
                 lheatstrg                                      , & !in
                 co2air ,o2air  ,solad  ,solai  ,cosz   ,igs    , & !in
                 eair   ,tbot   ,zsnso  ,zsoil  , & !in
                 elai   ,esai   ,fwet   ,foln   ,         & !in
                 fveg   ,pahv   ,pahg   ,pahb   ,                 & !in
                 qsnow  ,dzsnso ,lat    ,canliq ,canice ,iloc, jloc , & !in
                 z0wrf  ,                                         &
                 imelt  ,snicev ,snliqv ,epore  ,t2m    ,fsno   , & !out
                 sav    ,sag    ,qmelt  ,fsa    ,fsr    ,taux   , & !out
                 tauy   ,fira   ,fsh    ,fshx   ,fcev   ,fgev   ,fctr   , & !out
                 trad   ,psn    ,apar   ,ssoil  ,btrani ,btran  , & !out
                 ponding,ts     ,latheav , latheag , frozen_canopy,frozen_ground,                         & !out
                 tv     ,tg     ,stc    ,snowh  ,eah    ,tah    , & !inout
                 sneqvo ,sneqv  ,sh2o   ,smc    ,snice  ,snliq  , & !inout
                 albold ,cm     ,ch     ,dx     ,dz8w   ,q2     , & !inout
#ifdef CCPP
                 tauss  ,cpfac  ,errmsg ,errflg ,                 & !inout
#else
                 tauss  ,cpfac  ,                                 & !inout
#endif
!jref:start
                 qc     ,qsfc   ,psfc   , & !in 
                 t2mv   ,t2mb  ,fsrv   , &
                 fsrg   ,rssun   ,rssha ,bgap   ,wgap, tgv,tgb,&
                 q1     ,q2v    ,q2b    ,q2e    ,chv   ,chb     , & !out
                 emissi ,pah    ,                                 &
                 shg,shc,shb,evg,evb,ghv,ghb,irg,irc,irb,tr,evc,chleaf,chuc,chv2,chb2 )                                            !out
!jref:end
#ifdef CCPP
    if (errflg /= 0) return
#endif
    sice(:) = max(0.0, smc(:) - sh2o(:))   
    sneqvo  = sneqv

    qvap = max( fgev/latheag, 0.)       ! positive part of fgev; barlage change to ground v3.6
    qdew = abs( min(fgev/latheag, 0.))  ! negative part of fgev
    edir = qvap - qdew

! compute water budgets (water storages, et components, and runoff)

     call water (parameters,vegtyp ,nsnow  ,nsoil  ,imelt  ,dt     ,uu     , & !in
                 vv     ,fcev   ,fctr   ,qprecc ,qprecl ,elai   , & !in
                 esai   ,sfctmp ,qvap   ,qdew   ,zsoil  ,btrani , & !in
                 ficeold,ponding,tg     ,ist    ,fveg   ,iloc,jloc , smceq , & !in
                 bdfall ,fp     ,rain   ,snow   ,                 & !in  mb/an: v3.7
		 qsnow  ,qrain  ,snowhin,latheav,latheag,frozen_canopy,frozen_ground,  & !in  mb
                 isnow  ,canliq ,canice ,tv     ,snowh  ,sneqv  , & !inout
                 snice  ,snliq  ,stc    ,zsnso  ,sh2o   ,smc    , & !inout
                 sice   ,zwt    ,wa     ,wt     ,dzsnso ,wslake , & !inout
                 smcwtd ,deeprech,rech                          , & !inout
                 cmc    ,ecan   ,etran  ,fwet   ,runsrf ,runsub , & !out
                 qin    ,qdis   ,ponding1       ,ponding2,&
                 qsnbot ,esnow   )  !out

!     write(*,'(a20,10f15.5)') 'sflx:runoff=',runsrf*dt,runsub*dt,edir*dt

! compute carbon budgets (carbon storages and co2 & bvoc fluxes)

   if (dveg == 2 .or. dveg == 5) then
    call carbon (parameters,nsnow  ,nsoil  ,vegtyp ,dt     ,zsoil  , & !in
                 dzsnso ,stc    ,smc    ,tv     ,tg     ,psn    , & !in
                 foln   ,btran  ,apar   ,fveg   ,igs    , & !in
                 troot  ,ist    ,lat    ,iloc   ,jloc   , & !in
                 lfmass ,rtmass ,stmass ,wood   ,stblcp ,fastcp , & !inout
                 gpp    ,npp    ,nee    ,autors ,heters ,totsc  , & !out
                 totlb  ,lai    ,sai    )                   !out
   end if

! water and energy balance check

     call error (parameters,swdown ,fsa    ,fsr    ,fira   ,fsh   ,fcev   , & !in
                 fgev   ,fctr   ,ssoil  ,beg_wb ,canliq ,canice , & !in
                 sneqv  ,wa     ,smc    ,dzsnso ,prcp   ,ecan   , & !in
                 etran  ,edir   ,runsrf ,runsub ,dt     ,nsoil  , & !in
                 nsnow  ,ist    ,errwat ,iloc   , jloc  ,fveg   , &
                 sav    ,sag    ,fsrv   ,fsrg   ,zwt    ,pah    , &
#ifdef CCPP
                 pahv   ,pahg   ,pahb   ,errmsg, errflg)   !in ( except errwat [out] and errmsg, errflg [inout] )
#else
                 pahv   ,pahg   ,pahb   )   !in ( except errwat, which is out )
#endif

#ifdef CCPP
     if (errflg /= 0) return
#endif

! urban - jref
    qfx = etran + ecan + edir
    if ( parameters%urban_flag ) then
       qsfc = (qfx/rhoair*ch) + qair
       q2b = qsfc
    end if

    if(snowh <= 1.e-6 .or. sneqv <= 1.e-3) then
     snowh = 0.0
     sneqv = 0.0
    end if

    if(swdown.ne.0.) then
      albedo = fsr / swdown
    else
      albedo = -999.9
    end if
    

  end subroutine noahmp_sflx

!== begin atm ======================================================================================

!>\ingroup NoahMP_LSM
  subroutine atm (parameters,sfcprs  ,sfctmp   ,q2      ,                             &
                  prcpconv,prcpnonc ,prcpshcv,prcpsnow,prcpgrpl,prcphail , &
                  soldn   ,cosz     ,thair   ,qair    ,                    & 
                  eair    ,rhoair   ,qprecc  ,qprecl  ,solad   , solai   , &
		  swdown  ,bdfall   ,rain    ,snow    ,fp      , fpice   ,prcp )     
! --------------------------------------------------------------------------------------------------
! re-process atmospheric forcing
! ----------------------------------------------------------------------
  implicit none
! --------------------------------------------------------------------------------------------------
! inputs

  type (noahmp_parameters), intent(in) :: parameters
  real                          , intent(in)  :: sfcprs !pressure (pa)
  real                          , intent(in)  :: sfctmp !surface air temperature [k]
  real                          , intent(in)  :: q2     !mixing ratio (kg/kg)
  real                          , intent(in)  :: prcpconv ! convective precipitation entering  [mm/s]    ! mb/an : v3.7
  real                          , intent(in)  :: prcpnonc ! non-convective precipitation entering [mm/s] ! mb/an : v3.7
  real                          , intent(in)  :: prcpshcv ! shallow convective precip entering  [mm/s]   ! mb/an : v3.7
  real                          , intent(in)  :: prcpsnow ! snow entering land model [mm/s]              ! mb/an : v3.7
  real                          , intent(in)  :: prcpgrpl ! graupel entering land model [mm/s]           ! mb/an : v3.7
  real                          , intent(in)  :: prcphail ! hail entering land model [mm/s]              ! mb/an : v3.7
  real                          , intent(in)  :: soldn  !downward shortwave radiation (w/m2)
  real                          , intent(in)  :: cosz   !cosine solar zenith angle [0-1]

! outputs

  real                          , intent(out) :: thair  !potential temperature (k)
  real                          , intent(out) :: qair   !specific humidity (kg/kg) (q2/(1+q2))
  real                          , intent(out) :: eair   !vapor pressure air (pa)
  real                          , intent(out) :: rhoair !density air (kg/m3)
  real                          , intent(out) :: qprecc !convective precipitation (mm/s)
  real                          , intent(out) :: qprecl !large-scale precipitation (mm/s)
  real, dimension(       1:   2), intent(out) :: solad  !incoming direct solar radiation (w/m2)
  real, dimension(       1:   2), intent(out) :: solai  !incoming diffuse solar radiation (w/m2)
  real                          , intent(out) :: swdown !downward solar filtered by sun angle [w/m2]
  real                          , intent(out) :: bdfall  !!bulk density of snowfall (kg/m3) ajn
  real                          , intent(out) :: rain    !rainfall (mm/s) ajn
  real                          , intent(out) :: snow    !liquid equivalent snowfall (mm/s) ajn
  real                          , intent(out) :: fp      !fraction of area receiving precipitation  ajn
  real                          , intent(out) :: fpice   !fraction of ice                ajn
  real                          , intent(out) :: prcp    !total precipitation [mm/s]     ! mb/an : v3.7

!locals

  real                                        :: pair   !atm bottom level pressure (pa)
  real                                        :: prcp_frozen   !total frozen precipitation [mm/s] ! mb/an : v3.7
  real, parameter                             :: rho_grpl = 500.0  ! graupel bulk density [kg/m3] ! mb/an : v3.7
  real, parameter                             :: rho_hail = 917.0  ! hail bulk density [kg/m3]    ! mb/an : v3.7
! --------------------------------------------------------------------------------------------------

!jref: seems like pair should be p1000mb??
       pair   = sfcprs                   ! atm bottom level pressure (pa)
       thair  = sfctmp * (sfcprs/pair)**(rair/cpair) 

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

       prcp = prcpconv + prcpnonc + prcpshcv

!      if(opt_snf == 4) then
         qprecc = prcpconv + prcpshcv
	 qprecl = prcpnonc
!      else
!        qprecc = 0.10 * prcp          ! should be from the atmospheric model
!        qprecl = 0.90 * prcp          ! should be from the atmospheric model
!      end if

! fractional area that receives precipitation (see, niu et al. 2005)
   
    fp = 0.0
    if(qprecc + qprecl > 0.) & 
       fp = (qprecc + qprecl) / (10.*qprecc + qprecl)

! partition precipitation into rain and snow. moved from canwat mb/an: v3.7

! jordan (1991)

     if(opt_snf == 1) then
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

! hedstrom nr and jw pomeroy (1998), hydrol. processes, 12, 1611-1625
! fresh snow density

     bdfall = min(120.,67.92+51.25*exp((sfctmp-tfrz)/2.59))       !mb/an: change to min  
     if(opt_snf == 4) then
        prcp_frozen = prcpsnow + prcpgrpl + prcphail
        if(prcpnonc > 0. .and. prcp_frozen > 0.) then
	  fpice = min(1.0,prcp_frozen/prcp)
	  fpice = max(0.0,fpice)
	  bdfall = bdfall*(prcpsnow/prcp_frozen) + rho_grpl*(prcpgrpl/prcp_frozen) + &
	             rho_hail*(prcphail/prcp_frozen)
	else
	  fpice = 0.0
        endif
	
     endif

     rain   = prcp * (1.-fpice)
     snow   = prcp * fpice


  end subroutine atm

!== begin phenology ================================================================================

!>\ingroup NoahMP_LSM
  subroutine phenology (parameters,vegtyp , snowh  , tv     , lat   , yearlen , julian , & !in
                        lai    , sai    , troot  , elai    , esai   , igs)

! --------------------------------------------------------------------------------------------------
! vegetation phenology considering vegeation canopy being buries by snow and evolution in time
! --------------------------------------------------------------------------------------------------
  implicit none
! --------------------------------------------------------------------------------------------------
! inputs
  type (noahmp_parameters), intent(in) :: parameters
  integer                , intent(in   ) :: vegtyp !vegetation type 
  real                   , intent(in   ) :: snowh  !snow height [m]
  real                   , intent(in   ) :: tv     !vegetation temperature (k)
  real                   , intent(in   ) :: lat    !latitude (radians)
  integer                , intent(in   ) :: yearlen!number of days in the particular year
  real                   , intent(in   ) :: julian !julian day of year (fractional) ( 0 <= julian < yearlen )
  real                   , intent(in   ) :: troot  !root-zone averaged temperature (k)
  real                   , intent(inout) :: lai    !lai, unadjusted for burying by snow
  real                   , intent(inout) :: sai    !sai, unadjusted for burying by snow

! outputs
  real                   , intent(out  ) :: elai   !leaf area index, after burying by snow
  real                   , intent(out  ) :: esai   !stem area index, after burying by snow
  real                   , intent(out  ) :: igs    !growing season index (0=off, 1=on)

! locals

  real                                   :: db     !thickness of canopy buried by snow (m)
  real                                   :: fb     !fraction of canopy buried by snow
  real                                   :: snowhc !critical snow depth at which short vege
                                                   !is fully covered by snow

  integer                                :: k       !index
  integer                                :: it1,it2 !interpolation months
  real                                   :: day     !current day of year ( 0 <= day < yearlen )
  real                                   :: wt1,wt2 !interpolation weights
  real                                   :: t       !current month (1.00, ..., 12.00)
! --------------------------------------------------------------------------------------------------

  if ( dveg == 1 .or. dveg == 3 .or. dveg == 4 ) then

     if (lat >= 0.) then
        ! northern hemisphere
        day = julian
     else
        ! southern hemisphere.  day is shifted by 1/2 year.
        day = mod ( julian + ( 0.5 * yearlen ) , real(yearlen) )
     endif

     t = 12. * day / real(yearlen)
     it1 = t + 0.5
     it2 = it1 + 1
     wt1 = (it1+0.5) - t
     wt2 = 1.-wt1
     if (it1 .lt.  1) it1 = 12
     if (it2 .gt. 12) it2 = 1

     lai = wt1*parameters%laim(it1) + wt2*parameters%laim(it2)
     sai = wt1*parameters%saim(it1) + wt2*parameters%saim(it2)
  endif
  if (sai < 0.05) sai = 0.0                  ! mb: sai check, change to 0.05 v3.6
  if (lai < 0.05 .or. sai == 0.0) lai = 0.0  ! mb: lai check

  if ( ( vegtyp == parameters%iswater ) .or. ( vegtyp == parameters%isbarren ) .or. &
       ( vegtyp == parameters%isice   ) .or. ( parameters%urban_flag ) ) then
     lai  = 0.
     sai  = 0.
  endif

!buried by snow

     db = min( max(snowh - parameters%hvb,0.), parameters%hvt-parameters%hvb )
     fb = db / max(1.e-06,parameters%hvt-parameters%hvb)

     if(parameters%hvt> 0. .and. parameters%hvt <= 1.0) then          !mb: change to 1.0 and 0.2 to reflect
       snowhc = parameters%hvt*exp(-snowh/0.2)             !      changes to hvt in mptable
       fb     = min(snowh,snowhc)/snowhc
     endif

     elai =  lai*(1.-fb)
     esai =  sai*(1.-fb)
     if (esai < 0.05) esai = 0.0                   ! mb: esai check, change to 0.05 v3.6
     if (elai < 0.05 .or. esai == 0.0) elai = 0.0  ! mb: lai check

     if (tv .gt. parameters%tmin) then
         igs = 1.
     else
         igs = 0.
     endif

  end subroutine phenology

!== begin precip_heat ==============================================================================

!>\ingroup NoahMP_LSM
  subroutine precip_heat (parameters,iloc   ,jloc   ,vegtyp ,dt     ,uu     ,vv     , & !in
                          elai   ,esai   ,fveg   ,ist    ,                 & !in
                          bdfall ,rain   ,snow   ,fp     ,                 & !in
                          canliq ,canice ,tv     ,sfctmp ,tg     ,         & !in
                          qintr  ,qdripr ,qthror ,qints  ,qdrips ,qthros , & !out
			  pahv   ,pahg   ,pahb   ,qrain  ,qsnow  ,snowhin, & !out
			  fwet   ,cmc                                    )   !out

! ------------------------ code history ------------------------------
! michael barlage: oct 2013 - split canwater to calculate precip movement for 
!                             tracking of advected heat
! --------------------------------------------------------------------------------------------------
  implicit none
! ------------------------ input/output variables --------------------
! input
  type (noahmp_parameters), intent(in) :: parameters
  integer,intent(in)  :: iloc    !grid index
  integer,intent(in)  :: jloc    !grid index
  integer,intent(in)  :: vegtyp  !vegetation type
  integer,intent(in)  :: ist     !surface type 1-soil; 2-lake
  real,   intent(in)  :: dt      !main time step (s)
  real,   intent(in)  :: uu      !u-direction wind speed [m/s]
  real,   intent(in)  :: vv      !v-direction wind speed [m/s]
  real,   intent(in)  :: elai    !leaf area index, after burying by snow
  real,   intent(in)  :: esai    !stem area index, after burying by snow
  real,   intent(in)  :: fveg    !greeness vegetation fraction (-)
  real,   intent(in)  :: bdfall  !bulk density of snowfall (kg/m3)
  real,   intent(in)  :: rain    !rainfall (mm/s)
  real,   intent(in)  :: snow    !snowfall (mm/s)
  real,   intent(in)  :: fp      !fraction of the gridcell that receives precipitation
  real,   intent(in)  :: tv      !vegetation temperature (k)
  real,   intent(in)  :: sfctmp  !model-level temperature (k)
  real,   intent(in)  :: tg      !ground temperature (k)

! input & output
  real, intent(inout) :: canliq  !intercepted liquid water (mm)
  real, intent(inout) :: canice  !intercepted ice mass (mm)

! output
  real, intent(out)   :: qintr   !interception rate for rain (mm/s)
  real, intent(out)   :: qdripr  !drip rate for rain (mm/s)
  real, intent(out)   :: qthror  !throughfall for rain (mm/s)
  real, intent(out)   :: qints   !interception (loading) rate for snowfall (mm/s)
  real, intent(out)   :: qdrips  !drip (unloading) rate for intercepted snow (mm/s)
  real, intent(out)   :: qthros  !throughfall of snowfall (mm/s)
  real, intent(out)   :: pahv    !precipitation advected heat - vegetation net (w/m2)
  real, intent(out)   :: pahg    !precipitation advected heat - under canopy net (w/m2)
  real, intent(out)   :: pahb    !precipitation advected heat - bare ground net (w/m2)
  real, intent(out)   :: qrain   !rain at ground srf (mm/s) [+]
  real, intent(out)   :: qsnow   !snow at ground srf (mm/s) [+]
  real, intent(out)   :: snowhin !snow depth increasing rate (m/s)
  real, intent(out)   :: fwet    !wetted or snowed fraction of the canopy (-)
  real, intent(out)   :: cmc     !intercepted water (mm)
! --------------------------------------------------------------------

! ------------------------ local variables ---------------------------
  real                :: maxsno  !canopy capacity for snow interception (mm)
  real                :: maxliq  !canopy capacity for rain interception (mm)
  real                :: ft      !temperature factor for unloading rate
  real                :: fv      !wind factor for unloading rate
  real                :: pah_ac  !precipitation advected heat - air to canopy (w/m2)
  real                :: pah_cg  !precipitation advected heat - canopy to ground (w/m2)
  real                :: pah_ag  !precipitation advected heat - air to ground (w/m2)
  real                :: icedrip !canice unloading
! --------------------------------------------------------------------
! initialization

      qintr   = 0.
      qdripr  = 0.
      qthror  = 0.
      qintr   = 0.
      qints   = 0.
      qdrips  = 0.
      qthros  = 0.
      pah_ac  = 0.
      pah_cg  = 0.
      pah_ag  = 0.
      pahv    = 0.
      pahg    = 0.
      pahb    = 0.
      qrain   = 0.0
      qsnow   = 0.0
      snowhin = 0.0
      icedrip = 0.0
!      print*, "precip_heat begin canopy balance:",canliq+canice+(rain+snow)*dt
!      print*,  "precip_heat snow*3600.0:",snow*3600.0
!      print*,  "precip_heat rain*3600.0:",rain*3600.0
!      print*,  "precip_heat canice:",canice
!      print*,  "precip_heat canliq:",canliq

! --------------------------- liquid water ------------------------------
! maximum canopy water

      maxliq =  parameters%ch2op * (elai+ esai)

! average interception and throughfall

      if((elai+ esai).gt.0.) then
         qintr  = fveg * rain * fp  ! interception capability
         qintr  = min(qintr, (maxliq - canliq)/dt * (1.-exp(-rain*dt/maxliq)) )
         qintr  = max(qintr, 0.)
         qdripr = fveg * rain - qintr
         qthror = (1.-fveg) * rain
         canliq=max(0.,canliq+qintr*dt)
      else
         qintr  = 0.
         qdripr = 0.
         qthror = rain
	 if(canliq > 0.) then             ! for case of canopy getting buried
	   qdripr = qdripr + canliq/dt
	   canliq = 0.0
	 end if
      end if
      
! heat transported by liquid water

      pah_ac = fveg * rain * (cwat/1000.0) * (sfctmp - tv)
      pah_cg = qdripr * (cwat/1000.0) * (tv - tg)
      pah_ag = qthror * (cwat/1000.0) * (sfctmp - tg)
!      print*, "precip_heat pah_ac:",pah_ac
!      print*, "precip_heat pah_cg:",pah_cg
!      print*, "precip_heat pah_ag:",pah_ag

! --------------------------- canopy ice ------------------------------
! for canopy ice

      maxsno = 6.6*(0.27+46./bdfall) * (elai+ esai)

      if((elai+ esai).gt.0.) then
         qints = fveg * snow * fp
         qints = min(qints, (maxsno - canice)/dt * (1.-exp(-snow*dt/maxsno)) )
         qints = max(qints, 0.)
         ft = max(0.0,(tv - 270.15) / 1.87e5)
         fv = sqrt(uu*uu + vv*vv) / 1.56e5
	 ! mb: changed below to reflect the rain assumption that all precip gets intercepted 
	 icedrip = max(0.,canice) * (fv+ft)    !mb: removed /dt
         qdrips = (fveg * snow - qints) + icedrip
         qthros = (1.0-fveg) * snow
         canice= max(0.,canice + (qints - icedrip)*dt)
      else
         qints  = 0.
         qdrips = 0.
         qthros = snow
	 if(canice > 0.) then             ! for case of canopy getting buried
	   qdrips = qdrips + canice/dt
	   canice = 0.0
	 end if
      endif
!      print*, "precip_heat canopy through:",3600.0*(fveg * snow - qints)
!      print*, "precip_heat canopy drip:",3600.0*max(0.,canice) * (fv+ft)

! wetted fraction of canopy

      if(canice.gt.0.) then
           fwet = max(0.,canice) / max(maxsno,1.e-06)
      else
           fwet = max(0.,canliq) / max(maxliq,1.e-06)
      endif
      fwet = min(fwet, 1.) ** 0.667

! total canopy water

      cmc = canliq + canice

! heat transported by snow/ice

      pah_ac = pah_ac +  fveg * snow * (cice/1000.0) * (sfctmp - tv)
      pah_cg = pah_cg + qdrips * (cice/1000.0) * (tv - tg)
      pah_ag = pah_ag + qthros * (cice/1000.0) * (sfctmp - tg)
      
      pahv = pah_ac - pah_cg
      pahg = pah_cg
      pahb = pah_ag
      
      if (fveg > 0.0 .and. fveg < 1.0) then
        pahg = pahg / fveg         ! these will be multiplied by fraction later
	pahb = pahb / (1.0-fveg)
      elseif (fveg <= 0.0) then
        pahb = pahg + pahb         ! for case of canopy getting buried
        pahg = 0.0
	pahv = 0.0
      elseif (fveg >= 1.0) then
	pahb = 0.0
      end if
      
      pahv = max(pahv,-20.0)       ! put some artificial limits here for stability
      pahv = min(pahv,20.0)
      pahg = max(pahg,-20.0)
      pahg = min(pahg,20.0)
      pahb = max(pahb,-20.0)
      pahb = min(pahb,20.0)
      
!      print*, 'precip_heat sfctmp,tv,tg:',sfctmp,tv,tg
!      print*, 'precip_heat 3600.0*qints+qdrips+qthros:',3600.0*(qints+qdrips+qthros)
!      print*, "precip_heat maxsno:",maxsno
!      print*, "precip_heat pah_ac:",pah_ac
!      print*, "precip_heat pah_cg:",pah_cg
!      print*, "precip_heat pah_ag:",pah_ag
      
!      print*, "precip_heat pahv:",pahv
!      print*, "precip_heat pahg:",pahg
!      print*, "precip_heat pahb:",pahb
!      print*, "precip_heat fveg:",fveg
!      print*,  "precip_heat qints*3600.0:",qints*3600.0
!      print*,  "precip_heat qdrips*3600.0:",qdrips*3600.0
!      print*,  "precip_heat qthros*3600.0:",qthros*3600.0
      
! rain or snow on the ground

      qrain   = qdripr + qthror
      qsnow   = qdrips + qthros
      snowhin = qsnow/bdfall

      if (ist == 2 .and. tg > tfrz) then
         qsnow   = 0.
         snowhin = 0.
      end if
!      print*,  "precip_heat qsnow*3600.0:",qsnow*3600.0
!      print*,  "precip_heat qrain*3600.0:",qrain*3600.0
!      print*,  "precip_heat snowhin:",snowhin
!      print*,  "precip_heat canice:",canice
!      print*,  "precip_heat canliq:",canliq
!      print*, "precip_heat end canopy balance:",canliq+canice+(qrain+qsnow)*dt
      

  end subroutine precip_heat

!== begin error ====================================================================================

!>\ingroup NoahMP_LSM
  subroutine error (parameters,swdown ,fsa    ,fsr    ,fira   ,fsh    ,fcev   , &
                    fgev   ,fctr   ,ssoil  ,beg_wb ,canliq ,canice , &
                    sneqv  ,wa     ,smc    ,dzsnso ,prcp   ,ecan   , &
                    etran  ,edir   ,runsrf ,runsub ,dt     ,nsoil  , &
                    nsnow  ,ist    ,errwat, iloc   ,jloc   ,fveg   , &
                    sav    ,sag    ,fsrv   ,fsrg   ,zwt    ,pah    , &
#ifdef CCPP
                    pahv   ,pahg   ,pahb   ,errmsg, errflg)
#else
                    pahv   ,pahg   ,pahb   )
#endif
! --------------------------------------------------------------------------------------------------
! check surface energy balance and water balance
! --------------------------------------------------------------------------------------------------
  implicit none
! --------------------------------------------------------------------------------------------------
! inputs
  type (noahmp_parameters), intent(in) :: parameters
  integer                        , intent(in) :: nsnow  !maximum no. of snow layers        
  integer                        , intent(in) :: nsoil  !number of soil layers
  integer                        , intent(in) :: ist    !surface type 1->soil; 2->lake
  integer                        , intent(in) :: iloc   !grid index
  integer                        , intent(in) :: jloc   !grid index
  real                           , intent(in) :: swdown !downward solar filtered by sun angle [w/m2]
  real                           , intent(in) :: fsa    !total absorbed solar radiation (w/m2)
  real                           , intent(in) :: fsr    !total reflected solar radiation (w/m2)
  real                           , intent(in) :: fira   !total net longwave rad (w/m2)  [+ to atm]
  real                           , intent(in) :: fsh    !total sensible heat (w/m2)     [+ to atm]
  real                           , intent(in) :: fcev   !canopy evaporation heat (w/m2) [+ to atm]
  real                           , intent(in) :: fgev   !ground evaporation heat (w/m2) [+ to atm]
  real                           , intent(in) :: fctr   !transpiration heat flux (w/m2) [+ to atm]
  real                           , intent(in) :: ssoil  !ground heat flux (w/m2)        [+ to soil]
  real                           , intent(in) :: fveg
  real                           , intent(in) :: sav
  real                           , intent(in) :: sag
  real                           , intent(in) :: fsrv
  real                           , intent(in) :: fsrg
  real                           , intent(in) :: zwt

  real                           , intent(in) :: prcp   !precipitation rate (kg m-2 s-1)
  real                           , intent(in) :: ecan   !evaporation of intercepted water (mm/s)
  real                           , intent(in) :: etran  !transpiration rate (mm/s)
  real                           , intent(in) :: edir   !soil surface evaporation rate[mm/s]
  real                           , intent(in) :: runsrf !surface runoff [mm/s] 
  real                           , intent(in) :: runsub !baseflow (saturation excess) [mm/s]
  real                           , intent(in) :: canliq !intercepted liquid water (mm)
  real                           , intent(in) :: canice !intercepted ice mass (mm)
  real                           , intent(in) :: sneqv  !snow water eqv. [mm]
  real, dimension(       1:nsoil), intent(in) :: smc    !soil moisture (ice + liq.) [m3/m3]
  real, dimension(-nsnow+1:nsoil), intent(in) :: dzsnso !snow/soil layer thickness [m]
  real                           , intent(in) :: wa     !water storage in aquifer [mm]
  real                           , intent(in) :: dt     !time step [sec]
  real                           , intent(in) :: beg_wb !water storage at begin of a timesetp [mm]
  real                           , intent(out) :: errwat !error in water balance [mm/timestep]
  real, intent(in)   :: pah     !precipitation advected heat - total (w/m2)
  real, intent(in)   :: pahv    !precipitation advected heat - total (w/m2)
  real, intent(in)   :: pahg    !precipitation advected heat - total (w/m2)
  real, intent(in)   :: pahb    !precipitation advected heat - total (w/m2)

#ifdef CCPP
  character(len=*)               , intent(inout) :: errmsg
  integer                        , intent(inout) :: errflg
#endif

  integer                                     :: iz     !do-loop index
  real                                        :: end_wb !water storage at end of a timestep [mm]
  !kwm real                                        :: errwat !error in water balance [mm/timestep]
  real                                        :: erreng !error in surface energy balance [w/m2]
  real                                        :: errsw  !error in shortwave radiation balance [w/m2]
  real                                        :: fsrvg
  character(len=256)                          :: message
! --------------------------------------------------------------------------------------------------
!jref:start
   errsw   = swdown - (fsa + fsr)
!   errsw   = swdown - (sav+sag + fsrv+fsrg)
!   write(*,*) "errsw =",errsw
   if (abs(errsw) > 0.01) then            ! w/m2
   write(*,*) "vegetation!"
   write(*,*) "swdown*fveg =",swdown*fveg
   write(*,*) "fveg*(sav+sag) =",fveg*sav + sag
   write(*,*) "fveg*(fsrv +fsrg)=",fveg*fsrv + fsrg
   write(*,*) "ground!"
   write(*,*) "(1-.fveg)*swdown =",(1.-fveg)*swdown
   write(*,*) "(1.-fveg)*sag =",(1.-fveg)*sag
   write(*,*) "(1.-fveg)*fsrg=",(1.-fveg)*fsrg
   write(*,*) "fsrv   =",fsrv
   write(*,*) "fsrg   =",fsrg
   write(*,*) "fsr    =",fsr
   write(*,*) "sav    =",sav
   write(*,*) "sag    =",sag
   write(*,*) "fsa    =",fsa
!jref:end   
      write(message,*) 'errsw =',errsw
#ifdef CCPP
      errflg = 1
      errmsg = trim(message)//NEW_LINE('A')//"stop in noah-mp"
      return 
#else
      call wrf_message(trim(message))
      call wrf_error_fatal("stop in noah-mp")
#endif
   end if

   erreng = sav+sag-(fira+fsh+fcev+fgev+fctr+ssoil) +pah
!   erreng = fveg*sav+sag-(fira+fsh+fcev+fgev+fctr+ssoil)
   if(abs(erreng) > 0.01) then
      write(message,*) 'erreng =',erreng,' at i,j: ',iloc,jloc
#ifdef CCPP
      errmsg = trim(message)
#else
      call wrf_message(trim(message))
#endif
      write(message,'(a17,f10.4)') "net solar:       ",fsa
#ifdef CCPP
      errmsg = trim(errmsg)//NEW_LINE('A')//trim(message)
#else
      call wrf_message(trim(message))
#endif
      write(message,'(a17,f10.4)') "net longwave:    ",fira
#ifdef CCPP
      errmsg = trim(errmsg)//NEW_LINE('A')//trim(message)
#else
      call wrf_message(trim(message))
#endif
      write(message,'(a17,f10.4)') "total sensible:  ",fsh
#ifdef CCPP
      errmsg = trim(errmsg)//NEW_LINE('A')//trim(message)
#else
      call wrf_message(trim(message))
#endif
      write(message,'(a17,f10.4)') "canopy evap:     ",fcev
#ifdef CCPP
      errmsg = trim(errmsg)//NEW_LINE('A')//trim(message)
#else
      call wrf_message(trim(message))
#endif
      write(message,'(a17,f10.4)') "ground evap:     ",fgev
#ifdef CCPP
      errmsg = trim(errmsg)//NEW_LINE('A')//trim(message)
#else
      call wrf_message(trim(message))
#endif
      write(message,'(a17,f10.4)') "transpiration:   ",fctr
#ifdef CCPP
      errmsg = trim(errmsg)//NEW_LINE('A')//trim(message)
#else
      call wrf_message(trim(message))
#endif
      write(message,'(a17,f10.4)') "total ground:    ",ssoil
#ifdef CCPP
      errmsg = trim(errmsg)//NEW_LINE('A')//trim(message)
#else
      call wrf_message(trim(message))
#endif
      write(message,'(a17,4f10.4)') "precip advected: ",pah,pahv,pahg,pahb
#ifdef CCPP
      errmsg = trim(errmsg)//NEW_LINE('A')//trim(message)
#else
      call wrf_message(trim(message))
#endif
      write(message,'(a17,f10.4)') "precip: ",prcp
#ifdef CCPP
      errmsg = trim(errmsg)//NEW_LINE('A')//trim(message)
#else
      call wrf_message(trim(message))
#endif
      write(message,'(a17,f10.4)') "veg fraction: ",fveg
#ifdef CCPP
      errflg = 1
      errmsg = trim(errmsg)//NEW_LINE('A')//trim(message)//NEW_LINE('A')//"energy budget problem in noahmp lsm"
      return
#else
      call wrf_message(trim(message))
      call wrf_error_fatal("energy budget problem in noahmp lsm")
#endif
      
   end if

   if (ist == 1) then                                       !soil
        end_wb = canliq + canice + sneqv + wa
        do iz = 1,nsoil
          end_wb = end_wb + smc(iz) * dzsnso(iz) * 1000.
        end do
        errwat = end_wb-beg_wb-(prcp-ecan-etran-edir-runsrf-runsub)*dt

   else                 !kwm
      errwat = 0.0      !kwm
   endif

 end subroutine error

!== begin energy ===================================================================================

!>\ingroup NoahMP_LSM
  subroutine energy (parameters,ice    ,vegtyp ,ist    ,nsnow  ,nsoil  , & !in
                     isnow  ,dt     ,rhoair ,sfcprs ,qair   , & !in
                     sfctmp ,thair  ,lwdn   ,uu     ,vv     ,zref   , & !in
                     lheatstrg      ,  & !in
                     co2air ,o2air  ,solad  ,solai  ,cosz   ,igs    , & !in
                     eair   ,tbot   ,zsnso  ,zsoil  , & !in
                     elai   ,esai   ,fwet   ,foln   ,         & !in
                     fveg   ,pahv   ,pahg   ,pahb   ,                 & !in
                     qsnow  ,dzsnso ,lat    ,canliq ,canice ,iloc   , jloc, & !in
		     z0wrf  ,                                         &
                     imelt  ,snicev ,snliqv ,epore  ,t2m    ,fsno   , & !out
                     sav    ,sag    ,qmelt  ,fsa    ,fsr    ,taux   , & !out
                     tauy   ,fira   ,fsh    ,fshx   ,fcev   ,fgev   ,fctr   , & !out
                     trad   ,psn    ,apar   ,ssoil  ,btrani ,btran  , & !out
                     ponding,ts     ,latheav , latheag , frozen_canopy,frozen_ground,                       & !out
                     tv     ,tg     ,stc    ,snowh  ,eah    ,tah    , & !inout
                     sneqvo ,sneqv  ,sh2o   ,smc    ,snice  ,snliq  , & !inout
                     albold ,cm     ,ch     ,dx     ,dz8w   ,q2     , &   !inout
#ifdef CCPP
                     tauss  ,cpfac  ,errmsg ,errflg,                  & !inout
#else
                     tauss  ,cpfac  ,                                 & !inout
#endif
!jref:start
                     qc     ,qsfc   ,psfc   , & !in 
                     t2mv   ,t2mb   ,fsrv   , &
                     fsrg   ,rssun  ,rssha  ,bgap   ,wgap,tgv,tgb,&
                     q1     ,q2v    ,q2b    ,q2e    ,chv  ,chb, emissi,pah  ,&
		     shg,shc,shb,evg,evb,ghv,ghb,irg,irc,irb,tr,evc,chleaf,chuc,chv2,chb2 )   !out 
!jref:end                            

! --------------------------------------------------------------------------------------------------
! we use different approaches to deal with subgrid features of radiation transfer and turbulent
! transfer. we use 'tile' approach to compute turbulent fluxes, while we use modified two-
! stream to compute radiation transfer. tile approach, assemblying vegetation canopies together,
! may expose too much ground surfaces (either covered by snow or grass) to solar radiation. the
! modified two-stream assumes vegetation covers fully the gridcell but with gaps between tree
! crowns.
! --------------------------------------------------------------------------------------------------
! turbulence transfer : 'tile' approach to compute energy fluxes in vegetated fraction and
!                         bare fraction separately and then sum them up weighted by fraction
!                     --------------------------------------
!                    / o  o  o  o  o  o  o  o  /          / 
!                   /  |  |  |  |  |  |  |  | /          /
!                  / o  o  o  o  o  o  o  o  /          /
!                 /  |  |  |tile1|  |  |  | /  tile2   /
!                / o  o  o  o  o  o  o  o  /  bare    /
!               /  |  |  | vegetated |  | /          /
!              / o  o  o  o  o  o  o  o  /          /
!             /  |  |  |  |  |  |  |  | /          /
!            --------------------------------------
! --------------------------------------------------------------------------------------------------
! radiation transfer : modified two-stream (yang and friedl, 2003, jgr; niu ang yang, 2004, jgr)
!                     --------------------------------------  two-stream treats leaves as
!                    /   o   o   o   o   o   o   o   o    /  cloud over the entire grid-cell,
!                   /    |   |   |   |   |   |   |   |   / while the modified two-stream 
!                  /   o   o   o   o   o   o   o   o    / aggregates cloudy leaves into  
!                 /    |   |   |   |   |   |   |   |   / tree crowns with gaps (as shown in
!                /   o   o   o   o   o   o   o   o    / the left figure). we assume these
!               /    |   |   |   |   |   |   |   |   / tree crowns are evenly distributed
!              /   o   o   o   o   o   o   o   o    / within the gridcell with 100% veg
!             /    |   |   |   |   |   |   |   |   / fraction, but with gaps. the 'tile'
!            -------------------------------------- approach overlaps too much shadows.
! --------------------------------------------------------------------------------------------------
  implicit none
! --------------------------------------------------------------------------------------------------
! inputs
  type (noahmp_parameters), intent(in) :: parameters
  integer                           , intent(in)    :: iloc
  integer                           , intent(in)    :: jloc
  integer                           , intent(in)    :: ice    !ice (ice = 1)
  integer                           , intent(in)    :: vegtyp !vegetation physiology type
  integer                           , intent(in)    :: ist    !surface type: 1->soil; 2->lake
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
  real                              , intent(in)    :: thair  !potential temperature (k)
  real                              , intent(in)    :: lwdn   !downward longwave radiation (w/m2)
  real                              , intent(in)    :: uu     !wind speed in e-w dir (m/s)
  real                              , intent(in)    :: vv     !wind speed in n-s dir (m/s)
  real   , dimension(       1:    2), intent(in)    :: solad  !incoming direct solar rad. (w/m2)
  real   , dimension(       1:    2), intent(in)    :: solai  !incoming diffuse solar rad. (w/m2)
  real                              , intent(in)    :: cosz   !cosine solar zenith angle (0-1)
  real                              , intent(in)    :: elai   !lai adjusted for burying by snow
  real                              , intent(in)    :: esai   !lai adjusted for burying by snow
  real                              , intent(in)    :: fwet   !fraction of canopy that is wet [-]
  real                              , intent(in)    :: fveg   !greeness vegetation fraction (-)
  real                              , intent(in)    :: lat    !latitude (radians)
  real                              , intent(in)    :: canliq !canopy-intercepted liquid water (mm)
  real                              , intent(in)    :: canice !canopy-intercepted ice mass (mm)
  real                              , intent(in)    :: foln   !foliage nitrogen (%)
  real                              , intent(in)    :: co2air !atmospheric co2 concentration (pa)
  real                              , intent(in)    :: o2air  !atmospheric o2 concentration (pa)
  real                              , intent(in)    :: igs    !growing season index (0=off, 1=on)

  real                              , intent(in)    :: zref   !reference height (m)
  logical                           , intent(in)    :: lheatstrg ! flag for canopy heat storage parameterization       
  real                              , intent(in)    :: tbot   !bottom condition for soil temp. (k) 
  real   , dimension(-nsnow+1:nsoil), intent(in)    :: zsnso  !layer-bottom depth from snow surf [m]
  real   , dimension(       1:nsoil), intent(in)    :: zsoil  !layer-bottom depth from soil surf [m]
  real   , dimension(-nsnow+1:nsoil), intent(in)    :: dzsnso !depth of snow & soil layer-bottom [m]
  real, intent(in)   :: pahv    !precipitation advected heat - vegetation net (w/m2)
  real, intent(in)   :: pahg    !precipitation advected heat - under canopy net (w/m2)
  real, intent(in)   :: pahb    !precipitation advected heat - bare ground net (w/m2)

!jref:start; in 
  real                              , intent(in)    :: qc     !cloud water mixing ratio
  real                              , intent(inout) :: qsfc   !mixing ratio at lowest model layer
  real                              , intent(in)    :: psfc   !pressure at lowest model layer
  real                              , intent(in)    :: dx     !horisontal resolution
  real                              , intent(in)    :: dz8w   !thickness of lowest layer
  real                              , intent(in)    :: q2     !mixing ratio (kg/kg)
!jref:end

! outputs
  real                              , intent(out)   :: z0wrf  !combined z0 sent to coupled model
  integer, dimension(-nsnow+1:nsoil), intent(out)   :: imelt  !phase change index [1-melt; 2-freeze]
  real   , dimension(-nsnow+1:    0), intent(out)   :: snicev !partial volume ice [m3/m3]
  real   , dimension(-nsnow+1:    0), intent(out)   :: snliqv !partial volume liq. water [m3/m3]
  real   , dimension(-nsnow+1:    0), intent(out)   :: epore  !effective porosity [m3/m3]
  real                              , intent(out)   :: fsno   !snow cover fraction (-)
  real                              , intent(out)   :: qmelt  !snowmelt [mm/s]
  real                              , intent(out)   :: ponding!pounding at ground [mm]
  real                              , intent(out)   :: sav    !solar rad. absorbed by veg. (w/m2)
  real                              , intent(out)   :: sag    !solar rad. absorbed by ground (w/m2)
  real                              , intent(out)   :: fsa    !tot. absorbed solar radiation (w/m2)
  real                              , intent(out)   :: fsr    !tot. reflected solar radiation (w/m2)
  real                              , intent(out)   :: taux   !wind stress: e-w (n/m2)
  real                              , intent(out)   :: tauy   !wind stress: n-s (n/m2)
  real                              , intent(out)   :: fira   !total net lw. rad (w/m2)   [+ to atm]
  real                              , intent(out)   :: fsh    !total sensible heat (w/m2) [+ to atm]
  real                              , intent(out)   :: fshx   !total sensible heat (w/m2) [+ to atm]
  real                              , intent(out)   :: fcev   !canopy evaporation (w/m2)  [+ to atm]
  real                              , intent(out)   :: fgev   !ground evaporation (w/m2)  [+ to atm]
  real                              , intent(out)   :: fctr   !transpiration (w/m2)       [+ to atm]
  real                              , intent(out)   :: trad   !radiative temperature (k)
  real                              , intent(out)   :: t2m    !2 m height air temperature (k)
  real                              , intent(out)   :: psn    !total photosyn. (umolco2/m2/s) [+]
  real                              , intent(out)   :: apar   !total photosyn. active energy (w/m2)
  real                              , intent(out)   :: ssoil  !ground heat flux (w/m2)   [+ to soil]
  real   , dimension(       1:nsoil), intent(out)   :: btrani !soil water transpiration factor (0-1)
  real                              , intent(out)   :: btran  !soil water transpiration factor (0-1)
!  real                              , intent(out)   :: lathea !latent heat vap./sublimation (j/kg)
  real                              , intent(out)   :: latheav !latent heat vap./sublimation (j/kg)
  real                              , intent(out)   :: latheag !latent heat vap./sublimation (j/kg)
  logical                           , intent(out)   :: frozen_ground ! used to define latent heat pathway
  logical                           , intent(out)   :: frozen_canopy ! used to define latent heat pathway

!jref:start  
  real                              , intent(out)   :: fsrv    !veg. reflected solar radiation (w/m2)
  real                              , intent(out)   :: fsrg    !ground reflected solar radiation (w/m2)
  real, intent(out) :: rssun        !sunlit leaf stomatal resistance (s/m)
  real, intent(out) :: rssha        !shaded leaf stomatal resistance (s/m)
!jref:end - out for debug  

!jref:start; output
  real                              , intent(out)   :: t2mv   !2-m air temperature over vegetated part [k]
  real                              , intent(out)   :: t2mb   !2-m air temperature over bare ground part [k]
  real                              , intent(out)   :: bgap
  real                              , intent(out)   :: wgap
!jref:end

! input & output
  real                              , intent(inout) :: ts     !surface temperature (k)
  real                              , intent(inout) :: tv     !vegetation temperature (k)
  real                              , intent(inout) :: tg     !ground temperature (k)
  real   , dimension(-nsnow+1:nsoil), intent(inout) :: stc    !snow/soil temperature [k]
  real                              , intent(inout) :: snowh  !snow height [m]
  real                              , intent(inout) :: sneqv  !snow mass (mm)
  real                              , intent(inout) :: sneqvo !snow mass at last time step (mm)
  real   , dimension(       1:nsoil), intent(inout) :: sh2o   !liquid soil moisture [m3/m3]
  real   , dimension(       1:nsoil), intent(inout) :: smc    !soil moisture (ice + liq.) [m3/m3]
  real   , dimension(-nsnow+1:    0), intent(inout) :: snice  !snow ice mass (kg/m2)
  real   , dimension(-nsnow+1:    0), intent(inout) :: snliq  !snow liq mass (kg/m2)
  real                              , intent(inout) :: eah    !canopy air vapor pressure (pa)
  real                              , intent(inout) :: tah    !canopy air temperature (k)
  real                              , intent(inout) :: albold !snow albedo at last time step(class type)
  real                              , intent(inout) :: tauss  !non-dimensional snow age
  real                              , intent(inout) :: cpfac  !heat capacity enhancement factor due to heat storage
  real                              , intent(inout) :: cm     !momentum drag coefficient
  real                              , intent(inout) :: ch     !sensible heat exchange coefficient
  real                              , intent(inout) :: q1
#ifdef CCPP
  character(len=*)                  , intent(inout) :: errmsg
  integer                           , intent(inout) :: errflg
#endif
!  real                                              :: q2e
  real,                               intent(out)   :: emissi
  real,                               intent(out)   :: pah    !precipitation advected heat - total (w/m2)

! local
  integer                                           :: iz     !do-loop index
  logical                                           :: veg    !true if vegetated surface
  real                                              :: ur     !wind speed at height zlvl (m/s)
  real                                              :: zlvl   !reference height (m)
  real                                              :: fsun   !sunlit fraction of canopy [-]
  real                                              :: rb     !leaf boundary layer resistance (s/m)
  real                                              :: rsurf  !ground surface resistance (s/m)
  real                                              :: l_rsurf!dry-layer thickness for computing rsurf (sakaguchi and zeng, 2009)
  real                                              :: d_rsurf!reduced vapor diffusivity in soil for computing rsurf (sz09)
  real                                              :: bevap  !soil water evaporation factor (0- 1)
  real                                              :: mol    !monin-obukhov length (m)
  real                                              :: vai    !sum of lai  + stem area index [m2/m2]
  real                                              :: cwp    !canopy wind extinction parameter
  real                                              :: zpd    !zero plane displacement (m)
  real                                              :: z0m    !z0 momentum (m)
  real                                              :: zpdg   !zero plane displacement (m)
  real                                              :: z0mg   !z0 momentum, ground (m)
  real                                              :: emv    !vegetation emissivity
  real                                              :: emg    !ground emissivity
  real                                              :: fire   !emitted ir (w/m2)

  real                                              :: laisun !sunlit leaf area index (m2/m2)
  real                                              :: laisha !shaded leaf area index (m2/m2)
  real                                              :: psnsun !sunlit photosynthesis (umolco2/m2/s)
  real                                              :: psnsha !shaded photosynthesis (umolco2/m2/s)
!jref:start - for debug  
!  real                                              :: rssun  !sunlit stomatal resistance (s/m)
!  real                                              :: rssha  !shaded stomatal resistance (s/m)
!jref:end - for debug
  real                                              :: parsun !par absorbed per sunlit lai (w/m2)
  real                                              :: parsha !par absorbed per shaded lai (w/m2)

  real, dimension(-nsnow+1:nsoil)                   :: fact   !temporary used in phase change
  real, dimension(-nsnow+1:nsoil)                   :: df     !thermal conductivity [w/m/k]
  real, dimension(-nsnow+1:nsoil)                   :: hcpct  !heat capacity [j/m3/k]
  real                                              :: bdsno  !bulk density of snow (kg/m3)
  real                                              :: fmelt  !melting factor for snow cover frac
  real                                              :: gx     !temporary variable
  real, dimension(-nsnow+1:nsoil)                   :: phi    !light through water (w/m2)
!  real                                              :: gamma  !psychrometric constant (pa/k)
  real                                              :: gammav  !psychrometric constant (pa/k)
  real                                              :: gammag  !psychrometric constant (pa/k)
  real                                              :: psi    !surface layer soil matrix potential (m)
  real                                              :: rhsur  !raltive humidity in surface soil/snow air space (-)

! temperature and fluxes over vegetated fraction

  real                                              :: tauxv  !wind stress: e-w dir [n/m2]
  real                                              :: tauyv  !wind stress: n-s dir [n/m2]
  real,intent(out)                                              :: irc    !canopy net lw rad. [w/m2] [+ to atm]
  real,intent(out)                                              :: irg    !ground net lw rad. [w/m2] [+ to atm]
  real,intent(out)                                              :: shc    !canopy sen. heat [w/m2]   [+ to atm]
  real,intent(out)                                              :: shg    !ground sen. heat [w/m2]   [+ to atm]
!jref:start  
  real,intent(out)                                  :: q2v
  real,intent(out)                                  :: q2b
  real,intent(out)                                  :: q2e
!jref:end  
  real,intent(out)                                              :: evc    !canopy evap. heat [w/m2]  [+ to atm]
  real,intent(out)                                              :: evg    !ground evap. heat [w/m2]  [+ to atm]
  real,intent(out)                                              :: tr     !transpiration heat [w/m2] [+ to atm]
  real,intent(out)                                              :: ghv    !ground heat flux [w/m2]  [+ to soil]
  real,intent(out)                                  :: tgv    !ground surface temp. [k]
  real                                              :: cmv    !momentum drag coefficient
  real,intent(out)                                  :: chv    !sensible heat exchange coefficient

! temperature and fluxes over bare soil fraction

  real                                              :: tauxb  !wind stress: e-w dir [n/m2]
  real                                              :: tauyb  !wind stress: n-s dir [n/m2]
  real,intent(out)                                              :: irb    !net longwave rad. [w/m2] [+ to atm]
  real,intent(out)                                              :: shb    !sensible heat [w/m2]     [+ to atm]
  real,intent(out)                                              :: evb    !evaporation heat [w/m2]  [+ to atm]
  real,intent(out)                                              :: ghb    !ground heat flux [w/m2] [+ to soil]
  real,intent(out)                                  :: tgb    !ground surface temp. [k]
  real                                              :: cmb    !momentum drag coefficient
  real,intent(out)                                  :: chb    !sensible heat exchange coefficient
  real,intent(out)                                  :: chleaf !leaf exchange coefficient
  real,intent(out)                                  :: chuc   !under canopy exchange coefficient
!jref:start  
  real,intent(out)                                  :: chv2    !sensible heat conductance, canopy air to zlvl air (m/s)
  real,intent(out)                                  :: chb2    !sensible heat conductance, canopy air to zlvl air (m/s)
  real                                  :: noahmpres

!jref:end  

  real, parameter                   :: mpe    = 1.e-6
  real, parameter                   :: psiwlt = -150.  !metric potential for wilting point (m)
  real, parameter                   :: z0     = 0.01   ! bare-soil roughness length (m) (i.e., under the canopy)
!
! parameters for heat storage parametrization
!
  real, parameter :: z0min = 0.2 !minimum roughness length for heat storage
  real, parameter :: z0max = 1.0 !maximum roughness length for heat storage

! ---------------------------------------------------------------------------------------------------
! initialize fluxes from veg. fraction

    tauxv     = 0.    
    tauyv     = 0.
    irc       = 0.
    shc       = 0.
    irg       = 0.
    shg       = 0.
    evg       = 0.       
    evc       = 0.
    tr        = 0.
    ghv       = 0.       
    psnsun    = 0.
    psnsha    = 0.
    t2mv      = 0.
    q2v       = 0.
    chv       = 0.
    chleaf    = 0.
    chuc      = 0.
    chv2      = 0.

! wind speed at reference height: ur >= 1

    ur = max( sqrt(uu**2.+vv**2.), 1. )

! vegetated or non-vegetated

    vai = elai + esai
    veg = .false.
    if(vai > 0.) veg = .true.

! ground snow cover fraction [niu and yang, 2007, jgr]

     fsno = 0.
     if(snowh.gt.0.)  then
         bdsno    = sneqv / snowh
         fmelt    = (bdsno/100.)**parameters%mfsno
         fsno     = tanh( snowh /(2.5* z0 * fmelt))
     endif

! ground roughness length

     if(ist == 2) then
       if(tg .le. tfrz) then
         z0mg = 0.01 * (1.0-fsno) + fsno * parameters%z0sno
       else
         z0mg = 0.01  
       end if
     else
       z0mg = z0 * (1.0-fsno) + fsno * parameters%z0sno
     end if

! roughness length and displacement height

     zpdg  = snowh
     if(veg) then
        z0m  = parameters%z0mvt
        zpd  = 0.65 * parameters%hvt
        if(snowh.gt.zpd) zpd  = snowh
     else
        z0m  = z0mg
        zpd  = zpdg
     end if
!
!  compute heat capacity enhancement factor as a function of z0m to mimic heat storage
!
     if (lheatstrg .and. (.not. parameters%urban_flag) ) then
         cpfac = (z0m - z0min) / (z0max - z0min)
         cpfac = 1. + min(max(cpfac, 0.0), 1.0)
     endif

     zlvl = max(zpd,parameters%hvt) + zref
     if(zpdg >= zlvl) zlvl = zpdg + zref
!     ur   = ur*log(zlvl/z0m)/log(10./z0m)       !input ur is at 10m

! canopy wind absorption coeffcient

     cwp = parameters%cwpvt

! thermal properties of soil, snow, lake, and frozen soil

  call thermoprop (parameters,nsoil   ,nsnow   ,isnow   ,ist     ,dzsnso  , & !in
                   dt      ,snowh   ,snice   ,snliq   , & !in
                   smc     ,sh2o    ,tg      ,stc     ,ur      , & !in
                   lat     ,z0m     ,zlvl    ,vegtyp  , & !in
                   df      ,hcpct   ,snicev  ,snliqv  ,epore   , & !out
                   fact    )                              !out

! solar radiation: absorbed & reflected by the ground and canopy

  call  radiation (parameters,vegtyp  ,ist     ,ice     ,nsoil   , & !in 
                   sneqvo  ,sneqv   ,dt      ,cosz    ,snowh   , & !in
                   tg      ,tv      ,fsno    ,qsnow   ,fwet    , & !in
                   elai    ,esai    ,smc     ,solad   ,solai   , & !in
                   fveg    ,iloc    ,jloc    ,                   & !in
                   albold  ,tauss   ,                            & !inout
                   fsun    ,laisun  ,laisha  ,parsun  ,parsha  , & !out
                   sav     ,sag     ,fsr     ,fsa     ,fsrv    , & 
                   fsrg    ,bgap    ,wgap    )            !out

! vegetation and ground emissivity

     emv = 1. - exp(-(elai+esai)/1.0)
     if (ice == 1) then
       emg = 0.98*(1.-fsno) + 1.0*fsno
     else
       emg = parameters%eg(ist)*(1.-fsno) + 1.0*fsno
     end if

! soil moisture factor controlling stomatal resistance
   
     btran = 0.

     if(ist ==1 ) then
       do iz = 1, parameters%nroot
          if(opt_btr == 1) then                  ! noah
            gx    = (sh2o(iz)-parameters%smcwlt) / (parameters%smcref-parameters%smcwlt)
          end if
          if(opt_btr == 2) then                  ! clm
            psi   = max(psiwlt,-parameters%psisat*(max(0.01,sh2o(iz))/parameters%smcmax)**(-parameters%bexp) )
            gx    = (1.-psi/psiwlt)/(1.+parameters%psisat/psiwlt)
          end if
          if(opt_btr == 3) then                  ! ssib
            psi   = max(psiwlt,-parameters%psisat*(max(0.01,sh2o(iz))/parameters%smcmax)**(-parameters%bexp) )
            gx    = 1.-exp(-5.8*(log(psiwlt/psi))) 
          end if
       
          gx = min(1.,max(0.,gx))
          btrani(iz) = max(mpe,dzsnso(iz) / (-zsoil(parameters%nroot)) * gx)
          btran      = btran + btrani(iz)
       end do
       btran = max(mpe,btran)

       btrani(1:parameters%nroot) = btrani(1:parameters%nroot)/btran
     end if

! soil surface resistance for ground evap.

     bevap = max(0.0,sh2o(1)/parameters%smcmax)
     if(ist == 2) then
       rsurf = 1.          ! avoid being divided by 0
       rhsur = 1.0
     else

        ! rsurf based on sakaguchi and zeng, 2009
        ! taking the "residual water content" to be the wilting point, 
        ! and correcting the exponent on the d term (typo in sz09 ?)
        l_rsurf = (-zsoil(1)) * ( exp ( (1.0 - min(1.0,sh2o(1)/parameters%smcmax)) ** 5 ) - 1.0 ) / ( 2.71828 - 1.0 ) 
        d_rsurf = 2.2e-5 * parameters%smcmax * parameters%smcmax * ( 1.0 - parameters%smcwlt / parameters%smcmax ) ** (2.0+3.0/parameters%bexp)
        rsurf = l_rsurf / d_rsurf

        ! older rsurf computations:
        !    rsurf = fsno * 1. + (1.-fsno)* exp(8.25-4.225*bevap) !sellers (1992)
        !    rsurf = fsno * 1. + (1.-fsno)* exp(8.25-6.0  *bevap) !adjusted to decrease rsurf for wet soil

       if(sh2o(1) < 0.01 .and. snowh == 0.) rsurf = 1.e6
       psi   = -parameters%psisat*(max(0.01,sh2o(1))/parameters%smcmax)**(-parameters%bexp)   
       rhsur = fsno + (1.-fsno) * exp(psi*grav/(rw*tg)) 
     end if

! urban - jref 
     if (parameters%urban_flag .and. snowh == 0. ) then
        rsurf = 1.e6
     endif

! set psychrometric constant

     if (tv .gt. tfrz) then           ! barlage: add distinction between ground and 
        latheav = hvap                ! vegetation in v3.6
	frozen_canopy = .false.
     else
        latheav = hsub
	frozen_canopy = .true.
     end if
     gammav = cpair*cpfac*sfcprs/(0.622*latheav)

     if (tg .gt. tfrz) then
        latheag = hvap
	frozen_ground = .false.
     else
        latheag = hsub
	frozen_ground = .true.
     end if
     gammag = cpair*cpfac*sfcprs/(0.622*latheag)

!     if (sfctmp .gt. tfrz) then
!        lathea = hvap
!     else
!        lathea = hsub
!     end if
!     gamma = cpair*cpfac*sfcprs/(0.622*lathea)

! surface temperatures of the ground and canopy and energy fluxes

    if (veg .and. fveg > 0) then 
    tgv = tg
    cmv = cm
    chv = ch
! YRQ
!    write(*,*) 'cm,ch,tv,tgv, YRQ', cm,ch,tv,tgv
    call vege_flux (parameters,nsnow   ,nsoil   ,isnow   ,vegtyp  ,veg     , & !in
                    dt      ,sav     ,sag     ,lwdn    ,ur      , & !in
                    uu      ,vv      ,sfctmp  ,thair   ,qair    , & !in
                    eair    ,rhoair  ,snowh   ,vai     ,gammav  ,gammag    , & !in
                    fwet    ,laisun  ,laisha  ,cwp     ,dzsnso  , & !in
                    zlvl    ,cpfac   ,zpd     ,z0m     ,fveg    , & !in
                    z0mg    ,emv     ,emg     ,canliq  ,fsno, & !in
                    canice  ,stc     ,df      ,rssun   ,rssha   , & !in
                    rsurf   ,latheav ,latheag ,parsun  ,parsha  ,igs     , & !in
                    foln    ,co2air  ,o2air   ,btran   ,sfcprs  , & !in
                    rhsur   ,iloc    ,jloc    ,q2      ,pahv  ,pahg  , & !in
                    eah     ,tah     ,tv      ,tgv     ,cmv     , & !inout
#ifdef CCPP
                    chv     ,dx      ,dz8w    ,errmsg  ,errflg  , & !inout
#else
                    chv     ,dx      ,dz8w    ,                   & !inout
#endif
                    tauxv   ,tauyv   ,irg     ,irc     ,shg     , & !out
                    shc     ,evg     ,evc     ,tr      ,ghv     , & !out
                    t2mv    ,psnsun  ,psnsha  ,                   & !out
!jref:start
                    qc      ,qsfc    ,psfc    , & !in
                    q2v     ,chv2, chleaf, chuc)               !inout 
!jref:end
#ifdef CCPP
        if (errflg /= 0) return 
#endif                           
    end if

    tgb = tg
    cmb = cm
    chb = ch
    call bare_flux (parameters,nsnow   ,nsoil   ,isnow   ,dt      ,sag     , & !in
                    lwdn    ,ur      ,uu      ,vv      ,sfctmp  , & !in
                    thair   ,qair    ,eair    ,rhoair  ,snowh   , & !in
                    dzsnso  ,zlvl    ,zpdg    ,z0mg    ,fsno,             & !in
                    emg     ,stc     ,df      ,rsurf   ,latheag  , & !in
                    gammag   ,rhsur   ,iloc    ,jloc    ,q2      ,pahb  , & !in
#ifdef CCPP
                    tgb     ,cmb     ,chb     ,errmsg  ,errflg   , & !inout
#else
                    tgb     ,cmb     ,chb     ,                   & !inout
#endif
                    tauxb   ,tauyb   ,irb     ,shb     ,evb     , & !out
                    ghb     ,t2mb    ,dx      ,dz8w    ,vegtyp  , & !out
!jref:start
                    qc      ,qsfc    ,psfc    , & !in
                    sfcprs  ,q2b,   chb2)                          !in 
!jref:end
#ifdef CCPP
    if (errflg /= 0) return
#endif
!energy balance at vege canopy: sav          =(irc+shc+evc+tr)     *fveg  at   fveg 
!energy balance at vege ground: sag*    fveg =(irg+shg+evg+ghv)    *fveg  at   fveg
!energy balance at bare ground: sag*(1.-fveg)=(irb+shb+evb+ghb)*(1.-fveg) at 1-fveg

    if (veg .and. fveg > 0) then 
        taux  = fveg * tauxv     + (1.0 - fveg) * tauxb
        tauy  = fveg * tauyv     + (1.0 - fveg) * tauyb
        fira  = fveg * irg       + (1.0 - fveg) * irb       + irc
        fsh   = fveg * shg       + (1.0 - fveg) * shb       + shc
        fshx  = fveg * shg/cpfac + (1.0 - fveg) * shb + shc/cpfac
        fgev  = fveg * evg       + (1.0 - fveg) * evb
        ssoil = fveg * ghv       + (1.0 - fveg) * ghb
        fcev  = evc
        fctr  = tr
	pah   = fveg * pahg      + (1.0 - fveg) * pahb   + pahv
        tg    = fveg * tgv       + (1.0 - fveg) * tgb
        t2m   = fveg * t2mv      + (1.0 - fveg) * t2mb
        ts    = fveg * tv        + (1.0 - fveg) * tgb
        cm    = fveg * cmv       + (1.0 - fveg) * cmb      ! better way to average?
        ch    = fveg * chv       + (1.0 - fveg) * chb
        q1    = fveg * (eah*0.622/(sfcprs - 0.378*eah)) + (1.0 - fveg)*qsfc
        q2e   = fveg * q2v       + (1.0 - fveg) * q2b
	z0wrf = z0m
    else
        taux  = tauxb
        tauy  = tauyb
        fira  = irb
        fsh   = shb
        fshx  = shb
        fgev  = evb
        ssoil = ghb
        tg    = tgb
        t2m   = t2mb
        fcev  = 0.
        fctr  = 0.
	pah   = pahb
        ts    = tg
        cm    = cmb
        ch    = chb
        q1    = qsfc
        q2e   = q2b
        rssun = 0.0
        rssha = 0.0
        tgv   = tgb
        chv   = chb
	z0wrf = z0mg
    end if

    fire = lwdn + fira

    if(fire <=0.) then
       write(6,*) 'emitted longwave <0; skin t may be wrong due to inconsistent'
       write(6,*) 'input of shdfac with lai'
       write(6,*) iloc, jloc, 'shdfac=',fveg,'vai=',vai,'tv=',tv,'tg=',tg
       write(6,*) 'lwdn=',lwdn,'fira=',fira,'snowh=',snowh
#ifdef CCPP
      errflg = 1
      errmsg = "stop in noah-mp"
      return
#else
      call wrf_error_fatal("stop in noah-mp")
#endif
       
    end if

    ! compute a net emissivity
    emissi = fveg * ( emg*(1-emv) + emv + emv*(1-emv)*(1-emg) ) + &
         (1-fveg) * emg

    ! when we're computing a trad, subtract from the emitted ir the
    ! reflected portion of the incoming lwdn, so we're just
    ! considering the ir originating in the canopy/ground system.
    
    trad = ( ( fire - (1-emissi)*lwdn ) / (emissi*sb) ) ** 0.25

    ! old trad calculation not taking into account emissivity:
    ! trad = (fire/sb)**0.25

    apar = parsun*laisun + parsha*laisha
    psn  = psnsun*laisun + psnsha*laisha

! 3l snow & 4l soil temperatures

    call tsnosoi (parameters,ice     ,nsoil   ,nsnow   ,isnow   ,ist     , & !in
                  tbot    ,zsnso   ,ssoil   ,df      ,hcpct   , & !in
                  sag     ,dt      ,snowh   ,dzsnso  , & !in
                  tg      ,iloc    ,jloc    ,                   & !in
#ifdef CCPP
                  stc     ,errmsg  ,errflg     )                  !inout
#else
                  stc     )                                       !inout
#endif

#ifdef CCPP
    if (errflg /= 0) return
#endif

! adjusting snow surface temperature
     if(opt_stc == 2) then
      if (snowh > 0.05 .and. tg > tfrz) then
        tgv = tfrz
        tgb = tfrz
          if (veg .and. fveg > 0) then
             tg    = fveg * tgv       + (1.0 - fveg) * tgb
             ts    = fveg * tv        + (1.0 - fveg) * tgb
          else
             tg    = tgb
             ts    = tgb
          end if
      end if
     end if

! energy released or consumed by snow & frozen soil

 call phasechange (parameters,nsnow   ,nsoil   ,isnow   ,dt      ,fact    , & !in
                   dzsnso  ,hcpct   ,ist     ,iloc    ,jloc    , & !in
                   stc     ,snice   ,snliq   ,sneqv   ,snowh   , & !inout
#ifdef CCPP
                   smc     ,sh2o    ,errmsg  ,errflg  ,          & !inout
#else
                   smc     ,sh2o    ,                            & !inout
#endif
                   qmelt   ,imelt   ,ponding )                     !out
#ifdef CCPP
    if (errflg /= 0) return
#endif

  end subroutine energy

!== begin thermoprop ===============================================================================

!>\ingroup NoahMP_LSM
  subroutine thermoprop (parameters,nsoil   ,nsnow   ,isnow   ,ist     ,dzsnso  , & !in
                         dt      ,snowh   ,snice   ,snliq   , & !in
                         smc     ,sh2o    ,tg      ,stc     ,ur      , & !in
                         lat     ,z0m     ,zlvl    ,vegtyp  , & !in
                         df      ,hcpct   ,snicev  ,snliqv  ,epore   , & !out
                         fact    )                                       !out
! ------------------------------------------------------------------------------------------------- 
  implicit none
! --------------------------------------------------------------------------------------------------
! inputs
  type (noahmp_parameters), intent(in) :: parameters
  integer                        , intent(in)  :: nsoil   !number of soil layers
  integer                        , intent(in)  :: nsnow   !maximum no. of snow layers        
  integer                        , intent(in)  :: isnow   !actual no. of snow layers
  integer                        , intent(in)  :: ist     !surface type
  real                           , intent(in)  :: dt      !time step [s]
  real, dimension(-nsnow+1:    0), intent(in)  :: snice   !snow ice mass (kg/m2)
  real, dimension(-nsnow+1:    0), intent(in)  :: snliq   !snow liq mass (kg/m2)
  real, dimension(-nsnow+1:nsoil), intent(in)  :: dzsnso  !thickness of snow/soil layers [m]
  real, dimension(       1:nsoil), intent(in)  :: smc     !soil moisture (ice + liq.) [m3/m3]
  real, dimension(       1:nsoil), intent(in)  :: sh2o    !liquid soil moisture [m3/m3]
  real                           , intent(in)  :: snowh   !snow height [m]
  real,                            intent(in)  :: tg      !surface temperature (k)
  real, dimension(-nsnow+1:nsoil), intent(in)  :: stc     !snow/soil/lake temp. (k)
  real,                            intent(in)  :: ur      !wind speed at zlvl (m/s)
  real,                            intent(in)  :: lat     !latitude (radians)
  real,                            intent(in)  :: z0m     !roughness length (m)
  real,                            intent(in)  :: zlvl    !reference height (m)
  integer                        , intent(in)  :: vegtyp  !vegtyp type

! outputs
  real, dimension(-nsnow+1:nsoil), intent(out) :: df      !thermal conductivity [w/m/k]
  real, dimension(-nsnow+1:nsoil), intent(out) :: hcpct   !heat capacity [j/m3/k]
  real, dimension(-nsnow+1:    0), intent(out) :: snicev  !partial volume of ice [m3/m3]
  real, dimension(-nsnow+1:    0), intent(out) :: snliqv  !partial volume of liquid water [m3/m3]
  real, dimension(-nsnow+1:    0), intent(out) :: epore   !effective porosity [m3/m3]
  real, dimension(-nsnow+1:nsoil), intent(out) :: fact    !computing energy for phase change
! --------------------------------------------------------------------------------------------------
! locals

  integer :: iz
  real, dimension(-nsnow+1:    0)              :: cvsno   !volumetric specific heat (j/m3/k)
  real, dimension(-nsnow+1:    0)              :: tksno   !snow thermal conductivity (j/m3/k)
  real, dimension(       1:nsoil)              :: sice    !soil ice content
! --------------------------------------------------------------------------------------------------

! compute snow thermal conductivity and heat capacity

    call csnow (parameters,isnow   ,nsnow   ,nsoil   ,snice   ,snliq   ,dzsnso  , & !in
                tksno   ,cvsno   ,snicev  ,snliqv  ,epore   )   !out

    do iz = isnow+1, 0
      df   (iz) = tksno(iz)
      hcpct(iz) = cvsno(iz)
    end do

! compute soil thermal properties

    do  iz = 1, nsoil
       sice(iz)  = smc(iz) - sh2o(iz)
       hcpct(iz) = sh2o(iz)*cwat + (1.0-parameters%smcmax)*parameters%csoil &
                + (parameters%smcmax-smc(iz))*cpair + sice(iz)*cice
       call tdfcnd (parameters,df(iz), smc(iz), sh2o(iz))
    end do
       
    if ( parameters%urban_flag ) then
       do iz = 1,nsoil
         df(iz) = 3.24
       end do
    endif

! heat flux reduction effect from the overlying green canopy, adapted from 
! section 2.1.2 of peters-lidard et al. (1997, jgr, vol 102(d4)).
! not in use because of the separation of the canopy layer from the ground.
! but this may represent the effects of leaf litter (niu comments)
!       df1 = df1 * exp (sbeta * shdfac)

! compute lake thermal properties 
! (no consideration of turbulent mixing for this version)

    if(ist == 2) then
       do iz = 1, nsoil 
         if(stc(iz) > tfrz) then
            hcpct(iz) = cwat
            df(iz)    = tkwat  !+ keddy * cwat 
         else
            hcpct(iz) = cice
            df(iz)    = tkice 
         end if
       end do
    end if

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


  end subroutine thermoprop

!== begin csnow ====================================================================================

!>\ingroup NoahMP_LSM
  subroutine csnow (parameters,isnow   ,nsnow   ,nsoil   ,snice   ,snliq   ,dzsnso  , & !in
                    tksno   ,cvsno   ,snicev  ,snliqv  ,epore   )   !out
! --------------------------------------------------------------------------------------------------
! snow bulk density,volumetric capacity, and thermal conductivity
!---------------------------------------------------------------------------------------------------
  implicit none
!---------------------------------------------------------------------------------------------------
! inputs

  type (noahmp_parameters), intent(in) :: parameters
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

  end subroutine csnow

!== begin tdfcnd ===================================================================================

!>\ingroup NoahMP_LSM
  subroutine tdfcnd (parameters, df, smc, sh2o)
! --------------------------------------------------------------------------------------------------
! calculate thermal diffusivity and conductivity of the soil.
! peters-lidard approach (peters-lidard et al., 1998)
! --------------------------------------------------------------------------------------------------
! code history:
! june 2001 changes: frozen soil condition.
! --------------------------------------------------------------------------------------------------
    implicit none
  type (noahmp_parameters), intent(in) :: parameters
    real, intent(in)       :: smc    ! total soil water
    real, intent(in)       :: sh2o   ! liq. soil water
    real, intent(out)      :: df     ! thermal diffusivity

! local variables
    real  :: ake
    real  :: gammd
    real  :: thkdry
    real  :: thko     ! thermal conductivity for other soil components         
    real  :: thkqtz   ! thermal conductivity for quartz
    real  :: thksat   ! 
    real  :: thks     ! thermal conductivity for the solids
    real  :: thkw     ! water thermal conductivity
    real  :: satratio
    real  :: xu
    real  :: xunfroz
! --------------------------------------------------------------------------------------------------
! we now get quartz as an input argument (set in routine redprm):
!      data quartz /0.82, 0.10, 0.25, 0.60, 0.52,
!     &             0.35, 0.60, 0.40, 0.82/
! --------------------------------------------------------------------------------------------------
! if the soil has any moisture content compute a partial sum/product
! otherwise use a constant value which works well with most soils
! --------------------------------------------------------------------------------------------------
!  quartz ....quartz content (soil type dependent)
! --------------------------------------------------------------------------------------------------
! use as in peters-lidard, 1998 (modif. from johansen, 1975).

!                                  pablo grunmann, 08/17/98
! refs.:
!      farouki, o.t.,1986: thermal properties of soils. series on rock
!              and soil mechanics, vol. 11, trans tech, 136 pp.
!      johansen, o., 1975: thermal conductivity of soils. ph.d. thesis,
!              university of trondheim,
!      peters-lidard, c. d., et al., 1998: the effect of soil thermal
!              conductivity parameterization on surface energy fluxes
!              and temperatures. journal of the atmospheric sciences,
!              vol. 55, pp. 1209-1224.
! --------------------------------------------------------------------------------------------------
! needs parameters
! porosity(soil type):
!      poros = smcmax
! saturation ratio:
! parameters  w/(m.k)
    satratio = smc / parameters%smcmax
    thkw = 0.57
!      if (quartz .le. 0.2) thko = 3.0
    thko = 2.0
! solids' conductivity
! quartz' conductivity
    thkqtz = 7.7

! unfrozen fraction (from 1., i.e., 100%liquid, to 0. (100% frozen))
    thks = (thkqtz ** parameters%quartz)* (thko ** (1. - parameters%quartz))

! unfrozen volume for saturation (porosity*xunfroz)
    xunfroz = sh2o / smc
! saturated thermal conductivity
    xu = xunfroz * parameters%smcmax

! dry density in kg/m3
    thksat = thks ** (1. - parameters%smcmax)* tkice ** (parameters%smcmax - xu)* thkw **   &
         (xu)

! dry thermal conductivity in w.m-1.k-1
    gammd = (1. - parameters%smcmax)*2700.

    thkdry = (0.135* gammd+ 64.7)/ (2700. - 0.947* gammd)
! frozen
    if ( (sh2o + 0.0005) <  smc ) then
       ake = satratio
! unfrozen
! range of validity for the kersten number (ake)
    else

! kersten number (using "fine" formula, valid for soils containing at
! least 5% of particles with diameter less than 2.e-6 meters.)
! (for "coarse" formula, see peters-lidard et al., 1998).

       if ( satratio >  0.1 ) then

          ake = log10 (satratio) + 1.0

! use k = kdry
       else

          ake = 0.0
       end if
!  thermal conductivity

    end if

    df = ake * (thksat - thkdry) + thkdry


  end subroutine tdfcnd

!== begin radiation ================================================================================

!>\ingroup NoahMP_LSM
  subroutine radiation (parameters,vegtyp  ,ist     ,ice     ,nsoil   , & !in
                        sneqvo  ,sneqv   ,dt      ,cosz    ,snowh   , & !in
                        tg      ,tv      ,fsno    ,qsnow   ,fwet    , & !in
                        elai    ,esai    ,smc     ,solad   ,solai   , & !in
                        fveg    ,iloc    ,jloc    ,                   & !in
                        albold  ,tauss   ,                            & !inout
                        fsun    ,laisun  ,laisha  ,parsun  ,parsha  , & !out
                        sav     ,sag     ,fsr     ,fsa     ,fsrv    , &
                        fsrg    ,bgap    ,wgap)            !out
! --------------------------------------------------------------------------------------------------
  implicit none
! --------------------------------------------------------------------------------------------------
! input
  type (noahmp_parameters), intent(in) :: parameters
  integer, intent(in)                  :: iloc
  integer, intent(in)                  :: jloc
  integer, intent(in)                  :: vegtyp !vegetation type
  integer, intent(in)                  :: ist    !surface type
  integer, intent(in)                  :: ice    !ice (ice = 1)
  integer, intent(in)                  :: nsoil  !number of soil layers

  real, intent(in)                     :: dt     !time step [s]
  real, intent(in)                     :: qsnow  !snowfall (mm/s)
  real, intent(in)                     :: sneqvo !snow mass at last time step(mm)
  real, intent(in)                     :: sneqv  !snow mass (mm)
  real, intent(in)                     :: snowh  !snow height (mm)
  real, intent(in)                     :: cosz   !cosine solar zenith angle (0-1)
  real, intent(in)                     :: tg     !ground temperature (k)
  real, intent(in)                     :: tv     !vegetation temperature (k)
  real, intent(in)                     :: elai   !lai, one-sided, adjusted for burying by snow
  real, intent(in)                     :: esai   !sai, one-sided, adjusted for burying by snow
  real, intent(in)                     :: fwet   !fraction of canopy that is wet
  real, dimension(1:nsoil), intent(in) :: smc    !volumetric soil water [m3/m3]
  real, dimension(1:2)    , intent(in) :: solad  !incoming direct solar radiation (w/m2)
  real, dimension(1:2)    , intent(in) :: solai  !incoming diffuse solar radiation (w/m2)
  real, intent(in)                     :: fsno   !snow cover fraction (-)
  real, intent(in)                     :: fveg   !green vegetation fraction [0.0-1.0]

! inout
  real,                  intent(inout) :: albold !snow albedo at last time step (class type)
  real,                  intent(inout) :: tauss  !non-dimensional snow age.

! output
  real, intent(out)                    :: fsun   !sunlit fraction of canopy (-)
  real, intent(out)                    :: laisun !sunlit leaf area (-)
  real, intent(out)                    :: laisha !shaded leaf area (-)
  real, intent(out)                    :: parsun !average absorbed par for sunlit leaves (w/m2)
  real, intent(out)                    :: parsha !average absorbed par for shaded leaves (w/m2)
  real, intent(out)                    :: sav    !solar radiation absorbed by vegetation (w/m2)
  real, intent(out)                    :: sag    !solar radiation absorbed by ground (w/m2)
  real, intent(out)                    :: fsa    !total absorbed solar radiation (w/m2)
  real, intent(out)                    :: fsr    !total reflected solar radiation (w/m2)

!jref:start  
  real, intent(out)                    :: fsrv    !veg. reflected solar radiation (w/m2)
  real, intent(out)                    :: fsrg    !ground reflected solar radiation (w/m2)
  real, intent(out)                    :: bgap
  real, intent(out)                    :: wgap
!jref:end  

! local
  real                                 :: fage   !snow age function (0 - new snow)
  real, dimension(1:2)                 :: albgrd !ground albedo (direct)
  real, dimension(1:2)                 :: albgri !ground albedo (diffuse)
  real, dimension(1:2)                 :: albd   !surface albedo (direct)
  real, dimension(1:2)                 :: albi   !surface albedo (diffuse)
  real, dimension(1:2)                 :: fabd   !flux abs by veg (per unit direct flux)
  real, dimension(1:2)                 :: fabi   !flux abs by veg (per unit diffuse flux)
  real, dimension(1:2)                 :: ftdd   !down direct flux below veg (per unit dir flux)
  real, dimension(1:2)                 :: ftid   !down diffuse flux below veg (per unit dir flux)
  real, dimension(1:2)                 :: ftii   !down diffuse flux below veg (per unit dif flux)
!jref:start  
  real, dimension(1:2)                 :: frevi
  real, dimension(1:2)                 :: frevd
  real, dimension(1:2)                 :: fregi
  real, dimension(1:2)                 :: fregd
!jref:end

  real                                 :: fsha   !shaded fraction of canopy
  real                                 :: vai    !total lai + stem area index, one sided

  real,parameter :: mpe = 1.e-6
  logical veg  !true: vegetated for surface temperature calculation

! --------------------------------------------------------------------------------------------------

! surface abeldo

   call albedo (parameters,vegtyp ,ist    ,ice    ,nsoil  , & !in
                dt     ,cosz   ,fage   ,elai   ,esai   , & !in
                tg     ,tv     ,snowh  ,fsno   ,fwet   , & !in
                smc    ,sneqvo ,sneqv  ,qsnow  ,fveg   , & !in
                iloc   ,jloc   ,                         & !in
                albold ,tauss                          , & !inout
                albgrd ,albgri ,albd   ,albi   ,fabd   , & !out
                fabi   ,ftdd   ,ftid   ,ftii   ,fsun   , & !)   !out
                frevi  ,frevd   ,fregd ,fregi  ,bgap   , & !inout
                wgap)

! surface radiation

     fsha = 1.-fsun
     laisun = elai*fsun
     laisha = elai*fsha
     vai = elai+ esai
     if (vai .gt. 0.) then
        veg = .true.
     else
        veg = .false.
     end if

   call surrad (parameters,mpe    ,fsun   ,fsha   ,elai   ,vai    , & !in
                laisun ,laisha ,solad  ,solai  ,fabd   , & !in
                fabi   ,ftdd   ,ftid   ,ftii   ,albgrd , & !in
                albgri ,albd   ,albi   ,iloc   ,jloc   , & !in
                parsun ,parsha ,sav    ,sag    ,fsa    , & !out
                fsr    ,                                 & !out
                frevi  ,frevd  ,fregd  ,fregi  ,fsrv   , & !inout
                fsrg)

  end subroutine radiation

!== begin albedo ===================================================================================

!>\ingroup NoahMP_LSM
  subroutine albedo (parameters,vegtyp ,ist    ,ice    ,nsoil  , & !in
                     dt     ,cosz   ,fage   ,elai   ,esai   , & !in
                     tg     ,tv     ,snowh  ,fsno   ,fwet   , & !in
                     smc    ,sneqvo ,sneqv  ,qsnow  ,fveg   , & !in
                     iloc   ,jloc   ,                         & !in
                     albold ,tauss                          , & !inout
                     albgrd ,albgri ,albd   ,albi   ,fabd   , & !out
                     fabi   ,ftdd   ,ftid   ,ftii   ,fsun   , & !out
                     frevi  ,frevd  ,fregd  ,fregi  ,bgap   , & !out
                     wgap)

! --------------------------------------------------------------------------------------------------
! surface albedos. also fluxes (per unit incoming direct and diffuse
! radiation) reflected, transmitted, and absorbed by vegetation.
! also sunlit fraction of the canopy.
! --------------------------------------------------------------------------------------------------
  implicit none
! --------------------------------------------------------------------------------------------------
! input
  type (noahmp_parameters), intent(in) :: parameters
  integer,                  intent(in)  :: iloc
  integer,                  intent(in)  :: jloc
  integer,                  intent(in)  :: nsoil  !number of soil layers
  integer,                  intent(in)  :: vegtyp !vegetation type
  integer,                  intent(in)  :: ist    !surface type
  integer,                  intent(in)  :: ice    !ice (ice = 1)

  real,                     intent(in)  :: dt     !time step [sec]
  real,                     intent(in)  :: qsnow  !snowfall
  real,                     intent(in)  :: cosz   !cosine solar zenith angle for next time step
  real,                     intent(in)  :: snowh  !snow height (mm)
  real,                     intent(in)  :: tg     !ground temperature (k)
  real,                     intent(in)  :: tv     !vegetation temperature (k)
  real,                     intent(in)  :: elai   !lai, one-sided, adjusted for burying by snow
  real,                     intent(in)  :: esai   !sai, one-sided, adjusted for burying by snow
  real,                     intent(in)  :: fsno   !fraction of grid covered by snow
  real,                     intent(in)  :: fwet   !fraction of canopy that is wet
  real,                     intent(in)  :: sneqvo !snow mass at last time step(mm)
  real,                     intent(in)  :: sneqv  !snow mass (mm)
  real,                     intent(in)  :: fveg   !green vegetation fraction [0.0-1.0]
  real, dimension(1:nsoil), intent(in)  :: smc    !volumetric soil water (m3/m3)

! inout
  real,                  intent(inout)  :: albold !snow albedo at last time step (class type)
  real,                  intent(inout)  :: tauss  !non-dimensional snow age

! output
  real, dimension(1:    2), intent(out) :: albgrd !ground albedo (direct)
  real, dimension(1:    2), intent(out) :: albgri !ground albedo (diffuse)
  real, dimension(1:    2), intent(out) :: albd   !surface albedo (direct)
  real, dimension(1:    2), intent(out) :: albi   !surface albedo (diffuse)
  real, dimension(1:    2), intent(out) :: fabd   !flux abs by veg (per unit direct flux)
  real, dimension(1:    2), intent(out) :: fabi   !flux abs by veg (per unit diffuse flux)
  real, dimension(1:    2), intent(out) :: ftdd   !down direct flux below veg (per unit dir flux)
  real, dimension(1:    2), intent(out) :: ftid   !down diffuse flux below veg (per unit dir flux)
  real, dimension(1:    2), intent(out) :: ftii   !down diffuse flux below veg (per unit dif flux)
  real,                     intent(out) :: fsun   !sunlit fraction of canopy (-)
!jref:start
  real, dimension(1:    2), intent(out) :: frevd
  real, dimension(1:    2), intent(out) :: frevi
  real, dimension(1:    2), intent(out) :: fregd
  real, dimension(1:    2), intent(out) :: fregi
  real, intent(out) :: bgap
  real, intent(out) :: wgap
!jref:end

! ------------------------------------------------------------------------
! ------------------------ local variables -------------------------------
! local
  real                 :: fage     !snow age function
  real                 :: alb
  integer              :: ib       !indices
  integer              :: nband    !number of solar radiation wave bands
  integer              :: ic       !direct beam: ic=0; diffuse: ic=1

  real                 :: wl       !fraction of lai+sai that is lai
  real                 :: ws       !fraction of lai+sai that is sai
  real                 :: mpe      !prevents overflow for division by zero

  real, dimension(1:2) :: rho      !leaf/stem reflectance weighted by fraction lai and sai
  real, dimension(1:2) :: tau      !leaf/stem transmittance weighted by fraction lai and sai
  real, dimension(1:2) :: ftdi     !down direct flux below veg per unit dif flux = 0
  real, dimension(1:2) :: albsnd   !snow albedo (direct)
  real, dimension(1:2) :: albsni   !snow albedo (diffuse)

  real                 :: vai      !elai+esai
  real                 :: gdir     !average projected leaf/stem area in solar direction
  real                 :: ext      !optical depth direct beam per unit leaf + stem area

! --------------------------------------------------------------------------------------------------

  nband = 2
  mpe = 1.e-06
  bgap = 0.
  wgap = 0.

! initialize output because solar radiation only done if cosz > 0

  do ib = 1, nband
    albd(ib) = 0.
    albi(ib) = 0.
    albgrd(ib) = 0.
    albgri(ib) = 0.
    fabd(ib) = 0.
    fabi(ib) = 0.
    ftdd(ib) = 0.
    ftid(ib) = 0.
    ftii(ib) = 0.
    if (ib.eq.1) fsun = 0.
  end do

  if(cosz <= 0) goto 100

! weight reflectance/transmittance by lai and sai

  do ib = 1, nband
    vai = elai + esai
    wl  = elai / max(vai,mpe)
    ws  = esai / max(vai,mpe)
    rho(ib) = max(parameters%rhol(ib)*wl+parameters%rhos(ib)*ws, mpe)
    tau(ib) = max(parameters%taul(ib)*wl+parameters%taus(ib)*ws, mpe)
  end do

! snow age

   call snow_age (parameters,dt,tg,sneqvo,sneqv,tauss,fage)

! snow albedos: only if cosz > 0 and fsno > 0

  if(opt_alb == 1) &
     call snowalb_bats (parameters,nband, fsno,cosz,fage,albsnd,albsni)
  if(opt_alb == 2) then
     call snowalb_class (parameters,nband,qsnow,dt,alb,albold,albsnd,albsni,iloc,jloc)
     albold = alb
  end if

! ground surface albedo

  call groundalb (parameters,nsoil   ,nband   ,ice     ,ist     , & !in
                  fsno    ,smc     ,albsnd  ,albsni  ,cosz    , & !in
                  tg      ,iloc    ,jloc    ,                   & !in
                  albgrd  ,albgri  )                              !out

! loop over nband wavebands to calculate surface albedos and solar
! fluxes for unit incoming direct (ic=0) and diffuse flux (ic=1)

  do ib = 1, nband
      ic = 0      ! direct
      call twostream (parameters,ib     ,ic      ,vegtyp  ,cosz    ,vai    , & !in
                      fwet   ,tv      ,albgrd  ,albgri  ,rho    , & !in
                      tau    ,fveg    ,ist     ,iloc    ,jloc   , & !in
                      fabd   ,albd    ,ftdd    ,ftid    ,gdir   , &!)   !out
                      frevd  ,fregd   ,bgap    ,wgap)

      ic = 1      ! diffuse
      call twostream (parameters,ib     ,ic      ,vegtyp  ,cosz    ,vai    , & !in
                      fwet   ,tv      ,albgrd  ,albgri  ,rho    , & !in
                      tau    ,fveg    ,ist     ,iloc    ,jloc   , & !in
                      fabi   ,albi    ,ftdi    ,ftii    ,gdir   , & !)   !out
                      frevi  ,fregi   ,bgap    ,wgap)

  end do

! sunlit fraction of canopy. set fsun = 0 if fsun < 0.01.

  ext = gdir/cosz * sqrt(1.-rho(1)-tau(1))
  fsun = (1.-exp(-ext*vai)) / max(ext*vai,mpe)
  ext = fsun

  if (ext .lt. 0.01) then
     wl = 0.
  else
     wl = ext 
  end if
  fsun = wl

100 continue

  end subroutine albedo

!== begin surrad ===================================================================================

!>\ingroup NoahMP_LSM
  subroutine surrad (parameters,mpe     ,fsun    ,fsha    ,elai    ,vai     , & !in
                     laisun  ,laisha  ,solad   ,solai   ,fabd    , & !in
                     fabi    ,ftdd    ,ftid    ,ftii    ,albgrd  , & !in
                     albgri  ,albd    ,albi    ,iloc    ,jloc    , & !in
                     parsun  ,parsha  ,sav     ,sag     ,fsa     , & !out
                     fsr     , & !)                                       !out
                     frevi   ,frevd   ,fregd   ,fregi   ,fsrv    , &
                     fsrg) !inout

! --------------------------------------------------------------------------------------------------
  implicit none
! --------------------------------------------------------------------------------------------------
! input

  type (noahmp_parameters), intent(in) :: parameters
  integer, intent(in)              :: iloc
  integer, intent(in)              :: jloc
  real, intent(in)                 :: mpe     !prevents underflow errors if division by zero

  real, intent(in)                 :: fsun    !sunlit fraction of canopy
  real, intent(in)                 :: fsha    !shaded fraction of canopy
  real, intent(in)                 :: elai    !leaf area, one-sided
  real, intent(in)                 :: vai     !leaf + stem area, one-sided
  real, intent(in)                 :: laisun  !sunlit leaf area index, one-sided
  real, intent(in)                 :: laisha  !shaded leaf area index, one-sided

  real, dimension(1:2), intent(in) :: solad   !incoming direct solar radiation (w/m2)
  real, dimension(1:2), intent(in) :: solai   !incoming diffuse solar radiation (w/m2)
  real, dimension(1:2), intent(in) :: fabd    !flux abs by veg (per unit incoming direct flux)
  real, dimension(1:2), intent(in) :: fabi    !flux abs by veg (per unit incoming diffuse flux)
  real, dimension(1:2), intent(in) :: ftdd    !down dir flux below veg (per incoming dir flux)
  real, dimension(1:2), intent(in) :: ftid    !down dif flux below veg (per incoming dir flux)
  real, dimension(1:2), intent(in) :: ftii    !down dif flux below veg (per incoming dif flux)
  real, dimension(1:2), intent(in) :: albgrd  !ground albedo (direct)
  real, dimension(1:2), intent(in) :: albgri  !ground albedo (diffuse)
  real, dimension(1:2), intent(in) :: albd    !overall surface albedo (direct)
  real, dimension(1:2), intent(in) :: albi    !overall surface albedo (diffuse)

  real, dimension(1:2), intent(in) :: frevd    !overall surface albedo veg (direct)
  real, dimension(1:2), intent(in) :: frevi    !overall surface albedo veg (diffuse)
  real, dimension(1:2), intent(in) :: fregd    !overall surface albedo grd (direct)
  real, dimension(1:2), intent(in) :: fregi    !overall surface albedo grd (diffuse)

! output

  real, intent(out)                :: parsun  !average absorbed par for sunlit leaves (w/m2)
  real, intent(out)                :: parsha  !average absorbed par for shaded leaves (w/m2)
  real, intent(out)                :: sav     !solar radiation absorbed by vegetation (w/m2)
  real, intent(out)                :: sag     !solar radiation absorbed by ground (w/m2)
  real, intent(out)                :: fsa     !total absorbed solar radiation (w/m2)
  real, intent(out)                :: fsr     !total reflected solar radiation (w/m2)
  real, intent(out)                :: fsrv    !reflected solar radiation by vegetation
  real, intent(out)                :: fsrg    !reflected solar radiation by ground

! ------------------------ local variables ----------------------------------------------------
  integer                          :: ib      !waveband number (1=vis, 2=nir)
  integer                          :: nband   !number of solar radiation waveband classes

  real                             :: abs     !absorbed solar radiation (w/m2)
  real                             :: rnir    !reflected solar radiation [nir] (w/m2)
  real                             :: rvis    !reflected solar radiation [vis] (w/m2)
  real                             :: laifra  !leaf area fraction of canopy
  real                             :: trd     !transmitted solar radiation: direct (w/m2)
  real                             :: tri     !transmitted solar radiation: diffuse (w/m2)
  real, dimension(1:2)             :: cad     !direct beam absorbed by canopy (w/m2)
  real, dimension(1:2)             :: cai     !diffuse radiation absorbed by canopy (w/m2)
! ---------------------------------------------------------------------------------------------
   nband = 2

! zero summed solar fluxes

    sag = 0.
    sav = 0.
    fsa = 0.

! loop over nband wavebands

  do ib = 1, nband

! absorbed by canopy

    cad(ib) = solad(ib)*fabd(ib)    
    cai(ib) = solai(ib)*fabi(ib)
    sav     = sav + cad(ib) + cai(ib)
    fsa     = fsa + cad(ib) + cai(ib)
 
! transmitted solar fluxes incident on ground

    trd = solad(ib)*ftdd(ib)
    tri = solad(ib)*ftid(ib) + solai(ib)*ftii(ib)

! solar radiation absorbed by ground surface

    abs = trd*(1.-albgrd(ib)) + tri*(1.-albgri(ib))
    sag = sag + abs
    fsa = fsa + abs
  end do

! partition visible canopy absorption to sunlit and shaded fractions
! to get average absorbed par for sunlit and shaded leaves

     laifra = elai / max(vai,mpe)
     if (fsun .gt. 0.) then
        parsun = (cad(1)+fsun*cai(1)) * laifra / max(laisun,mpe)
        parsha = (fsha*cai(1))*laifra / max(laisha,mpe)
     else
        parsun = 0.
        parsha = (cad(1)+cai(1))*laifra /max(laisha,mpe)
     endif

! reflected solar radiation

     rvis = albd(1)*solad(1) + albi(1)*solai(1)
     rnir = albd(2)*solad(2) + albi(2)*solai(2)
     fsr  = rvis + rnir

! reflected solar radiation of veg. and ground (combined ground)
     fsrv = frevd(1)*solad(1)+frevi(1)*solai(1)+frevd(2)*solad(2)+frevi(2)*solai(2)
     fsrg = fregd(1)*solad(1)+fregi(1)*solai(1)+fregd(2)*solad(2)+fregi(2)*solai(2)


  end subroutine surrad

!== begin snow_age =================================================================================

!>\ingroup NoahMP_LSM
  subroutine snow_age (parameters,dt,tg,sneqvo,sneqv,tauss,fage)
! ----------------------------------------------------------------------
  implicit none
! ------------------------ code history ------------------------------------------------------------
! from bats
! ------------------------ input/output variables --------------------------------------------------
!input
  type (noahmp_parameters), intent(in) :: parameters
   real, intent(in) :: dt        !main time step (s)
   real, intent(in) :: tg        !ground temperature (k)
   real, intent(in) :: sneqvo    !snow mass at last time step(mm)
   real, intent(in) :: sneqv     !snow water per unit ground area (mm)

!output
   real, intent(out) :: fage     !snow age

!input/output
   real, intent(inout) :: tauss      !non-dimensional snow age
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
          dela0 = 1.e-6*dt
          arg   = 5.e3*(1./tfrz-1./tg)
          age1  = exp(arg)
          age2  = exp(amin1(0.,10.*arg))
          age3  = 0.3
          tage  = age1+age2+age3
          dela  = dela0*tage
          dels  = amax1(0.0,sneqv-sneqvo) / parameters%swemx
          sge   = (tauss+dela)*(1.0-dels)
          tauss = amax1(0.,sge)
   endif

   fage= tauss/(tauss+1.)

  end subroutine snow_age

!== begin snowalb_bats =============================================================================

!>\ingroup NoahMP_LSM
  subroutine snowalb_bats (parameters,nband,fsno,cosz,fage,albsnd,albsni)
! --------------------------------------------------------------------------------------------------
  implicit none
! --------------------------------------------------------------------------------------------------
! input

  type (noahmp_parameters), intent(in) :: parameters
  integer,intent(in) :: nband  !number of waveband classes

  real,intent(in) :: cosz    !cosine solar zenith angle
  real,intent(in) :: fsno    !snow cover fraction (-)
  real,intent(in) :: fage    !snow age correction

! output

  real, dimension(1:2),intent(out) :: albsnd !snow albedo for direct(1=vis, 2=nir)
  real, dimension(1:2),intent(out) :: albsni !snow albedo for diffuse
! ---------------------------------------------------------------------------------------------

! ------------------------ local variables ----------------------------------------------------
  integer :: ib          !waveband class

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

  end subroutine snowalb_bats

!== begin snowalb_class ============================================================================

!>\ingroup NoahMP_LSM
  subroutine snowalb_class (parameters,nband,qsnow,dt,alb,albold,albsnd,albsni,iloc,jloc)
! ----------------------------------------------------------------------
  implicit none
! --------------------------------------------------------------------------------------------------
! input

  type (noahmp_parameters), intent(in) :: parameters
  integer,intent(in) :: iloc !grid index
  integer,intent(in) :: jloc !grid index
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

! ------------------------ local variables ----------------------------------------------------
  integer :: ib          !waveband class

! ---------------------------------------------------------------------------------------------
! zero albedos for all points

        albsnd(1: nband) = 0.
        albsni(1: nband) = 0.

! when cosz > 0

         alb = 0.55 + (albold-0.55) * exp(-0.01*dt/3600.)

! 1 mm fresh snow(swe) -- 10mm snow depth, assumed the fresh snow density 100kg/m3
! here assume 1cm snow depth will fully cover the old snow

         if (qsnow > 0.) then
           alb = alb + min(qsnow,parameters%swemx/dt) * (0.84-alb)/(parameters%swemx/dt)
         endif

         albsni(1)= alb         ! vis diffuse
         albsni(2)= alb         ! nir diffuse
         albsnd(1)= alb         ! vis direct
         albsnd(2)= alb         ! nir direct

  end subroutine snowalb_class

!== begin groundalb ================================================================================

!>\ingroup NoahMP_LSM
  subroutine groundalb (parameters,nsoil   ,nband   ,ice     ,ist     , & !in
                        fsno    ,smc     ,albsnd  ,albsni  ,cosz    , & !in
                        tg      ,iloc    ,jloc    ,                   & !in
                        albgrd  ,albgri  )                              !out
! --------------------------------------------------------------------------------------------------
  implicit none
! --------------------------------------------------------------------------------------------------
!input

  type (noahmp_parameters), intent(in) :: parameters
  integer,                  intent(in)  :: iloc   !grid index
  integer,                  intent(in)  :: jloc   !grid index
  integer,                  intent(in)  :: nsoil  !number of soil layers
  integer,                  intent(in)  :: nband  !number of solar radiation waveband classes
  integer,                  intent(in)  :: ice    !value of ist for land ice
  integer,                  intent(in)  :: ist    !surface type
  real,                     intent(in)  :: fsno   !fraction of surface covered with snow (-)
  real,                     intent(in)  :: tg     !ground temperature (k)
  real,                     intent(in)  :: cosz   !cosine solar zenith angle (0-1)
  real, dimension(1:nsoil), intent(in)  :: smc    !volumetric soil water content (m3/m3)
  real, dimension(1:    2), intent(in)  :: albsnd !direct beam snow albedo (vis, nir)
  real, dimension(1:    2), intent(in)  :: albsni !diffuse snow albedo (vis, nir)

!output

  real, dimension(1:    2), intent(out) :: albgrd !ground albedo (direct beam: vis, nir)
  real, dimension(1:    2), intent(out) :: albgri !ground albedo (diffuse: vis, nir)

!local 

  integer                               :: ib     !waveband number (1=vis, 2=nir)
  real                                  :: inc    !soil water correction factor for soil albedo
  real                                  :: albsod !soil albedo (direct)
  real                                  :: albsoi !soil albedo (diffuse)
! --------------------------------------------------------------------------------------------------

  do ib = 1, nband
        inc = max(0.11-0.40*smc(1), 0.)
        if (ist .eq. 1)  then                     !soil
           albsod = min(parameters%albsat(ib)+inc,parameters%albdry(ib))
           albsoi = albsod
        else if (tg .gt. tfrz) then               !unfrozen lake, wetland
           albsod = 0.06/(max(0.01,cosz)**1.7 + 0.15)
           albsoi = 0.06
        else                                      !frozen lake, wetland
           albsod = parameters%alblak(ib)
           albsoi = albsod
        end if

! increase desert and semi-desert albedos

!        if (ist .eq. 1 .and. isc .eq. 9) then
!           albsod = albsod + 0.10
!           albsoi = albsoi + 0.10
!        end if

        albgrd(ib) = albsod*(1.-fsno) + albsnd(ib)*fsno
        albgri(ib) = albsoi*(1.-fsno) + albsni(ib)*fsno
  end do

  end subroutine groundalb

!== begin twostream ================================================================================

!>\ingroup NoahMP_LSM
  subroutine twostream (parameters,ib     ,ic      ,vegtyp  ,cosz    ,vai    , & !in
                        fwet   ,t       ,albgrd  ,albgri  ,rho    , & !in
                        tau    ,fveg    ,ist     ,iloc    ,jloc   , & !in
                        fab    ,fre     ,ftd     ,fti     ,gdir   , & !)   !out
                        frev   ,freg    ,bgap    ,wgap)

! --------------------------------------------------------------------------------------------------
! use two-stream approximation of dickinson (1983) adv geophysics
! 25:305-353 and sellers (1985) int j remote sensing 6:1335-1372
! to calculate fluxes absorbed by vegetation, reflected by vegetation,
! and transmitted through vegetation for unit incoming direct or diffuse
! flux given an underlying surface with known albedo.
! --------------------------------------------------------------------------------------------------
  implicit none
! --------------------------------------------------------------------------------------------------
! input

  type (noahmp_parameters), intent(in) :: parameters
   integer,              intent(in)  :: iloc    !grid index
   integer,              intent(in)  :: jloc    !grid index
   integer,              intent(in)  :: ist     !surface type
   integer,              intent(in)  :: ib      !waveband number
   integer,              intent(in)  :: ic      !0=unit incoming direct; 1=unit incoming diffuse
   integer,              intent(in)  :: vegtyp  !vegetation type

   real,                 intent(in)  :: cosz    !cosine of direct zenith angle (0-1)
   real,                 intent(in)  :: vai     !one-sided leaf+stem area index (m2/m2)
   real,                 intent(in)  :: fwet    !fraction of lai, sai that is wetted (-)
   real,                 intent(in)  :: t       !surface temperature (k)

   real, dimension(1:2), intent(in)  :: albgrd  !direct  albedo of underlying surface (-)
   real, dimension(1:2), intent(in)  :: albgri  !diffuse albedo of underlying surface (-)
   real, dimension(1:2), intent(in)  :: rho     !leaf+stem reflectance
   real, dimension(1:2), intent(in)  :: tau     !leaf+stem transmittance
   real,                 intent(in)  :: fveg    !green vegetation fraction [0.0-1.0]

! output

   real, dimension(1:2), intent(out) :: fab     !flux abs by veg layer (per unit incoming flux)
   real, dimension(1:2), intent(out) :: fre     !flux refl above veg layer (per unit incoming flux)
   real, dimension(1:2), intent(out) :: ftd     !down dir flux below veg layer (per unit in flux)
   real, dimension(1:2), intent(out) :: fti     !down dif flux below veg layer (per unit in flux)
   real,                 intent(out) :: gdir    !projected leaf+stem area in solar direction
   real, dimension(1:2), intent(out) :: frev    !flux reflected by veg layer   (per unit incoming flux) 
   real, dimension(1:2), intent(out) :: freg    !flux reflected by ground (per unit incoming flux)

! local
   real                              :: omega   !fraction of intercepted radiation that is scattered
   real                              :: omegal  !omega for leaves
   real                              :: betai   !upscatter parameter for diffuse radiation
   real                              :: betail  !betai for leaves
   real                              :: betad   !upscatter parameter for direct beam radiation
   real                              :: betadl  !betad for leaves
   real                              :: ext     !optical depth of direct beam per unit leaf area
   real                              :: avmu    !average diffuse optical depth

   real                              :: coszi   !0.001 <= cosz <= 1.000
   real                              :: asu     !single scattering albedo
   real                              :: chil    ! -0.4 <= xl <= 0.6

   real                              :: tmp0,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9
   real                              :: p1,p2,p3,p4,s1,s2,u1,u2,u3
   real                              :: b,c,d,d1,d2,f,h,h1,h2,h3,h4,h5,h6,h7,h8,h9,h10
   real                              :: phi1,phi2,sigma
   real                              :: ftds,ftis,fres
   real                              :: denfveg
   real                              :: vai_spread
!jref:start
   real                              :: freveg,frebar,ftdveg,ftiveg,ftdbar,ftibar
   real                              :: thetaz
!jref:end   

!  variables for the modified two-stream scheme
!  niu and yang (2004), jgr

   real, parameter :: pai = 3.14159265 
   real :: hd       !crown depth (m)
   real :: bb       !vertical crown radius (m)
   real :: thetap   !angle conversion from sza 
   real :: fa       !foliage volume density (m-1)
   real :: newvai   !effective lsai (-)

   real,intent(inout) :: bgap     !between canopy gap fraction for beam (-)
   real,intent(inout) :: wgap     !within canopy gap fraction for beam (-)

   real :: kopen    !gap fraction for diffue light (-)
   real :: gap      !total gap fraction for beam ( <=1-shafac )

! -----------------------------------------------------------------
! compute within and between gaps
     vai_spread = vai
     if(vai == 0.0) then
         gap     = 1.0
         kopen   = 1.0
     else
         if(opt_rad == 1) then
	   denfveg = -log(max(1.0-fveg,0.01))/(pai*parameters%rc**2)
           hd      = parameters%hvt - parameters%hvb
           bb      = 0.5 * hd           
           thetap  = atan(bb/parameters%rc * tan(acos(max(0.01,cosz))) )
           ! bgap    = exp(-parameters%den * pai * parameters%rc**2/cos(thetap) )
           bgap    = exp(-denfveg * pai * parameters%rc**2/cos(thetap) )
           fa      = vai/(1.33 * pai * parameters%rc**3.0 *(bb/parameters%rc)*denfveg)
           newvai  = hd*fa
           wgap    = (1.0-bgap) * exp(-0.5*newvai/cosz)
           gap     = min(1.0-fveg, bgap+wgap)

           kopen   = 0.05
         end if

         if(opt_rad == 2) then
           gap     = 0.0
           kopen   = 0.0
         end if

         if(opt_rad == 3) then
           gap     = 1.0-fveg
           kopen   = 1.0-fveg
         end if
     end if

! calculate two-stream parameters omega, betad, betai, avmu, gdir, ext.
! omega, betad, betai are adjusted for snow. values for omega*betad
! and omega*betai are calculated and then divided by the new omega
! because the product omega*betai, omega*betad is used in solution.
! also, the transmittances and reflectances (tau, rho) are linear
! weights of leaf and stem values.

     coszi  = max(0.001, cosz)
     chil   = min( max(parameters%xl, -0.4), 0.6)
     if (abs(chil) .le. 0.01) chil = 0.01
     phi1   = 0.5 - 0.633*chil - 0.330*chil*chil
     phi2   = 0.877 * (1.-2.*phi1)
     gdir   = phi1 + phi2*coszi
     ext    = gdir/coszi
     avmu   = ( 1. - phi1/phi2 * log((phi1+phi2)/phi1) ) / phi2
     omegal = rho(ib) + tau(ib)
     tmp0   = gdir + phi2*coszi
     tmp1   = phi1*coszi
     asu    = 0.5*omegal*gdir/tmp0 * ( 1.-tmp1/tmp0*log((tmp1+tmp0)/tmp1) )
     betadl = (1.+avmu*ext)/(omegal*avmu*ext)*asu
     betail = 0.5 * ( rho(ib)+tau(ib) + (rho(ib)-tau(ib))   &
            * ((1.+chil)/2.)**2 ) / omegal

! adjust omega, betad, and betai for intercepted snow

     if (t .gt. tfrz) then                                !no snow
        tmp0 = omegal
        tmp1 = betadl
        tmp2 = betail
     else
        tmp0 =   (1.-fwet)*omegal        + fwet*parameters%omegas(ib)
        tmp1 = ( (1.-fwet)*omegal*betadl + fwet*parameters%omegas(ib)*parameters%betads ) / tmp0
        tmp2 = ( (1.-fwet)*omegal*betail + fwet*parameters%omegas(ib)*parameters%betais ) / tmp0
     end if

     omega = tmp0
     betad = tmp1
     betai = tmp2

! absorbed, reflected, transmitted fluxes per unit incoming radiation

     b = 1. - omega + omega*betai
     c = omega*betai
     tmp0 = avmu*ext
     d = tmp0 * omega*betad
     f = tmp0 * omega*(1.-betad)
     tmp1 = b*b - c*c
     h = sqrt(tmp1) / avmu
     sigma = tmp0*tmp0 - tmp1
     if ( abs (sigma) < 1.e-6 ) sigma = sign(1.e-6,sigma)
     p1 = b + avmu*h
     p2 = b - avmu*h
     p3 = b + tmp0
     p4 = b - tmp0
     s1 = exp(-h*vai)
     s2 = exp(-ext*vai)
     if (ic .eq. 0) then
        u1 = b - c/albgrd(ib)
        u2 = b - c*albgrd(ib)
        u3 = f + c*albgrd(ib)
     else
        u1 = b - c/albgri(ib)
        u2 = b - c*albgri(ib)
        u3 = f + c*albgri(ib)
     end if
     tmp2 = u1 - avmu*h
     tmp3 = u1 + avmu*h
     d1 = p1*tmp2/s1 - p2*tmp3*s1
     tmp4 = u2 + avmu*h
     tmp5 = u2 - avmu*h
     d2 = tmp4/s1 - tmp5*s1
     h1 = -d*p4 - c*f
     tmp6 = d - h1*p3/sigma
     tmp7 = ( d - c - h1/sigma*(u1+tmp0) ) * s2
     h2 = ( tmp6*tmp2/s1 - p2*tmp7 ) / d1
     h3 = - ( tmp6*tmp3*s1 - p1*tmp7 ) / d1
     h4 = -f*p3 - c*d
     tmp8 = h4/sigma
     tmp9 = ( u3 - tmp8*(u2-tmp0) ) * s2
     h5 = - ( tmp8*tmp4/s1 + tmp9 ) / d2
     h6 = ( tmp8*tmp5*s1 + tmp9 ) / d2
     h7 = (c*tmp2) / (d1*s1)
     h8 = (-c*tmp3*s1) / d1
     h9 = tmp4 / (d2*s1)
     h10 = (-tmp5*s1) / d2

! downward direct and diffuse fluxes below vegetation
! niu and yang (2004), jgr.

     if (ic .eq. 0) then
        ftds = s2                           *(1.0-gap) + gap
        ftis = (h4*s2/sigma + h5*s1 + h6/s1)*(1.0-gap)
     else
        ftds = 0.
        ftis = (h9*s1 + h10/s1)*(1.0-kopen) + kopen
     end if
     ftd(ib) = ftds
     fti(ib) = ftis

! flux reflected by the surface (veg. and ground)

     if (ic .eq. 0) then
        fres   = (h1/sigma + h2 + h3)*(1.0-gap  ) + albgrd(ib)*gap        
        freveg = (h1/sigma + h2 + h3)*(1.0-gap  ) 
        frebar = albgrd(ib)*gap                   !jref - separate veg. and ground reflection
     else
        fres   = (h7 + h8) *(1.0-kopen) + albgri(ib)*kopen        
        freveg = (h7 + h8) *(1.0-kopen) + albgri(ib)*kopen
        frebar = 0                                !jref - separate veg. and ground reflection
     end if
     fre(ib) = fres

     frev(ib) = freveg 
     freg(ib) = frebar 

! flux absorbed by vegetation

     fab(ib) = 1. - fre(ib) - (1.-albgrd(ib))*ftd(ib) &
                            - (1.-albgri(ib))*fti(ib)

!if(iloc == 1.and.jloc ==  2) then
!  write(*,'(a7,2i2,5(a6,f8.4),2(a9,f8.4))') "ib,ic: ",ib,ic," gap: ",gap," ftd: ",ftd(ib)," fti: ",fti(ib)," fre: ", &
!         fre(ib)," fab: ",fab(ib)," albgrd: ",albgrd(ib)," albgri: ",albgri(ib)
!end if

  end subroutine twostream

!== begin vege_flux ================================================================================

!>\ingroup NoahMP_LSM
  subroutine vege_flux(parameters,nsnow   ,nsoil   ,isnow   ,vegtyp  ,veg     , & !in
                       dt      ,sav     ,sag     ,lwdn    ,ur      , & !in
                       uu      ,vv      ,sfctmp  ,thair   ,qair    , & !in
                       eair    ,rhoair  ,snowh   ,vai     ,gammav   ,gammag,  & !in
                       fwet    ,laisun  ,laisha  ,cwp     ,dzsnso  , & !in
                       zlvl    ,cpfac            , & !in
                       zpd     ,z0m     ,fveg    , & !in
                       z0mg    ,emv     ,emg     ,canliq  ,fsno,          & !in
                       canice  ,stc     ,df      ,rssun   ,rssha   , & !in
                       rsurf   ,latheav ,latheag  ,parsun  ,parsha  ,igs     , & !in
                       foln    ,co2air  ,o2air   ,btran   ,sfcprs  , & !in
                       rhsur   ,iloc    ,jloc    ,q2      ,pahv    ,pahg     , & !in
                       eah     ,tah     ,tv      ,tg      ,cm      , & !inout
#ifdef CCPP
                       ch      ,dx      ,dz8w    ,errmsg  ,errflg  , & !inout
#else
                       ch      ,dx      ,dz8w    ,                   & !inout
#endif
                       tauxv   ,tauyv   ,irg     ,irc     ,shg     , & !out
                       shc     ,evg     ,evc     ,tr      ,gh      , & !out
                       t2mv    ,psnsun  ,psnsha  ,                   & !out
                       qc      ,qsfc    ,psfc    ,                   & !in
                       q2v     ,cah2    ,chleaf  ,chuc    )            !inout 

! --------------------------------------------------------------------------------------------------
! use newton-raphson iteration to solve for vegetation (tv) and
! ground (tg) temperatures that balance the surface energy budgets

! vegetated:
! -sav + irc[tv] + shc[tv] + evc[tv] + tr[tv] = 0
! -sag + irg[tg] + shg[tg] + evg[tg] + gh[tg] = 0
! --------------------------------------------------------------------------------------------------
  implicit none
! --------------------------------------------------------------------------------------------------
! input
  type (noahmp_parameters), intent(in) :: parameters
  integer,                         intent(in) :: iloc   !grid index
  integer,                         intent(in) :: jloc   !grid index
  logical,                         intent(in) :: veg    !true if vegetated surface
  integer,                         intent(in) :: nsnow  !maximum no. of snow layers        
  integer,                         intent(in) :: nsoil  !number of soil layers
  integer,                         intent(in) :: isnow  !actual no. of snow layers
  integer,                         intent(in) :: vegtyp !vegetation physiology type
  real,                            intent(in) :: fveg   !greeness vegetation fraction (-)
  real,                            intent(in) :: sav    !solar rad absorbed by veg (w/m2)
  real,                            intent(in) :: sag    !solar rad absorbed by ground (w/m2)
  real,                            intent(in) :: lwdn   !atmospheric longwave radiation (w/m2)
  real,                            intent(in) :: ur     !wind speed at height zlvl (m/s)
  real,                            intent(in) :: uu     !wind speed in eastward dir (m/s)
  real,                            intent(in) :: vv     !wind speed in northward dir (m/s)
  real,                            intent(in) :: sfctmp !air temperature at reference height (k)
  real,                            intent(in) :: thair  !potential temp at reference height (k)
  real,                            intent(in) :: eair   !vapor pressure air at zlvl (pa)
  real,                            intent(in) :: qair   !specific humidity at zlvl (kg/kg)
  real,                            intent(in) :: rhoair !density air (kg/m**3)
  real,                            intent(in) :: dt     !time step (s)
  real,                            intent(in) :: fsno     !snow fraction

  real,                            intent(in) :: snowh  !actual snow depth [m]
  real,                            intent(in) :: fwet   !wetted fraction of canopy
  real,                            intent(in) :: cwp    !canopy wind parameter

  real,                            intent(in) :: vai    !total leaf area index + stem area index
  real,                            intent(in) :: laisun !sunlit leaf area index, one-sided (m2/m2)
  real,                            intent(in) :: laisha !shaded leaf area index, one-sided (m2/m2)
  real,                            intent(in) :: zlvl   !reference height (m)
  real,                            intent(in) :: cpfac  !heat capacity enhancement factor due to heat storage

  real,                            intent(in) :: zpd    !zero plane displacement (m)
  real,                            intent(in) :: z0m    !roughness length, momentum (m)
  real,                            intent(in) :: z0mg   !roughness length, momentum, ground (m)
  real,                            intent(in) :: emv    !vegetation emissivity
  real,                            intent(in) :: emg    !ground emissivity

  real, dimension(-nsnow+1:nsoil), intent(in) :: stc    !soil/snow temperature (k)
  real, dimension(-nsnow+1:nsoil), intent(in) :: df     !thermal conductivity of snow/soil (w/m/k)
  real, dimension(-nsnow+1:nsoil), intent(in) :: dzsnso !thinkness of snow/soil layers (m)
  real,                            intent(in) :: canliq !intercepted liquid water (mm)
  real,                            intent(in) :: canice !intercepted ice mass (mm)
  real,                            intent(in) :: rsurf  !ground surface resistance (s/m)
!  real,                            intent(in) :: gamma  !psychrometric constant (pa/k)
!  real,                            intent(in) :: lathea !latent heat of vaporization/subli (j/kg)
  real,                            intent(in) :: gammav  !psychrometric constant (pa/k)
  real,                            intent(in) :: latheav !latent heat of vaporization/subli (j/kg)
  real,                            intent(in) :: gammag  !psychrometric constant (pa/k)
  real,                            intent(in) :: latheag !latent heat of vaporization/subli (j/kg)
  real,                            intent(in) :: parsun !par absorbed per unit sunlit lai (w/m2)
  real,                            intent(in) :: parsha !par absorbed per unit shaded lai (w/m2)
  real,                            intent(in) :: foln   !foliage nitrogen (%)
  real,                            intent(in) :: co2air !atmospheric co2 concentration (pa)
  real,                            intent(in) :: o2air  !atmospheric o2 concentration (pa)
  real,                            intent(in) :: igs    !growing season index (0=off, 1=on)
  real,                            intent(in) :: sfcprs !pressure (pa)
  real,                            intent(in) :: btran  !soil water transpiration factor (0 to 1)
  real,                            intent(in) :: rhsur  !raltive humidity in surface soil/snow air space (-)

  real                           , intent(in) :: qc     !cloud water mixing ratio
  real                           , intent(in) :: psfc   !pressure at lowest model layer
  real                           , intent(in) :: dx     !grid spacing
  real                           , intent(in) :: q2     !mixing ratio (kg/kg)
  real                           , intent(in) :: dz8w   !thickness of lowest layer
  real                           , intent(inout) :: qsfc   !mixing ratio at lowest model layer
  real, intent(in)   :: pahv  !precipitation advected heat - canopy net in (w/m2)
  real, intent(in)   :: pahg  !precipitation advected heat - ground net in (w/m2)

! input/output
  real,                         intent(inout) :: eah    !canopy air vapor pressure (pa)
  real,                         intent(inout) :: tah    !canopy air temperature (k)
  real,                         intent(inout) :: tv     !vegetation temperature (k)
  real,                         intent(inout) :: tg     !ground temperature (k)
  real,                         intent(inout) :: cm     !momentum drag coefficient
  real,                         intent(inout) :: ch     !sensible heat exchange coefficient

#ifdef CCPP
  character(len=*),             intent(inout) :: errmsg
  integer,                      intent(inout) :: errflg
#endif

! output
! -fsa + fira + fsh + (fcev + fctr + fgev) + fcst + ssoil = 0
  real,                           intent(out) :: tauxv  !wind stress: e-w (n/m2)
  real,                           intent(out) :: tauyv  !wind stress: n-s (n/m2)
  real,                           intent(out) :: irc    !net longwave radiation (w/m2) [+= to atm]
  real,                           intent(out) :: shc    !sensible heat flux (w/m2)     [+= to atm]
  real,                           intent(out) :: evc    !evaporation heat flux (w/m2)  [+= to atm]
  real,                           intent(out) :: irg    !net longwave radiation (w/m2) [+= to atm]
  real,                           intent(out) :: shg    !sensible heat flux (w/m2)     [+= to atm]
  real,                           intent(out) :: evg    !evaporation heat flux (w/m2)  [+= to atm]
  real,                           intent(out) :: tr     !transpiration heat flux (w/m2)[+= to atm]
  real,                           intent(out) :: gh     !ground heat (w/m2) [+ = to soil]
  real,                           intent(out) :: t2mv   !2 m height air temperature (k)
  real,                           intent(out) :: psnsun !sunlit leaf photosynthesis (umolco2/m2/s)
  real,                           intent(out) :: psnsha !shaded leaf photosynthesis (umolco2/m2/s)
  real,                           intent(out) :: chleaf !leaf exchange coefficient
  real,                           intent(out) :: chuc   !under canopy exchange coefficient

  real,                           intent(out) :: q2v
  real :: cah    !sensible heat conductance, canopy air to zlvl air (m/s)
  real :: u10v    !10 m wind speed in eastward dir (m/s) 
  real :: v10v    !10 m wind speed in eastward dir (m/s) 
  real :: wspd

! ------------------------ local variables ----------------------------------------------------
  real :: cw           !water vapor exchange coefficient
  real :: fv           !friction velocity (m/s)
  real :: wstar        !friction velocity n vertical direction (m/s) (only for sfcdif2)
  real :: z0h          !roughness length, sensible heat (m)
  real :: z0hg         !roughness length, sensible heat (m)
  real :: rb           !bulk leaf boundary layer resistance (s/m)
  real :: ramc         !aerodynamic resistance for momentum (s/m)
  real :: rahc         !aerodynamic resistance for sensible heat (s/m)
  real :: rawc         !aerodynamic resistance for water vapor (s/m)
  real :: ramg         !aerodynamic resistance for momentum (s/m)
  real :: rahg         !aerodynamic resistance for sensible heat (s/m)
  real :: rawg         !aerodynamic resistance for water vapor (s/m)

  real, intent(out) :: rssun        !sunlit leaf stomatal resistance (s/m)
  real, intent(out) :: rssha        !shaded leaf stomatal resistance (s/m)

  real :: mol          !monin-obukhov length (m)
  real :: dtv          !change in tv, last iteration (k)
  real :: dtg          !change in tg, last iteration (k)

  real :: air,cir      !coefficients for ir as function of ts**4
  real :: csh          !coefficients for sh as function of ts
  real :: cev          !coefficients for ev as function of esat[ts]
  real :: cgh          !coefficients for st as function of ts
  real :: atr,ctr      !coefficients for tr as function of esat[ts]
  real :: ata,bta      !coefficients for tah as function of ts
  real :: aea,bea      !coefficients for eah as function of esat[ts]

  real :: estv         !saturation vapor pressure at tv (pa)
  real :: estg         !saturation vapor pressure at tg (pa)
  real :: destv        !d(es)/dt at ts (pa/k)
  real :: destg        !d(es)/dt at tg (pa/k)
  real :: esatw        !es for water
  real :: esati        !es for ice
  real :: dsatw        !d(es)/dt at tg (pa/k) for water
  real :: dsati        !d(es)/dt at tg (pa/k) for ice

  real :: fm           !momentum stability correction, weighted by prior iters
  real :: fh           !sen heat stability correction, weighted by prior iters
  real :: fhg          !sen heat stability correction, ground
  real :: hcan         !canopy height (m) [note: hcan >= z0mg]

  real :: a            !temporary calculation
  real :: b            !temporary calculation
  real :: cvh          !sensible heat conductance, leaf surface to canopy air (m/s)
  real :: caw          !latent heat conductance, canopy air zlvl air (m/s)
  real :: ctw          !transpiration conductance, leaf to canopy air (m/s)
  real :: cew          !evaporation conductance, leaf to canopy air (m/s)
  real :: cgw          !latent heat conductance, ground to canopy air (m/s)
  real :: cond         !sum of conductances (s/m)
  real :: uc           !wind speed at top of canopy (m/s)
  real :: kh           !turbulent transfer coefficient, sensible heat, (m2/s)
  real :: h            !temporary sensible heat flux (w/m2)
  real :: hg           !temporary sensible heat flux (w/m2)

  real :: moz          !monin-obukhov stability parameter
  real :: mozg         !monin-obukhov stability parameter
  real :: mozold       !monin-obukhov stability parameter from prior iteration
  real :: fm2          !monin-obukhov momentum adjustment at 2m
  real :: fh2          !monin-obukhov heat adjustment at 2m
  real :: ch2          !surface exchange at 2m
  real :: thstar          !surface exchange at 2m

  real :: thvair
  real :: thah 
  real :: rahc2        !aerodynamic resistance for sensible heat (s/m)
  real :: rawc2        !aerodynamic resistance for water vapor (s/m)
  real, intent(out):: cah2         !sensible heat conductance for diagnostics
  real :: ch2v         !exchange coefficient for 2m over vegetation. 
  real :: cq2v         !exchange coefficient for 2m over vegetation. 
  real :: eah2         !2m vapor pressure over canopy
  real :: qfx        !moisture flux
  real :: e1           


  real :: vaie         !total leaf area index + stem area index,effective
  real :: laisune      !sunlit leaf area index, one-sided (m2/m2),effective
  real :: laishae      !shaded leaf area index, one-sided (m2/m2),effective

  integer :: k         !index
  integer :: iter      !iteration index

!jref - niterc test from 5 to 20  
  integer, parameter :: niterc = 20   !number of iterations for surface temperature
!jref - niterg test from 3-5
  integer, parameter :: niterg = 5   !number of iterations for ground temperature
  integer :: mozsgn    !number of times moz changes sign
  real    :: mpe       !prevents overflow error if division by zero

  integer :: liter     !last iteration


  real :: t, tdc       !kelvin to degree celsius with limit -50 to +50

  character(len=80) ::  message

  tdc(t)   = min( 50., max(-50.,(t-tfrz)) )
! ---------------------------------------------------------------------------------------------

        mpe = 1e-6
        liter = 0
        fv = 0.1

! ---------------------------------------------------------------------------------------------
! initialization variables that do not depend on stability iteration
! ---------------------------------------------------------------------------------------------
        dtv = 0.
        dtg = 0.
        moz    = 0.
        mozsgn = 0
        mozold = 0.
        hg     = 0.
        h      = 0.
        qfx    = 0.

! YRQ
!       write(*,*) 'tv,tg,stc in input:YRQ', tv,tg,stc

! convert grid-cell lai to the fractional vegetated area (fveg)

        vaie    = min(6.,vai    / fveg)
        laisune = min(6.,laisun / fveg)
        laishae = min(6.,laisha / fveg)

! saturation vapor pressure at ground temperature

        t = tdc(tg)
        call esat(t, esatw, esati, dsatw, dsati)
        if (t .gt. 0.) then
           estg = esatw
        else
           estg = esati
        end if

!jref - consistent surface specific humidity for sfcdif3 and sfcdif4

        qsfc = 0.622*eair/(psfc-0.378*eair)  

! canopy height

        hcan = parameters%hvt
        uc = ur*log(hcan/z0m)/log(zlvl/z0m)
        uc = ur*log((hcan-zpd+z0m)/z0m)/log(zlvl/z0m)   ! mb: add zpd v3.7
        if((hcan-zpd) <= 0.) then
          write(message,*) "critical problem: hcan <= zpd"
#ifdef CCPP
          errmsg = trim(message)
#else
          call wrf_message ( message )
#endif
          write(message,*) 'i,j point=',iloc, jloc
#ifdef CCPP
          errmsg = trim(errmsg)//NEW_LINE('A')//trim(message)
#else
          call wrf_message ( message )
#endif
          write(message,*) 'hcan  =',hcan
#ifdef CCPP
          errmsg = trim(errmsg)//NEW_LINE('A')//trim(message)
#else
          call wrf_message ( message )
#endif
          write(message,*) 'zpd   =',zpd
#ifdef CCPP
          errmsg = trim(errmsg)//NEW_LINE('A')//trim(message)
#else
          call wrf_message ( message )
#endif
          write (message, *) 'snowh =',snowh
#ifdef CCPP
          errflg = 1
          errmsg = trim(errmsg)//NEW_LINE('A')//trim(message)//NEW_LINE('A')//"critical problem in module_sf_noahmplsm:vegeflux"
          return
#else
          call wrf_message ( message )
          call wrf_error_fatal ( "critical problem in module_sf_noahmplsm:vegeflux" )
#endif
          
        end if

! prepare for longwave rad.

        air = -emv*(1.+(1.-emv)*(1.-emg))*lwdn - emv*emg*sb*tg**4  
        cir = (2.-emv*(1.-emg))*emv*sb

! ---------------------------------------------------------------------------------------------
      loop1: do iter = 1, niterc    !  begin stability iteration

       if(iter == 1) then
            z0h  = z0m  
            z0hg = z0mg
       else
            z0h  = z0m    !* exp(-czil*0.4*258.2*sqrt(fv*z0m))
            z0hg = z0mg   !* exp(-czil*0.4*258.2*sqrt(fv*z0mg))
       end if

! aerodyn resistances between heights zlvl and d+z0v

       if(opt_sfc == 1) then
          call sfcdif1(parameters,iter   ,sfctmp ,rhoair ,h      ,qair   , & !in
                       zlvl   ,zpd    ,z0m    ,z0h    ,ur     , & !in
                       mpe    ,iloc   ,jloc   ,                 & !in
#ifdef CCPP
                       moz ,mozsgn ,fm ,fh ,fm2 ,fh2 ,errmsg ,errflg ,& !inout
#else
                       moz ,mozsgn ,fm ,fh ,fm2 ,fh2 ,           & !inout
#endif
                       cm     ,ch     ,fv     ,ch2     )          !out
#ifdef CCPP
          if (errflg /= 0) return
#endif
       endif
     
       if(opt_sfc == 2) then
          call sfcdif2(parameters,iter   ,z0m    ,tah    ,thair  ,ur     , & !in
                       zlvl   ,iloc   ,jloc   ,         & !in
                       cm     ,ch     ,moz    ,wstar  ,         & !in
                       fv     )                                   !out
          ! undo the multiplication by windspeed that sfcdif2 
          ! applies to exchange coefficients ch and cm:
          ch = ch / ur
          cm = cm / ur
       endif

       ramc = max(1.,1./(cm*ur))
       rahc = max(1.,1./(ch*ur))
       rawc = rahc

! aerodyn resistance between heights z0g and d+z0v, rag, and leaf
! boundary layer resistance, rb
       
       call ragrb(parameters,iter   ,vaie   ,rhoair ,hg     ,tah    , & !in
                  zpd    ,z0mg   ,z0hg   ,hcan   ,uc     , & !in
                  z0h    ,fv     ,cwp    ,vegtyp ,mpe    , & !in
                  tv     ,mozg   ,fhg    ,iloc   ,jloc   , & !inout
                  ramg   ,rahg   ,rawg   ,rb     )           !out

! es and d(es)/dt evaluated at tv

       t = tdc(tv)
       call esat(t, esatw, esati, dsatw, dsati)
       if (t .gt. 0.) then
          estv  = esatw
          destv = dsatw
       else
          estv  = esati
          destv = dsati
       end if

! stomatal resistance
        
     if(iter == 1) then
        if (opt_crs == 1) then  ! ball-berry
         call stomata (parameters,vegtyp,mpe   ,parsun ,foln  ,iloc  , jloc , & !in       
                       tv    ,estv  ,eah    ,sfctmp,sfcprs, & !in
                       o2air ,co2air,igs    ,btran ,rb    , & !in
                       rssun ,psnsun)                         !out

         call stomata (parameters,vegtyp,mpe   ,parsha ,foln  ,iloc  , jloc , & !in
                       tv    ,estv  ,eah    ,sfctmp,sfcprs, & !in
                       o2air ,co2air,igs    ,btran ,rb    , & !in
                       rssha ,psnsha)                         !out
        end if

        if (opt_crs == 2) then  ! jarvis
         call  canres (parameters,parsun,tv    ,btran ,eah    ,sfcprs, & !in
                       rssun ,psnsun,iloc  ,jloc   )          !out

         call  canres (parameters,parsha,tv    ,btran ,eah    ,sfcprs, & !in
                       rssha ,psnsha,iloc  ,jloc   )          !out
        end if
     end if

! prepare for sensible heat flux above veg.

        cah  = 1./rahc
        cvh  = 2.*vaie/rb
        cgh  = 1./rahg
        cond = cah + cvh + cgh
        ata  = (sfctmp*cah + tg*cgh) / cond
        bta  = cvh/cond
        csh  = (1.-bta)*rhoair*cpair*cpfac*cvh

! prepare for latent heat flux above veg.

        caw  = 1./rawc
        cew  = fwet*vaie/rb
        ctw  = (1.-fwet)*(laisune/(rb+rssun) + laishae/(rb+rssha))
        cgw  = 1./(rawg+rsurf)
        cond = caw + cew + ctw + cgw
        aea  = (eair*caw + estg*cgw) / cond
        bea  = (cew+ctw)/cond
        cev  = (1.-bea)*cew*rhoair*cpair*cpfac/gammav   ! barlage: change to vegetation v3.6
        ctr  = (1.-bea)*ctw*rhoair*cpair*cpfac/gammav

! evaluate surface fluxes with current temperature and solve for dts

        tah = ata + bta*tv               ! canopy air t.
        eah = aea + bea*estv             ! canopy air e

        irc = fveg*(air + cir*tv**4)
        shc = fveg*rhoair*cpair*cpfac*cvh * (  tv-tah)
        evc = fveg*rhoair*cpair*cpfac*cew * (estv-eah) / gammav ! barlage: change to v in v3.6
        tr  = fveg*rhoair*cpair*cpfac*ctw * (estv-eah) / gammav
	if (tv > tfrz) then
          evc = min(canliq*latheav/dt,evc)    ! barlage: add if block for canice in v3.6
	else
          evc = min(canice*latheav/dt,evc)
	end if

        b   = sav-irc-shc-evc-tr+pahv                          !additional w/m2
        a   = fveg*(4.*cir*tv**3 + csh + (cev+ctr)*destv) !volumetric heat capacity
        dtv = b/a

        irc = irc + fveg*4.*cir*tv**3*dtv
        shc = shc + fveg*csh*dtv
        evc = evc + fveg*cev*destv*dtv
        tr  = tr  + fveg*ctr*destv*dtv                               

! update vegetation surface temperature
        tv  = tv + dtv
!        tah = ata + bta*tv               ! canopy air t; update here for consistency

! for computing m-o length in the next iteration
        h  = rhoair*cpair*(tah - sfctmp) /rahc        
        hg = rhoair*cpair*(tg  - tah)   /rahg

! consistent specific humidity from canopy air vapor pressure
        qsfc = (0.622*eah)/(sfcprs-0.378*eah)

        if (liter == 1) then
           exit loop1 
        endif
        if (iter >= 5 .and. abs(dtv) <= 0.01 .and. liter == 0) then
           liter = 1
        endif

     end do loop1 ! end stability iteration

! under-canopy fluxes and tg

        air = - emg*(1.-emv)*lwdn - emg*emv*sb*tv**4
        cir = emg*sb
        csh = rhoair*cpair*cpfac/rahg
        cev = rhoair*cpair*cpfac / (gammag*(rawg+rsurf))  ! barlage: change to ground v3.6
        cgh = 2.*df(isnow+1)/dzsnso(isnow+1)
!        write(*,*)'inside tg=',tg,'stc(1)=',stc(1)

     loop2: do iter = 1, niterg

        t = tdc(tg)
        call esat(t, esatw, esati, dsatw, dsati)
        if (t .gt. 0.) then
            estg  = esatw
            destg = dsatw
        else
            estg  = esati
            destg = dsati
        end if

        irg = cir*tg**4 + air
        shg = csh * (tg         - tah         )
        evg = cev * (estg*rhsur - eah         )
        gh  = cgh * (tg         - stc(isnow+1))

        b = sag-irg-shg-evg-gh+pahg
        a = 4.*cir*tg**3+csh+cev*destg+cgh
        dtg = b/a

        irg = irg + 4.*cir*tg**3*dtg
        shg = shg + csh*dtg
        evg = evg + cev*destg*dtg
        gh  = gh  + cgh*dtg
        tg  = tg  + dtg

     end do loop2
     
!     tah = (cah*sfctmp + cvh*tv + cgh*tg)/(cah + cvh + cgh)

! if snow on ground and tg > tfrz: reset tg = tfrz. reevaluate ground fluxes.

     if(opt_stc == 1 .or. opt_stc == 3) then
     if (snowh > 0.05 .and. tg > tfrz) then
        tg  = tfrz
        if(opt_stc == 3) tg  = (1.-fsno)*tg + fsno*tfrz   ! mb: allow tg>0c during melt v3.7
        irg = cir*tg**4 - emg*(1.-emv)*lwdn - emg*emv*sb*tv**4
        shg = csh * (tg         - tah)
        evg = cev * (estg*rhsur - eah)
        gh  = sag+pahg - (irg+shg+evg)
     end if
     end if

! wind stresses

     tauxv = -rhoair*cm*ur*uu
     tauyv = -rhoair*cm*ur*vv

! consistent vegetation air temperature and vapor pressure since tg is not consistent with the tah/eah
! calculation.
!     tah = sfctmp + (shg+shc)/(rhoair*cpair*cpfac*cah) 
!     tah = sfctmp + (shg*fveg+shc)/(rhoair*cpair*cpfac*cah) ! ground flux need fveg
!     eah = eair + (evc+fveg*(tr+evg))/(rhoair*caw*cpair*cpfac/gammag )
!     qfx = (qsfc-qair)*rhoair*cpfac*caw !*cpair/gammag

! 2m temperature over vegetation ( corrected for low cq2v values )
   if (opt_sfc == 1 .or. opt_sfc == 2) then
!      cah2 = fv*1./vkc*log((2.+z0h)/z0h)
      cah2 = fv*vkc/log((2.+z0h)/z0h)
      cah2 = fv*vkc/(log((2.+z0h)/z0h)-fh2)
      cq2v = cah2
      if (cah2 .lt. 1.e-5 ) then
         t2mv = tah
!         q2v  = (eah*0.622/(sfcprs - 0.378*eah))
         q2v  = qsfc
      else
         t2mv = tah - (shg+shc/fveg)/(rhoair*cpair*cpfac) * 1./cah2
!         q2v = (eah*0.622/(sfcprs - 0.378*eah))- qfx/(rhoair*fv)* 1./vkc * log((2.+z0h)/z0h)
         q2v = qsfc - ((evc+tr)/fveg+evg)/(latheav*rhoair) * 1./cq2v
      endif
   endif

! update ch for output
     ch = cah
     chleaf = cvh
     chuc = 1./rahg

  end subroutine vege_flux

!== begin bare_flux ================================================================================

!>\ingroup NoahMP_LSM
  subroutine bare_flux (parameters,nsnow   ,nsoil   ,isnow   ,dt      ,sag     , & !in
                        lwdn    ,ur      ,uu      ,vv      ,sfctmp  , & !in
                        thair   ,qair    ,eair    ,rhoair  ,snowh   , & !in
                        dzsnso  ,zlvl    ,zpd     ,z0m     ,fsno    , & !in
                        emg     ,stc     ,df      ,rsurf   ,lathea  , & !in
                        gamma   ,rhsur   ,iloc    ,jloc    ,q2      ,pahb  , & !in
#ifdef CCPP
                        tgb     ,cm      ,ch      ,errmsg  ,errflg  , & !inout
#else
                        tgb     ,cm      ,ch      ,          & !inout
#endif
                        tauxb   ,tauyb   ,irb     ,shb     ,evb     , & !out
                        ghb     ,t2mb    ,dx      ,dz8w    ,ivgtyp  , & !out
                        qc      ,qsfc    ,psfc    ,                   & !in
                        sfcprs  ,q2b     ,ehb2    )                     !in 

! --------------------------------------------------------------------------------------------------
! use newton-raphson iteration to solve ground (tg) temperature
! that balances the surface energy budgets for bare soil fraction.

! bare soil:
! -sab + irb[tg] + shb[tg] + evb[tg] + ghb[tg] = 0
! ----------------------------------------------------------------------
  implicit none
! ----------------------------------------------------------------------
! input
  type (noahmp_parameters), intent(in) :: parameters
  integer                        , intent(in) :: iloc   !grid index
  integer                        , intent(in) :: jloc   !grid index
  integer,                         intent(in) :: nsnow  !maximum no. of snow layers
  integer,                         intent(in) :: nsoil  !number of soil layers
  integer,                         intent(in) :: isnow  !actual no. of snow layers
  real,                            intent(in) :: dt     !time step (s)
  real,                            intent(in) :: sag    !solar radiation absorbed by ground (w/m2)
  real,                            intent(in) :: lwdn   !atmospheric longwave radiation (w/m2)
  real,                            intent(in) :: ur     !wind speed at height zlvl (m/s)
  real,                            intent(in) :: uu     !wind speed in eastward dir (m/s)
  real,                            intent(in) :: vv     !wind speed in northward dir (m/s)
  real,                            intent(in) :: sfctmp !air temperature at reference height (k)
  real,                            intent(in) :: thair  !potential temperature at height zlvl (k)
  real,                            intent(in) :: qair   !specific humidity at height zlvl (kg/kg)
  real,                            intent(in) :: eair   !vapor pressure air at height (pa)
  real,                            intent(in) :: rhoair !density air (kg/m3)
  real,                            intent(in) :: snowh  !actual snow depth [m]
  real, dimension(-nsnow+1:nsoil), intent(in) :: dzsnso !thickness of snow/soil layers (m)
  real,                            intent(in) :: zlvl   !reference height (m)
  real,                            intent(in) :: zpd    !zero plane displacement (m)
  real,                            intent(in) :: z0m    !roughness length, momentum, ground (m)
  real,                            intent(in) :: emg    !ground emissivity
  real, dimension(-nsnow+1:nsoil), intent(in) :: stc    !soil/snow temperature (k)
  real, dimension(-nsnow+1:nsoil), intent(in) :: df     !thermal conductivity of snow/soil (w/m/k)
  real,                            intent(in) :: rsurf  !ground surface resistance (s/m)
  real,                            intent(in) :: lathea !latent heat of vaporization/subli (j/kg)
  real,                            intent(in) :: gamma  !psychrometric constant (pa/k)
  real,                            intent(in) :: rhsur  !raltive humidity in surface soil/snow air space (-)
  real,                            intent(in) :: fsno     !snow fraction

!jref:start; in 
  integer                        , intent(in) :: ivgtyp
  real                           , intent(in) :: qc     !cloud water mixing ratio
  real                           , intent(inout) :: qsfc   !mixing ratio at lowest model layer
  real                           , intent(in) :: psfc   !pressure at lowest model layer
  real                           , intent(in) :: sfcprs !pressure at lowest model layer
  real                           , intent(in) :: dx     !horisontal grid spacing
  real                           , intent(in) :: q2     !mixing ratio (kg/kg)
  real                           , intent(in) :: dz8w   !thickness of lowest layer
!jref:end
  real, intent(in)   :: pahb  !precipitation advected heat - ground net in (w/m2)

! input/output
  real,                         intent(inout) :: tgb    !ground temperature (k)
  real,                         intent(inout) :: cm     !momentum drag coefficient
  real,                         intent(inout) :: ch     !sensible heat exchange coefficient
#ifdef CCPP
  character(len=*),             intent(inout) :: errmsg
  integer,                      intent(inout) :: errflg
#endif

! output
! -sab + irb[tg] + shb[tg] + evb[tg] + ghb[tg] = 0

  real,                           intent(out) :: tauxb  !wind stress: e-w (n/m2)
  real,                           intent(out) :: tauyb  !wind stress: n-s (n/m2)
  real,                           intent(out) :: irb    !net longwave rad (w/m2)   [+ to atm]
  real,                           intent(out) :: shb    !sensible heat flux (w/m2) [+ to atm]
  real,                           intent(out) :: evb    !latent heat flux (w/m2)   [+ to atm]
  real,                           intent(out) :: ghb    !ground heat flux (w/m2)  [+ to soil]
  real,                           intent(out) :: t2mb   !2 m height air temperature (k)
!jref:start
  real,                           intent(out) :: q2b    !bare ground heat conductance
  real :: ehb    !bare ground heat conductance
  real :: u10b    !10 m wind speed in eastward dir (m/s)
  real :: v10b    !10 m wind speed in eastward dir (m/s)
  real :: wspd
!jref:end

! local variables 

  real :: taux       !wind stress: e-w (n/m2)
  real :: tauy       !wind stress: n-s (n/m2)
  real :: fira       !total net longwave rad (w/m2)      [+ to atm]
  real :: fsh        !total sensible heat flux (w/m2)    [+ to atm]
  real :: fgev       !ground evaporation heat flux (w/m2)[+ to atm]
  real :: ssoil      !soil heat flux (w/m2)             [+ to soil]
  real :: fire       !emitted ir (w/m2)
  real :: trad       !radiative temperature (k)
  real :: tah        !"surface" temperature at height z0h+zpd (k)

  real :: cw         !water vapor exchange coefficient
  real :: fv         !friction velocity (m/s)
  real :: wstar      !friction velocity n vertical direction (m/s) (only for sfcdif2)
  real :: z0h        !roughness length, sensible heat, ground (m)
  real :: rb         !bulk leaf boundary layer resistance (s/m)
  real :: ramb       !aerodynamic resistance for momentum (s/m)
  real :: rahb       !aerodynamic resistance for sensible heat (s/m)
  real :: rawb       !aerodynamic resistance for water vapor (s/m)
  real :: mol        !monin-obukhov length (m)
  real :: dtg        !change in tg, last iteration (k)

  real :: cir        !coefficients for ir as function of ts**4
  real :: csh        !coefficients for sh as function of ts
  real :: cev        !coefficients for ev as function of esat[ts]
  real :: cgh        !coefficients for st as function of ts

!jref:start
  real :: rahb2      !aerodynamic resistance for sensible heat 2m (s/m)
  real :: rawb2      !aerodynamic resistance for water vapor 2m (s/m)
  real,intent(out) :: ehb2       !sensible heat conductance for diagnostics
  real :: ch2b       !exchange coefficient for 2m temp.
  real :: cq2b       !exchange coefficient for 2m temp.
  real :: thvair     !virtual potential air temp
  real :: thgh       !potential ground temp
  real :: emb        !momentum conductance
  real :: qfx        !moisture flux
  real :: estg2      !saturation vapor pressure at 2m (pa)
  integer :: vegtyp     !vegetation type set to isbarren
  real :: e1
!jref:end

  real :: estg       !saturation vapor pressure at tg (pa)
  real :: destg      !d(es)/dt at tg (pa/k)
  real :: esatw      !es for water
  real :: esati      !es for ice
  real :: dsatw      !d(es)/dt at tg (pa/k) for water
  real :: dsati      !d(es)/dt at tg (pa/k) for ice

  real :: a          !temporary calculation
  real :: b          !temporary calculation
  real :: h          !temporary sensible heat flux (w/m2)
  real :: moz        !monin-obukhov stability parameter
  real :: mozold     !monin-obukhov stability parameter from prior iteration
  real :: fm         !momentum stability correction, weighted by prior iters
  real :: fh         !sen heat stability correction, weighted by prior iters
  integer :: mozsgn  !number of times moz changes sign
  real :: fm2          !monin-obukhov momentum adjustment at 2m
  real :: fh2          !monin-obukhov heat adjustment at 2m
  real :: ch2          !surface exchange at 2m

  integer :: iter    !iteration index
  integer :: niterb  !number of iterations for surface temperature
  real    :: mpe     !prevents overflow error if division by zero
!jref:start
!  data niterb /3/
  data niterb /5/
  save niterb
  real :: t, tdc     !kelvin to degree celsius with limit -50 to +50
  tdc(t)   = min( 50., max(-50.,(t-tfrz)) )

! -----------------------------------------------------------------
! initialization variables that do not depend on stability iteration
! -----------------------------------------------------------------
        mpe = 1e-6
        dtg = 0.
        moz    = 0.
        mozsgn = 0
        mozold = 0.
        h      = 0.
        qfx    = 0.
        fv     = 0.1

        cir = emg*sb
        cgh = 2.*df(isnow+1)/dzsnso(isnow+1)

! -----------------------------------------------------------------
      loop3: do iter = 1, niterb  ! begin stability iteration

        if(iter == 1) then
            z0h = z0m 
        else
            z0h = z0m !* exp(-czil*0.4*258.2*sqrt(fv*z0m))
        end if

        if(opt_sfc == 1) then
          call sfcdif1(parameters,iter   ,sfctmp ,rhoair ,h      ,qair   , & !in
                       zlvl   ,zpd    ,z0m    ,z0h    ,ur     , & !in
                       mpe    ,iloc   ,jloc   ,                 & !in
#ifdef CCPP
                       moz ,mozsgn ,fm ,fh ,fm2 ,fh2 ,errmsg ,errflg ,& !inout
#else
                       moz ,mozsgn ,fm ,fh ,fm2 ,fh2 ,           & !inout
#endif
                       cm     ,ch     ,fv     ,ch2     )          !out
#ifdef CCPP
          if (errflg /= 0) return
#endif
        endif

        if(opt_sfc == 2) then
          call sfcdif2(parameters,iter   ,z0m    ,tgb    ,thair  ,ur     , & !in
                       zlvl   ,iloc   ,jloc   ,         & !in
                       cm     ,ch     ,moz    ,wstar  ,         & !in
                       fv     )                                   !out
          ! undo the multiplication by windspeed that sfcdif2 
          ! applies to exchange coefficients ch and cm:
          ch = ch / ur
          cm = cm / ur
          if(snowh > 0.) then
             cm = min(0.01,cm)   ! cm & ch are too large, causing
             ch = min(0.01,ch)   ! computational instability
          end if

        endif

        ramb = max(1.,1./(cm*ur))
        rahb = max(1.,1./(ch*ur))
        rawb = rahb

!jref - variables for diagnostics         
        emb = 1./ramb
        ehb = 1./rahb

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

        b     = sag-irb-shb-evb-ghb+pahb
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
        qsfc = 0.622*(estg*rhsur)/(psfc-0.378*(estg*rhsur))

        qfx = (qsfc-qair)*cev*gamma/cpair

     end do loop3 ! end stability iteration
! -----------------------------------------------------------------

! if snow on ground and tg > tfrz: reset tg = tfrz. reevaluate ground fluxes.

     if(opt_stc == 1 .or. opt_stc == 3) then
     if (snowh > 0.05 .and. tgb > tfrz) then
          tgb = tfrz
          if(opt_stc == 3) tgb  = (1.-fsno)*tgb + fsno*tfrz  ! mb: allow tg>0c during melt v3.7
          irb = cir * tgb**4 - emg*lwdn
          shb = csh * (tgb        - sfctmp)
          evb = cev * (estg*rhsur - eair )          !estg reevaluate ?
          ghb = sag+pahb - (irb+shb+evb)
     end if
     end if

! wind stresses
         
     tauxb = -rhoair*cm*ur*uu
     tauyb = -rhoair*cm*ur*vv

!jref:start; errors in original equation corrected.
! 2m air temperature
     if(opt_sfc == 1 .or. opt_sfc ==2) then
       ehb2  = fv*vkc/log((2.+z0h)/z0h)
       ehb2  = fv*vkc/(log((2.+z0h)/z0h)-fh2)
       cq2b  = ehb2
       if (ehb2.lt.1.e-5 ) then
         t2mb  = tgb
         q2b   = qsfc
       else
         t2mb  = tgb - shb/(rhoair*cpair) * 1./ehb2
         q2b   = qsfc - evb/(lathea*rhoair)*(1./cq2b + rsurf)
       endif
       if (parameters%urban_flag) q2b = qsfc
     end if

! update ch 
     ch = ehb

  end subroutine bare_flux

!== begin ragrb ====================================================================================

!>\ingroup NoahMP_LSM
  subroutine ragrb(parameters,iter   ,vai    ,rhoair ,hg     ,tah    , & !in
                   zpd    ,z0mg   ,z0hg   ,hcan   ,uc     , & !in
                   z0h    ,fv     ,cwp    ,vegtyp ,mpe    , & !in
                   tv     ,mozg   ,fhg    ,iloc   ,jloc   , & !inout
                   ramg   ,rahg   ,rawg   ,rb     )           !out
! --------------------------------------------------------------------------------------------------
! compute under-canopy aerodynamic resistance rag and leaf boundary layer
! resistance rb
! --------------------------------------------------------------------------------------------------
  implicit none
! --------------------------------------------------------------------------------------------------
! inputs

  type (noahmp_parameters), intent(in) :: parameters
  integer,              intent(in) :: iloc   !grid index
  integer,              intent(in) :: jloc   !grid index
  integer,              intent(in) :: iter   !iteration index
  integer,              intent(in) :: vegtyp !vegetation physiology type
  real,                 intent(in) :: vai    !total lai + stem area index, one sided
  real,                 intent(in) :: rhoair !density air (kg/m3)
  real,                 intent(in) :: hg     !ground sensible heat flux (w/m2)
  real,                 intent(in) :: tv     !vegetation temperature (k)
  real,                 intent(in) :: tah    !air temperature at height z0h+zpd (k)
  real,                 intent(in) :: zpd    !zero plane displacement (m)
  real,                 intent(in) :: z0mg   !roughness length, momentum, ground (m)
  real,                 intent(in) :: hcan   !canopy height (m) [note: hcan >= z0mg]
  real,                 intent(in) :: uc     !wind speed at top of canopy (m/s)
  real,                 intent(in) :: z0h    !roughness length, sensible heat (m)
  real,                 intent(in) :: z0hg   !roughness length, sensible heat, ground (m)
  real,                 intent(in) :: fv     !friction velocity (m/s)
  real,                 intent(in) :: cwp    !canopy wind parameter
  real,                 intent(in) :: mpe    !prevents overflow error if division by zero

! in & out

  real,              intent(inout) :: mozg   !monin-obukhov stability parameter
  real,              intent(inout) :: fhg    !stability correction

! outputs
  real                             :: ramg   !aerodynamic resistance for momentum (s/m)
  real                             :: rahg   !aerodynamic resistance for sensible heat (s/m)
  real                             :: rawg   !aerodynamic resistance for water vapor (s/m)
  real                             :: rb     !bulk leaf boundary layer resistance (s/m)


  real :: kh           !turbulent transfer coefficient, sensible heat, (m2/s)
  real :: tmp1         !temporary calculation
  real :: tmp2         !temporary calculation
  real :: tmprah2      !temporary calculation for aerodynamic resistances
  real :: tmprb        !temporary calculation for rb
  real :: molg,fhgnew,cwpc
! --------------------------------------------------------------------------------------------------
! stability correction to below canopy resistance

       mozg = 0.
       molg = 0.

       if(iter > 1) then
        tmp1 = vkc * (grav/tah) * hg/(rhoair*cpair)
        if (abs(tmp1) .le. mpe) tmp1 = mpe
        molg = -1. * fv**3 / tmp1
        mozg = min( (zpd-z0mg)/molg, 1.)
       end if

       if (mozg < 0.) then
          fhgnew  = (1. - 15.*mozg)**(-0.25)
       else
          fhgnew  = 1.+ 4.7*mozg
       endif

       if (iter == 1) then
          fhg = fhgnew
       else
          fhg = 0.5 * (fhg+fhgnew)
       endif

       cwpc = (cwp * vai * hcan * fhg)**0.5
!       cwpc = (cwp*fhg)**0.5

       tmp1 = exp( -cwpc*z0hg/hcan )
       tmp2 = exp( -cwpc*(z0h+zpd)/hcan )
       tmprah2 = hcan*exp(cwpc) / cwpc * (tmp1-tmp2)

! aerodynamic resistances raw and rah between heights zpd+z0h and z0hg.

       kh  = max ( vkc*fv*(hcan-zpd), mpe )
       ramg = 0.
       rahg = tmprah2 / kh
       rawg = rahg

! leaf boundary layer resistance

       tmprb  = cwpc*50. / (1. - exp(-cwpc/2.))
       rb     = tmprb * sqrt(parameters%dleaf/uc)
!       rb = 200

  end subroutine ragrb

!== begin sfcdif1 ==================================================================================

!>\ingroup NoahMP_LSM
  subroutine sfcdif1(parameters,iter   ,sfctmp ,rhoair ,h      ,qair   , & !in
       &             zlvl   ,zpd    ,z0m    ,z0h    ,ur     , & !in
       &             mpe    ,iloc   ,jloc   ,                 & !in
#ifdef CCPP
       &             moz    ,mozsgn ,fm     ,fh     ,fm2,fh2,errmsg,errflg, & !inout
#else
       &             moz    ,mozsgn ,fm     ,fh     ,fm2,fh2, & !inout
#endif
       &             cm     ,ch     ,fv     ,ch2     )          !out
! -------------------------------------------------------------------------------------------------
! computing surface drag coefficient cm for momentum and ch for heat
! -------------------------------------------------------------------------------------------------
    implicit none
! -------------------------------------------------------------------------------------------------
! inputs
    
  type (noahmp_parameters), intent(in) :: parameters
    integer,              intent(in) :: iloc   !grid index
    integer,              intent(in) :: jloc   !grid index
    integer,              intent(in) :: iter   !iteration index
    real,                 intent(in) :: sfctmp !temperature at reference height (k)
    real,                 intent(in) :: rhoair !density air (kg/m**3)
    real,                 intent(in) :: h      !sensible heat flux (w/m2) [+ to atm]
    real,                 intent(in) :: qair   !specific humidity at reference height (kg/kg)
    real,                 intent(in) :: zlvl   !reference height  (m)
    real,                 intent(in) :: zpd    !zero plane displacement (m)
    real,                 intent(in) :: z0h    !roughness length, sensible heat, ground (m)
    real,                 intent(in) :: z0m    !roughness length, momentum, ground (m)
    real,                 intent(in) :: ur     !wind speed (m/s)
    real,                 intent(in) :: mpe    !prevents overflow error if division by zero
! in & out

    integer,           intent(inout) :: mozsgn !number of times moz changes sign
    real,              intent(inout) :: moz    !monin-obukhov stability (z/l)
    real,              intent(inout) :: fm     !momentum stability correction, weighted by prior iters
    real,              intent(inout) :: fh     !sen heat stability correction, weighted by prior iters
    real,              intent(inout) :: fm2    !sen heat stability correction, weighted by prior iters
    real,              intent(inout) :: fh2    !sen heat stability correction, weighted by prior iters
#ifdef CCPP
    character(len=*),  intent(inout) :: errmsg
    integer,           intent(inout) :: errflg
#endif

! outputs

    real,                intent(out) :: cm     !drag coefficient for momentum
    real,                intent(out) :: ch     !drag coefficient for heat
    real,                intent(out) :: fv     !friction velocity (m/s)
    real,                intent(out) :: ch2    !drag coefficient for heat

! locals
    real    :: mol                      !monin-obukhov length (m)
    real    :: tmpcm                    !temporary calculation for cm
    real    :: tmpch                    !temporary calculation for ch
    real    :: fmnew                    !stability correction factor, momentum, for current moz
    real    :: fhnew                    !stability correction factor, sen heat, for current moz
    real    :: mozold                   !monin-obukhov stability parameter from prior iteration
    real    :: tmp1,tmp2,tmp3,tmp4,tmp5 !temporary calculation
    real    :: tvir                     !temporary virtual temperature (k)
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
       write(*,*) 'critical problem: zlvl <= zpd; model stops'
#ifdef CCPP
       errflg = 1
       errmsg = "stop in noah-mp"
       return
#else
      call wrf_error_fatal("stop in noah-mp")
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

  end subroutine sfcdif1

!== begin sfcdif2 ==================================================================================

!>\ingroup NoahMP_LSM
  subroutine sfcdif2(parameters,iter   ,z0     ,thz0   ,thlm   ,sfcspd , & !in
                     zlm    ,iloc   ,jloc   ,         & !in
                     akms   ,akhs   ,rlmo   ,wstar2 ,         & !in
                     ustar  )                                   !out

! -------------------------------------------------------------------------------------------------
! subroutine sfcdif (renamed sfcdif_off to avoid clash with eta pbl)
! -------------------------------------------------------------------------------------------------
! calculate surface layer exchange coefficients via iterative process.
! see chen et al (1997, blm)
! -------------------------------------------------------------------------------------------------
    implicit none
  type (noahmp_parameters), intent(in) :: parameters
    integer, intent(in) :: iloc
    integer, intent(in) :: jloc
    integer, intent(in) :: iter
    real,    intent(in) :: zlm, z0, thz0, thlm, sfcspd
    real, intent(inout) :: akms
    real, intent(inout) :: akhs
    real, intent(inout) :: rlmo
    real, intent(inout) :: wstar2
    real,   intent(out) :: ustar

    real     zz, pslmu, pslms, pslhu, pslhs
    real     xx, pspmu, yy, pspms, psphu, psphs
    real     zilfc, zu, zt, rdz, cxch
    real     dthv, du2, btgh, zslu, zslt, rlogu, rlogt
    real     zetalt, zetalu, zetau, zetat, xlu4, xlt4, xu4, xt4

    real     xlu, xlt, xu, xt, psmz, simm, pshz, simh, ustark, rlmn,  &
         &         rlma

    integer  ilech, itr

    integer, parameter :: itrmx  = 5
    real,    parameter :: wwst   = 1.2
    real,    parameter :: wwst2  = wwst * wwst
    real,    parameter :: vkrm   = 0.40
    real,    parameter :: excm   = 0.001
    real,    parameter :: beta   = 1.0 / 270.0
    real,    parameter :: btg    = beta * grav
    real,    parameter :: elfc   = vkrm * btg
    real,    parameter :: wold   = 0.15
    real,    parameter :: wnew   = 1.0 - wold
    real,    parameter :: pihf   = 3.14159265 / 2.
    real,    parameter :: epsu2  = 1.e-4
    real,    parameter :: epsust = 0.07
    real,    parameter :: epsit  = 1.e-4
    real,    parameter :: epsa   = 1.e-8
    real,    parameter :: ztmin  = -5.0
    real,    parameter :: ztmax  = 1.0
    real,    parameter :: hpbl   = 1000.0
    real,    parameter :: sqvisc = 258.2
    real,    parameter :: ric    = 0.183
    real,    parameter :: rric   = 1.0 / ric
    real,    parameter :: fhneu  = 0.8
    real,    parameter :: rfc    = 0.191
    real,    parameter :: rfac   = ric / ( fhneu * rfc * rfc )

! ----------------------------------------------------------------------
! note: the two code blocks below define functions
! ----------------------------------------------------------------------
! lech's surface functions
    pslmu (zz)= -0.96* log (1.0-4.5* zz)
    pslms (zz)= zz * rric -2.076* (1. -1./ (zz +1.))
    pslhu (zz)= -0.96* log (1.0-4.5* zz)
    pslhs (zz)= zz * rfac -2.076* (1. -1./ (zz +1.))
! paulson's surface functions
    pspmu (xx)= -2.* log ( (xx +1.)*0.5) - log ( (xx * xx +1.)*0.5)   &
         &        +2.* atan (xx)                                            &
         &- pihf
    pspms (yy)= 5.* yy
    psphu (xx)= -2.* log ( (xx * xx +1.)*0.5)
    psphs (yy)= 5.* yy

! this routine sfcdif can handle both over open water (sea, ocean) and
! over solid surface (land, sea-ice).
! ----------------------------------------------------------------------
!     ztfc: ratio of zoh/zom  less or equal than 1
!     c......ztfc=0.1
!     czil: constant c in zilitinkevich, s. s.1995,:note about zt
! ----------------------------------------------------------------------
    ilech = 0

! ----------------------------------------------------------------------
    zilfc = - parameters%czil * vkrm * sqvisc
    zu = z0
    rdz = 1./ zlm
    cxch = excm * rdz
    dthv = thlm - thz0

! beljars correction of ustar
    du2 = max (sfcspd * sfcspd,epsu2)
    btgh = btg * hpbl

    if(iter == 1) then
        if (btgh * akhs * dthv .ne. 0.0) then
           wstar2 = wwst2* abs (btgh * akhs * dthv)** (2./3.)
        else
           wstar2 = 0.0
        end if
        ustar = max (sqrt (akms * sqrt (du2+ wstar2)),epsust)
        rlmo = elfc * akhs * dthv / ustar **3
    end if
 
! zilitinkevitch approach for zt
    zt = max(1.e-6,exp (zilfc * sqrt (ustar * z0))* z0)
    zslu = zlm + zu
    zslt = zlm + zt
    rlogu = log (zslu / zu)
    rlogt = log (zslt / zt)

! ----------------------------------------------------------------------
! 1./monin-obukkhov length-scale
! ----------------------------------------------------------------------
    zetalt = max (zslt * rlmo,ztmin)
    rlmo = zetalt / zslt
    zetalu = zslu * rlmo
    zetau = zu * rlmo
    zetat = zt * rlmo

    if (ilech .eq. 0) then
       if (rlmo .lt. 0.)then
          xlu4 = 1. -16.* zetalu
          xlt4 = 1. -16.* zetalt
          xu4  = 1. -16.* zetau
          xt4  = 1. -16.* zetat
          xlu  = sqrt (sqrt (xlu4))
          xlt  = sqrt (sqrt (xlt4))
          xu   = sqrt (sqrt (xu4))

          xt = sqrt (sqrt (xt4))
          psmz = pspmu (xu)
          simm = pspmu (xlu) - psmz + rlogu
          pshz = psphu (xt)
          simh = psphu (xlt) - pshz + rlogt
       else
          zetalu = min (zetalu,ztmax)
          zetalt = min (zetalt,ztmax)
          psmz = pspms (zetau)
          simm = pspms (zetalu) - psmz + rlogu
          pshz = psphs (zetat)
          simh = psphs (zetalt) - pshz + rlogt
       end if
! ----------------------------------------------------------------------
! lech's functions
! ----------------------------------------------------------------------
    else
       if (rlmo .lt. 0.)then
          psmz = pslmu (zetau)
          simm = pslmu (zetalu) - psmz + rlogu
          pshz = pslhu (zetat)
          simh = pslhu (zetalt) - pshz + rlogt
       else
          zetalu = min (zetalu,ztmax)
          zetalt = min (zetalt,ztmax)
          psmz = pslms (zetau)
          simm = pslms (zetalu) - psmz + rlogu
          pshz = pslhs (zetat)
          simh = pslhs (zetalt) - pshz + rlogt
       end if
! ----------------------------------------------------------------------
       end if

! ----------------------------------------------------------------------
! beljaars correction for ustar
! ----------------------------------------------------------------------
       ustar = max (sqrt (akms * sqrt (du2+ wstar2)),epsust)

! zilitinkevitch fix for zt
       zt = max(1.e-6,exp (zilfc * sqrt (ustar * z0))* z0)
       zslt = zlm + zt
!-----------------------------------------------------------------------
       rlogt = log (zslt / zt)
       ustark = ustar * vkrm
       akms = max (ustark / simm,cxch)
!-----------------------------------------------------------------------
! if statements to avoid tangent linear problems near zero
!-----------------------------------------------------------------------
       akhs = max (ustark / simh,cxch)

       if (btgh * akhs * dthv .ne. 0.0) then
          wstar2 = wwst2* abs (btgh * akhs * dthv)** (2./3.)
       else
          wstar2 = 0.0
       end if
!-----------------------------------------------------------------------
       rlmn = elfc * akhs * dthv / ustar **3
!-----------------------------------------------------------------------
!     if(abs((rlmn-rlmo)/rlma).lt.epsit)    go to 110
!-----------------------------------------------------------------------
       rlma = rlmo * wold+ rlmn * wnew
!-----------------------------------------------------------------------
       rlmo = rlma

!       write(*,'(a20,10f15.6)')'sfcdif: rlmo=',rlmo,rlmn,elfc , akhs , dthv , ustar
!    end do
! ----------------------------------------------------------------------
  end subroutine sfcdif2

!== begin esat =====================================================================================

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

!== begin stomata ==================================================================================

!>\ingroup NoahMP_LSM
  subroutine stomata (parameters,vegtyp  ,mpe     ,apar    ,foln    ,iloc    , jloc, & !in
                      tv      ,ei      ,ea      ,sfctmp  ,sfcprs  , & !in
                      o2      ,co2     ,igs     ,btran   ,rb      , & !in
                      rs      ,psn     )                              !out
! --------------------------------------------------------------------------------------------------
  implicit none
! --------------------------------------------------------------------------------------------------
! input
  type (noahmp_parameters), intent(in) :: parameters
      integer,intent(in)  :: iloc   !grid index
      integer,intent(in)  :: jloc   !grid index
      integer,intent(in)  :: vegtyp !vegetation physiology type

      real, intent(in)    :: igs    !growing season index (0=off, 1=on)
      real, intent(in)    :: mpe    !prevents division by zero errors

      real, intent(in)    :: tv     !foliage temperature (k)
      real, intent(in)    :: ei     !vapor pressure inside leaf (sat vapor press at tv) (pa)
      real, intent(in)    :: ea     !vapor pressure of canopy air (pa)
      real, intent(in)    :: apar   !par absorbed per unit lai (w/m2)
      real, intent(in)    :: o2     !atmospheric o2 concentration (pa)
      real, intent(in)    :: co2    !atmospheric co2 concentration (pa)
      real, intent(in)    :: sfcprs !air pressure at reference height (pa)
      real, intent(in)    :: sfctmp !air temperature at reference height (k)
      real, intent(in)    :: btran  !soil water transpiration factor (0 to 1)
      real, intent(in)    :: foln   !foliage nitrogen concentration (%)
      real, intent(in)    :: rb     !boundary layer resistance (s/m)

! output
      real, intent(out)   :: rs     !leaf stomatal resistance (s/m)
      real, intent(out)   :: psn    !foliage photosynthesis (umol co2 /m2/ s) [always +]

! in&out
      real                :: rlb    !boundary layer resistance (s m2 / umol)
! ---------------------------------------------------------------------------------------------

! ------------------------ local variables ----------------------------------------------------
      integer :: iter     !iteration index
      integer :: niter    !number of iterations

      data niter /3/
      save niter

      real :: ab          !used in statement functions
      real :: bc          !used in statement functions
      real :: f1          !generic temperature response (statement function)
      real :: f2          !generic temperature inhibition (statement function)
      real :: tc          !foliage temperature (degree celsius)
      real :: cs          !co2 concentration at leaf surface (pa)
      real :: kc          !co2 michaelis-menten constant (pa)
      real :: ko          !o2 michaelis-menten constant (pa)
      real :: a,b,c,q     !intermediate calculations for rs
      real :: r1,r2       !roots for rs
      real :: fnf         !foliage nitrogen adjustment factor (0 to 1)
      real :: ppf         !absorb photosynthetic photon flux (umol photons/m2/s)
      real :: wc          !rubisco limited photosynthesis (umol co2/m2/s)
      real :: wj          !light limited photosynthesis (umol co2/m2/s)
      real :: we          !export limited photosynthesis (umol co2/m2/s)
      real :: cp          !co2 compensation point (pa)
      real :: ci          !internal co2 (pa)
      real :: awc         !intermediate calculation for wc
      real :: vcmx        !maximum rate of carbonylation (umol co2/m2/s)
      real :: j           !electron transport (umol co2/m2/s)
      real :: cea         !constrain ea or else model blows up
      real :: cf          !s m2/umol -> s/m

      f1(ab,bc) = ab**((bc-25.)/10.)
      f2(ab) = 1. + exp((-2.2e05+710.*(ab+273.16))/(8.314*(ab+273.16)))
      real :: t
! ---------------------------------------------------------------------------------------------

! initialize rs=rsmax and psn=0 because will only do calculations
! for apar > 0, in which case rs <= rsmax and psn >= 0

         cf = sfcprs/(8.314*sfctmp)*1.e06
         rs = 1./parameters%bp * cf
         psn = 0.

         if (apar .le. 0.) return

         fnf = min( foln/max(mpe,parameters%folnmx), 1.0 )
         tc  = tv-tfrz
         ppf = 4.6*apar
         j   = ppf*parameters%qe25
         kc  = parameters%kc25 * f1(parameters%akc,tc)
         ko  = parameters%ko25 * f1(parameters%ako,tc)
         awc = kc * (1.+o2/ko)
         cp  = 0.5*kc/ko*o2*0.21
         vcmx = parameters%vcmx25 / f2(tc) * fnf * btran * f1(parameters%avcmx,tc)

! first guess ci

         ci = 0.7*co2*parameters%c3psn + 0.4*co2*(1.-parameters%c3psn)

! rb: s/m -> s m**2 / umol

         rlb = rb/cf

! constrain ea

         cea = max(0.25*ei*parameters%c3psn+0.40*ei*(1.-parameters%c3psn), min(ea,ei) )

! ci iteration
!jref: c3psn is equal to 1 for all veg types.
       do iter = 1, niter
            wj = max(ci-cp,0.)*j/(ci+2.*cp)*parameters%c3psn  + j*(1.-parameters%c3psn)
            wc = max(ci-cp,0.)*vcmx/(ci+awc)*parameters%c3psn + vcmx*(1.-parameters%c3psn)
            we = 0.5*vcmx*parameters%c3psn + 4000.*vcmx*ci/sfcprs*(1.-parameters%c3psn)
            psn = min(wj,wc,we) * igs

            cs = max( co2-1.37*rlb*sfcprs*psn, mpe )
            a = parameters%mp*psn*sfcprs*cea / (cs*ei) + parameters%bp
            b = ( parameters%mp*psn*sfcprs/cs + parameters%bp ) * rlb - 1.
            c = -rlb
            if (b .ge. 0.) then
               q = -0.5*( b + sqrt(b*b-4.*a*c) )
            else
               q = -0.5*( b - sqrt(b*b-4.*a*c) )
            end if
            r1 = q/a
            r2 = c/q
            rs = max(r1,r2)
            ci = max( cs-psn*sfcprs*1.65*rs, 0. )
       end do 

! rs, rb:  s m**2 / umol -> s/m

         rs = rs*cf

  end subroutine stomata

!== begin canres ===================================================================================

!>\ingroup NoahMP_LSM
  subroutine canres (parameters,par   ,sfctmp,rcsoil ,eah   ,sfcprs , & !in
                     rc    ,psn   ,iloc   ,jloc  )           !out

! --------------------------------------------------------------------------------------------------
! calculate canopy resistance which depends on incoming solar radiation,
! air temperature, atmospheric water vapor pressure deficit at the
! lowest model level, and soil moisture (preferably unfrozen soil
! moisture rather than total)
! --------------------------------------------------------------------------------------------------
! source:  jarvis (1976), noilhan and planton (1989, mwr), jacquemin and
! noilhan (1990, blm). chen et al (1996, jgr, vol 101(d3), 7251-7268), 
! eqns 12-14 and table 2 of sec. 3.1.2
! --------------------------------------------------------------------------------------------------
!niu    use module_noahlsm_utility
! --------------------------------------------------------------------------------------------------
    implicit none
! --------------------------------------------------------------------------------------------------
! inputs

  type (noahmp_parameters), intent(in) :: parameters
    integer,                  intent(in)  :: iloc   !grid index
    integer,                  intent(in)  :: jloc   !grid index
    real,                     intent(in)  :: par    !par absorbed per unit sunlit lai (w/m2)
    real,                     intent(in)  :: sfctmp !canopy air temperature
    real,                     intent(in)  :: sfcprs !surface pressure (pa)
    real,                     intent(in)  :: eah    !water vapor pressure (pa)
    real,                     intent(in)  :: rcsoil !soil moisture stress factor

!outputs

    real,                     intent(out) :: rc     !canopy resistance per unit lai
    real,                     intent(out) :: psn    !foliage photosynthesis (umolco2/m2/s)

!local

    real                                  :: rcq
    real                                  :: rcs
    real                                  :: rct
    real                                  :: ff
    real                                  :: q2     !water vapor mixing ratio (kg/kg)
    real                                  :: q2sat  !saturation q2
    real                                  :: dqsdt2 !d(q2sat)/d(t)

! rsmin, rsmax, topt, rgl, hs are canopy stress parameters set in redprm
! ----------------------------------------------------------------------
! initialize canopy resistance multiplier terms.
! ----------------------------------------------------------------------
    rc     = 0.0
    rcs    = 0.0
    rct    = 0.0
    rcq    = 0.0

!  compute q2 and q2sat

    q2 = 0.622 *  eah  / (sfcprs - 0.378 * eah) !specific humidity [kg/kg]
    q2 = q2 / (1.0 + q2)                        !mixing ratio [kg/kg]

    call calhum(parameters,sfctmp, sfcprs, q2sat, dqsdt2)

! contribution due to incoming solar radiation

    ff  = 2.0 * par / parameters%rgl                
    rcs = (ff + parameters%rsmin / parameters%rsmax) / (1.0+ ff)
    rcs = max (rcs,0.0001)

! contribution due to air temperature

    rct = 1.0- 0.0016* ( (parameters%topt - sfctmp)**2.0)
    rct = max (rct,0.0001)

! contribution due to vapor pressure deficit

    rcq = 1.0/ (1.0+ parameters%hs * max(0.,q2sat-q2))
    rcq = max (rcq,0.01)

! determine canopy resistance due to all factors

    rc  = parameters%rsmin / (rcs * rct * rcq * rcsoil)
    psn = -999.99       ! psn not applied for dynamic carbon

  end subroutine canres

!== begin calhum ===================================================================================

!>\ingroup NoahMP_LSM
        subroutine calhum(parameters,sfctmp, sfcprs, q2sat, dqsdt2)

        implicit none

  type (noahmp_parameters), intent(in) :: parameters
        real, intent(in)       :: sfctmp, sfcprs
        real, intent(out)      :: q2sat, dqsdt2
        real, parameter        :: a2=17.67,a3=273.15,a4=29.65, elwv=2.501e6,         &
                                  a23m4=a2*(a3-a4), e0=0.611, rv=461.0,             &
                                  epsilon=0.622
        real                   :: es, sfcprsx

! q2sat: saturated mixing ratio
        es = e0 * exp ( elwv/rv*(1./a3 - 1./sfctmp) )
! convert sfcprs from pa to kpa
        sfcprsx = sfcprs*1.e-3
        q2sat = epsilon * es / (sfcprsx-es)
! convert from  g/g to g/kg
        q2sat = q2sat * 1.e3
! q2sat is currently a 'mixing ratio'

! dqsdt2 is calculated assuming q2sat is a specific humidity
        dqsdt2=(q2sat/(1+q2sat))*a23m4/(sfctmp-a4)**2

! dg q2sat needs to be in g/g when returned for sflx
        q2sat = q2sat / 1.e3

        end subroutine calhum

!== begin tsnosoi ==================================================================================

!>\ingroup NoahMP_LSM
  subroutine tsnosoi (parameters,ice     ,nsoil   ,nsnow   ,isnow   ,ist     , & !in
                      tbot    ,zsnso   ,ssoil   ,df      ,hcpct   , & !in
                      sag     ,dt      ,snowh   ,dzsnso  , & !in
                      tg      ,iloc    ,jloc    ,                   & !in
#ifdef CCPP
                      stc     ,errmsg  ,errflg)                       !inout
#else
                      stc     )                                       !inout
#endif
! --------------------------------------------------------------------------------------------------
! compute snow (up to 3l) and soil (4l) temperature. note that snow temperatures
! during melting season may exceed melting point (tfrz) but later in phasechange
! subroutine the snow temperatures are reset to tfrz for melting snow.
! --------------------------------------------------------------------------------------------------
  implicit none
! --------------------------------------------------------------------------------------------------
!input

  type (noahmp_parameters), intent(in) :: parameters
    integer,                         intent(in)  :: iloc
    integer,                         intent(in)  :: jloc
    integer,                         intent(in)  :: ice    !
    integer,                         intent(in)  :: nsoil  !no of soil layers (4)
    integer,                         intent(in)  :: nsnow  !maximum no of snow layers (3)
    integer,                         intent(in)  :: isnow  !actual no of snow layers
    integer,                         intent(in)  :: ist    !surface type

    real,                            intent(in)  :: dt     !time step (s)
    real,                            intent(in)  :: tbot   !
    real,                            intent(in)  :: ssoil  !ground heat flux (w/m2)
    real,                            intent(in)  :: sag    !solar rad. absorbed by ground (w/m2)
    real,                            intent(in)  :: snowh  !snow depth (m)
    real,                            intent(in)  :: tg     !ground temperature (k)
    real, dimension(-nsnow+1:nsoil), intent(in)  :: zsnso  !layer-bot. depth from snow surf.(m)
    real, dimension(-nsnow+1:nsoil), intent(in)  :: dzsnso !snow/soil layer thickness (m)
    real, dimension(-nsnow+1:nsoil), intent(in)  :: df     !thermal conductivity
    real, dimension(-nsnow+1:nsoil), intent(in)  :: hcpct  !heat capacity (j/m3/k)

!input and output

    real, dimension(-nsnow+1:nsoil), intent(inout) :: stc
#ifdef CCPP
    character(len=*)               , intent(inout) :: errmsg
    integer                        , intent(inout) :: errflg
#endif

!local

    integer                                      :: iz
    real                                         :: zbotsno   !zbot from snow surface
    real, dimension(-nsnow+1:nsoil)              :: ai, bi, ci, rhsts
    real                                         :: eflxb !energy influx from soil bottom (w/m2)
    real, dimension(-nsnow+1:nsoil)              :: phi   !light through water (w/m2)

    real, dimension(-nsnow+1:nsoil) :: tbeg
    real                            :: err_est !heat storage error  (w/m2)
    real                            :: ssoil2  !ground heat flux (w/m2) (for energy check)
    real                            :: eflxb2  !heat flux from the bottom (w/m2) (for energy check)
    character(len=256)              :: message
! ----------------------------------------------------------------------
! compute solar penetration through water, needs more work

    phi(isnow+1:nsoil) = 0.

! adjust zbot from soil surface to zbotsno from snow surface

    zbotsno = parameters%zbot - snowh    !from snow surface

! snow/soil heat storage for energy balance check

    do iz = isnow+1, nsoil
       tbeg(iz) = stc(iz)
    enddo

! compute soil temperatures

      call hrt   (parameters,nsnow     ,nsoil     ,isnow     ,zsnso     , &
                  stc       ,tbot      ,zbotsno   ,dt        , &
                  df        ,hcpct     ,ssoil     ,phi       , &
                  ai        ,bi        ,ci        ,rhsts     , &
                  eflxb     )

      call hstep (parameters,nsnow     ,nsoil     ,isnow     ,dt        , &
                  ai        ,bi        ,ci        ,rhsts     , &
                  stc       ) 

! update ground heat flux just for energy check, but not for final output
! otherwise, it would break the surface energy balance

    if(opt_tbot == 1) then
       eflxb2  = 0.
    else if(opt_tbot == 2) then
       eflxb2  = df(nsoil)*(tbot-stc(nsoil)) / &
            (0.5*(zsnso(nsoil-1)+zsnso(nsoil)) - zbotsno)
    end if

    ! skip the energy balance check for now, until we can make it work
    ! right for small time steps.
    return

! energy balance check

    err_est = 0.0
    do iz = isnow+1, nsoil
       err_est = err_est + (stc(iz)-tbeg(iz)) * dzsnso(iz) * hcpct(iz) / dt
    enddo

    if (opt_stc == 1) then   ! semi-implicit
       err_est = err_est - (ssoil +eflxb)
    else                     ! full-implicit
       ssoil2 = df(isnow+1)*(tg-stc(isnow+1))/(0.5*dzsnso(isnow+1))   !m. barlage
       err_est = err_est - (ssoil2+eflxb2)
    endif

    if (abs(err_est) > 1.) then    ! w/m2
       write(message,*) 'tsnosoi is losing(-)/gaining(+) false energy',err_est,' w/m2'
#ifdef CCPP
       errmsg = trim(message)
#else
       call wrf_message(trim(message))
#endif
       write(message,'(i6,1x,i6,1x,i3,f18.13,5f20.12)') &
            iloc, jloc, ist,err_est,ssoil,snowh,tg,stc(isnow+1),eflxb
#ifdef CCPP
       errmsg = trim(errmsg)//NEW_LINE('A')//trim(message)
#else
       call wrf_message(trim(message))
#endif
       !niu      stop
    end if

  end subroutine tsnosoi

!== begin hrt ======================================================================================

!>\ingroup NoahMP_LSM
  subroutine hrt (parameters,nsnow     ,nsoil     ,isnow     ,zsnso     , &
                  stc       ,tbot      ,zbot      ,dt        , &
                  df        ,hcpct     ,ssoil     ,phi       , &
                  ai        ,bi        ,ci        ,rhsts     , &
                  botflx    )
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! calculate the right hand side of the time tendency term of the soil
! thermal diffusion equation.  also to compute ( prepare ) the matrix
! coefficients for the tri-diagonal matrix of the implicit time scheme.
! ----------------------------------------------------------------------
    implicit none
! ----------------------------------------------------------------------
! input

  type (noahmp_parameters), intent(in) :: parameters
    integer,                         intent(in)  :: nsoil  !no of soil layers (4)
    integer,                         intent(in)  :: nsnow  !maximum no of snow layers (3)
    integer,                         intent(in)  :: isnow  !actual no of snow layers
    real,                            intent(in)  :: tbot   !bottom soil temp. at zbot (k)
    real,                            intent(in)  :: zbot   !depth of lower boundary condition (m)
                                                           !from soil surface not snow surface
    real,                            intent(in)  :: dt     !time step (s)
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
    real, dimension(-nsnow+1:nsoil)              :: dz
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

  end subroutine hrt

!== begin hstep ====================================================================================

!>\ingroup NoahMP_LSM
  subroutine hstep (parameters,nsnow     ,nsoil     ,isnow     ,dt        ,  &
                    ai        ,bi        ,ci        ,rhsts     ,  &
                    stc       )  
! ----------------------------------------------------------------------
! calculate/update the soil temperature field.
! ----------------------------------------------------------------------
    implicit none
! ----------------------------------------------------------------------
! input

  type (noahmp_parameters), intent(in) :: parameters
    integer,                         intent(in)    :: nsoil
    integer,                         intent(in)    :: nsnow
    integer,                         intent(in)    :: isnow
    real,                            intent(in)    :: dt

! output & input
    real, dimension(-nsnow+1:nsoil), intent(inout) :: rhsts
    real, dimension(-nsnow+1:nsoil), intent(inout) :: ai
    real, dimension(-nsnow+1:nsoil), intent(inout) :: bi
    real, dimension(-nsnow+1:nsoil), intent(inout) :: ci
    real, dimension(-nsnow+1:nsoil), intent(inout) :: stc

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


    call rosr12 (ci,ai,bi,ciin,rhstsin,rhsts,isnow+1,nsoil,nsnow)

! update snow & soil temperature

    do k = isnow+1,nsoil
       stc (k) = stc (k) + ci (k)
    end do

  end subroutine hstep

!== begin rosr12 ===================================================================================

!>\ingroup NoahMP_LSM
  subroutine rosr12 (p,a,b,c,d,delta,ntop,nsoil,nsnow)
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
  end subroutine rosr12

!== begin phasechange ==============================================================================

!>\ingroup NoahMP_LSM
  subroutine phasechange (parameters,nsnow   ,nsoil   ,isnow   ,dt      ,fact    , & !in
                          dzsnso  ,hcpct   ,ist     ,iloc    ,jloc    , & !in
                          stc     ,snice   ,snliq   ,sneqv   ,snowh   , & !inout
#ifdef CCPP
                          smc     ,sh2o    ,errmsg  ,errflg  ,          & !inout
#else
                          smc     ,sh2o    ,                            & !inout
#endif
                          qmelt   ,imelt   ,ponding )                     !out
! ----------------------------------------------------------------------
! melting/freezing of snow water and soil water
! ----------------------------------------------------------------------
  implicit none
! ----------------------------------------------------------------------
! inputs

  type (noahmp_parameters), intent(in) :: parameters
  integer, intent(in)                             :: iloc   !grid index
  integer, intent(in)                             :: jloc   !grid index
  integer, intent(in)                             :: nsnow  !maximum no. of snow layers [=3]
  integer, intent(in)                             :: nsoil  !no. of soil layers [=4]
  integer, intent(in)                             :: isnow  !actual no. of snow layers [<=3]
  integer, intent(in)                             :: ist    !surface type: 1->soil; 2->lake
  real, intent(in)                                :: dt     !land model time step (sec)
  real, dimension(-nsnow+1:nsoil), intent(in)     :: fact   !temporary
  real, dimension(-nsnow+1:nsoil), intent(in)     :: dzsnso !snow/soil layer thickness [m]
  real, dimension(-nsnow+1:nsoil), intent(in)     :: hcpct  !heat capacity (j/m3/k)

! outputs
  integer, dimension(-nsnow+1:nsoil), intent(out) :: imelt  !phase change index
  real,                               intent(out) :: qmelt  !snowmelt rate [mm/s]
  real,                               intent(out) :: ponding!snowmelt when snow has no layer [mm]

! inputs and outputs

  real, intent(inout) :: sneqv
  real, intent(inout) :: snowh
  real, dimension(-nsnow+1:nsoil), intent(inout)  :: stc    !snow/soil layer temperature [k]
  real, dimension(       1:nsoil), intent(inout)  :: sh2o   !soil liquid water [m3/m3]
  real, dimension(       1:nsoil), intent(inout)  :: smc    !total soil water [m3/m3]
  real, dimension(-nsnow+1:0)    , intent(inout)  :: snice  !snow layer ice [mm]
  real, dimension(-nsnow+1:0)    , intent(inout)  :: snliq  !snow layer liquid water [mm]
#ifdef CCPP
  character(len=*)               , intent(inout)  :: errmsg
  integer                        , intent(inout)  :: errflg
#endif

! local

  integer                         :: j         !do loop index
  real, dimension(-nsnow+1:nsoil) :: hm        !energy residual [w/m2]
  real, dimension(-nsnow+1:nsoil) :: xm        !melting or freezing water [kg/m2]
  real, dimension(-nsnow+1:nsoil) :: wmass0
  real, dimension(-nsnow+1:nsoil) :: wice0 
  real, dimension(-nsnow+1:nsoil) :: wliq0 
  real, dimension(-nsnow+1:nsoil) :: mice      !soil/snow ice mass [mm]
  real, dimension(-nsnow+1:nsoil) :: mliq      !soil/snow liquid water mass [mm]
  real, dimension(-nsnow+1:nsoil) :: supercool !supercooled water in soil (kg/m2)
  real                            :: heatr     !energy residual or loss after melting/freezing
  real                            :: temp1     !temporary variables [kg/m2]
  real                            :: propor
  real                            :: smp       !frozen water potential (mm)
  real                            :: xmf       !total latent heat of phase change

! ----------------------------------------------------------------------
! initialization

    qmelt   = 0.
    ponding = 0.
    xmf     = 0.

    do j = -nsnow+1, nsoil
         supercool(j) = 0.0
    end do

    do j = isnow+1,0       ! all layers
         mice(j) = snice(j)
         mliq(j) = snliq(j)
    end do

    do j = 1, nsoil               ! soil
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

    if(ist == 1) then
      do j = 1,nsoil
         if (opt_frz == 1) then
            if(stc(j) < tfrz) then
               smp = hfus*(tfrz-stc(j))/(grav*stc(j))             !(m)
               supercool(j) = parameters%smcmax*(smp/parameters%psisat)**(-1./parameters%bexp)
               supercool(j) = supercool(j)*dzsnso(j)*1000.        !(mm)
            end if
         end if
         if (opt_frz == 2) then
#ifdef CCPP
               call frh2o (parameters,supercool(j),stc(j),smc(j),sh2o(j),errmsg,errflg)
               if (errflg /=0) return
#else
               call frh2o (parameters,supercool(j),stc(j),smc(j),sh2o(j))
#endif
               supercool(j) = supercool(j)*dzsnso(j)*1000.        !(mm)
         end if
      enddo
    end if

    do j = isnow+1,nsoil
         if (mice(j) > 0. .and. stc(j) >= tfrz) then  !melting 
             imelt(j) = 1
         endif
         if (mliq(j) > supercool(j) .and. stc(j) < tfrz) then
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
        heatr  = hm(1) - hfus*(temp1-sneqv)/dt  
        if (heatr > 0.) then
              xm(1) = heatr*dt/hfus             
              hm(1) = heatr                    
        else
              xm(1) = 0.
              hm(1) = 0.
        endif
        qmelt   = max(0.,(temp1-sneqv))/dt
        xmf     = hfus*qmelt
        ponding = temp1-sneqv
    endif

! the rate of melting and freezing for snow and soil

    do j = isnow+1,nsoil
      if (imelt(j) > 0 .and. abs(hm(j)) > 0.) then

         heatr = 0.
         if (xm(j) > 0.) then                            
            mice(j) = max(0., wice0(j)-xm(j))
            heatr = hm(j) - hfus*(wice0(j)-mice(j))/dt
         else if (xm(j) < 0.) then                      
            if (j <= 0) then                             ! snow
               mice(j) = min(wmass0(j), wice0(j)-xm(j))  
            else                                         ! soil
               if (wmass0(j) < supercool(j)) then
                  mice(j) = 0.
               else
                  mice(j) = min(wmass0(j) - supercool(j),wice0(j)-xm(j))
                  mice(j) = max(mice(j),0.0)
               endif
            endif
            heatr = hm(j) - hfus*(wice0(j)-mice(j))/dt
         endif

         mliq(j) = max(0.,wmass0(j)-mice(j))

         if (abs(heatr) > 0.) then
            stc(j) = stc(j) + fact(j)*heatr
            if (j <= 0) then                             ! snow
               if (mliq(j)*mice(j)>0.) stc(j) = tfrz
            end if
         endif

         xmf = xmf + hfus * (wice0(j)-mice(j))/dt

         if (j < 1) then
            qmelt = qmelt + max(0.,(wice0(j)-mice(j)))/dt
         endif
      endif
    enddo

    do j = isnow+1,0             ! snow
       snliq(j) = mliq(j)
       snice(j) = mice(j)
    end do

    do j = 1, nsoil              ! soil
       sh2o(j) =  mliq(j)            / (1000. * dzsnso(j))
       smc(j)  = (mliq(j) + mice(j)) / (1000. * dzsnso(j))
    end do
   
  end subroutine phasechange

!== begin frh2o ====================================================================================

!>\ingroup NoahMP_LSM
  subroutine frh2o (parameters,free,tkelv,smc,sh2o,&
#ifdef CCPP
     errmsg,errflg)
#else
     )
#endif

! ----------------------------------------------------------------------
! subroutine frh2o
! ----------------------------------------------------------------------
! calculate amount of supercooled liquid soil water content if
! temperature is below 273.15k (tfrz).  requires newton-type iteration
! to solve the nonlinear implicit equation given in eqn 17 of koren et al
! (1999, jgr, vol 104(d16), 19569-19585).
! ----------------------------------------------------------------------
! new version (june 2001): much faster and more accurate newton
! iteration achieved by first taking log of eqn cited above -- less than
! 4 (typically 1 or 2) iterations achieves convergence.  also, explicit
! 1-step solution option for special case of parameter ck=0, which
! reduces the original implicit equation to a simpler explicit form,
! known as the "flerchinger eqn". improved handling of solution in the
! limit of freezing point temperature tfrz.
! ----------------------------------------------------------------------
! input:

!   tkelv.........temperature (kelvin)
!   smc...........total soil moisture content (volumetric)
!   sh2o..........liquid soil moisture content (volumetric)
!   b.............soil type "b" parameter (from redprm)
!   psisat........saturated soil matric potential (from redprm)

! output:
!   free..........supercooled liquid water content [m3/m3]
! ----------------------------------------------------------------------
    implicit none
  type (noahmp_parameters), intent(in) :: parameters
    real, intent(in)     :: sh2o,smc,tkelv
    real, intent(out)    :: free
#ifdef CCPP
    character(len=*), intent(inout)  :: errmsg
    integer, intent(inout)           :: errflg
#endif
    real                 :: bx,denom,df,dswl,fk,swl,swlk
    integer              :: nlog,kcount
!      parameter(ck = 0.0)
    real, parameter      :: ck = 8.0, blim = 5.5, error = 0.005,       &
         dice = 920.0
    character(len=80)    :: message

! ----------------------------------------------------------------------
! limits on parameter b: b < 5.5  (use parameter blim)
! simulations showed if b > 5.5 unfrozen water content is
! non-realistically high at very low temperatures.
! ----------------------------------------------------------------------
    bx = parameters%bexp
! ----------------------------------------------------------------------
! initializing iterations counter and iterative solution flag.
! ----------------------------------------------------------------------

    if (parameters%bexp >  blim) bx = blim
    nlog = 0

! ----------------------------------------------------------------------
!  if temperature not significantly below freezing (tfrz), sh2o = smc
! ----------------------------------------------------------------------
    kcount = 0
    if (tkelv > (tfrz- 1.e-3)) then
       free = smc
    else

! ----------------------------------------------------------------------
! option 1: iterated solution in koren et al, jgr, 1999, eqn 17
! ----------------------------------------------------------------------
! initial guess for swl (frozen content)
! ----------------------------------------------------------------------
       if (ck /= 0.0) then
          swl = smc - sh2o
! ----------------------------------------------------------------------
! keep within bounds.
! ----------------------------------------------------------------------
          if (swl > (smc -0.02)) swl = smc -0.02
! ----------------------------------------------------------------------
!  start of iterations
! ----------------------------------------------------------------------
          if (swl < 0.) swl = 0.
1001      continue
          if (.not.( (nlog < 10) .and. (kcount == 0)))   goto 1002
          nlog = nlog +1
          df = alog ( ( parameters%psisat * grav / hfus ) * ( ( 1. + ck * swl )**2.) * &
               ( parameters%smcmax / (smc - swl) )** bx) - alog ( - (               &
               tkelv - tfrz)/ tkelv)
          denom = 2. * ck / ( 1. + ck * swl ) + bx / ( smc - swl )
          swlk = swl - df / denom
! ----------------------------------------------------------------------
! bounds useful for mathematical solution.
! ----------------------------------------------------------------------
          if (swlk > (smc -0.02)) swlk = smc - 0.02
          if (swlk < 0.) swlk = 0.

! ----------------------------------------------------------------------
! mathematical solution bounds applied.
! ----------------------------------------------------------------------
          dswl = abs (swlk - swl)
! if more than 10 iterations, use explicit method (ck=0 approx.)
! when dswl less or eq. error, no more iterations required.
! ----------------------------------------------------------------------
          swl = swlk
          if ( dswl <= error ) then
             kcount = kcount +1
          end if
! ----------------------------------------------------------------------
!  end of iterations
! ----------------------------------------------------------------------
! bounds applied within do-block are valid for physical solution.
! ----------------------------------------------------------------------
          goto 1001
1002      continue
          free = smc - swl
       end if
! ----------------------------------------------------------------------
! end option 1
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! option 2: explicit solution for flerchinger eq. i.e. ck=0
! in koren et al., jgr, 1999, eqn 17
! apply physical bounds to flerchinger solution
! ----------------------------------------------------------------------
       if (kcount == 0) then
          write(message, '("flerchinger used in new version. iterations=", i6)') nlog
#ifdef CCPP
          errmsg = trim(message)
#else
          call wrf_message(trim(message))
#endif
          fk = ( ( (hfus / (grav * ( - parameters%psisat)))*                    &
               ( (tkelv - tfrz)/ tkelv))** ( -1/ bx))* parameters%smcmax
          if (fk < 0.02) fk = 0.02
          free = min (fk, smc)
! ----------------------------------------------------------------------
! end option 2
! ----------------------------------------------------------------------
       end if
    end if
! ----------------------------------------------------------------------
  end subroutine frh2o
! ----------------------------------------------------------------------
! ==================================================================================================
! **********************end of energy subroutines***********************
! ==================================================================================================

!== begin water ====================================================================================

!>\ingroup NoahMP_LSM
  subroutine water (parameters,vegtyp ,nsnow  ,nsoil  ,imelt  ,dt     ,uu     , & !in
                    vv     ,fcev   ,fctr   ,qprecc ,qprecl ,elai   , & !in
                    esai   ,sfctmp ,qvap   ,qdew   ,zsoil  ,btrani , & !in
                    ficeold,ponding,tg     ,ist    ,fveg   ,iloc   ,jloc ,smceq , & !in
                    bdfall ,fp     ,rain   ,snow,                    & !in  mb/an: v3.7
		    qsnow  ,qrain  ,snowhin,latheav,latheag,frozen_canopy,frozen_ground,    & !in  mb
                    isnow  ,canliq ,canice ,tv     ,snowh  ,sneqv  , & !inout
                    snice  ,snliq  ,stc    ,zsnso  ,sh2o   ,smc    , & !inout
                    sice   ,zwt    ,wa     ,wt     ,dzsnso ,wslake , & !inout
                    smcwtd ,deeprech,rech                          , & !inout
                    cmc    ,ecan   ,etran  ,fwet   ,runsrf ,runsub , & !out
                    qin    ,qdis   ,ponding1       ,ponding2,        &
                    qsnbot ,esnow) 
! ----------------------------------------------------------------------  
! code history:
! initial code: guo-yue niu, oct. 2007
! ----------------------------------------------------------------------
  implicit none
! ----------------------------------------------------------------------
! input
  type (noahmp_parameters), intent(in) :: parameters
  integer,                         intent(in)    :: iloc    !grid index
  integer,                         intent(in)    :: jloc    !grid index
  integer,                         intent(in)    :: vegtyp  !vegetation type
  integer,                         intent(in)    :: nsnow   !maximum no. of snow layers
  integer                        , intent(in)    :: ist     !surface type 1-soil; 2-lake
  integer,                         intent(in)    :: nsoil   !no. of soil layers
  integer, dimension(-nsnow+1:0) , intent(in)    :: imelt   !melting state index [1-melt; 2-freeze]
  real,                            intent(in)    :: dt      !main time step (s)
  real,                            intent(in)    :: uu      !u-direction wind speed [m/s]
  real,                            intent(in)    :: vv      !v-direction wind speed [m/s]
  real,                            intent(in)    :: fcev    !canopy evaporation (w/m2) [+ to atm ]
  real,                            intent(in)    :: fctr    !transpiration (w/m2) [+ to atm]
  real,                            intent(in)    :: qprecc  !convective precipitation (mm/s)
  real,                            intent(in)    :: qprecl  !large-scale precipitation (mm/s)
  real,                            intent(in)    :: elai    !leaf area index, after burying by snow
  real,                            intent(in)    :: esai    !stem area index, after burying by snow
  real,                            intent(in)    :: sfctmp  !surface air temperature [k]
  real,                            intent(in)    :: qvap    !soil surface evaporation rate[mm/s]
  real,                            intent(in)    :: qdew    !soil surface dew rate[mm/s]
  real, dimension(       1:nsoil), intent(in)    :: zsoil   !depth of layer-bottom from soil surface
  real, dimension(       1:nsoil), intent(in)    :: btrani  !soil water stress factor (0 to 1)
  real, dimension(-nsnow+1:    0), intent(in)    :: ficeold !ice fraction at last timestep
!  real                           , intent(in)    :: ponding ![mm]
  real                           , intent(in)    :: tg      !ground temperature (k)
  real                           , intent(in)    :: fveg    !greeness vegetation fraction (-)
  real                           , intent(in)    :: bdfall   !bulk density of snowfall (kg/m3) ! mb/an: v3.7
  real                           , intent(in)    :: fp       !fraction of the gridcell that receives precipitation ! mb/an: v3.7
  real                           , intent(in)    :: rain     !rainfall (mm/s) ! mb/an: v3.7
  real                           , intent(in)    :: snow     !snowfall (mm/s) ! mb/an: v3.7
  real, dimension(       1:nsoil), intent(in)    :: smceq   !equilibrium soil water content [m3/m3] (used in m-m&f groundwater dynamics)
  real                           , intent(in)    :: qsnow   !snow at ground srf (mm/s) [+]
  real                           , intent(in)    :: qrain   !rain at ground srf (mm) [+]
  real                           , intent(in)    :: snowhin !snow depth increasing rate (m/s)

! input/output
  integer,                         intent(inout) :: isnow   !actual no. of snow layers
  real,                            intent(inout) :: canliq  !intercepted liquid water (mm)
  real,                            intent(inout) :: canice  !intercepted ice mass (mm)
  real,                            intent(inout) :: tv      !vegetation temperature (k)
  real,                            intent(inout) :: snowh   !snow height [m]
  real,                            intent(inout) :: sneqv   !snow water eqv. [mm]
  real, dimension(-nsnow+1:    0), intent(inout) :: snice   !snow layer ice [mm]
  real, dimension(-nsnow+1:    0), intent(inout) :: snliq   !snow layer liquid water [mm]
  real, dimension(-nsnow+1:nsoil), intent(inout) :: stc     !snow/soil layer temperature [k]
  real, dimension(-nsnow+1:nsoil), intent(inout) :: zsnso   !depth of snow/soil layer-bottom
  real, dimension(-nsnow+1:nsoil), intent(inout) :: dzsnso  !snow/soil layer thickness [m]
  real, dimension(       1:nsoil), intent(inout) :: sh2o    !soil liquid water content [m3/m3]
  real, dimension(       1:nsoil), intent(inout) :: sice    !soil ice content [m3/m3]
  real, dimension(       1:nsoil), intent(inout) :: smc     !total soil water content [m3/m3]
  real,                            intent(inout) :: zwt     !the depth to water table [m]
  real,                            intent(inout) :: wa      !water storage in aquifer [mm]
  real,                            intent(inout) :: wt      !water storage in aquifer 
                                                            !+ stuarated soil [mm]
  real,                            intent(inout) :: wslake  !water storage in lake (can be -) (mm)
  real                           , intent(inout) :: ponding ![mm]
  real,                            intent(inout) :: smcwtd !soil water content between bottom of the soil and water table [m3/m3]
  real,                            intent(inout) :: deeprech !recharge to or from the water table when deep [m]
  real,                            intent(inout) :: rech !recharge to or from the water table when shallow [m] (diagnostic)

! output
  real,                            intent(out)   :: cmc     !intercepted water per ground area (mm)
  real,                            intent(out)   :: ecan    !evap of intercepted water (mm/s) [+]
  real,                            intent(out)   :: etran   !transpiration rate (mm/s) [+]
  real,                            intent(out)   :: fwet    !wetted/snowed fraction of canopy (-)
  real,                            intent(out)   :: runsrf  !surface runoff [mm/s] 
  real,                            intent(out)   :: runsub  !baseflow (sturation excess) [mm/s]
  real,                            intent(out)   :: qin     !groundwater recharge [mm/s]
  real,                            intent(out)   :: qdis    !groundwater discharge [mm/s]
  real,                            intent(out)   :: ponding1
  real,                            intent(out)   :: ponding2
  real,                            intent(out)   :: esnow
  real,                            intent(out)   :: qsnbot  !melting water out of snow bottom [mm/s]
  real                              , intent(in)   :: latheav !latent heat vap./sublimation (j/kg)
  real                              , intent(in)   :: latheag !latent heat vap./sublimation (j/kg)
  logical                           , intent(in)   :: frozen_ground ! used to define latent heat pathway
  logical                           , intent(in)   :: frozen_canopy ! used to define latent heat pathway


! local
  integer                                        :: iz
  real                                           :: qinsur  !water input on soil surface [m/s]
  real                                           :: qseva   !soil surface evap rate [mm/s]
  real                                           :: qsdew   !soil surface dew rate [mm/s]
  real                                           :: qsnfro  !snow surface frost rate[mm/s]
  real                                           :: qsnsub  !snow surface sublimation rate [mm/s]
  real, dimension(       1:nsoil)                :: etrani  !transpiration rate (mm/s) [+]
  real, dimension(       1:nsoil)                :: wcnd   !hydraulic conductivity (m/s)
  real                                           :: qdrain  !soil-bottom free drainage [mm/s] 
  real                                           :: snoflow !glacier flow [mm/s]
  real                                           :: fcrmax !maximum of fcr (-)

  real, parameter ::  wslmax = 5000.      !maximum lake water storage (mm)


! ----------------------------------------------------------------------
! initialize

   etrani(1:nsoil) = 0.
   snoflow         = 0.
   runsub          = 0.
   qinsur          = 0.

! canopy-intercepted snowfall/rainfall, drips, and throughfall

   call canwater (parameters,vegtyp ,dt     , & !in
                  fcev   ,fctr   ,elai   , & !in
                  esai   ,tg     ,fveg   ,iloc   , jloc, & !in
                  bdfall ,frozen_canopy  , & !in     
                  canliq ,canice ,tv     ,                 & !inout
                  cmc    ,ecan   ,etran  , & !out
                  fwet      )                           !out

! sublimation, frost, evaporation, and dew

     qsnsub = 0.
     if (sneqv > 0.) then
       qsnsub = min(qvap, sneqv/dt)
     endif
     qseva = qvap-qsnsub
     esnow = qsnsub*2.83e+6

     qsnfro = 0.
     if (sneqv > 0.) then
        qsnfro = qdew
     endif
     qsdew = qdew - qsnfro

     call snowwater (parameters,nsnow  ,nsoil  ,imelt  ,dt     ,zsoil  , & !in
          &          sfctmp ,snowhin,qsnow  ,qsnfro ,qsnsub , & !in
          &          qrain  ,ficeold,iloc   ,jloc   ,         & !in
          &          isnow  ,snowh  ,sneqv  ,snice  ,snliq  , & !inout
          &          sh2o   ,sice   ,stc    ,zsnso  ,dzsnso , & !inout
          &          qsnbot ,snoflow,ponding1       ,ponding2)  !out

   if(frozen_ground) then
      sice(1) =  sice(1) + (qsdew-qseva)*dt/(dzsnso(1)*1000.)
      qsdew = 0.0
      qseva = 0.0
      if(sice(1) < 0.) then
         sh2o(1) = sh2o(1) + sice(1)
         sice(1) = 0.
      end if
   end if

! convert units (mm/s -> m/s)

    !ponding: melting water from snow when there is no layer
    qinsur = (ponding+ponding1+ponding2)/dt * 0.001
!    qinsur = ponding/dt * 0.001

    if(isnow == 0) then
       qinsur = qinsur+(qsnbot + qsdew + qrain) * 0.001
    else
       qinsur = qinsur+(qsnbot + qsdew) * 0.001
    endif

    qseva  = qseva * 0.001 

    do iz = 1, parameters%nroot
       etrani(iz) = etran * btrani(iz) * 0.001
    enddo


! lake/soil water balances

    if (ist == 2) then                                        ! lake
       runsrf = 0.
       if(wslake >= wslmax) runsrf = qinsur*1000.             !mm/s
       wslake = wslake + (qinsur-qseva)*1000.*dt -runsrf*dt   !mm
    else                                                      ! soil
       call      soilwater (parameters,nsoil  ,nsnow  ,dt     ,zsoil  ,dzsnso , & !in
                            qinsur ,qseva  ,etrani ,sice   ,iloc   , jloc , & !in
                            sh2o   ,smc    ,zwt    ,vegtyp , & !inout
                           smcwtd, deeprech                       , & !inout
                            runsrf ,qdrain ,runsub ,wcnd   ,fcrmax )   !out
 
       if(opt_run == 1) then 
          call groundwater (parameters,nsnow  ,nsoil  ,dt     ,sice   ,zsoil  , & !in
                            stc    ,wcnd   ,fcrmax ,iloc   ,jloc   , & !in
                            sh2o   ,zwt    ,wa     ,wt     ,         & !inout
                            qin    ,qdis   )                           !out
          runsub       = qdis          !mm/s
       end if

       if(opt_run == 3 .or. opt_run == 4) then 
          runsub       = runsub + qdrain        !mm/s
       end if

       do iz = 1,nsoil
           smc(iz) = sh2o(iz) + sice(iz)
       enddo
 
       if(opt_run == 5) then
          call shallowwatertable (parameters,nsnow  ,nsoil, zsoil, dt       , & !in
                         dzsnso ,smceq   ,iloc , jloc        , & !in
                         smc    ,zwt    ,smcwtd ,rech, qdrain  ) !inout

          sh2o(nsoil) = smc(nsoil) - sice(nsoil)
          runsub = runsub + qdrain !it really comes from subroutine watertable, which is not called with the same frequency as the soil routines here
          wa = 0.
       endif

    endif

    runsub       = runsub + snoflow         !mm/s

  end subroutine water

!== begin canwater =================================================================================

!>\ingroup NoahMP_LSM
  subroutine canwater (parameters,vegtyp ,dt     , & !in
                       fcev   ,fctr   ,elai   , & !in
                       esai   ,tg     ,fveg   ,iloc   , jloc , & !in
                       bdfall ,frozen_canopy  ,  & !in      
                       canliq ,canice ,tv     ,                 & !inout
                       cmc    ,ecan   ,etran  , & !out
                       fwet      )                           !out

! ------------------------ code history ------------------------------
! canopy hydrology
! --------------------------------------------------------------------
  implicit none
! ------------------------ input/output variables --------------------
! input
  type (noahmp_parameters), intent(in) :: parameters
  integer,intent(in)  :: iloc    !grid index
  integer,intent(in)  :: jloc    !grid index
  integer,intent(in)  :: vegtyp  !vegetation type
  real,   intent(in)  :: dt      !main time step (s)
  real,   intent(in)  :: fcev    !canopy evaporation (w/m2) [+ = to atm]
  real,   intent(in)  :: fctr    !transpiration (w/m2) [+ = to atm]
  real,   intent(in)  :: elai    !leaf area index, after burying by snow
  real,   intent(in)  :: esai    !stem area index, after burying by snow
  real,   intent(in)  :: tg      !ground temperature (k)
  real,   intent(in)  :: fveg    !greeness vegetation fraction (-)
  logical                           , intent(in)   :: frozen_canopy ! used to define latent heat pathway
  real                           , intent(in)    :: bdfall   !bulk density of snowfall (kg/m3) ! mb/an: v3.7

! input & output
  real, intent(inout) :: canliq  !intercepted liquid water (mm)
  real, intent(inout) :: canice  !intercepted ice mass (mm)
  real, intent(inout) :: tv      !vegetation temperature (k)

! output
  real, intent(out)   :: cmc     !intercepted water (mm)
  real, intent(out)   :: ecan    !evaporation of intercepted water (mm/s) [+]
  real, intent(out)   :: etran   !transpiration rate (mm/s) [+]
  real, intent(out)   :: fwet    !wetted or snowed fraction of the canopy (-)
! --------------------------------------------------------------------

! ------------------------ local variables ---------------------------
  real                :: maxsno  !canopy capacity for snow interception (mm)
  real                :: maxliq  !canopy capacity for rain interception (mm)
  real                :: qevac   !evaporation rate (mm/s)
  real                :: qdewc   !dew rate (mm/s)
  real                :: qfroc   !frost rate (mm/s)
  real                :: qsubc   !sublimation rate (mm/s)
  real                :: qmeltc  !melting rate of canopy snow (mm/s)
  real                :: qfrzc   !refreezing rate of canopy liquid water (mm/s)
  real                :: canmas  !total canopy mass (kg/m2)
! --------------------------------------------------------------------
! initialization

      ecan    = 0.0

! --------------------------- liquid water ------------------------------
! maximum canopy water

      maxliq =  parameters%ch2op * (elai+ esai)

! evaporation, transpiration, and dew

      if (.not.frozen_canopy) then             ! barlage: change to frozen_canopy
        etran = max( fctr/hvap, 0. )
        qevac = max( fcev/hvap, 0. )
        qdewc = abs( min( fcev/hvap, 0. ) )
        qsubc = 0.
        qfroc = 0.
      else
        etran = max( fctr/hsub, 0. )
        qevac = 0.
        qdewc = 0.
        qsubc = max( fcev/hsub, 0. )
        qfroc = abs( min( fcev/hsub, 0. ) )
      endif

! canopy water balance. for convenience allow dew to bring canliq above
! maxh2o or else would have to re-adjust drip

       qevac = min(canliq/dt,qevac)
       canliq=max(0.,canliq+(qdewc-qevac)*dt)
       if(canliq <= 1.e-06) canliq = 0.0

! --------------------------- canopy ice ------------------------------
! for canopy ice

      maxsno = 6.6*(0.27+46./bdfall) * (elai+ esai)

      qsubc = min(canice/dt,qsubc) 
      canice= max(0.,canice + (qfroc-qsubc)*dt)
      if(canice.le.1.e-6) canice = 0.
     
! wetted fraction of canopy

      if(canice.gt.0.) then
           fwet = max(0.,canice) / max(maxsno,1.e-06)
      else
           fwet = max(0.,canliq) / max(maxliq,1.e-06)
      endif
      fwet = min(fwet, 1.) ** 0.667

! phase change

      qmeltc = 0.
      qfrzc = 0.

      if(canice.gt.1.e-6.and.tv.gt.tfrz) then
         qmeltc = min(canice/dt,(tv-tfrz)*cice*canice/denice/(dt*hfus))
         canice = max(0.,canice - qmeltc*dt)
         canliq = max(0.,canliq + qmeltc*dt)
         tv     = fwet*tfrz + (1.-fwet)*tv
      endif

      if(canliq.gt.1.e-6.and.tv.lt.tfrz) then
         qfrzc  = min(canliq/dt,(tfrz-tv)*cwat*canliq/denh2o/(dt*hfus))
         canliq = max(0.,canliq - qfrzc*dt)
         canice = max(0.,canice + qfrzc*dt)
         tv     = fwet*tfrz + (1.-fwet)*tv
      endif

! total canopy water

      cmc = canliq + canice

! total canopy evaporation

      ecan = qevac + qsubc - qdewc - qfroc

  end subroutine canwater

!== begin snowwater ================================================================================

!>\ingroup NoahMP_LSM
  subroutine snowwater (parameters,nsnow  ,nsoil  ,imelt  ,dt     ,zsoil  , & !in
                        sfctmp ,snowhin,qsnow  ,qsnfro ,qsnsub , & !in
                        qrain  ,ficeold,iloc   ,jloc   ,         & !in
                        isnow  ,snowh  ,sneqv  ,snice  ,snliq  , & !inout
                        sh2o   ,sice   ,stc    ,zsnso  ,dzsnso , & !inout
                        qsnbot ,snoflow,ponding1       ,ponding2)  !out
! ----------------------------------------------------------------------
  implicit none
! ----------------------------------------------------------------------
! input
  type (noahmp_parameters), intent(in) :: parameters
  integer,                         intent(in)    :: iloc   !grid index
  integer,                         intent(in)    :: jloc   !grid index
  integer,                         intent(in)    :: nsnow  !maximum no. of snow layers
  integer,                         intent(in)    :: nsoil  !no. of soil layers
  integer, dimension(-nsnow+1:0) , intent(in)    :: imelt  !melting state index [0-no melt;1-melt]
  real,                            intent(in)    :: dt     !time step (s)
  real, dimension(       1:nsoil), intent(in)    :: zsoil  !depth of layer-bottom from soil surface
  real,                            intent(in)    :: sfctmp !surface air temperature [k]
  real,                            intent(in)    :: snowhin!snow depth increasing rate (m/s)
  real,                            intent(in)    :: qsnow  !snow at ground srf (mm/s) [+]
  real,                            intent(in)    :: qsnfro !snow surface frost rate[mm/s]
  real,                            intent(in)    :: qsnsub !snow surface sublimation rate[mm/s]
  real,                            intent(in)    :: qrain  !snow surface rain rate[mm/s]
  real, dimension(-nsnow+1:0)    , intent(in)    :: ficeold!ice fraction at last timestep

! input & output
  integer,                         intent(inout) :: isnow  !actual no. of snow layers
  real,                            intent(inout) :: snowh  !snow height [m]
  real,                            intent(inout) :: sneqv  !snow water eqv. [mm]
  real, dimension(-nsnow+1:    0), intent(inout) :: snice  !snow layer ice [mm]
  real, dimension(-nsnow+1:    0), intent(inout) :: snliq  !snow layer liquid water [mm]
  real, dimension(       1:nsoil), intent(inout) :: sh2o   !soil liquid moisture (m3/m3)
  real, dimension(       1:nsoil), intent(inout) :: sice   !soil ice moisture (m3/m3)
  real, dimension(-nsnow+1:nsoil), intent(inout) :: stc    !snow layer temperature [k]
  real, dimension(-nsnow+1:nsoil), intent(inout) :: zsnso  !depth of snow/soil layer-bottom
  real, dimension(-nsnow+1:nsoil), intent(inout) :: dzsnso !snow/soil layer thickness [m]

! output
  real,                              intent(out) :: qsnbot !melting water out of snow bottom [mm/s]
  real,                              intent(out) :: snoflow!glacier flow [mm]
  real,                              intent(out) :: ponding1
  real,                              intent(out) :: ponding2

! local
  integer :: iz,i
  real    :: bdsnow  !bulk density of snow (kg/m3)
! ----------------------------------------------------------------------
   snoflow = 0.0
   ponding1 = 0.0
   ponding2 = 0.0

   call snowfall (parameters,nsoil  ,nsnow  ,dt     ,qsnow  ,snowhin, & !in
                  sfctmp ,iloc   ,jloc   ,                 & !in
                  isnow  ,snowh  ,dzsnso ,stc    ,snice  , & !inout
                  snliq  ,sneqv  )                           !inout

! mb: do each if block separately

   if(isnow < 0) &        ! when multi-layer
   call  compact (parameters,nsnow  ,nsoil  ,dt     ,stc    ,snice  , & !in
                  snliq  ,zsoil  ,imelt  ,ficeold,iloc   , jloc ,& !in
                  isnow  ,dzsnso ,zsnso  )                   !inout

   if(isnow < 0) &        !when multi-layer
   call  combine (parameters,nsnow  ,nsoil  ,iloc   ,jloc   ,         & !in
                  isnow  ,sh2o   ,stc    ,snice  ,snliq  , & !inout
                  dzsnso ,sice   ,snowh  ,sneqv  ,         & !inout
                  ponding1       ,ponding2)                  !out

   if(isnow < 0) &        !when multi-layer
   call   divide (parameters,nsnow  ,nsoil  ,                         & !in
                  isnow  ,stc    ,snice  ,snliq  ,dzsnso )   !inout

   call  snowh2o (parameters,nsnow  ,nsoil  ,dt     ,qsnfro ,qsnsub , & !in 
                  qrain  ,iloc   ,jloc   ,                 & !in
                  isnow  ,dzsnso ,snowh  ,sneqv  ,snice  , & !inout
                  snliq  ,sh2o   ,sice   ,stc    ,         & !inout
                  qsnbot ,ponding1       ,ponding2)           !out

!set empty snow layers to zero

   do iz = -nsnow+1, isnow
        snice(iz) = 0.
        snliq(iz) = 0.
        stc(iz)   = 0.
        dzsnso(iz)= 0.
        zsnso(iz) = 0.
   enddo

!to obtain equilibrium state of snow in glacier region
       
   if(sneqv > 2000.) then   ! 2000 mm -> maximum water depth
      bdsnow      = snice(0) / dzsnso(0)
      snoflow     = (sneqv - 2000.)
      snice(0)    = snice(0)  - snoflow 
      dzsnso(0)   = dzsnso(0) - snoflow/bdsnow
      snoflow     = snoflow / dt
   end if

! sum up snow mass for layered snow

   if(isnow < 0) then  ! mb: only do for multi-layer
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

  end subroutine snowwater

!== begin snowfall =================================================================================

!>\ingroup NoahMP_LSM
  subroutine snowfall (parameters,nsoil  ,nsnow  ,dt     ,qsnow  ,snowhin , & !in
                       sfctmp ,iloc   ,jloc   ,                  & !in
                       isnow  ,snowh  ,dzsnso ,stc    ,snice   , & !inout
                       snliq  ,sneqv  )                            !inout
! ----------------------------------------------------------------------
! snow depth and density to account for the new snowfall.
! new values of snow depth & density returned.
! ----------------------------------------------------------------------
    implicit none
! ----------------------------------------------------------------------
! input

  type (noahmp_parameters), intent(in) :: parameters
  integer,                            intent(in) :: iloc   !grid index
  integer,                            intent(in) :: jloc   !grid index
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
 
    if(isnow == 0  .and. qsnow>0. .and. snowh >= 0.025) then !mb: change limit
!    if(isnow == 0  .and. qsnow>0. .and. snowh >= 0.05) then
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
  end subroutine snowfall

!== begin combine ==================================================================================

!>\ingroup NoahMP_LSM
  subroutine combine (parameters,nsnow  ,nsoil  ,iloc   ,jloc   ,         & !in
                      isnow  ,sh2o   ,stc    ,snice  ,snliq  , & !inout
                      dzsnso ,sice   ,snowh  ,sneqv  ,         & !inout
                      ponding1       ,ponding2)                  !out
! ----------------------------------------------------------------------
    implicit none
! ----------------------------------------------------------------------
! input

  type (noahmp_parameters), intent(in) :: parameters
    integer, intent(in)     :: iloc
    integer, intent(in)     :: jloc
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
    real,                            intent(out) :: ponding1
    real,                            intent(out) :: ponding2

! local variables:

    integer :: i,j,k,l               ! node indices
    integer :: isnow_old             ! number of top snow layer
    integer :: mssi                  ! node index
    integer :: neibor                ! adjacent node selected for combination
    real    :: zwice                 ! total ice mass in snow
    real    :: zwliq                 ! total liquid water in snow

    real    :: dzmin(3)              ! minimum of top snow layer
!    data dzmin /0.045, 0.05, 0.2/
    data dzmin /0.025, 0.025, 0.1/  ! mb: change limit
!-----------------------------------------------------------------------

       isnow_old = isnow

       do j = isnow_old+1,0
          if (snice(j) <= .1) then
             if(j /= 0) then
                snliq(j+1) = snliq(j+1) + snliq(j)
                snice(j+1) = snice(j+1) + snice(j)
             else
               if (isnow_old < -1) then    ! mb/km: change to isnow
                snliq(j-1) = snliq(j-1) + snliq(j)
                snice(j-1) = snice(j-1) + snice(j)
               else
	         if(snice(j) >= 0.) then
                  ponding1 = snliq(j)    ! isnow will get set to zero below; ponding1 will get 
                  sneqv = snice(j)       ! added to ponding from phasechange ponding should be
                  snowh = dzsnso(j)      ! zero here because it was calculated for thin snow
		 else   ! snice over-sublimated earlier
		  ponding1 = snliq(j) + snice(j)
		  if(ponding1 < 0.) then  ! if snice and snliq sublimates remove from soil
		   sice(1) = max(0.0,sice(1)+ponding1/(dzsnso(1)*1000.))
                   ponding1 = 0.0
		  end if
                  sneqv = 0.0
                  snowh = 0.0
		 end if
                 snliq(j) = 0.0
                 snice(j) = 0.0
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

       if (snowh < 0.025 .and. isnow < 0 ) then ! mb: change limit
!       if (snowh < 0.05 .and. isnow < 0 ) then
          isnow  = 0
          sneqv = zwice
          ponding2 = zwliq           ! limit of isnow < 0 means input ponding
          if(sneqv <= 0.) snowh = 0. ! should be zero; see above
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

                call combo (parameters,dzsnso(j), snliq(j), snice(j), &
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

  end subroutine combine

!== begin divide ===================================================================================

!>\ingroup NoahMP_LSM
  subroutine divide (parameters,nsnow  ,nsoil  ,                         & !in
                     isnow  ,stc    ,snice  ,snliq  ,dzsnso  )  !inout
! ----------------------------------------------------------------------
    implicit none
! ----------------------------------------------------------------------
! input

  type (noahmp_parameters), intent(in) :: parameters
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

             call combo (parameters,dz(2), swliq(2), swice(2), tsno(2), drr, &
                  zwliq, zwice, tsno(1))

             ! subdivide a new layer
             if (msno <= 2 .and. dz(2) > 0.20) then  ! mb: change limit
!             if (msno <= 2 .and. dz(2) > 0.10) then
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
             call combo (parameters,dz(3), swliq(3), swice(3), tsno(3), drr, &
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

  end subroutine divide

!== begin combo ====================================================================================

!>\ingroup NoahMP_LSM
  subroutine combo(parameters,dz,  wliq,  wice, t, dz2, wliq2, wice2, t2)
! ----------------------------------------------------------------------
    implicit none
! ----------------------------------------------------------------------

! ----------------------------------------------------------------------s
! input

  type (noahmp_parameters), intent(in) :: parameters
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

  end subroutine combo

!== begin compact ==================================================================================

!>\ingroup NoahMP_LSM
  subroutine compact (parameters,nsnow  ,nsoil  ,dt     ,stc    ,snice  , & !in
                      snliq  ,zsoil  ,imelt  ,ficeold,iloc   , jloc , & !in
                      isnow  ,dzsnso ,zsnso )                    !inout
! ----------------------------------------------------------------------
  implicit none
! ----------------------------------------------------------------------
! input
  type (noahmp_parameters), intent(in) :: parameters
   integer,                         intent(in)    :: iloc   !grid index
   integer,                         intent(in)    :: jloc   !grid index
   integer,                         intent(in)    :: nsoil  !no. of soil layers [ =4]
   integer,                         intent(in)    :: nsnow  !maximum no. of snow layers [ =3]
   integer, dimension(-nsnow+1:0) , intent(in)    :: imelt  !melting state index [0-no melt;1-melt]
   real,                            intent(in)    :: dt     !time step (sec)
   real, dimension(-nsnow+1:nsoil), intent(in)    :: stc    !snow layer temperature [k]
   real, dimension(-nsnow+1:    0), intent(in)    :: snice  !snow layer ice [mm]
   real, dimension(-nsnow+1:    0), intent(in)    :: snliq  !snow layer liquid water [mm]
   real, dimension(       1:nsoil), intent(in)    :: zsoil  !depth of layer-bottom from soil srf
   real, dimension(-nsnow+1:    0), intent(in)    :: ficeold!ice fraction at last timestep

! input and output
   integer,                         intent(inout) :: isnow  ! actual no. of snow layers
   real, dimension(-nsnow+1:nsoil), intent(inout) :: dzsnso ! snow layer thickness [m]
   real, dimension(-nsnow+1:nsoil), intent(inout) :: zsnso  ! depth of snow/soil layer-bottom

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

  end subroutine compact

!== begin snowh2o ==================================================================================

!>\ingroup NoahMP_LSM
  subroutine snowh2o (parameters,nsnow  ,nsoil  ,dt     ,qsnfro ,qsnsub , & !in 
                      qrain  ,iloc   ,jloc   ,                 & !in
                      isnow  ,dzsnso ,snowh  ,sneqv  ,snice  , & !inout
                      snliq  ,sh2o   ,sice   ,stc    ,         & !inout
                      qsnbot ,ponding1       ,ponding2)          !out
! ----------------------------------------------------------------------
! renew the mass of ice lens (snice) and liquid (snliq) of the
! surface snow layer resulting from sublimation (frost) / evaporation (dew)
! ----------------------------------------------------------------------
   implicit none
! ----------------------------------------------------------------------
! input

  type (noahmp_parameters), intent(in) :: parameters
   integer,                         intent(in)    :: iloc   !grid index
   integer,                         intent(in)    :: jloc   !grid index
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

! local variables:

   integer                     :: j         !do loop/array indices
   real                        :: qin       !water flow into the element (mm/s)
   real                        :: qout      !water flow out of the element (mm/s)
   real                        :: wgdif     !ice mass after minus sublimation
   real, dimension(-nsnow+1:0) :: vol_liq   !partial volume of liquid water in layer
   real, dimension(-nsnow+1:0) :: vol_ice   !partial volume of ice lens in layer
   real, dimension(-nsnow+1:0) :: epore     !effective porosity = porosity - vol_ice
   real :: propor, temp
   real :: ponding1, ponding2
! ----------------------------------------------------------------------

!for the case when sneqv becomes '0' after 'combine'

   if(sneqv == 0.) then
      sice(1) =  sice(1) + (qsnfro-qsnsub)*dt/(dzsnso(1)*1000.)  ! barlage: sh2o->sice v3.6
      if(sice(1) < 0.) then
         sh2o(1) = sh2o(1) + sice(1)
         sice(1) = 0.
      end if
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
         call  combine (parameters,nsnow  ,nsoil  ,iloc, jloc   , & !in
              isnow  ,sh2o   ,stc    ,snice  ,snliq  , & !inout
              dzsnso ,sice   ,snowh  ,sneqv  ,         & !inout
              ponding1, ponding2 )                       !out
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
               qout = max(0.,(vol_liq(j)-parameters%ssi*epore(j))*dzsnso(j))
               qout = min(qout,(1.-vol_ice(j+1)-vol_liq(j+1))*dzsnso(j+1))
            end if
         else
            qout = max(0.,(vol_liq(j) - parameters%ssi*epore(j))*dzsnso(j))
         end if
         qout = qout*1000.
         snliq(j) = snliq(j) - qout
         qin = qout
      end if
   end do

! liquid water from snow bottom to soil

   qsnbot = qout / dt           ! mm/s

  end subroutine snowh2o

!== begin soilwater ================================================================================

!>\ingroup NoahMP_LSM
  subroutine soilwater (parameters,nsoil  ,nsnow  ,dt     ,zsoil  ,dzsnso , & !in
                        qinsur ,qseva  ,etrani ,sice   ,iloc   , jloc, & !in
                        sh2o   ,smc    ,zwt    ,vegtyp ,& !inout
                        smcwtd, deeprech                       ,& !inout
                        runsrf ,qdrain ,runsub ,wcnd   ,fcrmax )   !out

! ----------------------------------------------------------------------
! calculate surface runoff and soil moisture.
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
  implicit none
! ----------------------------------------------------------------------
! input
  type (noahmp_parameters), intent(in) :: parameters
  integer,                     intent(in) :: iloc   !grid index
  integer,                     intent(in) :: jloc   !grid index
  integer,                     intent(in) :: nsoil  !no. of soil layers
  integer,                     intent(in) :: nsnow  !maximum no. of snow layers
  real,                        intent(in) :: dt     !time step (sec)
  real, intent(in)                        :: qinsur !water input on soil surface [mm/s]
  real, intent(in)                        :: qseva  !evap from soil surface [mm/s]
  real, dimension(1:nsoil),    intent(in) :: zsoil  !depth of soil layer-bottom [m]
  real, dimension(1:nsoil),    intent(in) :: etrani !evapotranspiration from soil layers [mm/s]
  real, dimension(-nsnow+1:nsoil), intent(in) :: dzsnso !snow/soil layer depth [m]
  real, dimension(1:nsoil), intent(in)   :: sice   !soil ice content [m3/m3]

  integer,                     intent(in) :: vegtyp

! input & output
  real, dimension(1:nsoil), intent(inout) :: sh2o   !soil liquid water content [m3/m3]
  real, dimension(1:nsoil), intent(inout) :: smc    !total soil water content [m3/m3]
  real, intent(inout)                     :: zwt    !water table depth [m]
  real,                     intent(inout) :: smcwtd !soil moisture between bottom of the soil and the water table [m3/m3]
  real                    , intent(inout) :: deeprech

! output
  real, intent(out)                       :: qdrain !soil-bottom free drainage [mm/s] 
  real, intent(out)                       :: runsrf !surface runoff [mm/s] 
  real, intent(out)                       :: runsub !subsurface runoff [mm/s] 
  real, intent(out)                       :: fcrmax !maximum of fcr (-)
  real, dimension(1:nsoil), intent(out)   :: wcnd   !hydraulic conductivity (m/s)

! local
  integer                                 :: k,iz   !do-loop index
  integer                                 :: iter   !iteration index
  real                                    :: dtfine !fine time step (s)
  real, dimension(1:nsoil)                :: rhstt  !right-hand side term of the matrix
  real, dimension(1:nsoil)                :: ai     !left-hand side term
  real, dimension(1:nsoil)                :: bi     !left-hand side term
  real, dimension(1:nsoil)                :: ci     !left-hand side term

  real                                    :: fff    !runoff decay factor (m-1)
  real                                    :: rsbmx  !baseflow coefficient [mm/s]
  real                                    :: pddum  !infiltration rate at surface (m/s)
  real                                    :: fice   !ice fraction in frozen soil
  real                                    :: wplus  !saturation excess of the total soil [m]
  real                                    :: rsat   !accumulation of wplus (saturation excess) [m]
  real                                    :: sicemax!maximum soil ice content (m3/m3)
  real                                    :: sh2omin!minimum soil liquid water content (m3/m3)
  real                                    :: wtsub  !sum of wcnd(k)*dzsnso(k)
  real                                    :: mh2o   !water mass removal (mm)
  real                                    :: fsat   !fractional saturated area (-)
  real, dimension(1:nsoil)                :: mliq   !
  real                                    :: xs     !
  real                                    :: watmin !
  real                                    :: qdrain_save !
  real                                    :: epore  !effective porosity [m3/m3]
  real, dimension(1:nsoil)                :: fcr    !impermeable fraction due to frozen soil
  integer                                 :: niter  !iteration times soil moisture (-)
  real                                    :: smctot !2-m averaged soil moisture (m3/m3)
  real                                    :: dztot  !2-m soil depth (m)
  real, parameter :: a = 4.0
! ----------------------------------------------------------------------
    runsrf = 0.0
    pddum  = 0.0
    rsat   = 0.0

! for the case when snowmelt water is too large

    do k = 1,nsoil
       epore   = max ( 1.e-4 , ( parameters%smcmax - sice(k) ) )
       rsat    = rsat + max(0.,sh2o(k)-epore)*dzsnso(k)  
       sh2o(k) = min(epore,sh2o(k))             
    end do

!impermeable fraction due to frozen soil

    do k = 1,nsoil
       fice    = min(1.0,sice(k)/parameters%smcmax)
       fcr(k)  = max(0.0,exp(-a*(1.-fice))- exp(-a)) /  &
                        (1.0              - exp(-a))
    end do

! maximum soil ice content and minimum liquid water of all layers

    sicemax = 0.0
    fcrmax  = 0.0
    sh2omin = parameters%smcmax
    do k = 1,nsoil
       if (sice(k) > sicemax) sicemax = sice(k)
       if (fcr(k)  > fcrmax)  fcrmax  = fcr(k)
       if (sh2o(k) < sh2omin) sh2omin = sh2o(k)
    end do

!subsurface runoff for runoff scheme option 2

    if(opt_run == 2) then 
        fff   = 2.0
        rsbmx = 4.0
        call zwteq (parameters,nsoil  ,nsnow  ,zsoil  ,dzsnso ,sh2o   ,zwt)
        runsub = (1.0-fcrmax) * rsbmx * exp(-parameters%timean) * exp(-fff*zwt)   ! mm/s
    end if

!surface runoff and infiltration rate using different schemes

!jref impermable surface at urban
    if ( parameters%urban_flag ) fcr(1)= 0.95

    if(opt_run == 1) then
       fff = 6.0
       fsat   = parameters%fsatmx*exp(-0.5*fff*(zwt-2.0))
       if(qinsur > 0.) then
         runsrf = qinsur * ( (1.0-fcr(1))*fsat + fcr(1) )
         pddum  = qinsur - runsrf                          ! m/s 
       end if
    end if

    if(opt_run == 5) then
       fff = 6.0
       fsat   = parameters%fsatmx*exp(-0.5*fff*max(-2.0-zwt,0.))
       if(qinsur > 0.) then
         runsrf = qinsur * ( (1.0-fcr(1))*fsat + fcr(1) )
         pddum  = qinsur - runsrf                          ! m/s
       end if
    end if

    if(opt_run == 2) then
       fff   = 2.0
       fsat   = parameters%fsatmx*exp(-0.5*fff*zwt)
       if(qinsur > 0.) then
         runsrf = qinsur * ( (1.0-fcr(1))*fsat + fcr(1) )
         pddum  = qinsur - runsrf                          ! m/s 
       end if
    end if

    if(opt_run == 3) then
       call infil (parameters,nsoil  ,dt     ,zsoil  ,sh2o   ,sice   , & !in
                   sicemax,qinsur ,                         & !in
                   pddum  ,runsrf )                           !out
    end if

    if(opt_run == 4) then
       smctot = 0.
       dztot  = 0.
       do k = 1,nsoil
          dztot   = dztot  + dzsnso(k)  
          smctot  = smctot + smc(k)*dzsnso(k)
          if(dztot >= 2.0) exit
       end do
       smctot = smctot/dztot
       fsat   = max(0.01,smctot/parameters%smcmax) ** 4.        !bats

       if(qinsur > 0.) then
         runsrf = qinsur * ((1.0-fcr(1))*fsat+fcr(1))  
         pddum  = qinsur - runsrf                       ! m/s
       end if
    end if

! determine iteration times and finer time step

    niter = 1

    if(opt_inf == 1) then    !opt_inf =2 may cause water imbalance
       niter = 3
       if (pddum*dt>dzsnso(1)*parameters%smcmax ) then
          niter = niter*2
       end if
    end if                 

    dtfine  = dt / niter

! solve soil moisture

    qdrain_save = 0.0
    do iter = 1, niter
       call srt   (parameters,nsoil  ,zsoil  ,dtfine ,pddum  ,etrani , & !in
                   qseva  ,sh2o   ,smc    ,zwt    ,fcr    , & !in
                   sicemax,fcrmax ,iloc   ,jloc   ,smcwtd ,         & !in
                   rhstt  ,ai     ,bi     ,ci     ,qdrain , & !out
                   wcnd   )                                   !out
  
       call sstep (parameters,nsoil  ,nsnow  ,dtfine ,zsoil  ,dzsnso , & !in
                   sice   ,iloc   ,jloc   ,zwt            ,                 & !in
                   sh2o   ,smc    ,ai     ,bi     ,ci     , & !inout
                   rhstt  ,smcwtd ,qdrain ,deeprech,                                 & !inout
                   wplus)                                     !out
       rsat =  rsat + wplus
       qdrain_save = qdrain_save + qdrain
    end do

    qdrain = qdrain_save/niter

    runsrf = runsrf * 1000. + rsat * 1000./dt  ! m/s -> mm/s
    qdrain = qdrain * 1000.

!wrf_hydro_djg...
!yw    infxsrt = runsrf * dt   !mm/s -> mm

! removal of soil water due to groundwater flow (option 2)

    if(opt_run == 2) then
         wtsub = 0.
         do k = 1, nsoil
           wtsub = wtsub + wcnd(k)*dzsnso(k)
         end do

         do k = 1, nsoil
           mh2o    = runsub*dt*(wcnd(k)*dzsnso(k))/wtsub       ! mm
           sh2o(k) = sh2o(k) - mh2o/(dzsnso(k)*1000.)
         end do
    end if

! limit mliq to be greater than or equal to watmin.
! get water needed to bring mliq equal watmin from lower layer.

   if(opt_run /= 1) then
      do iz = 1, nsoil
         mliq(iz) = sh2o(iz)*dzsnso(iz)*1000.
      end do

      watmin = 0.01           ! mm
      do iz = 1, nsoil-1
          if (mliq(iz) .lt. 0.) then
             xs = watmin-mliq(iz)
          else
             xs = 0.
          end if
          mliq(iz  ) = mliq(iz  ) + xs
          mliq(iz+1) = mliq(iz+1) - xs
      end do

        iz = nsoil
        if (mliq(iz) .lt. watmin) then
           xs = watmin-mliq(iz)
        else
           xs = 0.
        end if
        mliq(iz) = mliq(iz) + xs
        runsub   = runsub - xs/dt
        if(opt_run == 5)deeprech = deeprech - xs*1.e-3

      do iz = 1, nsoil
        sh2o(iz)     = mliq(iz) / (dzsnso(iz)*1000.)
      end do
   end if

  end subroutine soilwater

!== begin zwteq ====================================================================================

!>\ingroup NoahMP_LSM
  subroutine zwteq (parameters,nsoil  ,nsnow  ,zsoil  ,dzsnso ,sh2o   ,zwt)
! ----------------------------------------------------------------------
! calculate equilibrium water table depth (niu et al., 2005)
! ----------------------------------------------------------------------
  implicit none
! ----------------------------------------------------------------------
! input

  type (noahmp_parameters), intent(in) :: parameters
  integer,                         intent(in) :: nsoil  !no. of soil layers
  integer,                         intent(in) :: nsnow  !maximum no. of snow layers
  real, dimension(1:nsoil),        intent(in) :: zsoil  !depth of soil layer-bottom [m]
  real, dimension(-nsnow+1:nsoil), intent(in) :: dzsnso !snow/soil layer depth [m]
  real, dimension(1:nsoil),        intent(in) :: sh2o   !soil liquid water content [m3/m3]

! output

  real,                           intent(out) :: zwt    !water table depth [m]

! locals

  integer :: k                      !do-loop index
  integer, parameter :: nfine = 100 !no. of fine soil layers of 6m soil
  real    :: wd1                    !water deficit from coarse (4-l) soil moisture profile
  real    :: wd2                    !water deficit from fine (100-l) soil moisture profile
  real    :: dzfine                 !layer thickness of the 100-l soil layers to 6.0 m
  real    :: temp                   !temporary variable
  real, dimension(1:nfine) :: zfine !layer-bottom depth of the 100-l soil layers to 6.0 m
! ----------------------------------------------------------------------

   wd1 = 0.
   do k = 1,nsoil
     wd1 = wd1 + (parameters%smcmax-sh2o(k)) * dzsnso(k) ! [m]
   enddo

   dzfine = 3.0 * (-zsoil(nsoil)) / nfine  
   do k =1,nfine
      zfine(k) = float(k) * dzfine
   enddo

   zwt = -3.*zsoil(nsoil) - 0.001   ! initial value [m]

   wd2 = 0.
   do k = 1,nfine
     temp  = 1. + (zwt-zfine(k))/parameters%psisat
     wd2   = wd2 + parameters%smcmax*(1.-temp**(-1./parameters%bexp))*dzfine
     if(abs(wd2-wd1).le.0.01) then
        zwt = zfine(k)
        exit
     endif
   enddo

  end subroutine zwteq

!== begin infil ====================================================================================

!>\ingroup NoahMP_LSM
  subroutine infil (parameters,nsoil  ,dt     ,zsoil  ,sh2o   ,sice   , & !in
                    sicemax,qinsur ,                         & !in
                    pddum  ,runsrf )                           !out
! --------------------------------------------------------------------------------
! compute inflitration rate at soil surface and surface runoff
! --------------------------------------------------------------------------------
    implicit none
! --------------------------------------------------------------------------------
! inputs
  type (noahmp_parameters), intent(in) :: parameters
  integer,                  intent(in) :: nsoil  !no. of soil layers
  real,                     intent(in) :: dt     !time step (sec)
  real, dimension(1:nsoil), intent(in) :: zsoil  !depth of soil layer-bottom [m]
  real, dimension(1:nsoil), intent(in) :: sh2o   !soil liquid water content [m3/m3]
  real, dimension(1:nsoil), intent(in) :: sice   !soil ice content [m3/m3]
  real,                     intent(in) :: qinsur !water input on soil surface [mm/s]
  real,                     intent(in) :: sicemax!maximum soil ice content (m3/m3)

! outputs
  real,                    intent(out) :: runsrf !surface runoff [mm/s] 
  real,                    intent(out) :: pddum  !infiltration rate at surface

! locals
  integer :: ialp1, j, jj,  k
  real                     :: val
  real                     :: ddt
  real                     :: px
  real                     :: dt1, dd, dice
  real                     :: fcr
  real                     :: sum
  real                     :: acrt
  real                     :: wdf
  real                     :: wcnd
  real                     :: smcav
  real                     :: infmax
  real, dimension(1:nsoil) :: dmax
  integer, parameter       :: cvfrz = 3
! --------------------------------------------------------------------------------

    if (qinsur >  0.0) then
       dt1 = dt /86400.
       smcav = parameters%smcmax - parameters%smcwlt

! maximum infiltration rate

       dmax(1)= -zsoil(1) * smcav
       dice   = -zsoil(1) * sice(1)
       dmax(1)= dmax(1)* (1.0-(sh2o(1) + sice(1) - parameters%smcwlt)/smcav)

       dd = dmax(1)

       do k = 2,nsoil
          dice    = dice + (zsoil(k-1) - zsoil(k) ) * sice(k)
          dmax(k) = (zsoil(k-1) - zsoil(k)) * smcav
          dmax(k) = dmax(k) * (1.0-(sh2o(k) + sice(k) - parameters%smcwlt)/smcav)
          dd      = dd + dmax(k)
       end do

       val = (1. - exp ( - parameters%kdt * dt1))
       ddt = dd * val
       px  = max(0.,qinsur * dt)
       infmax = (px * (ddt / (px + ddt)))/ dt

! impermeable fraction due to frozen soil

       fcr = 1.
       if (dice >  1.e-2) then
          acrt = cvfrz * parameters%frzx / dice
          sum = 1.
          ialp1 = cvfrz - 1
          do j = 1,ialp1
             k = 1
             do jj = j +1,ialp1
                k = k * jj
             end do
             sum = sum + (acrt ** (cvfrz - j)) / float(k)
          end do
          fcr = 1. - exp (-acrt) * sum
       end if

! correction of infiltration limitation

       infmax = infmax * fcr

! jref for urban areas
!       if ( parameters%urban_flag ) infmax == infmax * 0.05

       call wdfcnd2 (parameters,wdf,wcnd,sh2o(1),sicemax)
       infmax = max (infmax,wcnd)
       infmax = min (infmax,px)

       runsrf= max(0., qinsur - infmax)
       pddum = qinsur - runsrf

    end if

  end subroutine infil

!== begin srt ======================================================================================

!>\ingroup NoahMP_LSM
  subroutine srt (parameters,nsoil  ,zsoil  ,dt     ,pddum  ,etrani , & !in
                  qseva  ,sh2o   ,smc    ,zwt    ,fcr    , & !in
                  sicemax,fcrmax ,iloc   ,jloc   ,smcwtd ,         & !in
                  rhstt  ,ai     ,bi     ,ci     ,qdrain , & !out
                  wcnd   )                                   !out
! ----------------------------------------------------------------------
! calculate the right hand side of the time tendency term of the soil
! water diffusion equation.  also to compute ( prepare ) the matrix
! coefficients for the tri-diagonal matrix of the implicit time scheme.
! ----------------------------------------------------------------------
    implicit none
! ----------------------------------------------------------------------
!input

  type (noahmp_parameters), intent(in) :: parameters
    integer,                  intent(in)  :: iloc   !grid index
    integer,                  intent(in)  :: jloc   !grid index
    integer,                  intent(in)  :: nsoil
    real, dimension(1:nsoil), intent(in)  :: zsoil
    real,                     intent(in)  :: dt
    real,                     intent(in)  :: pddum
    real,                     intent(in)  :: qseva
    real, dimension(1:nsoil), intent(in)  :: etrani
    real, dimension(1:nsoil), intent(in)  :: sh2o
    real, dimension(1:nsoil), intent(in)  :: smc
    real,                     intent(in)  :: zwt    ! water table depth [m]
    real, dimension(1:nsoil), intent(in)  :: fcr
    real, intent(in)                      :: fcrmax !maximum of fcr (-)
    real,                     intent(in)  :: sicemax!maximum soil ice content (m3/m3)
    real,                     intent(in)  :: smcwtd !soil moisture between bottom of the soil and the water table

! output

    real, dimension(1:nsoil), intent(out) :: rhstt
    real, dimension(1:nsoil), intent(out) :: ai
    real, dimension(1:nsoil), intent(out) :: bi
    real, dimension(1:nsoil), intent(out) :: ci
    real, dimension(1:nsoil), intent(out) :: wcnd    !hydraulic conductivity (m/s)
    real,                     intent(out) :: qdrain  !bottom drainage (m/s)

! local
    integer                               :: k
    real, dimension(1:nsoil)              :: ddz
    real, dimension(1:nsoil)              :: denom
    real, dimension(1:nsoil)              :: dsmdz
    real, dimension(1:nsoil)              :: wflux
    real, dimension(1:nsoil)              :: wdf
    real, dimension(1:nsoil)              :: smx
    real                                  :: temp1
    real                                  :: smxwtd !soil moisture between bottom of the soil and water table
    real                                  :: smxbot  !soil moisture below bottom to calculate flux

! niu and yang (2006), j. of hydrometeorology
! ----------------------------------------------------------------------

    if(opt_inf == 1) then
      do k = 1, nsoil
        call wdfcnd1 (parameters,wdf(k),wcnd(k),smc(k),fcr(k))
        smx(k) = smc(k)
      end do
        if(opt_run == 5)smxwtd=smcwtd
    end if

    if(opt_inf == 2) then
      do k = 1, nsoil
        call wdfcnd2 (parameters,wdf(k),wcnd(k),sh2o(k),sicemax)
        smx(k) = sh2o(k)
      end do
          if(opt_run == 5)smxwtd=smcwtd*sh2o(nsoil)/smc(nsoil)  !same liquid fraction as in the bottom layer
    end if

    do k = 1, nsoil
       if(k == 1) then
          denom(k) = - zsoil (k)
          temp1    = - zsoil (k+1)
          ddz(k)   = 2.0 / temp1
          dsmdz(k) = 2.0 * (smx(k) - smx(k+1)) / temp1
          wflux(k) = wdf(k) * dsmdz(k) + wcnd(k) - pddum + etrani(k) + qseva
       else if (k < nsoil) then
          denom(k) = (zsoil(k-1) - zsoil(k))
          temp1    = (zsoil(k-1) - zsoil(k+1))
          ddz(k)   = 2.0 / temp1
          dsmdz(k) = 2.0 * (smx(k) - smx(k+1)) / temp1
          wflux(k) = wdf(k  ) * dsmdz(k  ) + wcnd(k  )         &
                   - wdf(k-1) * dsmdz(k-1) - wcnd(k-1) + etrani(k)
       else
          denom(k) = (zsoil(k-1) - zsoil(k))
          if(opt_run == 1 .or. opt_run == 2) then
             qdrain   = 0.
          end if
          if(opt_run == 3) then
             qdrain   = parameters%slope*wcnd(k)
          end if
          if(opt_run == 4) then
             qdrain   = (1.0-fcrmax)*wcnd(k)
          end if
          if(opt_run == 5) then   !gmm new m-m&f water table dynamics formulation
             temp1    = 2.0 * denom(k)
             if(zwt < zsoil(nsoil)-denom(nsoil))then
!gmm interpolate from below, midway to the water table, to the middle of the auxiliary layer below the soil bottom
                smxbot = smx(k) - (smx(k)-smxwtd) *  denom(k) * 2./ (denom(k) + zsoil(k) - zwt)
             else
                smxbot = smxwtd
             endif
             dsmdz(k) = 2.0 * (smx(k) - smxbot) / temp1
             qdrain   = wdf(k  ) * dsmdz(k  ) + wcnd(k  )
          end if   
          wflux(k) = -(wdf(k-1)*dsmdz(k-1))-wcnd(k-1)+etrani(k) + qdrain
       end if
    end do

    do k = 1, nsoil
       if(k == 1) then
          ai(k)    =   0.0
          bi(k)    =   wdf(k  ) * ddz(k  ) / denom(k)
          ci(k)    = - bi (k)
       else if (k < nsoil) then
          ai(k)    = - wdf(k-1) * ddz(k-1) / denom(k)
          ci(k)    = - wdf(k  ) * ddz(k  ) / denom(k)
          bi(k)    = - ( ai (k) + ci (k) )
       else
          ai(k)    = - wdf(k-1) * ddz(k-1) / denom(k)
          ci(k)    = 0.0
          bi(k)    = - ( ai (k) + ci (k) )
       end if
          rhstt(k) = wflux(k) / (-denom(k))
    end do

! ----------------------------------------------------------------------
  end subroutine srt

!== begin sstep ====================================================================================

!>\ingroup NoahMP_LSM
  subroutine sstep (parameters,nsoil  ,nsnow  ,dt     ,zsoil  ,dzsnso , & !in
                    sice   ,iloc   ,jloc   ,zwt            ,                 & !in
                    sh2o   ,smc    ,ai     ,bi     ,ci     , & !inout
                    rhstt  ,smcwtd ,qdrain ,deeprech,                                 & !inout
                    wplus  )                                   !out

! ----------------------------------------------------------------------
! calculate/update soil moisture content values 
! ----------------------------------------------------------------------
    implicit none
! ----------------------------------------------------------------------
!input

  type (noahmp_parameters), intent(in) :: parameters
    integer,                         intent(in) :: iloc   !grid index
    integer,                         intent(in) :: jloc   !grid index
    integer,                         intent(in) :: nsoil  !
    integer,                         intent(in) :: nsnow  !
    real, intent(in)                            :: dt
    real, intent(in)                            :: zwt
    real, dimension(       1:nsoil), intent(in) :: zsoil
    real, dimension(       1:nsoil), intent(in) :: sice
    real, dimension(-nsnow+1:nsoil), intent(in) :: dzsnso ! snow/soil layer thickness [m]

!input and output
    real, dimension(1:nsoil), intent(inout) :: sh2o
    real, dimension(1:nsoil), intent(inout) :: smc
    real, dimension(1:nsoil), intent(inout) :: ai
    real, dimension(1:nsoil), intent(inout) :: bi
    real, dimension(1:nsoil), intent(inout) :: ci
    real, dimension(1:nsoil), intent(inout) :: rhstt
    real                    , intent(inout) :: smcwtd
    real                    , intent(inout) :: qdrain
    real                    , intent(inout) :: deeprech

!output
    real, intent(out)                       :: wplus     !saturation excess water (m)

!local
    integer                                 :: k
    real, dimension(1:nsoil)                :: rhsttin
    real, dimension(1:nsoil)                :: ciin
    real                                    :: stot
    real                                    :: epore
    real                                    :: wminus
! ----------------------------------------------------------------------
    wplus = 0.0

    do k = 1,nsoil
       rhstt (k) =   rhstt(k) * dt
       ai (k)    =      ai(k) * dt
       bi (k)    = 1. + bi(k) * dt
       ci (k)    =      ci(k) * dt
    end do

! copy values for input variables before calling rosr12

    do k = 1,nsoil
       rhsttin(k) = rhstt(k)
       ciin(k)    = ci(k)
    end do

! call rosr12 to solve the tri-diagonal matrix

    call rosr12 (ci,ai,bi,ciin,rhsttin,rhstt,1,nsoil,0)

    do k = 1,nsoil
        sh2o(k) = sh2o(k) + ci(k)
    enddo

!  excessive water above saturation in a layer is moved to
!  its unsaturated layer like in a bucket

!gmmwith opt_run=5 there is soil moisture below nsoil, to the water table
  if(opt_run == 5) then

!update smcwtd

     if(zwt < zsoil(nsoil)-dzsnso(nsoil))then
!accumulate qdrain to update deep water table and soil moisture later
        deeprech =  deeprech + dt * qdrain
     else
        smcwtd = smcwtd + dt * qdrain  / dzsnso(nsoil)
        wplus        = max((smcwtd-parameters%smcmax), 0.0) * dzsnso(nsoil)
        wminus       = max((1.e-4-smcwtd), 0.0) * dzsnso(nsoil)

        smcwtd = max( min(smcwtd,parameters%smcmax) , 1.e-4)
        sh2o(nsoil)    = sh2o(nsoil) + wplus/dzsnso(nsoil)

!reduce fluxes at the bottom boundaries accordingly
        qdrain = qdrain - wplus/dt
        deeprech = deeprech - wminus
     endif

  endif

    do k = nsoil,2,-1
      epore        = max ( 1.e-4 , ( parameters%smcmax - sice(k) ) )
      wplus        = max((sh2o(k)-epore), 0.0) * dzsnso(k)
      sh2o(k)      = min(epore,sh2o(k))
      sh2o(k-1)    = sh2o(k-1) + wplus/dzsnso(k-1)
    end do

    epore        = max ( 1.e-4 , ( parameters%smcmax - sice(1) ) )
    wplus        = max((sh2o(1)-epore), 0.0) * dzsnso(1) 
    sh2o(1)      = min(epore,sh2o(1))

  end subroutine sstep

!== begin wdfcnd1 ==================================================================================

!>\ingroup NoahMP_LSM
  subroutine wdfcnd1 (parameters,wdf,wcnd,smc,fcr)
! ----------------------------------------------------------------------
! calculate soil water diffusivity and soil hydraulic conductivity.
! ----------------------------------------------------------------------
    implicit none
! ----------------------------------------------------------------------
! input 
  type (noahmp_parameters), intent(in) :: parameters
    real,intent(in)  :: smc
    real,intent(in)  :: fcr

! output
    real,intent(out) :: wcnd
    real,intent(out) :: wdf

! local
    real :: expon
    real :: factr
    real :: vkwgt
! ----------------------------------------------------------------------

! soil water diffusivity

    factr = max(0.01, smc/parameters%smcmax)
    expon = parameters%bexp + 2.0
    wdf   = parameters%dwsat * factr ** expon
    wdf   = wdf * (1.0 - fcr)

! hydraulic conductivity

    expon = 2.0*parameters%bexp + 3.0
    wcnd  = parameters%dksat * factr ** expon
    wcnd  = wcnd * (1.0 - fcr)

  end subroutine wdfcnd1

!== begin wdfcnd2 ==================================================================================

!>\ingroup NoahMP_LSM
  subroutine wdfcnd2 (parameters,wdf,wcnd,smc,sice)
! ----------------------------------------------------------------------
! calculate soil water diffusivity and soil hydraulic conductivity.
! ----------------------------------------------------------------------
    implicit none
! ----------------------------------------------------------------------
! input
  type (noahmp_parameters), intent(in) :: parameters
    real,intent(in)  :: smc
    real,intent(in)  :: sice

! output
    real,intent(out) :: wcnd
    real,intent(out) :: wdf

! local
    real :: expon
    real :: factr
    real :: vkwgt
! ----------------------------------------------------------------------

! soil water diffusivity

    factr = max(0.01, smc/parameters%smcmax)
    expon = parameters%bexp + 2.0
    wdf   = parameters%dwsat * factr ** expon

    if (sice > 0.0) then
    vkwgt = 1./ (1. + (500.* sice)**3.)
    wdf   = vkwgt * wdf + (1.-vkwgt)*parameters%dwsat*(0.2/parameters%smcmax)**expon
    end if

! hydraulic conductivity

    expon = 2.0*parameters%bexp + 3.0
    wcnd  = parameters%dksat * factr ** expon

  end subroutine wdfcnd2

!== begin groundwater ==============================================================================

!>\ingroup NoahMP_LSM
  subroutine groundwater(parameters,nsnow  ,nsoil  ,dt     ,sice   ,zsoil  , & !in
                         stc    ,wcnd   ,fcrmax ,iloc   ,jloc   , & !in
                         sh2o   ,zwt    ,wa     ,wt     ,         & !inout
                         qin    ,qdis   )                           !out
! ----------------------------------------------------------------------
  implicit none
! ----------------------------------------------------------------------
! input
  type (noahmp_parameters), intent(in) :: parameters
  integer,                         intent(in) :: iloc  !grid index
  integer,                         intent(in) :: jloc  !grid index
  integer,                         intent(in) :: nsnow !maximum no. of snow layers
  integer,                         intent(in) :: nsoil !no. of soil layers
  real,                            intent(in) :: dt    !timestep [sec]
  real,                            intent(in) :: fcrmax!maximum fcr (-)
  real, dimension(       1:nsoil), intent(in) :: sice  !soil ice content [m3/m3]
  real, dimension(       1:nsoil), intent(in) :: zsoil !depth of soil layer-bottom [m]
  real, dimension(       1:nsoil), intent(in) :: wcnd  !hydraulic conductivity (m/s)
  real, dimension(-nsnow+1:nsoil), intent(in) :: stc   !snow/soil temperature (k)

! input and output
  real, dimension(    1:nsoil), intent(inout) :: sh2o  !liquid soil water [m3/m3]
  real,                         intent(inout) :: zwt   !the depth to water table [m]
  real,                         intent(inout) :: wa    !water storage in aquifer [mm]
  real,                         intent(inout) :: wt    !water storage in aquifer 
                                                           !+ saturated soil [mm]
! output
  real,                           intent(out) :: qin   !groundwater recharge [mm/s]
  real,                           intent(out) :: qdis  !groundwater discharge [mm/s]

! local
  real                                        :: fff   !runoff decay factor (m-1)
  real                                        :: rsbmx !baseflow coefficient [mm/s]
  integer                                     :: iz    !do-loop index
  integer                                     :: iwt   !layer index above water table layer
  real,  dimension(    1:nsoil)               :: dzmm  !layer thickness [mm]
  real,  dimension(    1:nsoil)               :: znode !node depth [m]
  real,  dimension(    1:nsoil)               :: mliq  !liquid water mass [kg/m2 or mm]
  real,  dimension(    1:nsoil)               :: epore !effective porosity [-]
  real,  dimension(    1:nsoil)               :: hk    !hydraulic conductivity [mm/s]
  real,  dimension(    1:nsoil)               :: smc   !total soil water  content [m3/m3]
  real(kind=8)                                :: s_node!degree of saturation of iwt layer
  real                                        :: dzsum !cumulative depth above water table [m]
  real                                        :: smpfz !matric potential (frozen effects) [mm]
  real                                        :: ka    !aquifer hydraulic conductivity [mm/s]
  real                                        :: wh_zwt!water head at water table [mm]
  real                                        :: wh    !water head at layer above zwt [mm]
  real                                        :: ws    !water used to fill air pore [mm]
  real                                        :: wtsub !sum of hk*dzmm
  real                                        :: watmin!minimum soil vol soil moisture [m3/m3]
  real                                        :: xs    !excessive water above saturation [mm]
  real, parameter                             :: rous = 0.2    !specific yield [-]
  real, parameter                             :: cmic = 0.20   !microprore content (0.0-1.0)
                                                               !0.0-close to free drainage
! -------------------------------------------------------------
      qdis      = 0.0
      qin       = 0.0

! derive layer-bottom depth in [mm]
!kwm:  derive layer thickness in mm

      dzmm(1) = -zsoil(1)*1.e3
      do iz = 2, nsoil
         dzmm(iz)  = 1.e3 * (zsoil(iz - 1) - zsoil(iz))
      enddo

! derive node (middle) depth in [m]
!kwm:  positive number, depth below ground surface in m
      znode(1) = -zsoil(1) / 2.
      do iz = 2, nsoil
         znode(iz)  = -zsoil(iz-1) + 0.5 * (zsoil(iz-1) - zsoil(iz))
      enddo

! convert volumetric soil moisture "sh2o" to mass

      do iz = 1, nsoil
         smc(iz)      = sh2o(iz) + sice(iz)
         mliq(iz)     = sh2o(iz) * dzmm(iz)
         epore(iz)    = max(0.01,parameters%smcmax - sice(iz))
         hk(iz)       = 1.e3*wcnd(iz)
      enddo

! the layer index of the first unsaturated layer,
! i.e., the layer right above the water table

      iwt = nsoil
      do iz = 2,nsoil
         if(zwt   .le. -zsoil(iz) ) then
            iwt = iz-1
            exit
         end if
      enddo

! groundwater discharge [mm/s]

      fff   = 6.0
      rsbmx = 5.0

      qdis = (1.0-fcrmax)*rsbmx*exp(-parameters%timean)*exp(-fff*(zwt-2.0))

! matric potential at the layer above the water table

      s_node = min(1.0,smc(iwt)/parameters%smcmax )
      s_node = max(s_node,real(0.01,kind=8))
      smpfz  = -parameters%psisat*1000.*s_node**(-parameters%bexp)   ! m --> mm
      smpfz  = max(-120000.0,cmic*smpfz)   

! recharge rate qin to groundwater

      ka  = hk(iwt)

      wh_zwt  = - zwt * 1.e3                          !(mm)
      wh      = smpfz  - znode(iwt)*1.e3              !(mm)
      qin     = - ka * (wh_zwt-wh)  /((zwt-znode(iwt))*1.e3)
      qin     = max(-10.0/dt,min(10./dt,qin))
     
! water storage in the aquifer + saturated soil

      wt  = wt + (qin - qdis) * dt     !(mm)

      if(iwt.eq.nsoil) then
         wa          = wa + (qin - qdis) * dt     !(mm)
         wt          = wa
         zwt         = (-zsoil(nsoil) + 25.) - wa/1000./rous      !(m)
         mliq(nsoil) = mliq(nsoil) - qin * dt        ! [mm]

         mliq(nsoil) = mliq(nsoil) + max(0.,(wa - 5000.))
         wa          = min(wa, 5000.)
      else
         
         if (iwt.eq.nsoil-1) then
            zwt = -zsoil(nsoil)                   &
                 - (wt-rous*1000*25.) / (epore(nsoil))/1000.
         else
            ws = 0.   ! water used to fill soil air pores
            do iz = iwt+2,nsoil
               ws = ws + epore(iz) * dzmm(iz)
            enddo
            zwt = -zsoil(iwt+1)                  &
                  - (wt-rous*1000.*25.-ws) /(epore(iwt+1))/1000.
         endif

         wtsub = 0.
         do iz = 1, nsoil
           wtsub = wtsub + hk(iz)*dzmm(iz)
         end do

         do iz = 1, nsoil           ! removing subsurface runoff
         mliq(iz) = mliq(iz) - qdis*dt*hk(iz)*dzmm(iz)/wtsub
         end do
      end if

      zwt = max(1.5,zwt)

!
! limit mliq to be greater than or equal to watmin.
! get water needed to bring mliq equal watmin from lower layer.
!
      watmin = 0.01
      do iz = 1, nsoil-1
          if (mliq(iz) .lt. 0.) then
             xs = watmin-mliq(iz)
          else
             xs = 0.
          end if
          mliq(iz  ) = mliq(iz  ) + xs
          mliq(iz+1) = mliq(iz+1) - xs
      end do

        iz = nsoil
        if (mliq(iz) .lt. watmin) then
           xs = watmin-mliq(iz)
        else
           xs = 0.
        end if
        mliq(iz) = mliq(iz) + xs
        wa       = wa - xs
        wt       = wt - xs

      do iz = 1, nsoil
        sh2o(iz)     = mliq(iz) / dzmm(iz)
      end do

  end subroutine groundwater

!== begin shallowwatertable ========================================================================

!>\ingroup NoahMP_LSM
  subroutine shallowwatertable (parameters,nsnow  ,nsoil  ,zsoil, dt    , & !in
                         dzsnso ,smceq ,iloc   ,jloc         , & !in
                         smc    ,wtd   ,smcwtd ,rech, qdrain  )  !inout
! ----------------------------------------------------------------------
!diagnoses water table depth and computes recharge when the water table is within the resolved soil layers,
!according to the miguez-macho&fan scheme
! ----------------------------------------------------------------------
  implicit none
! ----------------------------------------------------------------------
! input
  type (noahmp_parameters), intent(in) :: parameters
  integer,                         intent(in) :: nsnow !maximum no. of snow layers
  integer,                         intent(in) :: nsoil !no. of soil layers
  integer,                         intent(in) :: iloc,jloc
  real,                            intent(in) :: dt
  real, dimension(       1:nsoil), intent(in) :: zsoil !depth of soil layer-bottom [m]
  real, dimension(-nsnow+1:nsoil), intent(in) :: dzsnso ! snow/soil layer thickness [m]
  real,  dimension(      1:nsoil), intent(in) :: smceq  !equilibrium soil water  content [m3/m3]

! input and output
  real,  dimension(      1:nsoil), intent(inout) :: smc   !total soil water  content [m3/m3]
  real,                         intent(inout) :: wtd   !the depth to water table [m]
  real,                         intent(inout) :: smcwtd   !soil moisture between bottom of the soil and the water table [m3/m3]
  real,                         intent(out) :: rech ! groundwater recharge (net vertical flux across the water table), positive up
  real,                         intent(inout) :: qdrain
    
! local
  integer                                     :: iz    !do-loop index
  integer                                     :: iwtd   !layer index above water table layer
  integer                                     :: kwtd   !layer index where the water table layer is
  real                                        :: wtdold
  real                                        :: dzup
  real                                        :: smceqdeep
  real,  dimension(       0:nsoil)            :: zsoil0
! -------------------------------------------------------------


zsoil0(1:nsoil) = zsoil(1:nsoil)
zsoil0(0) = 0.         
 
!find the layer where the water table is
     do iz=nsoil,1,-1
        if(wtd + 1.e-6 < zsoil0(iz)) exit
     enddo
        iwtd=iz

        
        kwtd=iwtd+1  !layer where the water table is
        if(kwtd.le.nsoil)then    !wtd in the resolved layers
           wtdold=wtd
           if(smc(kwtd).gt.smceq(kwtd))then
        
               if(smc(kwtd).eq.parameters%smcmax)then !wtd went to the layer above
                      wtd=zsoil0(iwtd)
                      rech=-(wtdold-wtd) * (parameters%smcmax-smceq(kwtd))
                      iwtd=iwtd-1
                      kwtd=kwtd-1
                   if(kwtd.ge.1)then
                      if(smc(kwtd).gt.smceq(kwtd))then
                      wtdold=wtd
                      wtd = min( ( smc(kwtd)*dzsnso(kwtd) &
                        - smceq(kwtd)*zsoil0(iwtd) + parameters%smcmax*zsoil0(kwtd) ) / &
                        ( parameters%smcmax-smceq(kwtd) ), zsoil0(iwtd))
                      rech=rech-(wtdold-wtd) * (parameters%smcmax-smceq(kwtd))
                      endif
                   endif
               else  !wtd stays in the layer
                      wtd = min( ( smc(kwtd)*dzsnso(kwtd) &
                        - smceq(kwtd)*zsoil0(iwtd) + parameters%smcmax*zsoil0(kwtd) ) / &
                        ( parameters%smcmax-smceq(kwtd) ), zsoil0(iwtd))
                      rech=-(wtdold-wtd) * (parameters%smcmax-smceq(kwtd))
               endif
           
           else    !wtd has gone down to the layer below
               wtd=zsoil0(kwtd)
               rech=-(wtdold-wtd) * (parameters%smcmax-smceq(kwtd))
               kwtd=kwtd+1
               iwtd=iwtd+1
!wtd crossed to the layer below. now adjust it there
               if(kwtd.le.nsoil)then
                   wtdold=wtd
                   if(smc(kwtd).gt.smceq(kwtd))then
                   wtd = min( ( smc(kwtd)*dzsnso(kwtd) &
                   - smceq(kwtd)*zsoil0(iwtd) + parameters%smcmax*zsoil0(kwtd) ) / &
                       ( parameters%smcmax-smceq(kwtd) ) , zsoil0(iwtd) )
                   else
                   wtd=zsoil0(kwtd)
                   endif
                   rech = rech - (wtdold-wtd) * &
                                 (parameters%smcmax-smceq(kwtd))

                else
                   wtdold=wtd
!restore smoi to equilibrium value with water from the ficticious layer below
!                   smcwtd=smcwtd-(smceq(nsoil)-smc(nsoil))
!                   qdrain = qdrain - 1000 * (smceq(nsoil)-smc(nsoil)) * dzsnso(nsoil) / dt
!                   smc(nsoil)=smceq(nsoil)
!adjust wtd in the ficticious layer below
                   smceqdeep = parameters%smcmax * ( -parameters%psisat / ( -parameters%psisat - dzsnso(nsoil) ) ) ** (1./parameters%bexp)
                   wtd = min( ( smcwtd*dzsnso(nsoil) &
                   - smceqdeep*zsoil0(nsoil) + parameters%smcmax*(zsoil0(nsoil)-dzsnso(nsoil)) ) / &
                       ( parameters%smcmax-smceqdeep ) , zsoil0(nsoil) )
                   rech = rech - (wtdold-wtd) * &
                                 (parameters%smcmax-smceqdeep)
                endif
            
            endif
        elseif(wtd.ge.zsoil0(nsoil)-dzsnso(nsoil))then
!if wtd was already below the bottom of the resolved soil crust
           wtdold=wtd
           smceqdeep = parameters%smcmax * ( -parameters%psisat / ( -parameters%psisat - dzsnso(nsoil) ) ) ** (1./parameters%bexp)
           if(smcwtd.gt.smceqdeep)then
               wtd = min( ( smcwtd*dzsnso(nsoil) &
                 - smceqdeep*zsoil0(nsoil) + parameters%smcmax*(zsoil0(nsoil)-dzsnso(nsoil)) ) / &
                     ( parameters%smcmax-smceqdeep ) , zsoil0(nsoil) )
               rech = -(wtdold-wtd) * (parameters%smcmax-smceqdeep)
           else
               rech = -(wtdold-(zsoil0(nsoil)-dzsnso(nsoil))) * (parameters%smcmax-smceqdeep)
               wtdold=zsoil0(nsoil)-dzsnso(nsoil)
!and now even further down
               dzup=(smceqdeep-smcwtd)*dzsnso(nsoil)/(parameters%smcmax-smceqdeep)
               wtd=wtdold-dzup
               rech = rech - (parameters%smcmax-smceqdeep)*dzup
               smcwtd=smceqdeep
           endif

         
         endif

if(iwtd.lt.nsoil)smcwtd=parameters%smcmax

end  subroutine shallowwatertable

! ==================================================================================================
! ********************* end of water subroutines ******************************************
! ==================================================================================================

!== begin carbon ===================================================================================

!>\ingroup NoahMP_LSM
  subroutine carbon (parameters,nsnow  ,nsoil  ,vegtyp ,dt     ,zsoil  , & !in
                     dzsnso ,stc    ,smc    ,tv     ,tg     ,psn    , & !in
                     foln   ,btran  ,apar   ,fveg   ,igs    , & !in
                     troot  ,ist    ,lat    ,iloc   ,jloc   , & !in
                     lfmass ,rtmass ,stmass ,wood   ,stblcp ,fastcp , & !inout
                     gpp    ,npp    ,nee    ,autors ,heters ,totsc  , & !out
                     totlb  ,xlai   ,xsai   )                   !out
! ------------------------------------------------------------------------------------------
      implicit none
! ------------------------------------------------------------------------------------------
! inputs (carbon)

  type (noahmp_parameters), intent(in) :: parameters
  integer                        , intent(in) :: iloc   !grid index
  integer                        , intent(in) :: jloc   !grid index
  integer                        , intent(in) :: vegtyp !vegetation type 
  integer                        , intent(in) :: nsnow  !number of snow layers
  integer                        , intent(in) :: nsoil  !number of soil layers
  real                           , intent(in) :: lat    !latitude (radians)
  real                           , intent(in) :: dt     !time step (s)
  real, dimension(       1:nsoil), intent(in) :: zsoil  !depth of layer-bottom from soil surface
  real, dimension(-nsnow+1:nsoil), intent(in) :: dzsnso !snow/soil layer thickness [m]
  real, dimension(-nsnow+1:nsoil), intent(in) :: stc    !snow/soil temperature [k]
  real, dimension(       1:nsoil), intent(in) :: smc    !soil moisture (ice + liq.) [m3/m3]
  real                           , intent(in) :: tv     !vegetation temperature (k)
  real                           , intent(in) :: tg     !ground temperature (k)
  real                           , intent(in) :: foln   !foliage nitrogen (%)
  real                           , intent(in) :: btran  !soil water transpiration factor (0 to 1)
  real                           , intent(in) :: psn    !total leaf photosyn (umolco2/m2/s) [+]
  real                           , intent(in) :: apar   !par by canopy (w/m2)
  real                           , intent(in) :: igs    !growing season index (0=off, 1=on)
  real                           , intent(in) :: fveg   !vegetation greenness fraction
  real                           , intent(in) :: troot  !root-zone averaged temperature (k)
  integer                        , intent(in) :: ist    !surface type 1->soil; 2->lake

! input & output (carbon)

  real                        , intent(inout) :: lfmass !leaf mass [g/m2]
  real                        , intent(inout) :: rtmass !mass of fine roots [g/m2]
  real                        , intent(inout) :: stmass !stem mass [g/m2]
  real                        , intent(inout) :: wood   !mass of wood (incl. woody roots) [g/m2]
  real                        , intent(inout) :: stblcp !stable carbon in deep soil [g/m2]
  real                        , intent(inout) :: fastcp !short-lived carbon in shallow soil [g/m2]

! outputs: (carbon)

  real                          , intent(out) :: gpp    !net instantaneous assimilation [g/m2/s c]
  real                          , intent(out) :: npp    !net primary productivity [g/m2/s c]
  real                          , intent(out) :: nee    !net ecosystem exchange [g/m2/s co2]
  real                          , intent(out) :: autors !net ecosystem respiration [g/m2/s c]
  real                          , intent(out) :: heters !organic respiration [g/m2/s c]
  real                          , intent(out) :: totsc  !total soil carbon [g/m2 c]
  real                          , intent(out) :: totlb  !total living carbon ([g/m2 c]
  real                          , intent(out) :: xlai   !leaf area index [-]
  real                          , intent(out) :: xsai   !stem area index [-]
!  real                          , intent(out) :: vocflx(5) ! voc fluxes [ug c m-2 h-1]

! local variables

  integer :: j         !do-loop index
  real    :: wroot     !root zone soil water [-]
  real    :: wstres    !water stress coeficient [-]  (1. for wilting )
  real    :: lapm      !leaf area per unit mass [m2/g]
! ------------------------------------------------------------------------------------------

   if ( ( vegtyp == parameters%iswater ) .or. ( vegtyp == parameters%isbarren ) .or. &
        ( vegtyp == parameters%isice ) .or. (parameters%urban_flag) ) then
      xlai   = 0.
      xsai   = 0.
      gpp    = 0.
      npp    = 0.
      nee    = 0.
      autors = 0.
      heters = 0.
      totsc  = 0.
      totlb  = 0.
      lfmass = 0.
      rtmass = 0.
      stmass = 0.
      wood   = 0.
      stblcp = 0.
      fastcp = 0.

      return
   end if

      lapm       = parameters%sla / 1000.   ! m2/kg -> m2/g

! water stress

      wstres  = 1.- btran

      wroot  = 0.
      do j=1,parameters%nroot
        wroot = wroot + smc(j)/parameters%smcmax *  dzsnso(j) / (-zsoil(parameters%nroot))
      enddo

  call co2flux (parameters,nsnow  ,nsoil  ,vegtyp ,igs    ,dt     , & !in
                dzsnso ,stc    ,psn    ,troot  ,tv     , & !in
                wroot  ,wstres ,foln   ,lapm   ,         & !in
                lat    ,iloc   ,jloc   ,fveg   ,         & !in
                xlai   ,xsai   ,lfmass ,rtmass ,stmass , & !inout
                fastcp ,stblcp ,wood   ,                 & !inout
                gpp    ,npp    ,nee    ,autors ,heters , & !out
                totsc  ,totlb  )                           !out

!   call bvoc (parameters,vocflx,  vegtyp,  vegfac,   apar,   tv)
!   call ch4

  end subroutine carbon

!== begin co2flux ==================================================================================

!>\ingroup NoahMP_LSM
  subroutine co2flux (parameters,nsnow  ,nsoil  ,vegtyp ,igs    ,dt     , & !in
                      dzsnso ,stc    ,psn    ,troot  ,tv     , & !in
                      wroot  ,wstres ,foln   ,lapm   ,         & !in
                      lat    ,iloc   ,jloc   ,fveg   ,         & !in
                      xlai   ,xsai   ,lfmass ,rtmass ,stmass , & !inout
                      fastcp ,stblcp ,wood   ,                 & !inout
                      gpp    ,npp    ,nee    ,autors ,heters , & !out
                      totsc  ,totlb  )                           !out
! -----------------------------------------------------------------------------------------
! the original code is from re dickinson et al.(1998), modifed by guo-yue niu, 2004
! -----------------------------------------------------------------------------------------
  implicit none
! -----------------------------------------------------------------------------------------

! input

  type (noahmp_parameters), intent(in) :: parameters
  integer                        , intent(in) :: iloc   !grid index
  integer                        , intent(in) :: jloc   !grid index
  integer                        , intent(in) :: vegtyp !vegetation physiology type
  integer                        , intent(in) :: nsnow  !number of snow layers
  integer                        , intent(in) :: nsoil  !number of soil layers
  real                           , intent(in) :: dt     !time step (s)
  real                           , intent(in) :: lat    !latitude (radians)
  real                           , intent(in) :: igs    !growing season index (0=off, 1=on)
  real, dimension(-nsnow+1:nsoil), intent(in) :: dzsnso !snow/soil layer thickness [m]
  real, dimension(-nsnow+1:nsoil), intent(in) :: stc    !snow/soil temperature [k]
  real                           , intent(in) :: psn    !total leaf photosynthesis (umolco2/m2/s)
  real                           , intent(in) :: troot  !root-zone averaged temperature (k)
  real                           , intent(in) :: tv     !leaf temperature (k)
  real                           , intent(in) :: wroot  !root zone soil water
  real                           , intent(in) :: wstres !soil water stress
  real                           , intent(in) :: foln   !foliage nitrogen (%)
  real                           , intent(in) :: lapm   !leaf area per unit mass [m2/g]
  real                           , intent(in) :: fveg   !vegetation greenness fraction

! input and output

  real                        , intent(inout) :: xlai   !leaf  area index from leaf carbon [-]
  real                        , intent(inout) :: xsai   !stem area index from leaf carbon [-]
  real                        , intent(inout) :: lfmass !leaf mass [g/m2]
  real                        , intent(inout) :: rtmass !mass of fine roots [g/m2]
  real                        , intent(inout) :: stmass !stem mass [g/m2]
  real                        , intent(inout) :: fastcp !short lived carbon [g/m2]
  real                        , intent(inout) :: stblcp !stable carbon pool [g/m2]
  real                        , intent(inout) :: wood   !mass of wood (incl. woody roots) [g/m2]

! output

  real                          , intent(out) :: gpp    !net instantaneous assimilation [g/m2/s]
  real                          , intent(out) :: npp    !net primary productivity [g/m2]
  real                          , intent(out) :: nee    !net ecosystem exchange (autors+heters-gpp)
  real                          , intent(out) :: autors !net ecosystem resp. (maintance and growth)
  real                          , intent(out) :: heters !organic respiration
  real                          , intent(out) :: totsc  !total soil carbon (g/m2)
  real                          , intent(out) :: totlb  !total living carbon (g/m2)

! local

  real                   :: cflux    !carbon flux to atmosphere [g/m2/s]
  real                   :: lfmsmn   !minimum leaf mass [g/m2]
  real                   :: rswood   !wood respiration [g/m2]
  real                   :: rsleaf   !leaf maintenance respiration per timestep [g/m2]
  real                   :: rsroot   !fine root respiration per time step [g/m2]
  real                   :: nppl     !leaf net primary productivity [g/m2/s]
  real                   :: nppr     !root net primary productivity [g/m2/s]
  real                   :: nppw     !wood net primary productivity [g/m2/s]
  real                   :: npps     !wood net primary productivity [g/m2/s]
  real                   :: dielf    !death of leaf mass per time step [g/m2]

  real                   :: addnpplf !leaf assimil after resp. losses removed [g/m2]
  real                   :: addnppst !stem assimil after resp. losses removed [g/m2]
  real                   :: carbfx   !carbon assimilated per model step [g/m2]
  real                   :: grleaf   !growth respiration rate for leaf [g/m2/s]
  real                   :: grroot   !growth respiration rate for root [g/m2/s]
  real                   :: grwood   !growth respiration rate for wood [g/m2/s]
  real                   :: grstem   !growth respiration rate for stem [g/m2/s]
  real                   :: leafpt   !fraction of carbon allocated to leaves [-]
  real                   :: lfdel    !maximum  leaf mass  available to change [g/m2/s]
  real                   :: lftovr   !stem turnover per time step [g/m2]
  real                   :: sttovr   !stem turnover per time step [g/m2]
  real                   :: wdtovr   !wood turnover per time step [g/m2]
  real                   :: rssoil   !soil respiration per time step [g/m2]
  real                   :: rttovr   !root carbon loss per time step by turnover [g/m2]
  real                   :: stablc   !decay rate of fast carbon to slow carbon [g/m2/s]
  real                   :: woodf    !calculated wood to root ratio [-]
  real                   :: nonlef   !fraction of carbon to root and wood [-]
  real                   :: rootpt   !fraction of carbon flux to roots [-]
  real                   :: woodpt   !fraction of carbon flux to wood [-]
  real                   :: stempt   !fraction of carbon flux to stem [-]
  real                   :: resp     !leaf respiration [umol/m2/s]
  real                   :: rsstem   !stem respiration [g/m2/s]

  real                   :: fsw      !soil water factor for microbial respiration
  real                   :: fst      !soil temperature factor for microbial respiration
  real                   :: fnf      !foliage nitrogen adjustemt to respiration (<= 1)
  real                   :: tf       !temperature factor
  real                   :: rf       !respiration reduction factor (<= 1)
  real                   :: stdel
  real                   :: stmsmn
  real                   :: sapm     !stem area per unit mass (m2/g)
  real                   :: diest
! -------------------------- constants -------------------------------
  real                   :: bf       !parameter for present wood allocation [-]
  real                   :: rswoodc  !wood respiration coeficient [1/s]
  real                   :: stovrc   !stem turnover coefficient [1/s]
  real                   :: rsdryc   !degree of drying that reduces soil respiration [-]
  real                   :: rtovrc   !root turnover coefficient [1/s]
  real                   :: wstrc    !water stress coeficient [-]
  real                   :: laimin   !minimum leaf area index [m2/m2]
  real                   :: xsamin   !minimum leaf area index [m2/m2]
  real                   :: sc
  real                   :: sd
  real                   :: vegfrac

! respiration as a function of temperature

  real :: r,x
          r(x) = exp(0.08*(x-298.16))
! ---------------------------------------------------------------------------------

! constants
    rtovrc  = 2.0e-8        !original was 2.0e-8
    rsdryc  = 40.0          !original was 40.0
    rswoodc = 3.0e-10       !
    bf      = 0.90          !original was 0.90   ! carbon to roots
    wstrc   = 100.0
    laimin  = 0.05   
    xsamin  = 0.05     ! mb: change to prevent vegetation from not growing back in spring

    sapm    = 3.*0.001      ! m2/kg -->m2/g
    lfmsmn  = laimin/lapm
    stmsmn  = xsamin/sapm
! ---------------------------------------------------------------------------------

! respiration

     if(igs .eq. 0.) then
       rf = 0.5
     else
       rf = 1.0
     endif
            
     fnf     = min( foln/max(1.e-06,parameters%folnmx), 1.0 )
     tf      = parameters%arm**( (tv-298.16)/10. )
     resp    = parameters%rmf25 * tf * fnf * xlai * rf * (1.-wstres) ! umol/m2/s
     rsleaf  = min((lfmass-lfmsmn)/dt,resp*12.e-6)                         ! g/m2/s
     
     rsroot  = parameters%rmr25*(rtmass*1e-3)*tf *rf* 12.e-6         ! g/m2/s
     rsstem  = parameters%rms25*((stmass-stmsmn)*1e-3)*tf *rf* 12.e-6         ! g/m2/s
     rswood  = rswoodc * r(tv) * wood*parameters%wdpool

! carbon assimilation
! 1 mole -> 12 g carbon or 44 g co2; 1 umol -> 12.e-6 g carbon;

     carbfx  = psn * 12.e-6              ! umol co2 /m2/ s -> g/m2/s carbon

! fraction of carbon into leaf versus nonleaf

     leafpt = exp(0.01*(1.-exp(0.75*xlai))*xlai)
     if(vegtyp == parameters%eblforest) leafpt = exp(0.01*(1.-exp(0.50*xlai))*xlai)

     nonlef = 1.0 - leafpt
     stempt = xlai/10.0*leafpt
     leafpt = leafpt - stempt

!  fraction of carbon into wood versus root

     if(wood.gt.0) then
        woodf = (1.-exp(-bf*(parameters%wrrat*rtmass/wood))/bf)*parameters%wdpool
     else
        woodf = 0.
     endif

     rootpt = nonlef*(1.-woodf)
     woodpt = nonlef*woodf

! leaf and root turnover per time step

     lftovr = parameters%ltovrc*5.e-7*lfmass
     sttovr = parameters%ltovrc*5.e-7*stmass
     rttovr = rtovrc*rtmass
     wdtovr = 9.5e-10*wood

! seasonal leaf die rate dependent on temp and water stress
! water stress is set to 1 at permanent wilting point

     sc  = exp(-0.3*max(0.,tv-parameters%tdlef)) * (lfmass/120.) 
     sd  = exp((wstres-1.)*wstrc)
     dielf = lfmass*1.e-6*(parameters%dilefw * sd + parameters%dilefc*sc)
     diest = stmass*1.e-6*(parameters%dilefw * sd + parameters%dilefc*sc)

! calculate growth respiration for leaf, rtmass and wood

     grleaf = max(0.0,parameters%fragr*(leafpt*carbfx - rsleaf))
     grstem = max(0.0,parameters%fragr*(stempt*carbfx - rsstem))
     grroot = max(0.0,parameters%fragr*(rootpt*carbfx - rsroot))
     grwood = max(0.0,parameters%fragr*(woodpt*carbfx - rswood))

! impose lower t limit for photosynthesis

     addnpplf = max(0.,leafpt*carbfx - grleaf-rsleaf)
     addnppst = max(0.,stempt*carbfx - grstem-rsstem)
!     addnpplf = leafpt*carbfx - grleaf-rsleaf  ! mb: test kjetil 
!     addnppst = stempt*carbfx - grstem-rsstem  ! mb: test kjetil 
     if(tv.lt.parameters%tmin) addnpplf =0.
     if(tv.lt.parameters%tmin) addnppst =0.

! update leaf, root, and wood carbon
! avoid reducing leaf mass below its minimum value but conserve mass

     lfdel = (lfmass - lfmsmn)/dt
     stdel = (stmass - stmsmn)/dt
     dielf = min(dielf,lfdel+addnpplf-lftovr)
     diest = min(diest,stdel+addnppst-sttovr)

! net primary productivities

     nppl   = max(addnpplf,-lfdel)
     npps   = max(addnppst,-stdel)
     nppr   = rootpt*carbfx - rsroot - grroot
     nppw   = woodpt*carbfx - rswood - grwood

! masses of plant components

     lfmass = lfmass + (nppl-lftovr-dielf)*dt
     stmass = stmass + (npps-sttovr-diest)*dt   ! g/m2
     rtmass = rtmass + (nppr-rttovr)      *dt

     if(rtmass.lt.0.0) then
           rttovr = nppr
           rtmass = 0.0
     endif
     wood = (wood+(nppw-wdtovr)*dt)*parameters%wdpool

! soil carbon budgets

     fastcp = fastcp + (rttovr+lftovr+sttovr+wdtovr+dielf+diest)*dt  ! mb: add diest v3.7

     fst = 2.0**( (stc(1)-283.16)/10. )
     fsw = wroot / (0.20+wroot) * 0.23 / (0.23+wroot)
     rssoil = fsw * fst * parameters%mrp* max(0.,fastcp*1.e-3)*12.e-6

     stablc = 0.1*rssoil
     fastcp = fastcp - (rssoil + stablc)*dt
     stblcp = stblcp + stablc*dt

!  total carbon flux

     cflux  = - carbfx + rsleaf + rsroot + rswood + rsstem &  ! mb: add rsstem,grstem,0.9*rssoil v3.7
          + 0.9*rssoil + grleaf + grroot + grwood + grstem    ! g/m2/s

! for outputs

     gpp    = carbfx                                             !g/m2/s c
     npp    = nppl + nppw + nppr +npps                           !g/m2/s c
     autors = rsroot + rswood  + rsleaf + rsstem + &             !g/m2/s c  mb: add rsstem, grstem v3.7
              grleaf + grroot + grwood + grstem                  !g/m2/s c  mb: add 0.9* v3.7
     heters = 0.9*rssoil                                         !g/m2/s c
     nee    = (autors + heters - gpp)*44./12.                    !g/m2/s co2
     totsc  = fastcp + stblcp                                    !g/m2   c
     totlb  = lfmass + rtmass +stmass + wood                     !g/m2   c  mb: add stmass v3.7

! leaf area index and stem area index

     xlai    = max(lfmass*lapm,laimin)
     xsai    = max(stmass*sapm,xsamin)
    
  end subroutine co2flux

!== begin bvocflux =================================================================================

!  subroutine bvocflux(parameters,vocflx,  vegtyp,  vegfrac,  apar,   tv )
!
! ------------------------------------------------------------------------------------------
!      implicit none
! ------------------------------------------------------------------------------------------
!
! ------------------------ code history ---------------------------
! source file:       bvoc
! purpose:           bvoc emissions
! description:
! volatile organic compound emission 
! this code simulates volatile organic compound emissions
! following the algorithm presented in guenther, a., 1999: modeling
! biogenic volatile organic compound emissions to the atmosphere. in
! reactive hydrocarbons in the atmosphere, ch. 3
! this model relies on the assumption that 90% of isoprene and monoterpene
! emissions originate from canopy foliage:
!    e = epsilon * gamma * density * delta
! the factor delta (longterm activity factor) applies to isoprene emission
! from deciduous plants only. we neglect this factor at the present time.
! this factor is discussed in guenther (1997).
! subroutine written to operate at the patch level.
! in final implementation, remember:
! 1. may wish to call this routine only as freq. as rad. calculations
! 2. may wish to place epsilon values directly in pft-physiology file
! ------------------------ input/output variables -----------------
! input
!  integer                     ,intent(in) :: vegtyp  !vegetation type 
!  real                        ,intent(in) :: vegfrac !green vegetation fraction [0.0-1.0]
!  real                        ,intent(in) :: apar    !photosynthesis active energy by canopy (w/m2)
!  real                        ,intent(in) :: tv      !vegetation canopy temperature (k)
!
! output
!  real                        ,intent(out) :: vocflx(5) ! voc fluxes [ug c m-2 h-1]
!
! local variables
!
!  real, parameter :: r      = 8.314    ! univ. gas constant [j k-1 mol-1]
!  real, parameter :: alpha  = 0.0027   ! empirical coefficient
!  real, parameter :: cl1    = 1.066    ! empirical coefficient
!  real, parameter :: ct1    = 95000.0  ! empirical coefficient [j mol-1]
!  real, parameter :: ct2    = 230000.0 ! empirical coefficient [j mol-1]
!  real, parameter :: ct3    = 0.961    ! empirical coefficient
!  real, parameter :: tm     = 314.0    ! empirical coefficient [k]
!  real, parameter :: tstd   = 303.0    ! std temperature [k]
!  real, parameter :: bet    = 0.09     ! beta empirical coefficient [k-1]
!
!  integer ivoc        ! do-loop index
!  integer ityp        ! do-loop index
!  real epsilon(5)
!  real gamma(5)
!  real density
!  real elai
!  real par,cl,reciprod,ct
!
! epsilon :
!
!    do ivoc = 1, 5
!    epsilon(ivoc) = parameters%eps(vegtyp,ivoc)
!    end do
!
! gamma : activity factor. units [dimensionless]
!
!      reciprod = 1. / (r * tv * tstd)
!      ct = exp(ct1 * (tv - tstd) * reciprod) / &
!           (ct3 + exp(ct2 * (tv - tm) * reciprod))
!
!      par = apar * 4.6 ! (multiply w/m2 by 4.6 to get umol/m2/s)
!      cl  = alpha * cl1 * par * (1. + alpha * alpha * par * par)**(-0.5)
!
!   gamma(1) = cl * ct ! for isoprenes
!
!   do ivoc = 2, 5
!   gamma(ivoc) = exp(bet * (tv - tstd))
!   end do
!
! foliage density
!
! transform vegfrac to lai      
!
!   elai    = max(0.0,-6.5/2.5*alog((1.-vegfrac)))
!   density = elai / (parameters%slarea(vegtyp) * 0.5)
!
! calculate the voc flux
!
!   do ivoc = 1, 5
!   vocflx(ivoc) = epsilon(ivoc) * gamma(ivoc) * density
!   end do
!
!   end subroutine bvocflux
! ==================================================================================================
! ********************************* end of carbon subroutines *****************************
! ==================================================================================================

!== begin noahmp_options ===========================================================================

!>\ingroup NoahMP_LSM
  subroutine noahmp_options(idveg     ,iopt_crs  ,iopt_btr  ,iopt_run  ,iopt_sfc  ,iopt_frz , & 
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
  
  end subroutine noahmp_options
 

end module module_sf_noahmplsm


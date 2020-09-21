!>\file micro_mg3_0.F90
!! This file contains Morrison-Gettelman MP version 3.0 -
!! Update of MG microphysics with prognostic hail OR graupel.

!>\ingroup mg2mg3
!>\defgroup mg3_mp Morrison-Gettelman MP version 3.0
!> @{
!!---------------------------------------------------------------------------------
!! Purpose:
!!   MG microphysics version 3.0 - Update of MG microphysics with
!!                                 prognostic hail OR graupel.
!!
!! \authors Andrew Gettelman, Hugh Morrison
!!
!! \version 3 history: Sep 2016: development begun for hail, graupel
!!  This version:https://svn-ccsm-models.cgd.ucar.edu/cam1/branch_tags/mg3_tags/mg3_33_cam5_4_153/
!!
!! \version 2 history: Sep 2011: Development begun.
!!\n                    Feb 2013: Added of prognostic precipitation.
!!\n                    Aug 2015: Published and released version
!!
!! Contributions from:  Sean Santos, Peter Caldwell, Xiaohong Liu and Steve Ghan
!!
!! - Anning Cheng adopted mg2 for FV3GFS 9/29/2017
!!\n             add GMAO ice conversion and Liu et. al liquid water
!!\n              conversion in 10/12/2017
!!
!! - Anning showed promising results for FV3GFS on 10/15/2017
!! - S. Moorthi - Oct/Nov 2017 - optimized the MG2 code
!! - S. Moorthi - Nov 2017     - made the sedimentation quasi-implicit
!! - S. Moorthi - Feb 2018     - updated to MG3 - modified graupel sedimentation
!!                             other modifications to eliminate blowup.
!! - S. Moorthi - Mar 2018     - fixed a few bugs and added option to run as MG2
!! - S. Moorthi - Oct,29,2018  - change nlb from nlev/3 to levels with p/ps < 0.05 (nlball)
!!
!! invoked in CAM by specifying -microphys=mg3
!!
!! References:
!!
!!           Gettelman, A. and H. Morrison, Advanced Two-Moment Microphysics for Global Models.
!!           Part I: Off line tests and comparisons with other schemes.
!!           J. Climate, 28, 1268-1287. doi: 10.1175/JCLI-D-14-00102.1, 2015.
!!
!!           Gettelman, A., H. Morrison, S. Santos, P. Bogenschutz and P. H. Caldwell
!!           Advanced Two-Moment Microphysics for Global Models.
!!           Part II: Global model solutions and Aerosol-Cloud Interactions.
!!           J. Climate, 28, 1288-1307. doi:10.1175/JCLI-D-14-00103.1 , 2015.
!!
!! for questions contact Hugh Morrison, Andrew Gettelman
!! e-mail: morrison@ucar.edu, andrew@ucar.edu
!!---------------------------------------------------------------------------------
!!
!! NOTE: Modified to allow other microphysics packages (e.g. CARMA) to do ice
!! microphysics in cooperation with the MG liquid microphysics. This is
!! controlled by the do_cldice variable.
!!
!! If do_cldice is false, then MG microphysics should not update CLDICE or
!! NUMICE; it is assumed that the other microphysics scheme will have updated
!! CLDICE and NUMICE. The other microphysics should handle the following
!! processes that would have been done by MG:
!!   - Detrainment (liquid and ice)
!!   - Homogeneous ice nucleation
!!   - Heterogeneous ice nucleation
!!   - Bergeron process
!!   - Melting of ice
!!   - Freezing of cloud drops
!!   - Autoconversion (ice -> snow)
!!   - Growth/Sublimation of ice
!!   - Sedimentation of ice
!!
!! This option has not been updated since the introduction of prognostic
!! precipitation, and probably should be adjusted to cover snow as well.
!
!---------------------------------------------------------------------------------
!Version 3.O based on micro_mg2_0.F90 and WRF3.8.1 module_mp_morr_two_moment.F
!---------------------------------------------------------------------------------
! Based on micro_mg (restructuring of former cldwat2m_micro)
! Author: Andrew Gettelman, Hugh Morrison.
! Contributions from: Xiaohong Liu and Steve Ghan
! December 2005-May 2010
! Description in: Morrison and Gettelman, 2008. J. Climate (MG2008)
!                 Gettelman et al., 2010 J. Geophys. Res. - Atmospheres (G2010)
! for questions contact Hugh Morrison, Andrew Gettelman
! e-mail: morrison@ucar.edu, andrew@ucar.edu
!---------------------------------------------------------------------------------
! Code comments added by HM, 093011
! General code structure:
!
! Code is divided into two main subroutines:
!   subroutine micro_mg_init --> initializes microphysics routine, should be called
!                                  once at start of simulation
!   subroutine micro_mg_tend --> main microphysics routine to be called each time step
!                                this also calls several smaller subroutines to calculate
!                                microphysical processes and other utilities
!
! List of external functions:
!   qsat_water --> for calculating saturation vapor pressure with respect to liquid water
!   qsat_ice --> for calculating saturation vapor pressure with respect to ice
!   gamma   --> standard mathematical gamma function
! .........................................................................
! List of inputs through use statement in fortran90:
! Variable Name                      Description                Units
! .........................................................................
! gravit          acceleration due to gravity                    m s-2
! rair            dry air gas constant for air                  J kg-1 K-1
! tmelt           temperature of melting point for water          K
! cpair           specific heat at constant pressure for dry air J kg-1 K-1
! rh2o            gas constant for water vapor                  J kg-1 K-1
! latvap          latent heat of vaporization                   J kg-1
! latice          latent heat of fusion                         J kg-1
! qsat_water      external function for calculating liquid water
!                 saturation vapor pressure/humidity              -
! qsat_ice        external function for calculating ice
!                 saturation vapor pressure/humidity              pa
! rhmini          relative humidity threshold parameter for
!                 nucleating ice                                  -
! .........................................................................
! NOTE: List of all inputs/outputs passed through the call/subroutine statement
!       for micro_mg_tend is given below at the start of subroutine micro_mg_tend.
!---------------------------------------------------------------------------------

! Procedures required:
! 1) An implementation of the gamma function (if not intrinsic).
! 2) saturation vapor pressure and specific humidity over water
! 3) svp over ice

module micro_mg3_0

use machine,        only : r8    => kind_phys
use funcphys,       only : fpvsl, fpvsi

!use wv_sat_methods, only: &
!    qsat_water => wv_sat_qsat_water, &
!    qsat_ice => wv_sat_qsat_ice

! Parameters from the utilities module.
use micro_mg_utils, only : pi,   omsm,  qsmall, mincld, rhosn, rhoi, &
                           rhow, rhows, ac,     bc,     ai,    bi,   &
                           aj,   bj,    ar,     br,     as,    bs,   &
!++ag
                           ag,   bg,    ah,     bh,     rhog,  rhoh, &
!--ag
                           mi0,  rising_factorial

implicit none
private
save

public :: micro_mg_init, micro_mg_tend, qcvar

! Switches for specification rather than prediction of droplet and crystal number
! note: number will be adjusted as needed to keep mean size within bounds,
! even when specified droplet or ice number is used
!
! If constant cloud ice number is set (nicons = .true.),
! then all microphysical processes except mass transfer due to ice nucleation
! (mnuccd) are based on the fixed cloud ice number. Calculation of
! mnuccd follows from the prognosed ice crystal number ni.

logical :: nccons  !< nccons = .true. to specify constant cloud droplet number
logical :: nicons  !< nicons = .true. to specify constant cloud ice number
!++ag kt
logical :: ngcons  !< ngcons = .true. to specify constant graupel number
!--ag kt

! specified ice and droplet number concentrations
! note: these are local in-cloud values, not grid-mean
real(r8) :: ncnst  !< droplet num concentration when nccons=.true. (m-3)
real(r8) :: ninst  !< ice num concentration when nicons=.true. (m-3)
!++ag kt
real(r8) :: ngnst  !< graupel num concentration when ngcons=.true. (m-3)
!--ag kt

!=========================================================
! Private module parameters
!=========================================================

!> Range of cloudsat reflectivities (dBz) for analytic simulator
real(r8), parameter :: csmin = -30._r8
real(r8), parameter :: csmax = 26._r8
real(r8), parameter :: mindbz = -99._r8
real(r8), parameter :: minrefl = 1.26e-10_r8    ! minrefl = 10._r8**(mindbz/10._r8)

! autoconversion size threshold for cloud ice to snow (m)
real(r8)            :: dcs, ts_au, ts_au_min, qcvar

! minimum mass of new crystal due to freezing of cloud droplets done
! externally (kg)
real(r8), parameter :: mi0l_min = 4._r8/3._r8*pi*rhow*(4.e-6_r8)**3

! Ice number sublimation parameter. Assume some decrease in ice number with sublimation if non-zero. Else, no decrease in number with sublimation.
real(r8), parameter :: sublim_factor = 0.0_r8      !number sublimation factor.

real(r8), parameter :: zero=0.0_r8, one=1.0_r8,    two=2.0_r8,  three=3.0_r8, &
                       four=4.0_r8, five=5.0_r8,   six=6._r8,   half=0.5_r8,  &
                       ten=10.0_r8, forty=40.0_r8, oneo6=one/six

!=========================================================
! Constants set in initialization
!=========================================================

! Set using arguments to micro_mg_init
real(r8) :: g           !< gravity
real(r8) :: r           !< dry air gas constant
real(r8) :: rv          !< water vapor gas constant
real(r8) :: cpp         !< specific heat of dry air
real(r8) :: tmelt       !< freezing point of water (K)

! latent heats of:
real(r8) :: xxlv        !< vaporization
real(r8) :: xlf         !< freezing
real(r8) :: xxls        !< sublimation

real(r8) :: rhmini      !< Minimum rh for ice cloud fraction > 0.

! flags
logical  :: microp_uniform, do_cldice, use_hetfrz_classnuc,  &
!++ag
            do_hail,        do_graupel
!--ag

real(r8) :: rhosu       !< typical 850mn air density

real(r8) :: icenuct     !< ice nucleation temperature: currently -5 degrees C

real(r8) :: snowmelt    !< what temp to melt all snow: currently 2 degrees C
real(r8) :: rainfrze    !< what temp to freeze all rain: currently -5 degrees C

real(r8) :: rhogtmp     !< hail or graupel density (kg m-3)
real(r8) :: agtmp       !< tmp ag/ah parameter
real(r8) :: bgtmp       !< tmp fall speed parameter

! additional constants to help speed up code
real(r8) :: gamma_br_plus1, gamma_bs_plus1, gamma_bi_plus1, gamma_bj_plus1, gamma_bg_plus1
real(r8) :: gamma_br_plus4, gamma_bs_plus4, gamma_bi_plus4, gamma_bj_plus4, gamma_bg_plus4
real(r8) :: xxlv_squared,   xxls_squared
real(r8) :: omeps, epsqs

character(len=16)  :: micro_mg_precip_frac_method  !< type of precipitation fraction method
real(r8)           :: micro_mg_berg_eff_factor     !< berg efficiency factor

logical  :: allow_sed_supersat !< Allow supersaturated conditions after sedimentation loop
logical  :: do_sb_physics      !< do SB 2001 autoconversion or accretion physics
logical  :: do_ice_gmao
logical  :: do_liq_liu

!===============================================================================
contains
!===============================================================================

!>\ingroup mg3_mp
!! This subroutine initializes the microphysics
!! and needs to be called once at start of simulation.
!!\author Andrew Gettelman, Dec 2005
subroutine micro_mg_init(                                         &
     kind, gravit, rair, rh2o, cpair, eps,                        &
     tmelt_in, latvap, latice,                                    &
     rhmini_in, micro_mg_dcs,ts_auto, mg_qcvar,                   &
!++ag
     micro_mg_do_hail_in, micro_mg_do_graupel_in,                 &
!--ag
     microp_uniform_in, do_cldice_in, use_hetfrz_classnuc_in,     &
     micro_mg_precip_frac_method_in, micro_mg_berg_eff_factor_in, &
     allow_sed_supersat_in, do_sb_physics_in,                     &
     do_ice_gmao_in, do_liq_liu_in,                               &
     nccons_in, nicons_in, ncnst_in, ninst_in, ngcons_in, ngnst_in)
!    nccons_in, nicons_in, ncnst_in, ninst_in, ngcons_in, ngnst_in, errstring)

  use micro_mg_utils, only : micro_mg_utils_init
  use wv_saturation,  only : gestbl

  !-----------------------------------------------------------------------
  !
  ! Purpose:
  ! initialize constants for MG microphysics
  !
  ! Author: Andrew Gettelman Dec 2005
  !
  !-----------------------------------------------------------------------

  integer,  intent(in)  :: kind         ! Kind used for reals
  real(r8), intent(in)  :: gravit
  real(r8), intent(in)  :: rair
  real(r8), intent(in)  :: rh2o
  real(r8), intent(in)  :: cpair
  real(r8), intent(in)  :: eps
! real(r8), intent(in)  :: fv
  real(r8), intent(in)  :: tmelt_in     ! Freezing point of water (K)
  real(r8), intent(in)  :: latvap
  real(r8), intent(in)  :: latice
  real(r8), intent(in)  :: rhmini_in    ! Minimum rh for ice cloud fraction > 0.
  real(r8), intent(in)  :: micro_mg_dcs
  real(r8), intent(in)  :: ts_auto(2)
  real(r8), intent(in)  :: mg_qcvar

!++ag
!MG3 dense precipitating ice. Note, only 1 can be true, or both false.
  logical,  intent(in)  :: micro_mg_do_graupel_in  ! .true.  = configure with graupel
                                                   ! .false. = no graupel (hail possible)
  logical,  intent(in)  :: micro_mg_do_hail_in     ! .true.  = configure with hail
                                                   ! .false. = no hail (graupel possible)
!--ag

  logical,  intent(in)  :: microp_uniform_in ! .true.  = configure uniform for sub-columns
                                             ! .false. = use w/o sub-columns (standard)
  logical,  intent(in)  :: do_cldice_in      ! .true.  = do all processes (standard)
                                             ! .false. = skip all processes affecting cloud ice
  logical,  intent(in)  :: use_hetfrz_classnuc_in ! use heterogeneous freezing

  character(len=16),intent(in)  :: micro_mg_precip_frac_method_in  ! type of precipitation fraction method
  real(r8),         intent(in)  :: micro_mg_berg_eff_factor_in     ! berg efficiency factor
  logical,  intent(in)  ::  allow_sed_supersat_in ! allow supersaturated conditions after sedimentation loop
  logical,  intent(in)  ::  do_sb_physics_in ! do SB autoconversion and accretion physics
  logical,  intent(in)  ::  do_ice_gmao_in
  logical,  intent(in)  ::  do_liq_liu_in

  logical,  intent(in)  :: nccons_in, nicons_in, ngcons_in
  real(r8), intent(in)  :: ncnst_in,  ninst_in,  ngnst_in
  logical ip
  real(r8):: tmn, tmx, trice


! character(128), intent(out) :: errstring    ! Output status (non-blank for error return)

  !-----------------------------------------------------------------------

  dcs       = micro_mg_dcs * 1.0e-6_r8
  ts_au_min = ts_auto(1)
  ts_au     = ts_auto(2)
  qcvar     = mg_qcvar

  ! Initialize subordinate utilities module.
  call micro_mg_utils_init(kind, rair, rh2o, cpair, tmelt_in, latvap, latice, &
                           dcs)
!                          dcs, errstring)

! if (trim(errstring) /= "") return

  ! declarations for MG code (transforms variable names)

  g      = gravit              ! gravity
  r      = rair                ! dry air gas constant: note units(phys_constants are in J/K/kmol)
! write(0,*)' in micro_mg_utils_init=',' r=',r,' rair=',rair,' rh2o=',rh2o
  rv     = rh2o                ! water vapor gas constant
  cpp    = cpair               ! specific heat of dry air
  tmelt  = tmelt_in
  rhmini = rhmini_in
  micro_mg_precip_frac_method = micro_mg_precip_frac_method_in
  micro_mg_berg_eff_factor    = micro_mg_berg_eff_factor_in
  allow_sed_supersat          = allow_sed_supersat_in
  do_sb_physics               = do_sb_physics_in
  do_ice_gmao                 = do_ice_gmao_in
  do_liq_liu                  = do_liq_liu_in

  nccons = nccons_in
  nicons = nicons_in
  ncnst  = ncnst_in
  ninst  = ninst_in
!++ag
  ngcons = ngcons_in
  ngnst  = ngnst_in
!--ag

  ! latent heats

  xxlv = latvap         ! latent heat vaporization
  xlf  = latice         ! latent heat freezing
  xxls = xxlv + xlf     ! latent heat of sublimation

  ! flags
  microp_uniform      = microp_uniform_in
  do_cldice           = do_cldice_in
  use_hetfrz_classnuc = use_hetfrz_classnuc_in
!++ag
  do_hail             = micro_mg_do_hail_in
  do_graupel          = micro_mg_do_graupel_in
!
  if (do_hail) then
    agtmp   = ah
    bgtmp   = bh
    rhogtmp = rhoh
  elseif (do_graupel) then
    agtmp   = ag
    bgtmp   = bg
    rhogtmp = rhog
  else
    agtmp = zero
    bgtmp = zero
  endif
!--ag

  ! typical air density at 850 mb

  rhosu = 85000._r8 / (rair * tmelt)

  ! Maximum temperature at which snow is allowed to exist
  snowmelt = tmelt + two
  ! Minimum temperature at which rain is allowed to exist
  rainfrze = tmelt - forty

  ! Ice nucleation temperature
  icenuct  = tmelt - five

  ! Define constants to help speed up code (this limits calls to gamma function)
  gamma_br_plus1 = gamma(br+one)
  gamma_br_plus4 = gamma(br+four)
  gamma_bs_plus1 = gamma(bs+one)
  gamma_bs_plus4 = gamma(bs+four)
  gamma_bi_plus1 = gamma(bi+one)
  gamma_bi_plus4 = gamma(bi+four)
  gamma_bj_plus1 = gamma(bj+one)
  gamma_bj_plus4 = gamma(bj+four)
!
  gamma_bg_plus1 = gamma(bgtmp+one)
  gamma_bg_plus4 = gamma(bgtmp+four)

  xxlv_squared = xxlv * xxlv
  xxls_squared = xxls * xxls
  epsqs = eps
  omeps = one - epsqs
  tmn   = 173.16_r8
  tmx   = 375.16_r8
  trice = 35.00_r8
  ip    = .true.
!> - call gestbl()
  call gestbl(tmn ,tmx ,trice ,ip ,epsqs , latvap ,latice ,rh2o ,     &
              cpair ,tmelt_in )



end subroutine micro_mg_init

!===============================================================================
!microphysics routine for each timestep goes here...

!>\ingroup mg3_mp
!! This subroutine calculates the MG3 microphysical processes.
!>\authors Hugh Morrison, Andrew Gettelman, NCAR, Peter Caldwell, LLNL
!! e-mail: morrison@ucar.edu, andrew@ucar.edu
!!\section mg3_micro_mg_tend MG3 micro_mg_tend General Algorithm
!> @{
subroutine micro_mg_tend (                                       &
     mgncol,             nlev,               deltatin,           &
     t,                            q,                            &
     qcn,                          qin,                          &
     ncn,                          nin,                          &
     qrn,                          qsn,                          &
     nrn,                          nsn,                          &
!++ag
     qgr,                          ngr,                          &
!--ag
     relvar,                       accre_enhan_i,                &
     p,                            pdel,                         &
     cldn,    liqcldf,        icecldf,       qsatfac,            &
     qcsinksum_rate1ord,                                         &
     naai,                         npccnin,                      &
     rndst,                        nacon,                        &
     tlat,                         qvlat,                        &
     qctend,                       qitend,                       &
     nctend,                       nitend,                       &
     qrtend,                       qstend,                       &
     nrtend,                       nstend,                       &
!++ag
     qgtend,                       ngtend,                       &
!--ag
     effc,               effc_fn,            effi,               &
     sadice,                       sadsnow,                      &
     prect,                        preci,                        &
     nevapr,                       evapsnow,                     &
     am_evp_st,                                                  &
     prain,                        prodsnow,                     &
     cmeout,                       deffi,                        &
     pgamrad,                      lamcrad,                      &
     qsout,                        dsout,                        &
!++ag
     qgout,     ngout,             dgout,                        &
!--ag
     lflx,               iflx,                                   &
!++ag
     gflx,                                                       &
!--ag
     rflx,               sflx,               qrout,              &
!++ag
     reff_rain,          reff_snow,          reff_grau,          &
!--ag

     qcsevap,            qisevap,            qvres,              &
     cmeitot,            vtrmc,              vtrmi,              &
     umr,                          ums,                          &
!++ag
     umg,                qgsedten,                               &
!--ag
     qcsedten,                     qisedten,                     &
     qrsedten,                     qssedten,                     &
     pratot,                       prctot,                       &
     mnuccctot,          mnuccttot,          msacwitot,          &
     psacwstot,          bergstot,           bergtot,            &
     melttot,                      homotot,                      &
     qcrestot,           prcitot,            praitot,            &
!++ag
     qirestot,           mnuccrtot, mnuccritot, pracstot,        &
!--ag
     meltsdttot,         frzrdttot,          mnuccdtot,          &
!++ag
     pracgtot,           psacwgtot,          pgsacwtot,          &
     pgracstot,          prdgtot,                                &
     qmultgtot,          qmultrgtot,         psacrtot,           &
     npracgtot,          nscngtot,           ngracstot,          &
     nmultgtot,          nmultrgtot,         npsacwgtot,         &
!--ag
     nrout,                        nsout,                        &
     refl,               arefl,              areflz,             &
     frefl,              csrfl,              acsrfl,             &
     fcsrfl,                       rercld,                       &
     ncai,                         ncal,                         &
     qrout2,                       qsout2,                       &
     nrout2,                       nsout2,                       &
     drout2,                       dsout2,                       &
!++ag
     qgout2,             ngout2,   dgout2,   freqg,              &
!--ag
     freqs,                        freqr,                        &
     nfice,                        qcrat,                        &
     prer_evap, xlat, xlon, lprnt, iccn, nlball)

  ! Constituent properties.
  use micro_mg_utils, only: mg_liq_props,    &
                            mg_ice_props,    &
                            mg_rain_props,   &
!++ag
                            mg_graupel_props,&
!--ag
                            mg_snow_props

  ! Size calculation functions.
  use micro_mg_utils, only: size_dist_param_liq,   &
                            size_dist_param_basic, &
                            avg_diameter

  ! Microphysical processes.
  use micro_mg_utils, only: ice_deposition_sublimation,        &
                            sb2001v2_liq_autoconversion,       &
                            sb2001v2_accre_cld_water_rain,     &       
                            kk2000_liq_autoconversion,         &
                            ice_autoconversion,                &
                            immersion_freezing,                &
                            contact_freezing,                  &
                            snow_self_aggregation,             &
                            accrete_cloud_water_snow,          &
                            secondary_ice_production,          &
                            accrete_rain_snow,                 &
                            heterogeneous_rain_freezing,       &
                            accrete_cloud_water_rain,          &
                            self_collection_rain,              &
                            accrete_cloud_ice_snow,            &
                            evaporate_sublimate_precip,        &
                            bergeron_process_snow,             &
                            size_dist_param_ice,               &
!++ag
                            graupel_collecting_snow,           &
                            graupel_collecting_rain,           &
                            graupel_collecting_cld_water,      &
                            graupel_riming_liquid_snow,        &
                            graupel_rain_riming_snow,          &
                            graupel_rime_splintering,          &
                            evaporate_sublimate_precip_graupel,&
!                           graupel_sublimate_evap
!--ag
                            liu_liq_autoconversion,            &
                            gmao_ice_autoconversion

  !Authors: Hugh Morrison, Andrew Gettelman, NCAR, Peter Caldwell, LLNL
  ! e-mail: morrison@ucar.edu, andrew@ucar.edu

  ! input arguments
  integer,  intent(in) :: mgncol                  !< number of microphysics columns
  integer,  intent(in) :: nlev                    !< number of layers
  integer,  intent(in) :: nlball(mgncol)          !< sedimentation start level
  real(r8), intent(in) :: xlat,xlon               !< number of layers
  real(r8), intent(in) :: deltatin                !< time step (s)
  real(r8), intent(in) :: t(mgncol,nlev)          !< input temperature (K)
  real(r8), intent(in) :: q(mgncol,nlev)          !< input h20 vapor mixing ratio (kg/kg)

  ! note: all input cloud variables are grid-averaged
  real(r8), intent(in) :: qcn(mgncol,nlev)        !< cloud water mixing ratio (kg/kg)
  real(r8), intent(in) :: qin(mgncol,nlev)        !< cloud ice mixing ratio (kg/kg)
  real(r8), intent(in) :: ncn(mgncol,nlev)        !< cloud water number conc (1/kg)
  real(r8), intent(in) :: nin(mgncol,nlev)        !< cloud ice number conc (1/kg)

  real(r8), intent(in) :: qrn(mgncol,nlev)        !< rain mixing ratio (kg/kg)
  real(r8), intent(in) :: qsn(mgncol,nlev)        !< snow mixing ratio (kg/kg)
  real(r8), intent(in) :: nrn(mgncol,nlev)        !< rain number conc (1/kg)
  real(r8), intent(in) :: nsn(mgncol,nlev)        !< snow number conc (1/kg)
!++ag
  real(r8), intent(in) :: qgr(mgncol,nlev)        !< graupel/hail mixing ratio (kg/kg)
  real(r8), intent(in) :: ngr(mgncol,nlev)        !< graupel/hail number conc (1/kg)
!--ag

  real(r8)             :: relvar(mgncol,nlev)     !< cloud water relative variance (-)
  real(r8)             :: accre_enhan(mgncol,nlev)!< optional accretion
! real(r8), intent(in) :: relvar_i                !< cloud water relative variance (-)
  real(r8), intent(in) :: accre_enhan_i           !< optional accretion
                                                  !< enhancement factor (-)

  real(r8), intent(in) :: p(mgncol,nlev)          !< air pressure (pa)
  real(r8), intent(in) :: pdel(mgncol,nlev)       !< pressure difference across level (pa)

  real(r8), intent(in) :: cldn(mgncol,nlev)       !< cloud fraction (no units)
  real(r8), intent(in) :: liqcldf(mgncol,nlev)    !< liquid cloud fraction (no units)
  real(r8), intent(in) :: icecldf(mgncol,nlev)    !< ice cloud fraction (no units)
  real(r8), intent(in) :: qsatfac(mgncol,nlev)    !< subgrid cloud water saturation scaling factor (no units)
  logical, intent(in)  :: lprnt                   !< control flag for diagnostic print out 
  integer, intent(in)  :: iccn                    !< flag for IN and CCN forcing for Morrison-Gettelman microphysics 


  ! used for scavenging
  ! Inputs for aerosol activation
  real(r8), intent(inout) :: naai(mgncol,nlev)    !< ice nucleation number (from microp_aero_ts) (1/kg)
  real(r8), intent(in)    :: npccnin(mgncol,nlev) !< ccn activated number tendency (from microp_aero_ts) (1/kg*s)
! real(r8), intent(in)    :: npccn(mgncol,nlev)   !< ccn activated number tendency (from microp_aero_ts) (1/kg*s)
  real(r8)                :: npccn(mgncol,nlev)   !< ccn activated number tendency (from microp_aero_ts) (1/kg*s)

  ! Note that for these variables, the dust bin is assumed to be the last index.
  ! (For example, in CAM, the last dimension is always size 4.)
  real(r8), intent(in) :: rndst(mgncol,nlev,10)   !< radius of each dust bin, for contact freezing (from microp_aero_ts) (m)
  real(r8), intent(in) :: nacon(mgncol,nlev,10)   !< number in each dust bin, for contact freezing  (from microp_aero_ts) (1/m^3)
  
  ! output arguments

  real(r8), intent(out) :: qcsinksum_rate1ord(mgncol,nlev) !< 1st order rate for
                                                           !! direct cw to precip conversion
  real(r8), intent(out) :: tlat(mgncol,nlev)      !< latent heating rate       (W/kg)
  real(r8), intent(out) :: qvlat(mgncol,nlev)     !< microphysical tendency qv (1/s)
  real(r8), intent(out) :: qctend(mgncol,nlev)    !< microphysical tendency qc (1/s)
  real(r8), intent(out) :: qitend(mgncol,nlev)    !< microphysical tendency qi (1/s)
  real(r8), intent(out) :: nctend(mgncol,nlev)    !< microphysical tendency nc (1/(kg*s))
  real(r8), intent(out) :: nitend(mgncol,nlev)    !< microphysical tendency ni (1/(kg*s))

  real(r8), intent(out) :: qrtend(mgncol,nlev)    !< microphysical tendency qr (1/s)
  real(r8), intent(out) :: qstend(mgncol,nlev)    !< microphysical tendency qs (1/s)
  real(r8), intent(out) :: nrtend(mgncol,nlev)    !< microphysical tendency nr (1/(kg*s))
  real(r8), intent(out) :: nstend(mgncol,nlev)    !< microphysical tendency ns (1/(kg*s))
!++ag
  real(r8), intent(out) :: qgtend(mgncol,nlev)    !< microphysical tendency qg (1/s)
  real(r8), intent(out) :: ngtend(mgncol,nlev)    !< microphysical tendency ng (1/(kg*s))
!--ag
  real(r8), intent(out) :: effc(mgncol,nlev)      !< droplet effective radius (micron)
  real(r8), intent(out) :: effc_fn(mgncol,nlev)   !< droplet effective radius, assuming nc = 1.e8 kg-1
  real(r8), intent(out) :: effi(mgncol,nlev)      !< cloud ice effective radius (micron)
  real(r8), intent(out) :: sadice(mgncol,nlev)    !< cloud ice surface area density (cm2/cm3)
  real(r8), intent(out) :: sadsnow(mgncol,nlev)   !< cloud snow surface area density (cm2/cm3)
  real(r8), intent(out) :: prect(mgncol)          !< surface precip rate (m/s)
  real(r8), intent(out) :: preci(mgncol)          !< cloud ice/snow precip rate (m/s)
  real(r8), intent(out) :: nevapr(mgncol,nlev)    !< evaporation rate of rain + snow (1/s)
  real(r8), intent(out) :: evapsnow(mgncol,nlev)  !< sublimation rate of snow (1/s)
  real(r8), intent(out) :: am_evp_st(mgncol,nlev) !< stratiform evaporation area (frac)
  real(r8), intent(out) :: prain(mgncol,nlev)     !< production of rain + snow (1/s)
  real(r8), intent(out) :: prodsnow(mgncol,nlev)  !< production of snow (1/s)
  real(r8), intent(out) :: cmeout(mgncol,nlev)    !< evap/sub of cloud (1/s)
  real(r8), intent(out) :: deffi(mgncol,nlev)     !< ice effective diameter for optics (radiation) (micron)
  real(r8), intent(out) :: pgamrad(mgncol,nlev)   !< ice gamma parameter for optics (radiation) (no units)
  real(r8), intent(out) :: lamcrad(mgncol,nlev)   !< slope of droplet distribution for optics (radiation) (1/m)
  real(r8), intent(out) :: qsout(mgncol,nlev)     !< snow mixing ratio (kg/kg)
  real(r8), intent(out) :: dsout(mgncol,nlev)     !< snow diameter (m)
  real(r8), intent(out) :: lflx(mgncol,2:nlev+1)  !< grid-box average liquid condensate flux (kg m^-2 s^-1)
  real(r8), intent(out) :: iflx(mgncol,2:nlev+1)  !< grid-box average ice condensate flux (kg m^-2 s^-1)
  real(r8), intent(out) :: rflx(mgncol,2:nlev+1)  !< grid-box average rain flux (kg m^-2 s^-1)
  real(r8), intent(out) :: sflx(mgncol,2:nlev+1)  !< grid-box average snow flux (kg m^-2 s^-1)
!++ag
  real(r8), intent(out) :: gflx(mgncol,2:nlev+1)  !< grid-box average graupel/hail flux (kg m^-2 s^-1)
!--ag
  real(r8), intent(out) :: qrout(mgncol,nlev)     !< grid-box average rain mixing ratio (kg/kg)
  real(r8), intent(out) :: reff_rain(mgncol,nlev) !< rain effective radius (micron)
  real(r8), intent(out) :: reff_snow(mgncol,nlev) !< snow effective radius (micron)
!++ag
  real(r8), intent(out) :: reff_grau(mgncol,nlev) !< graupel effective radius (micron)
!--ag
  real(r8), intent(out) :: qcsevap(mgncol,nlev)   !< cloud water evaporation due to sedimentation (1/s)
  real(r8), intent(out) :: qisevap(mgncol,nlev)   !< cloud ice sublimation due to sedimentation (1/s)
  real(r8), intent(out) :: qvres(mgncol,nlev)     !< residual condensation term to ensure RH < 100% (1/s)
  real(r8), intent(out) :: cmeitot(mgncol,nlev)   !< grid-mean cloud ice sub/dep (1/s)
  real(r8), intent(out) :: vtrmc(mgncol,nlev)     !< mass-weighted cloud water fallspeed (m/s)
  real(r8), intent(out) :: vtrmi(mgncol,nlev)     !< mass-weighted cloud ice fallspeed (m/s)
  real(r8), intent(out) :: umr(mgncol,nlev)       !< mass weighted rain fallspeed (m/s)
  real(r8), intent(out) :: ums(mgncol,nlev)       !< mass weighted snow fallspeed (m/s)
!++ag
  real(r8), intent(out) :: umg(mgncol,nlev)       !< mass weighted graupel/hail fallspeed (m/s)
  real(r8), intent(out) :: qgsedten(mgncol,nlev)  !< qg sedimentation tendency (1/s)
!--ag

  real(r8), intent(out) :: qcsedten(mgncol,nlev)  !< qc sedimentation tendency (1/s)
  real(r8), intent(out) :: qisedten(mgncol,nlev)  !< qi sedimentation tendency (1/s)
  real(r8), intent(out) :: qrsedten(mgncol,nlev)  !< qr sedimentation tendency (1/s)
  real(r8), intent(out) :: qssedten(mgncol,nlev)  !< qs sedimentation tendency (1/s)

  ! microphysical process rates for output (mixing ratio tendencies) (all have units of 1/s)
  real(r8), intent(out) :: pratot(mgncol,nlev)    !< accretion of cloud by rain
  real(r8), intent(out) :: prctot(mgncol,nlev)    !< autoconversion of cloud to rain
  real(r8), intent(out) :: mnuccctot(mgncol,nlev) !< mixing ratio tend due to immersion freezing
  real(r8), intent(out) :: mnuccttot(mgncol,nlev) !< mixing ratio tend due to contact freezing
  real(r8), intent(out) :: msacwitot(mgncol,nlev) !< mixing ratio tend due to H-M splintering
  real(r8), intent(out) :: psacwstot(mgncol,nlev) !< collection of cloud water by snow
  real(r8), intent(out) :: bergstot(mgncol,nlev)  !< bergeron process on snow
  real(r8), intent(out) :: bergtot(mgncol,nlev)   !< bergeron process on cloud ice
  real(r8), intent(out) :: melttot(mgncol,nlev)   !< melting of cloud ice
  real(r8), intent(out) :: homotot(mgncol,nlev)   !< homogeneous freezing cloud water
  real(r8), intent(out) :: qcrestot(mgncol,nlev)  !< residual cloud condensation due to removal of excess supersat
  real(r8), intent(out) :: prcitot(mgncol,nlev)   !< autoconversion of cloud ice to snow
  real(r8), intent(out) :: praitot(mgncol,nlev)   !< accretion of cloud ice by snow
  real(r8), intent(out) :: qirestot(mgncol,nlev)  !< residual ice deposition due to removal of excess supersat
  real(r8), intent(out) :: mnuccrtot(mgncol,nlev) !< mixing ratio tendency due to heterogeneous freezing of rain to snow (1/s)
  real(r8), intent(out) :: mnuccritot(mgncol,nlev)!< mixing ratio tendency due to heterogeneous freezing of rain to ice (1/s)
  real(r8), intent(out) :: pracstot(mgncol,nlev)  !< mixing ratio tendency due to accretion of rain by snow (1/s)
  real(r8), intent(out) :: meltsdttot(mgncol,nlev)!< latent heating rate due to melting of snow  (W/kg)
  real(r8), intent(out) :: frzrdttot(mgncol,nlev) !< latent heating rate due to homogeneous freezing of rain (W/kg)
  real(r8), intent(out) :: mnuccdtot(mgncol,nlev) !< mass tendency from ice nucleation
!++ag Hail/Graupel Tendencies
  real(r8), intent(out) :: pracgtot(mgncol,nlev)  !< change in q collection rain by graupel  (precipf)
  real(r8), intent(out) :: psacwgtot(mgncol,nlev) !< change in q collection droplets by graupel (lcldm)
  real(r8), intent(out) :: pgsacwtot(mgncol,nlev) !< conversion q to graupel due to collection droplets by snow  (lcldm)
  real(r8), intent(out) :: pgracstot(mgncol,nlev) !< conversion q to graupel due to collection rain by snow (precipf)
  real(r8), intent(out) :: prdgtot(mgncol,nlev)   !< dep of graupel (precipf)
! real(r8), intent(out) :: eprdgtot(mgncol,nlev)  !< sub of graupel (precipf)
  real(r8), intent(out) :: qmultgtot(mgncol,nlev) !< change q due to ice mult droplets/graupel  (lcldm)
  real(r8), intent(out) :: qmultrgtot(mgncol,nlev)!< change q due to ice mult rain/graupel (precipf)
  real(r8), intent(out) :: psacrtot(mgncol,nlev)  !< conversion due to coll of snow by rain (precipf)
  real(r8), intent(out) :: npracgtot(mgncol,nlev) !< change n collection rain by graupel  (precipf)
  real(r8), intent(out) :: nscngtot(mgncol,nlev)  !< change n conversion to graupel due to collection droplets by snow (lcldm)
  real(r8), intent(out) :: ngracstot(mgncol,nlev) !< change n conversion to graupel due to collection rain by snow (precipf)
  real(r8), intent(out) :: nmultgtot(mgncol,nlev) !< ice mult due to acc droplets by graupel  (lcldm)
  real(r8), intent(out) :: nmultrgtot(mgncol,nlev)!< ice mult due to acc rain by graupel  (precipf)
  real(r8), intent(out) :: npsacwgtot(mgncol,nlev)!< change n collection droplets by graupel (lcldm?)
!--ag
  real(r8), intent(out) :: nrout(mgncol,nlev)     !< rain number concentration (1/m3)
  real(r8), intent(out) :: nsout(mgncol,nlev)     !< snow number concentration (1/m3)
  real(r8), intent(out) :: refl(mgncol,nlev)      !< analytic radar reflectivity
  real(r8), intent(out) :: arefl(mgncol,nlev)     !< average reflectivity will zero points outside valid range
  real(r8), intent(out) :: areflz(mgncol,nlev)    !< average reflectivity in z.
  real(r8), intent(out) :: frefl(mgncol,nlev)     !< fractional occurrence of radar reflectivity
  real(r8), intent(out) :: csrfl(mgncol,nlev)     !< cloudsat reflectivity
  real(r8), intent(out) :: acsrfl(mgncol,nlev)    !< cloudsat average
  real(r8), intent(out) :: fcsrfl(mgncol,nlev)    !< cloudsat fractional occurrence of radar reflectivity
  real(r8), intent(out) :: rercld(mgncol,nlev)    !< effective radius calculation for rain + cloud
  real(r8), intent(out) :: ncai(mgncol,nlev)      !< output number conc of ice nuclei available (1/m3)
  real(r8), intent(out) :: ncal(mgncol,nlev)      !< output number conc of CCN (1/m3)
  real(r8), intent(out) :: qrout2(mgncol,nlev)    !< copy of qrout as used to compute drout2
  real(r8), intent(out) :: qsout2(mgncol,nlev)    !< copy of qsout as used to compute dsout2
  real(r8), intent(out) :: nrout2(mgncol,nlev)    !< copy of nrout as used to compute drout2
  real(r8), intent(out) :: nsout2(mgncol,nlev)    !< copy of nsout as used to compute dsout2
  real(r8), intent(out) :: drout2(mgncol,nlev)    !< mean rain particle diameter (m)
  real(r8), intent(out) :: dsout2(mgncol,nlev)    !< mean snow particle diameter (m)
  real(r8), intent(out) :: freqs(mgncol,nlev)     !< fractional occurrence of snow
  real(r8), intent(out) :: freqr(mgncol,nlev)     !< fractional occurrence of rain
  real(r8), intent(out) :: nfice(mgncol,nlev)     !< fractional occurrence of ice
  real(r8), intent(out) :: qcrat(mgncol,nlev)     !< limiter for qc process rates (1=no limit --> 0. no qc)
!++ag
  real(r8), intent(out) :: qgout(mgncol,nlev)     !< graupel/hail mixing ratio (kg/kg)
  real(r8), intent(out) :: dgout(mgncol,nlev)     !< graupel/hail diameter (m)
  real(r8), intent(out) :: ngout(mgncol,nlev)     !< graupel/hail number concentration (1/m3)
!Not sure if these are needed since graupel/hail is prognostic?
  real(r8), intent(out) :: qgout2(mgncol,nlev)    !< copy of qgout as used to compute dgout2
  real(r8), intent(out) :: ngout2(mgncol,nlev)    !< copy of ngout as used to compute dgout2
  real(r8), intent(out) :: dgout2(mgncol,nlev)    !< mean graupel/hail particle diameter (m)
  real(r8), intent(out) :: freqg(mgncol,nlev)     !< fractional occurrence of graupel

!--ag

  real(r8), intent(out) :: prer_evap(mgncol,nlev)


  ! Tendencies calculated by external schemes that can replace MG's native
  ! process tendencies.

  ! Used with CARMA cirrus microphysics
  ! (or similar external microphysics model)
  ! real(r8), intent(in) :: tnd_qsnow(:,:)        !< snow mass tendency (kg/kg/s)
  ! real(r8), intent(in) :: tnd_nsnow(:,:)        !< snow number tendency (#/kg/s)
  ! real(r8), intent(in) :: re_ice(:,:)           !< ice effective radius (m)

  ! From external ice nucleation.
  !real(r8), intent(in) :: frzimm(:,:)            !< Number tendency due to immersion freezing (1/cm3)
  !real(r8), intent(in) :: frzcnt(:,:)            !< Number tendency due to contact freezing (1/cm3)
  !real(r8), intent(in) :: frzdep(:,:)            !< Number tendency due to deposition nucleation (1/cm3)

  ! local workspace
  ! all units mks unless otherwise stated

  ! local copies of input variables
  real(r8) :: qc(mgncol,nlev)         !< cloud liquid mixing ratio (kg/kg)
  real(r8) :: qi(mgncol,nlev)         !< cloud ice mixing ratio (kg/kg)
  real(r8) :: nc(mgncol,nlev)         !< cloud liquid number concentration (1/kg)
  real(r8) :: ni(mgncol,nlev)         !< cloud liquid number concentration (1/kg)
  real(r8) :: qr(mgncol,nlev)         !< rain mixing ratio (kg/kg)
  real(r8) :: qs(mgncol,nlev)         !< snow mixing ratio (kg/kg)
  real(r8) :: nr(mgncol,nlev)         !< rain number concentration (1/kg)
  real(r8) :: ns(mgncol,nlev)         !< snow number concentration (1/kg)
!++ag
  real(r8) :: qg(mgncol,nlev)         !< graupel mixing ratio (kg/kg)
  real(r8) :: ng(mgncol,nlev)         !< graupel number concentration (1/kg)
! real(r8) :: rhogtmp                 !< hail or graupel density (kg m-3)

!--ag

  ! general purpose variables
  real(r8) :: deltat                  !< sub-time step (s)
  real(r8) :: oneodt                  !< one / deltat
  real(r8) :: mtime                   !< the assumed ice nucleation timescale

  ! physical properties of the air at a given point
  real(r8) :: rho(mgncol,nlev)        ! density (kg m-3)
  real(r8) :: rhoinv(mgncol,nlev)     ! one / density (kg m-3)
  real(r8) :: dv(mgncol,nlev)         ! diffusivity of water vapor
  real(r8) :: mu(mgncol,nlev)         ! viscosity
  real(r8) :: sc(mgncol,nlev)         ! schmidt number
  real(r8) :: rhof(mgncol,nlev)       ! density correction factor for fallspeed

  ! cloud fractions
  real(r8) :: precip_frac(mgncol,nlev)! precip fraction assuming maximum overlap
  real(r8) :: cldm(mgncol,nlev)       ! cloud fraction
  real(r8) :: icldm(mgncol,nlev)      ! ice cloud fraction
  real(r8) :: lcldm(mgncol,nlev)      ! liq cloud fraction
  real(r8) :: qsfm(mgncol,nlev)       ! subgrid cloud water saturation scaling factor

  ! mass mixing ratios
  real(r8) :: qcic(mgncol,nlev)       ! in-cloud cloud liquid
  real(r8) :: qiic(mgncol,nlev)       ! in-cloud cloud ice
  real(r8) :: qsic(mgncol,nlev)       ! in-precip snow
  real(r8) :: qric(mgncol,nlev)       ! in-precip rain
!++ag
  real(r8) :: qgic(mgncol,nlev)       ! in-precip graupel/hail
!++ag


  ! number concentrations
  real(r8) :: ncic(mgncol,nlev)       ! in-cloud droplet
  real(r8) :: niic(mgncol,nlev)       ! in-cloud cloud ice
  real(r8) :: nsic(mgncol,nlev)       ! in-precip snow
  real(r8) :: nric(mgncol,nlev)       ! in-precip rain
!++ag
  real(r8) :: ngic(mgncol,nlev)       ! in-precip graupel/hail
!++ag

  ! maximum allowed ni value
  real(r8) :: nimax(mgncol,nlev)

  ! Size distribution parameters for:
  ! cloud ice
  real(r8) :: lami(mgncol,nlev)       ! slope
  real(r8) :: n0i(mgncol,nlev)        ! intercept
  ! cloud liquid
  real(r8) :: lamc(mgncol,nlev)       ! slope
  real(r8) :: pgam(mgncol,nlev)       ! spectral width parameter
  ! snow
  real(r8) :: lams(mgncol,nlev)       ! slope
  real(r8) :: n0s(mgncol,nlev)        ! intercept
  ! rain
  real(r8) :: lamr(mgncol,nlev)       ! slope
  real(r8) :: n0r(mgncol,nlev)        ! intercept
!++ag
  ! graupel/hail
  real(r8) :: lamg(mgncol,nlev)       ! slope
  real(r8) :: n0g(mgncol,nlev)        ! intercept
! real(r8) :: bgtmp                   ! tmp fall speed parameter
!--ag

  ! Rates/tendencies due to:

  ! Instantaneous snow melting
  real(r8) :: minstsm(mgncol,nlev)    ! mass mixing ratio
  real(r8) :: ninstsm(mgncol,nlev)    ! number concentration
!++ag
  ! Instantaneous graupel melting
  real(r8) :: minstgm(mgncol,nlev)    ! mass mixing ratio
  real(r8) :: ninstgm(mgncol,nlev)    ! number concentration
!--ag

  ! Instantaneous rain freezing
  real(r8) :: minstrf(mgncol,nlev)    ! mass mixing ratio
  real(r8) :: ninstrf(mgncol,nlev)    ! number concentration

  ! deposition of cloud ice
  real(r8) :: vap_dep(mgncol,nlev)    ! deposition from vapor to ice PMC 12/3/12
  ! sublimation of cloud ice
  real(r8) :: ice_sublim(mgncol,nlev) ! sublimation from ice to vapor PMC 12/3/12
  ! ice nucleation
  real(r8) :: nnuccd(mgncol,nlev)     ! number rate from deposition/cond.-freezing
  real(r8) :: mnuccd(mgncol,nlev)     ! mass mixing ratio
  ! freezing of cloud water
  real(r8) :: mnuccc(mgncol,nlev)     ! mass mixing ratio
  real(r8) :: nnuccc(mgncol,nlev)     ! number concentration
  ! contact freezing of cloud water
  real(r8) :: mnucct(mgncol,nlev)     ! mass mixing ratio
  real(r8) :: nnucct(mgncol,nlev)     ! number concentration
  ! deposition nucleation in mixed-phase clouds (from external scheme)
  real(r8) :: mnudep(mgncol,nlev)     ! mass mixing ratio
  real(r8) :: nnudep(mgncol,nlev)     ! number concentration
  ! ice multiplication
  real(r8) :: msacwi(mgncol,nlev)     ! mass mixing ratio
  real(r8) :: nsacwi(mgncol,nlev)     ! number concentration
  ! autoconversion of cloud droplets
  real(r8) :: prc(mgncol,nlev)        ! mass mixing ratio
  real(r8) :: nprc(mgncol,nlev)       ! number concentration (rain)
  real(r8) :: nprc1(mgncol,nlev)      ! number concentration (cloud droplets)
  ! self-aggregation of snow
  real(r8) :: nsagg(mgncol,nlev)      ! number concentration
  ! self-collection of rain
  real(r8) :: nragg(mgncol,nlev)      ! number concentration
  ! collection of droplets by snow
  real(r8) :: psacws(mgncol,nlev)     ! mass mixing ratio
  real(r8) :: npsacws(mgncol,nlev)    ! number concentration
  ! collection of rain by snow
  real(r8) :: pracs(mgncol,nlev)      ! mass mixing ratio
  real(r8) :: npracs(mgncol,nlev)     ! number concentration
  ! freezing of rain
  real(r8) :: mnuccr(mgncol,nlev)     ! mass mixing ratio
  real(r8) :: nnuccr(mgncol,nlev)     ! number concentration
  ! freezing of rain to form ice (mg add 4/26/13)
  real(r8) :: mnuccri(mgncol,nlev)    ! mass mixing ratio
  real(r8) :: nnuccri(mgncol,nlev)    ! number concentration
  ! accretion of droplets by rain
  real(r8) :: pra(mgncol,nlev)        ! mass mixing ratio
  real(r8) :: npra(mgncol,nlev)       ! number concentration
  ! autoconversion of cloud ice to snow
  real(r8) :: prci(mgncol,nlev)       ! mass mixing ratio
  real(r8) :: nprci(mgncol,nlev)      ! number concentration
  ! accretion of cloud ice by snow
  real(r8) :: prai(mgncol,nlev)       ! mass mixing ratio
  real(r8) :: nprai(mgncol,nlev)      ! number concentration
  ! evaporation of rain
  real(r8) :: pre(mgncol,nlev)        ! mass mixing ratio
  ! sublimation of snow
  real(r8) :: prds(mgncol,nlev)       ! mass mixing ratio
  ! number evaporation
  real(r8) :: nsubi(mgncol,nlev)      ! cloud ice
  real(r8) :: nsubc(mgncol,nlev)      ! droplet
  real(r8) :: nsubs(mgncol,nlev)      ! snow
  real(r8) :: nsubr(mgncol,nlev)      ! rain
  ! bergeron process
  real(r8) :: berg(mgncol,nlev)       ! mass mixing ratio (cloud ice)
  real(r8) :: bergs(mgncol,nlev)      ! mass mixing ratio (snow)

!++ag
  !graupel/hail processes
  real(r8) :: npracg(mgncol,nlev)     ! change n collection rain by graupel  (precipf)
  real(r8) :: nscng(mgncol,nlev)      ! change n conversion to graupel due to collection droplets by snow (lcldm)
  real(r8) :: ngracs(mgncol,nlev)     ! change n conversion to graupel due to collection rain by snow (precipf)
  real(r8) :: nmultg(mgncol,nlev)     ! ice mult due to acc droplets by graupel  (lcldm)
  real(r8) :: nmultrg(mgncol,nlev)    ! ice mult due to acc rain by graupel  (precipf)
  real(r8) :: npsacwg(mgncol,nlev)    ! change n collection droplets by graupel (lcldm)

  real(r8) :: psacr(mgncol,nlev)      ! conversion due to coll of snow by rain (precipf)
  real(r8) :: pracg(mgncol,nlev)      ! change in q collection rain by graupel  (precipf)
  real(r8) :: psacwg(mgncol,nlev)     ! change in q collection droplets by graupel (lcldm)
  real(r8) :: pgsacw(mgncol,nlev)     ! conversion q to graupel due to collection droplets by snow  (lcldm)
  real(r8) :: pgracs(mgncol,nlev)     ! conversion q to graupel due to collection rain by snow (precipf)
  real(r8) :: prdg(mgncol,nlev)       ! dep of graupel (precipf)
! real(r8) :: eprdg(mgncol,nlev)      ! evap/sub of graupel (precipf)
  real(r8) :: qmultg(mgncol,nlev)     ! change q due to ice mult droplets/graupel  (lcldm)
  real(r8) :: qmultrg(mgncol,nlev)    ! change q due to ice mult rain/graupel (precipf)
!--ag


  ! fallspeeds
  ! number-weighted
  real(r8) :: uns(mgncol,nlev)        ! snow
  real(r8) :: unr(mgncol,nlev)        ! rain
!++ag
  real(r8) :: ung(mgncol,nlev)        ! graupel/hail
!--ag
  ! air density corrected fallspeed parameters
  real(r8) :: arn(mgncol,nlev)        ! rain
  real(r8) :: asn(mgncol,nlev)        ! snow
!++a
  real(r8) :: agn(mgncol,nlev)        ! graupel
!--ag
  real(r8) :: acn(mgncol,nlev)        ! cloud droplet
  real(r8) :: ain(mgncol,nlev)        ! cloud ice
  real(r8) :: ajn(mgncol,nlev)        ! cloud small ice

  ! Mass of liquid droplets used with external heterogeneous freezing.
  real(r8) :: mi0l(mgncol)

  ! saturation vapor pressures
  real(r8) :: esl(mgncol,nlev)        ! liquid
  real(r8) :: esi(mgncol,nlev)        ! ice
  real(r8) :: esn                     ! checking for RH after rain evap

  ! saturation vapor mixing ratios
  real(r8) :: qvl(mgncol,nlev)        ! liquid
  real(r8) :: qvi(mgncol,nlev)        ! ice
  real(r8) :: qvn                     ! checking for RH after rain evap

  ! relative humidity
  real(r8) :: relhum(mgncol,nlev)

  ! parameters for cloud water and cloud ice sedimentation calculations
  real(r8) :: fc(mgncol,nlev)
  real(r8) :: fnc(mgncol,nlev)
  real(r8) :: fi(mgncol,nlev)
  real(r8) :: fni(mgncol,nlev)

!++ag
  real(r8) :: fg(mgncol,nlev)
  real(r8) :: fng(mgncol,nlev)
!--ag

  real(r8) :: fr(mgncol,nlev)
  real(r8) :: fnr(mgncol,nlev)
  real(r8) :: fs(mgncol,nlev)
  real(r8) :: fns(mgncol,nlev)

  real(r8) :: faloutc(nlev)
  real(r8) :: faloutnc(nlev)
  real(r8) :: falouti(nlev)
  real(r8) :: faloutni(nlev)

  real(r8) :: faloutr(nlev)
  real(r8) :: faloutnr(nlev)
  real(r8) :: falouts(nlev)
  real(r8) :: faloutns(nlev)

  real(r8) :: faltndc
  real(r8) :: faltndnc
  real(r8) :: faltndi
  real(r8) :: faltndni
  real(r8) :: faltndqie
  real(r8) :: faltndqce

  real(r8) :: faltndr
  real(r8) :: faltndnr
  real(r8) :: faltnds
  real(r8) :: faltndns

!++ag
  real(r8) :: faloutg(nlev)
  real(r8) :: faloutng(nlev)
  real(r8) :: faltndg
  real(r8) :: faltndng
!--ag

  real(r8) :: rainrt(mgncol,nlev)     ! rain rate for reflectivity calculation

  ! dummy variables
  real(r8) :: dum
  real(r8) :: dum1
  real(r8) :: dum2
!++ag
  real(r8) :: dum3
!--ag
  real(r8) :: dumni0
  real(r8) :: dumns0
  real(r8) :: tx1, tx2, tx3, tx4, tx5, tx6, tx7, grho
  ! dummies for checking RH
  real(r8) :: qtmp
  real(r8) :: ttmp
  ! dummies for conservation check
  real(r8) :: ratio
  real(r8) :: tmpfrz
  ! dummies for in-cloud variables
  real(r8) :: dumc(mgncol,nlev)   ! qc
  real(r8) :: dumnc(mgncol,nlev)  ! nc
  real(r8) :: dumi(mgncol,nlev)   ! qi
  real(r8) :: dumni(mgncol,nlev)  ! ni
  real(r8) :: dumr(mgncol,nlev)   ! rain mixing ratio
  real(r8) :: dumnr(mgncol,nlev)  ! rain number concentration
  real(r8) :: dums(mgncol,nlev)   ! snow mixing ratio
  real(r8) :: dumns(mgncol,nlev)  ! snow number concentration
!++ag
  real(r8) :: dumg(mgncol,nlev)   ! graupel mixing ratio
  real(r8) :: dumng(mgncol,nlev)  ! graupel number concentration
!--ag
  ! Array dummy variable
! real(r8) :: dum_2D(mgncol,nlev)
  real(r8) :: pdel_inv(mgncol,nlev)
  real(r8) :: ts_au_loc(mgncol)

  ! loop array variables
  ! "i" and "k" are column/level iterators for internal (MG) variables
  ! "n" is used for other looping (currently just sedimentation)
  integer i, k, n

  ! number of sub-steps for loops over "n" (for sedimentation)
  integer nstep, mdust, nlb, nstep_def

  ! Varaibles to scale fall velocity between small and regular ice regimes.
! real(r8) :: irad, ifrac, tsfac
  real(r8) :: irad, ifrac
! logical, parameter  :: do_ice_gmao=.false., do_liq_liu=.false.
! logical, parameter  :: do_ice_gmao=.false., do_liq_liu=.true.
! logical, parameter  :: do_ice_gmao=.true.,  do_liq_liu=.false.
! real(r8), parameter :: qimax=0.010, qimin=0.001, qiinv=one/(qimax-qimin), &
! real(r8), parameter :: qimax=0.010, qimin=0.001, qiinv=one/(qimax-qimin), &
  real(r8), parameter :: qimax=0.010_r8, qimin=0.005_r8, qiinv=one/(qimax-qimin)
!                        ts_au_min=180.0

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


  ! Process inputs

  !> - Assign variable deltat to deltatin
  deltat    = deltatin
  oneodt    = one / deltat
! nstep_def = max(1, nint(deltat/20))
  nstep_def = max(1, nint(deltat/5))
! tsfac     = log(ts_au/ts_au_min) * qiinv

  !> - Copies of input concentrations that may be changed internally.
  do k=1,nlev
    do i=1,mgncol
      qc(i,k) = qcn(i,k)
      nc(i,k) = ncn(i,k)
      qi(i,k) = qin(i,k)
      ni(i,k) = nin(i,k)
      qr(i,k) = qrn(i,k)
      nr(i,k) = nrn(i,k)
      qs(i,k) = qsn(i,k)
      ns(i,k) = nsn(i,k)
!++ag
      qg(i,k) = qgr(i,k)
      ng(i,k) = ngr(i,k)
    enddo
  enddo

  ! cldn: used to set cldm, unused for subcolumns
  ! liqcldf: used to set lcldm, unused for subcolumns
  ! icecldf: used to set icldm, unused for subcolumns
!> - Calculation liquid/ice cloud fraction
  if (microp_uniform) then
     ! subcolumns, set cloud fraction variables to one
     ! if cloud water or ice is present, if not present
     ! set to mincld (mincld used instead of zero, to prevent
     ! possible division by zero errors).

    do k=1,nlev
      do i=1,mgncol

        if (qc(i,k) >= qsmall) then
          lcldm(i,k) = one
        else
          lcldm(i,k) = mincld
        endif

        if (qi(i,k) >= qsmall) then
          icldm(i,k) = one
        else
          icldm(i,k) = mincld
        endif

        cldm(i,k) = max(icldm(i,k), lcldm(i,k))
!       qsfm(i,k) = one
        qsfm(i,k) = qsatfac(i,k)
      enddo
    enddo

  else         ! get cloud fraction, check for minimum
    do k=1,nlev
      do i=1,mgncol
        cldm(i,k)  = max(cldn(i,k),    mincld)
        lcldm(i,k) = max(liqcldf(i,k), mincld)
        icldm(i,k) = max(icecldf(i,k), mincld)
        qsfm(i,k)  = qsatfac(i,k)
      enddo
    enddo
  end if

! if (lprnt) write(0,*)' cldm=',cldm(1,nlev-20:nlev)
! if (lprnt) write(0,*)' liqcldf=',liqcldf(1,nlev-20:nlev)
! if (lprnt) write(0,*)' lcldm=',lcldm(1,nlev-20:nlev)
! if (lprnt) write(0,*)' icecldf=',icecldf(1,nlev-20:nlev)
! if (lprnt) write(0,*)' icldm=',icldm(1,nlev-20:nlev)
! if (lprnt) write(0,*)' qsfm=',qsfm(1,nlev-20:nlev)

  !> - Initialize local variables

  ! local physical properties

!    write(0,*)' in mg2 T=',t(1,:)
!    write(0,*)' in mg2 P=',p(1,:),' r=',r
  do k=1,nlev
    do i=1,mgncol
!     rho(i,k)    = p(i,k) / (r*t(i,k)*(one+fv*q(i,k)))
      rho(i,k)    = p(i,k) / (r*t(i,k))
      rhoinv(i,k) = one / rho(i,k)
      dv(i,k)     = 8.794E-5_r8 * t(i,k)**1.81_r8 / p(i,k)
      mu(i,k)     = 1.496E-6_r8 * t(i,k)*sqrt(t(i,k)) / (t(i,k) + 120._r8)
      sc(i,k)     = mu(i,k) / (rho(i,k)*dv(i,k))

  ! air density adjustment for fallspeed parameters
  ! includes air density correction factor to the
  ! power of 0.54 following Heymsfield and Bansemer 2007

      rhof(i,k) = (rhosu*rhoinv(i,k))**0.54_r8

      arn(i,k)  = ar    * rhof(i,k)
      asn(i,k)  = as    * rhof(i,k)
!++ag if do hail then agn = ah *rhof else ag*rhof
      agn(i,k)  = agtmp * rhof(i,k)
      acn(i,k)  = g*rhow/(18._r8*mu(i,k))
      tx1       = (rhosu*rhoinv(i,k))**0.35_r8
      ain(i,k)  = ai * tx1
      ajn(i,k)  = aj * tx1

  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Get humidity and saturation vapor pressures

! do k=1,nlev
!    do i=1,mgncol
!       relvar(i,k)      = relvar_i
        accre_enhan(i,k) = accre_enhan_i
!       call qsat_water(t(i,k), p(i,k), esl(i,k), qvl(i,k))
        esl(i,k)         = min(fpvsl(t(i,k)), p(i,k))
        qvl(i,k)         = epsqs*esl(i,k) / (p(i,k)-omeps*esl(i,k))


        ! make sure when above freezing that esi=esl, not active yet
        if (t(i,k) >= tmelt) then
           esi(i,k) = esl(i,k)
           qvi(i,k) = qvl(i,k)
        else
!          call qsat_ice(t(i,k), p(i,k), esi(i,k), qvi(i,k))
           esi(i,k) = min(fpvsi(t(i,k)), p(i,k))
           qvi(i,k) = epsqs*esi(i,k) / (p(i,k)-omeps*esi(i,k))
        end if

           ! Scale the water saturation values to reflect subgrid scale
           ! ice cloud fraction, where ice clouds begin forming at a
           ! gridbox average relative humidity of rhmini (not 1).
           !
           ! NOTE: For subcolumns and other non-subgrid clouds, qsfm will be 1.
        qvi(i,k) = qsfm(i,k) * qvi(i,k)
!       esi(i,k) = qsfm(i,k) * esi(i,k)
        qvl(i,k) = qsfm(i,k) * qvl(i,k)
!       esl(i,k) = qsfm(i,k) * esl(i,k)

        relhum(i,k) = max(zero, min(q(i,k)/max(qvl(i,k), qsmall), two))
     end do
  end do

  !===============================================

  ! set mtime here to avoid answer-changing
  mtime = deltat

  !> - initialize microphysics output
  do k=1,nlev
    do i=1,mgncol
      qcsevap(i,k)     = zero
      qisevap(i,k)     = zero
      qvres(i,k)       = zero
      cmeitot(i,k)     = zero
      vtrmc(i,k)       = zero
      vtrmi(i,k)       = zero
      qcsedten(i,k)    = zero
      qisedten(i,k)    = zero
      qrsedten(i,k)    = zero
      qssedten(i,k)    = zero
!++ag
      qgsedten(i,k)    = zero
!--ag


      pratot(i,k)      = zero
      prctot(i,k)      = zero
      mnuccctot(i,k)   = zero
      mnuccttot(i,k)   = zero
      msacwitot(i,k)   = zero
      psacwstot(i,k)   = zero
      bergstot(i,k)    = zero
      bergtot(i,k)     = zero
      melttot(i,k)     = zero
      homotot(i,k)     = zero
      qcrestot(i,k)    = zero
      prcitot(i,k)     = zero
      praitot(i,k)     = zero
      qirestot(i,k)    = zero
      mnuccrtot(i,k)   = zero
!++ag
      mnuccritot(i,k)  = zero
!--ag

      pracstot(i,k)    = zero
      meltsdttot(i,k)  = zero
      frzrdttot(i,k)   = zero
      mnuccdtot(i,k)   = zero

!++ag
      psacrtot(i,k)    = zero
      pracgtot(i,k)    = zero
      psacwgtot(i,k)   = zero
      pgsacwtot(i,k)   = zero
      pgracstot(i,k)   = zero
      prdgtot(i,k)     = zero
!     eprdgtot(i,k)    = zero
      qmultgtot(i,k)   = zero
      qmultrgtot(i,k)  = zero
      npracgtot(i,k)   = zero
      nscngtot(i,k)    = zero
      ngracstot(i,k)   = zero
      nmultgtot(i,k)   = zero
      nmultrgtot(i,k)  = zero
      npsacwgtot(i,k)  = zero
!need to zero these out to be totally switchable (for conservation)
      psacr(i,k)       = zero
      pracg(i,k)       = zero
      psacwg(i,k)      = zero
      pgsacw(i,k)      = zero
      pgracs(i,k)      = zero

      prdg(i,k)        = zero
!     eprdg(i,k)       = zero
      qmultg(i,k)      = zero
      qmultrg(i,k)     = zero
      npracg(i,k)      = zero
      nscng(i,k)       = zero
      ngracs(i,k)      = zero
      nmultg(i,k)      = zero
      nmultrg(i,k)     = zero
      npsacwg(i,k)     = zero
!--ag
      rflx(i,k+1)      = zero
      sflx(i,k+1)      = zero
      lflx(i,k+1)      = zero
      iflx(i,k+1)      = zero
!++ag
      gflx(i,k+1)      = zero
!--ag

  !> - initialize precip output

      qrout(i,k)       = zero
      qsout(i,k)       = zero
      nrout(i,k)       = zero
      nsout(i,k)       = zero
!++ag
      qgout(i,k)       = zero
      ngout(i,k)       = zero
      dgout(i,k)       = zero
!--ag

  ! for refl calc
      rainrt(i,k)      = zero

  !> - initialize rain size
      rercld(i,k)      = zero

      qcsinksum_rate1ord(i,k) = zero

  !> - initialize variables for trop_mozart
      nevapr(i,k)      = zero
      prer_evap(i,k)   = zero
      evapsnow(i,k)    = zero
      am_evp_st(i,k)   = zero
      prain(i,k)       = zero
      prodsnow(i,k)    = zero
      cmeout(i,k)      = zero

      precip_frac(i,k) = mincld

      lamc(i,k)        = zero

  !> - initialize microphysical tendencies

      tlat(i,k)   = zero
      qvlat(i,k)  = zero
      qctend(i,k) = zero
      qitend(i,k) = zero
      qstend(i,k) = zero
      qrtend(i,k) = zero
      nctend(i,k) = zero
      nitend(i,k) = zero
      nrtend(i,k) = zero
      nstend(i,k) = zero
!++ag
      qgtend(i,k) = zero
      ngtend(i,k) = zero
!--ag

  !> - initialize in-cloud and in-precip quantities to zero
      qcic(i,k)   = zero
      qiic(i,k)   = zero
      qsic(i,k)   = zero
      qric(i,k)   = zero
!++ag
      qgic(i,k)   = zero
!--ag


      ncic(i,k)   = zero
      niic(i,k)   = zero
      nsic(i,k)   = zero
      nric(i,k)   = zero
!++ag
      ngic(i,k)   = zero
!--ag
  !> - initialize precip fallspeeds to zero
      ums(i,k)    = zero 
      uns(i,k)    = zero
      umr(i,k)    = zero
      unr(i,k)    = zero
!++ag
      umg(i,k)    = zero
      ung(i,k)    = zero
!--ag

  !> - initialize limiter for output
      qcrat(i,k)  = one

  ! Many outputs have to be initialized here at the top to work around
  ! ifort problems, even if they are always overwritten later.
      effc(i,k)    = ten
      lamcrad(i,k) = zero
      pgamrad(i,k) = zero
      effc_fn(i,k) = ten
      effi(i,k)    = 25._r8
      sadice(i,k)  = zero
      sadsnow(i,k) = zero
      deffi(i,k)   = 50._r8

      qrout2(i,k) = zero
      nrout2(i,k) = zero
      drout2(i,k) = zero
      qsout2(i,k) = zero
      nsout2(i,k) = zero
      dsout(i,k)  = zero
      dsout2(i,k) = zero
!++ag
      qgout2(i,k) = zero
      ngout2(i,k) = zero
      freqg(i,k)  = zero
      dgout2(i,k) = zero
!--ag

      freqr(i,k) = zero
      freqs(i,k) = zero

      reff_rain(i,k) = zero
      reff_snow(i,k) = zero
!++ag
      reff_grau(i,k) = zero
      lamg(i,k)      = zero
      n0g(i,k)       = zero
!--ag

      refl(i,k)   = -9999._r8
      arefl(i,k)  = zero
      areflz(i,k) = zero
      frefl(i,k)  = zero
      csrfl(i,k)  = zero
      acsrfl(i,k) = zero
      fcsrfl(i,k) = zero

      ncal(i,k)   = zero
      ncai(i,k)   = zero

      nfice(i,k)  = zero
      npccn(i,k)  = zero
    enddo
  enddo
!> - initialize ccn activated number tendency (\p npccn)
  if (iccn == 1) then
    do k=1,nlev
      do i=1,mgncol
        npccn(i,k) = npccnin(i,k)
      enddo
    enddo
  else
    do k=1,nlev
      do i=1,mgncol
        npccn(i,k) = max((npccnin(i,k)*lcldm(i,k)-nc(i,k))*oneodt, zero)
      enddo
    enddo
  endif

  !> - initialize precip at surface

  do i=1,mgncol
      prect(i) = zero
      preci(i) = zero
  enddo

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! droplet activation
  ! get provisional droplet number after activation. This is used for
  ! all microphysical process calculations, for consistency with update of
  ! droplet mass before microphysics

  ! calculate potential for droplet activation if cloud water is present
  ! tendency from activation (npccn) is read in from companion routine

  ! output activated liquid and ice (convert from #/kg -> #/m3)
  !--------------------------------------------------
! where (qc >= qsmall .and. lcldm > mincld)
! where (qc >= qsmall)
!    npccn = max((npccnin*lcldm-nc)*oneodt, zero)
!    nc    = max(nc + npccn*deltat, zero)
!    ncal  = nc*rho/lcldm ! sghan minimum in #/cm3
! elsewhere
!    ncal = zero
! end where

! if (lprnt) write(0,*)' nc1=',nc(1,:)
  do k=1,nlev
    do i=1,mgncol
      if (qc(i,k) > qsmall .and. lcldm(i,k) >= mincld) then
        npccn(i,k) = max((npccnin(i,k)*lcldm(i,k)-nc(i,k))*oneodt, zero)
        nc(i,k)    = max(nc(i,k) + npccn(i,k)*deltat, zero)
        ncal(i,k)  = nc(i,k) * rho(i,k) / lcldm(i,k)
      else
        ncal(i,k) = 0.0
      endif
    enddo
  enddo

  if (iccn == 1) then
    do k=1,nlev
      do i=1,mgncol
        if (t(i,k) < icenuct) then
          ncai(i,k) = 0.005_r8*exp(0.304_r8*(273.15_r8-t(i,k))) * 1000._r8
!         ncai(i,k) = min(ncai(i,k), 208.9e3_r8)
          ncai(i,k) = min(ncai(i,k), 355.0e3_r8)
          naai(i,k) = (ncai(i,k)*rhoinv(i,k) + naai(i,k))  * half
          ncai(i,k) = naai(i,k)*rho(i,k)
        else
          naai(i,k) = zero
          ncai(i,k) = zero
        endif
      enddo
    enddo
  elseif (iccn == 2) then
    do k=1,nlev
      do i=1,mgncol
        if (t(i,k) < icenuct) then
          ncai(i,k) = naai(i,k)*rho(i,k)
          ncai(i,k) = min(ncai(i,k), 710.0e3_r8)
          naai(i,k) = ncai(i,k)*rhoinv(i,k)
        else
          naai(i,k) = zero
          ncai(i,k) = zero
        endif
      enddo
    enddo
  else
    do k=1,nlev
      do i=1,mgncol
        if (t(i,k) < icenuct) then
          ncai(i,k) = 0.005_r8*exp(0.304_r8*(273.15_r8-t(i,k))) * 1000._r8
          ncai(i,k) = min(ncai(i,k), 355.0e3_r8)
          naai(i,k) = ncai(i,k)*rhoinv(i,k) 
        else
          naai(i,k) = zero
          ncai(i,k) = zero
        endif
      enddo
    enddo

  endif


  !===============================================

  ! ice nucleation if activated nuclei exist at t<-5C AND rhmini + 5%
  !
  ! NOTE: If using gridbox average values, condensation will not occur until rh=1,
  ! so the threshold seems like it should be 1.05 and not rhmini + 0.05. For subgrid
  ! clouds (using rhmini and qsfacm), the relhum has already been adjusted, and thus
  ! the nucleation threshold should also be 1.05 and not rhmini + 0.05.

  !-------------------------------------------------------

  if (do_cldice) then
     where (naai > zero .and. t < icenuct .and. relhum*esl/esi > 1.05_r8)
!    where (naai > zero .and. t < icenuct .and. relhum*esl/esi > 1.05_r8 &
!                                         .and. icldm > mincld )

        !if NAAI > 0. then set numice = naai (as before)
        !note: this is gridbox averaged
        nnuccd = (naai-ni/icldm)/mtime*icldm
        nnuccd = max(nnuccd, zero)
        nimax  = naai*icldm

        !Calc mass of new particles using new crystal mass...
        !also this will be multiplied by mtime as nnuccd is...

        mnuccd = nnuccd * mi0

     elsewhere
        nnuccd = zero
        nimax  = zero
        mnuccd = zero
     end where

  end if


  !=============================================================================
  do k=1,nlev

     do i=1,mgncol

        ! calculate instantaneous precip processes (melting and homogeneous freezing)

        ! melting of snow at +2 C

        if (t(i,k) > snowmelt) then
           if (qs(i,k) > zero) then

              ! make sure melting snow doesn't reduce temperature below threshold
              dum = -(xlf/cpp) * qs(i,k)
              if (t(i,k)+dum < snowmelt) then
                 dum = min(one, max(zero, (cpp/xlf)*(t(i,k)-snowmelt)/qs(i,k)))
              else
                 dum = one
              end if

              minstsm(i,k) = dum*qs(i,k)
              ninstsm(i,k) = dum*ns(i,k)

              dum1 = - minstsm(i,k) * (xlf*oneodt)
              tlat(i,k)       = tlat(i,k)       + dum1
              meltsdttot(i,k) = meltsdttot(i,k) + dum1

!     if (lprnt .and. k >=40) write(0,*)' tlats=',tlat(i,k),' dum1=',dum1,&
!    ' minstsm=',minstsm(i,k),' qs=',qs(i,k),' xlf=',xlf,' oneodt=',oneodt, &
!    ' snowmelt=',snowmelt,' t=',t(i,k),' dum=',dum,' k=',k

              qs(i,k) = max(qs(i,k) - minstsm(i,k), zero)
              ns(i,k) = max(ns(i,k) - ninstsm(i,k), zero)
              qr(i,k) = max(qr(i,k) + minstsm(i,k), zero)
              nr(i,k) = max(nr(i,k) + ninstsm(i,k), zero)
           end if
        end if

     end do
  end do 
! if (lprnt) write(0,*)' tlat1=',tlat(1,:)*deltat
! if (lprnt) write(0,*)' qg1=',qg(1,:)

!++ag

  if (do_graupel .or. do_hail) then
! melting of graupel at +2 C

    do k=1,nlev
      do i=1,mgncol

        if (t(i,k) > snowmelt) then
          if (qg(i,k) > zero) then

! make sure melting graupel doesn't reduce temperature below threshold
            dum = -(xlf/cpp) * qg(i,k)
            if (t(i,k)+dum < snowmelt) then
              dum = max(zero, min(one, (cpp/xlf)*(t(i,k)-snowmelt)/qg(i,k)))
            else
              dum = one
            end if

            minstgm(i,k) = dum*qg(i,k)
            ninstgm(i,k) = dum*ng(i,k)

            dum1            = - minstgm(i,k) * (xlf*oneodt)
            tlat(i,k)       = dum1 + tlat(i,k)
            meltsdttot(i,k) = dum1 + meltsdttot(i,k)

!     if (lprnt .and. k >=40) write(0,*)' tlatg=',tlat(i,k),' dum1=',dum1,&
!    ' minstgm=',minstgm(i,k),' qg=',qg(i,k),' xlf=',xlf,' oneodt=',oneodt, &
!    ' snowmelt=',snowmelt,' t=',t(i,k),' k=',k,' cpp=',cpp

            qg(i,k) = max(qg(i,k) - minstgm(i,k), zero)
            ng(i,k) = max(ng(i,k) - ninstgm(i,k), zero)
            qr(i,k) = max(qr(i,k) + minstgm(i,k), zero)
            nr(i,k) = max(nr(i,k) + ninstgm(i,k), zero)
          end if
        end if

      end do
    end do
  endif

! if (lprnt) write(0,*)' tlat1g=',tlat(1,:)*deltat
!--ag

  do k=1,nlev
    do i=1,mgncol
        ! freezing of rain at -5 C

        if (t(i,k) < rainfrze) then

           if (qr(i,k) > zero) then

              ! make sure freezing rain doesn't increase temperature above threshold
              dum = (xlf/cpp) * qr(i,k)
              if (t(i,k)+dum > rainfrze) then
                 dum = -(t(i,k)-rainfrze) * (cpp/xlf)
                 dum = min(one, max(zero, dum/qr(i,k)))
              else
                 dum = one
              end if

              minstrf(i,k) = dum*qr(i,k)
              ninstrf(i,k) = dum*nr(i,k)

              ! heating tendency
              dum1 = minstrf(i,k) * (xlf*oneodt)
              tlat(i,k)      = tlat(i,k)      + dum1
              frzrdttot(i,k) = frzrdttot(i,k) + dum1

              qr(i,k) = max(qr(i,k) - minstrf(i,k), zero)
              nr(i,k) = max(nr(i,k) - ninstrf(i,k), zero)

!++ag
! freeze rain to graupel not snow.
              if(do_hail .or. do_graupel) then
                qg(i,k) = max(qg(i,k) + minstrf(i,k), zero)
                ng(i,k) = max(ng(i,k) + ninstrf(i,k), zero)
              else
                qs(i,k) = max(qs(i,k) + minstrf(i,k), zero)
                ns(i,k) = max(ns(i,k) + ninstrf(i,k), zero)
              end if
!--ag

           end if
        end if
     end do
  end do 

! if (lprnt) then
!  write(0,*)' tlat2=',tlat(1,:)*deltat
!  write(0,*)' lcldm=',lcldm(1,:)
!  write(0,*)' qc=',qc(1,:)
!  write(0,*)' nc=',nc(1,:)
!  write(0,*)' qg2=',qg(1,:)
! endif

  do k=1,nlev
    do i=1,mgncol
        ! obtain in-cloud values of cloud water/ice mixing ratios and number concentrations
        !-------------------------------------------------------
        ! for microphysical process calculations
        ! units are kg/kg for mixing ratio, 1/kg for number conc

!       if (qc(i,k) >= qsmall .and. lcldm(i,k) > mincld) then
        if (qc(i,k) >= qsmall) then
           ! limit in-cloud values to 0.005 kg/kg
           dum = one / lcldm(i,k)
!          qcic(i,k) = min(qc(i,k)*dum, 5.e-3_r8) ! limit in-cloud values to 0.005 kg/kg
           qcic(i,k) = min(qc(i,k)*dum, 0.05_r8)  ! limit in-cloud values to 0.05 kg/kg
           ncic(i,k) = max(nc(i,k)*dum, zero)

           ! specify droplet concentration
           if (nccons) then
              ncic(i,k) = ncnst * rhoinv(i,k)
           end if
        else
           qcic(i,k) = zero
           ncic(i,k) = zero
        end if

!       if (qi(i,k) >= qsmall .and. icldm(i,k) > mincld) then
        if (qi(i,k) >= qsmall) then
           dum = one / icldm(i,k)
!          qiic(i,k) = min(qi(i,k)*dum, 5.e-3_r8) ! limit in-cloud values to 0.005 kg/kg
           qiic(i,k) = min(qi(i,k)*dum, 0.05_r8)  ! limit in-cloud values to 0.05 kg/kg
           niic(i,k) = max(ni(i,k)*dum, zero)

           ! switch for specification of cloud ice number
           if (nicons) then
              niic(i,k) = ninst * rhoinv(i,k)
           end if
        else
           qiic(i,k) = zero
           niic(i,k) = zero
        end if

    end do
  end do

  !========================================================================

  ! for sub-columns cldm has already been set to 1 if cloud
  ! water or ice is present, so precip_frac will be correctly set below
  ! and nothing extra needs to be done here

  precip_frac = cldm

  micro_vert_loop: do k=1,nlev

     if (trim(micro_mg_precip_frac_method) == 'in_cloud') then

        if (k /= 1) then
           where (qc(:,k) < qsmall .and. qi(:,k) < qsmall)
              precip_frac(:,k) = precip_frac(:,k-1)
           end where
        endif

     else if (trim(micro_mg_precip_frac_method) == 'max_overlap') then

!++ag  add graupel to precip frac?
        ! calculate precip fraction based on maximum overlap assumption

        ! if rain or snow mix ratios are smaller than threshold,
        ! then leave precip_frac as cloud fraction at current level
        if (k /= 1) then
!++ag
!           where (qr(:,k-1) >= qsmall .or. qs(:,k-1) >= qsmall .or. qg(:,k-1) >= qsmall)
!--ag
           where (qr(:,k-1) >= qsmall .or. qs(:,k-1) >= qsmall)
              precip_frac(:,k) = max(precip_frac(:,k-1), precip_frac(:,k))
           end where
        end if

     endif


     !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
     ! get size distribution parameters based on in-cloud cloud water
     ! these calculations also ensure consistency between number and mixing ratio
     !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

     ! cloud liquid
     !-------------------------------------------

!     if (lprnt .and. k>=60 .and. k<=65) then
!     if (lprnt .and. k>=100) then
!     if (lprnt) then
!      write(0,*)' pgam=',pgam(1,k), ' qcic=',qcic(1,k),' ncic=',ncic(1,k),' rho=',rho(1,k),' k=',k
!     endif
     call size_dist_param_liq(mg_liq_props, qcic(:,k), ncic(:,k), rho(:,k), &
                              pgam(:,k), lamc(:,k), mgncol)
!     if (lprnt .and. k>=60 .and. k<=65) then
!     if (lprnt .and. k>=100) then
!     if (lprnt) then
!      write(0,*)' pgam=',pgam(1,k), ' lamc=',lamc(1,k),' k=',k
!     endif


     !========================================================================
     ! autoconversion of cloud liquid water to rain
     ! formula from Khrouditnov and Kogan (2000), modified for sub-grid distribution of qc
     ! minimum qc of 1 x 10^-8 prevents floating point error

     if (.not. do_sb_physics) then
       call kk2000_liq_autoconversion(microp_uniform, qcic(:,k), &
          ncic(:,k), rho(:,k), relvar(:,k), prc(:,k), nprc(:,k), nprc1(:,k), mgncol)
     endif

     ! assign qric based on prognostic qr, using assumed precip fraction
     ! note: this could be moved above for consistency with qcic and qiic calculations
     do i=1,mgncol
       if (precip_frac(i,k) > mincld) then
         dum = one / precip_frac(i,k)
       else
         dum = zero
       endif
!      qric(i,k) = min(qr(i,k)*dum, 0.01_r8) ! limit in-precip mixing ratios to 10 g/kg
       qric(i,k) = min(qr(i,k)*dum, 0.05_r8) ! limit in-precip mixing ratios to 50 g/kg
       nric(i,k) = nr(i,k) * dum


     ! add autoconversion to precip from above to get provisional rain mixing ratio
     ! and number concentration (qric and nric)

       if(qric(i,k) < qsmall) then
         qric(i,k) = zero
         nric(i,k) = zero
       endif

     ! make sure number concentration is a positive number to avoid
     ! taking root of negative later

       nric(i,k) = max(nric(i,k),zero)
     enddo 
     ! Get size distribution parameters for cloud ice

     call size_dist_param_ice(mg_ice_props, qiic(:,k), niic(:,k),   &
          lami(:,k), mgncol, n0=n0i(:,k))
  
     ! Alternative autoconversion 
     if (do_sb_physics) then
       if  (do_liq_liu) then 
         call liu_liq_autoconversion(pgam(:,k),qcic(:,k),ncic(:,k),        &
              qric(:,k),rho(:,k),relvar(:,k),prc(:,k),nprc(:,k),nprc1(:,k),mgncol)
       else
         call sb2001v2_liq_autoconversion(pgam(:,k),qcic(:,k),ncic(:,k),   &
              qric(:,k),rho(:,k),relvar(:,k),prc(:,k),nprc(:,k),nprc1(:,k), mgncol)     
       endif
     endif

     !.......................................................................
     ! Autoconversion of cloud ice to snow
     ! similar to Ferrier (1994)

     if (do_cldice) then
        do i=1,mgncol
          if (qiic(i,k) >= qimax) then
!         if (qi(i,k) >= qimax) then
            ts_au_loc(i) = ts_au_min
          elseif (qiic(i,k) <= qimin) then
!         elseif (qi(i,k) <= qimin) then
            ts_au_loc(i) = ts_au
          else
!           ts_au_loc(i) = (ts_au*(qimax-qi(i,k)) + ts_au_min*(qi(i,k)-qimin)) * qiinv
            ts_au_loc(i) = (ts_au*(qimax-qiic(i,k)) + ts_au_min*(qiic(i,k)-qimin)) * qiinv
!           ts_au_loc(i) = ts_au * exp(-tsfac*(qiic(i,k)-qimin))
          endif
!         if (ts_au_loc(i) > ts_au_min) ts_au_loc(i) = ts_au_loc(i)*min(five,sqrt(p(i,nlev)/p(i,k)))
        enddo
!       if (lprnt) write(0,*)' ts_au_loc=',ts_au_loc(1),' k=',k, ' qiic=',qiic(1,k),&
!       if (lprnt) write(0,*)' ts_au_loc=',ts_au_loc(1),' k=',k, ' qi=',qi(1,k),&
!        ' ts_au=',ts_au,' ts_au_min=',ts_au_min,' qimin=',qimin,' qimax=',qimax
!        ' ts_au=',ts_au,' ts_au_min=',ts_au_min,' tsfac=',tsfac,' qimin=',qimin,' qimax=',qimax

        if(do_ice_gmao) then
          call gmao_ice_autoconversion(t(:,k), qiic(:,k), niic(:,k), lami(:,k),     &
               n0i(:,k), dcs, ts_au_loc(:), prci(:,k), nprci(:,k), mgncol)
        else
          call ice_autoconversion(t(:,k), qiic(:,k), lami(:,k), n0i(:,k),           &
                         dcs, ts_au_loc(:), prci(:,k), nprci(:,k), mgncol)
        end if
     !else
        ! Add in the particles that we have already converted to snow, and
        ! don't do any further autoconversion of ice.
        !prci(:,k)  = tnd_qsnow(:,k) / cldm(:,k)
        !nprci(:,k) = tnd_nsnow(:,k) / cldm(:,k)
     end if

     ! note, currently we don't have this
     ! inside the do_cldice block, should be changed later
     ! assign qsic based on prognostic qs, using assumed precip fraction
     do i=1,mgncol
       if (precip_frac(i,k) > mincld) then
         dum = one / precip_frac(i,k)
       else
         dum = zero
       endif
!      qsic(i,k) = min(qs(i,k)*dum, 0.01_r8) ! limit in-precip mixing ratios to 10 g/kg
       qsic(i,k) = min(qs(i,k)*dum, 0.05_r8) ! limit in-precip mixing ratios to 50 g/kg
       nsic(i,k) = ns(i,k) * dum

     ! if precip mix ratio is zero so should number concentration

       if(qsic(i,k) < qsmall) then
         qsic(i,k) = zero
         nsic(i,k) = zero
       endif

     ! make sure number concentration is a positive number to avoid
     ! taking root of negative later

       nsic(i,k) = max(nsic(i,k), zero)

!++ also do this for graupel, which is assumed to be 'precip_frac'
       qgic(i,k) = min(qg(i,k)*dum, 0.01_r8) ! limit in-precip mixing ratios to 10 g/kg)
       ngic(i,k) = ng(i,k) * dum

     ! if precip mix ratio is zero so should number concentration
       if (qgic(i,k) < qsmall) then
         qgic(i,k) = zero
         ngic(i,k) = zero
       endif

     ! make sure number concentration is a positive number to avoid
     ! taking root of negative later

       ngic(i,k) = max(ngic(i,k), zero)
!--ag
     enddo

     !.......................................................................
     ! get size distribution parameters for precip
     !......................................................................
     ! rain

     call size_dist_param_basic(mg_rain_props, qric(:,k), nric(:,k), &
                                lamr(:,k), mgncol, n0=n0r(:,k))

     do i=1,mgncol
       if (lamr(i,k) >= qsmall) then
         dum  = arn(i,k) / lamr(i,k)**br
         dum1 = 9.1_r8*rhof(i,k)

        ! provisional rain number and mass weighted mean fallspeed (m/s)

         umr(i,k) = min(dum1, dum*gamma_br_plus4*oneo6)
         unr(i,k) = min(dum1, dum*gamma_br_plus1)
       else

         umr(i,k) = zero
         unr(i,k) = zero
       endif
     enddo

     !......................................................................
     ! snow

     call size_dist_param_basic(mg_snow_props, qsic(:,k), nsic(:,k), &
                                lams(:,k), mgncol, n0=n0s(:,k))

     do i=1,mgncol
       if (lams(i,k) >= qsmall) then

        ! provisional snow number and mass weighted mean fallspeed (m/s)

         dum  = asn(i,k) / lams(i,k)**bs
         dum1 = 1.2_r8*rhof(i,k)
         ums(i,k) = min(dum1, dum*gamma_bs_plus4*oneo6)
         uns(i,k) = min(dum1, dum*gamma_bs_plus1)

       else
         ums(i,k) = zero
         uns(i,k) = zero
       endif
     enddo

     if (do_graupel .or. do_hail) then
!++ag
!use correct bg or bh  (bgtmp=bg or bh)
     !......................................................................
     ! graupel/hail

!++AG SET rhog here and for mg_graupel_props?
! For now: rhog is constant. Set to same in micro_mg_utils.F90
! Ideally: find a method to set once. (Hail = 400, Graupel = 500 from M2005)

!mg,snow_props or mg_graupel props?

       call size_dist_param_basic(mg_graupel_props, qgic(:,k), ngic(:,k), &
                                  lamg(:,k), mgncol, n0=n0g(:,k))

       do i=1,mgncol
         if (lamg(i,k) >= qsmall) then

        ! provisional graupel/hail number and mass weighted mean fallspeed (m/s)

           dum  = agn(i,k) / lamg(i,k)**bgtmp
           dum1 = 20._r8*rhof(i,k)
           umg(i,k) = min(dum1, dum*gamma_bg_plus4*oneo6)
           ung(i,k) = min(dum1, dum*gamma_bg_plus1)
!          umg(i,k) = min(dum1, dum*gamma(four+bgtmp)*oneo6)
!          ung(i,k) = min(dum1, dum*gamma(one+bgtmp))

         else
           umg(i,k) = zero
           ung(i,k) = zero
         endif
       enddo
!--ag
     endif

     if (do_cldice) then
       if (.not. use_hetfrz_classnuc) then

         ! heterogeneous freezing of cloud water
         !----------------------------------------------

!        if (lprnt .and. k>=60 .and. k<=65) then
!        if (lprnt .and. k>=100) then
!        if (lprnt) then
!         write(0,*)' pgam=',pgam(1,k), ' lamc=',lamc(1,k),' qcic=',qcic(1,k),' ncic=',ncic(1,k),' t=',t(1,k),' k=',k,&
!                  ' relvar=',relvar(1,k)
!        endif

         call immersion_freezing(microp_uniform, t(:,k), pgam(:,k), lamc(:,k), &
                qcic(:,k), ncic(:,k), relvar(:,k), mnuccc(:,k), nnuccc(:,k), mgncol)

!        if (lprnt .and. k>=60 .and. k<=65) then
!        if (lprnt .and. k>=100) then
!        if (lprnt) then
!         write(0,*)' mnuccca=',mnuccc(1,k),' lcldm=',lcldm(1,k),' nnuccc=',nnuccc(1,k),' k=',k
!        endif

         ! make sure number of droplets frozen does not exceed available ice nuclei concentration
         ! this prevents 'runaway' droplet freezing

!        where (qcic(:,k) >= qsmall .and. t(:,k) < 269.15_r8 .and. lcldm(:,k) > mincld)
         where (qcic(:,k) >= qsmall .and. t(:,k) < 269.15_r8)
           where (nnuccc(:,k)*lcldm(:,k) > nnuccd(:,k))
              ! scale mixing ratio of droplet freezing with limit
              mnuccc(:,k) = mnuccc(:,k)*(nnuccd(:,k)/(nnuccc(:,k)*lcldm(:,k)))
              nnuccc(:,k) = nnuccd(:,k)/lcldm(:,k)
           end where
         end where

!     if (lprnt .and. k >= 60 .and. k <=65) write(0,*)' mnuccc=',mnuccc(1,60:65)
!     if (lprnt .and. k >= 100) write(0,*)' mnuccc=',mnuccc(1,k)
!     if (lprnt) write(0,*)' mnuccc=',mnuccc(1,k)

         mdust = size(rndst,3)
         call contact_freezing(microp_uniform, t(:,k), p(:,k), rndst(:,k,:),   &
              nacon(:,k,:), pgam(:,k), lamc(:,k), qcic(:,k), ncic(:,k), &
              relvar(:,k), mnucct(:,k), nnucct(:,k), mgncol, mdust)

!     if (lprnt .and. k >= 60 .and. k <=65) write(0,*)' mnucct=',mnucct(1,:)
!     if (lprnt .and. k >= 100 ) write(0,*)' mnucct=',mnucct(1,k)
!     if (lprnt) write(0,*)' mnucct=',mnucct(1,k)

         mnudep(:,k) = zero
         nnudep(:,k) = zero

      !else

           ! Mass of droplets frozen is the average droplet mass, except
           ! with two limiters: concentration must be at least 1/cm^3, and
           ! mass must be at least the minimum defined above.
           !mi0l = qcic(:,k)/max(ncic(:,k), 1.0e6_r8/rho(:,k))
           !mi0l = max(mi0l_min, mi0l)

           !where (qcic(:,k) >= qsmall)
              !nnuccc(:,k) = frzimm(:,k)*1.0e6_r8*rhoinv(:,k)
              !mnuccc(:,k) = nnuccc(:,k)*mi0l

              !nnucct(:,k) = frzcnt(:,k)*1.0e6_r8*rhoinv(:,k)
              !mnucct(:,k) = nnucct(:,k)*mi0l

              !nnudep(:,k) = frzdep(:,k)*1.0e6_r8*rhoinv(:,k)
              !mnudep(:,k) = nnudep(:,k)*mi0
           !elsewhere
              !nnuccc(:,k) = zero
              !mnuccc(:,k) = zero

              !nnucct(:,k) = zero
              !mnucct(:,k) = zero

              !nnudep(:,k) = zero
              !mnudep(:,k) = zero
           !end where

       end if

     else
       do i=1,mgncol
         mnuccc(i,k) = zero
         nnuccc(i,k) = zero
         mnucct(i,k) = zero
         nnucct(i,k) = zero
         mnudep(i,k) = zero
         nnudep(i,k) = zero
       enddo
     end if

     call snow_self_aggregation(t(:,k), rho(:,k), asn(:,k), rhosn, qsic(:,k), nsic(:,k), &
          nsagg(:,k), mgncol)

     call accrete_cloud_water_snow(t(:,k), rho(:,k), asn(:,k), uns(:,k), mu(:,k), &
          qcic(:,k), ncic(:,k), qsic(:,k), pgam(:,k), lamc(:,k), lams(:,k), n0s(:,k), &
          psacws(:,k), npsacws(:,k), mgncol)

     if (do_cldice) then
        call secondary_ice_production(t(:,k), psacws(:,k), msacwi(:,k), nsacwi(:,k), mgncol)
     else
        nsacwi(:,k) = zero
        msacwi(:,k) = zero
     end if

     call accrete_rain_snow(t(:,k), rho(:,k), umr(:,k), ums(:,k), unr(:,k), uns(:,k), &
          qric(:,k), qsic(:,k), lamr(:,k), n0r(:,k), lams(:,k), n0s(:,k), &
          pracs(:,k), npracs(:,k), mgncol)

     call heterogeneous_rain_freezing(t(:,k), qric(:,k), nric(:,k), lamr(:,k), &
          mnuccr(:,k), nnuccr(:,k), mgncol)

     if (do_sb_physics) then
       call sb2001v2_accre_cld_water_rain(qcic(:,k), ncic(:,k), qric(:,k), &
            rho(:,k), relvar(:,k), pra(:,k), npra(:,k), mgncol)     
     else
       call accrete_cloud_water_rain(microp_uniform, qric(:,k), qcic(:,k), &
            ncic(:,k), relvar(:,k), accre_enhan(:,k), pra(:,k), npra(:,k), mgncol)
     endif

     call self_collection_rain(rho(:,k), qric(:,k), nric(:,k), nragg(:,k), mgncol)

     if (do_cldice) then
        call accrete_cloud_ice_snow(t(:,k), rho(:,k), asn(:,k), qiic(:,k), niic(:,k), &
             qsic(:,k), lams(:,k), n0s(:,k), prai(:,k), nprai(:,k), mgncol)
     else
        prai(:,k)  = zero
        nprai(:,k) = zero
     end if

!++ag Moved below graupel conditional, now two different versions
!    if (.not. (do_hail .or. do_graupel))  then
!      call evaporate_sublimate_precip(t(:,k), rho(:,k),                           &
!          dv(:,k), mu(:,k), sc(:,k), q(:,k), qvl(:,k), qvi(:,k),                  &
!          lcldm(:,k), precip_frac(:,k), arn(:,k), asn(:,k), qcic(:,k), qiic(:,k), &
!          qric(:,k), qsic(:,k), lamr(:,k), n0r(:,k), lams(:,k), n0s(:,k),         &
!          pre(:,k), prds(:,k), am_evp_st(:,k), mgncol)
!    endif
!--ag

     call bergeron_process_snow(t(:,k), rho(:,k), dv(:,k), mu(:,k), sc(:,k),       &
          qvl(:,k), qvi(:,k), asn(:,k), qcic(:,k), qsic(:,k), lams(:,k), n0s(:,k), &
          bergs(:,k), mgncol)
!     if(lprnt) write(0,*)' bergs1=',bergs(1,k),' k=',k,' micro_mg_berg_eff_factor=',micro_mg_berg_eff_factor
!     if(lprnt) write(0,*)' t=',t(1,k),' rho=',rho(1,k),' dv=',dv(1,k),' mu=',mu(1,k),&
!        'qcic=',qcic(1,k),' qsic=',qsic(1,k),' qvl=',qvl(1,k),' qvi=',qvi(1,k),      &
!        ' mu=',mu(1,k),' sc=',sc(1,k),' asn=',asn(1,k),' lams=',lams(1,k),' n0s=',n0s(1,k),' ni=',ni(1,k)

     bergs(:,k) = bergs(:,k) * micro_mg_berg_eff_factor

     !+++PMC 12/3/12 - NEW VAPOR DEP/SUBLIMATION GOES HERE!!!
     if (do_cldice) then

        call ice_deposition_sublimation(t(:,k), q(:,k), qi(:,k), ni(:,k), &
             icldm(:,k), rho(:,k), dv(:,k), qvl(:,k), qvi(:,k),           &
             berg(:,k), vap_dep(:,k), ice_sublim(:,k), mgncol)

!     if(lprnt) write(0,*)' t=',t(1,k),' k=',k,' q=',q(1,k),' qi=',qi(1,k),&
!       ' ni=',ni(1,k),' icldm=',icldm(1,k),' rho=',rho(1,k),' dv=',dv(1,k),&
!       ' qvl=',qvl(1,k),' qvi=',qvi(1,k),' berg=',berg(1,k),' vap_dep=',&
!       vap_dep(1,k),' ice_sublim=',ice_sublim(1,k)
!     if(lprnt) write(0,*)' berg1=',berg(1,k),' k=',k,' micro_mg_berg_eff_factor=',micro_mg_berg_eff_factor
        do i=1,mgncol
! sublimation should not exceed available ice
          ice_sublim(i,k) = max(ice_sublim(i,k), -qi(i,k)*oneodt)
          berg(i,k) = berg(i,k) * micro_mg_berg_eff_factor
          if (ice_sublim(i,k) < zero .and. qi(i,k) > qsmall .and. icldm(i,k) > mincld) then
            nsubi(i,k) = sublim_factor * ice_sublim(i,k) * ni(i,k) / (qi(i,k) * icldm(i,k))
          else
            nsubi(i,k) = zero
          endif

        ! bergeron process should not reduce nc unless
        ! all ql is removed (which is handled elsewhere)
        !in fact, nothing in this entire file makes nsubc nonzero.
          nsubc(i,k) = zero
        end do

     end if !do_cldice
     !---PMC 12/3/12

!++ag Process rate calls for graupel here.
!   (Should this be in do_cldice loop?)
!===================================================================

     if(do_hail .or. do_graupel) then
       call graupel_collecting_snow(qsic(:,k),qric(:,k),umr(:,k),ums(:,k),        &
             rho(:,k),lamr(:,k),n0r(:,k),lams(:,k),n0s(:,k), psacr(:,k), mgncol)

       call graupel_collecting_cld_water(qgic(:,k),qcic(:,k),ncic(:,k),rho(:,k),  &
             n0g(:,k),lamg(:,k),bgtmp,agn(:,k), psacwg(:,k), npsacwg(:,k), mgncol)

       call graupel_riming_liquid_snow(psacws(:,k),qsic(:,k),qcic(:,k),nsic(:,k), &
             rho(:,k),rhosn,rhogtmp,asn(:,k),lams(:,k),n0s(:,k),deltat,           &
             pgsacw(:,k),nscng(:,k),mgncol)

!      if(lprnt .and. k >=100) then
!      if(lprnt) then
!        write(0,*)' k=',k,' qric=',qric(1,k),' qgic=',qgic(1,k),' umg=',umg(1,k),' umr=',umr(1,k),&
!        ' ung=',ung(1,k),' unr=',unr(1,k),' rho=',rho(1,k),' n0r=',n0r(1,k),' lamr=',lamr(1,k),&
!        ' n0g=',n0g(1,k),' lamg=',lamg(1,k),' pracg=',pracg(1,k)
!      endif
       call graupel_collecting_rain(qric(:,k),qgic(:,k),umg(:,k),                 &
             umr(:,k),ung(:,k),unr(:,k),rho(:,k),n0r(:,k),lamr(:,k),n0g(:,k),     &
             lamg(:,k), pracg(:,k),npracg(:,k),mgncol)
!      if(lprnt .and. k >=100) write(0,*)' k=',k,' pracg=',pracg(1,k),' npracg=',npracg(1,k)

!AG note: Graupel rain riming snow changes
!    pracs, npracs, (accretion of rain by snow)  psacr (collection of snow by rain)

!      if (lprnt .and. abs(k-81) <5) &
!               write(0,*)' k=',k,' pracs=',pracs(1,k),' npracs=',npracs(1,k),' psacr=',psacr(1,k),&
!       ' qsic=',qsic(1,k),' qric=',qric(1,k),' nric=',nric(1,k),' nsic=',nsic(1,k), &
!       ' n0s=',n0s(1,k),' lams=',lams(1,k),' n0r=',n0r(1,k),' lamr=',lamr(1,k),       &
!       ' pgracs=',pgracs(1,k),' ngracs=',ngracs(1,k)

       call graupel_rain_riming_snow(pracs(:,k),npracs(:,k),psacr(:,k),qsic(:,k), &
             qric(:,k),nric(:,k),nsic(:,k),n0s(:,k),lams(:,k),n0r(:,k),lamr(:,k), &
             deltat,pgracs(:,k),ngracs(:,k),mgncol)
!      if (lprnt .and. abs(k-81) <5) &
!               write(0,*)' k=',k,' pracs=',pracs(1,k),' npracs=',npracs(1,k),' psacr=',psacr(1,k),&
!       ' pgracs=',pgracs(1,k),' ngracs=',ngracs(1,k)

       call graupel_rime_splintering(t(:,k),qcic(:,k),qric(:,k),qgic(:,k),  &
             psacwg(:,k),pracg(:,k),qmultg(:,k),nmultg(:,k),qmultrg(:,k),   &
             nmultrg(:,k),mgncol)

!      if(lprnt .and. k >=100) write(0,*)' k=',k,' pracg2=',pracg(1,k)
!      if (lprnt .and. abs(k-81) <5) &
!      write(0,*)' k=',k,' pracg2=',pracg(1,k)

       call evaporate_sublimate_precip_graupel(t(:,k), rho(:,k),                &
             dv(:,k), mu(:,k), sc(:,k), q(:,k), qvl(:,k), qvi(:,k),             &
             lcldm(:,k), precip_frac(:,k), arn(:,k), asn(:,k), agn(:,k), bgtmp, &
             qcic(:,k), qiic(:,k), qric(:,k), qsic(:,k), qgic(:,k),             &
             lamr(:,k), n0r(:,k), lams(:,k), n0s(:,k), lamg(:,k), n0g(:,k),     &
             pre(:,k), prds(:,k), prdg(:,k), am_evp_st(:,k), mgncol)

!!Not used: part of above
!!        call graupel_sublimate_evap(t(:,k),q(:,k),qgic(:,k),rho(:,k),n0g(:,k), &
!!             lamg(:,k),qvi(:,k),dv(:,k),mu(:,k),sc(:,k),bgtmp,agn(:,k), &
!!             prdg(:,k),eprdg(:,k),mgncol)

!Checks for Debugging

!        if (minval(qmultg(:,k)).lt.0._r8) &
!             write(iulog,*) "OOPS, qmultg < 0 : min=",minval(qmultg(:,k))
!
!        if (minval(qmultrg(:,k)).lt.0._r8) &
!             write(iulog,*) "OOPS, qmultrg < 0 : min=",minval(qmultrg(:,k))
!
!        if (minval(pgracs(:,k)).lt.0._r8) &
!             write(iulog,*) "OOPS, pgracs < 0 : min=",minval(pgracs(:,k))
!
!        if (minval(psacwg(:,k)).lt.0._r8) &
!             write(iulog,*) "OOPS, psacwg < 0 : min=",minval(psacwg(:,k))
!
!        if (minval(npsacwg(:,k)).lt.0._r8) &
!             write(iulog,*) "OOPS, npsacwg < 0 : min=",minval(npsacwg(:,k))
!
!        if (minval(pracg(:,k)).lt.0._r8) &
!             write(iulog,*) "OOPS, pracg < 0 : min=",minval(pracg(:,k))
!
!        if (maxval(prdg(:,k)).gt.0._r8) &
!             write(iulog,*) "OOPS, prdg > 0 : max=",maxval(prdg(:,k))
!
!        if (minval(nmultg(:,k)).lt.0._r8) &
!             write(iulog,*) "OOPS, nmultg < 0 : min=",minval(nmultg(:,k))
!
!        if (minval(nmultrg(:,k)).lt.0._r8) &
!             write(iulog,*) "OOPS, nmultrg < 0 : min=",minval(nmultrg(:,k))
!
!        if (minval(ngracs(:,k)).lt.0._r8) &
!             write(iulog,*) "OOPS, ngracs < 0 : min=",minval(ngracs(:,k))
!
!        if (minval(psacr(:,k)).lt.0._r8) &
!             write(iulog,*) "OOPS, psacr < 0 : min=",minval(psacr(:,k))
!
!        if (minval(nscng(:,k)).lt.0._r8) &
!             write(iulog,*) "OOPS, nscng < 0 : min=",minval(nscng(:,k))

     else
! Routine without Graupel (original)
        call evaporate_sublimate_precip(t(:,k), rho(:,k),                         &
          dv(:,k), mu(:,k), sc(:,k), q(:,k), qvl(:,k), qvi(:,k),                  &
          lcldm(:,k), precip_frac(:,k), arn(:,k), asn(:,k), qcic(:,k), qiic(:,k), &
          qric(:,k), qsic(:,k), lamr(:,k), n0r(:,k), lams(:,k), n0s(:,k),         &
          pre(:,k), prds(:,k), am_evp_st(:,k), mgncol)


     end if ! end do_graupel/hail loop
!--ag

     do i=1,mgncol

        ! conservation to ensure no negative values of cloud water/precipitation
        ! in case microphysical process rates are large
        !===================================================================

        ! note: for check on conservation, processes are multiplied by omsm
        ! to prevent problems due to round off error

        ! conservation of qc
        !-------------------------------------------------------------------

!++ag Add graupel tendencies for qc to equation ON
!       dum = ((prc(i,k)+pra(i,k)+mnuccc(i,k)+mnucct(i,k)+msacwi(i,k)+ &
!             psacws(i,k)+bergs(i,k))*lcldm(i,k)+berg(i,k))*deltat
        dum = ( (prc(i,k)    + pra(i,k)   + mnuccc(i,k) + mnucct(i,k) + msacwi(i,k)             &
              +  psacws(i,k) + bergs(i,k) + qmultg(i,k) + psacwg(i,k) + pgsacw(i,k))*lcldm(i,k) &
              +  berg(i,k) ) * deltat
!--ag

        if (dum > qc(i,k) .and. abs(dum) > qsmall) then
!++ag
!          ratio = qc(i,k)/deltat/((prc(i,k)+pra(i,k)+mnuccc(i,k)+mnucct(i,k)+ &
!                  msacwi(i,k)+psacws(i,k)+bergs(i,k))*lcldm(i,k)+berg(i,k))*omsm

           ratio = qc(i,k) / dum * omsm

           qmultg(i,k) = ratio * qmultg(i,k)
           psacwg(i,k) = ratio * psacwg(i,k)
           pgsacw(i,k) = ratio * pgsacw(i,k)
!--ag
           prc(i,k)    = ratio * prc(i,k)
           pra(i,k)    = ratio * pra(i,k)
           mnuccc(i,k) = ratio * mnuccc(i,k)
           mnucct(i,k) = ratio * mnucct(i,k)
           msacwi(i,k) = ratio * msacwi(i,k)
           psacws(i,k) = ratio * psacws(i,k)
           bergs(i,k)  = ratio * bergs(i,k)
           berg(i,k)   = ratio * berg(i,k)
           qcrat(i,k)  = ratio
        else
           qcrat(i,k) = one
        end if

!     if(lprnt) write(0,*)' bergs2=',bergs(1,k),' k=',k,' ratio=',ratio

        !PMC 12/3/12: ratio is also frac of step w/ liquid.
        !thus we apply berg for "ratio" of timestep and vapor
        !deposition for the remaining frac of the timestep.
        if (qc(i,k) >= qsmall) then
           vap_dep(i,k) = vap_dep(i,k) * (one-qcrat(i,k))
        end if

     end do

     do i=1,mgncol

        !=================================================================
        ! apply limiter to ensure that ice/snow sublimation and rain evap
        ! don't push conditions into supersaturation, and ice deposition/nucleation don't
        ! push conditions into sub-saturation
        ! note this is done after qc conservation since we don't know how large
        ! vap_dep is before then
        ! estimates are only approximate since other process terms haven't been limited
        ! for conservation yet

        ! first limit ice deposition/nucleation vap_dep + mnuccd

        dum1 = vap_dep(i,k) + mnuccd(i,k)
        if (dum1 > 1.e-20_r8) then
           dum = (q(i,k)-qvi(i,k))/(one + xxls_squared*qvi(i,k)/(cpp*rv*t(i,k)*t(i,k)))*oneodt
           dum = max(dum, zero)
           if (dum1 > dum) then
              ! Allocate the limited "dum" tendency to mnuccd and vap_dep
              ! processes. Don't divide by cloud fraction; these are grid-
              ! mean rates.
              dum1 = mnuccd(i,k) / (vap_dep(i,k)+mnuccd(i,k))
              mnuccd(i,k)  = dum*dum1
              vap_dep(i,k) = dum - mnuccd(i,k)
           end if
        end if

     end do

     do i=1,mgncol

        !===================================================================
        ! conservation of nc
        !-------------------------------------------------------------------
!++ag NEW ONE ON
!       dum = (nprc1(i,k)+npra(i,k)+nnuccc(i,k)+nnucct(i,k)+ &
!            npsacws(i,k)-nsubc(i,k))*lcldm(i,k)*deltat
        dum = (nprc1(i,k)   + npra(i,k)  + nnuccc(i,k) + nnucct(i,k) &
             + npsacws(i,k) - nsubc(i,k) + npsacwg(i,k))*lcldm(i,k)*deltat
!--ag

        if (dum > nc(i,k) .and. abs(dum) > qsmall) then
           ratio = nc(i,k) / dum * omsm
!++ag
           npsacwg(i,k) = ratio * npsacwg(i,k)
!--ag

           nprc1(i,k)   = ratio * nprc1(i,k)
           npra(i,k)    = ratio * npra(i,k)
           nnuccc(i,k)  = ratio * nnuccc(i,k)
           nnucct(i,k)  = ratio * nnucct(i,k)
           npsacws(i,k) = ratio * npsacws(i,k)
           nsubc(i,k)   = ratio * nsubc(i,k)
        end if

        mnuccri(i,k) = zero
        nnuccri(i,k) = zero

        if (do_cldice) then

           ! freezing of rain to produce ice if mean rain size is smaller than Dcs
           if (lamr(i,k) > qsmall .and. one/lamr(i,k) < Dcs) then
              mnuccri(i,k) = mnuccr(i,k)
              nnuccri(i,k) = nnuccr(i,k)
              mnuccr(i,k)  = zero 
              nnuccr(i,k)  = zero
           end if
        end if

     end do

     do i=1,mgncol

        ! conservation of rain mixing ratio
        !-------------------------------------------------------------------
!++ag Implemented change for graupel
        dum1 = - pre(i,k)     + pracs(i,k) + mnuccr(i,k) + mnuccri(i,k)  &
               + qmultrg(i,k) + pracg(i,k) + pgracs(i,k)
        dum3 = dum1 * precip_frac(i,k)
        dum2 = (pra(i,k)+prc(i,k))*lcldm(i,k)
        dum = (dum3 - dum2) * deltat
!--ag

        ! note that qrtend is included below because of instantaneous freezing/melt
        if (dum > qr(i,k) .and. dum1 >= qsmall .and. abs(dum3) > qsmall) then
           ratio = (qr(i,k)*oneodt + dum2) / dum3 * omsm
!++ag
           qmultrg(i,k) = ratio * qmultrg(i,k)
           pracg(i,k)   = ratio * pracg(i,k)
           pgracs(i,k)  = ratio * pgracs(i,k)
!--ag
           pre(i,k)     = ratio * pre(i,k)
           pracs(i,k)   = ratio * pracs(i,k)
           mnuccr(i,k)  = ratio * mnuccr(i,k)
           mnuccri(i,k) = ratio * mnuccri(i,k)
        end if

     end do

     do i=1,mgncol

        ! conservation of rain number
        !-------------------------------------------------------------------

        ! Add evaporation of rain number.
        if (pre(i,k) < zero) then
           dum        = max(-one, pre(i,k)*deltat/qr(i,k))
           nsubr(i,k) = dum*nr(i,k) * oneodt
        else
           nsubr(i,k) = zero
        end if

     end do

     do i=1,mgncol

!++ag IMplemented change for graupel
!       dum1 = (-nsubr(i,k)+npracs(i,k)+nnuccr(i,k)+nnuccri(i,k)-nragg(i,k))*precip_frac(i,k)
!             nprc(i,k)*lcldm(i,k))*deltat

        dum1 = (-nsubr(i,k)+npracs(i,k)+nnuccr(i,k)+nnuccri(i,k)-nragg(i,k)  &
                +npracg(i,k)+ngracs(i,k))*precip_frac(i,k)
        dum2 = nprc(i,k)*lcldm(i,k)
        dum = (dum1 - dum2) * deltat
!--ag

        if (dum > nr(i,k) .and. abs(dum1) > qsmall) then
           ratio = (nr(i,k)*oneodt + dum2) / dum1 * omsm

!++ag
           npracg(i,k)  = ratio * npracg(i,k)
           ngracs(i,k)  = ratio * ngracs(i,k)
!--ag
           nragg(i,k)   = ratio * nragg(i,k)
           npracs(i,k)  = ratio * npracs(i,k)
           nnuccr(i,k)  = ratio * nnuccr(i,k)
           nsubr(i,k)   = ratio * nsubr(i,k)
           nnuccri(i,k) = ratio * nnuccri(i,k)
        end if

     end do

     if (do_cldice) then

        do i=1,mgncol

           ! conservation of qi
           !-------------------------------------------------------------------

!++ag

           dum1 = (prci(i,k)+prai(i,k))*icldm(i,k)-ice_sublim(i,k)
!          dum2 = vap_dep(i,k)+berg(i,k)+mnuccd(i,k)                           &
!               + (mnuccc(i,k)+mnucct(i,k)+mnudep(i,k)+msacwi(i,k))*lcldm(i,k) &
!               + mnuccri(i,k)*precip_frac(i,k)
           dum2 = vap_dep(i,k)+berg(i,k)+mnuccd(i,k)                                       &
                + (mnuccc(i,k)+mnucct(i,k)+mnudep(i,k)+msacwi(i,k)+qmultg(i,k))*lcldm(i,k) &
                + (qmultrg(i,k)+mnuccri(i,k))*precip_frac(i,k)
           dum  = (dum1 - dum2) * deltat
!-ag

           if (dum > qi(i,k) .and. abs(dum1) > qsmall) then
              ratio = (qi(i,k)*oneodt + dum2) / dum1 * omsm

!++ag
! Only sink terms are limited.
!             qmultg(i,k)     = ratio * qmultg(i,k)
!             qmultrg(i,k)    = ratio * qmultrg(i,k)
!--ag
              prci(i,k)       = ratio * prci(i,k)
              prai(i,k)       = ratio * prai(i,k)
              ice_sublim(i,k) = ratio * ice_sublim(i,k)
           end if

        end do

     end if

     if (do_cldice) then

        do i=1,mgncol

           ! conservation of ni
           !-------------------------------------------------------------------
           if (use_hetfrz_classnuc) then
              tmpfrz = nnuccc(i,k)
           else
              tmpfrz = zero
           end if
!++ag
           dum1 = (nprci(i,k)+nprai(i,k)-nsubi(i,k))*icldm(i,k)
!          dum2 = nnuccd(i,k)+(nnucct(i,k)+tmpfrz+nnudep(i,k)+nsacwi(i,k))*lcldm(i,k) &
!                            + nnuccri(i,k)*precip_frac(i,k)
           dum2 = nnuccd(i,k)+(nnucct(i,k)+tmpfrz+nnudep(i,k)+nsacwi(i,k)+nmultg(i,k))*lcldm(i,k) &
                             + (nmultrg(i,k)+nnuccri(i,k))*precip_frac(i,k)
!--ag
           dum  = (dum1 - dum2) * deltat

           if (dum > ni(i,k) .and. abs(dum1) > qsmall) then
              ratio = (ni(i,k)*oneodt + dum2) / dum1 * omsm

              nprci(i,k) = ratio * nprci(i,k)
              nprai(i,k) = ratio * nprai(i,k)
              nsubi(i,k) = ratio * nsubi(i,k)
           end if

        end do

     end if

     do i=1,mgncol

        ! conservation of snow mixing ratio
        !-------------------------------------------------------------------
!++ag
        if (do_hail .or. do_graupel) then
!NOTE: mnuccr is moved to graupel when active
!psacr is a positive value, but a loss for snow
!HM: psacr is positive in dum (two negatives)

          dum1 = (psacr(i,k) - prds(i,k)) * precip_frac(i,k)
          dum2 = pracs(i,k)*precip_frac(i,k)                                          &
               + (prai(i,k)+prci(i,k))*icldm(i,k) + (bergs(i,k)+psacws(i,k))*lcldm(i,k)
          dum  = (dum1 - dum2) * deltat
          if (dum > qs(i,k) .and. psacr(i,k)-prds(i,k) >= qsmall) then
            ratio = (qs(i,k)*oneodt + dum2) / dum1 * omsm
            psacr(i,k) = ratio * psacr(i,k)
            prds(i,k) = ratio * prds(i,k)
          endif
        else
          dum1 = - prds(i,k) * precip_frac(i,k)
          dum2 = (pracs(i,k)+mnuccr(i,k))*precip_frac(i,k)                            &
               + (prai(i,k)+prci(i,k))*icldm(i,k) + (bergs(i,k)+psacws(i,k))*lcldm(i,k)
          dum  = (dum1 - dum2) * deltat
          if (dum > qs(i,k) .and. -prds(i,k) >= qsmall) then
            ratio = (qs(i,k)*oneodt + dum2) / dum1 * omsm
            prds(i,k) = ratio * prds(i,k)
          endif
        endif

!--ag
!       dum1 = - prds(i,k) * precip_frac(i,k)
!       dum2 = (pracs(i,k)+mnuccr(i,k))*precip_frac(i,k)                            &
!            + (prai(i,k)+prci(i,k))*icldm(i,k) + (bergs(i,k)+psacws(i,k))*lcldm(i,k)

!       dum  = (dum1 - dum2) * deltat

!       if (dum > qs(i,k) .and. -prds(i,k) >= qsmall) then
!          ratio = (qs(i,k)*oneodt + dum2) / dum1 * omsm

!          prds(i,k) = ratio * prds(i,k)
!       end if

     end do

     do i=1,mgncol

       ! conservation of snow number
       !-------------------------------------------------------------------
       ! calculate loss of number due to sublimation
       ! for now neglect sublimation of ns
       nsubs(i,k) = zero

       ratio = one
!++ag Watch sign of nscng and ngracs. What is sign of nnuccr? Negative? Should be a source here.

       if (do_hail .or. do_graupel) then
!        dum1 = precip_frac(i,k)* (-nsubs(i,k)-nsagg(i,k)+ngracs(i,k))
!        dum2 = nprci(i,k)*icldm(i,k) + nscng(i,k)*lcldm(i,k)
!        dum  = (dum1 - dum2) * deltat
! check here -  this is slightly different from ag version - moorthi

         dum1 = precip_frac(i,k)* (-nsubs(i,k)-nsagg(i,k)+ngracs(i,k)) &
              - nscng(i,k)*lcldm(i,k)
         dum2 = nprci(i,k)*icldm(i,k)
         dum  = (dum1 - dum2) * deltat

         if (dum > ns(i,k) .and. abs(dum1) > qsmall) then
           ratio = (ns(i,k)*oneodt + dum2) / dum1 * omsm
           nscng(i,k)  = ratio * nscng(i,k)
           ngracs(i,k) = ratio * ngracs(i,k)
         end if

       else
         dum1 = precip_frac(i,k)* (-nsubs(i,k)-nsagg(i,k))
         dum2 = nnuccr(i,k)*precip_frac(i,k) + nprci(i,k)*icldm(i,k)
         dum  = (dum1 - dum2) * deltat

         if (dum > ns(i,k) .and. abs(dum1) > qsmall) then
           ratio = (ns(i,k)*oneodt + dum2) / dum1 * omsm
         end if
       endif
       nsubs(i,k) = ratio * nsubs(i,k)
       nsagg(i,k) = ratio * nsagg(i,k)

     end do

!++ag Graupel Conservation Checks
!-------------------------------------------------------------------
     if (do_hail .or. do_graupel) then
! conservation of graupel mass
!-------------------------------------------------------------------
       do i=1,mgncol

         dum1 = -prdg(i,k) * precip_frac(i,k)
         dum2 = (pracg(i,k)+pgracs(i,k)+psacr(i,k)+mnuccr(i,k))*precip_frac(i,k) &
              + (psacwg(i,k)+pgsacw(i,k))*lcldm(i,k)
         dum  = (dum1 - dum2) * deltat

         if (dum > qg(i,k) .and. abs(dum1) > qsmall) then

! hm added
! note: prdg is always negative (like prds), so it needs to be subtracted in ratio
           ratio = (qg(i,k)*oneodt + dum2) / dum1 * omsm

           prdg(i,k) = ratio * prdg(i,k)

         end if

       end do

! conservation of graupel number: not needed, no sinks
!-------------------------------------------------------------------
     end if
!--ag


     do i=1,mgncol

        ! next limit ice and snow sublimation and rain evaporation
        ! get estimate of q and t at end of time step
        ! don't include other microphysical processes since they haven't
        ! been limited via conservation checks yet

!++ag need to add graupel sublimation/evap here too (prdg)? May not need eprdg?
!++ag
        tx1 = pre(i,k)  * precip_frac(i,k)
        tx2 = prds(i,k) * precip_frac(i,k)
        tx6 = prdg(i,k) * precip_frac(i,k)
        tx5 = tx2 + tx6
        tx3 = tx1 + tx5 + ice_sublim(i,k)

        if (tx3 < -1.e-20_r8) then

           tx4 =  tx5 + ice_sublim(i,k) + vap_dep(i,k) + mnuccd(i,k)
           qtmp = q(i,k) - (tx1 + tx4) * deltat
           ttmp = t(i,k) + (tx1*xxlv + tx4*xxls) * (deltat/cpp)

           ! use rhw to allow ice supersaturation
           ! call qsat_water(ttmp, p(i,k), esn, qvn)
           esn = min(fpvsl(ttmp), p(i,k))
           qvn = epsqs*esn/(p(i,k)-omeps*esn) * qsfm(i,k)
!          qvn = epsqs*esn/(p(i,k)-omeps*esn)

           ! modify ice/precip evaporation rate if q > qsat
           if (qtmp > qvn) then

              tx4  = one / tx3
              dum1 = tx1 * tx4
              dum2 = tx2 * tx4
!++ag
              dum3 = tx6 * tx4
!--ag
              ! recalculate q and t after vap_dep and mnuccd but without evap or sublim
              tx5  = (vap_dep(i,k)+mnuccd(i,k)) * deltat
              qtmp = q(i,k) - tx5
              ttmp = t(i,k) + tx5 * (xxls/cpp)

              ! use rhw to allow ice supersaturation
              !call qsat_water(ttmp, p(i,k), esn, qvn)
              esn = min(fpvsl(ttmp), p(i,k))
              qvn = epsqs*esn / (p(i,k)-omeps*esn) * qsfm(i,k)
!             qvn = epsqs*esn / (p(i,k)-omeps*esn)

              dum = min(zero, (qtmp-qvn)/(one + xxlv_squared*qvn/(cpp*rv*ttmp*ttmp)))

              ! modify rates if needed, divide by precip_frac to get local (in-precip) value
              if (precip_frac(i,k) > mincld) then
                tx4 = oneodt / precip_frac(i,k)
              else
                tx4 = zero
              endif
              pre(i,k) = dum*dum1*tx4

              ! do separately using RHI for prds and ice_sublim
              !call qsat_ice(ttmp, p(i,k), esn, qvn)
              esn = min(fpvsi(ttmp), p(i,k))
              qvn = epsqs*esn / (p(i,k)-omeps*esn) * qsfm(i,k)
!             qvn = epsqs*esn / (p(i,k)-omeps*esn)


              dum = min(zero, (qtmp-qvn)/(one + xxls_squared*qvn/(cpp*rv*ttmp*ttmp)))

              ! modify rates if needed, divide by precip_frac to get local (in-precip) value
              prds(i,k) = dum*dum2*tx4
!++ag
              prdg(i,k) = dum*dum3*tx4
!--ag
!++ag
              ! don't divide ice_sublim by cloud fraction since it is grid-averaged
!             dum1 = one - dum1 - dum2
              dum1 = one - dum1 - dum2 - dum3
!--ag
              ice_sublim(i,k) = dum*dum1*oneodt
           end if
        end if

     end do

     ! Big "administration" loop enforces conservation, updates variables
     ! that accumulate over substeps, and sets output variables.

     do i=1,mgncol

        ! get tendencies due to microphysical conversion processes
        !==========================================================
        ! note: tendencies are multiplied by appropriate cloud/precip
        ! fraction to get grid-scale values
        ! note: vap_dep is already grid-average values

        ! The net tendencies need to be added to rather than overwritten,
        ! because they may have a value already set for instantaneous
        ! melting/freezing.

!++ag
!       qvlat(i,k) = qvlat(i,k) - (pre(i,k)+prds(i,k))*precip_frac(i,k)-&
!            vap_dep(i,k)-ice_sublim(i,k)-mnuccd(i,k)-mnudep(i,k)*lcldm(i,k)

        qvlat(i,k) = qvlat(i,k)-(pre(i,k)+prds(i,k))*precip_frac(i,k)                &
                    -vap_dep(i,k)-ice_sublim(i,k)-mnuccd(i,k)-mnudep(i,k)*lcldm(i,k) &
                    -prdg(i,k)*precip_frac(i,k)

!       tlat(i,k) = tlat(i,k) + ((pre(i,k)*precip_frac(i,k)) &
!            *xxlv+(prds(i,k)*precip_frac(i,k)+vap_dep(i,k)+ice_sublim(i,k)+mnuccd(i,k)+mnudep(i,k)*lcldm(i,k))*xxls+ &
!            ((bergs(i,k)+psacws(i,k)+mnuccc(i,k)+mnucct(i,k)+msacwi(i,k))*lcldm(i,k)+(mnuccr(i,k)+ &
!            pracs(i,k)+mnuccri(i,k))*precip_frac(i,k)+berg(i,k))*xlf)

!     if (lprnt .and. k >= 60 .and. k <=65) &
!     if (lprnt .and. k >= 100 ) &
!     if (lprnt .and. abs(k-81) <5) &
!     if (lprnt .and. k >= 60 ) &
!               write(0,*)' k=',k,' tlat=',tlat(i,k),' pre=',pre(i,k),' precip_frac=',precip_frac(i,k),&
!       ' prds=',prds(i,k),' prdg=',prdg(i,k),' vap_dep=',vap_dep(i,k),' ice_sublim=',ice_sublim(i,k), &
!       ' mnuccd=',mnuccd(i,k),' mnudep=',mnudep(i,k),' lcldm=',lcldm(i,k),' bergs=',bergs(i,k),       &
!       ' psacws=',psacws(i,k),' mnuccc=',mnuccc(i,k),' mnucct=',mnucct(i,k),' msacwi=',msacwi(i,k),   &
!       ' psacwg=',psacwg(i,k),' qmultg=',qmultg(i,k),' pgsacw=',pgsacw(i,k),' mnuccr=',mnuccr(i,k),   &
!       ' pracs=',pracs(i,k),' mnuccri=',mnuccri(i,k),' pracg=',pracg(i,k),' pgracs=',pgracs(i,k),     &
!       ' qmultrg=',qmultrg(i,k),' xlf=',xlf,' xxlv=',xxlv,' xxls=',xxls


         tlat(i,k) = tlat(i,k)+((pre(i,k)*precip_frac(i,k))*xxlv+ &
                   ((prds(i,k)+prdg(i,k))*precip_frac(i,k)+vap_dep(i,k)+ice_sublim(i,k)+ &
                     mnuccd(i,k)+mnudep(i,k)*lcldm(i,k))*xxls+ &
                   ((bergs(i,k)+psacws(i,k)+mnuccc(i,k)+mnucct(i,k)+msacwi(i,k)+psacwg(i,k)+ &
                     qmultg(i,k)+pgsacw(i,k))*lcldm(i,k)+ &
                   (mnuccr(i,k)+pracs(i,k)+mnuccri(i,k)+pracg(i,k)+pgracs(i,k)+qmultrg(i,k))*precip_frac(i,k)+ &
                     berg(i,k))*xlf)

!     if (lprnt .and. k >= 100 ) write(0,*)' k=',k,' tlat=',tlat(i,k)
!     if (lprnt)  write(0,*)' k=',k,' tlat=',tlat(i,k)
!     if (lprnt .and. k >= 60) write(0,*)' k=',k,' tlat=',tlat(i,k)

!      qctend(i,k) = qctend(i,k) + (-pra(i,k)-prc(i,k)-mnuccc(i,k)-mnucct(i,k)-msacwi(i,k)- &
!                                    psacws(i,k)-bergs(i,k))*lcldm(i,k)-berg(i,k)

       qctend(i,k) = qctend(i,k) +                                     &
             (-pra(i,k)-prc(i,k)-mnuccc(i,k)-mnucct(i,k)-msacwi(i,k) - &
             psacws(i,k)-bergs(i,k)-qmultg(i,k)-psacwg(i,k)-pgsacw(i,k))*lcldm(i,k)-berg(i,k)

        if (do_cldice) then
!          qitend(i,k) = qitend(i,k) + &
!               (mnuccc(i,k)+mnucct(i,k)+mnudep(i,k)+msacwi(i,k))*lcldm(i,k)+(-prci(i,k)- &
!               prai(i,k))*icldm(i,k)+vap_dep(i,k)+berg(i,k)+ice_sublim(i,k)+ &
!               mnuccd(i,k)+mnuccri(i,k)*precip_frac(i,k)

           qitend(i,k) = qitend(i,k)+                                                      &
                (mnuccc(i,k)+mnucct(i,k)+mnudep(i,k)+msacwi(i,k)+qmultg(i,k)) * lcldm(i,k) &
               + (-prci(i,k)-prai(i,k)) * icldm(i,k)                                       &
               + vap_dep(i,k)+berg(i,k)+ice_sublim(i,k)+mnuccd(i,k)                        &
               + (mnuccri(i,k)+qmultrg(i,k)) * precip_frac(i,k)

        end if

!       qrtend(i,k) = qrtend(i,k) + (pra(i,k)+prc(i,k))*lcldm(i,k)+(pre(i,k)-pracs(i,k)- &
!                                    mnuccr(i,k)-mnuccri(i,k))*precip_frac(i,k)

        qrtend(i,k) = qrtend(i,k) + (pra(i,k)+prc(i,k))*lcldm(i,k)+(pre(i,k)-pracs(i,k)- &
             mnuccr(i,k)-mnuccri(i,k)-qmultrg(i,k)-pracg(i,k)-pgracs(i,k))*precip_frac(i,k)


!       qstend(i,k) = qstend(i,k) + (prai(i,k)+prci(i,k))*icldm(i,k)+(psacws(i,k)+bergs(i,k))*lcldm(i,k) &
!                                 + (prds(i,k)+pracs(i,k)+mnuccr(i,k))*precip_frac(i,k)

        if (do_hail.or.do_graupel) then
           qgtend(i,k) = qgtend(i,k) + (pracg(i,k)+pgracs(i,k)+prdg(i,k)+psacr(i,k)+mnuccr(i,k))*precip_frac(i,k) &
                                     + (psacwg(i,k)+pgsacw(i,k))*lcldm(i,k)

           qstend(i,k) = qstend(i,k) + (prai(i,k)+prci(i,k))*icldm(i,k)+(psacws(i,k)+bergs(i,k))*lcldm(i,k) &
                                     + (prds(i,k)+pracs(i,k)-psacr(i,k))*precip_frac(i,k)

        else
           !necessary since mnuccr moved to graupel
           qstend(i,k) = qstend(i,k) + (prai(i,k)+prci(i,k))*icldm(i,k)+(psacws(i,k)+bergs(i,k))*lcldm(i,k) &
                                     + (prds(i,k)+pracs(i,k)+mnuccr(i,k))*precip_frac(i,k)

        end if
!--ag


        cmeout(i,k) = vap_dep(i,k) + ice_sublim(i,k) + mnuccd(i,k)

        ! add output for cmei (accumulate)
        cmeitot(i,k) = vap_dep(i,k) + ice_sublim(i,k) + mnuccd(i,k)

        ! assign variables for trop_mozart, these are grid-average
        !-------------------------------------------------------------------
        ! evaporation/sublimation is stored here as positive term

!++add evaporation/sublimation of graupel too? YES: After conservation checks.

!++ag
!ADD GRAUPEL to evapsnow: prdg. (sign? same as prds: negative, so this is a positive number)
!       evapsnow(i,k)  = -prds(i,k) * precip_frac(i,k)
        evapsnow(i,k)  = (-prds(i,k)-prdg(i,k)) * precip_frac(i,k)
!--ag
        nevapr(i,k)    = -pre(i,k)  * precip_frac(i,k)
        prer_evap(i,k) = -pre(i,k)  * precip_frac(i,k)

        ! change to make sure prain is positive: do not remove snow from
        ! prain used for wet deposition

!++AG NEED TO MAKE CONSISTENT WITH BUDGETS
        prain(i,k) = (pra(i,k)+prc(i,k))*lcldm(i,k)                                              &
                   - (pracs(i,k)+mnuccr(i,k)+mnuccri(i,k))*precip_frac(i,k)
        if (do_hail .or. do_graupel) then
! Subtract PSACR here or not? Ask Hugh
           prodsnow(i,k) = (prai(i,k)+prci(i,k))*icldm(i,k)+(psacws(i,k)+bergs(i,k))*lcldm(i,k)+ &
                            pracs(i,k)*precip_frac(i,k)
        else
           prodsnow(i,k) = (prai(i,k)+prci(i,k))*icldm(i,k)+(psacws(i,k)+bergs(i,k))*lcldm(i,k)+ &
                           (pracs(i,k)+mnuccr(i,k))*precip_frac(i,k)
        end if

        ! following are used to calculate 1st order conversion rate of cloud water
        !    to rain and snow (1/s), for later use in aerosol wet removal routine
        ! previously, wetdepa used (prain/qc) for this, and the qc in wetdepa may be smaller than the qc
        !    used to calculate pra, prc, ... in this routine
        ! qcsinksum_rate1ord = { rate of direct transfer of cloud water to rain & snow }
        !                      (no cloud ice or bergeron terms)

!++AG NEED TO MAKE CONSITANT: PGSACW, PSACWG (check budgets)? More sink terms? Check. No. Just loss to precip.
!Ask Hugh
!       qcsinksum_rate1ord(i,k) = (pra(i,k)+prc(i,k)+psacws(i,k))*lcldm(i,k)
        qcsinksum_rate1ord(i,k) = (pra(i,k)+prc(i,k)+psacws(i,k)+psacwg(i,k)+pgsacw(i,k))*lcldm(i,k)
!--ag
        ! Avoid zero/near-zero division.
        qcsinksum_rate1ord(i,k) = qcsinksum_rate1ord(i,k) / max(qc(i,k),1.0e-30_r8)


        ! microphysics output, note this is grid-averaged
        pratot(i,k)     = pra(i,k)     * lcldm(i,k)
        prctot(i,k)     = prc(i,k)     * lcldm(i,k)
        mnuccctot(i,k)  = mnuccc(i,k)  * lcldm(i,k)
        mnuccttot(i,k)  = mnucct(i,k)  * lcldm(i,k)
        msacwitot(i,k)  = msacwi(i,k)  * lcldm(i,k)
        psacwstot(i,k)  = psacws(i,k)  * lcldm(i,k)
        bergstot(i,k)   = bergs(i,k)   * lcldm(i,k)
        bergtot(i,k)    = berg(i,k)
        prcitot(i,k)    = prci(i,k)    * icldm(i,k)
        praitot(i,k)    = prai(i,k)    * icldm(i,k)
        mnuccdtot(i,k)  = mnuccd(i,k)  * icldm(i,k)

        pracstot(i,k)   = pracs(i,k)   * precip_frac(i,k)
        mnuccrtot(i,k)  = mnuccr(i,k)  * precip_frac(i,k)
!++ag
        mnuccritot(i,k) = mnuccri(i,k) * precip_frac(i,k)
!--ag

!++ag Hail/Graupel tendencies for output
        psacrtot(i,k)   = psacr(i,k)   * precip_frac(i,k)
        pracgtot(i,k)   = pracg(i,k)   * precip_frac(i,k)
        psacwgtot(i,k)  = psacwg(i,k)  * lcldm(i,k)
        pgsacwtot(i,k)  = pgsacw(i,k)  * lcldm(i,k)
        pgracstot(i,k)  = pgracs(i,k)  * precip_frac(i,k)
        prdgtot(i,k)    = prdg(i,k)    * precip_frac(i,k)
        qmultgtot(i,k)  = qmultg(i,k)  * lcldm(i,k)
        qmultrgtot(i,k) = qmultrg(i,k) * precip_frac(i,k)
        npracgtot(i,k)  = npracg(i,k)  * precip_frac(i,k)
        nscngtot(i,k)   = nscng(i,k)   * lcldm(i,k)
        ngracstot(i,k)  = ngracs(i,k)  * precip_frac(i,k)
        nmultgtot(i,k)  = nmultg(i,k)  * lcldm(i,k)
        nmultrgtot(i,k) = nmultrg(i,k) * precip_frac(i,k)
        npsacwgtot(i,k) = npsacwg(i,k) * lcldm(i,k)
!--ag

!++ag
!       nctend(i,k) = nctend(i,k) + (-nnuccc(i,k)-nnucct(i,k)-npsacws(i,k)+nsubc(i,k) &
!                                   - npra(i,k)-nprc1(i,k))*lcldm(i,k)

        nctend(i,k) = nctend(i,k) + (-nnuccc(i,k)-nnucct(i,k)-npsacws(i,k)+nsubc(i,k) &
                                     -npra(i,k)-nprc1(i,k)-npsacwg(i,k))*lcldm(i,k)

        if (do_cldice) then
           if (use_hetfrz_classnuc) then
              tmpfrz = nnuccc(i,k)
           else
              tmpfrz = zero
           end if
!          nitend(i,k) = nitend(i,k) + nnuccd(i,k)+ &
!               (nnucct(i,k)+tmpfrz+nnudep(i,k)+nsacwi(i,k))*lcldm(i,k)+(nsubi(i,k)-nprci(i,k)- &
!               nprai(i,k))*icldm(i,k)+nnuccri(i,k)*precip_frac(i,k)

          nitend(i,k) = nitend(i,k) + nnuccd(i,k)                                     &
               +  (nnucct(i,k)+tmpfrz+nnudep(i,k)+nsacwi(i,k)+nmultg(i,k))*lcldm(i,k) &
               + (nsubi(i,k)-nprci(i,k)-nprai(i,k))*icldm(i,k)                        &
               + (nnuccri(i,k)+nmultrg(i,k))*precip_frac(i,k)
        end if

        if(do_graupel.or.do_hail) then
!       nstend(i,k) = nstend(i,k) + (nsubs(i,k)+nsagg(i,k)+nnuccr(i,k))*precip_frac(i,k) &
!                                 + nprci(i,k)*icldm(i,k)

           nstend(i,k) = nstend(i,k) + (nsubs(i,k)+nsagg(i,k)-ngracs(i,k))*precip_frac(i,k) &
                                     + nprci(i,k)*icldm(i,k)-nscng(i,k)*lcldm(i,k)

           ngtend(i,k) = ngtend(i,k) + nscng(i,k)*lcldm(i,k)                                &
                                     + (ngracs(i,k)+nnuccr(i,k))*precip_frac(i,k)

        else
           !necessary since mnuccr moved to graupel
           nstend(i,k) = nstend(i,k) + (nsubs(i,k)+nsagg(i,k)+nnuccr(i,k))*precip_frac(i,k) &
                                     + nprci(i,k)*icldm(i,k)

        end if

!       nrtend(i,k) = nrtend(i,k) + nprc(i,k)*lcldm(i,k)+(nsubr(i,k)-npracs(i,k)-nnuccr(i,k) &
!                                 - nnuccri(i,k)+nragg(i,k))*precip_frac(i,k)

        nrtend(i,k) = nrtend(i,k)+ nprc(i,k)*lcldm(i,k)                                     &
                                 + (nsubr(i,k)-npracs(i,k)-nnuccr(i,k)                      &
                                   -nnuccri(i,k)+nragg(i,k)-npracg(i,k)-ngracs(i,k))*precip_frac(i,k)
!--ag

        ! make sure that ni at advanced time step does not exceed
        ! maximum (existing N + source terms*dt), which is possible if mtime < deltat
        ! note that currently mtime = deltat
        !================================================================

        if (do_cldice .and. nitend(i,k) > zero .and. ni(i,k)+nitend(i,k)*deltat > nimax(i,k)) then
           nitend(i,k) = max(zero, (nimax(i,k)-ni(i,k))*oneodt)
        end if

     end do

     ! End of "administration" loop

  end do micro_vert_loop ! end k loop

! if (lprnt) write(0,*)' tlat3=',tlat(1,:)*deltat
  !-----------------------------------------------------
  ! convert rain/snow q and N for output to history, note,
  ! output is for gridbox average

  do k=1,nlev
    do i=1,mgncol
      qrout(i,k) = qr(i,k)
      nrout(i,k) = nr(i,k) * rho(i,k)
      qsout(i,k) = qs(i,k)
      nsout(i,k) = ns(i,k) * rho(i,k)
!++ag
      qgout(i,k) = qg(i,k)
      ngout(i,k) = ng(i,k) * rho(i,k)
!--ag
    enddo
  enddo

  ! calculate n0r and lamr from rain mass and number
  ! divide by precip fraction to get in-precip (local) values of
  ! rain mass and number, divide by rhow to get rain number in kg^-1

  do k=1,nlev

     call size_dist_param_basic(mg_rain_props, qric(:,k), nric(:,k), lamr(:,k), mgncol, n0=n0r(:,k))

  enddo
     ! Calculate rercld

     ! calculate mean size of combined rain and cloud water

     call calc_rercld(lamr, n0r, lamc, pgam, qric, qcic, ncic, rercld, mgncol, nlev)


  ! Assign variables back to start-of-timestep values
  ! Some state variables are changed before the main microphysics loop
  ! to make "instantaneous" adjustments. Afterward, we must move those changes
  ! back into the tendencies.
  ! These processes:
  !  - Droplet activation (npccn, impacts nc)
  !  - Instantaneous snow melting  (minstsm/ninstsm, impacts qr/qs/nr/ns)
  !  - Instantaneous rain freezing (minstfr/ninstrf, impacts qr/qs/nr/ns)
  !================================================================================

  do k=1,nlev
    do i=1,mgncol
  ! Re-apply droplet activation tendency
      nc(i,k)     = ncn(i,k)
      nctend(i,k) = nctend(i,k) + npccn(i,k)

  ! Re-apply rain freezing and snow melting.
      qstend(i,k) = qstend(i,k) + (qs(i,k)-qsn(i,k)) * oneodt
      qs(i,k)     = qsn(i,k)

      nstend(i,k) = nstend(i,k) + (ns(i,k)-nsn(i,k)) * oneodt
      ns(i,k)     = nsn(i,k)

      qrtend(i,k) = qrtend(i,k) + (qr(i,k)-qrn(i,k)) * oneodt
      qr(i,k)     = qrn(i,k)

      nrtend(i,k) = nrtend(i,k) + (nr(i,k)-nrn(i,k)) * oneodt
      nr(i,k)     = nrn(i,k)

!++ag Re-apply graupel freezing/melting
      qgtend(i,k) = qgtend(i,k) + (qg(i,k)-qgr(i,k)) * oneodt
      qg(i,k)     = qgr(i,k)

!  if (maxval(dum_2D-qg).gt.0._r8) &
!       write(iulog,*) "OOPS, qg diff > 0 : max=",maxval(dum_2D-qg)
!  if (minval(dum_2D-qg).lt.0._r8) &
!       write(iulog,*) "OOPS, qg diff < 0 : min=",minval(dum_2D-qg)
!
!  write(iulog,*) "Max qgtend: 1st = ",maxval(qgtend)
!  write(iulog,*) "Min qgtend: 1st = ",minval(qgtend)
!  write(iulog,*) "Max qvtend: 1st = ",maxval(qvlat)
!  write(iulog,*) "Min qvtend: 1st = ",minval(qvlat)

      ngtend(i,k) = ngtend(i,k) + (ng(i,k)-ngr(i,k)) * oneodt
      ng(i,k)     = ngr(i,k)
!--ag

  !.............................................................................

  !================================================================================

  ! modify to include snow. in prain & evap (diagnostic here: for wet dep)
      nevapr(i,k) = nevapr(i,k) + evapsnow(i,k)
      prain(i,k)  = prain(i,k)  + prodsnow(i,k)


    enddo
  enddo

  do k=1,nlev

     do i=1,mgncol

        ! calculate sedimentation for cloud water and ice
!++ag   ! and Graupel (mg3)
        !================================================================================

        ! update in-cloud cloud mixing ratio and number concentration
        ! with microphysical tendencies to calculate sedimentation, assign to dummy vars
        ! note: these are in-cloud values***, hence we divide by cloud fraction

        if (lcldm(i,k) > mincld) then
          tx1        = one / lcldm(i,k)
          dumc(i,k)  = (qc(i,k)+qctend(i,k)*deltat) * tx1
          dumnc(i,k) = max((nc(i,k)+nctend(i,k)*deltat)*tx1, zero)
          ! switch for specification of droplet and crystal number
          if (nccons) then
             dumnc(i,k) = ncnst*rhoinv(i,k)
          end if
        else
          dumc(i,k)  = zero
          dumnc(i,k) = zero
        endif
        if (icldm(i,k) > mincld) then
          tx1        = one / icldm(i,k)
          dumi(i,k)  = (qi(i,k)+qitend(i,k)*deltat) * tx1
          dumni(i,k) = max((ni(i,k)+nitend(i,k)*deltat)*tx1, zero)
          ! switch for specification of cloud ice number
          if (nicons) then
             dumni(i,k) = ninst*rhoinv(i,k)
          end if
        else
          dumi(i,k)  = zero
          dumni(i,k) = zero
        endif
        if (precip_frac(i,k) > mincld) then
          tx1        = one / precip_frac(i,k)
          dumr(i,k)  = (qr(i,k)+qrtend(i,k)*deltat) * tx1
          dums(i,k)  = (qs(i,k)+qstend(i,k)*deltat) * tx1

          dumnr(i,k) = max((nr(i,k)+nrtend(i,k)*deltat)*tx1, zero)
          dumns(i,k) = max((ns(i,k)+nstend(i,k)*deltat)*tx1, zero)

!++ag Add graupel
          dumg(i,k)  = (qg(i,k)+qgtend(i,k)*deltat) * tx1
!  Moorthi testing
          if (dumg(i,k) > 0.01_r8) then
            tx2 = dumg(i,k) - 0.01_r8
            dumg(i,k) = 0.01_r8
            dums(i,k) = dums(i,k) + tx2
            qstend(i,k) = (dums(i,k)*precip_frac(i,k) - qs(i,k)) * oneodt
            qgtend(i,k) = (dumg(i,k)*precip_frac(i,k) - qg(i,k)) * oneodt
          endif
!  Moorthi testing

          dumng(i,k) = max((ng(i,k)+ngtend(i,k)*deltat)*tx1, zero)
          ! switch for specification of droplet and crystal number
          if (ngcons) then
             dumng(i,k) = ngnst*rhoinv(i,k)
          endif
!--ag
        else
          dumr(i,k)  = zero
          dumr(i,k)  = zero
          dums(i,k)  = zero
          dumns(i,k) = zero
!++ag Add graupel
          dumg(i,k)  = zero
          dumng(i,k) = zero
        endif
!--ag
     enddo
  enddo

  do k=1,nlev

! obtain new slope parameter to avoid possible singularity

     call size_dist_param_ice(mg_ice_props, dumi(:,k), dumni(:,k), &
                                lami(:,k), mgncol)

     call size_dist_param_liq(mg_liq_props, dumc(:,k), dumnc(:,k), rho(:,k), &
                              pgam(:,k), lamc(:,k), mgncol)

!    call size_dist_param_basic(mg_ice_props, dumi(:,k), dumni(:,k), &
!                               lami(:,k), mgncol)
! fallspeed for rain

     call size_dist_param_basic(mg_rain_props, dumr(:,k), dumnr(:,k), &
                                lamr(:,k), mgncol)
! fallspeed for snow
     call size_dist_param_basic(mg_snow_props, dums(:,k), dumns(:,k), &
                                lams(:,k), mgncol)
! fallspeed for graupel/hail
     if (do_graupel .or. do_hail) then
       call size_dist_param_basic(mg_graupel_props, dumg(:,k), dumng(:,k), &
                                  lamg(:,k), mgncol)
     endif
  enddo

  do k=1,nlev
     do i=1,mgncol

        ! calculate number and mass weighted fall velocity for droplets and cloud ice
        !-------------------------------------------------------------------

        grho = g*rho(i,k)

        if (dumc(i,k) >= qsmall) then

           tx1 = lamc(i,k)**bc
           vtrmc(i,k) = acn(i,k)*gamma(pgam(i,k)+four+bc)     &
                          / (tx1*gamma(pgam(i,k)+four))

           fc(i,k)  = grho * vtrmc(i,k)
           fnc(i,k) = grho * acn(i,k)*gamma(pgam(i,k)+one+bc) &
                               / (tx1*gamma(pgam(i,k)+one))
        else
           fc(i,k)  = zero
           fnc(i,k) = zero
        end if

        ! calculate number and mass weighted fall velocity for cloud ice

        if (dumi(i,k) >= qsmall) then

           tx3 = one / lami(i,k)
           tx1 = ain(i,k) * tx3**bi
           tx2 = 1.2_r8*rhof(i,k)
           vtrmi(i,k) = min(tx1*gamma_bi_plus4*oneo6, tx2)

           fi(i,k)  = grho * vtrmi(i,k)
           fni(i,k) = grho * min(tx1*gamma_bi_plus1, tx2)

           ! adjust the ice fall velocity for smaller (r < 20 um) ice
           ! particles (blend over 18-20 um)
           irad = (1.5_r8 * 1e6_r8) * tx3
           ifrac = min(one, max(zero, (irad-18._r8)*half))
 
           if (ifrac < one) then
              tx1 = ajn(i,k) / lami(i,k)**bj
              vtrmi(i,k) = ifrac*vtrmi(i,k) +  (one-ifrac) * min(tx1*gamma_bj_plus4*oneo6, tx2)

              fi(i,k)  = grho * vtrmi(i,k)
              fni(i,k) = ifrac * fni(i,k) + (one-ifrac) * grho * min(tx1*gamma_bj_plus1, tx2)
           end if
        else
           fi(i,k) = zero
           fni(i,k)= zero
        end if

        ! fallspeed for rain

!       if (lamr(i,k) >= qsmall) then
        if (dumr(i,k) >= qsmall) then

           ! 'final' values of number and mass weighted mean fallspeed for rain (m/s)

           tx1 = arn(i,k) / lamr(i,k)**br
           tx2 = 9.1_r8*rhof(i,k)
           umr(i,k) = min(tx1*gamma_br_plus4*oneo6, tx2)
           unr(i,k) = min(tx1*gamma_br_plus1,       tx2)

           fr(i,k)  = grho * umr(i,k)
           fnr(i,k) = grho * unr(i,k)

        else
           fr(i,k)  = zero
           fnr(i,k) = zero
        end if

        ! fallspeed for snow

!       if (lams(i,k) >= qsmall) then
        if (dums(i,k) >= qsmall) then

           ! 'final' values of number and mass weighted mean fallspeed for snow (m/s)
           tx1 = asn(i,k) / lams(i,k)**bs
           tx2 = 1.2_r8*rhof(i,k)
           ums(i,k) = min(tx1*gamma_bs_plus4*oneo6, tx2)
           uns(i,k) = min(tx1*gamma_bs_plus1,       tx2)

           fs(i,k)  = grho * ums(i,k)
           fns(i,k) = grho * uns(i,k)

        else
           fs(i,k)  = zero
           fns(i,k) = zero
        end if

        if (do_graupel .or. do_hail) then
!++ag
        ! fallspeed for graupel


!         if (lamg(i,k) >= qsmall) then
          if (dumg(i,k) >= qsmall) then

           ! 'final' values of number and mass weighted mean fallspeed for graupel (m/s)
            tx1 = agn(i,k) / lamg(i,k)**bgtmp
            tx2 = 20._r8 * rhof(i,k)
            umg(i,k) = min(tx1*gamma_bg_plus4*oneo6, tx2)
            ung(i,k) = min(tx1*gamma_bg_plus1,       tx2)

            fg(i,k)  = grho * umg(i,k)
            fng(i,k) = grho * ung(i,k)

          else
            fg(i,k)  = zero
            fng(i,k) = zero
          end if
        endif

        ! redefine dummy variables - sedimentation is calculated over grid-scale
        ! quantities to ensure conservation

        dumc(i,k)  = qc(i,k) + qctend(i,k)*deltat
        dumi(i,k)  = qi(i,k) + qitend(i,k)*deltat
        dumr(i,k)  = qr(i,k) + qrtend(i,k)*deltat
        dums(i,k)  = qs(i,k) + qstend(i,k)*deltat

        dumnc(i,k) = nc(i,k) + nctend(i,k)*deltat
        dumni(i,k) = ni(i,k) + nitend(i,k)*deltat
        dumnr(i,k) = nr(i,k) + nrtend(i,k)*deltat
        dumns(i,k) = ns(i,k) + nstend(i,k)*deltat
!++ag
        dumg(i,k)  = qg(i,k) + qgtend(i,k)*deltat
        dumng(i,k) = ng(i,k) + ngtend(i,k)*deltat
!--ag

        if (dumc(i,k) < qsmall) dumnc(i,k) = zero
        if (dumi(i,k) < qsmall) dumni(i,k) = zero
        if (dumr(i,k) < qsmall) dumnr(i,k) = zero
        if (dums(i,k) < qsmall) dumns(i,k) = zero
        if (dumg(i,k) < qsmall) dumng(i,k) = zero

     enddo
  end do       !!! vertical loop

  do k=1,nlev
     do i=1,mgncol
       pdel_inv(i,k) = one / pdel(i,k)
     enddo
  enddo
! if (lprnt) write(0,*)' bef sedimentation dumc=',dumc(i,nlev-10:nlev)

  ! initialize nstep for sedimentation sub-steps

  ! calculate number of split time steps to ensure courant stability criteria
  ! for sedimentation calculations
  !-------------------------------------------------------------------
  do i=1,mgncol
     nlb = nlball(i)
     nstep = 1 + nint(max( maxval( fi(i,nlb:nlev)*pdel_inv(i,nlb:nlev)), &
                           maxval(fni(i,nlb:nlev)*pdel_inv(i,nlb:nlev))) * deltat)
     nstep = min(nstep, nstep_def)


     ! loop over sedimentation sub-time step to ensure stability
     !==============================================================
     if (do_cldice) then
       tx2 = one / nstep
       tx1 = tx2 * deltat
       tx3 = tx2 / g

       do n = 1,nstep

        ! top of model

          k = 1

        ! add fallout terms to microphysical tendencies

          tx5           = dumi(i,k)
          tx7           = pdel_inv(i,k) * tx1
          dumi(i,k)     = tx5 / (one + fi(i,k)*tx7)
          tx6           = (dumi(i,k)-tx5) * oneodt
          qitend(i,k)   = qitend(i,k) + tx6
          tx5           = dumni(i,k)
          dumni(i,k)    = tx5 / (one + fni(i,k)*tx7)
          nitend(i,k)   = nitend(i,k) + (dumni(i,k)-tx5) * oneodt

        ! sedimentation tendency for output
          qisedten(i,k) = qisedten(i,k) + tx6

          falouti(k)    = fi(i,k)  * dumi(i,k)
          faloutni(k)   = fni(i,k) * dumni(i,k)

          iflx(i,k+1)   = iflx(i,k+1) + falouti(k) * tx3   ! Ice flux

          do k = 2,nlev

          ! for cloud liquid and ice, if cloud fraction increases with height
          ! then add flux from above to both vapor and cloud water of current level
          ! this means that flux entering clear portion of cell from above evaporates
          ! instantly

          ! note: this is not an issue with precip, since we assume max overlap

            if (icldm(i,k-1) > mincld) then
              dum1 = max(zero, min(one, icldm(i,k)/icldm(i,k-1)))
            else
              dum1 = one
            endif

            tx5           = dumi(i,k)
            tx7           = pdel_inv(i,k) * tx1
            dum2          = tx7 * dum1
            dumi(i,k)     = (tx5 + falouti(k-1)*dum2) / (one + fi(i,k)*tx7)
            tx6           = (dumi(i,k)-tx5) * oneodt
                                                                   ! add fallout terms to eulerian tendencies
            qitend(i,k)   = qitend(i,k) + tx6
            tx5           = dumni(i,k)
            dumni(i,k)    = (tx5 + faloutni(k-1)*dum2)  / (one + fni(i,k)*tx7)
            nitend(i,k)   = nitend(i,k) + (dumni(i,k)-tx5) * oneodt


            qisedten(i,k) = qisedten(i,k) + tx6                    ! sedimentation tendency for output


            falouti(k)    = fi(i,k)  * dumi(i,k)
            faloutni(k)   = fni(i,k) * dumni(i,k)

            dum2 = (one-dum1) * falouti(k-1) * pdel_inv(i,k) * tx2
            qvlat(i,k)    = qvlat(i,k)   + dum2                    ! add terms to evap/sub of cloud ice
            qisevap(i,k)  = qisevap(i,k) + dum2                    ! for output

            tlat(i,k)     = tlat(i,k)    - dum2 * xxls

            iflx(i,k+1)   = iflx(i,k+1) + falouti(k) * tx3         ! Ice flux
          end do

        ! units below are m/s
        ! sedimentation flux at surface is added to precip flux at surface
        ! to get total precip (cloud + precip water) rate

          prect(i) = prect(i) + falouti(nlev) * (tx3*0.001_r8)
          preci(i) = preci(i) + falouti(nlev) * (tx3*0.001_r8)

       end do
     end if

! if (lprnt) write(0,*)' tlat4=',tlat(1,:)*deltat
     ! calculate number of split time steps to ensure courant stability criteria
     ! for sedimentation calculations
     !-------------------------------------------------------------------
     nstep = 1 + nint(max( maxval( fc(i,nlb:nlev)*pdel_inv(i,nlb:nlev)), &
                           maxval(fnc(i,nlb:nlev)*pdel_inv(i,nlb:nlev))) * deltat)
     nstep = min(nstep, nstep_def)

     ! loop over sedimentation sub-time step to ensure stability
     !==============================================================
     tx2 = one / nstep
     tx1 = tx2 * deltat
     tx3 = tx2 / g

     do n = 1,nstep

        ! top of model
        k = 1

        tx5           = dumc(i,k)
        tx7           = pdel_inv(i,k) * tx1
        dumc(i,k)     = tx5 / (one + fc(i,k)*tx7)
        tx6           = (dumc(i,k)-tx5) * oneodt
        qctend(i,k)   = qctend(i,k) + tx6
        tx5           = dumnc(i,k)
        dumnc(i,k)    = tx5 / (one + fnc(i,k)*tx7)
        nctend(i,k)   = nctend(i,k) + (dumnc(i,k)-tx5) * oneodt


        ! sedimentation tendency for output
        qcsedten(i,k) = qcsedten(i,k) + tx6

        faloutc(k)    = fc(i,k)  * dumc(i,k)
        faloutnc(k)   = fnc(i,k) * dumnc(i,k)

        lflx(i,k+1)   = lflx(i,k+1) + faloutc(k) * tx3
        do k = 2,nlev

          if (lcldm(i,k-1) > mincld) then
            dum1 = max(zero, min(one, lcldm(i,k)/lcldm(i,k-1)))
          else
            dum1 = one
          endif

          tx5           = dumc(i,k)
          tx7           = pdel_inv(i,k) * tx1
          dum2          = tx7 * dum1
          dumc(i,k)     = (tx5 + faloutc(k-1)*dum2) / (one + fc(i,k)*tx7)
          tx6           = (dumc(i,k)-tx5) * oneodt
          qctend(i,k)   = qctend(i,k) + tx6
          tx5           = dumnc(i,k)
          dumnc(i,k)    = (tx5 + faloutnc(k-1)*dum2)  / (one + fnc(i,k)*tx7)
          nctend(i,k)   = nctend(i,k) + (dumnc(i,k)-tx5) * oneodt



          qcsedten(i,k) = qcsedten(i,k) + tx6                ! sedimentation tendency for output

          faloutc(k)    = fc(i,k)  * dumc(i,k)
          faloutnc(k)   = fnc(i,k) * dumnc(i,k)

          dum2 = (one-dum1) * faloutc(k-1) * pdel_inv(i,k) * tx2
          qvlat(i,k)    = qvlat(i,k)   + dum2                ! add terms to to evap/sub of cloud water
          qcsevap(i,k)  = qcsevap(i,k) + dum2                ! for output

          tlat(i,k)     = tlat(i,k)    - dum2 * xxlv

          lflx(i,k+1)   = lflx(i,k+1) + faloutc(k) * tx3     ! Liquid condensate flux here
        end do

        prect(i) = prect(i) + faloutc(nlev) * (tx3*0.001_r8)

     end do
! if (lprnt) write(0,*)' tlat5=',tlat(1,:)*deltat
! if (lprnt) write(0,*)' maxval=',maxval( fr(i,nlb:nlev)*pdel_inv(i,nlb:nlev))&
!                                ,maxval(fnr(i,nlb:nlev)*pdel_inv(i,nlb:nlev))

     ! calculate number of split time steps to ensure courant stability criteria
     ! for sedimentation calculations
     !-------------------------------------------------------------------
     nstep = 1 + nint(max( maxval( fr(i,nlb:nlev)*pdel_inv(i,nlb:nlev)), &
                           maxval(fnr(i,nlb:nlev)*pdel_inv(i,nlb:nlev))) * deltat)

     nstep = min(nstep, nstep_def)

     ! loop over sedimentation sub-time step to ensure stability
     !==============================================================
     tx2 = one / nstep
     tx1 = tx2 * deltat
     tx3 = tx2 / g

!    if(lprnt) then
!      write(0,*)' nstep=',nstep,' tx1=',tx1,' tx2=',tx2,' tx3=',tx3,' qsmall=',qsmall
!      write(0,*)' fr=',fr(i,:)
!      write(0,*)' dumr=',dumr(i,:)
!    endif

     do n = 1,nstep

        ! top of model
        k = 1

        ! add fallout terms to microphysical tendencies

        tx5           = dumr(i,k)
        tx7           = pdel_inv(i,k) * tx1
        dumr(i,k)     = tx5 / (one + fr(i,k)*tx7)
        tx6           = (dumr(i,k)-tx5) * oneodt
        qrtend(i,k)   = qrtend(i,k) + tx6
        tx5           = dumnr(i,k)
        dumnr(i,k)    = tx5 / (one + fnr(i,k)*tx7)
        nrtend(i,k)   = nrtend(i,k) + (dumnr(i,k)-tx5) * oneodt

        ! sedimentation tendency for output
        qrsedten(i,k) = qrsedten(i,k) + tx6

        faloutr(k)    = fr(i,k)  * dumr(i,k)
        faloutnr(k)   = fnr(i,k) * dumnr(i,k)

        rflx(i,k+1)   = rflx(i,k+1) + faloutr(k) * tx3

        do k = 2,nlev

          tx5           = dumr(i,k)
          tx7           = pdel_inv(i,k) * tx1
          dumr(i,k)     = (tx5 + faloutr(k-1)*tx7) / (one + fr(i,k)*tx7)
          tx6           = (dumr(i,k)-tx5) * oneodt
          qrtend(i,k)   = qrtend(i,k) + tx6
          tx5           = dumnr(i,k)
          dumnr(i,k)    = (tx5 + faloutnr(k-1)*tx7)  / (one + fnr(i,k)*tx7)
          nrtend(i,k)   = nrtend(i,k) + (dumnr(i,k)-tx5) * oneodt

          qrsedten(i,k) = qrsedten(i,k) + tx6            ! sedimentation tendency for output

          faloutr(k)    = fr(i,k)  * dumr(i,k)
          faloutnr(k)   = fnr(i,k) * dumnr(i,k)

          rflx(i,k+1)   = rflx(i,k+1) + faloutr(k) * tx3 ! Rain Flux
        end do

        prect(i) = prect(i) + faloutr(nlev) * (tx3*0.001_r8)

     end do

!    if (lprnt) write(0,*)' prectaftrain=',prect(i),' preci=',preci(i)

     ! calculate number of split time steps to ensure courant stability criteria
     ! for sedimentation calculations
     !-------------------------------------------------------------------
     nstep = 1 + nint(max( maxval( fs(i,nlb:nlev)*pdel_inv(i,nlb:nlev)), &
                           maxval(fns(i,nlb:nlev)*pdel_inv(i,nlb:nlev))) * deltat)
     nstep = min(nstep, nstep_def)


     ! loop over sedimentation sub-time step to ensure stability
     !==============================================================
     tx2 = one / nstep
     tx1 = tx2 * deltat
     tx3 = tx2 / g
     do n = 1,nstep

        ! top of model
        k = 1

        ! add fallout terms to microphysical tendencies

        tx5           = dums(i,k)
        tx7           = pdel_inv(i,k) * tx1
        dums(i,k)     = tx5 / (one + fs(i,k)*tx7)
        tx6           = (dums(i,k)-tx5) * oneodt
        qstend(i,k)   = qstend(i,k) + tx6
        tx5           = dumns(i,k)
        dumns(i,k)    = tx5 / (one + fns(i,k)*tx7)
        nstend(i,k)   = nstend(i,k) + (dumns(i,k)-tx5) * oneodt

        ! sedimentation tendency for output
        qssedten(i,k) = qssedten(i,k) + tx6

        falouts(k)    = fs(i,k)  * dums(i,k)
        faloutns(k)   = fns(i,k) * dumns(i,k)

        sflx(i,k+1)   = sflx(i,k+1) + falouts(k) * tx3

        do k = 2,nlev


           tx5           = dums(i,k)
           tx7           = pdel_inv(i,k) * tx1
           dums(i,k)     = (tx5 + falouts(k-1)*tx7) / (one + fs(i,k)*tx7)
           tx6           = (dums(i,k)-tx5) * oneodt
           qstend(i,k)   = qstend(i,k) + tx6
           tx5           = dumns(i,k)
           dumns(i,k)    = (tx5 + faloutns(k-1)*tx7)  / (one + fns(i,k)*tx7)
           nstend(i,k)   = nstend(i,k) + (dumns(i,k)-tx5) * oneodt


           qssedten(i,k) = qssedten(i,k) + tx6                      ! sedimentation tendency for output

           falouts(k)    = fs(i,k)  * dums(i,k)
           faloutns(k)   = fns(i,k) * dumns(i,k)

           sflx(i,k+1)   = sflx(i,k+1) + falouts(k) * tx3 ! Snow Flux
        end do   !! k loop

        prect(i) = prect(i) + falouts(nlev) * (tx3*0.001_r8)
        preci(i) = preci(i) + falouts(nlev) * (tx3*0.001_r8)

     end do   !! nstep loop

!    if (lprnt) write(0,*)' prectaftssno=',prect(i),' preci=',preci(i)
!    if (lprnt) write(0,*)' qgtnd1=',qgtend(1,:)

     if (do_graupel .or. do_hail) then
!++ag Graupel Sedimentation
    ! calculate number of split time steps to ensure courant stability criteria
     ! for sedimentation calculations
     !-------------------------------------------------------------------
       nstep = 1 + nint(max( maxval( fg(i,nlb:nlev)*pdel_inv(i,nlb:nlev)), &
                             maxval(fng(i,nlb:nlev)*pdel_inv(i,nlb:nlev))) * deltat)
       nstep = min(nstep, nstep_def)

       tx2 = one / nstep
       tx1 = tx2 * deltat
       tx3 = tx2 / g
     ! loop over sedimentation sub-time step to ensure stability
     !==============================================================
       do n = 1,nstep

        ! top of model
          k = 1

        ! add fallout terms to microphysical tendencies

          tx5         = dumg(i,k)
          tx7         = pdel_inv(i,k) * tx1
          dumg(i,k)   = tx5 / (one + fg(i,k)*tx7)
          tx6         = (dumg(i,k)-tx5) * oneodt
          qgtend(i,k) = qgtend(i,k) + tx6
          tx5         = dumng(i,k)
          dumng(i,k)  = tx5 / (one + fng(i,k)*tx7)
          ngtend(i,k) = ngtend(i,k) + (dumng(i,k)-tx5) * oneodt

       ! sedimentation tendency for output
          qgsedten(i,k) = qgsedten(i,k) + tx6

          faloutg(k)    = fg(i,k)  * dumg(i,k)
          faloutng(k)   = fng(i,k) * dumng(i,k)

          gflx(i,k+1)   = gflx(i,k+1) + faloutg(k) * tx3   ! Ice flux

          do k = 2,nlev

             tx5         = dumg(i,k)
             tx7         = pdel_inv(i,k) * tx1
             dumg(i,k)   = (tx5 + faloutg(k-1)*tx7) / (one + fg(i,k)*tx7)
             tx6         = (dumg(i,k)-tx5) * oneodt
                                                                   ! add fallout terms to eulerian tendencies
             qgtend(i,k) = qgtend(i,k) + tx6
             tx5         = dumng(i,k)
             dumng(i,k)  = (tx5 + faloutng(k-1)*tx7)  / (one + fng(i,k)*tx7)
             ngtend(i,k) = ngtend(i,k) + (dumng(i,k)-tx5) * oneodt


             qgsedten(i,k) = qgsedten(i,k) + tx6                      ! sedimentation tendency for output


             faloutg(k)    = fg(i,k)  * dumg(i,k)
             faloutng(k)   = fng(i,k) * dumng(i,k)

             gflx(i,k+1)   = gflx(i,k+1) + faloutg(k) * tx3           ! Ice flux
          end do

        ! units below are m/s
        ! sedimentation flux at surface is added to precip flux at surface
        ! to get total precip (cloud + precip water) rate

          prect(i) = prect(i) + faloutg(nlev) * (tx3*0.001_r8)
          preci(i) = preci(i) + faloutg(nlev) * (tx3*0.001_r8)

       end do   !! nstep loop
     endif
!  if (lprnt) write(0,*)' qgtnds=',qgtend(1,:)
!--ag
  enddo       ! end of i loop
  ! end sedimentation

! if (lprnt) write(0,*)' prectaftsed=',prect(i),' preci=',preci(i)

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  ! get new update for variables that includes sedimentation tendency
  ! note : here dum variables are grid-average, NOT in-cloud

  do k=1,nlev
     do i=1,mgncol
        dumc(i,k)  = max(qc(i,k)+qctend(i,k)*deltat, zero)
        dumi(i,k)  = max(qi(i,k)+qitend(i,k)*deltat, zero)
        dumnc(i,k) = max(nc(i,k)+nctend(i,k)*deltat, zero)
        dumni(i,k) = max(ni(i,k)+nitend(i,k)*deltat, zero)

        dumr(i,k)  = max(qr(i,k)+qrtend(i,k)*deltat, zero)
        dumnr(i,k) = max(nr(i,k)+nrtend(i,k)*deltat, zero)
        dums(i,k)  = max(qs(i,k)+qstend(i,k)*deltat, zero)
        dumns(i,k) = max(ns(i,k)+nstend(i,k)*deltat, zero)

!++ag
        dumg(i,k)  = max(qg(i,k)+qgtend(i,k)*deltat, zero)
!  Moorthi testing
        if (dumg(i,k) > 0.01_r8) then
          tx2 = dumg(i,k) - 0.01_r8
          dumg(i,k) = 0.01_r8
          dums(i,k) = dums(i,k) + tx2
          qstend(i,k) = (dums(i,k) - qs(i,k)) * oneodt
          qgtend(i,k) = (dumg(i,k) - qg(i,k)) * oneodt
        endif
!  Moorthi testing
        dumng(i,k) = max(ng(i,k)+ngtend(i,k)*deltat, zero)
!--ag

        ! switch for specification of droplet and crystal number
        if (nccons) then
           dumnc(i,k) = ncnst*rhoinv(i,k)*lcldm(i,k)
        end if

        ! switch for specification of cloud ice number
        if (nicons) then
           dumni(i,k) = ninst*rhoinv(i,k)*icldm(i,k)
        end if

!++ag
        ! switch for specification of graupel number
        if (ngcons) then
           dumng(i,k) = ngnst*rhoinv(i,k)*precip_frac(i,k)
        end if
!--ag

        if (dumc(i,k) < qsmall) dumnc(i,k) = zero
        if (dumi(i,k) < qsmall) dumni(i,k) = zero
        if (dumr(i,k) < qsmall) dumnr(i,k) = zero
        if (dums(i,k) < qsmall) dumns(i,k) = zero
!++ag
        if (dumg(i,k) < qsmall) dumng(i,k) = zero
!--ag

     enddo

  enddo

  ! calculate instantaneous processes (melting, homogeneous freezing)
  !====================================================================

  ! melting of snow at +2 C
  do k=1,nlev

     do i=1,mgncol

        tx1 = t(i,k) + tlat(i,k)*(deltat/cpp) - snowmelt
        if (tx1 > zero) then
           if (dums(i,k) > zero) then

              ! make sure melting snow doesn't reduce temperature below threshold
              dum = -(xlf/cpp) * dums(i,k)
              if (tx1+dum <  zero) then
                 dum = min(one, max(zero, -tx1/dum))
              else
                 dum = one
              end if

              tx2 = dum * oneodt
              qstend(i,k) = qstend(i,k) - tx2*dums(i,k)
              nstend(i,k) = nstend(i,k) - tx2*dumns(i,k)
              qrtend(i,k) = qrtend(i,k) + tx2*dums(i,k)
              nrtend(i,k) = nrtend(i,k) + tx2*dumns(i,k)

              dum1 = - xlf * tx2 * dums(i,k)
              tlat(i,k)       = dum1 + tlat(i,k)
              meltsdttot(i,k) = dum1 + meltsdttot(i,k)
           end if
        end if
     enddo
  enddo

  if (do_graupel .or. do_hail) then
!++ag

  ! melting of graupel at +2 C

    do k=1,nlev

      do i=1,mgncol

        tx1 = t(i,k) + tlat(i,k)*(deltat/cpp) - snowmelt
        if (tx1 > zero) then
          if (dumg(i,k) > zero) then

           ! make sure melting graupel doesn't reduce temperature below threshold
             dum = -(xlf/cpp) * dumg(i,k)
             if (tx1+dum <  zero) then
                dum = max(zero, min(one, -tx1/dum))
             else
                dum = one
             end if

             tx2 = dum * oneodt

             qgtend(i,k) = qgtend(i,k) - tx2*dumg(i,k)
             ngtend(i,k) = ngtend(i,k) - tx2*dumng(i,k)
             qrtend(i,k) = qrtend(i,k) + tx2*dumg(i,k)
             nrtend(i,k) = nrtend(i,k) + tx2*dumng(i,k)

             dum1 = - xlf*tx2*dumg(i,k)
             tlat(i,k)       = dum1 + tlat(i,k)
             meltsdttot(i,k) = dum1 + meltsdttot(i,k)
          end if
        end if
      enddo
    enddo

!--ag
  endif

  do k=1,nlev
     do i=1,mgncol

       ! freezing of rain at -5 C

       tx1 = t(i,k) + tlat(i,k) * (deltat/cpp) - rainfrze
       if (tx1 < zero) then

          if (dumr(i,k) > zero) then

             ! make sure freezing rain doesn't increase temperature above threshold
             dum = (xlf/cpp) * dumr(i,k)
             if (tx1+dum > zero) then
                dum = min(one, max(zero, -tx1/dum))
             else
                dum = one
             end if
             tx2 = dum * oneodt
             qrtend(i,k) = qrtend(i,k) - tx2 * dumr(i,k)
             nrtend(i,k) = nrtend(i,k) - tx2 * dumnr(i,k)

             ! get mean size of rain = 1/lamr, add frozen rain to either snow or cloud ice
             ! depending on mean rain size

             call size_dist_param_basic(mg_rain_props, dumr(i,k), dumnr(i,k), lamr(i,k))

             if (lamr(i,k) < one/Dcs) then
!++ag freeze rain to graupel
               if (do_hail .or. do_graupel) then
                  qgtend(i,k) = qgtend(i,k) + tx2 * dumr(i,k)
                  ngtend(i,k) = ngtend(i,k) + tx2 * dumnr(i,k)
               else
                  qstend(i,k) = qstend(i,k) + tx2 * dumr(i,k)
                  nstend(i,k) = nstend(i,k) + tx2 * dumnr(i,k)
               end if
!--ag
             else
               qitend(i,k) = qitend(i,k) + tx2 * dumr(i,k)
               nitend(i,k) = nitend(i,k) + tx2 * dumnr(i,k)
             end if
                                                         ! heating tendency
             dum1 = xlf*dum*dumr(i,k)*oneodt
             frzrdttot(i,k) = dum1 + frzrdttot(i,k)
             tlat(i,k)      = dum1 + tlat(i,k)

          end if
       end if

     enddo
  enddo
  if (do_cldice) then
     do k=1,nlev
        do i=1,mgncol
           tx1 = t(i,k) + tlat(i,k) * (deltat/cpp) - tmelt
           if (tx1 > zero) then
              if (dumi(i,k) > zero) then

                 ! limit so that melting does not push temperature below freezing
                 !-----------------------------------------------------------------
                 dum = -dumi(i,k)*xlf/cpp
                 if (tx1+dum < zero) then
                    dum = min(one, max(zero, -tx1/dum))
                 else
                    dum = one
                 end if

                 tx2 = dum * oneodt
                 qctend(i,k) = qctend(i,k) + tx2*dumi(i,k)

                 ! for output
                 melttot(i,k) = tx2*dumi(i,k)

                 ! assume melting ice produces droplet
                 ! mean volume radius of 8 micron

                 nctend(i,k) = nctend(i,k) + three*tx2*dumi(i,k)/(four*pi*5.12e-16_r8*rhow)

                 qitend(i,k) = ((one-dum)*dumi(i,k)-qi(i,k))  * oneodt
                 nitend(i,k) = ((one-dum)*dumni(i,k)-ni(i,k)) * oneodt
                 tlat(i,k)   = tlat(i,k) - xlf*tx2*dumi(i,k)
              end if
           end if
        enddo
     enddo

! if (lprnt) write(0,*)' tlat6=',tlat(1,:)*deltat
! if (lprnt) write(0,*)' qitend=',qitend(1,nlev-45:nlev)*deltat
! if (lprnt) write(0,*)' qctend=',qctend(1,nlev-45:nlev)*deltat

     ! homogeneously freeze droplets at -40 C
     !-----------------------------------------------------------------

     do k=1,nlev
        do i=1,mgncol
           tx1 = t(i,k) + tlat(i,k)*(deltat/cpp) - 233.15_r8
           if (tx1 < zero) then
              if (dumc(i,k) > zero) then

                 ! limit so that freezing does not push temperature above threshold
                 dum = (xlf/cpp) * dumc(i,k)
                 if (tx1+dum > zero) then
                    dum = min(one, max(zero, -tx1/dum))
                 else
                    dum = one
                 end if

                 tx2 = dum * oneodt * dumc(i,k)
                 qitend(i,k)  = tx2 + qitend(i,k)
                 homotot(i,k) = tx2                   ! for output

                 ! assume 25 micron mean volume radius of homogeneously frozen droplets
                 ! consistent with size of detrained ice in stratiform.F90

                 nitend(i,k) = nitend(i,k) + tx2*(three/(four*pi*1.563e-14_r8* 500._r8))
                 qctend(i,k) = ((one-dum)*dumc(i,k)-qc(i,k))  * oneodt
                 nctend(i,k) = ((one-dum)*dumnc(i,k)-nc(i,k)) * oneodt
                 tlat(i,k)   = tlat(i,k) + xlf*tx2
              end if
           end if
        enddo 
     enddo 
     ! remove any excess over-saturation, which is possible due to non-linearity when adding
     ! together all microphysical processes
     !-----------------------------------------------------------------
     ! follow code similar to old CAM scheme
     do k=1,nlev
        do i=1,mgncol

           qtmp = q(i,k) + qvlat(i,k) * deltat
           ttmp = t(i,k) + tlat(i,k)  * (deltat/cpp)

           ! use rhw to allow ice supersaturation
           !call qsat_water(ttmp, p(i,k), esn, qvn)
           esn = min(fpvsl(ttmp), p(i,k))
           qvn = epsqs*esn/(p(i,k)-omeps*esn) * qsfm(i,k)
!          qvn = epsqs*esn/(p(i,k)-omeps*esn)


           if (qtmp > qvn .and. qvn > zero .and. allow_sed_supersat) then
              ! expression below is approximate since there may be ice deposition
              dum = (qtmp-qvn)/(one+xxlv_squared*qvn/(cpp*rv*ttmp*ttmp)) * oneodt
              ! add to output cme
              cmeout(i,k) = cmeout(i,k) + dum
              ! now add to tendencies, partition between liquid and ice based on temperature
              if (ttmp > 268.15_r8) then
                 dum1 = zero
                 ! now add to tendencies, partition between liquid and ice based on te
                 !-------------------------------------------------------
              else if (ttmp < 238.15_r8) then
                 dum1 = one
              else
                 dum1 = (268.15_r8-ttmp)/30._r8
              end if

              tx1 = xxls*dum1 + xxlv*(one-dum1)
              dum = (qtmp-qvn)/(one+tx1*tx1*qvn/(cpp*rv*ttmp*ttmp)) * oneodt
              tx2 = dum*(one-dum1)
              qctend(i,k)   = qctend(i,k) + tx2
              qcrestot(i,k) = tx2                         ! for output
              qitend(i,k)   = qitend(i,k) + dum*dum1
              qirestot(i,k) = dum*dum1
              qvlat(i,k)    = qvlat(i,k) - dum
              ! for output
              qvres(i,k)    = -dum
              tlat(i,k)     = tlat(i,k) + dum*tx1
           end if
        enddo 
     enddo 
  end if

! if (lprnt) write(0,*)' tlat7=',tlat(1,:)*deltat

  ! calculate effective radius for pass to radiation code
  !=========================================================
  ! if no cloud water, default value is 10 micron for droplets,
  ! 25 micron for cloud ice

  ! update cloud variables after instantaneous processes to get effective radius
  ! variables are in-cloud to calculate size dist parameters
  do k=1,nlev
     do i=1,mgncol
        if (lcldm(i,k) > mincld) then
          tx1 = one / lcldm(i,k)
        else
          tx1 = zero
        endif
        if (icldm(i,k) > mincld) then
          tx2 = one / icldm(i,k)
        else
          tx2 = zero
        endif
        if (precip_frac(i,k) > mincld) then
          tx3 = one / precip_frac(i,k)
        else
          tx3 = zero
        endif
        dumc(i,k)  = max(qc(i,k)+qctend(i,k)*deltat, zero) * tx1
        dumi(i,k)  = max(qi(i,k)+qitend(i,k)*deltat, zero) * tx2
        dumnc(i,k) = max(nc(i,k)+nctend(i,k)*deltat, zero) * tx1
        dumni(i,k) = max(ni(i,k)+nitend(i,k)*deltat, zero) * tx2

        dumr(i,k)  = max(qr(i,k)+qrtend(i,k)*deltat, zero) * tx3
        dumnr(i,k) = max(nr(i,k)+nrtend(i,k)*deltat, zero) * tx3
        dums(i,k)  = max(qs(i,k)+qstend(i,k)*deltat, zero) * tx3
        dumns(i,k) = max(ns(i,k)+nstend(i,k)*deltat, zero) * tx3

!++ag
        dumg(i,k)  = max(qg(i,k)+qgtend(i,k)*deltat, zero) * tx3
        dumng(i,k) = max(ng(i,k)+ngtend(i,k)*deltat, zero) * tx3
!--ag

        ! switch for specification of droplet and crystal number
        if (nccons) then
           dumnc(i,k) = ncnst * rhoinv(i,k)
        end if

        ! switch for specification of cloud ice number
        if (nicons) then
           dumni(i,k) = ninst * rhoinv(i,k)
        end if

!++ag
        ! switch for specification of graupel number
        if (ngcons) then
           dumng(i,k) = ngnst*rhoinv(i,k)*precip_frac(i,k)
        end if
!--ag

        ! limit in-cloud mixing ratio to reasonable value of 5 g kg-1
!       dumc(i,k) = min(dumc(i,k), 5.e-3_r8)
!       dumi(i,k) = min(dumi(i,k), 5.e-3_r8)
        dumc(i,k) = min(dumc(i,k), 10.e-3_r8)
        dumi(i,k) = min(dumi(i,k), 10.e-3_r8)
        ! limit in-precip mixing ratios
        dumr(i,k) = min(dumr(i,k), 10.e-3_r8)
        dums(i,k) = min(dums(i,k), 10.e-3_r8)
!++ag
        dumg(i,k) = min(dumg(i,k), 10.e-3_r8)
!--ag
     enddo
  enddo
  ! cloud ice effective radius
  !-----------------------------------------------------------------

  if (do_cldice) then
     do k=1,nlev
        do i=1,mgncol
           if (dumi(i,k) >= qsmall) then

              tx1 = dumni(i,k)
              call size_dist_param_basic(mg_ice_props, dumi(i,k), dumni(i,k), &
                                         lami(i,k), dumni0)

              if (dumni(i,k) /= tx1) then
                 ! adjust number conc if needed to keep mean size in reasonable range
                 nitend(i,k) = (dumni(i,k)*icldm(i,k)-ni(i,k)) * oneodt
              end if

              tx1 = one / lami(i,k)
!             effi(i,k)   = (1.5_r8*1.e6_r8) * tx1
              effi(i,k)   = (three*1.e6_r8) * tx1
              sadice(i,k) = two*pi*(tx1*tx1*tx1)*dumni0*rho(i,k)*1.e-2_r8  ! m2/m3 -> cm2/cm3

           else
              effi(i,k)   = 50._r8
              sadice(i,k) = zero
           end if

           ! ice effective diameter for david mitchell's optics
           deffi(i,k) = effi(i,k) * (rhoi+rhoi)/rhows
        enddo
     enddo
  !else
     !do k=1,nlev
        !do i=1,mgncol
           ! NOTE: If CARMA is doing the ice microphysics, then the ice effective
           ! radius has already been determined from the size distribution.
           !effi(i,k) = re_ice(i,k) * 1.e6_r8      ! m -> um
           !deffi(i,k)=effi(i,k) * 2._r8
           !sadice(i,k) = 4._r8*pi*(effi(i,k)**2)*ni(i,k)*rho(i,k)*1e-2_r8
        !enddo
     !enddo
  end if

  ! cloud droplet effective radius
  !-----------------------------------------------------------------
  do k=1,nlev
     do i=1,mgncol
        if (dumc(i,k) >= qsmall) then


           ! switch for specification of droplet and crystal number
           if (nccons) then
              ! make sure nc is consistence with the constant N by adjusting tendency, need
              ! to multiply by cloud fraction
              ! note that nctend may be further adjusted below if mean droplet size is
              ! out of bounds

              nctend(i,k) = (ncnst*rhoinv(i,k)*lcldm(i,k)-nc(i,k)) * oneodt

           end if

           dum = dumnc(i,k)

           call size_dist_param_liq(mg_liq_props, dumc(i,k), dumnc(i,k), rho(i,k), &
                                    pgam(i,k), lamc(i,k))

           if (dum /= dumnc(i,k)) then
              ! adjust number conc if needed to keep mean size in reasonable range
              nctend(i,k) = (dumnc(i,k)*lcldm(i,k)-nc(i,k)) * oneodt
           end if

           effc(i,k) = (half*1.e6_r8) * (pgam(i,k)+three) / lamc(i,k)
           !assign output fields for shape here
           lamcrad(i,k) = lamc(i,k)
           pgamrad(i,k) = pgam(i,k)


           ! recalculate effective radius for constant number, in order to separate
           ! first and second indirect effects
           !======================================
           ! assume constant number of 10^8 kg-1

           dumnc(i,k) = 1.e8_r8

           ! Pass in "false" adjust flag to prevent number from being changed within
           ! size distribution subroutine.
           call size_dist_param_liq(mg_liq_props, dumc(i,k), dumnc(i,k), rho(i,k), &
                                    pgam(i,k), lamc(i,k))

           effc_fn(i,k) = (half*1.e6_r8) * (pgam(i,k)+three)/lamc(i,k)

        else
           effc(i,k)    = ten
           lamcrad(i,k) = zero
           pgamrad(i,k) = zero
           effc_fn(i,k) = ten
        end if
     enddo
  enddo
  ! recalculate 'final' rain size distribution parameters
  ! to ensure that rain size is in bounds, adjust rain number if needed
  do k=1,nlev
     do i=1,mgncol

        if (dumr(i,k) >= qsmall) then

           dum = dumnr(i,k)

           call size_dist_param_basic(mg_rain_props, dumr(i,k), dumnr(i,k), lamr(i,k))

           if (dum /= dumnr(i,k)) then
              ! adjust number conc if needed to keep mean size in reasonable range
              nrtend(i,k) = (dumnr(i,k)*precip_frac(i,k)-nr(i,k)) *oneodt
           end if

        end if
     enddo
  enddo
  ! recalculate 'final' snow size distribution parameters
  ! to ensure that snow size is in bounds, adjust snow number if needed
  do k=1,nlev
     do i=1,mgncol
        if (dums(i,k) >= qsmall) then

           dum = dumns(i,k)

           call size_dist_param_basic(mg_snow_props, dums(i,k), dumns(i,k), &
                                      lams(i,k), n0=dumns0)

           if (dum /= dumns(i,k)) then
              ! adjust number conc if needed to keep mean size in reasonable range
              nstend(i,k) = (dumns(i,k)*precip_frac(i,k)-ns(i,k)) * oneodt
           end if

           tx1 = (two*pi*1.e-2_r8) / (lams(i,k)*lams(i,k)*lams(i,k))
           sadsnow(i,k) = tx1*dumns0*rho(i,k)             ! m2/m3 -> cm2/cm3

        end if


     end do ! vertical k loop
  enddo
  do k=1,nlev
     do i=1,mgncol
        ! if updated q (after microphysics) is zero, then ensure updated n is also zero
        !=================================================================================
        if (qc(i,k)+qctend(i,k)*deltat < qsmall)                 nctend(i,k) = -nc(i,k) * oneodt
        if (do_cldice .and. qi(i,k)+qitend(i,k)*deltat < qsmall) nitend(i,k) = -ni(i,k) * oneodt
        if (qr(i,k)+qrtend(i,k)*deltat < qsmall)                 nrtend(i,k) = -nr(i,k) * oneodt
        if (qs(i,k)+qstend(i,k)*deltat < qsmall)                 nstend(i,k) = -ns(i,k) * oneodt
!++ag
        if (qg(i,k)+qgtend(i,k)*deltat < qsmall)                 ngtend(i,k) = -ng(i,k) * oneodt
!--ag

     end do

  end do

  ! DO STUFF FOR OUTPUT:
  !==================================================

  do k=1,nlev
    do i=1,mgncol

  ! qc and qi are only used for output calculations past here,
  ! so add qctend and qitend back in one more time
      qc(i,k) = qc(i,k) + qctend(i,k)*deltat
      qi(i,k) = qi(i,k) + qitend(i,k)*deltat

  ! averaging for snow and rain number and diameter
  !--------------------------------------------------

  ! drout2/dsout2:
  ! diameter of rain and snow
  ! dsout:
  ! scaled diameter of snow (passed to radiation in CAM)
  ! reff_rain/reff_snow:
  ! calculate effective radius of rain and snow in microns for COSP using Eq. 9 of COSP v1.3 manual

      if (qrout(i,k) > 1.e-7_r8 .and. nrout(i,k) > zero) then
         qrout2(i,k)    = qrout(i,k) * precip_frac(i,k)
         nrout2(i,k)    = nrout(i,k) * precip_frac(i,k)
     ! The avg_diameter call does the actual calculation; other diameter
     ! outputs are just drout2 times constants.
         drout2(i,k)    = avg_diameter(qrout(i,k), nrout(i,k), rho(i,k), rhow)
         freqr(i,k)     = precip_frac(i,k)

         reff_rain(i,k) = (1.e6_r8*1.5_r8) * drout2(i,k)
      else
         qrout2(i,k)    = zero
         nrout2(i,k)    = zero
         drout2(i,k)    = zero
         freqr(i,k)     = zero
         reff_rain(i,k) = zero
      endif

      if (qsout(i,k) > 1.e-7_r8 .and. nsout(i,k) > zero) then
         qsout2(i,k)    = qsout(i,k) * precip_frac(i,k)
         nsout2(i,k)    = nsout(i,k) * precip_frac(i,k)
     ! The avg_diameter call does the actual calculation; other diameter
     ! outputs are just dsout2 times constants.
         dsout2(i,k)    = avg_diameter(qsout(i,k), nsout(i,k), rho(i,k), rhosn)
         freqs(i,k)     = precip_frac(i,k)

         dsout(i,k)     = three*rhosn/rhows*dsout2(i,k)

         reff_snow(i,k) = (1.e6_r8*three) * dsout2(i,k)
      else
         dsout(i,k)     = zero
         qsout2(i,k)    = zero
         nsout2(i,k)    = zero 
         dsout2(i,k)    = zero
         freqs(i,k)     = zero
         reff_snow(i,k) = zero
      endif

    enddo
  enddo

  ! analytic radar reflectivity
  !--------------------------------------------------
  ! formulas from Matthew Shupe, NOAA/CERES
  ! *****note: radar reflectivity is local (in-precip average)
  ! units of mm^6/m^3

  do k=1,nlev
     do i = 1,mgncol
!       if (qc(i,k) >= qsmall .and. (nc(i,k)+nctend(i,k)*deltat) > ten .and. lcldm(i,k) > mincld) then
        if (qc(i,k) >= qsmall .and. (nc(i,k)+nctend(i,k)*deltat) > ten) then
           tx1 = rho(i,k) / lcldm(i,k)
           tx2 = 1000._r8 * qc(i,k) * tx1
           dum = tx2 * tx2 * lcldm(i,k)                            &
                /(0.109_r8*(nc(i,k)+nctend(i,k)*deltat)*tx1*1.e-6_r8*precip_frac(i,k))
!          dum = (qc(i,k)/lcldm(i,k)*rho(i,k)*1000._r8)**2 &
!               /(0.109_r8*(nc(i,k)+nctend(i,k)*deltat)/lcldm(i,k)*rho(i,k)/1.e6_r8)*lcldm(i,k)/precip_frac(i,k)
        else
           dum = zero
        end if
!       if (qi(i,k) >= qsmall .and. icldm(i,k) > mincld) then
        if (qi(i,k) >= qsmall) then
!          dum1 = (qi(i,k)*rho(i,k)/icldm(i,k)*1000._r8/0.1_r8)**(one/0.63_r8)*icldm(i,k)/precip_frac(i,k)
           dum1 = (qi(i,k)*rho(i,k)/icldm(i,k)*10000._r8)**(one/0.63_r8)*icldm(i,k)/precip_frac(i,k)
        else
           dum1 = zero
        end if

        if (qsout(i,k) >= qsmall) then
!          dum1 = dum1 + (qsout(i,k)*rho(i,k)*1000._r8/0.1_r8)**(one/0.63_r8)
           dum1 = dum1 + (qsout(i,k)*rho(i,k)*10000._r8)**(one/0.63_r8)
        end if

        refl(i,k) = dum + dum1

        ! add rain rate, but for 37 GHz formulation instead of 94 GHz
        ! formula approximated from data of Matrasov (2007)
        ! rainrt is the rain rate in mm/hr
        ! reflectivity (dum) is in DBz

        if (rainrt(i,k) >= 0.001_r8) then
           dum = rainrt(i,k) * rainrt(i,k)
           dum = log10(dum*dum*dum) + 16._r8

           ! convert from DBz to mm^6/m^3

           dum = ten**(dum/ten)
        else
           ! don't include rain rate in R calculation for values less than 0.001 mm/hr
           dum = zero
        end if

        ! add to refl

        refl(i,k) = refl(i,k) + dum

        !output reflectivity in Z.
        areflz(i,k) = refl(i,k) * precip_frac(i,k)

        ! convert back to DBz

        if (refl(i,k) > minrefl) then
           refl(i,k) = ten*log10(refl(i,k))
        else
           refl(i,k) = -9999._r8
        end if

        !set averaging flag
        if (refl(i,k) > mindbz) then
           arefl(i,k) = refl(i,k) * precip_frac(i,k)
           frefl(i,k) = precip_frac(i,k)
        else
           arefl(i,k)  = zero
           areflz(i,k) = zero
           frefl(i,k)  = zero
        end if

        ! bound cloudsat reflectivity

        csrfl(i,k) = min(csmax,refl(i,k))

        !set averaging flag
        if (csrfl(i,k) > csmin) then
           acsrfl(i,k) = refl(i,k) * precip_frac(i,k)
           fcsrfl(i,k) = precip_frac(i,k)
        else
           acsrfl(i,k) = zero
           fcsrfl(i,k) = zero
        end if

     end do
  end do

  do k=1,nlev
     do i = 1,mgncol
  !redefine fice here....
       tx2 = qsout(i,k) + qi(i,k)
       tx1 = tx2 + qrout(i,k) + qc(i,k)
       if ( tx2 > qsmall .and. tx1 > qsmall) then
         nfice(i,k) = min(tx2/tx1, one)
       else
         nfice(i,k) = zero
       endif
    enddo
  enddo

end subroutine micro_mg_tend
!> @}

!========================================================================
!OUTPUT CALCULATIONS
!========================================================================

!>\ingroup mg3_mp
!! This subroutine calculates effective radii for rain and cloud.
subroutine calc_rercld(lamr, n0r, lamc, pgam, qric, qcic, ncic, rercld, mgncol,nlev)
  integer, intent(in) :: mgncol, nlev                           ! horizontal and vertical dimension
  real(r8), dimension(mgncol,nlev), intent(in) :: lamr          ! rain size parameter (slope)
  real(r8), dimension(mgncol,nlev), intent(in) :: n0r           ! rain size parameter (intercept)
  real(r8), dimension(mgncol,nlev), intent(in) :: lamc          ! size distribution parameter (slope)
  real(r8), dimension(mgncol,nlev), intent(in) :: pgam          ! droplet size parameter
  real(r8), dimension(mgncol,nlev), intent(in) :: qric          ! in-cloud rain mass mixing ratio
  real(r8), dimension(mgncol,nlev), intent(in) :: qcic          ! in-cloud cloud liquid
  real(r8), dimension(mgncol,nlev), intent(in) :: ncic          ! in-cloud droplet number concentration

  real(r8), dimension(mgncol,nlev), intent(inout) :: rercld     ! effective radius calculation for rain + cloud

  ! combined size of precip & cloud drops
  real(r8) :: Atmp

  integer :: i, k

  do k=1,nlev
    do i=1,mgncol
      ! Rain drops
      if (lamr(i,k) > zero) then
        Atmp = n0r(i,k) * (half*pi) / (lamr(i,k)*lamr(i,k)*lamr(i,k))
      else
        Atmp = zero
      end if

      ! Add cloud drops
      if (lamc(i,k) > zero) then
        Atmp = Atmp + ncic(i,k) * pi * rising_factorial(pgam(i,k)+one, 2) &
                    / (four*lamc(i,k)*lamc(i,k))
      end if

      if (Atmp > zero) then
        rercld(i,k) = rercld(i,k) + three *(qric(i,k) + qcic(i,k)) / (four * rhow * Atmp)
      end if
    enddo
  enddo
end subroutine calc_rercld

!========================================================================

end module micro_mg3_0
!>@}

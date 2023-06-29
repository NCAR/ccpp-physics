
!=======================================================================
! Current Code Owner: DWD, Ulrich Schaettler
!  phone:  +49  69  8062 2739
!  fax:    +49  69  8236 1493
!  email:  uschaettler@dwd.d400.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.1        1998/03/11 Ulrich Schaettler
!  Initial release
! 1.8        1998/08/03 Ulrich Schaettler
!  Eliminated intgribf, intgribc, irealgrib, iwlength and put it to data_io.
! 1.10       1998/09/29 Ulrich Schaettler
!  Eliminated parameters for grid point and diagnostic calculations.
! !VERSION!  !DATE!     <Your name>
!  <Modification comments>
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!
! reorganize the FLake to module_FLake.F90 by Shaobo Zhang in 2016-7-13
! added a new layer for deep lakes by Shaobo Zhang in 2016-11-15
!
!=======================================================================

!------------------------------------------------------------------------------

MODULE data_parameters

!------------------------------------------------------------------------------
!
! Description:
!  Global parameters for the program are defined.
!  Actually, scratch that. We'll import them from machine.F instead.
!
  use machine, only: ireals=>kind_phys, iintegers=>kind_INTEGER

!=======================================================================

END MODULE data_parameters

!------------------------------------------------------------------------------

MODULE flake_albedo_ref

!------------------------------------------------------------------------------
!
! Description:
!
!  This module contains "reference" values of albedo 
!  for the lake water, lake ice and snow. 
!  As in "flake_paramoptic_ref", two ice categories, viz. white ice and blue ice,
!  and two snow categories, viz. dry snow and melting snow, are used.  
!

USE data_parameters, ONLY :      &
  ireals                       , & ! KIND-type parameter for real variables 
  iintegers                        ! KIND-type parameter for "normal" integer variables

use machine,               only: kind_phys
!==============================================================================

IMPLICIT NONE

!==============================================================================
!
! Declarations
 
!  Albedo for water, ice and snow.
!REAL (KIND = ireals), PARAMETER ::        &
!  albedo_water_ref       = 0.070  , & ! Water
!  albedo_whiteice_ref    = 0.600  , & ! White ice
!  albedo_blueice_ref     = 0.100  , & ! Blue ice
!  albedo_drysnow_ref     = 0.600  , & ! Dry snow 
!  albedo_meltingsnow_ref = 0.100      ! Melting snow 

!  Empirical parameters.
!REAL (KIND = ireals), PARAMETER :: &
!  c_albice_MR = 95.60          ! Constant in the interpolation formula for 
                                     ! the ice albedo (Mironov and Ritter 2004)
!  Albedo for water, ice and snow.
REAL (KIND = kind_phys), PARAMETER ::        &
  albedo_water_ref       = 0.07  , & ! Water
  albedo_whiteice_ref    = 0.60  , & ! White ice
  albedo_blueice_ref     = 0.10  , & ! Blue ice
!  albedo_drysnow_ref     = 0.60  , & ! Dry snow
  albedo_drysnow_ref     = 0.90  , & ! Dry snow
  albedo_meltingsnow_ref = 0.10      ! Melting snow

!  Empirical parameters.
REAL (KIND = kind_phys), PARAMETER :: &
  c_albice_MR = 95.6                 ! Constant in the interpolation formula for
                                     ! the ice albedo (Mironov and Ritter 2004)


!==============================================================================

END MODULE flake_albedo_ref

!------------------------------------------------------------------------------

MODULE flake_configure

!------------------------------------------------------------------------------
!
! Description:
!
!  Switches and reference values of parameters 
!  that configure the lake model FLake are set.
!

USE data_parameters , ONLY : &
  ireals                   , & ! KIND-type parameter for real variables 
  iintegers                    ! KIND-type parameter for "normal" integer variables

use machine,               only: kind_phys
!==============================================================================

IMPLICIT NONE

!==============================================================================
!
! Declarations
! changed by Shaobo Zhang
LOGICAL lflk_botsed_use
!LOGICAL, PARAMETER :: &
!  lflk_botsed_use   = .TRUE.        ! .TRUE. indicates that the bottom-sediment scheme is used
                                     ! to compute the depth penetrated by the thermal wave, 
                                     ! the temperature at this depth and the bottom heat flux.
                                     ! Otherwise, the heat flux at the water-bottom sediment interface
                                     ! is set to zero, the depth penetrated by the thermal wave 
                                     ! is set to a reference value defined below,
                                     ! and the temperature at this depth is set to 
                                     ! the temperature of maximum density of the fresh water.

!REAL (KIND = ireals), PARAMETER :: &
!  rflk_depth_bs_ref = 10.00    ! Reference value of the depth of the thermally active
                                     ! layer of bottom sediments [m].
                                     ! This value is used to (formally) define
                                     ! the depth penetrated by the thermal wave
                                     ! in case the bottom-sediment scheme is not used.

REAL (KIND = kind_phys), PARAMETER :: &
  rflk_depth_bs_ref = 10.0

!==============================================================================

END MODULE flake_configure

!------------------------------------------------------------------------------

MODULE flake_derivedtypes  

!------------------------------------------------------------------------------
!
! Description:
!
!  Derived type(s) is(are) defined.
!

USE data_parameters , ONLY : &
  ireals                   , & ! KIND-type parameter for real variables 
  iintegers                    ! KIND-type parameter for "normal" integer variables

!==============================================================================

IMPLICIT NONE

!==============================================================================
!
! Declarations

!  Maximum value of the wave-length bands 
!  in the exponential decay law for the radiation flux.
!  A storage for a ten-band approximation is allocated,
!  although a smaller number of bands is actually used.
INTEGER (KIND = iintegers), PARAMETER :: & 
  nband_optic_max = 10_iintegers

!  Define TYPE "opticpar_medium"
TYPE opticpar_medium
  INTEGER (KIND = iintegers)                        ::   & 
    nband_optic                                            ! Number of wave-length bands
  REAL (KIND = ireals), DIMENSION (nband_optic_max) ::   & 
    frac_optic                                         , & ! Fractions of total radiation flux 
    extincoef_optic                                        ! Extinction coefficients 
END TYPE opticpar_medium

!==============================================================================

END MODULE flake_derivedtypes  

!------------------------------------------------------------------------------

MODULE flake_paramoptic_ref

!------------------------------------------------------------------------------
!
! Description:
!
!  This module contains "reference" values of the optical characteristics
!  of the lake water, lake ice and snow. These reference values may be used 
!  if no information about the optical characteristics of the lake in question 
!  is available. An exponential decay law for the solar radiation flux is assumed.
!  In the simplest one-band approximation,
!  the extinction coefficient for water is set to a large value,
!  leading to the absorption of 95% of the incoming radiation 
!  within the uppermost 1 m of the lake water. 
!  The extinction coefficients for ice and snow are taken from 
!  Launiainen and Cheng (1998). The estimates for the ice correspond 
!  to the uppermost 0.1 m of the ice layer and to the clear sky conditions 
!  (see Table 2 in op. cit.).
!  Very large values of the extinction coefficients for ice and snow ("opaque")
!  can be used to prevent penetration of the solar radiation 
!  through the snow-ice cover.
!

USE data_parameters, ONLY :      &
  ireals                       , & ! KIND-type parameter for real variables 
  iintegers                        ! KIND-type parameter for "normal" integer variables

USE flake_derivedtypes, ONLY :   &
  nband_optic_max              , & ! Maximum value of the wave-length bands
  opticpar_medium                  ! Derived TYPE 

!==============================================================================

IMPLICIT NONE

!==============================================================================
!
! Declarations
 
INTEGER (KIND = iintegers), PRIVATE :: & ! Help variable(s)
  i                                      ! DO loop index

!  Optical characteristics for water, ice and snow.
!  The simplest one-band approximation is used as a reference.
!TYPE (opticpar_medium), PARAMETER ::                           & 
!  opticpar_water_ref = opticpar_medium(1,                      & ! Water (reference)
!    (/1.0, (0.0,i=2,nband_optic_max)/),            &
!    (/3.0, (1.E+100,i=2,nband_optic_max)/))      , &
!  opticpar_water_trans = opticpar_medium(2,                             & ! Transparent Water (two-band)
!    (/0.100, 0.900, (0.0,i=3,nband_optic_max)/),      &
!    (/2.00, 0.200, (1.E+100,i=3,nband_optic_max)/)) , &
!!_nu  opticpar_water_trans = opticpar_medium(1,                    & ! Transparent Water (one-band)
!!_nu    (/1.0, (0.0,i=2,nband_optic_max)/),            &
!!_nu    (/0.300, (1.E+100,i=2,nband_optic_max)/))    , &
!  opticpar_whiteice_ref = opticpar_medium(1,                   & ! White ice
!    (/1.0, (0.0,i=2,nband_optic_max)/),            &   
!    (/17.10, (1.E+100,i=2,nband_optic_max)/))    , &
!  opticpar_blueice_ref = opticpar_medium(1,                    & ! Blue ice
!    (/1.0, (0.0,i=2,nband_optic_max)/),            &
!    (/8.40, (1.E+100,i=2,nband_optic_max)/))     , &
!  opticpar_drysnow_ref = opticpar_medium(1,                    & ! Dry snow 
!    (/1.0, (0.0,i=2,nband_optic_max)/),            &
!    (/25.00, (1.E+100,i=2,nband_optic_max)/))    , &
!  opticpar_meltingsnow_ref = opticpar_medium(1,                & ! Melting snow 
!    (/1.0, (0.0,i=2,nband_optic_max)/),            &
!    (/15.00, (1.E+100,i=2,nband_optic_max)/))    , &
!  opticpar_ice_opaque = opticpar_medium(1,                     & ! Opaque ice
!    (/1.0, (0.0,i=2,nband_optic_max)/),            &
!    (/1.0E+070, (1.E+100,i=2,nband_optic_max)/)) , &
!  opticpar_snow_opaque = opticpar_medium(1,                    & ! Opaque snow
!    (/1.0, (0.0,i=2,nband_optic_max)/),            &
!    (/1.0E+070, (1.E+100,i=2,nband_optic_max)/)) 

TYPE (opticpar_medium), PARAMETER ::                           &
  opticpar_water_ref = opticpar_medium(1,                      & ! Water (reference)
    (/1., (0.,i=2,nband_optic_max)/),            &
    (/3., (1.E+10,i=2,nband_optic_max)/))      , &
  opticpar_water_trans = opticpar_medium(2,                             & ! Transparent Water (two-band)
    (/0.10, 0.90, (0.,i=3,nband_optic_max)/),      &
    (/2.0, 0.20, (1.E+10,i=3,nband_optic_max)/)) , &
  opticpar_whiteice_ref = opticpar_medium(1,                   & ! White ice
    (/1., (0.,i=2,nband_optic_max)/),            &
    (/17.1, (1.E+10,i=2,nband_optic_max)/))    , &
  opticpar_blueice_ref = opticpar_medium(1,                    & ! Blue ice
    (/1., (0.,i=2,nband_optic_max)/),            &
    (/8.4, (1.E+10,i=2,nband_optic_max)/))     , &
  opticpar_drysnow_ref = opticpar_medium(1,                    & ! Dry snow
    (/1., (0.,i=2,nband_optic_max)/),            &
    (/25.0, (1.E+10,i=2,nband_optic_max)/))    , &
  opticpar_meltingsnow_ref = opticpar_medium(1,                & ! Melting snow
    (/1., (0.,i=2,nband_optic_max)/),            &
    (/15.0, (1.E+10,i=2,nband_optic_max)/))    , &
  opticpar_ice_opaque = opticpar_medium(1,                     & ! Opaque ice
    (/1., (0.,i=2,nband_optic_max)/),            &
    (/1.0E+07, (1.E+10,i=2,nband_optic_max)/)) , &
  opticpar_snow_opaque = opticpar_medium(1,                    & ! Opaque snow
    (/1., (0.,i=2,nband_optic_max)/),            &
    (/1.0E+07, (1.E+10,i=2,nband_optic_max)/))


!==============================================================================

END MODULE flake_paramoptic_ref

!------------------------------------------------------------------------------

MODULE flake_parameters

!------------------------------------------------------------------------------
!
! Description:
!
!  Values of empirical constants of the lake model FLake 
!  and of several thermodynamic parameters are set.
!

USE data_parameters , ONLY : &
  ireals                   , & ! KIND-type parameter for real variables 
  iintegers                    ! KIND-type parameter for "normal" integer variables

use machine,               only: kind_phys

!==============================================================================

IMPLICIT NONE

!==============================================================================
!
! Declarations

!  Dimensionless constants 
!  in the equations for the mixed-layer depth 
!  and for the shape factor with respect to the temperature profile in the thermocline
REAL (KIND = kind_phys), PARAMETER ::         &
!  c_cbl_1       = 0.170            , & ! Constant in the CBL entrainment equation
!  c_cbl_2       = 1.0              , & ! Constant in the CBL entrainment equation
!  c_sbl_ZM_n    = 0.50             , & ! Constant in the ZM1996 equation for the equilibrium SBL depth
!  c_sbl_ZM_s    = 10.0             , & ! Constant in the ZM1996 equation for the equilibrium SBL depth
!  c_sbl_ZM_i    = 20.0             , & ! Constant in the ZM1996 equation for the equilibrium SBL depth
!  c_relax_h     = 0.0300           , & ! Constant in the relaxation equation for the SBL depth
!  c_relax_C     = 0.00300              ! Constant in the relaxation equation for the shape factor
                                             ! with respect to the temperature profile in the thermocline
  c_cbl_1       = 0.17            , & ! Constant in the CBL entrainment equation
  c_cbl_2       = 1.              , & ! Constant in the CBL entrainment equation
  c_sbl_ZM_n    = 0.5             , & ! Constant in the ZM1996 equation for the equilibrium SBL depth
  c_sbl_ZM_s    = 10.             , & ! Constant in the ZM1996 equation for the equilibrium SBL depth
  c_sbl_ZM_i    = 20.             , & ! Constant in the ZM1996 equation for the equilibrium SBL depth
  c_relax_h     = 0.030           , & ! Constant in the relaxation equation for the SBL depth
  c_relax_C     = 0.0030

!  Parameters of the shape functions 
!  Indices refer to T - thermocline, S - snow, I - ice,
!  B1 - upper layer of the bottom sediments, B2 - lower layer of the bottom sediments.
!  "pr0" and "pr1" denote zeta derivatives of the corresponding shape function 
!  at "zeta=0" ad "zeta=1", respectively.
REAL (KIND = kind_phys), PARAMETER ::         &
  C_T_min       = 0.5             , & ! Minimum value of the shape factor C_T (thermocline)
  C_T_max       = 0.8             , & ! Maximum value of the shape factor C_T (thermocline)
  Phi_T_pr0_1   = 40.0/3.0   , & ! Constant in the expression for the T shape-function derivative 
  Phi_T_pr0_2   = 20.0/3.0   , & ! Constant in the expression for the T shape-function derivative 
  C_TT_1        = 11.0/18.0  , & ! Constant in the expression for C_TT (thermocline)
  C_TT_2        = 7.0/45.0   , & ! Constant in the expression for C_TT (thermocline)
  C_B1          = 2.0/3.0    , & ! Shape factor (upper layer of bottom sediments)
  C_B2          = 3.0/5.0    , & ! Shape factor (lower layer of bottom sediments)
  Phi_B1_pr0    = 2.0              , & ! B1 shape-function derivative 
  C_S_lin       = 0.5             , & ! Shape factor (linear temperature profile in the snow layer)
  Phi_S_pr0_lin = 1.0              , & ! S shape-function derivative (linear profile) 
  C_I_lin       = 0.5             , & ! Shape factor (linear temperature profile in the ice layer)
  Phi_I_pr0_lin = 1.0              , & ! I shape-function derivative (linear profile) 
  Phi_I_pr1_lin = 1.0              , & ! I shape-function derivative (linear profile) 
  Phi_I_ast_MR  = 2.0              , & ! Constant in the MR2004 expression for I shape factor
  C_I_MR        = 1.0/12.0   , & ! Constant in the MR2004 expression for I shape factor
  H_Ice_max     = 3.0                  ! Maximum ice tickness in 
                                             ! the Mironov and Ritter (2004, MR2004) ice model [m] 

!  Security constants
REAL (KIND = kind_phys), PARAMETER ::         &
  h_Snow_min_flk = 1.0E-5         , & ! Minimum snow thickness [m]
  h_Ice_min_flk  = 1.0E-9         , & ! Minimum ice thickness [m]
  h_ML_min_flk   = 1.0E-2         , & ! Minimum mixed-layer depth [m]
  h_ML_max_flk   = 1.0E+3         , & ! Maximum mixed-layer depth [m]
  H_B1_min_flk   = 1.0E-3         , & ! Minimum thickness of the upper layer of bottom sediments [m]
  u_star_min_flk = 1.0E-6             ! Minimum value of the surface friction velocity [m s^{-1}]

!  Security constant(s)
REAL (KIND = kind_phys), PARAMETER ::         &
  c_small_flk    = 1.0E-10            ! A small number

!  Thermodynamic parameters
REAL (KIND = kind_phys), PARAMETER ::        &
  tpl_grav          = 9.81       , & ! Acceleration due to gravity [m s^{-2}]
  tpl_T_r           = 277.13     , & ! Temperature of maximum density of fresh water [K]
  tpl_T_f           = 273.15     , & ! Fresh water freezing point [K]
  tpl_a_T           = 1.6509E-05 , & ! Constant in the fresh-water equation of state [K^{-2}]
  tpl_rho_w_r       = 1.0E+03    , & ! Maximum density of fresh water [kg m^{-3}]
  tpl_rho_I         = 9.1E+02    , & ! Density of ice [kg m^{-3}]
  tpl_rho_S_min     = 1.0E+02    , & ! Minimum snow density [kg m^{-3}]
  tpl_rho_S_max     = 4.0E+02    , & ! Maximum snow density [kg m^{-3}]
  tpl_Gamma_rho_S   = 2.0E+02    , & ! Empirical parameter [kg m^{-4}]  
                                            ! in the expression for the snow density 
  tpl_L_f           = 3.3E+05    , & ! Latent heat of fusion [J kg^{-1}]
  tpl_c_w           = 4.2E+03    , & ! Specific heat of water [J kg^{-1} K^{-1}]
  tpl_c_I           = 2.1E+03    , & ! Specific heat of ice [J kg^{-1} K^{-1}]
  tpl_c_S           = 2.1E+03    , & ! Specific heat of snow [J kg^{-1} K^{-1}]
  tpl_kappa_w       = 5.46E-01   , & ! Molecular heat conductivity of water [J m^{-1} s^{-1} K^{-1}]
  tpl_kappa_I       = 2.29       , & ! Molecular heat conductivity of ice [J m^{-1} s^{-1} K^{-1}]
  tpl_kappa_S_min   = 0.2        , & ! Minimum molecular heat conductivity of snow [J m^{-1} s^{-1} K^{-1}]
  tpl_kappa_S_max   = 1.5        , & ! Maximum molecular heat conductivity of snow [J m^{-1} s^{-1} K^{-1}]
  tpl_Gamma_kappa_S = 1.3            ! Empirical parameter [J m^{-2} s^{-1} K^{-1}] 
                                     ! in the expression for the snow heat conductivity 

!==============================================================================

END MODULE flake_parameters

!------------------------------------------------------------------------------

MODULE flake

!------------------------------------------------------------------------------
!
! Description:
!
!  The main program unit of the lake model FLake,  
!  containing most of the FLake procedures.
!  Most FLake variables and local parameters are declared.
!
!  FLake (Fresh-water Lake) is a lake model capable of predicting the surface temperature 
!  in lakes of various depth on the time scales from a few hours to a year.
!  The model is based on a two-layer parametric representation of
!  the evolving temperature profile, where the structure of the stratified layer between the
!  upper mixed layer and the basin bottom, the lake thermocline,
!  is described using the concept of self-similarity of the temperature-depth curve.
!  The concept was put forward by Kitaigorodskii and Miropolsky (1970) 
!  to describe the vertical temperature structure of the oceanic seasonal thermocline.
!  It has since been successfully used in geophysical applications.
!  The concept of self-similarity of the evolving temperature profile
!  is also used to describe the vertical structure of the thermally active upper layer 
!  of bottom sediments and of the ice and snow cover.
!
!  The lake model incorporates the heat budget equations
!  for the four layers in question, viz., snow, ice, water and bottom sediments,
!  developed with due regard for the vertically distributed character
!  of solar radiation heating.
!  The entrainment equation that incorporates the Zilitinkevich (1975) spin-up term
!  is used to compute the depth of a convectively-mixed layer. 
!  A relaxation-type equation is used
!  to compute the wind-mixed layer depth in stable and neutral stratification,
!  where a multi-limit formulation for the equilibrium mixed-layer depth
!  proposed by Zilitinkevich and Mironov (1996)
!  accounts for the effects of the earth's rotation, of the surface buoyancy flux
!  and of the static stability in the thermocline.
!  The equations for the mixed-layer depth are developed with due regard for  
!  the volumetric character of the radiation heating.
!  Simple thermodynamic arguments are invoked to develop
!  the evolution equations for the ice thickness and for the snow thickness.
!  The heat flux through the water-bottom sediment interface is computed,
!  using a parameterization proposed by Golosov et al. (1998).
!  The heat flux trough the air-water interface 
!  (or through the air-ice or air-snow interface)
!  is provided by the driving atmospheric model.
!
!  Empirical constants and parameters of the lake model
!  are estimated, using independent empirical and numerical data.
!  They should not be re-evaluated when the model is applied to a particular lake.
!  The only lake-specific parameters are the lake depth,
!  the optical characteristics of lake water,
!  the temperature at the bottom of the thermally active layer
!  of bottom sediments and the depth of that layer.
!
!  A detailed description of the lake model is given in
!  Mironov, D. V., 2005:
!  Parameterization of Lakes in Numerical Weather Prediction.
!  Part 1: Description of a Lake Model.
!  Manuscript is available from the author.
!  Dmitrii Mironov 
!  German Weather Service, Kaiserleistr. 29/35, D-63067 Offenbach am Main, Germany.
!  dmitrii.mironov@dwd.de 
!
!  Lines embraced with "!_tmp" contain temporary parts of the code.
!  Lines embraced/marked with "!_dev" may be replaced
!  as improved parameterizations are developed and tested.
!  Lines embraced/marked with "!_dm" are DM's comments
!  that may be helpful to a user.
!  Lines embraced/marked with "!_dbg" are used
!  for debugging purposes only.
!

USE data_parameters , ONLY : &
  ireals                   , & ! KIND-type parameter for real variables
  iintegers                    ! KIND-type parameter for "normal" integer variables

use machine,               only: kind_phys
!==============================================================================

IMPLICIT NONE

!==============================================================================
!
! Declarations
!
!  The variables declared below
!  are accessible to all program units of the MODULE flake.
!  Some of them should be USEd by the driving routines that call flake routines.
!  These are basically the quantities computed by FLake.
!  All variables declared below have a suffix "flk".

!  FLake variables of type REAL

!  Temperatures at the previous time step ("p") and the updated temperatures ("n") 
REAL (KIND = kind_phys) ::           &
  T_mnw_p_flk, T_mnw_n_flk      , & ! Mean temperature of the water column [K] 
  T_snow_p_flk, T_snow_n_flk    , & ! Temperature at the air-snow interface [K] 
  T_ice_p_flk, T_ice_n_flk      , & ! Temperature at the snow-ice or air-ice interface [K] 
  T_wML_p_flk, T_wML_n_flk      , & ! Mixed-layer temperature [K] 
  T_bot_p_flk, T_bot_n_flk      , & ! Temperature at the water-bottom sediment interface [K] 
  T_B1_p_flk, T_B1_n_flk            ! Temperature at the bottom of the upper layer of the sediments [K] 

!  Thickness of various layers at the previous time step ("p") and the updated values ("n") 
REAL (KIND = kind_phys) ::           &
  h_snow_p_flk, h_snow_n_flk    , & ! Snow thickness [m]
  h_ice_p_flk, h_ice_n_flk      , & ! Ice thickness [m]
  h_ML_p_flk, h_ML_n_flk        , & ! Thickness of the mixed-layer [m] 
  H_B1_p_flk, H_B1_n_flk            ! Thickness of the upper layer of bottom sediments [m] 

!  The shape factor(s) at the previous time step ("p") and the updated value(s) ("n") 
REAL (KIND = kind_phys) ::           &
  C_T_p_flk, C_T_n_flk          , & ! Shape factor (thermocline)
  C_TT_flk                      , & ! Dimensionless parameter (thermocline)
  C_Q_flk                       , & ! Shape factor with respect to the heat flux (thermocline)
  C_I_flk                       , & ! Shape factor (ice)
  C_S_flk                           ! Shape factor (snow)

!  Derivatives of the shape functions
REAL (KIND = kind_phys) ::           &
  Phi_T_pr0_flk                 , & ! d\Phi_T(0)/d\zeta   (thermocline)
  Phi_I_pr0_flk                 , & ! d\Phi_I(0)/d\zeta_I (ice)
  Phi_I_pr1_flk                 , & ! d\Phi_I(1)/d\zeta_I (ice)
  Phi_S_pr0_flk                     ! d\Phi_S(0)/d\zeta_S (snow)

!  Heat and radiation fluxes
REAL (KIND = kind_phys) ::           &
  Q_snow_flk                    , & ! Heat flux through the air-snow interface [W m^{-2}]
  Q_ice_flk                     , & ! Heat flux through the snow-ice or air-ice interface [W m^{-2}]
  Q_w_flk                       , & ! Heat flux through the ice-water or air-water interface [W m^{-2}]
  Q_bot_flk                     , & ! Heat flux through the water-bottom sediment interface [W m^{-2}]
  I_atm_flk                     , & ! Radiation flux at the lower boundary of the atmosphere [W m^{-2}],
                                    ! i.e. the incident radiation flux with no regard for the surface albedo.
  I_snow_flk                    , & ! Radiation flux through the air-snow interface [W m^{-2}]
  I_ice_flk                     , & ! Radiation flux through the snow-ice or air-ice interface [W m^{-2}]
  I_w_flk                       , & ! Radiation flux through the ice-water or air-water interface [W m^{-2}]
  I_h_flk                       , & ! Radiation flux through the mixed-layer-thermocline interface [W m^{-2}]
  I_bot_flk                     , & ! Radiation flux through the water-bottom sediment interface [W m^{-2}]
  I_intm_0_h_flk                , & ! Mean radiation flux over the mixed layer [W m^{-1}]
  I_intm_h_D_flk                , & ! Mean radiation flux over the thermocline [W m^{-1}]
  I_intm_D_H_flk                , & ! Mean radiation flux over the deeper layer defined by Shaobo Zhang [W m^{-1}]
  I_HH_flk                      , & ! Radiation flux through the bottom of the deeper layer defined by Shaobo Zhang [W m^{-2}]
  Q_star_flk                        ! A generalized heat flux scale [W m^{-2}]

!  Velocity scales
REAL (KIND = kind_phys) ::           &
  u_star_w_flk                  , & ! Friction velocity in the surface layer of lake water [m s^{-1}]
  w_star_sfc_flk                    ! Convective velocity scale, 
                                    ! using a generalized heat flux scale [m s^{-1}]

!  The rate of snow accumulation
REAL (KIND = kind_phys) ::           &
  dMsnowdt_flk                      ! The rate of snow accumulation [kg m^{-2} s^{-1}]
!  The secondary layer temp
REAL (KIND = kind_phys) ::           &
  T_BOT_2_IN_FLK

!==============================================================================
! Procedures 
!==============================================================================

CONTAINS

!==============================================================================
!  The codes of the FLake procedures are stored in separate "*.incf" files
!  and are included below.
!------------------------------------------------------------------------------

!==============================================================================
! include 'flake_radflux.incf'
!------------------------------------------------------------------------------
! changed by Shaobo Zhang

SUBROUTINE flake_radflux ( depth_w, albedo_water, albedo_ice, albedo_snow, & 
                           opticpar_water, opticpar_ice, opticpar_snow,    &
                           depth_bs )       

!------------------------------------------------------------------------------
!
! Description:
!
!  Computes the radiation fluxes 
!  at the snow-ice, ice-water, air-water, 
!  mixed layer-thermocline and water column-bottom sediment interfaces,
!  the mean radiation flux over the mixed layer,
!  and the mean radiation flux over the thermocline.
!
!
! Declarations:
!
! Modules used:

!_dm Parameters are USEd in module "flake".
!_nu USE data_parameters , ONLY : &
!_nu     ireals,                  & ! KIND-type parameter for real variables
!_nu     iintegers                  ! KIND-type parameter for "normal" integer variables

USE flake_derivedtypes          ! Definitions of derived TYPEs

USE flake_parameters , ONLY : & 
  h_Snow_min_flk            , & ! Minimum snow thickness [m]
  h_Ice_min_flk             , & ! Minimum ice thickness [m]
  h_ML_min_flk                  ! Minimum mixed-layer depth [m]

use machine,               only: kind_phys
!==============================================================================

IMPLICIT NONE

!==============================================================================
!
! Declarations

!  Input (procedure arguments)

REAL (KIND = kind_phys), INTENT(IN) ::   &
  depth_w                           , & ! The lake depth [m]
  depth_bs                          , & ! The depth_bs added by Shaobo Zhang
  albedo_water                      , & ! Albedo of the water surface 
  albedo_ice                        , & ! Albedo of the ice surface
  albedo_snow                           ! Albedo of the snow surface

TYPE (opticpar_medium), INTENT(IN) :: & 
  opticpar_water                    , & ! Optical characteristics of water
  opticpar_ice                      , & ! Optical characteristics of ice
  opticpar_snow                         ! Optical characteristics of snow 


!  Local variables of type INTEGER
INTEGER (KIND = iintegers) :: & ! Help variable(s)
  i                             ! DO loop index

!==============================================================================
!  Start calculations
!------------------------------------------------------------------------------

  IF(h_ice_p_flk.GE.h_Ice_min_flk) THEN            ! Ice exists
    IF(h_snow_p_flk.GE.h_Snow_min_flk) THEN        ! There is snow above the ice
      I_snow_flk = I_atm_flk*(1.0-albedo_snow) 
      I_bot_flk = 0.0
      DO i=1, opticpar_snow%nband_optic
        I_bot_flk = I_bot_flk +                    & 
        opticpar_snow%frac_optic(i)*EXP(-opticpar_snow%extincoef_optic(i)*h_snow_p_flk) 
      END DO 
      I_ice_flk  = I_snow_flk*I_bot_flk
    ELSE                                           ! No snow above the ice 
      I_snow_flk = I_atm_flk  
      I_ice_flk  = I_atm_flk*(1.0-albedo_ice)
    END IF 
    I_bot_flk = 0.0
    DO i=1, opticpar_ice%nband_optic
      I_bot_flk = I_bot_flk +                      & 
      opticpar_ice%frac_optic(i)*EXP(-opticpar_ice%extincoef_optic(i)*h_ice_p_flk) 
    END DO 
    I_w_flk      = I_ice_flk*I_bot_flk
  ELSE                                             ! No ice-snow cover
    I_snow_flk   = I_atm_flk  
    I_ice_flk    = I_atm_flk
    I_w_flk      = I_atm_flk*(1.0-albedo_water)
  END IF 

  IF(h_ML_p_flk.GE.h_ML_min_flk) THEN           ! Radiation flux at the bottom of the mixed layer
    I_bot_flk = 0.0
    DO i=1, opticpar_water%nband_optic
      I_bot_flk = I_bot_flk +            & 
      opticpar_water%frac_optic(i)*EXP(-opticpar_water%extincoef_optic(i)*h_ML_p_flk) 
!      print*,'nband_optic=',opticpar_water%nband_optic
!      print*,'Extinction=',opticpar_water%extincoef_optic(i)
    END DO 
    I_h_flk = I_w_flk*I_bot_flk
  ELSE                                          ! Mixed-layer depth is less then a minimum value
    I_h_flk = I_w_flk
  END IF

  I_bot_flk = 0.0                         ! Radiation flux at the lake bottom
  DO i=1, opticpar_water%nband_optic
    I_bot_flk = I_bot_flk +              & 
    opticpar_water%frac_optic(i)*EXP(-opticpar_water%extincoef_optic(i)*depth_w) 
  END DO 
  I_bot_flk = I_w_flk*I_bot_flk

  IF(h_ML_p_flk.GE.h_ML_min_flk) THEN           ! Integral-mean radiation flux over the mixed layer
    I_intm_0_h_flk = 0.0
    DO i=1, opticpar_water%nband_optic
      I_intm_0_h_flk = I_intm_0_h_flk +                                &
      opticpar_water%frac_optic(i)/opticpar_water%extincoef_optic(i)*  &
      (1.0 - EXP(-opticpar_water%extincoef_optic(i)*h_ML_p_flk))
    END DO 
    I_intm_0_h_flk = I_w_flk*I_intm_0_h_flk/h_ML_p_flk
  ELSE
    I_intm_0_h_flk = I_h_flk
  END IF

  IF(h_ML_p_flk.LE.depth_w-h_ML_min_flk) THEN   ! Integral-mean radiation flux over the thermocline
    I_intm_h_D_flk = 0.0 
    DO i=1, opticpar_water%nband_optic
      I_intm_h_D_flk = I_intm_h_D_flk +                                &
      opticpar_water%frac_optic(i)/opticpar_water%extincoef_optic(i)*  &
      ( EXP(-opticpar_water%extincoef_optic(i)*h_ML_p_flk)             &
      - EXP(-opticpar_water%extincoef_optic(i)*depth_w) )
    END DO 
    I_intm_h_D_flk = I_w_flk*I_intm_h_D_flk/(depth_w-h_ML_p_flk)
  ELSE
    I_intm_h_D_flk = I_h_flk
  END IF

! Added by Shaobo Zhang

  IF(depth_bs.GE.h_ML_min_flk) THEN! Integral-mean radiation flux over the deeper layer defined by Shaobo Zhang
    I_intm_D_H_flk = 0.0 
    DO i=1, opticpar_water%nband_optic
      I_intm_D_H_flk = I_intm_D_H_flk +                                &
      opticpar_water%frac_optic(i)/opticpar_water%extincoef_optic(i)*  &
      ( EXP(-opticpar_water%extincoef_optic(i)*depth_w)             &
      - EXP(-opticpar_water%extincoef_optic(i)*(depth_w+depth_bs)) )
    END DO 
    I_intm_D_H_flk = I_w_flk*I_intm_D_H_flk/depth_bs
  ELSE
    I_intm_D_H_flk = I_bot_flk
  END IF

! Radiation flux at the bottom of the deeper layer defined by Shaobo Zhang
    I_HH_flk = 0.0
    DO i=1, opticpar_water%nband_optic
      I_HH_flk = I_HH_flk +            & 
      opticpar_water%frac_optic(i)*EXP(-opticpar_water%extincoef_optic(i)*(depth_w+depth_bs)) 
    END DO 
    I_HH_flk = I_w_flk*I_HH_flk

!------------------------------------------------------------------------------
!  End calculations
!==============================================================================

END SUBROUTINE flake_radflux

!==============================================================================

!==============================================================================
! include 'flake_main.incf'
!------------------------------------------------------------------------------

SUBROUTINE flake_main ( depthw, depthbs, T_bs, par_Coriolis,       &
                        extincoef_water_typ,                       &
                        del_time, T_sfc_p, T_sfc_n, T_bot_2_in,    &
                        T_bot_2_out  )         

!------------------------------------------------------------------------------
!
! Description:
!
!  The main driving routine of the lake model FLake 
!  where computations are performed.
!  Advances the surface temperature
!  and other FLake variables one time step.
!  At the moment, the Euler explicit scheme is used.
!
!  Lines embraced with "!_tmp" contain temporary parts of the code.
!  Lines embraced/marked with "!_dev" may be replaced
!  as improved parameterizations are developed and tested.
!  Lines embraced/marked with "!_dm" are DM's comments
!  that may be helpful to a user.
!  Lines embraced/marked with "!_dbg" are used 
!  for debugging purposes only.
!
! Declarations:
!
! Modules used:

!_dm Parameters are USEd in module "flake".
!_nu USE data_parameters , ONLY : &
!_nu     ireals,                  & ! KIND-type parameter for real variables
!_nu     iintegers                  ! KIND-type parameter for "normal" integer variables

USE flake_parameters            ! Thermodynamic parameters and dimensionless constants of FLake

USE flake_configure             ! Switches and parameters that configure FLake

use machine,               only: kind_phys
! ADDED by Shaobo Zhang
! USE mod_dynparam, only : lake_depth_max 

!==============================================================================

IMPLICIT NONE

!==============================================================================
!
! Declarations

!  Input (procedure arguments)

! changed by Shaobo Zhang
REAL (KIND = kind_phys), INTENT(IN) ::   &
  depthw                           , & ! The lake depth [m]
  depthbs                          , & ! Depth of the thermally active layer of bottom sediments [m]
  T_bs                              , & ! Temperature at the outer edge of 
                                        ! the thermally active layer of bottom sediments [K]
  par_Coriolis                      , & ! The Coriolis parameter [s^{-1}]
  extincoef_water_typ               , & ! "Typical" extinction coefficient of the lake water [m^{-1}],
                                        ! used to compute the equilibrium CBL depth
  del_time                          , & ! The model time step [s]
  T_sfc_p                           , & ! Surface temperature at the previous time step [K]  
  T_bot_2_in

REAL (KIND = kind_phys)             ::   &
  depth_w                           , & ! The lake depth [m]
  depth_bs                              ! Depth of the thermally active layer of bottom sediments [m]

!  Output (procedure arguments)

REAL (KIND = kind_phys), INTENT(OUT) ::  &
  T_sfc_n                              , &  ! Updated surface temperature [K] 
                                        ! (equal to the updated value of either T_ice, T_snow or T_wML)
  T_bot_2_out


!  Local variables of type LOGICAL
LOGICAL ::          &
  l_ice_create    , & ! Switch, .TRUE. = ice does not exist but should be created
  l_snow_exists   , & ! Switch, .TRUE. = there is snow above the ice
  l_ice_meltabove     ! Switch, .TRUE. = snow/ice melting from above takes place

!  Local variables of type INTEGER
INTEGER (KIND = iintegers) :: &
  i                             ! Loop index

!  Local variables of type REAL
REAL (KIND = kind_phys) ::    &
  d_T_mnw_dt             , & ! Time derivative of T_mnw [K s^{-1}] 
  d_T_ice_dt             , & ! Time derivative of T_ice [K s^{-1}] 
  d_T_bot_dt             , & ! Time derivative of T_bot [K s^{-1}] 
  d_T_B1_dt              , & ! Time derivative of T_B1 [K s^{-1}] 
  d_h_snow_dt            , & ! Time derivative of h_snow [m s^{-1}]
  d_h_ice_dt             , & ! Time derivative of h_ice [m s^{-1}]
  d_h_ML_dt              , & ! Time derivative of h_ML [m s^{-1}]
  d_H_B1_dt              , & ! Time derivative of H_B1 [m s^{-1}]
  d_h_D_dt               , & ! Time derivative of h_D, new defined by Shaobo Zhang
  d_T_H_dt               , & ! Time derivative of T_H, new defined by Shaobo Zhang
  d_C_T_dt                   ! Time derivative of C_T [s^{-1}]

!  Local variables of type REAL
REAL (KIND = kind_phys) ::    &
  N_T_mean               , & ! The mean buoyancy frequency in the thermocline [s^{-1}] 
  tmp                    , & ! temperary variable
  ZM_h_scale             , & ! The ZM96 equilibrium SBL depth scale [m] 
  conv_equil_h_scale         ! The equilibrium CBL depth scale [m]

!  Local variables of type REAL
REAL (KIND = kind_phys) :: &
  h_ice_threshold     , & ! If h_ice<h_ice_threshold, use quasi-equilibrium ice model 
  flk_str_1           , & ! Help storage variable
  flk_str_2           , & ! Help storage variable
  R_H_icesnow         , & ! Dimensionless ratio, used to store intermediate results
  R_rho_c_icesnow     , & ! Dimensionless ratio, used to store intermediate results
  R_TI_icesnow        , & ! Dimensionless ratio, used to store intermediate results
  R_Tstar_icesnow         ! Dimensionless ratio, used to store intermediate results

! ADDED by Shaobo Zhang
!REAL (KIND = kind_phys) :: T_bot_2_in, T_bot_2_out
REAL (KIND = kind_phys) :: CT, CTT, CQ
REAL (KIND = kind_phys) :: lake_depth_max

!==============================================================================
!  Start calculations
!------------------------------------------------------------------------------

!_dm 
! Security. Set time-rate-of-change of prognostic variables to zero.
! Set prognostic variables to their values at the previous time step.
! (This is to avoid spurious changes of prognostic variables 
! when FLake is used within a 3D model, e.g. to avoid spurious generation of ice 
! at the neighbouring lake points as noticed by Burkhardt Rockel.)
!_dm 

!  print*,'T_sfc_p=',T_sfc_p,' T_bot_2_in=',T_bot_2_in 
!  print*,'h_snow_p_flk=',h_snow_p_flk,' h_ice_p_flk=',h_ice_p_flk
!  print*,'T_snow_p_flk=',T_snow_p_flk, ' T_ice_p_flk=',T_ice_p_flk
!  print*,'T_wML_p_flk=',T_wML_p_flk, ' T_mnw_p_flk=',T_mnw_p_flk 
!  print*,'T_bot_p_flk=',T_bot_p_flk, ' T_B1_p_flk=',T_B1_p_flk
!  print*,'h_ML_p_flk=',h_ML_p_flk, ' H_B1_p_flk=',H_B1_p_flk
!  print*,'C_T_p_flk=',C_T_p_flk  
!  print*,'depthw= ', depthw
!  print*,'depthbs=',depthbs

lake_depth_max = 60.0
depth_w =depthw
depth_bs = depthbs
d_T_mnw_dt   = 0.0 
d_T_ice_dt   = 0.0 
d_T_bot_dt   = 0.0 
d_T_B1_dt    = 0.0 
d_h_snow_dt  = 0.0 
d_h_ice_dt   = 0.0 
d_h_ML_dt    = 0.0 
d_H_B1_dt    = 0.0 
d_C_T_dt     = 0.0 
T_snow_n_flk = T_snow_p_flk   
T_ice_n_flk  = T_ice_p_flk    
T_wML_n_flk  = T_wML_p_flk   
T_mnw_n_flk  = T_mnw_p_flk     
T_bot_n_flk  = T_bot_p_flk  
T_B1_n_flk   = T_B1_p_flk      
h_snow_n_flk = h_snow_p_flk 
h_ice_n_flk  = h_ice_p_flk   
h_ML_n_flk   = h_ML_p_flk    
H_B1_n_flk   = H_B1_p_flk   
C_T_n_flk    = C_T_p_flk    

!------------------------------------------------------------------------------
!  Compute fluxes, using variables from the previous time step.
!------------------------------------------------------------------------------

!_dm
! At this point, the heat and radiation fluxes, namely,
! Q_snow_flk, Q_ice_flk, Q_w_flk, 
! I_atm_flk, I_snow_flk, I_ice_flk, I_w_flk, I_h_flk, I_bot_flk,     
! the mean radiation flux over the mixed layer, I_intm_0_h_flk, 
! and the mean radiation flux over the thermocline, I_intm_h_D_flk, 
! should be known.
! They are computed within "flake_interface" (or within the driving model)
! and are available to "flake_main"
! through the above variables declared in the MODULE "flake".
! In case a lake is ice-covered, Q_w_flk is re-computed below.
!_dm

! Heat flux through the ice-water interface
IF(h_ice_p_flk.GE.h_Ice_min_flk) THEN    ! Ice exists 
  IF(h_ML_p_flk.LE.h_ML_min_flk) THEN    ! Mixed-layer depth is zero, compute flux 
    Q_w_flk = -tpl_kappa_w*(T_bot_p_flk-T_wML_p_flk)/depth_w  ! Flux with linear T(z) 
    Phi_T_pr0_flk = Phi_T_pr0_1*C_T_p_flk-Phi_T_pr0_2         ! d\Phi(0)/d\zeta (thermocline)
    Q_w_flk = Q_w_flk*MAX(Phi_T_pr0_flk, 1.0)           ! Account for an increased d\Phi(0)/d\zeta 
  ELSE                    
    Q_w_flk = 0.0          ! Mixed-layer depth is greater than zero, set flux to zero
  END IF   
ELSE                       ! If Ice doesn't exist, set flux to zero, YWu  2019
  Q_w_flk = 0.0
END IF   

! A generalized heat flux scale 
Q_star_flk = Q_w_flk + I_w_flk + I_h_flk - 2.0*I_intm_0_h_flk

! changed by Shaobo Zhang
! Heat flux through the water-bottom sediment interface
!IF(lflk_botsed_use) THEN
!  Q_bot_flk = -tpl_kappa_w*(T_B1_p_flk-T_bot_p_flk)/MAX(H_B1_p_flk, H_B1_min_flk)*Phi_B1_pr0
!ELSE  
!  Q_bot_flk = 0.0   ! The bottom-sediment scheme is not used
!END IF

IF(lflk_botsed_use) THEN   ! The bottom-sediment scheme is used
  Q_bot_flk = -tpl_kappa_w*(T_B1_p_flk-T_bot_p_flk)/MAX(H_B1_p_flk, H_B1_min_flk)*Phi_B1_pr0
ELSE  ! The scheme written by Shaobo Zhang for the deeper layer of a deep lake is used
  Q_bot_flk = -tpl_kappa_w*(T_bot_2_in-T_bot_p_flk)/MAX(depth_bs, H_B1_min_flk)*Phi_B1_pr0
END IF

!------------------------------------------------------------------------------
!  Check if ice exists or should be created.
!  If so, compute the thickness and the temperature of ice and snow.
!------------------------------------------------------------------------------

!_dm
! Notice that a quasi-equilibrium ice-snow model is used 
! to avoid numerical instability when the ice is thin.
! This is always the case when new ice is created.
!_dm

!_dev
! The dependence of snow density and of snow heat conductivity 
! on the snow thickness is accounted for parametrically.
! That is, the time derivatives of \rho_S and \kappa_S are neglected.
! The exception is the equation for the snow thickness 
! in case of snow accumulation and no melting, 
! where d\rho_S/dt is incorporated.
! Furthermore, some (presumably small) correction terms incorporating 
! the snow density and the snow heat conductivity are dropped out.
! Those terms may be included as better formulations 
! for \rho_S and \kappa_S are available.
!_dev

! Default values
l_ice_create    = .FALSE.  
l_ice_meltabove = .FALSE.  

Ice_exist: IF(h_ice_p_flk.LT.h_Ice_min_flk) THEN   ! Ice does not exist 

  l_ice_create = T_wML_p_flk.LE.(tpl_T_f+c_small_flk).AND.Q_w_flk.LT.0.0
  IF(l_ice_create) THEN                            ! Ice does not exist but should be created
    d_h_ice_dt = -Q_w_flk/tpl_rho_I/tpl_L_f                                  
    h_ice_n_flk = h_ice_p_flk + d_h_ice_dt*del_time                          ! Advance h_ice 
    T_ice_n_flk = tpl_T_f + h_ice_n_flk*Q_w_flk/tpl_kappa_I/Phi_I_pr0_lin    ! Ice temperature
    d_h_snow_dt = dMsnowdt_flk/tpl_rho_S_min 
    h_snow_n_flk = h_snow_p_flk + d_h_snow_dt*del_time                       ! Advance h_snow
    Phi_I_pr1_flk = Phi_I_pr1_lin                                    & 
                  + Phi_I_ast_MR*MIN(1.0, h_ice_n_flk/H_Ice_max)       ! d\Phi_I(1)/d\zeta_I (ice)
    R_H_icesnow = Phi_I_pr1_flk/Phi_S_pr0_lin*tpl_kappa_I/flake_snowheatconduct(h_snow_n_flk) &
                * h_snow_n_flk/MAX(h_ice_n_flk, h_Ice_min_flk)
    T_snow_n_flk = T_ice_n_flk + R_H_icesnow*(T_ice_n_flk-tpl_T_f)           ! Snow temperature
  END IF

ELSE Ice_exist                                     ! Ice exists

  l_snow_exists = h_snow_p_flk.GE.h_Snow_min_flk   ! Check if there is snow above the ice

  Melting: IF(T_snow_p_flk.GE.(tpl_T_f-c_small_flk)) THEN  ! T_sfc = T_f, check for melting from above
                                                           ! T_snow = T_ice if snow is absent 
    IF(l_snow_exists) THEN   ! There is snow above the ice
      flk_str_1 = Q_snow_flk + I_snow_flk - I_ice_flk        ! Atmospheric forcing
      IF(flk_str_1.GE.0.0) THEN  ! Melting of snow and ice from above
        l_ice_meltabove = .TRUE.
        d_h_snow_dt = (-flk_str_1/tpl_L_f+dMsnowdt_flk)/flake_snowdensity(h_snow_p_flk)
        d_h_ice_dt  = -(I_ice_flk - I_w_flk - Q_w_flk)/tpl_L_f/tpl_rho_I 
      END IF 
    ELSE                     ! No snow above the ice
      flk_str_1 = Q_ice_flk + I_ice_flk - I_w_flk - Q_w_flk  ! Atmospheric forcing + heating from the water
      IF(flk_str_1.GE.0.0) THEN  ! Melting of ice from above, snow accumulation may occur
        l_ice_meltabove = .TRUE.
        d_h_ice_dt  = -flk_str_1/tpl_L_f/tpl_rho_I 
        d_h_snow_dt = dMsnowdt_flk/tpl_rho_S_min
      END IF 
    END IF 
    IF(l_ice_meltabove) THEN  ! Melting from above takes place
      h_ice_n_flk  = h_ice_p_flk  + d_h_ice_dt *del_time  ! Advance h_ice
      h_snow_n_flk = h_snow_p_flk + d_h_snow_dt*del_time  ! Advance h_snow
      T_ice_n_flk  = tpl_T_f                              ! Set T_ice to the freezing point
      T_snow_n_flk = tpl_T_f                              ! Set T_snow to the freezing point
    END IF

  END IF Melting

  No_Melting: IF(.NOT.l_ice_meltabove) THEN                 ! No melting from above

    d_h_snow_dt = flake_snowdensity(h_snow_p_flk)  
    IF(d_h_snow_dt.LT.tpl_rho_S_max) THEN    ! Account for d\rho_S/dt
     flk_str_1 = h_snow_p_flk*tpl_Gamma_rho_S/tpl_rho_w_r
     flk_str_1 = flk_str_1/(1.0-flk_str_1)
    ELSE                                     ! Snow density is equal to its maximum value, d\rho_S/dt=0
     flk_str_1 = 0.0
    END IF
    d_h_snow_dt = dMsnowdt_flk/d_h_snow_dt/(1.0+flk_str_1)       ! Snow accumulation
    h_snow_n_flk = h_snow_p_flk + d_h_snow_dt*del_time                         ! Advance h_snow
    
    Phi_I_pr0_flk = h_ice_p_flk/H_Ice_max                              ! h_ice relative to its maximum value
    C_I_flk = C_I_lin - C_I_MR*(1.0+Phi_I_ast_MR)*Phi_I_pr0_flk  ! Shape factor (ice)
    Phi_I_pr1_flk = Phi_I_pr1_lin + Phi_I_ast_MR*Phi_I_pr0_flk         ! d\Phi_I(1)/d\zeta_I (ice)
    Phi_I_pr0_flk = Phi_I_pr0_lin - Phi_I_pr0_flk                      ! d\Phi_I(0)/d\zeta_I (ice)

    h_ice_threshold = MAX(1.0, 2.0*C_I_flk*tpl_c_I*(tpl_T_f-T_ice_p_flk)/tpl_L_f)
    h_ice_threshold = Phi_I_pr0_flk/C_I_flk*tpl_kappa_I/tpl_rho_I/tpl_c_I*h_ice_threshold
    h_ice_threshold = SQRT(h_ice_threshold*del_time)                   ! Threshold value of h_ice
    h_ice_threshold = MIN(0.9*H_Ice_max, MAX(h_ice_threshold, h_Ice_min_flk))
                                                                       ! h_ice(threshold) < 0.9*H_Ice_max

    IF(h_ice_p_flk.LT.h_ice_threshold) THEN  ! Use a quasi-equilibrium ice model

      IF(l_snow_exists) THEN   ! Use fluxes at the air-snow interface
        flk_str_1 = Q_snow_flk + I_snow_flk - I_w_flk
      ELSE                     ! Use fluxes at the air-ice interface
        flk_str_1 = Q_ice_flk + I_ice_flk - I_w_flk
      END IF
      d_h_ice_dt = -(flk_str_1-Q_w_flk)/tpl_L_f/tpl_rho_I
      h_ice_n_flk = h_ice_p_flk + d_h_ice_dt *del_time                         ! Advance h_ice
      T_ice_n_flk = tpl_T_f + h_ice_n_flk*flk_str_1/tpl_kappa_I/Phi_I_pr0_flk  ! Ice temperature

    ELSE                                     ! Use a complete ice model

      d_h_ice_dt = tpl_kappa_I*(tpl_T_f-T_ice_p_flk)/h_ice_p_flk*Phi_I_pr0_flk
      d_h_ice_dt = (Q_w_flk+d_h_ice_dt)/tpl_L_f/tpl_rho_I
      h_ice_n_flk = h_ice_p_flk  + d_h_ice_dt*del_time                         ! Advance h_ice

      R_TI_icesnow = tpl_c_I*(tpl_T_f-T_ice_p_flk)/tpl_L_f         ! Dimensionless parameter
      R_Tstar_icesnow = 1.0 - C_I_flk                        ! Dimensionless parameter
      IF(l_snow_exists) THEN  ! There is snow above the ice
        R_H_icesnow = Phi_I_pr1_flk/Phi_S_pr0_lin*tpl_kappa_I/flake_snowheatconduct(h_snow_p_flk) &
                    * h_snow_p_flk/h_ice_p_flk
        R_rho_c_icesnow = flake_snowdensity(h_snow_p_flk)*tpl_c_S/tpl_rho_I/tpl_c_I 
!_dev 
!_dm 
! These terms should be included as an improved understanding of the snow scheme is gained, 
! of the effect of snow density in particular. 
!_dm 
!_nu        R_Tstar_icesnow = R_Tstar_icesnow                                                           &
!_nu                        + (1.0+C_S_lin*h_snow_p_flk/h_ice_p_flk)*R_H_icesnow*R_rho_c_icesnow
!_dev

        R_Tstar_icesnow = R_Tstar_icesnow*R_TI_icesnow             ! Dimensionless parameter

!_dev
!_nu        R_Tstar_icesnow = R_Tstar_icesnow                                                         &
!_nu                        + (1.0-R_rho_c_icesnow)*tpl_c_I*T_ice_p_flk/tpl_L_f
!_dev
        flk_str_2 = Q_snow_flk+I_snow_flk-I_w_flk                  ! Atmospheric fluxes
        flk_str_1  = C_I_flk*h_ice_p_flk + (1.0+C_S_lin*R_H_icesnow)*R_rho_c_icesnow*h_snow_p_flk
        d_T_ice_dt = -(1.0-2.0*C_S_lin)*R_H_icesnow*(tpl_T_f-T_ice_p_flk)             & 
                   * tpl_c_S*dMsnowdt_flk                          ! Effect of snow accumulation
      ELSE                    ! No snow above the ice
        R_Tstar_icesnow = R_Tstar_icesnow*R_TI_icesnow             ! Dimensionless parameter
        flk_str_2 = Q_ice_flk+I_ice_flk-I_w_flk                    ! Atmospheric fluxes
        flk_str_1  = C_I_flk*h_ice_p_flk
        d_T_ice_dt = 0.0
      END IF 
      d_T_ice_dt = d_T_ice_dt + tpl_kappa_I*(tpl_T_f-T_ice_p_flk)/h_ice_p_flk*Phi_I_pr0_flk       &
                 * (1.0-R_Tstar_icesnow)                     ! Add flux due to heat conduction
      d_T_ice_dt = d_T_ice_dt - R_Tstar_icesnow*Q_w_flk            ! Add flux from water to ice
      d_T_ice_dt = d_T_ice_dt + flk_str_2                          ! Add atmospheric fluxes
      d_T_ice_dt = d_T_ice_dt/tpl_rho_I/tpl_c_I                    ! Total forcing
      d_T_ice_dt = d_T_ice_dt/flk_str_1                            ! dT_ice/dt 
      T_ice_n_flk = T_ice_p_flk + d_T_ice_dt*del_time                          ! Advance T_ice
    END IF

    Phi_I_pr1_flk = MIN(1.0, h_ice_n_flk/H_Ice_max)          ! h_ice relative to its maximum value
    Phi_I_pr1_flk = Phi_I_pr1_lin + Phi_I_ast_MR*Phi_I_pr1_flk     ! d\Phi_I(1)/d\zeta_I (ice)
    R_H_icesnow = Phi_I_pr1_flk/Phi_S_pr0_lin*tpl_kappa_I/flake_snowheatconduct(h_snow_n_flk) &
                 *h_snow_n_flk/MAX(h_ice_n_flk, h_Ice_min_flk)
    T_snow_n_flk = T_ice_n_flk + R_H_icesnow*(T_ice_n_flk-tpl_T_f)             ! Snow temperature

  END IF No_Melting

END IF Ice_exist   

! Security, limit h_ice by its maximum value
h_ice_n_flk = MIN(h_ice_n_flk, H_Ice_max)      

! Security, limit the ice and snow temperatures by the freezing point 
T_snow_n_flk = MIN(T_snow_n_flk, tpl_T_f)  
T_ice_n_flk =  MIN(T_ice_n_flk,  tpl_T_f)    

!_tmp
! Security, avoid too low values (these constraints are used for debugging purposes)
!  T_snow_n_flk = MAX(T_snow_n_flk, 73.15)  
!  T_ice_n_flk =  MAX(T_ice_n_flk,  73.15)  
!   The lowest natural temperature ever
!  directly recorded at ground level on Earth is -89.2 C (-128.6 F; 184.0 K)
!  at the Soviet Vostok Station in Antarctica on 21 July, 1983 by ground
!  measurements.
!  https://en.wikipedia.org/wiki/Lowest_temperature_recorded_on_Earth
 
  T_snow_n_flk = MAX(T_snow_n_flk, 184.0)
  T_ice_n_flk =  MAX(T_ice_n_flk,  184.0)
!_tmp

! Remove too thin ice and/or snow
IF(h_ice_n_flk.LT.h_Ice_min_flk)  THEN        ! Check ice
  h_ice_n_flk = 0.0       ! Ice is too thin, remove it, and
  T_ice_n_flk = tpl_T_f         ! set T_ice to the freezing point.
  h_snow_n_flk = 0.0      ! Remove snow when there is no ice, and
  T_snow_n_flk = tpl_T_f        ! set T_snow to the freezing point.
  l_ice_create = .FALSE.        ! "Exotic" case, ice has been created but proved to be too thin
ELSE IF(h_snow_n_flk.LT.h_Snow_min_flk) THEN  ! Ice exists, check snow
  h_snow_n_flk = 0.0      ! Snow is too thin, remove it, 
  T_snow_n_flk = T_ice_n_flk    ! and set the snow temperature equal to the ice temperature.
END IF


!------------------------------------------------------------------------------
!  Compute the mean temperature of the water column.
!------------------------------------------------------------------------------

IF(l_ice_create) Q_w_flk = 0.0     ! Ice has just been created, set Q_w to zero
d_T_mnw_dt = (Q_w_flk - Q_bot_flk + I_w_flk - I_bot_flk)/tpl_rho_w_r/tpl_c_w/depth_w

!print*,'d_T_mnw_dt= ',d_T_mnw_dt
!print*,'Q_w_flk=',Q_w_flk
!print*,'Q_bot_flk=',Q_bot_flk
!print*,'I_w_flk=',I_w_flk
!print*,'I_bot_flk=',I_bot_flk
!print*,'tpl_rho_w_r=',tpl_rho_w_r
!print*,'tpl_c_w=',tpl_c_w
!print*,'depth_w=',depth_w

T_mnw_n_flk = T_mnw_p_flk + d_T_mnw_dt*del_time   ! Advance T_mnw

T_mnw_n_flk = MAX(T_mnw_n_flk, tpl_T_f)           ! Limit T_mnw by the freezing point 


!------------------------------------------------------------------------------
!  Compute the mixed-layer depth, the mixed-layer temperature, 
!  the bottom temperature and the shape factor
!  with respect to the temperature profile in the thermocline. 
!  Different formulations are used, depending on the regime of mixing. 
!------------------------------------------------------------------------------

HTC_Water: IF(h_ice_n_flk.GE.h_Ice_min_flk) THEN    ! Ice exists

  T_mnw_n_flk = MIN(T_mnw_n_flk, tpl_T_r) ! Limit the mean temperature under the ice by T_r 

  T_wML_n_flk = tpl_T_f                   ! The mixed-layer temperature is equal to the freezing point 

  IF(l_ice_create) THEN                  ! Ice has just been created 
    IF(h_ML_p_flk.GE.depth_w-h_ML_min_flk) THEN    ! h_ML=D when ice is created 
      h_ML_n_flk = 0.0                 ! Set h_ML to zero 
      C_T_n_flk = C_T_min                    ! Set C_T to its minimum value 
    ELSE                                          ! h_ML<D when ice is created 
      h_ML_n_flk = h_ML_p_flk                ! h_ML remains unchanged 
      C_T_n_flk = C_T_p_flk                  ! C_T (thermocline) remains unchanged 
    END IF 
    T_bot_n_flk = T_wML_n_flk - (T_wML_n_flk-T_mnw_n_flk)/C_T_n_flk/(1.0-h_ML_n_flk/depth_w)
                                             ! Update the bottom temperature 

  ELSE IF(T_bot_p_flk.LT.tpl_T_r) THEN   ! Ice exists and T_bot < T_r, molecular heat transfer 
    h_ML_n_flk = h_ML_p_flk                  ! h_ML remains unchanged 
    C_T_n_flk = C_T_p_flk                    ! C_T (thermocline) remains unchanged 
    T_bot_n_flk = T_wML_n_flk - (T_wML_n_flk-T_mnw_n_flk)/C_T_n_flk/(1.0-h_ML_n_flk/depth_w)
                                             ! Update the bottom temperature 

  ELSE                                   ! Ice exists and T_bot = T_r, convection due to bottom heating 
    T_bot_n_flk = tpl_T_r                      ! T_bot is equal to the temperature of maximum density 
    IF(h_ML_p_flk.GE.c_small_flk) THEN   ! h_ML > 0 
      C_T_n_flk = C_T_p_flk                     ! C_T (thermocline) remains unchanged 
      h_ML_n_flk = depth_w*(1.0-(T_wML_n_flk-T_mnw_n_flk)/(T_wML_n_flk-T_bot_n_flk)/C_T_n_flk)
      h_ML_n_flk = MAX(h_ML_n_flk, 0.0)   ! Update the mixed-layer depth  
    ELSE                                 ! h_ML = 0 
      h_ML_n_flk = h_ML_p_flk                   ! h_ML remains unchanged 
      C_T_n_flk = (T_wML_n_flk-T_mnw_n_flk)/(T_wML_n_flk-T_bot_n_flk) 
      C_T_n_flk = MIN(C_T_max, MAX(C_T_n_flk, C_T_min)) ! Update the shape factor (thermocline)  
    END IF 
  END IF 

  T_bot_n_flk = MIN(T_bot_n_flk, tpl_T_r)    ! Security, limit the bottom temperature by T_r 

ELSE HTC_Water                                      ! Open water

! Generalised buoyancy flux scale and convective velocity scale
  flk_str_1 = flake_buoypar(T_wML_p_flk)*Q_star_flk/tpl_rho_w_r/tpl_c_w                    
  IF(flk_str_1.LT.0.0) THEN       
    w_star_sfc_flk = (-flk_str_1*h_ML_p_flk)**(1.0/3.0)  ! Convection     
  ELSE 
    w_star_sfc_flk = 0.0                                       ! Neutral or stable stratification
  END IF 

!_dm
! The equilibrium depth of the CBL due to surface cooling with the volumetric heating
! is not computed as a solution to the transcendental equation.
! Instead, an algebraic formula is used
! that interpolates between the two asymptotic limits.
!_dm
  conv_equil_h_scale = -Q_w_flk/MAX(I_w_flk, c_small_flk)
  IF(conv_equil_h_scale.GT.0.0 .AND. conv_equil_h_scale.LT.1.0  &
    .AND. T_wML_p_flk.GT.tpl_T_r) THEN   ! The equilibrium CBL depth scale is only used above T_r
    conv_equil_h_scale = SQRT(6.0*conv_equil_h_scale)                 &
                       + 2.0*conv_equil_h_scale/(1.0-conv_equil_h_scale)
    conv_equil_h_scale = MIN(depth_w, conv_equil_h_scale/extincoef_water_typ)
!    print*,'extincoef_water_typ=',extincoef_water_typ
  ELSE
    conv_equil_h_scale = 0.0       ! Set the equilibrium CBL depth to zero
  END IF

! Mean buoyancy frequency in the thermocline
  N_T_mean = flake_buoypar(0.5*(T_wML_p_flk+T_bot_p_flk))*(T_wML_p_flk-T_bot_p_flk)
  IF(h_ML_p_flk.LE.depth_w-h_ML_min_flk) THEN
    tmp = N_T_mean/(depth_w-h_ML_p_flk)
    N_T_mean = SQRT(abs(tmp))  ! Compute N                   
  ELSE 
    N_T_mean = 0.0                            ! h_ML=D, set N to zero
  END IF 


! The rate of change of C_T
  d_C_T_dt = MAX(w_star_sfc_flk, u_star_w_flk, u_star_min_flk)**2_iintegers
  d_C_T_dt = N_T_mean*(depth_w-h_ML_p_flk)**2_iintegers       &
           / c_relax_C/d_C_T_dt                               ! Relaxation time scale for C_T
  d_C_T_dt = (C_T_max-C_T_min)/MAX(d_C_T_dt, c_small_flk)     ! Rate-of-change of C_T 

! Compute the shape factor and the mixed-layer depth, 
! using different formulations for convection and wind mixing

  C_TT_flk = C_TT_1*C_T_p_flk-C_TT_2         ! C_TT, using C_T at the previous time step
  C_Q_flk = 2.0*C_TT_flk/C_T_p_flk     ! C_Q using C_T at the previous time step

  Mixing_regime: IF(flk_str_1.LT.0.0) THEN  ! Convective mixing 

    C_T_n_flk = C_T_p_flk + d_C_T_dt*del_time                        ! Update C_T, assuming dh_ML/dt>0
    C_T_n_flk = MIN(C_T_max, MAX(C_T_n_flk, C_T_min))                ! Limit C_T 
    d_C_T_dt = (C_T_n_flk-C_T_p_flk)/del_time                        ! Re-compute dC_T/dt

    IF(h_ML_p_flk.LE.depth_w-h_ML_min_flk) THEN       ! Compute dh_ML/dt
      IF(h_ML_p_flk.LE.h_ML_min_flk) THEN    ! Use a reduced entrainment equation (spin-up)
        d_h_ML_dt = c_cbl_1/c_cbl_2*MAX(w_star_sfc_flk, c_small_flk)

!_dbg
! print*, ' FLake: reduced entrainment eq. D_time*d_h_ML_dt  = ', d_h_ML_dt*del_time
! print*, '         w_*       = ', w_star_sfc_flk
! print*, '         \beta*Q_* = ', flk_str_1
!_dbg

      ELSE                                   ! Use a complete entrainment equation 
        R_H_icesnow     = depth_w/h_ML_p_flk
        R_rho_c_icesnow = R_H_icesnow-1.0
        R_TI_icesnow    = C_T_p_flk/C_TT_flk
        R_Tstar_icesnow = (R_TI_icesnow/2.0-1.0)*R_rho_c_icesnow + 1.0
        d_h_ML_dt = -Q_star_flk*(R_Tstar_icesnow*(1.0+c_cbl_1)-1.0) - Q_bot_flk
        d_h_ML_dt = d_h_ML_dt/tpl_rho_w_r/tpl_c_w                        ! Q_* and Q_b flux terms
        flk_str_2 = (depth_w-h_ML_p_flk)*(T_wML_p_flk-T_bot_p_flk)*C_TT_2/C_TT_flk*d_C_T_dt 
        d_h_ML_dt = d_h_ML_dt + flk_str_2                                 ! Add dC_T/dt term
        flk_str_2 = I_bot_flk + (R_TI_icesnow-1.0)*I_h_flk - R_TI_icesnow*I_intm_h_D_flk
        flk_str_2 = flk_str_2 + (R_TI_icesnow-2.0)*R_rho_c_icesnow*(I_h_flk-I_intm_0_h_flk)
        flk_str_2 = flk_str_2/tpl_rho_w_r/tpl_c_w
        d_h_ML_dt = d_h_ML_dt + flk_str_2                                 ! Add radiation terms
        flk_str_2 = -c_cbl_2*R_Tstar_icesnow*Q_star_flk/tpl_rho_w_r/tpl_c_w/MAX(w_star_sfc_flk, c_small_flk)
        flk_str_2 = flk_str_2 + C_T_p_flk*(T_wML_p_flk-T_bot_p_flk)
        d_h_ML_dt = d_h_ML_dt/flk_str_2                                   ! dh_ML/dt = r.h.s.
      END IF 
!_dm
! Notice that dh_ML/dt may appear to be negative  
! (e.g. due to buoyancy loss to bottom sediments and/or
! the effect of volumetric radiation heating),
! although a negative generalized buoyancy flux scale indicates 
! that the equilibrium CBL depth has not yet been reached
! and convective deepening of the mixed layer should take place.
! Physically, this situation reflects an approximate character of the lake model.
! Using the self-similar temperature profile in the thermocline, 
! there is always communication between the mixed layer, the thermocline 
! and the lake bottom. As a result, the rate of change of the CBL depth
! is always dependent on the bottom heat flux and the radiation heating of the thermocline.
! In reality, convective mixed-layer deepening may be completely decoupled
! from the processes underneath. In order to account for this fact,
! the rate of CBL deepening is set to a small value
! if dh_ML/dt proves to be negative.
! This is "double insurance" however, 
! as a negative dh_ML/dt is encountered very rarely.
!_dm

!_dbg
! IF(d_h_ML_dt.LT.0.0) THEN 
!   print*, 'FLake: negative d_h_ML_dt during convection, = ', d_h_ML_dt
!   print*, '                d_h_ML_dt*del_time = ', MAX(d_h_ML_dt, c_small_flk)*del_time
!   print*, '         u_*       = ', u_star_w_flk   
!   print*, '         w_*       = ', w_star_sfc_flk
!   print*, '         h_CBL_eqi = ', conv_equil_h_scale
!   print*, '         ZM scale  = ', ZM_h_scale
!   print*, '        h_ML_p_flk = ', h_ML_p_flk
! END IF
!   print*, 'FLake: Convection, = ', d_h_ML_dt
!   print*, '         Q_*       = ', Q_star_flk
!   print*, '         \beta*Q_* = ', flk_str_1
!_dbg

      d_h_ML_dt = MAX(d_h_ML_dt, c_small_flk)    
      h_ML_n_flk = h_ML_p_flk + d_h_ML_dt*del_time                       ! Update h_ML 
      h_ML_n_flk = MAX(h_ML_min_flk, MIN(h_ML_n_flk, depth_w))           ! Security, limit h_ML
    ELSE                                              ! Mixing down to the lake bottom
      h_ML_n_flk = depth_w
    END IF

  ELSE Mixing_regime                              ! Wind mixing

    d_h_ML_dt = MAX(u_star_w_flk, u_star_min_flk)                        ! The surface friction velocity
    ZM_h_scale = (ABS(par_Coriolis)/c_sbl_ZM_n + N_T_mean/c_sbl_ZM_i)*d_h_ML_dt**2_iintegers
    ZM_h_scale = ZM_h_scale + flk_str_1/c_sbl_ZM_s
    ZM_h_scale = MAX(ZM_h_scale, c_small_flk)
    ZM_h_scale = d_h_ML_dt**3_iintegers/ZM_h_scale 
    ZM_h_scale = MAX(h_ML_min_flk, MIN(ZM_h_scale, h_ML_max_flk))        ! The ZM96 SBL depth scale 
    ZM_h_scale = MAX(ZM_h_scale, conv_equil_h_scale)                     ! Equilibrium mixed-layer depth 

!_dm 
! In order to avoid numerical discretization problems,
! an analytical solution to the evolution equation 
! for the wind-mixed layer depth is used.
! That is, an exponential relaxation formula is applied
! over the time interval equal to the model time step.
!_dm 

    d_h_ML_dt = c_relax_h*d_h_ML_dt/ZM_h_scale*del_time
    h_ML_n_flk = ZM_h_scale - (ZM_h_scale-h_ML_p_flk)*EXP(-d_h_ML_dt)    ! Update h_ML 
    h_ML_n_flk = MAX(h_ML_min_flk, MIN(h_ML_n_flk, depth_w))             ! Limit h_ML 
    d_h_ML_dt = (h_ML_n_flk-h_ML_p_flk)/del_time                         ! Re-compute dh_ML/dt

    IF(h_ML_n_flk.LE.h_ML_p_flk)           &
      d_C_T_dt = -d_C_T_dt                 ! Mixed-layer retreat or stationary state, dC_T/dt<0
    C_T_n_flk = C_T_p_flk + d_C_T_dt*del_time                            ! Update C_T
    C_T_n_flk = MIN(C_T_max, MAX(C_T_n_flk, C_T_min))                    ! Limit C_T 
    d_C_T_dt = (C_T_n_flk-C_T_p_flk)/del_time                            ! Re-compute dC_T/dt

!_dbg
! print*, 'FLake: wind mixing: d_h_ML_dt*del_time = ', d_h_ML_dt*del_time
! print*, '              h_CBL_eqi = ', conv_equil_h_scale
! print*, '              ZM scale  = ', ZM_h_scale
! print*, '              w_*       = ', w_star_sfc_flk
! print*, '              u_*       = ', u_star_w_flk
! print*, '             h_ML_p_flk = ', h_ML_p_flk
!_dbg

  END IF Mixing_regime

! Compute the time-rate-of-change of the the bottom temperature, 
! depending on the sign of dh_ML/dt 
! Update the bottom temperature and the mixed-layer temperature

  IF(h_ML_n_flk.LE.depth_w-h_ML_min_flk) THEN       ! Mixing did not reach the bottom 

    IF(h_ML_n_flk.GT.h_ML_p_flk) THEN   ! Mixed-layer deepening 
      R_H_icesnow     = h_ML_p_flk/depth_w
      R_rho_c_icesnow = 1.0-R_H_icesnow 
      R_TI_icesnow    = 0.5*C_T_p_flk*R_rho_c_icesnow+C_TT_flk*(2.0*R_H_icesnow-1.0)
      R_Tstar_icesnow = (0.5+C_TT_flk-C_Q_flk)/R_TI_icesnow
      R_TI_icesnow    = (1.0-C_T_p_flk*R_rho_c_icesnow)/R_TI_icesnow
     
      d_T_bot_dt = (Q_w_flk-Q_bot_flk+I_w_flk-I_bot_flk)/tpl_rho_w_r/tpl_c_w
      d_T_bot_dt = d_T_bot_dt - C_T_p_flk*(T_wML_p_flk-T_bot_p_flk)*d_h_ML_dt
      d_T_bot_dt = d_T_bot_dt*R_Tstar_icesnow/depth_w                   ! Q+I fluxes and dh_ML/dt term

      flk_str_2 = I_intm_h_D_flk - (1.0-C_Q_flk)*I_h_flk - C_Q_flk*I_bot_flk
      flk_str_2 = flk_str_2*R_TI_icesnow/(depth_w-h_ML_p_flk)/tpl_rho_w_r/tpl_c_w
      d_T_bot_dt = d_T_bot_dt + flk_str_2                               ! Add radiation-flux term

      flk_str_2 = (1.0-C_TT_2*R_TI_icesnow)/C_T_p_flk
      flk_str_2 = flk_str_2*(T_wML_p_flk-T_bot_p_flk)*d_C_T_dt
      d_T_bot_dt = d_T_bot_dt + flk_str_2                               ! Add dC_T/dt term
      
    ELSE                                ! Mixed-layer retreat or stationary state
      d_T_bot_dt = 0.0                                            ! dT_bot/dt=0
    END IF

    T_bot_n_flk = T_bot_p_flk + d_T_bot_dt*del_time                      ! Update T_bot  
    T_bot_n_flk = MAX(T_bot_n_flk, tpl_T_f)           ! Security, limit T_bot by the freezing point
    flk_str_2 = (T_bot_n_flk-tpl_T_r)*flake_buoypar(T_mnw_n_flk)
    IF(flk_str_2.LT.0.0) T_bot_n_flk = tpl_T_r  ! Security, avoid T_r crossover 
    T_wML_n_flk = C_T_n_flk*(1.0-h_ML_n_flk/depth_w)
    T_wML_n_flk = (T_mnw_n_flk-T_bot_n_flk*T_wML_n_flk)/(1.0-T_wML_n_flk)
    T_wML_n_flk = MAX(T_wML_n_flk, tpl_T_f)           ! Security, limit T_wML by the freezing point
!    print*,'mark 2 T_wML_n_flk=',T_wML_n_flk

  ELSE                                              ! Mixing down to the lake bottom 

    h_ML_n_flk = depth_w
    T_wML_n_flk = T_mnw_n_flk
    T_bot_n_flk = T_mnw_n_flk
    C_T_n_flk = C_T_min
!    print*,'mark 3 T_wML_n_flk=',T_wML_n_flk

  END IF

END IF HTC_Water


!------------------------------------------------------------------------------
!  Compute the depth of the upper layer of bottom sediments
!  and the temperature at that depth.
!------------------------------------------------------------------------------

Use_sediment: IF(lflk_botsed_use) THEN   ! The bottom-sediment scheme is used
  
  IF(H_B1_p_flk.GE.depth_bs-H_B1_min_flk) THEN   ! No T(z) maximum (no thermal wave) 
    H_B1_p_flk = 0.0                       ! Set H_B1_p to zero
    T_B1_p_flk = T_bot_p_flk                     ! Set T_B1_p to the bottom temperature
  END IF 

  flk_str_1 = 2.0*Phi_B1_pr0/(1.0-C_B1)*tpl_kappa_w/tpl_rho_w_r/tpl_c_w*del_time
  h_ice_threshold = SQRT(flk_str_1)                              ! Threshold value of H_B1
  h_ice_threshold = MIN(0.9*depth_bs, h_ice_threshold)    ! Limit H_B1
  flk_str_2 = C_B2/(1.0-C_B2)*(T_bs-T_B1_p_flk)/(depth_bs-H_B1_p_flk)

  IF(H_B1_p_flk.LT.h_ice_threshold) THEN  ! Use a truncated equation for H_B1(t)
    H_B1_n_flk = SQRT(H_B1_p_flk**2_iintegers+flk_str_1)  ! Advance H_B1
    d_H_B1_dt = (H_B1_n_flk-H_B1_p_flk)/del_time          ! Re-compute dH_B1/dt 
  ELSE                                    ! Use a full equation for H_B1(t)
    flk_str_1 = (Q_bot_flk+I_bot_flk)/H_B1_p_flk/tpl_rho_w_r/tpl_c_w
    flk_str_1 = flk_str_1 - (1.0-C_B1)*(T_bot_n_flk-T_bot_p_flk)/del_time
    d_H_B1_dt = (1.0-C_B1)*(T_bot_p_flk-T_B1_p_flk)/H_B1_p_flk + C_B1*flk_str_2
    d_H_B1_dt = flk_str_1/d_H_B1_dt
    H_B1_n_flk = H_B1_p_flk + d_H_B1_dt*del_time          ! Advance H_B1
  END IF 
  d_T_B1_dt = flk_str_2*d_H_B1_dt
  T_B1_n_flk = T_B1_p_flk + d_T_B1_dt*del_time            ! Advance T_B1

!_dbg
! print*, 'BS module: '
! print*, '  Q_bot   = ', Q_bot_flk
! print*, '  d_H_B1_dt = ', d_H_B1_dt
! print*, '  d_T_B1_dt = ', d_T_B1_dt
! print*, '  H_B1    = ', H_B1_n_flk
! print*, '    T_bot = ', T_bot_n_flk
! print*, '  T_B1    = ', T_B1_n_flk
! print*, '    T_bs  = ',  T_bs
!_dbg

!_nu  
! Use a very simplistic procedure, where only the upper layer profile is used, 
! H_B1 is always set to depth_bs, and T_B1 is always set to T_bs.
! Then, the time derivatives are zero, and the sign of the bottom heat flux depends on 
! whether T_bot is smaller or greater than T_bs.
! This is, of course, an oversimplified scheme.
!_nu  d_H_B1_dt = 0.0
!_nu  d_T_B1_dt = 0.0
!_nu  H_B1_n_flk = H_B1_p_flk + d_H_B1_dt*del_time   ! Advance H_B1
!_nu  T_B1_n_flk = T_B1_p_flk + d_T_B1_dt*del_time   ! Advance T_B1
!_nu  

  l_snow_exists = H_B1_n_flk.GE.depth_bs-H_B1_min_flk                    & ! H_B1 reached depth_bs, or
             .OR. H_B1_n_flk.LT.H_B1_min_flk                             & ! H_B1 decreased to zero, or
             .OR.(T_bot_n_flk-T_B1_n_flk)*(T_bs-T_B1_n_flk).LE.0.0   ! there is no T(z) maximum
  IF(l_snow_exists) THEN      
    H_B1_n_flk = depth_bs                     ! Set H_B1 to the depth of the thermally active layer
    T_B1_n_flk = T_bs                         ! Set T_B1 to the climatological temperature 
  END IF

! changed by Shaobo Zhang

!ELSE Use_sediment                        ! The bottom-sediment scheme is not used
!
!  H_B1_n_flk = rflk_depth_bs_ref              ! H_B1 is set to a reference value 
!  T_B1_n_flk = tpl_T_r                        ! T_B1 is set to the temperature of maximum density

ELSE Use_sediment   ! The scheme written by Shaobo Zhang for the deeper layer of a deep lake is used

   if ( abs(T_bot_p_flk-T_bot_2_in) .lt. 0.01 ) then
    depth_bs  = depthbs
    depth_w   = depthw
    T_bot_2_out = T_bot_2_in
   else

   CT = 0.65
   CTT = 11.0/18.0 * CT - 7.0/45.0
   CQ = 2.0 * CTT / CT

! compute d_h_D_dt
    flk_str_1 = ( Q_bot_flk * CQ + I_bot_flk - I_intm_D_H_flk/depth_bs ) /tpl_rho_w_r/tpl_c_w
    flk_str_1 = flk_str_1 - depth_bs * (0.5-CTT) * (T_bot_n_flk-T_bot_p_flk)/del_time
    flk_str_1 = flk_str_1 - CTT/CT*( (Q_bot_flk+I_bot_flk-I_HH_flk)/tpl_rho_w_r/tpl_c_w - &
                depth_bs * ( 1.0 - CT ) * (T_bot_n_flk-T_bot_p_flk)/del_time )
    flk_str_2 = CTT * (T_bot_p_flk-T_bot_2_in)
    if(abs(flk_str_2)<0.01) then
       d_h_D_dt = 0.0
    else
       d_h_D_dt  = flk_str_1/flk_str_2
    endif

! compute d_T_H_dt
    flk_str_1 = (Q_bot_flk+I_bot_flk-I_HH_flk)/tpl_rho_w_r/tpl_c_w
    flk_str_1 = flk_str_1 - depth_bs*(1.0-CT)*(T_bot_n_flk-T_bot_p_flk)/del_time
    flk_str_1 = flk_str_1 - CT * (T_bot_p_flk-T_bot_2_in) * d_h_D_dt
    flk_str_2 = CT * depth_bs
    d_T_H_dt  = flk_str_1/flk_str_2

! update T_bot_2_out
    T_bot_2_out = T_bot_2_in + d_T_H_dt * del_time

    if ( (T_bot_2_out-tpl_T_r)*(T_bot_n_flk-tpl_T_r) .le. 0.0 ) then
      T_bot_2_out = tpl_T_r
    else if ( (T_bot_n_flk-tpl_T_r)/(T_bot_2_out-tpl_T_r) .le. 1.0 ) then
      T_bot_2_out = T_bot_n_flk
    end if

! update depth_w
    flk_str_2 = depth_w + depth_bs                     ! The total depth of the deep lake

    depth_w   = depth_w - d_h_D_dt * del_time          ! Advance depth_bs
!   depth_w   = max (lake_depth_max, min(depth_w, 0.4*flk_str_2))    ! Limit depth_bs
    depth_w   = max (max(lake_depth_max,0.3*flk_str_2), &
                min(depth_w, min(0.4*flk_str_2, 80.0)))    ! Limit depth_bs
    depth_bs  = flk_str_2  - depth_w                  ! Update depth_w
    end if

END IF Use_sediment


!------------------------------------------------------------------------------
!  Impose additional constraints.
!------------------------------------------------------------------------------

! In case of unstable stratification, force mixing down to the bottom
flk_str_2 = (T_wML_n_flk-T_bot_n_flk)*flake_buoypar(T_mnw_n_flk)
IF(flk_str_2.LT.0.0) THEN 

!_dbg
! print*, 'FLake: inverse (unstable) stratification !!! '
! print*, '       Mixing down to the bottom is forced.'
! print*, '  T_wML_p, T_wML_n ', T_wML_p_flk-tpl_T_f, T_wML_n_flk-tpl_T_f
! print*, '  T_mnw_p, T_mnw_n ', T_mnw_p_flk-tpl_T_f, T_mnw_n_flk-tpl_T_f
! print*, '  T_bot_p, T_bot_n ', T_bot_p_flk-tpl_T_f, T_bot_n_flk-tpl_T_f
! print*, '  h_ML_p,  h_ML_n  ', h_ML_p_flk,          h_ML_n_flk
! print*, '  C_T_p,   C_T_n   ', C_T_p_flk,           C_T_n_flk
!_dbg

  h_ML_n_flk = depth_w
  T_wML_n_flk = T_mnw_n_flk
  T_bot_n_flk = T_mnw_n_flk
  C_T_n_flk = C_T_min
!  print*,'mark 4 T_wML_n_flk=',T_wML_n_flk

END IF


!------------------------------------------------------------------------------
!  Update the surface temperature.
!------------------------------------------------------------------------------

IF(h_snow_n_flk.GE.h_Snow_min_flk) THEN   
  T_sfc_n = T_snow_n_flk                   ! Snow exists, use the snow temperature
ELSE IF(h_ice_n_flk.GE.h_Ice_min_flk) THEN
  T_sfc_n = T_ice_n_flk                    ! Ice exists but there is no snow, use the ice temperature
ELSE 
  T_sfc_n = T_wML_n_flk                    ! No ice-snow cover, use the mixed-layer temperature
END IF
if(T_sfc_n.lt.200.0 .or. T_sfc_n.gt.350.0) then
   print*,'h_snow_n_flk=',h_snow_n_flk,' h_Snow_min_flk=',h_Snow_min_flk
   print*,'h_ice_n_flk=',h_ice_n_flk, ' h_Ice_min_flk= ',h_Ice_min_flk
   print*,'T_wML_n_flk=',T_wML_n_flk
   print*,'T_sfc_n=',T_sfc_n
endif
!------------------------------------------------------------------------------
!  End calculations
!==============================================================================

END SUBROUTINE flake_main

!==============================================================================

!==============================================================================
! include 'flake_buoypar.incf'
!------------------------------------------------------------------------------

REAL (KIND = ireals) FUNCTION flake_buoypar (T_water)

!------------------------------------------------------------------------------
!
! Description:
!
!  Computes the buoyancy parameter,
!  using a quadratic equation of state for the fresh-water.
!  
! Declarations:
!
! Modules used:

!_dm Parameters are USEd in module "flake".
!_nu USE data_parameters , ONLY : &
!_nu     ireals,                  & ! KIND-type parameter for real variables
!_nu     iintegers                  ! KIND-type parameter for "normal" integer variables

USE flake_parameters , ONLY : &
  tpl_grav                  , & ! Acceleration due to gravity [m s^{-2}]
  tpl_T_r                   , & ! Temperature of maximum density of fresh water [K]
  tpl_a_T                       ! Constant in the fresh-water equation of state [K^{-2}]

use machine,               only: kind_phys
!==============================================================================

IMPLICIT NONE

!==============================================================================
!
! Declarations
 
!  Input (function argument) 
REAL (KIND = kind_phys), INTENT(IN) :: &
  T_water                             ! Water temperature [K]

!==============================================================================
!  Start calculations
!------------------------------------------------------------------------------

! Buoyancy parameter [m s^{-2} K^{-1}]

  flake_buoypar = tpl_grav*tpl_a_T*(T_water-tpl_T_r)

!------------------------------------------------------------------------------
!  End calculations
!==============================================================================

END FUNCTION flake_buoypar

!==============================================================================

!==============================================================================
! include 'flake_snowdensity.incf'
!------------------------------------------------------------------------------

REAL (KIND = ireals) FUNCTION flake_snowdensity (h_snow)

!------------------------------------------------------------------------------
!
! Description:
!
!  Computes the snow density,
!  using an empirical approximation from Heise et al. (2003).
!  
! Declarations:
!
! Modules used:

!_dm Parameters are USEd in module "flake".
!_nu USE data_parameters , ONLY : &
!_nu     ireals,                  & ! KIND-type parameter for real variables
!_nu     iintegers                  ! KIND-type parameter for "normal" integer variables

USE flake_parameters , ONLY : &
  tpl_rho_w_r               , & ! Maximum density of fresh water [kg m^{-3}]
  tpl_rho_S_min             , & ! Minimum snow density [kg m^{-3}]
  tpl_rho_S_max             , & ! Maximum snow density [kg m^{-3}]
  tpl_Gamma_rho_S           , & ! Empirical parameter [kg m^{-4}] in the expression for the snow density
  c_small_flk                   ! A small number

use machine,               only: kind_phys
!==============================================================================

IMPLICIT NONE

!==============================================================================
!
! Declarations
 
!  Input (function argument) 
REAL (KIND = kind_phys), INTENT(IN) :: &
  h_snow                              ! Snow thickness [m]

!==============================================================================
!  Start calculations
!------------------------------------------------------------------------------

! Snow density [kg m^{-3}]

!  Security. Ensure that the expression in () does not become negative at a very large h_snow.
  flake_snowdensity = MAX( c_small_flk, (1.0 - h_snow*tpl_Gamma_rho_S/tpl_rho_w_r) )
  flake_snowdensity = MIN( tpl_rho_S_max, tpl_rho_S_min/flake_snowdensity )

!------------------------------------------------------------------------------
!  End calculations
!==============================================================================

END FUNCTION flake_snowdensity 

!==============================================================================

!==============================================================================
! include 'flake_snowheatconduct.incf'
!------------------------------------------------------------------------------

REAL (KIND = ireals) FUNCTION flake_snowheatconduct (h_snow)

!------------------------------------------------------------------------------
!
! Description:
!
!  Computes the snow heat conductivity,
!  using an empirical approximation from Heise et al. (2003).
!  
! Declarations:
!
! Modules used:

!_dm Parameters are USEd in module "flake".
!_nu USE data_parameters , ONLY : &
!_nu     ireals,                  & ! KIND-type parameter for real variables
!_nu     iintegers                  ! KIND-type parameter for "normal" integer variables

USE flake_parameters , ONLY : &
  tpl_rho_w_r               , & ! Maximum density of fresh water [kg m^{-3}]
  tpl_kappa_S_min           , & ! Minimum molecular heat conductivity of snow [J m^{-1} s^{-1} K^{-1}]
  tpl_kappa_S_max           , & ! Maximum molecular heat conductivity of snow [J m^{-1} s^{-1} K^{-1}]
  tpl_Gamma_kappa_S             ! Empirical parameter [J m^{-2} s^{-1} K^{-1}] in the expression for kappa_S

use machine,               only: kind_phys
!==============================================================================

IMPLICIT NONE

!==============================================================================
!
! Declarations
 
!  Input (function argument) 
REAL (KIND = kind_phys), INTENT(IN) :: &
  h_snow                              ! Snow thickness [m]

!==============================================================================
!  Start calculations
!------------------------------------------------------------------------------

! Snow heat conductivity [J m^{-1} s^{-1} K^{-1} = kg m s^{-3} K^{-1}]

  flake_snowheatconduct = flake_snowdensity( h_snow )   ! Compute snow density
  flake_snowheatconduct = MIN( tpl_kappa_S_max, tpl_kappa_S_min                      &
                        + h_snow*tpl_Gamma_kappa_S*flake_snowheatconduct/tpl_rho_w_r )

!------------------------------------------------------------------------------
!  End calculations
!==============================================================================

END FUNCTION flake_snowheatconduct

!==============================================================================

END MODULE flake

!------------------------------------------------------------------------------

MODULE SfcFlx

!------------------------------------------------------------------------------
!
! Description:
!
!  The main program unit of 
!  the atmospheric surface-layer parameterization scheme "SfcFlx".
!  "SfcFlx" is used to compute fluxes 
!  of momentum and of sensible and latent heat over lakes.
!  The surface-layer scheme developed by Mironov (1991) was used as the starting point.
!  It was modified and further developed to incorporate new results as to 
!  the roughness lenghts for scalar quantities,
!  heat and mass transfer in free convection,
!  and the effect of limited fetch on the momentum transfer.
!  Apart from the momentum flux and sensible and latent heat fluxes,
!  the long-wave radiation flux from the water surface and
!  the long-wave radiation flux from the atmosphere can also be computed.
!  The atmospheric long-wave radiation flux is computed with simple empirical formulae,
!  where the atmospheric emissivity is taken to be dependent on 
!  the water vapour pressure and cloud fraction.
!
!  A description of SfcFlx is available from the author.
!  Dmitrii Mironov 
!  German Weather Service, Kaiserleistr. 29/35, D-63067 Offenbach am Main, Germany. 
!  dmitrii.mironov@dwd.de 
!
!  Lines embraced with "!_tmp" contain temporary parts of the code.
!  Lines embraced/marked with "!_dev" may be replaced
!  as improved parameterizations are developed and tested.
!  Lines embraced/marked with "!_dm" are DM's comments
!  that may be helpful to a user.
!  Lines embraced/marked with "!_dbg" are used
!  for debugging purposes only.
!


USE data_parameters  , ONLY :   &
  ireals                      , & ! KIND-type parameter for real variables
  iintegers                       ! KIND-type parameter for "normal" integer variables

USE flake_parameters , ONLY :   &
  tpl_grav                    , & ! Acceleration due to gravity [m s^{-2}]
  tpl_T_f                     , & ! Fresh water freezing point [K]
  tpl_rho_w_r                 , & ! Maximum density of fresh water [kg m^{-3}]
  tpl_c_w                     , & ! Specific heat of water [J kg^{-1} K^{-1}]
  tpl_L_f                     , & ! Latent heat of fusion [J kg^{-1}]
  h_Ice_min_flk                   ! Minimum ice thickness [m]

use machine,               only: kind_phys

!==============================================================================

IMPLICIT NONE

!==============================================================================
!
! Declarations

!  Dimensionless constants in the Monin-Obukhov surface-layer 
!  similarity relations and in the expressions for the roughness lengths.
REAL (KIND = kind_phys), PARAMETER ::   &
  c_Karman      = 0.40      , & ! The von Karman constant 
!  Pr_neutral    = 1.0       , & ! Turbulent Prandtl number at neutral static stability
  Pr_neutral    = 0.9       , & ! Turbulent Prandtl number at neutral static stability
  Sc_neutral    = 1.0       , & ! Turbulent Schmidt number at neutral static stability
  c_MO_u_stab   = 5.0       , & ! Constant of the MO theory (wind, stable stratification)
  c_MO_t_stab   = 5.0       , & ! Constant of the MO theory (temperature, stable stratification)
  c_MO_q_stab   = 5.0       , & ! Constant of the MO theory (humidity, stable stratification)
  c_MO_u_conv   = 15.0      , & ! Constant of the MO theory (wind, convection)
  c_MO_t_conv   = 15.0      , & ! Constant of the MO theory (temperature, convection)
  c_MO_q_conv   = 15.0      , & ! Constant of the MO theory (humidity, convection)
  c_MO_u_exp    = 0.25      , & ! Constant of the MO theory (wind, exponent)
  c_MO_t_exp    = 0.5       , & ! Constant of the MO theory (temperature, exponent)
  c_MO_q_exp    = 0.5       , & ! Constant of the MO theory (humidity, exponent)
  z0u_ice_rough = 1.0E-03   , & ! Aerodynamic roughness of the ice surface [m] (rough flow)
  c_z0u_smooth  = 0.1       , & ! Constant in the expression for z0u (smooth flow) 
  c_z0u_rough   = 1.23E-02  , & ! The Charnock constant in the expression for z0u (rough flow)
  c_z0u_rough_L = 1.00E-01  , & ! An increased Charnock constant (used as the upper limit)
  c_z0u_ftch_f  = 0.70      , & ! Factor in the expression for fetch-dependent Charnock parameter
  c_z0u_ftch_ex = 0.3333333 , & ! Exponent in the expression for fetch-dependent Charnock parameter
  c_z0t_rough_1 = 4.0       , & ! Constant in the expression for z0t (factor) 
  c_z0t_rough_2 = 3.2       , & ! Constant in the expression for z0t (factor)
  c_z0t_rough_3 = 0.5       , & ! Constant in the expression for z0t (exponent) 
  c_z0q_rough_1 = 4.0       , & ! Constant in the expression for z0q (factor)
  c_z0q_rough_2 = 4.2       , & ! Constant in the expression for z0q (factor)
  c_z0q_rough_3 = 0.5       , & ! Constant in the expression for z0q (exponent)
  c_z0t_ice_b0s = 1.250     , & ! Constant in the expression for z0t over ice
  c_z0t_ice_b0t = 0.149     , & ! Constant in the expression for z0t over ice
  c_z0t_ice_b1t = -0.550    , & ! Constant in the expression for z0t over ice
  c_z0t_ice_b0r = 0.317     , & ! Constant in the expression for z0t over ice
  c_z0t_ice_b1r = -0.565    , & ! Constant in the expression for z0t over ice
  c_z0t_ice_b2r = -0.183    , & ! Constant in the expression for z0t over ice
  c_z0q_ice_b0s = 1.610     , & ! Constant in the expression for z0q over ice
  c_z0q_ice_b0t = 0.351     , & ! Constant in the expression for z0q over ice
  c_z0q_ice_b1t = -0.628    , & ! Constant in the expression for z0q over ice
  c_z0q_ice_b0r = 0.396     , & ! Constant in the expression for z0q over ice
  c_z0q_ice_b1r = -0.512    , & ! Constant in the expression for z0q over ice
  c_z0q_ice_b2r = -0.180    , & ! Constant in the expression for z0q over ice
  Re_z0s_ice_t  = 2.5       , & ! Threshold value of the surface Reynolds number 
                                ! used to compute z0t and z0q over ice (Andreas 2002)
  Re_z0u_thresh = 0.1           ! Threshold value of the roughness Reynolds number 
                                ! [value from Zilitinkevich, Grachev, and Fairall (200),
                                ! currently not used] 

!  Dimensionless constants 
REAL (KIND = kind_phys), PARAMETER ::   &
  c_free_conv   = 0.14          ! Constant in the expressions for fluxes in free convection

!  Dimensionless constants 
REAL (KIND = kind_phys), PARAMETER ::   &
  c_lwrad_emis  = 0.99          ! Surface emissivity with respect to the long-wave radiation

!  Thermodynamic parameters
REAL (KIND = kind_phys), PARAMETER ::        &
  tpsf_C_StefBoltz = 5.67E-08    , & ! The Stefan-Boltzmann constant [W m^{-2} K^{-4}]
  tpsf_R_dryair    = 2.8705E+02  , & ! Gas constant for dry air [J kg^{-1} K^{-1}]
  tpsf_R_watvap    = 4.6151E+02  , & ! Gas constant for water vapour [J kg^{-1} K^{-1}]
  tpsf_c_a_p       = 1.005E+03   , & ! Specific heat of air at constant pressure [J kg^{-1} K^{-1}]
  tpsf_L_evap      = 2.501E+06   , & ! Specific heat of evaporation [J kg^{-1}]
  tpsf_nu_u_a      = 1.50E-05    , & ! Kinematic molecular viscosity of air [m^{2} s^{-1}]
  tpsf_kappa_t_a   = 2.20E-05    , & ! Molecular temperature conductivity of air [m^{2} s^{-1}]
  tpsf_kappa_q_a   = 2.40E-05        ! Molecular diffusivity of air for water vapour [m^{2} s^{-1}]

!  Derived thermodynamic parameters
REAL (KIND = kind_phys), PARAMETER ::                        &
  tpsf_Rd_o_Rv  = tpsf_R_dryair/tpsf_R_watvap           , & ! Ratio of gas constants (Rd/Rv)
  tpsf_alpha_q  = (1.0-tpsf_Rd_o_Rv)/tpsf_Rd_o_Rv     ! Diemsnionless ratio 

!  Thermodynamic parameters
REAL (KIND = kind_phys), PARAMETER ::     &
  P_a_ref             = 1.0E+05   ! Reference pressure [N m^{-2} = kg m^{-1} s^{-2}]


!  The variables declared below
!  are accessible to all program units of the MODULE "SfcFlx"
!  and to the driving routines that use "SfcFlx".
!  These are basically the quantities computed by SfcFlx.
!  Apart from these quantities, there a few local scalars 
!  used by SfcFlx routines mainly for security reasons.
!  All variables declared below have a suffix "sf".

!  SfcFlx variables of type REAL

!  Roughness lengths
REAL (KIND = kind_phys) ::    &
  z0u_sf                 , & ! Roughness length with respect to wind velocity [m]
  z0t_sf                 , & ! Roughness length with respect to potential temperature [m]
  z0q_sf                     ! Roughness length with respect to specific humidity [m]

!  Fluxes in the surface air layer
REAL (KIND = kind_phys) ::    &
  u_star_a_sf            , & ! Friction velocity [m s^{-1}]
  Q_mom_a_sf             , & ! Momentum flux [N m^{-2}]
  Q_sens_a_sf            , & ! Sensible heat flux [W m^{-2}]
  Q_lat_a_sf             , & ! Laten heat flux [W m^{-2}]
  Q_watvap_a_sf              ! Flux of water vapout [kg m^{-2} s^{-1}]

!  Security constants
REAL (KIND = kind_phys), PARAMETER ::   &
  u_wind_min_sf  = 1.0E-02  , & ! Minimum wind speed [m s^{-1}]
  u_star_min_sf  = 1.0E-04  , & ! Minimum value of friction velocity [m s^{-1}]
  c_accur_sf     = 1.0E-07  , & ! A small number (accuracy)
  c_small_sf     = 1.0E-04      ! A small number (used to compute fluxes)

!  Useful constants
REAL (KIND = kind_phys), PARAMETER ::     &
  num_1o3_sf = 1.0/3.0       ! 1/3

!==============================================================================
! Procedures 
!==============================================================================

CONTAINS

!==============================================================================
!  The codes of the SfcFlx procedures are stored in separate "*.incf" files
!  and are included below.
!------------------------------------------------------------------------------

!==============================================================================
! include 'SfcFlx_lwradatm.incf'
!------------------------------------------------------------------------------

REAL (KIND = kind_phys) FUNCTION SfcFlx_lwradatm (T_a, e_a, cl_tot, cl_low)

!------------------------------------------------------------------------------
!
! Description:
!
!  Computes the long-wave radiation flux from the atmosphere
!  as function of air temperature, water vapour pressure and cloud fraction. 
!
! Declarations:
!
! Modules used:

!_dm Parameters are USEd in module "SfcFlx".
!_nu USE data_parameters , ONLY : &
!_nu     ireals,                  & ! KIND-type parameter for real variables
!_nu     iintegers                  ! KIND-type parameter for "normal" integer variables

!==============================================================================

IMPLICIT NONE

!==============================================================================
!
! Declarations
 
!  Input (function argument) 
REAL (KIND = kind_phys), INTENT(IN) ::   &
  T_a                               , & ! Air temperature [K]
  e_a                               , & ! Water vapour pressure [N m^{-2} = kg m^{-1} s^{-2}]
  cl_tot                            , & ! Total cloud cover [0,1]
  cl_low                                ! Lowe-level cloud cover [0,1]


!  Local parameters  

!  Coefficients in the empirical formulation  
!  developed at the Main Geophysical Observatory (MGO), St. Petersburg, Russia.
REAL (KIND = kind_phys), PARAMETER ::   &
  c_lmMGO_1    = 43.057924  , & ! Empirical coefficient 
  c_lmMGO_2    = 540.795        ! Empirical coefficient 
!  Temperature-dependent cloud-correction coefficients in the MGO formula
INTEGER (KIND = iintegers), PARAMETER :: &
  nband_coef = 6_iintegers                 ! Number of temperature bands
REAL (KIND = kind_phys), PARAMETER, DIMENSION (nband_coef) ::      &
  corr_cl_tot     = (/0.70, 0.45, 0.32,    & 
                      0.23, 0.18, 0.13/) , & ! Total clouds
  corr_cl_low     = (/0.76, 0.49, 0.35,    & 
                      0.26, 0.20, 0.15/) , & ! Low-level clouds
  corr_cl_midhigh = (/0.46, 0.30, 0.21,    & 
                      0.15, 0.12, 0.09/)     ! Mid- and high-level clouds
REAL (KIND = kind_phys), PARAMETER ::   &
  T_low  = 253.15           , & ! Low-limit temperature in the interpolation formula [K]
  del_T  = 10.0                 ! Temperature step in the interpolation formula [K]

!  Coefficients in the empirical water-vapour correction function 
!  (see Fung et al. 1984, Zapadka and Wozniak 2000, Zapadka et al. 2001). 
REAL (KIND = kind_phys), PARAMETER ::     &
  c_watvap_corr_min = 0.6100  , & ! Empirical coefficient (minimum value of the correction function)
  c_watvap_corr_max = 0.7320  , & ! Empirical coefficient (maximum value of the correction function)
  c_watvap_corr_e   = 0.0050      ! Empirical coefficient [(N m^{-2})^{-1/2}]

!  Local variables of type INTEGER
INTEGER (KIND = iintegers) :: &
  i                             ! Loop index

!  Local variables of type REAL
REAL (KIND = kind_phys) ::    &
  c_cl_tot_corr          , & ! The MGO cloud correction coefficient, total clouds
  c_cl_low_corr          , & ! The MGO cloud correction coefficient, low-level clouds
  c_cl_midhigh_corr      , & ! The MGO cloud correction coefficient, mid- and high-level clouds
  T_corr                 , & ! Temperature used to compute the MGO cloud correction [K]
  f_wvpres_corr          , & ! Correction function with respect to water vapour
  f_cloud_corr               ! Correction function with respect to cloudiness
 
!==============================================================================
!  Start calculations
!------------------------------------------------------------------------------

! Water-vapour correction function
  f_wvpres_corr = c_watvap_corr_min + c_watvap_corr_e*SQRT(e_a) 
  f_wvpres_corr = MIN(f_wvpres_corr, c_watvap_corr_max)

! Cloud-correction coefficients using the MGO formulation with linear interpolation 
IF(T_a.LT.T_low) THEN 
  c_cl_tot_corr     = corr_cl_tot(1)   
  c_cl_low_corr     = corr_cl_low(1)
  c_cl_midhigh_corr = corr_cl_midhigh(1)
ELSE IF(T_a.GE.T_low+(nband_coef-1_iintegers)*del_T) THEN
  c_cl_tot_corr     = corr_cl_tot(nband_coef)   
  c_cl_low_corr     = corr_cl_low(nband_coef)
  c_cl_midhigh_corr = corr_cl_midhigh(nband_coef)
ELSE 
  T_corr = T_low
  DO i=1, nband_coef-1
    IF(T_a.GE.T_corr.AND.T_a.LT.T_corr+del_T) THEN 
      c_cl_tot_corr = (T_a-T_corr)/del_T
      c_cl_low_corr = corr_cl_low(i) + (corr_cl_low(i+1)-corr_cl_low(i))*c_cl_tot_corr
      c_cl_midhigh_corr = corr_cl_midhigh(i) + (corr_cl_midhigh(i+1)-corr_cl_midhigh(i))*c_cl_tot_corr
      c_cl_tot_corr = corr_cl_tot(i) + (corr_cl_tot(i+1)-corr_cl_tot(i))*c_cl_tot_corr
    END IF 
    T_corr = T_corr + del_T
  END DO
END IF
! Cloud correction function
IF(cl_low.LT.0.0) THEN  ! Total cloud cover only 
  f_cloud_corr = 1.0 + c_cl_tot_corr*cl_tot*cl_tot
ELSE                          ! Total and low-level cloud cover
  f_cloud_corr = (1.0 + c_cl_low_corr*cl_low*cl_low)  &
               * (1.0 + c_cl_midhigh_corr*(cl_tot*cl_tot-cl_low*cl_low))
END IF

! Long-wave radiation flux [W m^{-2}]

!  The MGO formulation  
!_nu The MGO formulation  
!_nu SfcFlx_lwradatm = -SfcFlx_lwradatm*c_lwrad_emis  &
!_nu                 * (c_lmMGO_1*SQRT(tpsf_C_StefBoltz*T_a**4_iintegers)-c_lmMGO_2)
!_nu 

!  "Conventional" formulation  
!  (see Fung et al. 1984, Zapadka and Wozniak 2000, Zapadka et al. 2001)  
SfcFlx_lwradatm = -c_lwrad_emis*tpsf_C_StefBoltz*T_a**4_iintegers  &
                * f_wvpres_corr*f_cloud_corr

!------------------------------------------------------------------------------
!  End calculations
!==============================================================================

END FUNCTION SfcFlx_lwradatm

!==============================================================================

!==============================================================================
! include 'SfcFlx_lwradwsfc.incf'
!------------------------------------------------------------------------------

REAL (KIND = kind_phys) FUNCTION SfcFlx_lwradwsfc (T)

!------------------------------------------------------------------------------
!
! Description:
!
!  Computes the surface long-wave radiation flux
!  as function of temperature. 
!  
! Declarations:
!
! Modules used:

!_dm Parameters are USEd in module "SfcFlx".
!_nu USE data_parameters , ONLY : &
!_nu     ireals,                  & ! KIND-type parameter for real variables
!_nu     iintegers                  ! KIND-type parameter for "normal" integer variables

!==============================================================================

IMPLICIT NONE

!==============================================================================
!
! Declarations
 
!  Input (function argument) 
REAL (KIND = kind_phys), INTENT(IN) ::   &
  T                                     ! Temperature [K]

!==============================================================================
!  Start calculations
!------------------------------------------------------------------------------

! Long-wave radiation flux [W m^{-2}]

SfcFlx_lwradwsfc = c_lwrad_emis*tpsf_C_StefBoltz*T**4_iintegers

!------------------------------------------------------------------------------
!  End calculations
!==============================================================================

END FUNCTION SfcFlx_lwradwsfc

!==============================================================================

!==============================================================================
! include 'SfcFlx_momsenlat.incf'
!------------------------------------------------------------------------------

SUBROUTINE SfcFlx_momsenlat ( height_u, height_tq, fetch,                &
                              U_a, T_a, q_a, T_s, P_a, h_ice,            &
                              Q_momentum, Q_sensible, Q_latent, Q_watvap,&
                              q_s,rho_a ) 

!------------------------------------------------------------------------------
!
! Description:
!
!  The SfcFlx routine 
!  where fluxes of momentum and of sensible and latent heat 
!  at the air-water or air-ice (air-snow) interface are computed. 
!
!  Lines embraced with "!_tmp" contain temporary parts of the code.
!  Lines embraced/marked with "!_dev" may be replaced
!  as improved parameterizations are developed and tested.
!  Lines embraced/marked with "!_dm" are DM's comments
!  that may be helpful to a user.
!  Lines embraced/marked with "!_dbg" are used 
!  for debugging purposes only.
!
! Declarations:
!
! Modules used:

!_dm Parameters are USEd in module "SfcFlx".
!_nu USE data_parameters , ONLY : &
!_nu     ireals,                  & ! KIND-type parameter for real variables
!_nu     iintegers                  ! KIND-type parameter for "normal" integer variables

use machine,               only: kind_phys
!==============================================================================

IMPLICIT NONE

!==============================================================================
!
! Declarations

!  Input (procedure arguments)

REAL (KIND = kind_phys), INTENT(IN) ::   &
  height_u                          , & ! Height where wind is measured [m]
  height_tq                         , & ! Height where temperature and humidity are measured [m]
  fetch                             , & ! Typical wind fetch [m]
  U_a                               , & ! Wind speed [m s^{-1}]
  T_a                               , & ! Air temperature [K]
  q_a                               , & ! Air specific humidity [-]
  T_s                               , & ! Surface temperature (water, ice or snow) [K]
  P_a                               , & ! Surface air pressure [N m^{-2} = kg m^{-1} s^{-2}]
  h_ice                                 ! Ice thickness [m]

!  Output (procedure arguments)

REAL (KIND = kind_phys), INTENT(OUT) ::   &
  Q_momentum                         , & ! Momentum flux [N m^{-2}]  
  Q_sensible                         , & ! Sensible heat flux [W m^{-2}]  
  Q_latent                           , & ! Laten heat flux [W m^{-2}]
  Q_watvap                               ! Flux of water vapout [kg m^{-2} s^{-1}]


!  Local parameters of type INTEGER
INTEGER (KIND = iintegers), PARAMETER ::  &
  n_iter_max     = 24                       ! Maximum number of iterations 

!  Local variables of type LOGICAL
LOGICAL ::          &
  l_conv_visc     , & ! Switch, TRUE = viscous free convection, the Nu=C Ra^(1/3) law is used
  l_conv_cbl          ! Switch, TRUE = CBL scale convective structures define surface fluxes 

!  Local variables of type INTEGER
INTEGER (KIND = iintegers) ::   &
  i                           , & ! Loop index
  n_iter                          ! Number of iterations performed 

!  Local variables of type REAL
REAL (KIND = kind_phys) ::    &
  rho_a                  , & ! Air density [kg m^{-3}]  
  wvpres_s               , & ! Saturation water vapour pressure at T=T_s [N m^{-2}]
  q_s                        ! Saturation specific humidity at T=T_s [-]

!  Local variables of type REAL
REAL (KIND = kind_phys) ::    &
  Q_mom_tur              , & ! Turbulent momentum flux [N m^{-2}]
  Q_sen_tur              , & ! Turbulent sensible heat flux [W m^{-2}]  
  Q_lat_tur              , & ! Turbulent laten heat flux [W m^{-2}]
  Q_mom_mol              , & ! Molecular momentum flux [N m^{-2}]
  Q_sen_mol              , & ! Molecular sensible heat flux [W m^{-2}]  
  Q_lat_mol              , & ! Molecular laten heat flux [W m^{-2}]
  Q_mom_con              , & ! Momentum flux in free convection [N m^{-2}]
  Q_sen_con              , & ! Sensible heat flux in free convection [W m^{-2}]  
  Q_lat_con                  ! Laten heat flux in free convection [W m^{-2}]

!  Local variables of type REAL
REAL (KIND = kind_phys) ::    &
  par_conv_visc          , & ! Viscous convection stability parameter
  par_conv_cbl           , & ! CBL convection stability parameter
  c_z0u_fetch            , & ! Fetch-dependent Charnock parameter
  U_a_thresh             , & ! Threshld value of the wind speed [m s^{-1}] 
  u_star_thresh          , & ! Threshld value of friction velocity [m s^{-1}]
  u_star_previter        , & ! Friction velocity from previous iteration [m s^{-1}]
  u_star_n               , & ! Friction velocity at neutral stratification [m s^{-1}]
  u_star_st              , & ! Friction velocity with due regard for stratification [m s^{-1}]
  ZoL                    , & ! The z/L ratio, z=height_u
  Ri                     , & ! Gradient Richardson number 
  Ri_cr                  , & ! Critical value of Ri 
  R_z                    , & ! Ratio of "height_tq" to "height_u"
  Fun                    , & ! A function of generic variable "x"
  Fun_prime              , & ! Derivative of "Fun" with respect to "x"
  Delta                  , & ! Relative error 
  psi_u                  , & ! The MO stability function for wind profile
  psi_t                  , & ! The MO stability function for temperature profile
  psi_q                      ! The MO stability function for specific humidity profile


!==============================================================================
!  Start calculations
!------------------------------------------------------------------------------

!write(*,*) 'Inside Flake Flux Module'
!write(*,*)   height_u,height_tq,fetch,U_a,T_a,q_a,T_s,P_a,h_ice

!_dm All fluxes are positive when directed upwards.

!------------------------------------------------------------------------------
!  Compute saturation specific humidity and the air density at T=T_s
!------------------------------------------------------------------------------

wvpres_s = SfcFlx_satwvpres(T_s, h_ice)  ! Saturation water vapour pressure at T=T_s
q_s = SfcFlx_spechum (wvpres_s, P_a)     ! Saturation specific humidity at T=T_s
rho_a = SfcFlx_rhoair(T_s, q_s, P_a)     ! Air density at T_s and q_s (surface values)

!------------------------------------------------------------------------------
!  Compute molecular fluxes of momentum and of sensible and latent heat
!------------------------------------------------------------------------------

!_dm The fluxes are in kinematic units
Q_mom_mol = -tpsf_nu_u_a*U_a/height_u 
Q_sen_mol = -tpsf_kappa_t_a*(T_a-T_s)/height_tq    
Q_lat_mol = -tpsf_kappa_q_a*(q_a-q_s)/height_tq  

!------------------------------------------------------------------------------
!  Compute fluxes in free convection
!------------------------------------------------------------------------------

par_conv_visc = (T_s-T_a)/T_s*SQRT(tpsf_kappa_t_a) + (q_s-q_a)*tpsf_alpha_q*SQRT(tpsf_kappa_q_a)
IF(par_conv_visc.GT.0.0) THEN   ! Viscous convection takes place
  l_conv_visc = .TRUE.
  par_conv_visc = (par_conv_visc*tpl_grav/tpsf_nu_u_a)**num_1o3_sf
  Q_sen_con = c_free_conv*SQRT(tpsf_kappa_t_a)*par_conv_visc  
  Q_sen_con = Q_sen_con*(T_s-T_a)
  Q_lat_con = c_free_conv*SQRT(tpsf_kappa_q_a)*par_conv_visc
  Q_lat_con = Q_lat_con*(q_s-q_a)
ELSE                                  ! No viscous convection, set fluxes to zero
  l_conv_visc = .FALSE.
  Q_sen_con = 0.0 
  Q_lat_con = 0.0
END IF
Q_mom_con = 0.0                 ! Momentum flux in free (viscous or CBL-scale) convection is zero  

!------------------------------------------------------------------------------
!  Compute turbulent fluxes
!------------------------------------------------------------------------------

R_z   = height_tq/height_u                        ! Ratio of "height_tq" to "height_u"
Ri_cr = c_MO_t_stab/c_MO_u_stab**2_iintegers*R_z  ! Critical Ri
Ri    = tpl_grav*((T_a-T_s)/T_s+tpsf_alpha_q*(q_a-q_s))/MAX(U_a,u_wind_min_sf)**2_iintegers
Ri    = Ri*height_u/Pr_neutral                    ! Gradient Richardson number

Turb_Fluxes: IF(U_a.LT.u_wind_min_sf.OR.Ri.GT.Ri_cr-c_small_sf) THEN  ! Low wind or Ri>Ri_cr 

u_star_st = 0.0                       ! Set turbulent fluxes to zero 
Q_mom_tur = 0.0                       
Q_sen_tur = 0.0   
Q_lat_tur = 0.0  

ELSE Turb_Fluxes                            ! Compute turbulent fluxes using MO similarity

! Compute z/L, where z=height_u
IF(Ri.GE.0.0) THEN   ! Stable stratification
  ZoL = SQRT(1.0-4.0*(c_MO_u_stab-R_z*c_MO_t_stab)*Ri)
  ZoL = ZoL - 1.0 + 2.0*c_MO_u_stab*Ri
  ZoL = ZoL/2.0/c_MO_u_stab/c_MO_u_stab/(Ri_cr-Ri)
ELSE                       ! Convection
  n_iter = 0_iintegers
  Delta = 1.0                ! Set initial error to a large value (as compared to the accuracy)
  u_star_previter = Ri*MAX(1.0, SQRT(R_z*c_MO_t_conv/c_MO_u_conv)) ! Initial guess for ZoL
  DO WHILE (Delta.GT.c_accur_sf.AND.n_iter.LT.n_iter_max) 
    Fun = u_star_previter**2_iintegers*(c_MO_u_conv*u_star_previter-1.0)  &
        + Ri**2_iintegers*(1.0-R_z*c_MO_t_conv*u_star_previter)
    Fun_prime = 3.0*c_MO_u_conv*u_star_previter**2_iintegers              &
              - 2.0*u_star_previter - R_z*c_MO_t_conv*Ri**2_iintegers
    ZoL = u_star_previter - Fun/Fun_prime
    Delta = ABS(ZoL-u_star_previter)/MAX(c_accur_sf, ABS(ZoL+u_star_previter))
    u_star_previter = ZoL
    n_iter = n_iter + 1_iintegers
  END DO 
!_dbg
!  IF(n_iter.GE.n_iter_max-1_iintegers)  & 
!    print*(*,*) 'ZoL: Max No. iters. exceeded (n_iter = ', n_iter, ')!'
!_dbg
END IF

!  Compute fetch-dependent Charnock parameter, use "u_star_min_sf"
CALL SfcFlx_roughness (fetch, U_a, u_star_min_sf, h_ice, c_z0u_fetch, u_star_thresh, z0u_sf, z0t_sf, z0q_sf)

!  Threshold value of wind speed 
u_star_st = u_star_thresh
CALL SfcFlx_roughness (fetch, U_a, u_star_st, h_ice, c_z0u_fetch, u_star_thresh, z0u_sf, z0t_sf, z0q_sf)
IF(ZoL.GT.0.0) THEN   ! MO function in stable stratification 
  psi_u = c_MO_u_stab*ZoL*(1.0-MIN(z0u_sf/height_u, 1.0))
ELSE                        ! MO function in convection
  psi_t = (1.0-c_MO_u_conv*ZoL)**c_MO_u_exp
  psi_q = (1.0-c_MO_u_conv*ZoL*MIN(z0u_sf/height_u, 1.0))**c_MO_u_exp
  psi_u = 2.0*(ATAN(psi_t)-ATAN(psi_q))                  &
        + 2.0*LOG((1.0+psi_q)/(1.0+psi_t))   &
        + LOG((1.0+psi_q*psi_q)/(1.0+psi_t*psi_t))   
END IF 
U_a_thresh = u_star_thresh/c_Karman*(LOG(height_u/z0u_sf)+psi_u)

!  Compute friction velocity 
n_iter = 0_iintegers
Delta = 1.0                ! Set initial error to a large value (as compared to the accuracy)
u_star_previter = u_star_thresh  ! Initial guess for friction velocity  
IF(U_a.LE.U_a_thresh) THEN  ! Smooth surface
  DO WHILE (Delta.GT.c_accur_sf.AND.n_iter.LT.n_iter_max) 
    CALL SfcFlx_roughness (fetch, U_a, MIN(u_star_thresh, u_star_previter), h_ice,   &
                           c_z0u_fetch, u_star_thresh, z0u_sf, z0t_sf, z0q_sf)
    IF(ZoL.GE.0.0) THEN  ! Stable stratification
      psi_u = c_MO_u_stab*ZoL*(1.0-MIN(z0u_sf/height_u, 1.0))
      Fun = LOG(height_u/z0u_sf) + psi_u
      Fun_prime = (Fun + 1.0 + c_MO_u_stab*ZoL*MIN(z0u_sf/height_u, 1.0))/c_Karman
      Fun = Fun*u_star_previter/c_Karman - U_a
    ELSE                       ! Convection 
      psi_t = (1.0-c_MO_u_conv*ZoL)**c_MO_u_exp
      psi_q = (1.0-c_MO_u_conv*ZoL*MIN(z0u_sf/height_u, 1.0))**c_MO_u_exp
      psi_u = 2.0*(ATAN(psi_t)-ATAN(psi_q))                  &
            + 2.0*LOG((1.0+psi_q)/(1.0+psi_t))   &
            + LOG((1.0+psi_q*psi_q)/(1.0+psi_t*psi_t))   
      Fun = LOG(height_u/z0u_sf) + psi_u
      Fun_prime = (Fun + 1.0/psi_q)/c_Karman
      Fun = Fun*u_star_previter/c_Karman - U_a
    END IF
    u_star_st = u_star_previter - Fun/Fun_prime
    Delta = ABS((u_star_st-u_star_previter)/(u_star_st+u_star_previter))
    u_star_previter = u_star_st
    n_iter = n_iter + 1_iintegers
  END DO 
ELSE                        ! Rough surface
  DO WHILE (Delta.GT.c_accur_sf.AND.n_iter.LT.n_iter_max) 
    CALL SfcFlx_roughness (fetch, U_a, MAX(u_star_thresh, u_star_previter), h_ice,   &
                           c_z0u_fetch, u_star_thresh, z0u_sf, z0t_sf, z0q_sf)
    IF(ZoL.GE.0.0) THEN  ! Stable stratification
      psi_u = c_MO_u_stab*ZoL*(1.0-MIN(z0u_sf/height_u, 1.0))
      Fun = LOG(height_u/z0u_sf) + psi_u
      Fun_prime = (Fun - 2.0 - 2.0*c_MO_u_stab*ZoL*MIN(z0u_sf/height_u, 1.0))/c_Karman
      Fun = Fun*u_star_previter/c_Karman - U_a
    ELSE                       ! Convection 
      psi_t = (1.0-c_MO_u_conv*ZoL)**c_MO_u_exp
      psi_q = (1.0-c_MO_u_conv*ZoL*MIN(z0u_sf/height_u, 1.0))**c_MO_u_exp
      psi_u = 2.0*(ATAN(psi_t)-ATAN(psi_q))                  &
            + 2.0*LOG((1.0+psi_q)/(1.0+psi_t))   &
            + LOG((1.0+psi_q*psi_q)/(1.0+psi_t*psi_t))   
      Fun = LOG(height_u/z0u_sf) + psi_u
      Fun_prime = (Fun - 2.0/psi_q)/c_Karman
      Fun = Fun*u_star_previter/c_Karman - U_a
    END IF
    IF(h_ice.GE.h_Ice_min_flk) THEN   ! No iteration is required for rough flow over ice
      u_star_st = c_Karman*U_a/MAX(c_small_sf, LOG(height_u/z0u_sf)+psi_u)
      u_star_previter = u_star_st
    ELSE                              ! Iterate in case of open water
      u_star_st = u_star_previter - Fun/Fun_prime
    END IF
    Delta = ABS((u_star_st-u_star_previter)/(u_star_st+u_star_previter))
    u_star_previter = u_star_st
    n_iter = n_iter + 1_iintegers
  END DO 
END IF

!_dbg
!  print*(*,*) 'MO stab. func. psi_u = ', psi_u, '   n_iter = ', n_iter
!  print*(*,*) '   Wind speed = ', U_a, '  u_* = ', u_star_st
!  print*(*,*) '   Fun = ', Fun
!_dbg

!_dbg
!  IF(n_iter.GE.n_iter_max-1_iintegers)  & 
!    print*(*,*) 'u_*: Max No. iters. exceeded (n_iter = ', n_iter, ')!'
!_dbg

!  Momentum flux
Q_mom_tur = -u_star_st*u_star_st

!  Temperature and specific humidity fluxes
CALL SfcFlx_roughness (fetch, U_a, u_star_st, h_ice, c_z0u_fetch, u_star_thresh, z0u_sf, z0t_sf, z0q_sf)
IF(ZoL.GE.0.0) THEN   ! Stable stratification 
  psi_t = c_MO_t_stab*R_z*ZoL*(1.0-MIN(z0t_sf/height_tq, 1.0))
  psi_q = c_MO_q_stab*R_z*ZoL*(1.0-MIN(z0q_sf/height_tq, 1.0))
!_dbg
!  print*(*,*) 'STAB: psi_t = ', psi_t, '   psi_q = ', psi_q
!_dbg
ELSE                        ! Convection 
  psi_u = (1.0-c_MO_t_conv*R_z*ZoL)**c_MO_t_exp
  psi_t = (1.0-c_MO_t_conv*R_z*ZoL*MIN(z0t_sf/height_tq, 1.0))**c_MO_t_exp
!  psi_t = 2.0*LOG((1.0+psi_t)/(1.0+psi_u))
  psi_t = abs(2.0*LOG((1.0+psi_t)/(1.0+psi_u)))
  psi_u = (1.0-c_MO_q_conv*R_z*ZoL)**c_MO_q_exp
  psi_q = (1.0-c_MO_q_conv*R_z*ZoL*MIN(z0q_sf/height_tq, 1.0))**c_MO_q_exp
!  psi_q = 2.0*LOG((1.0+psi_q)/(1.0+psi_u))
  psi_q = abs(2.0*LOG((1.0+psi_q)/(1.0+psi_u)))
!  write(0,*) 'psi_q= ',psi_q
!_dbg
!  print*(*,*) 'CONV: psi_t = ', psi_t, '   psi_q = ', psi_q
!_dbg
END IF 
Q_sen_tur = -(T_a-T_s)*u_star_st*c_Karman/Pr_neutral  &
          / MAX(c_small_sf, LOG(height_tq/z0t_sf)+psi_t)
if(MAX(c_small_sf, LOG(height_tq/z0t_sf)+psi_t) .lt. 10E-6) then
  write(0,*)'inside flake'
  write(0,*) Q_sen_tur, T_a, T_s, u_star_st, c_Karman, Pr_neutral 
  write(0,*) c_small_sf,height_tq,z0t_sf,psi_t
  write(0,*) 'nominator= ', (T_a-T_s)*u_star_st*c_Karman/Pr_neutral
  write(0,*) 'denominator= ',MAX(c_small_sf, LOG(height_tq/z0t_sf)+psi_t)
endif
Q_lat_tur = -(q_a-q_s)*u_star_st*c_Karman/Sc_neutral  &
          / MAX(c_small_sf, LOG(height_tq/z0q_sf)+psi_q)
if(Q_lat_tur .gt. 6.0E-4) then
  Q_lat_tur = -(q_a-q_s)*u_star_st*c_Karman/3.0  &
          / MAX(c_small_sf, LOG(height_tq/z0q_sf)+psi_q)
  write(0,*) 'Q_lat_tur= ',Q_lat_tur
  write(0,135) q_a,q_s,u_star_st,c_Karman
  write(0,136) MAX(c_small_sf,LOG(height_tq/z0q_sf)+psi_q),c_small_sf, LOG(height_tq/z0q_sf),psi_q 
endif
135   format(1x,4(f16.4))
136   format(1x,4(f16.4))

END IF Turb_Fluxes

!------------------------------------------------------------------------------
!  Decide between turbulent, molecular, and convective fluxes
!------------------------------------------------------------------------------

Q_momentum = MIN(Q_mom_tur, Q_mom_mol, Q_mom_con)  ! Momentum flux is negative          
IF(l_conv_visc) THEN    ! Convection, take fluxes that are maximal in magnitude 
  IF(ABS(Q_sen_tur).GE.ABS(Q_sen_con)) THEN
    Q_sensible = Q_sen_tur
  ELSE
    Q_sensible = Q_sen_con
  END IF
  IF(ABS(Q_sensible).LT.ABS(Q_sen_mol)) THEN
    Q_sensible = Q_sen_mol
  END IF
  IF(ABS(Q_lat_tur).GE.ABS(Q_lat_con)) THEN
    Q_latent = Q_lat_tur
  ELSE
    Q_latent = Q_lat_con
  END IF
  IF(ABS(Q_latent).LT.ABS(Q_lat_mol)) THEN
    Q_latent = Q_lat_mol
  END IF
ELSE                    ! Stable or neutral stratification, chose fluxes that are maximal in magnitude 
  IF(ABS(Q_sen_tur).GE.ABS(Q_sen_mol)) THEN 
    Q_sensible = Q_sen_tur
  ELSE 
    Q_sensible = Q_sen_mol    
  END IF
  IF(ABS(Q_lat_tur).GE.ABS(Q_lat_mol)) THEN 
    Q_latent = Q_lat_tur
  ELSE 
    Q_latent = Q_lat_mol  
  END IF
END IF

!------------------------------------------------------------------------------
!  Set output (notice that fluxes are no longer in kinematic units)
!------------------------------------------------------------------------------

Q_momentum = Q_momentum*rho_a 
!Q_sensible = Q_sensible*rho_a*tpsf_c_a_p
!write(0,*) 'Q_sensible= ',Q_sensible

Q_watvap   = Q_latent*rho_a

!Q_latent = tpsf_L_evap
IF(h_ice.GE.h_Ice_min_flk) Q_latent = Q_latent + tpl_L_f   ! Add latent heat of fusion over ice
!Q_latent = Q_watvap*Q_latent
Q_latent = Q_watvap*tpsf_L_evap
if(Q_latent .gt. 2000.00) then
   write(0,145) 'final Q_watvap= ',Q_watvap, 'tpsf_L_evap= ',tpsf_L_evap, 'Q_latent= ', Q_latent 
endif
!Q_latent = Q_watvap*Q_latent
145   format(A17,E12.5,1x,A13,1x,f10.2,1x,A10,1x,E12.4)
! Set "*_sf" variables to make fluxes accessible to driving routines that use "SfcFlx"
u_star_a_sf     = u_star_st 
Q_mom_a_sf      = Q_momentum  
Q_sens_a_sf     = Q_sensible 
Q_lat_a_sf      = Q_latent
Q_watvap_a_sf   = Q_watvap

!write(85,127) Q_sensible, Q_watvap, Q_latent
 127  format(1x, 3(f16.5,1x))

!------------------------------------------------------------------------------
!  End calculations
!==============================================================================

END SUBROUTINE SfcFlx_momsenlat

!==============================================================================

!==============================================================================
! include 'SfcFlx_rhoair.incf'
!------------------------------------------------------------------------------

REAL (KIND = kind_phys) FUNCTION SfcFlx_rhoair (T, q, P)

!------------------------------------------------------------------------------
!
! Description:
!
!  Computes the air density as function 
!  of temperature, specific humidity and pressure.
!  
! Declarations:
!
! Modules used:

!_dm Parameters are USEd in module "SfcFlx".
!_nu USE data_parameters , ONLY : &
!_nu     ireals,                  & ! KIND-type parameter for real variables
!_nu     iintegers                  ! KIND-type parameter for "normal" integer variables

use machine,               only: kind_phys
!==============================================================================

IMPLICIT NONE

!==============================================================================
!
! Declarations
 
!  Input (function argument) 
REAL (KIND = kind_phys), INTENT(IN) ::   &
  T                                 , & ! Temperature [K]
  q                                 , & ! Specific humidity 
  P                                     ! Pressure [N m^{-2} = kg m^{-1} s^{-2}]

!==============================================================================
!  Start calculations
!------------------------------------------------------------------------------

! Air density [kg m^{-3}] 

SfcFlx_rhoair = P/tpsf_R_dryair/T/(1.0+(1.0/tpsf_Rd_o_Rv-1.0)*q)

!------------------------------------------------------------------------------
!  End calculations
!==============================================================================

END FUNCTION SfcFlx_rhoair

!==============================================================================

!==============================================================================
! include 'SfcFlx_roughness.incf'
!------------------------------------------------------------------------------

SUBROUTINE SfcFlx_roughness (fetch, U_a, u_star, h_ice,   & 
                             c_z0u_fetch, u_star_thresh, z0u, z0t, z0q)

!------------------------------------------------------------------------------
!
! Description:
!
!  Computes the water-surface or the ice-surface roughness lengths
!  with respect to wind velocity, potential temperature and specific humidity.
!
!  The water-surface roughness lengths with respect to wind velocity is computed
!  from the Charnock formula when the surface is aerodynamically rough.
!  A simple empirical formulation is used to account for the dependence 
!  of the Charnock parameter on the wind fetch. 
!  When the flow is aerodynamically smooth, the roughness length with respect to 
!  wind velocity is proportional to the depth of the viscous sub-layer.
!  The water-surface roughness lengths for scalars are computed using the power-law 
!  formulations in terms of the roughness Reynolds number (Zilitinkevich et al. 2001).
!  The ice-surface aerodynamic roughness is taken to be constant.
!  The ice-surface roughness lengths for scalars 
!  are computed through the power-law formulations 
!  in terms of the roughness Reynolds number (Andreas 2002).
!
! Declarations:
!
! Modules used:

!_dm Parameters are USEd in module "SfcFlx".
!_nu USE data_parameters , ONLY :   &
!_nu   ireals                     , & ! KIND-type parameter for real variables
!_nu   iintegers                      ! KIND-type parameter for "normal" integer variables

use machine,               only: kind_phys
!==============================================================================

IMPLICIT NONE

!==============================================================================
!
! Declarations

!  Input (procedure arguments)
REAL (KIND = kind_phys), INTENT(IN) ::   &
  fetch                             , & ! Typical wind fetch [m]
  U_a                               , & ! Wind speed [m s^{-1}]
  u_star                            , & ! Friction velocity in the surface air layer [m s^{-1}]
  h_ice                                 ! Ice thickness [m]

!  Output (procedure arguments)
REAL (KIND = kind_phys), INTENT(OUT) ::   &
  c_z0u_fetch                        , & ! Fetch-dependent Charnock parameter
  u_star_thresh                      , & ! Threshold value of friction velocity [m s^{-1}]
  z0u                                , & ! Roughness length with respect to wind velocity [m]
  z0t                                , & ! Roughness length with respect to potential temperature [m]
  z0q                                    ! Roughness length with respect to specific humidity [m]

!  Local variables of type REAL
REAL (KIND = kind_phys) ::    &
  Re_s                   , & ! Surface Reynolds number 
  Re_s_thresh                ! Threshold value of Re_s

!==============================================================================
!  Start calculations
!------------------------------------------------------------------------------

Water_or_Ice: IF(h_ice.LT.h_Ice_min_flk) THEN  ! Water surface  

! The Charnock parameter as dependent on dimensionless fetch
  c_z0u_fetch = MAX(U_a, u_wind_min_sf)**2_iintegers/tpl_grav/fetch  ! Inverse dimensionless fetch
  c_z0u_fetch = c_z0u_rough + c_z0u_ftch_f*c_z0u_fetch**c_z0u_ftch_ex
  c_z0u_fetch = MIN(c_z0u_fetch, c_z0u_rough_L)                      ! Limit Charnock parameter

! Threshold value of friction velocity
  u_star_thresh = (c_z0u_smooth/c_z0u_fetch*tpl_grav*tpsf_nu_u_a)**num_1o3_sf

! Surface Reynolds number and its threshold value
  Re_s = u_star**3_iintegers/tpsf_nu_u_a/tpl_grav
  Re_s_thresh = c_z0u_smooth/c_z0u_fetch

! Aerodynamic roughness
  IF(Re_s.LE.Re_s_thresh) THEN                 
    z0u = c_z0u_smooth*tpsf_nu_u_a/u_star     ! Smooth flow
  ELSE
    z0u = c_z0u_fetch*u_star*u_star/tpl_grav  ! Rough flow
  END IF 
! Roughness for scalars  
  z0q = c_z0u_fetch*MAX(Re_s, Re_s_thresh)
  z0t = c_z0t_rough_1*z0q**c_z0t_rough_3 - c_z0t_rough_2
  z0q = c_z0q_rough_1*z0q**c_z0q_rough_3 - c_z0q_rough_2
  z0t = z0u*EXP(-c_Karman/Pr_neutral*z0t)
  z0q = z0u*EXP(-c_Karman/Sc_neutral*z0q) 

ELSE Water_or_Ice                              ! Ice surface

! The Charnock parameter is not used over ice, formally set "c_z0u_fetch" to its minimum value
  c_z0u_fetch = c_z0u_rough

! Threshold value of friction velocity
  u_star_thresh = c_z0u_smooth*tpsf_nu_u_a/z0u_ice_rough

! Aerodynamic roughness
  z0u = MAX(z0u_ice_rough, c_z0u_smooth*tpsf_nu_u_a/u_star)

! Roughness Reynolds number 
  Re_s = MAX(u_star*z0u/tpsf_nu_u_a, c_accur_sf)

! Roughness for scalars  
  IF(Re_s.LE.Re_z0s_ice_t) THEN 
    z0t = c_z0t_ice_b0t + c_z0t_ice_b1t*LOG(Re_s)
    z0t = MIN(z0t, c_z0t_ice_b0s)
    z0q = c_z0q_ice_b0t + c_z0q_ice_b1t*LOG(Re_s)
    z0q = MIN(z0q, c_z0q_ice_b0s)
  ELSE 
    z0t = c_z0t_ice_b0r + c_z0t_ice_b1r*LOG(Re_s) + c_z0t_ice_b2r*LOG(Re_s)**2_iintegers
    z0q = c_z0q_ice_b0r + c_z0q_ice_b1r*LOG(Re_s) + c_z0q_ice_b2r*LOG(Re_s)**2_iintegers
  END IF
  z0t = z0u*EXP(z0t)
  z0q = z0u*EXP(z0q)

END IF Water_or_Ice

!------------------------------------------------------------------------------
!  End calculations
!==============================================================================

END SUBROUTINE SfcFlx_roughness

!==============================================================================

!==============================================================================
! include 'SfcFlx_satwvpres.incf'
!------------------------------------------------------------------------------

REAL (KIND = kind_phys) FUNCTION SfcFlx_satwvpres (T, h_ice)

!------------------------------------------------------------------------------
!
! Description:
!
!  Computes saturation water vapour pressure 
!  over the water surface or over the ice surface
!  as function of temperature. 
!  
! Declarations:
!
! Modules used:

!_dm Parameters are USEd in module "SfcFlx".
!_nu USE data_parameters  , ONLY : &
!_nu     ireals,                   & ! KIND-type parameter for real variables
!_nu     iintegers                   ! KIND-type parameter for "normal" integer variables

!_dm The variable is USEd in module "SfcFlx".
!_nu USE flake_parameters , ONLY : &
!_nu   h_Ice_min_flk                 ! Minimum ice thickness [m]
use machine,               only: kind_phys

!==============================================================================

IMPLICIT NONE

!==============================================================================
!
! Declarations
 
!  Input (function argument) 
REAL (KIND = kind_phys), INTENT(IN) ::   &
  T                                 , & ! Temperature [K]
  h_ice                                 ! Ice thickness [m]

!  Local parameters
REAL (KIND = kind_phys), PARAMETER ::   &
   b1_vap   = 610.780        , & ! Coefficient [N m^{-2} = kg m^{-1} s^{-2}]
   b3_vap   = 273.160        , & ! Triple point [K]
   b2w_vap  = 17.26938820    , & ! Coefficient (water)
   b2i_vap  = 21.87455840    , & ! Coefficient (ice) 
   b4w_vap  = 35.860         , & ! Coefficient (temperature) [K]
   b4i_vap  = 7.660              ! Coefficient (temperature) [K]

!==============================================================================
!  Start calculations
!------------------------------------------------------------------------------

! Saturation water vapour pressure [N m^{-2} = kg m^{-1} s^{-2}]

IF(h_ice.LT.h_Ice_min_flk) THEN  ! Water surface
  SfcFlx_satwvpres = b1_vap*EXP(b2w_vap*(T-b3_vap)/(T-b4w_vap))
ELSE                             ! Ice surface
  SfcFlx_satwvpres = b1_vap*EXP(b2i_vap*(T-b3_vap)/(T-b4i_vap))
END IF 

!------------------------------------------------------------------------------
!  End calculations
!==============================================================================

END FUNCTION SfcFlx_satwvpres

!==============================================================================

!==============================================================================
! include 'SfcFlx_spechum.incf'
!------------------------------------------------------------------------------

REAL (KIND = kind_phys) FUNCTION SfcFlx_spechum (wvpres, P)

!------------------------------------------------------------------------------
!
! Description:
!
!  Computes specific humidity as function 
!  of water vapour pressure and air pressure. 
!  
! Declarations:
!
! Modules used:

!_dm Parameters are USEd in module "SfcFlx".
!_nu USE data_parameters , ONLY : &
!_nu     ireals,                  & ! KIND-type parameter for real variables
!_nu     iintegers                  ! KIND-type parameter for "normal" integer variables

use machine,               only: kind_phys
!==============================================================================

IMPLICIT NONE

!==============================================================================
!
! Declarations
 
!  Input (function argument) 
REAL (KIND = kind_phys), INTENT(IN) ::   &
  wvpres                            , & ! Water vapour pressure [N m^{-2} = kg m^{-1} s^{-2}]
  P                                     ! Air pressure [N m^{-2} = kg m^{-1} s^{-2}]

!==============================================================================
!  Start calculations
!------------------------------------------------------------------------------

! Specific humidity 

SfcFlx_spechum = tpsf_Rd_o_Rv*wvpres/(P-(1.0-tpsf_Rd_o_Rv)*wvpres)

!------------------------------------------------------------------------------
!  End calculations
!==============================================================================

END FUNCTION SfcFlx_spechum

!==============================================================================

!==============================================================================
! include 'SfcFlx_wvpreswetbulb.incf'
!------------------------------------------------------------------------------

REAL (KIND = ireals) FUNCTION SfcFlx_wvpreswetbulb (T_dry, T_wetbulb, satwvpres_bulb, P)             

!------------------------------------------------------------------------------
!
! Description:
!
!  Computes water vapour pressure as function of air temperature, 
!  wet bulb temperature, satururation vapour pressure at wet-bulb temperature,
!  and air pressure.
!  
! Declarations:
!
! Modules used:

!_dm Parameters are USEd in module "SfcFlx".
!_nu USE data_parameters , ONLY : &
!_nu     ireals,                  & ! KIND-type parameter for real variables
!_nu     iintegers                  ! KIND-type parameter for "normal" integer variables

use machine,               only: kind_phys
!==============================================================================

IMPLICIT NONE

!==============================================================================
!
! Declarations
 
!  Input (function argument) 
REAL (KIND = kind_phys), INTENT(IN) ::   &
  T_dry                             , & ! Dry air temperature [K]
  T_wetbulb                         , & ! Wet bulb temperature [K]
  satwvpres_bulb                    , & ! Satururation vapour pressure at wet-bulb temperature [N m^{-2}]
  P                                     ! Atmospheric pressure [N m^{-2}]

!==============================================================================
!  Start calculations
!------------------------------------------------------------------------------

! Water vapour pressure [N m^{-2} = kg m^{-1} s^{-2}]

SfcFlx_wvpreswetbulb = satwvpres_bulb & 
                     - tpsf_c_a_p*P/tpsf_L_evap/tpsf_Rd_o_Rv*(T_dry-T_wetbulb)


!------------------------------------------------------------------------------
!  End calculations
!==============================================================================

END FUNCTION SfcFlx_wvpreswetbulb

!==============================================================================

END MODULE SfcFlx


MODULE module_FLake
IMPLICIT NONE
CONTAINS

!------------------------------------------------------------------------------

SUBROUTINE flake_interface ( dMsnowdt_in, I_atm_in, Q_atm_lw_in, height_u_in, height_tq_in,        &
                             U_a_in, T_a_in, q_a_in, P_a_in,                                       &
                             
                             depth_w, fetch, depth_bs, T_bs, par_Coriolis, del_time,               &
                             T_snow_in,  T_ice_in,  T_mnw_in,  T_wML_in,  T_bot_in,  T_B1_in,      &
                             C_T_in,  h_snow_in,  h_ice_in,  h_ML_in,  H_B1_in, T_sfc_p,           &
                             ch, cm, albedo_water, water_extinc,   &
                             
                             T_snow_out, T_ice_out, T_mnw_out, T_wML_out, T_bot_out,               & 
                             T_B1_out, C_T_out, h_snow_out, h_ice_out, h_ML_out,                   &
                             H_B1_out, T_sfc_n, hflx_out, evap_out, gflx_out, lflx_out,            &
                             
                             T_bot_2_in, T_bot_2_out,ustar, q_sfc, chh, cmm )

!------------------------------------------------------------------------------
!
! Description:
!
!  The FLake interface is
!  a communication routine between "flake_main"
!  and a prediction system that uses FLake.
!  It assigns the FLake variables at the previous time step 
!  to their input values given by the driving model,
!  calls a number of routines to compute the heat and radiation fluxes,
!  calls "flake_main",
!  and returns the updated FLake variables to the driving model.
!  The "flake_interface" does not contain any Flake physics. 
!  It only serves as a convenient means to organize calls of "flake_main"
!  and of external routines that compute heat and radiation fluxes.
!  The interface may (should) be changed so that to provide 
!  the most convenient use of FLake.
!  Within a 3D atmospheric prediction system,
!  "flake_main" may be called in a DO loop within "flake_interface" 
!  for each grid-point where a lake is present.
!  In this way, the driving atmospheric model should call "flake_interface"
!  only once, passing the FLake variables to "flake_interface" as 2D fields. 
!
!  Lines embraced with "!_tmp" contain temporary parts of the code.
!  These should be removed prior to using FLake in applications.
!  Lines embraced/marked with "!_dev" may be replaced
!  as improved parameterizations are developed and tested.
!  Lines embraced/marked with "!_dm" are DM's comments
!  that may be helpful to a user.
!
use machine,               only: kind_phys

USE data_parameters , ONLY : &
    ireals,                  & ! KIND-type parameter for real variables
    iintegers                  ! KIND-type parameter for "normal" integer variables

USE flake_derivedtypes         ! Definitions of several derived TYPEs

!USE flake_parameters , ONLY :   &
!  tpl_T_f                     , & ! Fresh water freezing point [K]
!  tpl_rho_w_r                 , & ! Maximum density of fresh water [kg m^{-3}]
!  h_Snow_min_flk              , & ! Minimum snow thickness [m]
!  h_Ice_min_flk                   ! Minimum ice thickness [m]

USE flake_paramoptic_ref       ! Reference values of the optical characteristics
                               ! of the lake water, lake ice and snow 

USE flake_albedo_ref           ! Reference values the albedo for the lake water, lake ice and snow

USE flake           , ONLY :    &
  flake_main                , & ! Subroutine, FLake driver
  flake_radflux               , & ! Subroutine, computes radiation fluxes at various depths
                                  !
  T_snow_p_flk, T_snow_n_flk  , & ! Temperature at the air-snow interface [K]
  T_ice_p_flk, T_ice_n_flk    , & ! Temperature at the snow-ice or air-ice interface [K]
  T_mnw_p_flk, T_mnw_n_flk    , & ! Mean temperature of the water column [K]
  T_wML_p_flk, T_wML_n_flk    , & ! Mixed-layer temperature [K]
  T_bot_p_flk, T_bot_n_flk    , & ! Temperature at the water-bottom sediment interface [K]
  T_B1_p_flk, T_B1_n_flk      , & ! Temperature at the bottom of the upper layer of the sediments [K]
  C_T_p_flk, C_T_n_flk        , & ! Shape factor (thermocline)
  h_snow_p_flk, h_snow_n_flk  , & ! Snow thickness [m]
  h_ice_p_flk, h_ice_n_flk    , & ! Ice thickness [m]
  h_ML_p_flk, h_ML_n_flk      , & ! Thickness of the mixed-layer [m]
  H_B1_p_flk, H_B1_n_flk      , & ! Thickness of the upper layer of bottom sediments [m]
                                  !
  Q_snow_flk                  , & ! Heat flux through the air-snow interface [W m^{-2}]
  Q_ice_flk                   , & ! Heat flux through the snow-ice or air-ice interface [W m^{-2}]
  Q_w_flk                     , & ! Heat flux through the ice-water or air-water interface [W m^{-2}]
  Q_bot_flk                   , & ! Heat flux through the water-bottom sediment interface [W m^{-2}]
  I_atm_flk                   , & ! Radiation flux at the lower boundary of the atmosphere [W m^{-2}],
                                  ! i.e. the incident radiation flux with no regard for the surface albedo
  I_snow_flk                  , & ! Radiation flux through the air-snow interface [W m^{-2}]
  I_ice_flk                   , & ! Radiation flux through the snow-ice or air-ice interface [W m^{-2}]
  I_w_flk                     , & ! Radiation flux through the ice-water or air-water interface [W m^{-2}]
  I_h_flk                     , & ! Radiation flux through the mixed-layer-thermocline interface [W m^{-2}]
  I_bot_flk                   , & ! Radiation flux through the water-bottom sediment interface [W m^{-2}]
  I_intm_0_h_flk              , & ! Mean radiation flux over the mixed layer [W m^{-1}]
  I_intm_h_D_flk              , & ! Mean radiation flux over the thermocline [W m^{-1}]
  Q_star_flk                  , & ! A generalized heat flux scale [W m^{-2}]
  u_star_w_flk                , & ! Friction velocity in the surface layer of lake water [m s^{-1}]
  w_star_sfc_flk              , & ! Convective velocity scale, using a generalized heat flux scale [m s^{-1}]
  dMsnowdt_flk                , & ! The rate of snow accumulation [kg m^{-2} s^{-1}]
  T_bot_2_in_flk              


USE SfcFlx          , ONLY :    &
  SfcFlx_lwradwsfc            , & ! Function, returns the surface long-wave radiation flux
  SfcFlx_momsenlat                ! Subroutine, computes fluxes of momentum and of sensible and latent heat

!==============================================================================

IMPLICIT NONE

!==============================================================================
!
! Declarations

!  Input (procedure arguments)

REAL (KIND = kind_phys), INTENT(IN) ::   &
  dMsnowdt_in                       , & ! The rate of snow accumulation [kg m^{-2} s^{-1}]
  I_atm_in                          , & ! Solar radiation flux at the surface [W m^{-2}]
  Q_atm_lw_in                       , & ! Long-wave radiation flux from the atmosphere [W m^{-2}]
  height_u_in                       , & ! Height above the lake surface where the wind speed is measured [m]
  height_tq_in                      , & ! Height where temperature and humidity are measured [m]
  U_a_in                            , & ! Wind speed at z=height_u_in [m s^{-1}]
  T_a_in                            , & ! Air temperature at z=height_tq_in [K]
  q_a_in                            , & ! Air specific humidity at z=height_tq_in
  P_a_in                            , & ! Surface air pressure [N m^{-2} = kg m^{-1} s^{-2}]
  ch                                , &
  cm                                , &
  albedo_water,                       & ! Water surface albedo with respect to the solar radiation
  water_extinc

REAL (KIND = kind_phys), INTENT(IN) ::   &
  depth_w                           , & ! The lake depth [m]
  fetch                             , & ! Typical wind fetch [m]
  depth_bs                          , & ! Depth of the thermally active layer of the bottom sediments [m]
  T_bs                              , & ! Temperature at the outer edge of 
                                        ! the thermally active layer of the bottom sediments [K]
  par_Coriolis                      , & ! The Coriolis parameter [s^{-1}]
  del_time                              ! The model time step [s]

REAL (KIND = kind_phys), INTENT(IN)  :: &
  T_snow_in                        , & ! Temperature at the air-snow interface [K] 
  T_ice_in                         , & ! Temperature at the snow-ice or air-ice interface [K]
  T_mnw_in                         , & ! Mean temperature of the water column [K]
  T_wML_in                         , & ! Mixed-layer temperature [K]
  T_bot_in                         , & ! Temperature at the water-bottom sediment interface [K]
  T_B1_in                          , & ! Temperature at the bottom of the upper layer of the sediments [K]
  C_T_in                           , & ! Shape factor (thermocline)
  h_snow_in                        , & ! Snow thickness [m]
  h_ice_in                         , & ! Ice thickness [m]
  h_ML_in                          , & ! Thickness of the mixed-layer [m]
  H_B1_in                          , & ! Thickness of the upper layer of bottom sediments [m]
  T_sfc_p                          , & ! Surface temperature at the previous time step [K]  
  T_bot_2_in

!  Input/Output (procedure arguments)

!REAL (KIND = ireals), INTENT(INOUT)  :: &
REAL (KIND = kind_phys)                 :: &
  albedo_ice                          , & ! Ice surface albedo with respect to the solar radiation
  albedo_snow                             ! Snow surface albedo with respect to the solar radiation

!TYPE (opticpar_medium), INTENT(INOUT) :: & 
TYPE (opticpar_medium)                :: & 
  opticpar_water                       , & ! Optical characteristics of water
  opticpar_ice                         , & ! Optical characteristics of ice
  opticpar_snow                            ! Optical characteristics of snow 

!  Output (procedure arguments)

REAL (KIND = kind_phys), INTENT(OUT)  :: &
  T_snow_out                        , & ! Temperature at the air-snow interface [K] 
  T_ice_out                         , & ! Temperature at the snow-ice or air-ice interface [K]
  T_mnw_out                         , & ! Mean temperature of the water column [K]
  T_wML_out                         , & ! Mixed-layer temperature [K]
  T_bot_out                         , & ! Temperature at the water-bottom sediment interface [K]
  T_B1_out                          , & ! Temperature at the bottom of the upper layer of the sediments [K]
  C_T_out                           , & ! Shape factor (thermocline)
  h_snow_out                        , & ! Snow thickness [m]
  h_ice_out                         , & ! Ice thickness [m]
  h_ML_out                          , & ! Thickness of the mixed-layer [m]
  H_B1_out                          , & ! Thickness of the upper layer of bottom sediments [m]
  T_sfc_n                           , & ! Updated surface temperature [K]  
  hflx_out                          , & ! sensibl heat flux
  evap_out                          , & ! Latent heat flux
  gflx_out                          , & ! flux from to water
  lflx_out                          , & ! latent heat flux
  T_bot_2_out                       , & ! Bottom temperature
  ustar                             , &
  q_sfc                             , &
  chh                               , &
  cmm

!  Local variables of type REAL

REAL (KIND = kind_phys) ::    &
  Q_momentum             , & ! Momentum flux [N m^{-2}]
  Q_sensible             , & ! Sensible heat flux [W m^{-2}]
  Q_latent               , & ! Latent heat flux [W m^{-2}]
  Q_watvap               , & ! Flux of water vapour [kg m^{-2} s^{-1}]
  Q_w_flux               , & ! flux from ice to water
  rho_a

! ADDED by Shaobo Zhang
LOGICAL lflk_botsed_use
!REAL (KIND = kind_phys) :: T_bot_2_in, T_bot_2_out
REAL (KIND = kind_phys), parameter :: tpl_rho_w_r  = 1.0E+03 
REAL (KIND = kind_phys), parameter :: tpl_T_f      = 273.15
REAL (KIND = kind_phys), parameter :: h_Snow_min_flk = 1.0E-5
REAL (KIND = kind_phys), parameter :: h_Ice_min_flk  = 1.0E-9
!==============================================================================
!  Start calculations
!------------------------------------------------------------------------------
! lflk_botsed_use   = .TRUE.
 lflk_botsed_use   = .FALSE.
!------------------------------------------------------------------------------
!  Set albedos of the lake water, lake ice and snow
!------------------------------------------------------------------------------

! Use default value 
! albedo_water = albedo_water_ref
! Use empirical formulation proposed by Mironov and Ritter (2004) for GME 
!_nu albedo_ice   = albedo_whiteice_ref 
!albedo_ice   = EXP(-c_albice_MR*(tpl_T_f-T_sfc_p)/tpl_T_f)
!albedo_ice   = albedo_whiteice_ref*(1.0-albedo_ice) + albedo_blueice_ref*albedo_ice
! Snow is not considered
!albedo_snow  = albedo_ice  
albedo_ice   = albedo_whiteice_ref
!albedo_snow  = albedo_ice  
albedo_snow  = albedo_drysnow_ref
opticpar_water%extincoef_optic(1) = water_extinc
!write(0,*)'albedo= ',albedo_water,albedo_ice,albedo_snow

!------------------------------------------------------------------------------
!  Set optical characteristics of the lake water, lake ice and snow
!------------------------------------------------------------------------------

! Use default values
opticpar_water = opticpar_water_ref
opticpar_ice   = opticpar_ice_opaque   ! Opaque ice
opticpar_snow  = opticpar_snow_opaque  ! Opaque snow

!print*,'opticpar = ',opticpar_water, opticpar_ice,opticpar_snow

!------------------------------------------------------------------------------
!  Set initial values
!------------------------------------------------------------------------------
!print*,'Inter depth_w=',depth_w
!print*,'Inter depth_bs=',depth_bs

T_snow_p_flk = T_snow_in     
T_ice_p_flk  = T_ice_in         
T_mnw_p_flk  = T_mnw_in        
T_wML_p_flk  = T_wML_in       
T_bot_p_flk  = T_bot_in      
T_B1_p_flk   = T_B1_in       
C_T_p_flk    = C_T_in        
h_snow_p_flk = h_snow_in      
h_ice_p_flk  = h_ice_in       
h_ML_p_flk   = h_ML_in        
H_B1_p_flk   = H_B1_in       
T_bot_2_in_flk = T_bot_2_in

!write(71,120) T_sfc_p,T_mnw_in,T_wML_in,T_bot_in,T_B1_in,T_bot_2_in
 120  format(1x,6(f12.5,1x))
!------------------------------------------------------------------------------
!  Set the rate of snow accumulation
!------------------------------------------------------------------------------

dMsnowdt_flk = dMsnowdt_in  

!------------------------------------------------------------------------------
!  Compute solar radiation fluxes (positive downward)
!------------------------------------------------------------------------------

I_atm_flk = I_atm_in
CALL flake_radflux ( depth_w, albedo_water, albedo_ice, albedo_snow, &
                     opticpar_water, opticpar_ice, opticpar_snow,    &
                     depth_bs )

!------------------------------------------------------------------------------
!  Compute long-wave radiation fluxes (positive downward)
!------------------------------------------------------------------------------

Q_w_flk = Q_atm_lw_in                          ! Radiation of the atmosphere 
Q_w_flk = Q_w_flk - SfcFlx_lwradwsfc(T_sfc_p)  ! Radiation of the surface (notice the sign)

!------------------------------------------------------------------------------
!  Compute the surface friction velocity and fluxes of sensible and latent heat 
!------------------------------------------------------------------------------

CALL SfcFlx_momsenlat ( height_u_in, height_tq_in, fetch,                      &
                        U_a_in, T_a_in, q_a_in, T_sfc_p, P_a_in, h_ice_p_flk,  &
                        Q_momentum, Q_sensible, Q_latent, Q_watvap, q_sfc, rho_a )
!write(0,*)'tpl_rho_w_r= ', tpl_rho_w_r
!write(0,*) 'Q_momentum= ',Q_momentum
u_star_w_flk = SQRT(-Q_momentum/tpl_rho_w_r)
ustar = u_star_w_flk

!------------------------------------------------------------------------------
!  Compute heat fluxes Q_snow_flk, Q_ice_flk, Q_w_flk
!------------------------------------------------------------------------------

Q_w_flk = Q_w_flk - Q_sensible - Q_latent  ! Add sensible and latent heat fluxes (notice the signs)
IF(h_ice_p_flk.GE.h_Ice_min_flk) THEN            ! Ice exists
  IF(h_snow_p_flk.GE.h_Snow_min_flk) THEN        ! There is snow above the ice
    Q_snow_flk = Q_w_flk
    Q_ice_flk  = 0.0
    Q_w_flk    = 0.0
  ELSE                                           ! No snow above the ice
    Q_snow_flk = 0.0
    Q_ice_flk  = Q_w_flk
    Q_w_flk    = 0.0
  END IF
ELSE                                             ! No ice-snow cover
    Q_snow_flk = 0.0
    Q_ice_flk  = 0.0
END IF

!------------------------------------------------------------------------------
!  Advance FLake variables
!------------------------------------------------------------------------------

CALL flake_main ( depth_w, depth_bs, T_bs, par_Coriolis,         &
                  opticpar_water%extincoef_optic(1),             &
                  del_time, T_sfc_p, T_sfc_n, T_bot_2_in_flk,    &
                  T_bot_2_out )

!------------------------------------------------------------------------------
!  Set output values
!------------------------------------------------------------------------------

T_snow_out = T_snow_n_flk  
T_ice_out  = T_ice_n_flk      
T_mnw_out  = T_mnw_n_flk     
T_wML_out  = T_wML_n_flk    
T_bot_out  = T_bot_n_flk   
T_B1_out   = T_B1_n_flk    
C_T_out    = C_T_n_flk     
h_snow_out = h_snow_n_flk   
h_ice_out  = h_ice_n_flk    
h_ML_out   = h_ML_n_flk     
H_B1_out   = H_B1_n_flk    
hflx_out   = Q_sensible
evap_out   = Q_watvap
!evap_out   = Q_latent
gflx_out   = Q_w_flk
lflx_out   = Q_latent
chh        = ch * U_a_in * rho_a
cmm        = cm * U_a_in

!write(72,120) T_sfc_n,T_mnw_out,T_wML_out,T_bot_out,T_B1_out,T_bot_2_out 
!------------------------------------------------------------------------------
!  End calculations
!==============================================================================

END SUBROUTINE flake_interface

END MODULE module_FLake

module module_physics_driver

! # - define DHDEBUG

#include "debug_bitforbit_module_use.inc"

  use machine,               only: kind_phys
  use physcons,              only: con_cp, con_fvirt, con_g, con_rd,    &
                                   con_rv, con_hvap, con_hfus,          &
                                   con_rerth, con_pi, pa2mb, rlapse,    &
                                   con_eps, con_epsm1, con_cliq,        &
                                   con_cvap, con_t0c
#ifndef CCPP
  use physcons,              only: rhc_max, dxmin, dxinv,               &
                                   PQ0,A2A,A3,A4,RHmin
  use cs_conv,               only: cs_convr
  use ozne_def,              only: levozp,  oz_coeff, oz_pres
  use h2o_def,               only: levh2o, h2o_coeff, h2o_pres
  use gfs_fv3_needs,         only: get_prs_fv3, get_phi_fv3
  use module_nst_water_prop, only: get_dtzm_2d
#endif
  use GFS_typedefs,          only: GFS_statein_type, GFS_stateout_type, &
                                   GFS_sfcprop_type, GFS_coupling_type, &
                                   GFS_control_type, GFS_grid_type,     &
                                   GFS_tbd_type,     GFS_cldprop_type,  &
                                   GFS_radtend_type, GFS_diag_type
#ifdef CCPP
  use GFS_typedefs,          only: GFS_interstitial_type
#endif

#ifndef CCPP
  use gfdl_cloud_microphys_mod, only: gfdl_cloud_microphys_driver,      &
				   cloud_diagnosis
  use module_mp_thompson,    only: mp_gt_driver
#endif
  use module_mp_wsm6,        only: wsm6
  use funcphys,              only: ftdp
#ifndef CCPP
  use surface_perturbation,  only: cdfnor
#endif

#ifdef CCPP
  use ccpp_api,              only: ccpp_physics_run
  use CCPP_data,             only: cdata_block
#ifdef OPENMP
  use omp_lib
#endif
#endif


  implicit none


  !--- CONSTANT PARAMETERS
  real(kind=kind_phys), parameter :: hocp    = con_hvap/con_cp
  real(kind=kind_phys), parameter :: qmin    = 1.0e-10
  real(kind=kind_phys), parameter :: qsmall  = 1.0e-20
  real(kind=kind_phys), parameter :: rainmin = 1.0e-13
#ifndef CCPP
  real(kind=kind_phys), parameter :: p850    = 85000.0
#endif
  real(kind=kind_phys), parameter :: epsq    = 1.e-20
  real(kind=kind_phys), parameter :: hsub    = con_hvap+con_hfus
  real(kind=kind_phys), parameter :: czmin   = 0.0001      ! cos(89.994)
#ifndef CCPP
  real(kind=kind_phys), parameter :: one     = 1.0d0, onebg = one/con_g
  real(kind=kind_phys), parameter :: albdf   = 0.06
#endif
  real(kind=kind_phys), parameter :: tf=258.16, tcr=273.16, tcrf=1.0/(tcr-tf) !GF - once MG is CCPP-compliant, wrap this in #ifndef CCPP since moved to GFS_typedefs
  real(kind=kind_phys), parameter :: con_p001= 0.001d0
  real(kind=kind_phys), parameter :: con_d00 = 0.0d0
  real(kind=kind_phys), parameter :: con_day = 86400.d0


!> GFS Physics Implementation Layer
!> @brief Layer that invokes individual GFS physics routines
!> @{
!at tune step===========================================================!
!  description:                                                         !
!                                                                       !
!  usage:                                                               !
!                                                                       !
!    call GFS_physics_driver                                            !
!                                                                       !
!  ---  interface variables                                             !
!     type(GFS_control_type),         intent(in)    :: Model            !
!     type(GFS_statein_type),         intent(inout) :: Statein          !
!     type(GFS_stateout_type),        intent(inout) :: Stateout         !
!     type(GFS_sfcprop_type),         intent(inout) :: Sfcprop          !
!     type(GFS_coupling_type),        intent(inout) :: Coupling         !
!     type(GFS_grid_type),            intent(in)    :: Grid             !
!     type(GFS_tbd_type),             intent(inout  :: Tbd              !
!     type(GFS_cldprop_type),         intent(inout) :: Cldprop          !
!     type(GFS_radtend_type),         intent(inout) :: Radtend          !
!     type(GFS_diag_type),            intent(inout) :: Diag             !
#ifdef CCPP
!     type(GFS_interstitial_type),    intent(inout) :: Interstitial(:)  !
#endif
!                                                                       !
!  subprograms called:                                                  !
!                                                                       !
!     get_prs,  dcyc2t2_pre_rad (testing),    dcyc2t3,  sfc_diff,       !
!     sfc_ocean,sfc_drv,  sfc_land, sfc_sice, sfc_diag, moninp1,        !
!     moninp,   moninq1,  moninq,   satmedmfvdif,                       !
!     gwdps,    ozphys,   get_phi,                                      !
!     sascnv,   sascnvn,  samfdeepcnv, rascnv,   cs_convr, gwdc,        !
!     shalcvt3, shalcv,   samfshalcnv,                                  !
!     shalcnv,  cnvc90,   lrgscl,   gsmdrive, gscond,   precpd,         !
!     progt2.                                                           !
!                                                                       !
!                                                                       !
!  program history log:                                                 !
!           19xx  - ncep mrf/gfs                                        !
!           2002  - s. moorthi  modify and restructure and add Ferrier  !
!                               microphysics as an option               !
!           200x  - h. juang    modify (what?)                          !
!      nov  2004  - x. wu       modify sea-ice model                    !
!      may  2005  - s. moorthi  modify and restructure                  !
!           2005  - s. lu       modify to include noah lsm              !
!      oct  2006  - h. wei      modify lsm options to include both      !
!                               noah and osu lsms.                      !
!           2006  - s. moorthi  added a. johansson's convective gravity !
!                               wave parameterization code              !
!           2007  - s. moorthi  added j. han's modified pbl/sas options !
!      dec  2007  - xu li       modified the operational version for    !
!                               nst model                               !
!           2008  - s. moorthi  applied xu li's nst model to new gfs    !
!      mar  2008  - y.-t. hou   added sunshine duration var (suntim) as !
!                     an input/output argument.                         !
!           2008  - jun wang    added spfhmax/spfhmin as input/output.  !
!      apr  2008  - y.-t. hou   added lw sfc emissivity var (sfcemis),  !
!                     define the lw sfc dn/up fluxes in two forms: atmos!
!                     and ground. also changed sw sfc net flux direction!
!                     (positive) from ground -> atmos to the direction  !
!                     of atmos -> ground. recode the program and add    !
!                     program documentation block.
!           2008/ - s. moorthi and y.t. hou upgraded the code to more   !
!           2009    modern form and changed all the inputs to MKS units.!
!      feb  2009  - s. moorthi  upgraded to add Hochun's gocart changes !
!      jul  2009  - s. moorthi  added rqtk for sela's semi-lagrangian   !
!      aug  2009  - s. moorthi  added j. han and h. pan updated shallow !
!                               convection package                      !
!      sep  2009  - s. moorthi  updated for the mcica (rrtm3) radiation !
!      dec  2010  - sarah lu    lgocart added to input arg;             !
!                               compute dqdt_v if inline gocart is on   !
!      feb  2011  - sarah lu    add the option to update surface diag   !
!                               fields (t2m,q2m,u10m,v10m) at the end   !
!      Jun  2011  - s. moorthi and Xu Li - updated the nst model        !
!                               !
!      sep  2011  - sarah lu    correct dqdt_v calculations             !
!      apr  2012  - henry juang add idea                                !
!      sep  2012  - s. moorthi  merge with operational version          !
!      Mar  2013  - Jun Wang    set idea heating rate to tmp tendency   !
!      May  2013  - Jun Wang    tmp updated after idea phys             !
!      Jun  2013  - s. moorthi  corrected a bug in 3d diagnostics for T !
!      Aug  2013  - s. moorthi updating J. Whitekar's changes related   !
!                              to stochastic physics perturnbation      !
!      Oct  2013  - Xingren Wu  add dusfci/dvsfci                       !
!      Mar  2014  - Xingren Wu  add "_cpl" for coupling                 !
!      Mar  2014  - Xingren Wu  add "nir/vis beam and nir/vis diff"     !
!      Apr  2014  - Xingren Wu  add "NET LW/SW including nir/vis"       !
!      Jan  2014  - Jun Wang    merge Moorthi's gwdc change and H.Juang !
!                               and F. Yang's energy conversion from GWD!
!      jan  2014  - y-t hou     revised sw sfc spectral component fluxes!
!                     for coupled mdl, added estimation of ocean albedo !
!                     without ice contamination.                        !
!      Jun  2014  - Xingren Wu  update net SW fluxes over the ocean     !
!                               (no ice contamination)                  !
!      Jul  2014  - Xingren Wu  add Sea/Land/Ice Mask - slmsk_cpl       !
!      Jul  2014  - s. moorthi  merge with standalone gfs and cleanup   !
!      Aug  2014  - s. moorthi  add tracer fixer                        !
!      Sep  2014  - Sarah Lu    disable the option to compute tracer    !
!                               scavenging in GFS phys (set fscav=0.)   !
!      Dec  2014  - Jun Wang    add cnvqc_v for gocart                  !

!  ====================  defination of variables  ====================  !
!      ---  2014  - D. Dazlich  Added Chikira-Sugiyama (CS) convection  !
!                               as an option in opr GFS.                !
!      Apr  2015    S. Moorthi  Added CS scheme to NEMS/GSM             !
!      Jun  2015    S. Moorthi  Added SHOC  to NEMS/GSM                 !
!      Aug  2015  - Xu  Li      change nst_fcst to be nstf_name         !
!                               and introduce depth mean SST            !
!      Sep  2015  - Xingren Wu  remove oro_cpl & slmsk_cpl              !
!      Sep  2015  - Xingren Wu  add sfc_cice                            !
!      Sep  2015  - Xingren Wu  connect CICE output to sfc_cice         !
!      Jan  2016  - P. Tripp    NUOPC/GSM merge                         !
!      Mar  2016  - J. Han  -   add ncnvcld3d integer                   !
!                               for convective cloudiness enhancement   !
!      Mar  2016  - J. Han  -   change newsas & sashal to imfdeepcnv    !
!                            &  imfshalcnv, respectively                !
!      Mar  2016    F. Yang     add pgr to rayleigh damping call        !
!      Mar  2016    S. Moorthi  add ral_ts                              !
!      Mar  2016    Ann Cheng   add morrison 2m microphysics (gsfc)     !
!      May  2016    S. Moorthi  cleanup 2m microphysics implementation  !
!      Jun  2016    X. Li       change all nst_fld as inout             !
!      jul  2016    S. Moorthi  fix some bugs in shoc/2m microphysics   !
!      au-nv2016a   S. Moorthi  CS with AW and connect with shoc/2m     !
!      Dec  2016    Anning C.   Add prognostic rain and snow with 2M    !
!      Oct  2017    S. Moorthi  fix tracers to account for ice, snow etc!
!                               with this RAS and CSAW advect condensates!
!      Mar  2017    Ruiyu S.    Add Thompson's 2M aerosol MP            !
!      May  2017    Ruiyu S.    Add WSM6 MP                             !
!      Dec  2017    S. Moorthi  Merge/update Ruiyu's update on vertical !
!                               diffusion of tracers for all monins     !
!      Jan 04 2018  S. Moorthi  fix a bug in rhc for use in MG          !
!                               macrophysics and replace ntrac by nvdiff!
!                               in call to moninshoc                    !
!      Jun  2018    J. Han      Add scal-aware TKE-based moist EDMF     !
!                               vertical turbulent mixng scheme         !
!      Nov  2018    J. Han      Add canopy heat storage parameterization!
!      Feb  2019    Ruiyu S.    Add an alternate method to use          ! 
!				hydrometeors from GFDL MP in radiation  !
!
!  ====================    end of description    =====================
!  ====================  definition of variables  ====================  !

!> @details This subroutine is the suite driver for the GFS atmospheric physics and surface.
!! It is responsible for calculating and applying tendencies of the atmospheric state
!! variables due to the atmospheric physics and due to the surface layer scheme. In addition,
!! this routine applies radiative heating rates that were calculated during the
!! antecedent call to the radiation scheme. Code within this subroutine is executed on the
!! physics sub-timestep. The sub-timestep loop is executed in the subroutine gloopb.
!!
!!  \section general General Algorithm
!!  -# Prepare input variables for calling individual parameterizations.
!!  -# Using a two-iteration loop, calculate the state variable tendencies for the surface layer.
!!  -# Calculate the state variable tendencies due to the PBL (vertical diffusion) scheme.
!!  -# Calculate the state variable tendencies due to orographic gravity wave drag and Rayleigh damping.
!!  -# Apply tendencies to the state variables calculated so far:
!!   - for temperature: radiation, surface, PBL, oro. GWD, Rayleigh damping
!!   - for momentum: surface, PBL, oro. GWD, Rayleigh damping
!!   - for water vapor: surface, PBL
!!  -# Calculate and apply the tendency of ozone.
!!  -# Prepare input variables for physics routines that update the state variables within their subroutines.
!!  -# If SHOC is active and is supposed to be called before convection, call it and update the state variables within.
!!  -# Calculate and apply the state variable tendencies (within the subroutine) due to deep convection.
!!  -# Calculate the state variable tendencies due to convective gravity wave drag and apply them afterwards.
!!  -# Calculate and apply the state variable tendencies (within the subroutine) due to shallow convection.
!!  -# If SHOC is active and is supposed to be called after convection, call it and update the state variables within.
!!  -# Prepare for microphysics call by calculating preliminary variables.
!!  -# If necessary, call the moist convective adjustment subroutine and update the state temperature and moisture variable within.
!!  -# Calculate and apply the state variable tendencies (within the subroutine) due to microphysics.
!!  -# Determine the precipitation type and update land surface properties if necessary.
!!  -# Fill the output variables from the local variables as necessary and deallocate allocatable arrays.
!!  \section detailed Detailed Algorithm
!!  ## Prepare input variables for calling individual parameterizations.
!!  Before calling any parameterizations, there is a section at the beginning of the subroutine for
!!  preparing input arguments to the various schemes based on general input to the driver and initializing
!!  variables used throughout the driver.
!!  - General initialization:
!!   - set a flag for running in debug mode and the horizontal index of the column to print
!!   - calculate the pressure at layer centers, the exner function at layer centers and interfaces,
!!     geopotential at layer centers and interfaces, and the layer-centered pressure difference
!!   - calculate the ratio of dynamics time step to physics time step for applying tendencies
!!   - initialize local tendency arrays to zero
!!  - Radiation:
!!   - adjust radiative fluxes and heating rates to the shorter physics time step (from the longer radiation time step),
!!    unless idealized physics is true (lsidea) where radiative heating rates are set to 0
!!   - compute diagnostics from the radiation scheme needed for other schemes (e.g., downward longwave flux absorbed by the surface)
!!   - accumulate the upward and downward longwave fluxes at the surface
!!  - Surface:
!!   - set NOAH and OSU scheme variables from gbphys input arguments and initialize local soil moisture variables
!!   - set local sea ice variables from gbphys arguments
!!   - set up A/O/I coupling variables from gbphys arguments
!!  - PBL:
!!   - set the number of tracers that are diffused vertically
!!  - SHOC:
!!   - determine the index of TKE (ntk) in the convectively transported tracer array (clw)
!!   - allocate precipitation mixing ratio cloud droplet number concentration arrays
!!  - Deep Convection:
!!   - determine which tracers in the tracer input array undergo convective transport (valid for the RAS and Chikira-Sugiyama, and SAMF schemes) and allocate a local convective transported tracer array (clw)
!!   - apply an adjustment to the tracers from the dynamics
!!   - calculate horizontal grid-related parameters needed for some parameterizations
!!   - calculate the maxiumum cloud base updraft speed for the Chikira-Sugiyama scheme
!!   - allocate array for cloud water and cloud cover (for non-RAS and non-Chikira-Sugiyama deep convective schemes)
!!  - Shallow Convection:
!!   - when using the Tiedtke shallow convection scheme with the stratus modifications, find the lowest
!!     model level where a temperature inversion exists in the absence of CTEI
!!  - Microphysics:
!!   - for the Morrison (MGB) scheme, calculate 'FRLAND' if the grid point is over land
!!   - allocate arrays associated with the Morrison scheme
!!   - assign the local critical relative humidity variables from the gbphys arguments
!!  - Gravity Wave Drag:
!!   - calculate the deep convective cloud fraction at cloud top for the convective GWD scheme
!!  .
!!  ## Using a two-iteration loop, calculate the state variable tendencies for the surface layer.
!!   - Each iteration of the loop calls the following:
!!    - 'sfc_diff' to calculate surface exchange coefficients and near-surface wind
!!    - surface energy balances routines are called regardless of surface type; the surface type is checked within each to determine whether the routine is "active"
!!    - for the surface energy balance over the ocean, call 'sfc_nst' if NSST is on, otherwise, call 'sfc_ocean'
!!    - for the surface energy balance over the land, call 'sfc_drv' for the NOAH model and 'sfc_land' for the OSU model
!!    - for the surface energy balance over sea ice, call sfc_sice; if A/O/I coupling, call sfc_cice
!!   - The initial iteration has flag_guess = F unless wind < 2 m/s; flag_iter = T
!!   - After the initial iteration, flag_guess = F and flag_iter = F (unless wind < 2 m/s and over a land surface or an ocean surface with NSST on)
!!   - The following actions are performed after the iteration to calculate surface energy balance:
!!    - set surface output variables from their local values
!!    - call 'sfc_diag' to calculate state variable values at 2 and 10 m as appropriate from near-surface model levels and the surface exchange coefficients
!!    - if A/O/I coupling, set coupling variables from local variables and calculate the open water albedo
!!    - finally, accumulate surface-related diagnostics and calculate the max/min values of T and q at 2 m height.
!!  .
!!  ## Calculate the state variable tendencies due to the PBL (vertical diffusion) scheme.
!!   - Call the vertical diffusion scheme (PBL) based on the following logical flags: do_shoc, hybedmf, satmedmf, old_monin, mstrat
!!    - the PBL scheme is expected to return tendencies of the state variables
!!   - If A/O/I coupling and the surface is sea ice, overwrite some surface-related variables to their states before PBL was called
!!   - For diagnostics, do the following:
!!    - accumulate surface state variable tendencies and set the instantaneous values for output
!!    - accumulate the temperature tendency due to the PBL scheme in dt3dt(:,:,3), subtracting out the radiative heating rate if necessary
!!    - accumulate the u, v tendencies due to the PBL in du3dt(:,:,1:2) and dv3dt(:,:,1:2)
!!    - accumulate the water vapor tendency due to the PBL in dq3dt(:,:,1)
!!    - accumulate the ozone tendency in dq3dt(:,:,5)
!!  .
!!  ## Calculate the state variable tendencies due to orographic gravity wave drag and Rayleigh damping.
!!   - Based on the variable nmtvr, unpack orographic gravity wave varibles from the hprime array
!!   - Call 'gwdps' to calculate tendencies of u, v, T, and surface stress
!!   - Accumulate gravity wave drag surface stresses.
!!   - Accumulate change in u, v, and T due to oro. gravity wave drag in du3dt(:,:,2), dv3dt(:,:,2), and dt3dt(:,:,2)
!!   - Call 'rayleigh_damp' to calculate tendencies to u, v, and T due to Rayleigh friction
!!  .
!!  ## Apply tendencies to the state variables calculated so far.
!!  ## Calculate and apply the tendency of ozone.
!!   - Call the convective adjustment scheme for IDEA
!!   - Call 'ozphys_2015' or 'ozphys' depending on the value of pl_coeff, updating the ozone tracer within and outputing the tendency of ozone in dq3dt(:,:,6)
!!   - Call 'h2ophys' if necessary ("adaptation of NRL H2O phys for stratosphere and mesophere")
!!  .
!!  ## Prepare input variables for physics routines that update the state variables within their subroutines.
!!  - If diagnostics is active, save the updated values of the state variables in 'dudt', 'dvdt', 'dTdt', and 'dqdt(:,:,1)'
!!  - Call 'get_phi' to calculate geopotential from p, q, T
!!  - Initialize the cloud water and ice portions of the convectively transported tracer array (clw) and (if the deep convective scheme is not RAS or Chikira-Sugiyama) the convective cloud water and cloud cover.
!!  - If the dep convective scheme is RAS or Chikira-Sugiyama, fill the 'clw' array with tracers to be transported by convection
!!  - Initialize 'ktop' and 'kbot' (to be modified by all convective schemes except Chikira-Sugiyama)
!!  - Prepare for microphysics call (if cloud condensate is in the input tracer array):
!!   - all schemes: calculate critical relative humidity
!!   - Morrison et al. scheme (occasionally denoted MGB) (when ncld==2): set clw(:,:,1) to cloud ice and clw(:,:,2) to cloud liquid water
!!   - Ferrier scheme (num_p3d==3): set the cloud water variable and separate hydrometeors into cloud ice, cloud water, and rain; set clw(:,:,1) to cloud ice and clw(:,:,2) to cloud liquid water
!!   - Zhao-Carr scheme (num_p3d==4): calculate autoconversion coefficients from input constants and grid info; set set clw(:,:,1) to cloud liquid water
!!   - otherwise: set autoconversion parameters like in Zhao-Carr and set critical relative humidity to 1
!!  .
!!  ##  If SHOC is active and is supposed to be called before convection, call it and update the state variables within.
!!   - Prior to calling SHOC, prepare some microphysics variables:
!!    - if Morrison et al. scheme: set 'skip_macro', fill clw(:,:,1,2) with cloud ice, liquid from the tracer array, and fill cloud droplet number concentration arrays from the input tracer array
!!    - if Zhao-Carr scheme: initialize precip. mixing ratios to 0, fill clw(:,:,1,2) with cloud ice, liquid from the tracer array (as a function of temperature)
!!   - Call 'shoc' (modifies state variables within the subroutine)
!!   - Afterward, set updated cloud number concentrations in the tracer array from the updated 'ncpl' and 'ncpi'
!!  .
!!  ## Calculate and apply the state variable tendencies (within the subroutine) due to deep convection.
!!   - Call deep convective scheme according to the parameter 'imfdeepcnv', 'ras', and 'cscnv'.
!!    - if imfdeepcnv == 0, 1, or 2, no special processing is needed
!!    - if the Chikira-Sugiyama scheme (cscnv), convert rain rate to accumulated rain (rain1)
!!    - if RAS, initialize 'ccwfac', 'dlqfac', and revap before the call to 'rascnv'
!!   - Zero out 'cld1d' (cloud work function calculated in non-RAS, non-Chikira-Sugiyama schemes)
!!   - If 'lgocart', accumulate convective mass fluxes and convective cloud water
!!   - Update tracers in the tracer array (gq0) due to convective transport (RAS, CS only) from the 'clw' array
!!   - Calculate accumulated surface convective precip. for this physics time step (rainc)
!!   - If necessary, accumulate cloud work function, convective precipitation, and convective mass fluxes; accumulate dt3dt(:,:,4), dq3dt(:,:,2), du3dt(:,:,3), dv3dt(:,:,3) as change in state variables due to deep convection
!!   - If 'lgocart', repeat the accumulation of convective mass fluxes and convective cloud water; save convective tendency for water vapor in 'dqdt_v'
!!   - If PDF-based clouds are active and Zhao-Carr microphysics, save convective cloud cover and water in 'phy_f3d' array
!!    - otherwise, if non-PDF-based clouds and the "convective cloudiness enhancement" is active, save convective cloud water in 'phy_f3d' array
!!  .
!!  ## Calculate the state variable tendencies due to convective gravity wave drag and apply them afterwards.
!!   - Calculate the average deep convective heating rate in the column to pass into 'gwdc'
!!   - Call 'gwdc' to calculate tendencies of u, v due to convective GWD
!!   - For diagnostics, accumulate the vertically-integrated change in u, v due to conv. GWD; accumulate change in u, v, due to conv. GWD in du3dt(:,:,4) and dv3dt(:,:,4)
!!   - Calculate updated values of u, v, T using conv. GWD tendencies
!!  .
!!  ## Calculate and apply the state variable tendencies (within the subroutine) due to shallow convection.
!!   - If diagnostics are active, set 'dtdt' and 'dqdt' to updated values of T and q before shallow convection
!!   - If SHOC is not active, do the following:
!!    - for the mass-flux shallow convection scheme (imfshalcnv == 1), call 'shalcnv'
!!    - for the scale- and aerosol-aware scheme (imfshalcnv == 2), call 'samfshalcnv'
!!    - for either of the first two schemes, perform the following after the call:
!!     - if Zhao-Carr microphysics with PDF-based clouds, save convective cloud water an cover in 'phy_f3d'
!!     - if non-PDF-based clouds and convective cloudiness enhancement is active, save convective cloud water in 'phy_f3d'
!!     - calculate shallow convective precip. and add it to convective precip; accumulate convective precip.
!!    - for the Tiedtke scheme (imfshalcnv == 0), find the top level where shallow convection must stratosphere
!!     - if using Moorthi's approach to stratus, call 'shalcv'
!!     - otherwise, call 'shalcvt3'
!!    - for diagnostics, accumulate the change in water vapor due to shallow convection and save in dqdt_v if 'lgocart';
!!     - save the change in T and q due to shallow convection in dt3dt(:,:,5) and dq3dt(:,:,3); reset dtdt and dqdt to the updated values of T, q after shallow Convection
!!     - if 'clw' is not partitioned into ice/water, set 'clw(ice)' to zero
!!   - If SHOC is active (and shocaftcnv)
!!    - if Morrison et al. scheme: set 'skip_macro' and fill cloud droplet number concentration arrays from the input tracer array
!!    - initialize precip. mixing ratios to 0
!!    - call 'shoc' (modifies state variables within the subroutine)
!!    - afterward, set updated cloud number concentrations in the tracer array from the updated 'ncpl' and 'ncpi'
!!  .
!!  ## Prepare for microphysics call by calculating preliminary variables.
!!   - For Morrison et al. microphysics, set cloud water and ice arrays to the convecitvely transported values
!!   - For Ferrier microphysics, combine convectively transported cloud water and ice with column rain and save in cloud water array
!!    - calculate and save ice fraction and rain fraction in phy_f3d(1),(2)
!!   - For Zhao-Carr, combine convectively transported cloud water and ice into the cloud water array
!!   - Otherwise, combine convectively transported cloud water and ice into the convectively transported cloud water
!!   - Call 'cnvc90'; a "legacy routine which determines convective clouds"; outputs 'acv','acvb','acvt','cv','cvb','cvt'
!!  .
!!  ## If necessary, call the moist convective adjustment subroutine and update the state temperature and moisture variable within.
!!   - Updates T, q, 'rain1', cloud water array
!!   - Accumulate convective precip
!!   - For diagnostics, accumulate the change in T, q due to moist convective adjustment; reset 'dtdt' and 'dqdt' to updated T, q before call to microphysics
!!  .
!!  ## Calculate and apply the state variable tendencies (within the subroutine) due to microphysics.
!!   - If 'lgocart', calculate instantaneous moisture tendency in dqdt_v
!!   - If no cloud microphysics (ncld == 0), call 'lrgscl' to update T, q and output large scale precipitation and cloud water
!!   - Otherwise, a more advanced microphysics scheme is called (which scheme depends on values of 'num_p3d','do_shoc',and 'ncld')
!!   - Ferrier scheme (num_p3d == 3):
!!    - calculate droplet number concentration and minimum large ice fraction
!!    - call 'gsmdrive' (modifies T, q, cloud water, 'f_ice', 'f_rain', 'f_rimef', 'rain1')
!!   - Zhao-Carr-Sundqvist scheme (num_p3d == 4):
!!    - if non-PDF-based clouds:
!!     - if 'do_shoc', call 'precpd_shoc' (precpd modified for SHOC)
!!     - else, call 'gscond' (grid-scale condensation/evaporation); updates water vapor, cloud water, temperature
!!      - call 'precpd'; updates water vapor, cloud water, temperature and outputs precip., snow ratio, and rain water path
!!    - for PDF-based clouds:
!!     - call 'gscondp' followed by 'precpdp' (similar arguments to gscond, precpd above)
!!   - Morrison et al. scheme (ncld = 2):
!!    - if 'do_shoc', set clw(1),(2) from updated values; set phy_f3d(:,:,1) from phy_f3d(:,:,ntot3d-2)
!!    - else, set clw(1),(2) from updated values; set phy_f3d(:,:,1) to cloud cover from previous time step + convective cloud water from convective scheme
!!    - call 'm_micro_driver'; updates water vapor, temperature, droplet number concentrations, cloud cover
!!   - Combine large scale and convective precip.
!!   - For diagnostics, accumulate total surface precipitation and accumulate change in T and q due to microphysics in dt3dt(:,:,6) and dq3dt(:,:,4)
!!  .
!!  ## Determine the precipitation type and update land surface properties if necessary.
!!   - If 'cal_pre', diagnose the surface precipitation type
!!    - call 'calpreciptype'; set snow flag to 1 if snow or sleet, 0 otherwise
!!   - For rain/snow decision, calculate temperature at 850 hPa (\f$T_{850}\f$)
!!    - If not 'cal_pre', set snow flag to 1 if \f$T_{850}\f$ is below freezing
!!   - For coupling, accumulate rain if \f$T_{850}\f$ is above freezing, otherwise accumulate snow
!!   - If using the OSU land model, accumulate surface snow depth if \f$T_{850}\f$ is below freezing and not over an ocean surface
!!    - call 'progt2' (canopy and soil moisture?) and set the soil liquid water equal to soil total water
!!    - if 'lgocart', call 'sfc_diag' to update near-surface state variables (this "allows gocart to use filtered wind fields")
!!   - If necessary (lssav), update the 2m max/min values of T and q
!!   - If necessary (lssav), accumulate total runoff and surface runoff.
!!  .
!!  ## Fill the output variables from the local variables as necessary and deallocate allocatable arrays.
!!   - Set global sea ice thickness and concentration as well as the temperature of the sea ice
!!   - Set global soil moisture variables
!!   - Calculate precipitable water and water vapor mass change due to all physics for the column
!!   - Deallocate arrays for SHOC scheme, deep convective scheme, and Morrison et al. microphysics


  public GFS_physics_driver

  CONTAINS
!*******************************************************************************************

    subroutine GFS_physics_driver                                 &
                   (Model, Statein, Stateout, Sfcprop, Coupling,  &
#ifdef CCPP
                   Grid, Tbd, Cldprop, Radtend, Diag, Interstitial)
#else
                   Grid, Tbd, Cldprop, Radtend, Diag)
#endif

      implicit none
!
!  ---  interface variables
! DH* gfortran correctly throws an error if the intent() declarations
! for arguments differ between the actual routine (here) and the dummy
! interface routine (IPD_func0d_proc in IPD_typedefs.F90):
!
! Error: Interface mismatch in procedure pointer assignment at (1): INTENT mismatch in argument 'control'
!
! Since IPD_func0d_proc declares all arguments as intent(inout), we
! need to do the same here - however, this way we are loosing the
! valuable information on the actual intent to this routine. *DH
#ifdef __GFORTRAN__
      type(GFS_control_type),         intent(inout) :: Model
      type(GFS_statein_type),         intent(inout) :: Statein
      type(GFS_stateout_type),        intent(inout) :: Stateout
      type(GFS_sfcprop_type),         intent(inout) :: Sfcprop
      type(GFS_coupling_type),        intent(inout) :: Coupling
      type(GFS_grid_type),            intent(inout) :: Grid
      type(GFS_tbd_type),             intent(inout) :: Tbd
      type(GFS_cldprop_type),         intent(inout) :: Cldprop
      type(GFS_radtend_type),         intent(inout) :: Radtend
      type(GFS_diag_type),            intent(inout) :: Diag
#ifdef CCPP
      type(GFS_interstitial_type),    intent(inout) :: Interstitial(:)
#endif
#else
      type(GFS_control_type),         intent(in)    :: Model
      type(GFS_statein_type),         intent(inout) :: Statein
      type(GFS_stateout_type),        intent(inout) :: Stateout
      type(GFS_sfcprop_type),         intent(inout) :: Sfcprop
      type(GFS_coupling_type),        intent(inout) :: Coupling
      type(GFS_grid_type),            intent(in)    :: Grid
      type(GFS_tbd_type),             intent(inout) :: Tbd
      type(GFS_cldprop_type),         intent(inout) :: Cldprop
      type(GFS_radtend_type),         intent(inout) :: Radtend
      type(GFS_diag_type),            intent(inout) :: Diag
#ifdef CCPP
      type(GFS_interstitial_type),    intent(inout) :: Interstitial(:)
#endif
#endif
!
!  ---  local variables

!--- INTEGER VARIABLES
      integer :: me, ipr, ix, im, levs, ntrac, nvdiff, kdt,             &
                 ntoz, ntcw, ntiw, ncld,ntke,ntkev, ntlnc, ntinc, lsoil,&
                 ntrw, ntsw, ntrnc, ntsnc, ntot3d, ntgl, ntgnc, ntclamt,&
                 ims, ime, kms, kme, its, ite, kts, kte, imp_physics,   &
                 ntwa, ntia

#ifdef CCPP
      integer :: i, kk, ic, k, n, iter, levshcm, tracers,               &
                 tottracer, nsamftrac, num2, num3, ntk,                 &
                 nn, nncl, ntiwx
#else
      integer :: i, kk, ic, k1, k, n, iter, levshcm, tracers,           &
                 tottracer, nsamftrac, num2, num3, nshocm, nshoc, ntk,  &
                 nn, nncl, ntiwx, seconds
#endif

      integer, dimension(size(Grid%xlon,1)) ::                          &
           kbot, ktop, kcnv, soiltyp, vegtype, kpbl, slopetyp, kinver,  &
           levshc, islmsk,                                              &
!--- coupling inputs for physics
           islmsk_cice

!--- LOGICAL VARIABLES
      logical :: lprnt, revap, mg3_as_mg2, skip_macro

      logical, dimension(size(Grid%xlon,1)) ::                          &
           flag_iter, flag_guess, invrsn,                               &
!--- coupling inputs for physics
           flag_cice

      logical, dimension(Model%ntrac+1,2) :: otspt

      real(kind=kind_phys), dimension(Model%ntrac+2) :: trcmin

!--- REAL VARIABLES
#ifdef CCPP
      real(kind=kind_phys) ::                                           &
           dtf, dtp, rhbbot, rhbtop, rhpbl, frain, tem, tem1,           &
           dpshc!,                                                      &
!--- experimental for shoc sub-stepping
!           dtshoc
#else
      real(kind=kind_phys) ::                                           &
           dtf, dtp, rhbbot, rhbtop, rhpbl, frain, tem, tem1, tem2,     &
           xcosz_loc, zsea1, zsea2, eng0, eng1, dpshc,                  &
!--- experimental for shoc sub-stepping
           dtshoc
!--- GFDL Cloud microphysics
      real(kind=kind_phys) ::                                           &
           crain, csnow, total_precip
#endif

      real(kind=kind_phys), dimension(Model%ntrac-Model%ncld+2) ::      &
           fscav, fswtr

      real(kind=kind_phys), dimension(size(Grid%xlon,1))  ::            &
#ifdef CCPP
           ccwfac, garea, dlength, cumabs, gflx,                        &
#else
           ccwfac, garea, dlength, cumabs, cice, zice, tice, gflx,      &
#endif
           rain1, raincs, snowmt, cd, cdq, qss, dusfcg, dvsfcg, dusfc1, &
           dvsfc1,  dtsfc1, dqsfc1, rb, drain,  cld1d, evap, hflx,      &
           stress, t850, ep1d, gamt, gamq, sigmaf, oc, theta, gamma,    &
           sigma, elvmax, wind, work1, work2, runof, xmu, fm10, fh2,    &
           tsurf,  tx1, tx2, ctei_r, evbs, evcw, trans, sbsno, snowc,   &
#ifdef CCPP
           frland,                                                      &
#else
           frland, adjsfcdsw, adjsfcnsw, adjsfcdlw, adjsfculw,          &
#endif
           adjnirbmu, adjnirdfu, adjvisbmu, adjvisdfu, adjnirbmd,       &
           adjnirdfd, adjvisbmd, adjvisdfd, gabsbdlw, xcosz, tseal,     &
           snohf, dlqfac, work3, ctei_rml, cldf,                        &
#ifdef CCPP
           psautco_l, prautco_l, dtzm,                                  &
#else
           psautco_l, prautco_l, ocalnirbm_cpl, ocalnirdf_cpl,          &
           ocalvisbm_cpl, ocalvisdf_cpl, dtzm, temrain1,                &
#endif
           psaur_l, praur_l,                                            &
!--- coupling inputs for physics
           dtsfc_cice, dqsfc_cice, dusfc_cice, dvsfc_cice, ulwsfc_cice, &
           tisfc_cice, tsea_cice, hice_cice, fice_cice,                 &
!--- for CS-convection
           wcbmax

#ifndef CCPP
!--- for precipitation type algorithm only
      real(kind=kind_phys), dimension(size(Grid%xlon,1))  ::            &
           domr, domzr, domip, doms

!--- for GFDL MP only
      real(kind=kind_phys), dimension(size(Grid%xlon,1),1) ::           &
          area, land, rain0, snow0, ice0, graupel0
#endif

      real(kind=kind_phys), dimension(size(Grid%xlon,1),4) ::           &
           oa4, clx

#ifdef CCPP
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levs) ::  &
          del, rhc, dtdt, dudt, dvdt, gwdcu, gwdcv, dtdtc, rainp,       &
          ud_mf, dd_mf, dt_mf, prnum, sigmatot, sigmafrac
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levs-1) ::&
          dkt
#else
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levs) ::  &
          del, rhc, dtdt, dudt, dvdt, gwdcu, gwdcv, dtdtc, rainp,       &
          ud_mf, dd_mf, dt_mf, prnum, dkt, sigmatot, sigmafrac
#endif

!--- for isppt_deep
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levs) ::  &
          savet,saveq,saveu,savev

!--- GFDL modification for FV3
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levs+1) ::&
           del_gz
#ifndef CCPP
      real(kind=kind_phys), allocatable, dimension(:,:,:) ::            &
           delp, dz, uin, vin, pt, qv1, ql1, qr1, qg1, qa1, qn1, qi1,   &
           qs1, pt_dt, qa_dt, udt, vdt, w, qv_dt, ql_dt, qr_dt, qi_dt,  &
           qs_dt, qg_dt,p123,refl
      real(kind=kind_phys), allocatable, dimension(:,:) ::              &
           den
#endif
!
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levs,Model%ntrac) :: &
           dqdt

#ifndef CCPP
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levs,9) ::  &
           dq3dt_loc
#endif

!  mg, sfc perts
#ifdef CCPP
      real (kind=kind_phys), dimension(size(Grid%xlon,1)) :: &
         z01d, zt1d, bexp1d, xlai1d, vegf1d
#else
      real (kind=kind_phys), dimension(size(Grid%xlon,1)) :: &
         z01d, zt1d, bexp1d, xlai1d, alb1d, vegf1d
      real(kind=kind_phys) :: cdfz
#endif
!--- ALLOCATABLE ELEMENTS
      !--- in clw, the first two varaibles are cloud water and ice.
      !--- from third to ntrac are convective transportable tracers,
      !--- third being the ozone, when ntrac=3 (valid with ras, csaw, or samf)
      !--- Anning Cheng 9/21/2016 leave a hook here for diagnosed snow,
      !--- rain, and their numbers
      real(kind=kind_phys), allocatable ::                              &
           clw(:,:,:), qrn(:,:),  qsnw(:,:), ncpl(:,:), ncpi(:,:),      &
           ncpr(:,:),  ncps(:,:), cnvc(:,:), cnvw(:,:),                 &
           qgl(:,:),   ncgl(:,:)
!--- for 2 M microphysics
      real(kind=kind_phys), allocatable, dimension(:,:) ::              &
             qlcn, qicn, w_upi, cf_upi, CNV_MFD,           CNV_DQLDT,   &
             CLCN, CNV_FICE, CNV_NDROP, CNV_NICE
!#ifndef CCPP
       real(kind=kind_phys),parameter :: slope_mg = 0.02, slope_upmg = 0.04,  &
                          turnrhcrit = 0.900, turnrhcrit_upper = 0.150
!#endif
!
      !--- for 2 M Thmpson MP
      real(kind=kind_phys), allocatable, dimension(:,:,:) ::            &
            vdftra, dvdftra
      real(kind=kind_phys), allocatable, dimension(:,:)   ::            &
            ice00, liq0
!     real(kind=kind_phys), allocatable, dimension(:) ::  nwfa2d
#ifndef CCPP
      real(kind=kind_phys), parameter :: liqm = 4./3.*con_pi*1.e-12,    &
                              icem = 4./3.*con_pi*3.2768*1.e-14*890.
      real, allocatable, dimension(:) :: refd,REFD263K
      integer :: nsteps_per_reset
      integer :: kdtminus1
! For computing saturation vapor pressure and rh at 2m
! DH* real(kind_phys) ???
      real :: pshltr,QCQ,rh02
! *DH
#endif
#ifdef CCPP
      integer :: nb
      integer :: nt
      integer :: ierr
      character(len=512) :: errmsg
      integer            :: errflg

      errmsg = ''
      errflg = 0

      ! Retrieve block number and OpenMP thread number
      nb = Tbd%blkno
#ifdef OPENMP
      nt = OMP_GET_THREAD_NUM() + 1
#else
      nt = 1
#endif
      ! If this run uses non-uniform block sizes and the last block is being
      ! computed, then switch from the current Interstitial DDT element to the
      ! additional/special last element of the Interstitial DDT, which has
      ! the correct (smaller) blocksize for this block.
      if (Interstitial(nt)%non_uniform_blocks .and. nb==size(Model%blksz)) then
         nt = size(Interstitial(:))
      end if
#endif

      ! Initialize local variables (mainly for debugging purposes, because the
      ! corresponding variables Interstitial(nt)%... are reset to zero every time);
      ! these variables are only modified over parts of the entire domain (related
      ! to land surface mask etc.)
      snowmt = 0.
      gamq   = 0.
      gamt   = 0.
      gflx   = 0.
      hflx   = 0.

#ifdef CCPP
      ! Associate construct to replace local variables with CCPP variables wrapping the entire driver routine
      associate_ccpp_vars: associate(adjsfcdsw => Diag%dswsfci,               &
                                     adjsfcnsw => Diag%nswsfci,               &
                                     adjsfcdlw => Diag%dlwsfci,               &
                                     adjsfculw => Diag%ulwsfci,               &
                                     levozp    => Interstitial(nt)%levozp,    &
                                     oz_coeff  => Interstitial(nt)%oz_coeff,  &
                                     oz_pres   => Interstitial(nt)%oz_pres,   &
                                     levh2o    => Interstitial(nt)%levh2o,    &
                                     h2o_coeff => Interstitial(nt)%h2o_coeff, &
                                     h2o_pres  => Interstitial(nt)%h2o_pres,  &
                                     dxmin     => Model%dxmin,                &
                                     dxinv     => Model%dxinv,                &
                                     rhc_max   => Model%rhcmax                )

      if (Model%me==0) write(0,*) 'CCPP DEBUG: calling GFS_suite_interstitial_phys_reset through option B'
      ! Copy local variables from driver to appropriate interstitial variables
      !cdata_block(nb,nt)%errmsg = errmsg        ! intent(out)
      !cdata_block(nb,nt)%errflg = errflg        ! intent(out)
      call ccpp_physics_run(cdata_block(nb,nt), scheme_name="GFS_suite_interstitial_phys_reset", ierr=ierr)
      ! Copy intent(inout) and intent(out) interstitial variables to local variables in driver
      errmsg = trim(cdata_block(nb,nt)%errmsg)
      errflg = cdata_block(nb,nt)%errflg
      if (errflg/=0) then
          write(0,*) 'Error in call to GFS_suite_interstitial_phys_reset: ' // trim(errmsg)
          stop
      end if
      if (Model%me==0) write(0,*) 'CCPP DEBUG: calling GFS_suite_stateout_reset through option B'
      !Interstitial(nt)%im = im              ! intent(in) - set in Interstitial(nt)%create()
      !Model%levs                            ! intent(in)
      !Model%ntrac                           ! intent(in)
      !Statein%tgrs                          ! intent(in)
      !Statein%ugrs                          ! intent(in)
      !Statein%vgrs                          ! intent(in)
      !Statein%qgrs                          ! intent(in)
      !Stateout%gt0                          ! intent(out)
      !Stateout%gu0                          ! intent(out)
      !Stateout%gv0                          ! intent(out)
      !Stateout%gq0                          ! intent(out)
      !cdata_block(nb,nt)%errmsg = errmsg    ! intent(out)
      !cdata_block(nb,nt)%errflg = errflg    ! intent(out)
      call ccpp_physics_run(cdata_block(nb,nt), scheme_name="GFS_suite_stateout_reset", ierr=ierr)
      ! Copy intent(inout) and intent(out) interstitial variables to local variables in driver
      errmsg = trim(cdata_block(nb,nt)%errmsg)
      errflg = cdata_block(nb,nt)%errflg
      if (errflg/=0) then
          write(0,*) 'Error in call to GFS_suite_stateout_reset: ' // trim(errmsg)
          stop
      end if
#else
      ! DH Strictly speaking, this is not required. But when hunting for bit-for-bit differences,
      ! doing the same as in GFS_suite_stateout_reset makes life a lot easier.
      Stateout%gt0(:,:)   = Statein%tgrs(:,:)
      Stateout%gu0(:,:)   = Statein%ugrs(:,:)
      Stateout%gv0(:,:)   = Statein%vgrs(:,:)
      Stateout%gq0(:,:,:) = Statein%qgrs(:,:,:)
#endif

#include "debug_bitforbit_diagtoscreen.inc"

!
!===> ...  begin here

      me      = Model%me
#ifdef CCPP
      ix      = Interstitial(nt)%ix
      im      = Interstitial(nt)%im
      ipr     = Interstitial(nt)%ipr
#else
      ix      = size(Grid%xlon,1)
      im      = size(Grid%xlon,1)
      ipr     = min(im,10)
#endif
      levs    = Model%levs
      lsoil   = Model%lsoil
      ntrac   = Model%ntrac
      dtf     = Model%dtf
      dtp     = Model%dtp
      kdt     = Model%kdt
      lprnt   = Model%lprnt
      nvdiff  = ntrac           ! vertical diffusion of all tracers!
      ntcw    = Model%ntcw
      ntoz    = Model%ntoz
      ntiw    = Model%ntiw
      ncld    = Model%ncld
      ntke    = Model%ntke
!
      ntlnc   = Model%ntlnc
      ntinc   = Model%ntinc
      ntrw    = Model%ntrw
      ntsw    = Model%ntsw
      ntrnc   = Model%ntrnc
      ntsnc   = Model%ntsnc
      ntgl    = Model%ntgl
      ntgnc   = Model%ntgnc
      ntclamt = Model%ntclamt
      ntot3d  = Model%ntot3d
      ntwa    = Model%ntwa
      ntia    = Model%ntia

      imp_physics = Model%imp_physics

#ifdef CCPP
      nncl = Interstitial(nt)%nncl
      nvdiff = Interstitial(nt)%nvdiff
      mg3_as_mg2 = Interstitial(nt)%mg3_as_mg2
      ntkev = Interstitial(nt)%ntkev
#else
      nncl = ncld

      if (imp_physics == Model%imp_physics_thompson) then
        if (Model%ltaerosol) then
          nvdiff = 8
        else
          nvdiff = 5
        endif
        if (Model%satmedmf) then
          nvdiff = nvdiff + 1
        endif
        nncl = 5
      elseif (imp_physics == Model%imp_physics_wsm6) then
        nvdiff = ntrac -3
        nncl = 5
      elseif (ntclamt > 0) then             ! for GFDL MP don't diffuse cloud amount
        nvdiff = ntrac - 1
      endif

      if (imp_physics == Model%imp_physics_gfdl) then
        nncl = 5
      endif

      if (imp_physics == Model%imp_physics_mg) then
        if (abs(Model%fprcp) == 1) then
          nncl = 4                          ! MG2 with rain and snow
          mg3_as_mg2 = .false.
        elseif (Model%fprcp >= 2) then
          if(ntgl > 0 .and. (Model%mg_do_graupel .or. Model%mg_do_hail)) then
            nncl = 5                        ! MG3 with rain and snow and grapuel/hail
            mg3_as_mg2 = .false.
          else                              ! MG3 code run without graupel/hail i.e. as MG2
            nncl = 4
            mg3_as_mg2 = .true.
          endif
        endif
      endif

      ntkev = nvdiff
#endif

#ifndef CCPP
      nsteps_per_reset=nint(Model%avg_max_length/dtp)
      kdtminus1=kdt-1
#endif
!
!-------------------------------------------------------------------------------------------
      lprnt   = .false.

!     do i=1,im
!       lprnt = kdt >= 1 .and. abs(grid%xlon(i)*57.29578-119.78) < 0.101      &
!                        .and. abs(grid%xlat(i)*57.29578-19.49)  < 0.101
!       lprnt = kdt >= 250 .and. abs(grid%xlon(i)*57.29578-227.34) < 0.101      &
!                        .and. abs(grid%xlat(i)*57.29578-6.206)  < 0.101
!       lprnt = kdt >= 0 .and. abs(grid%xlon(i)*57.29578-90.9375) < 0.501      &
!                        .and. abs(grid%xlat(i)*57.29578-36.0)  < 0.501
!       lprnt = kdt >= 0 .and. abs(grid%xlon(i)*57.29578-285.938) < 0.501      &
!                        .and. abs(grid%xlat(i)*57.29578+46.286)  < 0.501
!       lprnt = kdt >= 0 .and. abs(grid%xlon(i)*57.29578-108.41) < 0.501      &
!                        .and. abs(grid%xlat(i)*57.29578-32.97)  < 0.501
!       if (kdt == 1) &
!         write(2000+me,*)' i=',i,' xlon=',grid%xlon(i)*57.29578,             &
!                       ' xlat=',grid%xlat(i)*57.29578,' me=',me
!       if (lprnt) then
!         ipr = i
!         write(0,*)' ipr=',ipr,'xlon=',grid%xlon(i)*57.29578,' xlat=',grid%xlat(i)*57.29578,' me=',me
!         exit
!       endif
!     enddo
!     lprnt = .false.
!     if (lprnt) write(0,*)' cloudsphysdriver=',Tbd%phy_f3d(ipr,:,1)*100,' kdt=',kdt
!-------------------------------------------------------------------------------------------
!
      skip_macro = .false.
#ifdef CCPP
      nn = Interstitial(nt)%nn
#else
      if (ntiw > 0) then
        if (ntclamt > 0) then
          nn = ntrac - 2
        else
          nn = ntrac - 1
        endif
      elseif (ntcw > 0) then
        nn = ntrac
      else
        nn = ntrac + 1
      endif
#endif
      allocate (clw(ix,levs,nn))

      num2 = Model%ncnvw
      num3 = Model%ncnvc

      if (Model%imfdeepcnv >= 0 .or.  Model%imfshalcnv > 0  .or. &
         (Model%npdf3d == 3     .and. Model%num_p3d   == 4) .or. &
         (Model%npdf3d == 0     .and. Model%ncnvcld3d == 1) ) then
        allocate (cnvc(ix,levs), cnvw(ix,levs))
        do k=1,levs
          do i=1,im
            cnvc(i,k)  = 0.0
            cnvw(i,k)  = 0.0
          enddo
        enddo
      endif
!
#ifndef CCPP
!  ---  set initial quantities for stochastic physics deltas
      if (Model%do_sppt) then
        Tbd%dtdtr = 0.0
        do i=1,im
          Tbd%drain_cpl(i) = Coupling%rain_cpl (i)
          Tbd%dsnow_cpl(i) = Coupling%snow_cpl (i)
        enddo
      endif

! mg, sfc-perts
!  ---  scale random patterns for surface perturbations with perturbation size
!  ---  turn vegetation fraction pattern into percentile pattern
      do i=1,im
         z01d(i)   = 0.
         zt1d(i)   = 0.
         bexp1d(i) = 0.
         xlai1d(i) = 0.
!        alb1d(i)  = 0.
         vegf1d(i) = 0.
      enddo
      if (Model%do_sfcperts) then
        if (Model%pertz0(1) > 0.) then
          z01d(:) = Model%pertz0(1) * Coupling%sfc_wts(:,1)
!          if (me == 0) print*,'Coupling%sfc_wts(:,1) min and max',minval(Coupling%sfc_wts(:,1)),maxval(Coupling%sfc_wts(:,1))
!          if (me == 0) print*,'z01d min and max ',minval(z01d),maxval(z01d)
        endif
        if (Model%pertzt(1) > 0.) then
          zt1d(:) = Model%pertzt(1) * Coupling%sfc_wts(:,2)
        endif
        if (Model%pertshc(1) > 0.) then
          bexp1d(:) = Model%pertshc(1) * Coupling%sfc_wts(:,3)
        endif
        if (Model%pertlai(1) > 0.) then
          xlai1d(:) = Model%pertlai(1) * Coupling%sfc_wts(:,4)
        endif
! --- do the albedo percentile calculation in GFS_radiation_driver instead --- !
!        if (Model%pertalb(1) > 0.) then
!          do i=1,im
!            call cdfnor(Coupling%sfc_wts(i,5),cdfz)
!            alb1d(i) = cdfz
!          enddo
!        endif
        if (Model%pertvegf(1) > 0.) then
          do i=1,im
            call cdfnor(Coupling%sfc_wts(i,6),cdfz)
            vegf1d(i) = cdfz
          enddo
        endif
      endif
#endif
!
      if (Model%do_shoc) then
        allocate (qrn(im,levs),   qsnw(im,levs), qgl(im,levs), &
                  ncpl(im,levs),  ncpi(im,levs))
        do k=1,levs
          do i=1,im
            ncpl(i,k) = 0.0
            ncpi(i,k) = 0.0
            qrn(i,k)  = 0.0
            qsnw(i,k) = 0.0
            qgl(i,k)  = 0.0
          enddo
        enddo
      endif
!
      if (imp_physics == Model%imp_physics_thompson) then
        if(Model%ltaerosol) then
          allocate(ice00(im,levs))
          allocate(liq0(im,levs))
!         allocate(nwfa2d(im))
        else
          allocate(ice00(im,levs))
        endif
      endif

      if (imp_physics == Model%imp_physics_mg) then         ! For MGB double moment microphysics
        allocate (qlcn(im,levs),      qicn(im,levs),    w_upi(im,levs),     &
                  cf_upi(im,levs),    CNV_MFD(im,levs),                     &
                  CNV_DQLDT(im,levs), clcn(im,levs),    cnv_fice(im,levs),  &
                  cnv_ndrop(im,levs), cnv_nice(im,levs))
        allocate (ncpr(im,levs), ncps(im,levs), ncgl(im,levs))
        if (.not. allocated(qrn))  allocate (qrn(im,levs))
        if (.not. allocated(qsnw)) allocate (qsnw(im,levs))
        if (.not. allocated(qgl))  allocate (qgl(im,levs))
        do k=1,levs
          do i=1,im
            qrn(i,k)  = 0.0
            qsnw(i,k) = 0.0
            qgl(i,k)  = 0.0
            ncpr(i,k) = 0.0
            ncps(i,k) = 0.0
            ncgl(i,k) = 0.0
          enddo
        enddo
!
      else
#ifdef CCPP
        ! Need to allocate those as for MG microphysics to be able to compile
        ! the CCPP hybrid code in DEBUG mode
        allocate (qlcn(im,levs),      qicn(im,levs),    w_upi(im,levs),     &
                  cf_upi(im,levs),    CNV_MFD(im,levs),                     &
                  CNV_DQLDT(im,levs), clcn(im,levs),    cnv_fice(im,levs),  &
                  cnv_ndrop(im,levs), cnv_nice(im,levs))
#else
        allocate (qlcn(1,1),    qicn(1,1),     w_upi(1,1),    cf_upi(1,1),  &
                  CNV_MFD(1,1),                CNV_DQLDT(1,1),              &
                  clcn(1,1),    cnv_fice(1,1), cnv_ndrop(1,1), cnv_nice(1,1))
#endif

#ifndef CCPP
        if (imp_physics == Model%imp_physics_gfdl) then       ! GFDL MP
          allocate (delp(im,1,levs),  dz(im,1,levs),    uin(im,1,levs),                    &
                    vin(im,1,levs),   pt(im,1,levs),    qv1(im,1,levs),   ql1(im,1,levs),  &
                    qr1(im,1,levs),   qg1(im,1,levs),   qa1(im,1,levs),   qn1(im,1,levs),  &
                    qi1(im,1,levs),   qs1(im,1,levs),   pt_dt(im,1,levs), qa_dt(im,1,levs),&
                    udt(im,1,levs),   vdt(im,1,levs),   w(im,1,levs),     qv_dt(im,1,levs),&
                    ql_dt(im,1,levs), qr_dt(im,1,levs), qi_dt(im,1,levs), qs_dt(im,1,levs),&
                    qg_dt(im,1,levs), p123(im,1,levs),  refl(im,1,levs))
          allocate (den(im,levs))
        endif
#endif
      endif

#ifdef GFS_HYDRO
      call get_prs(im, ix, levs, ntrac, Statein%tgrs, Statein%qgrs,     &
                   Model%thermodyn_id, Model%sfcpress_id,               &
                   Model%gen_coord_hybrid, Statein%prsi, Statein%prsik, &
                   Statein%prsl, Statein%prslk, Statein%phii, Statein%phil, del)
#else
!GFDL   Adjust the geopotential height hydrostatically in a way consistent with FV3 discretization
#ifdef CCPP
      if (Model%me==0) write(0,*) 'CCPP DEBUG: calling get_prs_fv3 through option B'
      ! Copy local variables from driver to appropriate interstitial variables
      !Interstitial(nt)%ix                    ! intent(in) - set in Interstitial(nt)%create()
      !Model%levs                             ! intent(in)
      !Statein%phii                           ! intent(in)
      !Statein%prsi                           ! intent(in)
      !Statein%tgrs                           ! intent(in)
      !Statein%qgrs(:,:,1)                    ! intent(in)
      Interstitial(nt)%del = del              ! intent(out)
      Interstitial(nt)%del_gz = del_gz        ! intent(out)
      !cdata_block(nb,nt)%errmsg = errmsg     ! intent(out)
      !cdata_block(nb,nt)%errflg = errflg     ! intent(out)
      !
      call ccpp_physics_run(cdata_block(nb,nt), scheme_name="get_prs_fv3", ierr=ierr)
      ! Copy back intent(inout) interstitial variables to local variables in driver
      del    = Interstitial(nt)%del       ! intent(out)
      del_gz = Interstitial(nt)%del_gz    ! intent(out)
      errmsg = trim(cdata_block(nb,nt)%errmsg)
      errflg = cdata_block(nb,nt)%errflg
      if (errflg/=0) then
        write(0,*) 'Error in call to get_prs_fv3: ' // trim(errmsg)
        stop
      end if
#else
      if (Model%me==0) write(0,*) 'CCPP DEBUG: calling non-CCPP compliant version of get_prs_fv3'
      call get_prs_fv3 (ix, levs, ntrac, Statein%phii, Statein%prsi, &
                        Statein%tgrs, Statein%qgrs, del, del_gz)
#endif
#endif
!

#ifdef CCPP
      if (Model%me==0) write(0,*) 'CCPP DEBUG: calling GFS_suite_interstitial_1 through option B'
      ! Copy local variables from driver to appropriate interstitial variables
      !Interstitial(nt)%im = im           ! intent(in) - set in Interstitial(nt)%create()
      !Model%levs                         ! intent(in)
      !Model%ntrac                        ! intent(in)
      !Model%crtrh                        ! intent(in)
      !Model%dtf                          ! intent(in)
      !Model%dtp                          ! intent(in)
      !Sfcprop%slmsk                      ! intent(in)
      !Grid%area                          ! intent(in)
      !Model%dxmin                        ! intent(in)
      !Model%dxinv                        ! intent(in)
      !Statein%pgr                        ! intent(in)
      !Interstitial(nt)%rhcbot = rhbbot   ! intent(out) - rhbbot uninitialized at this point
      !Interstitial(nt)%rhcpbl = rhpbl    ! intent(out) - rhpbl  uninitialized at this point
      !Interstitial(nt)%rhctop = rhbtop   ! intent(out) - rhbtop uninitialized at this point
      !Interstitial(nt)%frain  = frain    ! intent(out) - frain  uninitialized at this point
      !Interstitial(nt)%islmsk = islmsk    ! intent(out)
      !Interstitial(nt)%frland = frland    ! intent(out)
      !Interstitial(nt)%work1  = work1     ! intent(out)
      !Interstitial(nt)%work2  = work2     ! intent(out)
      !Diag%psurf                         ! intent(out)
      !Interstitial(nt)%dudt   = dudt      ! intent(out)
      !Interstitial(nt)%dvdt   = dvdt      ! intent(out)
      !Interstitial(nt)%dtdt   = dtdt      ! intent(out)
      !Interstitial(nt)%dtdtc  = dtdtc     ! intent(out)
      !Interstitial(nt)%dqdt   = dqdt      ! intent(out)
      !cdata_block(nb,nt)%errmsg = errmsg ! intent(out)
      !cdata_block(nb,nt)%errflg = errflg ! intent(out)
      !
      call ccpp_physics_run(cdata_block(nb,nt), scheme_name="GFS_suite_interstitial_1", ierr=ierr)
      ! Copy back intent(inout) interstitial variables to local variables in driver
      rhbbot = Interstitial(nt)%rhcbot
      rhpbl  = Interstitial(nt)%rhcpbl
      rhbtop = Interstitial(nt)%rhctop
      frain  = Interstitial(nt)%frain
      islmsk = Interstitial(nt)%islmsk
      frland = Interstitial(nt)%frland
      work1  = Interstitial(nt)%work1
      work2  = Interstitial(nt)%work2
      dudt   = Interstitial(nt)%dudt
      dvdt   = Interstitial(nt)%dvdt
      dtdt   = Interstitial(nt)%dtdt
      dtdtc  = Interstitial(nt)%dtdtc
      dqdt   = Interstitial(nt)%dqdt
      errmsg = trim(cdata_block(nb,nt)%errmsg)
      errflg = cdata_block(nb,nt)%errflg
      if (errflg/=0) then
          write(0,*) 'Error in call to GFS_suite_interstitial_1: ' // trim(errmsg)
          stop
      end if
#endif

#ifndef CCPP
      rhbbot = Model%crtrh(1)
      rhpbl  = Model%crtrh(2)
      rhbtop = Model%crtrh(3)
#endif
!
!  --- ...  frain=factor for centered difference scheme correction of rain amount.

#ifndef CCPP
      frain = dtf / dtp
#endif

      do i=1,im
#ifndef CCPP
        sigmaf(i) = max( Sfcprop%vfrac(i),0.01 )
        islmsk(i) = nint(Sfcprop%slmsk(i))

        if (islmsk(i) == 1) then
          frland(i) = 1.0
        else
          frland(i) = 0.
        endif

        if (islmsk(i) == 2) then
          if (Model%isot == 1) then
            soiltyp(i)  = 16
          else
            soiltyp(i)  = 9
          endif
          if (Model%ivegsrc == 1) then
            vegtype(i)  = 15
          elseif(Model%ivegsrc == 2) then
            vegtype(i)  = 13
          endif
          slopetyp(i) = 9
        else
          soiltyp(i)  = int( Sfcprop%stype(i)+0.5 )
          vegtype(i)  = int( Sfcprop%vtype(i)+0.5 )
          slopetyp(i) = int( Sfcprop%slope(i)+0.5 )    !! clu: slope -> slopetyp
        endif

!  --- ...  xw: transfer ice thickness & concentration from global to local variables
        zice(i) = Sfcprop%hice(i)
        cice(i) = Sfcprop%fice(i)
        tice(i) = Sfcprop%tisfc(i)
!
!GFDL   work1(i)   = (log(coslat(i) / (nlons(i)*latr)) - dxmin) * dxinv
!GFS         Moorthi thinks this should be area and not dx
!       work1(i)   = (log(Grid%dx(i)) - dxmin) * dxinv
        work1(i)   = (log(Grid%area(i)) - dxmin) * dxinv
        work1(i)   = max(0.0, min(1.0,work1(i)))
        work2(i)   = 1.0 - work1(i)
        Diag%psurf(i) = Statein%pgr(i)
        work3(i)   = Statein%prsik(i,1) / Statein%prslk(i,1)
#endif
!GFDL   tem1       = con_rerth * (con_pi+con_pi)*coslat(i)/nlons(i)
!GFDL   tem2       = con_rerth * con_pi / latr
!GFDL   garea(i)   = tem1 * tem2
        garea(i)   = Grid%area(i)

#ifndef CCPP
        tem1       = Grid%dx(i)
        tem2       = Grid%dx(i)
        dlength(i) = sqrt( tem1*tem1+tem2*tem2 )
        cldf(i)    = Model%cgwf(1)    * work1(i) + Model%cgwf(2)    * work2(i)
#endif
      enddo
!
! DH* note: this block is not yet in CCPP
      if (Model%cplflx) then
        do i=1,im
          islmsk_cice(i) = nint(Coupling%slimskin_cpl(i))
          flag_cice(i)   = (islmsk_cice(i) == 4)

          ulwsfc_cice(i) = Coupling%ulwsfcin_cpl(i)
          dusfc_cice(i)  = Coupling%dusfcin_cpl(i)
          dvsfc_cice(i)  = Coupling%dvsfcin_cpl(i)
          dtsfc_cice(i)  = Coupling%dtsfcin_cpl(i)
          dqsfc_cice(i)  = Coupling%dqsfcin_cpl(i)
          tisfc_cice(i)  = Sfcprop%tisfc(i)
          tsea_cice(i)   = Sfcprop%tsfc(i)
          fice_cice(i)   = Sfcprop%fice(i)
          hice_cice(i)   = Sfcprop%hice(i)
          if(flag_cice(i))Sfcprop%tsfc(i) = fice_cice(i)*tisfc_cice(i) + (1.0 - fice_cice(i))*tsea_cice(i)
        enddo
      endif
! *DH

#ifndef CCPP
      do k=1,levs
        do i=1,im
          dudt(i,k)  = 0.
          dvdt(i,k)  = 0.
          dtdt(i,k)  = 0.
          dtdtc(i,k) = 0.
        enddo
      enddo
      do n=1,ntrac
        do k=1,levs
          do i=1,im
            dqdt(i,k,n) = 0.
          enddo
        enddo
      enddo
#endif

!  --- ...  initialize dtdt with heating rate from dcyc2

!  --- ...  adjust mean radiation fluxes and heating rates to fit for
!           faster model time steps.
!      sw:  using cos of zenith angle as scaling factor
!      lw:  using surface air skin temperature as scaling factor

      if (Model%pre_rad) then
        call dcyc2t3_pre_rad                                                &
!  ---  inputs:
           ( Model%solhr, Model%slag, Model%sdec, Model%cdec, Grid%sinlat,  &
             Grid%coslat, Grid%xlon, Radtend%coszen, Sfcprop%tsfc,          &
             Statein%tgrs(1,1), Statein%tgrs(1,1), Coupling%sfcdsw,         &
             Coupling%sfcnsw, Coupling%sfcdlw, Radtend%htrsw, Radtend%htrlw,&
             Coupling%nirbmui, Coupling%nirdfui, Coupling%visbmui,          &
             Coupling%visdfui, Coupling%nirbmdi, Coupling%nirdfdi,          &
             Coupling%visbmdi, Coupling%visdfdi, ix, im, levs,              &
!  ---  input/output:
             dtdt,                                                          &
!  ---  outputs:
             adjsfcdsw, adjsfcnsw, adjsfcdlw, adjsfculw, xmu, xcosz,        &
             adjnirbmu, adjnirdfu, adjvisbmu, adjvisdfu,                    &
             adjnirbmd, adjnirdfd, adjvisbmd, adjvisdfd                     &
           )

      else
#ifdef CCPP
            if (Model%me==0) write(0,*) 'CCPP DEBUG: calling dcyc2t3 through option B'
            ! Copy local variables from driver to appropriate interstitial variables
            !Model%solhr                            ! intent(in)
            !Model%slag                             ! intent(in)
            !Model%sdec                             ! intent(in)
            !Model%cdec                             ! intent(in)
            !Grid%sinlat                            ! intent(in)
            !Grid%coslat                            ! intent(in)
            !Grid%xlon                              ! intent(in)
            !Radtend%coszen                         ! intent(in)
            !Sfcprop%tsfc                           ! intent(in)
            !Statein%tgrs(1,1)                      ! intent(in)
            !Radtend%tsflw                          ! intent(in)
            !Radtend%semis                          ! intent(in)
            !Coupling%sfcdsw                        ! intent(in)
            !Coupling%sfcnsw                        ! intent(in)
            !Coupling%sfcdlw                        ! intent(in)
            !Radtend%htrsw                          ! intent(in)
            !Radtend%swhc                           ! intent(in)
            !Radtend%htrlw                          ! intent(in)
            !Radtend%lwhc,&                         ! intent(in)
            !Coupling%nirbmui                       ! intent(in)
            !Coupling%nirdfui                       ! intent(in)
            !Coupling%visbmui                       ! intent(in)
            !Coupling%visdfui                       ! intent(in)
            !Coupling%nirbmdi                       ! intent(in)
            !Coupling%nirdfdi                       ! intent(in)
            !Coupling%visbmdi                       ! intent(in)
            !Coupling%visdfdi                       ! intent(in)
            !Interstitial(nt)%ix                    ! intent(in) - set in Interstitial(nt)%create()
            !Interstitial(nt)%im                    ! intent(in) - set in Interstitial(nt)%create()
            !Model%levs                             ! intent(in)
            !Model%dtf                              ! intent(in)
            Interstitial(nt)%dtdt       = dtdt      ! intent(inout)
            Interstitial(nt)%dtdtc      = dtdtc     ! intent(inout)
            !Diag%dswsfci                           ! intent(out)
            !Diag%nswsfci                           ! intent(out)
            !Diag%dlwsfci                           ! intent(out)
            !Diag%ulwsfci                           ! intent(out)
            Interstitial(nt)%xmu        = xmu       ! intent(out)
            Interstitial(nt)%xcosz      = xcosz     ! intent(out)
            Interstitial(nt)%adjnirbmu  = adjnirbmu ! intent(out)
            Interstitial(nt)%adjnirdfu  = adjnirdfu ! intent(out)
            Interstitial(nt)%adjvisbmu  = adjvisbmu ! intent(out)
            Interstitial(nt)%adjvisdfu  = adjvisdfu ! intent(out)
            Interstitial(nt)%adjnirbmd  = adjnirbmd ! intent(out)
            Interstitial(nt)%adjnirdfd  = adjnirdfd ! intent(out)
            Interstitial(nt)%adjvisbmd  = adjvisbmd ! intent(out)
            Interstitial(nt)%adjvisdfd  = adjvisdfd ! intent(out)
            !cdata_block(nb,nt)%errmsg                  ! intent(out)
            !cdata_block(nb,nt)%errflg                  ! intent(out)
            call ccpp_physics_run(cdata_block(nb,nt), scheme_name="dcyc2t3", ierr=ierr)
            ! Copy intent(inout) and intent(out) interstitial variables to local variables in driver
            dtdt      = Interstitial(nt)%dtdt
            dtdtc     = Interstitial(nt)%dtdtc
            xmu       = Interstitial(nt)%xmu
            xcosz     = Interstitial(nt)%xcosz
            adjnirbmu = Interstitial(nt)%adjnirbmu
            adjnirdfu = Interstitial(nt)%adjnirdfu
            adjvisbmu = Interstitial(nt)%adjvisbmu
            adjvisdfu = Interstitial(nt)%adjvisdfu
            adjnirbmd = Interstitial(nt)%adjnirbmd
            adjnirdfd = Interstitial(nt)%adjnirdfd
            adjvisbmd = Interstitial(nt)%adjvisbmd
            adjvisdfd = Interstitial(nt)%adjvisdfd
            errmsg    = trim(cdata_block(nb,nt)%errmsg)
            errflg    = cdata_block(nb,nt)%errflg
            if (errflg/=0) then
                write(0,*) 'Error in call to dcyc2t3: ' // trim(errmsg)
                stop
            end if
#else
        if (Model%me==0) write(0,*) 'CCPP DEBUG: calling non-CCPP compliant version of dcyc2t3'
        call dcyc2t3                                                        &
!  ---  inputs:
           ( Model%solhr, Model%slag, Model%sdec, Model%cdec, Grid%sinlat,  &
             Grid%coslat, Grid%xlon, Radtend%coszen, Sfcprop%tsfc,          &
             Statein%tgrs(1,1), Radtend%tsflw, Radtend%semis,               &
             Coupling%sfcdsw,  Coupling%sfcnsw, Coupling%sfcdlw,            &
             Radtend%htrsw,    Radtend%swhc,    Radtend%htrlw, Radtend%lwhc,&
             Coupling%nirbmui, Coupling%nirdfui, Coupling%visbmui,          &
             Coupling%visdfui, Coupling%nirbmdi, Coupling%nirdfdi,          &
             Coupling%visbmdi, Coupling%visdfdi, ix, im, levs, dtf,         &
!  ---  input/output:
             dtdt, dtdtc,                                                   &
!  ---  outputs:
             adjsfcdsw, adjsfcnsw, adjsfcdlw, adjsfculw, xmu, xcosz,        &
             adjnirbmu, adjnirdfu, adjvisbmu, adjvisdfu,                    &
             adjnirbmd, adjnirdfd, adjvisbmd, adjvisdfd                     &
           )
#endif
!
! save temp change due to radiation - need for sttp stochastic physics
!---------------------------------------------------------------------
      endif
!
      if (Model%lsidea) then                       !idea jw
        dtdt(:,:) = 0.
      endif

!  ---  convert lw fluxes for land/ocean/sea-ice models
!  note: for sw: adjsfcdsw and adjsfcnsw are zenith angle adjusted downward/net fluxes.
!        for lw: adjsfcdlw is (sfc temp adjusted) downward fluxe with no emiss effect.
!                adjsfculw is (sfc temp adjusted) upward fluxe including emiss effect.
!        one needs to be aware that that the absorbed downward lw flux (used by land/ocean
!        models as downward flux) is not the same as adjsfcdlw but a value reduced by
!        the factor of emissivity.  however, the net effects are the same when seeing
!        it either above the surface interface or below.
!
!   - flux above the interface used by atmosphere model:
!        down: adjsfcdlw;    up: adjsfculw = sfcemis*sigma*T**4 + (1-sfcemis)*adjsfcdlw
!        net = up - down = sfcemis * (sigma*T**4 - adjsfcdlw)
!   - flux below the interface used by lnd/oc/ice models:
!        down: sfcemis*adjsfcdlw;  up: sfcemis*sigma*T**4
!        net = up - down = sfcemis * (sigma*T**4 - adjsfcdlw)

!  --- ...  define the downward lw flux absorbed by ground

#ifdef CCPP
      if (Model%me==0) write(0,*) 'CCPP DEBUG: calling GFS_surface_generic_pre through option B'
      ! Copy local variables from driver to appropriate interstitial variables
      !Interstitial(nt)%im = im               ! intent(in) - set in Interstitial(nt)%create()
      !Model%levs                             ! intent(in)
      !Sfcprop%vfrac                          ! intent(in)
      Interstitial(nt)%islmsk = islmsk        ! intent(in)
      !Model%isot                             ! intent(in)
      !Model%ivegsrc                          ! intent(in)
      !Sfcprop%stype                          ! intent(in)
      !Sfcprop%vtype                          ! intent(in)
      !Sfcprop%slope                          ! intent(in)
      !Statein%prsik(:,1)                     ! intent(in)
      !Statein%prslk(:,1)                     ! intent(in)
      !Radtend%semis                          ! intent(in)
      !Diag%dlwsfci                           ! intent(in)
      !Sfcprop%tsfc                           ! intent(in)
      !Statein%phil                           ! intent(in)
      !con_g                                  ! intent(in) - physical constant in physcons.F90
      Interstitial(nt)%sigmaf    = sigmaf     ! intent(inout)
      Interstitial(nt)%soiltype  = soiltyp    ! intent(inout)
      Interstitial(nt)%vegtype   = vegtype    ! intent(inout)
      Interstitial(nt)%slopetype = slopetyp   ! intent(inout)
      Interstitial(nt)%work3     = work3      ! intent(inout)
      Interstitial(nt)%gabsbdlw  = gabsbdlw   ! intent(inout)
      Interstitial(nt)%tsurf     = tsurf      ! intent(inout)
      !Diag%zlvl                              ! intent(inout)
      !Model%do_sppt                          ! intent(in)
      !Tbd%dtdtr                              ! intent(out)
      !Tbd%drain_cpl                          ! intent(out)
      !Tbd%dsnow_cpl                          ! intent(out)
      !Tbd%rain_cpl                           ! intent(in)
      !Tbd%snow_cpl                           ! intent(in)
      !Model%do_sfcperts                      ! intent(in)
      !Model%nsfcpert                         ! intent(in)
      !Coupling%sfc_wts                       ! intent(in)
      !Model%pertz0                           ! intent(in)
      !Model%pertzt                           ! intent(in)
      !Model%pertshc                          ! intent(in)
      !Model%pertlai                          ! intent(in)
      !Model%pertvegf                         ! intent(in)
      !Interstitial(nt)%z01d                  ! intent(out) - overwritten within
      !Interstitial(nt)%zt1d                  ! intent(out) - overwritten within
      !Interstitial(nt)%bexp1d                ! intent(out) - overwritten within
      !Interstitial(nt)%xlai1d                ! intent(out) - overwritten within
      !Interstitial(nt)%vegf1d                ! intent(out) - overwritten within
      !cdata_block(nb,nt)%errmsg = errmsg     ! intent(out)
      !cdata_block(nb,nt)%errflg = errflg     ! intent(out)
      call ccpp_physics_run(cdata_block(nb,nt), scheme_name="GFS_surface_generic_pre", ierr=ierr)
      ! Copy intent(inout) and intent(out) interstitial variables to local variables in driver
      sigmaf   = Interstitial(nt)%sigmaf
      soiltyp  = Interstitial(nt)%soiltype
      vegtype  = Interstitial(nt)%vegtype
      slopetyp = Interstitial(nt)%slopetype
      work3    = Interstitial(nt)%work3
      gabsbdlw = Interstitial(nt)%gabsbdlw
      tsurf    = Interstitial(nt)%tsurf
      z01d     = Interstitial(nt)%z01d
      zt1d     = Interstitial(nt)%zt1d
      bexp1d   = Interstitial(nt)%bexp1d
      xlai1d   = Interstitial(nt)%xlai1d
      vegf1d   = Interstitial(nt)%vegf1d
      errmsg   = trim(cdata_block(nb,nt)%errmsg)
      errflg   = cdata_block(nb,nt)%errflg
      if (errflg/=0) then
          write(0,*) 'Error in call to GFS_surface_generic_pre: ' // trim(errmsg)
          stop
      end if
#endif

#ifndef CCPP
      gabsbdlw(:) = Radtend%semis(:) * adjsfcdlw(:)
#endif

      do i=1,im
        kcnv(i)   = 0
#ifndef CCPP
        kinver(i) = levs
        invrsn(i) = .false.
        tx1(i)    = 0.0
        tx2(i)    = 10.0
        ctei_r(i) = 10.0
#endif
      enddo
#ifdef CCPP
      kinver = Interstitial(nt)%kinver
#endif

#ifdef CCPP
      if (Model%me==0) write(0,*) 'CCPP DEBUG: calling GFS_suite_interstitial_2 through option B'
      ! Copy local variables from driver to appropriate interstitial variables
      !Interstitial(nt)%im = im                   ! intent(in) - set in Interstitial(nt)%create()
      !Model%levs                                 ! intent(in)
      !Model%lssav                                ! intent(in)
      !Model%ldiag3d                              ! intent(in)
      !Model%lsidea                               ! intent(in)
      !Model%cplflx                               ! intent(in)
      Interstitial(nt)%flag_cice = flag_cice      ! intent(in)
      !Model%shal_cnv                             ! intent(in)
      !Model%old_monin                            ! intent(in)
      !Model%mstrat                               ! intent(in)
      !Model%do_shoc                              ! intent(in)
      !Model%imfshalcnv                           ! intent(in)
      !Model%dtf                                  ! intent(in)
      Interstitial(nt)%xcosz = xcosz              ! intent(in)
      !Diag%dswsfci                               ! intent(in) - associated with adjsfcdsw
      !Diag%dlwsfci                               ! intent(in) - associated with adjsfcdlw
      !Statein%pgr                                ! intent(in)
      Interstitial(nt)%ulwsfc_cice = ulwsfc_cice  ! intent(in)
      !Radtend%lwhd                               ! intent(in)
      !Radtend%htrsw                              ! intent(in)
      !Radtend%htrlw                              ! intent(in)
      Interstitial(nt)%xmu = xmu                  ! intent(in)
      !Model%ctei_rm                              ! intent(in)
      Interstitial(nt)%work1 = work1              ! intent(in)
      Interstitial(nt)%work2 = work2              ! intent(in)
      !Statein%prsi                               ! intent(in)
      !Statein%tgrs                               ! intent(in)
      !Statein%prsl                               ! intent(in)
      !Statein%qgrs(:,:,1)                        ! intent(in)
      !Statein%qgrs(:,:,Model%ntcw)               ! intent(in)
      !con_cp                                     ! intent(in) - physical constant in physcons.f90
      !con_hvap                                   ! intent(in) - physical constant in physcons.f90
      !Statein%prslk                              ! intent(in)
      !Diag%suntim                                ! intent(inout)
      !Diag%ulwsfci                               ! intent(inout) - associated with adjsfculw
      !Diag%dlwsfc                                ! intent(inout)
      !Diag%ulwsfc                                ! intent(inout)
      !Diag%psmean                                ! intent(inout)
      !Diag%dt3dt(:,:,1)                          ! intent(inout)
      !Diag%dt3dt(:,:,2)                          ! intent(inout)
      !Diag%dt3dt(:,:,3)                          ! intent(inout)
      !Diag%dt3dt(:,:,4)                          ! intent(inout)
      !Diag%dt3dt(:,:,5)                          ! intent(inout)
      !Diag%dt3dt(:,:,6)                          ! intent(inout)
      Interstitial(nt)%ctei_rml = ctei_rml        ! intent(inout)
      Interstitial(nt)%ctei_r = ctei_r            ! intent(inout)
      Interstitial(nt)%kinver = kinver            ! intent(inout)
      !cdata_block(nb,nt)%errmsg = errmsg         ! intent(out)
      !cdata_block(nb,nt)%errflg = errflg         ! intent(out)
      call ccpp_physics_run(cdata_block(nb,nt), scheme_name="GFS_suite_interstitial_2", ierr=ierr)
      ! Copy intent(inout) and intent(out) interstitial variables to local variables in driver
      ctei_rml = Interstitial(nt)%ctei_rml
      ctei_r = Interstitial(nt)%ctei_r
      kinver = Interstitial(nt)%kinver
      errmsg   = trim(cdata_block(nb,nt)%errmsg)
      errflg   = cdata_block(nb,nt)%errflg
      if (errflg/=0) then
          write(0,*) 'Error in call to GFS_suite_interstitial_2: ' // trim(errmsg)
          stop
      end if
#else
      if (Model%lssav) then      !  --- ...  accumulate/save output variables

!  --- ...  sunshine duration time is defined as the length of time (in mdl output
!           interval) that solar radiation falling on a plane perpendicular to the
!           direction of the sun >= 120 w/m2

        do i=1,im
          if ( xcosz(i) >= czmin ) then   ! zenth angle > 89.994 deg
            tem1 = adjsfcdsw(i) / xcosz(i)
            if ( tem1 >= 120.0 ) then
              Diag%suntim(i) = Diag%suntim(i) + dtf
            endif
          endif
        enddo

!  --- ...  sfc lw fluxes used by atmospheric model are saved for output

        if (Model%cplflx) then
          do i=1,im
            if (flag_cice(i)) adjsfculw(i) = ulwsfc_cice(i)
          enddo
        endif
        do i=1,im
          Diag%dlwsfc(i) = Diag%dlwsfc(i) +   adjsfcdlw(i)*dtf
          Diag%ulwsfc(i) = Diag%ulwsfc(i) +   adjsfculw(i)*dtf
          Diag%psmean(i) = Diag%psmean(i) + Statein%pgr(i)*dtf        ! mean surface pressure
        enddo

        if (Model%ldiag3d) then
          if (Model%lsidea) then
            do k=1,levs
              do i=1,im
                Diag%dt3dt(i,k,1) = Diag%dt3dt(i,k,1) + Radtend%lwhd(i,k,1)*dtf
                Diag%dt3dt(i,k,2) = Diag%dt3dt(i,k,2) + Radtend%lwhd(i,k,2)*dtf
                Diag%dt3dt(i,k,3) = Diag%dt3dt(i,k,3) + Radtend%lwhd(i,k,3)*dtf
                Diag%dt3dt(i,k,4) = Diag%dt3dt(i,k,4) + Radtend%lwhd(i,k,4)*dtf
                Diag%dt3dt(i,k,5) = Diag%dt3dt(i,k,5) + Radtend%lwhd(i,k,5)*dtf
                Diag%dt3dt(i,k,6) = Diag%dt3dt(i,k,6) + Radtend%lwhd(i,k,6)*dtf
              enddo
            enddo
          else
            do k=1,levs
              do i=1,im
                Diag%dt3dt(i,k,1) = Diag%dt3dt(i,k,1) + Radtend%htrlw(i,k)*dtf
                Diag%dt3dt(i,k,2) = Diag%dt3dt(i,k,2) + Radtend%htrsw(i,k)*dtf*xmu(i)
              enddo
            enddo
          endif
        endif
      endif    ! end if_lssav_block


!    Only used for old shallow convection with mstrat=.true.

      if ((((Model%imfshalcnv == 0 .and. Model%shal_cnv) .or. Model%old_monin)        &
                                   .and. Model%mstrat)   .or. Model%do_shoc) then
        ctei_rml(:) = Model%ctei_rm(1)*work1(:) + Model%ctei_rm(2)*work2(:)
        do k=1,levs/2
          do i=1,im
            if (Statein%prsi(i,1)-Statein%prsi(i,k+1) < 0.35*Statein%prsi(i,1)       &
                .and. (.not. invrsn(i))) then
              tem = (Statein%tgrs(i,k+1) - Statein%tgrs(i,k))  &
                  / (Statein%prsl(i,k)   - Statein%prsl(i,k+1))

              if (((tem > 0.00010) .and. (tx1(i) < 0.0)) .or.  &
                  ((tem-abs(tx1(i)) > 0.0) .and. (tx2(i) < 0.0))) then
                invrsn(i) = .true.

                if (Statein%qgrs(i,k,1) > Statein%qgrs(i,k+1,1)) then
                  tem1 = Statein%tgrs(i,k+1) + hocp*max(Statein%qgrs(i,k+1,1),qmin)
                  tem2 = Statein%tgrs(i,k)   + hocp*max(Statein%qgrs(i,k,1),qmin)

                  tem1 = tem1 / Statein%prslk(i,k+1) - tem2 / Statein%prslk(i,k)

!  --- ...  (cp/l)(deltathetae)/(deltatwater) > ctei_rm -> conditon for CTEI
                  ctei_r(i) = (1.0/hocp)*tem1/(Statein%qgrs(i,k+1,1)-Statein%qgrs(i,k,1)  &
                            + Statein%qgrs(i,k+1,ntcw)-Statein%qgrs(i,k,ntcw))
                else
                  ctei_r(i) = 10
                endif

                if ( ctei_rml(i) > ctei_r(i) ) then
                  kinver(i) = k
                else
                  kinver(i) = levs
                endif
              endif

              tx2(i) = tx1(i)
              tx1(i) = tem
            endif
          enddo
        enddo
      endif
#endif

!  --- ...  lu: initialize flag_guess, flag_iter, tsurf

      do i=1,im
#ifndef CCPP
        tsurf(i)        = Sfcprop%tsfc(i)
#endif
        flag_guess(i)   = .false.
        flag_iter(i)    = .true.
        drain(i)        = 0.0
        ep1d(i)         = 0.0
        runof(i)        = 0.0
        hflx(i)         = 0.0
        evap(i)         = 0.0
        evbs(i)         = 0.0
        evcw(i)         = 0.0
        trans(i)        = 0.0
        sbsno(i)        = 0.0
        snowc(i)        = 0.0
        snohf(i)        = 0.0
#ifndef CCPP
        Diag%zlvl(i)    = Statein%phil(i,1) * onebg
#endif
        Diag%smcwlt2(i) = 0.0
        Diag%smcref2(i) = 0.0
      enddo

!  --- ...  lu: iter-loop over (sfc_diff,sfc_drv,sfc_ocean,sfc_sice)

      do iter=1,2

#include "debug_bitforbit_diagtoscreen.inc"

!  --- ...  surface exchange coefficients
!
!     if (lprnt) write(0,*)' tsfc=',Sfcprop%tsfc(ipr),' tsurf=',tsurf(ipr),iter
#ifdef CCPP
         if (Model%me==0) write(0,*) 'CCPP DEBUG: calling sfc_ex_coef through option B'
         ! Copy local variables from driver to appropriate interstitial variables
         !Interstitial(nt)%im = im              ! intent(in) - set in Interstitial(nt)%create()
         !Statein%pgr                           ! intent(in)
         !Statein%ugrs                          ! intent(in)
         !Statein%vgrs                          ! intent(in)
         !Statein%tgrs                          ! intent(in)
         !Statein%qgrs                          ! intent(in)
         !Diag%zlvl                             ! intent(in)
         !Sfcprop%snowd                         ! intent(in)
         !Sfcprop%tsfc                          ! intent(in)
         !Sfcprop%zorl                          ! intent(inout)
         Interstitial(nt)%cd = cd               ! intent(inout)
         Interstitial(nt)%cdq = cdq             ! intent(inout)
         Interstitial(nt)%rb = rb               ! intent(inout)
         !Statein%prsl(1,1)                     ! intent(in)
         Interstitial(nt)%work3 = work3         ! intent(in)
         Interstitial(nt)%islmsk = islmsk       ! intent(in)
         Interstitial(nt)%stress = stress       ! intent(inout)
         !Sfcprop%ffmm                          ! intent(inout)
         !Sfcprop%ffhh                          ! intent(inout)
         !Sfcprop%uustar                        ! intent(inout)
         Interstitial(nt)%wind = wind           ! intent(inout)
         !Tbd%phy_f2d(1,Model%num_p2d)          ! intent(in)
         Interstitial(nt)%fm10 = fm10           ! intent(inout)
         Interstitial(nt)%fh2 = fh2             ! intent(inout)
         Interstitial(nt)%sigmaf = sigmaf       ! intent(in)
         Interstitial(nt)%vegtype = vegtype     ! intent(in)
         !Sfcprop%shdmax                        ! intent(in)
         !Model%ivegsrc                         ! intent(in)
         Interstitial(nt)%z01d = z01d           ! intent(in)
         Interstitial(nt)%zt1d = zt1d           ! intent(in)
         Interstitial(nt)%tsurf = tsurf         ! intent(in)
         Interstitial(nt)%flag_iter = flag_iter ! intent(in)
         !Model%redrag                          ! intent(in)
         !cdata_block(nb,nt)%errmsg = errmsg    ! intent(out)
         !cdata_block(nb,nt)%errflg = errflg    ! intent(out)
         call ccpp_physics_run(cdata_block(nb,nt), scheme_name="sfc_ex_coef", ierr=ierr)
         ! Copy back intent(inout) interstitial variables to local variables in driver
         cd     = Interstitial(nt)%cd
         cdq    = Interstitial(nt)%cdq
         rb     = Interstitial(nt)%rb
         stress = Interstitial(nt)%stress
         wind   = Interstitial(nt)%wind
         fm10   = Interstitial(nt)%fm10
         fh2    = Interstitial(nt)%fh2
         errmsg = trim(cdata_block(nb,nt)%errmsg)
         errflg = cdata_block(nb,nt)%errflg
         if (errflg/=0) then
             write(0,*) 'Error in call to sfc_ex_coef: ' // trim(errmsg)
             stop
         end if
#else
         if (Model%me==0) write(0,*) 'CCPP DEBUG: calling non-CCPP compliant version of sfc_diff'
         call sfc_diff (im, Statein%pgr, Statein%ugrs, Statein%vgrs,      &
                       Statein%tgrs, Statein%qgrs, Diag%zlvl,             &
                       Sfcprop%snowd, Sfcprop%tsfc,  Sfcprop%zorl, cd,    &
                       cdq, rb, Statein%prsl(1,1), work3, islmsk, stress, &
                       Sfcprop%ffmm,  Sfcprop%ffhh, Sfcprop%uustar,       &
                       wind,  Tbd%phy_f2d(1,Model%num_p2d), fm10, fh2,    &
                       sigmaf, vegtype, Sfcprop%shdmax, Model%ivegsrc,    &
                       z01d, zt1d,                                        &  ! mg, sfc-perts
                       tsurf, flag_iter, Model%redrag)
#endif

!  --- ...  lu: update flag_guess

#ifdef CCPP
         if (Model%me==0) write(0,*) 'CCPP DEBUG: calling GFS_surface_loop_control_part1 through option B'
         ! Copy local variables from driver to appropriate interstitial variables
         !Interstitial(nt)%im = im                    ! intent(in) - set in Interstitial(nt)%create()
         ! Update the CCPP-internal loop counter - in the final CCPP version
         ! running a suite definition file without this GFS_physics_driver,
         ! this is done automatically using two subcycles for the surface
         ! processes.
         cdata_block(nb,nt)%loop_cnt     = iter       ! intent(in)
         Interstitial(nt)%wind           = wind       ! intent(in)
         Interstitial(nt)%flag_guess     = flag_guess ! intent(inout)
         !cdata_block(nb,nt)%errmsg = errmsg          ! intent(out)
         !cdata_block(nb,nt)%errflg = errflg          ! intent(out)
         call ccpp_physics_run(cdata_block(nb,nt), scheme_name="GFS_surface_loop_control_part1", ierr=ierr)
         ! Copy intent(inout) and intent(out) interstitial variables to local variables in driver
         flag_guess = Interstitial(nt)%flag_guess
         errmsg     = trim(cdata_block(nb,nt)%errmsg)
         errflg     = cdata_block(nb,nt)%errflg
         if (errflg/=0) then
             write(0,*) 'Error in call to GFS_surface_loop_control_part1: ' // trim(errmsg)
             stop
         end if
#else
         if (Model%me==0) write(0,*) 'CCPP DEBUG: calling non-CCPP compliant version of GFS_surface_loop_control_part1'
         do i=1,im
           if (iter == 1 .and. wind(i) < 2.0) then
             flag_guess(i) = .true.
           endif
         enddo
#endif

         if (Model%nstf_name(1) > 0) then

#ifdef CCPP
            if (Model%me==0) write(0,*) 'CCPP DEBUG: calling sfc_nst_pre through option B'
            ! Copy local variables from driver to appropriate interstitial variables
            !Interstitial(nt)%im = im                ! intent(in) - set in Interstitial(nt)%create()
            Interstitial(nt)%islmsk = islmsk         ! intent(in)
            Interstitial(nt)%tsurf = tsurf           ! intent(inout)
            Interstitial(nt)%tseal = tseal           ! intent(inout)
            !cdata_block(nb,nt)%errmsg = errmsg      ! intent(out)
            !cdata_block(nb,nt)%errflg = errflg      ! intent(out)
            call ccpp_physics_run(cdata_block(nb,nt), scheme_name="sfc_nst_pre", ierr=ierr)
            ! Copy intent(inout) and intent(out) interstitial variables to local variables in driver
            tsurf = Interstitial(nt)%tsurf
            tseal = Interstitial(nt)%tseal
            errmsg = trim(cdata_block(nb,nt)%errmsg)
            errflg = cdata_block(nb,nt)%errflg
            if (errflg/=0) then
                write(0,*) 'Error in call to sfc_nst_pre: ' // trim(errmsg)
                stop
            end if
            !
            if (Model%me==0) write(0,*) 'CCPP DEBUG: calling sfc_nst through option B'
            ! Copy local variables from driver to appropriate interstitial variables
            !Interstitial(nt)%im = im                ! intent(in) - set in Interstitial(nt)%create()
            Interstitial(nt)%islmsk = islmsk         ! intent(in)
            Interstitial(nt)%tsurf = tsurf           ! intent(inout)
            Interstitial(nt)%tseal = tseal           ! intent(inout)
            Interstitial(nt)%cd = cd                 ! intent(in)
            Interstitial(nt)%cdq = cdq               ! intent(in)
            Interstitial(nt)%work3 = work3           ! intent(in)
            Interstitial(nt)%stress = stress         ! intent(in)
            Interstitial(nt)%gabsbdlw = gabsbdlw     ! intent(in)
            !Diag%nswsfci                            ! intent(in)
            Interstitial(nt)%xcosz = xcosz           ! intent(in)
            Interstitial(nt)%flag_iter = flag_iter   ! intent(in)
            Interstitial(nt)%flag_guess = flag_guess ! intent(in)
            Interstitial(nt)%ipr = ipr               ! intent(in)
            Interstitial(nt)%qss = qss               ! intent(inout)
            Interstitial(nt)%gflx = gflx             ! intent(inout)
            Interstitial(nt)%evap = evap             ! intent(inout)
            Interstitial(nt)%hflx = hflx             ! intent(inout)
            Interstitial(nt)%ep1d = ep1d             ! intent(inout)
            !cdata_block(nb,nt)%errmsg = errmsg      ! intent(out)
            !cdata_block(nb,nt)%errflg = errflg      ! intent(out)
            call ccpp_physics_run(cdata_block(nb,nt), scheme_name="sfc_nst", ierr=ierr)
            ! Copy intent(inout) and intent(out) interstitial variables to local variables in driver
            tsurf = Interstitial(nt)%tsurf
            tseal = Interstitial(nt)%tseal
            qss = Interstitial(nt)%qss
            gflx = Interstitial(nt)%gflx
            evap = Interstitial(nt)%evap
            hflx = Interstitial(nt)%hflx
            ep1d = Interstitial(nt)%ep1d
            errmsg = trim(cdata_block(nb,nt)%errmsg)
            errflg = cdata_block(nb,nt)%errflg
            if (errflg/=0) then
                write(0,*) 'Error in call to sfc_nst: ' // trim(errmsg)
                stop
            end if
            !
            if (Model%me==0) write(0,*) 'CCPP DEBUG: calling sfc_nst_post through option B'
            ! Copy local variables from driver to appropriate interstitial variables
            !Interstitial(nt)%im = im                ! intent(in) - set in Interstitial(nt)%create()
            Interstitial(nt)%islmsk = islmsk         ! intent(in)
            Interstitial(nt)%tsurf = tsurf           ! intent(inout)
            Interstitial(nt)%dtzm = dtzm             ! intent(inout)
            !cdata_block(nb,nt)%errmsg = errmsg      ! intent(out)
            !cdata_block(nb,nt)%errflg = errflg      ! intent(out)
            call ccpp_physics_run(cdata_block(nb,nt), scheme_name="sfc_nst_post", ierr=ierr)
            ! Copy intent(inout) and intent(out) interstitial variables to local variables in driver
            tsurf = Interstitial(nt)%tsurf
            dtzm = Interstitial(nt)%dtzm
            errmsg = trim(cdata_block(nb,nt)%errmsg)
            errflg = cdata_block(nb,nt)%errflg
            if (errflg/=0) then
                write(0,*) 'Error in call to sfc_nst_post: ' // trim(errmsg)
                stop
            end if
#else
            if (Model%me==0) write(0,*) 'CCPP DEBUG: calling non-CCPP compliant version of sfc_nst'
            do i=1,im
              if (islmsk(i) == 0) then
                tem      = (Sfcprop%oro(i)-Sfcprop%oro_uf(i)) * rlapse
                tseal(i) = Sfcprop%tsfc(i) + tem
                tsurf(i) = tsurf(i)        + tem
              endif
            enddo

            call sfc_nst (im, Statein%pgr, Statein%ugrs,                     &
                          Statein%vgrs, Statein%tgrs, Statein%qgrs,          &
                          Sfcprop%tref, cd, cdq, Statein%prsl(1,1), work3,   &
                          islmsk, Grid%xlon, Grid%sinlat, stress,            &
                          Radtend%semis, gabsbdlw, adjsfcnsw, Sfcprop%tprcp, &
                          dtf, kdt, Model%solhr, xcosz,                      &
                          Tbd%phy_f2d(1,Model%num_p2d), flag_iter,           &
                          flag_guess, Model%nstf_name, lprnt, ipr,           &
!  --- Input/output
                          tseal, tsurf, Sfcprop%xt, Sfcprop%xs, Sfcprop%xu,  &
                          Sfcprop%xv, Sfcprop%xz, Sfcprop%zm, Sfcprop%xtts,  &
                          Sfcprop%xzts, Sfcprop%dt_cool, Sfcprop%z_c,        &
                          Sfcprop%c_0, Sfcprop%c_d, Sfcprop%w_0, Sfcprop%w_d,&
                          Sfcprop%d_conv, Sfcprop%ifd, Sfcprop%qrain,        &
!  ---  outputs:
                          qss, gflx, Diag%cmm, Diag%chh, evap, hflx, ep1d)

!         if (lprnt) write(0,*)' tseaz2=',tseal(ipr),' tref=', Sfcprop%tref(ipr),   &
!    &     ' dt_cool=',dt_cool(ipr),' dt_warm=',2.0*xt(ipr)/xz(ipr),       &
!    &     ' kdt=',kdt

            do i=1,im
              if ( islmsk(i) == 0 ) then
                tsurf(i) = tsurf(i) - (Sfcprop%oro(i)-Sfcprop%oro_uf(i)) * rlapse
              endif
            enddo

!  --- ...  run nsst model  ... ---

            if (Model%nstf_name(1) > 1) then
              zsea1 = 0.001*real(Model%nstf_name(4))
              zsea2 = 0.001*real(Model%nstf_name(5))
              call get_dtzm_2d (Sfcprop%xt,  Sfcprop%xz, Sfcprop%dt_cool,  &
                                Sfcprop%z_c, Sfcprop%slmsk, zsea1, zsea2,  &
                                im, 1, dtzm)
              do i=1,im
                if ( islmsk(i) == 0 ) then
                Sfcprop%tsfc(i) = max(271.2,Sfcprop%tref(i) + dtzm(i)) -   &
                                  (Sfcprop%oro(i)-Sfcprop%oro_uf(i))*rlapse
                endif
              enddo
            endif

!         if (lprnt) print *,' tseaz2=',Sfcprop%tsfc(ipr),' tref=',tref(ipr),   &
!    &    ' dt_cool=',dt_cool(ipr),' dt_warm=',dt_warm(ipr),' kdt=',kdt

#endif

         else

!  --- ...  surface energy balance over ocean
#ifdef CCPP
            if (Model%me==0) write(0,*) 'CCPP DEBUG: calling sfc_ocean through option B'
            ! Copy local variables from driver to appropriate interstitial variables
            !Interstitial(nt)%im = im                ! intent(in   ) - set in Interstitial(nt)%create()
            !Statein%pgr                             ! intent(in   )
            !Statein%ugrs                            ! intent(in   )
            !Statein%vgrs                            ! intent(in   )
            !Statein%tgrs                            ! intent(in   )
            !Statein%qgrs                            ! intent(in   )
            !Sfcprop%tsfc                            ! intent(in   )
            Interstitial(nt)%cd = cd                 ! intent(in   )
            Interstitial(nt)%cdq = cdq               ! intent(in   )
            !Statein%prsl(:,1)                       ! intent(in   )
            Interstitial(nt)%work3 = work3           ! intent(in   )
            Interstitial(nt)%islmsk = islmsk         ! intent(in   )
            !Tbd%phy_f2d(:,Model%num_p2d)            ! intent(in   )
            Interstitial(nt)%flag_iter = flag_iter   ! intent(in   )
            Interstitial(nt)%qss = qss               ! intent(inout)
            !Diag%cmm                                ! intent(inout)
            !Diag%chh                                ! intent(inout)
            Interstitial(nt)%gflx = gflx             ! intent(inout)
            Interstitial(nt)%evap = evap             ! intent(inout)
            Interstitial(nt)%hflx = hflx             ! intent(inout)
            Interstitial(nt)%ep1d = ep1d             ! intent(inout)
            !cdata_block(nb,nt)%errmsg = errmsg      ! intent(  out)
            !cdata_block(nb,nt)%errflg = errflg      ! intent(  out)
            call ccpp_physics_run(cdata_block(nb,nt), scheme_name="sfc_ocean", ierr=ierr)
            ! Copy intent(inout) and intent(out) interstitial variables to local variables in driver
            qss  = Interstitial(nt)%qss
            gflx = Interstitial(nt)%gflx
            evap = Interstitial(nt)%evap
            hflx = Interstitial(nt)%hflx
            ep1d = Interstitial(nt)%ep1d
            errmsg = trim(cdata_block(nb,nt)%errmsg)
            errflg = cdata_block(nb,nt)%errflg
            if (errflg/=0) then
                write(0,*) 'Error in call to sfc_ocean: ' // trim(errmsg)
                stop
            end if
#else
            if (Model%me==0) write(0,*) 'CCPP DEBUG: calling non-CCPP compliant version of sfc_ocean'
            call sfc_ocean                                                &
!  ---  inputs:
             (im, Statein%pgr, Statein%ugrs, Statein%vgrs, Statein%tgrs,  &
              Statein%qgrs, Sfcprop%tsfc, cd, cdq, Statein%prsl(1,1),     &
              work3, islmsk, Tbd%phy_f2d(1,Model%num_p2d), flag_iter,     &
!  ---  outputs:
               qss, Diag%cmm, Diag%chh, gflx, evap, hflx, ep1d)
#endif

         endif       ! if ( nstf_name(1) > 0 ) then

!       if (lprnt) write(0,*)' sfalb=',Radtend%sfalb(ipr),' ipr=',ipr          &
!     ,   ' weasd=',Sfcprop%weasd(ipr)             &
!     ,   ' tprcp=',Sfcprop%tprcp(ipr),' kdt=',kdt,' iter=',iter               &
!     ,' tseabefland=',Sfcprop%tsfc(ipr)

!  --- ...  surface energy balance over land
!
         if (Model%lsm == Model%lsm_noah) then                          ! noah lsm call

!     if (lprnt) write(0,*)' tsead=',tsea(ipr),' tsurf=',tsurf(ipr),iter &
!    &,' pgr=',pgr(ipr),' sfcemis=',sfcemis(ipr)
#ifdef CCPP
         if (Model%me==0) write(0,*) 'CCPP DEBUG: calling lsm_noah through option B'
         ! Copy local variables from driver to appropriate interstitial variables
         !Interstitial(nt)%im = im              ! intent(in) - set in Interstitial(nt)%create()
         !Model%lsoil                           ! intent(in)
         !Statein%pgr                           ! intent(in)
         !Statein%ugrs                          ! intent(in)
         !Statein%vgrs                          ! intent(in)
         !Statein%tgrs                          ! intent(in)
         !Statein%qgrs                          ! intent(in)
         Interstitial(nt)%soiltype= soiltyp     ! intent(in)
         Interstitial(nt)%vegtype = vegtype     ! intent(in)
         Interstitial(nt)%sigmaf = sigmaf       ! intent(in)
         !Radtend%semis                         ! intent(in)
         Interstitial(nt)%gabsbdlw = gabsbdlw   ! intent(in)
         !Diag%dswsfci                          ! intent(in)
         !Diag%nswsfci                          ! intent(in)
         !Model%dtf                             ! intent(in)
         !Sfcprop%tg3                           ! intent(in)
         Interstitial(nt)%cd = cd               ! intent(in)
         Interstitial(nt)%cdq = cdq             ! intent(in)
         !Statein%prsl(1,1)                     ! intent(in)
         Interstitial(nt)%work3 = work3         ! intent(in)
         !Diag%zlvl                             ! intent(in)
         Interstitial(nt)%islmsk = islmsk       ! intent(in)
         !Tbd%phy_f2d(1,Model%num_p2d)          ! intent(in)
         Interstitial(nt)%slopetype= slopetyp   ! intent(in)
         !Sfcprop%shdmin                        ! intent(in)
         !Sfcprop%shdmax                        ! intent(in)
         !Sfcprop%snoalb                        ! intent(in)
         !Radtend%sfalb                         ! intent(in)
         Interstitial(nt)%flag_iter = flag_iter ! intent(in)
         Interstitial(nt)%flag_guess= flag_guess! intent(in)
         !Model%lheatstrg                       ! intent(in)
         !Model%isot                            ! intent(in)
         !Model%ivegsrc                         ! intent(in)
         Interstitial(nt)%bexp1d = bexp1d       ! intent(in)
         Interstitial(nt)%xlai1d = xlai1d       ! intent(in)
         Interstitial(nt)%vegf1d = vegf1d       ! intent(in)
         !Model%pertvegf                        ! intent(in)
         !Sfcprop%weasd                         ! intent(inout)
         !Sfcprop%snowd                         ! intent(inout)
         !Sfcprop%tsfc                          ! intent(inout)
         !Sfcprop%tprcp                         ! intent(inout)
         !Sfcprop%srflag                        ! intent(inout)
         !Sfcprop%smc                           ! intent(inout)
         !Sfcprop%stc                           ! intent(inout)
         !Sfcprop%slc                           ! intent(inout)
         !Sfcprop%canopy                        ! intent(inout)
         Interstitial(nt)%trans  = trans        ! intent(inout)
         Interstitial(nt)%tsurf  = tsurf        ! intent(inout)
         !Sfcprop%zorl                          ! intent(inout)
         !Sfcprop%sncovr                        ! intent(out)
         Interstitial(nt)%qss    = qss          ! intent(out)
         Interstitial(nt)%gflx   = gflx         ! intent(out)
         Interstitial(nt)%drain  = drain        ! intent(out)
         Interstitial(nt)%evap   = evap         ! intent(out)
         Interstitial(nt)%hflx   = hflx         ! intent(out)
         Interstitial(nt)%ep1d   = ep1d         ! intent(out)
         Interstitial(nt)%runoff = runof        ! intent(out)
         !Diag%cmm                              ! intent(out)
         !Diag%chh                              ! intent(out)
         Interstitial(nt)%evbs   = evbs         ! intent(out)
         Interstitial(nt)%evcw   = evcw         ! intent(out)
         Interstitial(nt)%sbsno  = sbsno        ! intent(out)
         Interstitial(nt)%snowc  = snowc        ! intent(out)
         !Diag%soilm                            ! intent(out)
         Interstitial(nt)%snohf  = snohf        ! intent(out)
         !Diag%smcwlt2                          ! intent(out)
         !Diag%smcref2                          ! intent(out)
         !Sfcprop%wet1                          ! intent(out)
         !cdata_block(nb,nt)%errmsg = errmsg    ! intent(out)
         !cdata_block(nb,nt)%errflg = errflg    ! intent(out)
         ! Call NOAH LSM through CCPP
         call ccpp_physics_run(cdata_block(nb,nt), scheme_name="lsm_noah", ierr=ierr)
         ! Copy back intent(inout) interstitial variables to local variables in driver
         trans  = Interstitial(nt)%trans
         tsurf  = Interstitial(nt)%tsurf
         qss    = Interstitial(nt)%qss
         gflx   = Interstitial(nt)%gflx
         drain  = Interstitial(nt)%drain
         evap   = Interstitial(nt)%evap
         hflx   = Interstitial(nt)%hflx
         ep1d   = Interstitial(nt)%ep1d
         runof  = Interstitial(nt)%runoff
         evbs   = Interstitial(nt)%evbs
         evcw   = Interstitial(nt)%evcw
         sbsno  = Interstitial(nt)%sbsno
         snowc  = Interstitial(nt)%snowc
         snohf  = Interstitial(nt)%snohf
         errmsg = trim(cdata_block(nb,nt)%errmsg)
         errflg = cdata_block(nb,nt)%errflg
         if (errflg/=0) then
             write(0,*) 'Error in call to lsm_noah: ' // trim(errmsg)
             stop
         end if
#else
         if (Model%me==0) write(0,*) 'CCPP DEBUG: calling non-CCPP compliant version of sfc_drv'
         call sfc_drv                                                   &
!  ---  inputs:
           (im, lsoil, Statein%pgr, Statein%ugrs, Statein%vgrs,        &
            Statein%tgrs, Statein%qgrs, soiltyp, vegtype, sigmaf,      &
            Radtend%semis, gabsbdlw, adjsfcdsw, adjsfcnsw, dtf,        &
            Sfcprop%tg3, cd, cdq, Statein%prsl(1,1), work3, Diag%zlvl, &
            islmsk, Tbd%phy_f2d(1,Model%num_p2d), slopetyp,            &
            Sfcprop%shdmin, Sfcprop%shdmax, Sfcprop%snoalb,            &
            Radtend%sfalb, flag_iter, flag_guess, Model%lheatstrg,     &
            Model%isot,Model%ivegsrc,                                  &
            bexp1d, xlai1d, vegf1d, Model%pertvegf,                    &
!  ---  in/outs:
            Sfcprop%weasd, Sfcprop%snowd, Sfcprop%tsfc, Sfcprop%tprcp, &
            Sfcprop%srflag, Sfcprop%smc, Sfcprop%stc, Sfcprop%slc,     &
            Sfcprop%canopy, trans, tsurf, Sfcprop%zorl,                &
!  ---  outputs:
            Sfcprop%sncovr, qss, gflx, drain, evap, hflx, ep1d, runof, &
            Diag%cmm, Diag%chh, evbs, evcw, sbsno, snowc, Diag%soilm,  &
            snohf, Diag%smcwlt2, Diag%smcref2, Sfcprop%wet1)

#endif
!     if (lprnt) write(0,*)' tseae=',tsea(ipr),' tsurf=',tsurf(ipr),iter &
!    &,' phy_f2d=',phy_f2d(ipr,num_p2d)

        elseif (Model%lsm == Model%lsm_ruc) then
#ifdef CCPP  
       if (Model%me==0) write(0,*) 'CCPP DEBUG: calling lsm_ruc_run through option B' 
              !Model%me                              ! intent(in)
              !Model%kdt                             ! intent(in)
              ! Interstitial(nt)%im = im             ! intent(in)- set in Interstitial(nt)%create()
              !Model%levs                            ! intent(in)
              !Model%lsoil_lsm                       ! intent(in)
              !Model%zs                              ! intent(in)
              !Statein%ugrs(:,1,1)                   ! intent(in)
              !Statein%vgrs(:,1,1)                   ! intent(in)
              !Statein%tgrs(:,1,1)                   ! intent(in)
              !Statein%qgrs(:,1,1)                   ! intent(in)
              !Statein%qgrs(:,1,ntcw)                ! intent(in)
              Interstitial(nt)%soiltype= soiltyp     ! intent(in)
              Interstitial(nt)%vegtype = vegtype     ! intent(in)
              Interstitial(nt)%sigmaf  = sigmaf      ! intent(in)
              !Radtend%semis                         ! intent(in)
              !Diag%dlwsfci                          ! intent(in)
              !Diag%dswsfci                          ! intent(in)
              !Diag%nswsfci                          ! intent(in)
              !Model%dtf                             ! intent(in)
              !Sfcprop%tg3                           ! intent(in)
              Interstitial(nt)%cd = cd               ! intent(in)
              Interstitial(nt)%cdq = cdq             ! intent(in)
              !Statein%prsl(:,1)                     ! intent(in)
              !Diag%zlvl                             ! intent(in)
              Interstitial(nt)%islmsk = islmsk       ! intent(in)
              !Sfcprop%shdmin                        ! intent(in)
              !Sfcprop%shdmax                        ! intent(in)
              !Sfcprop%albedo                        ! intent(inout)
              !Sfcprop%snoalb                        ! intent(in)
              !Radtend%sfalb                         ! intent(inout)
              !Sfcprop%alvwf                         ! intent(in)
              !Sfcprop%alnwf                         ! intent(in)
              Interstitial(nt)%flag_iter = flag_iter ! intent(in)
              Interstitial(nt)%flag_guess= flag_guess! intent(in)
              !Model%isot                            ! intent(in)
              !Model%ivegsrc                         ! intent(in)
              !Model%fice                            ! intent(in)
              !Model%lsm_ruc                         ! intent(in)
              !Model%lsm                             ! intent(in)
              !Interstitial(nt)%con_cp    = con_cp   ! intent(in)
              !Interstitial(nt)%con_rv    = con_rv   ! intent(in)
              !Interstitial(nt)%con_rd    = con_rd   ! intent(in)
              !Interstitial(nt)%con_g     = con_g    ! intent(in)
              !Interstitial(nt)%con_pi    = con_pi   ! intent(in)
              !Interstitial(nt)%con_hvap  = con_hvap ! intent(in)
              !Interstitial(nt)%con_fvirt = con_fvirt! intent(in)
              !Sfcprop%tprcp                         ! intent(in)
              !Diag%rain                             ! intent(in)
              !Diag%rainc                            ! intent(in)
              !Diag%snow                             ! intent(in)
              !Diag%graupel                          ! intent(in)
              !Sfcprop%srflag                        ! intent(in)
              !Sfcprop%sr                            ! intent(in)
              !Sfcprop%smc                           ! intent(inout)
              !Sfcprop%stc                           ! intent(inout)
              !Sfcprop%slc                           ! intent(inout)
              !Sfcprop%weasd                         ! intent(inout)
              !Sfcprop%snowd                         ! intent(inout)
              !Sfcprop%tsfc                          ! intent(inout)
              !Sfcprop%smois                         ! intent(inout)
              !Sfcprop%tslb                          ! intent(inout)
              !Sfcprop%sh2o                          ! intent(inout)
              !Sfcprop%flag_frsoil                   ! intent(inout)
              !Sfcprop%keepsmfr                      ! intent(inout)
              !Sfcprop%canopy                        ! intent(inout)
              !Sfcprop%tsnow                         ! intent(inout)
              !Sfcprop%zorl                          ! intent(inout)
              !Sfcprop%clw_surf                      ! intent(inout)
              !Sfcprop%qwv_surf                      ! intent(inout)
              !Sfcprop%cndm_surf                     ! intent(out)
              !Sfcprop%sncovr                        ! intent(out)
              !Sfcprop%rhofr                         ! intent(out)
              Interstitial(nt)%trans  = trans        ! intent(inout)
              Interstitial(nt)%tsurf  = tsurf        ! intent(inout)
              !Sfcprop%tisfc(i)                      ! intent(inout)
              Interstitial(nt)%snowc  = snowc        ! intent(inout)
              Interstitial(nt)%qss    = qss          ! intent(out)
              Interstitial(nt)%gflx   = gflx         ! intent(out)
              Interstitial(nt)%drain  = drain        ! intent(out)
              Interstitial(nt)%evap   = evap         ! intent(out)
              Interstitial(nt)%hflx   = hflx         ! intent(out)
              Interstitial(nt)%runoff = runof        ! intent(out)
              Interstitial(nt)%evbs   = evbs         ! intent(out)
              Interstitial(nt)%evcw   = evcw         ! intent(out)
              Interstitial(nt)%sbsno  = sbsno        ! intent(out)
              Interstitial(nt)%snohf  = snohf        ! intent(out)
              !Diag%runoff                           ! intent(inout)
              !Diag%srunoff                          ! intent(inout)
              !Diag%soilm                            ! intent(out)
              !Sfcprop%wet1                          ! intent(out)
              !Diag%acsnow                           ! intent(out)
              !Diag%snowfallac                       ! intent(out)
              !cdata_block(nb,nt)%errmsg = errmsg    ! intent(out)
              !cdata_block(nb,nt)%errflg = errflg    ! intent(out)

              call ccpp_physics_run(cdata_block(nb,nt), scheme_name="lsm_ruc", ierr=ierr)

              ! Copy back intent(inout) interstitial variables to local variables in driver
              trans  = Interstitial(nt)%trans
              tsurf  = Interstitial(nt)%tsurf
              qss    = Interstitial(nt)%qss
              gflx   = Interstitial(nt)%gflx
              drain  = Interstitial(nt)%drain
              runof  = Interstitial(nt)%runoff
              evap   = Interstitial(nt)%evap
              hflx   = Interstitial(nt)%hflx
              evbs   = Interstitial(nt)%evbs
              evcw   = Interstitial(nt)%evcw
              sbsno  = Interstitial(nt)%sbsno
              snowc  = Interstitial(nt)%snowc
              snohf  = Interstitial(nt)%snohf
              errmsg = trim(cdata_block(nb,nt)%errmsg)
              errflg = cdata_block(nb,nt)%errflg
         if (errflg/=0) then
             write(0,*) 'Error in call to lsm_ruc: ' // trim(errmsg)
             stop
         end if
#else
        write (0,*) 'RUC LSM is available only in CCPP'
        stop
#endif
         end if

!       if (lprnt) write(0,*)' tseabeficemodel =',Sfcprop%tsfc(ipr),' me=',me   &
!    &,   ' kdt=',kdt

      if (Model%lsm == Model%lsm_noah) then  ! DH* Not for RUC - just for Noah LSM? Should we implement a namelist option for this?
!  --- ...  surface energy balance over seaice

         if (Model%cplflx) then
           do i=1,im
             if (flag_cice(i)) then
                islmsk (i) = islmsk_cice(i)
             endif
           enddo
         endif

#ifdef CCPP
         if (Model%me==0) write(0,*) 'CCPP DEBUG: calling sfc_sice through option B'
         ! Copy local variables from driver to appropriate interstitial variables
         !Interstitial(nt)%im = im              ! intent(in   ) - set in Interstitial(nt)%create()
         !Model%lsoil                           ! intent(in   )
         !Statein%pgr                           ! intent(in   )
         !Statein%ugrs                          ! intent(in   )
         !Statein%vgrs                          ! intent(in   )
         !Statein%tgrs                          ! intent(in   )
         !Statein%qgrs                          ! intent(in   )
         !Model%dtf                             ! intent(in   )
         !Radtend%semis                         ! intent(in   )
         Interstitial(nt)%gabsbdlw = gabsbdlw   ! intent(in   )
         !Diag%nswsfci                          ! intent(in   ) - associated with adjsfcnsw
         !Diag%dswsfci                          ! intent(in   ) - associated with adjsfcdsw
         !Sfcprop%srflag                        ! intent(in   )
         Interstitial(nt)%cd = cd               ! intent(in   )
         Interstitial(nt)%cdq = cdq             ! intent(in   )
         !Statein%prsl(:,1)                     ! intent(in   )
         Interstitial(nt)%work3 = work3         ! intent(in   )
         Interstitial(nt)%islmsk = islmsk       ! intent(in   )
         !Tbd%phy_f2d(1,Model%num_p2d)          ! intent(in   )
         Interstitial(nt)%flag_iter = flag_iter ! intent(in   )
         !Model%mom4ice                         ! intent(in   )
         !Model%lsm                             ! intent(in   )
         !Model%lprnt                           ! intent(in   )
         !Interstitial(nt)%ipr                  ! intent(in   ) - set in Interstitial(nt)%create()
         !Sfcprop%hice                          ! intent(inout)
         !Sfcprop%fice                          ! intent(inout)
         !Sfcprop%tisfc                         ! intent(inout)
         !Sfcprop%weasd                         ! intent(inout)
         !Sfcprop%tsfc                          ! intent(inout)
         !Sfcprop%tprcp                         ! intent(inout)
         !Sfcprop%stc                           ! intent(inout)
         Interstitial(nt)%ep1d = ep1d           ! intent(inout)
         !Sfcprop%snowd                         ! intent(inout)
         Interstitial(nt)%qss = qss             ! intent(inout)
         Interstitial(nt)%snowmt = snowmt       ! intent(inout)
         Interstitial(nt)%gflx = gflx           ! intent(inout)
         !Diag%cmm                              ! intent(inout)
         !Diag%chh                              ! intent(inout)
         Interstitial(nt)%evap = evap           ! intent(inout)
         Interstitial(nt)%hflx = hflx           ! intent(inout)
         !cdata_block(nb,nt)%errmsg = errmsg    ! intent(out)
         !cdata_block(nb,nt)%errflg = errflg    ! intent(out)
         call ccpp_physics_run(cdata_block(nb,nt), scheme_name="sfc_sice", ierr=ierr)
         ! Copy back intent(inout) interstitial variables to local variables in driver
         ep1d   = Interstitial(nt)%ep1d
         qss    = Interstitial(nt)%qss
         snowmt = Interstitial(nt)%snowmt
         gflx   = Interstitial(nt)%gflx
         evap   = Interstitial(nt)%evap
         hflx   = Interstitial(nt)%hflx
         errmsg = trim(cdata_block(nb,nt)%errmsg)
         errflg = cdata_block(nb,nt)%errflg
         if (errflg/=0) then
             write(0,*) 'Error in call to sfc_sice: ' // trim(errmsg)
             stop
         end if
#else
         if (Model%me==0) write(0,*) 'CCPP DEBUG: calling non-CCPP compliant version of sfc_sice'
         call sfc_sice                                                     &
!  ---  inputs:
              (im, lsoil, Statein%pgr, Statein%ugrs, Statein%vgrs,         &
               Statein%tgrs, Statein%qgrs, dtf, Radtend%semis, gabsbdlw,   &
               adjsfcnsw, adjsfcdsw, Sfcprop%srflag, cd, cdq,              &
               Statein%prsl(1,1), work3, islmsk,                           &
               Tbd%phy_f2d(1,Model%num_p2d), flag_iter, Model%mom4ice,     &
               Model%lsm, lprnt, ipr,                                      &
!  ---  input/outputs:
               zice, cice, tice, Sfcprop%weasd, Sfcprop%tsfc,              &
               Sfcprop%tprcp, Sfcprop%stc, ep1d,                           &
!  ---  outputs:
               Sfcprop%snowd, qss, snowmt, gflx, Diag%cmm, Diag%chh, evap, &
               hflx)
#endif

        if (Model%cplflx .or. Model%cplchm) then
          do i=1,im
            if (flag_cice(i)) then
               islmsk(i) = nint(Sfcprop%slmsk(i))
            endif
          enddo

           call sfc_cice                                                   &
!  ---  inputs:
              (im, Statein%ugrs, Statein%vgrs, Statein%tgrs, Statein%qgrs, &
               cd, cdq, Statein%prsl(1,1), work3, islmsk_cice,             &
               Tbd%phy_f2d(1,Model%num_p2d),flag_iter, dqsfc_cice,         &
               dtsfc_cice,                                                 &
!  ---  outputs:
             qss, Diag%cmm, Diag%chh, evap, hflx)
         endif
       endif ! lsm == lsm_noah for sea ice

!  --- ...  lu: update flag_iter and flag_guess
#ifdef CCPP
         if (Model%me==0) write(0,*) 'CCPP DEBUG: calling GFS_surface_loop_control_part2 through option B'
         ! Copy local variables from driver to appropriate interstitial variables
         !Interstitial(nt)%im = im                    ! intent(in) - set in Interstitial(nt)%create()
         cdata_block(nb,nt)%loop_cnt = iter           ! intent(in)
         Interstitial(nt)%wind       = wind           ! intent(in)
         Interstitial(nt)%flag_guess = flag_guess     ! intent(inout)
         Interstitial(nt)%flag_iter  = flag_iter      ! intent(inout)
         !Sfcprop%slmsk                               ! intent(in)
         !Model%nstf_name(1)                          ! intent(in)
         !cdata_block(nb,nt)%errmsg = errmsg          ! intent(out)
         !cdata_block(nb,nt)%errflg = errflg          ! intent(out)
         call ccpp_physics_run(cdata_block(nb,nt), scheme_name="GFS_surface_loop_control_part2", ierr=ierr)
         ! Copy intent(inout) and intent(out) interstitial variables to local variables in driver
         flag_guess = Interstitial(nt)%flag_guess
         flag_iter  = Interstitial(nt)%flag_iter
         errmsg     = trim(cdata_block(nb,nt)%errmsg)
         errflg     = cdata_block(nb,nt)%errflg
         if (errflg/=0) then
             write(0,*) 'Error in call to GFS_surface_loop_control_part2: ' // trim(errmsg)
             stop
         end if
#else
         if (Model%me==0) write(0,*) 'CCPP DEBUG: calling non-CCPP compliant version of GFS_surface_loop_control_part2'
         do i=1,im
           flag_iter(i)  = .false.
           flag_guess(i) = .false.

           if (iter == 1 .and. wind(i) < 2.0) then
             if (islmsk(i) == 1 .or. (islmsk(i) == 0 .and.           &
                                      Model%nstf_name(1) > 0)) then
               flag_iter(i) = .true.
             endif
           endif

!         if(islmsk(i) == 1 .and. iter == 1) then
!           if (wind(i) < 2.0) flag_iter(i) = .true.
!         elseif (islmsk(i) == 0 .and. iter == 1                        &
!    &                           .and. nstf_name(1) > 0) then
!           if (wind(i) < 2.0) flag_iter(i) = .true.
!         endif
         enddo
#endif

      enddo   ! end iter_loop

#include "debug_bitforbit_diagtoscreen.inc"

#ifdef CCPP
      if (Model%me==0) write(0,*) 'CCPP DEBUG: calling dcyc2t3_post through option B'
      ! Copy local variables from driver to appropriate interstitial variables
      !Interstitial(nt)%im                 ! intent(in) - set in Interstitial(nt)%create()
      !Diag%dswsfci                        ! intent(in)
      !Diag%nswsfci                        ! intent(in)
      !Diag%uswsfci                        ! intent(out)
      !cdata_block(nb,nt)%errmsg = errmsg  ! intent(out)
      !cdata_block(nb,nt)%errflg = errflg  ! intent(out)
      call ccpp_physics_run(cdata_block(nb,nt), scheme_name="dcyc2t3_post", ierr=ierr)
      ! Copy intent(inout) and intent(out) interstitial variables to local variables in driver
      errmsg     = trim(cdata_block(nb,nt)%errmsg)
      errflg     = cdata_block(nb,nt)%errflg
      if (errflg/=0) then
          write(0,*) 'Error in call to dcyc2t3_post: ' // trim(errmsg)
          stop
      end if
#endif

#ifndef CCPP
      do i=1,im
        Diag%epi(i)     = ep1d(i)
        Diag%dlwsfci(i) = adjsfcdlw(i)
        Diag%ulwsfci(i) = adjsfculw(i)
        Diag%uswsfci(i) = adjsfcdsw(i) - adjsfcnsw(i)
        Diag%dswsfci(i) = adjsfcdsw(i)
        Diag%gfluxi(i)  = gflx(i)
        Diag%t1(i)      = Statein%tgrs(i,1)
        Diag%q1(i)      = Statein%qgrs(i,1,1)
        Diag%u1(i)      = Statein%ugrs(i,1)
        Diag%v1(i)      = Statein%vgrs(i,1)
      enddo
#endif
!  --- ...  update near surface fields

#ifdef CCPP
      if (Model%me==0) write(0,*) 'CCPP DEBUG: calling sfc_diag through option B'
      ! Copy local variables from driver to appropriate interstitial variables
      !Interstitial(nt)%im                 ! intent(in ) - set in Interstitial(nt)%create()
      !con_g, con_cp, con_eps, con_epsm1   ! intent(in ) - physical constants from physcons.F90
      !Statein%pgr                         ! intent(in )
      !Stateout%gu0(:,1)                   ! intent(in )
      !Stateout%gv0(:,1)                   ! intent(in )
      !Stateout%gt0(:,1)                   ! intent(in )
      !Stateout%gq0(:,1,1)                 ! intent(in )
      !Sfcprop%tsfc                        ! intent(in )
      Interstitial(nt)%qss = qss           ! intent(in )
      !Sfcprop%f10m                        ! intent(out)
      !Diag%u10m                           ! intent(out)
      !Diag%v10m                           ! intent(out)
      !Sfcprop%t2m                         ! intent(out)
      !Sfcprop%q2m                         ! intent(out)
      Interstitial(nt)%work3 = work3       ! intent(in )
      Interstitial(nt)%evap  = evap        ! intent(in )
      !Sfcprop%ffmm                        ! intent(in )
      !Sfcprop%ffhh                        ! intent(in )
      Interstitial(nt)%fm10 = fm10         ! intent(in )
      Interstitial(nt)%fh2 = fh2           ! intent(in )
      !cdata_block(nb,nt)%errmsg = errmsg  ! intent(out)
      !cdata_block(nb,nt)%errflg = errflg  ! intent(out)
      call ccpp_physics_run(cdata_block(nb,nt), scheme_name="sfc_diag", ierr=ierr)
      ! Copy intent(inout) and intent(out) interstitial variables to local variables in driver
      errmsg     = trim(cdata_block(nb,nt)%errmsg)
      errflg     = cdata_block(nb,nt)%errflg
      if (errflg/=0) then
          write(0,*) 'Error in call to sfc_diag: ' // trim(errmsg)
          stop
      end if
#else
      if (Model%me==0) write(0,*) 'CCPP DEBUG: calling non-CCPP compliant version of sfc_diag'
      call sfc_diag (im, Statein%pgr, Statein%ugrs, Statein%vgrs,     &
                     Statein%tgrs, Statein%qgrs, Sfcprop%tsfc, qss,   &
                     Sfcprop%f10m, Diag%u10m,    Diag%v10m,           &
                     Sfcprop%t2m,  Sfcprop%q2m,  work3, evap,         &
                     Sfcprop%ffmm, Sfcprop%ffhh, fm10, fh2)
#endif

#ifndef CCPP
      Tbd%phy_f2d(:,Model%num_p2d) = 0.0
#endif

#ifdef CCPP
      if (Model%me==0) write(0,*) 'CCPP DEBUG: calling sfc_diag_post through option B'
      ! Copy local variables from driver to appropriate interstitial variables
      !Interstitial(nt)%im = im               ! intent(in) - set in Interstitial(nt)%create()
      !Model%lssav                            ! intent(in)
      !con_eps                                ! intent(in) - physical constant in physcons
      !con_epsm1                              ! intent(in) - physical constant in physcons
      !Sfcprop%t2m                            ! intent(in)
      !Sfcprop%q2m                            ! intent(in)
      !Sfcprop%u10m                           ! intent(in)
      !Sfcprop%v10m                           ! intent(in)
      !Diag%tmpmin                            ! intent(inout)
      !Diag%tmpmax                            ! intent(inout)
      !Diag%spfhmin                           ! intent(inout)
      !Diag%spfhmax                           ! intent(inout)
      !Diag%wind10mmax                        ! intent(inout)
      !Diag%u10mmax                           ! intent(inout)
      !Diag%v10mmax                           ! intent(inout)
      !Diag%dpt2m                             ! intent(inout)
      !cdata_block(nb,nt)%errmsg = errmsg     ! intent(out)
      !cdata_block(nb,nt)%errflg = errflg     ! intent(out)
      call ccpp_physics_run(cdata_block(nb,nt), scheme_name="sfc_diag_post", ierr=ierr)
      ! Copy intent(inout) and intent(out) interstitial variables to local variables in driver
      errmsg     = trim(cdata_block(nb,nt)%errmsg)
      errflg     = cdata_block(nb,nt)%errflg
      if (errflg/=0) then
          write(0,*) 'Error in call to sfc_diag_post: ' // trim(errmsg)
          stop
      end if
      !
      if (Model%me==0) write(0,*) 'CCPP DEBUG: calling GFS_surface_generic_post through option B'
      ! Copy local variables from driver to appropriate interstitial variables
      !Interstitial(nt)%im = im               ! intent(in) - set in Interstitial(nt)%create()
      !Model%cplflx                           ! intent(in)
      !Model%lssav                            ! intent(in)
      !Model%cplwav                           ! intent(in)
      Interstitial(nt)%islmsk = islmsk        ! intent(in)
      !Model%dtf                              ! intent(in)
      Interstitial(nt)%ep1d = ep1d            ! intent(in)
      Interstitial(nt)%gflx = gflx            ! intent(in)
      !Statein%tgrs(:,1)                      ! intent(in)
      !Statein%qgrs(:,1)                      ! intent(in)
      !Statein%ugrs(:,1)                      ! intent(in)
      !Statein%vgrs(:,1)                      ! intent(in)
      !Diag%dlwsfci                           ! intent(in) - associated with adjsfcdlw
      !Diag%dswsfci                           ! intent(in) - associated with adjsfcdsw
      Interstitial(nt)%adjnirbmd = adjnirbmd  ! intent(in)
      Interstitial(nt)%adjnirdfd = adjnirdfd  ! intent(in)
      Interstitial(nt)%adjvisbmd = adjvisbmd  ! intent(in)
      Interstitial(nt)%adjvisdfd = adjvisdfd  ! intent(in)
      !Diag%ulwsfci                           ! intent(in) - associated with adjsfculw
      Interstitial(nt)%adjnirbmu = adjnirbmu  ! intent(in)
      Interstitial(nt)%adjnirdfu = adjnirdfu  ! intent(in)
      Interstitial(nt)%adjvisbmu = adjvisbmu  ! intent(in)
      Interstitial(nt)%adjvisdfu = adjvisdfu  ! intent(in)
      !Sfcprop%t2m                            ! intent(in)
      !Sfcprop%q2m                            ! intent(in)
      !Sfcprop%u10m                           ! intent(in)
      !Sfcprop%v10m                           ! intent(in)
      !Sfcprop%tsfc                           ! intent(in)
      !Statein%pgr                            ! intent(in)
      Interstitial(nt)%xcosz = xcosz          ! intent(in)
      Interstitial(nt)%evbs = evbs            ! intent(in)
      Interstitial(nt)%evcw = evcw            ! intent(in)
      Interstitial(nt)%trans = trans          ! intent(in)
      Interstitial(nt)%sbsno = sbsno          ! intent(in)
      Interstitial(nt)%snowc = snowc          ! intent(in)
      Interstitial(nt)%snohf = snohf          ! intent(in)
      !Diag%epi                               ! intent(inout)
      !Diag%gfluxi                            ! intent(inout)
      !Diag%t1                                ! intent(inout)
      !Diag%q1                                ! intent(inout)
      !Diag%u1                                ! intent(inout)
      !Diag%v1                                ! intent(inout)
      !Coupling%dlwsfci_cpl                   ! intent(inout)
      !Coupling%dswsfci_cpl                   ! intent(inout)
      !Coupling%dlwsfc_cpl                    ! intent(inout)
      !Coupling%dswsfc_cpl                    ! intent(inout)
      !Coupling%dnirbmi_cpl                   ! intent(inout)
      !Coupling%dnirdfi_cpl                   ! intent(inout)
      !Coupling%dvisbmi_cpl                   ! intent(inout)
      !Coupling%dvisdfi_cpl                   ! intent(inout)
      !Coupling%dnirbm_cpl                    ! intent(inout)
      !Coupling%dnirdf_cpl                    ! intent(inout)
      !Coupling%dvisbm_cpl                    ! intent(inout)
      !Coupling%dvisdf_cpl                    ! intent(inout)
      !Coupling%nlwsfci_cpl                   ! intent(inout)
      !Coupling%nlwsfc_cpl                    ! intent(inout)
      !Coupling%t2mi_cpl                      ! intent(inout)
      !Coupling%q2mi_cpl                      ! intent(inout)
      !Coupling%u10mi_cpl                     ! intent(inout)
      !Coupling%v10mi_cpl                     ! intent(inout)
      !Coupling%tsfci_cpl                     ! intent(inout)
      !Coupling%psurfi_cpl                    ! intent(inout)
      !Coupling%nnirbmi_cpl                   ! intent(inout)
      !Coupling%nnirdfi_cpl                   ! intent(inout)
      !Coupling%nvisbmi_cpl                   ! intent(inout)
      !Coupling%nvisdfi_cpl                   ! intent(inout)
      !Coupling%swsfci_cpl                    ! intent(inout)
      !Coupling%nswsfc_cpl                    ! intent(inout)
      !Coupling%nnirbm_cpl                    ! intent(inout)
      !Coupling%nnirdf_cpl                    ! intent(inout)
      !Coupling%nvisbm_cpl                    ! intent(inout)
      !Coupling%nvisdf_cpl                    ! intent(inout)
      !Diag%gflux                             ! intent(inout)
      !Diag%evbsa                             ! intent(inout)
      !Diag%evcwa                             ! intent(inout)
      !Diag%transa                            ! intent(inout)
      !Diag%sbsnoa                            ! intent(inout)
      !Diag%snowca                            ! intent(inout)
      !Diag%snohfa                            ! intent(inout)
      !Diag%ep                                ! intent(inout)
      !Diag%runoff                            ! intent(inout)
      !Diag%srunoff                           ! intent(inout)
      Interstitial(nt)%runoff = runof         ! intent(in)
      Interstitial(nt)%drain  = drain         ! intent(in)
      !cdata_block(nb,nt)%errmsg = errmsg     ! intent(out)
      !cdata_block(nb,nt)%errflg = errflg     ! intent(out)
      call ccpp_physics_run(cdata_block(nb,nt), scheme_name="GFS_surface_generic_post", ierr=ierr)
      ! Copy intent(inout) and intent(out) interstitial variables to local variables in driver
      errmsg     = trim(cdata_block(nb,nt)%errmsg)
      errflg     = cdata_block(nb,nt)%errflg
      if (errflg/=0) then
          write(0,*) 'Error in call to GFS_surface_generic_post: ' // trim(errmsg)
          stop
      end if
#else
      if (Model%cplflx .or. Model%cplwav) then
        do i=1,im
          Coupling%u10mi_cpl   (i) = Diag%u10m(i)
          Coupling%v10mi_cpl   (i) = Diag%v10m(i)
        enddo
      endif 

      if (Model%cplflx) then
        do i=1,im
          Coupling%dlwsfci_cpl (i) = adjsfcdlw(i)
          Coupling%dswsfci_cpl (i) = adjsfcdsw(i)
          Coupling%dlwsfc_cpl  (i) = Coupling%dlwsfc_cpl(i) + adjsfcdlw(i)*dtf
          Coupling%dswsfc_cpl  (i) = Coupling%dswsfc_cpl(i) + adjsfcdsw(i)*dtf
          Coupling%dnirbmi_cpl (i) = adjnirbmd(i)
          Coupling%dnirdfi_cpl (i) = adjnirdfd(i)
          Coupling%dvisbmi_cpl (i) = adjvisbmd(i)
          Coupling%dvisdfi_cpl (i) = adjvisdfd(i)
          Coupling%dnirbm_cpl  (i) = Coupling%dnirbm_cpl(i) + adjnirbmd(i)*dtf
          Coupling%dnirdf_cpl  (i) = Coupling%dnirdf_cpl(i) + adjnirdfd(i)*dtf
          Coupling%dvisbm_cpl  (i) = Coupling%dvisbm_cpl(i) + adjvisbmd(i)*dtf
          Coupling%dvisdf_cpl  (i) = Coupling%dvisdf_cpl(i) + adjvisdfd(i)*dtf
          Coupling%nlwsfci_cpl (i) = adjsfcdlw(i)           - adjsfculw(i)
          Coupling%nlwsfc_cpl  (i) = Coupling%nlwsfc_cpl(i) + Coupling%nlwsfci_cpl(i)*dtf
          Coupling%t2mi_cpl    (i) = Sfcprop%t2m(i)
          Coupling%q2mi_cpl    (i) = Sfcprop%q2m(i)
          Coupling%tsfci_cpl   (i) = Sfcprop%tsfc(i)
          Coupling%psurfi_cpl  (i) = Statein%pgr(i)
        enddo

!  ---  estimate mean albedo for ocean point without ice cover and apply
!       them to net SW heat fluxes

        do i=1,im
          if (islmsk(i) /= 1) then  ! not a land point
!  ---  compute open water albedo
            xcosz_loc = max( 0.0, min( 1.0, xcosz(i) ))
            ocalnirdf_cpl(i) = 0.06
            ocalnirbm_cpl(i) = max(albdf, 0.026/(xcosz_loc**1.7+0.065)  &
     &                       + 0.15 * (xcosz_loc-0.1) * (xcosz_loc-0.5) &
     &                       * (xcosz_loc-1.0))
            ocalvisdf_cpl(i) = 0.06
            ocalvisbm_cpl(i) = ocalnirbm_cpl(i)

            Coupling%nnirbmi_cpl(i) = adjnirbmd(i)-adjnirbmd(i)*ocalnirbm_cpl(i)
            Coupling%nnirdfi_cpl(i) = adjnirdfd(i)-adjnirdfd(i)*ocalnirdf_cpl(i)
            Coupling%nvisbmi_cpl(i) = adjvisbmd(i)-adjvisbmd(i)*ocalvisbm_cpl(i)
            Coupling%nvisdfi_cpl(i) = adjvisdfd(i)-adjvisdfd(i)*ocalvisdf_cpl(i)
          else
            Coupling%nnirbmi_cpl(i) = adjnirbmd(i) - adjnirbmu(i)
            Coupling%nnirdfi_cpl(i) = adjnirdfd(i) - adjnirdfu(i)
            Coupling%nvisbmi_cpl(i) = adjvisbmd(i) - adjvisbmu(i)
            Coupling%nvisdfi_cpl(i) = adjvisdfd(i) - adjvisdfu(i)
          endif
          Coupling%nswsfci_cpl(i) = Coupling%nnirbmi_cpl(i) + Coupling%nnirdfi_cpl(i)   &
                                  + Coupling%nvisbmi_cpl(i) + Coupling%nvisdfi_cpl(i)
          Coupling%nswsfc_cpl(i)  = Coupling%nswsfc_cpl(i)  + Coupling%nswsfci_cpl(i)*dtf
          Coupling%nnirbm_cpl(i)  = Coupling%nnirbm_cpl(i)  + Coupling%nnirbmi_cpl(i)*dtf
          Coupling%nnirdf_cpl(i)  = Coupling%nnirdf_cpl(i)  + Coupling%nnirdfi_cpl(i)*dtf
          Coupling%nvisbm_cpl(i)  = Coupling%nvisbm_cpl(i)  + Coupling%nvisbmi_cpl(i)*dtf
          Coupling%nvisdf_cpl(i)  = Coupling%nvisdf_cpl(i)  + Coupling%nvisdfi_cpl(i)*dtf
        enddo
      endif

      if (Model%lssav) then
        do i=1,im
          Diag%gflux(i)   = Diag%gflux(i)  + gflx(i)  * dtf
          Diag%evbsa(i)   = Diag%evbsa(i)  + evbs(i)  * dtf
          Diag%evcwa(i)   = Diag%evcwa(i)  + evcw(i)  * dtf
          Diag%transa(i)  = Diag%transa(i) + trans(i) * dtf
          Diag%sbsnoa(i)  = Diag%sbsnoa(i) + sbsno(i) * dtf
          Diag%snowca(i)  = Diag%snowca(i) + snowc(i) * dtf
          Diag%snohfa(i)  = Diag%snohfa(i) + snohf(i) * dtf
          Diag%ep(i)      = Diag%ep(i)     + ep1d(i)  * dtf

          Diag%tmpmax(i)  = max(Diag%tmpmax(i), Sfcprop%t2m(i))
          Diag%tmpmin(i)  = min(Diag%tmpmin(i), Sfcprop%t2m(i))

          Diag%spfhmax(i) = max(Diag%spfhmax(i), Sfcprop%q2m(i))
          Diag%spfhmin(i) = min(Diag%spfhmin(i), Sfcprop%q2m(i))
        enddo

        do i=1, im
! find max wind speed then decompose
           tem = sqrt(Diag%u10m(i)*Diag%u10m(i) + Diag%v10m(i)*Diag%v10m(i))
           if (tem > Diag%wind10mmax(i)) then
              Diag%wind10mmax(i) = tem
              Diag%u10mmax(i)    = Diag%u10m(i)
              Diag%v10mmax(i)    = Diag%v10m(i)
           endif

! Compute dew point, first using vapor pressure
           tem = max(Statein%pgr(i) * Sfcprop%q2m(i) / ( con_eps - con_epsm1 * Sfcprop%q2m(i)), 1.e-8)
           Diag%dpt2m(i) = 243.5 / ( ( 17.67 / log(tem/611.2) ) - 1.) + 273.14
        enddo

      endif
#endif
!!!!!!!!!!!!!!!!!Commented by Moorthi on July 18, 2012 !!!!!!!!!!!!!!!!!!!
!     do i=1,im
!  --- ...  compute coefficient of evaporation in evapc
!
!       if (evapc(i) > 1.0e0) evapc(i) = 1.0e0
!  --- ...  over snow cover or ice or sea, coef of evap =1.0e0
!       if (weasd(i) > 0.0 .or. slmsk(i) /= 1.0) evapc(i) = 1.0e0
!     enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  --- ...  Boundary Layer and Free atmospheic turbulence parameterization

!     if (lprnt) write(0,*)' tsea3=',Sfcprop%tsfc(ipr),' slmsk=',Sfcprop%slmsk(ipr)     &
!    &, ' kdt=',kdt,' evap=',evap(ipr)
!     if (lprnt)  write(0,*)' dtdtb=',(dtdt(ipr,k),k=1,15)

!     do i=1,im
!       if (islmsk(i) == 0) then
!         oro_land(i) = 0.0
!       else
!         oro_land(i) = oro(i)
!       endif
!     enddo

!     write(0,*)' before monin clstp=',clstp,' kdt=',kdt,' lat=',lat
!  if (lprnt) write(0,*)'befmonshoc=',Statein%tgrs(ipr,:)
!  if (lprnt) write(0,*)'befmonshocdtdt=',dtdt(ipr,1:10)
!  if (lprnt) write(0,*)'befmonshoctkh=',Tbd%phy_f3d(ipr,1:10,ntot3d-1)
!  if (lprnt) write(0,*)'befmonshochflx=',hflx(ipr),' tsea=',Sfcprop%tsfc(ipr),&
!      ' evap=',evap(ipr)
!  if (lprnt) write(0,*)'befmonshoctke=',Statein%qgrs(ipr,:,ntke)

      if (nvdiff == ntrac .or. Model%do_ysu .or. Model%shinhong) then
!
#ifdef CCPP
        if (Model%me==0) write(0,*) 'CCPP DEBUG: calling PBL_generic_pre through option B'
        ! Copy local variables from driver to appropriate interstitial variables
        !Interstitial(nt)%im = im             ! intent(in) - set in Interstitial(nt)%create()
        !Model%levs                           ! intent(in)
        !Model%nvdiff                         ! intent(in)
        !Model%ntrac                          ! intent(in)
        !Model%imp_physics                    ! intent(in)
        !Model%imp_physics_gfdl               ! intent(in)
        !Model%imp_physics_thompson           ! intent(in)
        !Model%imp_physics_wsm6               ! intent(in)
        !Model%ltaerosol                      ! intent(in)
        !Statein%qgrs                         ! intent(in)
        !Statein%qgrs(:,:,1)                  ! intent(in)
        !Statein%qgrs(:,:,ntcw)               ! intent(in)
        !Statein%qgrs(:,:,ntiw)               ! intent(in)
        !Statein%qgrs(:,:,ntoz)               ! intent(in)
        !Statein%qgrs(:,:,ntlnc)              ! intent(in)
        !Statein%qgrs(:,:,ntinc)              ! intent(in)
        !Statein%qgrs(:,:,ntwa)               ! intent(in)
        !Statein%qgrs(:,:,ntia)               ! intent(in)
        !Statein%qgrs(:,:,ntrw)               ! intent(in)
        !Statein%qgrs(:,:,ntsw)               ! intent(in)
        !Statein%qgrs(:,:,ntgl)               ! intent(in)
        !Interstitial(nt)%vdftra = vdftra     ! intent(inout)
        !cdata_block(nb,nt)%errmsg = errmsg   ! intent(out)
        !cdata_block(nb,nt)%errflg = errflg   ! intent(out)
        call ccpp_physics_run(cdata_block(nb,nt), scheme_name="GFS_PBL_generic_pre", ierr=ierr)
        ! Copy back intent(inout) interstitial variables to local variables in driver
        vdftra = Interstitial(nt)%vdftra
        errmsg = trim(cdata_block(nb,nt)%errmsg)
        errflg = cdata_block(nb,nt)%errflg
        if (errflg/=0) then
            write(0,*) 'Error in call to GFS_PBL_generic_pre: ' // trim(errmsg)
            stop
        end if
        ! DH* update local copy of ntiwx with the value set in
        ! Interstitial(nt)%create(); only needed for hybrid mode
        ntiwx = Interstitial(nt)%ntiwx
        ! *DH
#else
        ntiwx = 0
#endif

        if (Model%do_shoc) then
#ifdef CCPP
          if (Model%me==0) write(0,*) 'CCPP DEBUG: calling moninshoc through option B'
          ! Copy local variables from driver to appropriate interstitial variables
          !Interstitial(nt)%ix = ix             ! intent(in) - set in Interstitial(nt)%create()
          !Interstitial(nt)%im = im             ! intent(in) - set in Interstitial(nt)%create()
          !Model%levs                           ! intent(in)
          Interstitial(nt)%nvdiff = nvdiff      ! intent(in)
          !Model%ntcw                           ! intent(in)
          Interstitial(nt)%nncl = nncl          ! intent(in)
          Interstitial(nt)%dvdt = dvdt          ! intent(inout)
          Interstitial(nt)%dudt = dudt          ! intent(inout)
          Interstitial(nt)%dtdt = dtdt          ! intent(inout)
          !Interstitial(nt)%dvdftra             ! intent(inout)
          !Statein%ugrs                         ! intent(in)
          !Statein%vgrs                         ! intent(in)
          !Statein%tgrs                         ! intent(in)
          !Interstitial(nt)%vdftra              ! intent(in)
          !Tbd%phy_f3d(1,1,ntot3d-1)            ! intent(in)
          Interstitial(nt)%prnum = prnum        ! intent(inout)
          !Model%ntke                           ! intent(in)
          !Statein%prsik                        ! intent(in)
          Interstitial(nt)%rb = rb              ! intent(in)
          !Sfcprop%zorl                         ! intent(in)
          !Diag%u10m                            ! intent(in)
          !Diag%v10m                            ! intent(in)
          !Sfcprop%ffmm                         ! intent(in)
          !Sfcprop%ffhh                         ! intent(in)
          !Sfcprop%tsfc                         ! intent(in)
          Interstitial(nt)%hflx = hflx          ! intent(in)
          Interstitial(nt)%evap = evap          ! intent(in)
          Interstitial(nt)%stress = stress      ! intent(in)
          Interstitial(nt)%wind = wind          ! intent(in)
          Interstitial(nt)%kpbl = kpbl          ! intent(out)
          !Statein%prsi                         ! intent(in)
          Interstitial(nt)%del = del            ! intent(in)
          !Statein%prsl                         ! intent(in)
          !Statein%prslk                        ! intent(in)
          !Statein%phii                         ! intent(in)
          !Statein%phil                         ! intent(in)
          !Model%dtp                            ! intent(in)
          Interstitial(nt)%dusfc1 = dusfc1      ! intent(out)
          Interstitial(nt)%dvsfc1 = dvsfc1      ! intent(out)
          Interstitial(nt)%dtsfc1 = dtsfc1      ! intent(out)
          Interstitial(nt)%dqsfc1 = dqsfc1      ! intent(out)
          Interstitial(nt)%dkt  = dkt           ! intent(out)
          !Diag%hpbl                            ! intent(out)
          Interstitial(nt)%kinver = kinver      ! intent(in)
          !Model%xkzm_m                         ! intent(in)
          !Model%xkzm_h                         ! intent(in)
          !Model%xkzm_s                         ! intent(in)
          !Model%lprnt                          ! intent(in)
          Interstitial(nt)%ipr = ipr            ! intent(in)
          !Model%me                             ! intent(in)
          !con_g, con_rd, con_cp, con_hvap, con_fvirt ! intent(in) - physical constant in physcons.F90
          !cdata_block(nb,nt)%errmsg = errmsg   ! intent(out)
          !cdata_block(nb,nt)%errflg = errflg   ! intent(out)
          call ccpp_physics_run(cdata_block(nb,nt), scheme_name="moninshoc", ierr=ierr)
          ! Copy back intent(inout) interstitial variables to local variables in driver
          dvdt   = Interstitial(nt)%dvdt
          dudt   = Interstitial(nt)%dudt
          dtdt   = Interstitial(nt)%dtdt
          prnum  = Interstitial(nt)%prnum
          kpbl   = Interstitial(nt)%kpbl
          dusfc1 = Interstitial(nt)%dusfc1
          dvsfc1 = Interstitial(nt)%dvsfc1
          dtsfc1 = Interstitial(nt)%dtsfc1
          dqsfc1 = Interstitial(nt)%dqsfc1
          dkt    = Interstitial(nt)%dkt
          errmsg = trim(cdata_block(nb,nt)%errmsg)
          errflg = cdata_block(nb,nt)%errflg
          if (errflg/=0) then
              write(0,*) 'Error in call to moninshoc: ' // trim(errmsg)
              stop
          end if
#else
          if (Model%me==0) write(0,*) 'CCPP DEBUG: calling non-CCPP compliant version of moninshoc'
          call moninshoc(ix, im, levs, nvdiff, ntcw, nncl, dvdt, dudt, dtdt, dqdt, &
                         Statein%ugrs, Statein%vgrs, Statein%tgrs, Statein%qgrs,   &
                         Tbd%phy_f3d(1,1,ntot3d-1), prnum, ntke,                   &
                         Statein%prsik(1,1), rb, Sfcprop%zorl, Diag%u10m,          &
                         Diag%v10m, Sfcprop%ffmm, Sfcprop%ffhh, Sfcprop%tsfc, hflx,&
                         evap, stress, wind, kpbl, Statein%prsi, del, Statein%prsl,&
                         Statein%prslk, Statein%phii, Statein%phil, dtp, dusfc1,   &
                         dvsfc1, dtsfc1, dqsfc1, dkt, Diag%hpbl, kinver,           &
                         Model%xkzm_m, Model%xkzm_h, Model%xkzm_s, lprnt, ipr, me)
#endif
!  if (lprnt) write(0,*)'aftmonshoc=',Statein%tgrs(ipr,:)
!  if (lprnt) write(0,*)'aftmonshoctke=',Statein%qgrs(ipr,:,ntke)
!  if (lprnt) write(0,*)'aftmonshocdtdt=',dtdt(ipr,1:10)
        else
          if (Model%satmedmf) then
#ifdef CCPP
            if (Model%me==0) write(0,*) 'CCPP DEBUG: calling satmedmfvdif through option B'
            ! Copy local variables from driver to appropriate interstitial variables
            !Interstitial(nt)%im = im             ! intent(in) - set in Interstitial(nt)%create()
            !Interstitial(nt)%ix = ix             ! intent(in) - set in Interstitial(nt)%create()
            !Model%levs                           ! intent(in)
            Interstitial(nt)%nvdiff = nvdiff      ! intent(in)
            !Model%ntcw                           ! intent(in)
            !Interstitial(nt)%ntkev = ntiwx       ! intent(in) - set in Interstitial(nt)%create()
            Interstitial(nt)%ntkev = ntkev        ! intent(in)
            !con_g
            !con_rd
            !con_cp
            !con_rv
            !con_hvap
            !con_hfus
            !con_fvirt
            !con_eps
            !con_epsm1
            Interstitial(nt)%dvdt = dvdt          ! intent(inout)
            Interstitial(nt)%dudt = dudt          ! intent(inout)
            Interstitial(nt)%dtdt = dtdt          ! intent(inout)
            !Interstitial(nt)%dvdftra             ! intent(inout)
            !Statein%ugrs                         ! intent(in)
            !Statein%vgrs                         ! intent(in)
            !Statein%tgrs                         ! intent(in)
            !Interstitial(nt)%vdftra              ! intent(in)
            !Radtend%htrsw                        ! intent(in)
            !Radtend%htrlw                        ! intent(in)
            Interstitial(nt)%xmu = xmu            ! intent(in)
            !Grid%area = garea                    ! intent(in)
            !Statein%prsik(1,1)                   ! intent(in)
            Interstitial(nt)%rb = rb              ! intent(in)
            !Sfcprop%zorl                         ! intent(in)
            !Diag%u10m                            ! intent(in)
            !Diag%v10m                            ! intent(in)
            !Sfcprop%ffmm                         ! intent(in)
            !Sfcprop%ffhh                         ! intent(in)
            !Sfcprop%tsfc                         ! intent(in)
            Interstitial(nt)%hflx = hflx          ! intent(in)
            Interstitial(nt)%evap = evap          ! intent(in)
            Interstitial(nt)%stress = stress      ! intent(in)
            Interstitial(nt)%wind = wind          ! intent(in)
            Interstitial(nt)%kpbl = kpbl          ! intent(out)
            !Statein%prsi                         ! intent(in)
            Interstitial(nt)%del = del            ! intent(in)
            !Statein%prsl                         ! intent(in)
            !Statein%prslk                        ! intent(in)
            !Statein%phii                         ! intent(in)
            !Statein%phil                         ! intent(in)
            !Model%dtp                            ! intent(in)
            !Model%dspheat                        ! intent(in)
            Interstitial(nt)%dusfc1 = dusfc1      ! intent(out)
            Interstitial(nt)%dvsfc1 = dvsfc1      ! intent(out)
            Interstitial(nt)%dtsfc1 = dtsfc1      ! intent(out)
            Interstitial(nt)%dqsfc1 = dqsfc1      ! intent(out)
            !Diag%hpbl                            ! intent(out)
            Interstitial(nt)%kinver = kinver      ! intent(in)
            !Model%xkzm_m                         ! intent(in)
            !Model%xkzm_h                         ! intent(in)
            !Model%xkzm_s                         ! intent(in)
            !cdata_block(nb,nt)%errmsg = errmsg   ! intent(out)
            !cdata_block(nb,nt)%errflg = errflg   ! intent(out)
            call ccpp_physics_run(cdata_block(nb,nt), scheme_name="satmedmfvdif", ierr=ierr)
            ! Copy back intent(inout) interstitial variables to local variables in driver
            dvdt   = Interstitial(nt)%dvdt
            dudt   = Interstitial(nt)%dudt
            dtdt   = Interstitial(nt)%dtdt
            !dvdftra = Interstitial(nt)%dvdftra
            kpbl   = Interstitial(nt)%kpbl
            dusfc1 = Interstitial(nt)%dusfc1
            dvsfc1 = Interstitial(nt)%dvsfc1
            dtsfc1 = Interstitial(nt)%dtsfc1
            dqsfc1 = Interstitial(nt)%dqsfc1
            errmsg = trim(cdata_block(nb,nt)%errmsg)
            errflg = cdata_block(nb,nt)%errflg
            if (errflg/=0) then
                write(0,*) 'Error in call to satmedmfvdif: ' // trim(errmsg)
                stop
            end if

#else
              if (Model%me==0) write(0,*) 'CCPP DEBUG: calling non-CCPP satmedmfvdif'
              call satmedmfvdif(ix, im, levs, nvdiff, ntcw, ntiw, ntke,             &
                       dvdt, dudt, dtdt, dqdt,                                      &
                       Statein%ugrs, Statein%vgrs, Statein%tgrs, Statein%qgrs,      &
                       Radtend%htrsw, Radtend%htrlw, xmu, garea,                    &
                       Statein%prsik(1,1), rb, Sfcprop%zorl, Diag%u10m, Diag%v10m,  &
                       Sfcprop%ffmm, Sfcprop%ffhh, Sfcprop%tsfc, hflx, evap,        &
                       stress, wind, kpbl, Statein%prsi, del, Statein%prsl,         &
                       Statein%prslk, Statein%phii, Statein%phil, dtp,              &
                       Model%dspheat, dusfc1, dvsfc1, dtsfc1, dqsfc1, Diag%hpbl,    &
                       kinver, Model%xkzm_m, Model%xkzm_h, Model%xkzm_s)
#endif
          elseif (Model%hybedmf) then
#ifdef CCPP
              if (Model%me==0) write(0,*) 'CCPP DEBUG: calling hedmf through option B'
              ! Copy local variables from driver to appropriate interstitial variables
              !Interstitial(nt)%im = im             ! intent(in) - set in Interstitial(nt)%create()
              !Interstitial(nt)%ix = ix             ! intent(in) - set in Interstitial(nt)%create()
              !Model%levs                           ! intent(in)
              Interstitial(nt)%nvdiff = nvdiff      ! intent(in)
              !Model%ntcw                           ! intent(in)
              Interstitial(nt)%dvdt = dvdt          ! intent(inout)
              Interstitial(nt)%dudt = dudt          ! intent(inout)
              Interstitial(nt)%dtdt = dtdt          ! intent(inout)
              !Interstitial(nt)%dvdftra             ! intent(inout)
              !Statein%ugrs                         ! intent(in)
              !Statein%vgrs                         ! intent(in)
              !Statein%tgrs                         ! intent(in)
              !Interstitial(nt)%vdftra              ! intent(in)
              !Radtend%htrsw                        ! intent(in)
              !Radtend%htrlw                        ! intent(in)
              Interstitial(nt)%xmu = xmu            ! intent(in)
              !Statein%prsik(1,1)                   ! intent(in)
              Interstitial(nt)%rb = rb              ! intent(in)
              !Sfcprop%zorl                         ! intent(in)
              !Diag%u10m                            ! intent(in)
              !Diag%v10m                            ! intent(in)
              !Sfcprop%ffmm                         ! intent(in)
              !Sfcprop%ffhh                         ! intent(in)
              !Sfcprop%tsfc                         ! intent(in)
              Interstitial(nt)%hflx = hflx          ! intent(in)
              Interstitial(nt)%evap = evap          ! intent(in)
              Interstitial(nt)%stress = stress      ! intent(in)
              Interstitial(nt)%wind = wind          ! intent(in)
              Interstitial(nt)%kpbl = kpbl          ! intent(out)
              !Statein%prsi                         ! intent(in)
              Interstitial(nt)%del = del            ! intent(in)
              !Statein%prsl                         ! intent(in)
              !Statein%prslk                        ! intent(in)
              !Statein%phii                         ! intent(in)
              !Statein%phil                         ! intent(in)
              !Model%dtp                            ! intent(in)
              !Model%dspheat                        ! intent(in)
              Interstitial(nt)%dusfc1 = dusfc1      ! intent(out)
              Interstitial(nt)%dvsfc1 = dvsfc1      ! intent(out)
              Interstitial(nt)%dtsfc1 = dtsfc1      ! intent(out)
              Interstitial(nt)%dqsfc1 = dqsfc1      ! intent(out)
              !Diag%hpbl                            ! intent(out)
              Interstitial(nt)%gamt = gamt          ! intent(inout)
              Interstitial(nt)%gamq = gamq          ! intent(inout)
              Interstitial(nt)%dkt  = dkt           ! intent(out)
              Interstitial(nt)%kinver = kinver      ! intent(in)
              !Model%xkzm_m                         ! intent(in)
              !Model%xkzm_h                         ! intent(in)
              !Model%xkzm_s                         ! intent(in)
              !Model%lprnt                          ! intent(in)
              Interstitial(nt)%ipr = ipr            ! intent(in)
              !Model%xkzminv                        ! intent(in)
              !Model%moninq_fac                     ! intent(in)
              !cdata_block(nb,nt)%errmsg = errmsg   ! intent(out)
              !cdata_block(nb,nt)%errflg = errflg   ! intent(out)
              call ccpp_physics_run(cdata_block(nb,nt), scheme_name="hedmf", ierr=ierr)
              ! Copy back intent(inout) interstitial variables to local variables in driver
              dvdt   = Interstitial(nt)%dvdt
              dudt   = Interstitial(nt)%dudt
              dtdt   = Interstitial(nt)%dtdt
              !dvdftra = Interstitial(nt)%dvdftra
              kpbl   = Interstitial(nt)%kpbl
              dusfc1 = Interstitial(nt)%dusfc1
              dvsfc1 = Interstitial(nt)%dvsfc1
              dtsfc1 = Interstitial(nt)%dtsfc1
              dqsfc1 = Interstitial(nt)%dqsfc1
              gamt   = Interstitial(nt)%gamt
              gamq   = Interstitial(nt)%gamq
              dkt    = Interstitial(nt)%dkt
              errmsg = trim(cdata_block(nb,nt)%errmsg)
              errflg = cdata_block(nb,nt)%errflg
              if (errflg/=0) then
                  write(0,*) 'Error in call to hedmf: ' // trim(errmsg)
                  stop
              end if
#else
              if (Model%me==0) write(0,*) 'CCPP DEBUG: calling non-CCPP compliant version of moninedmf'
              call moninedmf(ix, im, levs, nvdiff, ntcw, dvdt, dudt, dtdt, dqdt,    &
                           Statein%ugrs, Statein%vgrs, Statein%tgrs, Statein%qgrs,  &
                           Radtend%htrsw, Radtend%htrlw, xmu, Statein%prsik(1,1),   &
                           rb, Sfcprop%zorl, Diag%u10m, Diag%v10m, Sfcprop%ffmm,    &
                           Sfcprop%ffhh, Sfcprop%tsfc, qss, hflx, evap, stress,     &
                           wind, kpbl, Statein%prsi, del, Statein%prsl,             &
                           Statein%prslk, Statein%phii, Statein%phil, dtp,          &
                           Model%dspheat, dusfc1, dvsfc1, dtsfc1, dqsfc1, Diag%hpbl,&
                           gamt, gamq, dkt, kinver, Model%xkzm_m, Model%xkzm_h,     &
                           Model%xkzm_s, lprnt, ipr,                                &
                           Model%xkzminv, Model%moninq_fac)
!     if (lprnt)  write(0,*)' dtdtm=',(dtdt(ipr,k),k=1,15)
!     if (lprnt)  write(0,*)' dqdtm=',(dqdt(ipr,k,1),k=1,15)
#endif
          elseif (Model%do_ysu) then
#ifdef CCPP
            if (Model%me==0) write(0,*) 'CCPP DEBUG: calling ysuvdif through option B'
            ! Copy local variables from driver to appropriate interstitial variables
            !Interstitial(nt)%im = im             ! intent(in   ) - set in Interstitial(nt)%create()
            !Interstitial(nt)%ix = ix             ! intent(in   ) - set in Interstitial(nt)%create()
            !Model%levs                           ! intent(in   )
            !Statein%ugrs                         ! intent(in   )
            !Statein%vgrs                         ! intent(in   )
            !Statein%tgrs                         ! intent(in   )
            !Statein%qgrs                         ! intent(in   )
            !Statein%prsl                         ! intent(in   )
            !Statein%prsi                         ! intent(in   )
            !Statein%prslk                        ! intent(in   )
            Interstitial(nt)%dvdt     = dvdt      ! intent(inout)
            Interstitial(nt)%dudt     = dudt      ! intent(inout)
            Interstitial(nt)%dtdt     = dtdt      ! intent(inout)
            Interstitial(nt)%dqdt     = dqdt      ! intent(inout)
            !Tbd%htswc                            ! intent(in   )
            !Tbd%htlwc                            ! intent(in   )
            Interstitial(nt)%xmu      = xmu       ! intent(in   )
            !Model%ntrac                          ! intent(in   )
            Interstitial(nt)%nvdiff   = nvdiff    ! intent(in   )
            !Model%ntcw                           ! intent(in   )
            !Model%ntiw                           ! intent(in   )
            !Statein%phii                         ! intent(in   )
            !Statein%phil                         ! intent(in   )
            !Statein%pgr                          ! intent(in   )
            !Sfcprop%zorl                         ! intent(in   )
            Interstitial(nt)%stress   = stress    ! intent(in   )
            !Intdiag%hpbl                         ! intent(out  )
            !Sfcprop%ffmm                         ! intent(in   )
            !Sfcprop%ffhh                         ! intent(in   )
            Interstitial(nt)%islmsk   = islmsk    ! intent(in   )
            Interstitial(nt)%hflx     = hflx      ! intent(in   )
            Interstitial(nt)%evap     = evap      ! intent(in   )
            Interstitial(nt)%wind     = wind      ! intent(in   )
            Interstitial(nt)%rb       = rb        ! intent(in   )
            !con_g                                ! intent(in   )
            !con_rd                               ! intent(in   )
            !con_cp                               ! intent(in   )
            !con_rv                               ! intent(in   )
            !con_fvirt                            ! intent(in   )
            !con_eps                              ! intent(in   )
            !con_hvap                             ! intent(in   )
            Interstitial(nt)%dusfc1   = dusfc1    ! intent(out  )
            Interstitial(nt)%dvsfc1   = dvsfc1    ! intent(out  )
            Interstitial(nt)%dtsfc1   = dtsfc1    ! intent(out  )
            Interstitial(nt)%dqsfc1   = dqsfc1    ! intent(out  )
            !Model%dtp                            ! intent(in   )
            Interstitial(nt)%kpbl     = kpbl      ! intent(out  )
            !Intdiag%u10m                         ! intent(in   )
            !Intdiag%v10m                         ! intent(in   )
            !cdata_block(nb,nt)%errmsg = errmsg   ! intent(out  )
            !cdata_block(nb,nt)%errflg = errflg   ! intent(out  )
            call ccpp_physics_run(cdata_block(nb,nt), scheme_name="ysuvdif", ierr=ierr)
            ! Copy back intent(inout) interstitial variables to local variables in driver
            dvdt   = Interstitial(nt)%dvdt
            dudt   = Interstitial(nt)%dudt
            dtdt   = Interstitial(nt)%dtdt
            dqdt   = Interstitial(nt)%dqdt
            dusfc1 = Interstitial(nt)%dusfc1
            dvsfc1 = Interstitial(nt)%dvsfc1
            dtsfc1 = Interstitial(nt)%dtsfc1
            dqsfc1 = Interstitial(nt)%dqsfc1
            kpbl   = Interstitial(nt)%kpbl
            errmsg = trim(cdata_block(nb,nt)%errmsg)
            errflg = cdata_block(nb,nt)%errflg
            if (errflg/=0) then
                write(0,*) 'Error in call to ysuvdif: ' // trim(errmsg)
                stop
            end if
#else
            if (Model%me==0) then
                write(0,*) 'Error, ysuvdif only available through CCPP'
                stop
            end if
#endif
          elseif (Model%shinhong) then
#ifdef CCPP
            if (Model%me==0) write(0,*) 'CCPP DEBUG: calling shinhongvdif through option B'
            ! Copy local variables from driver to appropriate interstitial variables
            !Interstitial(nt)%im = im             ! intent(in   ) - set in Interstitial(nt)%create()
            !Interstitial(nt)%ix = ix             ! intent(in   ) - set in Interstitial(nt)%create()
            !Model%levs                           ! intent(in   )
            !Statein%ugrs                         ! intent(in   )
            !Statein%vgrs                         ! intent(in   )
            !Statein%tgrs                         ! intent(in   )
            !Statein%qgrs                         ! intent(in   )
            !Statein%prsl                         ! intent(in   )
            !Statein%prsi                         ! intent(in   )
            !Statein%prslk                        ! intent(in   )
            Interstitial(nt)%dvdt     = dvdt      ! intent(inout)
            Interstitial(nt)%dudt     = dudt      ! intent(inout)
            Interstitial(nt)%dtdt     = dtdt      ! intent(inout)
            Interstitial(nt)%dqdt     = dqdt      ! intent(inout)
            !Model%ntrac                          ! intent(in   )
            Interstitial(nt)%nvdiff   = nvdiff    ! intent(in   )
            !Model%ntcw                           ! intent(in   )
            !Model%ntiw                           ! intent(in   )
            !Statein%phii                         ! intent(in   )
            !Statein%phil                         ! intent(in   )
            !Statein%pgr                          ! intent(in   )
            !Sfcprop%zorl                         ! intent(in   )
            Interstitial(nt)%stress   = stress    ! intent(in   )
            !Intdiag%hpbl                         ! intent(out  )
            !Sfcprop%ffmm                         ! intent(in   )
            !Sfcprop%ffhh                         ! intent(in   )
            Interstitial(nt)%islmsk   = islmsk    ! intent(in   )
            Interstitial(nt)%hflx     = hflx      ! intent(in   )
            Interstitial(nt)%evap     = evap      ! intent(in   )
            Interstitial(nt)%wind     = wind      ! intent(in   )
            Interstitial(nt)%rb       = rb        ! intent(in   )
            !con_g                                ! intent(in   )
            !con_rd                               ! intent(in   )
            !con_cp                               ! intent(in   )
            !con_rv                               ! intent(in   )
            !con_fvirt                            ! intent(in   )
            !con_eps                              ! intent(in   )
            !con_hvap                             ! intent(in   )
            Interstitial(nt)%dusfc1   = dusfc1    ! intent(out  )
            Interstitial(nt)%dvsfc1   = dvsfc1    ! intent(out  )
            Interstitial(nt)%dtsfc1   = dtsfc1    ! intent(out  )
            Interstitial(nt)%dqsfc1   = dqsfc1    ! intent(out  )
            !Model%dtp                            ! intent(in   )
            Interstitial(nt)%kpbl     = kpbl      ! intent(out  )
            !Intdiag%u10m                         ! intent(in   )
            !Intdiag%v10m                         ! intent(in   )
            !Grid%dx                              ! intent(in   )
            !cdata_block(nb,nt)%errmsg = errmsg   ! intent(out  )
            !cdata_block(nb,nt)%errflg = errflg   ! intent(out  )
            call ccpp_physics_run(cdata_block(nb,nt), scheme_name="shinhongvdif", ierr=ierr)
            ! Copy back intent(inout) interstitial variables to local variables in driver
            dvdt   = Interstitial(nt)%dvdt
            dudt   = Interstitial(nt)%dudt
            dtdt   = Interstitial(nt)%dtdt
            dqdt   = Interstitial(nt)%dqdt
            dusfc1 = Interstitial(nt)%dusfc1
            dvsfc1 = Interstitial(nt)%dvsfc1
            dtsfc1 = Interstitial(nt)%dtsfc1
            dqsfc1 = Interstitial(nt)%dqsfc1
            kpbl   = Interstitial(nt)%kpbl
            errmsg = trim(cdata_block(nb,nt)%errmsg)
            errflg = cdata_block(nb,nt)%errflg
            if (errflg/=0) then
                write(0,*) 'Error in call to shinhongvdif: ' // trim(errmsg)
                stop
            end if
#else
            if (Model%me==0) then
                write(0,*) 'Error, shinhongvdif only available through CCPP'
                stop
            end if
#endif
          elseif (.not. Model%old_monin) then
            if (Model%me==0) write(0,*) 'CCPP DEBUG: calling non-CCPP compliant version of moninq'
            call moninq(ix, im, levs, nvdiff, ntcw, dvdt, dudt, dtdt, dqdt,         &
                        Statein%ugrs, Statein%vgrs, Statein%tgrs, Statein%qgrs,     &
                        Radtend%htrsw, Radtend%htrlw, xmu, Statein%prsik(1,1), rb,  &
                        Sfcprop%ffmm, Sfcprop%ffhh, Sfcprop%tsfc, qss, hflx, evap,  &
                        stress, wind, kpbl, Statein%prsi, del, Statein%prsl,        &
                        Statein%prslk, Statein%phii, Statein%phil, dtp,             &
                        Model%dspheat, dusfc1, dvsfc1, dtsfc1, dqsfc1, Diag%hpbl,   &
                        gamt, gamq, dkt, kinver, Model%xkzm_m, Model%xkzm_h,        &
                        Model%xkzm_s, lprnt, ipr,                                   &
                        Model%xkzminv, Model%moninq_fac, Model%rbcr)
          else
            if (Model%mstrat) then
              if (Model%me==0) write(0,*) 'CCPP DEBUG: calling non-CCPP compliant version of moninp1'
              call moninp1(ix, im, levs, nvdiff, dvdt, dudt, dtdt, dqdt,            &
                           Statein%ugrs, Statein%vgrs, Statein%tgrs, Statein%qgrs,  &
                           Statein%prsik(1,1), rb, Sfcprop%ffmm, Sfcprop%ffhh,      &
                           Sfcprop%tsfc, qss, hflx, evap, stress, wind, kpbl,       &
                           Statein%prsi, del, Statein%prsl, Statein%prslk,          &
                           Statein%phii, Statein%phil, dtp, dusfc1, dvsfc1,         &
                           dtsfc1, dqsfc1, Diag%hpbl, gamt, gamq, dkt, kinver,      &
                           Model%xkzm_m, Model%xkzm_h)
            else
              if (Model%me==0) write(0,*) 'CCPP DEBUG: calling non-CCPP compliant version of moninp'
              call moninp(ix, im, levs, nvdiff, dvdt, dudt, dtdt, dqdt,             &
                           Statein%ugrs, Statein%vgrs, Statein%tgrs, Statein%qgrs,  &
                           Statein%prsik(1,1), rb, Sfcprop%ffmm, Sfcprop%ffhh,      &
                           Sfcprop%tsfc, qss, hflx, evap, stress, wind, kpbl,       &
                           Statein%prsi, del, Statein%prsl, Statein%phii,           &
                           Statein%phil, dtp, dusfc1, dvsfc1, dtsfc1, dqsfc1,       &
                           Diag%hpbl, gamt, gamq, dkt, Model%xkzm_m, Model%xkzm_h)
            endif

          endif   ! end if_hybedmf
        endif     ! end if_do_shoc
      else
        allocate(vdftra(ix,levs,nvdiff), dvdftra(im,levs,nvdiff))
        dvdftra(:,:,:) = 0.0
#ifdef CCPP
        if (Model%me==0) write(0,*) 'CCPP DEBUG: calling PBL_generic_pre through option B'
        ! Copy local variables from driver to appropriate interstitial variables
        !Interstitial(nt)%im = im             ! intent(in) - set in Interstitial(nt)%create()
        !Model%levs                           ! intent(in)
        !Model%nvdiff                         ! intent(in)
        !Model%ntrac                          ! intent(in)
        !Model%imp_physics                    ! intent(in)
        !Model%imp_physics_gfdl               ! intent(in)
        !Model%imp_physics_thompson           ! intent(in)
        !Model%imp_physics_wsm6               ! intent(in)
        !Model%ltaerosol                      ! intent(in)
        !Statein%qgrs                         ! intent(in)
        !Statein%qgrs(:,:,1)                  ! intent(in)
        !Statein%qgrs(:,:,ntcw)               ! intent(in)
        !Statein%qgrs(:,:,ntiw)               ! intent(in)
        !Statein%qgrs(:,:,ntoz)               ! intent(in)
        !Statein%qgrs(:,:,ntlnc)              ! intent(in)
        !Statein%qgrs(:,:,ntinc)              ! intent(in)
        !Statein%qgrs(:,:,ntwa)               ! intent(in)
        !Statein%qgrs(:,:,ntia)               ! intent(in)
        !Statein%qgrs(:,:,ntrw)               ! intent(in)
        !Statein%qgrs(:,:,ntsw)               ! intent(in)
        !Statein%qgrs(:,:,ntgl)               ! intent(in)
        !Interstitial(nt)%vdftra = vdftra     ! intent(inout)
        !cdata_block(nb,nt)%errmsg = errmsg   ! intent(out)
        !cdata_block(nb,nt)%errflg = errflg   ! intent(out)
        call ccpp_physics_run(cdata_block(nb,nt), scheme_name="GFS_PBL_generic_pre", ierr=ierr)
        ! Copy back intent(inout) interstitial variables to local variables in driver
        vdftra = Interstitial(nt)%vdftra
        errmsg = trim(cdata_block(nb,nt)%errmsg)
        errflg = cdata_block(nb,nt)%errflg
        if (errflg/=0) then
            write(0,*) 'Error in call to GFS_PBL_generic_pre: ' // trim(errmsg)
            stop
        end if
        ! DH* update local copy of ntiwx with the value set in
        ! Interstitial(nt)%create(); only needed for hybrid mode
        ntiwx = Interstitial(nt)%ntiwx
        ! *DH
#else
        ntiwx = 0
!
        if (imp_physics == Model%imp_physics_wsm6) then
! WSM6
          do k=1,levs
            do i=1,im
              vdftra(i,k,1) = Statein%qgrs(i,k,1)
              vdftra(i,k,2) = Statein%qgrs(i,k,ntcw)
              vdftra(i,k,3) = Statein%qgrs(i,k,ntiw)
              vdftra(i,k,4) = Statein%qgrs(i,k,ntoz)
            enddo
          enddo
          ntiwx = 3
        elseif (imp_physics == Model%imp_physics_thompson) then
! Thompson
          if(Model%ltaerosol) then
            do k=1,levs
              do i=1,im
                vdftra(i,k,1) = Statein%qgrs(i,k,1)
                vdftra(i,k,2) = Statein%qgrs(i,k,ntcw)
                vdftra(i,k,3) = Statein%qgrs(i,k,ntiw)
                vdftra(i,k,4) = Statein%qgrs(i,k,ntlnc)
                vdftra(i,k,5) = Statein%qgrs(i,k,ntinc)
                vdftra(i,k,6) = Statein%qgrs(i,k,ntoz)
                vdftra(i,k,7) = Statein%qgrs(i,k,ntwa)
                vdftra(i,k,8) = Statein%qgrs(i,k,ntia)
              enddo
            enddo
            ntiwx = 3
          else
            do k=1,levs
              do i=1,im
                vdftra(i,k,1) = Statein%qgrs(i,k,1)
                vdftra(i,k,2) = Statein%qgrs(i,k,ntcw)
                vdftra(i,k,3) = Statein%qgrs(i,k,ntiw)
                vdftra(i,k,4) = Statein%qgrs(i,k,ntinc)
                vdftra(i,k,5) = Statein%qgrs(i,k,ntoz)
              enddo
            enddo
            ntiwx = 3
          endif
!
        elseif (imp_physics == Model%imp_physics_gfdl) then
! GFDL MP
          do k=1,levs
            do i=1,im
              vdftra(i,k,1) = Statein%qgrs(i,k,1)
              vdftra(i,k,2) = Statein%qgrs(i,k,ntcw)
              vdftra(i,k,3) = Statein%qgrs(i,k,ntiw)
              vdftra(i,k,4) = Statein%qgrs(i,k,ntrw)
              vdftra(i,k,5) = Statein%qgrs(i,k,ntsw)
              vdftra(i,k,6) = Statein%qgrs(i,k,ntgl)
              vdftra(i,k,7) = Statein%qgrs(i,k,ntoz)
            enddo
          enddo
          ntiwx = 3
        endif

        if (Model%satmedmf) then
          do k=1,levs
            do i=1,im
              vdftra(i,k,ntkev) = Statein%qgrs(i,k,ntke)
            enddo
          enddo
        endif

#endif
!
        if (Model%do_shoc) then
#ifdef CCPP
          if (Model%me==0) write(0,*) 'CCPP DEBUG: calling moninshoc through option B'
          ! Copy local variables from driver to appropriate interstitial variables
          !Interstitial(nt)%ix = ix             ! intent(in) - set in Interstitial(nt)%create()
          !Interstitial(nt)%im = im             ! intent(in) - set in Interstitial(nt)%create()
          !Model%levs                           ! intent(in)
          Interstitial(nt)%nvdiff = nvdiff      ! intent(in)
          !Model%ntcw                           ! intent(in)
          Interstitial(nt)%nncl = nncl          ! intent(in)
          Interstitial(nt)%dvdt = dvdt          ! intent(inout)
          Interstitial(nt)%dudt = dudt          ! intent(inout)
          Interstitial(nt)%dtdt = dtdt          ! intent(inout)
          !Interstitial(nt)%dvdftra = dvdftra   ! intent(inout)
          !Statein%ugrs                         ! intent(in)
          !Statein%vgrs                         ! intent(in)
          !Statein%tgrs                         ! intent(in)
          !Interstitial(nt)%vdftra = vdftra     ! intent(in)
          !Tbd%phy_f3d(1,1,ntot3d-1)            ! intent(in)
          Interstitial(nt)%prnum = prnum        ! intent(inout)
          !Model%ntke                           ! intent(in)
          !Statein%prsik                        ! intent(in)
          Interstitial(nt)%rb = rb              ! intent(in)
          !Sfcprop%zorl                         ! intent(in)
          !Diag%u10m                            ! intent(in)
          !Diag%v10m                            ! intent(in)
          !Sfcprop%ffmm                         ! intent(in)
          !Sfcprop%ffhh                         ! intent(in)
          !Sfcprop%tsfc                         ! intent(in)
          Interstitial(nt)%hflx = hflx          ! intent(in)
          Interstitial(nt)%evap = evap          ! intent(in)
          Interstitial(nt)%stress = stress      ! intent(in)
          Interstitial(nt)%wind = wind          ! intent(in)
          Interstitial(nt)%kpbl = kpbl          ! intent(out)
          !Statein%prsi                         ! intent(in)
          Interstitial(nt)%del = del            ! intent(in)
          !Statein%prsl                         ! intent(in)
          !Statein%prslk                        ! intent(in)
          !Statein%phii                         ! intent(in)
          !Statein%phil                         ! intent(in)
          !Model%dtp                            ! intent(in)
          Interstitial(nt)%dusfc1 = dusfc1      ! intent(out)
          Interstitial(nt)%dvsfc1 = dvsfc1      ! intent(out)
          Interstitial(nt)%dtsfc1 = dtsfc1      ! intent(out)
          Interstitial(nt)%dqsfc1 = dqsfc1      ! intent(out)
          Interstitial(nt)%dkt  = dkt           ! intent(out)
          !Diag%hpbl                            ! intent(out)
          Interstitial(nt)%kinver = kinver      ! intent(in)
          !Model%xkzm_m                         ! intent(in)
          !Model%xkzm_h                         ! intent(in)
          !Model%xkzm_s                         ! intent(in)
          !Model%lprnt                          ! intent(in)
          Interstitial(nt)%ipr = ipr            ! intent(in)
          !Model%me                             ! intent(in)
          !con_g, con_rd, con_cp, con_hvap, con_fvirt ! intent(in) - physical constant in physcons.F90
          !cdata_block(nb,nt)%errmsg = errmsg   ! intent(out)
          !cdata_block(nb,nt)%errflg = errflg   ! intent(out)
          call ccpp_physics_run(cdata_block(nb,nt), scheme_name="moninshoc", ierr=ierr)
          ! Copy back intent(inout) interstitial variables to local variables in driver
          dvdt   = Interstitial(nt)%dvdt
          dudt   = Interstitial(nt)%dudt
          dtdt   = Interstitial(nt)%dtdt
          !dvdftra = Interstitial(nt)%dvdftra
          prnum  = Interstitial(nt)%prnum
          kpbl   = Interstitial(nt)%kpbl
          dusfc1 = Interstitial(nt)%dusfc1
          dvsfc1 = Interstitial(nt)%dvsfc1
          dtsfc1 = Interstitial(nt)%dtsfc1
          dqsfc1 = Interstitial(nt)%dqsfc1
          dkt    = Interstitial(nt)%dkt
          errmsg = trim(cdata_block(nb,nt)%errmsg)
          errflg = cdata_block(nb,nt)%errflg
          if (errflg/=0) then
              write(0,*) 'Error in call to moninshoc: ' // trim(errmsg)
              stop
          end if
#else
          if (Model%me==0) write(0,*) 'CCPP DEBUG: calling non-CCPP compliant version of moninshoc'
          call moninshoc(ix, im, levs, nvdiff, ntcw, nncl, dvdt, dudt, dtdt, dvdftra, &
                         Statein%ugrs, Statein%vgrs, Statein%tgrs, vdftra,            &
                         Tbd%phy_f3d(1,1,ntot3d-1), prnum, ntke,                      &
                         Statein%prsik(1,1), rb, Sfcprop%zorl, Diag%u10m,             &
                         Diag%v10m, Sfcprop%ffmm, Sfcprop%ffhh, Sfcprop%tsfc, hflx,   &
                         evap, stress, wind, kpbl, Statein%prsi, del, Statein%prsl,   &
                         Statein%prslk, Statein%phii, Statein%phil, dtp, dusfc1,      &
                         dvsfc1, dtsfc1, dqsfc1, dkt, Diag%hpbl, kinver,              &
                         Model%xkzm_m, Model%xkzm_h, Model%xkzm_s, lprnt, ipr, me)
#endif
        else
          if (Model%satmedmf) then
#ifdef CCPP
              if (Model%me==0) write(0,*) 'CCPP DEBUG: calling satmedmfvdif through option B'
              ! Copy local variables from driver to appropriate interstitial variables
              !Interstitial(nt)%im = im             ! intent(in) - set in Interstitial(nt)%create()
              !Interstitial(nt)%ix = ix             ! intent(in) - set in Interstitial(nt)%create()
              !Model%levs                           ! intent(in)
              Interstitial(nt)%nvdiff = nvdiff      ! intent(in)
              !Model%ntcw                           ! intent(in)
              !Interstitial(nt)%ntkev = ntiwx       ! intent(in) - set in Interstitial(nt)%create()
              Interstitial(nt)%ntkev = ntkev        ! intent(in)
              !con_g
              !con_rd
              !con_cp
              !con_rv
              !con_hvap
              !con_hfus
              !con_fvirt
              !con_eps
              !con_epsm1
              Interstitial(nt)%dvdt = dvdt          ! intent(inout)
              Interstitial(nt)%dudt = dudt          ! intent(inout)
              Interstitial(nt)%dtdt = dtdt          ! intent(inout)
              !Interstitial(nt)%dvdftra             ! intent(inout)
              !Statein%ugrs                         ! intent(in)
              !Statein%vgrs                         ! intent(in)
              !Statein%tgrs                         ! intent(in)
              !Interstitial(nt)%vdftra              ! intent(in)
              !Radtend%htrsw                        ! intent(in)
              !Radtend%htrlw                        ! intent(in)
              Interstitial(nt)%xmu = xmu            ! intent(in)
              !Grid%area = garea                    ! intent(in)
              !Statein%prsik(1,1)                   ! intent(in)
              Interstitial(nt)%rb = rb              ! intent(in)
              !Sfcprop%zorl                         ! intent(in)
              !Diag%u10m                            ! intent(in)
              !Diag%v10m                            ! intent(in)
              !Sfcprop%ffmm                         ! intent(in)
              !Sfcprop%ffhh                         ! intent(in)
              !Sfcprop%tsfc                         ! intent(in)
              Interstitial(nt)%hflx = hflx          ! intent(in)
              Interstitial(nt)%evap = evap          ! intent(in)
              Interstitial(nt)%stress = stress      ! intent(in)
              Interstitial(nt)%wind = wind          ! intent(in)
              Interstitial(nt)%kpbl = kpbl          ! intent(out)
              !Statein%prsi                         ! intent(in)
              Interstitial(nt)%del = del            ! intent(in)
              !Statein%prsl                         ! intent(in)
              !Statein%prslk                        ! intent(in)
              !Statein%phii                         ! intent(in)
              !Statein%phil                         ! intent(in)
              !Model%dtp                            ! intent(in)
              !Model%dspheat                        ! intent(in)
              Interstitial(nt)%dusfc1 = dusfc1      ! intent(out)
              Interstitial(nt)%dvsfc1 = dvsfc1      ! intent(out)
              Interstitial(nt)%dtsfc1 = dtsfc1      ! intent(out)
              Interstitial(nt)%dqsfc1 = dqsfc1      ! intent(out)
              !Diag%hpbl                            ! intent(out)
              Interstitial(nt)%kinver = kinver      ! intent(in)
              !Model%xkzm_m                         ! intent(in)
              !Model%xkzm_h                         ! intent(in)
              !Model%xkzm_s                         ! intent(in)
              !cdata_block(nb,nt)%errmsg = errmsg   ! intent(out)
              !cdata_block(nb,nt)%errflg = errflg   ! intent(out)
              call ccpp_physics_run(cdata_block(nb,nt), scheme_name="satmedmfvdif", ierr=ierr)
              ! Copy back intent(inout) interstitial variables to local variables in driver
              dvdt   = Interstitial(nt)%dvdt
              dudt   = Interstitial(nt)%dudt
              dtdt   = Interstitial(nt)%dtdt
              !dvdftra = Interstitial(nt)%dvdftra
              kpbl   = Interstitial(nt)%kpbl
              dusfc1 = Interstitial(nt)%dusfc1
              dvsfc1 = Interstitial(nt)%dvsfc1
              dtsfc1 = Interstitial(nt)%dtsfc1
              dqsfc1 = Interstitial(nt)%dqsfc1
              errmsg = trim(cdata_block(nb,nt)%errmsg)
              errflg = cdata_block(nb,nt)%errflg
              if (errflg/=0) then
                  write(0,*) 'Error in call to satmedmfvdif: ' // trim(errmsg)
                  stop
              end if
#else
              if (Model%me==0) write(0,*) 'CCPP DEBUG: calling non-CCPP compliant version of satmedmfvdif'
              call satmedmfvdif(ix, im, levs, nvdiff, ntcw, ntiwx, ntkev,             &
                         dvdt, dudt, dtdt, dvdftra,                                   &
                         Statein%ugrs, Statein%vgrs, Statein%tgrs, vdftra,            &
                         Radtend%htrsw, Radtend%htrlw, xmu, garea,                    &
                         Statein%prsik(1,1), rb, Sfcprop%zorl, Diag%u10m, Diag%v10m,  &
                         Sfcprop%ffmm, Sfcprop%ffhh, Sfcprop%tsfc, hflx, evap,        &
                         stress, wind, kpbl, Statein%prsi, del, Statein%prsl,         &
                         Statein%prslk, Statein%phii, Statein%phil, dtp,              &
                         Model%dspheat, dusfc1, dvsfc1, dtsfc1, dqsfc1, Diag%hpbl,    &
                         kinver, Model%xkzm_m, Model%xkzm_h, Model%xkzm_s)
#endif
          elseif (Model%hybedmf) then
#ifdef CCPP
            if (Model%me==0) write(0,*) 'CCPP DEBUG: calling hedmf through option B'
            ! Copy local variables from driver to appropriate interstitial variables
            !Interstitial(nt)%im = im             ! intent(in) - set in Interstitial(nt)%create()
            !Interstitial(nt)%ix = ix             ! intent(in) - set in Interstitial(nt)%create()
            !Model%levs                           ! intent(in)
            Interstitial(nt)%nvdiff = nvdiff      ! intent(in)
            !Model%ntcw                           ! intent(in)
            Interstitial(nt)%dvdt = dvdt          ! intent(inout)
            Interstitial(nt)%dudt = dudt          ! intent(inout)
            Interstitial(nt)%dtdt = dtdt          ! intent(inout)
            !Interstitial(nt)%dvdftra = dvdftra   ! intent(inout)
            !Statein%ugrs                         ! intent(in)
            !Statein%vgrs                         ! intent(in)
            !Statein%tgrs                         ! intent(in)
            !Interstitial(nt)%vdftra = vdftra     ! intent(in)
            !Radtend%htrsw                        ! intent(in)
            !Radtend%htrlw                        ! intent(in)
            Interstitial(nt)%xmu = xmu            ! intent(in)
            !Statein%prsik(1,1)                   ! intent(in)
            Interstitial(nt)%rb = rb              ! intent(in)
            !Sfcprop%zorl                         ! intent(in)
            !Diag%u10m                            ! intent(in)
            !Diag%v10m                            ! intent(in)
            !Sfcprop%ffmm                         ! intent(in)
            !Sfcprop%ffhh                         ! intent(in)
            !Sfcprop%tsfc                         ! intent(in)
            Interstitial(nt)%hflx = hflx          ! intent(in)
            Interstitial(nt)%evap = evap          ! intent(in)
            Interstitial(nt)%stress = stress      ! intent(in)
            Interstitial(nt)%wind = wind          ! intent(in)
            Interstitial(nt)%kpbl = kpbl          ! intent(out)
            !Statein%prsi                         ! intent(in)
            Interstitial(nt)%del = del            ! intent(in)
            !Statein%prsl                         ! intent(in)
            !Statein%prslk                        ! intent(in)
            !Statein%phii                         ! intent(in)
            !Statein%phil                         ! intent(in)
            !Model%dtp                            ! intent(in)
            !Model%dspheat                        ! intent(in)
            Interstitial(nt)%dusfc1 = dusfc1      ! intent(out)
            Interstitial(nt)%dvsfc1 = dvsfc1      ! intent(out)
            Interstitial(nt)%dtsfc1 = dtsfc1      ! intent(out)
            Interstitial(nt)%dqsfc1 = dqsfc1      ! intent(out)
            !Diag%hpbl                            ! intent(out)
            Interstitial(nt)%gamt = gamt          ! intent(inout)
            Interstitial(nt)%gamq = gamq          ! intent(inout)
            Interstitial(nt)%dkt  = dkt           ! intent(out)
            Interstitial(nt)%kinver = kinver      ! intent(in)
            !Model%xkzm_m                         ! intent(in)
            !Model%xkzm_h                         ! intent(in)
            !Model%xkzm_s                         ! intent(in)
            !Model%lprnt                          ! intent(in)
            Interstitial(nt)%ipr = ipr            ! intent(in)
            !Model%xkzminv                        ! intent(in)
            !Model%moninq_fac                     ! intent(in)
            !cdata_block(nb,nt)%errmsg = errmsg   ! intent(out)
            !cdata_block(nb,nt)%errflg = errflg   ! intent(out)
            call ccpp_physics_run(cdata_block(nb,nt), scheme_name="hedmf", ierr=ierr)
            ! Copy back intent(inout) interstitial variables to local variables in driver
            dvdt   = Interstitial(nt)%dvdt
            dudt   = Interstitial(nt)%dudt
            dtdt   = Interstitial(nt)%dtdt
            !dvdftra = Interstitial(nt)%dvdftra
            kpbl   = Interstitial(nt)%kpbl
            dusfc1 = Interstitial(nt)%dusfc1
            dvsfc1 = Interstitial(nt)%dvsfc1
            dtsfc1 = Interstitial(nt)%dtsfc1
            dqsfc1 = Interstitial(nt)%dqsfc1
            gamt   = Interstitial(nt)%gamt
            gamq   = Interstitial(nt)%gamq
            dkt    = Interstitial(nt)%dkt
            errmsg = trim(cdata_block(nb,nt)%errmsg)
            errflg = cdata_block(nb,nt)%errflg
            if (errflg/=0) then
                write(0,*) 'Error in call to hedmf: ' // trim(errmsg)
                stop
            end if
#else
            if (Model%me==0) write(0,*) 'CCPP DEBUG: calling non-CCPP compliant version of moninedmf'
            call moninedmf(ix, im, levs, nvdiff, ntcw, dvdt, dudt, dtdt, dvdftra,     &
                           Statein%ugrs, Statein%vgrs, Statein%tgrs, vdftra,          &
                           Radtend%htrsw, Radtend%htrlw, xmu, Statein%prsik(1,1),     &
                           rb, Sfcprop%zorl, Diag%u10m, Diag%v10m, Sfcprop%ffmm,      &
                           Sfcprop%ffhh, Sfcprop%tsfc, qss, hflx, evap, stress,       &
                           wind, kpbl, Statein%prsi, del, Statein%prsl,               &
                           Statein%prslk, Statein%phii, Statein%phil, dtp,            &
                           Model%dspheat, dusfc1, dvsfc1, dtsfc1, dqsfc1, Diag%hpbl,  &
                           gamt, gamq, dkt, kinver, Model%xkzm_m, Model%xkzm_h,       &
                           Model%xkzm_s, lprnt, ipr,                                  &
                           Model%xkzminv, Model%moninq_fac)
#endif
          elseif (.not. Model%old_monin) then
            if (Model%me==0) write(0,*) 'CCPP DEBUG: calling non-CCPP compliant version of moninq'
            call moninq(ix, im, levs, nvdiff, ntcw, dvdt, dudt, dtdt, dvdftra,        &
                        Statein%ugrs, Statein%vgrs, Statein%tgrs, vdftra,             &
                        Radtend%htrsw, Radtend%htrlw, xmu, Statein%prsik(1,1), rb,    &
                        Sfcprop%ffmm, Sfcprop%ffhh, Sfcprop%tsfc, qss, hflx, evap,    &
                        stress, wind, kpbl, Statein%prsi, del, Statein%prsl,          &
                        Statein%prslk, Statein%phii, Statein%phil, dtp,               &
                        Model%dspheat, dusfc1, dvsfc1, dtsfc1, dqsfc1, Diag%hpbl,     &
                        gamt, gamq, dkt, kinver, Model%xkzm_m, Model%xkzm_h,          &
                        Model%xkzm_s, lprnt, ipr,                                     &
                        Model%xkzminv, Model%moninq_fac, Model%rbcr)
          else
            if (Model%mstrat) then
              if (Model%me==0) write(0,*) 'CCPP DEBUG: calling non-CCPP compliant version of moninp1'
              call moninp1(ix, im, levs, nvdiff, dvdt, dudt, dtdt, dvdftra,           &
                           Statein%ugrs, Statein%vgrs, Statein%tgrs, vdftra,          &
                           Statein%prsik(1,1), rb, Sfcprop%ffmm, Sfcprop%ffhh,        &
                           Sfcprop%tsfc, qss, hflx, evap, stress, wind, kpbl,         &
                           Statein%prsi, del, Statein%prsl, Statein%prslk,            &
                           Statein%phii, Statein%phil, dtp, dusfc1, dvsfc1,           &
                           dtsfc1, dqsfc1, Diag%hpbl, gamt, gamq, dkt, kinver,        &
                           Model%xkzm_m, Model%xkzm_h)
            else
              if (Model%me==0) write(0,*) 'CCPP DEBUG: calling non-CCPP compliant version of moninp'
              call moninp(ix, im, levs, nvdiff, dvdt, dudt, dtdt, dvdftra,            &
                          Statein%ugrs, Statein%vgrs, Statein%tgrs, vdftra,           &
                          Statein%prsik(1,1), rb, Sfcprop%ffmm, Sfcprop%ffhh,         &
                          Sfcprop%tsfc, qss, hflx, evap, stress, wind, kpbl,          &
                          Statein%prsi, del, Statein%prsl, Statein%phii,              &
                          Statein%phil, dtp, dusfc1, dvsfc1, dtsfc1, dqsfc1,          &
                          Diag%hpbl, gamt, gamq, dkt, Model%xkzm_m, Model%xkzm_h)
            endif

          endif   ! end if_hybedmf
        endif     ! end if_do_shoc
#ifndef CCPP
        if (imp_physics == Model%imp_physics_wsm6) then
! WSM6
          do k=1,levs
            do i=1,im
              dqdt(i,k,1)     = dvdftra(i,k,1)
              dqdt(i,k,ntcw)  = dvdftra(i,k,2)
              dqdt(i,k,ntiw)  = dvdftra(i,k,3)
              dqdt(i,k,ntoz)  = dvdftra(i,k,4)
            enddo
          enddo
        elseif (imp_physics == Model%imp_physics_thompson) then
! Thompson
          if(Model%ltaerosol) then
            do k=1,levs
              do i=1,im
                dqdt(i,k,1)     = dvdftra(i,k,1)
                dqdt(i,k,ntcw)  = dvdftra(i,k,2)
                dqdt(i,k,ntiw)  = dvdftra(i,k,3)
                dqdt(i,k,ntlnc) = dvdftra(i,k,4)
                dqdt(i,k,ntinc) = dvdftra(i,k,5)
                dqdt(i,k,ntoz)  = dvdftra(i,k,6)
                dqdt(i,k,ntwa)  = dvdftra(i,k,7)
                dqdt(i,k,ntia)  = dvdftra(i,k,8)
              enddo
            enddo
          else
            do k=1,levs
              do i=1,im
                dqdt(i,k,1)     = dvdftra(i,k,1)
                dqdt(i,k,ntcw)  = dvdftra(i,k,2)
                dqdt(i,k,ntiw)  = dvdftra(i,k,3)
                dqdt(i,k,ntinc) = dvdftra(i,k,4)
                dqdt(i,k,ntoz)  = dvdftra(i,k,5)
              enddo
            enddo
          endif
!
        elseif (imp_physics == Model%imp_physics_gfdl) then
! GFDL MP
          do k=1,levs
            do i=1,im
              dqdt(i,k,1)    = dvdftra(i,k,1)
              dqdt(i,k,ntcw) = dvdftra(i,k,2)
              dqdt(i,k,ntiw) = dvdftra(i,k,3)
              dqdt(i,k,ntrw) = dvdftra(i,k,4)
              dqdt(i,k,ntsw) = dvdftra(i,k,5)
              dqdt(i,k,ntgl) = dvdftra(i,k,6)
              dqdt(i,k,ntoz) = dvdftra(i,k,7)
            enddo
          enddo
        endif

        if (Model%satmedmf) then
          do k=1,levs
            do i=1,im
              dqdt(i,k,ntke)  = dvdftra(i,k,ntkev)
            enddo
          enddo
        endif

#endif
        deallocate(vdftra, dvdftra)
      endif

! DH* note: this block is not yet in CCPP
      if (Model%cplflx) then
        do i=1,im
          if (flag_cice(i)) then
#ifdef CCPP
            Sfcprop%fice(i) = fice_cice(i)
#else
                    cice(i) = fice_cice(i)
#endif
                  dusfc1(i) = dusfc_cice(i)
                  dvsfc1(i) = dvsfc_cice(i)
                  dqsfc1(i) = dqsfc_cice(i)
                  dtsfc1(i) = dtsfc_cice(i)
          endif
        enddo
      endif

      if (Model%cplchm) then
        do i = 1, im
          tem1 = max(Diag%q1(i), 1.e-8)
          tem  = Statein%prsl(i,1) / (con_rd*Diag%t1(i)*(1.0+con_fvirt*tem1))
          Coupling%ushfsfci(i) = -con_cp * tem * hflx(i) ! upward sensible heat flux
        enddo
#ifdef CCPP
        ! Interstitial(nt)%dkt has dimensions (1:im,1:Model%levs-1) as required by physics schemes
        Coupling%dkt     (:,1:Model%levs-1) = dkt (:,:)
#else
        Coupling%dkt     (:,:) = dkt (:,:)
#endif
      endif
! *DH

!     if (lprnt) then
!       write(0,*) ' dusfc1=',dusfc1(ipr),' kdt=',kdt,' lat=',lat
!       write(0,*)' dtsfc1=',dtsfc1(ipr)
!       write(0,*)' dqsfc1=',dqsfc1(ipr)
!       write(0,*)' dtdtc=',(dtdt(ipr,k),k=1,15)
!       write(0,*)' dqdtc=',(dqdt(ipr,k,1),k=1,15)
!       print *,' dudtm=',dudt(ipr,:)
!     endif

!  --- ...  coupling insertion
#ifdef CCPP
      if (Model%me==0) write(0,*) 'CCPP DEBUG: calling GFS_PBL_generic_post through option B'
      ! Copy local variables from driver to appropriate interstitial variables
      !Interstitial(nt)%im = im                           ! intent(in) - set in Interstitial(nt)%create()
      !Model%levs                                         ! intent(in)
      !Model%nvdiff                                       ! intent(in)
      !Model%ntrac                                        ! intent(in)
      !Model%ntoz                                         ! intent(in)
      !Model%imp_physics                                  ! intent(in)
      !Model%imp_physics_gfdl                             ! intent(in)
      !Model%imp_physics_thompson                         ! intent(in)
      !Model%imp_physics_wsm6                             ! intent(in)
      !Model%ltaerosol                                    ! intent(in)
      !Model%cplflx                                       ! intent(in)
      !Model%lssav                                        ! intent(in)
      !Model%ldiag3d                                      ! intent(in)
      !Model%lsidea                                       ! intent(in)
      !Model%hybedmf                                      ! intent(in)
      !Model%do_shoc                                      ! intent(in)
      !Model%satmedmf                                     ! intent(in)
      !Interstitial(nt)%dvdftra = dvdftra                 ! intent(in)
      Interstitial(nt)%dusfc1 = dusfc1                    ! intent(in)
      Interstitial(nt)%dvsfc1 = dvsfc1                    ! intent(in)
      Interstitial(nt)%dtsfc1 = dtsfc1                    ! intent(in)
      Interstitial(nt)%dqsfc1 = dqsfc1                    ! intent(in)
      !Model%dtf                                          ! intent(in)
      Interstitial(nt)%dudt = dudt                        ! intent(in)
      Interstitial(nt)%dvdt = dvdt                        ! intent(in)
      Interstitial(nt)%dtdt = dtdt                        ! intent(in)
      !Rattend%htrsw                                      ! intent(in)
      !Radtend%htrlw                                      ! intent(in)
      Interstitial(nt)%xmu = xmu                          ! intent(in)
      Interstitial(nt)%dqdt = dqdt                        ! intent(inout)
      !Interstitial(nt)%dqdt(:,:,1) = dqdt(:,:,1)          ! intent(inout)
      !Interstitial(nt)%dqdt(:,:,ntcw) = dqdt(:,:,ntcw)    ! intent(inout)
      !Interstitial(nt)%dqdt(:,:,ntiw) = dqdt(:,:,ntiw)    ! intent(inout)
      !Interstitial(nt)%dqdt(:,:,ntoz) = dqdt(:,:,ntoz)    ! intent(inout)
      !Interstitial(nt)%dqdt(:,:,ntlnc) = dqdt(:,:,ntlnc)  ! intent(inout)
      !Interstitial(nt)%dqdt(:,:,ntinc) = dqdt(:,:,ntinc)  ! intent(inout)
      !Interstitial(nt)%dqdt(:,:,ntwa) = dqdt(:,:,ntwa)    ! intent(inout)
      !Interstitial(nt)%dqdt(:,:,ntia) = dqdt(:,:,ntia)    ! intent(inout)
      !Interstitial(nt)%dqdt(:,:,ntrw) = dqdt(:,:,ntrw)    ! intent(inout)
      !Interstitial(nt)%dqdt(:,:,ntsw) = dqdt(:,:,ntsw)    ! intent(inout)
      !Interstitial(nt)%dqdt(:,:,ntgl) = dqdt(:,:,ntgl)    ! intent(inout)
      !Coupling%dusfc_cpl                                 ! intent(inout)
      !Coupling%dvsfc_cpl                                 ! intent(inout)
      !Coupling%dtsfc_cpl                                 ! intent(inout)
      !Coupling%dqsfc_cpl                                 ! intent(inout)
      !Coupling%dusfci_cpl                                ! intent(inout)
      !Coupling%dvsfci_cpl                                ! intent(inout)
      !Coupling%dtsfci_cpl                                ! intent(inout)
      !Coupling%dqsfci_cpl                                ! intent(inout)
      !Diag%dusfc                                         ! intent(inout)
      !Diag%dvsfc                                         ! intent(inout)
      !Diag%dtsfc                                         ! intent(inout)
      !Diag%dqsfc                                         ! intent(inout)
      !Diag%dusfci                                        ! intent(inout)
      !Diag%dvsfci                                        ! intent(inout)
      !Diag%dtsfci                                        ! intent(inout)
      !Diag%dqsfci                                        ! intent(inout)
      !Diag%dt3dt(:,:,3)                                  ! intent(inout)
      !Diag%du3dt(:,:,1)                                  ! intent(inout)
      !Diag%du3dt(:,:,2)                                  ! intent(inout)
      !Diag%dv3dt(:,:,1)                                  ! intent(inout)
      !Diag%dv3dt(:,:,2)                                  ! intent(inout)
      !Diag%dq3dt(:,:,1)                                  ! intent(inout)
      !Diag%dq3dt(:,:,5)                                  ! intent(inout)
      !cdata_block(nb,nt)%errmsg = errmsg                 ! intent(out)
      !cdata_block(nb,nt)%errflg = errflg                 ! intent(out)
      call ccpp_physics_run(cdata_block(nb,nt), scheme_name="GFS_PBL_generic_post", ierr=ierr)
      ! Copy back intent(inout) interstitial variables to local variables in driver
      dqdt = Interstitial(nt)%dqdt
      errmsg = trim(cdata_block(nb,nt)%errmsg)
      errflg = cdata_block(nb,nt)%errflg
      if (errflg/=0) then
          write(0,*) 'Error in call to GFS_PBL_generic_post: ' // trim(errmsg)
          stop
      end if
#else
      if (Model%cplflx) then
        do i=1,im
          Coupling%dusfc_cpl (i) = Coupling%dusfc_cpl(i) + dusfc1(i)*dtf
          Coupling%dvsfc_cpl (i) = Coupling%dvsfc_cpl(i) + dvsfc1(i)*dtf
          Coupling%dtsfc_cpl (i) = Coupling%dtsfc_cpl(i) + dtsfc1(i)*dtf
          Coupling%dqsfc_cpl (i) = Coupling%dqsfc_cpl(i) + dqsfc1(i)*dtf
          Coupling%dusfci_cpl(i) = dusfc1(i)
          Coupling%dvsfci_cpl(i) = dvsfc1(i)
          Coupling%dtsfci_cpl(i) = dtsfc1(i)
          Coupling%dqsfci_cpl(i) = dqsfc1(i)
        enddo
      endif
!-------------------------------------------------------lssav if loop ----------
      if (Model%lssav) then
        do i=1,im
          Diag%dusfc (i) = Diag%dusfc(i) + dusfc1(i)*dtf
          Diag%dvsfc (i) = Diag%dvsfc(i) + dvsfc1(i)*dtf
          Diag%dtsfc (i) = Diag%dtsfc(i) + dtsfc1(i)*dtf
          Diag%dqsfc (i) = Diag%dqsfc(i) + dqsfc1(i)*dtf
          Diag%dusfci(i) = dusfc1(i)
          Diag%dvsfci(i) = dvsfc1(i)
          Diag%dtsfci(i) = dtsfc1(i)
          Diag%dqsfci(i) = dqsfc1(i)
        enddo
!       if (lprnt) then
!         write(0,*)' dusfc=',dusfc(ipr),' dusfc1=',dusfc1(ipr),' dtf=',
!    &     dtf,' kdt=',kdt,' lat=',lat
!       endif

        if (Model%ldiag3d) then
          if (Model%lsidea) then
            Diag%dt3dt(1:im,:,3) = Diag%dt3dt(1:im,:,3) + dtdt(1:im,:)*dtf
          else
            do k=1,levs
              do i=1,im
                tem  = dtdt(i,k) - (Radtend%htrlw(i,k)+Radtend%htrsw(i,k)*xmu(i))
                Diag%dt3dt(i,k,3) = Diag%dt3dt(i,k,3) + tem*dtf
              enddo
            enddo
          endif
          do k=1,levs
            do i=1,im
              Diag%du3dt(i,k,1) = Diag%du3dt(i,k,1) + dudt(i,k) * dtf
              Diag%du3dt(i,k,2) = Diag%du3dt(i,k,2) - dudt(i,k) * dtf
              Diag%dv3dt(i,k,1) = Diag%dv3dt(i,k,1) + dvdt(i,k) * dtf
              Diag%dv3dt(i,k,2) = Diag%dv3dt(i,k,2) - dvdt(i,k) * dtf
            enddo
          enddo
! update dqdt_v to include moisture tendency due to vertical diffusion
!         if (lgocart) then
!           do k=1,levs
!             do i=1,im
!               dqdt_v(i,k)  = dqdt(i,k,1) * dtf
!             enddo
!           enddo
!         endif
!          do k=1,levs
!            do i=1,im
!              tem  = dqdt(i,k,1) * dtf
!              Diag%dq3dt(i,k,1) = Diag%dq3dt(i,k,1) + tem
!            enddo
!          enddo
!          if (ntoz > 0) then
!            do k=1,levs
!              do i=1,im
!                Diag%dq3dt(i,k,5) = Diag%dq3dt(i,k,5) + dqdt(i,k,ntoz) * dtf
!              enddo
!            enddo
!          endif
        endif

      endif   ! end if_lssav
#endif

!-------------------------------------------------------lssav if loop ----------
!
!            Orographic gravity wave drag parameterization
!            ---------------------------------------------
#ifdef CCPP
! OPTION B BEGIN
      if (Model%me==0) write(0,*) 'CCPP DEBUG: calling gwdps_pre through option B'
      ! Copy local variables from driver to appropriate interstitial variables
      !Interstitial(nt)%im       = im         ! intent(in) - set in Interstitial(nt)%create()
      !Model%nmtvr               = nmtvr      ! intent(in)
      Interstitial(nt)%hprime1  = Sfcprop%hprime(:,1)
      Interstitial(nt)%oc       = oc          ! intent(out)
      Interstitial(nt)%oa4      = oa4         ! intent(out)
      Interstitial(nt)%clx      = clx         ! intent(out)
      Interstitial(nt)%theta    = theta       ! intent(out)
      Interstitial(nt)%sigma    = sigma       ! intent(out)
      Interstitial(nt)%gamma    = gamma       ! intent(out)
      Interstitial(nt)%elvmax   = elvmax      ! intent(out)
      !Model%lssav                            ! intent(in)
      !Model%ldiag3d                          ! intent(in)
      Interstitial(nt)%dtdt = dtdt            ! intent(in)
      !Diag%dt3dt(:,:,7)                      ! intent(inout)
      !Model%dtf                              ! intent(in)
      !cdata_block(nb,nt)%errmsg = errmsg     ! intent(out)
      !cdata_block(nb,nt)%errflg = errflg     ! intent(out)
      call ccpp_physics_run(cdata_block(nb,nt), scheme_name="gwdps_pre", ierr=ierr)
      ! Copy back intent(inout) interstitial variables to local variables in driver
      Sfcprop%hprime(:,1) = Interstitial(nt)%hprime1
      oc     = Interstitial(nt)%oc
      oa4    = Interstitial(nt)%oa4
      clx    = Interstitial(nt)%clx
      theta  = Interstitial(nt)%theta
      sigma  = Interstitial(nt)%sigma
      gamma  = Interstitial(nt)%gamma
      elvmax = Interstitial(nt)%elvmax
      errmsg = trim(cdata_block(nb,nt)%errmsg)
      errflg = cdata_block(nb,nt)%errflg
! OPTION B END
      if (errflg/=0) then
          write(0,*) 'Error in call to gwdps_pre: ' // trim(errmsg)
          stop
      end if
#else
      if (Model%me==0) write(0,*) 'CCPP DEBUG: calling non-CCPP compliant version of gwdps_pre'

      if (Model%nmtvr == 14) then         ! current operational - as of 2014
        do i=1,im
          oc(i)     = Sfcprop%hprime(i,2)
          oa4(i,1)  = Sfcprop%hprime(i,3)
          oa4(i,2)  = Sfcprop%hprime(i,4)
          oa4(i,3)  = Sfcprop%hprime(i,5)
          oa4(i,4)  = Sfcprop%hprime(i,6)
          clx(i,1)  = Sfcprop%hprime(i,7)
          clx(i,2)  = Sfcprop%hprime(i,8)
          clx(i,3)  = Sfcprop%hprime(i,9)
          clx(i,4)  = Sfcprop%hprime(i,10)
          theta(i)  = Sfcprop%hprime(i,11)
          gamma(i)  = Sfcprop%hprime(i,12)
          sigma(i)  = Sfcprop%hprime(i,13)
          elvmax(i) = Sfcprop%hprime(i,14)
        enddo
      elseif (Model%nmtvr == 10) then
        do i=1,im
          oc(i)     = Sfcprop%hprime(i,2)
          oa4(i,1)  = Sfcprop%hprime(i,3)
          oa4(i,2)  = Sfcprop%hprime(i,4)
          oa4(i,3)  = Sfcprop%hprime(i,5)
          oa4(i,4)  = Sfcprop%hprime(i,6)
          clx(i,1)  = Sfcprop%hprime(i,7)
          clx(i,2)  = Sfcprop%hprime(i,8)
          clx(i,3)  = Sfcprop%hprime(i,9)
          clx(i,4)  = Sfcprop%hprime(i,10)
        enddo
      elseif (Model%nmtvr == 6) then
        do i=1,im
          oc(i)     = Sfcprop%hprime(i,2)
          oa4(i,1)  = Sfcprop%hprime(i,3)
          oa4(i,2)  = Sfcprop%hprime(i,4)
          oa4(i,3)  = Sfcprop%hprime(i,5)
          oa4(i,4)  = Sfcprop%hprime(i,6)
          clx(i,1)  = 0.0
          clx(i,2)  = 0.0
          clx(i,3)  = 0.0
          clx(i,4)  = 0.0
        enddo
      else
        oc = 0 ; oa4 = 0 ; clx = 0 ; theta = 0 ; gamma = 0 ; sigma = 0
        elvmax = 0

      endif   ! end if_nmtvr

      if (Model%lssav) then
        if (Model%ldiag3d) then
          do k=1,levs
            do i=1,im

#ifdef WORKAROUND_DT3DT7
              Diag%dt3dt(i,k,7) = Diag%dt3dt(i,k,7) - dtdt(i,k)
#else
              Diag%dt3dt(i,k,7) = Diag%dt3dt(i,k,7) - dtdt(i,k)*dtf
#endif
            enddo
         enddo
        endif
      endif

#endif

#ifdef CCPP
! OPTION B BEGIN
      if (Model%me==0) write(0,*) 'CCPP DEBUG: calling gwdps through option B'
      ! Copy local variables from driver to appropriate interstitial variables
      !Interstitial(nt)%im     = im          ! intent(in) - set in Interstitial(nt)%create()
      !Interstitial(nt)%ix     = ix          ! intent(in) - set in Interstitial(nt)%create()
      !Model%levs                            ! intent(in)
      Interstitial(nt)%dvdt   = dvdt         ! A: intent(inout)
      Interstitial(nt)%dudt   = dudt         ! B: intent(inout)
      Interstitial(nt)%dtdt   = dtdt         ! C: intent(inout)
      Interstitial(nt)%kpbl   = kpbl         ! intent(in)
      Interstitial(nt)%del    = del          ! intent(in)
      !Model%dtp                             ! intent(in)
      !Model%kdt                             ! intent(in)
      Interstitial(nt)%hprime1 = Sfcprop%hprime(:,1)  ! intent(in)
      Interstitial(nt)%oc     = oc           ! intent(in)
      Interstitial(nt)%oa4    = oa4          ! intent(in)
      Interstitial(nt)%clx    = clx          ! intent(in)
      Interstitial(nt)%theta  = theta        ! intent(in)
      Interstitial(nt)%sigma  = sigma        ! intent(in)
      Interstitial(nt)%gamma  = gamma        ! intent(in)
      Interstitial(nt)%elvmax = elvmax       ! intent(inout)
      Interstitial(nt)%dusfcg = dusfcg       ! intent(out)
      Interstitial(nt)%dvsfcg = dvsfcg       ! intent(out)
      !Model%me                              ! intent(in)
      !Model%lprnt                           ! intent(in)
      Interstitial(nt)%ipr      = ipr        ! intent(in)
      !cdata_block(nb,nt)%errmsg = errmsg    ! intent(out)
      !cdata_block(nb,nt)%errflg = errflg    ! intent(out)
      call ccpp_physics_run(cdata_block(nb,nt), scheme_name="gwdps", ierr=ierr)
      ! Copy back intent(inout) interstitial variables to local variables in driver
      dvdt   = Interstitial(nt)%dvdt
      dudt   = Interstitial(nt)%dudt
      dtdt   = Interstitial(nt)%dtdt
      elvmax = Interstitial(nt)%elvmax
      dusfcg = Interstitial(nt)%dusfcg
      dvsfcg = Interstitial(nt)%dvsfcg
      !rdxzb  = Diag%zmtnblck
      errmsg = trim(cdata_block(nb,nt)%errmsg)
      errflg = cdata_block(nb,nt)%errflg
! OPTION B END
      if (errflg/=0) then
          write(0,*) 'Error in call to gwdps: ' // trim(errmsg)
          stop
      end if
#else
      if (Model%me==0) write(0,*) 'CCPP DEBUG: calling non-CCPP compliant version of gwdps'

      call gwdps(im, ix, im, levs, dvdt, dudt, dtdt,        &
                 Statein%ugrs, Statein%vgrs, Statein%tgrs,  &
                 Statein%qgrs, kpbl, Statein%prsi, del,     &
                 Statein%prsl, Statein%prslk, Statein%phii, &
                 Statein%phil, dtp, kdt,                    &
                 Sfcprop%hprime(1,1), oc, oa4, clx, theta,  &
                 sigma, gamma, elvmax, dusfcg, dvsfcg,      &
                 con_g, con_cp, con_rd, con_rv, Model%lonr, &
                 Model%nmtvr, Model%cdmbgwd, me, lprnt,ipr, &
                 Diag%zmtnblck)
#endif

#ifdef CCPP
      if (Model%me==0) write(0,*) 'CCPP DEBUG: calling gwdps_post through option B'
      ! Copy local variables from driver to appropriate interstitial variables
      Interstitial(nt)%dusfcg   = dusfcg       ! intent(in)
      Interstitial(nt)%dvsfcg   = dvsfcg       ! intent(in)
      Interstitial(nt)%dudt     = dudt         ! intent(in)
      Interstitial(nt)%dvdt     = dvdt         ! intent(in)
      Interstitial(nt)%dtdt     = dtdt         ! intent(in)
      !cdata_block(nb,nt)%errmsg = errmsg      ! intent(out)
      !cdata_block(nb,nt)%errflg = errflg      ! intent(out)
      call ccpp_physics_run(cdata_block(nb,nt), scheme_name="gwdps_post", ierr=ierr)
      ! Copy back intent(inout) interstitial variables to local variables in driver
      errmsg = trim(cdata_block(nb,nt)%errmsg)
      errflg = cdata_block(nb,nt)%errflg
      if (errflg/=0) then
        write(0,*) 'Error in call to gwdps_post: ' // trim(errmsg)
        stop
      end if
#else
     if (Model%me==0) write(0,*) 'CCPP DEBUG: calling non-CCPP compliant version of gwdps_post'
      if (Model%lssav) then
        do i=1,im
          Diag%dugwd(i) = Diag%dugwd(i) + dusfcg(i)*dtf
          Diag%dvgwd(i) = Diag%dvgwd(i) + dvsfcg(i)*dtf
        enddo

        if (Model%ldiag3d) then
          do k=1,levs
            do i=1,im
              Diag%du3dt(i,k,2) = Diag%du3dt(i,k,2) + dudt(i,k) * dtf
              Diag%dv3dt(i,k,2) = Diag%dv3dt(i,k,2) + dvdt(i,k) * dtf
              Diag%dt3dt(i,k,7) = Diag%dt3dt(i,k,7) + dtdt(i,k) * dtf
            enddo
          enddo
        endif
      endif
#endif

!    Rayleigh damping  near the model top
      if( .not. Model%lsidea .and. Model%ral_ts > 0.0) then
#ifdef CCPP
        ! test for "if( .not. Model%lsidea .and. Model%ral_ts > 0.0) then" also in scheme
        if (Model%me==0) write(0,*) 'CCPP DEBUG: calling rayleigh_damp through option B'
        ! Copy local variables from driver to appropriate interstitial variables
        !Model%lsidea                          ! intent(in   )
        !Interstitial(nt)%im     = im          ! intent(in   ) - set in Interstitial(nt)%create()
        !Interstitial(nt)%ix     = ix          ! intent(in   ) - set in Interstitial(nt)%create()
        !Model%levs                            ! intent(in   )
        Interstitial(nt)%dvdt = dvdt           ! intent(inout)
        Interstitial(nt)%dudt = dudt           ! intent(inout)
        Interstitial(nt)%dtdt = dtdt           ! intent(inout)
        !Statein%ugrs                          ! intent(in   )
        !Statein%vgrs                          ! intent(in   )
        !Model%dtp                             ! intent(in   )
        !con_cp                                ! intent(in   ) - physical constant in physcons.F90
        !Model%levr                            ! intent(in   )
        !Statein%pgr                           ! intent(in   )
        !Statein%prsl                          ! intent(in   )
        !Model%prslrd0                         ! intent(in   )
        !Model%ral_ts                          ! intent(in   )
        !cdata_block(nb,nt)%errmsg = errmsg      ! intent(out)
        !cdata_block(nb,nt)%errflg = errflg      ! intent(out)
        call ccpp_physics_run(cdata_block(nb,nt), scheme_name="rayleigh_damp", ierr=ierr)
        ! Copy back intent(inout) interstitial variables to local variables in driver
        dvdt   = Interstitial(nt)%dvdt
        dudt   = Interstitial(nt)%dudt
        dtdt   = Interstitial(nt)%dtdt
        errmsg = trim(cdata_block(nb,nt)%errmsg)
        errflg = cdata_block(nb,nt)%errflg
        if (errflg/=0) then
          write(0,*) 'Error in call to rayleigh_damp: ' // trim(errmsg)
          stop
        end if
#else
        if (Model%me==0) write(0,*) 'CCPP DEBUG: calling non-CCPP compliant version of rayleigh_damp'
        call rayleigh_damp(im, ix, im, levs, dvdt, dudt, dtdt,      &
                           Statein%ugrs, Statein%vgrs, dtp, con_cp, &
                           Model%levr, Statein%pgr, Statein%prsl,   &
                           Model%prslrd0, Model%ral_ts)
#endif
      endif

!     if (lprnt) then
!       write(0,*)' tgrs1=',(Statein%tgrs(ipr,k),k=1,10)
!       write(0,*)' dtdt=',(dtdt(ipr,k),k=1,10)
!     endif

#ifdef CCPP
      if (Model%me==0) write(0,*) 'CCPP DEBUG: calling GFS_suite_stateout_update through option B'
      !Interstitial(nt)%im = im              ! intent(in) - set in Interstitial(nt)%create()
      !Model%levs                            ! intent(in)
      !Model%ntrac                           ! intent(in)
      !Model%dtp                             ! intent(in)
      !Statein%tgrs                          ! intent(in)
      !Statein%ugrs                          ! intent(in)
      !Statein%vgrs                          ! intent(in)
      !Statein%qgrs                          ! intent(in)
      Interstitial(nt)%dudt = dudt           ! intent(in)
      Interstitial(nt)%dvdt = dvdt           ! intent(in)
      Interstitial(nt)%dtdt = dtdt           ! intent(in)
      Interstitial(nt)%dqdt = dqdt           ! intent(in)
      !Stateout%gt0                          ! intent(out)
      !Stateout%gu0                          ! intent(out)
      !Stateout%gv0                          ! intent(out)
      !Stateout%gq0                          ! intent(out)
      !cdata_block(nb,nt)%errmsg = errmsg    ! intent(out)
      !cdata_block(nb,nt)%errflg = errflg    ! intent(out)
      call ccpp_physics_run(cdata_block(nb,nt), scheme_name="GFS_suite_stateout_update", ierr=ierr)
      ! Copy back intent(inout) interstitial variables to local variables in driver
      errmsg = trim(cdata_block(nb,nt)%errmsg)
      errflg = cdata_block(nb,nt)%errflg
      if (errflg/=0) then
        write(0,*) 'Error in call to GFS_suite_stateout_update: ' // trim(errmsg)
        stop
      end if
#else
      do k=1,levs
        do i=1,im
          Stateout%gt0(i,k)  = Statein%tgrs(i,k) + dtdt(i,k) * dtp
          Stateout%gu0(i,k)  = Statein%ugrs(i,k) + dudt(i,k) * dtp
          Stateout%gv0(i,k)  = Statein%vgrs(i,k) + dvdt(i,k) * dtp
        enddo
      enddo
      Stateout%gq0(1:im,:,:) = Statein%qgrs(1:im,:,:) + dqdt(1:im,:,:) * dtp
#endif

      if (Model%lsidea) then            ! idea convective adjustment
        call ideaca_up(Statein%prsi,Stateout%gt0,ix,im,levs+1)
      endif

!  --- ...  ozone physics
      if (ntoz > 0 .and. ntrac >= ntoz) then
        if (oz_coeff > 4) then
#ifdef CCPP
          if (Model%me==0) write(0,*) 'CCPP DEBUG: calling ozphys_2015 through option B'
          ! Copy local variables from driver to appropriate interstitial variables
          !Interstitial(nt)%ix = ix                        ! intent(in) - set in Interstitial(nt)%create()
          !Interstitial(nt)%im = im                        ! intent(in) - set in Interstitial(nt)%create()
          !Model%levs                                      ! intent(in)
          !Interstitial(nt)%levozp = levozp                    ! intent(in) - associated with levozp
          !Model%dtp                                       ! intent(in)
          !Stateout%gq0(:,:,1)                             ! intent(inout)
          !Stateout%gt0                                    ! intent(in)
          !Interstitial(nt)%oz_pres = oz_pres               ! intent(in) - associated with oz_pres
          !Statein%prsl                                    ! intent(in)
          !Tbd%ozpl = ozpl                                  ! intent(in)
          !Interstitial(nt)%oz_coeff = oz_coeff             ! intent(in) - associated with oz_coeff
          Interstitial(nt)%del = del                       ! intent(in)
          !Model%ldiag3d                                   ! intent(in)
          !Diag%dq3dt(:,:,6:6+Interstitial(nt)%oz_coeff-1) ! intent(inout)
          !con_g
          !Model%me                                        ! intent(in)
          !cdata_block(nb,nt)%errmsg = errmsg              ! intent(out)
          !cdata_block(nb,nt)%errflg = errflg              ! intent(out)
          call ccpp_physics_run(cdata_block(nb,nt), scheme_name="ozphys_2015",ierr=ierr)
          ! Copy back intent(inout) interstitial variables to local variables in
          ! driver
          errmsg = trim(cdata_block(nb,nt)%errmsg)
          errflg = cdata_block(nb,nt)%errflg
          if (errflg/=0) then
            write(0,*) 'Error in call to ozphys_2015: ' // trim(errmsg)
            stop
          end if

#else
          if (Model%me==0) write(0,*) 'CCPP DEBUG: calling non-CCPP compliant version of ozphys_2015'
          call ozphys_2015 (ix, im, levs, levozp, dtp,               &
                            Stateout%gq0(1,1,ntoz),                  &
                            Stateout%gq0(1,1,ntoz),                  &
                            Stateout%gt0, oz_pres, Statein%prsl,     &
                            Tbd%ozpl, oz_coeff, del, Model%ldiag3d,  &
                            dq3dt_loc(1,1,6), me)
!          if (Model%ldiag3d) then
!            do k=1,levs
!              do i=1,im
!                Diag%dq3dt(i,k,6) = dq3dt_loc(i,k,6)
!                Diag%dq3dt(i,k,7) = dq3dt_loc(i,k,7)
!                Diag%dq3dt(i,k,8) = dq3dt_loc(i,k,8)
!                Diag%dq3dt(i,k,9) = dq3dt_loc(i,k,9)
!              enddo
!            enddo
!          endif
#endif
        else
#ifdef CCPP
          if (Model%me==0) write(0,*) 'CCPP DEBUG: calling ozphys through option B'
          ! Copy local variables from driver to appropriate interstitial variables
          !Interstitial(nt)%ix = ix                        ! intent(in) - set in Interstitial(nt)%create()
          !Interstitial(nt)%im = im                        ! intent(in) - set in Interstitial(nt)%create()
          !Model%levs                                      ! intent(in)
          !Interstitial(nt)%levozp = levozp                         ! intent(in) - associated with levozp
          !Model%dtp                                       ! intent(in)
          !Stateout%gq0(:,:,1)                             ! intent(inout)
          !Stateout%gt0                                    ! intent(in)
          !Interstitial(nt)%oz_pres = oz_pres                        ! intent(in) - associated with oz_pres
          !Statein%prsl                                    ! intent(in)
          !Tbd%ozpl = ozpl                                       ! intent(in)
          !Interstitial(nt)%oz_coeff = oz_coeff             ! intent(in) - associated with oz_coeff
          Interstitial(nt)%del = del                       ! intent(in)
          !Model%ldiag3d                                   ! intent(in)
          !Diag%dq3dt(:,:,6:6+Interstitial(nt)%oz_coeff-1) ! intent(inout)
          !Model%me                                        ! intent(in)
          !cdata_block(nb,nt)%errmsg = errmsg              ! intent(out)
          !cdata_block(nb,nt)%errflg = errflg              ! intent(out)
          call ccpp_physics_run(cdata_block(nb,nt), scheme_name="ozphys", ierr=ierr)
          ! Copy back intent(inout) interstitial variables to local variables in driver
          errmsg = trim(cdata_block(nb,nt)%errmsg)
          errflg = cdata_block(nb,nt)%errflg
          if (errflg/=0) then
            write(0,*) 'Error in call to ozphys: ' // trim(errmsg)
            stop
          end if
#else
          if (Model%me==0) write(0,*) 'CCPP DEBUG: calling non-CCPP compliant version of ozphys'
          call ozphys (ix, im, levs, levozp, dtp,                 &
                       Stateout%gq0(1,1,ntoz),                    &
                       Stateout%gq0(1,1,ntoz),                    &
                       Stateout%gt0, oz_pres, Statein%prsl,       &
                       Tbd%ozpl, oz_coeff, del, Model%ldiag3d,    &
                       dq3dt_loc(1,1,6), me)
!          if (Model%ldiag3d) then
!            do k=1,levs
!              do i=1,im
!                Diag%dq3dt(i,k,6) = dq3dt_loc(i,k,6)
!                Diag%dq3dt(i,k,7) = dq3dt_loc(i,k,7)
!                Diag%dq3dt(i,k,8) = dq3dt_loc(i,k,8)
!                Diag%dq3dt(i,k,9) = dq3dt_loc(i,k,9)
!              enddo
!            enddo
!          endif
#endif
        endif
      endif

      if (Model%h2o_phys) then
#ifdef CCPP
! OPTION B BEGIN
        if (Model%me==0) write(0,*) 'CCPP DEBUG: calling h2ophys through option B'
        ! Copy local variables from driver to appropriate interstitial variables
        !Interstitial(nt)%ix = ix                   ! intent(in) - set in Interstitial(nt)%create()
        !Interstitial(nt)%im = im                   ! intent(in) - set in Interstitial(nt)%create()
        !Model%levs                                 ! intent(in)
        !Interstitial(nt)%levh2o                    ! intent(in) - associated with levh2o
        !Model%dtp                                  ! intent(in)
        !Stateout%gq0(:,:,1)                        ! intent(inout)
        !Interstitial(nt)%h2o_pres                  ! intent(in) - associated with h2o_pres
        !Statein%prsl                               ! intent(in)
        !Tbd%h2opl                                  ! intent(in)
        !Interstitial(nt)%h2o_coeff                 ! intent(in) - associated with h2o_coeff
        !Model%ldiag3d                              ! intent(in)
        !Model%me                                   ! intent(in)
        !cdata_block(nb,nt)%errmsg = errmsg         ! intent(out)
        !cdata_block(nb,nt)%errflg = errflg         ! intent(out)
        call ccpp_physics_run(cdata_block(nb,nt), scheme_name="h2ophys", ierr=ierr)
        ! Copy back intent(inout) interstitial variables to local variables in driver
        errmsg = trim(cdata_block(nb,nt)%errmsg)
        errflg = cdata_block(nb,nt)%errflg
! OPTION B END
        if (errflg/=0) then
          write(0,*) 'Error in call to h2phys: ' // trim(errmsg)
          stop
        end if
#else
        if (Model%me==0) write(0,*) 'CCPP DEBUG: calling non-CCPP compliant version of h2ophys'
        call h2ophys (ix, im, levs, levh2o, dtp, Stateout%gq0(1,1,1),  &
                      Stateout%gq0(1,1,1), h2o_pres, Statein%prsl,     &
                      Tbd%h2opl, h2o_coeff, Model%ldiag3d,             &
                      dq3dt_loc(1,1,1), me)
#endif
      endif

!  --- ...  to side-step the ozone physics

!      if (ntrac >= 2) then
!        do k = 1, levs
!          gq0(k,ntoz) = qgrs(k,ntoz)
!        enddo
!      endif

!     if (lprnt) then
!       write(0,*) ' levs=',levs,' jcap=',jcap,' dtp',dtp               &
!    &,  ' slmsk=',slmsk(ilon,ilat),' kdt=',kdt
!       print *,' rann=',rann,' ncld=',ncld,' iq=',iq,' lat=',lat
!       print *,' pgr=',pgr
!       print *,' del=',del(ipr,:)
!       print *,' prsl=',prsl(ipr,:)
!       print *,' prslk=',prslk(ipr,:)
!       print *,' rann=',rann(ipr,1)
!       write(0,*)' gt0=',Stateout%gt0(ipr,:)                    &
!    &,         ' kdt=',kdt,' xlon=',grid%xlon(ipr),' xlat=',grid%xlat(ipr)
!       print *,' dtdt=',dtdt(ipr,:)
!       print *,' gu0=',gu0(ipr,:)
!       print *,' gv0=',gv0(ipr,:)
!       write(0,*) ' gt0=',(gt0(ipr,k),k=1,levs),' kdt=',kdt
!       write(0,*)' gq0=',(gq0(ipr,k,1),k=1,levs),' lat=',lat
!       write(0,*)' gq0i2=',(gq0(ipr,k,ntiw),k=1,levs),' lat=',lat
!       write(0,*)' gq1=',(gq0(ipr,k,ntcw),k=1,levs)
!       print *,' vvel=',vvel
!     endif
!     if (lprnt) write(7000,*)' bef convection gu0=',gu0(ipr,:)
!    &,' lat=',lat,' kdt=',kdt,' me=',me
!     if (lprnt) write(7000,*)' bef convection gv0=',gv0(ipr,:)

#ifdef CCPP
      ! From this point on, dtdt, dudt, dvdt and dqdt have a different meaning:
      ! instead of holding tendencies, they are used as temporary storage for
      ! previous values of t, q, u, v. Reset them to zero here to avoid misguiding
      ! differences between hybrid CCPP and standalone CCPP runs.
      dudt = 0.
      dvdt = 0.
      dtdt = 0.
      dqdt = 0.

      if (Model%me==0) write(0,*) 'CCPP DEBUG: calling GFS_DCNV_generic_pre through option B'
      ! Copy local variables from driver to appropriate interstitial variables
      !Interstitial(nt)%im = im                   ! intent(in) - set in Interstitial(nt)%create()
      !Model%levs                                 ! intent(in)
      !Model%ldiag3d                              ! intent(in)
      !Model%cnvgwd                               ! intent(in)
      !Model%lgocart                              ! intent(in)
      !Stateout%gu0                               ! intent(in)
      !Stateout%gv0                               ! intent(in)
      !Stateout%gt0                               ! intent(in)
      !Stateout%gq0(:,:,1)                        ! intent(in)
      Interstitial(nt)%save_u   = dudt            ! intent(inout)
      Interstitial(nt)%save_v   = dvdt            ! intent(inout)
      Interstitial(nt)%save_t   = dtdt            ! intent(inout)
      Interstitial(nt)%save_q(:,:,1) = dqdt(:,:,1)! intent(inout)
      !cdata_block(nb,nt)%errmsg = errmsg         ! intent(out)
      !cdata_block(nb,nt)%errflg = errflg         ! intent(out)
      call ccpp_physics_run(cdata_block(nb,nt), scheme_name="GFS_DCNV_generic_pre", ierr=ierr)
      ! Copy back intent(inout/out) interstitial variables to local variables in driver
      dudt = Interstitial(nt)%save_u
      dvdt = Interstitial(nt)%save_v
      dtdt = Interstitial(nt)%save_t
      dqdt(:,:,1) = Interstitial(nt)%save_q(:,:,1)
      errmsg = trim(cdata_block(nb,nt)%errmsg)
      errflg = cdata_block(nb,nt)%errflg
      if (errflg/=0) then
          write(0,*) 'Error in call to GFS_DCNV_generic_pre: ' // trim(errmsg)
          stop
      end if
#else
      if (Model%me==0) write(0,*) 'CCPP DEBUG: calling non-CCPP compliant version of GFS_DCNV_generic_pre'
      if (Model%ldiag3d) then
        do k=1,levs
          do i=1,im
            dtdt(i,k) = Stateout%gt0(i,k)
            dudt(i,k) = Stateout%gu0(i,k)
            dvdt(i,k) = Stateout%gv0(i,k)
          enddo
        enddo
      elseif (Model%cnvgwd) then
        dtdt(1:im,:) = Stateout%gt0(1:im,:)
      endif   ! end if_ldiag3d/cnvgwd

      if (Model%ldiag3d .or. Model%lgocart) then
        dqdt(1:im,:,1) = Stateout%gq0(1:im,:,1)
      endif   ! end if_ldiag3d/lgocart
#endif

#ifdef GFS_HYDRO
      call get_phi(im, ix, levs, ntrac, Stateout%gt0, Stateout%gq0,    &
                   Model%thermodyn_id, Model%sfcpress_id,              &
                   Model%gen_coord_hybrid Statein%prsi, Statein%prsik, &
                   Statein%prsl, Statein%prslk, Statein%phii, Statein%phil)
#else
!GFDL   Adjust the height hydrostatically in a way consistent with FV3 discretization
#ifdef CCPP
      if (Model%me==0) write(0,*) 'CCPP DEBUG: calling get_phi_fv3 through option B'
      ! Copy local variables from driver to appropriate interstitial variables
      !Interstitial(nt)%ix                    ! intent(in) - set in Interstitial(nt)%create()
      !Model%levs                             ! intent(in)
      !Stateout%gt0                           ! intent(in)
      !Stateout%gq0(:,:,1)                    ! intent(in)
      Interstitial(nt)%del_gz = del_gz        ! intent(inout)
      !Statein%phii                           ! intent(out)
      !Statein%phil                           ! intent(out)
      !cdata_block(nb,nt)%errmsg = errmsg     ! intent(out)
      !cdata_block(nb,nt)%errflg = errflg     ! intent(out)
      !
      call ccpp_physics_run(cdata_block(nb,nt), scheme_name="get_phi_fv3", ierr=ierr)
      ! Copy back intent(inout) interstitial variables to local variables in driver
      del_gz = Interstitial(nt)%del_gz
      errmsg = trim(cdata_block(nb,nt)%errmsg)
      errflg = cdata_block(nb,nt)%errflg
      if (errflg/=0) then
          write(0,*) 'Error in call to get_phi_fv3: ' // trim(errmsg)
          stop
      end if
#else
      call get_phi_fv3 (ix, levs, ntrac, Stateout%gt0, Stateout%gq0, &
                        del_gz, Statein%phii, Statein%phil)
#endif
#endif

#ifdef CCPP
      if (Model%me==0) write(0,*) 'CCPP DEBUG: calling GFS_suite_interstitial_3 through option B'
      ! Copy local variables from driver to appropriate interstitial variables
      !Interstitial(nt)%im = im                   ! intent(in) - set in Interstitial(nt)%create()
      !Model%levs                                 ! intent(in)
      Interstitial(nt)%nn = nn                    ! intent(in)
      !Model%cscnv                                ! intent(in)
      !Model%satmedmf                             ! intent(in)
      !Model%trans_trac                           ! intent(in)
      !Model%do_shoc                              ! intent(in)
      !Model%ltaerosol                            ! intent(in)
      !Model%ntrac                                ! intent(in)
      !Model%ntcw                                 ! intent(in)
      !Model%ntiw                                 ! intent(in)
      !Model%ntclamt                              ! intent(in)
      !Model%ntrw                                 ! intent(in)
      !Model%ntsw                                 ! intent(in)
      !Model%ntrnc                                ! intent(in)
      !Model%ntsnc                                ! intent(in)
      !Model%ntgl                                 ! intent(in)
      !Model%ntgnc                                ! intent(in)
      !Stateout%gq0                               ! intent(in)
      !Model%imp_physics                          ! intent(in)
      !Model%imp_physics_mg                       ! intent(in)
      !Model%imp_physics_zhao_carr                ! intent(in)
      !Model%imp_physics_zhao_carr_pdf            ! intent(in)
      !Model%imp_physics_gfdl                     ! intent(in)
      !Model%imp_physics_thompson                 ! intent(in)
      !Model%imp_physics_wsm6                     ! intent(in)
      !Statein%prsi                               ! intent(in)
      !Statein%prsl                               ! intent(in)
      !Statein%prslk                              ! intent(in)
      Interstitial(nt)%rhcbot = rhbbot            ! intent(in)
      Interstitial(nt)%rhcpbl = rhpbl             ! intent(in)
      Interstitial(nt)%rhctop = rhbtop            ! intent(in)
      !Model%rhcmax                               ! intent(in)
      Interstitial(nt)%islmsk = islmsk            ! intent(in)
      Interstitial(nt)%work1  = work1             ! intent(in)
      Interstitial(nt)%work2  = work2             ! intent(in)
      Interstitial(nt)%kpbl   = kpbl              ! intent(in)
      !Interstitial(nt)%clw    = clw              ! intent(inout) - clw is only initialized in GFS_typedefs/interstitial_phys_reset at this point in the code (don't overwrite with uninitialized data from local variable)
      !Interstitial(nt)%rhc    = rhc              ! intent(inout) - rhc is only initialized in GFS_typedefs/interstitial_phys_reset at this point in the code (don't overwrite with uninitialized data from local variable)
      !Interstitial(nt)%save_q(:,:,Model%ntcw) = liq0  ! intent(inout) - liq0 is only initialized in GFS_typedefs/interstitial_phys_reset at this point in the code (don't overwrite with uninitialized data from local variable)
      !Interstitial(nt)%save_q(:,:,Model%ntiw) = ice00 ! intent(inout) - ice00 is only initialized in GFS_typedefs/interstitial_phys_reset at this point in the code (don't overwrite with uninitialized data from local variable)
      !cdata_block(nb,nt)%errmsg = errmsg     ! intent(out)
      !cdata_block(nb,nt)%errflg = errflg     ! intent(out)
      call ccpp_physics_run(cdata_block(nb,nt), scheme_name="GFS_suite_interstitial_3", ierr=ierr)
      ! Copy back intent(inout) interstitial variables to local variables in driver
      !GF* The following variables are calculated in GFS_typedefs/interstitial_setup_tracers; they are copied to local vars here,
      !    but are not used in GFS_suite_interstitial_3
      ntk = Interstitial(nt)%ntk
      tottracer = Interstitial(nt)%tracers_total
      otspt = Interstitial(nt)%otspt
      nsamftrac = Interstitial(nt)%nsamftrac
      !*GF
      !GF* The following variables are initialized in GFS_typedefs/interstitial_phys_reset; the are copied to local vars here,
      !    but are not used in GFS_suite_interstitial_3
      if (Model%imfdeepcnv >= 0 .or. Model%imfshalcnv > 0) then
        cnvc = Interstitial(nt)%cnvc
        cnvw = Interstitial(nt)%cnvw
      end if
      ktop = Interstitial(nt)%ktop
      kbot = Interstitial(nt)%kbot
      !*GF
      clw = Interstitial(nt)%clw
      rhc = Interstitial(nt)%rhc
      if (Model%imp_physics == Model%imp_physics_thompson) then
        if (Model%ltaerosol) then
          liq0 = Interstitial(nt)%save_q(:,:,Model%ntcw)
          ice00 = Interstitial(nt)%save_q(:,:,Model%ntiw)
        else
          ice00 = Interstitial(nt)%save_q(:,:,Model%ntiw)
        end if
      end if
      errmsg = trim(cdata_block(nb,nt)%errmsg)
      errflg = cdata_block(nb,nt)%errflg
      if (errflg/=0) then
          write(0,*) 'Error in call to GFS_suite_interstitial_3: ' // trim(errmsg)
          stop
      end if
#else
      do k=1,levs
        do i=1,im
          clw(i,k,1) = 0.0
          clw(i,k,2) = -999.9
        enddo
      enddo

      if(imp_physics == Model%imp_physics_thompson) then
        if(Model%ltaerosol) then
          ice00 (:,:) = 0.0
          liq0  (:,:) = 0.0
        else
          ice00 (:,:) = 0.0
        endif
      endif

!  --- ...  for convective tracer transport (while using ras, csaw, or samf)
!           (the code here implicitly assumes that ntiw=ntcw+1)
      ntk       = 0
      tottracer = 0
      if (Model%cscnv .or. Model%satmedmf .or. Model%trans_trac ) then
        otspt(:,:)   = .true.     ! otspt is used only for cscnv
        otspt(1:3,:) = .false.    ! this is for sp.hum, ice and liquid water
        tracers = 2
        do n=2,ntrac
          if ( n /= ntcw  .and. n /= ntiw  .and. n /= ntclamt .and. &
               n /= ntrw  .and. n /= ntsw  .and. n /= ntrnc   .and. &
               n /= ntsnc .and. n /= ntgl  .and. n /= ntgnc) then
            tracers = tracers + 1
            do k=1,levs
              do i=1,im
                clw(i,k,tracers) = Stateout%gq0(i,k,n)
              enddo
            enddo
            if (ntke  == n ) then
              otspt(tracers+1,1) = .false.
              ntk = tracers
            endif
            if (ntlnc == n .or. ntinc == n .or. ntrnc == n .or. ntsnc == n .or. ntgnc == n)    &
!           if (ntlnc == n .or. ntinc == n .or. ntrnc == n .or. ntsnc == n .or.&
!               ntrw  == n .or. ntsw  == n .or. ntgl  == n)                    &
                    otspt(tracers+1,1) = .false.
          endif
        enddo
        tottracer = tracers - 2
      endif   ! end if_ras or cfscnv or samf

!    if (kdt == 1 .and. me == 0)                                       &
!        write(0,*)' trans_trac=',Model%trans_trac,' tottracer=',      &
!    &               tottracer,' kdt=',kdt,' ntk=',ntk

      do i=1,im
        ktop(i) = 1
        kbot(i) = levs
      enddo

!  --- ...  calling condensation/precipitation processes
!           --------------------------------------------

      if (ntcw > 0) then
        if (imp_physics == Model%imp_physics_mg) then                           ! compute rhc for GMAO macro physics cloud pdf
          do i=1,im
            tx1(i) = 1.0 / Statein%prsi(i,1)
            tx2(i) = 1.0 - rhc_max * work1(i)-rhbbot*work2(i)
          enddo
          do k = 1, levs
            do i = 1, im
              kk   = max(2,kpbl(i))
              tem  = Statein%prsl(i,k) * tx1(i)
              tem1 = min(max((tem-Statein%prsi(i,kk)*tx1(i))/slope_mg, -20.0), 20.0)
              tem2 = min(max((0.3-0.2*abs(cos(Grid%xlat(i)))-tem)/slope_upmg, -20.0), 20.0) ! Anning
              if (islmsk(i) > 0) then
                tem1 = 1.0 / (1.0+exp(tem1+tem1))
              else
                tem1 = 2.0 / (1.0+exp(tem1+tem1))
              endif
              tem2 = 1.0 / (1.0+exp(tem2))

              rhc(i,k) = min(rhc_max, max(0.7, 1.0-tx2(i)*tem1*tem2))
!             rhc(i,k) = min(rhc_max, rhc_max*work1(i) + (1.0-tx2(i)*tem1*tem2)*work2(i))
            enddo
          enddo
        else
          do k=1,levs
            do i=1,im
              kk = max(10,kpbl(i))
              if (k < kk) then
                tem    = rhbbot - (rhbbot-rhpbl) * (1.0-Statein%prslk(i,k))                &
                                                 / (1.0-Statein%prslk(i,kk))
              else
                tem    = rhpbl - (rhpbl-rhbtop) * (Statein%prslk(i,kk)-Statein%prslk(i,k)) &
                                                / Statein%prslk(i,kk)
              endif
              tem      = rhc_max * work1(i) + tem * work2(i)
              rhc(i,k) = max(0.0, min(1.0,tem))
            enddo
          enddo
        endif
      endif      ! ntcw > 0

      if (imp_physics == 99 .or. imp_physics == 98) then   ! zhao-carr microphysics
        do i=1,im
          psautco_l(i) = Model%psautco(1)*work1(i) + Model%psautco(2)*work2(i)
          prautco_l(i) = Model%prautco(1)*work1(i) + Model%prautco(2)*work2(i)
        enddo
        do k=1,levs
          do i=1,im
            clw(i,k,1) = Stateout%gq0(i,k,ntcw)
          enddo
        enddo
      elseif (imp_physics == Model%imp_physics_gfdl) then
        clw(1:im,:,1) = Stateout%gq0(1:im,:,ntcw)
      elseif (imp_physics == Model%imp_physics_thompson) then
        do k=1,levs
          do i=1,im
            clw(i,k,1) = Stateout%gq0(i,k,ntiw)                    ! ice
            clw(i,k,2) = Stateout%gq0(i,k,ntcw)                    ! water
          enddo
        enddo
        if(Model%ltaerosol) then
          ice00(:,:) = clw(:,:,1)
          liq0(:,:)  = clw(:,:,2)
        else
          ice00(:,:) = clw(:,:,1)
        endif
      elseif (imp_physics == Model%imp_physics_wsm6 .or. imp_physics == Model%imp_physics_mg) then
        do k=1,levs
          do i=1,im
            clw(i,k,1) = Stateout%gq0(i,k,ntiw)                    ! ice
            clw(i,k,2) = Stateout%gq0(i,k,ntcw)                    ! water
          enddo
        enddo
      else       ! if_ntcw
        do i=1,im
          psautco_l(i) = Model%psautco(1)*work1(i) + Model%psautco(2)*work2(i)
          prautco_l(i) = Model%prautco(1)*work1(i) + Model%prautco(2)*work2(i)
        enddo
        rhc(:,:) = 1.0
      endif   ! end if_ntcw
#endif

!
!        Call SHOC if do_shoc is true and shocaftcnv is false
!
      if (Model%do_shoc .and. .not. Model%shocaftcnv) then
#ifdef CCPP
        if (Model%me==0) write(0,*) 'CCPP DEBUG: calling shoc through option B'
        ! Copy local variables from driver to appropriate interstitial variables
        !Interstitial(nt)%ix = ix                     ! intent(in) - set in Interstitial(nt)%create()
        !Interstitial(nt)%im = im                     ! intent(in) - set in Interstitial(nt)%create()
        !Model%levs                                   ! intent(in)
        !Model%do_shoc                                ! intent(in)
        !Model%shocaftcnv                             ! intent(in)
        !Model%mg3_as_mg2                             ! intent(in)
        !Model%imp_physics                            ! intent(in)
        !Model%imp_physics_gfdl                       ! intent(in)
        !Model%imp_physics_zhao_carr                  ! intent(in)
        !Model%imp_physics_zhao_carr_pdf              ! intent(in)
        !Model%imp_physics_mg                         ! intent(in)
        !Model%fprcp                                  ! intent(in)
        !Model%tcr                                    ! intent(in)
        !Model%tcrf                                   ! intent(in)
        !con_cp, con_g, con_hvap, con_hfus, con_rv, con_rd, con_pi, con_fvirt ! intent(in) - physical constant in physcons.F90
        !Stateout%gq0(:,:,ntiw)                       ! intent(in)
        !Stateout%gq0(:,:,ntrw)                       ! intent(in)
        !Stateout%gq0(:,:,ntsw)                       ! intent(in)
        !Stateout%gq0(:,:,ntgl)                       ! intent(in)
        !Model%dtp                                    ! intent(in)
        !Model%me                                     ! intent(in)
        !Statein%prsl                                 ! intent(in)
        !Statein%phii                                 ! intent(in)
        !Statein%phil                                 ! intent(in)
        !Stateout%gu0                                 ! intent(in)
        !Stateout%gv0                                 ! intent(in)
        !Statein%vvl                                  ! intent(in)
        Interstitial(nt)%rhc = rhc                    ! intent(in)
        !Model%sup                                    ! intent(in)
        !Model%shoc_parm(1)                           ! intent(in)
        !Model%shoc_parm(2)                           ! intent(in)
        !Model%shoc_parm(3)                           ! intent(in)
        !Model%shoc_parm(4)                           ! intent(in)
        !Model%shoc_parm(5)                           ! intent(in)
        Interstitial(nt)%hflx = hflx                  ! intent(in)
        Interstitial(nt)%evap = evap                  ! intent(in)
        Interstitial(nt)%prnum = prnum                ! intent(in)
        Interstitial(nt)%skip_macro = skip_macro      ! intent(inout)
        Interstitial(nt)%clw(:,:,1) = clw(:,:,1)      ! intent(inout)
        Interstitial(nt)%clw(:,:,2) = clw(:,:,2)      ! intent(inout)
        !Stateout%gq0(:,:,ntcw)                       ! intent(inout)
        !Stateout%gq0(:,:,ntlnc)                      ! intent(inout)
        !Stateout%gq0(:,:,ntinc)                      ! intent(inout)
        !Stateout%gt0                                 ! intent(inout)
        !Stateout%gq0(:,:,1)                          ! intent(inout)
        !Tbd%phy_f3d(:,:,ntot3d-2)                    ! intent(inout)
        Interstitial(nt)%clw(:,:,ntk) = clw(:,:,ntk)  ! intent(inout)
        !Tbd%phy_f3d(:,:,ntot3d-1)                    ! intent(inout)
        !Tbd%phy_f3d(:,:,ntot3d)                      ! intent(inout)
        !cdata_block(nb,nt)%errmsg = errmsg           ! intent(out)
        !cdata_block(nb,nt)%errflg = errflg           ! intent(out)
        call ccpp_physics_run(cdata_block(nb,nt), scheme_name="shoc", ierr=ierr)
        ! Copy back intent(inout) interstitial variables to local variables in driver
        skip_macro = Interstitial(nt)%skip_macro
        clw(:,:,1) = Interstitial(nt)%clw(:,:,1)
        clw(:,:,2) = Interstitial(nt)%clw(:,:,2)
        clw(:,:,ntk) = Interstitial(nt)%clw(:,:,ntk)
        errmsg = trim(cdata_block(nb,nt)%errmsg)
        errflg = cdata_block(nb,nt)%errflg
        if (errflg/=0) then
            write(0,*) 'Error in call to shoc: ' // trim(errmsg)
            stop
        end if
#else
        if (imp_physics == Model%imp_physics_mg) then
          skip_macro = Model%do_shoc
          do k=1,levs
            do i=1,im
              clw(i,k,1) = Stateout%gq0(i,k,ntiw)                    ! ice
              clw(i,k,2) = Stateout%gq0(i,k,ntcw)                    ! water
              ncpl(i,k)  = Stateout%gq0(i,k,ntlnc)
              ncpi(i,k)  = Stateout%gq0(i,k,ntinc)
            enddo
          enddo
          if (abs(Model%fprcp) == 1 .or. mg3_as_mg2) then
            do k=1,levs
              do i=1,im
                qrn(i,k)  = Stateout%gq0(i,k,ntrw)
                qsnw(i,k) = Stateout%gq0(i,k,ntsw)
              enddo
            enddo
          elseif (Model%fprcp > 1) then
            do k=1,levs
              do i=1,im
                qrn(i,k)  = Stateout%gq0(i,k,ntrw)
                qsnw(i,k) = Stateout%gq0(i,k,ntsw) + Stateout%gq0(i,k,ntgl)
                clw(i,k,1) = clw(i,k,1) + Stateout%gq0(i,k,ntgl)
              enddo
            enddo
          endif
        elseif (imp_physics == Model%imp_physics_gfdl) then  ! GFDL MP - needs modify for condensation
          do k=1,levs
            do i=1,im
              clw(i,k,1) = Stateout%gq0(i,k,ntiw)                    ! ice
              clw(i,k,2) = Stateout%gq0(i,k,ntcw)                    ! water
              qrn(i,k)   = Stateout%gq0(i,k,ntrw)
              qsnw(i,k)  = Stateout%gq0(i,k,ntsw)
            enddo
          enddo
        elseif (imp_physics == 99 .or. imp_physics == 98) then
          do k=1,levs
            do i=1,im
              if (abs(Stateout%gq0(i,k,ntcw)) < epsq) then
                Stateout%gq0(i,k,ntcw) = 0.0
              endif
              tem = Stateout%gq0(i,k,ntcw)              &
     &            * max(0.0, MIN(1.0, (TCR-Stateout%gt0(i,k))*TCRF))
              clw(i,k,1) = tem                              ! ice
              clw(i,k,2) = Stateout%gq0(i,k,ntcw) - tem              ! water
            enddo
          enddo
        endif

!  if (lprnt) write(0,*)'gt01=',Stateout%gt0(ipr,:)
!  if (lprnt) write(0,*)'gq01=',Stateout%gq0(ipr,1:60,1)
!  if (lprnt) write(0,*)'clwi=',clw(ipr,1:20,1)
!  if (lprnt) write(0,*)'clwl=',clw(ipr,1:10,2)
!       dtshoc = 60.0
!       dtshoc = 120.0
!       dtshoc = dtp
!       dtshoc = min(dtp, 300.0)
!       nshocm = max(1, nint(dtp/dtshoc))
!       dtshoc = dtp / nshocm
!       do nshoc=1,nshocm
!      if (lprnt) write(0,*)' before shoc tke=',clw(ipr,1:45,ntk), &
!    &' kdt=',kdt,'xlon=',grid%xlon(ipr),' xlat=',grid%xlat(ipr)

!     phy_f3d(1,1,ntot3d-2) - shoc determined sgs clouds
!     phy_f3d(1,1,ntot3d-1) - shoc determined diffusion coefficients
!     phy_f3d(1,1,ntot3d  ) - shoc determined  w'theta'
!
!     dqdt(1:im,:,1) = Stateout%gq0(1:im,:,1)
!     dqdt(1:im,:,2) = Stateout%gq0(1:im,:,ntiw)
!     dqdt(1:im,:,3) = Stateout%gq0(1:im,:,ntcw)
!GFDL lat has no meaning inside of shoc - changed to "1"
!GFDL          call shoc(ix, im, 1, levs, levs+1, dtp, me, lat,
!         call shoc (ix, im, 1, levs, levs+1, dtp, me, 1, Statein%prsl(1,1),  &
!         call shoc (ix, im, 1, levs, levs+1, dtshoc, me, 1, Statein%prsl(1,1),  &
          call shoc (ix, im, 1, levs, levs+1, dtp, me, 1, Statein%prsl(1,1),  &
                     Statein%phii(1,1), Statein%phil(1,1), Stateout%gu0(1,1), &
                     Stateout%gv0(1,1), Statein%vvl(1,1), Stateout%gt0(1,1),  &
                     Stateout%gq0(1,1,1), clw(1,1,1), clw(1,1,2),             &
                     qsnw, qrn, qgl,                                          &
                     rhc, Model%sup, Model%shoc_parm(1), Model%shoc_parm(2),  &
                     Model%shoc_parm(3), Model%shoc_parm(4),                  &
                     Model%shoc_parm(5), Tbd%phy_f3d(1,1,ntot3d-2),           &
                     clw(1,1,ntk), hflx, evap, prnum,                         &
                     Tbd%phy_f3d(1,1,ntot3d-1),                               &
                     Tbd%phy_f3d(1,1,ntot3d), lprnt, ipr, ncpl, ncpi)

!       enddo
          if (imp_physics == Model%imp_physics_mg .and. Model%fprcp > 1) then
            do k=1,levs
              do i=1,im
                clw(i,k,1) = clw(i,k,1) - Stateout%gq0(i,k,ntgl)
              enddo
            enddo
          endif
!     if (lprnt) write(0,*)'aftshocgt0=',Stateout%gt0(ipr,:)
!     if (lprnt) write(0,*)'aftshocgq0=',Stateout%gq0(ipr,1:60,1)
!     if (lprnt) write(0,*)' aft shoc tke=',clw(ipr,1:25,ntk), &
!    &' kdt=',kdt,'xlon=',grid%xlon(ipr),' xlat=',grid%xlat(ipr)
!     if (lprnt) write(0,*)' aftshoccld=',tbd%phy_f3d(ipr,:,ntot3d-2)*100
!     if (lprnt) write(0,*)' aftshocice=',clw(ipr,:,1)
!     if (lprnt) write(0,*)' aftshocwat=',clw(ipr,:,1)
!     write(1000+me,*)' at latitude = ',lat
!     rain1 = 0.0
!     call moist_bud(im,im,ix,levs,me,kdt,con_g,dtp,del,rain1
!    &,              dqdt(1,1,1), dqdt(1,1,2), dqdt(1,1,3)
!    &,              gq0(1,1,1),clw(1,1,2),clw(1,1,1),'shoc      ')
!     tem = 1000.0
!     call moist_bud(im,im,ix,levs,me,kdt,con_g,tem,del,rain1               &
!    &,              dqdt(1,1,1), dqdt(1,1,2), dqdt(1,1,3)                  &
!    &,              Stateout%gq0(1:ix,1:levs,1),clw(1,1,2),clw(1,1,1)      &
!    &,              '   shoc   ', grid%xlon(1:im), grid%xlat(1:im))

          if (ntlnc > 0 .and. ntinc > 0 .and. ncld >= 2) then
            do k=1,levs
              do i=1,im
                Stateout%gq0(i,k,ntlnc) = ncpl(i,k)
                Stateout%gq0(i,k,ntinc) = ncpi(i,k)
              enddo
            enddo
          endif
!       do k=1,levs
!         do i=1,im
!           sgs_cld(i,k) = sgs_cld(i,k) + shoc_cld(i,k)
!         enddo
!       enddo
!     if (lprnt) write(0,*)' gt03=',gt0(ipr,1:10)
!     if (lprnt) write(0,*)' tke=',clw(ipr,1:10,ntk)

!      if (lprnt) write(1000+me,*)' after shoc tke=',clw(1,:,ntk),
!    &' kdt=',kdt
!       enddo
!
!      do k=1,levs
!      write(1000+me,*)' maxcld=',maxval(sgs_cld(1:im,k)),
!      write(1000+me,*)' maxtkh=',maxval(phy_f3d(1:im,k,ntot3d-1)),
!    &' k=',k,' kdt=',kdt,' lat=',lat
!      enddo

!     write(0,*)' aft shoc gt0=',gt0(1,:),' lat=',lat
!     write(0,*)' aft shoc gq0=',gq0(1,:,1),' lat=',lat
!     write(0,*)' aft shoc gu0=',gu0(1,:),' lat=',lat
!
#endif
      endif   ! if(do_shoc)

!
!  --- ...  calling convective parameterization
!           -----------------------------------
      if (Model%do_deep) then

! For CCPP compliant physics, this code is in GFS_DCNV_generic_pre;
! for non-CCPP compliant physics (RAS only, this code is copied to
! before/after rascnv)
#ifndef CCPP
         if (Model%do_ca) then
           do k=1,levs
            do i=1,im
             Stateout%gq0(i,k,1) = Stateout%gq0(i,k,1)*(1.0 + Coupling%ca_deep(i)/500.)
            enddo
           enddo
         endif   

        if (Model%isppt_deep) then
           savet = Stateout%gt0
           saveq = Stateout%gq0(:,:,1)
           saveu = Stateout%gu0
           savev = Stateout%gv0
        endif
#endif

        if (.not. Model%ras .and. .not. Model%cscnv) then

          if (Model%imfdeepcnv == 1) then             ! no random cloud top
            call sascnvn (im, ix, levs, Model%jcap, dtp, del,                    &
                          Statein%prsl, Statein%pgr, Statein%phil, clw(:,:,1:2), &
                          Stateout%gq0(:,:,1), Stateout%gt0, Stateout%gu0,       &
                          Stateout%gv0, cld1d, rain1, kbot, ktop, kcnv,          &
                          islmsk, Statein%vvl, ncld, ud_mf, dd_mf,               &
                          dt_mf, cnvw, cnvc,                                     &
                          QLCN, QICN, w_upi,cf_upi, CNV_MFD,                     &
                          CNV_DQLDT,CLCN,CNV_FICE,CNV_NDROP,CNV_NICE,imp_physics,&
                          Model%clam_deep,   Model%c0s_deep,                     &
                          Model%c1_deep,     Model%betal_deep, Model%betas_deep, &
                          Model%evfact_deep, Model%evfactl_deep,                 &
                          Model%pgcon_deep)
          elseif (Model%imfdeepcnv == 2) then
#ifndef CCPP
            if(.not. Model%satmedmf .and. .not. Model%trans_trac) then
               nsamftrac = 0
            else
               nsamftrac = tottracer
            endif
#endif
#ifdef CCPP
            if (Model%me==0) write(0,*) 'CCPP DEBUG: calling samfdeepcnv through option B'
            ! Copy local variables from driver to appropriate interstitial variables
            !Interstitial(nt)%im = im             ! intent(in) - set in Interstitial(nt)%create()
            !Interstitial(nt)%ix = ix             ! intent(in) - set in Interstitial(nt)%create()
            !Model%levs                           ! intent(in)
            !physical constants                   ! intent(in) - physical constant in physcons.F90
            !Model%dtp                            ! intent(in)
            Interstitial(nt)%ntk = ntk            ! intent(in)
            Interstitial(nt)%nsamftrac = nsamftrac! intent(in)
            Interstitial(nt)%del = del            ! intent(in)
            !Statein%prsl                         ! intent(in)
            !Statein%pgr                          ! intent(in)
            !Statein%phil                         ! intent(in)
            Interstitial(nt)%clw = clw            ! intent(inout)
            !Stateout%gq0(:,:,1)                  ! intent(inout)
            !Stateout%gt0                         ! intent(inout)
            !Stateout%gu0                         ! intent(inout)
            !Stateout%gv0                         ! intent(inout)
            !Model%do_ca                          ! intent(in)
            !Coupling%ca_deep                     ! intent(in)
            Interstitial(nt)%cld1d = cld1d        ! intent(out)
            Interstitial(nt)%raincd = rain1       ! intent(out)
            Interstitial(nt)%kbot = kbot          ! intent(out)
            Interstitial(nt)%ktop = ktop          ! intent(out)
            Interstitial(nt)%kcnv = kcnv          ! intent(out)
            Interstitial(nt)%islmsk = islmsk      ! intent(in)
            !Grid%area                            ! intent(in)
            !Statein%vvl                          ! intent(in)
            !Model%ncld                           ! intent(in)
            Interstitial(nt)%ud_mf     = ud_mf    ! intent(out)
            Interstitial(nt)%dd_mf     = dd_mf    ! intent(out)
            Interstitial(nt)%dt_mf     = dt_mf    ! intent(out)
            Interstitial(nt)%cnvw      = cnvw     ! intent(out)
            Interstitial(nt)%cnvc      = cnvc     ! intent(out)
            Interstitial(nt)%qlcn      = qlcn     ! intent(inout)
            Interstitial(nt)%qicn      = qicn     ! intent(inout)
            Interstitial(nt)%w_upi     = w_upi    ! intent(inout)
            Interstitial(nt)%cf_upi    = cf_upi   ! intent(inout)
            Interstitial(nt)%cnv_mfd   = cnv_mfd  ! intent(inout)
            Interstitial(nt)%cnv_dqldt = cnv_dqldt! intent(inout)
            Interstitial(nt)%clcn      = clcn     ! intent(inout)
            Interstitial(nt)%cnv_fice  = cnv_fice ! intent(inout)
            Interstitial(nt)%cnv_ndrop = cnv_ndrop! intent(inout)
            Interstitial(nt)%cnv_nice  = cnv_nice ! intent(inout)
            !Model%imp_physics                    ! intent(in)
            !Model%imp_physics_mg                 ! intent(in)
            !Model%clam_deep                      ! intent(in)
            !Model%c0s_deep                       ! intent(in)
            !Model%c1_deep                        ! intent(in)
            !Model%betal_deep                     ! intent(in)
            !Model%betas_deep                     ! intent(in)
            !Model%evfact_deep                    ! intent(in)
            !Model%evfactl_deep                   ! intent(in)
            !Model%pgcon_deep                     ! intent(in)
            !Model%asolfac_deep                   ! intent(in)
            !cdata_block(nb,nt)%errmsg = errmsg   ! intent(out)
            !cdata_block(nb,nt)%errflg = errflg   ! intent(out)
            call ccpp_physics_run(cdata_block(nb,nt), scheme_name="samfdeepcnv", ierr=ierr)
            ! Copy back intent(inout) interstitial variables to local variables in driver
            clw       = Interstitial(nt)%clw
            cld1d     = Interstitial(nt)%cld1d
            rain1     = Interstitial(nt)%raincd
            kbot      = Interstitial(nt)%kbot
            ktop      = Interstitial(nt)%ktop
            kcnv      = Interstitial(nt)%kcnv
            ud_mf     = Interstitial(nt)%ud_mf
            dd_mf     = Interstitial(nt)%dd_mf
            dt_mf     = Interstitial(nt)%dt_mf
            cnvw      = Interstitial(nt)%cnvw
            cnvc      = Interstitial(nt)%cnvc
            qlcn      = Interstitial(nt)%qlcn
            qicn      = Interstitial(nt)%qicn
            w_upi     = Interstitial(nt)%w_upi
            cf_upi    = Interstitial(nt)%cf_upi
            cnv_mfd   = Interstitial(nt)%cnv_mfd
            cnv_dqldt = Interstitial(nt)%cnv_dqldt
            clcn      = Interstitial(nt)%clcn
            cnv_fice  = Interstitial(nt)%cnv_fice
            cnv_ndrop = Interstitial(nt)%cnv_ndrop
            cnv_nice  = Interstitial(nt)%cnv_nice
            errmsg = trim(cdata_block(nb,nt)%errmsg)
            errflg = cdata_block(nb,nt)%errflg
            if (errflg/=0) then
                write(0,*) 'Error in call to samfdeepcnv: ' // trim(errmsg)
                stop
            end if
#else
            if (Model%me==0) write(0,*) 'CCPP DEBUG: calling non-CCPP compliant version of samfdeepcnv'
            call samfdeepcnv(im, ix, levs, dtp, ntk, nsamftrac, del,             &
                             Statein%prsl, Statein%pgr, Statein%phil, clw,       &
                             Stateout%gq0(:,:,1), Stateout%gt0,                  &
                             Stateout%gu0, Stateout%gv0, Model%do_ca,            &
                             Coupling%ca_deep, cld1d, rain1, kbot, ktop, kcnv,   &
                             islmsk, garea,                                      &
                             Statein%vvl, ncld, ud_mf, dd_mf, dt_mf, cnvw, cnvc, &
                             QLCN, QICN, w_upi,cf_upi, CNV_MFD,                  &
                             CNV_DQLDT,CLCN,CNV_FICE,CNV_NDROP,CNV_NICE,         &
                             imp_physics,                                        &
                             Model%clam_deep,   Model%c0s_deep,                  &
                             Model%c1_deep,  Model%betal_deep, Model%betas_deep, &
                             Model%evfact_deep, Model%evfactl_deep,              &
                             Model%pgcon_deep,  Model%asolfac_deep)
!           if (lprnt) print *,' rain1=',rain1(ipr)
#endif
          elseif (Model%imfdeepcnv == 3) then
#ifdef CCPP
            if (Model%me==0) write(0,*) 'CCPP DEBUG: calling cu_gf_driver_pre through option B'
            ! Copy local variables from driver to appropriate interstitial variables
            !Model%kdt                              ! intent(in)
            !Model%fhour                            ! intent(in)
            !Model%dtp                              ! intent(in)
            !Statein(nb)%tgrs                       ! intent(in)
            !Statein(nb%qgrs(:,:,1)                 ! intent(in)
            !Tbd(nb%prevst                          ! intent(in)
            !Tbd(nb%prevsq                          ! intent(in)
            !Tbd(nb%forcet                          ! intent(out)
            !Tbd(nb%forceq                          ! intent(out)
            !Tbd(nb)%cactiv                         ! intent(out)
            !Sfcprop(nb)%conv_act                   ! intent(in)
            cdata_block(nb,nt)%errmsg = errmsg      ! intent(out)
            cdata_block(nb,nt)%errflg = errflg      ! intent(out)
            call ccpp_physics_run(cdata_block(nb,nt), scheme_name="cu_gf_driver_pre",ierr=ierr)
            ! Copy back intent(inout) interstitial variables to local variables in driver
            errmsg = trim(cdata_block(nb,nt)%errmsg)
            errflg = cdata_block(nb,nt)%errflg
            if (errflg/=0) then
                write(0,*) 'Error in call to cu_gf_driver_pre: '//trim(errmsg)
                stop
            end if
            if (Model%me==0) write(0,*) 'CCPP DEBUG: calling cu_gf_driver through option B'
            ! Copy local variables from driver to appropriate interstitial variables
            Interstitial(nt)%tracers_total=tottracer        ! intent(in) 
            !Grid%area                                      ! intent(in)
            Interstitial(nt)%im = im                        ! intent(in)
            Interstitial(nt)%ix = ix                        ! intent(in)
            !Model%levs                                     ! intent(in)
            !Model%dtp                                      ! intent(in)
            !Tbd%cactiv                                     ! intent(inout)
            !Tbd%forcet                                     ! intent(in)
            !Tbd%forceq                                     ! intent(in)
            !Statein%phil                                   ! intent(in)
            Interstitial(nt)%raincd = rain1                 ! intent(out)
            !Stateout%gq0                                   ! intent(inout)
            !Stateout%gt0                                   ! intent(inout)
            Interstitial(nt)%cld1d = cld1d                  ! intent(out)
            !Stateout%gu0                                   ! intent(inout)
            !Stateout%gv0                                   ! intent(inout)
            !Statein%tgrs                                   ! intent(in)
            !Statein%vvl                                    ! intent(in)
            !Statein%qgrs                                   ! intent(in)
            !Statein%prsl                                   ! intent(in)
            !Statein%pgr                                    ! intent(in)
            Interstitial(nt)%kbot  = kbot                   ! intent(out)
            Interstitial(nt)%ktop  = ktop                   ! intent(out)
            Interstitial(nt)%kcnv  = kcnv                   ! intent(out)
            Interstitial(nt)%islmsk= islmsk                 ! intent(in)
            Interstitial(nt)%hflx  = hflx                   ! intent(in)
            Interstitial(nt)%evap  = evap                   ! intent(in)
            Interstitial(nt)%clw   = clw                    ! intent(inout)
            !Diag%hpbl                                      ! intent(in)
            Interstitial(nt)%ud_mf = ud_mf                  ! intent(out)
            Interstitial(nt)%dd_mf = dd_mf                  ! intent(out)
            Interstitial(nt)%dt_mf = dt_mf                  ! intent(out)
            Interstitial(nt)%cnvw  = cnvw                   ! intent(out)
            Interstitial(nt)%cnvc  = cnvc                   ! intent(out)
            cdata_block(nb,nt)%errmsg = errmsg              ! intent(out)
            cdata_block(nb,nt)%errflg = errflg              ! intent(out)
            call ccpp_physics_run(cdata_block(nb,nt), scheme_name="cu_gf_driver",ierr=ierr)
            ! Copy back intent(inout) interstitial variables to local variables in driver
            cld1d  = Interstitial(nt)%cld1d
            rain1  = Interstitial(nt)%raincd
            kbot   = Interstitial(nt)%kbot
            ktop   = Interstitial(nt)%ktop
            kcnv   = Interstitial(nt)%kcnv
            clw    = Interstitial(nt)%clw
            ud_mf  = Interstitial(nt)%ud_mf
            dd_mf  = Interstitial(nt)%dd_mf
            dt_mf  = Interstitial(nt)%dt_mf
            cnvw   = Interstitial(nt)%cnvw
            cnvc   = Interstitial(nt)%cnvc
            errmsg = trim(cdata_block(nb,nt)%errmsg)
            errflg = cdata_block(nb,nt)%errflg
            if (errflg/=0) then
                write(0,*) 'Error in call to cu_gf_driver: '// trim(errmsg)
                stop
            end if
#else
            if (Model%me==0) then
                write(0,*) 'Error, GF convection scheme only available through CCPP'
                stop
            end if
#endif
          elseif (Model%imfdeepcnv == 4) then
#ifdef CCPP
            if (Model%me==0) write(0,*) 'CCPP DEBUG: calling cu_ntiedtke_pre through option B'
            ! Copy local variables from driver to appropriate interstitial variables
            !Model%kdt                              ! intent(in)
            !Model%fhour                            ! intent(in)
            !Model%dtp                              ! intent(in)
            !Statein(nb)%tgrs                       ! intent(in)
            !Statein(nb%qgrs(:,:,1)                 ! intent(in)
            !Tbd(nb%prevst                          ! intent(in)
            !Tbd(nb%prevsq                          ! intent(in)
            !Tbd(nb%forcet                          ! intent(out)
            !Tbd(nb%forceq                          ! intent(out)
            cdata_block(nb,nt)%errmsg = errmsg      ! intent(out)
            cdata_block(nb,nt)%errflg = errflg      ! intent(out)
            call ccpp_physics_run(cdata_block(nb,nt), scheme_name="cu_ntiedtke_pre",ierr=ierr)
            ! Copy back intent(inout) interstitial variables to local variables in driver
            errmsg = trim(cdata_block(nb,nt)%errmsg)
            errflg = cdata_block(nb,nt)%errflg
            if (errflg/=0) then
                write(0,*) 'Error in call to cu_ntiedtke_pre: '//trim(errmsg)
                stop
            end if
            if (Model%me==0) write(0,*) 'CCPP DEBUG: calling cu_ntiedtke through option B'
            ! Copy local variables from driver to appropriate interstitial variables
            !Stateout%gu0                                   ! intent(inout)
            !Stateout%gv0                                   ! intent(inout)
            !Stateout%gt0                                   ! intent(inout)
            !Stateout%gq0(:,:,Model%ntqv)                   ! intent(inout)
            !Tbd%forceq                                     ! intent(in   )
            !Tbd%forcet                                     ! intent(in   )
            Interstitial(nt)%clw           = clw            ! intent(inout)
            !Statein%phil                                   ! intent(in   )
            !Statein%phii                                   ! intent(in   )
            !Statein%prsl                                   ! intent(in   )
            !Statein%prsi                                   ! intent(in   )
            !Statein%vvl                                    ! intent(in   )
            Interstitial(nt)%evap          = evap           ! intent(in   )
            Interstitial(nt)%hflx          = hflx           ! intent(in   )
            Interstitial(nt)%raincd        = rain1          ! intent(out  )
            Interstitial(nt)%islmsk        = islmsk         ! intent(in   )
            Interstitial(nt)%im            = im             ! intent(in   )
            Interstitial(nt)%ix            = ix             ! intent(in   )
            !Model%levs                                     ! intent(in   )
            !Model%dtp                                      ! intent(in   )
            !Grid%dx                                        ! intent(in   )
            Interstitial(nt)%kbot          = kbot           ! intent(out  )
            Interstitial(nt)%ktop          = ktop           ! intent(out  )
            Interstitial(nt)%kcnv          = kcnv           ! intent(out  )
            Interstitial(nt)%tracers_total = tottracer      ! intent(in   )
            Interstitial(nt)%ud_mf         = ud_mf          ! intent(out  )
            Interstitial(nt)%dd_mf         = dd_mf          ! intent(out  )
            Interstitial(nt)%dt_mf         = dt_mf          ! intent(out  )
            Interstitial(nt)%cnvw          = cnvw           ! intent(out  )
            Interstitial(nt)%cnvc          = cnvc           ! intent(out  )
            cdata_block(nb,nt)%errmsg      = errmsg         ! intent(out  )
            cdata_block(nb,nt)%errflg      = errflg         ! intent(out  )
            call ccpp_physics_run(cdata_block(nb,nt), scheme_name="cu_ntiedtke",ierr=ierr)
            ! Copy back intent(inout) interstitial variables to local variables in driver
            clw    = Interstitial(nt)%clw
            rain1  = Interstitial(nt)%raincd
            kbot   = Interstitial(nt)%kbot
            ktop   = Interstitial(nt)%ktop
            kcnv   = Interstitial(nt)%kcnv
            ud_mf  = Interstitial(nt)%ud_mf
            dd_mf  = Interstitial(nt)%dd_mf
            dt_mf  = Interstitial(nt)%dt_mf
            cnvw   = Interstitial(nt)%cnvw
            cnvc   = Interstitial(nt)%cnvc
            errmsg = trim(cdata_block(nb,nt)%errmsg)
            errflg = cdata_block(nb,nt)%errflg
            if (errflg/=0) then
                write(0,*) 'Error in call to cu_ntiedtke: '// trim(errmsg)
                stop
            end if
#else
            if (Model%me==0) then
                write(0,*) 'Error, New Tiedtke convection scheme only available through CCPP'
                stop
            end if
#endif

          elseif (Model%imfdeepcnv == 0) then         ! random cloud top
            call sascnv (im, ix, levs, Model%jcap, dtp, del,                     &
                         Statein%prsl, Statein%pgr, Statein%phil, clw(:,:,1:2),  &
                         Stateout%gq0(:,:,1), Stateout%gt0, Stateout%gu0,        &
                         Stateout%gv0, cld1d, rain1, kbot, ktop, kcnv,           &
                         islmsk, Statein%vvl, Tbd%rann, ncld,                    &
                         ud_mf, dd_mf, dt_mf, cnvw, cnvc,                        &
                         QLCN, QICN, w_upi,cf_upi, CNV_MFD,                      &
                         CNV_DQLDT,CLCN,CNV_FICE,CNV_NDROP,CNV_NICE,imp_physics )
!           if (lprnt) print *,' rain1=',rain1(ipr),' rann=',rann(ipr,1)
          endif

!
#ifndef CCPP
          ! For CCPP, this is in GFS_DCNV_generic_post
          if(Model%do_ca) then
            Coupling%cape(:)=cld1d(:)
          endif
          if (Model%npdf3d == 3 .and. Model%num_p3d == 4) then
            do k=1,levs
              do i=1,im
                Tbd%phy_f3d(i,k,num2) = cnvw(i,k)
                Tbd%phy_f3d(i,k,num3) = cnvc(i,k)
                cnvw(i,k)             = 0.0
                cnvc(i,k)             = 0.0
              enddo
            enddo
          elseif (Model%npdf3d == 0 .and. Model%ncnvcld3d == 1) then
            do k=1,levs
              do i=1,im
                Tbd%phy_f3d(i,k,num2) = cnvw(i,k)
                cnvw(i,k)             = 0.0
              enddo
            enddo
          endif
#endif

        else        ! ras or cscnv
! cs_conv_pre and ras_pre
#ifndef CCPP
          fscav = 0
#else
          cld1d = 0
#endif
          if (Model%cscnv) then    ! Chikira-Sugiyama  convection scheme (via CSU)
#ifdef CCPP
! OPTION B - works with all compilers
            if (Model%me==0) write(0,*) 'CCPP DEBUG: calling cs_conv_pre through option B'
            ! Copy local variables from driver to appropriate interstitial variables
            !Interstitial(nt)%im              = im                ! intent(in) - set in Interstitial(nt)%create()
            !Model%levs                                           ! intent(in)
            !IPD_Contral%ntrac
            !Stateout%gq0(:,:,1)
            Interstitial(nt)%clw(:,:,1)                = clw(:,:,1)        ! intent(in)
            Interstitial(nt)%clw(:,:,2)                = clw(:,:,2)        ! intent(in)
            Interstitial(nt)%work1                     = work1             ! intent(in)
            Interstitial(nt)%work2                     = work2             ! intent(in)
            !Model%cs_parm(1)                                              ! intent(in)
            !Model%cs_parm(2)                                              ! intent(in)
            Interstitial(nt)%wcbmax                    = wcbmax            ! intent(out)
            Interstitial(nt)%fswtr                     = fswtr             ! intent(out)
            Interstitial(nt)%fscav                     = fscav             ! intent(out)
            Interstitial(nt)%save_q(:,:,1)             = dqdt(:,:,1)       ! intent(out)
            Interstitial(nt)%save_q(:,:,ntcw)          = dqdt(:,:,2)       ! intent(out) ntcw=2
            Interstitial(nt)%save_q(:,:,ntiw)          = dqdt(:,:,3)       ! intent(out) ntiw=3
            !cdata_block(nb,nt)%errmsg        = errmsg            ! intent(out)
            !cdata_block(nb,nt)%errflg        = errflg            ! intent(out)
            call ccpp_physics_run(cdata_block(nb,nt), scheme_name="cs_conv_pre", ierr=ierr)
            ! Copy intent(inout) and intent(out) interstitial variables to local variables in driver
            wcbmax        = Interstitial(nt)%wcbmax                   ! intent(out)
            fswtr         = Interstitial(nt)%fswtr                    ! intent(out)
            fscav         = Interstitial(nt)%fscav                    ! intent(out)
            dqdt(:,:,1)   = Interstitial(nt)%save_q(:,:,1)            ! intent(out)
            dqdt(:,:,2)   = Interstitial(nt)%save_q(:,:,ntcw)
            dqdt(:,:,3)   = Interstitial(nt)%save_q(:,:,ntiw)
            errmsg        = trim(cdata_block(nb,nt)%errmsg)
            errflg        = cdata_block(nb,nt)%errflg
            if (errflg/=0) then
                write(0,*) 'Error in call to cs_conv_pre_mp_cs_conv_pre: ' // trim(errmsg)
                stop
            end if

#else
         if (Model%me==0) write(0,*) 'CCPP DEBUG: calling non-CCPP compliant version of cs_conv_pre'
          fscav(:) = 0.0
          fswtr(:) = 0.0
          do i = 1,im
              wcbmax(i)  = Model%cs_parm(1) * work1(i) + Model%cs_parm(2) * work2(i)
          enddo
!     write(0,*)' bef cs_cconv phii=',phii(ipr,:)
!    &,' sizefsc=',size(fscav)
!     write(0,*)' bef cs_cconv otspt=',otspt,' kdt=',kdt,' me=',me
!           do k=1,levs
!             do i=1,im
!               dqdt(i,k,1) = Stateout%gq0(i,k,1)
!               dqdt(i,k,2) = clw(i,k,2)
!               dqdt(i,k,3) = clw(i,k,1)
!             enddo
!           enddo

#endif

#ifdef CCPP
            if (Model%me==0) write(0,*) 'CCPP DEBUG: calling cs_conv through option B'
            ! Copy local variables from driver to appropriate interstitial variables
            !Interstitial(nt)%ix       = ix                       ! intent(in) - set in Interstitial(nt)%create()
            !Interstitial(nt)%im       = im                       ! intent(in) - set in Interstitial(nt)%create()
            !Model%levs                                           ! intent(in)
            !Model%ntracp1 = Model%ntrac + 1                      ! intent(in) - set in Model%create()
            !Interstitial(nt)%nn = nn                             ! intent(in) - set in Interstitial(nt)%create()
            !Interstitial(nt)%ncstrac = tottracer+3               ! intent(in) - set in Interstitial(nt)%create()
            !Model%nctp                                           ! intent(in)
            Interstitial(nt)%otspt = otspt                        ! intent(in)
            ! lat = 1                                             ! intent(in)
            !Model%kdt                                            ! intent(in)
            !Stateout%gt0                                         ! intent(inout)
            !Stateout%gq0(:,:,1:1)                                ! intent(inout)
            Interstitial(nt)%raincd    = rain1                    ! intent(out)
            Interstitial(nt)%clw   = clw                          ! intent(inout)
            !Statein%phil                                         ! intent(in)
            !Statein%phii                                         ! intent(in)
            !Statein%prsl                                         ! intent(in)
            !Statein%prsi                                         ! intent(in)
            !Model%dtp                                            ! intent(in)
            !Model%dtf                                            ! intent(in)
!!
!! JLS NOTE:  The convective mass fluxes (dt_mf, dd_mf and ud_mf) passed in and out of cs_conv have not been multiplied by
!!            the timestep (kg/m2/sec) as they are in all other convective schemes.  EMC is aware of this problem, 
!!            and in the future will be fixing this discrepancy.  In the meantime, CCPP will use the same mass flux standard_name
!!            and long_name as the other convective schemes, where the units are in kg/m2. (Aug 2018)
!!
            Interstitial(nt)%ud_mf     = ud_mf                    ! intent(inout)
            Interstitial(nt)%dd_mf     = dd_mf                    ! intent(inout)
            Interstitial(nt)%dt_mf     = dt_mf                    ! intent(inout)
            !Stateout%gu0                                         ! intent(inout)
            !Stateout%gv0                                         ! intent(inout)
            Interstitial(nt)%fscav = fscav                        ! intent(in)
            Interstitial(nt)%fswtr = fswtr                        ! intent(in)
            !Tbd%phy_fctd                                         ! intent(inout)
            !IPD_Control%me                                       ! intent(in)
            Interstitial(nt)%wcbmax    = wcbmax                   ! intent(in) - set in cs_conv_pre
            !Model%cs_parm(3)                                     ! intent(in)
            !Model%cs_parm(4)                                     ! intent(in)
            !Model%cs_parm(9)                                     ! intent(in)
            Interstitial(nt)%sigmatot = sigmatot                  ! intent(out)
            !Model%do_aw                                          ! intent(in)
            !Model%do_awdd                                        ! intent(in)
            !Model%flx_form                                       ! intent(in)
            !IPD_Control%lprnt                                    ! intent(in)
            Interstitial(nt)%ipr       = ipr                      ! intent(in)
            Interstitial(nt)%kcnv      = kcnv                     ! intent(inout)
            Interstitial(nt)%qlcn      = qlcn                     ! intent(out)
            Interstitial(nt)%qicn      = qicn                     ! intent(out)
            Interstitial(nt)%w_upi     = w_upi                    ! intent(out)
            Interstitial(nt)%cf_upi    = cf_upi                   ! intent(out)
            Interstitial(nt)%cnv_mfd   = cnv_mfd                  ! intent(out)
            Interstitial(nt)%cnv_dqldt = cnv_dqldt                ! intent(out)
            Interstitial(nt)%clcn      = clcn                     ! intent(out)
            Interstitial(nt)%cnv_fice  = cnv_fice                 ! intent(out)
            Interstitial(nt)%cnv_ndrop = cnv_ndrop                ! intent(out)
            Interstitial(nt)%cnv_nice  = cnv_nice                 ! intent(out)
            !Model%imp_physics                                    ! intent(in)
            !cdata_block(nb,nt)%errmsg = errmsg                   ! intent(out)
            !cdata_block(nb,nt)%errflg = errflg                   ! intent(out)
            call ccpp_physics_run(cdata_block(nb,nt), scheme_name="cs_conv", ierr=ierr)
            ! Copy intent(inout) and intent(out) interstitial variables to local variables in driver
            rain1     = Interstitial(nt)%raincd
            clw       = Interstitial(nt)%clw
            ud_mf     = Interstitial(nt)%ud_mf
            dd_mf     = Interstitial(nt)%dd_mf
            dt_mf     = Interstitial(nt)%dt_mf
            sigmatot  = Interstitial(nt)%sigmatot
            kcnv      = Interstitial(nt)%kcnv
            qlcn      = Interstitial(nt)%qlcn
            qicn      = Interstitial(nt)%qicn
            w_upi     = Interstitial(nt)%w_upi
            cf_upi    = Interstitial(nt)%cf_upi
            cnv_mfd   = Interstitial(nt)%cnv_mfd
            cnv_dqldt = Interstitial(nt)%cnv_dqldt
            clcn      = Interstitial(nt)%clcn
            cnv_fice  = Interstitial(nt)%cnv_fice
            cnv_ndrop = Interstitial(nt)%cnv_ndrop
            cnv_nice  = Interstitial(nt)%cnv_nice
            errmsg    = trim(cdata_block(nb,nt)%errmsg)
            errflg    = cdata_block(nb,nt)%errflg
            if (errflg/=0) then
                write(0,*) 'Error in call to cs_conv_mp_cs_conv: ' // trim(errmsg)
                stop
            end if

#else
            if (Model%me==0) write(0,*) 'CCPP DEBUG: calling non-CCPP compliant version of cs_convr'

!           if (lprnt) write(0,*)'befcsgt0=',Stateout%gt0(ipr,:)
!           if (lprnt) write(0,*)'befcstke=',clw(ipr,1:25,ntk)

! NOTE:  The variable rain1 output from cs_convr (called prec inside the subroutine) is a precipitation flux (kg/m2/sec),
!         not meters LWE like the other schemes.  It is converted to m after the call to cs_convr.
              

            call cs_convr (ix, im, levs, Model%ntrac+1, nn,                 &
                           tottracer+3, Model%nctp,                         &
!                           otspt(1:tottracer+3,1:2), 1,                     &
                           otspt(1:Model%ntrac+1,1:2), 1,                   &
                           kdt, Stateout%gt0, Stateout%gq0(:,:,1:1), rain1, &
                           clw, Statein%phil, Statein%phii, Statein%prsl,   &
                           Statein%prsi, dtp, dtf, ud_mf, dd_mf, dt_mf,     &
                           Stateout%gu0, Stateout%gv0, fscav, fswtr,        &
                           Tbd%phy_fctd, me, wcbmax, Model%cs_parm(3),      &
                           Model%cs_parm(4), Model%cs_parm(9), sigmatot,    &
!                          Model%cs_parm(4), sigmai, sigmatot, vverti,      &
                           Model%do_aw, Model%do_awdd, Model%flx_form,      &
                           lprnt, ipr, kcnv, QLCN, QICN,                    &
                           w_upi, cf_upi, CNV_MFD,           CNV_DQLDT,     &
                           CLCN, CNV_FICE, CNV_NDROP, CNV_NICE, imp_physics)

!           if (lprnt) write(0,*)'aftcsgt0=',Stateout%gt0(ipr,:)
!           if (lprnt) write(0,*)'aftcstke=',clw(ipr,1:25,ntk)

!     write(1000+me,*)' at latitude = ',lat
!     call moist_bud(im,im,ix,levs,me,kdt,con_g,dtp,del,rain1
!    &,                    dqdt(1,1,1), dqdt(1,1,2), dqdt(1,1,3)
!    &,                    gq0(1,1,1),clw(1,1,2),clw(1,1,1),' cs_conv')
!     tem = 1000.0
!     call moist_bud(im,im,ix,levs,me,kdt,con_g,tem,del,rain1               &
!    &,              dqdt(1,1,1), dqdt(1,1,2), dqdt(1,1,3)                  &
!    &,              Stateout%gq0(1:ix,1:levs,1),clw(1,1,2),clw(1,1,1)      &
!    &,              '   cs_conv', grid%xlon(1:im), grid%xlat(1:im))

              rain1(:) = rain1(:) * (dtp*0.001)
#endif

#ifdef CCPP
            if (Model%me==0) write(0,*) 'CCPP DEBUG: calling cs_conv_post through option B'
            ! Copy local variables from driver to appropriate interstitial
            ! variables
            !Interstitial(nt)%im       = im                       ! intent(in) - set in Interstitial(nt)%create()
            !Model%levs                                           ! intent(in)
            !Model%do_aw                                          ! intent(in)
            Interstitial(nt)%sigmatot  = sigmatot                 ! intent(in)
            Interstitial(nt)%sigmafrac = sigmafrac                ! intent(out)
            !cdata_block(nb,nt)%errmsg = errmsg                   ! intent(out)
            !cdata_block(nb,nt)%errflg = errflg                   ! intent(out)
            call ccpp_physics_run(cdata_block(nb,nt), scheme_name="cs_conv_post", ierr=ierr)
            ! Copy intent(inout) and intent(out) interstitial variables to local
            ! variables in driver
            sigmafrac = Interstitial(nt)%sigmafrac
            errmsg    = trim(cdata_block(nb,nt)%errmsg)
            errflg    = cdata_block(nb,nt)%errflg
            if (errflg/=0) then
                write(0,*) 'Error in call to cs_conv_post_mp_cs_conv_post: ' // trim(errmsg)
                stop
            end if
#else
              if (Model%me==0) write(0,*) 'CCPP DEBUG: calling non-CCPP compliant version of cs_conv_post'
              if (Model%do_aw) then
                do k=1,levs
                  kk = min(k+1,levs)  ! assuming no cloud top reaches the model top
                  do i=1,im                                               !DD
                    sigmafrac(i,k) = 0.5 * (sigmatot(i,k)+sigmatot(i,kk))
                  enddo
                enddo
              endif

!             if (lprnt) then
!               write(0,*)' gt01=',stateout%gt0(ipr,:),' kdt=',kdt
!               write(0,*)' gq01=',stateout%gq0(ipr,:,1),' kdt=',kdt
!               write(0,*)' clw1=',clw(ipr,:,1),' kdt=',kdt
!               write(0,*)' clw2=',clw(ipr,:,1),' kdt=',kdt
!               write(0,*)' aft cs rain1=',rain1(ipr)*86400
!               write(0,*)' aft cs rain1=',rain1(ipr)
!            endif
! end of ifdef CCPP for cs_conv
#endif

          else      ! ras version 2

! For CCPP compliant physics, this code is in GFS_DCNV_generic_pre;
! for the non-CCPP compliant RAS, this must be executed before the
! call to RAS in CCPP (hybrid) mode
#ifdef CCPP
            if (Model%do_ca) then
              do k=1,levs
               do i=1,im
                Stateout%gq0(i,k,1) = Stateout%gq0(i,k,1)*(1.0 + Coupling%ca_deep(i)/500.)
               enddo
              enddo
            endif

            if (Model%isppt_deep) then
               savet = Stateout%gt0
               saveq = Stateout%gq0(:,:,1)
               saveu = Stateout%gu0
               savev = Stateout%gv0
            endif
#endif

            if (Model%ccwf(1) >= 0.0 .or. Model%ccwf(2) >= 0) then
              do i=1,im
                ccwfac(i)  = Model%ccwf(1)*work1(i)    + Model%ccwf(2)*work2(i)
                dlqfac(i)  = Model%dlqf(1)*work1(i)    + Model%dlqf(2)*work2(i)
                psaur_l(i) = Model%psauras(1)*work1(i) + Model%psauras(2)*work2(i)
                praur_l(i) = Model%prauras(1)*work1(i) + Model%prauras(2)*work2(i)
              enddo
            else
              do i=1,im
                ccwfac(i) = -999.0
                dlqfac(i) = 0.0
              enddo
            endif
!           if  (lprnt) write(0,*) ' calling ras for kdt=',kdt,' me=',me    &
!    &,                            ' lprnt=',lprnt,' ccwfac=',ccwfac(ipr)

!           do k=1,levs
!             do i=1,im
!               dqdt(i,k,1) = Stateout%gq0(i,k,1)
!               dqdt(i,k,2) = clw(i,k,2)
!               dqdt(i,k,3) = clw(i,k,1)
!             enddo
!           enddo

            revap = .true.
!           if (ncld ==2) revap = .false.
            trcmin(:)     = -999999.0
            if (ntk-2 > 0) trcmin(ntk-2) = 1.0e-4

!           if (lprnt) write(0,*)' gt04bras=',Stateout%gt0(ipr,1:60)
!           if (lprnt) write(0,*)' gq04bras=',Stateout%gq0(ipr,1:60,1)
!           if (lprnt) write(0,*)'befrastke=',clw(ipr,1:25,ntk)
!           if (lprnt) write(0,*)'trcmin=',trcmin(ntk-2),' ntk=',ntk

            call rascnv (im, ix, levs, dtp, dtf, Tbd%rann, Stateout%gt0,     &
                         Stateout%gq0, Stateout%gu0, Stateout%gv0, clw,      &
                         tottracer, fscav, Statein%prsi, Statein%prsl,       &
                         Statein%prsik, Statein%prslk, Statein%phil,         &
                         Statein%phii, kpbl, cd, rain1, kbot, ktop, kcnv,    &
                         Tbd%phy_f2d(1,Model%num_p2d), Model%flipv, pa2mb,   &
                         me, garea, ccwfac, Model%nrcm, rhc, ud_mf,          &
                         dd_mf, dt_mf, praur_l, Model%wminras(1),            &
                         psaur_l, Model%wminras(2), dlqfac,                  &
                         lprnt, ipr, kdt, revap, QLCN,                       &
                         QICN, w_upi, cf_upi, CNV_MFD,           CNV_DQLDT,  &
!                        QICN, w_upi, cf_upi, CNV_MFD, CNV_PRC3, CNV_DQLDT,  &
                         CLCN, CNV_FICE, CNV_NDROP, CNV_NICE, imp_physics,   &
!                        trcmin)
                         trcmin, ntk)

!          if (lprnt) write(0,*)' gt04=',Stateout%gt0(ipr,1:60)
!          if (lprnt) write(0,*)' gq04=',Stateout%gq0(ipr,1:60,1)
!          if (lprnt) write(0,*)'aftrastke=',clw(ipr,1:25,ntk)

! For CCPP compliant physics, this code is in GFS_DCNV_generic_post;
! for the non-CCPP compliant RAS, this must be executed after the
! call to RAS in CCPP (hybrid) mode
#ifdef CCPP
            cld1d = 0

!            if (Model%ldiag3d .or. Model%lgocart) then
!              do k=1,levs
!                do i=1,im
!                  Coupling%upd_mfi(i,k) = 0.
!                  Coupling%dwn_mfi(i,k) = 0.
!                  Coupling%det_mfi(i,k) = 0.
!                enddo
!              enddo
!            endif
!            if (Model%lgocart) then
!              do k=1,levs
!                do i=1,im
!                  Coupling%dqdti(i,k)  = 0.
!                  Coupling%cnvqci(i,k) = 0.
!                enddo
!              enddo
!            endif

!            if (Model%lgocart) then
!              do k=1,levs
!                do i=1,im
!                  Coupling%upd_mfi(i,k)  = Coupling%upd_mfi(i,k)  + ud_mf(i,k) * frain
!                  Coupling%dwn_mfi(i,k)  = Coupling%dwn_mfi(i,k)  + dd_mf(i,k) * frain
!                  Coupling%det_mfi(i,k)  = Coupling%det_mfi(i,k)  + dt_mf(i,k) * frain
!                  Coupling%cnvqci (i,k)  = Coupling%cnvqci (i,k)  + (clw(i,k,1)+clw(i,k,2) - &
!                                          Stateout%gq0(i,k,ntcw)) * frain
!                enddo
!              enddo
!            endif ! if (lgocart)

            if(Model%isppt_deep)then
               Coupling%tconvtend = Stateout%gt0 - savet
               Coupling%qconvtend = Stateout%gq0(:,:,1) - saveq
               Coupling%uconvtend = Stateout%gu0 - saveu
               Coupling%vconvtend = Stateout%gv0 - savev
            endif
#endif

          endif

!     write(1000+me,*)' at latitude = ',lat
!     tem = 1000.0
!     call moist_bud(im,im,ix,levs,me,kdt,con_g,tem,del,rain1               &
!    &,              dqdt(1,1,1), dqdt(1,1,2), dqdt(1,1,3)                  &
!    &,              Stateout%gq0(1:ix,1:levs,1),clw(1,1,2),clw(1,1,1)      &
!    &,              '  ras_conv', grid%xlon(1:im), grid%xlat(1:im))
!     if(lprnt) write(0,*)' after ras rain1=',rain1(ipr),' me=',me,' kdt=',kdt
!    &,' cnv_prc3sum=',sum(cnv_prc3(ipr,1:levs))
!     if (lprnt) write(0,*)' gt04=',gt0(ipr,1:10)

! For CCPP compliant physics, this code is partially in GFS_DCNV_generic_post;
! for non-CCPP compliant physics (RAS only, this code is copied to before/after rascnv)
#ifndef CCPP
          cld1d = 0

!          if (Model%ldiag3d .or. Model%lgocart) then
!            do k=1,levs
!              do i=1,im
!                Coupling%upd_mfi(i,k) = 0.
!                Coupling%dwn_mfi(i,k) = 0.
!                Coupling%det_mfi(i,k) = 0.
!              enddo
!            enddo
!          endif
!          if (Model%lgocart) then
!            do k=1,levs
!              do i=1,im
!                Coupling%dqdti(i,k)  = 0.
!                Coupling%cnvqci(i,k) = 0.
!              enddo
!            enddo
!          endif

!          if (Model%lgocart) then
!            do k=1,levs
!              do i=1,im
!                Coupling%upd_mfi(i,k)  = Coupling%upd_mfi(i,k)  + ud_mf(i,k) * frain
!                Coupling%dwn_mfi(i,k)  = Coupling%dwn_mfi(i,k)  + dd_mf(i,k) * frain
!                Coupling%det_mfi(i,k)  = Coupling%det_mfi(i,k)  + dt_mf(i,k) * frain
!                Coupling%cnvqci (i,k)  = Coupling%cnvqci (i,k)  + (clw(i,k,1)+clw(i,k,2) - &
!                                        Stateout%gq0(i,k,ntcw)) * frain
!              enddo
!            enddo
!          endif ! if (lgocart)
#endif

        endif   ! end if_not_ras

! For CCPP compliant physics, this code is partially in GFS_DCNV_generic_post;
! for non-CCPP compliant physics (RAS only, this code is copied to before/after rascnv)
#ifndef CCPP
        if(Model%isppt_deep)then
           Coupling%tconvtend = Stateout%gt0 - savet
           Coupling%qconvtend = Stateout%gq0(:,:,1) - saveq
           Coupling%uconvtend = Stateout%gu0 - saveu
           Coupling%vconvtend = Stateout%gv0 - savev
        endif
#endif

      else      ! no parameterized deep convection
        cld1d = 0.
        rain1 = 0.
        ud_mf = 0.
        dd_mf = 0.
        dt_mf = 0.
      endif

!     if (lprnt) then
!       write(0,*)' aftcnvgt0=',stateout%gt0(ipr,:),' kdt=',kdt
!       write(0,*)' aftcnvgq0=',(stateout%gq0(ipr,k,1),k=1,levs)
!       write(0,*)' gq0i2=',(stateout%gq0(ipr,k,ntiw),k=1,levs)
!       write(0,*)' aftcnvgq1=',(stateout%gq0(ipr,k,ntcw),k=1,levs)
!     endif
!
#ifdef CCPP
      if (Model%me==0) write(0,*) 'CCPP DEBUG: calling GFS_DCNV_generic_post through option B'
      ! Copy local variables from driver to appropriate interstitial variables
      !Interstitial(nt)%im = im                   ! intent(in) - set in Interstitial(nt)%create()
      !Model%levs                                 ! intent(in)
      !Model%lssav                                ! intent(in)
      !Model%ldiag3d                              ! intent(in)
      !Model%lgocart                              ! intent(in)
      !Model%ras                                  ! intent(in)
      !Model%cscnv                                ! intent(in)
      Interstitial(nt)%frain = frain              ! intent(in)
      Interstitial(nt)%raincd = rain1             ! intent(in)
      !Model%dtf                                  ! intent(in)
      Interstitial(nt)%cld1d = cld1d              ! intent(in)
      Interstitial(nt)%save_u = dudt              ! intent(in)
      Interstitial(nt)%save_v = dvdt              ! intent(in)
      Interstitial(nt)%save_t = dtdt              ! intent(in)
      Interstitial(nt)%save_q(:,:,1) = dqdt(:,:,1)! intent(in)
      !Stateout%gu0                               ! intent(in)
      !Stateout%gv0                               ! intent(in)
      !Stateout%gt0                               ! intent(in)
      !Stateout%gq0(:,:,1)                        ! intent(in)
      Interstitial(nt)%ud_mf = ud_mf              ! intent(in)
      Interstitial(nt)%dd_mf = dd_mf              ! intent(in)
      Interstitial(nt)%dt_mf = dt_mf              ! intent(in)
      !con_g                                      ! intent(in)
      Interstitial(nt)%clw(:,:,1) = clw(:,:,1)    ! intent(in)
      Interstitial(nt)%clw(:,:,2) = clw(:,:,2)    ! intent(in)
      !Model%npdf3d                               ! intent(in)
      !Model%num_p3d                              ! intent(in)
      !Model%ncnvcld3d                            ! intent(in)
      !Diag%rainc                                 ! intent(inout)
      !Diag%cldwrk                                ! intent(inout)
      !Diag%cnvprcp                               ! intent(inout)
      !Diag%cnvprcpb                              ! intent(inout)
      !Diag%dt3dt(:,:,4)                          ! intent(inout)
      !Diag%dq3dt(:,:,2)                          ! intent(inout)
      !Diag%du3dt(:,:,3)                          ! intent(inout)
      !Diag%dv3dt(:,:,3)                          ! intent(inout)
      !Diag%upd_mf                                ! intent(inout)
      !Diag%dwn_mf                                ! intent(inout)
      !Diag%det_mf                                ! intent(inout)
      !Coupling%dqdti                             ! intent(inout)
      !Coupling%cnvqci                            ! intent(inout)
      !Coupling%upd_mfi                           ! intent(inout)
      !Coupling%dwn_mfi                           ! intent(inout)
      !Coupling%det_mfi                           ! intent(inout)
      if (allocated(cnvc) .and. allocated(cnvw)) then
        Interstitial(nt)%cnvw = cnvw              ! intent(inout)
        Interstitial(nt)%cnvc = cnvc              ! intent(inout)
      end if
      !Tbd%phy_f3d(:,:,Model%ncnvw)               ! intent(inout)
      !Tbd%phy_f3d(:,:,Model%ncnvw+1)             ! intent(inout)
      !cdata_block(nb,nt)%errmsg = errmsg         ! intent(out)
      !cdata_block(nb,nt)%errflg = errflg         ! intent(out)
      call ccpp_physics_run(cdata_block(nb,nt), scheme_name="GFS_DCNV_generic_post", ierr=ierr)
      ! Copy back intent(inout/out) interstitial variables to local variables in driver
      if (allocated(cnvc) .and. allocated(cnvw)) then
        cnvw = Interstitial(nt)%cnvw
        cnvc = Interstitial(nt)%cnvc
      end if
      errmsg = trim(cdata_block(nb,nt)%errmsg)
      errflg = cdata_block(nb,nt)%errflg
      if (errflg/=0) then
          write(0,*) 'Error in call to GFS_DCNV_generic_post: ' // trim(errmsg)
          stop
      end if
#else
      if (Model%me==0) write(0,*) 'CCPP DEBUG: calling non-CCPP compliant version of GFS_DCNV_generic_post'
      do i=1,im
        Diag%rainc(i) = frain * rain1(i)
      enddo
!
      if (Model%lssav) then
        do i=1,im
          Diag%cldwrk (i)  = Diag%cldwrk (i)  + cld1d(i) * dtf
          Diag%cnvprcp(i)  = Diag%cnvprcp(i)  + Diag%rainc(i)
          Diag%cnvprcpb(i) = Diag%cnvprcpb(i) + Diag%rainc(i)
        enddo

        if (Model%ldiag3d) then
          do k=1,levs
            do i=1,im
              Diag%dt3dt(i,k,4) = Diag%dt3dt(i,k,4) + (Stateout%gt0(i,k)-dtdt(i,k)) * frain
!              Diag%dq3dt(i,k,2) = Diag%dq3dt(i,k,2) + (Stateout%gq0(i,k,1)-dqdt(i,k,1)) * frain
              Diag%du3dt(i,k,3) = Diag%du3dt(i,k,3) + (Stateout%gu0(i,k)-dudt(i,k)) * frain
              Diag%dv3dt(i,k,3) = Diag%dv3dt(i,k,3) + (Stateout%gv0(i,k)-dvdt(i,k)) * frain

!              Diag%upd_mf(i,k)  = Diag%upd_mf(i,k)  + ud_mf(i,k) * (con_g*frain)
!              Diag%dwn_mf(i,k)  = Diag%dwn_mf(i,k)  + dd_mf(i,k) * (con_g*frain)
!              Diag%det_mf(i,k)  = Diag%det_mf(i,k)  + dt_mf(i,k) * (con_g*frain)
            enddo
          enddo
        endif ! if (ldiag3d)

      endif   ! end if_lssav
!
!       update dqdt_v to include moisture tendency due to deep convection
!      if (Model%lgocart) then
!        do k=1,levs
!          do i=1,im
!            Coupling%dqdti  (i,k) = (Stateout%gq0(i,k,1)  - dqdt(i,k,1)) * frain
!            Coupling%upd_mfi(i,k) = Coupling%upd_mfi(i,k) + ud_mf(i,k)   * frain
!            Coupling%dwn_mfi(i,k) = Coupling%dwn_mfi(i,k) + dd_mf(i,k)   * frain
!            Coupling%det_mfi(i,k) = Coupling%det_mfi(i,k) + dt_mf(i,k)   * frain
!            Coupling%cnvqci (i,k) = Coupling%cnvqci (i,k) + (clw(i,k,1)+clw(i,k,2))*frain
!          enddo
!        enddo
!      endif ! if (lgocart)
#endif

!     if (lprnt) write(7000,*)' bef cnvgwd gu0=',gu0(ipr,:)
!    &,' lat=',lat,' kdt=',kdt,' me=',me
!     if (lprnt) write(7000,*)' bef cnvgwd gv0=',gv0(ipr,:)
!
!----------------Convective gravity wave drag parameterization starting --------

      if (Model%cnvgwd) then         !        call convective gravity wave drag

!  --- ...  calculate maximum convective heating rate
!           cuhr = temperature change due to deep convection


#ifdef CCPP
! OPTION B BEGIN
      if (Model%me==0) write(0,*) 'CCPP DEBUG: calling gwdc_pre through option B'
      ! Copy local variables from driver to appropriate interstitial variables
      !Interstitial(nt)%im = im            ! intent(in) - set in Interstitial(nt)%create()
      Interstitial(nt)%work1 = work1       ! intent(in)
      Interstitial(nt)%work2 = work2       ! intent(in)
      Interstitial(nt)%dlength = dlength   ! intent(out)
      Interstitial(nt)%cldf = cldf         ! intent(out)
      Interstitial(nt)%kbot = kbot         ! intent(in)
      Interstitial(nt)%ktop = ktop         ! intent(in)
      Interstitial(nt)%save_t = dtdt       ! intent(in)
      Interstitial(nt)%del = del           ! intent(in)
      Interstitial(nt)%cumabs = cumabs     ! intent(out)
      !cdata_block(nb,nt)%errmsg = errmsg  ! intent(out)
      !cdata_block(nb,nt)%errflg = errflg  ! intent(out)
      call ccpp_physics_run(cdata_block(nb,nt), scheme_name="gwdc_pre", ierr=ierr)
      ! Copy back intent(inout) interstitial variables to local variables in driver
      dlength = Interstitial(nt)%dlength
      cldf = Interstitial(nt)%cldf
      cumabs = Interstitial(nt)%cumabs
      errmsg = trim(cdata_block(nb,nt)%errmsg)
      errflg = cdata_block(nb,nt)%errflg
! OPTION B END
      if (errflg/=0) then
          write(0,*) 'Error in call to gwdc_pre: ' // trim(errmsg)
          stop
      end if
#else
      if (Model%me==0) write(0,*) 'CCPP DEBUG: calling non-CCPP compliant version of gwdc_pre'

        do i=1,im
          cumabs(i) = 0.0
          work3 (i)  = 0.0
        enddo
        do k=1,levs
          do i=1,im
            if (k >= kbot(i) .and. k <= ktop(i)) then
              cumabs(i) = cumabs(i) + (Stateout%gt0(i,k)-dtdt(i,k)) * del(i,k)
              work3(i)  = work3(i)  + del(i,k)
            endif
          enddo
        enddo
        do i=1,im
          if (work3(i) > 0.0) cumabs(i) = cumabs(i) / (dtp*work3(i))
        enddo

#endif
! DH* 20180817 - note: the above non-CCPP code modifies work3, which until then was defined
! as the ratio of the exner function between midlayer and interface at lowest model layer:
!    work3(i) = Statein%prsik(i,1) / Statein%prslk(i,1)
! This does not happen for the CCPP code, because gwdc_pre uses an internal array
! work3 (maybe not a good name, given that we have work1/2/3 in GFS_physics_driver and
! in the IPD_Interstitial DDT). Therefore, work3 is different from here on until the end
! of GFS_physics_driver. This is ok as long as Model%lgocart is set to .false. - if
! Model%lgocart is set to .true., sfc_diag is called again, which uses work3 as input.
! This work3 used in sfc_diag should be the ratio of the exner function, not the modified
! value derived in the non-CCPP code above. If we get different results for the surface
! diagnstics with Model%lgocart=.true., then the CCPP code is correct! *DH 20180817

!       do i = 1, im
!         do k = kbot(i), ktop(i)
!           do k1 = kbot(i), k
!             cumchr(i,k) = cuhr(i,k1) + cumchr(i,k)
!           enddo
!           cumchr(i,k) = cumchr(i,k) / cumabs(i)
!         enddo
!       enddo

!  --- ...  begin check print ******************************************

!       if (lprnt) then
!         if (kbot(ipr) <= ktop(ipr)) then
!           write(*,*) 'kbot <= ktop     for (lat,lon) = ',             &
!    &            xlon(ipr)*57.29578,xlat(ipr)*57.29578
!           write(*,*) 'kcnv kbot ktop dlength  ',kcnv(ipr),       &
!    &            kbot(ipr),ktop(ipr),dlength(ipr)
!           write(*,9000) kdt
!9000       format(/,3x,'k',5x,'cuhr(k)',4x,'cumchr(k)',5x,             &
!    &            'at kdt = ',i4,/)

!           do k = ktop(ipr), kbot(ipr),-1
!             write(*,9010) k,(86400.*cuhr(ipr,k)),(100.*cumchr(ipr,k))
!9010         format(2x,i2,2x,f8.2,5x,f6.0)
!           enddo
!         endif

!         if (fhour >= fhourpr) then
!           print *,' before gwdc in gbphys start print'
!           write(*,*) 'fhour ix im levs = ',fhour,ix,im,levs
!           print *,'dtp  dtf  = ',dtp,dtf

!           write(*,9100)
!9100       format(//,14x,'pressure levels',//                          &
!    &             ' ilev',7x,'prsi',8x,'prsl',8x,'delp',/)

!           k = levs + 1
!           write(*,9110) k,(10.*prsi(ipr,k))
!9110       format(i4,2x,f10.3)

!           do k = levs, 1, -1
!             write(*,9120) k,(10.*prsl(ipr,k)),(10.*del(ipr,k))
!             write(*,9110) k,(10.*prsi(ipr,k))
!           enddo
!9120       format(i4,12x,2(2x,f10.3))

!           write(*,9130)
!9130       format(//,10x,'before gwdc in gbphys',//,' ilev',6x,        &
!    &             'ugrs',9x,'gu0',8x,'vgrs',9x,'gv0',8x,               &
!    &             'tgrs',9x,'gt0',8x,'gt0b',8x,'dudt',8x,'dvdt',/)

!           do k = levs, 1, -1
!             write(*,9140) k,ugrs(ipr,k),gu0(ipr,k),                   &
!    &                        vgrs(ipr,k),gv0(ipr,k),                   &
!    &                        tgrs(ipr,k),gt0(ipr,k),dtdt(ipr,k),       &
!    &                        dudt(ipr,k),dvdt(ipr,k)
!           enddo
!9140       format(i4,9(2x,f10.3))

!           print *,' before gwdc in gbphys end print'
!         endif
!       endif   ! end if_lprnt

!  --- ...  end check print ********************************************

#ifdef CCPP
! OPTION B BEGIN
      if (Model%me==0) write(0,*) 'CCPP DEBUG: calling gwdc through option B'
      ! Copy local variables from driver to appropriate interstitial variables
      !Interstitial(nt)%im       = im            ! intent(in) - set in Interstitial(nt)%create()
      !Interstitial(nt)%ix       = ix            ! intent(in) - set in Interstitial(nt)%create()
      Interstitial(nt)%del      = del            ! intent(in)
      Interstitial(nt)%cumabs   = cumabs         ! intent(in)
      Interstitial(nt)%ktop     = ktop           ! intent(in)
      Interstitial(nt)%kbot     = kbot           ! intent(in)
      Interstitial(nt)%kcnv     = kcnv           ! intent(in)
      Interstitial(nt)%cldf     = cldf           ! intent(in)
      Interstitial(nt)%dlength  = dlength        ! intent(in)
      Interstitial(nt)%ipr      = ipr            ! intent(in)
      Interstitial(nt)%gwdcu    = gwdcu          ! intent(out)
      Interstitial(nt)%gwdcv    = gwdcv          ! intent(out)
      Interstitial(nt)%dusfcg   = dusfcg         ! intent(out)
      Interstitial(nt)%dvsfcg   = dvsfcg         ! intent(out)
      !cdata_block(nb,nt)%errmsg = errmsg        ! intent(out)
      !cdata_block(nb,nt)%errflg = errflg        ! intent(out)
      call ccpp_physics_run(cdata_block(nb,nt), scheme_name="gwdc", ierr=ierr)
      ! Copy back intent(inout) interstitial variables to local variables in driver
      gwdcu  = Interstitial(nt)%gwdcu
      gwdcv  = Interstitial(nt)%gwdcv
      dusfcg = Interstitial(nt)%dusfcg
      dvsfcg = Interstitial(nt)%dvsfcg
      errmsg = trim(cdata_block(nb,nt)%errmsg)
      errflg = cdata_block(nb,nt)%errflg
! OPTION B END
      if (errflg/=0) then
          write(0,*) 'Error in call to gwdc: ' // trim(errmsg)
          stop
      end if
#else
      if (Model%me==0) write(0,*) 'CCPP DEBUG: calling non-CCPP compliant version of gwdc'


!GFDL replacing lat with "1"
!       call gwdc(im, ix, im, levs, lat, gu0, gv0, gt0, gq0, dtp,       &
       call gwdc (im, ix, im, levs, 1, Statein%ugrs, Statein%vgrs,     &
                   Statein%tgrs, Statein%qgrs, dtp, Statein%prsl,       &
                   Statein%prsi, del, cumabs, ktop, kbot, kcnv, cldf,   &
                   con_g, con_cp, con_rd, con_fvirt, con_pi, dlength,   &
                   lprnt, ipr, Model%fhour, gwdcu, gwdcv, dusfcg, dvsfcg)
#endif

!       if (lprnt) then
!         if (fhour >= fhourpr) then
!           print *,' after gwdc in gbphys start print'

!           write(*,9131)
!9131       format(//,10x,'after gwdc in gbphys',//,' ilev',6x,         &
!    &             'ugrs',9x,'gu0',8x,'vgrs',9x,'gv0',8x,               &
!    &             'tgrs',9x,'gt0',8x,'gt0b',7x,'gwdcu',7x,'gwdcv',/)

!           do k = levs, 1, -1
!             write(*,9141) k,ugrs(ipr,k),gu0(ipr,k),                   &
!    &                        vgrs(ipr,k),gv0(ipr,k),                   &
!    &                        tgrs(ipr,k),gt0(ipr,k),dtdt(ipr,k),       &
!    &                        gwdcu(ipr,k),gwdcv(ipr,k)
!           enddo
!9141       format(i4,9(2x,f10.3))

!           print *,' after gwdc in gbphys end print'
!         endif
!       endif

!  --- ...  write out cloud top stress and wind tendencies

#ifdef CCPP
! OPTION B BEGIN
      if (Model%me==0) write(0,*) 'CCPP DEBUG: calling gwdc_post through option B'
      ! Copy local variables from driver to appropriate interstitial variables
      !Interstitial(nt)%im       = im        ! intent(in) - set in Interstitial(nt)%create()
      Interstitial(nt)%dusfcg   = dusfcg     ! intent(in)
      Interstitial(nt)%dvsfcg   = dvsfcg     ! intent(in)
      Interstitial(nt)%gwdcu    = gwdcu      ! intent(in)
      Interstitial(nt)%gwdcv    = gwdcv      ! intent(in)
      !cdata_block(nb,nt)%errmsg = errmsg    ! intent(out)
      !cdata_block(nb,nt)%errflg = errflg    ! intent(out)
      call ccpp_physics_run(cdata_block(nb,nt), scheme_name="gwdc_post", ierr=ierr)
      ! Copy back intent(inout) interstitial variables to local variables in driver
      errmsg = trim(cdata_block(nb,nt)%errmsg)
      errflg = cdata_block(nb,nt)%errflg
! OPTION B END
      if (errflg/=0) then
          write(0,*) 'Error in call to gwdc_post: ' // trim(errmsg)
          stop
      end if
#else
      if (Model%me==0) write(0,*) 'CCPP DEBUG: calling non-CCPP compliant version of gwdc_post'

        if (Model%lssav) then
          do i=1,im
            Diag%dugwd(i) = Diag%dugwd(i) + dusfcg(i)*dtf
            Diag%dvgwd(i) = Diag%dvgwd(i) + dvsfcg(i)*dtf
          enddo

          if (Model%ldiag3d) then
            do k=1,levs
              do i=1,im
                Diag%du3dt(i,k,4) = Diag%du3dt(i,k,4) + gwdcu(i,k)  * dtf
                Diag%dv3dt(i,k,4) = Diag%dv3dt(i,k,4) + gwdcv(i,k)  * dtf
              enddo
            enddo
          endif
        endif   ! end if_lssav

!  --- ...  update the wind components with  gwdc tendencies

        do k=1,levs
          do i=1,im
            eng0               = 0.5*(Stateout%gu0(i,k)*Stateout%gu0(i,k)+Stateout%gv0(i,k)*Stateout%gv0(i,k))
            Stateout%gu0(i,k)  = Stateout%gu0(i,k) + gwdcu(i,k) * dtp
            Stateout%gv0(i,k)  = Stateout%gv0(i,k) + gwdcv(i,k) * dtp
            eng1               = 0.5*(Stateout%gu0(i,k)*Stateout%gu0(i,k)+Stateout%gv0(i,k)*Stateout%gv0(i,k))
            Stateout%gt0(i,k)  = Stateout%gt0(i,k) + (eng0-eng1)/(dtp*con_cp)
          enddo
!         if (lprnt) write(7000,*)' gu0=',gu0(ipr,k),' gwdcu=',
!    &gwdcu(ipr,k), ' gv0=', gv0(ipr,k),' gwdcv=',gwdcv(ipr,k)
!    &,' k=',k
        enddo
#endif

!       if (lprnt) then
!         if (fhour >= fhourpr) then
!           print *,' after tendency gwdc in gbphys start print'

!           write(*,9132)
!9132       format(//,10x,'after tendency gwdc in gbphys',//,' ilev',6x,&
!    &             'ugrs',9x,'gu0',8x,'vgrs',9x,'gv0',8x,               &
!    &             'tgrs',9x,'gt0',8x,'gt0b',7x,'gwdcu',7x,'gwdcv',/)

!           do k = levs, 1, -1
!             write(*,9142) k,ugrs(ipr,k),gu0(ipr,k),vgrs(ipr,k),       &
!    &              gv0(ipr,k),tgrs(ipr,k),gt0(ipr,k),dtdt(ipr,k),      &
!    &              gwdcu(ipr,k),gwdcv(ipr,k)
!           enddo
!9142       format(i4,9(2x,f10.3))

!           print *,' after tendency gwdc in gbphys end print'
!         endif
!       endif

      endif   ! end if_cnvgwd (convective gravity wave drag)

!     if (lprnt) write(7000,*)' aft cnvgwd gu0=',gu0(ipr,:)
!     if (lprnt) write(7000,*)' aft cnvgwd gv0=',gv0(ipr,:)
!    &,' lat=',lat,' kdt=',kdt,' me=',me
!----------------Convective gravity wave drag parameterization over --------

#ifdef CCPP
      if (Model%me==0) write(0,*) 'CCPP DEBUG: calling GFS_SCNV_generic_pre through option B'
      ! Copy local variables from driver to appropriate interstitial variables
      !Interstitial(nt)%im = im                   ! intent(in) - set in Interstitial(nt)%create()
      !Model%levs                                 ! intent(in)
      !Model%ldiag3d                              ! intent(in)
      !Model%lgocart                              ! intent(in)
      !Stateout%gt0                               ! intent(in)
      !Stateout%gq0(:,:,1)                        ! intent(in)
      Interstitial(nt)%save_t = dtdt              ! intent(inout)
      Interstitial(nt)%save_q(:,:,1) = dqdt(:,:,1)! intent(inout)
      cdata_block(nb,nt)%errmsg   = errmsg        ! intent(out)
      cdata_block(nb,nt)%errflg   = errflg        ! intent(out)
      call ccpp_physics_run(cdata_block(nb,nt), scheme_name="GFS_SCNV_generic_pre", ierr=ierr)
      ! Copy back intent(inout) interstitial variables to local variables in driver
      dtdt = Interstitial(nt)%save_t
      dqdt(:,:,1) = Interstitial(nt)%save_q(:,:,1)
      errmsg = trim(cdata_block(nb,nt)%errmsg)
      errflg = cdata_block(nb,nt)%errflg
      if (errflg/=0) then
          write(0,*) 'Error in call to GFS_SCNV_generic_pre: ' // trim(errmsg)
          stop
      end if
#else
      if (Model%me==0) write(0,*) 'CCPP DEBUG: calling non-CCPP compliant version of GFS_SCNV_generic_pre'
      if (Model%ldiag3d) then
        do k=1,levs
          do i=1,im
            dtdt(i,k)   = Stateout%gt0(i,k)
          enddo
        enddo
      endif
!      if (Model%ldiag3d .or. Model%lgocart) then
!        do k=1,levs
!          do i=1,im
!            dqdt(i,k,1) = Stateout%gq0(i,k,1)
!          enddo
!        enddo
!      endif
#endif


!     write(0,*)' before do_shoc shal clstp=',clstp,' kdt=',kdt,
!    &         ' lat=',lat
!           if (lprnt) then
!             write(0,*)' prsl=',prsl(ipr,:)
!             write(0,*) ' del=',del(ipr,:)
!             write(0,*) ' befshalgt0=',gt0(ipr,:),' kdt=',kdt
!             write(0,*) ' befshalgq0=',gt0(ipr,:),' kdt=',kdt
!             write(0,*) ' befshalgq0=',gq0(ipr,:,1),' kdt=',kdt
!             write(0,*) ' befshalgqw=',gq0(ipr,:,3),' kdt=',kdt
!           endif

      if (.not. Model%do_shoc) then

        if (Model%shal_cnv) then               ! Shallow convection parameterizations
!                                               --------------------------------------
          if (Model%imfshalcnv == 1) then      ! opr option now at 2014
                                               !-----------------------
            call shalcnv (im, ix, levs, Model%jcap, dtp, del, Statein%prsl, &
                          Statein%pgr, Statein%phil, clw, Stateout%gq0,     &
                          Stateout%gt0, Stateout%gu0, Stateout%gv0, rain1,  &
                          kbot, ktop, kcnv, islmsk, Statein%vvl, ncld,      &
                          Diag%hpbl, hflx, evap, ud_mf, dt_mf, cnvw, cnvc,  &
                          Model%clam_shal, Model%c0s_shal, Model%c1_shal,   &
                          Model%pgcon_shal)

            do i=1,im
              raincs(i)     = frain * rain1(i)
              Diag%rainc(i) = Diag%rainc(i) + raincs(i)
            enddo
            if (Model%lssav) then
              do i=1,im
                Diag%cnvprcp(i)  = Diag%cnvprcp(i)  + raincs(i)
                Diag%cnvprcpb(i) = Diag%cnvprcpb(i) + raincs(i)
              enddo
            endif
! in shalcnv,  'cnvw' and 'cnvc' are not set to zero
            if (Model%shcnvcw .and. Model%num_p3d == 4 .and. Model%npdf3d == 3) then
              do k=1,levs
                do i=1,im
                  Tbd%phy_f3d(i,k,num2) = Tbd%phy_f3d(i,k,num2) + cnvw(i,k)
                  Tbd%phy_f3d(i,k,num3) = Tbd%phy_f3d(i,k,num3) + cnvc(i,k)
                enddo
              enddo
            elseif (Model%npdf3d == 0 .and. Model%ncnvcld3d == 1) then
              do k=1,levs
                do i=1,im
                  Tbd%phy_f3d(i,k,num2) = Tbd%phy_f3d(i,k,num2) + cnvw(i,k)
                enddo
              enddo
            endif

          elseif (Model%imfshalcnv == 2) then
#ifndef CCPP
            if(.not. Model%satmedmf .and. .not. Model%trans_trac) then
               nsamftrac = 0
            else
               nsamftrac = tottracer
            endif
#endif
#ifdef CCPP
            if (Model%me==0) write(0,*) 'CCPP DEBUG: calling samfshalcnv through option B'
            ! Copy local variables from driver to appropriate interstitial variables
            !Interstitial(nt)%im = im             ! intent(in) - set in Interstitial(nt)%create()
            !Interstitial(nt)%ix = ix             ! intent(in) - set in Interstitial(nt)%create()
            !Model%levs                           ! intent(in)
            !physical constants                   ! intent(in) - physical constant in physcons.F90
            !Model%dtp                            ! intent(in)
            Interstitial(nt)%ntk = ntk            ! intent(in)
            Interstitial(nt)%nsamftrac = nsamftrac! intent(in)
            Interstitial(nt)%del = del            ! intent(in)
            !Statein%prsl                         ! intent(in)
            !Statein%pgr                          ! intent(in)
            !Statein%phil                         ! intent(in)
            Interstitial(nt)%clw = clw            ! intent(inout)
            !Stateout%gq0(:,:,1)                  ! intent(inout)
            !Stateout%gt0                         ! intent(inout)
            !Stateout%gu0                         ! intent(inout)
            !Stateout%gv0                         ! intent(inout)
            Interstitial(nt)%raincs = rain1       ! intent(out)
            Interstitial(nt)%kbot = kbot          ! intent(out)
            Interstitial(nt)%ktop = ktop          ! intent(out)
            Interstitial(nt)%kcnv = kcnv          ! intent(out)
            Interstitial(nt)%islmsk = islmsk      ! intent(in)
            !Grid%area                            ! intent(in)
            !Statein%vvl                          ! intent(in)
            !Model%ncld                           ! intent(in)
            !Diag%hpbl                            ! intent(in)
            Interstitial(nt)%ud_mf = ud_mf        ! intent(out)
            Interstitial(nt)%dt_mf = dt_mf        ! intent(out)
            Interstitial(nt)%cnvw = cnvw          ! intent(out)
            Interstitial(nt)%cnvc = cnvc          ! intent(out)
            !Model%clam_shal                      ! intent(in)
            !Model%c0s_shal                       ! intent(in)
            !Model%c1_shal                        ! intent(in)
            !Model%pgcon_shal                     ! intent(in)
            !Model%asolfac_shal                   ! intent(in)
            !cdata_block(nb,nt)%errmsg = errmsg   ! intent(out)
            !cdata_block(nb,nt)%errflg = errflg   ! intent(out)
            call ccpp_physics_run(cdata_block(nb,nt), scheme_name="samfshalcnv", ierr=ierr)
            ! Copy back intent(inout) interstitial variables to local variables in driver
            clw    = Interstitial(nt)%clw
            rain1  = Interstitial(nt)%raincs
            kbot   = Interstitial(nt)%kbot
            ktop   = Interstitial(nt)%ktop
            kcnv   = Interstitial(nt)%kcnv
            ud_mf  = Interstitial(nt)%ud_mf
            dt_mf  = Interstitial(nt)%dt_mf
            cnvw   = Interstitial(nt)%cnvw
            cnvc   = Interstitial(nt)%cnvc
            errmsg = trim(cdata_block(nb,nt)%errmsg)
            errflg = cdata_block(nb,nt)%errflg
            if (errflg/=0) then
                write(0,*) 'Error in call to samfshalcnv: ' // trim(errmsg)
                stop
            end if
#else
            if (Model%me==0) write(0,*) 'CCPP DEBUG: calling non-CCPP compliant version of samfshalcnv'
            call samfshalcnv (im, ix, levs, dtp, ntk, nsamftrac, del,          &
                              Statein%prsl, Statein%pgr, Statein%phil, clw,    &
                              Stateout%gq0(:,:,1), Stateout%gt0,               &
                              Stateout%gu0, Stateout%gv0,                      &
                              rain1, kbot, ktop, kcnv, islmsk, garea,          &
                              Statein%vvl, ncld, Diag%hpbl, ud_mf,             &
                              dt_mf, cnvw, cnvc,                               &
                              Model%clam_shal,  Model%c0s_shal, Model%c1_shal, &
                              Model%pgcon_shal, Model%asolfac_shal)
#endif
#ifdef CCPP
            if (Model%me==0) write(0,*) 'CCPP DEBUG: calling samfshalcnv_post through option B'
            ! Copy local variables from driver to appropriate interstitial variables
            !Interstitial(nt)%im = im                   ! intent(in) - set in Interstitial(nt)%create()
            !Model%levs                                 ! intent(in)
            !Model%lssav                                ! intent(in)
            !Model%shcnvcw                              ! intent(in)
            Interstitial(nt)%frain = frain              ! intent(in)
            Interstitial(nt)%raincs = rain1             ! intent(in)
            !Model%npdf3d                               ! intent(in)
            !Model%num_p3d                              ! intent(in)
            !Model%ncnvcld3d                            ! intent(in)
            Interstitial(nt)%cnvc = cnvc                ! intent(in)
            Interstitial(nt)%cnvw = cnvw                ! intent(in)
            !Diag%rainc                                 ! intent(inout)
            !Diag%cnvprcp                               ! intent(inout)
            !Diag%cnvprcpb                              ! intent(inout)
            !Tbd%phy_f3d(:,:,Model%ncnvw)               ! intent(inout)
            !Tbd%phy_f3d(:,:,Model%ncnvw+1)             ! intent(inout)
            !cdata_block(nb,nt)%errmsg = errmsg         ! intent(out)
            !cdata_block(nb,nt)%errflg = errflg         ! intent(out)
            call ccpp_physics_run(cdata_block(nb,nt), scheme_name="samfshalcnv_post", ierr=ierr)
            ! Copy back intent(inout) interstitial variables to local variables in driver
            errmsg = trim(cdata_block(nb,nt)%errmsg)
            errflg = cdata_block(nb,nt)%errflg
            if (errflg/=0) then
                write(0,*) 'Error in call to samfshalcnv_post: ' // trim(errmsg)
                stop
            end if
#else
            if (Model%me==0) write(0,*) 'CCPP DEBUG: calling non-CCPP compliant version of samfshalcnv_post'

            do i=1,im
              raincs(i)     = frain * rain1(i)
              Diag%rainc(i) = Diag%rainc(i) + raincs(i)
            enddo
            if (Model%lssav) then
              do i=1,im
                Diag%cnvprcp(i)  = Diag%cnvprcp(i)  + raincs(i)
                Diag%cnvprcpb(i) = Diag%cnvprcpb(i) + raincs(i)
              enddo
            endif
! in  mfshalcnv,  'cnvw' and 'cnvc' are set to zero before computation starts:
            if (Model%shcnvcw .and. Model%num_p3d == 4 .and. Model%npdf3d == 3) then
              do k=1,levs
                do i=1,im
                  Tbd%phy_f3d(i,k,num2) = Tbd%phy_f3d(i,k,num2) + cnvw(i,k)
                  Tbd%phy_f3d(i,k,num3) = Tbd%phy_f3d(i,k,num3) + cnvc(i,k)
                enddo
              enddo
            elseif (Model%npdf3d == 0 .and. Model%ncnvcld3d == 1) then
              do k=1,levs
                do i=1,im
                  Tbd%phy_f3d(i,k,num2) = Tbd%phy_f3d(i,k,num2) +  cnvw(i,k)
                enddo
              enddo
            endif
#endif

          !elseif (Model%imfshalcnv == 3) then
          !if (Model%me==0) write(0,*) "CCPP DEBUG: shallow convection of GF is called in gf_driver"

          !elseif (Model%imfshalcnv == 4) then
          !if (Model%me==0) write(0,*) "CCPP DEBUG: shallow convection of New Tiedtke is called in cu_tiedtke"

          elseif (Model%imfshalcnv == 0) then    ! modified Tiedtke Shallow convecton
                                                 !-----------------------------------
            levshc(:) = 0
            do k=2,levs
              do i=1,im
                dpshc = 0.3 * Statein%prsi(i,1)
                if (Statein%prsi(i,1)-Statein%prsi(i,k) <= dpshc) levshc(i) = k
              enddo
            enddo
            levshcm = 1
            do i=1,im
              levshcm = max(levshcm, levshc(i))
            enddo

!           if (lprnt) print *,' levshcm=',levshcm,' gt0sh=',gt0(ipr,:)
!    &,    ' lat=',lat

            if (Model%mstrat) then             !  As in CFSv2
              call shalcv (im, ix, levshcm, dtp, del, Statein%prsi,        &
                           Statein%prsl, Statein%prslk,kcnv, Stateout%gq0, &
                           Stateout%gt0, levshc, Statein%phil, kinver,     &
                           ctei_r, ctei_rml, lprnt, ipr)
            else
              call shalcvt3 (im, ix, levshcm, dtp, del, Statein%prsi, &
                             Statein%prsl, Statein%prslk, kcnv,       &
                             Stateout%gq0, Stateout%gt0)
            endif
!           if (lprnt) print *,' levshcm=',levshcm,' gt0sha=',gt0(ipr,:)

          endif   ! end if_imfshalcnv
        endif     ! end if_shal_cnv

#ifdef CCPP
      if (Model%me==0) write(0,*) 'CCPP DEBUG: calling GFS_SCNV_generic_post through option B'
      ! Copy local variables from driver to appropriate interstitial variables
      !Interstitial(nt)%im = im                   ! intent(in) - set in Interstitial(nt)%create()
      !Model%levs                                 ! intent(in)
      Interstitial(nt)%nn = nn                    ! intent(in)
      !Model%lssav                                ! intent(in)
      !Model%ldiag3d                              ! intent(in)
      !Model%lgocart                              ! intent(in)
      Interstitial(nt)%frain = frain              ! intent(in)
      !Stateout%gt0                               ! intent(in)
      !Stateout%gq0(:,:,1)                        ! intent(in)
      Interstitial(nt)%save_t = dtdt              ! intent(in)
      Interstitial(nt)%save_q(:,:,1) = dqdt(:,:,1)! intent(in)
      !Coupling%dqdti                             ! intent(inout)
      !Diag%dt3dt(:,:,5)                          ! intent(inout)
      !Diag%dq3dt(:,:,3)                          ! intent(inout)
      Interstitial(nt)%clw  = clw                 ! intent(inout)
      cdata_block(nb,nt)%errmsg   = errmsg        ! intent(out)
      cdata_block(nb,nt)%errflg   = errflg        ! intent(out)
      call ccpp_physics_run(cdata_block(nb,nt), scheme_name="GFS_SCNV_generic_post", ierr=ierr)
      ! Copy back intent(inout) interstitial variables to local variables in driver
      clw = Interstitial(nt)%clw
      errmsg = trim(cdata_block(nb,nt)%errmsg)
      errflg = cdata_block(nb,nt)%errflg
      if (errflg/=0) then
          write(0,*) 'Error in call to GFS_SCNV_generic_post: ' // trim(errmsg)
          stop
      end if
#else
        if (Model%lssav) then
        !          update dqdt_v to include moisture tendency due to shallow convection
          if (Model%lgocart) then
            do k=1,levs
              do i=1,im
                tem  = (Stateout%gq0(i,k,1)-dqdt(i,k,1)) * frain
                Coupling%dqdti(i,k) = Coupling%dqdti(i,k)  + tem
              enddo
            enddo
          endif
          if (Model%ldiag3d) then
            do k=1,levs
              do i=1,im
                Diag%dt3dt(i,k,5) = Diag%dt3dt(i,k,5) + (Stateout%gt0(i,k)  -dtdt(i,k))   * frain
!                Diag%dq3dt(i,k,3) = Diag%dq3dt(i,k,3) + (Stateout%gq0(i,k,1)-dqdt(i,k,1)) * frain
              enddo
            enddo
          endif
        endif   ! end if_lssav

        do k=1,levs
          do i=1,im
            if (clw(i,k,2) <= -999.0) clw(i,k,2) = 0.0
          enddo
        enddo
#endif


!       if (lprnt) then
!         write(0,*)' prsl=',prsl(ipr,:)
!         write(0,*) ' del=',del(ipr,:)
!         write(0,*) ' befshgt0=',gt0(ipr,:)
!         write(0,*) ' befshgq0=',gq0(ipr,:,1)
!       endif

      elseif (Model%shocaftcnv) then ! if do_shoc is true and shocaftcnv is true call shoc
#ifdef CCPP
        ! DH*
        write(0,*) "Error, calling CCPP-SHOC after convection has not been tested (argument list differs?)"
        stop
        ! *DH
        if (Model%me==0) write(0,*) 'CCPP DEBUG: calling shoc through option B'
        ! Copy local variables from driver to appropriate interstitial variables
        !Interstitial(nt)%ix = ix                     ! intent(in) - set in Interstitial(nt)%create()
        !Interstitial(nt)%im = im                     ! intent(in) - set in Interstitial(nt)%create()
        !Model%levs                                   ! intent(in)
        !Model%do_shoc                                ! intent(in)
        !Model%shocaftcnv                             ! intent(in)
        !Model%mg3_as_mg2                             ! intent(in)
        !Model%imp_physics                            ! intent(in)
        !Model%imp_physics_gfdl                       ! intent(in)
        !Model%imp_physics_zhao_carr                  ! intent(in)
        !Model%imp_physics_zhao_carr_pdf              ! intent(in)
        !Model%imp_physics_mg                         ! intent(in)
        !Model%fprcp                                  ! intent(in)
        !Model%tcr                                    ! intent(in)
        !Model%tcrf                                   ! intent(in)
        !con_cp, con_g, con_hvap, con_hfus, con_rv, con_rd, con_pi, con_fvirt ! intent(in) - physical constant in physcons.F90
        !Stateout%gq0(:,:,ntiw)                       ! intent(in)
        !Stateout%gq0(:,:,ntrw)                       ! intent(in)
        !Stateout%gq0(:,:,ntsw)                       ! intent(in)
        !Stateout%gq0(:,:,ntgl)                       ! intent(in)
        !Model%dtp                                    ! intent(in)
        !Model%me                                     ! intent(in)
        !Statein%prsl                                 ! intent(in)
        !Statein%phii                                 ! intent(in)
        !Statein%phil                                 ! intent(in)
        !Stateout%gu0                                 ! intent(in)
        !Stateout%gv0                                 ! intent(in)
        !Statein%vvl                                  ! intent(in)
        Interstitial(nt)%rhc = rhc                    ! intent(in)
        !Model%sup                                    ! intent(in)
        !Model%shoc_parm(1)                           ! intent(in)
        !Model%shoc_parm(2)                           ! intent(in)
        !Model%shoc_parm(3)                           ! intent(in)
        !Model%shoc_parm(4)                           ! intent(in)
        !Model%shoc_parm(5)                           ! intent(in)
        Interstitial(nt)%hflx = hflx                  ! intent(in)
        Interstitial(nt)%evap = evap                  ! intent(in)
        Interstitial(nt)%prnum = prnum                ! intent(in)
        Interstitial(nt)%skip_macro = skip_macro      ! intent(inout)
        Interstitial(nt)%clw(:,:,1) = clw(:,:,1)      ! intent(inout)
        Interstitial(nt)%clw(:,:,2) = clw(:,:,2)      ! intent(inout)
        !Stateout%gq0(:,:,ntcw)                       ! intent(inout)
        !Stateout%gq0(:,:,ntlnc)                      ! intent(inout)
        !Stateout%gq0(:,:,ntinc)                      ! intent(inout)
        !Stateout%gt0                                 ! intent(inout)
        !Stateout%gq0(:,:,1)                          ! intent(inout)
        !Tbd%phy_f3d(:,:,ntot3d-2)                    ! intent(inout)
        Interstitial(nt)%clw(:,:,ntk) = clw(:,:,ntk)  ! intent(inout)
        !Tbd%phy_f3d(:,:,ntot3d-1)                    ! intent(inout)
        !Tbd%phy_f3d(:,:,ntot3d)                      ! intent(inout)
        !cdata_block(nb,nt)%errmsg = errmsg           ! intent(out)
        !cdata_block(nb,nt)%errflg = errflg           ! intent(out)
        call ccpp_physics_run(cdata_block(nb,nt), scheme_name="shoc", ierr=ierr)
        ! Copy back intent(inout) interstitial variables to local variables in driver
        skip_macro = Interstitial(nt)%skip_macro
        clw(:,:,1) = Interstitial(nt)%clw(:,:,1)
        clw(:,:,2) = Interstitial(nt)%clw(:,:,2)
        clw(:,:,ntk) = Interstitial(nt)%clw(:,:,ntk)
        errmsg = trim(cdata_block(nb,nt)%errmsg)
        errflg = cdata_block(nb,nt)%errflg
        if (errflg/=0) then
            write(0,*) 'Error in call to shoc: ' // trim(errmsg)
            stop
        end if
#else
        if (imp_physics == Model%imp_physics_mg) then
          skip_macro = Model%do_shoc
          do k=1,levs
            do i=1,im
              ncpl(i,k)  = Stateout%gq0(i,k,ntlnc)
              ncpi(i,k)  = Stateout%gq0(i,k,ntinc)
            enddo
          enddo

!       else
!         if (clw(1,1,2) < -999.0) then ! if clw is not partitioned to ice and water
!           do k=1,levs
!             do i=1,im
!               tem = gq0(i,k,ntcw)                                     &
!    &              * max(0.0, MIN(1.0, (TCR-gt0(i,k))*TCRF))
!               clw(i,k,1) = tem                              ! ice
!               clw(i,k,2) = gq0(i,k,ntcw) - tem              ! water
!             enddo
!           enddo
!         endif     ! Anning ncld ==2
          if (abs(Model%fprcp) == 1 .or. mg3_as_mg2) then
            do k=1,levs
              do i=1,im
                qrn(i,k)  = Stateout%gq0(i,k,ntrw)
                qsnw(i,k) = Stateout%gq0(i,k,ntsw)
              enddo
            enddo
          elseif (Model%fprcp > 1) then
            do k=1,levs
              do i=1,im
                qrn(i,k)  = Stateout%gq0(i,k,ntrw)
                qsnw(i,k) = Stateout%gq0(i,k,ntsw) + Stateout%gq0(i,k,ntgl)
              enddo
            enddo
          endif
        endif

!       dtshoc = 60.0
!       dtshoc = min(dtp, 300.0)
!       nshocm = max(1, nint(dtp/dtshoc))
!       dtshoc = dtp / nshocm
!       do nshoc=1,nshocm
!       call shoc(im, 1, levs, levs+1, dtp, me, lat,        &
!!       call shoc(im, 1, levs, levs+1, dtshoc, me, lat, &
!    &                       prsl(1:im,:), phii (1:im,:),  phil(1:im,:),&
!    &          gu0(1:im,:),gv0(1:im,:), vvl(1:im,:), gt0(1:im,:),     &
!    &                                                   gq0(1:im,:,1), &
!    &          clw(1:im,:,1), clw(1:im,:,2), qsnw, qrn, sgs_cld(1:im,:)&
!    &,         gq0(1:im,:,ntke),                                       &
!    &          phy_f3d(1:im,:,ntot3d-1), phy_f3d(1:im,:,ntot3d),       &
!    &          lprnt, ipr,                                             &
!    &          con_cp, con_g, con_hvap, con_hfus, con_hvap+con_hfus,   &
!    &          con_rv, con_rd, con_pi, con_fvirt)

!GFDL  replace lat with "1:
!       call shoc(ix, im, 1, levs, levs+1, dtshoc, me, lat,             &
        call shoc (ix, im, 1, levs, levs+1, dtp, me, 1, Statein%prsl(1,1),    &
                   Statein%phii(1,1), Statein%phil(1,1), Stateout%gu0(1,1),   &
                   Stateout%gv0(1,1), Statein%vvl(1,1), Stateout%gt0(1,1),    &
                   Stateout%gq0(1,1,1), clw(1,1,1), clw(1,1,2),               &
                   qsnw, qrn, qgl,                                            &
                   rhc, Model%sup, Model%shoc_parm(1), Model%shoc_parm(2),    &
                   Model%shoc_parm(3), Model%shoc_parm(4),                    &
                   Model%shoc_parm(5), Tbd%phy_f3d(1,1,ntot3d-2),             &
                   Stateout%gq0(1,1,ntke), hflx, evap, prnum,                 &
                   Tbd%phy_f3d(1,1,ntot3d-1), Tbd%phy_f3d(1,1,ntot3d),        &
                   lprnt, ipr, ncpl, ncpi)
!       enddo

        if (ntlnc > 0 .and. ntinc > 0 .and. ncld >= 2) then
          do k=1,levs
            do i=1,im
              Stateout%gq0(i,k,ntlnc) = ncpl(i,k)
              Stateout%gq0(i,k,ntinc) = ncpi(i,k)
            enddo
          enddo
        endif

!
!      do k=1,levs
!      write(1000+me,*)' maxtkh=',maxval(phy_f3d(1:im,k,ntot3d-1)), &
!     ' k=',k,' kdt=',kdt,' lat=',lat
!      enddo

!     write(0,*)' aft shoc gt0=',gt0(1,:),' lat=',lat
!     write(0,*)' aft shoc gq0=',gq0(1,:,1),' lat=',lat
!     write(0,*)' aft shoc gu0=',gu0(1,:),' lat=',lat
!
#endif
      endif   ! if( .not. do_shoc)
!
!       if (lprnt) then
!         write(0,*)' prsl=',prsl(ipr,:)
!         write(0,*) ' del=',del(ipr,:)
!         write(0,*) ' aftshgt0=',gt0(ipr,:)
!         write(0,*) ' aftshgq0=',gq0(ipr,:,1)
!       endif
!
#ifdef CCPP
! OPTION B - works with all compilers
      if (Model%me==0) write(0,*) 'CCPP DEBUG: calling GFS_suite_interstitial_4 through option B'
      ! Copy local variables from driver to appropriate interstitial variables
      !Interstitial(nt)%im = im                   ! intent(in) - set in Interstitial(nt)%create()
      !Model%levs                                 ! intent(in)
      !Model%ltaerosol                            ! intent(in)
      Interstitial(nt)%tracers_total = tottracer  ! intent(in)
      !Model%ntrac                                ! intent(in)
      !Model%ntcw                                 ! intent(in)
      !Model%ntiw                                 ! intent(in)
      !Model%ntclamt                              ! intent(in)
      !Model%ntrw                                 ! intent(in)
      !Model%ntsw                                 ! intent(in)
      !Model%ntrnc                                ! intent(in)
      !Model%ntsnc                                ! intent(in)
      !Model%ntgl                                 ! intent(in)
      !Model%ntgnc                                ! intent(in)
      !Model%ntlnc                                ! intent(in)
      !Model%ntinc                                ! intent(in)
      Interstitial(nt)%nn = nn                    ! intent(in)
      !Model%imp_physics                          ! intent(in)
      !Model%imp_physics_gfdl                     ! intent(in)
      !Model%imp_physics_thompson                 ! intent(in)
      !Model%imp_physics_zhao_carr                ! intent(in)
      !Model%imp_physics_zhao_carr_pdf            ! intent(in)
      if (Model%imp_physics == Model%imp_physics_thompson) then
        if (Model%ltaerosol) then
          Interstitial(nt)%save_q(:,:,Model%ntcw) = liq0 ! intent(in)
          Interstitial(nt)%save_q(:,:,Model%ntiw) = ice00            ! intent(in)
        else
          Interstitial(nt)%save_q(:,:,Model%ntiw) = ice00            ! intent(in)
        end if
      end if
      !con_pi                                     ! intent(in) - from physcons
      !Stateout%gq0                               ! intent(inout)
      Interstitial(nt)%clw = clw                  ! intent(inout)
      cdata_block(nb,nt)%errmsg   = errmsg        ! intent(out)
      cdata_block(nb,nt)%errflg   = errflg        ! intent(out)
      call ccpp_physics_run(cdata_block(nb,nt), scheme_name="GFS_suite_interstitial_4", ierr=ierr)
      ! Copy back intent(inout) interstitial variables to local variables in driver
      clw = Interstitial(nt)%clw
      errmsg = trim(cdata_block(nb,nt)%errmsg)
      errflg = cdata_block(nb,nt)%errflg
      if (errflg/=0) then
          write(0,*) 'Error in call to GFS_suite_interstitial_4: ' // trim(errmsg)
          stop
      end if
#else
      if (Model%me==0) write(0,*) 'CCPP DEBUG: calling non-CCPP compliant version of GFS_suite_interstitial_4'
!------------------------------------------------------------------------------
!  --- update the tracers due to deep & shallow cumulus convective transport
!           (except for suspended water and ice)
!
      if (tottracer > 0) then
        tracers = 2
        do n=2,ntrac
!         if ( n /= ntcw .and. n /= ntiw .and. n /= ntclamt) then
          if ( n /= ntcw  .and. n /= ntiw  .and. n /= ntclamt .and. &
               n /= ntrw  .and. n /= ntsw  .and. n /= ntrnc   .and. &
               n /= ntsnc .and. n /= ntgl  .and. n /= ntgnc ) then
              tracers = tracers + 1
            do k=1,levs
              do i=1,im
                Stateout%gq0(i,k,n) = clw(i,k,tracers)
              enddo
            enddo
          endif
        enddo
      endif
!-------------------------------------------------------------------------------
!
      if (ntcw > 0) then

!  for microphysics
        if (imp_physics == 99 .or. imp_physics == 98    &
                              .or. imp_physics == Model%imp_physics_gfdl) then
           Stateout%gq0(1:im,:,ntcw) = clw(1:im,:,1) + clw(1:im,:,2)
        elseif (ntiw > 0) then
          do k=1,levs
            do i=1,im
              Stateout%gq0(i,k,ntiw) = clw(i,k,1)                     ! ice
              Stateout%gq0(i,k,ntcw) = clw(i,k,2)                     ! water
            enddo
          enddo
          if (imp_physics == Model%imp_physics_thompson) then
            if (Model%ltaerosol) then
              do k=1,levs
                do i=1,im
                  Stateout%gq0(i,k,ntlnc) = Stateout%gq0(i,k,ntlnc)  &
                           +  max(0.0, (clw(i,k,2)-liq0(i,k))) / liqm
                  Stateout%gq0(i,k,ntinc) = Stateout%gq0(i,k,ntinc)  &
                           +  max(0.0, (clw(i,k,1)-ice00(i,k))) / icem
                enddo
              enddo
            else
              do k=1,levs
                do i=1,im
                  Stateout%gq0(i,k,ntinc) = Stateout%gq0(i,k,ntinc)  &
                           +  max(0.0, (clw(i,k,1)-ice00(i,k))) / icem
                enddo
              enddo
            endif
          endif
        else
          do k=1,levs
            do i=1,im
              Stateout%gq0(i,k,ntcw) = clw(i,k,1) + clw(i,k,2)
            enddo
          enddo
        endif   ! end if_ntiw
      else
        do k=1,levs
          do i=1,im
            clw(i,k,1) = clw(i,k,1) + clw(i,k,2)
          enddo
        enddo
      endif   ! end if_ntcw
#endif

!  Legacy routine which determines convectve clouds - should be removed at some point

#ifdef CCPP
      if (Model%me==0) write(0,*) 'CCPP DEBUG: calling cnvc90 through option B'
      ! Copy local variables from driver to appropriate interstitial variables
      !Interstitial(nt)%im = im            ! intent(in) - set in Interstitial(nt)%create()
      !Interstitial(nt)%ix = ix            ! intent(in) - set in Interstitial(nt)%create()
      !Diag%rainc                          ! intent(in)
      Interstitial(nt)%kbot = kbot         ! intent(in)
      Interstitial(nt)%ktop = ktop         ! intent(in)
      !Model%levs = levs                   ! intent(in)
      !Statein%prsi                        ! intent(in)
      !Tbd%acv                             ! intent(inout)
      !Tbd%acvb                            ! intent(inout)
      !Tbd%acvt                            ! intent(inout)
      !Cldprop%cv                          ! intent(inout)
      !Cldprop%cvb                         ! intent(inout)
      !Cldprop%cvt                         ! intent(inout)
      !cdata_block(nb,nt)%errmsg = errmsg  ! intent(out)
      !cdata_block(nb,nt)%errflg = errflg  ! intent(out)
      call ccpp_physics_run(cdata_block(nb,nt), scheme_name="cnvc90", ierr=ierr)
      ! Copy intent(inout) and intent(out) interstitial variables to local variables in driver
      errmsg = trim(cdata_block(nb,nt)%errmsg)
      errflg = cdata_block(nb,nt)%errflg
      if (errflg/=0) then
          write(0,*) 'Error in call to cnvc90: ' // trim(errmsg)
          stop
      end if
#else
      if (Model%me==0) write(0,*) 'CCPP DEBUG: calling non-CCPP compliant version of cnvc90'
      call cnvc90 (Model%clstp, im, ix, Diag%rainc, kbot, ktop, levs, Statein%prsi,  &
                   Tbd%acv, Tbd%acvb, Tbd%acvt, Cldprop%cv, Cldprop%cvb, Cldprop%cvt)
#endif

! DH* - this block is not yet in CCPP
      if (Model%moist_adj) then       ! moist convective adjustment
!                                     ---------------------------
!
!       To call moist convective adjustment
!
!       if (lprnt) then
!         print *,' prsl=',prsl(ipr,:)
!         print *,' del=',del(ipr,:)
!         print *,' gt0b=',gt0(ipr,:)
!         print *,' gq0b=',gq0(ipr,:,1)
!       endif

        call mstcnv (im, ix, levs, dtp, Stateout%gt0, Stateout%gq0, &
                     Statein%prsl,del, Statein%prslk, rain1,        &
                     Stateout%gq0(1,1,ntcw), rhc, lprnt, ipr)

!       if (lprnt) then
!         print *,' rain1=',rain1(ipr),' rainc=',rainc(ipr)
!         print *,' gt0a=',gt0(ipr,:)
!         print *,' gq0a=',gq0(ipr,:,1)
!       endif
        do i=1,im
          Diag%rainc(i) = Diag%rainc(i) + frain * rain1(i)
        enddo
        if(Model%lssav) then
          do i=1,im
            Diag%cnvprcp(i)  = Diag%cnvprcp(i)  + rain1(i) * frain
            Diag%cnvprcpb(i) = Diag%cnvprcpb(i) + rain1(i) * frain
          enddo

! update dqdt_v to include moisture tendency due to surface processes
! dqdt_v : instaneous moisture tendency (kg/kg/sec)
!          if (lgocart) then
!            do k=1,levs
!              do i=1,im
!                tem = (gq0(i,k,1)-dqdt(i,k,1)) * frain
!                dqdt_v(i,k) = dqdt_v(i,k) + tem
!                dqdt_v(i,k) = dqdt_v(i,k) / dtf
!              enddo
!            enddo
!          endif
          if (Model%ldiag3d) then
            do k=1,levs
              do i=1,im
! DH* dt3dt has last dim 1:7 only (GFS_typedefs.F90)
!                 Diag%dt3dt(i,k,8) = Diag%dt3dt(i,k,8) + (Stateout%gt0(i,k)  -dtdt(i,k)  ) * frain
! *DH
!                Diag%dq3dt(i,k,2) = Diag%dq3dt(i,k,2) + (Stateout%gq0(i,k,1)-dqdt(i,k,1)) * frain
              enddo
            enddo
          endif
        endif
      endif               !       moist convective adjustment over
! *DH
!
#ifdef CCPP
! OPTION B - works with all compilers
      if (Model%me==0) write(0,*) 'CCPP DEBUG: calling GFS_MP_generic_pre through option B'
      ! Copy local variables from driver to appropriate interstitial variables
      !Interstitial(nt)%im = im                   ! intent(in) - set in Interstitial(nt)%create()
      !Model%levs                                 ! intent(in)
      !Model%ldiag3d                              ! intent(in)
      !Model%do_aw                                ! intent(in)
      !Model%ntcw                                 ! intent(in)
      Interstitial(nt)%nncl = nncl                ! intent(in)
      !Model%ntrac                                ! intent(in)
      !Stateout%gt0                               ! intent(in)
      !Stateout%gq0                               ! intent(in)
      Interstitial(nt)%save_t = dtdt              ! intent(inout)
      Interstitial(nt)%save_q = dqdt              ! intent(inout)
      !cdata_block(nb,nt)%errmsg = errmsg         ! intent(out)
      !cdata_block(nb,nt)%errflg = errflg         ! intent(out)
      call ccpp_physics_run(cdata_block(nb,nt), scheme_name="GFS_MP_generic_pre", ierr=ierr)
      ! Copy back intent(inout) interstitial variables to local variables in driver
      dtdt = Interstitial(nt)%save_t
      dqdt = Interstitial(nt)%save_q
      errmsg = trim(cdata_block(nb,nt)%errmsg)
      errflg = cdata_block(nb,nt)%errflg
      if (errflg/=0) then
          write(0,*) 'Error in call to GFS_MP_generic_pre: ' // trim(errmsg)
          stop
      end if
#else
      if (Model%ldiag3d .or. Model%do_aw) then
        do k=1,levs
          do i=1,im
            dtdt(i,k)   = Stateout%gt0(i,k)
            dqdt(i,k,1) = Stateout%gq0(i,k,1)
          enddo
        enddo
        do n=ntcw,ntcw+nncl-1
          dqdt(1:im,:,n) = Stateout%gq0(1:im,:,n)
        enddo
      endif
#endif

! dqdt_v : instaneous moisture tendency (kg/kg/sec)
#ifndef CCPP
      !GF* the following code is executed in GFS_suite_interstitial_4 (relevant for shallow and deep convection)
      if (Model%lgocart) then
        do k=1,levs
          do i=1,im
            Coupling%dqdti(i,k) = Coupling%dqdti(i,k) * (1.0 / dtf)
          enddo
        enddo
      endif
      !*GF
#endif
!
!     grid-scale condensation/precipitations and microphysics parameterization
!     ------------------------------------------------------------------------

      if (ncld == 0) then                   ! no cloud microphysics

        call lrgscl (ix, im, levs, dtp, Stateout%gt0, Stateout%gq0, &
                     Statein%prsl, del, Statein%prslk, rain1, clw)

      else                                  ! all microphysics
        if (imp_physics == 99) then         ! call zhao/carr/sundqvist microphysics
                                            ! ------------

!         if (lprnt) then
!           write(0,*)' prsl=',prsl(ipr,:)
!           write(0,*) ' del=',del(ipr,:)
!           write(0,*) ' beflsgt0=',gt0(ipr,:),' kdt=',kdt
!           write(0,*) ' beflsgq0=',gq0(ipr,:,1),' kdt=',kdt
!           write(0,*) ' beflsgw0=',gq0(ipr,:,3),' kdt=',kdt
!         endif
                                              ! ------------------
          if (Model%do_shoc) then
            call precpd_shoc (im, ix, levs, dtp, del, Statein%prsl,            &
                              Stateout%gq0(1,1,1), Stateout%gq0(1,1,ntcw),     &
                              Stateout%gt0, rain1, Sfcprop%sr, rainp, rhc,     &
                              psautco_l, prautco_l, Model%evpco, Model%wminco, &
                              Tbd%phy_f3d(1,1,ntot3d-2), lprnt, ipr)
          else
#ifdef CCPP
            if (Model%me==0) write(0,*) 'CCPP DEBUG: calling zhaocarr_gscond through option B'
            ! Copy local variables from driver to appropriate interstitial variables
            !Interstitial(nt)%im     = im                           ! intent(in) - set in Interstitial(nt)%create()
            !Interstitial(nt)%ix     = ix                           ! intent(in) - set in Interstitial(nt)%create()
            Interstitial(nt)%clw(:,:,1) = clw(:,:,1)                ! intent(in)
            Interstitial(nt)%clw(:,:,2) = clw(:,:,2)                ! intent(in)
            Interstitial(nt)%rhc    = rhc                           ! intent(in)
            Interstitial(nt)%ipr    = ipr                           ! intent(in)
            !cdata_block(nb,nt)%errmsg = errmsg                     ! intent(out)
            !cdata_block(nb,nt)%errflg = errflg                     ! intent(out)
            call ccpp_physics_run(cdata_block(nb,nt), scheme_name="zhaocarr_gscond", ierr=ierr)
            ! Copy back intent(inout) interstitial variables to local variables in driver
            errmsg = trim(cdata_block(nb,nt)%errmsg)
            errflg = cdata_block(nb,nt)%errflg
            if (errflg/=0) then
              write(0,*) 'Error in call to zhaocarr_gscond: ' // trim(errmsg)
              stop
            end if
            !
            if (Model%me==0) write(0,*) 'CCPP DEBUG: calling zhaocarr_precpd through option B'
            ! Copy local variables from driver to appropriate interstitial variables
            !Interstitial(nt)%im     = im                           ! intent(in   ) - set in Interstitial(nt)%create()
            !Interstitial(nt)%ix     = ix                           ! intent(in   ) - set in Interstitial(nt)%create()
            !Model%levs                                             ! intent(in   )
            !Model%dtp                                              ! intent(in   )
            Interstitial(nt)%del = del                              ! intent(in   )
            !Statein%prsl                                           ! intent(in   )
            !Stateout%gq0(:,:,1)                                    ! intent(inout)
            !Stateout%gq0(:,:,Model%ntcw)                           ! intent(inout)
            !Stateout%gt0                                           ! intent(inout)
            Interstitial(nt)%prcpmp = rain1                         ! intent(out  )
            !Sfcprop%sr                                             ! intent(out  )
            Interstitial(nt)%rainp = rainp                          ! intent(out  )
            Interstitial(nt)%rhc = rhc                              ! intent(in   )
            !Model%psautco                                          ! intent(in   )
            !Model%prautco                                          ! intent(in   )
            !Model%evpco                                            ! intent(in   )
            !Model%wminco                                           ! intent(in   )
            !Interstitial(nt)%work1                                 ! intent(in   ) - set in GFS_suite_interstitial_1
            !Model%lprnt                                            ! intent(in   )
            !Model%ipr                                              ! intent(in   )
            !cdata_block(nb,nt)%errmsg = errmsg                     ! intent(out  )
            !cdata_block(nb,nt)%errflg = errflg                     ! intent(out  )
            call ccpp_physics_run(cdata_block(nb,nt), scheme_name="zhaocarr_precpd", ierr=ierr)
            ! Copy back intent(inout) interstitial variables to local variables in driver
            rain1  = Interstitial(nt)%prcpmp
            rainp  = Interstitial(nt)%rainp
            errmsg = trim(cdata_block(nb,nt)%errmsg)
            errflg = cdata_block(nb,nt)%errflg
            if (errflg/=0) then
              write(0,*) 'Error in call to zhaocarr_precpd: ' // trim(errmsg)
              stop
            end if
#else
            if (Model%me==0) write(0,*) 'CCPP DEBUG: calling non-CCPP compliant version of gscond'
            call gscond (im, ix, levs, dtp, dtf, Statein%prsl, Statein%pgr,    &
                         Stateout%gq0(1,1,1), Stateout%gq0(1,1,ntcw),          &
                         Stateout%gt0, Tbd%phy_f3d(1,1,1), Tbd%phy_f3d(1,1,2), &
                         Tbd%phy_f2d(1,1), Tbd%phy_f3d(1,1,3),                 &
                         Tbd%phy_f3d(1,1,4), Tbd%phy_f2d(1,2), rhc,lprnt, ipr)
            if (Model%me==0) write(0,*) 'CCPP DEBUG: calling non-CCPP compliant version of precpd'
            call precpd (im, ix, levs, dtp, del, Statein%prsl,                 &
                        Stateout%gq0(1,1,1), Stateout%gq0(1,1,ntcw),           &
                        Stateout%gt0, rain1, Sfcprop%sr, rainp, rhc, psautco_l,&
                        prautco_l, Model%evpco, Model%wminco, lprnt, ipr)
#endif
          endif
!         if (lprnt) then
!           write(0,*)' prsl=',prsl(ipr,:)
!           write(0,*) ' del=',del(ipr,:)
!           write(0,*) ' aftlsgt0=',gt0(ipr,:),' kdt=',kdt
!           write(0,*) ' aftlsgq0=',gq0(ipr,:,1),' kdt=',kdt
!           write(0,*) ' aftlsgw0=',gq0(ipr,:,3),' kdt=',kdt
!           write(0,*)' aft precpd rain1=',rain1(1:3),' lat=',lat
!           endif

        elseif (imp_physics == 98) then       ! with pdf clouds
            call gscondp (im, ix, levs, dtp, dtf, Statein%prsl,        &
                          Statein%pgr, Stateout%gq0(1,1,1),            &
                          Stateout%gq0(1,1,ntcw), Stateout%gt0,        &
                          Tbd%phy_f3d(1,1,1), Tbd%phy_f3d(1,1,2),      &
                          Tbd%phy_f2d(1,1),   Tbd%phy_f3d(1,1,3),      &
                          Tbd%phy_f3d(1,1,4), Tbd%phy_f2d(1,2), rhc,   &
                          Tbd%phy_f3d(1,1,Model%num_p3d+1), Model%sup, &
                          lprnt, ipr, kdt)

            call precpdp (im, ix, levs,  dtp, del, Statein%prsl,       &
                          Statein%pgr, Stateout%gq0(1,1,1),            &
                          Stateout%gq0(1,1,ntcw), Stateout%gt0,        &
                          rain1, Sfcprop%sr, rainp, rhc,               &
                          Tbd%phy_f3d(1,1,Model%num_p3d+1), psautco_l, &
                          prautco_l, Model%evpco, Model%wminco, lprnt, ipr)
!          endif   ! end of grid-scale precip/microphysics options
!        endif     ! end if_num_p3d

!     if (lprnt) write(0,*) ' rain1=',rain1(ipr),' rainc=',rainc(ipr),' lat=',lat

        elseif (imp_physics == Model%imp_physics_thompson) then      !  Thompson MP
                                            ! ------------
          ims = 1 ; ime = ix ; kms = 1 ; kme = levs ; its = 1 ; ite = ix ; kts = 1 ; kte = levs

#ifndef CCPP
          if (Model%ltaerosol) then
            print*,'aerosol version of the Thompson scheme is not included'

!           call mp_gt_driver(ims,ime,kms,kme,its,ite,kts,kte,                             &
!              Stateout%gq0(1:im,1:levs,1),                                                &
!              Stateout%gq0(1:im,1:levs,Model%ntcw), Stateout%gq0(1:im,1:levs,Model%ntrw), &
!              Stateout%gq0(1:im,1:levs,Model%ntiw), Stateout%gq0(1:im,1:levs,Model%ntsw), &
!              Stateout%gq0(1:im,1:levs,Model%ntgl), Stateout%gq0(1:im,1:levs,Model%ntinc),&
!              Stateout%gq0(1:im,1:im,Model%ntrnc),                                        &
!              Stateout%gt0, Statein%prsl, Statein%vvl, del, dtp, kdt,                     &
!              rain1,                                                                      &
!              Sfcprop%sr,                                                                 &
!!             Diag%refl_10cm, Model%lradar,                                               &
!!             Tbd%phy_f3d(:,:,1),Tbd%phy_f3d(:,:,2),Tbd%phy_f3d(:,:,3),                   & !has_reqc, has_reqi, has_reqs,
!!             ims,ime,kms,kme,its,ite,kts,kte)
!              Tbd%phy_f3d(:,:,1),Tbd%phy_f3d(:,:,2),Tbd%phy_f3d(:,:,3),me,                &
!              nc=Stateout%gq0(1:im,1:levs,Model%ntlnc),                                   &
!              nwfa=Stateout%gq0(1:im,1:levs,Model%ntwa),                                  &
!              nifa=Stateout%gq0(1:im,1:levs,Model%ntia),                                  &
!!             nwfa2d=Sfcprop%nwfa2d(1:im))
!              nwfa2d=Coupling%nwfa2d(1:im))
          else
            call mp_gt_driver(ims,ime,kms,kme,its,ite,kts,kte,                             &
               Stateout%gq0(1:im,1:levs,1),                                                &
               Stateout%gq0(1:im,1:levs,Model%ntcw), Stateout%gq0(1:im,1:levs,Model%ntrw), &
               Stateout%gq0(1:im,1:levs,Model%ntiw), Stateout%gq0(1:im,1:levs,Model%ntsw), &
               Stateout%gq0(1:im,1:levs,Model%ntgl), Stateout%gq0(1:im,1:levs,Model%ntinc),&
               Stateout%gq0(1:im,1:levs,Model%ntrnc),                                      &
!2014v         Stateout%gt0, Statein%prsl, Statein%vvl, del, dtp, kdt,                     &
               Stateout%gt0, Statein%prsl, del, dtp, kdt,                                  &
               rain1,                                                                      &
               Sfcprop%sr,                                                                 &
               islmsk,                                                                     &
               Diag%refl_10cm, Model%lradar,                                               &
               Tbd%phy_f3d(:,:,1),Tbd%phy_f3d(:,:,2),Tbd%phy_f3d(:,:,3),me,Statein%phii)
          endif
#else
          if (Model%me==0) write(0,*) 'CCPP DEBUG: calling mp_thompson_pre_run through option B'
          !Interstitial(nt)%im = im                             ! intent(in) - set in Interstitial(nt)%create()
          !Model%levs                                           ! intent(in)
          !Model%kdt                                            ! intent(in)
          !con_g                                                ! intent(in)
          !con_rd                                               ! intent(in)
          !Model%ltaerosol                                      ! intent(in)
          !Stateout%gq0(:,:,IPD_Control%ntwa)                   ! intent(inout)
          !Stateout%gq0(:,:,IPD_Control%ntia)                   ! intent(inout)
          !Coupling%nwfa2d                                      ! intent(inout)
          !Coupling%nifa2d                                      ! intent(inout)
          !Stateout%gt0                                         ! intent(in)
          !Interstitial(nt)%save_t                              ! intent(  out)
          !Statein%prsl                                         ! intent(in)
          !Statein%phil                                         ! intent(in)
          !Model%me                                             ! intent(in)
          !Model%master                                         ! intent(in)
          !Model%communicator                                   ! intent(in)
          !cdata_block(nb,nt)%errmsg = errmsg                   ! intent(  out)
          !cdata_block(nb,nt)%errflg = errflg                   ! intent(  out)
          !
          call ccpp_physics_run(cdata_block(nb,nt), scheme_name="mp_thompson_pre", ierr=ierr)
          ! Copy back intent(inout) and intent(out) interstitial variables to local variables in driver
          errmsg = trim(cdata_block(nb,nt)%errmsg)
          errflg = cdata_block(nb,nt)%errflg
          !
          if (errflg/=0) then
              write(0,*) 'Error in call to mp_thompson_mp_mp_thompson_run: ' // trim(errmsg)
              stop
          end if
          if (Model%me==0) write(0,*) 'CCPP DEBUG: calling mp_thompson_run through option B'
          ! Copy local variables from driver to appropriate interstitial variables
          !Interstitial(nt)%im = im                             ! intent(in) - set in Interstitial(nt)%create()
          !Model%levs                                           ! intent(in)
          !con_g                                                ! intent(in)
          !con_rd                                               ! intent(in)
          !Stateout%gq0(:,:,1)                                  ! intent(inout)
          !Stateout%gq0(:,:,Model%ntcw)                         ! intent(inout)
          !Stateout%gq0(:,:,Model%ntrw)                         ! intent(inout)
          !Stateout%gq0(:,:,Model%ntiw)                         ! intent(inout)
          !Stateout%gq0(:,:,Model%ntsw)                         ! intent(inout)
          !Stateout%gq0(:,:,Model%ntgl)                         ! intent(inout)
          !Stateout%gq0(:,:,Model%ntinc)                        ! intent(inout)
          !Stateout%gq0(:,:,Model%ntrnc)                        ! intent(inout)
          !Model%ltaerosol                                      ! intent(in)
          !Stateout%gq0(:,:,Model%ntlnc)                        ! intent(inout)
          !Stateout%gq0(:,:,Model%ntwa)                         ! intent(inout)
          !Stateout%gq0(:,:,Model%ntia)                         ! intent(inout)
          !Coupling%nwfa2d                                      ! intent(in)
          !Coupling%nifa2d                                      ! intent(in)
          !Stateout%gt0                                         ! intent(inout)
          !Statein%prsl                                         ! intent(in)
          !Statein%phii                                         ! intent(in)
          !Statein%vvl                                          ! intent(in)
          !Model%dtp                                            ! intent(in)
          !Interstitial(nt)%prcpmp                              ! intent(  out)
          !Interstitial(nt)%rainmp                              ! intent(  out)
          !Interstitial(nt)%graupelmp                           ! intent(  out)
          !Interstitial(nt)%icemp                               ! intent(  out)
          !Interstitial(nt)%snowmp                              ! intent(  out)
          !Sfcprop%sr                                           ! intent(  out)
          !Diag%refl_10cm                                       ! intent(  out)
          !Model%lradar                                         ! intent(in   )
          !Tbd%phy_f3d(:,:,IPD_Control%nleffr)                  ! intent(inout)
          !Tbd%phy_f3d(:,:,IPD_Control%nieffr)                  ! intent(inout)
          !Tbd%phy_f3d(:,:,IPD_Control%nseffr)                  ! intent(inout)
          !Model%me                                             ! intent(in)
          !Model%master                                         ! intent(in)
          !Model%communicator                                   ! intent(in)
          !cdata_block(nb,nt)%errmsg = errmsg                   ! intent(  out)
          !cdata_block(nb,nt)%errflg = errflg                   ! intent(  out)
          !
          call ccpp_physics_run(cdata_block(nb,nt), scheme_name="mp_thompson", ierr=ierr)
          ! Copy back intent(inout) and intent(out) interstitial variables to local variables in driver
          rain1 = Interstitial(nt)%prcpmp
          errmsg = trim(cdata_block(nb,nt)%errmsg)
          errflg = cdata_block(nb,nt)%errflg
          !
          if (errflg/=0) then
              write(0,*) 'Error in call to mp_thompson: ' // trim(errmsg)
              stop
          end if
          !
          if (Model%me==0) write(0,*) 'CCPP DEBUG: calling mp_thompson_post_run through option B'
          ! Copy local variables from driver to appropriate interstitial variables
          !Interstitial(nt)%im = im                             ! intent(in) - set in Interstitial(nt)%create()
          !Model%levs                                           ! intent(in)
          !Interstitial(nt)%save_t                              ! intent(in) - coming straight from mp_thompson_pre_run
          !Stateout%gt0                                         ! intent(inout)
          !Statein%prslk                                        ! intent(in)
          !Model%dtp                                            ! intent(in)
          !cdata_block(nb,nt)%errmsg = errmsg                   ! intent(  out)
          !cdata_block(nb,nt)%errflg = errflg                   ! intent(  out)
          !
          call ccpp_physics_run(cdata_block(nb,nt), scheme_name="mp_thompson_post", ierr=ierr)
          ! Copy back intent(inout) and intent(out) interstitial variables to local variables in driver
          errmsg = trim(cdata_block(nb,nt)%errmsg)
          errflg = cdata_block(nb,nt)%errflg
          !
          if (errflg/=0) then
              write(0,*) 'Error in call to mp_thompson_mp_mp_thompson_post_run: ' // trim(errmsg)
              stop
          end if
#endif
        elseif (imp_physics == Model%imp_physics_wsm6) then      ! WSM6
                                            ! -----
          ims = 1 ; ime = ix ; kms = 1 ; kme = levs ; its = 1 ; ite = ix ; kts = 1 ; kte = levs

           call wsm6(Stateout%gt0, Statein%phii(1:im,1:levs+1),                                 &
                                Stateout%gq0(1:im,1:levs,1),                                    &
                                Stateout%gq0(1:im,1:levs,Model%ntcw),                           &
                                Stateout%gq0(1:im,1:levs,Model%ntrw),                           &
                                Stateout%gq0(1:im,1:levs,Model%ntiw),                           &
                                Stateout%gq0(1:im,1:levs,Model%ntsw),                           &
                                Stateout%gq0(1:im,1:levs,Model%ntgl),                           &
                                Statein%prsl, del, dtp, rain1,                                  &
                                Sfcprop%sr,                                                     &
                                islmsk,                                                         &
                                Tbd%phy_f3d(:,:,1),Tbd%phy_f3d(:,:,2),Tbd%phy_f3d(:,:,3),       &
                                ims,ime, kms,kme,                                               &
                                its,ite, kts,kte)
!
      elseif (imp_physics == Model%imp_physics_mg) then       ! MGB double-moment microphysics
                                                              ! ------------------------------
#ifdef CCPP
        if (Model%me==0) write(0,*) 'CCPP DEBUG: calling m_micro_pre_run through option B'
        !Interstitial(nt)%im = im                             ! intent(in) - set in Interstitial(nt)%create()
        !Model%levs                                           ! intent(in)
        !Model%do_shoc                                        ! intent(in)
        !Model%fprcp                                          ! intent(in)
        !Model%mg3_as_mg2                                     ! intent(in)
        !Model%imfdeepcnv                                     ! intent(in)
        !Model%imfshalcnv                                     ! intent(in)
        !Stateout%gq0(:,:,ntiw)                               ! intent(in)
        !Stateout%gq0(:,:,ntcw)                               ! intent(in)
        !Stateout%gq0(:,:,ntrw)                               ! intent(in)
        !Stateout%gq0(:,:,ntsw)                               ! intent(in)
        !Stateout%gq0(:,:,ntgl)                               ! intent(in)
        !Stateout%gq0(:,:,ntrnc)                              ! intent(in)
        !Stateout%gq0(:,:,ntsnc)                              ! intent(in)
        !Stateout%gq0(:,:,ntgnc)                              ! intent(in)
        !Tbd%phy_f3d(:,:,ntot3d-2)                            ! intent(in)
        if (allocated(cnvc) .and. allocated(cnvw)) then
          Interstitial(nt)%cnvc = cnvc                        ! intent(in)
          Interstitial(nt)%cnvw = cnvw                        ! intent(in)
        end if
        !Model%tcr                                            ! intent(in)
        !Model%tcrf                                           ! intent(in)
        !Stateout%gt0                                         ! intent(in)
        Interstitial(nt)%qrn = qrn                            ! intent(inout)
        Interstitial(nt)%qsnw = qsnw                          ! intent(inout)
        Interstitial(nt)%qgl = qgl                            ! intent(inout)
        Interstitial(nt)%ncpr = ncpr                          ! intent(inout)
        Interstitial(nt)%ncps = ncps                          ! intent(inout)
        Interstitial(nt)%ncgl = ncgl                          ! intent(inout)
        !Tbd%phy_f3d(:,:,1)                                   ! intent(inout)
        Interstitial(nt)%qlcn = qlcn                          ! intent(inout)
        Interstitial(nt)%qicn = qicn                          ! intent(inout)
        Interstitial(nt)%cf_upi = cf_upi                      ! intent(inout)
        !Interstitial(nt)%clw(:,:,2)                          ! intent(out)
        !Interstitial(nt)%clw(:,:,1)                          ! intent(out)
        Interstitial(nt)%clcn = clcn                          ! intent(in)
        !cdata_block(nb,nt)%errmsg = errmsg                   ! intent(  out)
        !cdata_block(nb,nt)%errflg = errflg                   ! intent(  out)
        call ccpp_physics_run(cdata_block(nb,nt), scheme_name="m_micro_pre", ierr=ierr)
        ! Copy back intent(inout) and intent(out) interstitial variables to local variables in driver
        qrn = Interstitial(nt)%qrn
        qsnw = Interstitial(nt)%qsnw
        qgl = Interstitial(nt)%qgl
        ncpr = Interstitial(nt)%ncpr
        ncps = Interstitial(nt)%ncps
        ncgl = Interstitial(nt)%ncgl
        qlcn = Interstitial(nt)%qlcn
        qicn = Interstitial(nt)%qicn
        cf_upi = Interstitial(nt)%cf_upi
        clw(:,:,2) = Interstitial(nt)%clw(:,:,2)
        clw(:,:,1) = Interstitial(nt)%clw(:,:,1)
        errmsg = trim(cdata_block(nb,nt)%errmsg)
        errflg = cdata_block(nb,nt)%errflg
        if (errflg/=0) then
            write(0,*) 'Error in call to m_micro_pre_run: ' // trim(errmsg)
            stop
        end if
#else

        kk = 5
        if (Model%fprcp >= 2) kk = 6

!       Acheng used clw here for other code to run smoothly and minimum change
!       to make the code work. However, the nc and clw should be treated
!       in other procceses too.  August 28/2015; Hope that can be done next
!       year. I believe this will make the physical interaction more reasonable
!       Anning 12/5/2015 changed ntcw hold liquid only
        if (Model%do_shoc) then
          if (Model%fprcp == 0) then
            do k=1,levs
              do i=1,im
                clw(i,k,1) = Stateout%gq0(i,k,ntiw)             ! ice
                clw(i,k,2) = Stateout%gq0(i,k,ntcw)             ! water
                Tbd%phy_f3d(i,k,1) = Tbd%phy_f3d(i,k,ntot3d-2) ! clouds from shoc
              enddo
            enddo
          elseif (abs(Model%fprcp) == 1 .or. mg3_as_mg2) then
            do k=1,levs
              do i=1,im
                clw(i,k,1) = Stateout%gq0(i,k,ntiw)             ! ice
                clw(i,k,2) = Stateout%gq0(i,k,ntcw)             ! water
                qrn(i,k)   = Stateout%gq0(i,k,ntrw)
                qsnw(i,k)  = Stateout%gq0(i,k,ntsw)
                ncpr(i,k)  = Stateout%gq0(i,k,ntrnc)
                ncps(i,k)  = Stateout%gq0(i,k,ntsnc)
                Tbd%phy_f3d(i,k,1) = Tbd%phy_f3d(i,k,ntot3d-2) ! clouds from shoc
              enddo
            enddo
          else
            do k=1,levs
              do i=1,im
                clw(i,k,1) = Stateout%gq0(i,k,ntiw)             ! ice
                clw(i,k,2) = Stateout%gq0(i,k,ntcw)             ! water
                qrn(i,k)   = Stateout%gq0(i,k,ntrw)
                qsnw(i,k)  = Stateout%gq0(i,k,ntsw)
                qgl(i,k)   = Stateout%gq0(i,k,ntgl)
                ncpr(i,k)  = Stateout%gq0(i,k,ntrnc)
                ncps(i,k)  = Stateout%gq0(i,k,ntsnc)
                ncgl(i,k)  = Stateout%gq0(i,k,ntgnc)
                Tbd%phy_f3d(i,k,1) = Tbd%phy_f3d(i,k,ntot3d-2) ! clouds from shoc
              enddo
            enddo

          endif

        else
                                                     ! clouds from t-dt and cnvc
          if (Model%fprcp == 0 ) then
            do k=1,levs
              do i=1,im
                clw(i,k,1) = Stateout%gq0(i,k,ntiw)             ! ice
                clw(i,k,2) = Stateout%gq0(i,k,ntcw)             ! water
              enddo
            enddo
          elseif (abs(Model%fprcp) == 1 .or. mg3_as_mg2) then
            do k=1,levs
              do i=1,im
                clw(i,k,1) = Stateout%gq0(i,k,ntiw)             ! ice
                clw(i,k,2) = Stateout%gq0(i,k,ntcw)             ! water
                qrn(i,k)   = Stateout%gq0(i,k,ntrw)
                qsnw(i,k)  = Stateout%gq0(i,k,ntsw)
                ncpr(i,k)  = Stateout%gq0(i,k,ntrnc)
                ncps(i,k)  = Stateout%gq0(i,k,ntsnc)
              enddo
            enddo
          else
            do k=1,levs
              do i=1,im
                clw(i,k,1) = Stateout%gq0(i,k,ntiw)             ! ice
                clw(i,k,2) = Stateout%gq0(i,k,ntcw)             ! water
                qrn(i,k)   = Stateout%gq0(i,k,ntrw)
                qsnw(i,k)  = Stateout%gq0(i,k,ntsw)
                qgl(i,k)   = Stateout%gq0(i,k,ntgl)
                ncpr(i,k)  = Stateout%gq0(i,k,ntrnc)
                ncps(i,k)  = Stateout%gq0(i,k,ntsnc)
                ncgl(i,k)  = Stateout%gq0(i,k,ntgnc)
              enddo
            enddo
          endif
        endif

! add convective cloud fraction
        do k = 1,levs
          do i = 1,im
            Tbd%phy_f3d(i,k,1) = min(1.0, Tbd%phy_f3d(i,k,1) + clcn(i,k))
          enddo
        enddo

#endif

!       notice clw ix instead of im
!       call m_micro_driver(im,ix,levs,flipv,del,dtp,prsl,prsi,
!    &    prslk,prsik,pgr,vvl,clw(1,1,2), QLCN, clw(1,1,1),QICN,
!       if (lprnt) write(0,*)' cnv_mfdbef=',cnv_mfd(ipr,:),' flipv=',flipv
!       if(lprnt) write(0,*) ' befgt0=',Stateout%gt0(ipr,:),' kdt=',kdt
!       if(lprnt) write(0,*) ' befgq0=',Stateout%gq0(ipr,:,1),' kdt=',kdt
!       if(lprnt) write(0,*) ' befntlnc=',Stateout%gq0(ipr,:,ntlnc),' kdt=',kdt
!       if (lprnt) write(0,*)' clw1bef=',clw(ipr,:,1),' kdt=',kdt
!       if (lprnt) write(0,*)' clw2bef=',clw(ipr,:,2),' kdt=',kdt
!       if (lprnt) write(0,*)' qrnb=',qrn(ipr,:),' kdt=',kdt
!       if (lprnt) write(0,*)' qsnwb=',qsnw(ipr,:),' kdt=',kdt
!       if (lprnt) write(0,*)' qglb=',qgl(ipr,:),' kdt=',kdt
!       if (lprnt) write(0,*)' rhc=',rhc(ipr,:),' kdt=',kdt,' kk=',kk
!       if (lprnt) write(0,*)' cloudsb=',tbd%phy_f3d(ipr,:,1)*100,' kdt=',kdt
!       if (lprnt) write(0,*)' cloudsb=',Tbd%phy_f3d(ipr,:,1)*100,' kdt=',kdt
!       if (lprnt) write(0,*)' clcn=',clcn(ipr,:)*100,' kdt=',kdt
!       txa(:,:) = Stateout%gq0(:,:,1)
!       do k=1,levs
!       write(1000+me,*)' maxwatncb=',maxval(Stateout%gq0(1:im,k,ntlnc)),' k=',k,' kdt',kdt
!       enddo

#ifdef CCPP
        if (Model%me==0) write(0,*) 'CCPP DEBUG: calling m_micro_run through option B'
        !Interstitial(nt)%im = im                             ! intent(in) - set in Interstitial(nt)%create()
        !Interstitial(nt)%ix = ix                             ! intent(in) - set in Interstitial(nt)%create()
        !Model%levs                                           ! intent(in)
        !Model%flipv                                          ! intent(in)
        !Model%dtp                                            ! intent(in)
        !Statein%prsl                                         ! intent(in)
        !Statein%prsi                                         ! intent(in)
        !Statein%phil                                         ! intent(in)
        !Statein%phii                                         ! intent(in)
        !Statein%vvl                                          ! intent(in)
        Interstitial(nt)%clw(:,:,2) = clw(:,:,2)              ! intent(in)
        Interstitial(nt)%qlcn = qlcn                          ! intent(in)
        Interstitial(nt)%clw(:,:,1) = clw(:,:,1)              ! intent(in)
        Interstitial(nt)%qicn = qicn                          ! intent(in)
        !Radtend%htrlw                                        ! intent(in)
        !Radtend%htrsw                                        ! intent(in)
        Interstitial(nt)%w_upi = w_upi                        ! intent(in)
        Interstitial(nt)%cf_upi = cf_upi                      ! intent(in)
        Interstitial(nt)%frland = frland                      ! intent(in)
        !Diag%hpbl                                            ! intent(in)
        Interstitial(nt)%cnv_mfd = cnv_mfd                    ! intent(in)
        Interstitial(nt)%cnv_dqldt = cnv_dqldt                ! intent(in)
        Interstitial(nt)%clcn = clcn                          ! intent(in)
        !Stateout%gu0                                         ! intent(in)
        !Stateout%gv0                                         ! intent(in)
        !Diag%dusfc                                           ! intent(in)
        !Diag%dvsfc                                           ! intent(in)
        Interstitial(nt)%dusfc1 = dusfc1                      ! intent(in)
        Interstitial(nt)%dvsfc1 = dvsfc1                      ! intent(in)
        Interstitial(nt)%cnv_fice = cnv_fice                  ! intent(in)
        Interstitial(nt)%cnv_ndrop = cnv_ndrop                ! intent(in)
        Interstitial(nt)%cnv_nice = cnv_nice                  ! intent(in)
        !Stateout%gq0(:,:,1)                                  ! intent(inout)
        !Stateout%gq0(:,:,ntcw)                               ! intent(out)
        !Stateout%gq0(:,:,ntiw)                               ! intent(out)
        !Stateout%gt0                                         ! intent(inout)
        !Interstitial(nt)%prcpmp                              ! intent(out)
        !Sfcprop%sr                                           ! intent(out)
        !Stateout%gq0(:,:,ntlnc)                              ! intent(inout)
        !Stateout%gq0(:,:,ntinc)                              ! intent(inout)
        !Model%fprcp                                          ! intent(in)
        Interstitial(nt)%qrn = qrn                            ! intent(inout)
        Interstitial(nt)%qsnw = qsnw                          ! intent(inout)
        Interstitial(nt)%qgl = qgl                            ! intent(inout)
        Interstitial(nt)%ncpr = ncpr                          ! intent(inout)
        Interstitial(nt)%ncps = ncps                          ! intent(inout)
        Interstitial(nt)%ncgl = ncgl                          ! intent(inout)
        !Tbd%phy_f3d(1,1,1)                                   ! intent(inout)
        Interstitial(nt)%kbot = kbot                          ! intent(inout)
        !Tbd%phy_f3d(1,1,2)                                   ! intent(out)
        !Tbd%phy_f3d(1,1,3)                                   ! intent(out)
        !Tbd%phy_f3d(1,1,4)                                   ! intent(out)
        !Tbd%phy_f3d(1,1,5)                                   ! intent(out)
        !Tbd%phy_f3d(1,1,Model%ngeffr)                        ! intent(out)
        !Tbd%aer_nm                                           ! intent(in)
        !Model%aero_in                                        ! intent(in)
        !Tbd%in_nm                                            ! intent(in)
        !Tbd%ccn_nm                                           ! intent(in)
        !Model%iccn                                           ! intent(in)
        Interstitial(nt)%skip_macro = skip_macro              ! intent(in)
        !Model%lprnt                                          ! intent(in)
        !Model%mg_alf                                         ! intent(in)
        !Model%mg_qcmin                                       ! intent(in)
        !Model%pdfflag                                        ! intent(in)
        Interstitial(nt)%ipr = ipr                            ! intent(in)
        !Model%kdt                                            ! intent(in)
        !Grid%xlat                                            ! intent(in)
        !Grid%xlon                                            ! intent(in)
        Interstitial(nt)%rhc = rhc                            ! intent(in)
        !cdata_block(nb,nt)%errmsg = errmsg                   ! intent(  out)
        !cdata_block(nb,nt)%errflg = errflg                   ! intent(  out)
        call ccpp_physics_run(cdata_block(nb,nt), scheme_name="m_micro", ierr=ierr)
        ! Copy back intent(inout) and intent(out) interstitial variables to local variables in driver
        rain1 = Interstitial(nt)%prcpmp
        qrn = Interstitial(nt)%qrn
        qsnw = Interstitial(nt)%qsnw
        qgl = Interstitial(nt)%qgl
        ncpr = Interstitial(nt)%ncpr
        ncps = Interstitial(nt)%ncps
        ncgl = Interstitial(nt)%ncgl
        kbot = Interstitial(nt)%kbot
        errmsg = trim(cdata_block(nb,nt)%errmsg)
        errflg = cdata_block(nb,nt)%errflg
        if (errflg/=0) then
            write(0,*) 'Error in call to m_micro_run: ' // trim(errmsg)
            stop
        end if
#else

        call m_micro_driver (im, ix, levs, Model%flipv, dtp,  Statein%prsl,      &
                             Statein%prsi, Statein%phil, Statein%phii,           &
                             Statein%vvl, clw(1,1,2), QLCN, clw(1,1,1), QICN,    &
                             Radtend%htrlw, Radtend%htrsw, w_upi, cf_upi,        &
                             FRLAND, Diag%HPBL, CNV_MFD,           CNV_DQLDT,    &
                             CLCN, Stateout%gu0, Stateout%gv0, Diag%dusfc,       &
                             Diag%dvsfc, dusfc1, dvsfc1, dusfc1, dvsfc1,         &
                             CNV_FICE, CNV_NDROP, CNV_NICE, Stateout%gq0(1,1,1), &
                             Stateout%gq0(1,1,ntcw),                             &
                             Stateout%gq0(1,1,ntiw), Stateout%gt0, rain1,        &
                             Sfcprop%sr, Stateout%gq0(1,1,ntlnc),                &
                             Stateout%gq0(1,1,ntinc), Model%fprcp, qrn,          &
                             qsnw, qgl, ncpr, ncps, ncgl,                        &
                             Tbd%phy_f3d(1,1,1),  kbot,                          &
                             Tbd%phy_f3d(1,1,2),  Tbd%phy_f3d(1,1,3),            &
                             Tbd%phy_f3d(1,1,4),  Tbd%phy_f3d(1,1,5),            &
                             Tbd%phy_f3d(1,1,kk), Tbd%aer_nm,                    &
                             Model%aero_in, Tbd%in_nm, Tbd%ccn_nm, Model%iccn,   &
                             skip_macro,                 lprnt,                  &
                             Model%mg_alf, Model%mg_qcmin, Model%pdfflag,        &
                             ipr, kdt, Grid%xlat, Grid%xlon, rhc)
#endif

!     do k=1,levs
!     write(1000+me,*)' maxwatnca=',maxval(Stateout%gq0(1:im,k,ntlnc)),' k=',k,' kdt=',kdt
!     enddo
!     write(1000+me,*)' at kdt = ',kdt
!     tem = 1000.0

!     call moist_bud2(im,ix,ix,levs,me,kdt,con_g,tem,del,rain1 &
!    &,               txa, clw(1,1,2), clw(1,1,1)         &
!    &,           Stateout%gq0(1:ix,1:levs,ntrw),Stateout%gq0(1:ix,1:levs,ntsw)&
!    &,           Stateout%gq0(1:ix,1:levs,ntgl)                       &
!    &,           Stateout%gq0(1:ix,1:levs,1),Stateout%gq0(1:ix,1:levs,ntcw)   &
!    &,           Stateout%gq0(1:ix,1:levs,ntiw)                       &
!    &,           qrn, qsnw, qgl, ' m_micro  ', grid%xlon(1:im), grid%xlat(1:im))

!       if (lprnt) write(0,*) ' rain1=',rain1(ipr)*86400.0, &
!    &' rainc=',diag%rainc(ipr)*86400.0                        &
!    &,' cn_prc=',cn_prc(ipr),' cn_snr=',cn_snr(ipr),' kdt=',kdt
!       if(lprnt) write(0,*) ' aftgt0=',Stateout%gt0(ipr,:),' kdt=',kdt
!       if (lprnt) write(0,*) ' aftlsgq0=',stateout%gq0(ipr,:,1),' kdt=',kdt
!       if (lprnt) write(0,*)' clw1aft=',stateout%gq0(ipr,:,ntiw),' kdt=',kdt
!       if (ntgl > 0 .and. lprnt)  &
!                  write(0,*)' cgw1aft=',stateout%gq0(ipr,:,ntgl),' kdt=',kdt
!       if (lprnt) write(0,*)' cloudsm=',tbd%phy_f3d(ipr,:,1)*100,' kdt=',kdt
!       if (lprnt) write(0,*)' clw2aft=',stateout%gq0(ipr,:,ntcw),' kdt=',kdt
!       if (lprnt) write(0,*)' qrna=',qrn(ipr,:),' kdt=',kdt
!       if (lprnt) write(0,*)' qsnwa=',qsnw(ipr,:),' kdt=',kdt
!       if (lprnt) write(0,*)' qglba',qgl(ipr,:),' kdt=',kdt

#ifdef CCPP
        if (Model%me==0) write(0,*) 'CCPP DEBUG: calling m_micro_post_run through option B'
        !Interstitial(nt)%im = im                             ! intent(in) - set in Interstitial(nt)%create()
        !Model%levs                                           ! intent(in)
        !Model%fprcp                                          ! intent(in)
        !Model%mg3_as_mg2                                     ! intent(in)
        Interstitial(nt)%ncpr = ncpr                          ! intent(in)
        Interstitial(nt)%ncps = ncps                          ! intent(in)
        Interstitial(nt)%ncgl = ncgl                          ! intent(in)
        Interstitial(nt)%qrn = qrn                            ! intent(inout)
        Interstitial(nt)%qsnw = qsnw                          ! intent(inout)
        Interstitial(nt)%qgl = qgl                            ! intent(inout)
        !Stateout%gq0(:,:,ntrw)                               ! intent(out)
        !Stateout%gq0(:,:,ntsw)                               ! intent(out)
        !Stateout%gq0(:,:,ntgl)                               ! intent(out)
        !Stateout%gq0(:,:,ntrnc)                              ! intent(out)
        !Stateout%gq0(:,:,ntsnc)                              ! intent(out)
        !Stateout%gq0(:,:,ntgnc)                              ! intent(out)
        !cdata_block(nb,nt)%errmsg = errmsg                   ! intent(  out)
        !cdata_block(nb,nt)%errflg = errflg                   ! intent(  out)
        call ccpp_physics_run(cdata_block(nb,nt), scheme_name="m_micro_post", ierr=ierr)
        ! Copy back intent(inout) and intent(out) interstitial variables to local variables in driver
        qrn = Interstitial(nt)%qrn
        qsnw = Interstitial(nt)%qsnw
        qgl = Interstitial(nt)%qgl
        errmsg = trim(cdata_block(nb,nt)%errmsg)
        errflg = cdata_block(nb,nt)%errflg
        if (errflg/=0) then
            write(0,*) 'Error in call to m_micro_post_run: ' // trim(errmsg)
            stop
        end if
#else
        if (abs(Model%fprcp) == 1 .or. mg3_as_mg2) then
          do k=1,levs
            do i=1,im
              if (abs(qrn(i,k))  < qsmall) qrn(i,k)  = 0.0
              if (abs(qsnw(i,k)) < qsmall) qsnw(i,k) = 0.0
              Stateout%gq0(i,k,ntrw)  = qrn(i,k)
              Stateout%gq0(i,k,ntsw)  = qsnw(i,k)
              Stateout%gq0(i,k,ntrnc) = ncpr(i,k)
              Stateout%gq0(i,k,ntsnc) = ncps(i,k)
            enddo
          enddo
        elseif (Model%fprcp > 1) then
          do k=1,levs
            do i=1,im
              if (abs(qrn(i,k))  < qsmall) qrn(i,k)  = 0.0
              if (abs(qsnw(i,k)) < qsmall) qsnw(i,k) = 0.0
              if (abs(qgl(i,k))  < qsmall) qgl(i,k)  = 0.0
              Stateout%gq0(i,k,ntrw)  = qrn(i,k)
              Stateout%gq0(i,k,ntsw)  = qsnw(i,k)
              Stateout%gq0(i,k,ntgl)  = qgl(i,k)
              Stateout%gq0(i,k,ntrnc) = ncpr(i,k)
              Stateout%gq0(i,k,ntsnc) = ncps(i,k)
              Stateout%gq0(i,k,ntgnc) = ncgl(i,k)
            enddo
          enddo

        endif
#endif

!       if (lprnt) write(0,*)' cloudsm=',tbd%phy_f3d(ipr,:,1)*100,' kdt=',kdt
!       if (lprnt) write(0,*)' clw2aft=',stateout%gq0(ipr,:,ntcw),' kdt=',kdt
!       if (lprnt) write(0,*)' qrna=',qrn(ipr,:),' kdt=',kdt
!       if (lprnt) write(0,*)' qsnwa=',qsnw(ipr,:),' kdt=',kdt
!       if (lprnt) write(0,*)' qglba',qgl(ipr,:),' kdt=',kdt
!

        elseif (imp_physics == Model%imp_physics_gfdl) then      ! GFDL MP

#ifdef CCPP
          ! OPTION B - works with all compilers
          if (Model%me==0) write(0,*) 'CCPP DEBUG: calling gfdl_cloud_microphys through option B'
          ! Copy local variables from driver to appropriate interstitial variables
          !Model%levs                                  ! intent(in)
          !Interstitial(nt)%im                         ! intent(in) - set in Interstitial(nt)%create
          !con_g                                       ! intent(in) - physical constant in physcons.F90
          !con_fvirt                                   ! intent(in) - physical constant in physcons.F90
          !con_rd                                      ! intent(in) - physical constant in physcons.F90
          !Interstitial(nt)%frland                     ! intent(in) - set in GFS_suite_interstitital_1
          !Grid%area                                   ! intent(in)
          !Stateout%gq0(:,:,1      )                   ! intent(inout)
          !Stateout%gq0(:,:,ntcw   )                   ! intent(inout)
          !Stateout%gq0(:,:,ntrw   )                   ! intent(inout)
          !Stateout%gq0(:,:,ntiw   )                   ! intent(inout)
          !Stateout%gq0(:,:,ntsw   )                   ! intent(inout)
          !Stateout%gq0(:,:,ntgl   )                   ! intent(inout)
          !Stateout%gq0(:,:,ntclamt)                   ! intent(inout)
          !Stateout%gt0                                ! intent(inout)
          !Stateout%gu0                                ! intent(inout)
          !Stateout%gv0                                ! intent(inout)
          !Statein%vvl                                 ! intent(in)
          !Statein%prsl                                ! intent(in)
          !Statein%phii                                ! intent(in)
          Interstitial(nt)%del = del                   ! intent(in)
          !Interstitial(nt)%rainmp                     ! intent(out)
          !Interstitial(nt)%icemp                      ! intent(out)
          !Interstitial(nt)%snowmp                     ! intent(out)
          !Interstitial(nt)%graupelmp                  ! intent(out)
          !Interstitial(nt)%prcpmp                     ! intent(out)
          !Sfcprop%sr                                  ! intent(out)
          !Model%dtp                                   ! intent(in)
          !Model%hydrostatic                           ! intent(in)
          !Interstitial(nt)%phys_hydrostatic           ! intent(in) - set in Interstitial(nt)%create
          !Model%lradar                                ! intent(in)
          !Diag%refl_10cm                              ! intent(inout)
          !Model%effr_in                               ! intent(in)
          !Tbd%phy_f3d(1:im,1:levs,1)                  ! intent(inout)
          !Tbd%phy_f3d(1:im,1:levs,2)                  ! intent(inout)
          !Tbd%phy_f3d(1:im,1:levs,3)                  ! intent(inout)
          !Tbd%phy_f3d(1:im,1:levs,4)                  ! intent(inout)
          !Tbd%phy_f3d(1:im,1:levs,5)                  ! intent(inout)
          !cdata_block(nb,nt)%errmsg = errmsg          ! intent(out)
          !cdata_block(nb,nt)%errflg = errflg          ! intent(out)
          !
          call ccpp_physics_run(cdata_block(nb,nt), scheme_name="gfdl_cloud_microphys", ierr=ierr)
          ! Copy back intent(inout) interstitial variables to local variables
          rain1 = Interstitial(nt)%prcpmp              ! explicit precipitation (rain, ice, snow, graupel) on physics timestep
          errmsg = trim(cdata_block(nb,nt)%errmsg)
          errflg = cdata_block(nb,nt)%errflg
          !
          if (errflg/=0) then
              write(0,*) 'Error in call to gfdl_cloud_microphys: ' // trim(errmsg)
              stop
          end if
#else
          if (Model%me==0) write(0,*) 'CCPP DEBUG: calling non-CCPP compliant version of gfdl_cloud_microphys_driver'
          do i = 1, im
            land     (i,1)   = frland(i)
            area     (i,1)   = Grid%area(i)
            rain0    (i,1)   = 0.0
            snow0    (i,1)   = 0.0
            ice0     (i,1)   = 0.0
            graupel0 (i,1)   = 0.0
          enddo
          do k = 1, levs
            kk = levs-k+1
            do i = 1, im
              qn1  (i,1,k) = 0.0
              qv_dt(i,1,k) = 0.0
              ql_dt(i,1,k) = 0.0
              qr_dt(i,1,k) = 0.0
              qi_dt(i,1,k) = 0.0
              qs_dt(i,1,k) = 0.0
              qg_dt(i,1,k) = 0.0
              qa_dt(i,1,k) = 0.0
              pt_dt(i,1,k) = 0.0
              udt  (i,1,k) = 0.0
              vdt  (i,1,k) = 0.0
!
              qv1  (i,1,k) =  Stateout%gq0(i,kk,1)
              ql1  (i,1,k) =  Stateout%gq0(i,kk,ntcw)
              qr1  (i,1,k) =  Stateout%gq0(i,kk,ntrw)
              qi1  (i,1,k) =  Stateout%gq0(i,kk,ntiw)
              qs1  (i,1,k) =  Stateout%gq0(i,kk,ntsw)
              qg1  (i,1,k) =  Stateout%gq0(i,kk,ntgl)
              qa1  (i,1,k) =  Stateout%gq0(i,kk,ntclamt)
              pt   (i,1,k) =  Stateout%gt0(i,kk)
              w    (i,1,k) = -Statein%vvl(i,kk)*(one+con_fvirt*qv1(i,1,k))   &
                           *  Stateout%gt0(i,kk) / Statein%prsl(i,kk) * (con_rd*onebg)
              uin  (i,1,k) =  Stateout%gu0(i,kk)
              vin  (i,1,k) =  Stateout%gv0(i,kk)
              delp (i,1,k) =  del(i,kk)
              dz   (i,1,k) = (Statein%phii(i,kk)-Statein%phii(i,kk+1)) * onebg
              p123 (i,1,k) = Statein%prsl(i,kk)
              refl (i,1,k) = Diag%refl_10cm(i,kk)
            enddo
          enddo

          call gfdl_cloud_microphys_driver(qv1, ql1, qr1, qi1, qs1, qg1, qa1, &
                                           qn1, qv_dt, ql_dt, qr_dt, qi_dt,   &
                                           qs_dt, qg_dt, qa_dt, pt_dt, pt, w, &
                                           uin, vin, udt, vdt, dz, delp,      &
                                           area, dtp, land, rain0, snow0,     &
                                           ice0, graupel0, .false., .true.,   &
                                           1, im, 1, 1, 1, levs, 1, levs,     &
                                           seconds,p123,Model%lradar,refl,kdt,&
                                           nsteps_per_reset)
          tem = dtp * con_p001 / con_day
          do i = 1, im
!            rain0(i,1) = max(con_d00, rain0(i,1))
!            snow0(i,1) = max(con_d00, snow0(i,1))
!            ice0(i,1)  = max(con_d00, ice0(i,1))
!            graupel0(i,1)  = max(con_d00, graupel0(i,1))
            if(rain0(i,1)*tem < rainmin) then
              rain0(i,1) = 0.0
            endif
            if(ice0(i,1)*tem < rainmin) then
              ice0(i,1) = 0.0
            endif
            if(snow0(i,1)*tem < rainmin) then
              snow0(i,1) = 0.0
            endif
            if(graupel0(i,1)*tem < rainmin) then
              graupel0(i,1) = 0.0
            endif

            rain1(i)        = (rain0(i,1)+snow0(i,1)+ice0(i,1)+graupel0(i,1)) * tem
            Diag%ice(i)     = ice0    (i,1) * tem
            Diag%snow(i)    = snow0   (i,1) * tem
            Diag%graupel(i) = graupel0(i,1) * tem
            if ( rain1(i) > rainmin ) then
              Sfcprop%sr(i) = (snow0(i,1) + ice0(i,1)  + graupel0(i,1)) &
                            / (rain0(i,1) + snow0(i,1) + ice0(i,1) + graupel0(i,1))
            else
              Sfcprop%sr(i) = 0.0
            endif
          enddo
#ifndef WORKAROUND_SRFLAG
          ! DH* Convert rain0, ice0, graupel0 and snow0 from mm/day to m/physics-timestep
          ! for later use (approx. lines 7970, calculation of srflag)
          rain0 = tem*rain0
          ice0  = tem*ice0
          snow0 = tem*snow0
          graupel0 = tem*graupel0
          ! *DH
#endif
          do k = 1, levs
            kk = levs-k+1
            do i=1,im
              Stateout%gq0(i,k,1   )    = qv1(i,1,kk)       + qv_dt(i,1,kk) * dtp
              Stateout%gq0(i,k,ntcw)    = ql1(i,1,kk)       + ql_dt(i,1,kk) * dtp
              Stateout%gq0(i,k,ntrw)    = qr1(i,1,kk)       + qr_dt(i,1,kk) * dtp
              Stateout%gq0(i,k,ntiw)    = qi1(i,1,kk)       + qi_dt(i,1,kk) * dtp
              Stateout%gq0(i,k,ntsw)    = qs1(i,1,kk)       + qs_dt(i,1,kk) * dtp
              Stateout%gq0(i,k,ntgl)    = qg1(i,1,kk)       + qg_dt(i,1,kk) * dtp
              Stateout%gq0(i,k,ntclamt) = qa1(i,1,kk)       + qa_dt(i,1,kk) * dtp
              Stateout%gt0(i,k)         = Stateout%gt0(i,k) + pt_dt(i,1,kk) * dtp
              Stateout%gu0(i,k)         = Stateout%gu0(i,k) + udt  (i,1,kk) * dtp
              Stateout%gv0(i,k)         = Stateout%gv0(i,k) + vdt  (i,1,kk) * dtp
              Diag%refl_10cm(i,k)       = refl(i,1,kk)
            enddo

            if(Model%effr_in) then 
              do i =1, im
                den(i,k)=0.622*Statein%prsl(i,k)/ &
                      (con_rd*Stateout%gt0(i,k)*(Stateout%gq0(i,k,1)+0.622))
              enddo
            endif 
          enddo
! DH* moved this down to the other maximum hourly diagnpstics calculations
!!Calculate hourly max 1-km agl and -10C reflectivity
!          if(Model%lradar .and. (imp_physics == 11 .or. imp_physics == 8))then
!            allocate(refd(im))
!            allocate(refd263k(im))
!            call max_fields(Statein%phil,Diag%refl_10cm,con_g,im,levs,refd,Stateout%gt0,refd263k)
!            do i=1,im
!              if(mod(kdtminus1,nsteps_per_reset)==0)then
!                Diag%refdmax(I) = -35.
!                Diag%refdmax263k(I) = -35.
!              endif
!              Diag%refdmax(i) = max(Diag%refdmax(i),refd(i))
!              Diag%refdmax263k(i) = max(Diag%refdmax263k(i),refd263k(i))
!            enddo
!            deallocate (refd) 
!            deallocate (refd263k)
!          endif
! *DH
!
          if(Model%effr_in) then 
            call cloud_diagnosis (1, im, 1, levs, den(1:im,1:levs),             & 
               Stateout%gq0(1:im,1:levs,ntcw), Stateout%gq0(1:im,1:levs,ntiw),  & 
               Stateout%gq0(1:im,1:levs,ntrw), Stateout%gq0(1:im,1:levs,ntsw),  & 
               Stateout%gq0(1:im,1:levs,ntgl), Stateout%gt0(1:im,1:levs),       &
               Tbd%phy_f3d(1:im,1:levs,1),     Tbd%phy_f3d(1:im,1:levs,2),      & 
               Tbd%phy_f3d(1:im,1:levs,3),     Tbd%phy_f3d(1:im,1:levs,4),      & 
               Tbd%phy_f3d(1:im,1:levs,5))

!            do k = 1, levs
!              do i=1,im
!
!                if(Model%me==0) then
!		  if(Tbd%phy_f3d(i,k,1) > 5.) then 
!                    write(6,*) 'phy driver:cloud radii:',Model%kdt, i,k,        &
!				Tbd%phy_f3d(i,k,1)
!                  endif 
!		  if(Tbd%phy_f3d(i,k,3)> 0.0) then 
!                    write(6,*) 'phy driver:rain radii:',Model%kdt, i,k,         & 
!				Tbd%phy_f3d(i,k,3)
!                  endif 
!
!                endif 
!              enddo 
!            enddo 

          endif 
#endif
        endif  ! end of if(Model%imp_physics)
      endif    ! end if_ncld

!     if (lprnt) write(0,*)' rain1 after ls=',rain1(ipr)
!
#ifdef CCPP
      if (Model%cscnv .and. Model%do_aw) then
         if (Model%me==0) write(0,*) 'CCPP DEBUG: calling cs_conv_aw_adj through option B'
         ! Copy local variables from driver to appropriate interstitial variables
         !Interstitial(nt)%im                         ! intent(in) - set in Interstitial(nt)%create
         !Model%levs                                  ! intent(in)
         !Model%cscnv                                 ! intent(in)
         !Model%do_aw                                 ! intent(in)
         !Model%do_shoc                               ! intent(in)
         !Model%ntrac                                 ! intent(in)
         !Model%ncld                                  ! intent(in)
         !Model%ntcw                                  ! intent(in)
         !Model%ntclamt                               ! intent(in)
         !Interstitial(nt)%nncl                       ! intent(in) - set in Interstitial(nt)%create
         !con_g                                       ! intent(in) - physical constant in physcons.F90
         Interstitial(nt)%sigmafrac = sigmafrac       ! intent(in)
         !Stateout%gt0                                ! intent(inout)
         !Stateout%gq0                                ! intent(inout)
         Interstitial(nt)%save_t = dtdt               ! intent(in)
         Interstitial(nt)%save_q = dqdt               ! intent(in)
         !Statein(nb)%prsi                            ! intent(in)
         !Tbd(nb)%phy_f3d(:,:,1)                      ! intent(inout)
         !Tbd(nb)%phy_f3d(:,:,IPD_Control%ntot3d-2)   ! intent(inout)
         Interstitial(nt)%prcpmp = rain1              ! intent(inout)
         !cdata_block(nb,nt)%errmsg = errmsg          ! intent(out)
         !cdata_block(nb,nt)%errflg = errflg          ! intent(out)
         !
         call ccpp_physics_run(cdata_block(nb,nt), scheme_name="cs_conv_aw_adj", ierr=ierr)
         ! Copy back intent(inout) interstitial variables to local variables
         rain1  = Interstitial(nt)%prcpmp
         errmsg = trim(cdata_block(nb,nt)%errmsg)
         errflg = cdata_block(nb,nt)%errflg
         !
         if (errflg/=0) then
            write(0,*) 'Error in call to cs_conv_aw_adj: ' // trim(errmsg)
            stop
         end if
      end if
#else
      if (Model%cscnv .and. Model%do_aw) then
!  Arakawa-Wu adjustment of large-scale microphysics tendencies:
!  reduce by factor of (1-sigma)
!  these are microphysics increments. We want to keep (1-sigma) of the increment,
!  we will remove sigma*increment from final values
!         fsigma = 0.  ! don't apply any AW correction, in addition comment next line
!         fsigma = sigmafrac

!  adjust sfc rainrate for conservation
!  vertically integrate reduction of water increments, reduce precip by that amount

        temrain1(:) = 0.0
        do k = 1,levs
          do i = 1,im
            tem1               = sigmafrac(i,k)
            Stateout%gt0(i,k)  = Stateout%gt0(i,k) - tem1 * (Stateout%gt0(i,k)-dtdt(i,k))
            tem2                = tem1 * (Stateout%gq0(i,k,1)-dqdt(i,k,1))
            Stateout%gq0(i,k,1) = Stateout%gq0(i,k,1) - tem2
            temrain1(i) = temrain1(i) - (Statein%prsi(i,k)-Statein%prsi(i,k+1)) &
                                      * tem2 * onebg
          enddo
        enddo

!     if (lprnt) write(0,*)' gt0aftpraw=',Stateout%gt0(ipr,:),' kdt=',kdt,'me=',me
        do n=ntcw,ntcw+nncl-1
          do k = 1,levs
            do i = 1,im
              tem1                = sigmafrac(i,k) * (Stateout%gq0(i,k,n)-dqdt(i,k,n))
              Stateout%gq0(i,k,n) = Stateout%gq0(i,k,n) - tem1
              temrain1(i)         = temrain1(i) - (Statein%prsi(i,k)-Statein%prsi(i,k+1)) &
                                        * tem1 * onebg
            enddo
          enddo
        enddo
!     write(1000+me,*)' rain1=',rain1(4),' temrain1=',temrain1(i)*0.001
        do i = 1,im
          rain1(i) = max(rain1(i) - temrain1(i)*0.001, 0.0_kind_phys)
        enddo
      endif
#endif

#ifdef CCPP
      ! OPTION B - works with all compilers
      if (Model%me==0) write(0,*) 'CCPP DEBUG: calling GFS_MP_generic_post through option B'
      !Interstitial(nt)%im = im                   ! intent(in) - set in Interstitial(nt)%create()
      !Interstitial(nt)%ix = im                   ! intent(in) - set in Interstitial(nt)%create()
      !Model%levs                                 ! intent(in)
      !Model%kdt                                  ! intent(in)
      !Model%nrcm                                 ! intent(in)
      !Model%ncld                                 ! intent(in)
      !Model%nncl                                 ! intent(in)
      !Model%ntcw                                 ! intent(in)
      !Model%ntrac                                ! intent(in)
      !Model%imp_physics                          ! intent(in)
      !Model%imp_physics_gfdl                     ! intent(in)
      !Model%imp_physics_thompson                 ! intent(in)
      !Model%cal_pre                              ! intent(in)
      !Model%lssav                                ! intent(in)
      !Model%ldiag3d                              ! intent(in)
      !Model%cplflx                               ! intent(in)
      !Model%cplchm                               ! intent(in)
      !con_g                                      ! intent(in) - in physcons
      !Model%dtf                                  ! intent(in)
      Interstitial(nt)%frain = frain              ! intent(in)
      !Diag%rainc                                 ! intent(in)
      Interstitial(nt)%prcpmp = rain1             ! intent(in)
      !Tbd%rann                                   ! intent(in)
      !Grid%xlat                                  ! intent(in)
      !Grid%xlon                                  ! intent(in)
      !Stateout%gt0                               ! intent(in)
      !Stateout%gq0                               ! intent(in)
      !Statein%prsl                               ! intent(in)
      !Statein%prsi                               ! intent(in)
      !Statein%phii                               ! intent(in)
      !Sfcprop%tsfc                               ! intent(in)
      !Diag%ice                                   ! intent(in)
      !Diag%snow                                  ! intent(in)
      !Diag%graupel                               ! intent(in)
      Interstitial(nt)%save_t = dtdt              ! intent(in)
      Interstitial(nt)%save_q(:,:,1) = dqdt(:,:,1)! intent(in)
      !if (imp_physics == 11) then
      !  Interstitial(nt)%icemp = ice0(:,1)           ! intent(in) ice0 is not allocated if CCPP is being used - coming straight from gfdl_cloud_microphys
      !  Interstitial(nt)%snowmp = snow0(:,1)         ! intent(in) snow0 is not allocated if CCPP is being used - coming straight from gfdl_cloud_microphys
      !  Interstitial(nt)%graupelmp = graupel0(:,1)   ! intent(in) graupel0 is not allocated if CCPP is being used - coming straight from gfdl_cloud_microphys
      !end if
      Interstitial(nt)%del = del                  ! intent(in)
      !Diag%rain                                  ! intent(inout)
      !Diag%tdomr                                 ! intent(inout)
      !Diag%tdomzr                                ! intent(inout)
      !Diag%tdomip                                ! intent(inout)
      !Diag%tdoms                                 ! intent(inout)
      !Sfcprop%tprcp                              ! intent(inout)
      !Sfcprop%srflag                             ! intent(inout)
      !Diag%totprcp                               ! intent(inout)
      !Diag%totice                                ! intent(inout)
      !Diag%totsnw                                ! intent(inout)
      !Diag%totgrp                                ! intent(inout)
      !Diag%totprcpb                              ! intent(inout)
      !Diag%toticeb                               ! intent(inout)
      !Diag%totsnwb                               ! intent(inout)
      !Diag%totgrpb                               ! intent(inout)
      !Diag%dt3dt(:,:,6)                          ! intent(inout)
      !Diag%dq3dt(:,:,4)                          ! intent(inout)
      !Coupling%rain_cpl                          ! intent(inout)
      !Coupling%rainc_cpl                         ! intent(inout)
      !Coupling%snow_cpl                          ! intent(inout)
      !Diag%pwat                                  ! intent(inout)
      !Model%do_sppt                              ! intent(in)
      !Tbd%dtdtr                                  ! intent(inout)
      Interstitial(nt)%dtdtc = dtdtc              ! intent(in)
      !Tbd%drain_cpl                              ! intent(inout)
      !Tbd%dsnow_cpl                              ! intent(inout)
      !cdata_block(nb,nt)%errmsg = errmsg         ! intent(out)
      !cdata_block(nb,nt)%errflg = errflg         ! intent(out)
      call ccpp_physics_run(cdata_block(nb,nt), scheme_name="GFS_MP_generic_post", ierr=ierr)
      ! Copy back intent(inout) interstitial variables to local variables
      errmsg = trim(cdata_block(nb,nt)%errmsg)
      errflg = cdata_block(nb,nt)%errflg
      if (errflg/=0) then
          write(0,*) 'Error in call to GFS_MP_generic_post: ' // trim(errmsg)
          stop
      end if
#else
      Diag%rain(:) = Diag%rainc(:) + frain * rain1(:)

      if (Model%cal_pre) then       ! hchuang: add dominant precipitation type algorithm
!
        call calpreciptype (kdt, Model%nrcm, im, ix, levs, levs+1,           &
                            Tbd%rann, Grid%xlat, Grid%xlon, Stateout%gt0,    &
                            Stateout%gq0, Statein%prsl, Statein%prsi,        &
                            Diag%rain, Statein%phii, Sfcprop%tsfc,           &  !input
                            domr, domzr, domip, doms)                           ! output
!
!        if (lprnt) print*,'debug calpreciptype: DOMR,DOMZR,DOMIP,DOMS '
!     &,DOMR(ipr),DOMZR(ipr),DOMIP(ipr),DOMS(ipr)
!        do i=1,im
!         if (abs(xlon(i)*57.29578-114.0) .lt. 0.2  .and.
!     &    abs(xlat(i)*57.29578-40.0) .lt. 0.2)
!     &    print*,'debug calpreciptype: DOMR,DOMZR,DOMIP,DOMS ',
!     &    DOMR(i),DOMZR(i),DOMIP(i),DOMS(i)
!       end do
!       HCHUANG: use new precipitation type to decide snow flag for LSM snow accumulation

        if (Model%imp_physics /= Model%imp_physics_gfdl) then
          do i=1,im
            Sfcprop%tprcp(i)  = max(0.0, Diag%rain(i) )
            if(doms(i) > 0.0 .or. domip(i) > 0.0) then
              Sfcprop%srflag(i) = 1.
            else
              Sfcprop%srflag(i) = 0.
            end if
          enddo
        endif

      endif

      if (Model%lssav) then
!        if (Model%me == 0) print *,'in phys drive, kdt=',Model%kdt, &
!          'totprcpb=', Diag%totprcpb(1),'totprcp=',Diag%totprcp(1), &
!          'rain=',Diag%rain(1)
        do i=1,im
          Diag%totprcp (i) = Diag%totprcp (i) + Diag%rain(i)
          Diag%totice  (i) = Diag%totice  (i) + Diag%ice(i)
          Diag%totsnw  (i) = Diag%totsnw  (i) + Diag%snow(i)
          Diag%totgrp  (i) = Diag%totgrp  (i) + Diag%graupel(i)
          Diag%totprcpb(i) = Diag%totprcpb(i) + Diag%rain(i)
          Diag%toticeb (i) = Diag%toticeb (i) + Diag%ice(i)
          Diag%totsnwb (i) = Diag%totsnwb (i) + Diag%snow(i)
          Diag%totgrpb (i) = Diag%totgrpb (i) + Diag%graupel(i)
!
          if (Model%cal_pre) then
            Diag%tdomr(i)  = Diag%tdomr(i)  + domr(i)  * dtf
            Diag%tdomzr(i) = Diag%tdomzr(i) + domzr(i) * dtf
            Diag%tdomip(i) = Diag%tdomip(i) + domip(i) * dtf
            Diag%tdoms(i)  = Diag%tdoms(i)  + doms(i)  * dtf
          endif
        enddo

        if (Model%ldiag3d) then
          do k=1,levs
            do i=1,im
              Diag%dt3dt(i,k,6) = Diag%dt3dt(i,k,6) + (Stateout%gt0(i,k)-dtdt(i,k)) * frain
!              Diag%dq3dt(i,k,4) = Diag%dq3dt(i,k,4) + (Stateout%gq0(i,k,1)-dqdt(i,k,1)) * frain
            enddo
          enddo
        endif
      endif

!  --- ...  estimate t850 for rain-snow decision

      t850(1:im) = Stateout%gt0(1:im,1)

      do k = 1, levs-1
        do i = 1, im
          if (Statein%prsl(i,k) > p850 .and. Statein%prsl(i,k+1) <= p850) then
            t850(i) = Stateout%gt0(i,k) - (Statein%prsl(i,k)-p850) / &
                      (Statein%prsl(i,k)-Statein%prsl(i,k+1)) *      &
                      (Stateout%gt0(i,k)-Stateout%gt0(i,k+1))
          endif
        enddo
      enddo

      if (Model%imp_physics == Model%imp_physics_gfdl) then
! determine convective rain/snow by surface temperature
! determine large-scale rain/snow by rain/snow coming out directly from MP
        tem = dtp * con_p001 / con_day
        do i = 1, im
          Sfcprop%tprcp(i)  = max(0.0, Diag%rain(i) )! clu: rain -> tprcp
          Sfcprop%srflag(i) = 0.                     ! clu: default srflag as 'rain' (i.e. 0)
          if (Sfcprop%tsfc(i) >= 273.15) then
            crain = Diag%rainc(i)
            csnow = 0.0
          else
            crain = 0.0
            csnow = Diag%rainc(i)
          endif
!         if (snow0(i,1)+ice0(i,1)+graupel0(i,1)+csnow > rain0(i,1)+crain) then
!          if (snow0(i,1)+ice0(i,1)+graupel0(i,1)+csnow > 0.0) then
!            Sfcprop%srflag(i) = 1.                   ! clu: set srflag to 'snow' (i.e. 1)
!          endif
! compute fractional srflag
          total_precip = snow0(i,1)+ice0(i,1)+graupel0(i,1)+rain0(i,1)+Diag%rainc(i)
          if (total_precip*tem > rainmin) then
            Sfcprop%srflag(i) = (snow0(i,1)+ice0(i,1)+graupel0(i,1)+csnow)/total_precip
          endif
        enddo
      elseif( .not. Model%cal_pre) then
        do i = 1, im
          Sfcprop%tprcp(i)  = max(0.0, Diag%rain(i) )! clu: rain -> tprcp
          Sfcprop%srflag(i) = 0.                     ! clu: default srflag as 'rain' (i.e. 0)
          if (t850(i) <= 273.16) then
            Sfcprop%srflag(i) = 1.                   ! clu: set srflag to 'snow' (i.e. 1)
          endif
        enddo
      endif

!  --- ...  coupling insertion

      if (Model%cplflx) then
        do i = 1, im
          if (t850(i) > 273.16) then
            Coupling%rain_cpl(i) = Coupling%rain_cpl(i) + Diag%rain(i)
          else
            Coupling%snow_cpl(i) = Coupling%snow_cpl(i) + Diag%rain(i)
          endif
        enddo
      endif

      if (Model%cplchm.and. .not. Model%cplflx) then
        do i = 1, im
          Coupling%rain_cpl(i)  = Coupling%rain_cpl(i)  + Diag%rain(i)
          Coupling%rainc_cpl(i) = Coupling%rainc_cpl(i) + Diag%rainc(i)
        enddo
      endif

!  --- ...  calculate column precipitable water "pwat"
      Diag%pwat(:) = 0.0
      do k = 1, levs
        do i=1,im
          work1(i) = 0.0
        enddo
        if (ncld > 0) then
          do ic = ntcw, ntcw+nncl-1
            do i=1,im
              work1(i) = work1(i) + Stateout%gq0(i,k,ic)
            enddo
          enddo
        endif
        do i=1,im
          Diag%pwat(i) = Diag%pwat(i) + del(i,k)*(Stateout%gq0(i,k,1)+work1(i))
        enddo
!     if (lprnt .and. i == ipr) write(0,*)' gq0=',
!    &gq0(i,k,1),' qgrs=',qgrs(i,k,1),' work2=',work2(i),' k=',k
      enddo
      do i=1,im
        Diag%pwat(i) = Diag%pwat(i) * onebg
      enddo
#endif
!  --- ...  end coupling insertion

!!! update surface diagnosis fields at the end of phys package
!!! this change allows gocart to use filtered wind fields
!!!
      if (Model%lgocart) then
#ifdef CCPP
        if (Model%me==0) write(0,*) 'CCPP DEBUG: calling sfc_diag through option B'
        ! Copy local variables from driver to appropriate interstitial variables
        !Interstitial(nt)%im                 ! intent(in ) - set in Interstitial(nt)%create()
        !con_g, con_cp, con_eps, con_epsm1   ! intent(in ) - physical constants from physcons.F90
        !Statein%pgr                         ! intent(in )
        !Stateout%gu0(:,1)                   ! intent(in )
        !Stateout%gv0(:,1)                   ! intent(in )
        !Stateout%gt0(:,1)                   ! intent(in )
        !Stateout%gq0(:,1,1)                 ! intent(in )
        !Sfcprop%tsfc                        ! intent(in )
        Interstitial(nt)%qss = qss           ! intent(in )
        !Sfcprop%f10m                        ! intent(out)
        !Diag%u10m                           ! intent(out)
        !Diag%v10m                           ! intent(out)
        !Sfcprop%t2m                         ! intent(out)
        !Sfcprop%q2m                         ! intent(out)
        Interstitial(nt)%work3 = work3       ! intent(in )
        Interstitial(nt)%evap  = evap        ! intent(in )
        !Sfcprop%ffmm                        ! intent(in )
        !Sfcprop%ffhh                        ! intent(in )
        Interstitial(nt)%fm10 = fm10         ! intent(in )
        Interstitial(nt)%fh2 = fh2           ! intent(in )
        !cdata_block(nb,nt)%errmsg = errmsg  ! intent(out)
        !cdata_block(nb,nt)%errflg = errflg  ! intent(out)
        call ccpp_physics_run(cdata_block(nb,nt), scheme_name="sfc_diag", ierr=ierr)
        ! Copy intent(inout) and intent(out) interstitial variables to local variables in driver
        errmsg     = trim(cdata_block(nb,nt)%errmsg)
        errflg     = cdata_block(nb,nt)%errflg
        if (errflg/=0) then
            write(0,*) 'Error in call to sfc_diag: ' // trim(errmsg)
            stop
        end if
        if (Model%me==0) write(0,*) 'CCPP DEBUG: calling sfc_diag_post through option B'
        ! Copy local variables from driver to appropriate interstitial variables
        !Interstitial(nt)%im = im               ! intent(in) - set in Interstitial(nt)%create()
        !Model%lssav                            ! intent(in)
        !con_eps                                ! intent(in) - physical constant in physcons
        !con_epsm1                              ! intent(in) - physical constant in physcons
        !Sfcprop%t2m                            ! intent(in)
        !Sfcprop%q2m                            ! intent(in)
        !Sfcprop%u10m                           ! intent(in)
        !Sfcprop%v10m                           ! intent(in)
        !Diag%tmpmin                            ! intent(inout)
        !Diag%tmpmax                            ! intent(inout)
        !Diag%spfhmin                           ! intent(inout)
        !Diag%spfhmax                           ! intent(inout)
        !Diag%wind10mmax                        ! intent(inout)
        !Diag%u10mmax                           ! intent(inout)
        !Diag%v10mmax                           ! intent(inout)
        !Diag%dpt2m                             ! intent(inout)
        !cdata_block(nb,nt)%errmsg = errmsg     ! intent(out)
        !cdata_block(nb,nt)%errflg = errflg     ! intent(out)
        call ccpp_physics_run(cdata_block(nb,nt), scheme_name="sfc_diag_post", ierr=ierr)
        ! Copy intent(inout) and intent(out) interstitial variables to local variables in driver
        errmsg     = trim(cdata_block(nb,nt)%errmsg)
        errflg     = cdata_block(nb,nt)%errflg
        if (errflg/=0) then
            write(0,*) 'Error in call to sfc_diag_post: ' // trim(errmsg)
            stop
        end if
#else
!        ! DH* 20180817 - see my comment further up (around gwdc_pre) that
        ! work3 for the non-CCPP code is incorrect (CCPP code is correct) *DH
        if (Model%me==0) write(0,*) 'CCPP DEBUG: calling non-CCPP compliant version of sfc_diag'
        call sfc_diag (im, Statein%pgr, Stateout%gu0, Stateout%gv0,     &
                       Stateout%gt0, Stateout%gq0, Sfcprop%tsfc, qss,   &
                       Sfcprop%f10m, Diag%u10m, Diag%v10m, Sfcprop%t2m, &
                       Sfcprop%q2m,  work3, evap, Sfcprop%ffmm,         &
                       Sfcprop%ffhh, fm10, fh2)

        if (Model%lssav) then
          do i=1,im
            Diag%tmpmax (i) = max(Diag%tmpmax (i),Sfcprop%t2m(i))
            Diag%tmpmin (i) = min(Diag%tmpmin (i),Sfcprop%t2m(i))
            Diag%spfhmax(i) = max(Diag%spfhmax(i),Sfcprop%q2m(i))
            Diag%spfhmin(i) = min(Diag%spfhmin(i),Sfcprop%q2m(i))
          enddo
          !find max wind speed then decompose
          do i=1, im
             tem = sqrt(Diag%u10m(i)**2 + Diag%v10m(i)**2 )
             if (tem > Diag%wind10mmax(i)) then
                Diag%wind10mmax(i) = tem
                Diag%u10mmax(i)    = Diag%u10m(i)
                Diag%v10mmax(i)    = Diag%v10m(i)
             endif
           !Compute dew point, first using vapor pressure
           tem = max(Statein%pgr(i) * Sfcprop%q2m(i) / ( con_eps - con_epsm1 * Sfcprop%q2m(i)), 1.e-8)
           Diag%dpt2m(i) = 243.5 / ( ( 17.67 / log(tem/611.2) ) - 1.) + 273.14
          enddo
        endif
#endif
      endif

#ifndef CCPP
      ! CCPP: this code is now in GFS_surface_generic_post
!  --- ...  total runoff is composed of drainage into water table and
!           runoff at the surface and is accumulated in unit of meters
      if (Model%lssav) then
        tem = dtf * 0.001
        do i=1,im
          Diag%runoff(i)  = Diag%runoff(i)  + (drain(i)+runof(i)) * tem
          Diag%srunoff(i) = Diag%srunoff(i) + runof(i) * tem
        enddo
      endif
#endif

      if (Model%lsm == Model%lsm_noah) then  ! DH* Not for RUC - just for Noah LSM? Should we implement a namelist option for this?
!  --- ...  xw: return updated ice thickness & concentration to global array
#ifdef CCPP
        if (Model%me==0) write(0,*) 'CCPP DEBUG: calling sfc_sice_post through option B'
        ! Copy local variables from driver to appropriate interstitial variables
        !Interstitial(nt)%im = im             ! intent(in) - set in Interstitial(nt)%create
        Interstitial(nt)%islmsk = islmsk      ! intent(in)
        !Sfcprop%tsfc                         ! intent(in)
        !Sfcprop%fice                         ! intent(inout)
        !Sfcprop%hice                         ! intent(inout)
        !Sfcprop%tisfc                        ! intent(inout)
        !cdata_block(nb,nt)%errmsg = errmsg   ! intent(out)
        !cdata_block(nb,nt)%errflg = errflg   ! intent(out)
        call ccpp_physics_run(cdata_block(nb,nt), scheme_name="sfc_sice_post", ierr=ierr)
        ! Copy back intent(inout) interstitial variables to local variables in driver
        errmsg = trim(cdata_block(nb,nt)%errmsg)
        errflg = cdata_block(nb,nt)%errflg
        if (errflg/=0) then
          write(0,*) 'Error in call to sfc_sice_post: ' //trim(errmsg)
          stop
        end if
#else
        do i = 1, im
          if (islmsk(i) == 2) then
            Sfcprop%hice(i)  = zice(i)
            Sfcprop%fice(i)  = cice(i)
            Sfcprop%tisfc(i) = tice(i)
          else
            Sfcprop%hice(i)  = 0.0
            Sfcprop%fice(i)  = 0.0
            Sfcprop%tisfc(i) = Sfcprop%tsfc(i)
          endif
        enddo
#endif
       endif ! lsm == lsm_noah for sea ice

!     tem = dtf * 0.03456 / 86400.0
!       write(1000+me,*)' pwat=',pwat(i),'i=',i,',
!    &' rain=',rain(i)*1000.0,' dqsfc1=',dqsfc1(i)*tem,' kdt=',kdt
!    &,' e-p=',dqsfc1(i)*tem-rain(i)*1000.0
!     if (lprnt) write(0,*)' pwat=',pwat(ipr),',
!    &' rain=',rain(ipr)*1000.0,' dqsfc1=',dqsfc1(ipr)*tem,' kdt=',kdt
!    &,' e-p=',dqsfc1(ipr)*tem-rain(ipr)*1000.0

!
!     if (lprnt .and. rain(ipr) > 5) call mpi_quit(5678)
!     if (lat == 45) write(1000+me,*)' pwat=',pwat(1),' kdt=',kdt
!       if (lprnt) then
!         write(7000,*) ' endgu0=',gu0(ipr,:),' kdt=',kdt
!         write(7000,*) ' endgv0=',gv0(ipr,:),' kdt=',kdt,' nnp=',nnp
!         write(0,*) ' endgt0=',Stateout%gt0(ipr,:),' kdt=',kdt
!         write(0,*) ' endgq0=',Stateout%gq0(ipr,:,1),' kdt=',kdt
!         write(0,*) ' endgw0=',gq0(ipr,:,3),' kdt=',kdt,' lat=',lat
!       endif
#ifdef CCPP
      if (Model%imfdeepcnv == 3) then
        if (Model%me==0) write(0,*) 'CCPP DEBUG: calling cu_gf_driver_post through option B'
        ! Copy local variables from driver to appropriate interstitial variables
        !Interstitial(nt)%im               ! intent(in) - set in Interstitial(nt)%create
        !Stateout(nb)%gt0                  ! intent(in)
        !Stateout(nb)%gq0(:,:,1)           ! intent(in)
        !Tbd(nb)%prevst                    ! intent(out)
        !Tbd(nb)%prevsq                    ! intent(out)
        !Tbd(nb)%cactiv                    ! intent(in)
        !Sfcprop(nb)%conv_act              ! intent(out)
        cdata_block(nb,nt)%errmsg = errmsg ! intent(out)
        cdata_block(nb,nt)%errflg = errflg ! intent(out)
        call ccpp_physics_run(cdata_block(nb,nt), scheme_name="cu_gf_driver_post",ierr=ierr)
       ! Copy back intent(inout) interstitial variables to local variables in driver
        errmsg = trim(cdata_block(nb,nt)%errmsg)
        errflg = cdata_block(nb,nt)%errflg
        if (errflg/=0) then
            write(0,*) 'Error in call to cu_gf_driver_post: '//trim(errmsg)
            stop
        end if
      else if (Model%imfdeepcnv == 4) then
        if (Model%me==0) write(0,*) 'CCPP DEBUG: calling cu_ntiedtke_post through option B'
        ! Copy local variables from driver to appropriate interstitial variables
        !Interstitial(nt)%im               ! intent(in) - set in Interstitial(nt)%create
        !Stateout(nb)%gt0                  ! intent(in)
        !Stateout(nb)%gq0(:,:,1)           ! intent(in)
        !Tbd(nb)%prevst                    ! intent(out)
        !Tbd(nb)%prevsq                    ! intent(out)
        cdata_block(nb,nt)%errmsg = errmsg ! intent(out)
        cdata_block(nb,nt)%errflg = errflg ! intent(out)
        call ccpp_physics_run(cdata_block(nb,nt), scheme_name="cu_ntiedtke_post",ierr=ierr)
       ! Copy back intent(inout) interstitial variables to local variables in driver
        errmsg = trim(cdata_block(nb,nt)%errmsg)
        errflg = cdata_block(nb,nt)%errflg
        if (errflg/=0) then
            write(0,*) 'Error in call to cu_ntiedtke_post: '//trim(errmsg)
            stop
        end if
      end if
#endif

#ifndef CCPP
      if (Model%do_sppt) then
!--- radiation heating rate
        Tbd%dtdtr(1:im,:) = Tbd%dtdtr(1:im,:) + dtdtc(1:im,:)*dtf
        do i = 1, im
          if (t850(i) > 273.16) then
!--- change in change in rain precip
             Tbd%drain_cpl(i) = Diag%rain(i) - Tbd%drain_cpl(i)
          else
!--- change in change in snow precip
             Tbd%dsnow_cpl(i) = Diag%rain(i) - Tbd%dsnow_cpl(i)
          endif
        enddo
      endif
#endif

                            deallocate (clw)
      if (allocated(cnvc))  deallocate(cnvc)
      if (allocated(cnvw))  deallocate(cnvw)
      if (allocated(qrn))   deallocate(qrn)
      if (allocated(qsnw))  deallocate(qsnw)
      if (allocated(qgl))   deallocate(qgl)
      if (allocated(ncpl))  deallocate(ncpl)
      if (allocated(ncpi))  deallocate(ncpi)
      if (allocated(ncpr))  deallocate(ncpr)
      if (allocated(ncps))  deallocate(ncps)
      if (allocated(ncgl))  deallocate(ncgl)

      if (allocated(liq0))  deallocate(liq0)
      if (allocated(ice00)) deallocate(ice00)


!     deallocate (fscav, fswtr)
!
!     if (lprnt) write(0,*)' end of gbphys maxu=',
!    &maxval(gu0(1:im,1:levs)),' minu=',minval(gu0(1:im,1:levs))
!    &,' maxv=',maxval(gv0(1:im,1:levs)),' minv=',
!    & minval(gv0(1:im,1:levs)),' kdt=',kdt,' lat=',lat,' nnp=',nnp
!     if (lprnt) write(0,*)' end of gbphys gv0=',gv0(:,120:128)
!     if (lprnt) write(0,*)' end of gbphys at kdt=',kdt,
!    &' rain=',rain(ipr),' rainc=',rainc(ipr)
!     if (lprnt) call mpi_quit(7)
!     if (kdt > 2 ) call mpi_quit(70)
!    if (lprnt) write(0,*)'qt0out=',Stateout%gt0(ipr,:)    &
!    if (lprnt) write(0,*)'gq0outtke=',Stateout%gq0(ipr,1:25,ntke)    &
!      ,'xlon=',grid%xlon(ipr)*57.29578,' xlat=',grid%xlat(ipr)*57.29578
!     if (lprnt) write(0,*)' clouddriverend=',Tbd%phy_f3d(ipr,:,1)*100,' kdt=',kdt

      deallocate (qlcn, qicn, w_upi, cf_upi, CNV_MFD,           &
                  CNV_DQLDT, clcn, cnv_fice, cnv_ndrop, cnv_nice)
#ifndef CCPP
      if (imp_physics == Model%imp_physics_gfdl) then
        deallocate (delp,  dz,    uin,   vin,   pt,    qv1,   ql1, qr1,        &
                    qg1,   qa1,   qn1,   qi1,   qs1,   pt_dt, qa_dt, udt, vdt, &
                    w,     qv_dt, ql_dt, qr_dt, qi_dt, qs_dt, qg_dt,p123,refl)
        deallocate (den)
      endif
#endif

#ifdef CCPP
         if (Model%me==0) write(0,*) 'CCPP DEBUG: calling maximum_hourly_diagnostics through option B'
         ! Copy local variables from driver to appropriate interstitial variables
         !Interstitial(nt)%im               ! intent(in) - set in Interstitial(nt)%create
         !Model%levs                        ! intent(in)
         !Model%kdt                         ! intent(in)
         !Interstitial(nt)%nsteps_per_reset ! intent(in) - set in Interstitial(nt)%create
         !Model%lradar                      ! intent(in)
         !Model%imp_physics                 ! intent(in)
         !Model%imp_physics_gfdl            ! intent(in)
         !Model%imp_physics_thompson        ! intent(in)
         !con_g                             ! intent(in)
         !Statein%phil                      ! intent(in)
         !Stateout%gt0                      ! intent(in)
         !Diag%refl_10cm                    ! intent(in)
         !Diag%refdmax                      ! intent(inout)
         !Diag%refdmax263k                  ! intent(inout)
         !Diag%u10m                         ! intent(in)
         !Diag%v10m                         ! intent(in)
         !Diag%u10max                       ! intent(inout)
         !Diag%v10max                       ! intent(inout)
         !Diag%spd10max                     ! intent(inout)
         !Statein%pgr                       ! intent(in)
         !Sfcprop%t2m                       ! intent(in)
         !Sfcprop%q2m                       ! intent(in)
         !Diag%t02max                       ! intent(inout)
         !Diag%t02min                       ! intent(inout)
         !Diag%rh02max                      ! intent(inout)
         !Diag%rh02min                      ! intent(inout)
         cdata_block(nb,nt)%errmsg = errmsg ! intent(out)
         cdata_block(nb,nt)%errflg = errflg ! intent(out)
         call ccpp_physics_run(cdata_block(nb,nt), scheme_name="maximum_hourly_diagnostics",ierr=ierr)
        ! Copy back intent(inout) interstitial variables to local variables in driver
         errmsg = trim(cdata_block(nb,nt)%errmsg)
         errflg = cdata_block(nb,nt)%errflg
         if (errflg/=0) then
             write(0,*) 'Error in call to maximum_hourly_diagnostics: '//trim(errmsg)
             stop
         end if
#else
         if (Model%me==0) write(0,*) 'CCPP DEBUG: calling non-CCPP compliant version of maximum_hourly_diagnostics'
! DH* moved here from inside "if (imp_physics == Model%imp_physics_gfdl)" block
!Calculate hourly max 1-km agl and -10C reflectivity
         if(Model%lradar .and. (imp_physics == Model%imp_physics_gfdl .or. imp_physics == Model%imp_physics_thompson))then
           allocate(refd(im))
           allocate(refd263k(im))
           call max_fields(Statein%phil,Diag%refl_10cm,con_g,im,levs,refd,Stateout%gt0,refd263k)
           do i=1,im
             if(mod(kdtminus1,nsteps_per_reset)==0)then
               Diag%refdmax(I) = -35.
               Diag%refdmax263k(I) = -35.
             endif
             Diag%refdmax(i) = max(Diag%refdmax(i),refd(i))
             Diag%refdmax263k(i) = max(Diag%refdmax263k(i),refd263k(i))
           enddo
           deallocate (refd) 
           deallocate (refd263k)
         endif
! *DH
         do i=1, im
! find max hourly wind speed then decompose
            tem = sqrt(Diag%u10m(i)*Diag%u10m(i) + Diag%v10m(i)*Diag%v10m(i))
            if(mod(kdtminus1,nsteps_per_reset)==0)then
               Diag%spd10max(i) = -999.
               Diag%u10max(i)    = -999.
               Diag%v10max(i)    = -999.
               Diag%t02max(i)    = -999.
               Diag%t02min(i)    = 999.
               Diag%rh02max(i)    = -999.
               Diag%rh02min(i)    = 999.
            endif
            if (tem > Diag%spd10max(i)) then
               Diag%spd10max(i) = tem
               Diag%u10max(i)    = Diag%u10m(i)
               Diag%v10max(i)    = Diag%v10m(i)
            endif
            pshltr=Statein%pgr(i)*exp(-0.068283/Stateout%gt0(i,1))
            QCQ=PQ0/pshltr*EXP(A2A*(Sfcprop%t2m(i)-A3)/(Sfcprop%t2m(i)-A4))
            rh02=Sfcprop%q2m(i)/QCQ
            IF (rh02.GT.1.0) THEN
               rh02=1.0
            ENDIF
            IF (rh02.LT.RHmin) THEN  !use smaller RH limit for stratosphere
               rh02=RHmin
            ENDIF
            Diag%rh02max(i)=max(Diag%rh02max(i),rh02)
            Diag%rh02min(i)=min(Diag%rh02min(i),rh02)
            Diag%T02MAX(I)=MAX(Diag%T02MAX(I),Sfcprop%t2m(i))  !<--- Hourly max 2m T
            Diag%T02MIN(I)=MIN(Diag%T02MIN(I),Sfcprop%t2m(i))  !<--- Hourly min 2m T
         enddo
#endif

#ifdef CCPP
      ! End associate construct to replace local variables with CCPP variables
      end associate associate_ccpp_vars
#endif

      return
!...................................
      end subroutine GFS_physics_driver
!-----------------------------------


#ifndef CCPP
      subroutine max_fields(phil,ref3D,grav,im,levs,refd,tk,refd263k)
      use machine, only : kind_phys
      integer, intent(in)               :: im,levs
      real (kind=kind_phys), intent(in) :: grav
      real (kind=kind_phys), intent(in),dimension(im,levs)  :: phil,ref3D,tk
      integer               :: i,k,ll,ipt,kpt
      real :: dbz1avg,zmidp1,zmidloc,refl,fact
      real, dimension(im,levs) :: z
      real, dimension(im) :: zintsfc
      real, dimension(im), intent(inout) :: refd,refd263k
      REAL :: dbz1(2),dbzk,dbzk1
      logical counter
      do i=1,im
         do k=1,levs
            z(i,k)=phil(i,k)/grav
         enddo
      enddo
      do i=1,im
         refd(I) = -35.
  vloop:  do k=1,levs-1
            if ( (z(i,k+1)) .ge. 1000.     &
             .and.(z(i,k))   .le. 1000.)  then
               zmidp1=z(i,k+1)
               zmidLOC=z(i,k)
               dbz1(1)=ref3d(i,k+1)   !- dBZ (not Z) values
               dbz1(2)=ref3d(i,k) !- dBZ values
               exit vloop
            endif
         enddo vloop

!!! Initial curefl value without reduction above freezing level
!
!         curefl=0.
!         if (cprate(i,j)>0.) then
!           cuprate=rdtphs*cprate(i,j)
!           curefl=cu_a*cuprate**cu_b
!         endif
         do ll=1,2
           refl=0.
           if (dbz1(ll)>-35.) refl=10.**(0.1*dbz1(ll))
!           dbz1(l)=curefl+refl    !- in Z units
             dbz1(ll)=refl
         enddo
!-- Vertical interpolation of Z (units of mm**6/m**3)
         fact=(1000.-zmidloc)/(zmidloc-zmidp1)
         dbz1avg=dbz1(2)+(dbz1(2)-dbz1(1))*fact
!-- Convert to dBZ (10*logZ) as the last step
         if (dbz1avg>0.01) then
           dbz1avg=10.*alog10(dbz1avg)
         else
           dbz1avg=-35.
         endif
         refd(I)=max(refd(I),dbz1avg)
      enddo

!-- refl at -10C
      do i=1,im
         dbz1(1) = -35.
         dbz1(2) = -35.
  vloopm10:  do k=1,levs-1
            if (tk(i,k+1) .le. 263.15 .and. tk(i,k) .ge. 263.15)  then     
               dbz1(1)=ref3d(i,k+1)   !- dBZ (not Z) values
               dbz1(2)=ref3d(i,k) !- dBZ values
               exit vloopm10
            endif
         enddo vloopm10
         
         do ll=1,2
           refl=0.
           if (dbz1(ll)>-35.) refl=10.**(0.1*dbz1(ll))
!           dbz1(l)=curefl+refl    !- in Z units
             dbz1(ll)=refl
         enddo
!-- Take max of bounding reflectivity values 
         dbz1avg=maxval(dbz1)
!-- Convert to dBZ (10*logZ) as the last step
         if (dbz1avg>0.01) then
           dbz1avg=10.*alog10(dbz1avg)
         else
           dbz1avg=-35.
         endif
         refd263K(I)=dbz1avg
      enddo
      end subroutine max_fields
#endif

 subroutine moist_bud(im,ix,ix2,levs,me,kdt,grav,dtp,delp,rain, &
                           qv0,ql0,qi0,qv1,ql1,qi1,comp, xlon, xlat)
!  nov 2016 - S. Moorthi - routine to compute local moisture budget
      use machine, only : kind_phys
      implicit none
      character*10          :: comp
      integer               :: im,ix,ix2,levs,me,kdt
      real (kind=kind_phys) :: grav, rain(im), dtp, xlon(im), xlat(im)
      real (kind=kind_phys), dimension(ix,levs)  :: qv0,ql0,qi0,delp
      real (kind=kind_phys), dimension(ix2,levs) :: qv1,ql1,qi1
      REAL (kind=kind_phys), dimension(im) :: sumq, sumqv, sumql, sumqi
      integer               :: i, k
!
      do i=1,im
        sumqv(i) = 0.0
        sumql(i) = 0.0
        sumqi(i) = 0.0
        sumq (i) = 0.0
      enddo
      do k=1,levs
        do i=1,im
          sumqv(i) = sumqv(i) + (qv1(i,k) - qv0(i,k)) * delp(i,k)
          sumql(i) = sumql(i) + (ql1(i,k) - ql0(i,k)) * delp(i,k)
          sumqi(i) = sumqi(i) + (qi1(i,k) - qi0(i,k)) * delp(i,k)
        enddo
      enddo
      do i=1,im
        sumqv(i) = - sumqv(i) * (1.0/grav)
        sumql(i) = - sumql(i) * (1.0/grav)
        sumqi(i) = - sumqi(i) * (1.0/grav)
        sumq (i) =  sumqv(i) + sumql(i) + sumqi(i)
      enddo
      do i=1,im
        write(2000+me,*)' in moist_bud:',' i=',i,' sumq=',sumq(i), &
       ' sumqv=',sumqv(i),' sumql=',sumql(i),' sumqi=',sumqi(i),   &
       ' rain=',rain(i)*dtp,' kdt=',kdt,' component=',trim(comp),  &
       ' qv:=',qv1(i,1),qv0(i,1),' ql=',ql1(i,1),ql0(i,1),         &
       ' qi=',qi1(i,1), qi0(i,1),' xlon=',xlon(i),' xlat=',xlat(i)
      enddo
      return

      end subroutine moist_bud


      subroutine moist_bud2(im,ix,ix2,levs,me,kdt,grav,dtp,delp,rain, &
                            qv0,ql0,qi0,qr0,qs0,qg0,                  &
                            qv1,ql1,qi1,qr1,qs1,qg1,comp,xlon,xlat)
!  aug 2018 - S. Moorthi - routine to compute local moisture budget
      use machine, only : kind_phys
      implicit none
      character*10          :: comp
      integer               :: im,ix,ix2,levs,me,kdt
      real (kind=kind_phys) :: grav, rain(im), dtp, oneog, xlon(im), xlat(im)
      real (kind=kind_phys), dimension(ix,levs)  :: qv0,ql0,qi0,delp, &
                                                    qr0,qs0,qg0
      real (kind=kind_phys), dimension(ix2,levs) :: qv1,ql1,qi1,      &
                                                    qr1,qs1,qg1
      REAL (kind=kind_phys), dimension(im) :: sumq, sumqv, sumql, sumqi, &
                                              sumqr, sumqs, sumqg
      integer               :: i, k
!
      do i=1,im
        sumqv(i) = 0.0
        sumql(i) = 0.0
        sumqi(i) = 0.0
        sumqr(i) = 0.0
        sumqs(i) = 0.0
        sumqg(i) = 0.0
        sumq (i) = 0.0
      enddo
      do k=1,levs
        do i=1,im
          sumqv(i) = sumqv(i) + (qv1(i,k) - qv0(i,k)) * delp(i,k)
          sumql(i) = sumql(i) + (ql1(i,k) - ql0(i,k)) * delp(i,k)
          sumqi(i) = sumqi(i) + (qi1(i,k) - qi0(i,k)) * delp(i,k)
          sumqr(i) = sumqr(i) + (qr1(i,k) - qr0(i,k)) * delp(i,k)
          sumqs(i) = sumqs(i) + (qs1(i,k) - qs0(i,k)) * delp(i,k)
          sumqg(i) = sumqg(i) + (qg1(i,k) - qg0(i,k)) * delp(i,k)
        enddo
      enddo
      oneog = 1.0 / grav
      do i=1,im
        sumqv(i) = - sumqv(i) * oneog
        sumql(i) = - sumql(i) * oneog
        sumqi(i) = - sumqi(i) * oneog
        sumqr(i) = - sumqr(i) * oneog
        sumqs(i) = - sumqs(i) * oneog
        sumqg(i) = - sumqg(i) * oneog
        sumq (i) =  sumqv(i) + sumql(i) + sumqi(i) + sumqr(i) &
                 +  sumqs(i) + sumqg(i)
      enddo
      do i=1,im
        write(1000+me,*)' in moist_bud:',' i=',i,' sumq=',sumq(i), &
       ' sumqv=',sumqv(i),' sumql=',sumql(i),' sumqi=',sumqi(i),   &
       ' sumqr=',sumqr(i),' sumqs=',sumqs(i),' sumqg=',sumqg(i),   &
       ' rain=',rain(i)*dtp,' kdt=',kdt,' component=',trim(comp),  &
       ' qv:=',qv1(i,1),qv0(i,1),' ql=',ql1(i,1),ql0(i,1),         &
       ' qi=',qi1(i,1), qi0(i,1),' qr=',qr1(i,1),qr0(i,1),         &
       ' qs=',qs1(i,1), qs0(i,1),' qg=',qg1(i,1),qg0(i,1),         &
       ' xlon=',xlon(i),' xlat=',xlat(i)
      enddo
      return

      end subroutine moist_bud2

!> @}

end module module_physics_driver


module module_physics_driver

  use machine,               only: kind_phys
  use physcons,              only: con_cp, con_fvirt, con_g, con_rd, &
                                   con_rv, con_hvap, con_hfus,       &
                                   con_rerth, con_pi, rhc_max, dxmin,&
                                   dxinv, pa2mb, rlapse
  use cs_conv,               only: cs_convr
  use ozne_def,              only: levozp,  oz_coeff, oz_pres
  use h2o_def,               only: levh2o, h2o_coeff, h2o_pres
  use gfs_fv3_needs,         only: get_prs_fv3, get_phi_fv3
  use sasas_deep,            only: sasasdeep_run
  use sfc_nst,               only: sfc_nst_run
  use sfc_nst_pre,           only: sfc_nst_pre_run
  use sfc_nst_post,          only: sfc_nst_post_run
  use sasas_shal,            only: sasasshal_run
  use sasas_shal_post,       only: sasasshal_post_run
  use GFS_typedefs,          only: GFS_statein_type, GFS_stateout_type, &
                                   GFS_sfcprop_type, GFS_coupling_type, &
                                   GFS_control_type, GFS_grid_type,     &
                                   GFS_tbd_type,     GFS_cldprop_type,  &
                                   GFS_radtend_type, GFS_diag_type
  use edmf,                  only: edmf_run
  use GFS_PBL_generic_pre,   only: GFS_PBL_generic_pre_run
  use GFS_PBL_generic_post,  only: GFS_PBL_generic_post_run
!  use sasas_deep,             only: sasasdeep_run
  use GFS_DCNV_generic_pre,   only: GFS_DCNV_generic_pre_run
  use GFS_DCNV_generic_post,  only: GFS_DCNV_generic_post_run
  use GFS_SCNV_generic_pre,   only: GFS_SCNV_generic_pre_run
  use GFS_SCNV_generic_post,  only: GFS_SCNV_generic_post_run
  use GFS_suite_interstitial_1, only: GFS_suite_interstitial_1_run
  use GFS_suite_interstitial_2, only: GFS_suite_interstitial_2_run
  use GFS_suite_interstitial_3, only: GFS_suite_interstitial_3_run
  use GFS_suite_update_stateout, only: GFS_suite_update_stateout_run
  use GFS_suite_interstitial_4, only: GFS_suite_interstitial_4_run
  use GFS_suite_interstitial_5, only: GFS_suite_interstitial_5_run

  use GFS_zhaocarr_gscond,       only: gscond_run
  use GFS_zhaocarr_precpd,       only: precpd_run
  use GFS_calpreciptype,         only: GFS_calpreciptype_run
  use GFS_MP_generic_post,       only: GFS_MP_generic_post_run
  use GFS_MP_generic_pre,        only: GFS_MP_generic_pre_run
  use GFS_zhao_carr_pre,         only: GFS_zhao_carr_pre_run

  use GFS_surface_generic_pre,   only: GFS_surface_generic_pre_run
  use GFS_surface_generic_post,  only: GFS_surface_generic_post_run
  implicit none


  !--- CONSTANT PARAMETERS
  real(kind=kind_phys), parameter :: hocp    = con_hvap/con_cp
  real(kind=kind_phys), parameter :: qmin    = 1.0e-10
  real(kind=kind_phys), parameter :: p850    = 85000.0
  real(kind=kind_phys), parameter :: epsq    = 1.e-20
  real(kind=kind_phys), parameter :: hsub    = con_hvap+con_hfus
  real(kind=kind_phys), parameter :: czmin   = 0.0001      ! cos(89.994)
  real(kind=kind_phys), parameter :: onebg   = 1.0/con_g
  real(kind=kind_phys), parameter :: albdf   = 0.06
  real(kind=kind_phys) tf, tcr, tcrf
  parameter (tf=258.16, tcr=273.16, tcrf=1.0/(tcr-tf))


!> GFS Physics Implementation Layer
!> @brief Layer that invokes individual GFS physics routines
!> @{
!at tune step===========================================================!
!  description:                                                         !
!                                                                       !
!  usage:                                                               !
!                                                                       !
!    call gbphys                                                        !
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
!                                                                       !
!  subprograms called:                                                  !
!                                                                       !
!     get_prs,  dcyc2t2_pre_rad (testing),    dcyc2t3,  sfc_diff,       !
!     sfc_ocean,sfc_drv,  sfc_land, sfc_sice, sfc_diag, moninp1,        !
!     moninp,   moninq1,  moninq,   gwdps,    ozphys,   get_phi,        !
!     sascnv,   sascnvn,  rascnv,   cs_convr, gwdc,     shalcvt3,shalcv,!
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
!!   - determine which tracers in the tracer input array undergo convective transport (valid only for the RAS and Chikira-Sugiyama schemes) and allocate a local convective transported tracer array (clw)
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
!!   - Call the vertical diffusion scheme (PBL) based on the following logical flags: do_shoc, hybedmf, old_monin, mstrat
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
!!   - Call 'h20phys' if necessary ("adaptation of NRL H20 phys for stratosphere and mesophere")
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
!!    - if RAS, initialize 'ccwfac', 'dlqfac', 'lmh', and revap before the call to 'rascnv'
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
!!    - for the mass-flux shallow convection scheme (imfdeepcnv == 1), call 'shalcnv'
!!    - for the scale- and aerosol-aware scheme (imfshalcnv == 2), call 'mfshalcnv'
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

    subroutine GFS_physics_driver                         &
         (Model, Statein, Stateout, Sfcprop, Coupling,  &
          Grid, Tbd, Cldprop, Radtend, Diag)

      implicit none
!
!  ---  interface variables
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
!
!  ---  local variables

      !--- INTEGER VARIABLES
      integer :: me, lprint, ipr, ix, im, levs, ntrac, nvdiff, kdt
      integer :: i, kk, ic, k, n, k1, iter, levshcm, tracers,           &
                 trc_shft, tottracer, num2, num3, nshocm, nshoc, ntk

      integer, dimension(size(Grid%xlon,1)) ::                          &
           kbot, ktop, kcnv, soiltyp, vegtype, kpbl, slopetyp, kinver,  &
           lmh, levshc, islmsk,                                         &
           !--- coupling inputs for physics
           islmsk_cice

      !--- LOGICAL VARIABLES
      logical :: lprnt, revap, do_awdd

      logical, dimension(size(Grid%xlon,1)) ::                          &
           flag_iter, flag_guess, invrsn, skip_macro,                   &
           !--- coupling inputs for physics
           flag_cice

      logical, dimension(Model%ntrac-Model%ncld+2,2) ::                 &
           otspt

      !--- REAL VARIABLES
      real(kind=kind_phys) ::                                           &
           dtf, dtp, rhbbot, rhbtop, rhpbl, frain, tem, tem1, tem2,     &
           xcosz_loc, zsea1, zsea2, eng0, eng1, dpshc,                  &
           !--- experimental for shoc sub-stepping
           dtshoc

      real(kind=kind_phys), dimension(size(Grid%xlon,1))  ::            &
           ccwfac, garea, dlength, cumabs, cice, zice, tice, gflx,      &
           rain1, raincs, snowmt, cd, cdq, qss, dusfcg, dvsfcg, dusfc1, &
           dvsfc1,  dtsfc1, dqsfc1, rb, drain,  cld1d, evap, hflx,      &
           stress, t850, ep1d, gamt, gamq, sigmaf, oc, theta, gamma,    &
           sigma, elvmax, wind, work1, work2, runof, xmu, fm10, fh2,    &
           tsurf,  tx1, tx2, ctei_r, evbs, evcw, trans, sbsno, snowc,   &
           frland, adjsfcdsw, adjsfcnsw, adjsfcdlw, adjsfculw,          &
           adjnirbmu, adjnirdfu, adjvisbmu, adjvisdfu, adjnirbmd,       &
           adjnirdfd, adjvisbmd, adjvisdfd, gabsbdlw, xcosz, tseal,     &
           snohf, dlqfac, work3, ctei_rml, cldf, domr, domzr, domip,    &
           doms, psautco_l, prautco_l, ocalnirbm_cpl, ocalnirdf_cpl,    &
           ocalvisbm_cpl, ocalvisdf_cpl, dtzm, temrain1, fscav, fswtr,  &
           !--- coupling inputs for physics
           dtsfc_cice, dqsfc_cice, dusfc_cice, dvsfc_cice, ulwsfc_cice, &
           tisfc_cice, tsea_cice, hice_cice, fice_cice,                 &
           !--- for CS-convection
           wcbmax

      real(kind=kind_phys), dimension(size(Grid%xlon,1),4) ::           &
           oa4, clx

      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%lsoil) :: &
          smsoil, stsoil, slsoil

      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levs) ::  &
          del, rhc, dtdt, dudt, dvdt, gwdcu, gwdcv, dtdtc, rainp,       &
          ud_mf, dd_mf, dt_mf, prnum, dkt, sigmatot, sigmafrac

      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levs) ::  &
          initial_u, initial_v, initial_t, initial_qv

      !--- GFDL modification for FV3
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levs+1) ::&
           del_gz

      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levs,Model%ntrac) ::  &
           dqdt

      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levs,Model%nctp) ::  &
           sigmai, vverti

      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levs,oz_coeff+5) ::  &
           dq3dt_loc

      !--- ALLOCATABLE ELEMENTS
      !--- in clw, the first two varaibles are cloud water and ice.
      !--- from third to ntrac are convective transportable tracers,
      !--- third being the ozone, when ntrac=3 (valid only with ras)
      !--- Anning Cheng 9/21/2016 leave a hook here for diagnosed snow,
      !--- rain, and their number
      real(kind=kind_phys), allocatable ::                              &
           clw(:,:,:), qpl(:,:),  qpi(:,:), ncpl(:,:), ncpi(:,:),       &
           qrn(:,:), qsnw(:,:), ncpr(:,:), ncps(:,:), cnvc(:,:),        &
           cnvw(:,:)
      !--- for 2 M microphysics
      real(kind=kind_phys), allocatable, dimension(:) ::                &
             cn_prc, cn_snr
      real(kind=kind_phys), allocatable, dimension(:,:) ::              &
             qlcn, qicn, w_upi, cf_upi, CNV_MFD, CNV_PRC3, CNV_DQLDT,   &
             CLCN, CNV_FICE, CNV_NDROP, CNV_NICE
!CCPP
!      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levs) ::  &
!          initial_t, initial_qv

!
!
!===> ...  begin here

      me     = Model%me
      ix     = size(Grid%xlon,1)
      im     = size(Grid%xlon,1)
      levs   = Model%levs
      ntrac  = Model%ntrac
      dtf    = Model%dtf
      dtp    = Model%dtp
      kdt    = Model%kdt
      lprnt  = Model%lprnt
      nvdiff = ntrac           ! vertical diffusion of all tracers!
      ipr    = min(im,10)

      do i = 1, im
        if(nint(Sfcprop%slmsk(i)) == 1) then
          frland(i) = 1.0
        else
          frland(i) = 0.
        endif
      enddo
!
!  --- ...                       figure out number of extra tracers
!
      ! tottracer = 0            ! no convective transport of tracers
      ! if (Model%trans_trac .or. Model%cscnv) then
      !   if (Model%ntcw > 0) then
      !     if (Model%ntoz < Model%ntcw) then
      !       trc_shft = Model%ntcw + Model%ncld - 1
      !     else
      !       trc_shft = Model%ntoz
      !     endif
      !   elseif (Model%ntoz > 0) then
      !     trc_shft = Model%ntoz
      !   else
      !     trc_shft = 1
      !   endif
      !
      !   tracers   = Model%ntrac - trc_shft
      !   tottracer = tracers
      !   if (Model%ntoz > 0) tottracer = tottracer + 1  ! ozone is added separately
      ! endif
      ! if (Model%ntke > 0) ntk = Model%ntke - trc_shft + 3


!     if (lprnt) write(0,*)' trans_trac=',trans_trac,' tottracer=',     &
!                write(0,*)' trans_trac=',trans_trac,' tottracer=',     &
!    &                   tottracer,' trc_shft=',trc_shft,' kdt=',kdt
!    &,                  ntrac-ncld+2,' clstp=',clstp,' kdt=',kdt
!    &,' ntk=',ntk,' lat=',lat

      ! skip_macro = .false.
      !
      ! allocate ( clw(ix,levs,tottracer+2) )
      ! if (Model%imfdeepcnv >= 0 .or. Model%imfshalcnv > 0) then
      !   allocate (cnvc(ix,levs), cnvw(ix,levs))
      ! endif

      call GFS_suite_interstitial_1_run (Model, Grid, tottracer, trc_shft, tracers, ntk, skip_macro, clw, cnvc, cnvw)
!
!  ---  set initial quantities for stochastic physics deltas
      if (Model%do_sppt) then
        Tbd%dtdtr     = 0.0
        Tbd%dtotprcp  (:) = Diag%rain    (:)
        Tbd%dcnvprcp  (:) = Diag%rainc   (:)
        Tbd%drain_cpl (:) = Coupling%rain_cpl (:)
        Tbd%dsnow_cpl (:) = Coupling%snow_cpl (:)
      endif

      if (Model%do_shoc) then
        allocate (qpl(im,levs),  qpi(im,levs), ncpl(im,levs), ncpi(im,levs))
        do k=1,levs
          do i=1,im
            ncpl(i,k) = 0.0
            ncpi(i,k) = 0.0
          enddo
        enddo
      endif

      if (Model%ncld == 2) then         ! For MGB double moment microphysics
        allocate (qlcn(im,levs), qicn(im,levs), w_upi(im,levs),         &
                 cf_upi(im,levs), CNV_MFD(im,levs), CNV_PRC3(im,levs),  &
                 CNV_DQLDT(im,levs), clcn(im,levs),  cnv_fice(im,levs), &
                 cnv_ndrop(im,levs), cnv_nice(im,levs))
        allocate (cn_prc(im),   cn_snr(im))
        allocate (qrn(im,levs), qsnw(im,levs), ncpr(im,levs), ncps(im,levs))
      else
        allocate (qlcn(1,1), qicn(1,1), w_upi(1,1), cf_upi(1,1),  &
                  CNV_MFD(1,1), CNV_PRC3(1,1), CNV_DQLDT(1,1),    &
                  clcn(1,1), cnv_fice(1,1), cnv_ndrop(1,1), cnv_nice(1,1))
      endif


#ifdef GFS_HYDRO
      call get_prs(im, ix, levs, ntrac, Statein%tgrs, Statein%qgrs,     &
                   Model%thermodyn_id, Model%sfcpress_id,               &
                   Model%gen_coord_hybrid, Statein%prsi, Statein%prsik, &
                   Statein%prsl, Statein%prslk, Statein%phii, Statein%phil, del)
#else
!GFDL   Adjust the geopotential height hydrostatically in a way consistent with FV3 discretization
      call get_prs_fv3 (ix, levs, ntrac, Statein%phii, Statein%prsi, &
                        Statein%tgrs, Statein%qgrs, del, del_gz)
#endif

      call GFS_suite_interstitial_2_run (Model, Grid, Sfcprop, Statein, &
        Diag, rhbbot, rhpbl, rhbtop, frain, islmsk, work1, work2, &
        dudt, dvdt, dtdt, dtdtc, dqdt )
!
!
!  --- ...  frain=factor for centered difference scheme correction of rain amount.

      ! frain = dtf / dtp

      do i = 1, im
        !  --- ...  xw: transfer ice thickness & concentration from global to local variables
        zice(i) = Sfcprop%hice(i)
        cice(i) = Sfcprop%fice(i)
        tice(i) = Sfcprop%tisfc(i)
!
!GFDL        work1(i)   = (log(coslat(i) / (nlons(i)*latr)) - dxmin) * dxinv
!       work1(i)   = (log(Grid%dx(i)) - dxmin) * dxinv
        ! work1(i)   = (log(Grid%area(i)) - dxmin) * dxinv
        ! work1(i)   = max(0.0, min(1.0,work1(i)))
        ! work2(i)   = 1.0 - work1(i)
        ! Diag%psurf(i)   = Statein%pgr(i)
        !GFDL        tem1       = con_rerth * (con_pi+con_pi)*coslat(i)/nlons(i)
!GFDL        tem2       = con_rerth * con_pi / latr
!GFDL        garea(i)   = tem1 * tem2
        tem1       = Grid%dx(i)
        tem2       = Grid%dx(i)
        garea(i)   = Grid%area(i)
        dlength(i) = sqrt( tem1*tem1+tem2*tem2 )
        cldf(i)    = Model%cgwf(1)*work1(i) + Model%cgwf(2)*work2(i)
        wcbmax(i)  = Model%cs_parm(1)*work1(i) + Model%cs_parm(2)*work2(i)
      enddo
!
      if (Model%cplflx) then
        do i = 1, im
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
        enddo
      endif

!  --- ...  transfer soil moisture and temperature from global to local variables
      smsoil(:,:) = Sfcprop%smc(:,:)
      slsoil(:,:) = Sfcprop%slc(:,:)          !! clu: slc -> slsoil
      ! dudt(:,:)  = 0.
      ! dvdt(:,:)  = 0.
      ! dtdt(:,:)  = 0.
      ! dtdtc(:,:) = 0.
      ! dqdt(:,:,:) = 0.

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

        call dcyc2t3                                                        &
!  ---  inputs:
           ( Model%solhr, Model%slag, Model%sdec, Model%cdec, Grid%sinlat,  &
             Grid%coslat, Grid%xlon, Radtend%coszen, Sfcprop%tsfc,          &
             Statein%tgrs(1,1), Radtend%tsflw, Radtend%semis,               &
             Coupling%sfcdsw, Coupling%sfcnsw, Coupling%sfcdlw,             &
             Radtend%htrsw, Radtend%swhc, Radtend%htrlw, Radtend%lwhc,      &
             Coupling%nirbmui, Coupling%nirdfui, Coupling%visbmui,          &
             Coupling%visdfui, Coupling%nirbmdi, Coupling%nirdfdi,          &
             Coupling%visbmdi, Coupling%visdfdi, ix, im, levs,              &
!  ---  input/output:
             dtdt, dtdtc,                                                   &
!  ---  outputs:
             adjsfcdsw, adjsfcnsw, adjsfcdlw, adjsfculw, xmu, xcosz,        &
             adjnirbmu, adjnirdfu, adjvisbmu, adjvisdfu,                    &
             adjnirbmd, adjnirdfd, adjvisbmd, adjvisdfd                     &
           )

!
! save temp change due to radiation - need for sttp stochastic physics
!---------------------------------------------------------------------
      endif
!
      if (Model%lsidea) then                       !idea jw
        dtdt(:,:) = 0.
      endif

!       if (Model%lssav) then      !  --- ...  accumulate/save output variables
!
! !  --- ...  sunshine duration time is defined as the length of time (in mdl output
! !           interval) that solar radiation falling on a plane perpendicular to the
! !           direction of the sun >= 120 w/m2
!
!         do i = 1, im
!           if ( xcosz(i) >= czmin ) then   ! zenth angle > 89.994 deg
!             tem1 = adjsfcdsw(i) / xcosz(i)
!             if ( tem1 >= 120.0 ) then
!               Diag%suntim(i) = Diag%suntim(i) + dtf
!             endif
!           endif
!         enddo
!
! !  --- ...  sfc lw fluxes used by atmospheric model are saved for output
!
!         if (Model%cplflx) then
!           do i = 1, im
!             if (flag_cice(i)) adjsfculw(i) = ulwsfc_cice(i)
!           enddo
!         endif
!         Diag%dlwsfc(:) = Diag%dlwsfc(:) +   adjsfcdlw(:)*dtf
!         Diag%ulwsfc(:) = Diag%ulwsfc(:) +   adjsfculw(:)*dtf
!         Diag%psmean(:) = Diag%psmean(:) + Statein%pgr(:)*dtf        ! mean surface pressure
!
!         if (Model%ldiag3d) then
!           if (Model%lsidea) then
!             Diag%dt3dt(:,:,1) = Diag%dt3dt(:,:,1) + Radtend%lwhd(:,:,1)*dtf
!             Diag%dt3dt(:,:,2) = Diag%dt3dt(:,:,2) + Radtend%lwhd(:,:,2)*dtf
!             Diag%dt3dt(:,:,3) = Diag%dt3dt(:,:,3) + Radtend%lwhd(:,:,3)*dtf
!             Diag%dt3dt(:,:,4) = Diag%dt3dt(:,:,4) + Radtend%lwhd(:,:,4)*dtf
!             Diag%dt3dt(:,:,5) = Diag%dt3dt(:,:,5) + Radtend%lwhd(:,:,5)*dtf
!             Diag%dt3dt(:,:,6) = Diag%dt3dt(:,:,6) + Radtend%lwhd(:,:,6)*dtf
!           else
!             do k = 1, levs
!               Diag%dt3dt(:,k,1) = Diag%dt3dt(:,k,1) + Radtend%htrlw(:,k)*dtf
!               Diag%dt3dt(:,k,2) = Diag%dt3dt(:,k,2) + Radtend%htrsw(:,k)*dtf*xmu(:)
!             enddo
!           endif
!         endif
!       endif    ! end if_lssav_block
      call GFS_suite_interstitial_3_run (Model, Grid, Statein, Radtend, xcosz, &
        adjsfcdsw, adjsfcdlw, adjsfculw, xmu, Diag, kcnv, hflx, evap)
      call GFS_surface_generic_pre_run (Model, Grid, Sfcprop, Radtend, Statein,&
        adjsfcdlw, Diag sigmaf, islmsk, soiltyp, vegtype, slopetyp, work3,     &
        stsoil, gabsbdlw, tsurf, flag_guess, flag_iter, ep1d)
      call GFS_PBL_generic_pre_run (im, levs, kinver)

      !kcnv(:)   = 0
      !kinver(:) = levs
      invrsn(:) = .false.
      tx1(:)    = 0.0
      tx2(:)    = 10.0
      ctei_r(:) = 10.0

!    Only used for old shallow convection with mstrat=.true.

      if (((Model%imfshalcnv == 0 .and. Model%shal_cnv) .or. Model%old_monin)        &
                                           .and. Model%mstrat) then
        ctei_rml(:) = Model%ctei_rm(1)*work1(:) + Model%ctei_rm(2)*work2(:)
        do k = 1, levs/2
          do i = 1, im
            if (Statein%prsi(i,1)-Statein%prsi(i,k+1) < 0.35*Statein%prsi(i,1)       &
                .and. (.not. invrsn(i))) then
              tem = (Statein%tgrs(i,k+1)-Statein%tgrs(i,k)) / (Statein%prsl(i,k)-Statein%prsl(i,k+1))

              if (((tem > 0.00010) .and. (tx1(i) < 0.0)) .or.  &
                  ((tem-abs(tx1(i)) > 0.0) .and. (tx2(i) < 0.0))) then
                invrsn(i) = .true.

                if (Statein%qgrs(i,k,1) > Statein%qgrs(i,k+1,1)) then
                  tem1 = Statein%tgrs(i,k+1) + hocp*max(Statein%qgrs(i,k+1,1),qmin)
                  tem2 = Statein%tgrs(i,k)   + hocp*max(Statein%qgrs(i,k,1),qmin)

                  tem1 = tem1 / Statein%prslk(i,k+1) - tem2 / Statein%prslk(i,k)

!  --- ...  (cp/l)(deltathetae)/(deltatwater) > ctei_rm -> conditon for CTEI
                  ctei_r(i) = (1.0/hocp)*tem1/(Statein%qgrs(i,k+1,1)-Statein%qgrs(i,k,1)  &
                            + Statein%qgrs(i,k+1,Model%ntcw)-Statein%qgrs(i,k,Model%ntcw))
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

!  --- ...  lu: initialize flag_guess, flag_iter, tsurf


      drain(:)      = 0.0
      runof(:)      = 0.0
      !hflx(:)       = 0.0
      !evap(:)       = 0.0
      evbs(:)       = 0.0
      evcw(:)       = 0.0
      trans(:)      = 0.0
      sbsno(:)      = 0.0
      snowc(:)      = 0.0
      snohf(:)      = 0.0
      Diag%smcwlt2(:) = 0.0
      Diag%smcref2(:) = 0.0

!  --- ...  lu: iter-loop over (sfc_diff,sfc_drv,sfc_ocean,sfc_sice)

      do iter = 1, 2

!  --- ...  surface exchange coefficients
!
!     if (lprnt) write(0,*)' tsea=',tsea(ipr),' tsurf=',tsurf(ipr),iter

        call sfc_diff (im,Statein%pgr, Statein%ugrs, Statein%vgrs,        &
                       Statein%tgrs, Statein%qgrs, Diag%zlvl,             &
                       Sfcprop%snowd, Sfcprop%tsfc,  Sfcprop%zorl, cd,    &
                       cdq, rb, Statein%prsl(1,1), work3, islmsk, stress, &
                       Sfcprop%ffmm,  Sfcprop%ffhh, Sfcprop%uustar,       &
                       wind,  Tbd%phy_f2d(1,Model%num_p2d), fm10, fh2,    &
                       sigmaf, vegtype, Sfcprop%shdmax, Model%ivegsrc,    &
                       tsurf, flag_iter, Model%redrag)

!  --- ...  lu: update flag_guess

        do i = 1, im
          if (iter == 1 .and. wind(i) < 2.0) then
            flag_guess(i) = .true.
          endif
        enddo

        if (Model%nstf_name(1) > 0) then

          call sfc_nst_pre_run(im, islmsk, Sfcprop%oro, Sfcprop%oro_uf,        &
                               Sfcprop%tsfc,                                   &
!  --- Input/output
                               tsurf,                                          &
!  ---  outputs:
                               tseal)

          call sfc_nst_run (im, Model%lsoil, Statein%pgr, Statein%ugrs(:,1),   &
                            Statein%vgrs(:,1), Statein%tgrs(:,1),              &
                            Statein%qgrs(:,1,1), Sfcprop%tref, cd, cdq,        &
                            Statein%prsl(:,1), work3, islmsk, Grid%xlon,       &
                            Grid%sinlat, stress, Radtend%semis, gabsbdlw,      &
                            adjsfcnsw, Sfcprop%tprcp, dtf, kdt, Model%solhr,   &
                            xcosz, Tbd%phy_f2d(:,Model%num_p2d), flag_iter,    &
                            flag_guess, Model%nstf_name(1), Model%nstf_name(4),&
                            Model%nstf_name(5), lprnt, ipr,                    &
!  --- Input/output
                            tseal, tsurf, Sfcprop%xt, Sfcprop%xs, Sfcprop%xu,  &
                            Sfcprop%xv, Sfcprop%xz, Sfcprop%zm, Sfcprop%xtts,  &
                            Sfcprop%xzts, Sfcprop%dt_cool, Sfcprop%z_c,        &
                            Sfcprop%c_0, Sfcprop%c_d, Sfcprop%w_0, Sfcprop%w_d,&
                            Sfcprop%d_conv, Sfcprop%ifd, Sfcprop%qrain,        &
!  ---  outputs:
                            qss, gflx, Diag%cmm, Diag%chh, evap, hflx, ep1d)


          call sfc_nst_post_run(im, islmsk, Sfcprop%oro, Sfcprop%oro_uf,       &
                                Model%nstf_name(1), Model%nstf_name(4),        &
                                Model%nstf_name(5), Sfcprop%xt, Sfcprop%xz,    &
                                Sfcprop%dt_cool, Sfcprop%z_c, Sfcprop%slmsk,   &
                                Sfcprop%tref, Grid%xlon,                       &
!  --- Input/output
                                tsurf,                                         &
!  ---  outputs:
                                dtzm, Sfcprop%tsfc)

        else

!  --- ...  surface energy balance over ocean

          call sfc_ocean                                                &
!  ---  inputs:
           (im, Statein%pgr, Statein%ugrs, Statein%vgrs, Statein%tgrs,  &
            Statein%qgrs, Sfcprop%tsfc, cd, cdq, Statein%prsl(1,1),     &
            work3, islmsk, Tbd%phy_f2d(1,Model%num_p2d), flag_iter,     &
!  ---  outputs:
             qss, Diag%cmm, Diag%chh, gflx, evap, hflx, ep1d)

        endif       ! if ( nstf_name(1) > 0 ) then

!       if (lprnt) write(0,*)' sfalb=',sfalb(ipr),' ipr=',ipr          &
!    &,   ' weasd=',weasd(ipr),' snwdph=',snwdph(ipr)                  &
!    &,   ' tprcp=',tprcp(ipr),' kdt=',kdt,' iter=',iter               &
!    &,' tseabefland=',tsea(ipr)

!  --- ...  surface energy balance over land
!
        if (Model%lsm == 1) then                          ! noah lsm call

!     if (lprnt) write(0,*)' tsead=',tsea(ipr),' tsurf=',tsurf(ipr),iter
!    &,' pgr=',pgr(ipr),' sfcemis=',sfcemis(ipr)

          call sfc_drv                                                 &
!  ---  inputs:
           (im, Model%lsoil, Statein%pgr, Statein%ugrs, Statein%vgrs,  &
            Statein%tgrs, Statein%qgrs, soiltyp, vegtype, sigmaf,      &
            Radtend%semis, gabsbdlw, adjsfcdsw, adjsfcnsw, dtf,        &
            Sfcprop%tg3, cd, cdq, Statein%prsl(1,1), work3, DIag%zlvl, &
            islmsk, Tbd%phy_f2d(1,Model%num_p2d), slopetyp,            &
            Sfcprop%shdmin, Sfcprop%shdmax, Sfcprop%snoalb,            &
            Radtend%sfalb, flag_iter, flag_guess, Model%isot,          &
            Model%ivegsrc,                                             &
!  ---  in/outs:
            Sfcprop%weasd, Sfcprop%snowd, Sfcprop%tsfc, Sfcprop%tprcp, &
            Sfcprop%srflag, smsoil, stsoil, slsoil, Sfcprop%canopy,    &
            trans, tsurf, Sfcprop%zorl,                                &
!  ---  outputs:
            Sfcprop%sncovr, qss, gflx, drain, evap, hflx, ep1d, runof, &
            Diag%cmm, Diag%chh, evbs, evcw, sbsno, snowc, Diag%soilm,  &
            snohf, Diag%smcwlt2, Diag%smcref2, Diag%wet1)

!     if (lprnt) write(0,*)' tseae=',tsea(ipr),' tsurf=',tsurf(ipr),iter
!    &,' phy_f2d=',phy_f2d(ipr,num_p2d)

        endif

!       if (lprnt) write(0,*)' tseabeficemodel =',tsea(ipr),' me=',me   &
!    &,   ' kdt=',kdt

!  --- ...  surface energy balance over seaice

        if (Model%cplflx) then
          do i = 1, im
            if (flag_cice(i)) then
               islmsk (i) = islmsk_cice(i)
            endif
          enddo
        endif

        call sfc_sice                                                   &
!  ---  inputs:
           (im, Model%lsoil, Statein%pgr, Statein%ugrs, Statein%vgrs,   &
            Statein%tgrs, Statein%qgrs, dtf, Radtend%semis, gabsbdlw,   &
            adjsfcnsw, adjsfcdsw, Sfcprop%srflag, cd, cdq,              &
            Statein%prsl(1,1), work3, islmsk,                           &
            Tbd%phy_f2d(1,Model%num_p2d), flag_iter, Model%mom4ice,     &
            Model%lsm, lprnt, ipr,                                      &
!  ---  input/outputs:
            zice, cice, tice, Sfcprop%weasd, Sfcprop%tsfc,              &
            Sfcprop%tprcp, stsoil, ep1d,                                &
!  ---  outputs:
            Sfcprop%snowd, qss, snowmt, gflx, Diag%cmm, Diag%chh, evap, &
            hflx)

        if (Model%cplflx) then
          do i = 1, im
            if (flag_cice(i)) then
               islmsk(i) = nint(Sfcprop%slmsk(i))
            endif
          enddo

          call sfc_cice                                                 &
!  ---     inputs:
           (im, Statein%ugrs, Statein%vgrs, Statein%tgrs, Statein%qgrs, &
            cd, cdq, Statein%prsl(1,1), work3, islmsk_cice,             &
            Tbd%phy_f2d(1,Model%num_p2d),flag_iter, dqsfc_cice,         &
            dtsfc_cice,                                                 &
!  ---     outputs:
            qss, Diag%cmm, Diag%chh, evap, hflx)
        endif

!  --- ...  lu: update flag_iter and flag_guess

        do i = 1, im
          flag_iter(i)  = .false.
          flag_guess(i) = .false.

          if (iter == 1 .and. wind(i) < 2.0) then
            if ((islmsk(i) == 1) .or. ((islmsk(i) == 0) .and.           &
                                       (Model%nstf_name(1) > 0))) then
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

      enddo   ! end iter_loop


      Diag%dlwsfci(:) = adjsfcdlw(:)
      Diag%ulwsfci(:) = adjsfculw(:)
      Diag%uswsfci(:) = adjsfcdsw(:) - adjsfcnsw(:)
      Diag%dswsfci(:) = adjsfcdsw(:)

      ! Diag%t1(:)      = Statein%tgrs(:,1)
      ! Diag%q1(:)      = Statein%qgrs(:,1,1)
      ! Diag%u1(:)      = Statein%ugrs(:,1)
      ! Diag%v1(:)      = Statein%vgrs(:,1)

!  --- ...  update near surface fields

      call sfc_diag (im, Statein%pgr, Statein%ugrs, Statein%vgrs,     &
                     Statein%tgrs, Statein%qgrs, Sfcprop%tsfc, qss,   &
                     Sfcprop%f10m, Diag%u10m, Diag%v10m,        &
                     Sfcprop%t2m, Sfcprop%q2m, work3, evap,           &
                     Sfcprop%ffmm, Sfcprop%ffhh, fm10, fh2)

      Tbd%phy_f2d(:,Model%num_p2d) = 0.0

      if (Model%cplflx) then
        Coupling%dlwsfci_cpl (:) = adjsfcdlw(:)
        Coupling%dswsfci_cpl (:) = adjsfcdsw(:)
        Coupling%dlwsfc_cpl  (:) = Coupling%dlwsfc_cpl(:) + adjsfcdlw(:)*dtf
        Coupling%dswsfc_cpl  (:) = Coupling%dswsfc_cpl(:) + adjsfcdsw(:)*dtf
        Coupling%dnirbmi_cpl (:) = adjnirbmd(:)
        Coupling%dnirdfi_cpl (:) = adjnirdfd(:)
        Coupling%dvisbmi_cpl (:) = adjvisbmd(:)
        Coupling%dvisdfi_cpl (:) = adjvisdfd(:)
        Coupling%dnirbm_cpl  (:) = Coupling%dnirbm_cpl(:) + adjnirbmd(:)*dtf
        Coupling%dnirdf_cpl  (:) = Coupling%dnirdf_cpl(:) + adjnirdfd(:)*dtf
        Coupling%dvisbm_cpl  (:) = Coupling%dvisbm_cpl(:) + adjvisbmd(:)*dtf
        Coupling%dvisdf_cpl  (:) = Coupling%dvisdf_cpl(:) + adjvisdfd(:)*dtf
        Coupling%nlwsfci_cpl (:) = adjsfcdlw(:)  - adjsfculw(:)
        Coupling%nlwsfc_cpl  (:) = Coupling%nlwsfc_cpl(:) + Coupling%nlwsfci_cpl(:)*dtf
        Coupling%t2mi_cpl    (:) = Sfcprop%t2m(:)
        Coupling%q2mi_cpl    (:) = Sfcprop%q2m(:)
        Coupling%u10mi_cpl   (:) = Diag%u10m(:)
        Coupling%v10mi_cpl   (:) = Diag%v10m(:)
        Coupling%tsfci_cpl   (:) = Sfcprop%tsfc(:)
        Coupling%psurfi_cpl  (:) = Statein%pgr(:)

!  ---  estimate mean albedo for ocean point without ice cover and apply
!       them to net SW heat fluxes

        do i = 1, im
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
          Coupling%nswsfci_cpl(i) = Coupling%nnirbmi_cpl(i) + Coupling%nnirdfi_cpl(i) + &
                                    Coupling%nvisbmi_cpl(i) + Coupling%nvisdfi_cpl(i)
          Coupling%nswsfc_cpl(i)  = Coupling%nswsfc_cpl(i)  + Coupling%nswsfci_cpl(i)*dtf
          Coupling%nnirbm_cpl(i)  = Coupling%nnirbm_cpl(i)  + Coupling%nnirbmi_cpl(i)*dtf
          Coupling%nnirdf_cpl(i)  = Coupling%nnirdf_cpl(i)  + Coupling%nnirdfi_cpl(i)*dtf
          Coupling%nvisbm_cpl(i)  = Coupling%nvisbm_cpl(i)  + Coupling%nvisbmi_cpl(i)*dtf
          Coupling%nvisdf_cpl(i)  = Coupling%nvisdf_cpl(i)  + Coupling%nvisdfi_cpl(i)*dtf
        enddo
      endif

      call GFS_surface_generic_post_run (Model, Grid, ep1d, gflx, evbs, evcw, &
        trans, sbsno, snowc, snohf, stsoil, Diag, Sfcprop)

!!!!!!!!!!!!!!!!!Commented by Moorthi on July 18, 2012 !!!!!!!!!!!!!!!!!!!
!     do i = 1, im
!  --- ...  compute coefficient of evaporation in evapc
!
!       if (evapc(i) > 1.0e0) evapc(i) = 1.0e0
!  --- ...  over snow cover or ice or sea, coef of evap =1.0e0
!       if (weasd(i) > 0.0 .or. slmsk(i) /= 1.0) evapc(i) = 1.0e0
!     enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  --- ...  Boundary Layer and Free atmospheic turbulence parameterization

!     if (lprnt) write(0,*)' tsea3=',tsea(ipr),' slmsk=',slmsk(ipr)     &
!    &, ' kdt=',kdt,' evap=',evap(ipr)
!     if (lprnt)  write(0,*)' dtdtb=',(dtdt(ipr,k),k=1,15)

!     do i = 1, im
!       if (islmsk(i) == 0) then
!         oro_land(i) = 0.0
!       else
!         oro_land(i) = oro(i)
!       endif
!     enddo

!     write(0,*)' before monin clstp=',clstp,' kdt=',kdt,' lat=',lat

      if (Model%do_shoc) then
        call moninshoc(ix, im, levs, ntrac, Model%ntcw, dvdt, dudt, dtdt, dqdt,  &
                       Statein%ugrs, Statein%vgrs, Statein%tgrs, Statein%qgrs,   &
                       Tbd%phy_f3d(1,1,Model%ntot3d-1), prnum, Model%ntke,       &
                       Statein%prsik(1,1), rb, Sfcprop%zorl, Diag%u10m,          &
                       Diag%v10m, Sfcprop%ffmm, Sfcprop%ffhh, Sfcprop%tsfc, hflx,&
                       evap, stress, wind, kpbl, Statein%prsi, del, Statein%prsl,&
                       Statein%prslk, Statein%phii, Statein%phil, dtp, dusfc1,   &
                       dvsfc1, dtsfc1, dqsfc1, dkt, Diag%hpbl, kinver,           &
                       Model%xkzm_m, Model%xkzm_h, Model%xkzm_s, lprnt, ipr, me)
      else
        if (Model%hybedmf) then
          call edmf_run (ix, im, levs, nvdiff, Model%ntcw, dvdt, dudt, dtdt, dqdt,&
                         Statein%ugrs, Statein%vgrs, Statein%tgrs, Statein%qgrs,  &
                         Radtend%htrsw, Radtend%htrlw, xmu, Statein%prsik(:,1),   &
                         rb, Sfcprop%zorl, Diag%u10m, Diag%v10m, Sfcprop%ffmm,    &
                         Sfcprop%ffhh, Sfcprop%tsfc, hflx, evap, stress,     &
                         wind, kpbl, Statein%prsi, del, Statein%prsl,             &
                         Statein%prslk, Statein%phii, Statein%phil, dtp,          &
                         Model%dspheat, dusfc1, dvsfc1, dtsfc1, dqsfc1, Diag%hpbl,&
                         gamt, gamq, dkt, kinver, Model%xkzm_m, Model%xkzm_h,     &
                         Model%xkzm_s, lprnt, ipr)
!     if (lprnt)  write(0,*)' dtdtm=',(dtdt(ipr,k),k=1,15)
!     if (lprnt)  write(0,*)' dqdtm=',(dqdt(ipr,k,1),k=1,15)
        elseif (.not. Model%old_monin) then
          call moninq(ix, im, levs, nvdiff, Model%ntcw, dvdt, dudt, dtdt, dqdt, &
                      Statein%ugrs, Statein%vgrs, Statein%tgrs, Statein%qgrs,   &
                      Radtend%htrsw, Radtend%htrlw, xmu, Statein%prsik(1,1), rb,&
                      Sfcprop%ffmm, Sfcprop%ffhh, Sfcprop%tsfc, qss, hflx, evap,&
                      stress, wind, kpbl, Statein%prsi, del, Statein%prsl,      &
                      Statein%prslk, Statein%phii, Statein%phil, dtp,           &
                      Model%dspheat, dusfc1, dvsfc1, dtsfc1, dqsfc1, Diag%hpbl, &
                      gamt, gamq, dkt, kinver, Model%xkzm_m, Model%xkzm_h,      &
                      Model%xkzm_s, lprnt, ipr)
        else
          if (Model%mstrat) then
            call moninp1(ix, im, levs, nvdiff, dvdt, dudt, dtdt, dqdt,          &
                         Statein%ugrs, Statein%vgrs, Statein%tgrs, Statein%qgrs,&
                         Statein%prsik(1,1), rb, Sfcprop%ffmm, Sfcprop%ffhh,    &
                         Sfcprop%tsfc, qss, hflx, evap, stress, wind, kpbl,     &
                         Statein%prsi, del, Statein%prsl, Statein%prslk,        &
                         Statein%phii, Statein%phil, dtp, dusfc1, dvsfc1,       &
                         dtsfc1, dqsfc1, Diag%hpbl, gamt, gamq, dkt, kinver,    &
                         Model%xkzm_m, Model%xkzm_h)
          else
            call moninp(ix, im, levs, nvdiff, dvdt, dudt, dtdt, dqdt,            &
                         Statein%ugrs, Statein%vgrs, Statein%tgrs, Statein%qgrs, &
                         Statein%prsik(1,1), rb, Sfcprop%ffmm, Sfcprop%ffhh,     &
                         Sfcprop%tsfc, qss, hflx, evap, stress, wind, kpbl,      &
                         Statein%prsi, del, Statein%prsl, Statein%phii,          &
                         Statein%phil, dtp, dusfc1, dvsfc1, dtsfc1, dqsfc1,      &
                         Diag%hpbl, gamt, gamq, dkt, Model%xkzm_m, Model%xkzm_h)
          endif

        endif   ! end if_hybedmf
      endif   ! end if_do_shoc

      if (Model%cplflx) then
        do i = 1, im
          if (flag_cice(i)) then
                  cice(i)   = fice_cice(i)
            Sfcprop%tsfc(i)   = tsea_cice(i)
                  dusfc1(i) = dusfc_cice(i)
                  dvsfc1(i) = dvsfc_cice(i)
                  dqsfc1(i) = dqsfc_cice(i)
                  dtsfc1(i) = dtsfc_cice(i)
          endif
        enddo
      endif

!     if (lprnt) then
!       write(0,*) ' dusfc1=',dusfc1(ipr),' kdt=',kdt,' lat=',lat
!       write(0,*)' dtsfc1=',dtsfc1(ipr)
!       write(0,*)' dqsfc1=',dqsfc1(ipr)
!       write(0,*)' dtdtc=',(dtdt(ipr,k),k=1,15)
!       write(0,*)' dqdtc=',(dqdt(ipr,k,1),k=1,15)
!       print *,' dudtm=',dudt(ipr,:)
!     endif

!  --- ...  coupling insertion

      if (Model%cplflx) then
        Coupling%dusfc_cpl (:) = Coupling%dusfc_cpl(:) + dusfc1(:)*dtf
        Coupling%dvsfc_cpl (:) = Coupling%dvsfc_cpl(:) + dvsfc1(:)*dtf
        Coupling%dtsfc_cpl (:) = Coupling%dtsfc_cpl(:) + dtsfc1(:)*dtf
        Coupling%dqsfc_cpl (:) = Coupling%dqsfc_cpl(:) + dqsfc1(:)*dtf
        Coupling%dusfci_cpl(:) = dusfc1(:)
        Coupling%dvsfci_cpl(:) = dvsfc1(:)
        Coupling%dtsfci_cpl(:) = dtsfc1(:)
        Coupling%dqsfci_cpl(:) = dqsfc1(:)
      endif

      call GFS_PBL_generic_post_run (Grid, Model, Radtend, dusfc1, dvsfc1,   &
        dtsfc1, dqsfc1, dudt, dvdt, dtdt, dqdt, xmu, Diag)
! !-------------------------------------------------------lssav if loop ----------
!       if (Model%lssav) then
!         Diag%dusfc (:) = Diag%dusfc(:) + dusfc1(:)*dtf
!         Diag%dvsfc (:) = Diag%dvsfc(:) + dvsfc1(:)*dtf
!         Diag%dtsfc (:) = Diag%dtsfc(:) + dtsfc1(:)*dtf
!         Diag%dqsfc (:) = Diag%dqsfc(:) + dqsfc1(:)*dtf
!         Diag%dusfci(:) = dusfc1(:)
!         Diag%dvsfci(:) = dvsfc1(:)
!         Diag%dtsfci(:) = dtsfc1(:)
!         Diag%dqsfci(:) = dqsfc1(:)
! !       if (lprnt) then
! !         write(0,*)' dusfc=',dusfc(ipr),' dusfc1=',dusfc1(ipr),' dtf=',
! !    &     dtf,' kdt=',kdt,' lat=',lat
! !       endif
!
!         if (Model%ldiag3d) then
!           if (Model%lsidea) then
!             Diag%dt3dt(:,:,3) = Diag%dt3dt(:,:,3) + dtdt(:,:)*dtf
!           else
!             do k = 1, levs
!               do i = 1, im
!                 tem          = dtdt(i,k) - (Radtend%htrlw(i,k)+Radtend%htrsw(i,k)*xmu(i))
!                 Diag%dt3dt(i,k,3) = Diag%dt3dt(i,k,3) + tem*dtf
!               enddo
!             enddo
!           endif
!           Diag%du3dt(:,:,1) = Diag%du3dt(:,:,1) + dudt(:,:) * dtf
!           Diag%du3dt(:,:,2) = Diag%du3dt(:,:,2) - dudt(:,:) * dtf
!           Diag%dv3dt(:,:,1) = Diag%dv3dt(:,:,1) + dvdt(:,:) * dtf
!           Diag%dv3dt(:,:,2) = Diag%dv3dt(:,:,2) - dvdt(:,:) * dtf
! ! update dqdt_v to include moisture tendency due to vertical diffusion
! !         if (lgocart) then
! !           do k = 1, levs
! !             do i = 1, im
! !               dqdt_v(i,k)  = dqdt(i,k,1) * dtf
! !             enddo
! !           enddo
! !         endif
!           do k = 1, levs
!             do i = 1, im
!               tem  = dqdt(i,k,1) * dtf
!               Diag%dq3dt(i,k,1) = Diag%dq3dt(i,k,1) + tem
!             enddo
!           enddo
!           if (Model%ntoz > 0) then
!             Diag%dq3dt(:,:,5) = Diag%dq3dt(:,:,5) + dqdt(i,k,Model%ntoz) * dtf
!           endif
!         endif
!
!       endif   ! end if_lssav
!-------------------------------------------------------lssav if loop ----------
!
!            Orographic gravity wave drag parameterization
!            ---------------------------------------------

      if (Model%nmtvr == 14) then         ! current operational - as of 2014
        oc(:)     = Sfcprop%hprime(:,2)
        oa4(:,1)  = Sfcprop%hprime(:,3)
        oa4(:,2)  = Sfcprop%hprime(:,4)
        oa4(:,3)  = Sfcprop%hprime(:,5)
        oa4(:,4)  = Sfcprop%hprime(:,6)
        clx(:,1)  = Sfcprop%hprime(:,7)
        clx(:,2)  = Sfcprop%hprime(:,8)
        clx(:,3)  = Sfcprop%hprime(:,9)
        clx(:,4)  = Sfcprop%hprime(:,10)
        theta(:)  = Sfcprop%hprime(:,11)
        gamma(:)  = Sfcprop%hprime(:,12)
        sigma(:)  = Sfcprop%hprime(:,13)
        elvmax(:) = Sfcprop%hprime(:,14)
      elseif (Model%nmtvr == 10) then
        oc(:)     = Sfcprop%hprime(:,2)
        oa4(:,1)  = Sfcprop%hprime(:,3)
        oa4(:,2)  = Sfcprop%hprime(:,4)
        oa4(:,3)  = Sfcprop%hprime(:,5)
        oa4(:,4)  = Sfcprop%hprime(:,6)
        clx(:,1)  = Sfcprop%hprime(:,7)
        clx(:,2)  = Sfcprop%hprime(:,8)
        clx(:,3)  = Sfcprop%hprime(:,9)
        clx(:,4)  = Sfcprop%hprime(:,10)
      elseif (Model%nmtvr == 6) then
        oc(:)     = Sfcprop%hprime(:,2)
        oa4(:,1)  = Sfcprop%hprime(:,3)
        oa4(:,2)  = Sfcprop%hprime(:,4)
        oa4(:,3)  = Sfcprop%hprime(:,5)
        oa4(:,4)  = Sfcprop%hprime(:,6)
        clx(:,1)  = 0.0
        clx(:,2)  = 0.0
        clx(:,3)  = 0.0
        clx(:,4)  = 0.0
      else
        oc = 0 ; oa4 = 0 ; clx = 0 ; theta = 0 ; gamma = 0 ; sigma = 0
        elvmax = 0

      endif   ! end if_nmtvr

!     write(0,*)' before gwd clstp=',clstp,' kdt=',kdt,' lat=',lat
      call gwdps(im, ix, im, levs, dvdt, dudt, dtdt,        &
                 Statein%ugrs, Statein%vgrs, Statein%tgrs,  &
                 Statein%qgrs, kpbl, Statein%prsi, del,     &
                 Statein%prsl, Statein%prslk, Statein%phii, &
                 Statein%phil, dtp, kdt,                    &
                 Sfcprop%hprime(1,1), oc, oa4, clx, theta,  &
                 sigma, gamma, elvmax, dusfcg, dvsfcg,      &
                 con_g, con_cp, con_rd, con_rv, Model%lonr, &
                 Model%nmtvr, Model%cdmbgwd, me, lprnt,ipr)

!     if (lprnt)  print *,' dudtg=',dudt(ipr,:)

      if (Model%lssav) then
        Diag%dugwd(:) = Diag%dugwd(:) + dusfcg(:)*dtf
        Diag%dvgwd(:) = Diag%dvgwd(:) + dvsfcg(:)*dtf

!       if (lprnt) print *,' dugwd=',dugwd(ipr),' dusfcg=',dusfcg(ipr)
!       if (lprnt) print *,' dvgwd=',dvgwd(ipr),' dvsfcg=',dvsfcg(ipr)

        if (Model%ldiag3d) then
          Diag%du3dt(:,:,2) = Diag%du3dt(:,:,2) + dudt(:,:) * dtf
          Diag%dv3dt(:,:,2) = Diag%dv3dt(:,:,2) + dvdt(:,:) * dtf
          Diag%dt3dt(:,:,2) = Diag%dt3dt(:,:,2) + dtdt(:,:) * dtf
        endif
      endif

!    Rayleigh damping  near the model top
      if( .not. Model%lsidea .and. Model%ral_ts > 0.0) then
        call rayleigh_damp(im, ix, im, levs, dvdt, dudt, dtdt,      &
                           Statein%ugrs, Statein%vgrs, dtp, con_cp, &
                           Model%levr, Statein%pgr, Statein%prsl,   &
                           Model%prslrd0, Model%ral_ts)
      endif

!     if (lprnt) then
!       write(0,*)' tgrs1=',(tgrs(ipr,ik),k=1,10)
!       write(0,*)' dtdt=',(dtdt(ipr,ik),k=1,10)
!     endif

      ! Stateout%gt0(:,:)   = Statein%tgrs(:,:) + dtdt(:,:) * dtp
      ! Stateout%gu0(:,:)   = Statein%ugrs(:,:) + dudt(:,:) * dtp
      ! Stateout%gv0(:,:)   = Statein%vgrs(:,:) + dvdt(:,:) * dtp
      ! Stateout%gq0(:,:,:) = Statein%qgrs(:,:,:) + dqdt(:,:,:) * dtp

      call GFS_suite_update_stateout_run (Statein, Model, Grid, dudt, dvdt, dtdt, dqdt, Stateout)

!     if (lprnt) then
!                write(7000,*)' ugrs=',ugrs(ipr,:)
!    &,' lat=',lat,' kdt=',kdt,' me=',me
!                write(7000,*)' dudt*dtp=',dudt(ipr,:)*dtp
!                write(7000,*)' vgrs=',vgrs(ipr,:)
!                write(7000,*)' dvdt*dtp ',dvdt(ipr,:)*dtp
!     endif
!     if(lprnt) write(1000+me,*)' gq0w=',gq0(ipr,:,ntcw)
!     if(lprnt) write(0,*)' gq0i=',gq0(ipr,:,ntiw)

      if (Model%lsidea) then            ! idea convective adjustment
        call ideaca_up(Statein%prsi,Stateout%gt0,ix,im,levs+1)
      endif

!  --- ...  ozone physics

      if ((Model%ntoz > 0) .and. (ntrac >= Model%ntoz)) then
        if (oz_coeff > 4) then
          call ozphys_2015 (ix, im, levs, levozp, dtp,               &
                            Stateout%gq0(1,1,Model%ntoz),            &
                            Stateout%gq0(1,1,Model%ntoz),            &
                            Stateout%gt0, oz_pres, Statein%prsl,     &
                            Tbd%ozpl, oz_coeff, del, Model%ldiag3d,  &
                            dq3dt_loc(1,1,6), me)
          if (Model%ldiag3d) then
            Diag%dq3dt(:,:,6) = dq3dt_loc(:,:,6)
            Diag%dq3dt(:,:,7) = dq3dt_loc(:,:,7)
            Diag%dq3dt(:,:,8) = dq3dt_loc(:,:,8)
            Diag%dq3dt(:,:,9) = dq3dt_loc(:,:,9)
          endif
        else
          call ozphys (ix, im, levs, levozp, dtp,                 &
                       Stateout%gq0(1,1,Model%ntoz),              &
                       Stateout%gq0(1,1,Model%ntoz),              &
                       Stateout%gt0, oz_pres, Statein%prsl,       &
                       Tbd%ozpl, oz_coeff, del, Model%ldiag3d,    &
                       dq3dt_loc(1,1,6), me)
          if (Model%ldiag3d) then
            Diag%dq3dt(:,:,6) = dq3dt_loc(:,:,6)
            Diag%dq3dt(:,:,7) = dq3dt_loc(:,:,7)
            Diag%dq3dt(:,:,8) = dq3dt_loc(:,:,8)
            Diag%dq3dt(:,:,9) = dq3dt_loc(:,:,9)
          endif
        endif
      endif

      if (Model%h2o_phys) then
        call h2ophys (ix, im, levs, levh2o, dtp, Stateout%gq0(1,1,1),  &
                      Stateout%gq0(1,1,1), h2o_pres, Statein%prsl,     &
                      Tbd%h2opl, h2o_coeff, Model%ldiag3d,             &
                      dq3dt_loc(1,1,1), me)
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
!       write(0,*)' gt0=',gt0(ipr,:)                                    &
!    &,         ' kdt=',kdt,' xlon=',xlon(ipr),' xlat=',xlat(ipr)
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

      ! if (Model%ldiag3d) then
      !   dtdt(:,:) = Stateout%gt0(:,:)
      !   dudt(:,:) = Stateout%gu0(:,:)
      !   dvdt(:,:) = Stateout%gv0(:,:)
      ! elseif (Model%cnvgwd) then
      !   dtdt(:,:) = Stateout%gt0(:,:)
      ! endif   ! end if_ldiag3d/cnvgwd
      !
      ! if (Model%ldiag3d .or. Model%lgocart) then
      !   dqdt(:,:,1) = Stateout%gq0(:,:,1)
      ! endif   ! end if_ldiag3d/lgocart

      call GFS_DCNV_generic_pre_run (Model, Stateout, Grid, initial_u, initial_v, initial_t, initial_qv)

#ifdef GFS_HYDRO
      call get_phi(im, ix, levs, ntrac, Stateout%gt0, Stateout%gq0,    &
                   Model%thermodyn_id, Model%sfcpress_id,              &
                   Model%gen_coord_hybrid Statein%prsi, Statein%prsik, &
                   Statein%prsl, Statein%prslk, Statein%phii, Statein%phil)
#else
!GFDL   Adjust the height hydrostatically in a way consistent with FV3 discretization
      call get_phi_fv3 (ix, levs, ntrac, Stateout%gt0, Stateout%gq0, &
                        del_gz, Statein%phii, Statein%phil)
#endif

!     if (lprnt) then
!       print *,' phii2=',phii(ipr,k=1,levs)
!       print *,' phil2=',phil(ipr,:)
!     endif
      call GFS_suite_interstitial_4_run (Model, Grid, Statein, rhbbot, &
        rhbtop, work1, work2, clw, cnvc, cnvw, ktop, kbot, rhc)
      ! clw(:,:,1) = 0.0
      ! clw(:,:,2) = -999.9
      ! if ((Model%imfdeepcnv >= 0) .or. (Model%imfshalcnv > 0)) then
      !   cnvc(:,:)  = 0.0
      !   cnvw(:,:)  = 0.0
      ! endif

!     write(0,*)' before cnv clstp=',clstp,' kdt=',kdt,' lat=',lat

!  --- ...  for convective tracer transport (while using ras)

      if (Model%ras .or. Model%cscnv) then
        if (tottracer > 0) then
          if (Model%ntoz > 0) then
            clw(:,:,3) = Stateout%gq0(:,:,Model%ntoz)
            if (tracers > 0) then
              do n=1,tracers
                clw(:,:,3+n) = Stateout%gq0(:,:,n+trc_shft)
              enddo
            endif
          else
            do n=1,tracers
              clw(:,:,2+n) = Stateout%gq0(:,:,n+trc_shft)
            enddo
          endif
        endif
      endif   ! end if_ras or cfscnv

      ktop(:)  = 1
      kbot(:)  = levs

!  --- ...  calling condensation/precipitation processes
!           --------------------------------------------

      if (Model%ntcw > 0) then
        if (Model%ncld == 2) then
          clw(:,:,1) = Stateout%gq0(:,:,Model%ntiw)                    ! ice
          clw(:,:,2) = Stateout%gq0(:,:,Model%ntcw)                    ! water
        else
          if (Model%num_p3d == 4) then   ! zhao-carr microphysics
!zhang: precpd interstitial
!            psautco_l(:) = Model%psautco(1)*work1(:) + Model%psautco(2)*work2(:)
!            prautco_l(:) = Model%prautco(1)*work1(:) + Model%prautco(2)*work2(:)
!zhang: zhao_carr_pre
!            clw(:,:,1) = Stateout%gq0(:,:,Model%ntcw)
             call GFS_zhao_carr_pre_run (im,ix, levs,Stateout%gq0(:,:,Model%ntcw),clw(:,:,1))
          endif  ! end if_num_p3d
        endif    ! end if (ncld == 2)
      else    ! if_ntcw
        psautco_l(:) = Model%psautco(1)*work1(:) + Model%psautco(2)*work2(:)
        prautco_l(:) = Model%prautco(1)*work1(:) + Model%prautco(2)*work2(:)
        rhc(:,:) = 1.0
      endif   ! end if_ntcw
!
!        Call SHOC if do_shoc is true and shocaftcnv is false
!
      if (Model%do_shoc .and. .not. Model%shocaftcnv) then
        if (Model%ncld == 2) then
          skip_macro = Model%do_shoc
          clw(:,:,1) = Stateout%gq0(:,:,Model%ntiw)                    ! ice
          clw(:,:,2) = Stateout%gq0(:,:,Model%ntcw)                    ! water
          ncpl(:,:)  = Stateout%gq0(:,:,Model%ntlnc)
          ncpi(:,:)  = Stateout%gq0(:,:,Model%ntinc)
        elseif (Model%num_p3d == 4) then
          do k=1,levs
            do i=1,im
              qpl(i,k)   = 0.0
              qpi(i,k)   = 0.0
              if (abs(Stateout%gq0(i,k,Model%ntcw)) < epsq) then
                Stateout%gq0(i,k,Model%ntcw) = 0.0
              endif
              tem = Stateout%gq0(i,k,Model%ntcw)              &
     &            * max(0.0, MIN(1.0, (TCR-Stateout%gt0(i,k))*TCRF))
              clw(i,k,1) = tem                              ! ice
              clw(i,k,2) = Stateout%gq0(i,k,Model%ntcw) - tem              ! water
            enddo
          enddo
        endif

!       dtshoc = 60.0
!       dtshoc = 120.0
!       dtshoc = dtp
!       nshocm = (dtp/dtshoc) + 0.001
!       dtshoc = dtp / nshocm
!       do nshoc=1,nshocm
!      if (lprnt) write(1000+me,*)' before shoc tke=',clw(ipr,:,ntk),
!    &' kdt=',kdt,' lat=',lat,'xlon=',xlon(ipr),' xlat=',xlat(ipr)

!     phy_f3d(1,1,ntot3d-2) - shoc determined sgs clouds
!     phy_f3d(1,1,ntot3d-1) - shoc determined diffusion coefficients
!     phy_f3d(1,1,ntot3d  ) - shoc determined  w'theta'
!
!     dqdt(1:im,:,1) = gq0(1:im,:,1)
!     dqdt(1:im,:,2) = gq0(1:im,:,ntiw)
!     dqdt(1:im,:,3) = gq0(1:im,:,ntcw)
!GFDL lat has no meaning inside of shoc - changed to "1"
!GFDL          call shoc(ix, im, 1, levs, levs+1, dtp, me, lat,
          call shoc (ix, im, 1, levs, levs+1, dtp, me, 1, Statein%prsl(1,1),  &
                     Statein%phii(1,1), Statein%phil(1,1), Stateout%gu0(1,1), &
                     Stateout%gv0(1,1), Statein%vvl(1,1), Stateout%gt0(1,1),  &
                     Stateout%gq0(1,1,1), clw(1,1,1), clw(1,1,2), qpi, qpl,   &
                     rhc, Model%sup, Tbd%phy_f3d(1,1,Model%ntot3d-2),         &
                     clw(1,1,ntk), hflx, evap, prnum,                         &
                     Tbd%phy_f3d(1,1,Model%ntot3d-1),                         &
                     Tbd%phy_f3d(1,1,Model%ntot3d), lprnt, ipr, ncpl, ncpi, kdt)

!     if (lprnt) write(0,*)' aftshoccld=',phy_f3d(ipr,:,ntot3d-2)*100
!     if (lprnt) write(0,*)' aftshocice=',clw(ipr,:,1)
!     if (lprnt) write(0,*)' aftshocwat=',clw(ipr,:,1)
!     write(1000+me,*)' at latitude = ',lat
!     rain1 = 0.0
!     call moist_bud(im,im,ix,levs,me,kdt,con_g,dtp,del,rain1
!    &,                    dqdt(1,1,1), dqdt(1,1,2), dqdt(1,1,3)
!    &,                    gq0(1,1,1),clw(1,1,2),clw(1,1,1),'shoc      ')

          if ((Model%ntlnc > 0) .and. (Model%ntinc > 0) .and. (Model%ncld >= 2)) then
            Stateout%gq0(:,:,Model%ntlnc) = ncpl(:,:)
            Stateout%gq0(:,:,Model%ntinc) = ncpi(:,:)
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
      endif   ! if(do_shoc)

!  --- ...  calling convective parameterization
!
      if (.not. Model%ras .and. .not. Model%cscnv) then

        if (Model%imfdeepcnv == 1) then             ! no random cloud top
          call sascnvn (im, ix, levs, Model%jcap, dtp, del,             &
                        Statein%prsl, Statein%pgr, Statein%phil, clw(:,:,1:2),   &
                        Stateout%gq0, Stateout%gt0, Stateout%gu0,       &
                        Stateout%gv0, cld1d, rain1, kbot, ktop, kcnv,   &
                        islmsk, Statein%vvl, Model%ncld, ud_mf, dd_mf,  &
                        dt_mf, cnvw, cnvc)
        elseif (Model%imfdeepcnv == 2) then
          call sasasdeep_run (im, ix, levs, dtp, del, Statein%prsl,     &
                          Statein%pgr, Statein%phil, clw(:,:,1),        &
                          clw(:,:,2), Stateout%gq0(:,:,1),              &
                          Stateout%gt0, Stateout%gu0, Stateout%gv0,     &
                          cld1d, rain1, kbot, ktop, kcnv, islmsk,       &
                          garea, Statein%vvl, Model%ncld, ud_mf, dd_mf, &
                          dt_mf, cnvw, cnvc)
!         if (lprnt) print *,' rain1=',rain1(ipr)
        elseif (Model%imfdeepcnv == 0) then         ! random cloud top
          call sascnv (im, ix, levs, Model%jcap, dtp, del,              &
                       Statein%prsl, Statein%pgr, Statein%phil, clw(:,:,1:2),    &
                       Stateout%gq0, Stateout%gt0, Stateout%gu0,        &
                       Stateout%gv0, cld1d, rain1, kbot, ktop, kcnv,    &
                       islmsk, Statein%vvl, Tbd%rann, Model%ncld,       &
                       ud_mf, dd_mf, dt_mf, cnvw, cnvc)
!         if (lprnt) print *,' rain1=',rain1(ipr),' rann=',rann(ipr,1)
        endif
      else        ! ras or cscnv
        if (Model%cscnv) then    ! Chikira-Sugiyama  convection scheme (via CSU)
          otspt(:,:)   = .true.
          otspt(1:3,:) = .false.
          if (Model%ntke > 0) then
            otspt(Model%ntke-trc_shft+4,1)  = .false.
          endif
          if (Model%ncld == 2) then
            otspt(Model%ntlnc-trc_shft+4,1) = .false.
            otspt(Model%ntinc-trc_shft+4,1) = .false.
          endif

         fscav(:) = 0.0
         fswtr(:) = 0.0
!     write(0,*)' bef cs_cconv phii=',phii(ipr,:)
!    &,' sizefsc=',size(fscav)
!     write(0,*)' bef cs_cconv otspt=',otspt,' kdt=',kdt,' me=',me
          dqdt(:,:,1) = Stateout%gq0(:,:,1)
          dqdt(:,:,2) = max(0.0,clw(:,:,2))
          dqdt(:,:,3) = max(0.0,clw(:,:,1))
!     if (lprnt) write(0,*)' gq0bfcs=',gq0(ipr,1:35,1)
!     if (lprnt) write(0,*)' gq0bfcs3=',gq0(ipr,1:35,3)
!     if (lprnt) write(0,*)' gq0bfcs4=',gq0(ipr,1:35,4)

          do_awdd = ((Model%do_aw) .and. (Model%cs_parm(6) > 0.0))
!     if (lprnt) write(0,*)' do_awdd=',do_awdd
!GFDL  again lat replaced with "1"
!GFDL     &                  otspt, lat, kdt     ,                     &
          call cs_convr (ix, im, levs, tottracer+3, Model%nctp, otspt, 1, &
                         kdt, Stateout%gt0, Stateout%gq0(1,1,1:1), rain1, &
                         clw, Statein%phil, Statein%phii, Statein%prsl,   &
                         Statein%prsi, dtp, dtf, ud_mf, dd_mf, dt_mf,     &
                         Stateout%gu0, Stateout%gv0, fscav, fswtr,        &
                         Tbd%phy_fctd, me, wcbmax, Model%cs_parm(3),      &
                         Model%cs_parm(4), sigmai, sigmatot, vverti,      &
                         Model%do_aw, do_awdd, lprnt, ipr, QLCN, QICN,    &
                         w_upi, cf_upi, CNV_MFD, CNV_PRC3, CNV_DQLDT,     &
                         CLCN, CNV_FICE, CNV_NDROP, CNV_NICE, Model%ncld)

!     if (lprnt) write(0,*)' gq0afcs=',gq0(ipr,1:35,1)
!     if (lprnt) write(0,*)' gq0afcs3=',gq0(ipr,1:35,3)
!     if (lprnt) write(0,*)' gq0afcs4=',gq0(ipr,1:35,4)
!     write(1000+me,*)' at latitude = ',lat
!     call moist_bud(im,im,ix,levs,me,kdt,con_g,dtp,del,rain1
!    &,                    dqdt(1,1,1), dqdt(1,1,2), dqdt(1,1,3)
!    &,                    gq0(1,1,1),clw(1,1,2),clw(1,1,1),' cs_conv')

          rain1(:) = rain1(:) * (dtp*0.001)
          if (Model%do_aw) then
            do k=1,levs
              kk = min(k+1,levs)  ! assuming no cloud top reaches the model top
              do i = 1,im                                               !DD
                sigmafrac(i,k) = 0.5 * (sigmatot(i,k)+sigmatot(i,kk))
              enddo
            enddo
          endif

!     if (lprnt) then
!       write(0,*)' gt01=',gt0(ipr,:),' kdt=',kdt
!       write(0,*)' gq01=',gq0(ipr,:,1),' kdt=',kdt
!       write(0,*)' clw1=',clw(ipr,:,1),' kdt=',kdt
!       write(0,*)' clw2=',clw(ipr,:,1),' kdt=',kdt
!       write(0,*)' aft cs rain1=',rain1(ipr)*86400
!       write(0,*)' aft cs rain1=',rain1(ipr)
!     endif

        else      ! ras version 2

          if ((Model%ccwf(1) >= 0.0) .or. (Model%ccwf(2) >= 0)) then
            ccwfac(:) = Model%ccwf(1)*work1(:) + Model%ccwf(2)*work2(:)
            dlqfac(:) = Model%dlqf(1)*work1(:) + Model%dlqf(2)*work2(:)
            lmh   (:) = levs
          else
            ccwfac(:) = -999.0
            dlqfac(:) = 0.0
            lmh   (:) = levs
          endif
!         if  (lprnt) write(0,*) ' calling ras for kdt=',kdt,' me=',me    &
!    &,                       ' lprnt=',lprnt,' ccwfac=',ccwfac(ipr)

!         do k=1,levs
!           do i=1,im
!             dqdt(i,k,1) = gq0(i,k,1)
!             dqdt(i,k,2) = max(0.0,clw(i,k,2))
!             dqdt(i,k,3) = max(0.0,clw(i,k,1))
!           enddo
!         enddo
!     if (lat == 64 .and. kdt == 1) write(0,*)' qliq=',clw(1,:,1)
!     if (lat == 64 .and. kdt == 1) write(0,*)' qice=',clw(1,:,2)

          revap = .true.
!         if (ncld ==2) revap = .false.
          call rascnv (im, ix, levs, dtp, dtf, Tbd%rann, Stateout%gt0,     &
                       Stateout%gq0, Stateout%gu0, Stateout%gv0, clw,      &
                       tottracer, fscav, Statein%prsi, Statein%prsl,       &
                       Statein%prsik, Statein%prslk, Statein%phil,         &
                       Statein%phii, kpbl, cd, rain1, kbot, ktop, kcnv,    &
                       Tbd%phy_f2d(1,Model%num_p2d), Model%flipv, pa2mb,   &
                       me, garea, lmh, ccwfac, Model%nrcm, rhc, ud_mf,     &
                       dd_mf, dt_mf, dlqfac, lprnt, ipr, kdt, revap, QLCN, &
                       QICN, w_upi, cf_upi, CNV_MFD, CNV_PRC3, CNV_DQLDT,  &
                       CLCN, CNV_FICE, CNV_NDROP, CNV_NICE, Model%ncld )
        endif

!     write(1000+me,*)' at latitude = ',lat
!     tx1 = 1000.0
!     call moist_bud(im,im,ix,levs,me,kdt,con_g,tx1,del,rain1
!    &,                    dqdt(1,1,1), dqdt(1,1,2), dqdt(1,1,3)
!    &,                    gq0(1,1,1),clw(1,1,2),clw(1,1,1),' ras_conv')
!     if(lprnt) write(0,*)' after ras rain1=',rain1(ipr)
!    &,' cnv_prc3sum=',sum(cnv_prc3(ipr,1:levs))
!     if (lprnt) write(0,*)' gt04=',gt0(ipr,1:10)
!     if (lprnt) write(0,*)' gq04=',gq0(ipr,:,1)

        cld1d = 0

        if (Model%ldiag3d .or. Model%lgocart) then
          Coupling%upd_mfi(:,:) = 0.
          Coupling%dwn_mfi(:,:) = 0.
          Coupling%det_mfi(:,:) = 0.
        endif
        if (Model%lgocart) then
          Coupling%dqdti(:,:)  = 0.
          Coupling%cnvqci(:,:) = 0.
        endif

        if (Model%lgocart) then
          Coupling%upd_mfi(:,:)  = Coupling%upd_mfi(:,:)  + ud_mf(:,:) * frain
          Coupling%dwn_mfi(:,:)  = Coupling%dwn_mfi(:,:)  + dd_mf(:,:) * frain
          Coupling%det_mfi(:,:)  = Coupling%det_mfi(:,:)  + dt_mf(:,:) * frain
          Coupling%cnvqci (:,:)  = Coupling%cnvqci (:,:)  + (clw(:,:,1)+clw(:,:,2) - &
                                      Stateout%gq0(:,:,Model%ntcw)) * frain
        endif ! if (lgocart)

!  --- ...  update the tracers due to convective transport

        if (tottracer > 0) then
          if (Model%ntoz > 0) then                         ! for ozone
            Stateout%gq0(:,:,Model%ntoz) = clw(:,:,3)

            if (tracers > 0) then                    ! for other tracers
              do n=1,tracers
                Stateout%gq0(:,:,n+trc_shft) = clw(:,:,3+n)
              enddo
            endif
          else
            do n=1,tracers
              Stateout%gq0(:,:,n+trc_shft) = clw(:,:,2+n)
            enddo
          endif
        endif
      endif   ! end if_not_ras

!     if (lprnt) then
!       write(0,*)' aftcnvgt0=',gt0(ipr,:),' kdt=',kdt,' lat=',lat
!       write(0,*)' aftcnvgq0=',(gq0(ipr,k,1),k=1,levs),' lat=',lat
!       write(0,*)' gq0i2=',(gq0(ipr,k,ntiw),k=1,levs),' lat=',lat
!       write(0,*)' aftcnvgq1=',(gq0(ipr,k,ntcw),k=1,levs)
!     endif
!
!       do i = 1, im
!         Diag%rainc(:) = frain * rain1(:)
!       enddo
! !
!       if (Model%lssav) then
!         Diag%cldwrk (:) = Diag%cldwrk (:) + cld1d(:) * dtf
!         Diag%cnvprcp(:) = Diag%cnvprcp(:) + Diag%rainc(:)
!
!         if (Model%ldiag3d) then
!           Diag%dt3dt(:,:,4) = Diag%dt3dt(:,:,4) + (Stateout%gt0(:,:)-dtdt(:,:)) * frain
!           Diag%dq3dt(:,:,2) = Diag%dq3dt(:,:,2) + (Stateout%gq0(:,:,1)-dqdt(:,:,1)) * frain
!           Diag%du3dt(:,:,3) = Diag%du3dt(:,:,3) + (Stateout%gu0(:,:)-dudt(:,:)) * frain
!           Diag%dv3dt(:,:,3) = Diag%dv3dt(:,:,3) + (Stateout%gv0(:,:)-dvdt(:,:)) * frain
!
!           Diag%upd_mf(:,:)  = Diag%upd_mf(:,:)  + ud_mf(:,:) * (con_g*frain)
!           Diag%dwn_mf(:,:)  = Diag%dwn_mf(:,:)  + dd_mf(:,:) * (con_g*frain)
!           Diag%det_mf(:,:)  = Diag%det_mf(:,:)  + dt_mf(:,:) * (con_g*frain)
!         endif ! if (ldiag3d)
!
!       endif   ! end if_lssav

      call GFS_DCNV_generic_post_run (Grid, Model, Stateout, frain, rain1, cld1d, &
        initial_u, initial_v, initial_t, initial_qv, ud_mf, dd_mf, dt_mf, cnvw, cnvc, Diag, Tbd)
!
!       update dqdt_v to include moisture tendency due to deep convection
      if (Model%lgocart) then
        Coupling%dqdti  (:,:)  = (Stateout%gq0(:,:,1) - initial_qv(:,:))  * frain
        Coupling%upd_mfi(:,:)  = Coupling%upd_mfi(:,:) + ud_mf(:,:) * frain
        Coupling%dwn_mfi(:,:)  = Coupling%dwn_mfi(:,:) + dd_mf(:,:) * frain
        Coupling%det_mfi(:,:)  = Coupling%det_mfi(:,:) + dt_mf(:,:) * frain
        Coupling%cnvqci (:,:)  = Coupling%cnvqci (:,:) + (clw(:,:,1)+clw(:,:,2))*frain
      endif ! if (lgocart)
!
      ! if ((Model%npdf3d == 3) .and. (Model%num_p3d == 4)) then
      !   num2 = Model%num_p3d + 2
      !   num3 = num2 + 1
      !   Tbd%phy_f3d(:,:,num2) = cnvw(:,:)
      !   Tbd%phy_f3d(:,:,num3) = cnvc(:,:)
      ! elseif ((Model%npdf3d == 0) .and. (Model%ncnvcld3d == 1)) then
      !   num2 = Model%num_p3d + 1
      !   Tbd%phy_f3d(:,:,num2) = cnvw(:,:)
      ! endif

!     if (lprnt) write(7000,*)' bef cnvgwd gu0=',gu0(ipr,:)
!    &,' lat=',lat,' kdt=',kdt,' me=',me
!     if (lprnt) write(7000,*)' bef cnvgwd gv0=',gv0(ipr,:)
!
!----------------Convective gravity wave drag parameterization starting --------

      if (Model%cnvgwd) then         !        call convective gravity wave drag

!  --- ...  calculate maximum convective heating rate
!           cuhr = temperature change due to deep convection

        cumabs(:) = 0.0
        work3 (:)  = 0.0
        do k = 1, levs
          do i = 1, im
            if (k >= kbot(i) .and. k <= ktop(i)) then
              cumabs(i) = cumabs(i) + (Stateout%gt0(i,k)-initial_t(i,k)) * del(i,k)
              work3(i)  = work3(i)  + del(i,k)
            endif
          enddo
        enddo
        do i=1,im
          if (work3(i) > 0.0) cumabs(i) = cumabs(i) / (dtp*work3(i))
        enddo

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

!GFDL replacing lat with "1"
!       call gwdc(im, ix, im, levs, lat, gu0, gv0, gt0, gq0, dtp,       &
        call gwdc (im, ix, im, levs, 1, Statein%ugrs, Statein%vgrs,     &
                   Statein%tgrs, Statein%qgrs, dtp, Statein%prsl,       &
                   Statein%prsi, del, cumabs, ktop, kbot, kcnv, cldf,   &
                   con_g, con_cp, con_rd, con_fvirt, con_pi, dlength,   &
                   lprnt, ipr, Model%fhour, gwdcu, gwdcv, dusfcg, dvsfcg)

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

        if (Model%lssav) then
          Diag%dugwd(:) = Diag%dugwd(:) + dusfcg(:)*dtf
          Diag%dvgwd(:) = Diag%dvgwd(:) + dvsfcg(:)*dtf

          if (Model%ldiag3d) then
            Diag%du3dt(:,:,4) = Diag%du3dt(:,:,4) + gwdcu(:,:)  * dtf
            Diag%dv3dt(:,:,4) = Diag%dv3dt(:,:,4) + gwdcv(:,:)  * dtf
          endif
        endif   ! end if_lssav

!  --- ...  update the wind components with  gwdc tendencies

        do k = 1, levs
          do i = 1, im
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

      ! if (Model%ldiag3d) then
      !   initial_t(:,:)   = Stateout%gt0(:,:)
      ! endif
      ! if (Model%ldiag3d .or. Model%lgocart) then
      !   initial_qv(:,:) = Stateout%gq0(:,:,1)
      ! endif

      call GFS_SCNV_generic_pre_run (Model, Stateout, Grid, initial_t, initial_qv)

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
                         kbot, ktop, kcnv, islmsk, Statein%vvl, Model%ncld,&
                         Diag%hpbl, hflx, evap, ud_mf, dt_mf, cnvw, cnvc)

            raincs(:)     = frain * rain1(:)
            Diag%rainc(:) = Diag%rainc(:) + raincs(:)
            if (Model%lssav) then
              Diag%cnvprcp(:) = Diag%cnvprcp(:) + raincs(:)
            endif
            if ((Model%shcnvcw) .and. (Model%num_p3d == 4) .and. (Model%npdf3d == 3)) then
              Tbd%phy_f3d(:,:,num2) = cnvw(:,:)
              Tbd%phy_f3d(:,:,num3) = cnvc(:,:)
            elseif ((Model%npdf3d == 0) .and. (Model%ncnvcld3d == 1)) then
              num2 = Model%num_p3d + 1
              Tbd%phy_f3d(:,:,num2) = cnvw(:,:)
            endif

          elseif (Model%imfshalcnv == 2) then
            call sasasshal_run (im, ix, levs, dtp, del, Statein%prsl,     &
                            Statein%pgr, Statein%phil, clw(:,:,1),        &
                            clw(:,:,2), Stateout%gq0(:,:,1),              &
                            Stateout%gt0, Stateout%gu0, Stateout%gv0,     &
                            rain1, kbot, ktop, kcnv, islmsk, garea,       &
                            Statein%vvl, Model%ncld, DIag%hpbl, ud_mf,    &
                            dt_mf, cnvw, cnvc)

            call sasasshal_post_run (frain, rain1, cnvc, cnvw, Model, Grid, Diag, Tbd)

          elseif (Model%imfshalcnv == 0) then    ! modified Tiedtke Shallow convecton
                                                 !-----------------------------------
            levshc(:) = 0
            do k = 2, levs
              do i = 1, im
                dpshc = 0.3 * Statein%prsi(i,1)
                if (Statein%prsi(i,1)-Statein%prsi(i,k) <= dpshc) levshc(i) = k
              enddo
            enddo
            levshcm = 1
            do i = 1, im
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

!         if (Model%lssav) then
! !          update dqdt_v to include moisture tendency due to shallow convection
!           if (Model%lgocart) then
!             do k = 1, levs
!               do i = 1, im
!                 tem  = (Stateout%gq0(i,k,1)-initial_qv(i,k)) * frain
!                 Coupling%dqdti(i,k) = Coupling%dqdti(i,k)  + tem
!               enddo
!             enddo
!           endif
!           if (Model%ldiag3d) then
!             Diag%dt3dt(:,:,5) = Diag%dt3dt(:,:,5) + (Stateout%gt0(:,:)-initial_t(:,:)) * frain
!             Diag%dq3dt(:,:,3) = Diag%dq3dt(:,:,3) + (Stateout%gq0(:,:,1)-initial_qv(:,:)) * frain
!           endif
!         endif   ! end if_lssav
! !
!         do k = 1, levs
!           do i = 1, im
!             if (clw(i,k,2) <= -999.0) clw(i,k,2) = 0.0
!           enddo
!         enddo

        call GFS_SCNV_generic_post_run (Model, Stateout, Grid, initial_t, initial_qv, frain, Diag, clw)

!       if (lprnt) then
!         write(0,*)' prsl=',prsl(ipr,:)
!         write(0,*) ' del=',del(ipr,:)
!         write(0,*) ' befshgt0=',gt0(ipr,:)
!         write(0,*) ' befshgq0=',gq0(ipr,:,1)
!       endif

      elseif (Model%shocaftcnv) then ! if do_shoc is true and shocaftcnv is true call shoc
        if (Model%ncld == 2) then
          skip_macro = Model%do_shoc
           ncpl(:,:)  = Stateout%gq0(:,:,Model%ntlnc)
           ncpi(:,:)  = Stateout%gq0(:,:,Model%ntinc)

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
        endif
        qpl(:,:) = 0.0
        qpi(:,:) = 0.0
!       dtshoc = 60.0
!       nshocm = (dtp/dtshoc) + 0.001
!       dtshoc = dtp / nshocm
!       do nshoc=1,nshocm
!       call shoc(im, 1, levs, levs+1, dtp, me, lat,        &
!!       call shoc(im, 1, levs, levs+1, dtshoc, me, lat, &
!    &                       prsl(1:im,:), phii (1:im,:),  phil(1:im,:),&
!    &          gu0(1:im,:),gv0(1:im,:), vvl(1:im,:), gt0(1:im,:),     &
!    &                                                   gq0(1:im,:,1), &
!    &          clw(1:im,:,1), clw(1:im,:,2), qpi, qpl,  sgs_cld(1:im,:)&
!    &,         gq0(1:im,:,ntke),                                       &
!    &          phy_f3d(1:im,:,ntot3d-1), phy_f3d(1:im,:,ntot3d),       &
!    &          lprnt, ipr,                                             &
!    &          con_cp, con_g, con_hvap, con_hfus, con_hvap+con_hfus,   &
!    &          con_rv, con_rd, con_pi, con_fvirt)

!GFDL  replace lat with "1:
!       call shoc(ix, im, 1, levs, levs+1, dtshoc, me, lat,             &
        call shoc (ix, im, 1, levs, levs+1, dtp, me, 1, Statein%prsl(1,1),        &
                   Statein%phii(1,1), Statein%phil(1,1), Stateout%gu0(1,1),       &
                   Stateout%gv0(1,1), Statein%vvl(1,1), Stateout%gt0(1,1),        &
                   Stateout%gq0(1,1,1), clw(1,1,1), clw(1,1,2), qpi, qpl,rhc,     &
                   Model%sup, Tbd%phy_f3d(1,1,Model%ntot3d-2),                    &
                   Stateout%gq0(1,1,Model%ntke), hflx, evap, prnum,               &
                   Tbd%phy_f3d(1,1,Model%ntot3d-1), Tbd%phy_f3d(1,1,Model%ntot3d),&
                   lprnt, ipr, ncpl, ncpi, kdt)

        if ((Model%ntlnc > 0) .and. (Model%ntinc > 0) .and. (Model%ncld >= 2)) then
          Stateout%gq0(:,:,Model%ntlnc) = ncpl(:,:)
          Stateout%gq0(:,:,Model%ntinc) = ncpi(:,:)
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
      endif   ! if( .not. do_shoc)
!
!       if (lprnt) then
!         write(0,*)' prsl=',prsl(ipr,:)
!         write(0,*) ' del=',del(ipr,:)
!         write(0,*) ' aftshgt0=',gt0(ipr,:)
!         write(0,*) ' aftshgq0=',gq0(ipr,:,1)
!       endif

      if (Model%ntcw > 0) then

!  for microphysics
        if (Model%ncld == 2) then           ! morrison microphysics
          Stateout%gq0(:,:,Model%ntiw) = clw(:,:,1)                     ! ice
          Stateout%gq0(:,:,Model%ntcw) = clw(:,:,2)                     ! water
        elseif (Model%num_p3d == 4) then    ! if_num_p3d
! 2275 zhang: Z-C_pre
!          Stateout%gq0(:,:,Model%ntcw) = clw(:,:,1) + clw(:,:,2)
        endif   ! end if_num_p3d

      else    ! if_ntcw

        clw(:,:,1) = clw(:,:,1) + clw(:,:,2)


      endif   ! end if_ntcw

!  Legacy routine which determines convectve clouds - should be removed at some point

      call cnvc90 (Model%clstp, im, ix, Diag%rainc, kbot, ktop, levs, Statein%prsi,  &
                   Tbd%acv, Tbd%acvb, Tbd%acvt, Cldprop%cv, Cldprop%cvb, Cldprop%cvt)

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
                     Stateout%gq0(1,1,Model%ntcw), rhc, lprnt, ipr)

!       if (lprnt) then
!         print *,' rain1=',rain1(ipr),' rainc=',rainc(ipr)
!         print *,' gt0a=',gt0(ipr,:)
!         print *,' gq0a=',gq0(ipr,:,1)
!       endif
        Diag%rainc(:) = Diag%rainc(:) + frain * rain1(:)
        if(Model%lssav) then
          Diag%cnvprcp(:) = Diag%cnvprcp(:) + rain1(:) * frain

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
            Diag%dt3dt(:,:,4) = Diag%dt3dt(:,:,4) + (Stateout%gt0(:,:)  -initial_t(:,:)  ) * frain
            Diag%dq3dt(:,:,2) = Diag%dq3dt(:,:,2) + (Stateout%gq0(:,:,1)-initial_qv(:,:)) * frain
          endif
         endif
      endif               !       moist convective adjustment over
!
       call GFS_MP_generic_pre_run (im, ix,levs,clw(:,:,1),clw(:,:,2),            &
                             Model%ldiag3d, Model%ntcw, Model%ncld,               &
                             Model%num_p3d, Stateout%gt0,Stateout%gq0(:,:,1),     &
                             dtdt,dqdt(:,:,1),dqdt(:,:,3) )

! dqdt_v : instaneous moisture tendency (kg/kg/sec)
      if (Model%lgocart) then
        Coupling%dqdti(:,:) = Coupling%dqdti(:,:) / dtf
      endif
!
!     grid-scale condensation/precipitations and microphysics parameterization
!     ------------------------------------------------------------------------

      if (Model%ncld == 0) then           ! no cloud microphysics

        call lrgscl (ix, im, levs, dtp, Stateout%gt0, Stateout%gq0, &
                     Statein%prsl, del, Statein%prslk, rain1, clw)

      elseif (Model%ncld == 1) then       ! microphysics with single cloud condensate

        if (Model%num_p3d == 4) then  ! call zhao/carr/sundqvist microphysics

          if (Model%npdf3d /= 3) then               ! without pdf clouds

!           if (lprnt) then
!             write(0,*)' prsl=',prsl(ipr,:)
!             write(0,*) ' del=',del(ipr,:)
!             write(0,*) ' beflsgt0=',gt0(ipr,:),' kdt=',kdt
!             write(0,*) ' beflsgq0=',gq0(ipr,:,1),' kdt=',kdt
!             write(0,*) ' beflsgw0=',gq0(ipr,:,3),' kdt=',kdt
!           endif
                                              ! ------------------
            if (Model%do_shoc) then
              call precpd_shoc (im, ix, levs, dtp, del, Statein%prsl,              &
                                Stateout%gq0(1,1,1), Stateout%gq0(1,1,Model%ntcw), &
                                Stateout%gt0, rain1, Diag%sr, rainp, rhc,          &
                                psautco_l, prautco_l, Model%evpco, Model%wminco,   &
                                Tbd%phy_f3d(1,1,Model%ntot3d-2), lprnt, ipr)
            else
              call gscond_run (im, ix, levs, dtp, dtf, Statein%prsl, Statein%pgr,&
                           Stateout%gq0(:,:,1), clw(:,:,1), clw(:,:,2),          &
                           Stateout%gq0(:,:,Model%ntcw),                         &
                           Stateout%gt0, Tbd%phy_f3d(:,:,1), Tbd%phy_f3d(:,:,2), &
                           Tbd%phy_f2d(:,1), Tbd%phy_f3d(:,:,3),                 &
                           Tbd%phy_f3d(:,:,4), Tbd%phy_f2d(:,2), rhc,lprnt, ipr)

              call precpd_run (im, ix, levs, dtp, del, Statein%prsl,           &
                          Stateout%gq0(:,:,1), Stateout%gq0(:,:,Model%ntcw),   &
                          Stateout%gt0, rain1, Diag%sr, rainp, rhc,            &
                          Model%psautco, Model%prautco, Model%evpco,           &
                          Model%wminco, work1, lprnt, ipr)
            endif
!           if (lprnt) then
!             write(0,*)' prsl=',prsl(ipr,:)
!             write(0,*) ' del=',del(ipr,:)
!             write(0,*) ' aftlsgt0=',gt0(ipr,:),' kdt=',kdt
!             write(0,*) ' aftlsgq0=',gq0(ipr,:,1),' kdt=',kdt
!             write(0,*) ' aftlsgw0=',gq0(ipr,:,3),' kdt=',kdt
!             write(0,*)' aft precpd rain1=',rain1(1:3),' lat=',lat
!           endif
          else                                ! with pdf clouds
                                              ! ---------------
            call gscondp (im, ix, levs, dtp, dtf, Statein%prsl,        &
                          Statein%pgr, Stateout%gq0(1,1,1),            &
                          Stateout%gq0(1,1,Model%ntcw), Stateout%gt0,  &
                          Tbd%phy_f3d(1,1,1), Tbd%phy_f3d(1,1,2),      &
                          Tbd%phy_f2d(1,1), Tbd%phy_f3d(1,1,3),        &
                          Tbd%phy_f3d(1,1,4), Tbd%phy_f2d(1,2), rhc,   &
                          Tbd%phy_f3d(1,1,Model%num_p3d+1), Model%sup, &
                          lprnt, ipr, kdt)

            call precpdp (im, ix, levs,  dtp, del, Statein%prsl,       &
                          Statein%pgr, Stateout%gq0(1,1,1),            &
                          Stateout%gq0(1,1,Model%ntcw), Stateout%gt0,  &
                          rain1, Diag%sr, rainp, rhc,                  &
                          Tbd%phy_f3d(1,1,Model%num_p3d+1), psautco_l, &
                          prautco_l, Model%evpco, Model%wminco, lprnt, ipr)
          endif   ! end of grid-scale precip/microphysics options
        endif     ! end if_num_p3d

!     if (lprnt) write(0,*) ' rain1=',rain1(ipr),' rainc=',rainc(ipr),' lat=',lat

      elseif (Model%ncld == 2) then       ! MGB double-moment microphysics
!       Acheng used clw here for other code to run smoothly and minimum change
!       to make the code work. However, the nc and clw should be treated
!       in other procceses too.  August 28/2015; Hope that can be done next
!       year. I believe this will make the physical interaction more reasonable
!       Anning 12/5/2015 changed ntcw hold liquid only
        if (Model%do_shoc) then
          if (Model%fprcp == 0) then
            clw(:,:,1) = Stateout%gq0(:,:,Model%ntiw)             ! ice
            clw(:,:,2) = Stateout%gq0(:,:,Model%ntcw)             ! water
            qrn(:,:)   = 0.
            qsnw(:,:)  = 0.
            ncpr(:,:)  = 0.
            ncps(:,:)  = 0.
            Tbd%phy_f3d(:,:,1) = Tbd%phy_f3d(:,:,Model%ntot3d-2) ! clouds from shoc
          else
            clw(:,:,1) = Stateout%gq0(:,:,Model%ntiw)             ! ice
            clw(:,:,2) = Stateout%gq0(:,:,Model%ntcw)             ! water
            qrn(:,:)   = Stateout%gq0(:,:,Model%ntrw)
            qsnw(:,:)  = Stateout%gq0(:,:,Model%ntsw)
            ncpr(:,:)  = Stateout%gq0(:,:,Model%ntrnc)
            ncps(:,:)  = Stateout%gq0(:,:,Model%ntsnc)
            Tbd%phy_f3d(:,:,1) = Tbd%phy_f3d(:,:,Model%ntot3d-2) ! clouds from shoc
          end if
        elseif ((Model%imfdeepcnv >= 0) .or. (Model%imfshalcnv > 0)) then
          if (Model%fprcp == 0) then
            clw(:,:,1) = Stateout%gq0(:,:,Model%ntiw)             ! ice
            clw(:,:,2) = Stateout%gq0(:,:,Model%ntcw)             ! water
            Tbd%phy_f3d(:,:,1) = max(0.0, min(1.0,Tbd%phy_f3d(:,:,1)+cnvc(:,:)))
                                                       ! clouds from t-dt and cnvc
            qrn(:,:)   = 0.
            qsnw(:,:)  = 0.
            ncpr(:,:)  = 0.
            ncps(:,:)  = 0.
          else
            clw(:,:,1) = Stateout%gq0(:,:,Model%ntiw)             ! ice
            clw(:,:,2) = Stateout%gq0(:,:,Model%ntcw)             ! water
            Tbd%phy_f3d(:,:,1) = max(0.0, min(1.0,Tbd%phy_f3d(:,:,1)+cnvc(:,:)))
                                                       ! clouds from t-dt and cnvc
            qrn(:,:)   = Stateout%gq0(:,:,Model%ntrw)
            qsnw(:,:)  = Stateout%gq0(:,:,Model%ntsw)
            ncpr(:,:)  = Stateout%gq0(:,:,Model%ntrnc)
            ncps(:,:)  = Stateout%gq0(:,:,Model%ntsnc)
          endif
        else
                                                     ! clouds from t-dt and cnvc
          if (Model%fprcp == 0 ) then
            clw(:,:,1) = Stateout%gq0(:,:,Model%ntiw)             ! ice
            clw(:,:,2) = Stateout%gq0(:,:,Model%ntcw)             ! water
            qrn(:,:)   = 0.
            qsnw(:,:)  = 0.
            ncpr(:,:)  = 0.
            ncps(:,:)  = 0.
            Tbd%phy_f3d(:,:,1) = min(1.0, Tbd%phy_f3d(:,:,1)+cnvc(:,:))
          else
            clw(:,:,1) = Stateout%gq0(:,:,Model%ntiw)             ! ice
            clw(:,:,2) = Stateout%gq0(:,:,Model%ntcw)             ! water
            qrn(:,:)   = Stateout%gq0(:,:,Model%ntrw)
            qsnw(:,:)  = Stateout%gq0(:,:,Model%ntsw)
            ncpr(:,:)  = Stateout%gq0(:,:,Model%ntrnc)
            ncps(:,:)  = Stateout%gq0(:,:,Model%ntsnc)
            Tbd%phy_f3d(:,:,1) = min(1.0, Tbd%phy_f3d(:,:,1)+cnvc(:,:))
          endif
        endif
!       notice clw ix instead of im
!       call m_micro_driver(im,ix,levs,flipv,del,dtp,prsl,prsi,
!    &    prslk,prsik,pgr,vvl,clw(1,1,2), QLCN, clw(1,1,1),QICN,
!     if (lprnt) write(0,*)' cnv_mfdbef=',cnv_mfd(ipr,:),' flipv=',flipv
!     if(lprnt) write(0,*) ' befgq0=',gq0(ipr,:,1),' kdt=',kdt
!     if (lprnt) write(0,*)' clw1bef=',clw(ipr,:,1),' kdt=',kdt
!     if (lprnt) write(0,*)' clw2bef=',clw(ipr,:,2),' kdt=',kdt
!     if (lprnt) write(0,*)' cloudsb=',phy_f3d(ipr,:,1)*100,' kdt=',kdt
!       txa(:,:) = gq0(:,:,1)
        call m_micro_driver (im, ix, levs, Model%flipv, dtp,  Statein%prsl,      &
                             Statein%prsi, Statein%prslk, Statein%prsik,         &
                             Statein%vvl, clw(1,1,2), QLCN, clw(1,1,1), QICN,    &
                             Radtend%htrlw, Radtend%htrsw, w_upi, cf_upi,        &
                             FRLAND, Diag%HPBL, CNV_MFD, CNV_PRC3, CNV_DQLDT,    &
                             CLCN, Stateout%gu0, Stateout%gv0, Diag%dusfc,       &
                             Diag%dvsfc, dusfc1, dvsfc1, dusfc1, dvsfc1,         &
                             CNV_FICE, CNV_NDROP, CNV_NICE, Stateout%gq0(1,1,1), &
                             Stateout%gq0(1,1,Model%ntcw),                       &
                             Stateout%gq0(1,1,Model%ntiw), Stateout%gt0, rain1,  &
                             Diag%sr, Stateout%gq0(1,1,Model%ntlnc),             &
                             Stateout%gq0(1,1,Model%ntinc), Model%fprcp, qrn,    &
                             qsnw, ncpr, ncps, Tbd%phy_f3d(1,1,1), kbot,         &
                             Model%aero_in, skip_macro, cn_prc, cn_snr, lprnt,   &
                             ipr, kdt, Grid%xlat, Grid%xlon)

!     write(1000+me,*)' at latitude = ',lat
!     tx1 = 1000.0
!     call moist_bud(im,ix,ix,levs,me,kdt,con_g,tx1,del,rain1
!    &,                    txa, clw(1,1,2), clw(1,1,1)
!    &,           gq0(1,1,1),gq0(1,1,ntcw),gq0(1,1,ntcw+1),' m_micro  ')

!     if (lprnt) write(0,*) ' rain1=',rain1(ipr)*86400.0,
!    &' rainc=',rainc(ipr)*86400.0
!    &,' cn_prc=',cn_prc(ipr),' cn_snr=',cn_snr(ipr)
!     if (lprnt) write(0,*) ' aftlsgq0=',gq0(ipr,:,1),' kdt=',kdt
!     if (lprnt) write(0,*)' clw1aft=',gq0(ipr,:,ntiw),' kdt=',kdt
!     if (lprnt) write(0,*)' cloudsm=',phy_f3d(ipr,:,1)*100,' kdt=',kdt
!     if (lprnt) write(0,*)' clw2aft=',gq0(ipr,:,ntcw),' kdt=',kdt

        if (Model%fprcp == 1) then
          Stateout%gq0(:,:,Model%ntrw)  = qrn(:,:)
          Stateout%gq0(:,:,Model%ntsw)  = qsnw(:,:)
          Stateout%gq0(:,:,Model%ntrnc) = ncpr(:,:)
          Stateout%gq0(:,:,Model%ntsnc) = ncps(:,:)
        endif
      endif       ! end if_ncld
!     if (lprnt) write(0,*)' rain1 after ls=',rain1(ipr)
!
      if (Model%do_aw) then
!  Arakawa-Wu adjustment of large-scale microphysics tendencies:
!   reduce by factor of (1-sigma)
! these are microphysics increments. We want to keep (1-sigma) of the increment,
!   we will remove sigma*increment from final values
!         fsigma = 0.  ! don't apply any AW correction, in addition comment next line
!         fsigma = sigmafrac

!  adjust sfc rainrate for conservation
!  vertically integrate reduction of water increments, reduce precip by that amount

        temrain1(:) = 0.0
        do k = 1,levs
          do i = 1,im
            tem1        = sigmafrac(i,k)
            Stateout%gt0(i,k)    = Stateout%gt0(i,k) - tem1 * (Stateout%gt0(i,k)-initial_t(i,k))
            tem2        = tem1 * (Stateout%gq0(i,k,1)-initial_qv(i,k))
            Stateout%gq0(i,k,1)  = Stateout%gq0(i,k,1) - tem2
            temrain1(i) = temrain1(i) - (Statein%prsi(i,k)-Statein%prsi(i,k+1)) &
                                      * tem2 * onebg
          enddo
        enddo
        do n=Model%ntcw,Model%ntcw+Model%ncld-1
          do k = 1,levs
            do i = 1,im
              tem1        = sigmafrac(i,k) * (Stateout%gq0(i,k,n)-dqdt(i,k,n))
              Stateout%gq0(i,k,n)  = Stateout%gq0(i,k,n) - tem1
              temrain1(i) = temrain1(i) - (Statein%prsi(i,k)-Statein%prsi(i,k+1)) &
                                        * tem1 * onebg
            enddo
          enddo
        enddo
!     write(1000+me,*)' rain1=',rain1(4),' temrain1=',temrain1(i)*0.001
        rain1(:) = max(rain1(:) - temrain1(:)*0.001, 0.0_kind_phys)
      endif

!      Diag%rain(:)  = Diag%rainc(:) + frain * rain1(:)

      call GFS_calpreciptype_run (kdt, Model%nrcm, im, ix, levs, levs+1,  &
                            Tbd%rann, Model%cal_pre, Stateout%gt0,        &
                            Stateout%gq0, Statein%prsl, Statein%prsi,     &
                            Diag%rainc,frain,rain1, Statein%phii, Model%num_p3d,       &
                            Sfcprop%tsfc, Diag%sr, Tbd%phy_f3d(:,:,3),    &   ! input !zhang:Tbd%phy_f3d(:,:,3) comes from gscond_run
                            Diag%rain, domr, domzr, domip, doms, Sfcprop%srflag,     &   ! output
                            Sfcprop%tprcp)

      call GFS_MP_generic_post_run (im, ix, levs, dtf, del,               &
                         Model%lssav, Model%ldiag3d, Diag%rain,frain,     &
                         Model%ntcw, Model%ncld,                          &
                         Stateout%gq0(:,:,Model%ntcw),                    &
                         Stateout%gt0, Stateout%gq0(:,:,1),               &
                         dtdt,dqdt(:,:,1),Diag%totprcp, Diag%dt3dt(:,:,6),&
                         Diag%dq3dt(:,:,6), Diag%pwat    )

!      if (Model%cal_pre) then       ! hchuang: add dominant precipitation type algorithm
!        i = min(3,Model%num_p3d)
!        call calpreciptype_run (kdt, Model%nrcm, im, ix, levs, levs+1,        &
!                            Tbd%rann, Grid%xlat, Grid%xlon, Stateout%gt0, &
!                            Stateout%gq0, Statein%prsl, Statein%prsi,     &
!                            Diag%rain, Statein%phii, Model%num_p3d,       &
!                            Sfcprop%tsfc, Diag%sr, Tbd%phy_f3d(:,:,i),    &     ! input
!                            domr, domzr, domip, doms)                           ! output
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

!        do i=1,im
!          if(doms(i) > 0.0 .or. domip(i) > 0.0) then
!            Sfcprop%srflag(i) = 1.
!          else
!            Sfcprop%srflag(i) = 0.
!          end if
!        enddo
!      endif

!      if (Model%lssav) then
!        Diag%totprcp(:) = Diag%totprcp(:) + Diag%rain(:)

!  --- ...  estimate t850 for rain-snow decision

!      t850(:) = Stateout%gt0(:,1)

!      do k = 1, levs-1
!        do i = 1, im
!          if (Statein%prsl(i,k) > p850 .and. Statein%prsl(i,k+1) <= p850) then
!            t850(i) = Stateout%gt0(i,k) - (Statein%prsl(i,k)-p850) / &
!                      (Statein%prsl(i,k)-Statein%prsl(i,k+1)) *      &
!                      (Stateout%gt0(i,k)-Stateout%gt0(i,k+1))
!          endif
!        enddo
!      enddo

!  --- ...  lu: snow-rain detection is performed in land/sice module

!      if (Model%cal_pre) then ! hchuang: new precip type algorithm defines srflag
!        Sfcprop%tprcp(:) = max(0.0, Diag%rain(:))  ! clu: rain -> tprcp
!      else
!        do i = 1, im
!          Sfcprop%tprcp(i)  = max(0.0, Diag%rain(i) )! clu: rain -> tprcp
!          Sfcprop%srflag(i) = 0.                     ! clu: default srflag as 'rain' (i.e. 0)
!          if (t850(i) <= 273.16) then
!            Sfcprop%srflag(i) = 1.                   ! clu: set srflag to 'snow' (i.e. 1)
!          endif
!        enddo
!      endif

!  --- ...  coupling insertion

      if (Model%cplflx .or. Model%do_sppt) then
        do i = 1, im
          if (t850(i) > 273.16) then
             Coupling%rain_cpl(i) = Coupling%rain_cpl(i) + Diag%rain(i)
          else
             Coupling%snow_cpl(i) = Coupling%snow_cpl(i) + Diag%rain(i)
          endif
        enddo
      endif

!  --- ...  end coupling insertion

!!! update surface diagnosis fields at the end of phys package
!!! this change allows gocart to use filtered wind fields
!!!
      if (Model%lgocart) then
        call sfc_diag (im, Statein%pgr, Stateout%gu0, Stateout%gv0,     &
                       Stateout%gt0, Stateout%gq0, Sfcprop%tsfc, qss,   &
                       Sfcprop%f10m, Diag%u10m, Diag%v10m, Sfcprop%t2m, &
                       Sfcprop%q2m, work3, evap, Sfcprop%ffmm,          &
                       Sfcprop%ffhh, fm10, fh2)

        if (Model%lssav) then
          Diag%tmpmax (:) = max(Diag%tmpmax (:),Sfcprop%t2m(:))
          Diag%tmpmin (:) = min(Diag%tmpmin (:),Sfcprop%t2m(:))
          Diag%spfhmax(:) = max(Diag%spfhmax(:),Sfcprop%q2m(:))
          Diag%spfhmin(:) = min(Diag%spfhmin(:),Sfcprop%q2m(:))
        endif
      endif

!  --- ...  total runoff is composed of drainage into water table and
!           runoff at the surface and is accumulated in unit of meters
      if (Model%lssav) then
        tem = dtf * 0.001
        Diag%runoff(:)  = Diag%runoff(:)  + (drain(:)+runof(:)) * tem
        Diag%srunoff(:) = Diag%srunoff(:) + runof(:) * tem
      endif

!  --- ...  xw: return updated ice thickness & concentration to global array
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

!  --- ...  return updated smsoil and stsoil to global arrays
      Sfcprop%smc(:,:) = smsoil(:,:)
      Sfcprop%slc(:,:) = slsoil(:,:)

!  --- ...  calculate column precipitable water "pwat"
!      Diag%pwat(:) = 0.0
!      tem = dtf * 0.03456 / 86400.0
!      do k = 1, levs
!        work1(:) = 0.0
!        if (Model%ncld > 0) then
!          do ic = Model%ntcw, Model%ntcw+Model%ncld-1
!            work1(:) = work1(:) +  Stateout%gq0(:,k,ic)
!          enddo
!        endif
!        Diag%pwat(:) = Diag%pwat(:) + del(:,k)*(Stateout%gq0(:,k,1)+work1(:))
!     if (lprnt .and. i == ipr) write(0,*)' gq0=',
!    &gq0(i,k,1),' qgrs=',qgrs(i,k,1),' work2=',work2(i),' k=',k
!      enddo
!      Diag%pwat(:) = Diag%pwat(:) * onebg

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
!         write(0,*) ' endgt0=',gt0(ipr,:),' kdt=',kdt
!         write(0,*) ' endgq0=',gq0(ipr,:,1),' kdt=',kdt,' lat=',lat
!         write(0,*) ' endgw0=',gq0(ipr,:,3),' kdt=',kdt,' lat=',lat
!       endif

      if (Model%do_sppt) then
        !--- radiation heating rate
        Tbd%dtdtr(:,:) = Tbd%dtdtr(:,:) + dtdtc(:,:)*dtf
        !--- change in total precip
        Tbd%dtotprcp  (:) = Diag%rain  (:) - Tbd%dtotprcp(:)
        !--- change in convective precip
        Tbd%dcnvprcp  (:) = Diag%rainc (:) - Tbd%dcnvprcp(:)
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

      call GFS_suite_interstitial_5_run (clw, cnvc, cnvw)
      !deallocate (clw)
      if (Model%do_shoc) then
        deallocate (qpl, qpi, ncpl, ncpi)
      endif
      ! if (allocated(cnvc)) deallocate(cnvc)
      ! if (allocated(cnvw)) deallocate(cnvw)

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
      if (Model%ncld == 2) then         ! For MGB double moment microphysics

        deallocate (qlcn, qicn, w_upi, cf_upi, CNV_MFD, CNV_PRC3, &
                    CNV_DQLDT, clcn, cnv_fice, cnv_ndrop, cnv_nice)
        deallocate (qrn, qsnw, ncpr, ncps)
      endif

      return
!...................................
      end subroutine GFS_physics_driver
!-----------------------------------


      subroutine moist_bud(im,ix,ix2,levs,me,kdt,grav,dtp,delp,rain, &
                           qv0,ql0,qi0,qv1,ql1,qi1,comp)
!  nov 2016 - S. Moorthi - routine to compute local moisture budget
      use machine, only : kind_phys
      implicit none
      character*10          :: comp
      integer               :: im,ix,ix2,levs,me,kdt
      real (kind=kind_phys) :: grav, rain(im), dtp
      real (kind=kind_phys), dimension(ix,levs)  :: qv0,ql0,qi0,delp
      real (kind=kind_phys), dimension(ix2,levs) :: qv1,ql1,qi1
      REAL (kind=kind_phys), dimension(im) :: sumq, sumqv, sumql, sumqi
      integer               :: i, k
!
      sumqv(:) = 0.0
      sumql(:) = 0.0
      sumqi(:) = 0.0
      sumq (:) = 0.0
      do i=1,im
        sumqv(:) = sumqv(:) + (qv1(:,k) - qv0(:,k)) * delp(:,k)
        sumql(:) = sumql(:) + (ql1(:,k) - ql0(:,k)) * delp(:,k)
        sumqi(:) = sumqi(:) + (qi1(:,k) - qi0(:,k)) * delp(:,k)
      enddo
      sumqv(:) = - sumqv(:) * (1.0/grav)
      sumql(:) = - sumql(:) * (1.0/grav)
      sumqi(:) = - sumqi(:) * (1.0/grav)
      sumq (:)   =  sumqv(:) + sumql(:) + sumqi(:)
      do i=1,im
        write(1000+me,*)' in moist_bud:',' i=',i,' sumq=',sumq(i), &
       ' sumqv=',sumqv(i),' sumql=',sumql(i),' sumqi=',sumqi(i),   &
       ' rain=',rain(i)*dtp,' kdt=',kdt,' component=',comp,        &
       ' qv:=',qv1(i,1),qv0(i,1),' ql=',ql1(i,1),ql0(i,1),         &
       ' qi=',qi1(i,1), qi0(i,1)
!      if(sumq(i) > 100) then
!       write(1000+me,*)' i=',i,'  sumq=',sumq(i)
!       write(1000+me,*)' qv1=',(qv1(i,k),k=1,levs)
!       write(1000+me,*)' qv0=',(qv0(i,k),k=1,levs)
!       write(1000+me,*)' ql1=',(ql1(i,k),k=1,levs)
!       write(1000+me,*)' ql0=',(ql0(i,k),k=1,levs)
!       write(1000+me,*)' qi1=',(qi1(i,k),k=1,levs)
!       write(1000+me,*)' qi0=',(qi0(i,k),k=1,levs)
!       endif
      enddo
      return

    end subroutine moist_bud
!> @}

end module module_physics_driver

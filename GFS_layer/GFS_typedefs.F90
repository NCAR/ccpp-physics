module GFS_typedefs

       use machine,                  only: kind_phys, kind_evod
       use module_radsw_parameters,  only: topfsw_type, sfcfsw_type
       use module_radlw_parameters,  only: topflw_type, sfcflw_type
       use ozne_def,                 only: levozp, oz_coeff
       use h2o_def,                  only: levh2o, h2o_coeff
       ! Required for surface-cycling
       use sfccyc_module,             only: sfccycle_clima_type


       implicit none

       !--- parameter constants used for default initializations
       real(kind=kind_phys), parameter :: zero      = 0.0_kind_phys
       real(kind=kind_phys), parameter :: huge      = 9.9999D15
       real(kind=kind_phys), parameter :: clear_val = zero
      !real(kind=kind_phys), parameter :: clear_val = -9.9999e80
       real(kind=kind_phys), parameter :: rann_init = 0.6_kind_phys
       real(kind=kind_phys), parameter :: cn_one    = 1._kind_phys
       real(kind=kind_phys), parameter :: cn_100    = 100._kind_phys
       real(kind=kind_phys), parameter :: cn_th     = 1000._kind_phys
       real(kind=kind_phys), parameter :: cn_hr     = 3600._kind_phys


!----------------
! Data Containers
!----------------
!    !--- GFS external initialization type
!    GFS_init_type
!    !--- GFS Derived Data Types (DDTs)
!    GFS_statein_type        !< prognostic state data in from dycore
!    GFS_stateout_type       !< prognostic state or tendencies return to dycore
!    GFS_sfcprop_type        !< surface fields
!    GFS_sfccycle_type       !< surface cycling/update fields
!    GFS_coupling_type       !< fields to/from coupling with other components (e.g. land/ice/ocean/etc.)
!    !---GFS specific containers
!    GFS_control_type        !< model control parameters 
!    GFS_grid_type           !< grid and interpolation related data
!    GFS_tbd_type            !< to be determined data that doesn't fit in any one container
!    GFS_clprop_type         !< cloud fields needed by radiation from physics
!    GFS_radtend_type        !< radiation tendencies needed in physics
!    GFS_diag_type           !< fields targetted for diagnostic output

!--------------------------------------------------------------------------------
! GFS_init_type
!--------------------------------------------------------------------------------
!   This container is the minimum set of data required from the dycore/atmosphere
!   component to allow proper initialization of the GFS physics
!--------------------------------------------------------------------------------
  type GFS_init_type
    integer :: me                                !< my MPI-rank
    integer :: master                            !< master MPI-rank
    integer :: isc                               !< starting i-index for this MPI-domain
    integer :: jsc                               !< starting j-index for this MPI-domain
    integer :: nx                                !< number of points in i-dir for this MPI rank
    integer :: ny                                !< number of points in j-dir for this MPI rank
    integer :: levs                              !< number of vertical levels
    integer :: cnx                               !< number of points in i-dir for this cubed-sphere face
                                                 !< equal to gnx for lat-lon grids
    integer :: cny                               !< number of points in j-dir for this cubed-sphere face
                                                 !< equal to gny for lat-lon grids
    integer :: gnx                               !< number of global points in x-dir (i) along the equator
    integer :: gny                               !< number of global points in y-dir (j) along any meridian
    integer :: nlunit                            !< fortran unit number for file opens
    integer :: logunit                           !< fortran unit number for writing logfile
    integer :: bdat(8)                           !< model begin date in GFS format   (same as idat)
    integer :: cdat(8)                           !< model current date in GFS format (same as jdat)
    real(kind=kind_phys) :: dt_dycore            !< dynamics time step in seconds
    real(kind=kind_phys) :: dt_phys              !< physics  time step in seconds
    !--- blocking data
    integer, pointer :: blksz(:)                 !< for explicit data blocking
                                                 !< default blksz(1)=[nx*ny]
    !--- ak/bk for pressure level calculations
    real(kind=kind_phys), pointer :: ak(:)                    !< from surface (k=1) to TOA (k=levs)
    real(kind=kind_phys), pointer :: bk(:)                    !< from surface (k=1) to TOA (k=levs)
    !--- grid metrics
    real(kind=kind_phys), pointer :: xlon(:,:)   !< column longitude for MPI rank
    real(kind=kind_phys), pointer :: xlat(:,:)   !< column latitude  for MPI rank
    real(kind=kind_phys), pointer :: area(:,:)   !< column area for length scale calculations

    character(len=32), pointer :: tracer_names(:) !< tracers names to dereference tracer id
                                                  !< based on name location in array
    character(len=65) :: fn_nml                   !< namelist filename
  end type GFS_init_type


!----------------------------------------------------------------
! GFS_statein_type
!   prognostic state variables with layer and level specific data
!----------------------------------------------------------------
  type GFS_statein_type

    !--- level geopotential and pressures 
    real (kind=kind_phys), pointer :: phii  (:,:) => null()   !< interface geopotential height
    real (kind=kind_phys), pointer :: prsi  (:,:) => null()   !< model level pressure in Pa
    real (kind=kind_phys), pointer :: prsik (:,:) => null()   !< Exner function at interface

    !--- layer geopotential and pressures
    real (kind=kind_phys), pointer :: phil  (:,:) => null()   !< layer geopotential height
    real (kind=kind_phys), pointer :: prsl  (:,:) => null()   !< model layer mean pressure Pa
    real (kind=kind_phys), pointer :: prslk (:,:) => null()   !< exner function = (p/p0)**rocp

    !--- prognostic variables
    real (kind=kind_phys), pointer :: pgr  (:)     => null()  !< surface pressure (Pa) real
    real (kind=kind_phys), pointer :: ugrs (:,:)   => null()  !< u component of layer wind
    real (kind=kind_phys), pointer :: vgrs (:,:)   => null()  !< v component of layer wind
    real (kind=kind_phys), pointer :: vvl  (:,:)   => null()  !< layer mean vertical velocity in pa/sec
    real (kind=kind_phys), pointer :: tgrs (:,:)   => null()  !< model layer mean temperature in k
    real (kind=kind_phys), pointer :: qgrs (:,:,:) => null()  !< layer mean tracer concentration

    contains
      procedure :: create  => statein_create  !<   allocate array data
  end type GFS_statein_type


!------------------------------------------------------------------
! GFS_stateout_type
!   prognostic state or tendencies after physical parameterizations
!------------------------------------------------------------------
  type GFS_stateout_type

    !-- Out (physics only)
    real (kind=kind_phys), pointer :: gu0 (:,:)   => null()  !< updated zonal wind
    real (kind=kind_phys), pointer :: gv0 (:,:)   => null()  !< updated meridional wind
    real (kind=kind_phys), pointer :: gt0 (:,:)   => null()  !< updated temperature
    real (kind=kind_phys), pointer :: gq0 (:,:,:) => null()  !< updated tracers

    contains
      procedure :: create  => stateout_create  !<   allocate array data
  end type GFS_stateout_type


!---------------------------------------------------------------------------------------
! GFS_sfcprop_type
!   surface properties that may be read in and/or updated by climatology or observations
!---------------------------------------------------------------------------------------
  type GFS_sfcprop_type

    !--- In (radiation and physics)
    real (kind=kind_phys), pointer :: slmsk  (:)   => null()  !< sea/land mask array (sea:0,land:1,sea-ice:2)
    real (kind=kind_phys), pointer :: tsfc   (:)   => null()  !< surface temperature in k 
                                                              !< [tsea in gbphys.f]
    real (kind=kind_phys), pointer :: tisfc  (:)   => null()  !< surface temperature over ice fraction 
    real (kind=kind_phys), pointer :: snowd  (:)   => null()  !< snow depth water equivalent in mm ; same as snwdph
    real (kind=kind_phys), pointer :: zorl   (:)   => null()  !< surface roughness in cm 
    real (kind=kind_phys), pointer :: fice   (:)   => null()  !< ice fraction over open water grid 
    real (kind=kind_phys), pointer :: hprim  (:)   => null()  !< topographic standard deviation in m            !
    real (kind=kind_phys), pointer :: hprime (:,:) => null()  !< orographic metrics

    !--- In (radiation only)
    real (kind=kind_phys), pointer :: sncovr (:)   => null()  !< snow cover in fraction
    real (kind=kind_phys), pointer :: snoalb (:)   => null()  !< maximum snow albedo in fraction
    real (kind=kind_phys), pointer :: alvsf  (:)   => null()  !< mean vis albedo with strong cosz dependency
    real (kind=kind_phys), pointer :: alnsf  (:)   => null()  !< mean nir albedo with strong cosz dependency
    real (kind=kind_phys), pointer :: alvwf  (:)   => null()  !< mean vis albedo with weak cosz dependency
    real (kind=kind_phys), pointer :: alnwf  (:)   => null()  !< mean nir albedo with weak cosz dependency
    real (kind=kind_phys), pointer :: facsf  (:)   => null()  !< fractional coverage with strong cosz dependency
    real (kind=kind_phys), pointer :: facwf  (:)   => null()  !< fractional coverage with   weak cosz dependency

    !--- In (physics only)
    real (kind=kind_phys), pointer :: slope  (:)   => null()  !< sfc slope type for lsm
    real (kind=kind_phys), pointer :: shdmin (:)   => null()  !< min fractional coverage of green veg
    real (kind=kind_phys), pointer :: shdmax (:)   => null()  !< max fractnl cover of green veg (not used)
    real (kind=kind_phys), pointer :: tg3    (:)   => null()  !< deep soil temperature
    real (kind=kind_phys), pointer :: vfrac  (:)   => null()  !< vegetation fraction
    real (kind=kind_phys), pointer :: vtype  (:)   => null()  !< vegetation type
    real (kind=kind_phys), pointer :: stype  (:)   => null()  !< soil type
    real (kind=kind_phys), pointer :: uustar (:)   => null()  !< boundary layer parameter
    real (kind=kind_phys), pointer :: oro    (:)   => null()  !< orography 
    real (kind=kind_phys), pointer :: oro_uf (:)   => null()  !< unfiltered orography

    !-- In/Out
    real (kind=kind_phys), pointer :: hice   (:)   => null()  !< sea ice thickness
    real (kind=kind_phys), pointer :: weasd  (:)   => null()  !< water equiv of accumulated snow depth (kg/m**2)
                                                              !< over land and sea ice
    real (kind=kind_phys), pointer :: canopy (:)   => null()  !< canopy water
    real (kind=kind_phys), pointer :: ffmm   (:)   => null()  !< fm parameter from PBL scheme
    real (kind=kind_phys), pointer :: ffhh   (:)   => null()  !< fh parameter from PBL scheme
    real (kind=kind_phys), pointer :: f10m   (:)   => null()  !< fm at 10m - Ratio of sigma level 1 wind and 10m wind
    real (kind=kind_phys), pointer :: tprcp  (:)   => null()  !< sfc_fld%tprcp - total precipitation
    real (kind=kind_phys), pointer :: srflag (:)   => null()  !< sfc_fld%srflag - snow/rain flag for precipitation
    real (kind=kind_phys), pointer :: slc    (:,:) => null()  !< liquid soil moisture
    real (kind=kind_phys), pointer :: smc    (:,:) => null()  !< total soil moisture
    real (kind=kind_phys), pointer :: stc    (:,:) => null()  !< soil temperature

    !--- Out
    real (kind=kind_phys), pointer :: t2m    (:)   => null()  !< 2 meter temperature
    real (kind=kind_phys), pointer :: q2m    (:)   => null()  !< 2 meter humidity

    !--- NSSTM variables  (only allocated when [Model%nstf_name(1) > 0])
    real (kind=kind_phys), pointer :: tref   (:)   => null()  !< nst_fld%Tref - Reference Temperature
    real (kind=kind_phys), pointer :: z_c    (:)   => null()  !< nst_fld%z_c - Sub layer cooling thickness
    real (kind=kind_phys), pointer :: c_0    (:)   => null()  !< nst_fld%c_0 - coefficient1 to calculate d(Tz)/d(Ts)
    real (kind=kind_phys), pointer :: c_d    (:)   => null()  !< nst_fld%c_d - coefficient2 to calculate d(Tz)/d(Ts)
    real (kind=kind_phys), pointer :: w_0    (:)   => null()  !< nst_fld%w_0 - coefficient3 to calculate d(Tz)/d(Ts)
    real (kind=kind_phys), pointer :: w_d    (:)   => null()  !< nst_fld%w_d - coefficient4 to calculate d(Tz)/d(Ts)
    real (kind=kind_phys), pointer :: xt     (:)   => null()  !< nst_fld%xt      heat content in DTL
    real (kind=kind_phys), pointer :: xs     (:)   => null()  !< nst_fld%xs      salinity  content in DTL
    real (kind=kind_phys), pointer :: xu     (:)   => null()  !< nst_fld%xu      u current content in DTL
    real (kind=kind_phys), pointer :: xv     (:)   => null()  !< nst_fld%xv      v current content in DTL
    real (kind=kind_phys), pointer :: xz     (:)   => null()  !< nst_fld%xz      DTL thickness
    real (kind=kind_phys), pointer :: zm     (:)   => null()  !< nst_fld%zm      MXL thickness
    real (kind=kind_phys), pointer :: xtts   (:)   => null()  !< nst_fld%xtts    d(xt)/d(ts)
    real (kind=kind_phys), pointer :: xzts   (:)   => null()  !< nst_fld%xzts    d(xz)/d(ts)
    real (kind=kind_phys), pointer :: d_conv (:)   => null()  !< nst_fld%d_conv  thickness of Free Convection Layer (FCL)
    real (kind=kind_phys), pointer :: ifd    (:)   => null()  !< nst_fld%ifd     index to start DTM run or not
    real (kind=kind_phys), pointer :: dt_cool(:)   => null()  !< nst_fld%dt_cool Sub layer cooling amount
    real (kind=kind_phys), pointer :: qrain  (:)   => null()  !< nst_fld%qrain   sensible heat flux due to rainfall (watts)

    contains
      procedure :: create  => sfcprop_create  !<   allocate array data
  end type GFS_sfcprop_type


  !---------------------------------------------------------------------
  ! GFS_sfccycle_type
  !   fields required for sfccyle routine called through gcycle
  !---------------------------------------------------------------------
  type GFS_sfccycle_type

    integer                             :: ifp                     !<
    real (kind=kind_phys), pointer      :: glacir (:)   => null()  !<
    real (kind=kind_phys), pointer      :: amxice (:)   => null()  !<
    real (kind=kind_phys), pointer      :: tsfcl0 (:)   => null()  !<
    type (sfccycle_clima_type)          :: clima                   !<

    contains
      procedure :: create  => sfccycle_create  !<   allocate array data
  end type GFS_sfccycle_type


!---------------------------------------------------------------------
! GFS_coupling_type
!   fields to/from other coupled components (e.g. land/ice/ocean/etc.)
!---------------------------------------------------------------------
  type GFS_coupling_type

    !--- Out (radiation only)
    real (kind=kind_phys), pointer :: nirbmdi(:)     => null()   !< sfc nir beam sw downward flux (w/m2) 
    real (kind=kind_phys), pointer :: nirdfdi(:)     => null()   !< sfc nir diff sw downward flux (w/m2)
    real (kind=kind_phys), pointer :: visbmdi(:)     => null()   !< sfc uv+vis beam sw downward flux (w/m2)
    real (kind=kind_phys), pointer :: visdfdi(:)     => null()   !< sfc uv+vis diff sw downward flux (w/m2)
    real (kind=kind_phys), pointer :: nirbmui(:)     => null()   !< sfc nir beam sw upward flux (w/m2)
    real (kind=kind_phys), pointer :: nirdfui(:)     => null()   !< sfc nir diff sw upward flux (w/m2)
    real (kind=kind_phys), pointer :: visbmui(:)     => null()   !< sfc uv+vis beam sw upward flux (w/m2)
    real (kind=kind_phys), pointer :: visdfui(:)     => null()   !< sfc uv+vis diff sw upward flux (w/m2)

    !--- In (physics only)
    real (kind=kind_phys), pointer :: sfcdsw(:)      => null()   !< total sky sfc downward sw flux ( w/m**2 )
                                                                 !< GFS_radtend_type%sfcfsw%dnfxc
    real (kind=kind_phys), pointer :: sfcnsw(:)      => null()   !< total sky sfc netsw flx into ground(w/m**2)
                                                                 !< difference of dnfxc & upfxc from GFS_radtend_type%sfcfsw
    real (kind=kind_phys), pointer :: sfcdlw(:)      => null()   !< total sky sfc downward lw flux ( w/m**2 )
                                                                 !< GFS_radtend_type%sfclsw%dnfxc

    !--- incoming quantities
    real (kind=kind_phys), pointer :: dusfcin_cpl(:) => null()   !< aoi_fld%dusfcin(item,lan)
    real (kind=kind_phys), pointer :: dvsfcin_cpl(:) => null()   !< aoi_fld%dvsfcin(item,lan)
    real (kind=kind_phys), pointer :: dtsfcin_cpl(:) => null()   !< aoi_fld%dtsfcin(item,lan)
    real (kind=kind_phys), pointer :: dqsfcin_cpl(:) => null()   !< aoi_fld%dqsfcin(item,lan)
    real (kind=kind_phys), pointer :: ulwsfcin_cpl(:)=> null()   !< aoi_fld%ulwsfcin(item,lan)
    !--- only variable needed for cplwav=.TRUE.
    real (kind=kind_phys), pointer :: slimskin_cpl(:)=> null()   !< aoi_fld%slimskin(item,lan)

    !--- outgoing accumulated quantities
    real (kind=kind_phys), pointer :: rain_cpl  (:)  => null()   !< total rain precipitation
    real (kind=kind_phys), pointer :: snow_cpl  (:)  => null()   !< total snow precipitation  
    real (kind=kind_phys), pointer :: dusfc_cpl (:)  => null()   !< sfc u momentum flux
    real (kind=kind_phys), pointer :: dvsfc_cpl (:)  => null()   !< sfc v momentum flux
    real (kind=kind_phys), pointer :: dtsfc_cpl (:)  => null()   !< sfc sensible heat flux
    real (kind=kind_phys), pointer :: dqsfc_cpl (:)  => null()   !< sfc   latent heat flux
    real (kind=kind_phys), pointer :: dlwsfc_cpl(:)  => null()   !< sfc downward lw flux (w/m**2)
    real (kind=kind_phys), pointer :: dswsfc_cpl(:)  => null()   !< sfc downward sw flux (w/m**2)
    real (kind=kind_phys), pointer :: dnirbm_cpl(:)  => null()   !< sfc nir beam downward sw flux (w/m**2) 
    real (kind=kind_phys), pointer :: dnirdf_cpl(:)  => null()   !< sfc nir diff downward sw flux (w/m**2) 
    real (kind=kind_phys), pointer :: dvisbm_cpl(:)  => null()   !< sfc uv+vis beam dnwd sw flux (w/m**2) 
    real (kind=kind_phys), pointer :: dvisdf_cpl(:)  => null()   !< sfc uv+vis diff dnwd sw flux (w/m**2) 
    real (kind=kind_phys), pointer :: nlwsfc_cpl(:)  => null()   !< net downward lw flux (w/m**2)
    real (kind=kind_phys), pointer :: nswsfc_cpl(:)  => null()   !< net downward sw flux (w/m**2)
    real (kind=kind_phys), pointer :: nnirbm_cpl(:)  => null()   !< net nir beam downward sw flux (w/m**2)
    real (kind=kind_phys), pointer :: nnirdf_cpl(:)  => null()   !< net nir diff downward sw flux (w/m**2)
    real (kind=kind_phys), pointer :: nvisbm_cpl(:)  => null()   !< net uv+vis beam downward sw rad flux (w/m**2)
    real (kind=kind_phys), pointer :: nvisdf_cpl(:)  => null()   !< net uv+vis diff downward sw rad flux (w/m**2)

    !--- outgoing instantaneous quantities
    real (kind=kind_phys), pointer :: dusfci_cpl (:) => null()   !< instantaneous sfc u momentum flux
    real (kind=kind_phys), pointer :: dvsfci_cpl (:) => null()   !< instantaneous sfc v momentum flux
    real (kind=kind_phys), pointer :: dtsfci_cpl (:) => null()   !< instantaneous sfc sensible heat flux
    real (kind=kind_phys), pointer :: dqsfci_cpl (:) => null()   !< instantaneous sfc   latent heat flux
    real (kind=kind_phys), pointer :: dlwsfci_cpl(:) => null()   !< instantaneous sfc downward lw flux
    real (kind=kind_phys), pointer :: dswsfci_cpl(:) => null()   !< instantaneous sfc downward sw flux
    real (kind=kind_phys), pointer :: dnirbmi_cpl(:) => null()   !< instantaneous sfc nir beam downward sw flux
    real (kind=kind_phys), pointer :: dnirdfi_cpl(:) => null()   !< instantaneous sfc nir diff downward sw flux
    real (kind=kind_phys), pointer :: dvisbmi_cpl(:) => null()   !< instantaneous sfc uv+vis beam downward sw flux
    real (kind=kind_phys), pointer :: dvisdfi_cpl(:) => null()   !< instantaneous sfc uv+vis diff downward sw flux
    real (kind=kind_phys), pointer :: nlwsfci_cpl(:) => null()   !< instantaneous net sfc downward lw flux
    real (kind=kind_phys), pointer :: nswsfci_cpl(:) => null()   !< instantaneous net sfc downward sw flux
    real (kind=kind_phys), pointer :: nnirbmi_cpl(:) => null()   !< instantaneous net nir beam sfc downward sw flux
    real (kind=kind_phys), pointer :: nnirdfi_cpl(:) => null()   !< instantaneous net nir diff sfc downward sw flux
    real (kind=kind_phys), pointer :: nvisbmi_cpl(:) => null()   !< instantaneous net uv+vis beam downward sw flux
    real (kind=kind_phys), pointer :: nvisdfi_cpl(:) => null()   !< instantaneous net uv+vis diff downward sw flux
    real (kind=kind_phys), pointer :: t2mi_cpl   (:) => null()   !< instantaneous T2m
    real (kind=kind_phys), pointer :: q2mi_cpl   (:) => null()   !< instantaneous Q2m
    real (kind=kind_phys), pointer :: u10mi_cpl  (:) => null()   !< instantaneous U10m
    real (kind=kind_phys), pointer :: v10mi_cpl  (:) => null()   !< instantaneous V10m
    real (kind=kind_phys), pointer :: tsfci_cpl  (:) => null()   !< instantaneous sfc temperature
    real (kind=kind_phys), pointer :: psurfi_cpl (:) => null()   !< instantaneous sfc pressure

    !--- topography-based information for the coupling system
    real (kind=kind_phys), pointer :: oro_cpl    (:) => null()   !< orography          (  oro from GFS_sfcprop_type)
    real (kind=kind_phys), pointer :: slmsk_cpl  (:) => null()   !< Land/Sea/Ice mask  (slmsk from GFS_sfcprop_type)

    !--- stochastic physics
    real (kind=kind_phys), pointer :: shum_wts  (:,:)   => null()  !
    real (kind=kind_phys), pointer :: sppt_wts  (:,:)   => null()  !
    real (kind=kind_phys), pointer :: skebu_wts (:,:)   => null()  !
    real (kind=kind_phys), pointer :: skebv_wts (:,:)   => null()  !
    real (kind=kind_phys), pointer :: vcu_wts   (:,:)   => null()  !
    real (kind=kind_phys), pointer :: vcv_wts   (:,:)   => null()  !

    !--- instantaneous quantities for GoCart and will be accumulated for 3D diagnostics
    real (kind=kind_phys), pointer :: dqdti   (:,:)   => null()  !< instantaneous total moisture tendency (kg/kg/s)
    real (kind=kind_phys), pointer :: cnvqci  (:,:)   => null()  !< instantaneous total convective conensate (kg/kg)
    real (kind=kind_phys), pointer :: upd_mfi (:,:)   => null()  !< instantaneous convective updraft mass flux
    real (kind=kind_phys), pointer :: dwn_mfi (:,:)   => null()  !< instantaneous convective downdraft mass flux
    real (kind=kind_phys), pointer :: det_mfi (:,:)   => null()  !< instantaneous convective detrainment mass flux
    real (kind=kind_phys), pointer :: cldcovi (:,:)   => null()  !< instantaneous 3D cloud fraction

    contains
      procedure :: create  => coupling_create  !<   allocate array data
  end type GFS_coupling_type


!----------------------------------------------------------------------------------
! GFS_control_type
!   model control parameters input from a namelist and/or derived from others
!   list of those that can be modified during the run are at the bottom of the list 
!----------------------------------------------------------------------------------
  type GFS_control_type

    integer              :: me              !< MPI rank designator
    integer              :: master          !< MPI rank of master atmosphere processor
    integer              :: nlunit          !< unit for namelist
    character(len=64)    :: fn_nml          !< namelist filename for surface data cycling
    real(kind=kind_phys) :: fhzero          !< seconds between clearing of diagnostic buckets
    logical              :: ldiag3d         !< flag for 3d diagnostic fields
    logical              :: lssav           !< logical flag for storing diagnostics
    real(kind=kind_phys) :: fhcyc           !< frequency for surface data cycling (secs)
    logical              :: lgocart         !< flag for 3d diagnostic fields for gocart 1
    real(kind=kind_phys) :: fhgoc3d         !< seconds between calls to gocart
    integer              :: thermodyn_id    !< valid for GFS only for get_prs/phi
    integer              :: sfcpress_id     !< valid for GFS only for get_prs/phi
    logical              :: gen_coord_hybrid!< for Henry's gen coord

    !--- set some grid extent parameters
    integer              :: isc             !< starting i-index for this MPI-domain
    integer              :: jsc             !< starting j-index for this MPI-domain
    integer              :: nx              !< number of points in the i-dir for this MPI-domain
    integer              :: ny              !< number of points in the j-dir for this MPI-domain
    integer              :: levs            !< number of vertical levels
    integer              :: cnx             !< number of points in the i-dir for this cubed-sphere face
    integer              :: cny             !< number of points in the j-dir for this cubed-sphere face
    integer              :: lonr            !< number of global points in x-dir (i) along the equator
    integer              :: latr            !< number of global points in y-dir (j) along any meridian

    !--- coupling parameters
    logical              :: cplflx          !< default no cplflx collection
    logical              :: cplwav          !< default no cplwav collection

    !--- integrated dynamics through earth's atmosphere
    logical              :: lsidea         

    !--- calendars and time parameters and activation triggers
    real(kind=kind_phys) :: dtp             !< physics timestep in seconds
    real(kind=kind_phys) :: dtf             !< dynamics timestep in seconds
    integer              :: nscyc           !< trigger for surface data cycling
    integer              :: nszero          !< trigger for zeroing diagnostic buckets
    integer              :: idat(1:8)       !< initialization date and time
                                            !< (yr, mon, day, t-zone, hr, min, sec, mil-sec)
    integer              :: idate(4)        !< initial date with different size and ordering
                                            !< (hr, mon, day, yr)
    !--- radiation control parameters
    real(kind=kind_phys) :: fhswr           !< frequency for shortwave radiation (secs)
    real(kind=kind_phys) :: fhlwr           !< frequency for longwave radiation (secs)
    integer              :: nsswr           !< integer trigger for shortwave radiation
    integer              :: nslwr           !< integer trigger for longwave  radiation
    integer              :: levr            !< number of vertical levels for radiation calculations
    integer              :: nfxr            !< second dimension for fluxr diagnostic variable (radiation)
    logical              :: aero_in         !< aerosol flag for gbphys
    logical              :: lmfshal         !< parameter for radiation
    logical              :: lmfdeep2        !< parameter for radiation
    integer              :: nrcm            !< second dimension of random number stream for RAS
    integer              :: iflip           !< iflip - is not the same as flipv
    integer              :: isol            !< use prescribed solar constant
    integer              :: ico2            !< prescribed global mean value (old opernl)
    integer              :: ialb            !< use climatology alb, based on sfc type
                                            !< 1 => use modis based alb
    integer              :: iems            !< use fixed value of 1.0
    integer              :: iaer            !< default aerosol effect in sw only
    integer              :: iovr_sw         !< sw: max-random overlap clouds
    integer              :: iovr_lw         !< lw: max-random overlap clouds
    integer              :: ictm            !< ictm=0 => use data at initial cond time, if not
                                            !<           available; use latest; no extrapolation.
                                            !< ictm=1 => use data at the forecast time, if not
                                            !<           available; use latest; do extrapolation.
                                            !< ictm=yyyy0 => use yyyy data for the forecast time;
                                            !<           no extrapolation.
                                            !< ictm=yyyy1 = > use yyyy data for the fcst. If needed, 
                                            !<           do extrapolation to match the fcst time.
                                            !< ictm=-1 => use user provided external data for
                                            !<           the fcst time; no extrapolation.
                                            !< ictm=-2 => same as ictm=0, but add seasonal cycle
                                            !<           from climatology; no extrapolation.
    integer              :: isubc_sw        !< sw clouds without sub-grid approximation
    integer              :: isubc_lw        !< lw clouds without sub-grid approximation
                                            !< =1 => sub-grid cloud with prescribed seeds
                                            !< =2 => sub-grid cloud with randomly generated
                                            !< seeds
    logical              :: crick_proof     !< CRICK-Proof cloud water
    logical              :: ccnorm          !< Cloud condensate normalized by cloud cover 
    logical              :: norad_precip    !< radiation precip flag for Ferrier/Moorthi
    logical              :: lwhtr           !< flag to output lw heating rate (Radtend%lwhc)
    logical              :: swhtr           !< flag to output sw heating rate (Radtend%swhc)

    !--- microphysical switch
    integer              :: ncld            !< cnoice of cloud scheme
    !--- Z-C microphysical parameters
    logical              :: zhao_mic        !< flag for Zhao-Carr microphysics
    real(kind=kind_phys) :: psautco(2)      !< [in] auto conversion coeff from ice to snow
    real(kind=kind_phys) :: prautco(2)      !< [in] auto conversion coeff from cloud to rain
    real(kind=kind_phys) :: evpco           !< [in] coeff for evaporation of largescale rain
    real(kind=kind_phys) :: wminco(2)       !< [in] water and ice minimum threshold for Zhao

    !--- M-G microphysical parameters
    integer              :: fprcp           !< no prognostic rain and snow (MG)
    real(kind=kind_phys) :: mg_dcs          !< Morrison-Gettleman microphysics parameters
    real(kind=kind_phys) :: mg_qcvar      
    real(kind=kind_phys) :: mg_ts_auto_ice  !< ice auto conversion time scale

    !--- land/surface model parameters
    integer              :: lsm             !< flag for land surface model lsm=1 for noah lsm
    integer              :: lsoil           !< number of soil layers
    integer              :: ivegsrc         !< ivegsrc = 0   => USGS, 
                                            !< ivegsrc = 1   => IGBP (20 category)
                                            !< ivegsrc = 2   => UMD  (13 category)
    integer              :: isot            !< isot = 0   => Zobler soil type  ( 9 category)
                                            !< isot = 1   => STATSGO soil type (19 category)
    logical              :: mom4ice         !< flag controls mom4 sea ice
    logical              :: use_ufo         !< flag for gcycle surface option

    !--- tuning parameters for physical parameterizations
    logical              :: ras             !< flag for ras convection scheme
    logical              :: flipv           !< flag for vertical direction flip (ras)
                                            !< .true. implies surface at k=1
    logical              :: trans_trac      !< flag for convective transport of tracers (RAS only)
    logical              :: old_monin       !< flag for diff monin schemes
    logical              :: cnvgwd          !< flag for conv gravity wave drag
    logical              :: mstrat          !< flag for moorthi approach for stratus
    logical              :: moist_adj       !< flag for moist convective adjustment
    logical              :: cscnv           !< flag for Chikira-Sugiyama convection
    logical              :: cal_pre         !< flag controls precip type algorithm
    logical              :: do_aw           !< AW scale-aware option in cs convection
    logical              :: do_shoc         !< flag for SHOC
    logical              :: shocaftcnv      !< flag for SHOC
    logical              :: shoc_cld        !< flag for clouds
    logical              :: uni_cld         !< flag for clouds in grrad
    logical              :: h2o_phys        !< flag for stratosphere h2o
    logical              :: pdfcld          !< flag for pdfcld
    logical              :: shcnvcw         !< flag for shallow convective cloud
    logical              :: redrag          !< flag for reduced drag coeff. over sea
    logical              :: hybedmf         !< flag for hybrid edmf pbl scheme
    logical              :: dspheat         !< flag for tke dissipative heating
    logical              :: cnvcld        
    logical              :: random_clds     !< flag controls whether clouds are random
    logical              :: shal_cnv        !< flag for calling shallow convection
    integer              :: imfshalcnv      !< flag for mass-flux shallow convection scheme
                                            !<     1: July 2010 version of mass-flux shallow conv scheme
                                            !<         current operational version as of 2016
                                            !<     2: scale- & aerosol-aware mass-flux shallow conv scheme (2017)
                                            !<     0: modified Tiedtke's eddy-diffusion shallow conv scheme
                                            !<    -1: no shallow convection used
    integer              :: imfdeepcnv      !< flag for mass-flux deep convection scheme
                                            !<     1: July 2010 version of SAS conv scheme
                                            !<           current operational version as of 2016
                                            !<     2: scale- & aerosol-aware mass-flux deep conv scheme (2017)
                                            !<     0: old SAS Convection scheme before July 2010
    integer              :: nmtvr           !< number of topographic variables such as variance etc
                                            !< used in the GWD parameterization
    integer              :: jcap            !< number of spectral wave trancation used only by sascnv shalcnv
    real(kind=kind_phys) :: cs_parm(10)     !< tunable parameters for Chikira-Sugiyama convection
    real(kind=kind_phys) :: flgmin(2)       !< [in] ice fraction bounds
    real(kind=kind_phys) :: cgwf(2)         !< multiplication factor for convective GWD
    real(kind=kind_phys) :: ccwf(2)         !< multiplication factor for critical cloud
                                            !< workfunction for RAS
    real(kind=kind_phys) :: cdmbgwd(2)      !< multiplication factors for cdmb and gwd
    real(kind=kind_phys) :: sup             !< supersaturation in pdf cloud when t is very low
    real(kind=kind_phys) :: ctei_rm(2)      !< critical cloud top entrainment instability criteria 
                                            !< (used if mstrat=.true.)
    real(kind=kind_phys) :: crtrh(3)        !< critical relative humidity at the surface
                                            !< PBL top and at the top of the atmosphere
    real(kind=kind_phys) :: dlqf(2)         !< factor for cloud condensate detrainment 
                                            !< from cloud edges for RAS
    integer              :: seed0           !< random seed for radiation

    !--- Rayleigh friction
    real(kind=kind_phys) :: prslrd0         !< pressure level from which Rayleigh Damping is applied
    real(kind=kind_phys) :: ral_ts          !< time scale for Rayleigh damping in days

    !--- near surface temperature model
    logical              :: nst_anl         !< flag for NSSTM analysis in gcycle/sfcsub
    integer              :: lsea           
    real(kind=kind_phys) :: xkzm_m          !< [in] bkgd_vdif_m  background vertical diffusion for momentum  
    real(kind=kind_phys) :: xkzm_h          !< [in] bkgd_vdif_h  background vertical diffusion for heat q  
    real(kind=kind_phys) :: xkzm_s          !< [in] bkgd_vdif_s  sigma threshold for background mom. diffusion  
    integer              :: nstf_name(5)    !< flag 0 for no nst  1 for uncoupled nst  and 2 for coupled NST
                                            !< nstf_name contains the NSST related parameters
                                            !< nstf_name(1) : 0 = NSSTM off, 1 = NSSTM on but uncoupled, 2 =
                                            !< nstf_name(2) : 1 = NSSTM spin up on, 0 = NSSTM spin up off
                                            !< nstf_name(3) : 1 = NSST analysis on, 0 = NSSTM analysis off
                                            !< nstf_name(4) : zsea1 in mm
                                            !< nstf_name(5) : zsea2 in mm
     
    !--- stochastic physics control parameters
    logical              :: do_sppt
    logical              :: do_shum
    logical              :: do_skeb
    logical              :: do_vc  
    real(kind=kind_phys) :: sppt(5)         !< stochastic physics tendency amplitude
    real(kind=kind_phys) :: shum(5)         !< stochastic boundary layer spf hum amp
    real(kind=kind_phys) :: skeb(5)         !< stochastic KE backscatter amplitude
    real(kind=kind_phys) :: vcamp(5)        !< stochastic vorticity confinment amp
    real(kind=kind_phys) :: vc              !< deterministic vorticity confinement parameter.

    !--- tracer handling
    character(len=32), pointer :: tracer_names(:) !< array of initialized tracers from dynamic core
    integer              :: ntrac           !< number of tracers
    integer              :: ntoz            !< tracer index for ozone mixing ratio
    integer              :: ntcw            !< tracer index for cloud condensate (or liquid water)
    integer              :: ntiw            !< tracer index for  ice water
    integer              :: ntrw            !< tracer index for rain water
    integer              :: ntsw            !< tracer index for snow water
    integer              :: ntgl            !< tracer index for graupel
    integer              :: ntlnc           !< tracer index for liquid number concentration
    integer              :: ntinc           !< tracer index for ice    number concentration
    integer              :: ntrnc           !< tracer index for rain   number concentration
    integer              :: ntsnc           !< tracer index for snow   number concentration
    integer              :: ntke            !< tracer index for kinetic energy
    integer              :: nto             !< tracer index for oxygen ion
    integer              :: nto2            !< tracer index for oxygen
 
    !--- derived totals for phy_f*d
    integer              :: ntot2d          !< total number of variables for phyf2d
    integer              :: ntot3d          !< total number of variables for phyf3d
    integer              :: num_p2d         !< number of 2D arrays needed for microphysics
    integer              :: num_p3d         !< number of 3D arrays needed for microphysics
    integer              :: nshoc_2d        !< number of 2d fields for SHOC
    integer              :: nshoc_3d        !< number of 3d fields for SHOC
    integer              :: ncnvcld3d       !< number of convective 3d clouds fields
    integer              :: npdf3d          !< number of 3d arrays associated with pdf based clouds/microphysics
    integer              :: nctp            !< number of cloud types in Chikira-Sugiyama scheme

    !--- debug flag
    logical              :: debug         
    logical              :: pre_rad         !< flag for testing purpose

    !--- variables modified at each time step
    integer              :: ipt             !< index for diagnostic printout point
    logical              :: lprnt           !< control flag for diagnostic print out
    logical              :: lsswr           !< logical flags for sw radiation calls
    logical              :: lslwr           !< logical flags for lw radiation calls
    real(kind=kind_phys) :: solhr           !< hour time after 00z at the t-step
    real(kind=kind_phys) :: solcon          !< solar constant (sun-earth distant adjusted)  [set via radupdate]
    real(kind=kind_phys) :: slag            !< equation of time ( radian )                  [set via radupdate]
    real(kind=kind_phys) :: sdec            !< sin of the solar declination angle           [set via radupdate]
    real(kind=kind_phys) :: cdec            !< cos of the solar declination angle           [set via radupdate]
    real(kind=kind_phys) :: clstp           !< index used by cnvc90 (for convective clouds) 
                                            !< legacy stuff - does not affect forecast
    real(kind=kind_phys) :: phour           !< previous forecast hour
    real(kind=kind_phys) :: fhour           !< curent forecast hour
    real(kind=kind_phys) :: zhour           !< previous hour diagnostic buckets emptied
    integer              :: kdt             !< current forecast iteration
    integer              :: jdat(1:8)       !< current forecast date and time
                                            !< (yr, mon, day, t-zone, hr, min, sec, mil-sec)

    contains
      procedure :: init  => control_initialize
      procedure :: print => control_print
  end type GFS_control_type


!--------------------------------------------------------------------
! GFS_grid_type
!   grid data needed for interpolations and length-scale calculations
!--------------------------------------------------------------------
  type GFS_grid_type
    
    real (kind=kind_phys), pointer :: xlon   (:)    => null()   !< grid longitude in radians, ok for both 0->2pi  
                                                                !! or -pi -> +pi ranges   
    real (kind=kind_phys), pointer :: xlat   (:)    => null()   !< grid latitude in radians, default to pi/2 ->  
                                                                !! -pi/2 range, otherwise adj in subr called   
    real (kind=kind_phys), pointer :: xlat_d (:)    => null()   !< grid latitude in degrees, default to 90 -> 
                                                                !! -90 range, otherwise adj in subr called   
    real (kind=kind_phys), pointer :: sinlat (:)    => null()   !< sine of the grids corresponding latitudes   
    real (kind=kind_phys), pointer :: coslat (:)    => null()   !< cosine of the grids corresponding latitudes   
    real (kind=kind_phys), pointer :: area   (:)    => null()   !< area of the grid cell   
    real (kind=kind_phys), pointer :: dx     (:)    => null()   !< relative dx for the grid cell

    !--- grid-related interpolation data for prognostic ozone
    real (kind=kind_phys), pointer :: ddy_o3    (:) => null()   !< interpolation     weight for ozone
    integer,               pointer :: jindx1_o3 (:) => null()   !< interpolation  low index for ozone
    integer,               pointer :: jindx2_o3 (:) => null()   !< interpolation high index for ozone

    !--- grid-related interpolation data for stratosphere water
    real (kind=kind_phys), pointer :: ddy_h     (:) => null()   !< interpolation     weight for h2o
    integer,               pointer :: jindx1_h  (:) => null()   !< interpolation  low index for h2o
    integer,               pointer :: jindx2_h  (:) => null()   !< interpolation high index for h2o
    contains
      procedure :: create   => grid_create   !<   allocate array data
  end type GFS_grid_type


!-----------------------------------------------
! GFS_tbd_type
!   data not yet assigned to a defined container
!-----------------------------------------------
  type GFS_tbd_type

    !--- radiation random seeds
    integer,               pointer :: icsdsw   (:)     => null()  !< (rad. only) auxiliary cloud control arrays passed to main
    integer,               pointer :: icsdlw   (:)     => null()  !< (rad. only) radiations. if isubcsw/isubclw (input to init) 
                                                                  !< (rad. only) are set to 2, the arrays contains provided     
                                                                  !< (rad. only) random seeds for sub-column clouds generators

    !--- In
    real (kind=kind_phys), pointer :: ozpl     (:,:,:) => null()  !< ozone forcing data
    real (kind=kind_phys), pointer :: h2opl    (:,:,:) => null()  !< water forcing data

    !--- active when ((.not. newsas .or. cal_pre) .and. random_clds)
    real (kind=kind_phys), pointer :: rann     (:,:)   => null()  !< random number array (0-1)

    !--- In/Out
    real (kind=kind_phys), pointer :: acv      (:)     => null()  !< array containing accumulated convective clouds
    real (kind=kind_phys), pointer :: acvb     (:)     => null()  !< arrays used by cnvc90 bottom
    real (kind=kind_phys), pointer :: acvt     (:)     => null()  !< arrays used by cnvc90 top (cnvc90.f)

    !--- Stochastic physics properties calculated in physics_driver
    real (kind=kind_phys), pointer :: dtdtr     (:,:)   => null()  !< temperature change due to radiative heating per time step (K)
    real (kind=kind_phys), pointer :: dtotprcp  (:)     => null()  !< change in totprcp  (diag_type)
    real (kind=kind_phys), pointer :: dcnvprcp  (:)     => null()  !< change in cnvprcp  (diag_type)
    real (kind=kind_phys), pointer :: drain_cpl (:)     => null()  !< change in rain_cpl (coupling_type)
    real (kind=kind_phys), pointer :: dsnow_cpl (:)     => null()  !< change in show_cpl (coupling_type)

    !--- phy_f*d variables needed for seamless restarts and moving data between grrad and gbphys
    real (kind=kind_phys), pointer :: phy_fctd (:,:)   => null()  !< For CS convection
    real (kind=kind_phys), pointer :: phy_f2d  (:,:)   => null()  !< 2d arrays saved for restart
    real (kind=kind_phys), pointer :: phy_f3d  (:,:,:) => null()  !< 3d arrays saved for restart

    integer                        :: blkno                       !< for explicit data blocking: block number of this block
    integer,               pointer :: blksz(:)         => null()  !< for explicit data blocking: horizontal block sizes of all blocks

    contains
      procedure :: create  => tbd_create  !<   allocate array data
  end type GFS_tbd_type


!------------------------------------------------------------------
! GFS_cldprop_type
!  cloud properties and tendencies needed by radiation from physics 
!------------------------------------------------------------------
  type GFS_cldprop_type

    !--- In     (radiation)
    !--- In/Out (physics)
    real (kind=kind_phys), pointer :: cv  (:)     => null()  !< fraction of convective cloud ; phys
    real (kind=kind_phys), pointer :: cvt (:)     => null()  !< convective cloud top pressure in pa ; phys
    real (kind=kind_phys), pointer :: cvb (:)     => null()  !< convective cloud bottom pressure in pa ; phys, cnvc90

    contains
      procedure :: create  => cldprop_create  !<   allocate array data
  end type GFS_cldprop_type


!-----------------------------------------
! GFS_radtend_type
!   radiation tendencies needed by physics
!-----------------------------------------
  type GFS_radtend_type

    type (sfcfsw_type),    pointer :: sfcfsw(:)   => null()   !< sw radiation fluxes at sfc
                                                              !< [dim(im): created in grrad.f], components:        
                                                              !!     (check module_radsw_parameters for definition)   
                                                              !!\n   %upfxc - total sky upward sw flux at sfc (w/m**2)     
                                                              !!\n   %upfx0 - clear sky upward sw flux at sfc (w/m**2)     
                                                              !!\n   %dnfxc - total sky downward sw flux at sfc (w/m**2)   
                                                              !!\n   %dnfx0 - clear sky downward sw flux at sfc (w/m**2)   

    type (sfcflw_type),    pointer :: sfcflw(:)    => null()  !< lw radiation fluxes at sfc
                                                              !< [dim(im): created in grrad.f], components:         
                                                              !!     (check module_radlw_paramters for definition)  
                                                              !!\n   %upfxc - total sky upward lw flux at sfc (w/m**2)     
                                                              !!\n   %upfx0 - clear sky upward lw flux at sfc (w/m**2)     
                                                              !!\n   %dnfxc - total sky downward lw flux at sfc (w/m**2)   
                                                              !!\n   %dnfx0 - clear sky downward lw flux at sfc (w/m**2)   

    !--- Out (radiation only)
    real (kind=kind_phys), pointer :: htrsw (:,:)  => null()  !< swh  total sky sw heating rate in k/sec 
    real (kind=kind_phys), pointer :: htrlw (:,:)  => null()  !< hlw  total sky lw heating rate in k/sec
    real (kind=kind_phys), pointer :: sfalb (:)    => null()  !< mean surface diffused sw albedo 

    real (kind=kind_phys), pointer :: coszen(:)    => null()  !< mean cos of zenith angle over rad call period 
    real (kind=kind_phys), pointer :: tsflw (:)    => null()  !< surface air temp during lw calculation in k 
    real (kind=kind_phys), pointer :: semis (:)    => null()  !< surface lw emissivity in fraction

    !--- In/Out (???) (radiaition only)
    real (kind=kind_phys), pointer :: coszdg(:)    => null()  !< daytime mean cosz over rad call period

    !--- In/Out (???) (physics only)
    real (kind=kind_phys), pointer :: swhc (:,:)   => null()  !< clear sky sw heating rates ( k/s ) 
    real (kind=kind_phys), pointer :: lwhc (:,:)   => null()  !< clear sky lw heating rates ( k/s ) 
    real (kind=kind_phys), pointer :: lwhd (:,:,:) => null()  !< idea sky lw heating rates ( k/s ) 

    contains
      procedure :: create  => radtend_create   !<   allocate array data
  end type GFS_radtend_type

!----------------------------------------------------------------
! GFS_diag_type
!  internal diagnostic type used as arguments to gbphys and grrad 
!----------------------------------------------------------------
  type GFS_diag_type

    !! Input/Output only in radiation
    real (kind=kind_phys), pointer :: fluxr (:,:)   => null()   !< to save time accumulated 2-d fields defined as:!
                                                                !< hardcoded field indices, opt. includes aerosols!
    type (topfsw_type),    pointer :: topfsw(:)     => null()   !< sw radiation fluxes at toa, components:        
                                               !       %upfxc    - total sky upward sw flux at toa (w/m**2)     
                                               !       %dnfxc    - total sky downward sw flux at toa (w/m**2)   
                                               !       %upfx0    - clear sky upward sw flux at toa (w/m**2)     
    type (topflw_type),    pointer :: topflw(:)     => null()   !< lw radiation fluxes at top, component:
                                               !       %upfxc    - total sky upward lw flux at toa (w/m**2)
                                               !       %upfx0    - clear sky upward lw flux at toa (w/m**2)

    ! Input/output - used by physics
    real (kind=kind_phys), pointer :: srunoff(:)    => null()   !< surface water runoff (from lsm)
    real (kind=kind_phys), pointer :: evbsa  (:)    => null()   !< noah lsm diagnostics
    real (kind=kind_phys), pointer :: evcwa  (:)    => null()   !< noah lsm diagnostics
    real (kind=kind_phys), pointer :: snohfa (:)    => null()   !< noah lsm diagnostics
    real (kind=kind_phys), pointer :: transa (:)    => null()   !< noah lsm diagnostics
    real (kind=kind_phys), pointer :: sbsnoa (:)    => null()   !< noah lsm diagnostics
    real (kind=kind_phys), pointer :: snowca (:)    => null()   !< noah lsm diagnostics
    real (kind=kind_phys), pointer :: soilm  (:)    => null()   !< soil moisture
    real (kind=kind_phys), pointer :: tmpmin (:)    => null()   !< min temperature at 2m height (k)
    real (kind=kind_phys), pointer :: tmpmax (:)    => null()   !< max temperature at 2m height (k)
    real (kind=kind_phys), pointer :: dusfc  (:)    => null()   !< u component of surface stress
    real (kind=kind_phys), pointer :: dvsfc  (:)    => null()   !< v component of surface stress
    real (kind=kind_phys), pointer :: dtsfc  (:)    => null()   !< sensible heat flux (w/m2)
    real (kind=kind_phys), pointer :: dqsfc  (:)    => null()   !< latent heat flux (w/m2)
    real (kind=kind_phys), pointer :: totprcp(:)    => null()   !< accumulated total precipitation (kg/m2)
    real (kind=kind_phys), pointer :: gflux  (:)    => null()   !< groud conductive heat flux
    real (kind=kind_phys), pointer :: dlwsfc (:)    => null()   !< time accumulated sfc dn lw flux ( w/m**2 )
    real (kind=kind_phys), pointer :: ulwsfc (:)    => null()   !< time accumulated sfc up lw flux ( w/m**2 )
    real (kind=kind_phys), pointer :: suntim (:)    => null()   !< sunshine duration time (s)
    real (kind=kind_phys), pointer :: runoff (:)    => null()   !< total water runoff
    real (kind=kind_phys), pointer :: ep     (:)    => null()   !< potential evaporation
    real (kind=kind_phys), pointer :: cldwrk (:)    => null()   !< cloud workfunction (valid only with sas)
    real (kind=kind_phys), pointer :: dugwd  (:)    => null()   !< vertically integrated u change by OGWD
    real (kind=kind_phys), pointer :: dvgwd  (:)    => null()   !< vertically integrated v change by OGWD
    real (kind=kind_phys), pointer :: psmean (:)    => null()   !< surface pressure (kPa)
    real (kind=kind_phys), pointer :: cnvprcp(:)    => null()   !< accumulated convective precipitation (kg/m2)
    real (kind=kind_phys), pointer :: spfhmin(:)    => null()   !< minimum specific humidity
    real (kind=kind_phys), pointer :: spfhmax(:)    => null()   !< maximum specific humidity
    real (kind=kind_phys), pointer :: rain   (:)    => null()   !< total rain at this time step
    real (kind=kind_phys), pointer :: rainc  (:)    => null()   !< convective rain at this time step
    real (kind=kind_phys), pointer :: ice    (:)    => null()   !< ice fall at this time step
    real (kind=kind_phys), pointer :: snow   (:)    => null()   !< snow fall at this time step
    real (kind=kind_phys), pointer :: graupel(:)    => null()   !< graupel fall at this time step
    real (kind=kind_phys), pointer :: totice (:)    => null()   !< accumulated ice precipitation (kg/m2)
    real (kind=kind_phys), pointer :: totsnw (:)    => null()   !< accumulated snow precipitation (kg/m2)
    real (kind=kind_phys), pointer :: totgrp (:)    => null()   !< accumulated graupel precipitation (kg/m2)

    ! Output - only in physics
    real (kind=kind_phys), pointer :: u10m   (:)    => null()   !< 10 meater u/v wind speed
    real (kind=kind_phys), pointer :: v10m   (:)    => null()   !< 10 meater u/v wind speed
    real (kind=kind_phys), pointer :: zlvl   (:)    => null()   !< layer 1 height (m)
    real (kind=kind_phys), pointer :: psurf  (:)    => null()   !< surface pressure (Pa)
    real (kind=kind_phys), pointer :: hpbl   (:)    => null()   !< pbl height (m)
    real (kind=kind_phys), pointer :: pwat   (:)    => null()   !< precipitable water
    real (kind=kind_phys), pointer :: t1     (:)    => null()   !< layer 1 temperature (K)
    real (kind=kind_phys), pointer :: q1     (:)    => null()   !< layer 1 specific humidity (kg/kg)
    real (kind=kind_phys), pointer :: u1     (:)    => null()   !< layer 1 zonal wind (m/s)
    real (kind=kind_phys), pointer :: v1     (:)    => null()   !< layer 1 merdional wind (m/s)
    real (kind=kind_phys), pointer :: chh    (:)    => null()   !< thermal exchange coefficient
    real (kind=kind_phys), pointer :: cmm    (:)    => null()   !< momentum exchange coefficient
    real (kind=kind_phys), pointer :: dlwsfci(:)    => null()   !< instantaneous sfc dnwd lw flux ( w/m**2 )
    real (kind=kind_phys), pointer :: ulwsfci(:)    => null()   !< instantaneous sfc upwd lw flux ( w/m**2 )
    real (kind=kind_phys), pointer :: dswsfci(:)    => null()   !< instantaneous sfc dnwd sw flux ( w/m**2 )
    real (kind=kind_phys), pointer :: uswsfci(:)    => null()   !< instantaneous sfc upwd sw flux ( w/m**2 )
    real (kind=kind_phys), pointer :: dusfci (:)    => null()   !< instantaneous u component of surface stress
    real (kind=kind_phys), pointer :: dvsfci (:)    => null()   !< instantaneous v component of surface stress
    real (kind=kind_phys), pointer :: dtsfci (:)    => null()   !< instantaneous sfc sensible heat flux
    real (kind=kind_phys), pointer :: dqsfci (:)    => null()   !< instantaneous sfc latent heat flux
    real (kind=kind_phys), pointer :: gfluxi (:)    => null()   !< instantaneous sfc ground heat flux
    real (kind=kind_phys), pointer :: epi    (:)    => null()   !< instantaneous sfc potential evaporation
    real (kind=kind_phys), pointer :: smcwlt2(:)    => null()   !< wilting point (volumetric)
    real (kind=kind_phys), pointer :: smcref2(:)    => null()   !< soil moisture threshold (volumetric)
    real (kind=kind_phys), pointer :: wet1   (:)    => null()   !< normalized soil wetness
    real (kind=kind_phys), pointer :: sr     (:)    => null()   !< snow ratio : ratio of snow to total precipitation

    !--- accumulated quantities for 3D diagnostics
    real (kind=kind_phys), pointer :: du3dt (:,:,:) => null()   !< u momentum change due to physics
    real (kind=kind_phys), pointer :: dv3dt (:,:,:) => null()   !< v momentum change due to physics
    real (kind=kind_phys), pointer :: dt3dt (:,:,:) => null()   !< temperature change due to physics
    real (kind=kind_phys), pointer :: dq3dt (:,:,:) => null()   !< moisture change due to physics
 
    !--- accumulated quantities for 3D diagnostics
    real (kind=kind_phys), pointer :: upd_mf (:,:)   => null()  !< instantaneous convective updraft mass flux
    real (kind=kind_phys), pointer :: dwn_mf (:,:)   => null()  !< instantaneous convective downdraft mass flux
    real (kind=kind_phys), pointer :: det_mf (:,:)   => null()  !< instantaneous convective detrainment mass flux
    real (kind=kind_phys), pointer :: cldcov (:,:)   => null()  !< instantaneous 3D cloud fraction

    contains
      procedure create    => diag_create
      procedure rad_zero  => diag_rad_zero
      procedure phys_zero => diag_phys_zero
  end type GFS_diag_type

!----------------
! PUBLIC ENTITIES
!----------------
  public GFS_init_type
  public GFS_statein_type,  GFS_stateout_type, GFS_sfcprop_type, &
         GFS_coupling_type, GFS_sfccycle_type
  public GFS_control_type,  GFS_grid_type,     GFS_tbd_type, &
         GFS_cldprop_type,  GFS_radtend_type,  GFS_diag_type

!*******************************************************************************************
  CONTAINS

!------------------------
! GFS_statein_type%create
!------------------------
  subroutine statein_create (Statein, IM, Model) 
    implicit none

    class(GFS_statein_type)             :: Statein
    integer,                 intent(in) :: IM
    type(GFS_control_type),  intent(in) :: Model

    !--- level geopotential and pressures
    allocate (Statein%phii  (IM,Model%levs+1))
    allocate (Statein%prsi  (IM,Model%levs+1))
    allocate (Statein%prsik (IM,Model%levs+1))

    Statein%phii  = clear_val
    Statein%prsi  = clear_val
    Statein%prsik = clear_val

    !--- layer geopotential and pressures
    allocate (Statein%phil  (IM,Model%levs))
    allocate (Statein%prsl  (IM,Model%levs))
    allocate (Statein%prslk (IM,Model%levs))

    Statein%phil  = clear_val
    Statein%prsl  = clear_val
    Statein%prslk = clear_val

    !--- shared radiation and physics variables
    allocate (Statein%vvl  (IM,Model%levs))
    allocate (Statein%tgrs (IM,Model%levs))

    Statein%vvl  = clear_val
    Statein%tgrs = clear_val

    !--- physics only variables
    allocate (Statein%pgr    (IM))
    allocate (Statein%ugrs   (IM,Model%levs))
    allocate (Statein%vgrs   (IM,Model%levs))
    allocate (Statein%qgrs   (IM,Model%levs,Model%ntrac))

    Statein%qgrs   = clear_val
    Statein%pgr    = clear_val
    Statein%ugrs   = clear_val
    Statein%vgrs   = clear_val

  end subroutine statein_create


!-------------------------
! GFS_stateout_type%create
!-------------------------
  subroutine stateout_create (Stateout, IM, Model)

    implicit none

    class(GFS_stateout_type)           :: Stateout
    integer,                intent(in) :: IM
    type(GFS_control_type), intent(in) :: Model

    allocate (Stateout%gu0 (IM,Model%levs))
    allocate (Stateout%gv0 (IM,Model%levs))
    allocate (Stateout%gt0 (IM,Model%levs))
    allocate (Stateout%gq0 (IM,Model%levs,Model%ntrac))

    Stateout%gu0 = clear_val
    Stateout%gv0 = clear_val
    Stateout%gt0 = clear_val
    Stateout%gq0 = clear_val

 end subroutine stateout_create


!------------------------
! GFS_sfcprop_type%create
!------------------------
  subroutine sfcprop_create (Sfcprop, IM, Model)
                                
    implicit none

    class(GFS_sfcprop_type)            :: Sfcprop
    integer,                intent(in) :: IM
    type(GFS_control_type), intent(in) :: Model

    !--- physics and radiation
    allocate (Sfcprop%slmsk  (IM))
    allocate (Sfcprop%tsfc   (IM))
    allocate (Sfcprop%tisfc  (IM))
    allocate (Sfcprop%snowd  (IM))
    allocate (Sfcprop%zorl   (IM))
    allocate (Sfcprop%fice   (IM))
    allocate (Sfcprop%hprim  (IM))
    allocate (Sfcprop%hprime (IM,Model%nmtvr))

    Sfcprop%slmsk   = clear_val
    Sfcprop%tsfc    = clear_val
    Sfcprop%tisfc   = clear_val
    Sfcprop%snowd   = clear_val
    Sfcprop%zorl    = clear_val
    Sfcprop%fice    = clear_val
    Sfcprop%hprim   = clear_val
    Sfcprop%hprime  = clear_val

    !--- In (radiation only)
    allocate (Sfcprop%sncovr (IM))
    allocate (Sfcprop%snoalb (IM))
    allocate (Sfcprop%alvsf  (IM))
    allocate (Sfcprop%alnsf  (IM))
    allocate (Sfcprop%alvwf  (IM))
    allocate (Sfcprop%alnwf  (IM))
    allocate (Sfcprop%facsf  (IM))
    allocate (Sfcprop%facwf  (IM))

    Sfcprop%sncovr = clear_val
    Sfcprop%snoalb = clear_val
    Sfcprop%alvsf  = clear_val
    Sfcprop%alnsf  = clear_val
    Sfcprop%alvwf  = clear_val
    Sfcprop%alnwf  = clear_val
    Sfcprop%facsf  = clear_val
    Sfcprop%facwf  = clear_val

    !--- physics surface props
    !--- In
    allocate (Sfcprop%slope   (IM))
    allocate (Sfcprop%shdmin  (IM))
    allocate (Sfcprop%shdmax  (IM))
    allocate (Sfcprop%snoalb  (IM))
    allocate (Sfcprop%tg3     (IM))
    allocate (Sfcprop%vfrac   (IM))
    allocate (Sfcprop%vtype   (IM))
    allocate (Sfcprop%stype   (IM))
    allocate (Sfcprop%uustar  (IM))
    allocate (Sfcprop%oro     (IM))
    allocate (Sfcprop%oro_uf  (IM))

    Sfcprop%slope   = clear_val
    Sfcprop%shdmin  = clear_val
    Sfcprop%shdmax  = clear_val
    Sfcprop%snoalb  = clear_val
    Sfcprop%tg3     = clear_val
    Sfcprop%vfrac   = clear_val
    Sfcprop%vtype   = clear_val
    Sfcprop%stype   = clear_val
    Sfcprop%uustar  = clear_val
    Sfcprop%oro     = clear_val
    Sfcprop%oro_uf  = clear_val

    !--- In/Out
    allocate (Sfcprop%hice   (IM))
    allocate (Sfcprop%weasd  (IM))
    allocate (Sfcprop%sncovr (IM))
    allocate (Sfcprop%canopy (IM))
    allocate (Sfcprop%ffmm   (IM))
    allocate (Sfcprop%ffhh   (IM))
    allocate (Sfcprop%f10m   (IM))
    allocate (Sfcprop%tprcp  (IM))
    allocate (Sfcprop%srflag (IM))
    allocate (Sfcprop%slc    (IM,Model%lsoil))
    allocate (Sfcprop%smc    (IM,Model%lsoil))
    allocate (Sfcprop%stc    (IM,Model%lsoil))

    Sfcprop%hice   = clear_val
    Sfcprop%weasd  = clear_val
    Sfcprop%sncovr = clear_val
    Sfcprop%canopy = clear_val
    Sfcprop%ffmm   = clear_val
    Sfcprop%ffhh   = clear_val
    Sfcprop%f10m   = clear_val
    Sfcprop%tprcp  = clear_val
    Sfcprop%srflag = clear_val
    Sfcprop%slc    = clear_val
    Sfcprop%smc    = clear_val
    Sfcprop%stc    = clear_val

    !--- Out
    allocate (Sfcprop%t2m (IM))
    allocate (Sfcprop%q2m (IM))

    Sfcprop%t2m = clear_val
    Sfcprop%q2m = clear_val

    if (Model%nstf_name(1) > 0) then
      allocate (Sfcprop%tref   (IM))
      allocate (Sfcprop%z_c    (IM))
      allocate (Sfcprop%c_0    (IM))
      allocate (Sfcprop%c_d    (IM))
      allocate (Sfcprop%w_0    (IM))
      allocate (Sfcprop%w_d    (IM))
      allocate (Sfcprop%xt     (IM))
      allocate (Sfcprop%xs     (IM))
      allocate (Sfcprop%xu     (IM))
      allocate (Sfcprop%xv     (IM))
      allocate (Sfcprop%xz     (IM))
      allocate (Sfcprop%zm     (IM))
      allocate (Sfcprop%xtts   (IM))
      allocate (Sfcprop%xzts   (IM))
      allocate (Sfcprop%d_conv (IM))
      allocate (Sfcprop%ifd    (IM))
      allocate (Sfcprop%dt_cool(IM))
      allocate (Sfcprop%qrain  (IM))

      Sfcprop%tref    = zero
      Sfcprop%z_c     = zero
      Sfcprop%c_0     = zero
      Sfcprop%c_d     = zero
      Sfcprop%w_0     = zero
      Sfcprop%w_d     = zero
      Sfcprop%xt      = zero
      Sfcprop%xs      = zero
      Sfcprop%xu      = zero
      Sfcprop%xv      = zero
      Sfcprop%xz      = zero
      Sfcprop%zm      = zero
      Sfcprop%xtts    = zero
      Sfcprop%xzts    = zero
      Sfcprop%d_conv  = zero
      Sfcprop%ifd     = zero
      Sfcprop%dt_cool = zero
      Sfcprop%qrain   = zero
    endif

  end subroutine sfcprop_create


  !-------------------------
  ! GFS_sfccycle_type%create
  !-------------------------
  subroutine sfccycle_create (Sfccycle, IM, Model)

    implicit none

    class(GFS_sfccycle_type)           :: Sfccycle
    integer,                intent(in) :: IM
    type(GFS_control_type), intent(in) :: Model

    Sfccycle%ifp = 0

    allocate (Sfccycle%glacir  (IM))
    allocate (Sfccycle%amxice  (IM))
    allocate (Sfccycle%tsfcl0  (IM))

    Sfccycle%glacir  = clear_val
    Sfccycle%amxice  = clear_val
    Sfccycle%tsfcl0  = clear_val

    call Sfccycle%clima%create (IM, Model%lsoil)

  end subroutine sfccycle_create


!-------------------------
! GFS_coupling_type%create
!-------------------------
  subroutine coupling_create (Coupling, IM, Model)

    implicit none

    class(GFS_coupling_type)           :: Coupling
    integer,                intent(in) :: IM
    type(GFS_control_type), intent(in) :: Model

    !--- radiation out
    !--- physics in
    allocate (Coupling%nirbmdi  (IM))
    allocate (Coupling%nirdfdi  (IM))
    allocate (Coupling%visbmdi  (IM))   
    allocate (Coupling%visdfdi  (IM))   
    allocate (Coupling%nirbmui  (IM))   
    allocate (Coupling%nirdfui  (IM))   
    allocate (Coupling%visbmui  (IM))   
    allocate (Coupling%visdfui  (IM))   

    Coupling%nirbmdi = clear_val
    Coupling%nirdfdi = clear_val
    Coupling%visbmdi = clear_val
    Coupling%visdfdi = clear_val
    Coupling%nirbmui = clear_val
    Coupling%nirdfui = clear_val
    Coupling%visbmui = clear_val
    Coupling%visdfui = clear_val

    allocate (Coupling%sfcdsw    (IM))
    allocate (Coupling%sfcnsw    (IM))
    allocate (Coupling%sfcdlw    (IM))

    Coupling%sfcdsw    = clear_val
    Coupling%sfcnsw    = clear_val
    Coupling%sfcdlw    = clear_val

    if (Model%cplflx .or. Model%do_sppt) then
      allocate (Coupling%rain_cpl     (IM))
      allocate (Coupling%snow_cpl     (IM))

      Coupling%rain_cpl     = clear_val
      Coupling%snow_cpl     = clear_val
    endif

    if (Model%cplflx) then
      !--- incoming quantities
      allocate (Coupling%slimskin_cpl (IM))
      allocate (Coupling%dusfcin_cpl  (IM))
      allocate (Coupling%dvsfcin_cpl  (IM))
      allocate (Coupling%dtsfcin_cpl  (IM))
      allocate (Coupling%dqsfcin_cpl  (IM))
      allocate (Coupling%ulwsfcin_cpl (IM))

      Coupling%slimskin_cpl = clear_val
      Coupling%dusfcin_cpl  = clear_val
      Coupling%dvsfcin_cpl  = clear_val
      Coupling%dtsfcin_cpl  = clear_val
      Coupling%dqsfcin_cpl  = clear_val
      Coupling%ulwsfcin_cpl = clear_val

      !--- accumulated quantities
      allocate (Coupling%dusfc_cpl    (IM))
      allocate (Coupling%dvsfc_cpl    (IM))
      allocate (Coupling%dtsfc_cpl    (IM))
      allocate (Coupling%dqsfc_cpl    (IM))
      allocate (Coupling%dlwsfc_cpl   (IM))
      allocate (Coupling%dswsfc_cpl   (IM))
      allocate (Coupling%dnirbm_cpl   (IM))
      allocate (Coupling%dnirdf_cpl   (IM))
      allocate (Coupling%dvisbm_cpl   (IM))
      allocate (Coupling%dvisdf_cpl   (IM))
      allocate (Coupling%nlwsfc_cpl   (IM))
      allocate (Coupling%nswsfc_cpl   (IM))
      allocate (Coupling%nnirbm_cpl   (IM))
      allocate (Coupling%nnirdf_cpl   (IM))
      allocate (Coupling%nvisbm_cpl   (IM))
      allocate (Coupling%nvisdf_cpl   (IM))

      Coupling%dusfc_cpl    = clear_val
      Coupling%dvsfc_cpl    = clear_val
      Coupling%dtsfc_cpl    = clear_val
      Coupling%dqsfc_cpl    = clear_val
      Coupling%dlwsfc_cpl   = clear_val
      Coupling%dswsfc_cpl   = clear_val
      Coupling%dnirbm_cpl   = clear_val
      Coupling%dnirdf_cpl   = clear_val
      Coupling%dvisbm_cpl   = clear_val
      Coupling%dvisdf_cpl   = clear_val
      Coupling%nlwsfc_cpl   = clear_val
      Coupling%nswsfc_cpl   = clear_val
      Coupling%nnirbm_cpl   = clear_val
      Coupling%nnirdf_cpl   = clear_val
      Coupling%nvisbm_cpl   = clear_val
      Coupling%nvisdf_cpl   = clear_val

      !--- instantaneous quantities
      allocate (Coupling%dusfci_cpl  (IM))
      allocate (Coupling%dvsfci_cpl  (IM))
      allocate (Coupling%dtsfci_cpl  (IM))
      allocate (Coupling%dqsfci_cpl  (IM))
      allocate (Coupling%dlwsfci_cpl (IM))
      allocate (Coupling%dswsfci_cpl (IM))
      allocate (Coupling%dnirbmi_cpl (IM))
      allocate (Coupling%dnirdfi_cpl (IM))
      allocate (Coupling%dvisbmi_cpl (IM))
      allocate (Coupling%dvisdfi_cpl (IM))
      allocate (Coupling%nlwsfci_cpl (IM))
      allocate (Coupling%nswsfci_cpl (IM))
      allocate (Coupling%nnirbmi_cpl (IM))
      allocate (Coupling%nnirdfi_cpl (IM))
      allocate (Coupling%nvisbmi_cpl (IM))
      allocate (Coupling%nvisdfi_cpl (IM))
      allocate (Coupling%t2mi_cpl    (IM))
      allocate (Coupling%q2mi_cpl    (IM))
      allocate (Coupling%u10mi_cpl   (IM))
      allocate (Coupling%v10mi_cpl   (IM))
      allocate (Coupling%tsfci_cpl   (IM))
      allocate (Coupling%psurfi_cpl  (IM))
      allocate (Coupling%oro_cpl     (IM))
      allocate (Coupling%slmsk_cpl   (IM))

      Coupling%dusfci_cpl  = clear_val
      Coupling%dvsfci_cpl  = clear_val
      Coupling%dtsfci_cpl  = clear_val
      Coupling%dqsfci_cpl  = clear_val
      Coupling%dlwsfci_cpl = clear_val
      Coupling%dswsfci_cpl = clear_val
      Coupling%dnirbmi_cpl = clear_val
      Coupling%dnirdfi_cpl = clear_val
      Coupling%dvisbmi_cpl = clear_val
      Coupling%dvisdfi_cpl = clear_val
      Coupling%nlwsfci_cpl = clear_val
      Coupling%nswsfci_cpl = clear_val
      Coupling%nnirbmi_cpl = clear_val
      Coupling%nnirdfi_cpl = clear_val
      Coupling%nvisbmi_cpl = clear_val
      Coupling%nvisdfi_cpl = clear_val
      Coupling%t2mi_cpl    = clear_val
      Coupling%q2mi_cpl    = clear_val
      Coupling%u10mi_cpl   = clear_val
      Coupling%v10mi_cpl   = clear_val
      Coupling%tsfci_cpl   = clear_val
      Coupling%psurfi_cpl  = clear_val
!!    Coupling%oro_cpl     = clear_val  !< pointer to sfcprop%oro
!!    Coupling%slmsk_cpl   = clear_val  !< pointer to sfcprop%slmsk
    endif

    !--- stochastic physics option
    if (Model%do_sppt) then
      allocate (Coupling%sppt_wts  (IM,Model%levs))
      Coupling%sppt_wts = clear_val
    endif

    !--- stochastic shum option
    if (Model%do_shum) then
      allocate (Coupling%shum_wts  (IM,Model%levs))
      Coupling%shum_wts = clear_val
    endif

    !--- stochastic skeb option
    if (Model%do_skeb) then
      allocate (Coupling%skebu_wts (IM,Model%levs))
      allocate (Coupling%skebv_wts (IM,Model%levs))

      Coupling%skebu_wts = clear_val
      Coupling%skebv_wts = clear_val
    endif

    !--- stochastic vc option
    if (Model%do_vc) then
      allocate (Coupling%vcu_wts (IM,Model%levs))
      allocate (Coupling%vcv_wts (IM,Model%levs))

      Coupling%vcu_wts = clear_val
      Coupling%vcv_wts = clear_val
    endif

    !--- needed for either GoCart or 3D diagnostics
    if (Model%lgocart .or. Model%ldiag3d) then
      allocate (Coupling%dqdti   (IM,Model%levs))
      allocate (Coupling%cnvqci  (IM,Model%levs))
      allocate (Coupling%upd_mfi (IM,Model%levs))
      allocate (Coupling%dwn_mfi (IM,Model%levs))
      allocate (Coupling%det_mfi (IM,Model%levs))
      allocate (Coupling%cldcovi (IM,Model%levs))

      Coupling%dqdti    = clear_val
      Coupling%cnvqci   = clear_val
      Coupling%upd_mfi  = clear_val
      Coupling%dwn_mfi  = clear_val
      Coupling%det_mfi  = clear_val
      Coupling%cldcovi  = clear_val

    elseif (.not.Model%uni_cld) then

      allocate (Coupling%cldcovi (IM,Model%levs))
      Coupling%cldcovi  = clear_val

    endif

  end subroutine coupling_create


!----------------------
! GFS_control_type%init
!----------------------
  subroutine control_initialize (Model, nlunit, fn_nml, me, master, &
                                 logunit, isc, jsc, nx, ny, levs,   &
                                 cnx, cny, gnx, gny, dt_dycore,     &
                                 dt_phys, idat, jdat, tracer_names)

    !--- modules
    use physcons,         only: max_lon, max_lat, min_lon, min_lat, &
                                dxmax, dxmin, dxinv, con_rerth, con_pi
    use mersenne_twister, only: random_setseed, random_number
    use module_ras,       only: nrcmax
    use parse_tracers,    only: get_tracer_index
    use wam_f107_kp_mod,  only: f107_kp_size, f107_kp_interval, &
                                f107_kp_skip_size, f107_kp_data_size
    implicit none

    !--- interface variables
    class(GFS_control_type)            :: Model
    integer,                intent(in) :: nlunit
    character(len=64),      intent(in) :: fn_nml
    integer,                intent(in) :: me
    integer,                intent(in) :: master
    integer,                intent(in) :: logunit
    integer,                intent(in) :: isc
    integer,                intent(in) :: jsc
    integer,                intent(in) :: nx
    integer,                intent(in) :: ny
    integer,                intent(in) :: levs
    integer,                intent(in) :: cnx
    integer,                intent(in) :: cny
    integer,                intent(in) :: gnx
    integer,                intent(in) :: gny
    real(kind=kind_phys),   intent(in) :: dt_dycore
    real(kind=kind_phys),   intent(in) :: dt_phys
    integer,                intent(in) :: idat(8)
    integer,                intent(in) :: jdat(8)
    character(len=32),      intent(in) :: tracer_names(:)
    !--- local variables
    integer :: n
    integer :: ios
    integer :: seed0
    logical :: exists
    real(kind=kind_phys) :: tem
    real(kind=kind_phys) :: rinc(5)
    real(kind=kind_evod) :: wrk(1)
    real(kind=kind_phys), parameter :: con_hr = 3600.

    !--- BEGIN NAMELIST VARIABLES
    real(kind=kind_phys) :: fhzero         = 0.0             !< seconds between clearing of diagnostic buckets
    logical              :: ldiag3d        = .false.         !< flag for 3d diagnostic fields
    logical              :: lssav          = .false.         !< logical flag for storing diagnostics
    real(kind=kind_phys) :: fhcyc          = 0.              !< frequency for surface data cycling (secs)
    logical              :: lgocart        = .false.         !< flag for 3d diagnostic fields for gocart 1
    real(kind=kind_phys) :: fhgoc3d        = 0.0             !< seconds between calls to gocart
    integer              :: thermodyn_id   =  1              !< valid for GFS only for get_prs/phi
    integer              :: sfcpress_id    =  1              !< valid for GFS only for get_prs/phi

    !--- coupling parameters
    logical              :: cplflx         = .false.         !< default no cplflx collection
    logical              :: cplwav         = .false.         !< default no cplwav collection

    !--- integrated dynamics through earth's atmosphere
    logical              :: lsidea         = .false.

    !--- radiation parameters
    real(kind=kind_phys) :: fhswr          = 3600.           !< frequency for shortwave radiation (secs)
    real(kind=kind_phys) :: fhlwr          = 3600.           !< frequency for longwave radiation (secs)
    integer              :: levr           = -99             !< number of vertical levels for radiation calculations
    integer              :: nfxr           = 39              !< second dimension of input/output array fluxr   
    logical              :: aero_in        = .false.         !< flag for initializing aero data 
    integer              :: iflip          =  1              !< iflip - is not the same as flipv
    integer              :: isol           =  0              !< use prescribed solar constant
    integer              :: ico2           =  0              !< prescribed global mean value (old opernl)
    integer              :: ialb           =  0              !< use climatology alb, based on sfc type
                                                             !< 1 => use modis based alb
    integer              :: iems           =  0              !< use fixed value of 1.0
    integer              :: iaer           =  1              !< default aerosol effect in sw only
    integer              :: iovr_sw        =  1              !< sw: max-random overlap clouds
    integer              :: iovr_lw        =  1              !< lw: max-random overlap clouds
    integer              :: ictm           =  1              !< ictm=0 => use data at initial cond time, if not
                                                             !<           available; use latest; no extrapolation.
                                                             !< ictm=1 => use data at the forecast time, if not
                                                             !<           available; use latest; do extrapolation.
                                                             !< ictm=yyyy0 => use yyyy data for the forecast time;
                                                             !<           no extrapolation.
                                                             !< ictm=yyyy1 = > use yyyy data for the fcst. If needed, 
                                                             !<           do extrapolation to match the fcst time.
                                                             !< ictm=-1 => use user provided external data for
                                                             !<           the fcst time; no extrapolation.
                                                             !< ictm=-2 => same as ictm=0, but add seasonal cycle
                                                             !<           from climatology; no extrapolation.
    integer              :: isubc_sw       =  0              !< sw clouds without sub-grid approximation
    integer              :: isubc_lw       =  0              !< lw clouds without sub-grid approximation
                                                             !< =1 => sub-grid cloud with prescribed seeds
                                                             !< =2 => sub-grid cloud with randomly generated
                                                             !< seeds
    logical              :: crick_proof    = .false.         !< CRICK-Proof cloud water
    logical              :: ccnorm         = .false.         !< Cloud condensate normalized by cloud cover 
    logical              :: norad_precip   = .false.         !< radiation precip flag for Ferrier/Moorthi
    logical              :: lwhtr          = .true.          !< flag to output lw heating rate (Radtend%lwhc)
    logical              :: swhtr          = .true.          !< flag to output sw heating rate (Radtend%swhc)

    !--- Z-C microphysical parameters
    integer              :: ncld           =  1                 !< cnoice of cloud scheme
    logical              :: zhao_mic       = .false.            !< flag for Zhao-Carr microphysics
    real(kind=kind_phys) :: psautco(2)     = (/6.0d-4,3.0d-4/)  !< [in] auto conversion coeff from ice to snow
    real(kind=kind_phys) :: prautco(2)     = (/1.0d-4,1.0d-4/)  !< [in] auto conversion coeff from cloud to rain
    real(kind=kind_phys) :: evpco          = 2.0d-5             !< [in] coeff for evaporation of largescale rain
    real(kind=kind_phys) :: wminco(2)      = (/1.0d-5,1.0d-5/)  !< [in] water and ice minimum threshold for Zhao

    !--- M-G microphysical parameters
    integer              :: fprcp          =  0                 !< no prognostic rain and snow (MG)
    real(kind=kind_phys) :: mg_dcs         = 350.0              !< Morrison-Gettleman microphysics parameters
    real(kind=kind_phys) :: mg_qcvar       = 2.0
    real(kind=kind_phys) :: mg_ts_auto_ice = 3600.0             !< ice auto conversion time scale

    !--- land/surface model parameters
    integer              :: lsm            =  1              !< flag for land surface model to use =0  for osu lsm; =1  for noah lsm
    integer              :: lsoil          =  4              !< number of soil layers
    integer              :: ivegsrc        =  2              !< ivegsrc = 0   => USGS,
                                                             !< ivegsrc = 1   => IGBP (20 category)
                                                             !< ivegsrc = 2   => UMD  (13 category)
    integer              :: isot           =  0              !< isot = 0   => Zobler soil type  ( 9 category)
                                                             !< isot = 1   => STATSGO soil type (19 category)
    logical              :: mom4ice        = .false.         !< flag controls mom4 sea ice
    logical              :: use_ufo        = .false.         !< flag for gcycle surface option

    !--- tuning parameters for physical parameterizations
    logical              :: ras            = .false.                  !< flag for ras convection scheme
    logical              :: flipv          = .true.                   !< flag for vertical direction flip (ras)
                                                                      !< .true. implies surface at k=1
    logical              :: trans_trac     = .false.                  !< flag for convective transport of tracers (RAS only)
    logical              :: old_monin      = .false.                  !< flag for diff monin schemes
    logical              :: cnvgwd         = .false.                  !< flag for conv gravity wave drag
    logical              :: mstrat         = .false.                  !< flag for moorthi approach for stratus
    logical              :: moist_adj      = .false.                  !< flag for moist convective adjustment
    logical              :: cscnv          = .false.                  !< flag for Chikira-Sugiyama convection
    logical              :: cal_pre        = .false.                  !< flag controls precip type algorithm
    logical              :: do_aw          = .false.                  !< AW scale-aware option in cs convection
    logical              :: do_shoc        = .false.                  !< flag for SHOC
    logical              :: shocaftcnv     = .false.                  !< flag for SHOC
    logical              :: shoc_cld       = .false.                  !< flag for SHOC in grrad
    logical              :: h2o_phys       = .false.                  !< flag for stratosphere h2o
    logical              :: pdfcld         = .false.                  !< flag for pdfcld
    logical              :: shcnvcw        = .false.                  !< flag for shallow convective cloud
    logical              :: redrag         = .false.                  !< flag for reduced drag coeff. over sea
    logical              :: hybedmf        = .false.                  !< flag for hybrid edmf pbl scheme
    logical              :: dspheat        = .false.                  !< flag for tke dissipative heating
    logical              :: cnvcld         = .false.
    logical              :: random_clds    = .false.                  !< flag controls whether clouds are random
    logical              :: shal_cnv       = .false.                  !< flag for calling shallow convection
    integer              :: imfshalcnv     =  1                       !< flag for mass-flux shallow convection scheme
                                                                      !<     1: July 2010 version of mass-flux shallow conv scheme
                                                                      !<         current operational version as of 2016
                                                                      !<     2: scale- & aerosol-aware mass-flux shallow conv scheme (2017)
                                                                      !<     0: modified Tiedtke's eddy-diffusion shallow conv scheme
                                                                      !<    -1: no shallow convection used
    integer              :: imfdeepcnv     =  1                       !< flag for mass-flux deep convection scheme
                                                                      !<     1: July 2010 version of SAS conv scheme
                                                                      !<           current operational version as of 2016
                                                                      !<     2: scale- & aerosol-aware mass-flux deep conv scheme (2017)
    integer              :: nmtvr          = 14                       !< number of topographic variables such as variance etc
                                                                      !< used in the GWD parameterization
    integer              :: jcap           =  1              !< number of spectral wave trancation used only by sascnv shalcnv
    real(kind=kind_phys) :: cs_parm(10) = (/5.0,2.5,1.0e3,3.0e3,20.0,-999.,-999.,0.,0.,0./)
    real(kind=kind_phys) :: flgmin(2)      = (/0.180,0.220/)          !< [in] ice fraction bounds
    real(kind=kind_phys) :: cgwf(2)        = (/0.5d0,0.05d0/)         !< multiplication factor for convective GWD
    real(kind=kind_phys) :: ccwf(2)        = (/1.0d0,1.0d0/)          !< multiplication factor for critical cloud
                                                                      !< workfunction for RAS
    real(kind=kind_phys) :: cdmbgwd(2)     = (/2.0d0,0.25d0/)         !< multiplication factors for cdmb and gwd
    real(kind=kind_phys) :: sup            = 1.1                      !< supersaturation in pdf cloud when t is very low
    real(kind=kind_phys) :: ctei_rm(2)     = (/10.0d0,10.0d0/)        !< critical cloud top entrainment instability criteria 
                                                                      !< (used if mstrat=.true.)
    real(kind=kind_phys) :: crtrh(3)       = (/0.90d0,0.90d0,0.90d0/) !< critical relative humidity at the surface
                                                                      !< PBL top and at the top of the atmosphere
    real(kind=kind_phys) :: dlqf(2)        = (/0.0d0,0.0d0/)          !< factor for cloud condensate detrainment 
                                                                      !< from cloud edges for RAS

    !--- Rayleigh friction
    real(kind=kind_phys) :: prslrd0        = 0.0d0           !< pressure level from which Rayleigh Damping is applied
    real(kind=kind_phys) :: ral_ts         = 0.0d0           !< time scale for Rayleigh damping in days

    !--- near surface temperature model
    logical              :: nst_anl        = .false.         !< flag for NSSTM analysis in gcycle/sfcsub
    integer              :: lsea           = 0 
    real(kind=kind_phys) :: xkzm_m         = 1.0d0           !< [in] bkgd_vdif_m  background vertical diffusion for momentum  
    real(kind=kind_phys) :: xkzm_h         = 1.0d0           !< [in] bkgd_vdif_h  background vertical diffusion for heat q  
    real(kind=kind_phys) :: xkzm_s         = 1.0d0           !< [in] bkgd_vdif_s  sigma threshold for background mom. diffusion  
    integer              :: nstf_name(5)   = (/0,0,1,0,5/)   !< flag 0 for no nst  1 for uncoupled nst  and 2 for coupled NST
                                                             !< nstf_name contains the NSSTM related parameters
                                                             !< nstf_name(1) : 0 = NSSTM off, 1 = NSSTM on but uncoupled
                                                             !<                2 = NSSTM on and coupled
                                                             !< nstf_name(2) : 1 = NSSTM spin up on, 0 = NSSTM spin up off
                                                             !< nstf_name(3) : 1 = NSSTM analysis on, 0 = NSSTM analysis off
                                                             !< nstf_name(4) : zsea1 in mm
                                                             !< nstf_name(5) : zsea2 in mm
     
    !--- stochastic physics options
    real(kind=kind_phys) :: sppt(5)        = -999.           !< stochastic physics tendency amplitude
    real(kind=kind_phys) :: shum(5)        = -999.           !< stochastic boundary layer spf hum amp
    real(kind=kind_phys) :: skeb(5)        = -999.           !< stochastic KE backscatter amplitude
    real(kind=kind_phys) :: vcamp(5)       = -999.           !< stochastic vorticity confinment amp
    real(kind=kind_phys) :: vc             = 0.              !< deterministic vorticity confinement parameter.

    !--- debug flag
    logical              :: debug          = .false.
    logical              :: pre_rad        = .false.         !< flag for testing purpose
    !--- END NAMELIST VARIABLES

    NAMELIST /gfs_physics_nml/                                                              &
                          !--- general parameters
                               fhzero, ldiag3d, lssav, fhcyc, lgocart, fhgoc3d,             &
                               thermodyn_id, sfcpress_id,                                   &
                          !--- coupling parameters
                               cplflx, cplwav, lsidea,                                      &
                          !--- radiation parameters
                               fhswr, fhlwr, levr, nfxr, aero_in, iflip, isol, ico2, ialb,  &
                               isot, iems,  iaer, iovr_sw, iovr_lw, ictm, isubc_sw,         &
                               isubc_lw, crick_proof, ccnorm, lwhtr, swhtr,                 &
                          !--- microphysical parameterizations
                               ncld, zhao_mic, psautco, prautco, evpco, wminco,             &
                               fprcp, mg_dcs, mg_qcvar, mg_ts_auto_ice,                     &
                          !--- land/surface model control
                               lsm, lsoil, nmtvr, ivegsrc, mom4ice, use_ufo,                &
                          !--- physical parameterizations
                               ras, trans_trac, old_monin, cnvgwd, mstrat, moist_adj,       &
                               cscnv, cal_pre, do_aw, do_shoc, shocaftcnv, shoc_cld,        &
                               h2o_phys, pdfcld, shcnvcw, redrag, hybedmf, dspheat, cnvcld, &
                               random_clds, shal_cnv, imfshalcnv, imfdeepcnv, jcap,         &
                               cs_parm, flgmin, cgwf, ccwf, cdmbgwd, sup, ctei_rm, crtrh,   &
                               dlqf,                                                        &
                          !--- Rayleigh friction
                               prslrd0, ral_ts,                                             &
                          !--- near surface temperature model
                               nst_anl, lsea, xkzm_m, xkzm_h, xkzm_s, nstf_name,            &
                          !--- stochastic physics
                               sppt, shum, skeb, vcamp, vc,                                 &
                          !--- debug options
                               debug, pre_rad

    !--- other parameters 
    integer :: nctp    =  0                !< number of cloud types in CS scheme
    logical :: gen_coord_hybrid = .false.  !< for Henry's gen coord

    !--- SHOC parameters
    integer :: nshoc_2d  = 0  !< number of 2d fields for SHOC
    integer :: nshoc_3d  = 0  !< number of 3d fields for SHOC

    !--- convective clouds
    integer :: ncnvcld3d = 0       !< number of convective 3d clouds fields

    !--- stochastic physics control parameters
    logical :: do_sppt   = .false.
    logical :: do_shum   = .false.
    logical :: do_skeb   = .false.
    logical :: do_vc     = .false.

    !--- read in the namelist
    inquire (file=trim(fn_nml), exist=exists)
    if (.not. exists) then
      write(6,*) 'GFS_namelist_read:: namelist file: ',trim(fn_nml),' does not exist'
      stop
    else
      open (unit=nlunit, file=fn_nml, READONLY, status='OLD', iostat=ios)
    endif
    rewind(nlunit)
    read (nlunit, nml=gfs_physics_nml)
    close (nlunit)
    !--- write version number and namelist to log file ---
    if (me == master) write(logunit, nml=gfs_physics_nml)

    !--- MPI parameters
    Model%me               = me
    Model%master           = master
    Model%nlunit           = nlunit
    Model%fn_nml           = fn_nml
    Model%fhzero           = fhzero
    Model%ldiag3d          = ldiag3d
    Model%lssav            = lssav
    Model%fhcyc            = fhcyc
    Model%lgocart          = lgocart
    Model%fhgoc3d          = fhgoc3d
    Model%thermodyn_id     = thermodyn_id
    Model%sfcpress_id      = sfcpress_id
    Model%gen_coord_hybrid = gen_coord_hybrid

    !--- set some grid extent parameters
    Model%isc              = isc
    Model%jsc              = jsc
    Model%nx               = nx
    Model%ny               = ny
    Model%levs             = levs
    Model%cnx              = cnx
    Model%cny              = cny
    Model%lonr             = gnx
    Model%latr             = gny

    !--- coupling parameters
    Model%cplflx           = cplflx
    Model%cplwav           = cplwav

    !--- integrated dynamics through earth's atmosphere
    Model%lsidea           = lsidea

    !--- calendars and time parameters and activation triggers
    Model%dtp              = dt_phys
    Model%dtf              = dt_dycore
    Model%nscyc            = nint(fhcyc*3600./Model%dtp)
    Model%nszero           = nint(Model%fhzero*con_hr/Model%dtp)
    Model%idat(1:8)        = idat(1:8)
    Model%idate            = 0
    Model%idate(1)         = Model%idat(5)
    Model%idate(2)         = Model%idat(2)
    Model%idate(3)         = Model%idat(3)
    Model%idate(4)         = Model%idat(1)

    !--- radiation control parameters
    Model%fhswr            = fhswr
    Model%fhlwr            = fhlwr
    Model%nsswr            = nint(fhswr/Model%dtp)
    Model%nslwr            = nint(fhlwr/Model%dtp)
    if (levr < 0) then
      Model%levr           = levs
    else
      Model%levr           = levr
    endif
    Model%nfxr             = nfxr
    Model%aero_in          = aero_in
    Model%iflip            = iflip
    Model%isol             = isol
    Model%ico2             = ico2
    Model%ialb             = ialb
    Model%iems             = iems
    Model%iaer             = iaer
    Model%iovr_sw          = iovr_sw
    Model%iovr_lw          = iovr_lw
    Model%ictm             = ictm
    Model%isubc_sw         = isubc_sw
    Model%isubc_lw         = isubc_lw
    Model%crick_proof      = crick_proof
    Model%ccnorm           = ccnorm
    Model%lwhtr            = lwhtr
    Model%swhtr            = swhtr

    !--- microphysical switch
    Model%ncld             = ncld
    !--- Zhao-Carr MP parameters
    Model%zhao_mic         = zhao_mic
    Model%psautco          = psautco
    Model%prautco          = prautco
    Model%evpco            = evpco
    Model%wminco           = wminco
    !--- Morroson-Gettleman MP parameters
    Model%fprcp            = fprcp
    Model%mg_dcs           = mg_dcs
    Model%mg_qcvar         = mg_qcvar
    Model%mg_ts_auto_ice   = mg_ts_auto_ice

    !--- land/surface model parameters
    Model%lsm              = lsm
    Model%lsoil            = lsoil
    Model%ivegsrc          = ivegsrc
    Model%isot             = isot
    Model%mom4ice          = mom4ice
    Model%use_ufo          = use_ufo

    !--- tuning parameters for physical parameterizations
    Model%ras              = ras
    Model%flipv            = flipv
    Model%trans_trac       = trans_trac
    Model%old_monin        = old_monin
    Model%cnvgwd           = cnvgwd
    Model%mstrat           = mstrat
    Model%moist_adj        = moist_adj
    Model%cscnv            = cscnv
    Model%cal_pre          = cal_pre
    Model%do_aw            = do_aw
    Model%do_shoc          = do_shoc
    Model%shocaftcnv       = shocaftcnv
    Model%shoc_cld         = shoc_cld
    Model%h2o_phys         = h2o_phys
    Model%pdfcld           = pdfcld
    Model%shcnvcw          = shcnvcw
    Model%redrag           = redrag
    Model%hybedmf          = hybedmf
    Model%dspheat          = dspheat
    Model%cnvcld           = cnvcld
    Model%random_clds      = random_clds
    Model%shal_cnv         = shal_cnv
    Model%imfshalcnv       = imfshalcnv
    Model%imfdeepcnv       = imfdeepcnv
    Model%nmtvr            = nmtvr
    Model%jcap             = jcap
    Model%cs_parm          = cs_parm
    Model%flgmin           = flgmin
    Model%cs_parm          = cs_parm
    Model%cgwf             = cgwf
    Model%ccwf             = ccwf
    Model%cdmbgwd          = cdmbgwd
    Model%sup              = sup
    Model%ctei_rm          = ctei_rm
    Model%crtrh            = crtrh
    Model%dlqf             = dlqf

    !--- Rayleigh friction
    Model%prslrd0          = prslrd0
    Model%ral_ts           = ral_ts

    !--- near surface temperature model
    Model%nst_anl          = nst_anl
    Model%lsea             = lsea
    Model%xkzm_m           = xkzm_m
    Model%xkzm_h           = xkzm_h
    Model%xkzm_s           = xkzm_s
    Model%nstf_name        = nstf_name

    !--- stochastic physics options
    Model%sppt             = sppt
    Model%shum             = shum
    Model%skeb             = skeb
    Model%vcamp            = vcamp
    Model%vc               = vc
    Model%do_sppt          = do_sppt
    Model%do_shum          = do_shum
    Model%do_skeb          = do_skeb
    Model%do_vc            = do_vc

    !--- tracer handling
    Model%ntrac            = size(tracer_names)
    allocate (Model%tracer_names(Model%ntrac))
    Model%tracer_names(:)  = tracer_names(:)
    Model%ntoz             = get_tracer_index(Model%tracer_names, 'o3mr',     Model%me, Model%master, Model%debug)
    Model%ntcw             = get_tracer_index(Model%tracer_names, 'liq_wat',  Model%me, Model%master, Model%debug)
    Model%ntiw             = get_tracer_index(Model%tracer_names, 'ice_wat',  Model%me, Model%master, Model%debug)
    Model%ntrw             = get_tracer_index(Model%tracer_names, 'rainwat',  Model%me, Model%master, Model%debug)
    Model%ntsw             = get_tracer_index(Model%tracer_names, 'snowwat',  Model%me, Model%master, Model%debug)
    Model%ntgl             = get_tracer_index(Model%tracer_names, 'graupel',  Model%me, Model%master, Model%debug)
    Model%ntlnc            = get_tracer_index(Model%tracer_names, 'water_nc', Model%me, Model%master, Model%debug)
    Model%ntinc            = get_tracer_index(Model%tracer_names, 'ice_nc',   Model%me, Model%master, Model%debug)
    Model%ntrnc            = get_tracer_index(Model%tracer_names, 'rain_nc',  Model%me, Model%master, Model%debug)
    Model%ntsnc            = get_tracer_index(Model%tracer_names, 'snow_nc',  Model%me, Model%master, Model%debug)
    Model%ntke             = get_tracer_index(Model%tracer_names, 'sgs_tke',  Model%me, Model%master, Model%debug)

    !--- quantities to be used to derive phy_f*d totals
    Model%nshoc_2d         = nshoc_2d
    Model%nshoc_3d         = nshoc_3d
    Model%ncnvcld3d        = ncnvcld3d
    Model%nctp             = nctp

    !--- debug flag
    Model%debug            = debug
    Model%pre_rad          = pre_rad

    !--- set initial values for time varying properties
    Model%ipt              = 1
    Model%lprnt            = .false.
    Model%lsswr            = .false.
    Model%lslwr            = .false.
    Model%solhr            = -9999.
    Model%solcon           = -9999.
    Model%slag             = -9999.
    Model%sdec             = -9999.
    Model%cdec             = -9999.
    Model%clstp            = -9999
    rinc(1:5)              = 0 
    call w3difdat(jdat,idat,4,rinc)
    Model%phour            = rinc(4)/con_hr
    Model%fhour            = (rinc(4) + Model%dtp)/con_hr
    Model%zhour            = mod(Model%phour,Model%fhzero)
    Model%kdt              = 0
    Model%jdat(1:8)        = jdat(1:8)

    !--- stored in wam_f107_kp module
    f107_kp_size      = 56
    f107_kp_skip_size = 0
    f107_kp_data_size = 56
    f107_kp_interval  = 10800

    !--- BEGIN CODE FROM GFS_PHYSICS_INITIALIZE
    !--- define physcons module variables
    tem   = con_rerth*con_rerth*(con_pi+con_pi)*con_pi
    dxmax = log(tem/(max_lon*max_lat))
    dxmin = log(tem/(min_lon*min_lat))
    dxinv = 1.0d0 / (dxmax-dxmin)
    if (Model%me == Model%master) write(0,*)' dxmax=',dxmax,' dxmin=',dxmin,' dxinv=',dxinv

    !--- set nrcm 
    if (Model%ras) then
      Model%nrcm = min(nrcmax, Model%levs-1) * (Model%dtp/1200.d0) + 0.10001d0
    else
      Model%nrcm = 2
    endif

    !--- cal_pre
    if (Model%cal_pre) then
      Model%random_clds = .true.
    endif
    !--- END CODE FROM GFS_PHYSICS_INITIALIZE


    !--- BEGIN CODE FROM COMPNS_PHYSICS
    !--- shoc scheme
    if (do_shoc) then
      Model%nshoc_3d   = 3
      Model%nshoc_2d   = 0
      Model%shal_cnv   = .false.
      Model%imfshalcnv = -1
      Model%hybedmf    = .false.
      if (Model%me == Model%master) print *,' Simplified Higher Order Closure Model used for', &
                                            ' Boundary layer and Shallow Convection',          &
                                            ' nshoc_3d=',Model%nshoc_3d,                       &
                                            ' nshoc_2d=',Model%nshoc_2d,                       &
                                            ' ntke=',Model%ntke
    endif

    !--- set number of cloud types
    if (Model%cscnv) then
      Model%nctp = nint(Model%cs_parm(5))
      Model%nctp = max(Model%nctp,10)
      if (Model%cs_parm(7) < 0.0) Model%cs_parm(7) = Model%dtp
    endif
    Model%nctp = max(Model%nctp,1)

    !--- output information about the run
    if (Model%me == Model%master) then
      if (Model%lsm == 1) then
        print *,' NOAH Land Surface Model used'
      elseif (Model%lsm == 0) then
        print *,' OSU no longer supported - job aborted'
        stop
      else
        print *,' Unsupported LSM type - job aborted - lsm=',Model%lsm
        stop
      endif
      print *,' nst_anl=',Model%nst_anl,' use_ufo=',Model%use_ufo
      if (Model%nstf_name(1) > 0 ) then 
        print *,' NSSTM is active '
        print *,' nstf_name(1)=',Model%nstf_name(1)
        print *,' nstf_name(2)=',Model%nstf_name(2)
        print *,' nstf_name(3)=',Model%nstf_name(3)
        print *,' nstf_name(4)=',Model%nstf_name(4)
        print *,' nstf_name(5)=',Model%nstf_name(5)
      endif
      if (.not. Model%cscnv) then
        if (Model%ras) then
          print *,' RAS Convection scheme used with ccwf=',Model%ccwf
          Model%imfdeepcnv = -1
        else
          if (Model%imfdeepcnv == 0) then
            print *,' old SAS Convection scheme before July 2010 used'
          elseif(Model%imfdeepcnv == 1) then
            print *,' July 2010 version of SAS conv scheme used'
          elseif(Model%imfdeepcnv == 2) then
          print *,' scale & aerosol-aware mass-flux deep conv scheme'
          endif
        endif
      else
        if (Model%do_aw) then
          print *,'Chikira-Sugiyama convection scheme with Arakawa-Wu'&
     &,              ' unified parameterization used'
        else
            print *,'Chikira-Sugiyama convection scheme used'
        endif
        print *,' cs_parm=',Model%cs_parm,' nctp=',Model%nctp
      endif
      if (.not. Model%old_monin .and. .not. Model%do_shoc) print *,' New PBL scheme used'
      if (.not. Model%shal_cnv) then
        Model%imfshalcnv = -1
        print *,' No shallow convection used'
      else
        if (Model%imfshalcnv == 0) then
          print *,' modified Tiedtke eddy-diffusion shallow conv scheme used'
        elseif (Model%imfshalcnv == 1) then
          print *,' July 2010 version of mass-flux shallow conv scheme used'
        elseif (Model%imfshalcnv == 2) then
          print *,' scale- & aerosol-aware mass-flux shallow conv scheme (2017)'
        else
          print *,' unknown mass-flux scheme in use - defaulting to no shallow convection'
          Model%imfshalcnv = -1
        endif
      endif
      if (Model%cnvgwd)      print *,' Convective GWD parameterization used'
      if (Model%crick_proof) print *,' CRICK-Proof cloud water used in radiation '
      if (Model%ccnorm)      print *,' Cloud condensate normalized by cloud cover for radiation'

      print *,' Radiative heating calculated at',Model%levr, ' layers'
      if (Model%iovr_sw == 0) then
        print *,' random cloud overlap for Shortwave IOVR_SW=',Model%iovr_sw
      else
        print *,' max-random cloud overlap for Shortwave IOVR_SW=',Model%iovr_sw
      endif
      if (Model%iovr_lw == 0) then
        print *,' random cloud overlap for Longwave IOVR_LW=',Model%iovr_lw
      else
        print *,' max-random cloud overlap for Longwave IOVR_LW=',Model%iovr_lw
      endif
      if (Model%isubc_sw == 0) then
        print *,' no sub-grid cloud for Shortwave ISUBC_SW=',Model%isubc_sw
      else
        print *,' sub-grid cloud for Shortwave ISUBC_SW=',Model%isubc_sw
      endif
      if (Model%isubc_lw == 0) then
        print *,' no sub-grid cloud for Longwave ISUBC_LW=',Model%isubc_lw
      else
        print *,' sub-grid cloud for Longwave ISUBC_LW=',Model%isubc_lw
      endif
    endif

    !--- set up cloud schemes and tracer elements
    if (Model%ncld <= 1) then
      if (Model%zhao_mic) then        ! default setup for Zhao Microphysics
        Model%num_p3d = 4
        Model%num_p2d = 3
        if (Model%pdfcld) then
          Model%npdf3d = 3
        else
          Model%shcnvcw = .false.
        endif
        if (Model%me == Model%master) print *,' Using Zhao/Carr/Sundqvist Microphysics'
      else
        print *,' Ferrier Microphysics scheme has been deprecated - job aborted'
        stop
      endif
    elseif (Model%ncld == 2) then
      Model%num_p3d = 1
      Model%num_p2d = 1
      Model%pdfcld  = .false.
      Model%shcnvcw = .false.
      if (Model%me == Model%master) print *,' Using Morrison-Gettelman double moment', &
                                            ' microphysics',' aero_in=',Model%aero_in, &
                                            ' mg_dcs=',Model%mg_dcs,' mg_qcvar=',Model%mg_qcvar, &
                                            ' mg_ts_auto_ice=',Model%mg_ts_auto_ice
    endif

    Model%uni_cld = .false.
    if ((Model%shoc_cld) .or. (Model%ncld == 2)) then
      Model%uni_cld = .true.
    endif

    if(Model%ras .or. Model%cscnv) Model%cnvcld = .false.
    if(Model%do_shoc .or. Model%pdfcld) Model%cnvcld = .false.
    if(Model%cnvcld) Model%ncnvcld3d = 1
 
    !--- derived totals for phy_f*d
    Model%ntot2d = Model%num_p2d + Model%nshoc_2d
    Model%ntot3d = Model%num_p3d + Model%nshoc_3d + Model%npdf3d + Model%ncnvcld3d
    if (me == Model%master) print *,' num_p3d=',Model%num_p3d,' num_p2d=',Model%num_p2d,    &
                                    ' crtrh=',Model%crtrh,' npdf3d=',Model%npdf3d,          &
                                    ' pdfcld=',Model%pdfcld,' shcnvcw=',Model%shcnvcw,      &
                                    ' cnvcld=',Model%cnvcld,' ncnvcld3d=',Model%ncnvcld3d,  &
                                    ' do_shoc=',Model%do_shoc,' nshoc3d=',Model%nshoc_3d,   &
                                    ' nshoc_2d=',Model%nshoc_2d,' shoc_cld=',Model%shoc_cld,& 
                                    ' ntot3d=',Model%ntot3d,' ntot2d=',Model%ntot2d,        &
                                    ' shocaftcnv=',Model%shocaftcnv

    !--- stochastic physics
    if (Model%sppt(1) > 0 ) Model%do_sppt = .true.
    if (Model%shum(1) > 0 ) Model%do_shum = .true.
    if (Model%skeb(1) > 0 ) Model%do_skeb = .true.
    if (Model%vc > tiny(Model%vc) .or. Model%vcamp(1) > 0 ) Model%do_vc = .true.
    if (Model%me == Model%master) write(0,*)' in compns_physics do_sppt=',Model%do_sppt,         &
                                            ' do_shum=',Model%do_shum,' do_skeb=',Model%do_skeb, &
                                            ' do_vc=',Model%do_vc,' nctp=',Model%nctp
    !--- END CODE FROM COMPNS_PHYSICS


    !--- BEGIN CODE FROM GLOOPR
    !--- set up parameters for Xu & Randell's cloudiness computation (Radiation)
    Model%lmfshal  = (Model%shal_cnv .and. (Model%imfshalcnv > 0))
    Model%lmfdeep2 = (Model%imfdeepcnv == 2)
    !--- END CODE FROM GLOOPR

    !--- BEGIN CODE FROM GLOOPB
    !--- set up random number seed needed for RAS and old SAS and when cal_pre=.true.
    if ((Model%imfdeepcnv <= 0) .or. (Model%cal_pre)) then
      if (Model%random_clds) then
        seed0 = Model%idate(1) + Model%idate(2) + Model%idate(3) + Model%idate(4)
        call random_setseed(seed0)
        call random_number(wrk)
        Model%seed0 = seed0 + nint(wrk(1)*1000.0d0)
      endif
    endif
    !--- END CODE FROM GLOOPB

    call Model%print ()

  end subroutine control_initialize


!------------------
! GFS_control%print
!------------------
  subroutine control_print(Model)

    implicit none

    !--- interface variables
    class(GFS_control_type) :: Model
 
    if (Model%me == Model%master) then
      print *, ' '
      print *, 'basic control parameters'
      print *, ' me                : ', Model%me
      print *, ' master            : ', Model%master
      print *, ' nlunit            : ', Model%nlunit
      print *, ' fn_nml            : ', trim(Model%fn_nml)
      print *, ' fhzero            : ', Model%fhzero
      print *, ' ldiag3d           : ', Model%ldiag3d
      print *, ' lssav             : ', Model%lssav
      print *, ' fhcyc             : ', Model%fhcyc
      print *, ' lgocart           : ', Model%lgocart
      print *, ' fhgoc3d           : ', Model%fhgoc3d
      print *, ' thermodyn_id      : ', Model%thermodyn_id
      print *, ' sfcpress_id       : ', Model%sfcpress_id
      print *, ' gen_coord_hybrid  : ', Model%gen_coord_hybrid
      print *, ' '
      print *, 'grid extent parameters'
      print *, ' isc               : ', Model%isc
      print *, ' jsc               : ', Model%jsc
      print *, ' nx                : ', Model%nx
      print *, ' ny                : ', Model%ny
      print *, ' levs              : ', Model%levs
      print *, ' cnx               : ', Model%cnx
      print *, ' cny               : ', Model%cny
      print *, ' lonr              : ', Model%lonr
      print *, ' latr              : ', Model%latr
      print *, ' '
      print *, 'coupling parameters'
      print *, ' cplflx            : ', Model%cplflx
      print *, ' cplwav            : ', Model%cplwav
      print *, ' '
      print *, 'integrated dynamics through earth atmosphere'
      print *, ' lsidea            : ', Model%lsidea
      print *, ' '
      print *, 'calendars and time parameters and activation triggers'
      print *, ' dtp               : ', Model%dtp
      print *, ' dtf               : ', Model%dtf
      print *, ' nscyc             : ', Model%nscyc
      print *, ' nszero            : ', Model%nszero
      print *, ' idat              : ', Model%idat
      print *, ' idate             : ', Model%idate
      print *, ' '
      print *, 'radiation control parameters'
      print *, ' fhswr             : ', Model%fhswr
      print *, ' fhlwr             : ', Model%fhlwr
      print *, ' nsswr             : ', Model%nsswr
      print *, ' nslwr             : ', Model%nslwr
      print *, ' levr              : ', Model%levr
      print *, ' nfxr              : ', Model%nfxr
      print *, ' aero_in           : ', Model%aero_in
      print *, ' lmfshal           : ', Model%lmfshal
      print *, ' lmfdeep2          : ', Model%lmfdeep2
      print *, ' nrcm              : ', Model%nrcm
      print *, ' iflip             : ', Model%iflip
      print *, ' isol              : ', Model%isol
      print *, ' ico2              : ', Model%ico2
      print *, ' ialb              : ', Model%ialb
      print *, ' iems              : ', Model%iems
      print *, ' iaer              : ', Model%iaer
      print *, ' iovr_sw           : ', Model%iovr_sw
      print *, ' iovr_lw           : ', Model%iovr_lw
      print *, ' ictm              : ', Model%ictm
      print *, ' isubc_sw          : ', Model%isubc_sw
      print *, ' isubc_lw          : ', Model%isubc_lw
      print *, ' crick_proof       : ', Model%crick_proof
      print *, ' ccnorm            : ', Model%ccnorm
      print *, ' norad_precip      : ', Model%norad_precip
      print *, ' lwhtr             : ', Model%lwhtr
      print *, ' swhtr             : ', Model%swhtr
      print *, ' '
      print *, 'microphysical switch'
      print *, ' ncld              : ', Model%ncld
      print *, ' Z-C microphysical parameters'
      print *, ' zhao_mic          : ', Model%zhao_mic
      print *, ' psautco           : ', Model%psautco
      print *, ' prautco           : ', Model%prautco
      print *, ' evpco             : ', Model%evpco
      print *, ' wminco            : ', Model%wminco
      print *, ' M-G microphysical parameters'
      print *, ' fprcp             : ', Model%fprcp
      print *, ' mg_dcs            : ', Model%mg_dcs
      print *, ' mg_qcvar          : ', Model%mg_qcvar
      print *, ' mg_ts_auto_ice    : ', Model%mg_ts_auto_ice
      print *, ' '
      print *, 'land/surface model parameters'
      print *, ' lsm               : ', Model%lsm
      print *, ' lsoil             : ', Model%lsoil
      print *, ' ivegsrc           : ', Model%ivegsrc
      print *, ' isot              : ', Model%isot
      print *, ' mom4ice           : ', Model%mom4ice
      print *, ' use_ufo           : ', Model%use_ufo
      print *, ' '
      print *, 'tuning parameters for physical parameterizations'
      print *, ' ras               : ', Model%ras
      print *, ' flipv             : ', Model%flipv
      print *, ' trans_trac        : ', Model%trans_trac
      print *, ' old_monin         : ', Model%old_monin
      print *, ' cnvgwd            : ', Model%cnvgwd
      print *, ' mstrat            : ', Model%mstrat
      print *, ' moist_adj         : ', Model%moist_adj
      print *, ' cscnv             : ', Model%cscnv
      print *, ' cal_pre           : ', Model%cal_pre
      print *, ' do_aw             : ', Model%do_aw
      print *, ' do_shoc           : ', Model%do_shoc
      print *, ' shocaftcnv        : ', Model%shocaftcnv
      print *, ' shoc_cld          : ', Model%shoc_cld
      print *, ' uni_cld           : ', Model%uni_cld
      print *, ' h2o_phys          : ', Model%h2o_phys
      print *, ' pdfcld            : ', Model%pdfcld
      print *, ' shcnvcw           : ', Model%shcnvcw
      print *, ' redrag            : ', Model%redrag
      print *, ' hybedmf           : ', Model%hybedmf
      print *, ' dspheat           : ', Model%dspheat
      print *, ' cnvcld            : ', Model%cnvcld
      print *, ' random_clds       : ', Model%random_clds
      print *, ' shal_cnv          : ', Model%shal_cnv
      print *, ' imfshalcnv        : ', Model%imfshalcnv
      print *, ' imfdeepcnv        : ', Model%imfdeepcnv
      print *, ' nmtvr             : ', Model%nmtvr
      print *, ' jcap              : ', Model%jcap
      print *, ' cs_parm           : ', Model%cs_parm
      print *, ' flgmin            : ', Model%flgmin
      print *, ' cgwf              : ', Model%cgwf
      print *, ' ccwf              : ', Model%ccwf
      print *, ' cdmbgwd           : ', Model%cdmbgwd
      print *, ' sup               : ', Model%sup
      print *, ' ctei_rm           : ', Model%ctei_rm
      print *, ' crtrh             : ', Model%crtrh
      print *, ' dlqf              : ', Model%dlqf
      print *, ' seed0             : ', Model%seed0
      print *, ' '
      print *, 'Rayleigh friction'
      print *, ' prslrd0           : ', Model%prslrd0
      print *, ' ral_ts            : ', Model%ral_ts
      print *, ' '
      print *, 'near surface temperature model'
      print *, ' nst_anl           : ', Model%nst_anl
      print *, ' lsea              : ', Model%lsea
      print *, ' xkzm_m            : ', Model%xkzm_m
      print *, ' xkzm_h            : ', Model%xkzm_h
      print *, ' xkzm_s            : ', Model%xkzm_s
      print *, ' nstf_name         : ', Model%nstf_name
      print *, ' '
      print *, 'stochastic physics'
      print *, ' do_sppt           : ', Model%do_sppt
      print *, ' do_shum           : ', Model%do_shum
      print *, ' do_skeb           : ', Model%do_skeb
      print *, ' do_vc             : ', Model%do_vc
      print *, ' sppt              : ', Model%sppt
      print *, ' shum              : ', Model%shum
      print *, ' skeb              : ', Model%skeb
      print *, ' vcamp             : ', Model%vcamp
      print *, ' vc                : ', Model%vc
      print *, ' '
      print *, 'tracers'
      print *, ' tracer_names      : ', Model%tracer_names
      print *, ' ntrac             : ', Model%ntrac
      print *, ' ntoz              : ', Model%ntoz
      print *, ' ntcw              : ', Model%ntcw
      print *, ' ntiw              : ', Model%ntiw
      print *, ' ntrw              : ', Model%ntrw
      print *, ' ntsw              : ', Model%ntsw
      print *, ' ntgl              : ', Model%ntgl
      print *, ' ntlnc             : ', Model%ntlnc
      print *, ' ntinc             : ', Model%ntinc
      print *, ' ntrnc             : ', Model%ntrnc
      print *, ' ntsnc             : ', Model%ntsnc
      print *, ' ntke              : ', Model%ntke
      print *, ' nto               : ', Model%nto
      print *, ' nto2              : ', Model%nto2
      print *, ' '
      print *, 'derived totals for phy_f*d'
      print *, ' ntot2d            : ', Model%ntot2d
      print *, ' ntot3d            : ', Model%ntot3d
      print *, ' num_p2d           : ', Model%num_p2d
      print *, ' num_p3d           : ', Model%num_p3d
      print *, ' nshoc_2d          : ', Model%nshoc_2d
      print *, ' nshoc_3d          : ', Model%nshoc_3d
      print *, ' ncnvcld3d         : ', Model%ncnvcld3d
      print *, ' npdf3d            : ', Model%npdf3d
      print *, ' nctp              : ', Model%nctp
      print *, ' '
      print *, 'debug flags'
      print *, ' debug             : ', Model%debug 
      print *, ' pre_rad           : ', Model%pre_rad
      print *, ' '
      print *, 'variables modified at each time step'
      print *, ' ipt               : ', Model%ipt
      print *, ' lprnt             : ', Model%lprnt
      print *, ' lsswr             : ', Model%lsswr
      print *, ' lslwr             : ', Model%lslwr
      print *, ' solhr             : ', Model%solhr
      print *, ' solcon            : ', Model%solcon
      print *, ' slag              : ', Model%slag
      print *, ' sdec              : ', Model%sdec
      print *, ' cdec              : ', Model%cdec
      print *, ' clstp             : ', Model%clstp
      print *, ' phour             : ', Model%phour
      print *, ' fhour             : ', Model%fhour
      print *, ' zhour             : ', Model%zhour
      print *, ' kdt               : ', Model%kdt
      print *, ' jdat              : ', Model%jdat
    endif

  end subroutine control_print


!----------------
! GFS_grid%create
!----------------
  subroutine grid_create (Grid, IM, Model)

    implicit none

    class(GFS_grid_type)              :: Grid
    integer,                intent(in) :: IM
    type(GFS_control_type), intent(in) :: Model

    allocate (Grid%xlon   (IM))
    allocate (Grid%xlat   (IM))
    allocate (Grid%xlat_d (IM))
    allocate (Grid%sinlat (IM))
    allocate (Grid%coslat (IM))
    allocate (Grid%area   (IM))
    allocate (Grid%dx     (IM))

    Grid%xlon   = clear_val
    Grid%xlat   = clear_val
    Grid%xlat_d = clear_val
    Grid%sinlat = clear_val
    Grid%coslat = clear_val
    Grid%area   = clear_val
    Grid%dx     = clear_val

    !--- ozone active
    if ( Model%ntoz > 0 ) then
      allocate (Grid%ddy_o3    (IM))
      allocate (Grid%jindx1_o3 (IM))
      allocate (Grid%jindx2_o3 (IM))
    endif

    !--- stratosphere h2o active
    if ( Model%h2o_phys ) then
      allocate (Grid%ddy_h    (IM))
      allocate (Grid%jindx1_h (IM))
      allocate (Grid%jindx2_h (IM))
    endif
 end subroutine grid_create


!--------------------
! GFS_tbd_type%create
!--------------------
  subroutine tbd_create (Tbd, IM, BLKSZ, BLKNO, Model)

    implicit none

    class(GFS_tbd_type)                :: Tbd
    integer,                intent(in) :: IM
    integer, dimension(:),  intent(in) :: BLKSZ
    integer,                intent(in) :: BLKNO
    type(GFS_control_type), intent(in) :: Model

    !--- In
    !--- sub-grid cloud radiation
    if ( Model%isubc_lw == 2 .or. Model%isubc_sw == 2 ) then
      allocate (Tbd%icsdsw (IM))
      allocate (Tbd%icsdlw (IM))
    endif

    !--- ozone and stratosphere h2o needs
    allocate (Tbd%ozpl  (IM,levozp,oz_coeff))
    allocate (Tbd%h2opl (IM,levh2o,h2o_coeff))
    Tbd%ozpl  = clear_val
    Tbd%h2opl = clear_val

    allocate (Tbd%rann (IM,Model%nrcm))
    Tbd%rann = rann_init

    !--- In/Out
    allocate (Tbd%acv  (IM))
    allocate (Tbd%acvb (IM))
    allocate (Tbd%acvt (IM))

    Tbd%acv  = clear_val
    Tbd%acvb = clear_val
    Tbd%acvt = clear_val

    if (Model%do_sppt) then
      allocate (Tbd%dtdtr     (IM,Model%levs))
      allocate (Tbd%dtotprcp  (IM))
      allocate (Tbd%dcnvprcp  (IM))
      allocate (Tbd%drain_cpl (IM))
      allocate (Tbd%dsnow_cpl (IM))

      Tbd%dtdtr     = clear_val
      Tbd%dtotprcp  = clear_val
      Tbd%dcnvprcp  = clear_val
      Tbd%drain_cpl = clear_val
      Tbd%dsnow_cpl = clear_val
    endif

    allocate (Tbd%phy_fctd (IM,Model%nctp))
    allocate (Tbd%phy_f2d  (IM,Model%ntot2d))
    allocate (Tbd%phy_f3d  (IM,Model%levs,Model%ntot3d))

    allocate (Tbd%blksz (size(BLKSZ)))

    Tbd%phy_fctd = clear_val
    Tbd%phy_f2d  = clear_val
    Tbd%phy_f3d  = clear_val

    Tbd%blkno = BLKNO
    Tbd%blksz = BLKSZ

  end subroutine tbd_create


!------------------------
! GFS_cldprop_type%create
!------------------------
  subroutine cldprop_create (Cldprop, IM, Model)

    implicit none

    class(GFS_cldprop_type)            :: Cldprop
    integer,                intent(in) :: IM
    type(GFS_control_type), intent(in) :: Model

    allocate (Cldprop%cv  (IM))
    allocate (Cldprop%cvt (IM)) 
    allocate (Cldprop%cvb (IM))
    
    Cldprop%cv  = clear_val
    Cldprop%cvt = clear_val
    Cldprop%cvb = clear_val
  
  end subroutine cldprop_create


!******************************************
! GFS_radtend_type%create
!******************************************
  subroutine radtend_create (Radtend, IM, Model)
                               
    implicit none
       
    class(GFS_radtend_type)            :: Radtend
    integer,                intent(in) :: IM
    type(GFS_control_type), intent(in) :: Model

    !--- Out (radiation only) 
    allocate (Radtend%sfcfsw (IM))
    allocate (Radtend%sfcflw (IM))

    Radtend%sfcfsw%upfxc = clear_val
    Radtend%sfcfsw%upfx0 = clear_val
    Radtend%sfcfsw%dnfxc = clear_val
    Radtend%sfcfsw%dnfx0 = clear_val
    Radtend%sfcflw%upfxc = clear_val
    Radtend%sfcflw%upfx0 = clear_val
    Radtend%sfcflw%dnfxc = clear_val
    Radtend%sfcflw%dnfx0 = clear_val
         
    allocate (Radtend%htrsw  (IM,Model%levs))
    allocate (Radtend%htrlw  (IM,Model%levs))
    allocate (Radtend%sfalb  (IM))
    allocate (Radtend%coszen (IM))
    allocate (Radtend%tsflw  (IM))
    allocate (Radtend%semis  (IM))

    Radtend%htrsw  = clear_val
    Radtend%htrlw  = clear_val
    Radtend%sfalb  = clear_val
    Radtend%coszen = clear_val
    Radtend%tsflw  = clear_val
    Radtend%semis  = clear_val
             
    !--- In/Out (???) (radiation only)
    allocate (Radtend%coszdg (IM))

    Radtend%coszdg = clear_val
             
    !--- In/Out (???) (physics only)
    allocate (Radtend%swhc  (IM,Model%levs))
    allocate (Radtend%lwhc  (IM,Model%levs))
    allocate (Radtend%lwhd  (IM,Model%levs,6))

    Radtend%lwhd  = clear_val
    Radtend%lwhc  = clear_val
    Radtend%swhc  = clear_val

  end subroutine radtend_create


!----------------
! GFS_diag%create
!----------------
  subroutine diag_create (Diag, IM, Model)
    class(GFS_diag_type)               :: Diag
    integer,                intent(in) :: IM
    type(GFS_control_type), intent(in) :: Model

    !--- Radiation
    allocate (Diag%fluxr   (IM,Model%nfxr))
    allocate (Diag%topfsw  (IM))
    allocate (Diag%topflw  (IM))
    !--- Physics
    !--- In/Out
    allocate (Diag%srunoff (IM))
    allocate (Diag%evbsa   (IM))
    allocate (Diag%evcwa   (IM))
    allocate (Diag%snohfa  (IM))
    allocate (Diag%transa  (IM))
    allocate (Diag%sbsnoa  (IM))
    allocate (Diag%snowca  (IM))
    allocate (Diag%soilm   (IM))
    allocate (Diag%tmpmin  (IM))
    allocate (Diag%tmpmax  (IM))
    allocate (Diag%dusfc   (IM))
    allocate (Diag%dvsfc   (IM))
    allocate (Diag%dtsfc   (IM))
    allocate (Diag%dqsfc   (IM))
    allocate (Diag%totprcp (IM))
    allocate (Diag%gflux   (IM))
    allocate (Diag%dlwsfc  (IM))
    allocate (Diag%ulwsfc  (IM))
    allocate (Diag%suntim  (IM))
    allocate (Diag%runoff  (IM))
    allocate (Diag%ep      (IM))
    allocate (Diag%cldwrk  (IM))
    allocate (Diag%dugwd   (IM))
    allocate (Diag%dvgwd   (IM))
    allocate (Diag%psmean  (IM))
    allocate (Diag%cnvprcp (IM))
    allocate (Diag%spfhmin (IM))
    allocate (Diag%spfhmax (IM))
    allocate (Diag%rain    (IM))
    allocate (Diag%rainc   (IM))
    allocate (Diag%ice     (IM))
    allocate (Diag%snow    (IM))
    allocate (Diag%graupel (IM))
    allocate (Diag%totice  (IM))
    allocate (Diag%totsnw  (IM))
    allocate (Diag%totgrp  (IM))
    allocate (Diag%u10m    (IM))
    allocate (Diag%v10m    (IM))
    allocate (Diag%zlvl    (IM))
    allocate (Diag%psurf   (IM))
    allocate (Diag%hpbl    (IM))
    allocate (Diag%pwat    (IM))
    allocate (Diag%t1      (IM))
    allocate (Diag%q1      (IM))
    allocate (Diag%u1      (IM))
    allocate (Diag%v1      (IM))
    allocate (Diag%chh     (IM))
    allocate (Diag%cmm     (IM))
    allocate (Diag%dlwsfci (IM))
    allocate (Diag%ulwsfci (IM))
    allocate (Diag%dswsfci (IM))
    allocate (Diag%uswsfci (IM))
    allocate (Diag%dusfci  (IM))
    allocate (Diag%dvsfci  (IM))
    allocate (Diag%dtsfci  (IM))
    allocate (Diag%dqsfci  (IM))
    allocate (Diag%gfluxi  (IM))
    allocate (Diag%epi     (IM))
    allocate (Diag%smcwlt2 (IM))
    allocate (Diag%smcref2 (IM))
    allocate (Diag%wet1    (IM))
    allocate (Diag%sr      (IM))
    !--- 3D diagnostics
    if (Model%ldiag3d) then
      allocate (Diag%du3dt  (IM,Model%levs,4))
      allocate (Diag%dv3dt  (IM,Model%levs,4))
      allocate (Diag%dt3dt  (IM,Model%levs,6))
      allocate (Diag%dq3dt  (IM,Model%levs,oz_coeff+5))
      !--- needed to allocate GoCart coupling fields
      allocate (Diag%upd_mf (IM,Model%levs))
      allocate (Diag%dwn_mf (IM,Model%levs))
      allocate (Diag%det_mf (IM,Model%levs))
      allocate (Diag%cldcov (IM,Model%levs))
    endif

    call Diag%rad_zero  (Model)
    call Diag%phys_zero (Model)

  end subroutine diag_create

!-----------------------
! GFS_diag%rad_zero
!-----------------------
  subroutine diag_rad_zero(Diag, Model)
    class(GFS_diag_type)               :: Diag
    type(GFS_control_type), intent(in) :: Model

    Diag%fluxr        = zero
    Diag%topfsw%upfxc = zero
    Diag%topfsw%dnfxc = zero
    Diag%topfsw%upfx0 = zero
    Diag%topflw%upfxc = zero
    Diag%topflw%upfx0 = zero
    if (Model%ldiag3d) then
      Diag%cldcov     = zero
    endif


  end subroutine diag_rad_zero

!------------------------
! GFS_diag%phys_zero
!------------------------
  subroutine diag_phys_zero (Diag, Model)
    class(GFS_diag_type)               :: Diag
    type(GFS_control_type), intent(in) :: Model

    !--- In/Out
    Diag%srunoff = zero
    Diag%evbsa   = zero
    Diag%evcwa   = zero
    Diag%snohfa  = zero
    Diag%transa  = zero
    Diag%sbsnoa  = zero
    Diag%snowca  = zero
    Diag%soilm   = zero
    Diag%tmpmin  = huge
    Diag%tmpmax  = zero
    Diag%dusfc   = zero
    Diag%dvsfc   = zero
    Diag%dtsfc   = zero
    Diag%dqsfc   = zero
    Diag%totprcp = zero
    Diag%gflux   = zero
    Diag%dlwsfc  = zero
    Diag%ulwsfc  = zero
    Diag%suntim  = zero
    Diag%runoff  = zero
    Diag%ep      = zero
    Diag%cldwrk  = zero
    Diag%dugwd   = zero
    Diag%dvgwd   = zero
    Diag%psmean  = zero
    Diag%cnvprcp = zero
    Diag%spfhmin = huge
    Diag%spfhmax = zero
    Diag%rain    = zero
    Diag%rainc   = zero
    Diag%ice     = zero
    Diag%snow    = zero
    Diag%graupel = zero
    Diag%totice  = zero
    Diag%totsnw  = zero
    Diag%totgrp  = zero

    !--- Out
    Diag%u10m    = zero
    Diag%v10m    = zero
    Diag%zlvl    = zero
    Diag%psurf   = zero
    Diag%hpbl    = zero
    Diag%pwat    = zero
    Diag%t1      = zero
    Diag%q1      = zero
    Diag%u1      = zero
    Diag%v1      = zero
    Diag%chh     = zero
    Diag%cmm     = zero
    Diag%dlwsfci = zero
    Diag%ulwsfci = zero
    Diag%dswsfci = zero
    Diag%uswsfci = zero
    Diag%dusfci  = zero
    Diag%dvsfci  = zero
    Diag%dtsfci  = zero
    Diag%dqsfci  = zero
    Diag%gfluxi  = zero
    Diag%epi     = zero
    Diag%smcwlt2 = zero
    Diag%smcref2 = zero
    Diag%wet1    = zero
    Diag%sr      = zero

    if (Model%ldiag3d) then
      Diag%du3dt   = zero
      Diag%dv3dt   = zero
      Diag%dt3dt   = zero
      Diag%dq3dt   = zero
      Diag%upd_mf  = zero
      Diag%dwn_mf  = zero
      Diag%det_mf  = zero
    endif

  end subroutine diag_phys_zero

end module GFS_typedefs

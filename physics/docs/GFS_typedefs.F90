!>\file GFS_typedefs.F90
!! This file contains DDTs on the host model side.

#undef MULTI_GASES

!>\defgroup fv3gfs_typedefs GFS typedefs
module GFS_typedefs

       use machine,                  only: kind_phys
#ifdef CCPP
       use physcons,                 only: con_cp, con_fvirt, con_g, &
                                           con_hvap, con_hfus, con_pi, con_rd, con_rv, &
                                           con_t0c, con_cvap, con_cliq, con_eps, &
                                           con_epsm1, con_ttp
       use module_radsw_parameters,  only: topfsw_type, sfcfsw_type, cmpfsw_type, NBDSW
       use module_radlw_parameters,  only: topflw_type, sfcflw_type, NBDLW
#else
       use module_radsw_parameters,  only: topfsw_type, sfcfsw_type
       use module_radlw_parameters,  only: topflw_type, sfcflw_type
       use ozne_def,                 only: levozp, oz_coeff
       use h2o_def,                  only: levh2o, h2o_coeff
       use aerclm_def,               only: ntrcaer, ntrcaerm
#endif


       implicit none

#ifdef CCPP
      ! To ensure that these values match what's in the physics,
      ! array sizes are compared during model init in GFS_rrtmg_setup_init()
      private :: NF_AESW, NF_AELW, NSPC, NSPC1, NF_CLDS, NF_VGAS, NF_ALBD, ntrcaerm
      ! from module_radiation_aerosols
      integer, parameter :: NF_AESW = 3
      integer, parameter :: NF_AELW = 3
      integer, parameter :: NSPC    = 5
      integer, parameter :: NSPC1   = NSPC + 1
      ! from module_radiation_clouds
      integer, parameter :: NF_CLDS = 9
      ! from module_radiation_gases
      integer, parameter :: NF_VGAS = 10
      ! from module_radiation_surface
      integer, parameter :: NF_ALBD = 4
      ! from aerclm_def
      integer, parameter :: ntrcaerm = 15

      ! These will be set later in GFS_Control%initialize,
      ! since they depend on the runtime config (e.g. Model%ntoz, Model%h2o_phys, Model%aero_in)
      private :: levozp, oz_coeff, levh2o, h2o_coeff, ntrcaer
      integer :: levozp, oz_coeff, levh2o, h2o_coeff, ntrcaer
#endif

#if 0
!> \section arg_table_GFS_typedefs
!! | local_name                      | standard_name                                            | long_name                                               | units         | rank | type                  |    kind   | intent | optional |
!! |---------------------------------|----------------------------------------------------------|---------------------------------------------------------|---------------|------|-----------------------|-----------|--------|----------|
!! | GFS_cldprop_type                | GFS_cldprop_type                                         | definition of type GFS_cldprop_type                     | DDT           |    0 | GFS_cldprop_type      |           | none   | F        |
!! | GFS_control_type                | GFS_control_type                                         | definition of type GFS_control_type                     | DDT           |    0 | GFS_control_type      |           | none   | F        |
!! | GFS_coupling_type               | GFS_coupling_type                                        | definition of type GFS_coupling_type                    | DDT           |    0 | GFS_coupling_type     |           | none   | F        |
!! | GFS_data_type                   | GFS_data_type                                            | definition of type GFS_data_type                        | DDT           |    0 | GFS_data_type         |           | none   | F        |
!! | GFS_diag_type                   | GFS_diag_type                                            | definition of type GFS_diag_type                        | DDT           |    0 | GFS_diag_type         |           | none   | F        |
!! | GFS_grid_type                   | GFS_grid_type                                            | definition of type GFS_grid_type                        | DDT           |    0 | GFS_grid_type         |           | none   | F        |
!! | GFS_interstitial_type           | GFS_interstitial_type                                    | definition of type GFS_interstitial_type                | DDT           |    0 | GFS_interstitial_type |           | none   | F        |
!! | GFS_radtend_type                | GFS_radtend_type                                         | definition of type GFS_radtend_type                     | DDT           |    0 | GFS_radtend_type      |           | none   | F        |
!! | GFS_sfcprop_type                | GFS_sfcprop_type                                         | definition of type GFS_sfcprop_type                     | DDT           |    0 | GFS_sfcprop_type      |           | none   | F        |
!! | GFS_statein_type                | GFS_statein_type                                         | definition of type GFS_statein_type                     | DDT           |    0 | GFS_statein_type      |           | none   | F        |
!! | GFS_stateout_type               | GFS_stateout_type                                        | definition of type GFS_stateout_type                    | DDT           |    0 | GFS_stateout_type     |           | none   | F        |
!! | GFS_tbd_type                    | GFS_tbd_type                                             | definition of type GFS_tbd_type                         | DDT           |    0 | GFS_tbd_type          |           | none   | F        |
!! | cmpfsw_type                     | cmpfsw_type                                              | definition of type cmpfsw_type                          | DDT           |    0 | cmpfsw_type           |           | none   | F        |
!! | sfcflw_type                     | sfcflw_type                                              | definition of type sfcflw_type                          | DDT           |    0 | sfcflw_type           |           | none   | F        |
!! | sfcfsw_type                     | sfcfsw_type                                              | definition of type sfcfsw_type                          | DDT           |    0 | sfcfsw_type           |           | none   | F        |
!! | topflw_type                     | topflw_type                                              | definition of type topflw_type                          | DDT           |    0 | topflw_type           |           | none   | F        |
!! | topfsw_type                     | topfsw_type                                              | definition of type topfsw_type                          | DDT           |    0 | topfsw_type           |           | none   | F        |
!! | LTP                             | extra_top_layer                                          | extra top layer for radiation                           | none          |    0 | integer               |           | none   | F        |
!! | con_cliq                        | specific_heat_of_liquid_water_at_constant_pressure       | specific heat of liquid water at constant pressure      | J kg-1 K-1    |    0 | real                  | kind_phys | none   | F        |
!! | con_cp                          | specific_heat_of_dry_air_at_constant_pressure            | specific heat of dry air at constant pressure           | J kg-1 K-1    |    0 | real                  | kind_phys | none   | F        |
!! | con_cvap                        | specific_heat_of_water_vapor_at_constant_pressure        | specific heat of water vapor at constant pressure       | J kg-1 K-1    |    0 | real                  | kind_phys | none   | F        |
!! | con_eps                         | ratio_of_dry_air_to_water_vapor_gas_constants            | rd/rv                                                   | none          |    0 | real                  | kind_phys | none   | F        |
!! | con_epsm1                       | ratio_of_dry_air_to_water_vapor_gas_constants_minus_one  | (rd/rv) - 1                                             | none          |    0 | real                  | kind_phys | none   | F        |
!! | con_fvirt                       | ratio_of_vapor_to_dry_air_gas_constants_minus_one        | (rv/rd) - 1 (rv = ideal gas constant for water vapor)   | none          |    0 | real                  | kind_phys | none   | F        |
!! | con_g                           | gravitational_acceleration                               | gravitational acceleration                              | m s-2         |    0 | real                  | kind_phys | none   | F        |
!! | con_hvap                        | latent_heat_of_vaporization_of_water_at_0C               | latent heat of evaporation/sublimation                  | J kg-1        |    0 | real                  | kind_phys | none   | F        |
!! | con_hfus                        | latent_heat_of_fusion_of_water_at_0C                     | latent heat of fusion                                   | J kg-1        |    0 | real                  | kind_phys | none   | F        |
!! | con_pi                          | pi                                                       | ratio of a circle's circumference to its diameter       | radians       |    0 | real                  | kind_phys | none   | F        |
!! | con_rd                          | gas_constant_dry_air                                     | ideal gas constant for dry air                          | J kg-1 K-1    |    0 | real                  | kind_phys | none   | F        |
!! | con_rv                          | gas_constant_water_vapor                                 | ideal gas constant for water vapor                      | J kg-1 K-1    |    0 | real                  | kind_phys | none   | F        |
!! | con_t0c                         | temperature_at_zero_celsius                              | temperature at 0 degrees Celsius                        | K             |    0 | real                  | kind_phys | none   | F        |
!! | con_ttp                         | triple_point_temperature_of_water                        | triple point temperature of water                       | K             |    0 | real                  | kind_phys | none   | F        |
!!
#endif

       !--- version of physics
       character(len=64) :: phys_version = 'v2018 FV3GFS BETA VERSION PHYSICS'

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
#ifdef CCPP
       !> optional extra top layer on top of low ceiling models
       !! this parameter was originally defined in the radiation driver
       !! (and is still for standard non-CCPP builds), but is required
       !! here for CCPP to allocate arrays used for the interstitial
       !! calculations previously in GFS_{physics,radiation}_driver.F90
       !! LTP=0: no extra top layer
       integer, parameter :: LTP = 0   !< no extra top layer
       !integer, parameter :: LTP = 1   ! add an extra top layer
#endif

!----------------
! Data Containers
!----------------
!    !--- GFS external initialization type
!    GFS_init_type
!    !--- GFS Derived Data Types (DDTs)
!    GFS_statein_type        !< prognostic state data in from dycore
!    GFS_stateout_type       !< prognostic state or tendencies return to dycore
!    GFS_sfcprop_type        !< surface fields
!    GFS_coupling_type       !< fields to/from coupling with other components (e.g. land/ice/ocean/etc.)
!    !---GFS specific containers
!    GFS_control_type        !< model control parameters
!    GFS_grid_type           !< grid and interpolation related data
!    GFS_tbd_type            !< to be determined data that doesn't fit in any one container
!    GFS_cldprop_type        !< cloud fields needed by radiation from physics
!    GFS_radtend_type        !< radiation tendencies needed in physics
!    GFS_diag_type           !< fields targetted for diagnostic output
#ifdef CCPP
!    GFS_interstitial_type   !< fields required to replace interstitial code in GFS_{physics,radiation}_driver.F90 in CCPP
!    GFS_data_type           !< combined type of all of the above except GFS_control_type and GFS_interstitial_type
#endif

!--------------------------------------------------------------------------------
! GFS_init_type
!--------------------------------------------------------------------------------
!>   This container is the minimum set of data required from the dycore/atmosphere
!!   component to allow proper initialization of the GFS physics
!--------------------------------------------------------------------------------
#if 0
!! \section arg_table_GFS_init_type
!! | local_name     | standard_name                                          | long_name                                               | units         | rank | type     |    kind   | intent | optional |
!! |----------------|--------------------------------------------------------|---------------------------------------------------------|---------------|------|----------|-----------|--------|----------|
!! | me             |                                                        | current MPI-rank                                        | none          |    0 | integer  |           | none   | F        |
!! | master         |                                                        | master MPI-rank                                         | none          |    0 | integer  |           | none   | F        |
!! | tile_num       |                                                        | tile number                                             | none          |    0 | integer  |           | none   | F        |
!! | isc            |                                                        | starting i-index for this MPI-domain                    | index         |    0 | integer  |           | none   | F        |
!! | jsc            |                                                        | starting j-index for this MPI-domain                    | index         |    0 | integer  |           | none   | F        |
!! | nx             |                                                        | number of points in i-dir for this MPI rank             | count         |    0 | integer  |           | none   | F        |
!! | ny             |                                                        | number of points in j-dir for this MPI rank             | count         |    0 | integer  |           | none   | F        |
!! | levs           |                                                        | number of vertical levels                               | count         |    0 | integer  |           | none   | F        |
!! | cnx            |                                                        | number of points in i-dir for this cubed-sphere face    | count         |    0 | integer  |           | none   | F        |
!! | cny            |                                                        | number of points in j-dir for this cubed-sphere face    | count         |    0 | integer  |           | none   | F        |
!! | gnx            |                                                        | number of global points in x-dir (i) along the equator  | count         |    0 | integer  |           | none   | F        |
!! | gny            |                                                        | number of global points in y-dir (j) along any meridian | count         |    0 | integer  |           | none   | F        |
!! | nlunit         |                                                        | fortran unit number for file opens                      | none          |    0 | integer  |           | none   | F        |
!! | logunit        |                                                        | fortran unit number for writing logfile                 | none          |    0 | integer  |           | none   | F        |
!! | bdat           |                                                        | model begin date in GFS format   (same as idat)         | none          |    0 | integer  |           | none   | F        |
!! | cdat           |                                                        | model current date in GFS format (same as jdat)         | none          |    0 | integer  |           | none   | F        |
!! | dt_dycore      |                                                        | dynamics time step in seconds                           | s             |    0 | real     | kind_phys | none   | F        |
!! | dt_phys        |                                                        | physics  time step in seconds                           | s             |    0 | real     | kind_phys | none   | F        |
!! | restart        |                                                        | flag for restart (warmstart) or coldstart               | flag          |    0 | logical  |           | none   | F        |
!! | hydrostatic    |                                                        | flag for hydrostatic solver from dynamics               | flag          |    0 | logical  |           | none   | F        |
!! | blksz          |                                                        | for explicit data blocking                              | count         |    1 | integer  |           | none   | F        |
!! | ak             |                                                        | a parameter for sigma pressure level calculations       | Pa            |    1 | real     | kind_phys | none   | F        |
!! | bk             |                                                        | b parameter for sigma pressure level calculations       | none          |    1 | real     | kind_phys | none   | F        |
!! | xlon           |                                                        | column longitude for MPI rank                           | radians       |    2 | real     | kind_phys | none   | F        |
!! | xlat           |                                                        | column latitude  for MPI rank                           | radians       |    2 | real     | kind_phys | none   | F        |
!! | area           |                                                        | column area for length scale calculations               | m2            |    2 | real     | kind_phys | none   | F        |
!! | tracer_names   |                                                        | tracers names to dereference tracer id                  | none          |    1 | charater | len=32    | none   | F        |
!! | fn_nml         |                                                        | namelist filename                                       | none          |    0 | charater | len=65    | none   | F        |
!! | input_nml_file |                                                        | namelist filename for internal file reads               | none          |    0 | charater | len=256   | none   | F        |
!!
#endif
  type GFS_init_type
    integer :: me                                !< my MPI-rank
    integer :: master                            !< master MPI-rank
    integer :: tile_num                          !< tile number for this MPI rank
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
#ifdef CCPP
!--- restart information
    logical :: restart                           !< flag whether this is a coldstart (.false.) or a warmstart/restart (.true.)
!--- hydrostatic/non-hydrostatic flag
    logical :: hydrostatic                       !< flag whether this is a hydrostatic or non-hydrostatic run
#endif
!--- blocking data
    integer, pointer :: blksz(:)                 !< for explicit data blocking
                                                 !< default blksz(1)=[nx*ny]
!--- ak/bk for pressure level calculations
    real(kind=kind_phys), pointer :: ak(:)       !< from surface (k=1) to TOA (k=levs)
    real(kind=kind_phys), pointer :: bk(:)       !< from surface (k=1) to TOA (k=levs)
!--- grid metrics
    real(kind=kind_phys), pointer :: xlon(:,:)   !< column longitude for MPI rank
    real(kind=kind_phys), pointer :: xlat(:,:)   !< column latitude  for MPI rank
    real(kind=kind_phys), pointer :: area(:,:)   !< column area for length scale calculations

    character(len=32), pointer :: tracer_names(:) !< tracers names to dereference tracer id
                                                  !< based on name location in array
    character(len=65) :: fn_nml                   !< namelist filename
    character(len=256), pointer :: input_nml_file(:) !< character string containing full namelist
                                                     !< for use with internal file reads
  end type GFS_init_type


!----------------------------------------------------------------
! GFS_statein_type
!   prognostic state variables with layer and level specific data
!----------------------------------------------------------------
#if 0
!! \section arg_table_GFS_statein_type
!! | local_name                                                | standard_name                                             | long_name                                                                           | units         | rank | type    |    kind   | intent | optional |
!! |-----------------------------------------------------------|-----------------------------------------------------------|-------------------------------------------------------------------------------------|---------------|------|---------|-----------|--------|----------|
!! | GFS_Data(cdata%blk_no)%Statein%phii                       | geopotential_at_interface                                 | geopotential at model layer interfaces                                              | m2 s-2        |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Statein%prsi                       | air_pressure_at_interface                                 | air pressure at model layer interfaces                                              | Pa            |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Statein%prsik                      | dimensionless_exner_function_at_model_interfaces          | dimensionless Exner function at model layer interfaces                              | none          |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Statein%prsik(:,1)                 | dimensionless_exner_function_at_lowest_model_interface    | dimensionless Exner function at lowest model interface                              | none          |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Statein%phil                       | geopotential                                              | geopotential at model layer centers                                                 | m2 s-2        |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Statein%prsl                       | air_pressure                                              | mean layer pressure                                                                 | Pa            |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Statein%prsl(:,1)                  | air_pressure_at_lowest_model_layer                        | mean pressure at lowest model layer                                                 | Pa            |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Statein%prslk                      | dimensionless_exner_function_at_model_layers              | dimensionless Exner function at model layer centers                                 | none          |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Statein%prslk(:,1)                 | dimensionless_exner_function_at_lowest_model_layer        | dimensionless Exner function at lowest model layer                                  | none          |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Statein%pgr                        | surface_air_pressure                                      | surface pressure                                                                    | Pa            |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Statein%ugrs                       | x_wind                                                    | zonal wind                                                                          | m s-1         |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Statein%ugrs(:,1)                  | x_wind_at_lowest_model_layer                              | zonal wind at lowest model layer                                                    | m s-1         |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Statein%vgrs                       | y_wind                                                    | meridional wind                                                                     | m s-1         |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Statein%vgrs(:,1)                  | y_wind_at_lowest_model_layer                              | meridional wind at lowest model layer                                               | m s-1         |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Statein%vvl                        | omega                                                     | layer mean vertical velocity                                                        | Pa s-1        |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Statein%tgrs                       | air_temperature                                           | model layer mean temperature                                                        | K             |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Statein%tgrs(:,1)                  | air_temperature_at_lowest_model_layer                     | mean temperature at lowest model layer                                              | K             |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Statein%qgrs                       | tracer_concentration                                      | model layer mean tracer concentration                                               | kg kg-1       |    3 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Statein%qgrs(:,:,GFS_Control%ntqv) | water_vapor_specific_humidity                             | water vapor specific humidity                                                       | kg kg-1       |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Statein%qgrs(:,1,GFS_Control%ntqv) | water_vapor_specific_humidity_at_lowest_model_layer       | water vapor specific humidity at lowest model layer                                 | kg kg-1       |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Statein%qgrs(:,:,GFS_Control%ntcw) | cloud_condensed_water_mixing_ratio                        | moist (dry+vapor, no condensates) mixing ratio of cloud water (condensate)          | kg kg-1       |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Statein%qgrs(:,1,GFS_Control%ntcw) | cloud_condensed_water_mixing_ratio_at_lowest_model_layer  | moist (dry+vapor, no condensates) mixing ratio of cloud water at lowest model layer | kg kg-1       |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Statein%qgrs(:,:,GFS_Control%ntiw) | ice_water_mixing_ratio                                    | moist (dry+vapor, no condensates) mixing ratio of ice water                         | kg kg-1       |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Statein%qgrs(:,:,GFS_Control%ntrw) | rain_water_mixing_ratio                                   | moist (dry+vapor, no condensates) mixing ratio of rain water                        | kg kg-1       |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Statein%qgrs(:,:,GFS_Control%ntsw) | snow_water_mixing_ratio                                   | moist (dry+vapor, no condensates) mixing ratio of snow water                        | kg kg-1       |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Statein%qgrs(:,:,GFS_Control%ntgl) | graupel_mixing_ratio                                      | moist (dry+vapor, no condensates) mixing ratio of graupel                           | kg kg-1       |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Statein%qgrs(:,:,GFS_Control%ntoz) | ozone_mixing_ratio                                        | ozone mixing ratio                                                                  | kg kg-1       |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Statein%qgrs(:,:,GFS_Control%ntwa) | water_friendly_aerosol_number_concentration               | number concentration of water-friendly aerosols                                     | kg-1          |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Statein%qgrs(:,:,GFS_Control%ntia) | ice_friendly_aerosol_number_concentration                 | number concentration of ice-friendly aerosols                                       | kg-1          |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Statein%qgrs(:,:,GFS_Control%ntlnc)| cloud_droplet_number_concentration                        | number concentration of cloud droplets (liquid)                                     | kg-1          |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Statein%qgrs(:,:,GFS_Control%ntinc)| ice_number_concentration                                  | number concentration of ice                                                         | kg-1          |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Statein%qgrs(:,:,GFS_Control%ntrnc)| rain_number_concentration                                 | number concentration of rain                                                        | kg-1          |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Statein%qgrs(:,:,GFS_Control%ntsnc)| snow_number_concentration                                 | number concentration of snow                                                        | kg-1          |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Statein%qgrs(:,:,GFS_Control%ntgnc)| graupel_number_concentration                              | number concentration of graupel                                                     | kg-1          |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Statein%qgrs(:,:,GFS_Control%ntke) | turbulent_kinetic_energy                                  | turbulent kinetic energy                                                            | J             |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Statein%diss_est                   | dissipation_estimate_of_air_temperature_at_model_layers   | dissipation estimate model layer mean temperature                                   | K             |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Statein%smc                        |                                                           | total soil moisture                                                                 | frac          |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Statein%stc                        |                                                           | soil temperature                                                                    | K             |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Statein%slc                        |                                                           | liquid soil moisture                                                                | frac          |    2 | real    | kind_phys | none   | F        |
!!
#endif
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
! dissipation estimate
    real (kind=kind_phys), pointer :: diss_est(:,:)   => null()  !< model layer mean temperature in k
    ! soil state variables - for soil SPPT - sfc-perts, mgehne
    real (kind=kind_phys), pointer :: smc (:,:)   => null()  !< soil moisture content
    real (kind=kind_phys), pointer :: stc (:,:)   => null()  !< soil temperature content
    real (kind=kind_phys), pointer :: slc (:,:)   => null()  !< soil liquid water content

    contains
      procedure :: create  => statein_create  !<   allocate array data
  end type GFS_statein_type


!------------------------------------------------------------------
! GFS_stateout_type
!   prognostic state or tendencies after physical parameterizations
!------------------------------------------------------------------
#if 0
!! \section arg_table_GFS_stateout_type
!! | local_name                                                   | standard_name                                                          | long_name                                                                                  | units   | rank | type    |    kind   | intent | optional |
!! |--------------------------------------------------------------|------------------------------------------------------------------------|--------------------------------------------------------------------------------------------|---------|------|---------|-----------|--------|----------|
!! | GFS_Data(cdata%blk_no)%Stateout%gu0                          | x_wind_updated_by_physics                                              | zonal wind updated by physics                                                              | m s-1   |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Stateout%gu0(:,1)                     | x_wind_at_lowest_model_layer_updated_by_physics                        | zonal wind at lowest model level updated by physics                                        | m s-1   |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Stateout%gv0                          | y_wind_updated_by_physics                                              | meridional wind updated by physics                                                         | m s-1   |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Stateout%gv0(:,1)                     | y_wind_at_lowest_model_layer_updated_by_physics                        | meridional wind at lowest model level updated by physics                                   | m s-1   |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Stateout%gt0                          | air_temperature_updated_by_physics                                     | temperature updated by physics                                                             | K       |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Stateout%gt0(:,1)                     | air_temperature_at_lowest_model_layer_updated_by_physics               | temperature at lowest model layer updated by physics                                       | K       |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Stateout%gq0                          | tracer_concentration_updated_by_physics                                | tracer concentration updated by physics                                                    | kg kg-1 |    3 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Stateout%gq0(:,:,GFS_Control%ntqv)    | water_vapor_specific_humidity_updated_by_physics                       | water vapor specific humidity updated by physics                                           | kg kg-1 |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Stateout%gq0(:,1,GFS_Control%ntqv)    | water_vapor_specific_humidity_at_lowest_model_layer_updated_by_physics | water vapor specific humidity at lowest model layer updated by physics                     | kg kg-1 |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Stateout%gq0(:,:,GFS_Control%ntoz)    | ozone_concentration_updated_by_physics                                 | ozone concentration updated by physics                                                     | kg kg-1 |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Stateout%gq0(:,:,GFS_Control%ntcw)    | cloud_condensed_water_mixing_ratio_updated_by_physics                  | moist (dry+vapor, no condensates) mixing ratio of cloud condensed water updated by physics | kg kg-1 |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Stateout%gq0(:,:,GFS_Control%ntiw)    | ice_water_mixing_ratio_updated_by_physics                              | moist (dry+vapor, no condensates) mixing ratio of ice water updated by physics             | kg kg-1 |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Stateout%gq0(:,:,GFS_Control%ntrw)    | rain_water_mixing_ratio_updated_by_physics                             | moist (dry+vapor, no condensates) mixing ratio of rain water updated by physics            | kg kg-1 |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Stateout%gq0(:,:,GFS_Control%ntsw)    | snow_water_mixing_ratio_updated_by_physics                             | moist (dry+vapor, no condensates) mixing ratio of snow water updated by physics            | kg kg-1 |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Stateout%gq0(:,:,GFS_Control%ntgl)    | graupel_mixing_ratio_updated_by_physics                                | moist (dry+vapor, no condensates) mixing ratio of graupel updated by physics               | kg kg-1 |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Stateout%gq0(:,:,GFS_Control%ntwa)    | water_friendly_aerosol_number_concentration_updated_by_physics         | number concentration of water-friendly aerosols updated by physics                         | kg-1    |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Stateout%gq0(:,:,GFS_Control%ntia)    | ice_friendly_aerosol_number_concentration_updated_by_physics           | number concentration of ice-friendly aerosols updated by physics                           | kg-1    |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Stateout%gq0(:,:,GFS_Control%ntlnc)   | cloud_droplet_number_concentration_updated_by_physics                  | number concentration of cloud droplets updated by physics                                  | kg-1    |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Stateout%gq0(:,:,GFS_Control%ntinc)   | ice_number_concentration_updated_by_physics                            | number concentration of ice updated by physics                                             | kg-1    |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Stateout%gq0(:,:,GFS_Control%ntrnc)   | rain_number_concentration_updated_by_physics                           | number concentration of rain updated by physics                                            | kg-1    |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Stateout%gq0(:,:,GFS_Control%ntsnc)   | snow_number_concentration_updated_by_physics                           | number concentration of snow updated by physics                                            | kg-1    |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Stateout%gq0(:,:,GFS_Control%ntgnc)   | graupel_number_concentration_updated_by_physics                        | number concentration of graupel updated by physics                                         | kg-1    |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Stateout%gq0(:,:,GFS_Control%ntclamt) | cloud_fraction_updated_by_physics                                      | cloud fraction updated by physics                                                          | frac    |    2 | real    | kind_phys | none   | F        |
!!
#endif
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
!> This container defines surface properties that may be read in 
!! and/or updated by climatology or observations
!---------------------------------------------------------------------------------------
#if 0
!! \section arg_table_GFS_sfcprop_type
!! | local_name                                 | standard_name                                                          | long_name                                              | units         | rank | type    |    kind   | intent | optional |
!! |--------------------------------------------|------------------------------------------------------------------------|--------------------------------------------------------|---------------|------|---------|-----------|--------|----------|
!! | GFS_Data(cdata%blk_no)%Sfcprop%slmsk       | sea_land_ice_mask_real                                                 | landmask: sea/land/ice=0/1/2                           | flag          |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%lakemsk     | lake_mask_real                                                         | lake mask: non-lake/lake=0/1                           | flag          |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%tsfc        | surface_skin_temperature                                               | surface skin temperature                               | K             |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%tisfc       | sea_ice_temperature                                                    | sea uce surface skin temperature                       | K             |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%snowd       | surface_snow_thickness_water_equivalent                                | water equivalent snow depth over land                  | mm            |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%zorl        | surface_roughness_length                                               | surface roughness length                               | cm            |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%fice        | sea_ice_concentration                                                  | ice fraction over open water                           | frac          |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%hprim       |                                                                        | topographic standard deviation                         | m             |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%hprime      | statistical_measures_of_subgrid_orography                              | orographic metrics                                     | various       |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%sncovr      | surface_snow_area_fraction_for_diagnostics                             | surface snow area fraction                             | frac          |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%snoalb      | upper_bound_on_max_albedo_over_deep_snow                               | maximum snow albedo                                    | frac          |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%alvsf       |                                                                        | mean vis albedo with strong cosz dependency            | frac          |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%alnsf       |                                                                        | mean nir albedo with strong cosz dependency            | frac          |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%alvwf       | mean_vis_albedo_with_weak_cosz_dependency                              | mean vis albedo with weak cosz dependency              | frac          |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%alnwf       | mean_nir_albedo_with_weak_cosz_dependency                              | mean nir albedo with weak cosz dependency              | frac          |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%facsf       |                                                                        | fractional coverage with strong cosz dependency        | frac          |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%facwf       |                                                                        | fractional coverage with weak cosz dependency          | frac          |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%slope       | surface_slope_classification_real                                      | sfc slope type for lsm                                 | index         |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%shdmin      | minimum_vegetation_area_fraction                                       | min fractional coverage of green vegetation            | frac          |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%shdmax      | maximum_vegetation_area_fraction                                       | max fractional coverage of green vegetation            | frac          |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%tg3         | deep_soil_temperature                                                  | deep soil temperature                                  | K             |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%vfrac       | vegetation_area_fraction                                               | areal fractional cover of green vegetation             | frac          |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%vtype       | vegetation_type_classification_real                                    | vegetation type for lsm                                | index         |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%stype       | soil_type_classification_real                                          | soil type for lsm                                      | index         |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%uustar      | surface_friction_velocity                                              | boundary layer parameter                               | m s-1         |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%oro         | orography                                                              | orography                                              | m             |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%oro_uf      | orography_unfiltered                                                   | unfiltered orography                                   | m             |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%conv_act    | gf_memory_counter                                                      | Memory counter for GF                                  | none          |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%hice        | sea_ice_thickness                                                      | sea ice thickness                                      | m             |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%weasd       | water_equivalent_accumulated_snow_depth                                | water equiv of acc snow depth over land and sea ice    | mm            |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%canopy      | canopy_water_amount                                                    | canopy water amount                                    | kg m-2        |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%ffmm        | Monin-Obukhov_similarity_function_for_momentum                         | Monin-Obukhov similarity function for momentum         | none          |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%ffhh        | Monin-Obukhov_similarity_function_for_heat                             | Monin-Obukhov similarity function for heat             | none          |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%f10m        | ratio_of_wind_at_lowest_model_layer_and_wind_at_10m                    | ratio of sigma level 1 wind and 10m wind               | ratio         |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%tprcp       | nonnegative_lwe_thickness_of_precipitation_amount_on_dynamics_timestep | total precipitation amount in each time step           | m             |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%srflag      | flag_for_precipitation_type                                            | snow/rain flag for precipitation                       | flag          |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%sr          | ratio_of_snowfall_to_rainfall                                          | snow ratio: ratio of snow to total precipitation       | frac          |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%slc         | volume_fraction_of_unfrozen_soil_moisture                              | liquid soil moisture                                   | frac          |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%smc         | volume_fraction_of_soil_moisture                                       | total soil moisture                                    | frac          |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%stc         | soil_temperature                                                       | soil temperature                                       | K             |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%wet1        | normalized_soil_wetness                                                | normalized soil wetness                                | frac          |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%t2m         | temperature_at_2m                                                      | 2 meter temperature                                    | K             |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%th2m        | potential_temperature_at_2m                                            | 2 meter potential temperature                          | K             |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%q2m         | specific_humidity_at_2m                                                | 2 meter specific humidity                              | kg kg-1       |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%tref        | sea_surface_reference_temperature                                      | sea surface reference temperature                      | K             |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%z_c         | sub-layer_cooling_thickness                                            | sub-layer cooling thickness                            | m             |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%c_0         | coefficient_c_0                                                        | coefficient 1 to calculate d(Tz)/d(Ts)                 | none          |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%c_d         | coefficient_c_d                                                        | coefficient 2 to calculate d(Tz)/d(Ts)                 | none          |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%w_0         | coefficient_w_0                                                        | coefficient 3 to calculate d(Tz)/d(Ts)                 | none          |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%w_d         | coefficient_w_d                                                        | coefficient 4 to calculate d(Tz)/d(Ts)                 | none          |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%xt          | diurnal_thermocline_layer_heat_content                                 | heat content in diurnal thermocline layer              | K m           |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%xs          | sea_water_salinity                                                     | salinity  content in diurnal thermocline layer         | ppt m         |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%xu          | diurnal_thermocline_layer_x_current                                    | u-current content in diurnal thermocline layer         | m2 s-1        |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%xv          | diurnal_thermocline_layer_y_current                                    | v-current content in diurnal thermocline layer         | m2 s-1        |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%xz          | diurnal_thermocline_layer_thickness                                    | diurnal thermocline layer thickness                    | m             |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%zm          | ocean_mixed_layer_thickness                                            | mixed layer thickness                                  | m             |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%xtts        | sensitivity_of_dtl_heat_content_to_surface_temperature                 | d(xt)/d(ts)                                            | m             |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%xzts        | sensitivity_of_dtl_thickness_to_surface_temperature                    | d(xz)/d(ts)                                            | m K-1         |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%d_conv      | free_convection_layer_thickness                                        | thickness of free convection layer (FCL)               | m             |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%ifd         | index_of_dtlm_start                                                    | index to start dtlm run or not                         | index         |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%dt_cool     | sub-layer_cooling_amount                                               | sub-layer cooling amount                               | K             |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%qrain       | sensible_heat_flux_due_to_rainfall                                     | sensible heat flux due to rainfall                     | W             |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%sh2o        | volume_fraction_of_unfrozen_soil_moisture_for_land_surface_model       | volume fraction of unfrozen soil moisture for lsm      | frac          |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%keepsmfr    | volume_fraction_of_frozen_soil_moisture_for_land_surface_model         | volume fraction of frozen soil moisture for lsm        | frac          |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%smois       | volume_fraction_of_soil_moisture_for_land_surface_model                | volumetric fraction of soil moisture for lsm           | frac          |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%tslb        | soil_temperature_for_land_surface_model                                | soil temperature for land surface model                | K             |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%zs          | depth_of_soil_levels_for_land_surface_model                            | depth of soil levels for land surface model            | m             |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%clw_surf    | cloud_condensed_water_mixing_ratio_at_surface                          | moist cloud water mixing ratio at surface              | kg kg-1       |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%qwv_surf    | water_vapor_mixing_ratio_at_surface                                    | water vapor mixing ratio at surface                    | kg kg-1       |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%cndm_surf   | surface_condensation_mass                                              | surface condensation mass                              | kg m-2        |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%flag_frsoil | flag_for_frozen_soil_physics                                           | flag for frozen soil physics (RUC)                     | flag          |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%rhofr       | density_of_frozen_precipitation                                        | density of frozen precipitation                        | kg m-3        |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%tsnow       | snow_temperature_bottom_first_layer                                    | snow temperature at the bottom of the first snow layer | K             |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%ustm        | surface_friction_velocity_drag                                         | friction velocity isolated for momentum only           | m s-1         |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%zol         | surface_stability_parameter                                            | monin obukhov surface stability parameter              | none          |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%mol         | theta_star                                                             | temperature flux divided by ustar (temperature scale)  | K             |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%rmol        | reciprocal_of_obukhov_length                                           | one over obukhov length                                | m-1           |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%flhc        | surface_exchange_coefficient_for_heat                                  | surface exchange coefficient for heat                  | W m-2 K-1     |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%flqc        | surface_exchange_coefficient_for_moisture                              | surface exchange coefficient for moisture              | kg m-2 s-1    |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%chs2        | surface_exchange_coefficient_for_heat_at_2m                            | exchange coefficient for heat at 2 meters              | m s-1         |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%cqs2        | surface_exchange_coefficient_for_moisture_at_2m                        | exchange coefficient for moisture at 2 meters          | m s-1         |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Sfcprop%lh          | surface_latent_heat                                                    | latent heating at the surface (pos = up)               | W m-2         |    1 | real    | kind_phys | none   | F        |
!!
#endif
  type GFS_sfcprop_type

!--- In (radiation and physics)
    real (kind=kind_phys), pointer :: slmsk  (:)   => null()  !< sea/land mask array (sea:0,land:1,sea-ice:2)
    real (kind=kind_phys), pointer :: lakemsk(:)   => null()  !< lake mask array (lake:1, non-lake:0)
    real (kind=kind_phys), pointer :: tsfc   (:)   => null()  !< surface temperature in k
                                                              !< [tsea in gbphys.f]
    real (kind=kind_phys), pointer :: tisfc  (:)   => null()  !< surface temperature over ice fraction
    real (kind=kind_phys), pointer :: snowd  (:)   => null()  !< snow depth water equivalent in mm ; same as snwdph
    real (kind=kind_phys), pointer :: zorl   (:)   => null()  !< surface roughness in cm
    real (kind=kind_phys), pointer :: fice   (:)   => null()  !< ice fraction over open water grid
    real (kind=kind_phys), pointer :: hprim  (:)   => null()  !< topographic standard deviation in m
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
#ifdef CCPP
    real (kind=kind_phys), pointer :: conv_act(:)  => null()  !< convective activity counter hli 09/2017
#endif
    real (kind=kind_phys), pointer :: hice   (:)   => null()  !< sea ice thickness
    real (kind=kind_phys), pointer :: weasd  (:)   => null()  !< water equiv of accumulated snow depth (kg/m**2)
                                                              !< over land and sea ice
    real (kind=kind_phys), pointer :: canopy (:)   => null()  !< canopy water
    real (kind=kind_phys), pointer :: ffmm   (:)   => null()  !< fm parameter from PBL scheme
    real (kind=kind_phys), pointer :: ffhh   (:)   => null()  !< fh parameter from PBL scheme
    real (kind=kind_phys), pointer :: f10m   (:)   => null()  !< fm at 10m - Ratio of sigma level 1 wind and 10m wind
    real (kind=kind_phys), pointer :: tprcp  (:)   => null()  !< sfc_fld%tprcp - total precipitation
    real (kind=kind_phys), pointer :: srflag (:)   => null()  !< sfc_fld%srflag - snow/rain flag for precipitation
    real (kind=kind_phys), pointer :: sr     (:)   => null()  !< snow ratio : ratio of snow to total precipitation
    real (kind=kind_phys), pointer :: slc    (:,:) => null()  !< liquid soil moisture
    real (kind=kind_phys), pointer :: smc    (:,:) => null()  !< total soil moisture
    real (kind=kind_phys), pointer :: stc    (:,:) => null()  !< soil temperature
    real (kind=kind_phys), pointer :: wet1   (:)   => null()  !< normalized soil wetness

!--- Out
    real (kind=kind_phys), pointer :: t2m    (:)   => null()  !< 2 meter temperature
#ifdef CCPP
    real (kind=kind_phys), pointer :: th2m   (:)   => null()  !< 2 meter potential temperature
#endif
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

#ifdef CCPP
    ! Soil properties for land-surface model (if number of levels different from NOAH 4-layer model)
    real (kind=kind_phys), pointer :: sh2o(:,:)        => null()  !< volume fraction of unfrozen soil moisture for lsm
    real (kind=kind_phys), pointer :: keepsmfr(:,:)    => null()  !< RUC LSM: frozen moisture in soil
    real (kind=kind_phys), pointer :: smois(:,:)       => null()  !< volumetric fraction of soil moisture for lsm
    real (kind=kind_phys), pointer :: tslb(:,:)        => null()  !< soil temperature for land surface model
    real (kind=kind_phys), pointer :: flag_frsoil(:,:) => null()  !< RUC LSM: flag for frozen soil physics
    !
    real (kind=kind_phys), pointer :: zs(:)            => null()  !< depth of soil levels for land surface model
    real (kind=kind_phys), pointer :: clw_surf(:)      => null()  !< RUC LSM: moist cloud water mixing ratio at surface
    real (kind=kind_phys), pointer :: qwv_surf(:)      => null()  !< RUC LSM: water vapor mixing ratio at surface
    real (kind=kind_phys), pointer :: cndm_surf(:)     => null()  !< RUC LSM: surface condensation mass
    real (kind=kind_phys), pointer :: rhofr(:)         => null()  !< RUC LSM: density of frozen precipitation
    real (kind=kind_phys), pointer :: tsnow(:)         => null()  !< RUC LSM: snow temperature at the bottom of the first soil layer
    !  MYNN surface layer
    real (kind=kind_phys), pointer :: ustm (:)         => null()  !u* including drag
    real (kind=kind_phys), pointer :: zol(:)           => null()  !surface stability parameter
    real (kind=kind_phys), pointer :: mol(:)           => null()  !theta star
    real (kind=kind_phys), pointer :: rmol(:)          => null()  !reciprocal of obukhov length
    real (kind=kind_phys), pointer :: flhc(:)          => null()  !drag coeff for heat
    real (kind=kind_phys), pointer :: flqc(:)          => null()  !drag coeff for moisture
    real (kind=kind_phys), pointer :: chs2(:)          => null()  !exch coeff for heat at 2m
    real (kind=kind_phys), pointer :: cqs2(:)          => null()  !exch coeff for moisture at 2m
    real (kind=kind_phys), pointer :: lh(:)            => null()  !latent heating at the surface
#endif

    contains
      procedure :: create  => sfcprop_create  !<   allocate array data
  end type GFS_sfcprop_type


!---------------------------------------------------------------------
! GFS_coupling_type
!   fields to/from other coupled components (e.g. land/ice/ocean/etc.)
!---------------------------------------------------------------------
#if 0
!! \section arg_table_GFS_coupling_type
!! | local_name                           | standard_name                                                                             | long_name                                            | units         | rank | type    |    kind   | intent | optional |
!! |--------------------------------------|-------------------------------------------------------------------------------------------|------------------------------------------------------|---------------|------|---------|-----------|--------|----------|
!! | GFS_Data(cdata%blk_no)%Coupling%nirbmdi        | surface_downwelling_direct_near_infrared_shortwave_flux_on_radiation_time_step            | sfc nir beam sw downward flux                        | W m-2         |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%nirdfdi        | surface_downwelling_diffuse_near_infrared_shortwave_flux_on_radiation_time_step           | sfc nir diff sw downward flux                        | W m-2         |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%visbmdi        | surface_downwelling_direct_ultraviolet_and_visible_shortwave_flux_on_radiation_time_step  | sfc uv+vis beam sw downward flux                     | W m-2         |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%visdfdi        | surface_downwelling_diffuse_ultraviolet_and_visible_shortwave_flux_on_radiation_time_step | sfc uv+vis diff sw downward flux                     | W m-2         |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%nirbmui        | surface_upwelling_direct_near_infrared_shortwave_flux_on_radiation_time_step              | sfc nir beam sw upward flux                          | W m-2         |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%nirdfui        | surface_upwelling_diffuse_near_infrared_shortwave_flux_on_radiation_time_step             | sfc nir diff sw upward flux                          | W m-2         |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%visbmui        | surface_upwelling_direct_ultraviolet_and_visible_shortwave_flux_on_radiation_time_step    | sfc uv+vis beam sw upward flux                       | W m-2         |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%visdfui        | surface_upwelling_diffuse_ultraviolet_and_visible_shortwave_flux_on_radiation_time_step   | sfc uv+vis diff sw upward flux                       | W m-2         |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%sfcdsw         | surface_downwelling_shortwave_flux_on_radiation_time_step                                 | total sky sfc downward sw flux                       | W m-2         |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%sfcnsw         | surface_net_downwelling_shortwave_flux_on_radiation_time_step                             | total sky sfc netsw flx into ground                  | W m-2         |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%sfcdlw         | surface_downwelling_longwave_flux_on_radiation_time_step                                  | total sky sfc downward lw flux                       | W m-2         |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%dusfcin_cpl    |                                                                                           | aoi_fld%dusfcin(item,lan)                            |               |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%dvsfcin_cpl    |                                                                                           | aoi_fld%dvsfcin(item,lan)                            |               |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%dtsfcin_cpl    |                                                                                           | aoi_fld%dtsfcin(item,lan)                            |               |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%dqsfcin_cpl    |                                                                                           | aoi_fld%dqsfcin(item,lan)                            |               |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%ulwsfcin_cpl   |                                                                                           | aoi_fld%ulwsfcin(item,lan)                           |               |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%tseain_cpl     |                                                                                           | aoi_fld%tseain(item,lan)                             |               |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%tisfcin_cpl    |                                                                                           | aoi_fld%tisfcin(item,lan)                            |               |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%ficein_cpl     |                                                                                           | aoi_fld%ficein(item,lan)                             |               |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%hicein_cpl     |                                                                                           | aoi_fld%hicein(item,lan)                             |               |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%hsnoin_cpl     |                                                                                           | aoi_fld%hsnoin(item,lan)                             |               |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%slimskin_cpl   |                                                                                           | aoi_fld%slimskin(item,lan)                           |               |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%rain_cpl       | lwe_thickness_of_precipitation_amount_for_coupling                                        | total rain precipitation                             | m             |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%rainc_cpl      | lwe_thickness_of_convective_precipitation_amount_for_coupling                             | total convective precipitation                       | m             |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%snow_cpl       | lwe_thickness_of_snow_amount_for_coupling                                                 | total snow precipitation                             | m             |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%dusfc_cpl      | cumulative_surface_x_momentum_flux_for_coupling_multiplied_by_timestep                    | cumulative sfc x momentum flux multiplied by timestep| Pa s          |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%dvsfc_cpl      | cumulative_surface_y_momentum_flux_for_coupling_multiplied_by_timestep                    | cumulative sfc y momentum flux multiplied by timestep| Pa s          |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%dtsfc_cpl      | cumulative_surface_upward_sensible_heat_flux_for_coupling_multiplied_by_timestep          | cumulative sfc sensible heat flux multiplied by timestep | W m-2 s   |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%dqsfc_cpl      | cumulative_surface_upward_latent_heat_flux_for_coupling_multiplied_by_timestep            | cumulative sfc latent heat flux multiplied by timestep | W m-2 s     |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%dlwsfc_cpl     | cumulative_surface_downwelling_longwave_flux_for_coupling_multiplied_by_timestep          | cumulative sfc downward lw flux mulitplied by timestep | W m-2 s     |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%dswsfc_cpl     | cumulative_surface_downwelling_shortwave_flux_for_coupling_multiplied_by_timestep         | cumulative sfc downward sw flux multiplied by timestep | W m-2 s     |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%dnirbm_cpl     | cumulative_surface_downwelling_direct_near_infrared_shortwave_flux_for_coupling_multiplied_by_timestep | cumulative sfc nir beam downward sw flux multiplied by timestep            | W m-2 s       |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%dnirdf_cpl     | cumulative_surface_downwelling_diffuse_near_infrared_shortwave_flux_for_coupling_multiplied_by_timestep | cumulative sfc nir diff downward sw flux multiplied by timestep                       | W m-2 s       |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%dvisbm_cpl     | cumulative_surface_downwelling_direct_ultraviolet_and_visible_shortwave_flux_for_coupling_multiplied_by_timestep | cumulative sfc uv+vis beam dnwd sw flux multiplied by timestep                        | W m-2 s       |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%dvisdf_cpl     | cumulative_surface_downwelling_diffuse_ultraviolet_and_visible_shortwave_flux_for_coupling_multiplied_by_timestep | cumulative sfc uv+vis diff dnwd sw flux multiplied by timestep                        | W m-2 s       |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%nlwsfc_cpl     | cumulative_surface_net_downward_longwave_flux_for_coupling_multiplied_by_timestep         | cumulative net downward lw flux multiplied by timestep                     | W m-2 s       |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%nswsfc_cpl     | cumulative_surface_net_downward_shortwave_flux_for_coupling_multiplied_by_timestep        | cumulative net downward sw flux multiplied by timestep                      | W m-2 s       |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%nnirbm_cpl     | cumulative_surface_net_downward_direct_near_infrared_shortwave_flux_for_coupling_multiplied_by_timestep | cumulative net nir beam downward sw flux multiplied by timestep | W m-2 s       |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%nnirdf_cpl     | cumulative_surface_net_downward_diffuse_near_infrared_shortwave_flux_for_coupling_multiplied_by_timestep | cumulative net nir diff downward sw flux multiplied by timestep | W m-2 s       |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%nvisbm_cpl     | cumulative_surface_net_downward_direct_ultraviolet_and_visible_shortwave_flux_for_coupling_multiplied_by_timestep | cumulative net uv+vis beam downward sw rad flux multiplied by timestep | W m-2 s       |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%nvisdf_cpl     | cumulative_surface_net_downward_diffuse_ultraviolet_and_visible_shortwave_flux_for_coupling_multiplied_by_timestep | cumulative net uv+vis diff downward sw rad flux multiplied by timestep | W m-2 s       |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%dusfci_cpl     | instantaneous_surface_x_momentum_flux_for_coupling                                        | instantaneous sfc x momentum flux                    | Pa            |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%dvsfci_cpl     | instantaneous_surface_y_momentum_flux_for_coupling                                        | instantaneous sfc y momentum flux                    | Pa            |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%dtsfci_cpl     | instantaneous_surface_upward_sensible_heat_flux_for_coupling                              | instantaneous sfc sensible heat flux                 | W m-2         |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%dqsfci_cpl     | instantaneous_surface_upward_latent_heat_flux_for_coupling                                | instantaneous sfc latent heat flux                   | W m-2         |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%dlwsfci_cpl    | instantaneous_surface_downwelling_longwave_flux_for_coupling                              | instantaneous sfc downward lw flux                   | W m-2         |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%dswsfci_cpl    | instantaneous_surface_downwelling_shortwave_flux_for_coupling                             | instantaneous sfc downward sw flux                   | W m-2         |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%dnirbmi_cpl    | instantaneous_surface_downwelling_direct_near_infrared_shortwave_flux_for_coupling        | instantaneous sfc nir beam downward sw flux          | W m-2         |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%dnirdfi_cpl    | instantaneous_surface_downwelling_diffuse_near_infrared_shortwave_flux_for_coupling       | instantaneous sfc nir diff downward sw flux          | W m-2         |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%dvisbmi_cpl    | instantaneous_surface_downwelling_direct_ultraviolet_and_visible_shortwave_flux_for_coupling | instantaneous sfc uv+vis beam downward sw flux       | W m-2         |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%dvisdfi_cpl    | instantaneous_surface_downwelling_diffuse_ultraviolet_and_visible_shortwave_flux_for_coupling | instantaneous sfc uv+vis diff downward sw flux       | W m-2         |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%nlwsfci_cpl    | instantaneous_surface_net_downward_longwave_flux_for_coupling                             | instantaneous net sfc downward lw flux               | W m-2         |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%nswsfci_cpl    | instantaneous_surface_net_downward_shortwave_flux_for_coupling                            | instantaneous net sfc downward sw flux               | W m-2         |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%nnirbmi_cpl    | instantaneous_surface_net_downward_direct_near_infrared_shortwave_flux_for_coupling       | instantaneous net nir beam sfc downward sw flux      | W m-2         |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%nnirdfi_cpl    | instantaneous_surface_net_downward_diffuse_near_infrared_shortwave_flux_for_coupling      | instantaneous net nir diff sfc downward sw flux      | W m-2         |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%nvisbmi_cpl    | instantaneous_surface_net_downward_direct_ultraviolet_and_visible_shortwave_flux_for_coupling  | instantaneous net uv+vis beam downward sw flux         | W m-2         |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%nvisdfi_cpl    | instantaneous_surface_net_downward_diffuse_ultraviolet_and_visible_shortwave_flux_for_coupling | instantaneous net uv+vis diff downward sw flux         | W m-2         |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%t2mi_cpl       | instantaneous_temperature_at_2m_for_coupling                                                   | instantaneous T2m                                      | K             |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%q2mi_cpl       | instantaneous_specific_humidity_at_2m_for_coupling                                             | instantaneous Q2m                                      | kg kg-1       |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%u10mi_cpl      | instantaneous_x_wind_at_10m_for_coupling                                                       | instantaneous U10m                                     | m s-1         |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%v10mi_cpl      | instantaneous_y_wind_at_10m_for_coupling                                                       | instantaneous V10m                                     | m s-1         |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%tsfci_cpl      | instantaneous_surface_skin_temperature_for_coupling                                            | instantaneous sfc temperature                          | K             |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%psurfi_cpl     | instantaneous_surface_air_pressure_for_coupling                                                | instantaneous sfc pressure                             | Pa            |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%oro_cpl        |                                                                                                | orography                                              |               |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%slmsk_cpl      |                                                                                                | land/sea/ice mask                                      |               |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%tconvtend      | tendency_of_air_temperature_due_to_deep_convection_for_coupling_on_physics_timestep            | tendency of air temperature due to deep convection     | K             |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%qconvtend      | tendency_of_water_vapor_specific_humidity_due_to_deep_convection_for_coupling_on_physics_timestep | tendency of specific humidity due to deep convection| kg kg-1       |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%uconvtend      | tendency_of_x_wind_due_to_deep_convection_for_coupling_on_physics_timestep                     | tendency_of_x_wind_due_to_deep_convection              | m s-1         |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%vconvtend      | tendency_of_y_wind_due_to_deep_convection_for_coupling_on_physics_timestep                     | tendency_of_y_wind_due_to_deep_convection              | m s-1         |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%ca_out         |                                                                                                |                                                        |               |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%ca_deep        | fraction_of_cellular_automata_for_deep_convection                                              | fraction of cellular automata for deep convection      | frac          |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%ca_turb        |                                                                                                |                                                        |               |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%ca_shal        |                                                                                                |                                                        |               |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%ca_rad         |                                                                                                |                                                        |               |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%ca_micro       |                                                                                                |                                                        |               |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%cape           | convective_available_potential_energy_for_coupling                                             | convective available potential energy for coupling DH* CHECK THIS DOESN'T MAKE SENSE!!! *DH | m2 s-2        |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%shum_wts       | weights_for_stochastic_shum_perturbation                                                       | weights for stochastic shum perturbation               | none          |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%sppt_wts       | weights_for_stochastic_sppt_perturbation                                                       | weights for stochastic sppt perturbation               | none          |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%skebu_wts      | weights_for_stochastic_skeb_perturbation_of_x_wind                                             | weights for stochastic skeb perturbation of x wind     | none          |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%skebv_wts      | weights_for_stochastic_skeb_perturbation_of_y_wind                                             | weights for stochastic skeb perturbation of y wind     | none          |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%sfc_wts        | weights_for_stochastic_surface_physics_perturbation                                            | weights for stochastic surface physics perturbation    | none          |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%nsfcpert       |                                                                                                | number of surface perturbations (redundant? Model%...) |               |    0 | integer |           | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%vcu_wts        |                                                                                                |                                                        |               |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%vcv_wts        |                                                                                                |                                                        |               |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%dqdti          | instantaneous_water_vapor_specific_humidity_tendency_due_to_convection                         | instantaneous moisture tendency due to convection      | kg kg-1 s-1   |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%cnvqci         | instantaneous_deep_convective_cloud_condensate_mixing_ratio_on_dynamics_time_step              | instantaneous total convective condensate mixing ratio | kg kg-1       |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%upd_mfi        | instantaneous_atmosphere_updraft_convective_mass_flux_on_dynamics_timestep                     | (updraft mass flux) * delt                             | kg m-2        |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%dwn_mfi        | instantaneous_atmosphere_downdraft_convective_mass_flux_on_dynamics_timestep                   | (downdraft mass flux) * delt                           | kg m-2        |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%det_mfi        | instantaneous_atmosphere_detrainment_convective_mass_flux_on_dynamics_timestep                 | (detrainment mass flux) * delt                         | kg m-2        |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%cldcovi        |                                                                                                | instantaneous 3D cloud fraction                        |               |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%nwfa2d         | tendency_of_water_friendly_aerosols_at_surface                                                 | instantaneous water-friendly sfc aerosol source        | kg-1 s-1      |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%nifa2d         | tendency_of_ice_friendly_aerosols_at_surface                                                   | instantaneous ice-friendly sfc aerosol source          | kg-1 s-1      |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%ushfsfci       | instantaneous_upward_sensible_heat_flux                                                        | instantaneous upward sensible heat flux                | W m-2         |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Coupling%dkt            | instantaneous_atmosphere_heat_diffusivity                                                      | instantaneous atmospheric heat diffusivity             | m2 s-1        |    2 | real    | kind_phys | none   | F        |
!!
#endif
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
    real (kind=kind_phys), pointer :: tseain_cpl(:)  => null()   !< aoi_fld%tseain(item,lan)
    real (kind=kind_phys), pointer :: tisfcin_cpl(:) => null()   !< aoi_fld%tisfcin(item,lan)
    real (kind=kind_phys), pointer :: ficein_cpl(:)  => null()   !< aoi_fld%ficein(item,lan)
    real (kind=kind_phys), pointer :: hicein_cpl(:)  => null()   !< aoi_fld%hicein(item,lan)
    real (kind=kind_phys), pointer :: hsnoin_cpl(:)  => null()   !< aoi_fld%hsnoin(item,lan)
    !--- only variable needed for cplwav=.TRUE.
    !--- also needed for ice/ocn coupling - Xingren
    real (kind=kind_phys), pointer :: slimskin_cpl(:)=> null()   !< aoi_fld%slimskin(item,lan)

!--- outgoing accumulated quantities
    real (kind=kind_phys), pointer :: rain_cpl  (:)  => null()   !< total rain precipitation
    real (kind=kind_phys), pointer :: rainc_cpl (:)  => null()   !< convective rain precipitation
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

    !--- cellular automata
    real (kind=kind_phys), pointer :: tconvtend(:,:) => null()
    real (kind=kind_phys), pointer :: qconvtend(:,:) => null()
    real (kind=kind_phys), pointer :: uconvtend(:,:) => null()
    real (kind=kind_phys), pointer :: vconvtend(:,:) => null()
    real (kind=kind_phys), pointer :: ca_out   (:)   => null()  !
    real (kind=kind_phys), pointer :: ca_deep  (:)   => null()  !
    real (kind=kind_phys), pointer :: ca_turb  (:)   => null()  !
    real (kind=kind_phys), pointer :: ca_shal  (:)   => null()  !
    real (kind=kind_phys), pointer :: ca_rad   (:)   => null()  !
    real (kind=kind_phys), pointer :: ca_micro (:)   => null()  !
    real (kind=kind_phys), pointer :: cape     (:)   => null()  !
    
    !--- stochastic physics
    real (kind=kind_phys), pointer :: shum_wts  (:,:)   => null()  !
    real (kind=kind_phys), pointer :: sppt_wts  (:,:)   => null()  !
    real (kind=kind_phys), pointer :: skebu_wts (:,:)   => null()  !
    real (kind=kind_phys), pointer :: skebv_wts (:,:)   => null()  !
    real (kind=kind_phys), pointer :: sfc_wts   (:,:)   => null()  ! mg, sfc-perts
    integer                        :: nsfcpert=6                   !< number of sfc perturbations

!--- instantaneous quantities for GoCart and will be accumulated for 3D diagnostics
    real (kind=kind_phys), pointer :: dqdti   (:,:)   => null()  !< instantaneous total moisture tendency (kg/kg/s)
    real (kind=kind_phys), pointer :: cnvqci  (:,:)   => null()  !< instantaneous total convective conensate (kg/kg)
    real (kind=kind_phys), pointer :: upd_mfi (:,:)   => null()  !< instantaneous convective updraft mass flux
    real (kind=kind_phys), pointer :: dwn_mfi (:,:)   => null()  !< instantaneous convective downdraft mass flux
    real (kind=kind_phys), pointer :: det_mfi (:,:)   => null()  !< instantaneous convective detrainment mass flux
    real (kind=kind_phys), pointer :: cldcovi (:,:)   => null()  !< instantaneous 3D cloud fraction
    real (kind=kind_phys), pointer :: nwfa2d  (:)     => null()  !< instantaneous water-friendly sfc aerosol source
    real (kind=kind_phys), pointer :: nifa2d  (:)     => null()  !< instantaneous ice-friendly sfc aerosol source

    !--- instantaneous quantities for GSDCHEM coupling
    real (kind=kind_phys), pointer :: ushfsfci(:)     => null()  !< instantaneous upward sensible heat flux (w/m**2)
    real (kind=kind_phys), pointer :: dkt     (:,:)   => null()  !< instantaneous dkt diffusion coefficient for temperature (m**2/s)

    contains
      procedure :: create  => coupling_create  !<   allocate array data
  end type GFS_coupling_type


!----------------------------------------------------------------------------------
! GFS_control_type
!   model control parameters input from a namelist and/or derived from others
!   list of those that can be modified during the run are at the bottom of the list
!----------------------------------------------------------------------------------
#if 0
!! \section arg_table_GFS_control_type
!! | local_name                           | standard_name                                                                 | long_name                                               | units         | rank | type      |    kind   | intent | optional |
!! |--------------------------------------|-------------------------------------------------------------------------------|---------------------------------------------------------|---------------|------|-----------|-----------|--------|----------|
!! | GFS_Control%me                       | mpi_rank                                                                      | current MPI-rank                                        | index         |    0 | integer   |           | none   | F        |
!! | GFS_Control%master                   | mpi_root                                                                      | master MPI-rank                                         | index         |    0 | integer   |           | none   | F        |
!! | GFS_Control%communicator             | mpi_comm                                                                      | MPI communicator                                        | index         |    0 | integer   |           | none   | F        |
!! | GFS_Control%ntasks                   | mpi_size                                                                      | number of MPI tasks in communicator                     | count         |    0 | integer   |           | none   | F        |
!! | GFS_Control%nthreads                 | omp_threads                                                                   | number of OpenMP threads available for physics schemes  | count         |    0 | integer   |           | none   | F        |
!! | GFS_Control%nlunit                   | iounit_namelist                                                               | fortran unit number for file opens                      | none          |    0 | integer   |           | none   | F        |
!! | GFS_Control%fn_nml                   | namelist_filename                                                             | namelist filename                                       | none          |    0 | character | len=64    | none   | F        |
!! | GFS_Control%input_nml_file           | namelist_filename_for_internal_file_reads                                     | namelist filename for internal file reads               | none          |    1 | character | len=256   | none   | F        |
!! | GFS_Control%logunit                  | iounit_log                                                                    | fortran unit number for logfile                         | none          |    0 | integer   |           | none   | F        |
!! | GFS_Control%fhzero                   |                                                                               | seconds between clearing of diagnostic buckets          | s             |    0 | real      | kind_phys | none   | F        |
!! | GFS_Control%ldiag3d                  | flag_diagnostics_3D                                                           | flag for 3d diagnostic fields                           | flag          |    0 | logical   |           | none   | F        |
!! | GFS_Control%lssav                    | flag_diagnostics                                                              | logical flag for storing diagnostics                    | flag          |    0 | logical   |           | none   | F        |
!! | GFS_Control%fhcyc                    |                                                                               | frequency for surface data cycling (secs)               | s             |    0 | real      | kind_phys | none   | F        |
!! | GFS_Control%lgocart                  | flag_gocart                                                                   | flag for 3d diagnostic fields for gocart 1              | flag          |    0 | logical   |           | none   | F        |
!! | GFS_Control%fhgoc3d                  |                                                                               | seconds between calls to gocart                         | s             |    0 | real      | kind_phys | none   | F        |
!! | GFS_Control%thermodyn_id             |                                                                               | valid for GFS only for get_prs/phi                      | index         |    0 | integer   |           | none   | F        |
!! | GFS_Control%sfcpress_id              |                                                                               | valid for GFS only for get_prs/phi                      | index         |    0 | integer   |           | none   | F        |
!! | GFS_Control%gen_coord_hybrid         |                                                                               | flag for Henry's gen coord                              | flag          |    0 | logical   |           | none   | F        |
!! | GFS_Control%isc                      |                                                                               | starting i-index for this MPI-domain                    | index         |    0 | integer   |           | none   | F        |
!! | GFS_Control%jsc                      |                                                                               | starting j-index for this MPI-domain                    | index         |    0 | integer   |           | none   | F        |
!! | GFS_Control%nx                       |                                                                               | number of points in i-dir for this MPI rank             | count         |    0 | integer   |           | none   | F        |
!! | GFS_Control%ny                       |                                                                               | number of points in j-dir for this MPI rank             | count         |    0 | integer   |           | none   | F        |
!! | GFS_Control%levs                     | vertical_dimension                                                            | number of vertical levels                               | count         |    0 | integer   |           | none   | F        |
!! | GFS_Control%ak                       |                                                                               | a parameter for sigma pressure level calculations       | Pa            |    1 | real      | kind_phys | none   | F        |
!! | GFS_Control%bk                       |                                                                               | b parameter for sigma pressure level calculations       | none          |    1 | real      | kind_phys | none   | F        |
!! | GFS_Control%cnx                      |                                                                               | number of points in i-dir for this cubed-sphere face    | count         |    0 | integer   |           | none   | F        |
!! | GFS_Control%cny                      |                                                                               | number of points in j-dir for this cubed-sphere face    | count         |    0 | integer   |           | none   | F        |
!! | GFS_Control%lonr                     | number_of_equatorial_longitude_points                                         | number of global points in x-dir (i) along the equator  | count         |    0 | integer   |           | none   | F        |
!! | GFS_Control%latr                     |                                                                               | number of global points in y-dir (j) along the meridian | count         |    0 | integer   |           | none   | F        |
!! | GFS_Control%blksz                    | horizontal_block_size                                                         | for explicit data blocking: block sizes of all blocks   | count         |    1 | integer   |           | none   | F        |
!! | GFS_Control%blksz(cdata%blk_no)      | horizontal_loop_extent                                                        | horizontal loop extent                                  | count         |    0 | integer   |           | none   | F        |
!! | GFS_Control%blksz(cdata%blk_no)      | horizontal_dimension                                                          | horizontal dimension                                    | count         |    0 | integer   |           | none   | F        |
!! | GFS_Control%tile_num                 | number_of_tile                                                                | tile number                                             | none          |    0 | integer   |           | none   | F        |
!! | GFS_Control%cplflx                   | flag_for_flux_coupling                                                        | flag controlling cplflx collection (default off)        | flag          |    0 | logical   |           | none   | F        |
!! | GFS_Control%cplwav                   | flag_for_wave_coupling                                                        | flag controlling cplwav collection (default off)        | flag          |    0 | logical   |           | none   | F        |
!! | GFS_Control%cplchm                   | flag_for_chemistry_coupling                                                   | flag controlling cplchm collection (default off)        | flag          |    0 | logical   |           | none   | F        |
!! | GFS_Control%lsidea                   | flag_idealized_physics                                                        | flag for idealized physics                              | flag          |    0 | logical   |           | none   | F        |
!! | GFS_Control%dtp                      | time_step_for_physics                                                         | physics timestep                                        | s             |    0 | real      | kind_phys | none   | F        |
!! | GFS_Control%dtf                      | time_step_for_dynamics                                                        | dynamics timestep                                       | s             |    0 | real      | kind_phys | none   | F        |
!! | GFS_Control%nscyc                    |                                                                               | trigger for surface data cycling                        |               |    0 | integer   |           | none   | F        |
!! | GFS_Control%nszero                   |                                                                               | trigger for zeroing diagnostic buckets                  |               |    0 | integer   |           | none   | F        |
!! | GFS_Control%idat                     | date_and_time_at_model_initialization                                         | initialization date and time                            | none          |    1 | integer   |           | none   | F        |
!! | GFS_Control%idate                    | date_and_time_at_model_initialization_reordered                               | initial date with different size and ordering           | none          |    1 | integer   |           | none   | F        |
!! | GFS_Control%fhswr                    | frequency_for_shortwave_radiation                                             | frequency for shortwave radiation                       | s             |    0 | real      | kind_phys | none   | F        |
!! | GFS_Control%fhlwr                    | frequency_for_longwave_radiation                                              | frequency for longwave radiation                        | s             |    0 | real      | kind_phys | none   | F        |
!! | GFS_Control%nsswr                    |                                                                               | integer trigger for shortwave radiation                 |               |    0 | integer   |           | none   | F        |
!! | GFS_Control%nslwr                    |                                                                               | integer trigger for longwave  radiation                 |               |    0 | integer   |           | none   | F        |
!! | GFS_Control%levr                     | number_of_vertical_layers_for_radiation_calculations                          | number of vertical levels for radiation calculations    | count         |    0 | integer   |           | none   | F        |
!! | GFS_Control%nfxr                     |                                                                               | second dimension for fluxr diagnostic variable          |               |    0 | integer   |           | none   | F        |
!! | GFS_Control%aero_in                  | flag_for_aerosol_input_MG                                                     | flag for using aerosols in Morrison-Gettelman microphysics | flag       |    0 | logical   |           | none   | F        |
!! | GFS_Control%lmfshal                  |                                                                               | parameter for radiation                                 |               |    0 | logical   |           | none   | F        |
!! | GFS_Control%lmfdeep2                 |                                                                               | parameter for radiation                                 |               |    0 | logical   |           | none   | F        |
!! | GFS_Control%nrcm                     | array_dimension_of_random_number                                              | second dimension of random number stream for RAS        | count         |    0 | integer   |           | none   | F        |
!! | GFS_Control%iflip                    | flag_for_vertical_index_direction_control                                     | iflip - is not the same as flipv                        | flag          |    0 | integer   |           | none   | F        |
!! | GFS_Control%isol                     | flag_for_solar_constant                                                       | use prescribed solar constant                           | flag          |    0 | integer   |           | none   | F        |
!! | GFS_Control%ico2                     | flag_for_using_prescribed_global_mean_co2_value                               | prescribed global mean value (old opernl)               | flag          |    0 | integer   |           | none   | F        |
!! | GFS_Control%ialb                     | flag_for_using_climatology_albedo                                             | flag for using climatology alb, based on sfc type       | flag          |    0 | integer   |           | none   | F        |
!! | GFS_Control%iems                     | flag_for_surface_emissivity_control                                           | surface emissivity control flag, use fixed value of 1   | flag          |    0 | integer   |           | none   | F        |
!! | GFS_Control%iaer                     | flag_for_default_aerosol_effect_in_shortwave_radiation                        | default aerosol effect in sw only                       | flag          |    0 | integer   |           | none   | F        |
!! | GFS_Control%icliq_sw                 | flag_for_optical_property_for_liquid_clouds_for_shortwave_radiation           | sw optical property for liquid clouds                   | flag          |    0 | integer   |           | none   | F        |
!! | GFS_Control%iovr_sw                  | flag_for_max-random_overlap_clouds_for_shortwave_radiation                    | sw: max-random overlap clouds                           | flag          |    0 | integer   |           | none   | F        |
!! | GFS_Control%iovr_lw                  | flag_for_max-random_overlap_clouds_for_longwave_radiation                     | lw: max-random overlap clouds                           | flag          |    0 | integer   |           | none   | F        |
!! | GFS_Control%ictm                     | flag_for_initial_time-date_control                                            | flag for initial conditions and forcing                 | flag          |    0 | integer   |           | none   | F        |
!! | GFS_Control%isubc_sw                 | flag_for_sw_clouds_without_sub-grid_approximation                             | flag for sw clouds without sub-grid approximation       | flag          |    0 | integer   |           | none   | F        |
!! | GFS_Control%isubc_lw                 | flag_for_lw_clouds_without_sub-grid_approximation                             | flag for lw clouds without sub-grid approximation       | flag          |    0 | integer   |           | none   | F        |
!! | GFS_Control%crick_proof              | flag_for_CRICK-proof_cloud_water                                              | flag for CRICK-Proof cloud water                        | flag          |    0 | logical   |           | none   | F        |
!! | GFS_Control%ccnorm                   | flag_for_cloud_condensate_normalized_by_cloud_cover                           | flag for cloud condensate normalized by cloud cover     | flag          |    0 | logical   |           | none   | F        |
!! | GFS_Control%norad_precip             | flag_for_precipitation_effect_on_radiation                                    | radiation precip flag for Ferrier/Moorthi               | flag          |    0 | logical   |           | none   | F        |
!! | GFS_Control%lwhtr                    | flag_for_output_of_longwave_heating_rate                                      | flag to output lw heating rate (Radtend%lwhc)           | flag          |    0 | logical   |           | none   | F        |
!! | GFS_Control%swhtr                    | flag_for_output_of_shortwave_heating_rate                                     | flag to output sw heating rate (Radtend%swhc)           | flag          |    0 | logical   |           | none   | F        |
!! | GFS_Control%ncld                     | number_of_hydrometeors                                                        | choice of cloud scheme / number of hydrometeors         | count         |    0 | integer   |           | none   | F        |
!! | GFS_Control%imp_physics              | flag_for_microphysics_scheme                                                  | choice of microphysics scheme                           | flag          |    0 | integer   |           | none   | F        |
!! | GFS_Control%imp_physics_gfdl         | flag_for_gfdl_microphysics_scheme                                             | choice of GFDL microphysics scheme                      | flag          |    0 | integer   |           | none   | F        |
!! | GFS_Control%imp_physics_thompson     | flag_for_thompson_microphysics_scheme                                         | choice of Thompson microphysics scheme                  | flag          |    0 | integer   |           | none   | F        |
!! | GFS_Control%imp_physics_wsm6         | flag_for_wsm6_microphysics_scheme                                             | choice of WSM6 microphysics scheme                      | flag          |    0 | integer   |           | none   | F        |
!! | GFS_Control%imp_physics_zhao_carr    | flag_for_zhao_carr_microphysics_scheme                                        | choice of Zhao-Carr microphysics scheme                 | flag          |    0 | integer   |           | none   | F        |
!! | GFS_Control%imp_physics_zhao_carr_pdf| flag_for_zhao_carr_pdf_microphysics_scheme                                    | choice of Zhao-Carr microphysics scheme with PDF clouds | flag          |    0 | integer   |           | none   | F        |
!! | GFS_Control%imp_physics_mg           | flag_for_morrison_gettelman_microphysics_scheme                               | choice of Morrison-Gettelman rmicrophysics scheme       | flag          |    0 | integer   |           | none   | F        |
!! | GFS_Control%psautco                  | coefficient_from_cloud_ice_to_snow                                            | auto conversion coeff from ice to snow                  | none          |    1 | real      | kind_phys | none   | F        |
!! | GFS_Control%prautco                  | coefficient_from_cloud_water_to_rain                                          | auto conversion coeff from cloud to rain                | none          |    1 | real      | kind_phys | none   | F        |
!! | GFS_Control%evpco                    | coefficient_for_evaporation_of_rainfall                                       | coeff for evaporation of largescale rain                | none          |    0 | real      | kind_phys | none   | F        |
!! | GFS_Control%wminco                   | cloud_condensed_water_conversion_threshold                                    | water and ice minimum threshold for Zhao                | none          |    1 | real      | kind_phys | none   | F        |
!! | GFS_Control%avg_max_length           | time_interval_for_maximum_hourly_fields                                       | reset time interval for maximum hourly fields           | s             |    0 | real      | kind_phys | none   | F        |
!! | GFS_Control%fprcp                    | number_of_frozen_precipitation_species                                        | number of frozen precipitation species                  | count         |    0 | integer   |           | none   | F        |
!! | GFS_Control%pdfflag                  | flag_for_pdf_for_morrison_gettelman_microphysics_scheme                       | pdf flag for MG macrophysics                            | flag          |    0 | integer   |           | none   | F        |
!! | GFS_Control%mg_dcs                   | mg_autoconversion_size_threshold_ice_snow                                     | autoconversion size threshold for cloud ice to snow for MG microphysics       | um      |    0 | real       | kind_phys | none   | F        |
!! | GFS_Control%mg_qcvar                 | mg_cloud_water_variance                                                       | cloud water relative variance for MG microphysics                             |         |    0 | real       | kind_phys | none   | F        |
!! | GFS_Control%mg_ts_auto_ice           | mg_time_scale_for_autoconversion_of_ice                                       | autoconversion time scale for ice for MG microphysics                         | s       |    1 | real       | kind_phys | none   | F        |
!! | GFS_Control%mg_rhmini                | mg_minimum_rh_for_ice                                                         | relative humidity threshold parameter for nucleating ice for MG microphysics  | none    |    0 | real       | kind_phys | none   | F        |
!! | GFS_Control%mg_ncnst                 | mg_drop_concentration_constant                                                | droplet concentration constant for MG microphysics                            | m-3     |    0 | real       | kind_phys | none   | F        |
!! | GFS_Control%mg_ninst                 | mg_ice_concentration_constant                                                 | ice concentration constant for MG microphysics                                | m-3     |    0 | real       | kind_phys | none   | F        |
!! | GFS_Control%mg_ngnst                 | mg_graupel_concentration_constant                                             | graupel concentration constant for MG microphysics                            | m-3     |    0 | real       | kind_phys | none   | F        |
!! | GFS_Control%mg_berg_eff_factor       | mg_bergeron_efficiency_factor                                                 | bergeron efficiency factor for MG microphysics                                | frac    |    0 | real       | kind_phys | none   | F        |
!! | GFS_Control%mg_alf                   | mg_tuning_factor_for_alphas                                                   | tuning factor for alphas (alpha = 1 - critical relative humidity)             | none    |    0 | real       | kind_phys | none   | F        |
!! | GFS_Control%mg_qcmin                 | mg_minimum_cloud_condensed_water_and_ice_mixing_ratio                         | minimum cloud condensed water and ice mixing ratio in MG macro clouds         | kg kg-1 |    1 | real       | kind_phys | none   | F        |
!! | GFS_Control%mg_qcmin(1)              | mg_minimum_cloud_condensed_water_mixing_ratio                                 | minimum cloud condensed water mixing ratio in MG macro clouds                 | kg kg-1 |    0 | real       | kind_phys | none   | F        |
!! | GFS_Control%mg_qcmin(2)              | mg_minimum_ice_mixing_ratio                                                   | minimum ice mixing ratio in MG macro clouds                                   | kg kg-1 |    0 | real       | kind_phys | none   | F        |
!! | GFS_Control%mg_precip_frac_method    | mg_type_of_precip_fraction_method                                             | type of precip fraction method for MG microphysics (in_cloud or max_overlap)  | none    |    0 | character  | len=16    | none   | F        |
!! | GFS_Control%tf                       | frozen_cloud_threshold_temperature                                            | threshold temperature below which all cloud is ice                            | K       |    0 | real       | kind_phys | none   | F        |
!! | GFS_Control%tcr                      | cloud_phase_transition_threshold_temperature                                  | threshold temperature below which cloud starts to freeze                      | K       |    0 | real       | kind_phys | none   | F        |
!! | GFS_Control%tcrf                     | cloud_phase_transition_denominator                                            | denominator in cloud phase transition = 1/(tcr-tf)                            | K-1     |    0 | real       | kind_phys | none   | F        |
!! | GFS_Control%effr_in                  | flag_for_cloud_effective_radii                                                | flag for cloud effective radii calculations in microphysics                   |         |    0 | logical    |           | none   | F        |
!! | GFS_Control%microp_uniform           | mg_flag_for_uniform_subcolumns                                                | flag for uniform subcolumns for MG microphysics                               | flag    |    0 | logical    |           | none   | F        |
!! | GFS_Control%do_cldliq                |                                                                               |                                                                               |         |    0 | logical    |           | none   | F        |
!! | GFS_Control%do_cldice                | mg_flag_for_cloud_ice_processes                                               | flag for cloud ice processes for MG microphysics                              | flag    |    0 | logical    |           | none   | F        |
!! | GFS_Control%hetfrz_classnuc          | mg_flag_for_heterogeneous_freezing                                            | flag for heterogeneous freezing for MG microphysics                           | flag    |    0 | logical    |           | none   | F        |
!! | GFS_Control%mg_nccons                | mg_flag_drop_concentration_constant                                           | flag for constant droplet concentration for MG microphysics                   | flag    |    0 | logical    |           | none   | F        |
!! | GFS_Control%mg_nicons                | mg_flag_ice_concentration_constant                                            | flag for constant ice concentration for MG microphysics                       | flag    |    0 | logical    |           | none   | F        |
!! | GFS_Control%mg_ngcons                | mg_flag_graupel_concentration_constant                                        | flag for constant graupel concentration for MG microphysics                   | flag    |    0 | logical    |           | none   | F        |
!! | GFS_Control%sed_supersat             | mg_allow_supersat_after_sed                                                   | allow supersaturation after sedimentation for MG microphysics                 | flag    |    0 | logical    |           | none   | F        |
!! | GFS_Control%do_sb_physics            | mg_flag_for_sb2001_autoconversion                                             | flag for SB 2001 autoconversion or accretion for MG microphysics              | flag    |    0 | logical    |           | none   | F        |
!! | GFS_Control%mg_do_graupel            | mg_flag_for_graupel                                                           | flag for graupel for MG microphysics (hail possible if false)                 | flag    |    0 | logical    |           | none   | F        |
!! | GFS_Control%mg_do_hail               | mg_flag_for_hail                                                              | flag for hail for MG microphysics (graupel possible if false)                 | flag    |    0 | logical    |           | none   | F        |
!! | GFS_Control%mg_do_ice_gmao           | mg_flag_for_gmao_ice_formulation                                              | flag for gmao ice formulation                                                 | flag    |    0 | logical    |           | none   | F        |
!! | GFS_Control%mg_do_liq_liu            | mg_flag_for_liu_liquid_treatment                                              | flag for liu liquid treatment                                                 | flag    |    0 | logical    |           | none   | F        |
!! | GFS_Control%shoc_parm(1)             | shoc_tke_dissipatation_pressure_threshold                                     | pressure below which extra TKE diss. is applied in SHOC | Pa            |    0 | real      | kind_phys | none   | F        |
!! | GFS_Control%shoc_parm(2)             | shoc_tke_dissipation_tunable_parameter                                        | mult. tuning parameter for TKE diss. in SHOC            | none          |    0 | real      | kind_phys | none   | F        |
!! | GFS_Control%shoc_parm(3)             | shoc_tke_dissipation_tunable_parameter_near_surface                           | mult. tuning parameter for TKE diss. at surface in SHOC | none          |    0 | real      | kind_phys | none   | F        |
!! | GFS_Control%shoc_parm(4)             | shoc_implicit_TKE_integration_uncentering_term                                | uncentering term for TKE integration in SHOC            | none          |    0 | real      | kind_phys | none   | F        |
!! | GFS_Control%shoc_parm(5)             | shoc_flag_for_optional_surface_TKE_dissipation                                | flag for alt. TKE diss. near surface in SHOC (>0 = ON)  | none          |    0 | real      | kind_phys | none   | F        |
!! | GFS_Control%ncnd                     | number_of_cloud_condensate_types                                              | number of cloud condensate types                        | count         |    0 | integer   |           | none   | F        |
!! | GFS_Control%ltaerosol                | flag_for_aerosol_physics                                                      | flag for aerosol physics                                | flag          |    0 | logical   |           | none   | F        |
!! | GFS_Control%lradar                   | flag_for_radar_reflectivity                                                   | flag for radar reflectivity                             | flag          |    0 | logical   |           | none   | F        |
!! | GFS_Control%ttendlim                 | limit_for_temperature_tendency_for_microphysics                               | temperature tendency limiter per physics time step      | K s-1         |    0 | real      | kind_phys | none   | F        |
!! | GFS_Control%lgfdlmprad               |                                                                               | flag for GFDL mp scheme and radiation consistency       |               |    0 | logical   |           | none   | F        |
!! | GFS_Control%lsm                      | flag_for_land_surface_scheme                                                  | flag for land surface model                             | flag          |    0 | integer   |           | none   | F        |
!! | GFS_Control%lsm_noah                 | flag_for_noah_land_surface_scheme                                             | flag for NOAH land surface model                        | flag          |    0 | integer   |           | none   | F        |
!! | GFS_Control%lsm_ruc                  | flag_for_ruc_land_surface_scheme                                              | flag for RUC land surface model                         | flag          |    0 | integer   |           | none   | F        |
!! | GFS_Control%lsoil                    | soil_vertical_dimension                                                       | number of soil layers                                   | count         |    0 | integer   |           | none   | F        |
!! | GFS_Control%lsoil_lsm                | soil_vertical_dimension_for_land_surface_model                                | number of soil layers for land surface model            | count         |    0 | integer   |           | none   | F        |
!! | GFS_Control%ivegsrc                  | vegetation_type_dataset_choice                                                | land use dataset choice                                 | index         |    0 | integer   |           | none   | F        |
!! | GFS_Control%isot                     | soil_type_dataset_choice                                                      | soil type dataset choice                                | index         |    0 | integer   |           | none   | F        |
!! | GFS_Control%mom4ice                  | flag_for_mom4_coupling                                                        | flag controls mom4 sea ice                              | flag          |    0 | logical   |           | none   | F        |
!! | GFS_Control%use_ufo                  |                                                                               | flag for gcycle surface option                          |               |    0 | logical   |           | none   | F        |
!! | GFS_Control%ras                      | flag_for_ras_deep_convection                                                  | flag for ras convection scheme                          | flag          |    0 | logical   |           | none   | F        |
!! | GFS_Control%flipv                    | flag_flip                                                                     | vertical flip logical                                   | flag          |    0 | logical   |           | none   | F        |
!! | GFS_Control%trans_trac               | flag_for_convective_transport_of_tracers                                      | flag for convective transport of tracers                | flag          |    0 | logical   |           | none   | F        |
!! | GFS_Control%old_monin                | flag_for_old_PBL_scheme                                                       | flag for using old PBL schemes                          | flag          |    0 | logical   |           | none   | F        |
!! | GFS_Control%cnvgwd                   | flag_convective_gravity_wave_drag                                             | flag for conv gravity wave drag                         | flag          |    0 | logical   |           | none   | F        |
!! | GFS_Control%mstrat                   | flag_for_moorthi_stratus                                                      | flag for moorthi approach for stratus                   | flag          |    0 | logical   |           | none   | F        |
!! | GFS_Control%moist_adj                |                                                                               | flag for moist convective adjustment                    |               |    0 | logical   |           | none   | F        |
!! | GFS_Control%cscnv                    | flag_for_Chikira_Sugiyama_deep_convection                                     | flag for Chikira-Sugiyama convection                    | flag          |    0 | logical   |           | none   | F        |
!! | GFS_Control%satmedmf                 | flag_for_scale_aware_TKE_moist_EDMF_PBL                                       | flag for scale-aware TKE moist EDMF PBL scheme          | flag          |    0 | logical   |           | none   | F        |
!! | GFS_Control%shinhong                 | flag_for_scale_aware_Shinhong_PBL                                             | flag for scale-aware Shinhong PBL scheme                | flag          |    0 | logical   |           | none   | F        |
!! | GFS_Control%do_ysu                   | flag_for_ysu                                                                  | flag for YSU PBL scheme                                 | flag          |    0 | logical   |           | none   | F        |
!! | GFS_Control%cal_pre                  | flag_for_precipitation_type_algorithm                                         | flag controls precip type algorithm                     | flag          |    0 | logical   |           | none   | F        |
!! | GFS_Control%do_aw                    | flag_for_Arakawa_Wu_adjustment                                                | flag for Arakawa Wu scale-aware adjustment              | flag          |    0 | logical   |           | none   | F        |
!! | GFS_Control%do_awdd                  | flag_arakawa_wu_downdraft                                                     | AW scale-aware option in cs convection downdraft        | flag          |    0 | logical   |           | none   | F        |
!! | GFS_Control%flx_form                 | flag_flux_form_CS                                                             | enable use of flux form of equations in CS scheme       | flag          |    0 | logical   |           | none   | F        |
!! | GFS_Control%do_shoc                  | flag_for_shoc                                                                 | flag for SHOC                                           | flag          |    0 | logical   |           | none   | F        |
!! | GFS_Control%shocaftcnv               | flag_for_shoc_after_convection                                                | flag to execute SHOC after convection                   | flag          |    0 | logical   |           | none   | F        |
!! | GFS_Control%shoc_cld                 |                                                                               | flag for clouds                                         |               |    0 | logical   |           | none   | F        |
!! | GFS_Control%uni_cld                  |                                                                               | flag for clouds in grrad                                |               |    0 | logical   |           | none   | F        |
!! | GFS_Control%oz_phys                  | flag_for_ozone_physics                                                        | flag for old (2006) ozone physics                       | flag          |    0 | logical   |           | none   | F        |
!! | GFS_Control%oz_phys_2015             | flag_for_2015_ozone_physics                                                   | flag for new (2015) ozone physics                       | flag          |    0 | logical   |           | none   | F        |
!! | GFS_Control%h2o_phys                 |                                                                               | flag for stratosphere h2o                               |               |    0 | logical   |           | none   | F        |
!! | GFS_Control%pdfcld                   |                                                                               | flag for pdfcld                                         |               |    0 | logical   |           | none   | F        |
!! | GFS_Control%shcnvcw                  | flag_shallow_convective_cloud                                                 | flag for shallow convective cloud                       |               |    0 | logical   |           | none   | F        |
!! | GFS_Control%redrag                   | flag_for_reduced_drag_coefficient_over_sea                                    | flag for reduced drag coeff. over sea                   | flag          |    0 | logical   |           | none   | F        |
!! | GFS_Control%hybedmf                  | flag_for_hedmf                                                                | flag for hybrid edmf pbl scheme (moninedmf)             | flag          |    0 | logical   |           | none   | F        |
!! | GFS_Control%dspheat                  | flag_TKE_dissipation_heating                                                  | flag for tke dissipative heating                        | flag          |    0 | logical   |           | none   | F        |
!! | GFS_Control%lheatstrg                | flag_for_canopy_heat_storage                                                  | flag for canopy heat storage parameterization           | flag          |    0 | logical   |           | none   | F        |
!! | GFS_Control%cnvcld                   |                                                                               |                                                         |               |    0 | logical   |           | none   | F        |
!! | GFS_Control%random_clds              |                                                                               | flag controls whether clouds are random                 |               |    0 | logical   |           | none   | F        |
!! | GFS_Control%shal_cnv                 | flag_for_shallow_convection                                                   | flag for calling shallow convection                     | flag          |    0 | logical   |           | none   | F        |
!! | GFS_Control%do_deep                  |                                                                               | whether to do deep convection                           |               |    0 | logical   |           | none   | F        |
!! | GFS_Control%imfshalcnv               | flag_for_mass_flux_shallow_convection_scheme                                  | flag for mass-flux shallow convection scheme            | flag          |    0 | integer   |           | none   | F        |
!! | GFS_Control%imfdeepcnv               | flag_for_mass_flux_deep_convection_scheme                                     | flag for mass-flux deep convection scheme               | flag          |    0 | integer   |           | none   | F        |
!! | GFS_Control%nmtvr                    | number_of_statistical_measures_of_subgrid_orography                           | number of topographic variables in GWD                  | count         |    0 | integer   |           | none   | F        |
!! | GFS_Control%jcap                     |                                                                               | number of spectral wave trancation                      |               |    0 | integer   |           | none   | F        |
!! | GFS_Control%cs_parm                  |                                                                               | tunable parameters for Chikira-Sugiyama convection      |               |    1 | real      | kind_phys | none   | F        |
!! | GFS_Control%cs_parm(1)               | updraft_velocity_tunable_parameter_1_CS                                       | tunable parameter 1 for Chikira-Sugiyama convection     | m s-1         |    0 | real      | kind_phys | none   | F        |
!! | GFS_Control%cs_parm(2)               | updraft_velocity_tunable_parameter_2_CS                                       | tunable parameter 2 for Chikira-Sugiyama convection     | m s-1         |    0 | real      | kind_phys | none   | F        |
!! | GFS_Control%cs_parm(3)               | detrainment_and_precipitation_tunable_parameter_3_CS   | partition water between detrainment and precipitation (decrease for more precipitation) | m    |    0 | real      | kind_phys | none   | F        |
!! | GFS_Control%cs_parm(4)               | detrainment_and_precipitation_tunable_parameter_4_CS   | partition water between detrainment and precipitation (decrease for more precipitation) | m    |    0 | real      | kind_phys | none   | F        |
!! | GFS_Control%cs_parm(9)               | entrainment_efficiency_tunable_parameter_9_CS                                 | entrainment efficiency                                  | none          |    0 | real      | kind_phys | none   | F        |
!! | GFS_Control%flgmin                   |                                                                               | [in] ice fraction bounds                                |               |    1 | real      | kind_phys | none   | F        |
!! | GFS_Control%cgwf                     | multiplication_factors_for_convective_gravity_wave_drag                       | multiplication factor for convective GWD                | none          |    1 | real      | kind_phys | none   | F        |
!! | GFS_Control%ccwf                     |                                                                               | multiplication factor for critical cloud workfunction   | none          |    1 | real      | kind_phys | none   | F        |
!! | GFS_Control%cdmbgwd                  | multiplication_factors_for_mountain_blocking_and_orographic_gravity_wave_drag | multiplication factors for cdmb and gwd                 | none          |    1 | real      | kind_phys | none   | F        |
!! | GFS_Control%sup                      | ice_supersaturation_threshold                                                 | ice supersaturation parameter for PDF clouds            | none          |    0 | real      | kind_phys | none   | F        |
!! | GFS_Control%ctei_rm                  | critical_cloud_top_entrainment_instability_criteria                           | critical cloud top entrainment instability criteria     | none          |    1 | real      | kind_phys | none   | F        |
!! | GFS_Control%crtrh                    | critical_relative_humidity_at_sfc_pbltop_toa                                  | critical relative humidity at SFC, PBL top and TOA      | frac          |    1 | real      | kind_phys | none   | F        |
!! | GFS_Control%dlqf                     |                                                                               | factor for cloud condensate detrainment from cloud edges|               |    1 | real      | kind_phys | none   | F        |
!! | GFS_Control%psauras                  |                                                                               | auto conversion coeff from ice to snow in ras           |               |    1 | real      | kind_phys | none   | F        |
!! | GFS_Control%prauras                  |                                                                               | auto conversion coeff from cloud to rain in ras         |               |    1 | real      | kind_phys | none   | F        |
!! | GFS_Control%wminras                  |                                                                               | water and ice minimum threshold for ras                 |               |    1 | real      | kind_phys | none   | F        |
!! | GFS_Control%seed0                    |                                                                               | random seed for radiation                               |               |    0 | integer   |           | none   | F        |
!! | GFS_Control%rbcr                     |                                                                               | Critical Richardson Number in the PBL scheme            |               |    0 | real      | kind_phys | none   | F        |
!! | GFS_Control%prslrd0                  | pressure_cutoff_for_rayleigh_damping                                          | pressure level from which Rayleigh Damping is applied   | Pa            |    0 | real      | kind_phys | none   | F        |
!! | GFS_Control%ral_ts                   | time_scale_for_rayleigh_damping                                               | time scale for Rayleigh damping in days                 | d             |    0 | real      | kind_phys | none   | F        |
!! | GFS_Control%clam_deep                | entrainment_rate_coefficient_deep_convection                                  | entrainment rate coefficient for deep conv.                                          | none    |    0 | real      | kind_phys | none   | F        |
!! | GFS_Control%c0s_deep                 | rain_conversion_parameter_deep_convection                                     | convective rain conversion parameter for deep conv.                                  | m-1     |    0 | real      | kind_phys | none   | F        |
!! | GFS_Control%c1_deep                  | detrainment_conversion_parameter_deep_convection                              | convective detrainment conversion parameter for deep conv.                           | m-1     |    0 | real      | kind_phys | none   | F        |
!! | GFS_Control%betal_deep               | downdraft_fraction_reaching_surface_over_land_deep_convection                 | downdraft fraction reaching surface over land for deep conv.                         | frac    |    0 | real      | kind_phys | none   | F        |
!! | GFS_Control%betas_deep               | downdraft_fraction_reaching_surface_over_ocean_deep_convection                | downdraft fraction reaching surface over ocean for deep conv.                        | frac    |    0 | real      | kind_phys | none   | F        |
!! | GFS_Control%evfact_deep              | rain_evaporation_coefficient_deep_convection                                  | convective rain evaporation coefficient for deep conv.                               | frac    |    0 | real      | kind_phys | none   | F        |
!! | GFS_Control%evfactl_deep             | rain_evaporation_coefficient_over_land_deep_convection                        | convective rain evaporation coefficient over land for deep conv.                     | frac    |    0 | real      | kind_phys | none   | F        |
!! | GFS_Control%pgcon_deep               | momentum_transport_reduction_factor_pgf_deep_convection                       | reduction factor in momentum transport due to deep conv. induced pressure gradient force  | frac    |    0 | real      | kind_phys | none   | F        |
!! | GFS_Control%asolfac_deep             | aerosol_aware_parameter_deep_convection                                       | aerosol-aware parameter inversely proportional to CCN number concentraion from Lim (2011) for deep conv. | none    |    0 | real      | kind_phys | none   | F        |
!! | GFS_Control%clam_shal                | entrainment_rate_coefficient_shallow_convection                               | entrainment rate coefficient for shal conv.                                                              | none    |    0 | real      | kind_phys | none   | F        |
!! | GFS_Control%c0s_shal                 | rain_conversion_parameter_shallow_convection                                  | convective rain conversion parameter for shal conv.                                                      | m-1     |    0 | real      | kind_phys | none   | F        |
!! | GFS_Control%c1_shal                  | detrainment_conversion_parameter_shallow_convection                           | convective detrainment conversion parameter for shal conv.                                               | m-1     |    0 | real      | kind_phys | none   | F        |
!! | GFS_Control%pgcon_shal               | momentum_transport_reduction_factor_pgf_shallow_convection                    | reduction factor in momentum transport due to shal conv. induced pressure gradient force                 | frac    |    0 | real      | kind_phys | none   | F        |
!! | GFS_Control%asolfac_shal             | aerosol_aware_parameter_shallow_convection                                    | aerosol-aware parameter inversely proportional to CCN number concentraion from Lim (2011) for shal conv. | none    |    0 | real      | kind_phys | none   | F        |
!! | GFS_Control%nst_anl                  |                                                                               | flag for NSSTM analysis in gcycle/sfcsub                             |               |    0 | logical   |           | none   | F        |
!! | GFS_Control%lsea                     |                                                                               |                                                                      |               |    0 | integer   |           | none   | F        |
!! | GFS_Control%xkzm_m                   | atmosphere_momentum_diffusivity_background                                    | background vertical diffusion for momentum                           | m2 s-1        |    0 | real      | kind_phys | none   | F        |
!! | GFS_Control%xkzm_h                   | atmosphere_heat_diffusivity_background                                        | background vertical diffusion for heat q                             | m2 s-1        |    0 | real      | kind_phys | none   | F        |
!! | GFS_Control%xkzm_s                   | diffusivity_background_sigma_level                                            | sigma threshold for background mom. diffusion                        | none          |    0 | real      | kind_phys | none   | F        |
!! | GFS_Control%nstf_name                |                                                                               |                                                                      |               |    1 | integer   |           | none   | F        |
!! | GFS_Control%nstf_name(1)             | flag_for_nsstm_run                                                            | NSSTM flag: off/uncoupled/coupled=0/1/2                              | flag          |    0 | integer   |           | none   | F        |
!! | GFS_Control%nstf_name(4)             | vertical_temperature_average_range_lower_bound                                | zsea1 in mm                                                          | mm            |    0 | integer   |           | none   | F        |
!! | GFS_Control%nstf_name(5)             | vertical_temperature_average_range_upper_bound                                | zsea2 in mm                                                          | mm            |    0 | integer   |           | none   | F        |
!! | GFS_Control%xkzminv                  | atmosphere_heat_diffusivity_background_maximum                                | maximum background value of heat diffusivity                         | m2 s-1        |    0 | real      | kind_phys | none   | F        |
!! | GFS_Control%moninq_fac               | atmosphere_diffusivity_coefficient_factor                                     | multiplicative constant for atmospheric diffusivities                | none          |    0 | real      | kind_phys | none   | F        |
!! | GFS_Control%nca                      | number_of_independent_cellular_automata                                       | number of independent cellular automata                              | count         |    0 | integer   |           | none   | F        |
!! | GFS_Control%nlives                   | cellular_automata_lifetime                                                    | cellular automata lifetime                                           | count         |    0 | integer   |           | none   | F        |
!! | GFS_Control%ncells                   | cellular_automata_finer_grid                                                  | cellular automata finer grid                                         | count         |    0 | integer   |           | none   | F        |
!! | GFS_Control%nfracseed                | cellular_automata_seed_probability                                            | cellular automata seed probability                                   | fraction      |    0 | real      | kind_phys | none   | F        |
!! | GFS_Control%nseed                    | cellular_automata_seed_frequency                                              | cellular automata seed frequency in units of time steps              | count         |    0 | integer   |           | none   | F        |
!! | GFS_Control%do_ca                    | flag_for_cellular_automata                                                    | cellular automata main switch                                        | flag          |    0 | logical   |           | none   | F        |
!! | GFS_Control%ca_sgs                   | flag_for_sgs_cellular_automata                                                | switch for sgs ca                                                    | flag          |    0 | logical   |           | none   | F        |
!! | GFS_Control%ca_global                | flag_for_global_cellular_automata                                             | switch for global ca                                                 | flag          |    0 | logical   |           | none   | F        |
!! | GFS_Control%ca_smooth                | flag_for_gaussian_spatial_filter                                              | switch for gaussian spatial filter                                   | flag          |    0 | logical   |           | none   | F        |
!! | GFS_Control%isppt_deep               | flag_for_combination_of_sppt_with_isppt_deep                                  | switch for combination with isppt_deep.                              | flag          |    0 | logical   |           | none   | F        |
!! | GFS_Control%iseed_ca                 | seed_for_random_number_generation_in_cellular_automata_scheme                 | seed for random number generation in ca scheme                       | none          |    0 | integer   |           | none   | F        |
!! | GFS_Control%nspinup                  | number_of_iterations_to_spin_up_cellular_automata                             | number of iterations to spin up the ca                               | count         |    0 | integer   |           | none   | F        |
!! | GFS_Control%nthresh                  | threshold_for_perturbed_vertical_velocity                                     | threshold used for perturbed vertical velocity                       | m s-1         |    0 | real      | kind_phys | none   | F        |
!! | GFS_Control%do_sppt                  | flag_for_stochastic_surface_physics_perturbations                             | flag for stochastic surface physics perturbations                    | flag          |    0 | logical   |           | none   | F        |
!! | GFS_Control%use_zmtnblck             | flag_for_mountain_blocking                                                    | flag for mountain blocking                                           | flag          |    0 | logical   |           | none   | F        |
!! | GFS_Control%do_shum                  | flag_for_stochastic_shum_option                                               | flag for stochastic shum option                                      | flag          |    0 | logical   |           | none   | F        |
!! | GFS_Control%do_skeb                  | flag_for_stochastic_skeb_option                                               | flag for stochastic skeb option                                      | flag          |    0 | logical   |           | none   | F        |
!! | GFS_Control%skeb_npass               |                                                                               |                                                                      |               |    0 | integer   |           | none   | F        |
!! | GFS_Control%do_sfcperts              | flag_for_stochastic_surface_perturbations                                     | flag for stochastic surface perturbations option                     | flag          |    0 | logical   |           | none   | F        |
!! | GFS_Control%nsfcpert                 | number_of_surface_perturbations                                               | number of surface perturbations                                      | count         |    0 | integer   |           | none   | F        |
!! | GFS_Control%pertz0                   | magnitude_of_perturbation_of_momentum_roughness_length                        | magnitude of perturbation of momentum roughness length               | frac          |    1 | real      | kind_phys | none   | F        |
!! | GFS_Control%pertzt                   | magnitude_of_perturbation_of_heat_to_momentum_roughness_length_ratio          | magnitude of perturbation of heat to momentum roughness length ratio | frac          |    1 | real      | kind_phys | none   | F        |
!! | GFS_Control%pertshc                  | magnitude_of_perturbation_of_soil_type_b_parameter                            | magnitude of perturbation of soil type b parameter                   | frac          |    1 | real      | kind_phys | none   | F        |
!! | GFS_Control%pertlai                  | magnitude_of_perturbation_of_leaf_area_index                                  | magnitude of perturbation of leaf area index                         | frac          |    1 | real      | kind_phys | none   | F        |
!! | GFS_Control%pertalb                  | magnitude_of_surface_albedo_perturbation                                      | magnitude of surface albedo perturbation                             | frac          |    1 | real      | kind_phys | none   | F        |
!! | GFS_Control%pertvegf                 | magnitude_of_perturbation_of_vegetation_fraction                              | magnitude of perturbation of vegetation fraction                     | frac          |    1 | real      | kind_phys | none   | F        |
!! | GFS_Control%tracer_names             |                                                                               | array of initialized tracers from dynamic core                       |               |    1 | character |           | none   | F        |
!! | GFS_Control%ntrac                    | number_of_tracers                                                             | number of tracers                                                    | count         |    0 | integer   |           | none   | F        |
!! | GFS_Control%ntracp1                  | number_of_tracers_plus_one                                                    | number of tracers plus one                                           | count         |    0 | integer   |           | none   | F        |
!! | GFS_Control%ntqv                     | index_for_water_vapor                                                         | tracer index for water vapor (specific humidity)                     | index         |    0 | integer   |           | none   | F        |
!! | GFS_Control%ntoz                     | index_for_ozone                                                               | tracer index for ozone mixing ratio                                  | index         |    0 | integer   |           | none   | F        |
!! | GFS_Control%ntcw                     | index_for_liquid_cloud_condensate                                             | tracer index for cloud condensate (or liquid water)                  | index         |    0 | integer   |           | none   | F        |
!! | GFS_Control%ntiw                     | index_for_ice_cloud_condensate                                                | tracer index for  ice water                                          | index         |    0 | integer   |           | none   | F        |
!! | GFS_Control%ntrw                     | index_for_rain_water                                                          | tracer index for rain water                                          | index         |    0 | integer   |           | none   | F        |
!! | GFS_Control%ntsw                     | index_for_snow_water                                                          | tracer index for snow water                                          | index         |    0 | integer   |           | none   | F        |
!! | GFS_Control%ntgl                     | index_for_graupel                                                             | tracer index for graupel                                             | index         |    0 | integer   |           | none   | F        |
!! | GFS_Control%ntclamt                  | index_for_cloud_amount                                                        | tracer index for cloud amount integer                                | index         |    0 | integer   |           | none   | F        |
!! | GFS_Control%ntlnc                    | index_for_liquid_cloud_number_concentration                                   | tracer index for liquid number concentration                         | index         |    0 | integer   |           | none   | F        |
!! | GFS_Control%ntinc                    | index_for_ice_cloud_number_concentration                                      | tracer index for ice    number concentration                         | index         |    0 | integer   |           | none   | F        |
!! | GFS_Control%ntrnc                    | index_for_rain_number_concentration                                           | tracer index for rain   number concentration                         | index         |    0 | integer   |           | none   | F        |
!! | GFS_Control%ntsnc                    | index_for_snow_number_concentration                                           | tracer index for snow   number concentration                         | index         |    0 | integer   |           | none   | F        |
!! | GFS_Control%ntgnc                    | index_for_graupel_number_concentration                                        | tracer index for graupel number concentration                        | index         |    0 | integer   |           | none   | F        |
!! | GFS_Control%ntke                     | index_for_turbulent_kinetic_energy                                            | tracer index for turbulent kinetic energy                            | index         |    0 | integer   |           | none   | F        |
!! | GFS_Control%nto                      |                                                                               | tracer index for oxygen ion                                          |               |    0 | integer   |           | none   | F        |
!! | GFS_Control%nto2                     |                                                                               | tracer index for oxygen                                              |               |    0 | integer   |           | none   | F        |
!! | GFS_Control%ntwa                     | index_for_water_friendly_aerosols                                             | tracer index for water friendly aerosol                              | index         |    0 | integer   |           | none   | F        |
!! | GFS_Control%ntia                     | index_for_ice_friendly_aerosols                                               | tracer index for ice friendly aerosol                                | index         |    0 | integer   |           | none   | F        |
!! | GFS_Control%ntchm                    | number_of_chemical_tracers                                                    | number of chemical tracers                                           | count         |    0 | integer   |           | none   | F        |
!! | GFS_Control%ntchs                    | index_for_first_chemical_tracer                                               | tracer index for first chemical tracer                               | index         |    0 | integer   |           | none   | F        |
!! | GFS_Control%ntdiag                   | diagnostics_control_for_chemical_tracers                                      | array to control diagnostics for chemical tracers                    | flag          |    1 | logical   |           | none   | F        |
!! | GFS_Control%ntot2d                   |                                                                               | total number of variables for phyf2d                                 |               |    0 | integer   |           | none   | F        |
!! | GFS_Control%ntot3d                   |                                                                               | total number of variables for phyf3d                                 |               |    0 | integer   |           | none   | F        |
!! | GFS_Control%indcld                   | index_for_cloud_fraction_in_3d_arrays_for_microphysics                        | index of cloud fraction in phyf3d (used only for SHOC or MG)         | index         |    0 | integer   |           | none   | F        |
!! | GFS_Control%num_p2d                  | array_dimension_of_2d_arrays_for_microphysics                                 | number of 2D arrays needed for microphysics                          | count         |    0 | integer   |           | none   | F        |
!! | GFS_Control%num_p3d                  | array_dimension_of_3d_arrays_for_microphysics                                 | number of 3D arrays needed for microphysics                          | count         |    0 | integer   |           | none   | F        |
!! | GFS_Control%nshoc_2d                 |                                                                               | number of 2d fields for SHOC                                         |               |    0 | integer   |           | none   | F        |
!! | GFS_Control%nshoc_3d                 |                                                                               | number of 3d fields for SHOC                                         |               |    0 | integer   |           | none   | F        |
!! | GFS_Control%ncnvcld3d                | number_of_convective_3d_cloud_fields                                          | number of convective 3d clouds fields                                | count         |    0 | integer   |           | none   | F        |
!! | GFS_Control%npdf3d                   | number_of_3d_arrays_associated_with_pdf-based_clouds                          | number of 3d arrays associated with pdf based clouds/mp              | count         |    0 | integer   |           | none   | F        |
!! | GFS_Control%nctp                     | number_of_cloud_types_CS                                                      | number of cloud types in Chikira-Sugiyama scheme                     | count         |    0 | integer   |           | none   | F        |
!! | GFS_Control%ncnvw                    |                                                                               | the index of cnvw in phy_f3d                                         |               |    0 | integer   |           | none   | F        |
!! | GFS_Control%ncnvc                    |                                                                               | the index of cnvc in phy_f3d                                         |               |    0 | integer   |           | none   | F        |
!! | GFS_Control%nleffr                   | index_for_cloud_liquid_water_effective_radius                                 | the index of cloud liquid water effective radius in phy_f3d          |               |    0 | integer   |           | none   | F        |
!! | GFS_Control%nieffr                   | index_for_ice_effective_radius                                                | the index of ice effective radius in phy_f3d                         |               |    0 | integer   |           | none   | F        |
!! | GFS_Control%nreffr                   | index_for_rain_effective_radius                                               | the index of rain effective radius in phy_f3d                        |               |    0 | integer   |           | none   | F        |
!! | GFS_Control%nseffr                   | index_for_snow_effective_radius                                               | the index of snow effective radius in phy_f3d                        |               |    0 | integer   |           | none   | F        |
!! | GFS_Control%ngeffr                   | index_for_graupel_effective_radius                                            | the index of graupel effective radius in phy_f3d                     |               |    0 | integer   |           | none   | F        |
!! | GFS_Control%debug                    |                                                                               | debug flag                                                           |               |    0 | logical   |           | none   | F        |
!! | GFS_Control%pre_rad                  |                                                                               | flag for testing purpose                                             |               |    0 | logical   |           | none   | F        |
!! | GFS_Control%ipt                      |                                                                               | index for diagnostic printout point                                  |               |    0 | integer   |           | none   | F        |
!! | GFS_Control%lprnt                    | flag_print                                                                    | control flag for diagnostic print out                                | flag          |    0 | logical   |           | none   | F        |
!! | GFS_Control%lsswr                    | flag_to_calc_sw                                                               | logical flags for sw radiation calls                                 | flag          |    0 | logical   |           | none   | F        |
!! | GFS_Control%lslwr                    | flag_to_calc_lw                                                               | logical flags for lw radiation calls                                 | flag          |    0 | logical   |           | none   | F        |
!! | GFS_Control%solhr                    | forecast_hour                                                                 | hour time after 00z at the t-step                                    | h             |    0 | real      | kind_phys | none   | F        |
!! | GFS_Control%solcon                   | solar_constant                                                                | solar constant (sun-earth distant adjusted)                          | W m-2         |    0 | real      | kind_phys | none   | F        |
!! | GFS_Control%slag                     | equation_of_time                                                              | equation of time (radian)                                            | radians       |    0 | real      | kind_phys | none   | F        |
!! | GFS_Control%sdec                     | sine_of_solar_declination_angle                                               | sin of the solar declination angle                                   | none          |    0 | real      | kind_phys | none   | F        |
!! | GFS_Control%cdec                     | cosine_of_solar_declination_angle                                             | cos of the solar declination angle                                   | none          |    0 | real      | kind_phys | none   | F        |
!! | GFS_Control%clstp                    | convective_cloud_switch                                                       | index used by cnvc90 (for convective clouds)                         | none          |    0 | real      | kind_phys | none   | F        |
!! | GFS_Control%phour                    |                                                                               | previous forecast time                                               | h             |    0 | real      | kind_phys | none   | F        |
!! | GFS_Control%fhour                    | forecast_time                                                                 | curent forecast time                                                 | h             |    0 | real      | kind_phys | none   | F        |
!! | GFS_Control%zhour                    |                                                                               | previous hour diagnostic buckets emptied                             | h             |    0 | real      | kind_phys | none   | F        |
!! | GFS_Control%kdt                      | index_of_time_step                                                            | current forecast iteration                                           | index         |    0 | integer   |           | none   | F        |
!! | GFS_Control%first_time_step          | flag_for_first_time_step                                                      | flag for first time step for time integration loop (cold/warmstart)  | flag          |    0 | logical   |           | none   | F        |
!! | GFS_Control%restart                  | flag_for_restart                                                              | flag for restart (warmstart) or coldstart                            | flag          |    0 | logical   |           | none   | F        |
!! | GFS_Control%hydrostatic              | flag_for_hydrostatic_solver                                                   | flag for hydrostatic solver from dynamics                            | flag          |    0 | logical   |           | none   | F        |
!! | GFS_Control%jdat                     | forecast_date_and_time                                                        | current forecast date and time                                       | none          |    1 | integer   |           | none   | F        |
!! | GFS_Control%iccn                     | flag_for_in_ccn_forcing_for_morrison_gettelman_microphysics                   | flag for IN and CCN forcing for morrison gettelman microphysics      | flag          |    0 | logical   |           | none   | F        |
!! | GFS_Control%sec                      | seconds_elapsed_since_model_initialization                                    | seconds elapsed since model initialization                           | s             |    0 | real      | kind_phys | none   | F        |
!! | GFS_Control%si                       | vertical_sigma_coordinate_for_radiation_initialization                        | vertical sigma coordinate for radiation initialization               | none          |    1 | real      | kind_phys | none   | F        |
!! | GFS_Control%iau_delthrs              |                                                                               | iau time interval (to scale increments) in hours                     |               |    0 | real      | kind_phys | none   | F        |
!! | GFS_Control%iau_inc_files            |                                                                               | list of increment files                                              |               |    1 | character | len=240   | none   | F        |
!! | GFS_Control%iaufhrs                  |                                                                               | forecast hours associated with increment files                       |               |    1 | real      | kind_phys | none   | F        |
!! | GFS_Control%iau_filter_increments    |                                                                               |                                                                      |               |    0 | logical   |           | none   | F        |
!! | GFS_Control%dxinv                    | inverse_scaling_factor_for_critical_relative_humidity                         | inverse scaling factor for critical relative humidity                | rad2 m-2      |    0 | real      | kind_phys | none   | F        |
!! | GFS_Control%dxmax                    | maximum_scaling_factor_for_critical_relative_humidity                         | maximum scaling factor for critical relative humidity                | m2 rad-2      |    0 | real      | kind_phys | none   | F        |
!! | GFS_Control%dxmin                    | minimum_scaling_factor_for_critical_relative_humidity                         | minimum scaling factor for critical relative humidity                | m2 rad-2      |    0 | real      | kind_phys | none   | F        |
!! | GFS_Control%rhcmax                   | maximum_critical_relative_humidity                                            | maximum critical relative humidity                                   | frac          |    0 | real      | kind_phys | none   | F        |
!! | GFS_Control%do_mynnedmf              | do_mynnedmf                                                                   | flag to activate MYNN-EDMF                                           | flag          |    0 | logical   |           | none   | F        |
!! | GFS_Control%do_mynnsfclay            | do_mynnsfclay                                                                 | flag to activate MYNN surface layer                                  | flag          |    0 | logical   |           | none   | F        |
!! | GFS_Control%grav_settling            | grav_settling                                                                 | flag to activate gravitational setting of fog                        | flag          |    0 | integer   |           | none   | F        |
!! | GFS_Control%bl_mynn_tkebudget        | tke_budget                                                                    | flag for activating TKE budget                                       | flag          |    0 | integer   |           | none   | F        |
!! | GFS_Control%bl_mynn_tkeadvect        | tke_advect                                                                    | flag for activating TKE advection                                    | flag          |    0 | logical   |           | none   | F        |
!! | GFS_Control%bl_mynn_cloudpdf         | cloudpdf                                                                      | flag to determine which cloud PDF to use                             | flag          |    0 | integer   |           | none   | F        |
!! | GFS_Control%bl_mynn_mixlength        | mixing_length_flag                                                            | flag to determine which mixing length form to use                    | flag          |    0 | integer   |           | none   | F        |
!! | GFS_Control%bl_mynn_edmf             | edmf_flag                                                                     | flag to activate the mass-flux scheme                                | flag          |    0 | integer   |           | none   | F        |
!! | GFS_Control%bl_mynn_edmf_mom         | edmf_momentum_transport_flag                                                  | flag to activate the transport of momentum                           | flag          |    0 | integer   |           | none   | F        |
!! | GFS_Control%bl_mynn_edmf_tke         | edmf_tke_transport_flag                                                       | flag to activate the transport of TKE                                | flag          |    0 | integer   |           | none   | F        |
!! | GFS_Control%bl_mynn_edmf_part        | edmf_partition_flag                                                           | flag to partitioning og the MF and ED areas                          | flag          |    0 | integer   |           | none   | F        |
!! | GFS_Control%bl_mynn_cloudmix         | cloud_specie_mix_flag                                                         | flag to activate mixing of cloud species                             | flag          |    0 | integer   |           | none   | F        |
!! | GFS_Control%bl_mynn_mixqt            | mix_total_water_flag                                                          | flag to mix total water or individual species                        | flag          |    0 | integer   |           | none   | F        |
!! | GFS_Control%icloud_bl                | couple_sgs_clouds_to_radiation_flag                                           | flag for coupling sgs clouds to radiation                            | flag          |    0 | integer   |           | none   | F        |
!!
#endif
  type GFS_control_type

    integer              :: me              !< MPI rank designator
    integer              :: master          !< MPI rank of master atmosphere processor
#ifdef CCPP
    integer              :: communicator    !< MPI communicator
    integer              :: ntasks          !< MPI size in communicator
    integer              :: nthreads        !< OpenMP threads available for physics
#endif
    integer              :: nlunit          !< unit for namelist
    character(len=64)    :: fn_nml          !< namelist filename for surface data cycling
    character(len=256), pointer :: input_nml_file(:) !< character string containing full namelist
                                                   !< for use with internal file reads
#ifdef CCPP
    integer              :: logunit
#endif
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
#ifdef CCPP
    !--- ak/bk for pressure level calculations
    real(kind=kind_phys), pointer :: ak(:)  !< from surface (k=1) to TOA (k=levs)
    real(kind=kind_phys), pointer :: bk(:)  !< from surface (k=1) to TOA (k=levs)
#endif
    integer              :: cnx             !< number of points in the i-dir for this cubed-sphere face
    integer              :: cny             !< number of points in the j-dir for this cubed-sphere face
    integer              :: lonr            !< number of global points in x-dir (i) along the equator
    integer              :: latr            !< number of global points in y-dir (j) along any meridian
    integer              :: tile_num
#ifdef CCPP
    integer,     pointer :: blksz(:)        !< for explicit data blocking: block sizes of all blocks
#endif

!--- coupling parameters
    logical              :: cplflx          !< default no cplflx collection
    logical              :: cplwav          !< default no cplwav collection
    logical              :: cplchm          !< default no cplchm collection

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
    integer              :: ico2            !< \f$CO_2\f$ data source control flag
                                            !!\brief 0: prescribed value (380 ppmv)
                                            !!\n 1: yearly global averaged annual mean from observations
                                            !!\n 2: monthly 15 degree horizontal resolution from observations
    integer              :: ialb            !< SW surface albedo control flag
                                            !!\brief 0: using climatology surface albedo scheme for SW
                                            !!\n 1: using MODIS based land surface albedo for SW
    integer              :: iems            !< LW global surface emissivity control flag
                                            !!\brief 0: black-body emissivity
                                            !!\n 1: surface type based climatology in \f$1^0\f$ horizontal resolution
    integer              :: iaer            !< aerosol effect control flag; 3-digit flag "abc" (volcanic, LW, SW)
                                            !!\brief a: stratospheric volcanic aerosols
                                            !!\n b: tropospheric aerosols for LW
                                            !!\n c: tropospheric aerosols for SW
                                            !!\n 0: aerosol effect is not included; 1: aerosol effect is included
    integer              :: icliq_sw        !< sw optical property for liquid clouds
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
    integer              :: isubc_sw        !< subgrid cloud approximation control flag in SW radiation
                                            !!\brief 0: no McICA approximation in SW radiation
                                            !!\n 1: use McICA with prescribed permutation seeds (test mode)
                                            !!\n 2: use McICA with randomly generated permutation seeds 
    integer              :: isubc_lw        !< subgrid cloud approximation in LW radiation
                                            !!\brief 0: no McICA approximation in LW radiation
                                            !!\n 1: use McICA with prescribed permutation seeds (test mode)
                                            !!\n 2: use McICA with randomly generated permutation seeds
    logical              :: crick_proof     !< CRICK-Proof cloud water
    logical              :: ccnorm          !< Cloud condensate normalized by cloud cover
    logical              :: norad_precip    !< radiation precip flag for Ferrier/Moorthi
    logical              :: lwhtr           !< flag to output lw heating rate (Radtend%lwhc)
    logical              :: swhtr           !< flag to output sw heating rate (Radtend%swhc)

!--- microphysical switch
    integer              :: ncld            !< number of hydrometeors
    !--- new microphysical switch
    integer              :: imp_physics                    !< choice of microphysics scheme
                                                           !!\brief 11: choice of GFDL microphysics scheme
                                                           !!\n 8: choice of Thompson microphysics scheme
                                                           !!\n 6: choice of WSMG microphysics scheme
                                                           !!\n 99: choice of Zhao-Carr microphysics scheme
                                                           !!\n 98: choice of Zhao-Carr microphysics scheme with PDF clouds
                                                           !!\n 10: choice of Morrison-Gettelman microphysics scheme
    integer              :: imp_physics_gfdl = 11          !< choice of GFDL     microphysics scheme
    integer              :: imp_physics_thompson = 8       !< choice of Thompson microphysics scheme
    integer              :: imp_physics_wsm6 = 6           !< choice of WSMG     microphysics scheme
    integer              :: imp_physics_zhao_carr = 99     !< choice of Zhao-Carr microphysics scheme
    integer              :: imp_physics_zhao_carr_pdf = 98 !< choice of Zhao-Carr microphysics scheme with PDF clouds
    integer              :: imp_physics_mg = 10            !< choice of Morrison-Gettelman microphysics scheme
    !--- Z-C microphysical parameters
    real(kind=kind_phys) :: psautco(2)         !< [in] auto conversion coeff from ice to snow
    real(kind=kind_phys) :: prautco(2)         !< [in] auto conversion coeff from cloud to rain
    real(kind=kind_phys) :: evpco              !< [in] coeff for evaporation of largescale rain
    real(kind=kind_phys) :: wminco(2)          !< [in] water and ice minimum threshold for Zhao
    real(kind=kind_phys) :: avg_max_length     !< reset time in seconds for max hourly fields
    !--- M-G microphysical parameters
    integer              :: fprcp              !< no prognostic rain and snow (MG)
    integer              :: pdfflag            !< pdf flag for MG macrophysics
    real(kind=kind_phys) :: mg_dcs             !< autoconversion size threshold for cloud ice to snow for MG micrphysics
    real(kind=kind_phys) :: mg_qcvar           !< cloud water relative variance for MG microphysics
    real(kind=kind_phys) :: mg_ts_auto_ice(2)  !< autoconversion time scale for ice in MG microphysics
#ifdef CCPP
    real(kind=kind_phys) :: mg_rhmini          !< relative humidity threshold parameter for nucleating ice
#endif

    real(kind=kind_phys) :: mg_ncnst           !< constant droplet num concentration (m-3)
    real(kind=kind_phys) :: mg_ninst           !< constant ice num concentration (m-3)
    real(kind=kind_phys) :: mg_ngnst           !< constant graupel/hail num concentration (m-3)
    real(kind=kind_phys) :: mg_berg_eff_factor !< berg efficiency factor
    real(kind=kind_phys) :: mg_alf             !< tuning factor for alphas (alpha = 1-critial relative humidity) in MG microphysics
    real(kind=kind_phys) :: mg_qcmin(2)        !< min liquid and ice mixing ratio in Mg macro clouds
    character(len=16)    :: mg_precip_frac_method ! type of precipitation fraction method
#ifdef CCPP
    real(kind=kind_phys) :: tf
    real(kind=kind_phys) :: tcr
    real(kind=kind_phys) :: tcrf
#endif
!
    logical              :: effr_in            !< eg to turn on ffective radii for MG
    logical              :: microp_uniform     !< flag for uniform subcolumns for MG microphysics
    logical              :: do_cldliq
    logical              :: do_cldice
    logical              :: hetfrz_classnuc

    logical              :: mg_nccons
    logical              :: mg_nicons
    logical              :: mg_ngcons
    logical              :: sed_supersat
    logical              :: do_sb_physics
    logical              :: mg_do_graupel
    logical              :: mg_do_hail
    logical              :: mg_do_ice_gmao
    logical              :: mg_do_liq_liu

    real(kind=kind_phys) :: shoc_parm(5)    !< critical pressure in Pa for tke dissipation in shoc
    integer              :: ncnd            !< number of cloud condensate types

    !--- Thompson's microphysical paramters
    logical              :: ltaerosol       !< flag for aerosol version
    logical              :: lradar          !< flag for radar reflectivity
    real(kind=kind_phys) :: ttendlim        !< temperature tendency limiter per time step in K/s

    !--- GFDL microphysical paramters
    logical              :: lgfdlmprad      !< flag for GFDL mp scheme and radiation consistency

    !--- land/surface model parameters
    integer              :: lsm             !< flag for land surface model lsm=1 for noah lsm
    integer              :: lsm_noah=1      !< flag for NOAH land surface model
    integer              :: lsm_ruc=2       !< flag for RUC land surface model
    integer              :: lsoil           !< number of soil layers
#ifdef CCPP
    integer              :: lsoil_lsm       !< number of soil layers internal to land surface model
#endif
    integer              :: ivegsrc         !< ivegsrc = 0   => USGS,
                                            !< ivegsrc = 1   => IGBP (20 category)
                                            !< ivegsrc = 2   => UMD  (13 category)
    integer              :: isot            !< isot = 0   => Zobler soil type  ( 9 category)
                                            !< isot = 1   => STATSGO soil type (19 category)
    logical              :: mom4ice         !< flag controls mom4 sea ice
    logical              :: use_ufo         !< flag for using unfiltered orography surface option

!--- tuning parameters for physical parameterizations
    logical              :: ras             !< flag for ras convection scheme
    logical              :: flipv           !< flag for vertical direction flip (ras)
                                            !< .true. implies surface at k=1
    logical              :: trans_trac      !< flag for convective transport of tracers (RAS, CS, or SAMF)
    logical              :: old_monin       !< flag for diff monin schemes
    logical              :: cnvgwd          !< flag for conv gravity wave drag
    logical              :: mstrat          !< flag for moorthi approach for stratus
    logical              :: moist_adj       !< flag for moist convective adjustment
    logical              :: cscnv           !< flag for Chikira-Sugiyama convection
    logical              :: cal_pre         !< flag controls precipitation type algorithm
    logical              :: do_aw           !< AW scale-aware option in cs convection
    logical              :: do_awdd         !< AW scale-aware option in cs convection
    logical              :: flx_form        !< AW scale-aware option in cs convection
    logical              :: do_shoc         !< flag for SHOC
    logical              :: shocaftcnv      !< flag for SHOC
    logical              :: shoc_cld        !< flag for clouds
    logical              :: uni_cld         !< flag for clouds in grrad
#ifdef CCPP
    logical              :: oz_phys         !< flag for old (2006) ozone physics
    logical              :: oz_phys_2015    !< flag for new (2015) ozone physics
#endif
    logical              :: h2o_phys        !< flag for stratosphere h2o
    logical              :: pdfcld          !< flag for pdfcld
    logical              :: shcnvcw         !< flag for shallow convective cloud
    logical              :: redrag          !< flag for reduced drag coefficient for high wind over sea
    logical              :: hybedmf         !< flag for hybrid edmf pbl scheme
    logical              :: satmedmf        !< flag for scale-aware TKE-based moist edmf
                                            !< vertical turbulent mixing scheme
    logical              :: shinhong        !< flag for scale-aware Shinhong vertical turbulent mixing scheme
    logical              :: do_ysu          !< flag for YSU turbulent mixing scheme
    logical              :: dspheat         !< flag for TKE dissipative heating
    logical              :: lheatstrg       !< flag for canopy heat storage parameterization
    logical              :: cnvcld          !< flag for convective clouds
    logical              :: random_clds     !< flag controls whether clouds are random
    logical              :: shal_cnv        !< flag for calling shallow convection
    logical              :: do_deep         !< whether to do deep convection
    integer              :: imfshalcnv      !< flag for mass-flux shallow convection scheme
                                            !!\n  1: July 2010 version of mass-flux shallow conv scheme
                                            !! (current operational version as of 2016)
                                            !!\n  2: scale- & aerosol-aware mass-flux shallow conv scheme (2017)
                                            !!\n  0: modified Tiedtke's eddy-diffusion shallow conv scheme
                                            !!\n -1: no shallow convection used
    integer              :: imfdeepcnv      !< flag for mass-flux deep convection scheme
                                            !!\n 1: July 2010 version of SAS conv scheme
                                            !!    (current operational version as of 2016)
                                            !!\n 2: scale- & aerosol-aware mass-flux deep conv scheme (2017)
                                            !!\n 0: old SAS Convection scheme before July 2010
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
    real(kind=kind_phys) :: crtrh(3)        !< critical relative humidity at the surface,
                                            !! PBL top and at the top of the atmosphere
    real(kind=kind_phys) :: dlqf(2)         !< factor for cloud condensate detrainment
                                            !< from cloud edges for RAS
    real(kind=kind_phys) :: psauras(2)      !< [in] auto conversion coeff from ice to snow in ras
    real(kind=kind_phys) :: prauras(2)      !< [in] auto conversion coeff from cloud to rain in ras
    real(kind=kind_phys) :: wminras(2)      !< [in] water and ice minimum threshold for ras

    integer              :: seed0           !< random seed for radiation

    real(kind=kind_phys) :: rbcr            !< Critical Richardson Number in the PBL scheme
#ifdef CCPP
    !--- MYNN parameters/switches
    logical              :: do_mynnedmf
    logical              :: do_mynnsfclay
    integer              :: grav_settling      !< flag for initalizing fist time step
    integer              :: bl_mynn_tkebudget  !< flag for activating TKE budget
    logical              :: bl_mynn_tkeadvect  !< activate computation of TKE advection (not yet in use for FV3)
    integer              :: bl_mynn_cloudpdf   !< flag to determine which cloud PDF to use
    integer              :: bl_mynn_mixlength  !< flag for different version of mixing length formulation
    integer              :: bl_mynn_edmf       !< flag to activate the mass-flux scheme
    integer              :: bl_mynn_edmf_mom   !< flag to activate the transport of momentum
    integer              :: bl_mynn_edmf_tke   !< flag to activate the transport of TKE
    integer              :: bl_mynn_edmf_part  !< flag to partitioning og the MF and ED areas
    integer              :: bl_mynn_cloudmix   !< flag to activate mixing of cloud species
    integer              :: bl_mynn_mixqt      !< flag to mix total water or individual species
    integer              :: icloud_bl          !< flag for coupling sgs clouds to radiation
#endif

!--- Rayleigh friction
    real(kind=kind_phys) :: prslrd0         !< pressure level from which Rayleigh Damping is applied
    real(kind=kind_phys) :: ral_ts          !< time scale for Rayleigh damping in days

!--- mass flux deep convection
    real(kind=kind_phys) :: clam_deep       !< c_e for deep convection (Han and Pan, 2011, eq(6))
    real(kind=kind_phys) :: c0s_deep        !< convective rain conversion parameter
    real(kind=kind_phys) :: c1_deep         !< conversion parameter of detrainment from liquid water into grid-scale cloud water
    real(kind=kind_phys) :: betal_deep      !< fraction factor of downdraft air mass reaching ground surface over land
    real(kind=kind_phys) :: betas_deep      !< fraction factor of downdraft air mass reaching ground surface over sea
    real(kind=kind_phys) :: evfact_deep     !< evaporation factor from convective rain
    real(kind=kind_phys) :: evfactl_deep    !< evaporation factor from convective rain over land
    real(kind=kind_phys) :: pgcon_deep      !< reduction factor in momentum transport due to convection induced pressure gradient force
                                            !< 0.7 : Gregory et al. (1997, QJRMS)
                                            !< 0.55: Zhang & Wu (2003, JAS)
    real(kind=kind_phys) :: asolfac_deep    !< aerosol-aware parameter based on Lim (2011)
                                            !< asolfac= cx / c0s(=.002)
                                            !< cx = min([-0.7 ln(Nccn) + 24]*1.e-4, c0s)
                                            !< Nccn: CCN number concentration in cm^(-3)
                                            !< Until a realistic Nccn is provided, Nccns are assumed
                                            !< as Nccn=100 for sea and Nccn=1000 for land

!--- mass flux shallow convection
    real(kind=kind_phys) :: clam_shal       !< c_e for shallow convection (Han and Pan, 2011, eq(6))
    real(kind=kind_phys) :: c0s_shal        !< convective rain conversion parameter
    real(kind=kind_phys) :: c1_shal         !< conversion parameter of detrainment from liquid water into grid-scale cloud water
    real(kind=kind_phys) :: pgcon_shal      !< reduction factor in momentum transport due to convection induced pressure gradient force
                                            !< 0.7 : Gregory et al. (1997, QJRMS)
                                            !< 0.55: Zhang & Wu (2003, JAS)
    real(kind=kind_phys) :: asolfac_shal    !< aerosol-aware parameter based on Lim (2011)
                                            !< asolfac= cx / c0s(=.002)
                                            !< cx = min([-0.7 ln(Nccn) + 24]*1.e-4, c0s)
                                            !< Nccn: CCN number concentration in cm^(-3)
                                            !< Until a realistic Nccn is provided, Nccns are assumed
                                            !< as Nccn=100 for sea and Nccn=1000 for land

!--- near surface temperature model
    logical              :: nst_anl         !< flag for NSSTM analysis in gcycle/sfcsub
    integer              :: lsea           
    integer              :: nstf_name(5)    !< flag 0 for no nst; 1 for uncoupled nst; and 2 for coupled NST
                                            !!\n nstf_name contains the NSST related parameters:
                                            !!\n nstf_name(1) : 0 = NSSTM off, 1 = NSSTM on but uncoupled,
                                            !!                2 = NSSTM on and coupled
                                            !!\n nstf_name(2) : 1 = NSSTM spin up on, 0 = NSSTM spin up off
                                            !!\n nstf_name(3) : 1 = NSST analysis on, 0 = NSSTM analysis off
                                            !!\n nstf_name(4) : zsea1 in mm
                                            !!\n nstf_name(5) : zsea2 in mm
!--- background vertical diffusion
    real(kind=kind_phys) :: xkzm_m          !< [in] bkgd_vdif_m  background vertical diffusion for momentum  
    real(kind=kind_phys) :: xkzm_h          !< [in] bkgd_vdif_h  background vertical diffusion for heat q  
    real(kind=kind_phys) :: xkzm_s          !< [in] bkgd_vdif_s  sigma threshold for background mom. diffusion  
    real(kind=kind_phys) :: xkzminv         !< diffusivity in inversion layers
    real(kind=kind_phys) :: moninq_fac      !< turbulence diffusion coefficient factor

 !---cellular automata control parameters
    integer              :: nca             !< number of independent cellular automata 
    integer              :: nlives          !< cellular automata lifetime
    integer              :: ncells          !< cellular automata finer grid
    real(kind=kind_phys) :: nfracseed       !< cellular automata seed probability 
    integer              :: nseed           !< cellular automata seed frequency
    logical              :: do_ca           !< cellular automata main switch
    logical              :: ca_sgs          !< switch for sgs ca
    logical              :: ca_global       !< switch for global ca
    logical              :: ca_smooth       !< switch for gaussian spatial filter
    logical              :: isppt_deep      !< switch for combination with isppt_deep. OBS! Switches off SPPT on other tendencies!
    integer              :: iseed_ca        !< seed for random number generation in ca scheme
    integer              :: nspinup         !< number of iterations to spin up the ca
    real(kind=kind_phys) :: nthresh         !< threshold used for perturbed vertical velocity
    
!--- stochastic physics control parameters
    logical              :: do_sppt
    logical              :: use_zmtnblck
    logical              :: do_shum
    logical              :: do_skeb
    integer              :: skeb_npass
    logical              :: do_sfcperts
    integer              :: nsfcpert=6
    real(kind=kind_phys) :: pertz0(5)          ! mg, sfc-perts
    real(kind=kind_phys) :: pertzt(5)          ! mg, sfc-perts
    real(kind=kind_phys) :: pertshc(5)         ! mg, sfc-perts
    real(kind=kind_phys) :: pertlai(5)         ! mg, sfc-perts
    real(kind=kind_phys) :: pertalb(5)         ! mg, sfc-perts
    real(kind=kind_phys) :: pertvegf(5)        ! mg, sfc-perts
!--- tracer handling
    character(len=32), pointer :: tracer_names(:) !< array of initialized tracers from dynamic core
    integer              :: ntrac           !< number of tracers
#ifdef CCPP
    integer              :: ntracp1         !< number of tracers plus one
    integer              :: ntqv            !< tracer index for water vapor (specific humidity)
#endif
    integer              :: ntoz            !< tracer index for ozone mixing ratio
    integer              :: ntcw            !< tracer index for cloud condensate (or liquid water)
    integer              :: ntiw            !< tracer index for ice water
    integer              :: ntrw            !< tracer index for rain water
    integer              :: ntsw            !< tracer index for snow water
    integer              :: ntgl            !< tracer index for graupel
    integer              :: ntclamt         !< tracer index for cloud amount
    integer              :: ntlnc           !< tracer index for liquid number concentration
    integer              :: ntinc           !< tracer index for ice    number concentration
    integer              :: ntrnc           !< tracer index for rain   number concentration
    integer              :: ntsnc           !< tracer index for snow   number concentration
    integer              :: ntgnc           !< tracer index for graupel number concentration
    integer              :: ntke            !< tracer index for kinetic energy
    integer              :: nto             !< tracer index for oxygen ion
    integer              :: nto2            !< tracer index for oxygen
    integer              :: ntwa            !< tracer index for water friendly aerosol
    integer              :: ntia            !< tracer index for ice friendly aerosol
    integer              :: ntchm           !< number of chemical tracers
    integer              :: ntchs           !< tracer index for first chemical tracer
    logical, pointer     :: ntdiag(:) => null() !< array to control diagnostics for chemical tracers
 
    !--- derived totals for phy_f*d
    integer              :: ntot2d          !< total number of variables for phyf2d
    integer              :: ntot3d          !< total number of variables for phyf3d
    integer              :: indcld          !< location of cloud fraction in phyf3d (used only for SHOC or MG)
    integer              :: num_p2d         !< number of 2D arrays needed for microphysics
    integer              :: num_p3d         !< number of 3D arrays needed for microphysics
    integer              :: nshoc_2d        !< number of 2d fields for SHOC
    integer              :: nshoc_3d        !< number of 3d fields for SHOC
    integer              :: ncnvcld3d       !< number of convective 3d clouds fields
    integer              :: npdf3d          !< number of 3d arrays associated with pdf based clouds/microphysics
    integer              :: nctp            !< number of cloud types in Chikira-Sugiyama scheme
    integer              :: ncnvw           !< the index of cnvw in phy_f3d
    integer              :: ncnvc           !< the index of cnvc in phy_f3d
    integer              :: nleffr          !< the index of cloud liquid water effective radius in phy_f3d
    integer              :: nieffr          !< the index of ice effective radius in phy_f3d
    integer              :: nreffr          !< the index of rain effective radius in phy_f3d
    integer              :: nseffr          !< the index of snow effective radius in phy_f3d
    integer              :: ngeffr          !< the index of graupel effective radius in phy_f3d

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
#ifdef CCPP
    logical              :: first_time_step !< flag signaling first time step for time integration routine
    logical              :: restart         !< flag whether this is a coldstart (.false.) or a warmstart/restart (.true.)
    logical              :: hydrostatic     !< flag whether this is a hydrostatic or non-hydrostatic run
#endif
    integer              :: jdat(1:8)       !< current forecast date and time
                                            !< (yr, mon, day, t-zone, hr, min, sec, mil-sec)
    logical              :: iccn            !< using IN CCN forcing for MG2/3
#ifdef CCPP
    real(kind=kind_phys)          :: sec    !< seconds since model initialization
    real(kind=kind_phys), pointer :: si(:)  !< vertical sigma coordinate for model initialization
#endif

!--- IAU
    real(kind=kind_phys) :: iau_delthrs     ! iau time interval (to scale increments) in hours
    character(len=240)   :: iau_inc_files(7)! list of increment files
    real(kind=kind_phys) :: iaufhrs(7)      ! forecast hours associated with increment files
    logical              :: iau_filter_increments

#ifdef CCPP
    ! From physcons.F90, updated/set in control_initialize
    real(kind=kind_phys) :: dxinv           ! inverse scaling factor for critical relative humidity, replaces dxinv in physcons.F90
    real(kind=kind_phys) :: dxmax           ! maximum scaling factor for critical relative humidity, replaces dxmax in physcons.F90
    real(kind=kind_phys) :: dxmin           ! minimum scaling factor for critical relative humidity, replaces dxmin in physcons.F90
    real(kind=kind_phys) :: rhcmax          ! maximum critical relative humidity, replaces rhc_max in physcons.F90
#endif

    contains
      procedure :: init  => control_initialize
      procedure :: print => control_print
  end type GFS_control_type


!--------------------------------------------------------------------
! GFS_grid_type
!   grid data needed for interpolations and length-scale calculations
!--------------------------------------------------------------------
#if 0
!! \section arg_table_GFS_grid_type
!! | local_name                           | standard_name                                                   | long_name                           | units         | rank | type     |    kind   | intent | optional |
!! |--------------------------------------|-----------------------------------------------------------------|-------------------------------------|---------------|------|----------|-----------|--------|----------|
!! | GFS_Data(cdata%blk_no)%Grid%area     | cell_area                                                       | area of the grid cell               | m2            |    1 | real     | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Grid%dx       | cell_size                                                       | relative dx for the grid cell       | m             |    1 | real     | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Grid%xlat     | latitude                                                        | latitude                            | radians       |    1 | real     | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Grid%xlon     | longitude                                                       | longitude                           | radians       |    1 | real     | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Grid%coslat   | cosine_of_latitude                                              | cosine of latitude                  | none          |    1 | real     | kind_phys | in     | F        |
!! | GFS_Data(cdata%blk_no)%Grid%sinlat   | sine_of_latitude                                                | sine of latitude                    | none          |    1 | real     | kind_phys | in     | F        |
!!
#endif
  type GFS_grid_type

    real (kind=kind_phys), pointer :: xlon   (:)    => null()   !< grid longitude in radians, ok for both 0->2pi
                                                                !! or -pi -> +pi ranges
    real (kind=kind_phys), pointer :: xlat   (:)    => null()   !< grid latitude in radians, default to pi/2 ->
                                                                !! -pi/2 range, otherwise adj in subr called
    real (kind=kind_phys), pointer :: xlat_d (:)    => null()   !< grid latitude in degrees, default to 90 ->
                                                                !! -90 range, otherwise adj in subr called
    real (kind=kind_phys), pointer :: xlon_d (:)    => null()   !< grid longitude in degrees, default to 0 ->
                                                                !! 360 range, otherwise adj in subr called
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

!--- grid-related interpolation data for prognostic iccn
    real (kind=kind_phys), pointer :: ddy_ci    (:) => null()   !< interpolation     weight for iccn
    integer,               pointer :: jindx1_ci (:) => null()   !< interpolation  low index for iccn
    integer,               pointer :: jindx2_ci (:) => null()   !< interpolation high index for iccn
    real (kind=kind_phys), pointer :: ddx_ci    (:) => null()   !< interpolation     weight for iccn
    integer,               pointer :: iindx1_ci (:) => null()   !< interpolation  low index for iccn
    integer,               pointer :: iindx2_ci (:) => null()   !< interpolation high index for iccn

!--- grid-related interpolation data for prescribed aerosols
    real (kind=kind_phys), pointer :: ddy_aer    (:) => null()   !< interpolation     weight for iaerclm
    integer,               pointer :: jindx1_aer (:) => null()   !< interpolation  low index for iaerclm
    integer,               pointer :: jindx2_aer (:) => null()   !< interpolation high index for iaerclm
    real (kind=kind_phys), pointer :: ddx_aer    (:) => null()   !< interpolation     weight for iaerclm
    integer,               pointer :: iindx1_aer (:) => null()   !< interpolation  low index for iaerclm
    integer,               pointer :: iindx2_aer (:) => null()   !< interpolation high index for iaerclm
    contains
      procedure :: create   => grid_create   !<   allocate array data
  end type GFS_grid_type


!-----------------------------------------------
! GFS_tbd_type
!   data not yet assigned to a defined container
!-----------------------------------------------
#if 0
!! \section arg_table_GFS_tbd_type
!! | local_name                                                  | standard_name                                                                                  | long_name                                               | units         | rank | type    |    kind   | intent | optional |
!! |-------------------------------------------------------------|------------------------------------------------------------------------------------------------|---------------------------------------------------------|---------------|------|---------|-----------|--------|----------|
!! | GFS_Data(cdata%blk_no)%Tbd%icsdsw                           | seed_random_numbers_sw                                                                         | random seeds for sub-column cloud generators sw         | none          |    1 | integer |           | none   | F        |
!! | GFS_Data(cdata%blk_no)%Tbd%icsdlw                           | seed_random_numbers_lw                                                                         | random seeds for sub-column cloud generators lw         | none          |    1 | integer |           | none   | F        |
!! | GFS_Data(cdata%blk_no)%Tbd%ozpl                             | ozone_forcing                                                                                  | ozone forcing data                                      | various       |    3 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Tbd%h2opl                            | h2o_forcing                                                                                    | water forcing data                                      | various       |    3 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Tbd%in_nm                            | in_number_concentration                                                                        | IN number concentration                                 | kg-1?         |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Tbd%ccn_nm                           | ccn_number_concentration                                                                       | CCN number concentration                                | kg-1?         |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Tbd%aer_nm                           | aerosol_number_concentration_from_gocart_aerosol_climatology                                   | GOCART aerosol climatology number concentration         | kg-1?         |    3 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Tbd%imap                             | map_of_block_column_number_to_global_i_index                                                   | map of local index ix to global index i for this block  | none          |    1 | integer |           | none   | F        |
!! | GFS_Data(cdata%blk_no)%Tbd%jmap                             | map_of_block_column_number_to_global_j_index                                                   | map of local index ix to global index j for this block  | none          |    1 | integer |           | none   | F        |
!! | GFS_Data(cdata%blk_no)%Tbd%rann                             | random_number_array                                                                            | random number array (0-1)                               | none          |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Tbd%acv                              | accumulated_lwe_thickness_of_convective_precipitation_amount_cnvc90                            | accumulated convective rainfall amount for cnvc90 only  | m             |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Tbd%acvb                             | smallest_cloud_base_vertical_index_encountered_thus_far                                        | smallest cloud base vertical index encountered thus far | index         |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Tbd%acvt                             | largest_cloud_top_vertical_index_encountered_thus_far                                          | largest cloud top vertical index encountered thus far   | index         |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Tbd%dtdtr                            | tendency_of_air_temperature_due_to_radiative_heating_on_physics_time_step                      | temp. change due to radiative heating per time step     | K             |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Tbd%dtotprcp                         |                                                                                                | change in totprcp  (diag_type)                          |               |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Tbd%dcnvprcp                         |                                                                                                | change in cnvprcp  (diag_type)                          |               |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Tbd%drain_cpl                        | tendency_of_lwe_thickness_of_precipitation_amount_for_coupling                                 | change in rain_cpl (coupling_type)                      | m             |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Tbd%dsnow_cpl                        | tendency_of_lwe_thickness_of_snow_amount_for_coupling                                          | change in show_cpl (coupling_type)                      | m             |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Tbd%phy_fctd                         | cloud_base_mass_flux                                                                           | cloud base mass flux for CS convection                  | kg m-2 s-1    |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Tbd%phy_f2d                          |                                                                                                | 2d arrays saved for restart                             |               |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Tbd%phy_f2d(:,1)                     | surface_air_pressure_two_time_steps_back                                                       | surface air pressure two time steps back                | Pa            |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Tbd%phy_f2d(:,2)                     | surface_air_pressure_at_previous_time_step                                                     | surface air pressure at previous time step              | Pa            |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Tbd%phy_f2d(:,GFS_Control%num_p2d)   | surface_wind_enhancement_due_to_convection                                                     | surface wind enhancement due to convection              | m s-1         |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Tbd%phy_f3d                          |                                                                                                | 3d arrays saved for restart                             |               |    3 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Tbd%phy_f3d(:,:,1)                   | air_temperature_two_time_steps_back                                                            | air temperature two time steps back                     | K             |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Tbd%phy_f3d(:,:,2)                   | water_vapor_specific_humidity_two_time_steps_back                                              | water vapor specific humidity two time steps back       | kg kg-1       |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Tbd%phy_f3d(:,:,3)                   | air_temperature_at_previous_time_step                                                          | air temperature at previous time step                   | K             |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Tbd%phy_f3d(:,:,4)                   | water_vapor_specific_humidity_at_previous_time_step                                            | water vapor specific humidity at previous time step     | kg kg-1       |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Tbd%phy_f3d(:,:,GFS_Control%ncnvw)   | convective_cloud_water_mixing_ratio_in_phy_f3d                                                 | convective cloud water mixing ratio in the phy_f3d array| kg kg-1       |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Tbd%phy_f3d(:,:,GFS_Control%ncnvc)   | convective_cloud_cover_in_phy_f3d                                                              | convective cloud cover in the phy_f3d array             | frac          |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Tbd%phy_f3d(:,:,GFS_Control%ntot3d)  | kinematic_buoyancy_flux_from_shoc                                                              | upward kinematic buoyancy flux from the SHOC scheme     | K m s-1       |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Tbd%phy_f3d(:,:,GFS_Control%ntot3d-1)| atmosphere_heat_diffusivity_from_shoc                                                          | diffusivity for heat from the SHOC scheme               | m2 s-1        |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Tbd%phy_f3d(:,:,GFS_Control%ntot3d-2)| subgrid_scale_cloud_fraction_from_shoc                                                         | subgrid-scale cloud fraction from the SHOC scheme       | frac          |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Tbd%phy_f3d(:,:,GFS_Control%indcld)  | cloud_fraction_for_MG                                                                          | cloud fraction used by Morrison-Gettelman MP            | frac          |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Tbd%phy_f3d(:,:,GFS_Control%nleffr)  | effective_radius_of_stratiform_cloud_liquid_water_particle_in_um                               | eff. radius of cloud liquid water particle in micrometer| um            |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Tbd%phy_f3d(:,:,GFS_Control%nieffr)  | effective_radius_of_stratiform_cloud_ice_particle_in_um                                        | eff. radius of cloud ice water particle in micrometer   | um            |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Tbd%phy_f3d(:,:,GFS_Control%nreffr)  | effective_radius_of_stratiform_cloud_rain_particle_in_um                                       | effective radius of cloud rain particle in micrometers  | um            |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Tbd%phy_f3d(:,:,GFS_Control%nseffr)  | effective_radius_of_stratiform_cloud_snow_particle_in_um                                       | effective radius of cloud snow particle in micrometers  | um            |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Tbd%phy_f3d(:,:,GFS_Control%ngeffr)  | effective_radius_of_stratiform_cloud_graupel_particle_in_um                                    | eff. radius of cloud graupel particle in micrometer     | um            |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Tbd%htlwc                            | tendency_of_air_temperature_due_to_longwave_heating_on_radiation_time_step                     | total sky heating rate due to longwave radiation        | K s-1         |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Tbd%htlw0                            | tendency_of_air_temperature_due_to_longwave_heating_assuming_clear_sky_on_radiation_time_step  | clear sky heating rate due to longwave radiation        | K s-1         |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Tbd%htswc                            | tendency_of_air_temperature_due_to_shortwave_heating_on_radiation_time_step                    | total sky heating rate due to shortwave radiation       | K s-1         |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Tbd%htsw0                            | tendency_of_air_temperature_due_to_shortwave_heating_assuming_clear_sky_on_radiation_time_step | clear sky heating rates due to shortwave radiation      | K s-1         |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Tbd%forcet                           | temperature_tendency_due_to_dynamics                                                           | temperature tendency due to dynamics only               | K s-1         |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Tbd%forceq                           | moisture_tendency_due_to_dynamics                                                              | moisture tendency due to dynamics only                  | kg kg-1 s-1   |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Tbd%prevst                           | temperature_from_previous_timestep                                                             | temperature from previous time step                     | K             |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Tbd%prevsq                           | moisture_from_previous_timestep                                                                | moisture from previous time step                        | kg kg-1       |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Tbd%cactiv                           | conv_activity_counter                                                                          | convective activity memory                              | none          |    1 | integer |           | none   | F        |
!! | GFS_Data(cdata%blk_no)%Tbd%raincprv                         | lwe_thickness_of_explicit_rainfall_amount_from_previous_timestep                               | explicit rainfall from previous timestep                | m             |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Tbd%rainncprv                        | lwe_thickness_of_convective_precipitation_amount_from_previous_timestep                        | convective_precipitation_amount from previous timestep  | m             |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Tbd%iceprv                           | lwe_thickness_of_ice_amount_from_previous_timestep                                             | ice amount from previous timestep                       | m             |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Tbd%snowprv                          | lwe_thickness_of_snow_amount_from_previous_timestep                                            | snow amount from previous timestep                      | m             |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Tbd%graupelprv                       | lwe_thickness_of_graupel_amount_from_previous_timestep                                         | graupel amount from previous timestep                   | m             |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Tbd%CLDFRA_BL                        | subgrid_cloud_fraction_pbl                                                                     | subgrid cloud fraction from PBL scheme                  | frac          |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Tbd%QC_BL                            | subgrid_cloud_mixing_ratio_pbl                                                                 | subgrid cloud cloud mixing ratio from PBL scheme        | kg kg-1       |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Tbd%el_pbl                           | mixing_length                                                                                  | mixing length in meters                                 | m             |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Tbd%Sh3D                             | stability_function_for_heat                                                                    | stability function for heat                             | none          |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Tbd%qke                              | tke_at_mass_points                                                                             | 2 x tke at mass points                                  | m2 s-2        |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Tbd%tsq                              | t_prime_squared                                                                                | temperature fluctuation squared                         | K2            |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Tbd%qsq                              | q_prime_squared                                                                                | water vapor fluctuation squared                         | kg2 kg-2      |    2 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Tbd%cov                              | t_prime_q_prime                                                                                | covariance of temperature and moisture                  | K kg kg-1     |    2 | real    | kind_phys | none   | F        |
!!
#endif
  type GFS_tbd_type

!--- radiation random seeds
    integer,               pointer :: icsdsw   (:)     => null()  !< (rad. only) auxiliary cloud control arrays passed to main
    integer,               pointer :: icsdlw   (:)     => null()  !< (rad. only) radiations. if isubcsw/isubclw (input to init)
                                                                  !< (rad. only) are set to 2, the arrays contains provided
                                                                  !< (rad. only) random seeds for sub-column clouds generators

!--- In
    real (kind=kind_phys), pointer :: ozpl     (:,:,:) => null()  !< ozone forcing data
    real (kind=kind_phys), pointer :: h2opl    (:,:,:) => null()  !< water forcing data
    real (kind=kind_phys), pointer :: in_nm    (:,:)   => null()  !< IN number concentration
    real (kind=kind_phys), pointer :: ccn_nm   (:,:)   => null()  !< CCN number concentration
    real (kind=kind_phys), pointer :: aer_nm   (:,:,:) => null()  !< GOCART aerosol climo

    !--- active when ((.not. newsas .or. cal_pre) .and. random_clds)
#ifdef CCPP
    integer,               pointer :: imap     (:)     => null()  !< map of local index ix to global index i for this block
    integer,               pointer :: jmap     (:)     => null()  !< map of local index ix to global index j for this block
#endif
    real (kind=kind_phys), pointer :: rann     (:,:)   => null()  !< random number array (0-1)

!--- In/Out
    real (kind=kind_phys), pointer :: acv      (:)     => null()  !< array containing accumulated convective clouds
    real (kind=kind_phys), pointer :: acvb     (:)     => null()  !< arrays used by cnvc90 bottom
    real (kind=kind_phys), pointer :: acvt     (:)     => null()  !< arrays used by cnvc90 top (cnvc90.f)

!--- Stochastic physics properties calculated in physics_driver
    real (kind=kind_phys), pointer :: dtdtr     (:,:)  => null()  !< temperature change due to radiative heating per time step (K)
    real (kind=kind_phys), pointer :: dtotprcp  (:)    => null()  !< change in totprcp  (diag_type)
    real (kind=kind_phys), pointer :: dcnvprcp  (:)    => null()  !< change in cnvprcp  (diag_type)
    real (kind=kind_phys), pointer :: drain_cpl (:)    => null()  !< change in rain_cpl (coupling_type)
    real (kind=kind_phys), pointer :: dsnow_cpl (:)    => null()  !< change in show_cpl (coupling_type)

!--- phy_f*d variables needed for seamless restarts and moving data between grrad and gbphys
    real (kind=kind_phys), pointer :: phy_fctd (:,:)   => null()  !< cloud base mass flux for CS convection
    real (kind=kind_phys), pointer :: phy_f2d  (:,:)   => null()  !< 2d arrays saved for restart
    real (kind=kind_phys), pointer :: phy_f3d  (:,:,:) => null()  !< 3d arrays saved for restart

#if !defined(CCPP) || defined(HYBRID)
!--- for explicit data blocking
    integer                        :: blkno                       !< block number of this block
#endif

#ifdef CCPP
    !--- radiation variables that need to be carried over from radiation to physics
    real (kind=kind_phys), pointer :: htlwc(:,:)       => null()  !<
    real (kind=kind_phys), pointer :: htlw0(:,:)       => null()  !<
    real (kind=kind_phys), pointer :: htswc(:,:)       => null()  !<
    real (kind=kind_phys), pointer :: htsw0(:,:)       => null()  !<

    !--- dynamical forcing variables for Grell-Freitas convection
    real (kind=kind_phys), pointer :: forcet (:,:)     => null()  !<
    real (kind=kind_phys), pointer :: forceq (:,:)     => null()  !<
    real (kind=kind_phys), pointer :: prevst (:,:)     => null()  !<
    real (kind=kind_phys), pointer :: prevsq (:,:)     => null()  !<
    integer,               pointer :: cactiv   (:)     => null()  !< convective activity memory contour
    
    !---- precipitation amounts from previous time step for RUC LSM
    real (kind=kind_phys), pointer :: raincprv  (:)    => null()  !< explicit rainfall from previous timestep
    real (kind=kind_phys), pointer :: rainncprv (:)    => null()  !< convective_precipitation_amount from previous timestep
    real (kind=kind_phys), pointer :: iceprv    (:)    => null()  !< ice amount from previous timestep
    real (kind=kind_phys), pointer :: snowprv   (:)    => null()  !< snow amount from previous timestep
    real (kind=kind_phys), pointer :: graupelprv(:)    => null()  !< graupel amount from previous timestep

    !--- MYNN prognostic variables that can't be in the Intdiag or Interstitial DDTs
    real (kind=kind_phys), pointer :: CLDFRA_BL  (:,:)   => null()  !
    real (kind=kind_phys), pointer :: QC_BL      (:,:)   => null()  !
    real (kind=kind_phys), pointer :: el_pbl     (:,:)   => null()  !
    real (kind=kind_phys), pointer :: Sh3D       (:,:)   => null()  !
    real (kind=kind_phys), pointer :: qke        (:,:)   => null()  !
    real (kind=kind_phys), pointer :: tsq        (:,:)   => null()  !
    real (kind=kind_phys), pointer :: qsq        (:,:)   => null()  !
    real (kind=kind_phys), pointer :: cov        (:,:)   => null()  !
#endif

    contains
      procedure :: create  => tbd_create  !<   allocate array data
  end type GFS_tbd_type


!------------------------------------------------------------------
! GFS_cldprop_type
!  cloud properties and tendencies needed by radiation from physics
!------------------------------------------------------------------
#if 0
!! \section arg_table_GFS_cldprop_type
!! | local_name                                | standard_name                                           | long_name                                               | units         | rank | type    |    kind   | intent | optional |
!! |-------------------------------------------|---------------------------------------------------------|---------------------------------------------------------|---------------|------|---------|-----------|--------|----------|
!! | GFS_Data(cdata%blk_no)%Cldprop%cv         | fraction_of_convective_cloud                            | fraction of convective cloud                            | frac          |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Cldprop%cvt        | pressure_at_top_of_convective_cloud                     | convective cloud top pressure                           | Pa            |    1 | real    | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Cldprop%cvb        | pressure_at_bottom_of_convective_cloud                  | convective cloud bottom pressure                        | Pa            |    1 | real    | kind_phys | none   | F        |
!!
#endif
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
#if 0
!! \section arg_table_GFS_radtend_type
!! | local_name                                | standard_name                                                                                 | long_name                                               | units         | rank | type        |    kind   | intent | optional |
!! |-------------------------------------------|-----------------------------------------------------------------------------------------------|---------------------------------------------------------|---------------|------|-------------|-----------|--------|----------|
!! | GFS_Data(cdata%blk_no)%Radtend%sfcfsw     | sw_fluxes_sfc                                                                                 | sw radiation fluxes at sfc                              | W m-2         |    1 | sfcfsw_type |           | none   | F        |
!! | GFS_Data(cdata%blk_no)%Radtend%sfcflw     | lw_fluxes_sfc                                                                                 | lw radiation fluxes at sfc                              | W m-2         |    1 | sfcflw_type |           | none   | F        |
!! | GFS_Data(cdata%blk_no)%Radtend%htrsw      | tendency_of_air_temperature_due_to_shortwave_heating_on_radiation_timestep                    | total sky sw heating rate                               | K s-1         |    2 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Radtend%htrlw      | tendency_of_air_temperature_due_to_longwave_heating_on_radiation_timestep                     | total sky lw heating rate                               | K s-1         |    2 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Radtend%sfalb      | surface_diffused_shortwave_albedo                                                             | mean surface diffused sw albedo                         | frac          |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Radtend%coszen     | cosine_of_zenith_angle                                                                        | mean cos of zenith angle over rad call period           | none          |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Radtend%tsflw      | surface_midlayer_air_temperature_in_longwave_radiation                                        | surface air temp during lw calculation                  | K             |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Radtend%semis      | surface_longwave_emissivity                                                                   | surface lw emissivity in fraction                       | frac          |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Radtend%coszdg     |                                                                                               | daytime mean cosz over rad call period                  | none          |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Radtend%swhc       | tendency_of_air_temperature_due_to_shortwave_heating_assuming_clear_sky_on_radiation_timestep | clear sky sw heating rates                              | K s-1         |    2 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Radtend%lwhc       | tendency_of_air_temperature_due_to_longwave_heating_assuming_clear_sky_on_radiation_timestep  | clear sky lw heating rates                              | K s-1         |    2 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Radtend%lwhd       | tendency_of_air_temperature_due_to_longwave_heating_for_idea                                  | idea sky lw heating rates                               | K s-1         |    3 | real        | kind_phys | none   | F        |
!!
#endif
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
#if 0
!! \section arg_table_GFS_diag_type
!! | local_name                                          | standard_name                                                           | long_name                                                       | units         | rank | type        |    kind   | intent | optional |
!! |-----------------------------------------------------|-------------------------------------------------------------------------|-----------------------------------------------------------------|---------------|------|-------------|-----------|--------|----------|
!! | GFS_Data(cdata%blk_no)%Intdiag%fluxr                |                                                                         | accumulated 2-d fields, opt. includes aerosols                  |               |    2 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%topfsw               | sw_fluxes_top_atmosphere                                                | sw radiation fluxes at toa                                      | W m-2         |    1 | topfsw_type |           | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%topflw               | lw_fluxes_top_atmosphere                                                | lw radiation fluxes at top                                      | W m-2         |    1 | topflw_type |           | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%srunoff              | surface_runoff                                                          | surface water runoff (from lsm)                                 | kg m-2        |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%evbsa                | cumulative_soil_upward_latent_heat_flux_multiplied_by_timestep          | cumulative soil upward latent heat flux multiplied by timestep  | W m-2 s       |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%evcwa                | cumulative_canopy_upward_latent_heat_flu_multiplied_by_timestep         | cumulative canopy upward latent heat flux multiplied by timestep| W m-2 s       |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%snohfa               | cumulative_snow_freezing_rain_upward_latent_heat_flux_multiplied_by_timestep | cumulative latent heat flux due to snow and frz rain multiplied by timestep | W m-2 s         |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%transa               | cumulative_transpiration_flux_multiplied_by_timestep                    | cumulative total plant transpiration rate multiplied by timestep | kg m-2        |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%sbsnoa               | cumulative_snow_deposition_sublimation_upward_latent_heat_flux_multiplied_by_timestep | cumulative latent heat flux from snow depo/subl multiplied by timestep | W m-2 s       |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%snowca               | cumulative_surface_snow_area_fraction_multiplied_by_timestep            | cumulative surface snow area fraction multiplied by timestep    | s             |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%soilm                | soil_moisture_content                                                   | soil moisture                                                   | kg m-2        |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%tmpmin               | minimum_temperature_at_2m                                               | min temperature at 2m height                                    | K             |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%tmpmax               | maximum_temperature_at_2m                                               | max temperature at 2m height                                    | K             |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%dusfc                | cumulative_surface_x_momentum_flux_for_diag_multiplied_by_timestep           | cumulative sfc x momentum flux multiplied by timestep      | Pa s          |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%dvsfc                | cumulative_surface_y_momentum_flux_for_diag_multiplied_by_timestep           | cumulative sfc y momentum flux multiplied by timestep      | Pa s          |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%dtsfc                | cumulative_surface_upward_sensible_heat_flux_for_diag_multiplied_by_timestep | cumulative sfc sensible heat flux multiplied by timestep   | W m-2 s       |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%dqsfc                | cumulative_surface_upward_latent_heat_flux_for_diag_multiplied_by_timestep   | cumulative sfc latent heat flux multiplied by timestep     | W m-2 s       |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%totprcp              | accumulated_lwe_thickness_of_precipitation_amount                       | accumulated total precipitation                                 | m             |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%totice               | accumulated_lwe_thickness_of_ice_amount                                 | accumulated ice precipitation                                   | kg m-2        |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%totsnw               | accumulated_lwe_thickness_of_snow_amount                                | accumulated snow precipitation                                  | kg m-2        |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%totgrp               | accumulated_lwe_thickness_of_graupel_amount                             | accumulated graupel precipitation                               | kg m-2        |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%totprcpb             | accumulated_lwe_thickness_of_precipitation_amount_in_bucket             | accumulated total precipitation in bucket                       | m             |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%toticeb              | accumulated_lwe_thickness_of_ice_amount_in_bucket                       | accumulated ice precipitation in bucket                         | kg m-2        |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%totsnwb              | accumulated_lwe_thickness_of_snow_amount_in_bucket                      | accumulated snow precipitation in bucket                        | kg m-2        |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%totgrpb              | accumulated_lwe_thickness_of_graupel_amount_in_bucket                   | accumulated graupel precipitation in bucket                     | kg m-2        |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%gflux                | cumulative_surface_ground_heat_flux_multiplied_by_timestep              | cumulative groud conductive heat flux multiplied by timestep    | W m-2 s       |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%dlwsfc               | cumulative_surface_downwelling_longwave_flux_multiplied_by_timestep     | cumulative surface downwelling LW flux multiplied by timestep   | W m-2 s       |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%ulwsfc               | cumulative_surface_upwelling_longwave_flux_multiplied_by_timestep       | cumulative surface upwelling LW flux multiplied by timestep     | W m-2 s       |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%suntim               | duration_of_sunshine                                                    | sunshine duration time                                          | s             |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%runoff               | total_runoff                                                            | total water runoff                                              | kg m-2        |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%ep                   | cumulative_surface_upward_potential_latent_heat_flux_multiplied_by_timestep | cumulative surface upward potential latent heat flux multiplied by timestep | W m-2 s        |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%cldwrk               | cumulative_cloud_work_function                                          | cumulative cloud work function (valid only with sas)            | m2 s-1        |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%dugwd                | time_integral_of_x_stress_due_to_gravity_wave_drag                      | vertically integrated u change by OGWD                          | Pa s          |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%dvgwd                | time_integral_of_y_stress_due_to_gravity_wave_drag                      | vertically integrated v change by OGWD                          | Pa s          |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%psmean               | cumulative_surface_pressure_multiplied_by_timestep                      | cumulative surface pressure multiplied by timestep              | Pa s          |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%cnvprcp              | cumulative_lwe_thickness_of_convective_precipitation_amount             | cumulative convective precipitation                             | m             |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%cnvprcpb             | cumulative_lwe_thickness_of_convective_precipitation_amount_in_bucket   | cumulative convective precipitation in bucket                   | m             |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%spfhmin              | minimum_specific_humidity_at_2m                                         | minimum specific humidity at 2m height                          | kg kg-1       |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%spfhmax              | maximum_specific_humidity_at_2m                                         | maximum specific humidity at 2m height                          | kg kg-1       |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%u10mmax              | maximum_x_wind_at_10m                                                   | maximum x wind at 10 m                                          | m s-1         |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%v10mmax              | maximum_y_wind_at_10m                                                   | maximum y wind at 10 m                                          | m s-1         |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%wind10mmax           | maximum_wind_at_10m                                                     | maximum wind speed at 10 m                                      | m s-1         |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%u10max               | maximum_u_wind_at_10m_over_maximum_hourly_time_interval                 | maximum u wind at 10m over maximum hourly time interval         | m s-1         |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%v10max               | maximum_v_wind_at_10m_over_maximum_hourly_time_interval                 | maximum v wind at 10m over maximum hourly time interval         | m s-1         |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%spd10max             | maximum_wind_at_10m_over_maximum_hourly_time_interval                   | maximum wind at 10m over maximum hourly time interval           | m s-1         |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%rain                 | lwe_thickness_of_precipitation_amount_on_dynamics_timestep              | total rain at this time step                                    | m             |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%rainc                | lwe_thickness_of_convective_precipitation_amount_on_dynamics_timestep   | convective rain at this time step                               | m             |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%ice                  | lwe_thickness_of_ice_amount_on_dynamics_timestep                        | ice fall at this time step                                      | m             |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%snow                 | lwe_thickness_of_snow_amount_on_dynamics_timestep                       | snow fall at this time step                                     | m             |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%graupel              | lwe_thickness_of_graupel_amount_on_dynamics_timestep                    | graupel fall at this time step                                  | m             |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%u10m                 | x_wind_at_10m                                                           | 10 meter u wind speed                                           | m s-1         |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%v10m                 | y_wind_at_10m                                                           | 10 meter v wind speed                                           | m s-1         |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%dpt2m                | dewpoint_temperature_at_2m                                              | 2 meter dewpoint temperature                                    | K             |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%zlvl                 | height_above_ground_at_lowest_model_layer                               | layer 1 height above ground (not MSL)                           | m             |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%psurf                | surface_air_pressure_diag                                               | surface air pressure diagnostic                                 | Pa            |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%hpbl                 | atmosphere_boundary_layer_thickness                                     | pbl height                                                      | m             |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%pwat                 | column_precipitable_water                                               | precipitable water                                              | kg m-2        |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%t1                   | air_temperature_at_lowest_model_layer_for_diag                          | layer 1 temperature for diag                                    | K             |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%q1                   | water_vapor_specific_humidity_at_lowest_model_layer_for_diag            | layer 1 specific humidity for diag                              | kg kg-1       |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%u1                   | x_wind_at_lowest_model_layer_for_diag                                   | layer 1 x wind for diag                                         | m s-1         |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%v1                   | y_wind_at_lowest_model_layer_for_diag                                   | layer 1 y wind for diag                                         | m s-1         |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%chh                  | surface_drag_mass_flux_for_heat_and_moisture_in_air                     | thermal exchange coefficient                                    | kg m-2 s-1    |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%cmm                  | surface_drag_wind_speed_for_momentum_in_air                             | momentum exchange coefficient                                   | m s-1         |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%dlwsfci              | surface_downwelling_longwave_flux                                       | surface downwelling longwave flux at current time               | W m-2         |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%ulwsfci              | surface_upwelling_longwave_flux                                         | surface upwelling longwave flux at current time                 | W m-2         |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%dswsfci              | surface_downwelling_shortwave_flux                                      | surface downwelling shortwave flux at current time              | W m-2         |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%nswsfci              | surface_net_downwelling_shortwave_flux                                  | surface net downwelling shortwave flux at current time          | W m-2         |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%uswsfci              | surface_upwelling_shortwave_flux                                        | surface upwelling shortwave flux at current time                | W m-2         |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%dusfci               | instantaneous_surface_x_momentum_flux_for_diag                          | instantaneous sfc x momentum flux multiplied by timestep        | Pa            |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%dvsfci               | instantaneous_surface_y_momentum_flux_for_diag                          | instantaneous sfc y momentum flux multiplied by timestep        | Pa            |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%dtsfci               | instantaneous_surface_upward_sensible_heat_flux_for_diag                | instantaneous sfc sensible heat flux multiplied by timestep     | W m-2         |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%dqsfci               | instantaneous_surface_upward_latent_heat_flux_for_diag                  | instantaneous sfc latent heat flux multiplied by timestep       | W m-2         |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%gfluxi               | instantaneous_surface_ground_heat_flux                                  | instantaneous sfc ground heat flux                              | W m-2         |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%epi                  | instantaneous_surface_potential_evaporation                             | instantaneous sfc potential evaporation                         | W m-2         |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%smcwlt2              | volume_fraction_of_condensed_water_in_soil_at_wilting_point             | wilting point (volumetric)                                      | frac          |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%smcref2              | threshold_volume_fraction_of_condensed_water_in_soil                    | soil moisture threshold (volumetric)                            | frac          |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%tdomr                | dominant_rain_type                                                      | dominant rain type                                              | none          |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%tdomzr               | dominant_freezing_rain_type                                             | dominant freezing rain type                                     | none          |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%tdomip               | dominant_sleet_type                                                     | dominant sleet type                                             | none          |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%tdoms                | dominant_snow_type                                                      | dominant snow type                                              | none          |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%snowfallac           | total_accumulated_snowfall                                              | run-total snow accumulation on the ground                       | kg m-2        |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%acsnow               | accumulated_water_equivalent_of_frozen_precip                           | snow water equivalent of run-total frozen precip                | kg m-2        |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%ca_out               |                                                                         | cellular automata fraction                                      |               |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%ca_deep              |                                                                         | cellular automata fraction                                      |               |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%ca_turb              |                                                                         | cellular automata fraction                                      |               |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%ca_shal              |                                                                         | cellular automata fraction                                      |               |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%ca_rad               |                                                                         | cellular automata fraction                                      |               |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%ca_micro             |                                                                         | cellular automata fraction                                      |               |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%skebu_wts            | weights_for_stochastic_skeb_perturbation_of_x_wind_flipped              | weights for stochastic skeb perturbation of x wind, flipped     | none          |    2 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%skebv_wts            | weights_for_stochastic_skeb_perturbation_of_y_wind_flipped              | weights for stochastic skeb perturbation of y wind, flipped     | none          |    2 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%sppt_wts             | weights_for_stochastic_sppt_perturbation_flipped                        | weights for stochastic sppt perturbation, flipped               | none          |    2 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%shum_wts             | weights_for_stochastic_shum_perturbation_flipped                        | weights for stochastic shum perturbation, flipped               | none          |    2 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%zmtnblck             | level_of_dividing_streamline                                            | level of the dividing streamline                                | none          |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%du3dt                |                                                                         | u momentum change due to physics                                |               |    3 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%du3dt(:,:,1)         | cumulative_change_in_x_wind_due_to_PBL                                  | cumulative change in x wind due to PBL                          | m s-1         |    2 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%du3dt(:,:,2)         | cumulative_change_in_x_wind_due_to_orographic_gravity_wave_drag         | cumulative change in x wind due to orographic gravity wave drag | m s-1         |    2 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%du3dt(:,:,3)         | cumulative_change_in_x_wind_due_to_deep_convection                      | cumulative change in x wind due to deep convection              | m s-1         |    2 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%du3dt(:,:,4)         | cumulative_change_in_x_wind_due_to_convective_gravity_wave_drag         | cumulative change in x wind due to convective gravity wave drag | m s-1         |    2 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%dv3dt                |                                                                         | v momentum change due to physics                                |               |    3 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%dv3dt(:,:,1)         | cumulative_change_in_y_wind_due_to_PBL                                  | cumulative change in y wind due to PBL                          | m s-1         |    2 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%dv3dt(:,:,2)         | cumulative_change_in_y_wind_due_to_orographic_gravity_wave_drag         | cumulative change in y wind due to orographic gravity wave drag | m s-1         |    2 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%dv3dt(:,:,3)         | cumulative_change_in_y_wind_due_to_deep_convection                      | cumulative change in y wind due to deep convection              | m s-1         |    2 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%dv3dt(:,:,4)         | cumulative_change_in_y_wind_due_to_convective_gravity_wave_drag         | cumulative change in y wind due to convective gravity wave drag | m s-1         |    2 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%dt3dt                |                                                                         | temperature change due to physics                               |               |    3 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%dt3dt(:,:,1)         | cumulative_change_in_temperature_due_to_longwave_radiation              | cumulative change in temperature due to longwave radiation      | K             |    2 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%dt3dt(:,:,2)         | cumulative_change_in_temperature_due_to_shortwave_radiation             | cumulative change in temperature due to shortwave radiation     | K             |    2 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%dt3dt(:,:,3)         | cumulative_change_in_temperature_due_to_PBL                               | cumulative change in temperature due to PBL                            | K             |    2 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%dt3dt(:,:,4)         | cumulative_change_in_temperature_due_to_deep_convection                   | cumulative change in temperature due to deep conv.                     | K             |    2 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%dt3dt(:,:,5)         | cumulative_change_in_temperature_due_to_shal_convection                   | cumulative change in temperature due to shal conv.                     | K             |    2 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%dt3dt(:,:,6)         | cumulative_change_in_temperature_due_to_microphysics                      | cumulative change in temperature due to microphysics                   | K             |    2 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%dt3dt(:,:,7)         | cumulative_change_in_temperature_due_to_orographic_gravity_wave_drag      | cumulative change in temperature due to orographic gravity wave drag   | K             |    2 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%dq3dt                | cumulative_change_in_water_vapor_specific_humidity_due_to_physics         | cumulative change in water vapor specific humidity due to physics      | kg kg-1       |    3 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%dq3dt(:,:,1)         | cumulative_change_in_water_vapor_specific_humidity_due_to_PBL             | cumulative change in water vapor specific humidity due to PBL          | kg kg-1       |    2 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%dq3dt(:,:,2)         | cumulative_change_in_water_vapor_specific_humidity_due_to_deep_convection | cumulative change in water vapor specific humidity due to deep conv.   | kg kg-1       |    2 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%dq3dt(:,:,3)         | cumulative_change_in_water_vapor_specific_humidity_due_to_shal_convection | cumulative change in water vapor specific humidity due to shal conv.   | kg kg-1       |    2 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%dq3dt(:,:,4)         | cumulative_change_in_water_vapor_specific_humidity_due_to_microphysics    | cumulative change in water vapor specific humidity due to microphysics | kg kg-1       |    2 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%dq3dt(:,:,5)         | cumulative_change_in_ozone_mixing_ratio_due_to_PBL                        | cumulative change in ozone mixing ratio due to PBL                     | kg kg-1       |    2 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%dq3dt(:,:,6)         | cumulative_change_in_ozone_concentration_due_to_production_and_loss_rate  | cumulative change in ozone concentration due to production and loss rate | kg kg-1     |    2 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%dq3dt(:,:,7)         | cumulative_change_in_ozone_concentration_due_to_ozone_mixing_ratio        | cumulative change in ozone concentration due to ozone mixing ratio     | kg kg-1       |    2 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%dq3dt(:,:,8)         | cumulative_change_in_ozone_concentration_due_to_temperature               | cumulative change in ozone concentration due to temperature            | kg kg-1       |    2 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%dq3dt(:,:,9)         | cumulative_change_in_ozone_concentration_due_to_overhead_ozone_column     | cumulative change in ozone concentration due to overhead ozone column  | kg kg-1       |    2 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%refdmax              | maximum_reflectivity_at_1km_agl_over_maximum_hourly_time_interval         | maximum reflectivity at 1km agl over maximum hourly time interval      | dBZ           |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%refdmax263k          | maximum_reflectivity_at_minus10c_over_maximum_hourly_time_interval        | maximum reflectivity at minus10c over maximum hourly time interval     | dBZ           |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%t02max               | maximum_temperature_at_2m_over_maximum_hourly_time_interval               | maximum temperature at 2m over maximum hourly time interval            | K             |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%t02min               | minimum_temperature_at_2m_over_maximum_hourly_time_interval               | minumum temperature at 2m over maximum hourly time interval            | K             |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%rh02max              | maximum_relative_humidity_at_2m_over_maximum_hourly_time_interval         | maximum relative humidity at 2m over maximum hourly time interval      | %             |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%rh02min              | minimum_relative_humidity_at_2m_over_maximum_hourly_time_interval         | minumum relative humidity at 2m over maximum hourly time interval      | %             |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%upd_mf               | cumulative_atmosphere_updraft_convective_mass_flux                        | cumulative updraft mass flux                                           | Pa            |    2 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%dwn_mf               | cumulative_atmosphere_downdraft_convective_mass_flux                      | cumulative downdraft mass flux                                         | Pa            |    2 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%det_mf               | cumulative_atmosphere_detrainment_convective_mass_flux                    | cumulative detrainment mass flux                                       | Pa            |    2 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%cldcov               |                                                                           | instantaneous 3D cloud fraction                                        |               |    2 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%refl_10cm            | radar_reflectivity_10cm                                                   | instantaneous refl_10cm                                                | dBZ           |    2 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%duem                 | instantaneous_dust_emission_flux                                          | instantaneous dust emission flux                                       | kg m-2 s-1    |    2 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%ssem                 | instantaneous_seasalt_emission_flux                                       | instantaneous sea salt emission flux                                   | kg m-2 s-1    |    2 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%sedim                | instantaneous_sedimentation                                               | instantaneous sedimentation                                            | kg m-2 s-1    |    2 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%drydep               | instantaneous_dry_deposition                                              | instantaneous dry deposition                                           | kg m-2 s-1    |    2 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%wetdpl               | instantaneous_large-scale_wet_deposition                                  | instantaneous large-scale wet deposition                               | kg m-2 s-1    |    2 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%wetdpc               | instantaneous_convective-scale_wet_deposition                             | instantaneous convective-scale wet deposition                          | kg m-2 s-1    |    2 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%abem                 | instantaneous_anthopogenic_and_biomass_burning_emissions | instantaneous anthopogenic and biomass burning emissions for black carbon, organic carbon, and sulfur dioxide | ug m-2 s-1 | 2 | real | kind_phys | none | F |
!! | GFS_Data(cdata%blk_no)%Intdiag%aecm                 | instantaneous_aerosol_column_mass_densities              | instantaneous aerosol column mass densities for pm2.5, black carbon, organic carbon, sulfate, dust, sea salt  | g m-2      | 2 | real | kind_phys | none | F |
!! | GFS_Data(cdata%blk_no)%Intdiag%edmf_a               | emdf_updraft_area                                                         | updraft area from mass flux scheme                                     | frac          |    2 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%edmf_w               | emdf_updraft_vertical_velocity                                            | updraft vertical velocity from mass flux scheme                        | m s-1         |    2 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%edmf_qt              | emdf_updraft_total_water                                                  | updraft total water from mass flux scheme                              | kg kg-1       |    2 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%edmf_thl             | emdf_updraft_theta_l                                                      | updraft theta-l from mass flux scheme                                  | K             |    2 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%edmf_ent             | emdf_updraft_entrainment_rate                                             | updraft entranment rate from mass flux scheme                          | s-1           |    2 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%edmf_qc              | emdf_updraft_cloud_water                                                  | updraft cloud water from mass flux scheme                              | kg kg-1       |    2 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%nupdraft             | number_of_plumes                                                          | number of plumes per grid column                                       | count         |    1 | integer     |           | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%maxMF                | maximum_mass_flux                                                         | maximum mass flux within a column                                      | m s-1         |    1 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%ktop_shallow         | k_level_of_highest_reaching_plume                                         | k-level of highest reaching plume                                      | count         |    1 | integer     |           | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%exch_h               | atmosphere_heat_diffusivity_for_mynnpbl                                   | diffusivity for heat for MYNN PBL (defined for all mass levels)        | m2 s-1        |    2 | real        | kind_phys | none   | F        |
!! | GFS_Data(cdata%blk_no)%Intdiag%exch_m               | atmosphere_momentum_diffusivity_for_mynnpbl                               | diffusivity for momentum for MYNN PBL (defined for all mass levels)    | m2 s-1        |    2 | real        | kind_phys | none   | F        |
!!
#endif
  type GFS_diag_type

!! Input/Output only in radiation
    real (kind=kind_phys), pointer :: fluxr(:,:)     => null()   !< to save time accumulated 2-d fields defined as:!
                                                                 !< hardcoded field indices, opt. includes aerosols!
    type (topfsw_type),    pointer :: topfsw(:)      => null()   !< sw radiation fluxes at toa, components:
                                               !       %upfxc    - total sky upward sw flux at toa (w/m**2)
                                               !       %dnfxc    - total sky downward sw flux at toa (w/m**2)
                                               !       %upfx0    - clear sky upward sw flux at toa (w/m**2)
    type (topflw_type),    pointer :: topflw(:)      => null()   !< lw radiation fluxes at top, component:
                                               !        %upfxc    - total sky upward lw flux at toa (w/m**2)
                                               !        %upfx0    - clear sky upward lw flux at toa (w/m**2)

! Input/output - used by physics
    real (kind=kind_phys), pointer :: srunoff(:)     => null()   !< surface water runoff (from lsm)
    real (kind=kind_phys), pointer :: evbsa  (:)     => null()   !< noah lsm diagnostics
    real (kind=kind_phys), pointer :: evcwa  (:)     => null()   !< noah lsm diagnostics
    real (kind=kind_phys), pointer :: snohfa (:)     => null()   !< noah lsm diagnostics
    real (kind=kind_phys), pointer :: transa (:)     => null()   !< noah lsm diagnostics
    real (kind=kind_phys), pointer :: sbsnoa (:)     => null()   !< noah lsm diagnostics
    real (kind=kind_phys), pointer :: snowca (:)     => null()   !< noah lsm diagnostics
#ifdef CCPP
    real (kind=kind_phys), pointer :: snowfallac (:) => null()   !< ruc lsm diagnostics
    real (kind=kind_phys), pointer :: acsnow     (:) => null()   !< ruc lsm diagnostics
#endif
    real (kind=kind_phys), pointer :: soilm  (:)     => null()   !< soil moisture
    real (kind=kind_phys), pointer :: tmpmin (:)     => null()   !< min temperature at 2m height (k)
    real (kind=kind_phys), pointer :: tmpmax (:)     => null()   !< max temperature at 2m height (k)
    real (kind=kind_phys), pointer :: dusfc  (:)     => null()   !< u component of surface stress
    real (kind=kind_phys), pointer :: dvsfc  (:)     => null()   !< v component of surface stress
    real (kind=kind_phys), pointer :: dtsfc  (:)     => null()   !< sensible heat flux (w/m2)
    real (kind=kind_phys), pointer :: dqsfc  (:)     => null()   !< latent heat flux (w/m2)
    real (kind=kind_phys), pointer :: totprcp(:)     => null()   !< accumulated total precipitation (kg/m2)
    real (kind=kind_phys), pointer :: totprcpb(:)    => null()   !< accumulated total precipitation in bucket(kg/m2)
    real (kind=kind_phys), pointer :: gflux  (:)     => null()   !< groud conductive heat flux
    real (kind=kind_phys), pointer :: dlwsfc (:)     => null()   !< time accumulated sfc dn lw flux ( w/m**2 )
    real (kind=kind_phys), pointer :: ulwsfc (:)     => null()   !< time accumulated sfc up lw flux ( w/m**2 )
    real (kind=kind_phys), pointer :: suntim (:)     => null()   !< sunshine duration time (s)
    real (kind=kind_phys), pointer :: runoff (:)     => null()   !< total water runoff
    real (kind=kind_phys), pointer :: ep     (:)     => null()   !< potential evaporation
    real (kind=kind_phys), pointer :: cldwrk (:)     => null()   !< cloud workfunction (valid only with sas)
    real (kind=kind_phys), pointer :: dugwd  (:)     => null()   !< vertically integrated u change by OGWD
    real (kind=kind_phys), pointer :: dvgwd  (:)     => null()   !< vertically integrated v change by OGWD
    real (kind=kind_phys), pointer :: psmean (:)     => null()   !< surface pressure (kPa)
    real (kind=kind_phys), pointer :: cnvprcp(:)     => null()   !< accumulated convective precipitation (kg/m2)
    real (kind=kind_phys), pointer :: cnvprcpb(:)    => null()   !< accumulated convective precipitation in bucket (kg/m2)
    real (kind=kind_phys), pointer :: spfhmin(:)     => null()   !< minimum specific humidity
    real (kind=kind_phys), pointer :: spfhmax(:)     => null()   !< maximum specific humidity
    real (kind=kind_phys), pointer :: u10mmax(:)     => null()   !< maximum u-wind
    real (kind=kind_phys), pointer :: v10mmax(:)     => null()   !< maximum v-wind
    real (kind=kind_phys), pointer :: wind10mmax(:)  => null()   !< maximum wind speed
    real (kind=kind_phys), pointer :: u10max(:)      => null()   !< maximum u-wind used with avg_max_length
    real (kind=kind_phys), pointer :: v10max(:)      => null()   !< maximum v-wind used with avg_max_length
    real (kind=kind_phys), pointer :: spd10max(:)    => null()   !< maximum wind speed used with avg_max_length
    real (kind=kind_phys), pointer :: rain   (:)     => null()   !< total rain at this time step
    real (kind=kind_phys), pointer :: rainc  (:)     => null()   !< convective rain at this time step
    real (kind=kind_phys), pointer :: ice    (:)     => null()   !< ice fall at this time step
    real (kind=kind_phys), pointer :: snow   (:)     => null()   !< snow fall at this time step
    real (kind=kind_phys), pointer :: graupel(:)     => null()   !< graupel fall at this time step
    real (kind=kind_phys), pointer :: totice (:)     => null()   !< accumulated ice precipitation (kg/m2)
    real (kind=kind_phys), pointer :: totsnw (:)     => null()   !< accumulated snow precipitation (kg/m2)
    real (kind=kind_phys), pointer :: totgrp (:)     => null()   !< accumulated graupel precipitation (kg/m2)
    real (kind=kind_phys), pointer :: toticeb(:)     => null()   !< accumulated ice precipitation in bucket (kg/m2)
    real (kind=kind_phys), pointer :: totsnwb(:)     => null()   !< accumulated snow precipitation in bucket (kg/m2)
    real (kind=kind_phys), pointer :: totgrpb(:)     => null()   !< accumulated graupel precipitation in bucket (kg/m2)

#ifdef CCPP
    !--- MYNN variables                                              
    real (kind=kind_phys), pointer :: edmf_a     (:,:)   => null()  !
    real (kind=kind_phys), pointer :: edmf_w     (:,:)   => null()  !
    real (kind=kind_phys), pointer :: edmf_qt    (:,:)   => null()  !
    real (kind=kind_phys), pointer :: edmf_thl   (:,:)   => null()  !
    real (kind=kind_phys), pointer :: edmf_ent   (:,:)   => null()  !
    real (kind=kind_phys), pointer :: edmf_qc    (:,:)   => null()  !
    real (kind=kind_phys), pointer :: maxMF       (:)    => null()  !
    integer, pointer               :: nupdraft    (:)    => null()  !
    integer, pointer               :: ktop_shallow (:)   => null()  !
    real (kind=kind_phys), pointer :: exch_h     (:,:)   => null()  !
    real (kind=kind_phys), pointer :: exch_m     (:,:)   => null()  !
#endif

! Output - only in physics
    real (kind=kind_phys), pointer :: u10m   (:)     => null()   !< 10 meter u/v wind speed
    real (kind=kind_phys), pointer :: v10m   (:)     => null()   !< 10 meter u/v wind speed
    real (kind=kind_phys), pointer :: dpt2m  (:)     => null()   !< 2 meter dew point temperature
    real (kind=kind_phys), pointer :: zlvl   (:)     => null()   !< layer 1 height (m)
    real (kind=kind_phys), pointer :: psurf  (:)     => null()   !< surface pressure (Pa)
    real (kind=kind_phys), pointer :: hpbl   (:)     => null()   !< pbl height (m)
    real (kind=kind_phys), pointer :: pwat   (:)     => null()   !< precipitable water
    real (kind=kind_phys), pointer :: t1     (:)     => null()   !< layer 1 temperature (K)
    real (kind=kind_phys), pointer :: q1     (:)     => null()   !< layer 1 specific humidity (kg/kg)
    real (kind=kind_phys), pointer :: u1     (:)     => null()   !< layer 1 zonal wind (m/s)
    real (kind=kind_phys), pointer :: v1     (:)     => null()   !< layer 1 merdional wind (m/s)
    real (kind=kind_phys), pointer :: chh    (:)     => null()   !< thermal exchange coefficient
    real (kind=kind_phys), pointer :: cmm    (:)     => null()   !< momentum exchange coefficient
    real (kind=kind_phys), pointer :: dlwsfci(:)     => null()   !< instantaneous sfc dnwd lw flux ( w/m**2 )
    real (kind=kind_phys), pointer :: ulwsfci(:)     => null()   !< instantaneous sfc upwd lw flux ( w/m**2 )
    real (kind=kind_phys), pointer :: dswsfci(:)     => null()   !< instantaneous sfc dnwd sw flux ( w/m**2 )
#ifdef CCPP
    real (kind=kind_phys), pointer :: nswsfci(:)     => null()   !< instantaneous sfc net dnwd sw flux ( w/m**2 )
#endif
    real (kind=kind_phys), pointer :: uswsfci(:)     => null()   !< instantaneous sfc upwd sw flux ( w/m**2 )
    real (kind=kind_phys), pointer :: dusfci (:)     => null()   !< instantaneous u component of surface stress
    real (kind=kind_phys), pointer :: dvsfci (:)     => null()   !< instantaneous v component of surface stress
    real (kind=kind_phys), pointer :: dtsfci (:)     => null()   !< instantaneous sfc sensible heat flux
    real (kind=kind_phys), pointer :: dqsfci (:)     => null()   !< instantaneous sfc latent heat flux
    real (kind=kind_phys), pointer :: gfluxi (:)     => null()   !< instantaneous sfc ground heat flux
    real (kind=kind_phys), pointer :: epi    (:)     => null()   !< instantaneous sfc potential evaporation
    real (kind=kind_phys), pointer :: smcwlt2(:)     => null()   !< wilting point (volumetric)
    real (kind=kind_phys), pointer :: smcref2(:)     => null()   !< soil moisture threshold (volumetric)
    real (kind=kind_phys), pointer :: tdomr  (:)     => null()   !< dominant accumulated rain type
    real (kind=kind_phys), pointer :: tdomzr (:)     => null()   !< dominant accumulated freezing rain type
    real (kind=kind_phys), pointer :: tdomip (:)     => null()   !< dominant accumulated sleet type
    real (kind=kind_phys), pointer :: tdoms  (:)     => null()   !< dominant accumulated snow type

    real (kind=kind_phys), pointer :: ca_out  (:)    => null()   !< cellular automata fraction
    real (kind=kind_phys), pointer :: ca_deep (:)    => null()   !< cellular automata fraction
    real (kind=kind_phys), pointer :: ca_turb (:)    => null()   !< cellular automata fraction
    real (kind=kind_phys), pointer :: ca_shal (:)    => null()   !< cellular automata fraction
    real (kind=kind_phys), pointer :: ca_rad  (:)    => null()   !< cellular automata fraction
    real (kind=kind_phys), pointer :: ca_micro(:)    => null()   !< cellular automata fraction

    real (kind=kind_phys), pointer :: skebu_wts(:,:) => null()   !< 10 meter u wind speed
    real (kind=kind_phys), pointer :: skebv_wts(:,:) => null()   !< 10 meter v wind speed
    real (kind=kind_phys), pointer :: sppt_wts(:,:)  => null()   !<
    real (kind=kind_phys), pointer :: shum_wts(:,:)  => null()   !<
!--- accumulated quantities for 3D diagnostics
    real (kind=kind_phys), pointer :: zmtnblck(:)    => null()   !< mountain blocking evel
    real (kind=kind_phys), pointer :: du3dt (:,:,:)  => null()   !< u momentum change due to physics
    real (kind=kind_phys), pointer :: dv3dt (:,:,:)  => null()   !< v momentum change due to physics
    real (kind=kind_phys), pointer :: dt3dt (:,:,:)  => null()   !< temperature change due to physics
    real (kind=kind_phys), pointer :: dq3dt (:,:,:)  => null()   !< moisture change due to physics
    real (kind=kind_phys), pointer :: refdmax (:)    => null()   !< max hourly 1-km agl reflectivity
    real (kind=kind_phys), pointer :: refdmax263k (:)=> null()   !< max hourly -10C reflectivity
    real (kind=kind_phys), pointer :: t02max (:)     => null()   !< max hourly 2m T
    real (kind=kind_phys), pointer :: t02min (:)     => null()   !< min hourly 2m T
    real (kind=kind_phys), pointer :: rh02max (:)    => null()   !< max hourly 2m RH
    real (kind=kind_phys), pointer :: rh02min (:)    => null()   !< min hourly 2m RH
!--- accumulated quantities for 3D diagnostics
    real (kind=kind_phys), pointer :: upd_mf (:,:)   => null()   !< instantaneous convective updraft mass flux
    real (kind=kind_phys), pointer :: dwn_mf (:,:)   => null()   !< instantaneous convective downdraft mass flux
    real (kind=kind_phys), pointer :: det_mf (:,:)   => null()   !< instantaneous convective detrainment mass flux
    real (kind=kind_phys), pointer :: cldcov (:,:)   => null()   !< instantaneous 3D cloud fraction

    !--- MP quantities for 3D diagnositics
    real (kind=kind_phys), pointer :: refl_10cm(:,:) => null()   !< instantaneous refl_10cm

    !--- Output diagnostics for coupled chemistry
    real (kind=kind_phys), pointer :: duem  (:,:) => null()    !< instantaneous dust emission flux                             ( kg/m**2/s )
    real (kind=kind_phys), pointer :: ssem  (:,:) => null()    !< instantaneous sea salt emission flux                         ( kg/m**2/s )
    real (kind=kind_phys), pointer :: sedim (:,:) => null()    !< instantaneous sedimentation                                  ( kg/m**2/s )
    real (kind=kind_phys), pointer :: drydep(:,:) => null()    !< instantaneous dry deposition                                 ( kg/m**2/s )
    real (kind=kind_phys), pointer :: wetdpl(:,:) => null()    !< instantaneous large-scale wet deposition                     ( kg/m**2/s )
    real (kind=kind_phys), pointer :: wetdpc(:,:) => null()    !< instantaneous convective-scale wet deposition                ( kg/m**2/s )
    real (kind=kind_phys), pointer :: abem  (:,:) => null()    !< instantaneous anthopogenic and biomass burning emissions
                                                               !< for black carbon, organic carbon, and sulfur dioxide         ( ug/m**2/s )
    real (kind=kind_phys), pointer :: aecm  (:,:) => null()    !< instantaneous aerosol column mass densities for
                                                               !< pm2.5, black carbon, organic carbon, sulfate, dust, sea salt ( g/m**2 )
    contains
      procedure :: create    => diag_create
      procedure :: rad_zero  => diag_rad_zero
      procedure :: phys_zero => diag_phys_zero
      procedure :: chem_init => diag_chem_init
  end type GFS_diag_type

#ifdef CCPP
!---------------------------------------------------------------------
! GFS_interstitial_type
!> This container defines fields required for interstitial code in CCPP schemes, previously
!   in GFS_{physics,radiation}_driver.F90
!---------------------------------------------------------------------
#if 0
!! \section arg_table_GFS_interstitial_type
!! | local_name                                                    | standard_name                                                                                  | long_name                                                                           | units         | rank | type        |    kind   | intent | optional |
!! |---------------------------------------------------------------|------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------|---------------|------|-------------|-----------|--------|----------|
!! | GFS_Interstitial(cdata%thrd_no)%adjnirbmd                     | surface_downwelling_direct_near_infrared_shortwave_flux                                        | surface downwelling beam near-infrared shortwave flux at current time               | W m-2         |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%adjnirbmu                     | surface_upwelling_direct_near_infrared_shortwave_flux                                          | surface upwelling beam near-infrared shortwave flux at current time                 | W m-2         |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%adjnirdfd                     | surface_downwelling_diffuse_near_infrared_shortwave_flux                                       | surface downwelling diffuse near-infrared shortwave flux at current time            | W m-2         |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%adjnirdfu                     | surface_upwelling_diffuse_near_infrared_shortwave_flux                                         | surface upwelling diffuse near-infrared shortwave flux at current time              | W m-2         |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%adjvisbmd                     | surface_downwelling_direct_ultraviolet_and_visible_shortwave_flux                              | surface downwelling beam ultraviolet plus visible shortwave flux at current time    | W m-2         |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%adjvisbmu                     | surface_upwelling_direct_ultraviolet_and_visible_shortwave_flux                                | surface upwelling beam ultraviolet plus visible shortwave flux at current time      | W m-2         |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%adjvisdfu                     | surface_upwelling_diffuse_ultraviolet_and_visible_shortwave_flux                               | surface upwelling diffuse ultraviolet plus visible shortwave flux at current time   | W m-2         |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%adjvisdfd                     | surface_downwelling_diffuse_ultraviolet_and_visible_shortwave_flux                             | surface downwelling diffuse ultraviolet plus visible shortwave flux at current time | W m-2         |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%aerodp                        | atmosphere_optical_thickness_due_to_ambient_aerosol_particles                                  | vertical integrated optical depth for various aerosol species                       | none          |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%alb1d                         | surface_albedo_perturbation                                                                    | surface albedo perturbation                                                         | frac          |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%bexp1d                        | perturbation_of_soil_type_b_parameter                                                          | perturbation of soil type "b" parameter                                             | frac          |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%cd                            | surface_drag_coefficient_for_momentum_in_air                                                   | surface exchange coeff for momentum                                                 | none          |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%cdq                           | surface_drag_coefficient_for_heat_and_moisture_in_air                                          | surface exchange coeff heat & moisture                                              | none          |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%cf_upi                        | convective_cloud_fraction_for_microphysics                                                     | convective cloud fraction for microphysics                                          | frac          |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%clcn                          | convective_cloud_volume_fraction                                                               | convective cloud volume fraction                                                    | frac          |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%cldf                          | cloud_area_fraction                                                                            | fraction of grid box area in which updrafts occur                                   | frac          |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%cldsa                         | cloud_area_fraction_for_radiation                                                              | fraction of clouds for low, middle, high, total and BL                              | frac          |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%cldtaulw                      | cloud_optical_depth_layers_at_10mu_band                                                        | approx 10mu band layer cloud optical depth                                          | none          |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%cldtausw                      | cloud_optical_depth_layers_at_0.55mu_band                                                      | approx .55mu band layer cloud optical depth                                         | none          |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%cld1d                         | cloud_work_function                                                                            | cloud work function                                                                 | m2 s-2        |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%clouds                        |                                                                                                |                                                                                     | various       |    3 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%clouds(:,:,1)                 | total_cloud_fraction                                                                           | layer total cloud fraction                                                          | frac          |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%clouds(:,:,2)                 | cloud_liquid_water_path                                                                        | layer cloud liquid water path                                                       | g m-2         |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%clouds(:,:,3)                 | mean_effective_radius_for_liquid_cloud                                                         | mean effective radius for liquid cloud                                              | micron        |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%clouds(:,:,4)                 | cloud_ice_water_path                                                                           | layer cloud ice water path                                                          | g m-2         |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%clouds(:,:,5)                 | mean_effective_radius_for_ice_cloud                                                            | mean effective radius for ice cloud                                                 | micron        |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%clouds(:,:,6)                 | cloud_rain_water_path                                                                          | cloud rain water path                                                               | g m-2         |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%clouds(:,:,7)                 | mean_effective_radius_for_rain_drop                                                            | mean effective radius for rain drop                                                 | micron        |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%clouds(:,:,8)                 | cloud_snow_water_path                                                                          | cloud snow water path                                                               | g m-2         |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%clouds(:,:,9)                 | mean_effective_radius_for_snow_flake                                                           | mean effective radius for snow flake                                                | micron        |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%clw                           | convective_transportable_tracers                                                               | array to contain cloud water and other convective trans. tracers                    | kg kg-1       |    3 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%clw(:,:,1)                    | ice_water_mixing_ratio_convective_transport_tracer                            | moist (dry+vapor, no condensates) mixing ratio of ice water in the convectively transported tracer array | kg kg-1   |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%clw(:,:,2)                    | cloud_condensed_water_mixing_ratio_convective_transport_tracer | moist (dry+vapor, no condensates) mixing ratio of cloud water (condensate) in the convectively transported tracer array | kg kg-1   |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%clw(:,:,GFS_Interstitial(cdata%thrd_no)%ntk) | turbulent_kinetic_energy_convective_transport_tracer                            | turbulent kinetic energy in the convectively transported tracer array               | m2 s-2        |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%clx                           | fraction_of_grid_box_with_subgrid_orography_higher_than_critical_height                        | frac. of grid box with by subgrid orography higher than critical height             | frac          |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%cnv_dqldt                     | tendency_of_cloud_water_due_to_convective_microphysics                                         | tendency of cloud water due to convective microphysics                              | kg m-2 s-1    |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%cnv_fice                      | ice_fraction_in_convective_tower                                                               | ice fraction in convective tower                                                    | frac          |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%cnv_mfd                       | detrained_mass_flux                                                                            | detrained mass flux                                                                 | kg m-2 s-1    |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%cnv_ndrop                     | number_concentration_of_cloud_liquid_water_particles_for_detrainment                           | droplet number concentration in convective detrainment                              | m-3           |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%cnv_nice                      | number_concentration_of_ice_crystals_for_detrainment                                           | crystal number concentration in convective detrainment                              | m-3           |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%cnvc                          | convective_cloud_cover                                                                         | convective cloud cover                                                              | frac          |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%cnvw                          | convective_cloud_water_mixing_ratio                                                            | moist convective cloud water mixing ratio                                           | kg kg-1       |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%ctei_r                        | cloud_top_entrainment_instability_value                                                        | cloud top entrainment instability value                                             | none          |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%ctei_rml                      | grid_sensitive_critical_cloud_top_entrainment_instability_criteria                             | grid sensitive critical cloud top entrainment instability criteria                  | none          |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%cumabs                        | maximum_column_heating_rate                                                                    | maximum heating rate in column                                                      | K s-1         |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%dd_mf                         | instantaneous_atmosphere_downdraft_convective_mass_flux                                        | (downdraft mass flux) * delt                                                        | kg m-2        |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%de_lgth                       | cloud_decorrelation_length                                                                     | cloud decorrelation length                                                          | km            |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%del                           | air_pressure_difference_between_midlayers                                                      | air pressure difference between midlayers                                           | Pa            |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%del_gz                        | geopotential_difference_between_midlayers_divided_by_midlayer_virtual_temperature              | difference between mid-layer geopotentials divided by mid-layer virtual temperature | m2 s-2 K-1    |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%delr                          | layer_pressure_thickness_for_radiation                                                         | layer pressure thickness on radiation levels                                        | hPa           |    2 | real        | kind_phys | none   | F        | 
!! | GFS_Interstitial(cdata%thrd_no)%dkt                           | atmosphere_heat_diffusivity                                                                    | diffusivity for heat                                                                | m2 s-1        |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%dlength                       | characteristic_grid_length_scale                                                               | representative horizontal length scale of grid box                                  | m             |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%dqdt                          | tendency_of_tracers_due_to_model_physics                                                       | updated tendency of the tracers due to model physics                                | kg kg-1 s-1   |    3 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%dqdt(:,:,GFS_Control%ntqv)    | tendency_of_water_vapor_specific_humidity_due_to_model_physics                                 | water vapor specific humidity tendency due to model physics                         | kg kg-1 s-1   |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%dqdt(:,:,GFS_Control%ntcw)    | tendency_of_liquid_cloud_water_mixing_ratio_due_to_model_physics                               | cloud condensed water mixing ratio tendency due to model physics                    | kg kg-1 s-1   |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%dqdt(:,:,GFS_Control%ntiw)    | tendency_of_ice_cloud_water_mixing_ratio_due_to_model_physics                                  | cloud condensed water mixing ratio tendency due to model physics                    | kg kg-1 s-1   |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%dqdt(:,:,GFS_Control%ntoz)    | tendency_of_ozone_mixing_ratio_due_to_model_physics                                            | ozone mixing ratio tendency due to model physics                                    | kg kg-1 s-1   |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%dqdt(:,:,GFS_Control%ntlnc)   | tendency_of_cloud_droplet_number_concentration_due_to_model_physics                            | number concentration of cloud droplets (liquid) tendency due to model physics       | kg-1 s-1      |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%dqdt(:,:,GFS_Control%ntinc)   | tendency_of_ice_number_concentration_due_to_model_physics                                      | number concentration of ice tendency due to model physics                           | kg-1 s-1      |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%dqdt(:,:,GFS_Control%ntwa)    | tendency_of_water_friendly_aerosol_number_concentration_due_to_model_physics                   | number concentration of water-friendly aerosols tendency due to model physics       | kg-1 s-1      |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%dqdt(:,:,GFS_Control%ntia)    | tendency_of_ice_friendly_aerosol_number_concentration_due_to_model_physics                     | number concentration of ice-friendly aerosols tendency due to model physics         | kg-1 s-1      |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%dqdt(:,:,GFS_Control%ntrw)    | tendency_of_rain_water_mixing_ratio_due_to_model_physics                                       | moist (dry+vapor, no condensates) mixing ratio of rain water tendency due to model physics | kg kg-1 s-1   |    2 | real | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%dqdt(:,:,GFS_Control%ntsw)    | tendency_of_snow_water_mixing_ratio_due_to_model_physics                                       | moist (dry+vapor, no condensates) mixing ratio of snow water tendency due to model physics | kg kg-1 s-1   |    2 | real | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%dqdt(:,:,GFS_Control%ntgl)    | tendency_of_graupel_mixing_ratio_due_to_model_physics                                          | moist (dry+vapor, no condensates) mixing ratio of graupel tendency due to model physics    | kg kg-1 s-1   |    2 | real | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%dqdt(:,:,GFS_Control%ntke)    | tendency_of_turbulent_kinetic_energy_due_to_model_physics                                      | turbulent kinetic energy tendency due to model physics                              | J s-1         |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%dqsfc1                        | instantaneous_surface_upward_latent_heat_flux                                                  | surface upward latent heat flux                                                     | W m-2         |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%drain                         | subsurface_runoff_flux                                                                         | subsurface runoff flux                                                              | g m-2 s-1     |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%dtdt                          | tendency_of_air_temperature_due_to_model_physics                                               | air temperature tendency due to model physics                                       | K s-1         |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%dtdtc                         | tendency_of_air_temperature_due_to_radiative_heating_assuming_clear_sky                        | clear sky radiative (shortwave + longwave) heating rate at current time             | K s-1         |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%dtsfc1                        | instantaneous_surface_upward_sensible_heat_flux                                                | surface upward sensible heat flux                                                   | W m-2         |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%dtzm                          | mean_change_over_depth_in_sea_water_temperature                                                | mean of dT(z)  (zsea1 to zsea2)                                                     | K             |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%dt_mf                         | instantaneous_atmosphere_detrainment_convective_mass_flux                                      | (detrainment mass flux) * delt                                                      | kg m-2        |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%dudt                          | tendency_of_x_wind_due_to_model_physics                                                        | zonal wind tendency due to model physics                                            | m s-2         |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%dusfcg                        | instantaneous_x_stress_due_to_gravity_wave_drag                                                | zonal surface stress due to orographic gravity wave drag                            | Pa            |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%dusfc1                        | instantaneous_surface_x_momentum_flux                                                          | x momentum flux                                                                     | Pa            |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%dvdftra                       | tendency_of_vertically_diffused_tracer_concentration                                           | updated tendency of the tracers due to vertical diffusion in PBL scheme             | kg kg-1 s-1   |    3 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%dvdt                          | tendency_of_y_wind_due_to_model_physics                                                        | meridional wind tendency due to model physics                                       | m s-2         |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%dvsfcg                        | instantaneous_y_stress_due_to_gravity_wave_drag                                                | meridional surface stress due to orographic gravity wave drag                       | Pa            |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%dvsfc1                        | instantaneous_surface_y_momentum_flux                                                          | y momentum flux                                                                     | Pa            |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%dzlyr                         | layer_thickness_for_radiation                                                                  | layer thickness on radiation levels                                                 | km            |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%elvmax                        | maximum_subgrid_orography                                                                      | maximum of subgrid orography                                                        | m             |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%ep1d                          | surface_upward_potential_latent_heat_flux                                                      | surface upward potential latent heat flux                                           | W m-2         |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%evap                          | kinematic_surface_upward_latent_heat_flux                                                      | kinematic surface upward latent heat flux                                           | kg kg-1 m s-1 |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%evbs                          | soil_upward_latent_heat_flux                                                                   | soil upward latent heat flux                                                        | W m-2         |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%evcw                          | canopy_upward_latent_heat_flux                                                                 | canopy upward latent heat flux                                                      | W m-2         |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%faerlw                        | aerosol_optical_properties_for_longwave_bands_01-16                                            | aerosol optical properties for longwave bands 01-16                                 | various       |    4 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%faerlw(:,:,:,1)               | aerosol_optical_depth_for_longwave_bands_01-16                                                 | aerosol optical depth for longwave bands 01-16                                      | none          |    3 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%faerlw(:,:,:,2)               | aerosol_single_scattering_albedo_for_longwave_bands_01-16                                      | aerosol single scattering albedo for longwave bands 01-16                           | frac          |    3 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%faerlw(:,:,:,3)               | aerosol_asymmetry_parameter_for_longwave_bands_01-16                                           | aerosol asymmetry parameter for longwave bands 01-16                                | none          |    3 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%faersw                        | aerosol_optical_properties_for_shortwave_bands_01-16                                           | aerosol optical properties for shortwave bands 01-16                                | various       |    4 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%faersw(:,:,:,1)               | aerosol_optical_depth_for_shortwave_bands_01-16                                                | aerosol optical depth for shortwave bands 01-16                                     | none          |    3 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%faersw(:,:,:,2)               | aerosol_single_scattering_albedo_for_shortwave_bands_01-16                                     | aerosol single scattering albedo for shortwave bands 01-16                          | frac          |    3 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%faersw(:,:,:,3)               | aerosol_asymmetry_parameter_for_shortwave_bands_01-16                                          | aerosol asymmetry parameter for shortwave bands 01-16                               | none          |    3 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%fh2                           | Monin-Obukhov_similarity_function_for_heat_at_2m                                               | Monin-Obukhov similarity parameter for heat at 2m                                   | none          |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%flag_cice                     | flag_for_cice                                                                                  | flag for cice                                                                       | flag          |    1 | logical     |           | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%flag_guess                    | flag_for_guess_run                                                                             | flag for guess run                                                                  | flag          |    1 | logical     |           | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%flag_iter                     | flag_for_iteration                                                                             | flag for iteration                                                                  | flag          |    1 | logical     |           | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%fm10                          | Monin-Obukhov_similarity_function_for_momentum_at_10m                                          | Monin-Obukhov similarity parameter for momentum at 10m                              | none          |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%frain                         | dynamics_to_physics_timestep_ratio                                                             | ratio of dynamics timestep to physics timestep                                      | none          |    0 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%frland                        | land_area_fraction                                                                             | land area fraction                                                                  | frac          |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%fscav                         | fraction_of_tracer_scavenged                                                                   | fraction of the tracer (aerosols) that is scavenged by convection                   | km-1          |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%fswtr                         | fraction_of_cloud_top_water_scavenged                                                          | fraction of the tracer (cloud top water) that is scavenged by convection            | km-1          |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%gabsbdlw                      | surface_downwelling_longwave_flux_absorbed_by_ground                                           | total sky surface downward longwave flux absorbed by the ground                     | W m-2         |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%gamma                         | anisotropy_of_subgrid_orography                                                                | anisotropy of subgrid orography                                                     | none          |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%gamq                          | countergradient_mixing_term_for_water_vapor                                                    | countergradient mixing term for water vapor                                         | kg kg-1       |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%gamt                          | countergradient_mixing_term_for_temperature                                                    | countergradient mixing term for temperature                                         | K             |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%gasvmr                        |                                                                                                | gas volume mixing ratios                                                            | kg kg-1       |    3 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%gasvmr(:,:,1)                 | volume_mixing_ratio_co2                                                                        | volume mixing ratio co2                                                             | kg kg-1       |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%gasvmr(:,:,2)                 | volume_mixing_ratio_n2o                                                                        | volume mixing ratio no2                                                             | kg kg-1       |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%gasvmr(:,:,3)                 | volume_mixing_ratio_ch4                                                                        | volume mixing ratio ch4                                                             | kg kg-1       |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%gasvmr(:,:,4)                 | volume_mixing_ratio_o2                                                                         | volume mixing ratio o2                                                              | kg kg-1       |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%gasvmr(:,:,5)                 | volume_mixing_ratio_co                                                                         | volume mixing ratio co                                                              | kg kg-1       |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%gasvmr(:,:,6)                 | volume_mixing_ratio_cfc11                                                                      | volume mixing ratio cfc11                                                           | kg kg-1       |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%gasvmr(:,:,7)                 | volume_mixing_ratio_cfc12                                                                      | volume mixing ratio cfc12                                                           | kg kg-1       |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%gasvmr(:,:,8)                 | volume_mixing_ratio_cfc22                                                                      | volume mixing ratio cfc22                                                           | kg kg-1       |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%gasvmr(:,:,9)                 | volume_mixing_ratio_ccl4                                                                       | volume mixing ratio ccl4                                                            | kg kg-1       |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%gasvmr(:,:,10)                | volume_mixing_ratio_cfc113                                                                     | volume mixing ratio cfc113                                                          | kg kg-1       |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%gflx                          | upward_heat_flux_in_soil                                                                       | soil heat flux                                                                      | W m-2         |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%graupelmp                     | lwe_thickness_of_graupel_amount                                                                | explicit graupel fall on physics timestep                                           | m             |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%gwdcu                         | tendency_of_x_wind_due_to_convective_gravity_wave_drag                                         | zonal wind tendency due to convective gravity wave drag                             | m s-2         |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%gwdcv                         | tendency_of_y_wind_due_to_convective_gravity_wave_drag                                         | meridional wind tendency due to convective gravity wave drag                        | m s-2         |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%h2o_coeff                     | number_of_coefficients_in_h2o_forcing_data                                                     | number of coefficients in h2o forcing data                                          | index         |    0 | integer     |           | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%h2o_pres                      | natural_log_of_h2o_forcing_data_pressure_levels                                                | natural log of h2o forcing data pressure levels                                     | log(Pa)       |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%hflx                          | kinematic_surface_upward_sensible_heat_flux                                                    | kinematic surface upward sensible heat flux                                         | K m s-1       |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%hprime1                       | standard_deviation_of_subgrid_orography                                                        | standard deviation of subgrid orography                                             | m             |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%icemp                         | lwe_thickness_of_ice_amount                                                                    | explicit ice fall on physics timestep                                               | m             |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%idxday                        | daytime_points                                                                                 | daytime points                                                                      | index         |    1 | integer     |           | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%im                            | horizontal_loop_extent_OLD_NOW_MODEL_PERCENT_BLKSZ_OF_NB                                       | horizontal loop extent                                                              | count         |    0 | integer     |           | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%ipr                           | horizontal_index_of_printed_column                                                             | horizontal index of printed column                                                  | index         |    0 | integer     |           | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%islmsk                        | sea_land_ice_mask                                                                              | sea/land/ice mask (=0/1/2)                                                          | flag          |    1 | integer     |           | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%ix                            | horizontal_dimension_OLD_NOW_MODEL_PERCENT_BLKSZ_OF_NB                                         | horizontal dimension                                                                | count         |    0 | integer     |           | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%kb                            | vertical_index_difference_between_layer_and_lower_bound                                        | vertical index difference between layer and lower bound                             | index         |    0 | integer     |           | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%kbot                          | vertical_index_at_cloud_base                                                                   | vertical index at cloud base                                                        | index         |    1 | integer     |           | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%kcnv                          | flag_deep_convection                                                                           | flag indicating whether convection occurs in column (0 or 1)                        | flag          |    1 | integer     |           | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%kd                            | vertical_index_difference_between_inout_and_local                                              | vertical index difference between in/out and local                                  | index         |    0 | integer     |           | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%kinver                        | index_of_highest_temperature_inversion                                                         | index of highest temperature inversion                                              | index         |    1 | integer     |           | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%kpbl                          | vertical_index_at_top_of_atmosphere_boundary_layer                                             | vertical index at top atmospheric boundary layer                                    | index         |    1 | integer     |           | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%kt                            | vertical_index_difference_between_layer_and_upper_bound                                        | vertical index difference between layer and upper bound                             | index         |    0 | integer     |           | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%ktop                          | vertical_index_at_cloud_top                                                                    | vertical index at cloud top                                                         | index         |    1 | integer     |           | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%latidxprnt                    | latitude_index_in_debug_printouts                                                              | latitude index in debug printouts                                                   | index         |    0 | integer     |           | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%levi                          | vertical_interface_dimension                                                                   | vertical interface dimension                                                        | count         |    0 | integer     |           | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%levh2o                        | vertical_dimension_of_h2o_forcing_data                                                         | number of vertical layers in h2o forcing data                                       | count         |    0 | integer     |           | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%levozp                        | vertical_dimension_of_ozone_forcing_data                                                       | number of vertical layers in ozone forcing data                                     | count         |    0 | integer     |           | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%lm                            | vertical_layer_dimension_for_radiation                                                         | number of vertical layers for radiation                                             | count         |    0 | integer     |           | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%lmk                           | adjusted_vertical_layer_dimension_for_radiation                                                | adjusted number of vertical layers for radiation                                    | count         |    0 | integer     |           | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%lmp                           | adjusted_vertical_level_dimension_for_radiation                                                | adjusted number of vertical levels for radiation                                    | count         |    0 | integer     |           | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%mbota                         | model_layer_number_at_cloud_base                                                               | vertical indices for low, middle and high cloud bases                               | index         |    2 | integer     |           | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%mg3_as_mg2                    | flag_mg3_as_mg2                                                                                | flag for controlling prep for Morrison-Gettelman microphysics                       | flag          |    0 | logical     |           | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%mtopa                         | model_layer_number_at_cloud_top                                                                | vertical indices for low, middle and high cloud tops                                | index         |    2 | integer     |           | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%ncstrac                       | number_of_tracers_for_CS                                                                       | number of convectively transported tracers in Chikira-Sugiyama deep conv. scheme    | count         |    0 | integer     |           | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%nday                          | daytime_points_dimension                                                                       | daytime points dimension                                                            | count         |    0 | integer     |           | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%nn                            | number_of_tracers_for_convective_transport                                                     | number of tracers for convective transport                                          | count         |    0 | integer     |           | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%nncl                          | number_of_tracers_for_cloud_condensate                                                         | number of tracers for cloud condensate                                              | count         |    0 | integer     |           | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%nsamftrac                     | number_of_tracers_for_samf                                                                     | number of tracers for scale-aware mass flux schemes                                 | count         |    0 | integer     |           | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%ncgl                          | local_graupel_number_concentration                                                             | number concentration of graupel local to physics                                    | kg-1          |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%ncpi                          | local_ice_number_concentration                                                                 | number concentration of ice local to physics                                        | kg-1          |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%ncpl                          | local_condesed_water_number_concentration                                                      | number concentration of condensed water local to physics                            | kg-1          |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%ncpr                          | local_rain_number_concentration                                                                | number concentration of rain local to physics                                       | kg-1          |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%ncps                          | local_snow_number_concentration                                                                | number concentration of snow local to physics                                       | kg-1          |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%nsteps_per_reset              | number_of_time_steps_per_maximum_hourly_time_interval                                          | number_of_time_steps_per_maximum_hourly_time_interval                               | count         |    0 | integer     |           | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%ntiwx                         | index_for_ice_cloud_condensate_vertical_diffusion_tracer                                       | index for ice cloud condensate n the vertically diffused tracer array               | index         |    0 | integer     |           | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%ntk                           | index_for_turbulent_kinetic_energy_convective_transport_tracer                                 | index for turbulent kinetic energy in the convectively transported tracer array     | index         |    0 | integer     |           | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%ntkev                         | index_for_turbulent_kinetic_energy_vertical_diffusion_tracer                                   | index for turbulent kinetic energy in the vertically diffused tracer array          | index         |    0 | integer     |           | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%nvdiff                        | number_of_vertical_diffusion_tracers                                                           | number of tracers to diffuse vertically                                             | count         |    0 | integer     |           | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%oa4                           | asymmetry_of_subgrid_orography                                                                 | asymmetry of subgrid orography                                                      | none          |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%oc                            | convexity_of_subgrid_orography                                                                 | convexity of subgrid orography                                                      | none          |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%olyr                          | ozone_concentration_at_layer_for_radiation                                                     | ozone concentration layer                                                           | kg kg-1       |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%otspt                         | flag_convective_tracer_transport                                                               | flag to enable tracer transport by updrafts/downdrafts[(:,1)] or subsidence [(:,2)] | flag          |    2 | logical     |           | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%oz_coeff                      | number_of_coefficients_in_ozone_forcing_data                                                   | number of coefficients in ozone forcing data                                        | index         |    0 | integer     |           | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%oz_pres                       | natural_log_of_ozone_forcing_data_pressure_levels                                              | natural log of ozone forcing data pressure levels                                   | log(Pa)       |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%phys_hydrostatic              | flag_for_hydrostatic_heating_from_physics                                                      | flag for use of hydrostatic heating in physics                                      | flag          |    0 | logical     |           | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%plvl                          | air_pressure_at_interface_for_radiation_in_hPa                                                 | air pressure at vertical interface for radiation calculation                        | hPa           |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%plyr                          | air_pressure_at_layer_for_radiation_in_hPa                                                     | air pressure at vertical layer for radiation calculation                            | hPa           |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%prnum                         | prandtl_number                                                                                 | turbulent Prandtl number                                                            | none          |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%qgl                           | local_graupel_mixing_ratio                                                                     | moist (dry+vapor, no condensates) mixing ratio of graupel local to physics          | kg kg-1       |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%qicn                          | mass_fraction_of_convective_cloud_ice                                                          | mass fraction of convective cloud ice water                                         | kg kg-1       |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%qlcn                          | mass_fraction_of_convective_cloud_liquid_water                                                 | mass fraction of convective cloud liquid water                                      | kg kg-1       |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%qlyr                          | water_vapor_specific_humidity_at_layer_for_radiation                                           | specific humidity layer                                                             | kg kg-1       |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%qrn                           | local_rain_water_mixing_ratio                                                                  | moist (dry+vapor, no condensates) mixing ratio of rain water local to physics       | kg kg-1       |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%qsnw                          | local_snow_water_mixing_ratio                                                                  | moist (dry+vapor, no condensates) mixing ratio of snow water local to physics       | kg kg-1       |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%prcpmp                        | lwe_thickness_of_explicit_precipitation_amount                                                 | explicit precipitation (rain, ice, snow, graupel, ...) on physics timestep          | m             |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%qss                           | surface_specific_humidity                                                                      | surface air saturation specific humidity                                            | kg kg-1       |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%raddt                         | time_step_for_radiation                                                                        | radiation time step                                                                 | s             |    0 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%raincd                        | lwe_thickness_of_deep_convective_precipitation_amount                                          | deep convective rainfall amount on physics timestep                                 | m             |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%raincs                        | lwe_thickness_of_shallow_convective_precipitation_amount                                       | shallow convective rainfall amount on physics timestep                              | m             |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%rainmcadj                     | lwe_thickness_of_moist_convective_adj_precipitation_amount                                     | adjusted moist convective rainfall amount on physics timestep                       | m             |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%rainmp                        | lwe_thickness_of_explicit_rain_amount                                                          | explicit rain on physics timestep                                                   | m             |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%rainp                         | tendency_of_rain_water_mixing_ratio_due_to_microphysics                                        | tendency of rain water mixing ratio due to microphysics                             | kg kg-1 s-1   |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%rb                            | bulk_richardson_number_at_lowest_model_level                                                   | bulk Richardson number at the surface                                               | none          |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%rhc                           | critical_relative_humidity                                                                     | critical relative humidity                                                          | frac          |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%rhcbot                        | critical_relative_humidity_at_surface                                                          | critical relative humidity at the surface                                           | frac          |    0 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%rhcpbl                        | critical_relative_humidity_at_PBL_top                                                          | critical relative humidity at the PBL top                                           | frac          |    0 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%rhctop                        | critical_relative_humidity_at_top_of_atmosphere                                                | critical relative humidity at the top of atmosphere                                 | frac          |    0 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%runoff                        | surface_runoff_flux                                                                            | surface runoff flux                                                                 | g m-2 s-1     |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%save_q(:,:,GFS_Control%ntcw)  | cloud_condensed_water_mixing_ratio_save                                | moist (dry+vapor, no condensates) mixing ratio of cloud water (condensate) before entering a physics scheme | kg kg-1       |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%save_q(:,:,GFS_Control%ntiw)  | ice_water_mixing_ratio_save                                                                    | cloud ice water mixing ratio before entering a physics scheme                       | kg kg-1       |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%save_q(:,:,1)                 | water_vapor_specific_humidity_save                                                             | water vapor specific humidity before entering a physics scheme                      | kg kg-1       |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%save_q                        | tracer_concentration_save                                                                      | tracer concentration before entering a physics scheme                               | kg kg-1       |    3 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%save_t                        | air_temperature_save                                                                           | air temperature before entering a physics scheme                                    | K             |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%save_u                        | x_wind_save                                                                                    | x-wind before entering a physics scheme                                             | m s-1         |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%save_v                        | y_wind_save                                                                                    | y-wind before entering a physics scheme                                             | m s-1         |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%sbsno                         | snow_deposition_sublimation_upward_latent_heat_flux                                            | latent heat flux from snow depo/subl                                                | W m-2         |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%scmpsw                        | components_of_surface_downward_shortwave_fluxes                                                | derived type for special components of surface downward shortwave fluxes            | W m-2         |    1 | cmpfsw_type |           | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%sfcalb                        |                                                                                                |                                                                                     | frac          |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%sfcalb(:,1)                   | surface_albedo_due_to_near_IR_direct                                                           | surface albedo due to near IR direct beam                                           | frac          |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%sfcalb(:,2)                   | surface_albedo_due_to_near_IR_diffused                                                         | surface albedo due to near IR diffused beam                                         | frac          |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%sfcalb(:,3)                   | surface_albedo_due_to_UV_and_VIS_direct                                                        | surface albedo due to UV+VIS direct beam                                            | frac          |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%sfcalb(:,4)                   | surface_albedo_due_to_UV_and_VIS_diffused                                                      | surface albedo due to UV+VIS diffused beam                                          | frac          |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%sigma                         | slope_of_subgrid_orography                                                                     | slope of subgrid orography                                                          | none          |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%sigmaf                        | bounded_vegetation_area_fraction                                                               | areal fractional cover of green vegetation bounded on the bottom                    | frac          |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%sigmafrac                     | convective_updraft_area_fraction                                                               | convective updraft area fraction                                                    | frac          |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%sigmatot                      | convective_updraft_area_fraction_at_model_interfaces                                           | convective updraft area fraction at model interfaces                                | frac          |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%skip_macro                    | flag_skip_macro                                                                                | flag to skip cloud macrophysics in Morrison scheme                                  | flag          |    0 | logical     |           | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%slopetype                     | surface_slope_classification                                                                   | surface slope type at each grid cell                                                | index         |    1 | integer     |           | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%snowc                         | surface_snow_area_fraction                                                                     | surface snow area fraction                                                          | frac          |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%snohf                         | snow_freezing_rain_upward_latent_heat_flux                                                     | latent heat flux due to snow and frz rain                                           | W m-2         |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%snowmp                        | lwe_thickness_of_snow_amount                                                                   | explicit snow fall on physics timestep                                              | m             |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%snowmt                        | surface_snow_melt                                                                              | snow melt during timestep                                                           | m             |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%soiltype                      | soil_type_classification                                                                       | soil type at each grid cell                                                         | index         |    1 | integer     |           | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%stress                        | surface_wind_stress                                                                            | surface wind stress                                                                 | m2 s-2        |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%theta                         | angle_from_east_of_maximum_subgrid_orographic_variations                                       | angle with_respect to east of maximum subgrid orographic variations                 | degrees       |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%tlvl                          | air_temperature_at_interface_for_radiation                                                     | air temperature at vertical interface for radiation calculation                     | K             |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%tlyr                          | air_temperature_at_layer_for_radiation                                                         | air temperature at vertical layer for radiation calculation                         | K             |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%tracers_start_index           | start_index_of_other_tracers                                                                   | beginning index of the non-water tracer species                                     | index         |    0 | integer     |           | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%tracers_total                 | number_of_total_tracers                                                                        | total number of tracers                                                             | count         |    0 | integer     |           | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%trans                         | transpiration_flux                                                                             | total plant transpiration rate                                                      | kg m-2 s-1    |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%tseal                         | surface_skin_temperature_for_nsst                                                              | ocean surface skin temperature                                                      | K             |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%tsfa                          | surface_air_temperature_for_radiation                                                          | lowest model layer air temperature for radiation                                    | K             |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%tsfg                          | surface_ground_temperature_for_radiation                                                       | surface ground temperature for radiation                                            | K             |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%tsurf                         | surface_skin_temperature_after_iteration                                                       | surface skin temperature after iteration                                            | K             |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%tracers_water                 | number_of_water_tracers                                                                        | number of water-related tracers                                                     | count         |    0 | integer     |           | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%ud_mf                         | instantaneous_atmosphere_updraft_convective_mass_flux                                          | (updraft mass flux) * delt                                                          | kg m-2        |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%ulwsfc_cice                   | surface_upwelling_longwave_flux_for_cice                                                       | surface upwelling longwave flux for cice                                            | W m-2         |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%vdftra                        | vertically_diffused_tracer_concentration                                                       | tracer concentration diffused by PBL scheme                                         | kg kg-1       |    3 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%vegf1d                        | perturbation_of_vegetation_fraction                                                            | perturbation of vegetation fraction                                                 | frac          |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%vegtype                       | vegetation_type_classification                                                                 | vegetation type at each grid cell                                                   | index         |    1 | integer     |           | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%w_upi                         | vertical_velocity_for_updraft                                                                  | vertical velocity for updraft                                                       | m s-1         |    2 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%wcbmax                        | maximum_updraft_velocity_at_cloud_base                                                         | maximum updraft velocity at cloud base                                              | m s-1         |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%wind                          | wind_speed_at_lowest_model_layer                                                               | wind speed at lowest model level                                                    | m s-1         |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%work1                         | grid_size_related_coefficient_used_in_scale-sensitive_schemes                                  | grid size related coefficient used in scale-sensitive schemes                       | none          |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%work2                         | grid_size_related_coefficient_used_in_scale-sensitive_schemes_complement                       | complement to work1                                                                 | none          |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%work3                         | ratio_of_exner_function_between_midlayer_and_interface_at_lowest_model_layer                   | Exner function ratio bt midlayer and interface at 1st layer                         | ratio         |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%xcosz                         | instantaneous_cosine_of_zenith_angle                                                           | cosine of zenith angle at current time                                              | none          |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%xlai1d                        | perturbation_of_leaf_area_index                                                                | perturbation of leaf area index                                                     | frac          |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%xmu                           | zenith_angle_temporal_adjustment_factor_for_shortwave_fluxes                                   | zenith angle temporal adjustment factor for shortwave                               | none          |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%z01d                          | perturbation_of_momentum_roughness_length                                                      | perturbation of momentum roughness length                                           | frac          |    1 | real        | kind_phys | none   | F        |
!! | GFS_Interstitial(cdata%thrd_no)%zt1d                          | perturbation_of_heat_to_momentum_roughness_length_ratio                                        | perturbation of heat to momentum roughness length ratio                             | frac          |    1 | real        | kind_phys | none   | F        |
!!
#endif
  type GFS_interstitial_type

    real (kind=kind_phys), pointer      :: adjnirbmd(:)     => null()  !<
    real (kind=kind_phys), pointer      :: adjnirbmu(:)     => null()  !<
    real (kind=kind_phys), pointer      :: adjnirdfd(:)     => null()  !<
    real (kind=kind_phys), pointer      :: adjnirdfu(:)     => null()  !<
    real (kind=kind_phys), pointer      :: adjvisbmd(:)     => null()  !<
    real (kind=kind_phys), pointer      :: adjvisbmu(:)     => null()  !<
    real (kind=kind_phys), pointer      :: adjvisdfu(:)     => null()  !<
    real (kind=kind_phys), pointer      :: adjvisdfd(:)     => null()  !<
    real (kind=kind_phys), pointer      :: aerodp(:,:)      => null()  !<
    real (kind=kind_phys), pointer      :: alb1d(:)         => null()  !<
    real (kind=kind_phys), pointer      :: bexp1d(:)        => null()  !<
    real (kind=kind_phys), pointer      :: cd(:)            => null()  !<
    real (kind=kind_phys), pointer      :: cdq(:)           => null()  !<
    real (kind=kind_phys), pointer      :: cf_upi(:,:)      => null()  !<
    real (kind=kind_phys), pointer      :: clcn(:,:)        => null()  !<
    real (kind=kind_phys), pointer      :: cldf(:)          => null()  !<
    real (kind=kind_phys), pointer      :: cldsa(:,:)       => null()  !<
    real (kind=kind_phys), pointer      :: cldtaulw(:,:)    => null()  !<
    real (kind=kind_phys), pointer      :: cldtausw(:,:)    => null()  !<
    real (kind=kind_phys), pointer      :: cld1d(:)         => null()  !<
    real (kind=kind_phys), pointer      :: clouds(:,:,:)    => null()  !<
    real (kind=kind_phys), pointer      :: clw(:,:,:)       => null()  !<
    real (kind=kind_phys), pointer      :: clw_surf(:)      => null()  !<
    real (kind=kind_phys), pointer      :: clx(:,:)         => null()  !<
    real (kind=kind_phys), pointer      :: cndm_surf(:)     => null()  !<
    real (kind=kind_phys), pointer      :: cnv_dqldt(:,:)   => null()  !<
    real (kind=kind_phys), pointer      :: cnv_fice(:,:)    => null()  !<
    real (kind=kind_phys), pointer      :: cnv_mfd(:,:)     => null()  !<
    real (kind=kind_phys), pointer      :: cnv_ndrop(:,:)   => null()  !<
    real (kind=kind_phys), pointer      :: cnv_nice(:,:)    => null()  !<
    real (kind=kind_phys), pointer      :: cnvc(:,:)        => null()  !<
    real (kind=kind_phys), pointer      :: cnvw(:,:)        => null()  !<
    real (kind=kind_phys), pointer      :: ctei_r(:)        => null()  !<
    real (kind=kind_phys), pointer      :: ctei_rml(:)      => null()  !<
    real (kind=kind_phys), pointer      :: cumabs(:)        => null()  !<
    real (kind=kind_phys), pointer      :: dd_mf(:,:)       => null()  !<
    real (kind=kind_phys), pointer      :: de_lgth(:)       => null()  !<
    real (kind=kind_phys), pointer      :: del(:,:)         => null()  !<
    real (kind=kind_phys), pointer      :: del_gz(:,:)      => null()  !<
    real (kind=kind_phys), pointer      :: delr(:,:)        => null()  !<
    real (kind=kind_phys), pointer      :: dkt(:,:)         => null()  !<
    real (kind=kind_phys), pointer      :: dlength(:)       => null()  !<
    real (kind=kind_phys), pointer      :: dqdt(:,:,:)      => null()  !<
    real (kind=kind_phys), pointer      :: dqsfc1(:)        => null()  !<
    real (kind=kind_phys), pointer      :: drain(:)         => null()  !<
    real (kind=kind_phys), pointer      :: dtdt(:,:)        => null()  !<
    real (kind=kind_phys), pointer      :: dtdtc(:,:)       => null()  !<
    real (kind=kind_phys), pointer      :: dtsfc1(:)        => null()  !<
    real (kind=kind_phys), pointer      :: dtzm(:)          => null()  !<
    real (kind=kind_phys), pointer      :: dt_mf(:,:)       => null()  !<
    real (kind=kind_phys), pointer      :: dudt(:,:)        => null()  !<
    real (kind=kind_phys), pointer      :: dusfcg(:)        => null()  !<
    real (kind=kind_phys), pointer      :: dusfc1(:)        => null()  !<
    real (kind=kind_phys), pointer      :: dvdftra(:,:,:)   => null()  !<
    real (kind=kind_phys), pointer      :: dvdt(:,:)        => null()  !<
    real (kind=kind_phys), pointer      :: dvsfcg(:)        => null()  !<
    real (kind=kind_phys), pointer      :: dvsfc1(:)        => null()  !<
    real (kind=kind_phys), pointer      :: dzlyr(:,:)       => null()  !<
    real (kind=kind_phys), pointer      :: elvmax(:)        => null()  !<
    real (kind=kind_phys), pointer      :: ep1d(:)          => null()  !<
    real (kind=kind_phys), pointer      :: evap(:)          => null()  !<
    real (kind=kind_phys), pointer      :: evbs(:)          => null()  !<
    real (kind=kind_phys), pointer      :: evcw(:)          => null()  !<
    real (kind=kind_phys), pointer      :: faerlw(:,:,:,:)  => null()  !<
    real (kind=kind_phys), pointer      :: faersw(:,:,:,:)  => null()  !<
    real (kind=kind_phys), pointer      :: fh2(:)           => null()  !<
    logical,               pointer      :: flag_cice(:)     => null()  !<
    logical,               pointer      :: flag_guess(:)    => null()  !<
    logical,               pointer      :: flag_iter(:)     => null()  !<
    real (kind=kind_phys), pointer      :: fm10(:)          => null()  !<
    real (kind=kind_phys)               :: frain                       !<
    real (kind=kind_phys), pointer      :: frland(:)        => null()  !<
    real (kind=kind_phys), pointer      :: fscav(:)         => null()  !<
    real (kind=kind_phys), pointer      :: fswtr(:)         => null()  !<
    real (kind=kind_phys), pointer      :: gabsbdlw(:)      => null()  !<
    real (kind=kind_phys), pointer      :: gamma(:)         => null()  !<
    real (kind=kind_phys), pointer      :: gamq(:)          => null()  !<
    real (kind=kind_phys), pointer      :: gamt(:)          => null()  !<
    real (kind=kind_phys), pointer      :: gasvmr(:,:,:)    => null()  !<
    real (kind=kind_phys), pointer      :: gflx(:)          => null()  !<
    real (kind=kind_phys), pointer      :: graupelmp(:)     => null()  !<
    real (kind=kind_phys), pointer      :: gwdcu(:,:)       => null()  !<
    real (kind=kind_phys), pointer      :: gwdcv(:,:)       => null()  !<
    integer                             :: h2o_coeff                   !<
    real (kind=kind_phys), pointer      :: h2o_pres(:)      => null()  !<
    real (kind=kind_phys), pointer      :: hflx(:)          => null()  !<
    real (kind=kind_phys), pointer      :: hprime1(:)       => null()  !<
    real (kind=kind_phys), pointer      :: icemp(:)         => null()  !<
    integer,               pointer      :: idxday(:)        => null()  !<
    integer                             :: im                          !<
    integer                             :: ipr                         !<
    integer,               pointer      :: islmsk(:)        => null()  !<
    integer                             :: ix                          !<
    integer                             :: kb                          !<
    integer,               pointer      :: kbot(:)          => null()  !<
    integer,               pointer      :: kcnv(:)          => null()  !<
    integer                             :: kd                          !<
    integer,               pointer      :: kinver(:)        => null()  !<
    integer,               pointer      :: kpbl(:)          => null()  !<
    integer                             :: kt                          !<
    integer,               pointer      :: ktop(:)          => null()  !<
    integer                             :: latidxprnt                  !<
    integer                             :: levi                        !<
    integer                             :: levh2o                      !<
    integer                             :: levozp                      !<
    integer                             :: lm                          !<
    integer                             :: lmk                         !<
    integer                             :: lmp                         !<
    integer,               pointer      :: mbota(:,:)       => null()  !<
    logical                             :: mg3_as_mg2                  !<
    integer,               pointer      :: mtopa(:,:)       => null()  !<
    real (kind=kind_phys), pointer      :: ncgl(:,:)        => null()  !<
    real (kind=kind_phys), pointer      :: ncpi(:,:)        => null()  !<
    real (kind=kind_phys), pointer      :: ncpl(:,:)        => null()  !<
    real (kind=kind_phys), pointer      :: ncpr(:,:)        => null()  !<
    real (kind=kind_phys), pointer      :: ncps(:,:)        => null()  !<
    integer                             :: ncstrac                     !<
    integer                             :: nday                        !<
    integer                             :: nn                          !<
    integer                             :: nncl                        !<
    integer                             :: nsamftrac                   !<
    integer                             :: nsteps_per_reset            !<
    integer                             :: ntiwx                       !<
    integer                             :: ntk                         !<
    integer                             :: ntkev                       !<
    integer                             :: nvdiff                      !<
    real (kind=kind_phys), pointer      :: oa4(:,:)         => null()  !<
    real (kind=kind_phys), pointer      :: oc(:)            => null()  !<
    real (kind=kind_phys), pointer      :: olyr(:,:)        => null()  !<
    logical              , pointer      :: otspt(:,:)       => null()  !<
    integer                             :: oz_coeff                    !<
    real (kind=kind_phys), pointer      :: oz_pres(:)       => null()  !<
    logical                             :: phys_hydrostatic            !<
    real (kind=kind_phys), pointer      :: plvl(:,:)        => null()  !<
    real (kind=kind_phys), pointer      :: plyr(:,:)        => null()  !<
    real (kind=kind_phys), pointer      :: prcpmp(:)        => null()  !<
    real (kind=kind_phys), pointer      :: prnum(:,:)       => null()  !<
    real (kind=kind_phys), pointer      :: qgl(:,:)         => null()  !<
    real (kind=kind_phys), pointer      :: qicn(:,:)        => null()  !<
    real (kind=kind_phys), pointer      :: qlcn(:,:)        => null()  !<
    real (kind=kind_phys), pointer      :: qlyr(:,:)        => null()  !<
    real (kind=kind_phys), pointer      :: qrn(:,:)         => null()  !<
    real (kind=kind_phys), pointer      :: qsnw(:,:)        => null()  !<
    real (kind=kind_phys), pointer      :: qss(:)           => null()  !<
    real (kind=kind_phys)               :: raddt                       !<
    real (kind=kind_phys), pointer      :: rainmp(:)        => null()  !<
    real (kind=kind_phys), pointer      :: raincd(:)        => null()  !<
    real (kind=kind_phys), pointer      :: raincs(:)        => null()  !<
    real (kind=kind_phys), pointer      :: rainmcadj(:)     => null()  !<
    real (kind=kind_phys), pointer      :: rainp(:,:)       => null()  !<
    real (kind=kind_phys), pointer      :: rb(:)            => null()  !<
    real (kind=kind_phys), pointer      :: rhc(:,:)         => null()  !<
    real (kind=kind_phys)               :: rhcbot                      !<
    real (kind=kind_phys)               :: rhcpbl                      !<
    real (kind=kind_phys)               :: rhctop                      !<
    real (kind=kind_phys), pointer      :: runoff(:)        => null()  !<
    real (kind=kind_phys), pointer      :: save_q(:,:,:)    => null()  !<
    real (kind=kind_phys), pointer      :: save_t(:,:)      => null()  !<
    real (kind=kind_phys), pointer      :: save_u(:,:)      => null()  !<
    real (kind=kind_phys), pointer      :: save_v(:,:)      => null()  !<
    real (kind=kind_phys), pointer      :: sbsno(:)         => null()  !<
    type (cmpfsw_type),    pointer      :: scmpsw(:)        => null()  !<
    real (kind=kind_phys), pointer      :: sfcalb(:,:)      => null()  !<
    real (kind=kind_phys), pointer      :: sigma(:)         => null()  !<
    real (kind=kind_phys), pointer      :: sigmaf(:)        => null()  !<
    real (kind=kind_phys), pointer      :: sigmafrac(:,:)   => null()  !<
    real (kind=kind_phys), pointer      :: sigmatot(:,:)    => null()  !<
    logical                             :: skip_macro                  !<
    integer, pointer                    :: slopetype(:)     => null()  !<
    real (kind=kind_phys), pointer      :: snowc(:)         => null()  !<
    real (kind=kind_phys), pointer      :: snohf(:)         => null()  !<
    real (kind=kind_phys), pointer      :: snowmp(:)        => null()  !<
    real (kind=kind_phys), pointer      :: snowmt(:)        => null()  !<
    integer, pointer                    :: soiltype(:)      => null()  !<
    real (kind=kind_phys), pointer      :: stress(:)        => null()  !<
    real (kind=kind_phys), pointer      :: theta(:)         => null()  !<
    real (kind=kind_phys), pointer      :: tlvl(:,:)        => null()  !<
    real (kind=kind_phys), pointer      :: tlyr(:,:)        => null()  !<
    integer                             :: tracers_start_index         !<
    integer                             :: tracers_total               !<
    integer                             :: tracers_water               !<
    real (kind=kind_phys), pointer      :: trans(:)         => null()  !<
    real (kind=kind_phys), pointer      :: tseal(:)         => null()  !<
    real (kind=kind_phys), pointer      :: tsfa(:)          => null()  !<
    real (kind=kind_phys), pointer      :: tsfg(:)          => null()  !<
    real (kind=kind_phys), pointer      :: tsnow(:)         => null()  !<
    real (kind=kind_phys), pointer      :: tsurf(:)         => null()  !<
    real (kind=kind_phys), pointer      :: ud_mf(:,:)       => null()  !<
    real (kind=kind_phys), pointer      :: ulwsfc_cice(:)   => null()  !<
    real (kind=kind_phys), pointer      :: vdftra(:,:,:)    => null()  !<
    real (kind=kind_phys), pointer      :: vegf1d(:)        => null()  !<
    integer, pointer                    :: vegtype(:)       => null()  !<
    real (kind=kind_phys), pointer      :: w_upi(:,:)       => null()  !<
    real (kind=kind_phys), pointer      :: wcbmax(:)        => null()  !<
    real (kind=kind_phys), pointer      :: wind(:)          => null()  !<
    real (kind=kind_phys), pointer      :: work1(:)         => null()  !<
    real (kind=kind_phys), pointer      :: work2(:)         => null()  !<
    real (kind=kind_phys), pointer      :: work3(:)         => null()  !<
    real (kind=kind_phys), pointer      :: xcosz(:)         => null()  !<
    real (kind=kind_phys), pointer      :: xlai1d(:)        => null()  !<
    real (kind=kind_phys), pointer      :: xmu(:)           => null()  !<
    real (kind=kind_phys), pointer      :: z01d(:)          => null()  !<
    real (kind=kind_phys), pointer      :: zt1d(:)          => null()  !<
#ifdef HYBRID
    logical                             :: non_uniform_blocks          !< flag to indicate non-uniform blocks are used, only for CCPP hybrid
#endif

    contains
      procedure :: create      => interstitial_create     !<   allocate array data
      procedure :: rad_reset   => interstitial_rad_reset  !<   reset array data for radiation
      procedure :: phys_reset  => interstitial_phys_reset !<   reset array data for physics
      procedure :: mprint      => interstitial_print      !<   print array data

  end type GFS_interstitial_type
#endif

!-------------------------
! GFS sub-containers
!-------------------------
#ifdef CCPP
  type GFS_data_type
     type(GFS_statein_type)  :: Statein
     type(GFS_stateout_type) :: Stateout
     type(GFS_sfcprop_type)  :: Sfcprop
     type(GFS_coupling_type) :: Coupling
     type(GFS_grid_type)     :: Grid
     type(GFS_tbd_type)      :: Tbd
     type(GFS_cldprop_type)  :: Cldprop
     type(GFS_radtend_type)  :: Radtend
     type(GFS_diag_type)     :: Intdiag
  end type GFS_data_type
#endif

!----------------
! PUBLIC ENTITIES
!----------------
  public GFS_init_type
  public GFS_statein_type,  GFS_stateout_type, GFS_sfcprop_type, &
         GFS_coupling_type
  public GFS_control_type,  GFS_grid_type,     GFS_tbd_type, &
         GFS_cldprop_type,  GFS_radtend_type,  GFS_diag_type
#ifdef CCPP
  public GFS_interstitial_type
#endif

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
! stochastic physics SKEB variable
    allocate (Statein%diss_est(IM,Model%levs))
    Statein%diss_est= clear_val
    !--- physics only variables
    allocate (Statein%pgr    (IM))
    allocate (Statein%ugrs   (IM,Model%levs))
    allocate (Statein%vgrs   (IM,Model%levs))
    allocate (Statein%qgrs   (IM,Model%levs,Model%ntrac))

    Statein%qgrs   = clear_val
    Statein%pgr    = clear_val
    Statein%ugrs   = clear_val
    Statein%vgrs   = clear_val

    !--- soil state variables - for soil SPPT - sfc-perts, mgehne
    allocate (Statein%smc  (IM,Model%lsoil))
    allocate (Statein%stc  (IM,Model%lsoil))
    allocate (Statein%slc  (IM,Model%lsoil))

    Statein%smc   = clear_val
    Statein%stc   = clear_val
    Statein%slc   = clear_val

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
    allocate (Sfcprop%lakemsk(IM))
    allocate (Sfcprop%tsfc   (IM))
    allocate (Sfcprop%tisfc  (IM))
    allocate (Sfcprop%snowd  (IM))
    allocate (Sfcprop%zorl   (IM))
    allocate (Sfcprop%fice   (IM))
    allocate (Sfcprop%hprim  (IM))
    allocate (Sfcprop%hprime (IM,Model%nmtvr))

    Sfcprop%slmsk   = clear_val
    Sfcprop%lakemsk = clear_val
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
    allocate (Sfcprop%sr     (IM))
    allocate (Sfcprop%slc    (IM,Model%lsoil))
    allocate (Sfcprop%smc    (IM,Model%lsoil))
    allocate (Sfcprop%stc    (IM,Model%lsoil))
    allocate (Sfcprop%wet1   (IM))

    Sfcprop%hice   = clear_val
    Sfcprop%weasd  = clear_val
    Sfcprop%sncovr = clear_val
    Sfcprop%canopy = clear_val
    Sfcprop%ffmm   = clear_val
    Sfcprop%ffhh   = clear_val
    Sfcprop%f10m   = clear_val
    Sfcprop%tprcp  = clear_val
    Sfcprop%srflag = clear_val
    Sfcprop%sr     = clear_val
    Sfcprop%slc    = clear_val
    Sfcprop%smc    = clear_val
    Sfcprop%stc    = clear_val
    Sfcprop%wet1   = clear_val

!--- Out
    allocate (Sfcprop%t2m (IM))
#ifdef CCPP
    allocate (Sfcprop%th2m(IM))
#endif
    allocate (Sfcprop%q2m (IM))

    Sfcprop%t2m  = clear_val
#ifdef CCPP
    Sfcprop%th2m = clear_val
#endif
    Sfcprop%q2m  = clear_val

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

#ifdef CCPP
    if (Model%lsm == Model%lsm_ruc) then
       ! For land surface models with different numbers of levels than the four NOAH levels
       allocate (Sfcprop%sh2o        (IM,Model%lsoil_lsm))
       allocate (Sfcprop%keepsmfr    (IM,Model%lsoil_lsm))
       allocate (Sfcprop%smois       (IM,Model%lsoil_lsm))
       allocate (Sfcprop%tslb        (IM,Model%lsoil_lsm))
       allocate (Sfcprop%flag_frsoil (IM,Model%lsoil_lsm))
       allocate (Sfcprop%zs          (Model%lsoil_lsm))
       allocate (Sfcprop%clw_surf    (IM))
       allocate (Sfcprop%qwv_surf    (IM))
       allocate (Sfcprop%cndm_surf   (IM))
       allocate (Sfcprop%rhofr       (IM))
       allocate (Sfcprop%tsnow       (IM))
       !
       Sfcprop%sh2o        = clear_val
       Sfcprop%keepsmfr    = clear_val
       Sfcprop%smois       = clear_val
       Sfcprop%tslb        = clear_val
       Sfcprop%zs          = clear_val
       Sfcprop%clw_surf    = clear_val
       Sfcprop%qwv_surf    = clear_val
       Sfcprop%cndm_surf   = clear_val
       Sfcprop%flag_frsoil = clear_val
       Sfcprop%rhofr       = clear_val
       Sfcprop%tsnow       = clear_val
    end if
    if (Model%do_mynnsfclay) then
    ! For MYNN surface layer scheme
       !print*,"Allocating all MYNN-sfclay variables"
       allocate (Sfcprop%ustm   (IM ))
       allocate (Sfcprop%zol    (IM ))
       allocate (Sfcprop%mol    (IM ))
       allocate (Sfcprop%rmol   (IM ))
       allocate (Sfcprop%flhc   (IM ))
       allocate (Sfcprop%flqc   (IM ))
       allocate (Sfcprop%chs2   (IM ))
       allocate (Sfcprop%cqs2   (IM ))
       allocate (Sfcprop%lh     (IM ))
       !
       !print*,"Initializing all MYNN-SfcLay variables with ",clear_val
       Sfcprop%ustm        = clear_val
       Sfcprop%zol         = clear_val
       Sfcprop%mol         = clear_val
       Sfcprop%rmol        = clear_val
       Sfcprop%flhc        = clear_val
       Sfcprop%flqc        = clear_val
       Sfcprop%chs2        = clear_val
       Sfcprop%cqs2        = clear_val
       Sfcprop%lh          = clear_val
    end if
    if (Model%imfdeepcnv == 3) then
        allocate (Sfcprop%conv_act(IM))
        Sfcprop%conv_act = zero
    end if
#endif

  end subroutine sfcprop_create


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

    if (Model%cplflx .or. Model%do_sppt .or. Model%cplchm) then
      allocate (Coupling%rain_cpl     (IM))
      allocate (Coupling%snow_cpl     (IM))
      Coupling%rain_cpl     = clear_val
      Coupling%snow_cpl     = clear_val
    endif

    if (Model%cplflx .or. Model%cplwav) then
      !--- instantaneous quantities 
      allocate (Coupling%u10mi_cpl   (IM))
      allocate (Coupling%v10mi_cpl   (IM))

      Coupling%u10mi_cpl   = clear_val
      Coupling%v10mi_cpl   = clear_val
    endif 

    if (Model%cplflx) then
      !--- incoming quantities
      allocate (Coupling%slimskin_cpl (IM))
      allocate (Coupling%dusfcin_cpl  (IM))
      allocate (Coupling%dvsfcin_cpl  (IM))
      allocate (Coupling%dtsfcin_cpl  (IM))
      allocate (Coupling%dqsfcin_cpl  (IM))
      allocate (Coupling%ulwsfcin_cpl (IM))
      allocate (Coupling%tseain_cpl   (IM))
      allocate (Coupling%tisfcin_cpl  (IM))
      allocate (Coupling%ficein_cpl   (IM))
      allocate (Coupling%hicein_cpl   (IM))
      allocate (Coupling%hsnoin_cpl   (IM))

      Coupling%slimskin_cpl = clear_val
      Coupling%dusfcin_cpl  = clear_val
      Coupling%dvsfcin_cpl  = clear_val
      Coupling%dtsfcin_cpl  = clear_val
      Coupling%dqsfcin_cpl  = clear_val
      Coupling%ulwsfcin_cpl = clear_val
      Coupling%tseain_cpl   = clear_val
      Coupling%tisfcin_cpl  = clear_val
      Coupling%ficein_cpl   = clear_val
      Coupling%hicein_cpl   = clear_val
      Coupling%hsnoin_cpl   = clear_val

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
      Coupling%tsfci_cpl   = clear_val
      Coupling%psurfi_cpl  = clear_val
      Coupling%oro_cpl     = clear_val  !< pointer to sfcprop%oro
      Coupling%slmsk_cpl   = clear_val  !< pointer to sfcprop%slmsk
    endif

   !-- cellular automata
    allocate (Coupling%tconvtend (IM,Model%levs))
    allocate (Coupling%qconvtend (IM,Model%levs))
    allocate (Coupling%uconvtend (IM,Model%levs))
    allocate (Coupling%vconvtend (IM,Model%levs))
    allocate (Coupling%cape (IM))
    allocate (Coupling%ca_out (IM))
    allocate (Coupling%ca_deep (IM))
    allocate (Coupling%ca_turb (IM))
    allocate (Coupling%ca_shal (IM))
    allocate (Coupling%ca_rad (IM))
    allocate (Coupling%ca_micro (IM))
    Coupling%ca_out = clear_val
    Coupling%ca_deep = clear_val
    Coupling%ca_turb = clear_val
    Coupling%ca_shal = clear_val
    Coupling%ca_rad =clear_val
    Coupling%ca_micro = clear_val   
    Coupling%cape = clear_val
    Coupling%tconvtend = clear_val
    Coupling%qconvtend = clear_val
    Coupling%uconvtend = clear_val
    Coupling%vconvtend = clear_val

    ! -- GSDCHEM coupling options
    if (Model%cplchm) then
      !--- outgoing instantaneous quantities
      allocate (Coupling%ushfsfci  (IM))
      allocate (Coupling%dkt       (IM,Model%levs))
      !--- accumulated convective rainfall
      allocate (Coupling%rainc_cpl (IM))

      Coupling%rainc_cpl = clear_val
      Coupling%ushfsfci  = clear_val
      Coupling%dkt       = clear_val
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

    !--- stochastic physics option
    if (Model%do_sfcperts) then
      allocate (Coupling%sfc_wts  (IM,Model%nsfcpert))
      Coupling%sfc_wts = clear_val
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
    endif

    !--- needed for Thompson's aerosol option
    if(Model%imp_physics == Model%imp_physics_thompson) then
      allocate (Coupling%nwfa2d (IM))
      allocate (Coupling%nifa2d (IM))
      Coupling%nwfa2d   = clear_val
      Coupling%nifa2d   = clear_val
    endif

  end subroutine coupling_create


!----------------------
! GFS_control_type%init
!----------------------
  subroutine control_initialize (Model, nlunit, fn_nml, me, master, &
                                 logunit, isc, jsc, nx, ny, levs,   &
                                 cnx, cny, gnx, gny, dt_dycore,     &
                                 dt_phys, idat, jdat, tracer_names, &
                                 input_nml_file, tile_num           &
#ifdef CCPP
                                ,ak, bk, blksz,                     &
                                 restart, hydrostatic,              &
                                 communicator, ntasks, nthreads     &
#endif
                                 )

!--- modules
#ifdef CCPP
    use physcons,         only: con_rerth, con_pi
#else
    use physcons,         only: dxmax, dxmin, dxinv, con_rerth, con_pi, rhc_max
#endif
    use mersenne_twister, only: random_setseed, random_number
#if !defined(CCPP) || defined(HYBRID)
    use module_ras,       only: nrcmax
#endif
    use parse_tracers,    only: get_tracer_index
#if !defined(CCPP) || defined(HYBRID)
    use wam_f107_kp_mod,  only: f107_kp_size, f107_kp_interval,     &
                                f107_kp_skip_size, f107_kp_data_size
#endif
    implicit none

!--- interface variables
    class(GFS_control_type)            :: Model
    integer,                intent(in) :: nlunit
    character(len=64),      intent(in) :: fn_nml
    integer,                intent(in) :: me
    integer,                intent(in) :: master
    integer,                intent(in) :: logunit
    integer,                intent(in) :: tile_num
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
    character(len=256),     intent(in), pointer :: input_nml_file(:)
#ifdef CCPP
    real(kind=kind_phys), dimension(:), intent(in) :: ak
    real(kind=kind_phys), dimension(:), intent(in) :: bk
    integer,                intent(in) :: blksz(:)
    logical,                intent(in) :: restart
    logical,                intent(in) :: hydrostatic
    integer,                intent(in) :: communicator
    integer,                intent(in) :: ntasks
    integer,                intent(in) :: nthreads
#endif
    !--- local variables
    integer :: n
    integer :: ios
    integer :: seed0
    logical :: exists
    real(kind=kind_phys) :: tem
    real(kind=kind_phys) :: rinc(5)
    real(kind=kind_phys) :: wrk(1)
    real(kind=kind_phys), parameter :: con_hr = 3600.
#ifdef CCPP
    real(kind=kind_phys), parameter   :: p_ref = 101325.0d0
#endif

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
    logical              :: cplchm         = .false.         !< default no cplchm collection

!--- integrated dynamics through earth's atmosphere
    logical              :: lsidea         = .false.

!--- radiation parameters
    real(kind=kind_phys) :: fhswr          = 3600.           !< frequency for shortwave radiation (secs)
    real(kind=kind_phys) :: fhlwr          = 3600.           !< frequency for longwave radiation (secs)
    integer              :: levr           = -99             !< number of vertical levels for radiation calculations
    integer              :: nfxr           = 39+6            !< second dimension of input/output array fluxr   
    logical              :: aero_in        = .false.         !< flag for initializing aero data 
    logical              :: iccn           = .false.         !< logical to use IN CCN forcing for MG2/3
    integer              :: iflip          =  1              !< iflip - is not the same as flipv
    integer              :: isol           =  0              !< solar constant scheme control flag
                                                             !!\brief 0: fixed value = 1366.0 \f$W m^{-2}\f$ (old standard)
                                                             !!\n 10: fixed value = 1360.8 \f$W m^{-2}\f$ (new standard)
                                                             !!\n 1: NOAA ABS-scale TSI table (yearly) with 11-yr cycle approximation
                                                             !!\n 2: NOAA TIM-scale TSI table (yearly) with 11-yr cycle approximation
                                                             !!\n 3: CMIP5 TIM-scale TSI table (yearly) with 11-yr cycle approximation
                                                             !!\n 4: CMIP5 TIM-scale TSI table (monthly) with 11-yr cycle approximation
    integer              :: ico2           =  0              !< prescribed global mean value (old opernl)
    integer              :: ialb           =  0              !< use climatology alb, based on sfc type
                                                             !< 1 => use modis based alb
    integer              :: iems           =  0              !< use fixed value of 1.0
    integer              :: iaer           =  1              !< default aerosol effect in sw only
    integer              :: icliq_sw       =  1              !< sw optical property for liquid clouds
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
    logical              :: lwhtr          = .true.          !< flag to output LW heating rate (Radtend%lwhc)
    logical              :: swhtr          = .true.          !< flag to output SW heating rate (Radtend%swhc)

!--- Z-C microphysical parameters
    integer              :: ncld           =  1                 !< choice of cloud scheme
    integer              :: imp_physics    =  99                !< choice of cloud scheme
    real(kind=kind_phys) :: psautco(2)     = (/6.0d-4,3.0d-4/)  !< [in] auto conversion coeff from ice to snow
    real(kind=kind_phys) :: prautco(2)     = (/1.0d-4,1.0d-4/)  !< [in] auto conversion coeff from cloud to rain
    real(kind=kind_phys) :: evpco          = 2.0d-5             !< [in] coeff for evaporation of largescale rain
    real(kind=kind_phys) :: wminco(2)      = (/1.0d-5,1.0d-5/)  !< [in] water and ice minimum threshold for Zhao
!---Max hourly
    real(kind=kind_phys) :: avg_max_length = 3600.              !< reset value in seconds for max hourly.

!--- M-G microphysical parameters
    integer              :: fprcp             =  0                 !< no prognostic rain and snow (MG)
    integer              :: pdfflag           =  4                 !< pdf flag for MG macro physics
    real(kind=kind_phys) :: mg_dcs            = 200.0              !< Morrison-Gettelman microphysics parameters
    real(kind=kind_phys) :: mg_qcvar          = 1.0
    real(kind=kind_phys) :: mg_ts_auto_ice(2) = (/180.0,180.0/)    !< ice auto conversion time scale
#ifdef CCPP
    real(kind=kind_phys) :: mg_rhmini         = 1.01               !< relative humidity threshold parameter for nucleating ice
#endif
    real(kind=kind_phys) :: mg_ncnst          = 100.e6             !< constant droplet num concentration (m-3)
    real(kind=kind_phys) :: mg_ninst          = 0.15e6             !< constant ice num concentration (m-3)
    real(kind=kind_phys) :: mg_ngnst          = 0.10e6             !< constant graupel/hail num concentration (m-3) = 0.1e6_r8
    real(kind=kind_phys) :: mg_alf            = 1.0                !< tuning factor for alphs in MG macrophysics
    real(kind=kind_phys) :: mg_qcmin(2)       = (/1.0d-9,1.0d-9/)  !< min liquid and ice mixing ratio in Mg macro clouds
    real(kind=kind_phys) :: mg_berg_eff_factor = 2.0               !< berg efficiency factor
    character(len=16)    :: mg_precip_frac_method = 'max_overlap'  !< type of precipitation fraction method
#ifdef CCPP
    real(kind=kind_phys) :: tf              = 258.16
    real(kind=kind_phys) :: tcr             = 273.16
#endif
!
    logical              :: effr_in         = .false.           !< flag to use effective radii of cloud species in radiation
    logical              :: microp_uniform  = .true.
    logical              :: do_cldliq       = .true.
    logical              :: do_cldice       = .true.
    logical              :: hetfrz_classnuc = .false.
    logical              :: mg_nccons       = .false.           !< set .true. to specify constant cloud droplet number
    logical              :: mg_nicons       = .false.           !< set .true. to specify constant cloud ice number
    logical              :: mg_ngcons       = .false.           !< set .true. to specify constant graupel/hail number
    logical              :: sed_supersat    = .true.
    logical              :: do_sb_physics   = .true.
    logical              :: mg_do_graupel   = .true.            !< set .true. to turn on prognostic grapuel (with fprcp=2)
    logical              :: mg_do_hail      = .false.           !< set .true. to turn on prognostic hail (with fprcp=2)
    logical              :: mg_do_ice_gmao  = .false.           !< set .true. to turn on gmao ice autoconversion for MG microphysics
    logical              :: mg_do_liq_liu   = .true.            !< set .true. to turn on liu liquid treatment for MG microphysics

    !--- Thompson microphysical parameters
    logical              :: ltaerosol      = .false.            !< flag for aerosol version
    logical              :: lradar         = .false.            !< flag for radar reflectivity
    real(kind=kind_phys) :: ttendlim       = -999.0             !< temperature tendency limiter, set to <0 to deactivate

    !--- GFDL microphysical parameters
    logical              :: lgfdlmprad     = .false.            !< flag for GFDLMP radiation interaction

    !--- land/surface model parameters
    integer              :: lsm            =  1              !< flag for land surface model to use =0  for osu lsm; =1  for noah lsm; =2  for RUC lsm
    integer              :: lsoil          =  4              !< number of soil layers
#ifdef CCPP
    integer              :: lsoil_lsm      =  -1             !< number of soil layers internal to land surface model; -1 use lsoil
#endif
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
    logical              :: trans_trac     = .false.                  !< flag for convective transport of tracers (RAS, CS, or SAMF)
    logical              :: old_monin      = .false.                  !< flag for diff monin schemes
    logical              :: cnvgwd         = .false.                  !< flag for convective gravity wave drag
    logical              :: mstrat         = .false.                  !< flag for moorthi approach for stratus
    logical              :: moist_adj      = .false.                  !< flag for moist convective adjustment
    logical              :: cscnv          = .false.                  !< flag for Chikira-Sugiyama convection
    logical              :: cal_pre        = .false.                  !< flag controls precip type algorithm
    logical              :: do_aw          = .false.                  !< AW scale-aware option in cs convection
    logical              :: do_awdd        = .false.                  !< AW scale-aware option in cs convection
    logical              :: flx_form       = .false.                  !< AW scale-aware option in cs convection
    logical              :: do_shoc        = .false.                  !< flag for SHOC
    logical              :: shocaftcnv     = .false.                  !< flag for SHOC
    logical              :: shoc_cld       = .false.                  !< flag for SHOC in grrad
#ifdef CCPP
    logical              :: oz_phys        = .true.                   !< flag for old (2006) ozone physics
    logical              :: oz_phys_2015   = .false.                  !< flag for new (2015) ozone physics
#endif
    logical              :: h2o_phys       = .false.                  !< flag for stratosphere h2o
    logical              :: pdfcld         = .false.                  !< flag for pdfcld
    logical              :: shcnvcw        = .false.                  !< flag for shallow convective cloud
    logical              :: redrag         = .false.                  !< flag for reduced drag coeff. over sea
    logical              :: hybedmf        = .false.                  !< flag for hybrid EDMF PBL scheme
    logical              :: satmedmf       = .false.                  !< flag for scale-aware TKE-based moist edmf
                                                                      !< vertical turbulent mixing scheme
    logical              :: shinhong       = .false.                  !< flag for scale-aware Shinhong vertical turbulent mixing scheme
    logical              :: do_ysu         = .false.                  !< flag for YSU vertical turbulent mixing scheme
    logical              :: dspheat        = .false.                  !< flag for tke dissipative heating
    logical              :: lheatstrg      = .false.                  !< flag for canopy heat storage parameterization
    logical              :: cnvcld         = .false.
    logical              :: random_clds    = .false.                  !< flag controls whether clouds are random
    logical              :: shal_cnv       = .false.                  !< flag for calling shallow convection
    integer              :: imfshalcnv     =  1                       !< flag for mass-flux shallow convection scheme
                                                                      !<     1: July 2010 version of mass-flux shallow conv scheme
                                                                      !<         current operational version as of 2016
                                                                      !<     2: scale- & aerosol-aware mass-flux shallow conv scheme (2017)
                                                                      !<     3: scale- & aerosol-aware Grell-Freitas scheme (GSD)
                                                                      !<     4: New Tiedtke scheme (CAPS)
                                                                      !<     0: modified Tiedtke's eddy-diffusion shallow conv scheme
                                                                      !<    -1: no shallow convection used
    integer              :: imfdeepcnv     =  1                       !< flag for mass-flux deep convection scheme
                                                                      !<     1: July 2010 version of SAS conv scheme
                                                                      !<           current operational version as of 2016
                                                                      !<     2: scale- & aerosol-aware mass-flux deep conv scheme (2017)
                                                                      !<     3: scale- & aerosol-aware Grell-Freitas scheme (GSD)
                                                                      !<     4: New Tiedtke scheme (CAPS)
    logical              :: do_deep        = .true.                   !< whether to do deep convection
#ifdef CCPP
    logical              :: do_mynnedmf       = .false.               !< flag for MYNN-EDMF
    logical              :: do_mynnsfclay     = .false.               !< flag for MYNN Surface Layer Scheme
    integer              :: grav_settling     = 0
    integer              :: bl_mynn_tkebudget = 0
    logical              :: bl_mynn_tkeadvect = .false.
    integer              :: bl_mynn_cloudpdf  = 2
    integer              :: bl_mynn_mixlength = 2
    integer              :: bl_mynn_edmf      = 0
    integer              :: bl_mynn_edmf_mom  = 1
    integer              :: bl_mynn_edmf_tke  = 0
    integer              :: bl_mynn_edmf_part = 0
    integer              :: bl_mynn_cloudmix  = 1
    integer              :: bl_mynn_mixqt     = 0
    integer              :: icloud_bl         = 1
#endif
    integer              :: nmtvr          = 14                       !< number of topographic variables such as variance etc
                                                                      !< used in the GWD parameterization
    integer              :: jcap           =  1              !< number of spectral wave trancation used only by sascnv shalcnv
!   real(kind=kind_phys) :: cs_parm(10) = (/5.0,2.5,1.0e3,3.0e3,20.0,-999.,-999.,0.,0.,0./)
    real(kind=kind_phys) :: cs_parm(10) = (/8.0,4.0,1.0e3,3.5e3,20.0,1.0,-999.,1.,0.6,0./)
    real(kind=kind_phys) :: flgmin(2)      = (/0.180,0.220/)          !< [in] ice fraction bounds
    real(kind=kind_phys) :: cgwf(2)        = (/0.5d0,0.05d0/)         !< multiplication factor for convective GWD
    real(kind=kind_phys) :: ccwf(2)        = (/1.0d0,1.0d0/)          !< multiplication factor for critical cloud
                                                                      !< workfunction for RAS
    real(kind=kind_phys) :: cdmbgwd(2)     = (/2.0d0,0.25d0/)         !< multiplication factors for cdmb and gwd
    real(kind=kind_phys) :: sup            = 1.0                      !< supersaturation in pdf cloud (IMP_physics=98) when t is very low
                                                                      !< or ice super saturation in SHOC (when do_shoc=.true.)
    real(kind=kind_phys) :: ctei_rm(2)     = (/10.0d0,10.0d0/)        !< critical cloud top entrainment instability criteria
                                                                      !< (used if mstrat=.true.)
    real(kind=kind_phys) :: crtrh(3)       = (/0.90d0,0.90d0,0.90d0/) !< critical relative humidity at the surface
                                                                      !< PBL top and at the top of the atmosphere
    real(kind=kind_phys) :: dlqf(2)        = (/0.0d0,0.0d0/)          !< factor for cloud condensate detrainment
                                                                      !< from cloud edges for RAS
    real(kind=kind_phys) :: psauras(2)     = (/1.0d-3,1.0d-3/)        !< [in] auto conversion coeff from ice to snow in ras
    real(kind=kind_phys) :: prauras(2)     = (/2.0d-3,2.0d-3/)        !< [in] auto conversion coeff from cloud to rain in ras
    real(kind=kind_phys) :: wminras(2)     = (/1.0d-5,1.0d-5/)        !< [in] water and ice minimum threshold for ras

    real(kind=kind_phys) :: rbcr           = 0.25                     !< Critical Richardson Number in PBL scheme
    real(kind=kind_phys) :: shoc_parm(5)   = (/7000.0,1.0,4.2857143,0.7,-999.0/)  !< some tunable parameters for shoc

!--- Rayleigh friction
    real(kind=kind_phys) :: prslrd0        = 0.0d0           !< pressure level from which Rayleigh Damping is applied
    real(kind=kind_phys) :: ral_ts         = 0.0d0           !< time scale for Rayleigh damping in days

!--- mass flux deep convection
    real(kind=kind_phys) :: clam_deep      = 0.1             !< c_e for deep convection (Han and Pan, 2011, eq(6))
    real(kind=kind_phys) :: c0s_deep       = 0.002           !< convective rain conversion parameter
    real(kind=kind_phys) :: c1_deep        = 0.002           !< conversion parameter of detrainment from liquid water into grid-scale cloud water
    real(kind=kind_phys) :: betal_deep     = 0.05            !< fraction factor of downdraft air mass reaching ground surface over land
    real(kind=kind_phys) :: betas_deep     = 0.05            !< fraction factor of downdraft air mass reaching ground surface over sea
    real(kind=kind_phys) :: evfact_deep    = 0.3             !< evaporation factor from convective rain
    real(kind=kind_phys) :: evfactl_deep   = 0.3             !< evaporation factor from convective rain over land
    real(kind=kind_phys) :: pgcon_deep     = 0.55            !< reduction factor in momentum transport due to convection induced pressure gradient force
                                                             !< 0.7 : Gregory et al. (1997, QJRMS)
                                                             !< 0.55: Zhang & Wu (2003, JAS)
    real(kind=kind_phys) :: asolfac_deep   = 0.958           !< aerosol-aware parameter based on Lim (2011)
                                                             !< asolfac= cx / c0s(=.002)
                                                             !< cx = min([-0.7 ln(Nccn) + 24]*1.e-4, c0s)
                                                             !< Nccn: CCN number concentration in cm^(-3)
                                                             !< Until a realistic Nccn is provided, Nccns are assumed
                                                             !< as Nccn=100 for sea and Nccn=1000 for land

!--- mass flux shallow convection
    real(kind=kind_phys) :: clam_shal      = 0.3             !< c_e for shallow convection (Han and Pan, 2011, eq(6))
    real(kind=kind_phys) :: c0s_shal       = 0.002           !< conversion parameter of detrainment from liquid water into convetive precipitaiton
    real(kind=kind_phys) :: c1_shal        = 5.e-4           !< conversion parameter of detrainment from liquid water into grid-scale cloud water
    real(kind=kind_phys) :: pgcon_shal     = 0.55            !< reduction factor in momentum transport due to convection induced pressure gradient force
                                                             !< 0.7 : Gregory et al. (1997, QJRMS)
                                                             !< 0.55: Zhang & Wu (2003, JAS)
    real(kind=kind_phys) :: asolfac_shal   = 0.958           !< aerosol-aware parameter based on Lim (2011)
                                                             !< asolfac= cx / c0s(=.002)
                                                             !< cx = min([-0.7 ln(Nccn) + 24]*1.e-4, c0s)
                                                             !< Nccn: CCN number concentration in cm^(-3)
                                                             !< Until a realistic Nccn is provided, Nccns are assumed
                                                             !< as Nccn=100 for sea and Nccn=1000 for land

!--- near surface sea temperature model
    logical              :: nst_anl        = .false.         !< flag for NSSTM analysis in gcycle/sfcsub
    integer              :: lsea           = 0
    integer              :: nstf_name(5)   = (/0,0,1,0,5/)   !< flag 0 for no nst  1 for uncoupled nst  and 2 for coupled NST
                                                             !< nstf_name(1) : 0 = NSSTM off, 1 = NSSTM on but uncoupled
                                                             !<                2 = NSSTM on and coupled
                                                             !< nstf_name(2) : 1 = NSSTM spin up on, 0 = NSSTM spin up off
                                                             !< nstf_name(3) : 1 = NSSTM analysis on, 0 = NSSTM analysis off
                                                             !< nstf_name(4) : zsea1 in mm
                                                             !< nstf_name(5) : zsea2 in mm
!--- background vertical diffusion
    real(kind=kind_phys) :: xkzm_m         = 1.0d0           !< [in] bkgd_vdif_m  background vertical diffusion for momentum
    real(kind=kind_phys) :: xkzm_h         = 1.0d0           !< [in] bkgd_vdif_h  background vertical diffusion for heat q
    real(kind=kind_phys) :: xkzm_s         = 1.0d0           !< [in] bkgd_vdif_s  sigma threshold for background mom. diffusion
    real(kind=kind_phys) :: xkzminv        = 0.3             !< diffusivity in inversion layers
    real(kind=kind_phys) :: moninq_fac     = 1.0             !< turbulence diffusion coefficient factor

!---Cellular automaton options
    integer              :: nca            = 1
    integer              :: ncells         = 5
    integer              :: nlives         = 10
    real(kind=kind_phys) :: nfracseed      = 0.5
    integer              :: nseed          = 100000
    integer              :: iseed_ca       = 0
    integer              :: nspinup        = 1
    logical              :: do_ca          = .false.
    logical              :: ca_sgs         = .false. 
    logical              :: ca_global      = .false.
    logical              :: ca_smooth      = .false.
    logical              :: isppt_deep     = .false.
    real(kind=kind_phys) :: nthresh        = 0.0
  
    !--- IAU options
    real(kind=kind_phys)  :: iau_delthrs = 6                 ! iau time interval (to scale increments)
    character(len=240)    :: iau_inc_files(7)=''             ! list of increment files
    real(kind=kind_phys)  :: iaufhrs(7)=-1                   ! forecast hours associated with increment files
    logical               :: iau_filter_increments = .false. ! filter IAU increments

!--- debug flag
    logical              :: debug          = .false.
    logical              :: pre_rad        = .false.         !< flag for testing purpose

!  max and min lon and lat for critical relative humidity
    integer :: max_lon=5000, max_lat=2000, min_lon=192, min_lat=94
    real(kind=kind_phys) :: rhcmax = 0.9999999               !< max critical rel. hum.

!--- stochastic physics control parameters
    logical :: do_sppt      = .false.
    logical :: use_zmtnblck = .false.
    logical :: do_shum      = .false.
    logical :: do_skeb      = .false.
    integer :: skeb_npass = 11
    logical :: do_sfcperts = .false.   ! mg, sfc-perts
    integer :: nsfcpert    =  6        ! mg, sfc-perts
    real(kind=kind_phys) :: pertz0   = -999.
    real(kind=kind_phys) :: pertzt   = -999.
    real(kind=kind_phys) :: pertshc  = -999.
    real(kind=kind_phys) :: pertlai  = -999.
    real(kind=kind_phys) :: pertalb  = -999.
    real(kind=kind_phys) :: pertvegf = -999.
!--- END NAMELIST VARIABLES

    NAMELIST /gfs_physics_nml/                                                              &
                          !--- general parameters
                               fhzero, ldiag3d, lssav, fhcyc, lgocart, fhgoc3d,             &
                               thermodyn_id, sfcpress_id,                                   &
                          !--- coupling parameters
                               cplflx, cplwav, cplchm, lsidea,                              &
                          !--- radiation parameters
                               fhswr, fhlwr, levr, nfxr, aero_in, iflip, isol, ico2, ialb,  &
                               isot, iems, iaer, icliq_sw, iovr_sw, iovr_lw, ictm, isubc_sw,&
                               isubc_lw, crick_proof, ccnorm, lwhtr, swhtr,                 &
                          ! IN CCN forcing
                               iccn,                                                        &
                          !--- microphysical parameterizations
                               ncld, imp_physics, psautco, prautco, evpco, wminco,          &
#ifdef CCPP
                               fprcp, pdfflag, mg_dcs, mg_qcvar, mg_ts_auto_ice, mg_rhmini, &
                               effr_in, tf, tcr,                                            &
#else
                               fprcp, pdfflag, mg_dcs, mg_qcvar, mg_ts_auto_ice, effr_in,   &
#endif
                               microp_uniform, do_cldice, hetfrz_classnuc,                  &
                               mg_do_graupel, mg_do_hail, mg_nccons, mg_nicons, mg_ngcons,  &
                               mg_ncnst, mg_ninst, mg_ngnst, sed_supersat, do_sb_physics,   &
                               mg_alf,   mg_qcmin, mg_do_ice_gmao, mg_do_liq_liu,           &
                               ltaerosol, lradar, ttendlim, lgfdlmprad,                     &
                          !--- max hourly
                               avg_max_length,                                              &
                          !--- land/surface model control
#ifdef CCPP
                               lsm, lsoil, lsoil_lsm, nmtvr, ivegsrc, mom4ice, use_ufo,     &
#else
                               lsm, lsoil, nmtvr, ivegsrc, mom4ice, use_ufo,                &
#endif
                          !--- physical parameterizations
                               ras, trans_trac, old_monin, cnvgwd, mstrat, moist_adj,       &
                               cscnv, cal_pre, do_aw, do_shoc, shocaftcnv, shoc_cld,        &
#ifdef CCPP
                               oz_phys, oz_phys_2015,                                       &
                               do_mynnedmf, do_mynnsfclay,                                  &
                               bl_mynn_cloudpdf, bl_mynn_edmf, bl_mynn_edmf_mom,            &
                               bl_mynn_edmf_tke, bl_mynn_edmf_part, bl_mynn_cloudmix,       &
                               bl_mynn_mixqt, icloud_bl, bl_mynn_tkeadvect,                 &
#endif
                               h2o_phys, pdfcld, shcnvcw, redrag, hybedmf, satmedmf,        &
                               shinhong, do_ysu, dspheat, dspheat, lheatstrg, cnvcld,       &
                               random_clds, shal_cnv, imfshalcnv, imfdeepcnv, do_deep, jcap,&
                               cs_parm, flgmin, cgwf, ccwf, cdmbgwd, sup, ctei_rm, crtrh,   &
                               dlqf, rbcr, shoc_parm, psauras, prauras, wminras,            &
#ifdef CCPP
                               do_sppt, do_shum, do_skeb, do_sfcperts,                      &
#endif
                          !--- Rayleigh friction
                               prslrd0, ral_ts,                                             &
                          !--- mass flux deep convection
                               clam_deep, c0s_deep, c1_deep, betal_deep,                    &
                               betas_deep, evfact_deep, evfactl_deep, pgcon_deep,           &
                               asolfac_deep,                                                &
                          !--- mass flux shallow convection
                               clam_shal, c0s_shal, c1_shal, pgcon_shal, asolfac_shal,      &
                          !--- near surface sea temperature model
                               nst_anl, lsea, nstf_name,                                    &
                          !    background vertical diffusion
                               xkzm_m, xkzm_h, xkzm_s, xkzminv, moninq_fac,                 &
                          !--- cellular automata                         
                               nca, ncells, nlives, nfracseed,nseed, nthresh, do_ca,        &
                               ca_sgs, ca_global,iseed_ca,ca_smooth,isppt_deep,nspinup,     &
                          !--- IAU
                               iau_delthrs,iaufhrs,iau_inc_files,iau_filter_increments,     &
                          !--- debug options
                               debug, pre_rad,                                              &
                          !--- parameter range for critical relative humidity
                               max_lon, max_lat, min_lon, min_lat, rhcmax,                  &
                               phys_version

!--- other parameters
    integer :: nctp    =  0                !< number of cloud types in CS scheme
    logical :: gen_coord_hybrid = .false.  !< for Henry's gen coord

!--- SHOC parameters
    integer :: nshoc_2d  = 0  !< number of 2d fields for SHOC
    integer :: nshoc_3d  = 0  !< number of 3d fields for SHOC

!--- convective clouds
    integer :: ncnvcld3d = 0       !< number of convective 3d clouds fields

!--- read in the namelist
#ifdef INTERNAL_FILE_NML
    Model%input_nml_file => input_nml_file
    read(Model%input_nml_file, nml=gfs_physics_nml)
#else
    inquire (file=trim(fn_nml), exist=exists)
    if (.not. exists) then
      write(6,*) 'GFS_namelist_read:: namelist file: ',trim(fn_nml),' does not exist'
      stop
    else
      open (unit=nlunit, file=fn_nml, action='READ', status='OLD', iostat=ios)
    endif
    rewind(nlunit)
    read (nlunit, nml=gfs_physics_nml)
    close (nlunit)
#endif
!--- write version number and namelist to log file ---
    if (me == master) then
      write(logunit, '(a80)') '================================================================================'
      write(logunit, '(a64)') phys_version
      write(logunit, nml=gfs_physics_nml)
    endif

!--- MPI parameters
    Model%me               = me
    Model%master           = master
#ifdef CCPP
    Model%communicator     = communicator
    Model%ntasks           = ntasks
    Model%nthreads         = nthreads
#endif
    Model%nlunit           = nlunit
    Model%fn_nml           = fn_nml
#ifdef CCPP
    Model%logunit          = logunit
#endif
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
    Model%tile_num         = tile_num
    Model%isc              = isc
    Model%jsc              = jsc
    Model%nx               = nx
    Model%ny               = ny
    Model%levs             = levs
#ifdef CCPP
    allocate(Model%ak(1:size(ak)))
    allocate(Model%bk(1:size(bk)))
    Model%ak               = ak
    Model%bk               = bk
#endif
    Model%cnx              = cnx
    Model%cny              = cny
    Model%lonr             = gnx         ! number longitudinal points
    Model%latr             = gny         ! number of latitudinal points from pole to pole
#ifdef CCPP
    allocate(Model%blksz(1:size(blksz)))
    Model%blksz            = blksz
#endif

!--- coupling parameters
    Model%cplflx           = cplflx
    Model%cplwav           = cplwav
    Model%cplchm           = cplchm

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
    if (Model%aero_in) then
      ntrcaer = ntrcaerm
    else
      ntrcaer = 1
    endif
    Model%iccn             = iccn
    if (Model%aero_in) Model%iccn = .false.
    ! further down: set Model%iccn to .false.
    ! for all microphysics schemes except
    ! MG2/3 (these are the only ones using ICCN)
    Model%iflip            = iflip
    Model%isol             = isol
    Model%ico2             = ico2
    Model%ialb             = ialb
    Model%iems             = iems
    Model%iaer             = iaer
    Model%icliq_sw         = icliq_sw
    Model%iovr_sw          = iovr_sw
    Model%iovr_lw          = iovr_lw
    Model%ictm             = ictm
    Model%isubc_sw         = isubc_sw
    Model%isubc_lw         = isubc_lw
    Model%crick_proof      = crick_proof
    Model%ccnorm           = ccnorm
    Model%lwhtr            = lwhtr
    Model%swhtr            = swhtr
#ifdef CCPP
    ! The CCPP versions of the RRTMG lw/sw schemes are configured
    ! such that lw and sw heating rate are output, i.e. they rely
    ! on the corresponding arrays to be allocated.
    if (.not.lwhtr .or. .not.swhtr) then
      write(0,*) "Logic error, the CCPP version of RRTMG lwrad/swrad require the output" // &
             " of the lw/sw heating rates to be turned on (namelist options lwhtr and swhtr)"
      stop
    end if
#endif

!--- microphysical switch
    Model%ncld             = ncld
    Model%imp_physics      = imp_physics
    ! DH*
    ! turn off ICCN interpolation when MG2/3 are not used
    if (.not. Model%imp_physics==Model%imp_physics_mg) Model%iccn = .false.
    ! *DH
!--- Zhao-Carr MP parameters
    Model%psautco          = psautco
    Model%prautco          = prautco
    Model%evpco            = evpco
    Model%wminco           = wminco
!--- Max hourly
    Model%avg_max_length   = avg_max_length
!--- Morrison-Gettelman MP parameters
    Model%fprcp            = fprcp
    Model%pdfflag          = pdfflag
    Model%mg_dcs           = mg_dcs
    Model%mg_qcvar         = mg_qcvar
    Model%mg_ts_auto_ice   = mg_ts_auto_ice
#ifdef CCPP
    Model%mg_rhmini        = mg_rhmini
#endif
    Model%mg_alf           = mg_alf
    Model%mg_qcmin         = mg_qcmin
    Model%effr_in          = effr_in
    Model%microp_uniform   = microp_uniform
    Model%do_cldice        = do_cldice
    Model%hetfrz_classnuc  = hetfrz_classnuc
    Model%mg_do_graupel    = mg_do_graupel
    Model%mg_do_hail       = mg_do_hail
    Model%mg_do_ice_gmao   = mg_do_ice_gmao
    Model%mg_do_liq_liu    = mg_do_liq_liu
    Model%mg_nccons        = mg_nccons
    Model%mg_nicons        = mg_nicons
    Model%mg_ngcons        = mg_ngcons
    Model%mg_ncnst         = mg_ncnst
    Model%mg_ninst         = mg_ninst
    Model%mg_ngnst         = mg_ngnst
    Model%sed_supersat     = sed_supersat
    Model%do_sb_physics    = do_sb_physics
    Model%mg_precip_frac_method  = mg_precip_frac_method
    Model%mg_berg_eff_factor     = mg_berg_eff_factor
#ifdef CCPP
    Model%tf               = tf
    Model%tcr              = tcr
    Model%tcrf             = 1.0/(tcr-tf)
#endif


!--- Thompson MP parameters
    Model%ltaerosol        = ltaerosol
    Model%lradar           = lradar
    Model%ttendlim         = ttendlim

!--- gfdl  MP parameters
    Model%lgfdlmprad       = lgfdlmprad

!--- land/surface model parameters
    Model%lsm              = lsm
    Model%lsoil            = lsoil
#ifdef CCPP
    if (lsoil_lsm==-1) then
      Model%lsoil_lsm      = lsoil
    else
      Model%lsoil_lsm      = lsoil_lsm
    end if
#endif
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
    Model%cs_parm          = cs_parm
    Model%do_shoc          = do_shoc
    Model%shoc_parm        = shoc_parm
    Model%shocaftcnv       = shocaftcnv
    Model%shoc_cld         = shoc_cld
#ifdef CCPP
    if (oz_phys .and. oz_phys_2015) then
       write(*,*) 'Logic error: can only use one ozone physics option (oz_phys or oz_phys_2015), not both. Exiting.'
       stop
    end if
    Model%oz_phys          = oz_phys
    Model%oz_phys_2015     = oz_phys_2015
#endif
    Model%h2o_phys         = h2o_phys
#ifdef CCPP
    ! To ensure that these values match what's in the physics,
    ! array sizes are compared during model init in GFS_phys_time_vary_init()
    !
    ! from module h2ointerp
    if (h2o_phys) then
       levh2o    = 72
       h2o_coeff = 3
    else
       levh2o    = 1
       h2o_coeff = 1
    end if
#endif
    Model%pdfcld           = pdfcld
    Model%shcnvcw          = shcnvcw
    Model%redrag           = redrag
    Model%hybedmf          = hybedmf
    Model%satmedmf         = satmedmf
    Model%shinhong         = shinhong
    Model%do_ysu           = do_ysu
    Model%dspheat          = dspheat
    Model%lheatstrg        = lheatstrg
    Model%cnvcld           = cnvcld
    Model%random_clds      = random_clds
    Model%shal_cnv         = shal_cnv
    Model%imfshalcnv       = imfshalcnv
    Model%imfdeepcnv       = imfdeepcnv
    Model%do_deep          = do_deep
    Model%nmtvr            = nmtvr
    Model%jcap             = jcap
    Model%flgmin           = flgmin
    Model%cgwf             = cgwf
    Model%ccwf             = ccwf
    Model%cdmbgwd          = cdmbgwd
    Model%sup              = sup
    Model%ctei_rm          = ctei_rm
    Model%crtrh            = crtrh
    Model%dlqf             = dlqf
    Model%psauras          = psauras
    Model%prauras          = prauras
    Model%wminras          = wminras
    Model%rbcr             = rbcr
#ifdef CCPP
    Model%do_mynnedmf       = do_mynnedmf
    Model%do_mynnsfclay     = do_mynnsfclay
    Model%bl_mynn_cloudpdf  = bl_mynn_cloudpdf
    Model%bl_mynn_mixlength = bl_mynn_mixlength
    Model%bl_mynn_edmf      = bl_mynn_edmf
    Model%bl_mynn_edmf_mom  = bl_mynn_edmf_mom
    Model%bl_mynn_edmf_tke  = bl_mynn_edmf_tke
    Model%bl_mynn_cloudmix  = bl_mynn_cloudmix
    Model%bl_mynn_mixqt     = bl_mynn_mixqt
    Model%bl_mynn_edmf_part = bl_mynn_edmf_part
    Model%bl_mynn_tkeadvect = bl_mynn_tkeadvect
    Model%grav_settling     = grav_settling
    Model%icloud_bl         = icloud_bl
#endif

!--- Rayleigh friction
    Model%prslrd0          = prslrd0
    Model%ral_ts           = ral_ts

!--- mass flux deep convection
    Model%clam_deep        = clam_deep
    Model%c0s_deep         = c0s_deep
    Model%c1_deep          = c1_deep
    Model%betal_deep       = betal_deep
    Model%betas_deep       = betas_deep
    Model%evfact_deep      = evfact_deep
    Model%evfactl_deep     = evfactl_deep
    Model%pgcon_deep       = pgcon_deep
    Model%asolfac_deep     = asolfac_deep

!--- mass flux shallow convection
    Model%clam_shal        = clam_shal
    Model%c0s_shal         = c0s_shal
    Model%c1_shal          = c1_shal
    Model%pgcon_shal       = pgcon_shal
    Model%asolfac_shal     = asolfac_shal

!--- near surface sea temperature model
    Model%nst_anl          = nst_anl
    Model%lsea             = lsea
    Model%nstf_name        = nstf_name

!--- backgroud vertical diffusion
    Model%xkzm_m           = xkzm_m
    Model%xkzm_h           = xkzm_h
    Model%xkzm_s           = xkzm_s
    Model%xkzminv          = xkzminv
    Model%moninq_fac       = moninq_fac

!--- stochastic physics options
    ! DH* 20180730
    ! For the standard/non-CCPP build, do_sppt, do_shum, do_skeb and do_sfcperts
    ! are set to false here and updated later as part of init_stochastic_physics,
    ! depending on values of other namelist parameters. Since these values are used
    ! as conditionals for allocating components of Coupling and other DDTs, the call
    ! to init_stochastic_physics is placed between creating the "Model" DDT and the
    ! other DDTs.
    ! This is confusing and not compatible with CCPP, because the CCPP physics init
    ! can happen only after ALL DDTs are created and added to the CCPP data structure
    ! (cdata). Hence, for CCPP do_sppt, do_shum, do_skeb and do_sfcperts are additional
    ! namelist variables in group physics that are parsed here and then compared in
    ! stochastic_physics_init (the CCPP version of init_stochastic_physics) to the
    ! stochastic physics namelist parameters to ensure consistency.
    Model%do_sppt          = do_sppt
    Model%use_zmtnblck     = use_zmtnblck
    Model%do_shum          = do_shum
    Model%do_skeb          = do_skeb
    Model%do_sfcperts      = do_sfcperts ! mg, sfc-perts
    Model%nsfcpert         = nsfcpert    ! mg, sfc-perts
    Model%pertz0           = pertz0
    Model%pertzt           = pertzt
    Model%pertshc          = pertshc
    Model%pertlai          = pertlai
    Model%pertalb          = pertalb
    Model%pertvegf         = pertvegf
    ! *DH 20180730

    !--- cellular automata options
    Model%nca              = nca
    Model%ncells           = ncells
    Model%nlives           = nlives
    Model%nfracseed        = nfracseed
    Model%nseed            = nseed
    Model%ca_global        = ca_global
    Model%do_ca            = do_ca
    Model%ca_sgs           = ca_sgs
    Model%iseed_ca         = iseed_ca
    Model%ca_smooth        = ca_smooth
    Model%isppt_deep       = isppt_deep
    Model%nspinup          = nspinup  
    Model%nthresh          = nthresh 

    ! IAU flags
    !--- iau parameters
    Model%iaufhrs         = iaufhrs
    Model%iau_inc_files   = iau_inc_files
    Model%iau_delthrs     = iau_delthrs
    Model%iau_filter_increments = iau_filter_increments

!--- tracer handling
    Model%ntrac            = size(tracer_names)
#ifdef CCPP
    Model%ntracp1          = Model%ntrac + 1
#endif
    allocate (Model%tracer_names(Model%ntrac))
    Model%tracer_names(:)  = tracer_names(:)
#ifdef CCPP
    Model%ntqv             = 1
#endif
#ifdef MULTI_GASES
    Model%nto              = get_tracer_index(Model%tracer_names, 'spfo',        Model%me, Model%master, Model%debug)
    Model%nto2             = get_tracer_index(Model%tracer_names, 'spfo2',       Model%me, Model%master, Model%debug)
    Model%ntoz             = get_tracer_index(Model%tracer_names, 'spfo3',       Model%me, Model%master, Model%debug)
#else
    Model%ntoz             = get_tracer_index(Model%tracer_names, 'o3mr',       Model%me, Model%master, Model%debug)
#endif
    Model%ntcw             = get_tracer_index(Model%tracer_names, 'liq_wat',    Model%me, Model%master, Model%debug)
    Model%ntiw             = get_tracer_index(Model%tracer_names, 'ice_wat',    Model%me, Model%master, Model%debug)
    Model%ntrw             = get_tracer_index(Model%tracer_names, 'rainwat',    Model%me, Model%master, Model%debug)
    Model%ntsw             = get_tracer_index(Model%tracer_names, 'snowwat',    Model%me, Model%master, Model%debug)
    Model%ntgl             = get_tracer_index(Model%tracer_names, 'graupel',    Model%me, Model%master, Model%debug)
    Model%ntclamt          = get_tracer_index(Model%tracer_names, 'cld_amt',    Model%me, Model%master, Model%debug)
    Model%ntlnc            = get_tracer_index(Model%tracer_names, 'water_nc',   Model%me, Model%master, Model%debug)
    Model%ntinc            = get_tracer_index(Model%tracer_names, 'ice_nc',     Model%me, Model%master, Model%debug)
    Model%ntrnc            = get_tracer_index(Model%tracer_names, 'rain_nc',    Model%me, Model%master, Model%debug)
    Model%ntsnc            = get_tracer_index(Model%tracer_names, 'snow_nc',    Model%me, Model%master, Model%debug)
    Model%ntgnc            = get_tracer_index(Model%tracer_names, 'graupel_nc', Model%me, Model%master, Model%debug)
    Model%ntke             = get_tracer_index(Model%tracer_names, 'sgs_tke',    Model%me, Model%master, Model%debug)
    Model%ntwa             = get_tracer_index(Model%tracer_names, 'liq_aero',   Model%me, Model%master, Model%debug)
    Model%ntia             = get_tracer_index(Model%tracer_names, 'ice_aero',   Model%me, Model%master, Model%debug)
    Model%ntchm            = 0
    Model%ntchs            = get_tracer_index(Model%tracer_names, 'so2',        Model%me, Model%master, Model%debug)
    if (Model%ntchs > 0) then
      Model%ntchm          = get_tracer_index(Model%tracer_names, 'pp10',       Model%me, Model%master, Model%debug)
      if (Model%ntchm > 0) then
        Model%ntchm = Model%ntchm - Model%ntchs + 1
        allocate(Model%ntdiag(Model%ntchm))
        ! -- turn on all tracer diagnostics to .true. by default, except for so2
        Model%ntdiag(1)  = .false.
        Model%ntdiag(2:) = .true.
        ! -- turn off diagnostics for DMS
        n = get_tracer_index(Model%tracer_names, 'DMS', Model%me, Model%master, Model%debug) - Model%ntchs + 1
        if (n > 0) Model%ntdiag(n) = .false.
        ! -- turn off diagnostics for msa
        n = get_tracer_index(Model%tracer_names, 'msa', Model%me, Model%master, Model%debug) - Model%ntchs + 1
        if (n > 0) Model%ntdiag(n) = .false.
      endif
    endif

#ifdef CCPP
    ! To ensure that these values match what's in the physics,
    ! array sizes are compared during model init in GFS_phys_time_vary_init()
    !
    ! from module ozinterp
    if (Model%ntoz>0) then
       if (Model%oz_phys) then
          levozp   = 80
          oz_coeff = 4
       else if (Model%oz_phys_2015) then
          levozp   = 53
          oz_coeff = 6
       else
          write(*,*) 'Logic error, ntoz>0 but no ozone physics selected'
          stop
       end if
    else
       if (Model%oz_phys .or. Model%oz_phys_2015) then
          write(*,*) 'Logic error, ozone physics are selected, but ntoz<=0'
          stop
       else
          levozp   = 1
          oz_coeff = 0
       end if
    end if
#endif

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
#ifdef CCPP
    Model%first_time_step  = .true.
    Model%restart          = restart
    Model%hydrostatic      = hydrostatic
#endif
    Model%jdat(1:8)        = jdat(1:8)
#ifdef CCPP
    Model%sec              = 0
    allocate(Model%si(Model%levr+1))
    !--- Define sigma level for radiation initialization
    !--- The formula converting hybrid sigma pressure coefficients to sigma coefficients follows Eckermann (2009, MWR)
    !--- ps is replaced with p0. The value of p0 uses that in http://www.emc.ncep.noaa.gov/officenotes/newernotes/on461.pdf
    !--- ak/bk have been flipped from their original FV3 orientation and are defined sfc -> toa
    Model%si = (ak + bk * p_ref - ak(Model%levr+1)) / (p_ref - ak(Model%levr+1))
#endif

#if !defined(CCPP) || defined(HYBRID)
    ! Beware! The values set here reside in wam_f107_kp_mod and determine sizes of arrays
    ! inside that module. These arrays get used later in modules idea_tracer.f, idea_ion.f,
    ! idea_solar_heating.f, efield.f, and idea_composition.f.
    ! Since in wam_f107_kp_mod no default values are assigned to the four integers below, not
    ! setting them here can lead to memory corruption that are hard to detect.
!--- stored in wam_f107_kp module
    f107_kp_size      = 56
    f107_kp_skip_size = 0
    f107_kp_data_size = 56
    f107_kp_interval  = 10800
#endif

!--- BEGIN CODE FROM GFS_PHYSICS_INITIALIZE
!--- define physcons module variables
    tem     = con_rerth*con_rerth*(con_pi+con_pi)*con_pi
#ifdef CCPP
    Model%dxmax  = log(tem/(max_lon*max_lat))
    Model%dxmin  = log(tem/(min_lon*min_lat))
    Model%dxinv  = 1.0d0 / (Model%dxmax-Model%dxmin)
    Model%rhcmax = rhcmax
    if (Model%me == Model%master) write(*,*)' dxmax=',Model%dxmax,' dxmin=',Model%dxmin,' dxinv=',Model%dxinv, &
       'max_lon=',max_lon,' max_lat=',max_lat,' min_lon=',min_lon,' min_lat=',min_lat,       &
       ' rhc_max=',Model%rhcmax
#else
    dxmax   = log(tem/(max_lon*max_lat))
    dxmin   = log(tem/(min_lon*min_lat))
    dxinv   = 1.0d0 / (dxmax-dxmin)
    rhc_max = rhcmax
    if (Model%me == Model%master) write(*,*)' dxmax=',dxmax,' dxmin=',dxmin,' dxinv=',dxinv, &
       'max_lon=',max_lon,' max_lat=',max_lat,' min_lon=',min_lon,' min_lat=',min_lat,       &
       ' rhc_max=',rhc_max
#endif

!--- set nrcm

#if !defined(CCPP) || defined(HYBRID)
    if (Model%ras) then
      Model%nrcm = min(nrcmax, Model%levs-1) * (Model%dtp/1200.d0) + 0.10001d0
    else
      Model%nrcm = 2
    endif
#else
    Model%nrcm = 2
#endif

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
      Model%satmedmf   = .false.
      if (Model%me == Model%master) print *,' Simplified Higher Order Closure Model used for', &
                                            ' Boundary layer and Shallow Convection',          &
                                            ' nshoc_3d=',Model%nshoc_3d,                       &
                                            ' nshoc_2d=',Model%nshoc_2d,                       &
                                            ' ntke=',Model%ntke,' shoc_parm=',shoc_parm
    endif

#ifdef CCPP
    !--- mynn-edmf scheme
    if (Model%bl_mynn_edmf > 0) then
      Model%do_shoc    = .false.
!      Model%shal_cnv   = .false.
!      Model%imfshalcnv = -1
      Model%hybedmf    = .false.
      Model%satmedmf   = .false.
      if (Model%me == Model%master) print *,' MYNN-EDMF scheme is used for both',                &
                                            ' boundary layer turbulence and shallow convection', &
                                            ' bl_mynn_cloudpdf=',Model%bl_mynn_cloudpdf,         &
                                            ' bl_mynn_mixlength=',Model%bl_mynn_mixlength,       &
                                            ' bl_mynn_edmf=',Model%bl_mynn_edmf
    endif
#endif

!--- set number of cloud types
    if (Model%cscnv) then
      Model%nctp = nint(Model%cs_parm(5))
      Model%nctp = max(Model%nctp,10)
      if (Model%cs_parm(7) < 0.0) Model%cs_parm(7) = Model%dtp
      Model%do_awdd  = Model%do_aw .and. Model%cs_parm(6) > 0.0
!     Model%flx_form = Model%do_aw .and. Model%cs_parm(8) > 0.0
      Model%flx_form = Model%cs_parm(8) > 0.0
    endif
    Model%nctp = max(Model%nctp,1)

!--- output information about the run
    if (Model%me == Model%master) then
      if (Model%lsm == 1) then
        print *,' NOAH Land Surface Model used'
#ifdef CCPP
      elseif (Model%lsm == Model%lsm_ruc) then
        print *,' RUC Land Surface Model used'
#else
      elseif (Model%lsm == Model%lsm_ruc) then
        print *,' RUC Land Surface Model only available through CCPP - job aborted'
        stop
#endif
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
      if (Model%do_deep) then
#ifdef CCPP
        ! Consistency check for NTDK convection: deep and shallow convection are bundled
        ! and cannot be combined with any other deep or shallow convection scheme
        if ( (Model%imfdeepcnv == 4 .or. Model%imfshalcnv == 4) .and. &
            .not. (Model%imfdeepcnv == 4 .and. Model%imfshalcnv == 4) ) then
            write(0,*) "Logic error: if NTDK deep convection is used, must also use NTDK shallow convection (and vice versa)"
            stop
        end if
#else
        if (Model%imfdeepcnv == 3 .or. Model%imfshalcnv == 3) then
            write(0,*) "Error, GF convection scheme only available through CCPP"
            stop
        else if (Model%imfdeepcnv == 4 .or. Model%imfshalcnv == 4) then
            write(0,*) "Error, NTDK convection scheme only available through CCPP"
            stop
        end if
#endif
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
            elseif(Model%imfdeepcnv == 3) then
               print *,' Grell-Freitas scale & aerosol-aware mass-flux deep conv scheme'
            elseif(Model%imfdeepcnv == 4) then
               print *,' New Tiedtke cumulus scheme'
            endif
          endif
        else
          if (Model%do_aw) then
            print *,'Chikira-Sugiyama convection scheme with Arakawa-Wu'&
     &,                ' unified parameterization used'
          else
              print *,'Chikira-Sugiyama convection scheme used'
          endif
          print *,' cs_parm=',Model%cs_parm,' nctp=',Model%nctp
        endif
      else
        print*, ' Deep convection scheme disabled'
      endif
#ifdef CPPP
      if (.not. Model%old_monin .and. .not. Model%do_shoc .and. .not. Model%do_mynnedmf) print *,' New PBL scheme used'
#else
      if (.not. Model%old_monin .and. .not. Model%do_shoc) print *,' New PBL scheme used'
#endif
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
        elseif (Model%imfshalcnv == 3) then
          print *,' Grell-Freitas scale- & aerosol-aware mass-flux shallow conv scheme (2013)'
        elseif (Model%imfshalcnv == 4) then
          print *,' New Tiedtke cumulus scheme'
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
    Model%nleffr = -999
    Model%nieffr = -999
    Model%nreffr = -999
    Model%nseffr = -999
    Model%ngeffr = -999
    if (Model%imp_physics == Model%imp_physics_zhao_carr) then
      Model%npdf3d  = 0
      Model%num_p3d = 4
      Model%num_p2d = 3
      Model%shcnvcw = .false.
      Model%ncnd    = 1                   ! ncnd is the number of cloud condensate types
      if (Model%me == Model%master) print *,' Using Zhao/Carr/Sundqvist Microphysics'

    elseif (Model%imp_physics == Model%imp_physics_zhao_carr_pdf) then !Zhao Microphysics with PDF cloud
      Model%npdf3d  = 3
      Model%num_p3d = 4
      Model%num_p2d = 3
      Model%ncnd    = 1
      if (Model%me == Model%master) print *,'Using Zhao/Carr/Sundqvist Microphysics with PDF Cloud'

    else if (Model%imp_physics == 5) then        ! F-A goes here
      print *,' Ferrier Microphysics scheme has been deprecated - job aborted'
      stop

    elseif (Model%imp_physics == Model%imp_physics_wsm6) then !WSM6 microphysics
      Model%npdf3d  = 0
      Model%num_p3d = 3
      Model%num_p2d = 1
      Model%pdfcld  = .false.
      Model%shcnvcw = .false.
      Model%ncnd    = 5
      Model%nleffr = 1
      Model%nieffr = 2
      Model%nseffr = 3
      if (Model%me == Model%master) print *,' Using wsm6 microphysics'

    elseif (Model%imp_physics == Model%imp_physics_thompson) then !Thompson microphysics
      Model%npdf3d  = 0
      Model%num_p3d = 3
      Model%num_p2d = 1
      Model%pdfcld  = .false.
      Model%shcnvcw = .false.
      Model%ncnd    = 5
      Model%nleffr = 1
      Model%nieffr = 2
      Model%nseffr = 3
      if (Model%me == Model%master) print *,' Using Thompson double moment', &
                                          ' microphysics',' ltaerosol = ',Model%ltaerosol, &
                                          ' lradar =',Model%lradar,' ttendlim =',Model%ttendlim, &
                                          Model%num_p3d,Model%num_p2d

    else if (Model%imp_physics == Model%imp_physics_mg) then        ! Morrison-Gettelman Microphysics
      Model%npdf3d  = 0
      Model%num_p3d = 5
      Model%num_p2d = 1
      Model%pdfcld  = .false.
      Model%shcnvcw = .false.
      Model%ncnd    = 2
      Model%nleffr = 2
      Model%nieffr = 3
      Model%nreffr = 4
      Model%nseffr = 5
      if (abs(Model%fprcp) == 1) then
        Model%ncnd = 4
      elseif (Model%fprcp >= 2) then
        Model%ncnd = 4
        if (Model%mg_do_graupel .or. Model%mg_do_hail) then
          Model%ncnd = 5
        endif
        Model%num_p3d = 6
        Model%ngeffr = 6
      endif
      if (Model%me == Model%master)                                                                 &
         print *,' Using Morrison-Gettelman double moment microphysics',                            &
                 ' aero_in=',         Model%aero_in,         ' iccn=',          Model%iccn,         &
                 ' mg_dcs=',          Model%mg_dcs,          ' mg_qcvar=',      Model%mg_qcvar,     &
                 ' mg_ts_auto_ice=',  Model%mg_ts_auto_ice,  ' pdfflag=',       Model%pdfflag,      &
                 ' mg_do_graupel=',   Model%mg_do_graupel,   ' mg_do_hail=',    Model%mg_do_hail,   &
                 ' mg_nccons=',       Model%mg_nccons,       ' mg_nicon=',      Model%mg_nicons,    &
                 ' mg_ngcons=',       Model%mg_ngcons ,      ' mg_ncnst=',      Model%mg_ncnst,     &
                 ' mg_ninst=',        Model%mg_ninst ,       ' mg_ngnst=',      Model%mg_ngnst,     &
                 ' sed_supersat=',    Model%sed_supersat ,   ' do_sb_physics=', Model%do_sb_physics,&
                 ' microp_uniform=',  Model%microp_uniform,  ' do_cldice=',     Model%do_cldice,    &
                 ' hetfrz_classnuc=', Model%hetfrz_classnuc, ' ncnd=',          Model%ncnd,         &
                 ' mg_alf=',          Model%mg_alf,          ' mg_qcmin=',      Model%mg_qcmin,     &
                 ' mg_do_ice_gmao=',  Model%mg_do_ice_gmao,  ' mg_do_liq_liu=', Model%mg_do_liq_liu

    elseif (Model%imp_physics == Model%imp_physics_gfdl) then !GFDL microphysics
      Model%npdf3d  = 0
      Model%num_p3d = 1
      if(Model%effr_in) then
         Model%num_p3d = 5
      endif
      Model%nleffr = 1
      Model%nieffr = 2
      Model%nreffr = 3
      Model%nseffr = 4
      Model%ngeffr = 5
      Model%num_p2d = 1
      Model%pdfcld  = .false.
      Model%shcnvcw = .false.
      Model%ncnd    = 5
      if (Model%me == Model%master) print *,' avg_max_length=',Model%avg_max_length
      if (Model%me == Model%master) print *,' Using GFDL Cloud Microphysics'
    else
      if (Model%me == Model%master) print *,'Wrong imp_physics value. Job abort.'
      stop
    endif

    if(Model%ras     .or. Model%cscnv)  Model%cnvcld = .false.
#ifdef CCPP
    if(Model%do_shoc .or. Model%pdfcld .or. Model%do_mynnedmf) Model%cnvcld = .false.
#else
    if(Model%do_shoc .or. Model%pdfcld) Model%cnvcld = .false.
#endif
    if(Model%cnvcld) Model%ncnvcld3d = 1

!--- get cnvw and cnvc indices in phy_f3d
    Model%ncnvw = -999
    Model%ncnvc = -999
    if ((Model%npdf3d == 3) .and. (Model%num_p3d == 4)) then
      Model%ncnvw = Model%num_p3d + 2
      Model%ncnvc = Model%ncnvw + 1
    elseif ((Model%npdf3d == 0) .and. (Model%ncnvcld3d == 1)) then
      Model%ncnvw = Model%num_p3d + 1
    endif

!--- derived totals for phy_f*d
    Model%ntot2d = Model%num_p2d + Model%nshoc_2d
    Model%ntot3d = Model%num_p3d + Model%nshoc_3d + Model%npdf3d + Model%ncnvcld3d
!
!   Unified cloud for SHOC and/or MG3
    Model%uni_cld = .false.
    Model%indcld  = -1
!   if (Model%shoc_cld .or. Model%ncld == 2 .or. Model%ntclamt > 0) then
    if (Model%imp_physics == Model%imp_physics_mg) then
      Model%uni_cld = .true.
      Model%indcld  = 1
    elseif (Model%shoc_cld) then
      Model%uni_cld = .true.
      Model%indcld  = Model%ntot3d - 2
    endif

    if (me == Model%master)                                                     &
      write(0,*) ' num_p3d=',   Model%num_p3d,   ' num_p2d=',  Model%num_p2d,   &
                 ' crtrh=',     Model%crtrh,     ' npdf3d=',   Model%npdf3d,    &
                 ' pdfcld=',    Model%pdfcld,    ' shcnvcw=',  Model%shcnvcw,   &
                 ' cnvcld=',    Model%cnvcld,    ' ncnvcld3d=',Model%ncnvcld3d, &
                 ' do_shoc=',   Model%do_shoc,   ' nshoc3d=',  Model%nshoc_3d,  &
                 ' nshoc_2d=',  Model%nshoc_2d,  ' shoc_cld=', Model%shoc_cld,  & 
                 ' uni_cld=',   Model%uni_cld,                                  &
                 ' ntot3d=',    Model%ntot3d,    ' ntot2d=',   Model%ntot2d,    &
                 ' shocaftcnv=',Model%shocaftcnv,' indcld=',   Model%indcld,    &
                 ' shoc_parm=', Model%shoc_parm,                                &
                 ' ncnvw=',     Model%ncnvw,     ' ncnvc=',     Model%ncnvc

!--- END CODE FROM COMPNS_PHYSICS


!--- BEGIN CODE FROM GLOOPR
!--- set up parameters for Xu & Randell's cloudiness computation (Radiation)

    Model%lmfshal  = (Model%shal_cnv .and. Model%imfshalcnv > 0)
#ifdef CCPP
    Model%lmfdeep2 = (Model%imfdeepcnv == 2 .or. Model%imfdeepcnv == 3 .or. Model%imfdeepcnv == 4)
#else
    Model%lmfdeep2 = (Model%imfdeepcnv == 2)
#endif
!--- END CODE FROM GLOOPR

!--- BEGIN CODE FROM GLOOPB
!--- set up random number seed needed for RAS and old SAS and when cal_pre=.true.
!    Model%imfdeepcnv < 0 when Model%ras = .true.

    if (Model%imfdeepcnv <= 0 .or. Model%cal_pre) then
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
#ifdef CCPP
      print *, ' communicator      : ', Model%communicator
#endif
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
#ifdef CCPP
      print *, ' blksz(1)          : ', Model%blksz(1)
      print *, ' blksz(nblks)      : ', Model%blksz(size(Model%blksz))
#endif
      print *, ' '
      print *, 'coupling parameters'
      print *, ' cplflx            : ', Model%cplflx
      print *, ' cplwav            : ', Model%cplwav
      print *, ' cplchm            : ', Model%cplchm
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
      print *, ' icliq_sw          : ', Model%icliq_sw
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
      print *, ' imp_physics       : ', Model%imp_physics
      print *, ' '

      if (Model%imp_physics == Model%imp_physics_zhao_carr .or. Model%imp_physics == Model%imp_physics_zhao_carr_pdf) then
        print *, ' Z-C microphysical parameters'
        print *, ' psautco           : ', Model%psautco
        print *, ' prautco           : ', Model%prautco
        print *, ' evpco             : ', Model%evpco
        print *, ' wminco            : ', Model%wminco
        print *, ' '
      endif
      if (Model%imp_physics == Model%imp_physics_wsm6 .or. Model%imp_physics == Model%imp_physics_thompson) then
        print *, ' Thompson microphysical parameters'
        print *, ' ltaerosol         : ', Model%ltaerosol
        print *, ' lradar            : ', Model%lradar
        print *, ' ttendlim          : ', Model%ttendlim
        print *, ' '
      endif
      if (Model%imp_physics == Model%imp_physics_mg) then
        print *, ' M-G microphysical parameters'
        print *, ' fprcp             : ', Model%fprcp
        print *, ' mg_dcs            : ', Model%mg_dcs
        print *, ' mg_qcvar          : ', Model%mg_qcvar
        print *, ' mg_ts_auto_ice    : ', Model%mg_ts_auto_ice
        print *, ' mg_alf            : ', Model%mg_alf
        print *, ' mg_qcmin          : ', Model%mg_qcmin
        print *, ' pdfflag           : ', Model%pdfflag
        print *, ' '
      endif
      if (Model%imp_physics == Model%imp_physics_gfdl) then
        print *, ' GFDL microphysical parameters'
        print *, ' GFDL MP radiation inter: ', Model%lgfdlmprad
        print *, ' '
      endif

      print *, 'land/surface model parameters'
      print *, ' lsm               : ', Model%lsm
      print *, ' lsoil             : ', Model%lsoil
#ifdef CCPP
      print *, ' lsoil_lsm         : ', Model%lsoil_lsm
#endif
      print *, ' ivegsrc           : ', Model%ivegsrc
      print *, ' isot              : ', Model%isot
      print *, ' mom4ice           : ', Model%mom4ice
      print *, ' use_ufo           : ', Model%use_ufo
      print *, ' '
      print *, 'tuning parameters for physical parameterizations'
      print *, ' ras               : ', Model%ras
      if (Model%ras) then
        print *, ' psauras           : ', Model%psauras
        print *, ' prauras           : ', Model%prauras
        print *, ' wminras           : ', Model%wminras
      endif
      print *, ' flipv             : ', Model%flipv
      print *, ' trans_trac        : ', Model%trans_trac
      print *, ' old_monin         : ', Model%old_monin
      print *, ' cnvgwd            : ', Model%cnvgwd
      print *, ' mstrat            : ', Model%mstrat
      print *, ' moist_adj         : ', Model%moist_adj
      print *, ' cscnv             : ', Model%cscnv
      print *, ' cal_pre           : ', Model%cal_pre
      print *, ' do_aw             : ', Model%do_aw
      print *, ' flx_form          : ', Model%flx_form
      print *, ' do_shoc           : ', Model%do_shoc
      print *, ' shoc_parm         : ', Model%shoc_parm
      print *, ' shocaftcnv        : ', Model%shocaftcnv
      print *, ' shoc_cld          : ', Model%shoc_cld
      print *, ' uni_cld           : ', Model%uni_cld
      print *, ' h2o_phys          : ', Model%h2o_phys
      print *, ' pdfcld            : ', Model%pdfcld
      print *, ' shcnvcw           : ', Model%shcnvcw
      print *, ' redrag            : ', Model%redrag
      print *, ' hybedmf           : ', Model%hybedmf
      print *, ' satmedmf          : ', Model%satmedmf
      print *, ' shinhong          : ', Model%shinhong
      print *, ' do_ysu            : ', Model%do_ysu
      print *, ' dspheat           : ', Model%dspheat
      print *, ' lheatstrg         : ', Model%lheatstrg
      print *, ' cnvcld            : ', Model%cnvcld
      print *, ' random_clds       : ', Model%random_clds
      print *, ' shal_cnv          : ', Model%shal_cnv
      print *, ' imfshalcnv        : ', Model%imfshalcnv
      print *, ' imfdeepcnv        : ', Model%imfdeepcnv
      print *, ' do_deep           : ', Model%do_deep
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
      print *, ' rbcr              : ', Model%rbcr
#ifdef CCPP
      print *, ' do_mynnedmf       : ', Model%do_mynnedmf
      print *, ' do_mynnsfclay     : ', Model%do_mynnsfclay
#endif
      print *, ' '
      print *, 'Rayleigh friction'
      print *, ' prslrd0           : ', Model%prslrd0
      print *, ' ral_ts            : ', Model%ral_ts
      print *, ' '
      if (Model%imfdeepcnv >= 0) then
        print *, 'mass flux deep convection'
        print *, ' clam_deep         : ', Model%clam_deep
        print *, ' c0s_deep          : ', Model%c0s_deep
        print *, ' c1_deep           : ', Model%c1_deep
        print *, ' betal_deep        : ', Model%betal_deep
        print *, ' betas_deep        : ', Model%betas_deep
        print *, ' evfact_deep       : ', Model%evfact_deep
        print *, ' evfactl_deep      : ', Model%evfactl_deep
        print *, ' pgcon_deep        : ', Model%pgcon_deep
        print *, ' asolfac_deep      : ', Model%asolfac_deep
        print *, ' '
      endif
      if (Model%imfshalcnv >= 0) then
        print *, 'mass flux shallow convection'
        print *, ' clam_shal         : ', Model%clam_shal
        print *, ' c0s_shal          : ', Model%c0s_shal
        print *, ' c1_shal           : ', Model%c1_shal
        print *, ' pgcon_shal        : ', Model%pgcon_shal
        print *, ' asolfac_shal      : ', Model%asolfac_shal
      endif
      print *, ' '
      print *, 'near surface sea temperature model'
      print *, ' nst_anl           : ', Model%nst_anl
      print *, ' nstf_name         : ', Model%nstf_name
      print *, ' lsea              : ', Model%lsea
      print *, ' '
      print *, 'background vertical diffusion coefficients'
      print *, ' xkzm_m            : ', Model%xkzm_m
      print *, ' xkzm_h            : ', Model%xkzm_h
      print *, ' xkzm_s            : ', Model%xkzm_s
      print *, ' xkzminv           : ', Model%xkzminv
      print *, ' moninq_fac        : ', Model%moninq_fac
      print *, ' '
      print *, 'stochastic physics'
      print *, ' do_sppt           : ', Model%do_sppt
      print *, ' do_shum           : ', Model%do_shum
      print *, ' do_skeb           : ', Model%do_skeb
      print *, ' do_sfcperts       : ', Model%do_sfcperts
      print *, ' '
      print *, 'cellular automata'
      print *, ' nca               : ', Model%ncells
      print *, ' ncells            : ', Model%ncells
      print *, ' nlives            : ', Model%nlives
      print *, ' nfracseed         : ', Model%nfracseed
      print *, ' nseed             : ', Model%nseed
      print *, ' ca_global         : ', Model%ca_global
      print *, ' ca_sgs            : ', Model%ca_sgs
      print *, ' do_ca             : ', Model%do_ca
      print *, ' iseed_ca          : ', Model%iseed_ca
      print *, ' ca_smooth         : ', Model%ca_smooth
      print *, ' isppt_deep        : ', Model%isppt_deep
      print *, ' nspinup           : ', Model%nspinup
      print *, ' nthresh           : ', Model%nthresh
      print *, ' '
      print *, 'tracers'
      print *, ' tracer_names      : ', Model%tracer_names
      print *, ' ntrac             : ', Model%ntrac
#ifdef CCPP
      print *, ' ntqv              : ', Model%ntqv
#endif
      print *, ' ntoz              : ', Model%ntoz
      print *, ' ntcw              : ', Model%ntcw
      print *, ' ntiw              : ', Model%ntiw
      print *, ' ntrw              : ', Model%ntrw
      print *, ' ntsw              : ', Model%ntsw
      print *, ' ntgl              : ', Model%ntgl
      print *, ' ntclamt           : ', Model%ntclamt
      print *, ' ntlnc             : ', Model%ntlnc
      print *, ' ntinc             : ', Model%ntinc
      print *, ' ntrnc             : ', Model%ntrnc
      print *, ' ntsnc             : ', Model%ntsnc
      print *, ' ntgnc             : ', Model%ntgnc
      print *, ' ntke              : ', Model%ntke
      print *, ' nto               : ', Model%nto
      print *, ' nto2              : ', Model%nto2
      print *, ' ntwa              : ', Model%ntwa
      print *, ' ntia              : ', Model%ntia
      print *, ' ntchm             : ', Model%ntchm
      print *, ' ntchs             : ', Model%ntchs
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
#ifdef CCPP
      print *, ' sec               : ', Model%sec
      print *, ' si                : ', Model%si
      print *, ' first_time_step   : ', Model%first_time_step
      print *, ' restart           : ', Model%restart
      print *, ' hydrostatic       : ', Model%hydrostatic
#endif
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
    allocate (Grid%xlon_d (IM))
    allocate (Grid%sinlat (IM))
    allocate (Grid%coslat (IM))
    allocate (Grid%area   (IM))
    allocate (Grid%dx     (IM))

    Grid%xlon   = clear_val
    Grid%xlat   = clear_val
    Grid%xlat_d = clear_val
    Grid%xlon_d = clear_val
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

!--- iccn active
    if ( Model%iccn ) then
      allocate (Grid%ddy_ci    (IM))
      allocate (Grid%jindx1_ci (IM))
      allocate (Grid%jindx2_ci (IM))
      allocate (Grid%ddx_ci    (IM))
      allocate (Grid%iindx1_ci (IM))
      allocate (Grid%iindx2_ci (IM))
    endif

!--- iaerclm active
    if ( Model%aero_in ) then
      allocate (Grid%ddy_aer   (IM))
      allocate (Grid%jindx1_aer(IM))
      allocate (Grid%jindx2_aer(IM))
      allocate (Grid%ddx_aer   (IM))
      allocate (Grid%iindx1_aer(IM))
      allocate (Grid%iindx2_aer(IM))
    endif
 end subroutine grid_create


!--------------------
! GFS_tbd_type%create
!--------------------
#if !defined(CCPP) || defined(HYBRID)
  subroutine tbd_create (Tbd, IM, BLKNO, Model)
#else
  subroutine tbd_create (Tbd, IM, Model)
#endif

    implicit none

    class(GFS_tbd_type)                :: Tbd
    integer,                intent(in) :: IM
#if !defined(CCPP) || defined(HYBRID)
    integer,                intent(in) :: BLKNO
#endif
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

!--- ccn and in needs
    allocate (Tbd%in_nm  (IM,Model%levs))
    allocate (Tbd%ccn_nm (IM,Model%levs))
    Tbd%in_nm  = clear_val
    Tbd%ccn_nm = clear_val

!--- aerosol fields
    allocate (Tbd%aer_nm  (IM,Model%levs,ntrcaer))
    Tbd%aer_nm = clear_val

#ifdef CCPP
! DH* TODO - MOVE THIS TO a block-vector dependent structure in GFS_control?
! e.g. GFS_Control%imap(blk), GFS_Control%jmap(blk), or ii instead if imap etc? *DH
!--- maps of local index ix to global indices i and j for this block
    allocate (Tbd%imap (IM))
    allocate (Tbd%jmap (IM))
    Tbd%imap = 0
    Tbd%jmap = 0
#endif

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

    Tbd%phy_fctd = clear_val
    Tbd%phy_f2d  = clear_val
    Tbd%phy_f3d  = clear_val
!   if (Model%do_shoc) Tbd%phy_f3d(:,1,Model%ntot3d-1) = 3.0
!   if (Model%do_shoc) Tbd%phy_f3d(:,:,Model%ntot3d-1) = 1.0

#if !defined(CCPP) || defined(HYBRID)
    Tbd%blkno = BLKNO
#endif

#ifdef CCPP
    allocate (Tbd%htlwc (IM,Model%levr+LTP))
    allocate (Tbd%htlw0 (IM,Model%levr+LTP))
    allocate (Tbd%htswc (IM,Model%levr+LTP))
    allocate (Tbd%htsw0 (IM,Model%levr+LTP))

    Tbd%htlwc = clear_val
    Tbd%htlw0 = clear_val
    Tbd%htswc = clear_val
    Tbd%htsw0 = clear_val

    if (Model%imfdeepcnv == 3 .or. Model%imfdeepcnv == 4) then
       allocate(Tbd%forcet(IM, Model%levs))
       allocate(Tbd%forceq(IM, Model%levs))
       allocate(Tbd%prevst(IM, Model%levs))
       allocate(Tbd%prevsq(IM, Model%levs))
       Tbd%forcet = clear_val
       Tbd%forceq = clear_val
       Tbd%prevst = clear_val
       Tbd%prevsq = clear_val
   end if

   if (Model%imfdeepcnv == 3) then
      allocate(Tbd%cactiv(IM))
      Tbd%cactiv = zero
  end if

   if (Model%lsm == Model%lsm_ruc) then
       allocate(Tbd%raincprv  (IM))
       allocate(Tbd%rainncprv (IM))
       allocate(Tbd%iceprv    (IM))
       allocate(Tbd%snowprv   (IM))
       allocate(Tbd%graupelprv(IM))
       Tbd%raincprv   = clear_val
       Tbd%rainncprv  = clear_val
       Tbd%iceprv     = clear_val
       Tbd%snowprv    = clear_val
       Tbd%graupelprv = clear_val
    end if

    !--- MYNN variables:
    if (Model%do_mynnedmf) then
       !print*,"Allocating all MYNN-EDMF variables:"
       allocate (Tbd%cldfra_bl (IM,Model%levs))
       allocate (Tbd%qc_bl     (IM,Model%levs))
       allocate (Tbd%el_pbl    (IM,Model%levs))
       allocate (Tbd%sh3d      (IM,Model%levs))
       allocate (Tbd%qke       (IM,Model%levs))
       allocate (Tbd%tsq       (IM,Model%levs))
       allocate (Tbd%qsq       (IM,Model%levs))
       allocate (Tbd%cov       (IM,Model%levs))
       !print*,"Allocating all MYNN-EDMF variables:"
       Tbd%cldfra_bl     = clear_val
       Tbd%qc_bl         = clear_val
       Tbd%el_pbl        = clear_val
       Tbd%sh3d          = clear_val
       Tbd%qke           = zero
       Tbd%tsq           = clear_val
       Tbd%qsq           = clear_val
       Tbd%cov           = clear_val
    end if
#endif

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

!
    logical, save :: linit

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
#ifdef CCPP
    allocate (Diag%snowfallac  (IM))
    allocate (Diag%acsnow      (IM))
#endif
    allocate (Diag%soilm   (IM))
    allocate (Diag%tmpmin  (IM))
    allocate (Diag%tmpmax  (IM))
    allocate (Diag%dusfc   (IM))
    allocate (Diag%dvsfc   (IM))
    allocate (Diag%dtsfc   (IM))
    allocate (Diag%dqsfc   (IM))
    allocate (Diag%totprcp (IM))
    allocate (Diag%totprcpb(IM))
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
    allocate (Diag%cnvprcpb(IM))
    allocate (Diag%spfhmin (IM))
    allocate (Diag%spfhmax (IM))
    allocate (Diag%u10mmax (IM))
    allocate (Diag%v10mmax (IM))
    allocate (Diag%wind10mmax (IM))
    allocate (Diag%u10max  (IM))
    allocate (Diag%v10max  (IM))
    allocate (Diag%spd10max (IM))
    allocate (Diag%rain    (IM))
    allocate (Diag%rainc   (IM))
    allocate (Diag%ice     (IM))
    allocate (Diag%snow    (IM))
    allocate (Diag%graupel (IM))
    allocate (Diag%totice  (IM))
    allocate (Diag%totsnw  (IM))
    allocate (Diag%totgrp  (IM))
    allocate (Diag%toticeb (IM))
    allocate (Diag%totsnwb (IM))
    allocate (Diag%totgrpb (IM))
    allocate (Diag%u10m    (IM))
    allocate (Diag%v10m    (IM))
    allocate (Diag%dpt2m   (IM))
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
#ifdef CCPP
    allocate (Diag%nswsfci (IM))
#endif
    allocate (Diag%uswsfci (IM))
    allocate (Diag%dusfci  (IM))
    allocate (Diag%dvsfci  (IM))
    allocate (Diag%dtsfci  (IM))
    allocate (Diag%dqsfci  (IM))
    allocate (Diag%gfluxi  (IM))
    allocate (Diag%epi     (IM))
    allocate (Diag%smcwlt2 (IM))
    allocate (Diag%smcref2 (IM))
    allocate (Diag%tdomr   (IM))
    allocate (Diag%tdomzr  (IM))
    allocate (Diag%tdomip  (IM))
    allocate (Diag%tdoms   (IM))
    allocate (Diag%skebu_wts(IM,Model%levs))
    allocate (Diag%skebv_wts(IM,Model%levs))
    allocate (Diag%sppt_wts(IM,Model%levs))
    allocate (Diag%shum_wts(IM,Model%levs))
!--- 3D diagnostics
    allocate (Diag%zmtnblck(IM))

    allocate (Diag%ca_out  (IM))
    allocate (Diag%ca_deep  (IM))
    allocate (Diag%ca_turb  (IM))
    allocate (Diag%ca_shal  (IM))
    allocate (Diag%ca_rad (IM))
    allocate (Diag%ca_micro  (IM))

#ifdef CCPP
    ! need to allocate these arrays to avoid Fortran runtime errors when
    ! trying to add slices of non-allocated fields to the CCPP data structure
#else
    if (Model%ldiag3d) then
#endif
      allocate (Diag%du3dt  (IM,Model%levs,4))
      allocate (Diag%dv3dt  (IM,Model%levs,4))
      allocate (Diag%dt3dt  (IM,Model%levs,7))
      ! DH*
      !if (Model%me==0) then
      !   write(0,*) "DH WARNING: TEMPORARY ALLOCATE Diag%dq3dt with size (IM,Model%levs,9) to avoid crash on MacOSX/GNU (and others?) in PROD mode"
      !end if
      allocate (Diag%dq3dt  (IM,Model%levs,9))
      ! *DH
!      allocate (Diag%dq3dt  (IM,Model%levs,oz_coeff+5))
!--- needed to allocate GoCart coupling fields
!      allocate (Diag%upd_mf (IM,Model%levs))
!      allocate (Diag%dwn_mf (IM,Model%levs))
!      allocate (Diag%det_mf (IM,Model%levs))
!      allocate (Diag%cldcov (IM,Model%levs))
#ifndef CCPP
    endif
#endif
    !--- 3D diagnostics for Thompson MP / GFDL MP
    allocate (Diag%refl_10cm(IM,Model%levs))

    !--  New max hourly diag.
    allocate (Diag%refdmax(IM))
    allocate (Diag%refdmax263k(IM))
    allocate (Diag%t02max(IM))
    allocate (Diag%t02min(IM))
    allocate (Diag%rh02max(IM))
    allocate (Diag%rh02min(IM))

#ifdef CCPP
    !--- MYNN variables:
    if (Model%do_mynnedmf) then
      !print*,"Allocating all MYNN-EDMF variables:"
      allocate (Diag%edmf_a    (IM,Model%levs))
      allocate (Diag%edmf_w    (IM,Model%levs))
      allocate (Diag%edmf_qt   (IM,Model%levs))
      allocate (Diag%edmf_thl  (IM,Model%levs))
      allocate (Diag%edmf_ent  (IM,Model%levs))
      allocate (Diag%edmf_qc   (IM,Model%levs))
      allocate (Diag%nupdraft  (IM))
      allocate (Diag%maxmf     (IM))
      allocate (Diag%ktop_shallow(IM))
      allocate (Diag%exch_h    (IM,Model%levs))
      allocate (Diag%exch_m    (IM,Model%levs))
      !print*,"Initializing all MYNN-EDMF variables with ",clear_val
      Diag%edmf_a        = clear_val
      Diag%edmf_w        = clear_val
      Diag%edmf_qt       = clear_val
      Diag%edmf_thl      = clear_val
      Diag%edmf_ent      = clear_val
      Diag%edmf_qc       = clear_val
      Diag%nupdraft      = 0
      Diag%maxmf         = clear_val
      Diag%ktop_shallow  = 0
      Diag%exch_h        = clear_val
      Diag%exch_m        = clear_val
    endif
#endif

    !--- diagnostics for coupled chemistry
    if (Model%cplchm) call Diag%chem_init(IM,Model)

    call Diag%rad_zero  (Model)
!    print *,'in diag_create, call phys_zero'
    linit = .true.
    call Diag%phys_zero (Model, linit=linit)
    linit = .false.

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
  subroutine diag_phys_zero (Diag, Model, linit)
    class(GFS_diag_type)               :: Diag
    type(GFS_control_type), intent(in) :: Model
    logical,optional, intent(in)       :: linit

    !--- In/Out
    Diag%srunoff    = zero
    Diag%evbsa      = zero
    Diag%evcwa      = zero
    Diag%snohfa     = zero
    Diag%transa     = zero
    Diag%sbsnoa     = zero
    Diag%snowca     = zero
#ifdef CCPP
    Diag%snowfallac = zero
    Diag%acsnow     = zero
#endif
    Diag%soilm      = zero
    Diag%tmpmin     = huge
    Diag%tmpmax     = zero
    Diag%dusfc      = zero
    Diag%dvsfc      = zero
    Diag%dtsfc      = zero
    Diag%dqsfc      = zero
    Diag%gflux      = zero
    Diag%dlwsfc     = zero
    Diag%ulwsfc     = zero
    Diag%suntim     = zero
    Diag%runoff     = zero
    Diag%ep         = zero
    Diag%cldwrk     = zero
    Diag%dugwd      = zero
    Diag%dvgwd      = zero
    Diag%psmean     = zero
    Diag%spfhmin    = huge
    Diag%spfhmax    = zero
    Diag%u10mmax    = zero
    Diag%v10mmax    = zero
    Diag%wind10mmax = zero
    Diag%u10max     = zero
    Diag%v10max     = zero
    Diag%spd10max   = zero
    Diag%rain       = zero
    Diag%rainc      = zero
    Diag%ice        = zero
    Diag%snow       = zero
    Diag%graupel    = zero

    !--- Out
    Diag%u10m       = zero
    Diag%v10m       = zero
    Diag%dpt2m      = zero
    Diag%zlvl       = zero
    Diag%psurf      = zero
    Diag%hpbl       = zero
    Diag%pwat       = zero
    Diag%t1         = zero
    Diag%q1         = zero
    Diag%u1         = zero
    Diag%v1         = zero
    Diag%chh        = zero
    Diag%cmm        = zero
    Diag%dlwsfci    = zero
    Diag%ulwsfci    = zero
    Diag%dswsfci    = zero
#ifdef CCPP
    Diag%nswsfci    = zero
#endif
    Diag%uswsfci    = zero
    Diag%dusfci     = zero
    Diag%dvsfci     = zero
    Diag%dtsfci     = zero
    Diag%dqsfci     = zero
    Diag%gfluxi     = zero
    Diag%epi        = zero
    Diag%smcwlt2    = zero
    Diag%smcref2    = zero
    Diag%tdomr      = zero
    Diag%tdomzr     = zero
    Diag%tdomip     = zero
    Diag%tdoms      = zero
    Diag%skebu_wts  = zero
    Diag%skebv_wts  = zero
    Diag%sppt_wts   = zero
    Diag%shum_wts   = zero
    Diag%zmtnblck   = zero
    Diag%totprcpb   = zero
    Diag%cnvprcpb   = zero
    Diag%toticeb    = zero
    Diag%totsnwb    = zero
    Diag%totgrpb    = zero
    Diag%ca_out     = zero
    Diag%ca_deep    = zero
    Diag%ca_turb    = zero
    Diag%ca_shal    = zero
    Diag%ca_rad     = zero
    Diag%ca_micro   = zero
!    if(Model%me == Model%master) print *,'in diag_phys_zero, totprcpb set to 0,kdt=',Model%kdt

    if (Model%ldiag3d) then
      Diag%du3dt   = zero
      Diag%dv3dt   = zero
      Diag%dt3dt   = zero
!      Diag%dq3dt   = zero
!      Diag%upd_mf  = zero
!      Diag%dwn_mf  = zero
!      Diag%det_mf  = zero
    endif

    Diag%refl_10cm = zero

! max hourly diagnostics
    Diag%refl_10cm = zero
    Diag%refdmax = -35.
    Diag%refdmax263k = -35.
    Diag%t02max = -999.
    Diag%t02min = 999.
    Diag%rh02max = -999.
    Diag%rh02min = 999.

    if (present(linit)) then
      if (linit) then
        Diag%totprcp = zero
        Diag%cnvprcp = zero
        Diag%totice  = zero
        Diag%totsnw  = zero
        Diag%totgrp  = zero
!       if(Model%me == Model%master) print *,'in diag_phys_zero, called in init step,set precip diag variable to zero',&
!                                            'size(Diag%totprcp)=',size(Diag%totprcp),'me=',Model%me,'kdt=',Model%kdt
      endif
    endif
  end subroutine diag_phys_zero

!-----------------------
! GFS_diag%chem_init
!-----------------------
  subroutine diag_chem_init(Diag, IM, Model)

    use parse_tracers,    only: get_tracer_index, NO_TRACER

    class(GFS_diag_type)               :: Diag
    integer,                intent(in) :: IM
    type(GFS_control_type), intent(in) :: Model

    ! -- local variables
    integer :: n

    ! -- initialize diagnostic variables depending on
    ! -- specific chemical tracers
    if (Model%ntchm > 0) then
      ! -- retrieve number of dust bins
      n = get_number_bins('dust')
      if (n > 0) then
        allocate (Diag%duem(IM,n))
        Diag%duem = zero
      end if

      ! -- retrieve number of sea salt bins
      n = get_number_bins('seas')
      if (n > 0) then
        allocate (Diag%ssem(IM,n))
        Diag%ssem = zero
      end if
    end if

    ! -- sedimentation and dry/wet deposition diagnostics
    if (associated(Model%ntdiag)) then
      ! -- get number of tracers with enabled diagnostics
      n = count(Model%ntdiag)

      ! -- initialize sedimentation
      allocate (Diag%sedim(IM,n))
      Diag%sedim = zero

      ! -- initialize dry deposition
      allocate (Diag%drydep(IM,n))
      Diag%drydep = zero

      ! -- initialize large-scale wet deposition
      allocate (Diag%wetdpl(IM,n))
      Diag%wetdpl = zero

      ! -- initialize convective-scale wet deposition
      allocate (Diag%wetdpc(IM,n))
      Diag%wetdpc = zero
    end if

    ! -- initialize anthropogenic and biomass
    ! -- burning emission diagnostics for
    ! -- (in order): black carbon,
    ! -- organic carbon, and sulfur dioxide
    allocate (Diag%abem(IM,6))
    Diag%abem = zero

    ! -- initialize column burden diagnostics
    ! -- for aerosol species (in order): pm2.5
    ! -- black carbon, organic carbon, sulfate,
    ! -- dust, sea salt
    allocate (Diag%aecm(IM,6))
    Diag%aecm = zero

  contains

    integer function get_number_bins(tracer_type)
      character(len=*), intent(in) :: tracer_type

      logical :: next
      integer :: n
      character(len=5) :: name

      get_number_bins = 0

      n = 0
      next = .true.
      do while (next)
        n = n + 1
        write(name,'(a,i1)') tracer_type, n + 1
        next = get_tracer_index(Model%tracer_names, name, &
          Model%me, Model%master, Model%debug) /= NO_TRACER
      end do

      get_number_bins = n

    end function get_number_bins

  end subroutine diag_chem_init

#ifdef CCPP
  !-------------------------
  ! GFS_interstitial_type%create
  !-------------------------
#ifdef HYBRID
  subroutine interstitial_create (Interstitial, IM, non_uniform_blocks, Model)
#else
  subroutine interstitial_create (Interstitial, IM, Model)
#endif
    !
    implicit none
    !
    class(GFS_interstitial_type)       :: Interstitial
    integer,                intent(in) :: IM
#ifdef HYBRID
    logical,                intent(in) :: non_uniform_blocks
#endif
    type(GFS_control_type), intent(in) :: Model
    !
    allocate (Interstitial%otspt      (Model%ntracp1,2))
    ! Set up numbers of tracers for PBL, convection, etc: sets
    ! Interstitial%{nncl,nvdiff,mg3_as_mg2,nn,tracers_total,ntiwx,ntk,ntkev,otspt,nsamftrac,ncstrac}
    call interstitial_setup_tracers(Interstitial, Model)
    ! Allocate arrays
    allocate (Interstitial%adjnirbmd  (IM))
    allocate (Interstitial%adjnirbmu  (IM))
    allocate (Interstitial%adjnirdfd  (IM))
    allocate (Interstitial%adjnirdfu  (IM))
    allocate (Interstitial%adjvisbmd  (IM))
    allocate (Interstitial%adjvisbmu  (IM))
    allocate (Interstitial%adjvisdfu  (IM))
    allocate (Interstitial%adjvisdfd  (IM))
    allocate (Interstitial%aerodp     (IM,NSPC1))
    allocate (Interstitial%alb1d      (IM))
    allocate (Interstitial%bexp1d     (IM))
    allocate (Interstitial%cd         (IM))
    allocate (Interstitial%cdq        (IM))
    allocate (Interstitial%cf_upi     (IM,Model%levs))
    allocate (Interstitial%clcn       (IM,Model%levs))
    allocate (Interstitial%cldf       (IM))
    allocate (Interstitial%cldsa      (IM,5))
    allocate (Interstitial%cldtaulw   (IM,Model%levr+LTP))
    allocate (Interstitial%cldtausw   (IM,Model%levr+LTP))
    allocate (Interstitial%cld1d      (IM))
    allocate (Interstitial%clouds     (IM,Model%levr+LTP,NF_CLDS))
    allocate (Interstitial%clw        (IM,Model%levs,Interstitial%nn))
    allocate (Interstitial%clx        (IM,4))
    allocate (Interstitial%cnv_dqldt  (IM,Model%levs))
    allocate (Interstitial%cnv_fice   (IM,Model%levs))
    allocate (Interstitial%cnv_mfd    (IM,Model%levs))
    allocate (Interstitial%cnv_ndrop  (IM,Model%levs))
    allocate (Interstitial%cnv_nice   (IM,Model%levs))
    allocate (Interstitial%cnvc       (IM,Model%levs))
    allocate (Interstitial%cnvw       (IM,Model%levs))
    allocate (Interstitial%ctei_r     (IM))
    allocate (Interstitial%ctei_rml   (IM))
    allocate (Interstitial%cumabs     (IM))
    allocate (Interstitial%dd_mf      (IM,Model%levs))
    allocate (Interstitial%de_lgth    (IM))
    allocate (Interstitial%del        (IM,Model%levs))
    allocate (Interstitial%del_gz     (IM,Model%levs+1))
    allocate (Interstitial%delr       (IM,Model%levr+LTP))
    allocate (Interstitial%dkt        (IM,Model%levs-1))
    allocate (Interstitial%dlength    (IM))
    allocate (Interstitial%dqdt       (IM,Model%levs,Model%ntrac))
    allocate (Interstitial%dqsfc1     (IM))
    allocate (Interstitial%drain      (IM))
    allocate (Interstitial%dtdt       (IM,Model%levs))
    allocate (Interstitial%dtdtc      (IM,Model%levs))
    allocate (Interstitial%dtsfc1     (IM))
    allocate (Interstitial%dt_mf      (IM,Model%levs))
    allocate (Interstitial%dtzm       (IM))
    allocate (Interstitial%dudt       (IM,Model%levs))
    allocate (Interstitial%dusfcg     (IM))
    allocate (Interstitial%dusfc1     (IM))
    allocate (Interstitial%dvdt       (IM,Model%levs))
    allocate (Interstitial%dvsfcg     (IM))
    allocate (Interstitial%dvsfc1     (IM))
    allocate (Interstitial%dvdftra    (IM,Model%levs,Interstitial%nvdiff))
    allocate (Interstitial%dzlyr      (IM,Model%levr+LTP))
    allocate (Interstitial%elvmax     (IM))
    allocate (Interstitial%ep1d       (IM))
    allocate (Interstitial%evap       (IM))
    allocate (Interstitial%evbs       (IM))
    allocate (Interstitial%evcw       (IM))
    allocate (Interstitial%faerlw     (IM,Model%levr+LTP,NBDLW,NF_AELW))
    allocate (Interstitial%faersw     (IM,Model%levr+LTP,NBDSW,NF_AESW))
    allocate (Interstitial%fh2        (IM))
    allocate (Interstitial%flag_cice  (IM))
    allocate (Interstitial%flag_guess (IM))
    allocate (Interstitial%flag_iter  (IM))
    allocate (Interstitial%fm10       (IM))
    allocate (Interstitial%frland     (IM))
    allocate (Interstitial%fscav      (Model%ntrac-Model%ncld+2))
    allocate (Interstitial%fswtr      (Model%ntrac-Model%ncld+2))
    allocate (Interstitial%gabsbdlw   (IM))
    allocate (Interstitial%gamma      (IM))
    allocate (Interstitial%gamq       (IM))
    allocate (Interstitial%gamt       (IM))
    allocate (Interstitial%gasvmr     (IM,Model%levr+LTP,NF_VGAS))
    allocate (Interstitial%gflx       (IM))
    allocate (Interstitial%gwdcu      (IM,Model%levs))
    allocate (Interstitial%gwdcv      (IM,Model%levs))
    allocate (Interstitial%h2o_pres   (levh2o))
    allocate (Interstitial%hflx       (IM))
    allocate (Interstitial%hprime1    (IM))
    allocate (Interstitial%idxday     (IM))
    allocate (Interstitial%islmsk     (IM))
    allocate (Interstitial%kbot       (IM))
    allocate (Interstitial%kcnv       (IM))
    allocate (Interstitial%kinver     (IM))
    allocate (Interstitial%kpbl       (IM))
    allocate (Interstitial%ktop       (IM))
    allocate (Interstitial%mbota      (IM,3))
    allocate (Interstitial%mtopa      (IM,3))
    allocate (Interstitial%oa4        (IM,4))
    allocate (Interstitial%oc         (IM))
    allocate (Interstitial%olyr       (IM,Model%levr+LTP))
    allocate (Interstitial%oz_pres    (levozp))
    allocate (Interstitial%plvl       (IM,Model%levr+1+LTP))
    allocate (Interstitial%plyr       (IM,Model%levr+LTP))
    allocate (Interstitial%prnum      (IM,Model%levs))
    allocate (Interstitial%qicn       (IM,Model%levs))
    allocate (Interstitial%qlcn       (IM,Model%levs))
    allocate (Interstitial%qlyr       (IM,Model%levr+LTP))
    allocate (Interstitial%prcpmp     (IM))
    allocate (Interstitial%qss        (IM))
    allocate (Interstitial%raincd     (IM))
    allocate (Interstitial%raincs     (IM))
    allocate (Interstitial%rainmcadj  (IM))
    allocate (Interstitial%rainp      (IM,Model%levs))
    allocate (Interstitial%rb         (IM))
    allocate (Interstitial%rhc        (IM,Model%levs))
    allocate (Interstitial%runoff     (IM))
    allocate (Interstitial%save_q     (IM,Model%levs,Model%ntrac))
    allocate (Interstitial%save_t     (IM,Model%levs))
    allocate (Interstitial%save_u     (IM,Model%levs))
    allocate (Interstitial%save_v     (IM,Model%levs))
    allocate (Interstitial%sbsno      (IM))
    allocate (Interstitial%scmpsw     (IM))
    allocate (Interstitial%sfcalb     (IM,NF_ALBD))
    allocate (Interstitial%sigma      (IM))
    allocate (Interstitial%sigmaf     (IM))
    allocate (Interstitial%sigmafrac  (IM,Model%levs))
    allocate (Interstitial%sigmatot   (IM,Model%levs))
    allocate (Interstitial%slopetype  (IM))
    allocate (Interstitial%snowc      (IM))
    allocate (Interstitial%snohf      (IM))
    allocate (Interstitial%snowmt     (IM))
    allocate (Interstitial%soiltype   (IM))
    allocate (Interstitial%stress     (IM))
    allocate (Interstitial%theta      (IM))
    allocate (Interstitial%tlvl       (IM,Model%levr+1+LTP))
    allocate (Interstitial%tlyr       (IM,Model%levr+LTP))
    allocate (Interstitial%trans      (IM))
    allocate (Interstitial%tseal      (IM))
    allocate (Interstitial%tsfa       (IM))
    allocate (Interstitial%tsfg       (IM))
    allocate (Interstitial%tsurf      (IM))
    allocate (Interstitial%ud_mf      (IM,Model%levs))
    allocate (Interstitial%ulwsfc_cice(IM))
    allocate (Interstitial%vdftra     (IM,Model%levs,Interstitial%nvdiff))  !GJF first dimension was set as 'IX' in GFS_physics_driver
    allocate (Interstitial%vegf1d     (IM))
    allocate (Interstitial%vegtype    (IM))
    allocate (Interstitial%w_upi      (IM,Model%levs))
    allocate (Interstitial%wcbmax     (IM))
    allocate (Interstitial%wind       (IM))
    allocate (Interstitial%work1      (IM))
    allocate (Interstitial%work2      (IM))
    allocate (Interstitial%work3      (IM))
    allocate (Interstitial%xcosz      (IM))
    allocate (Interstitial%xlai1d     (IM))
    allocate (Interstitial%xmu        (IM))
    allocate (Interstitial%z01d       (IM))
    allocate (Interstitial%zt1d       (IM))
    ! Allocate arrays that are conditional on physics choices
    if (Model%imp_physics == Model%imp_physics_gfdl .or. Model%imp_physics == Model%imp_physics_thompson) then
       allocate (Interstitial%graupelmp  (IM))
       allocate (Interstitial%icemp      (IM))
       allocate (Interstitial%rainmp     (IM))
       allocate (Interstitial%snowmp     (IM))
    else if (Model%imp_physics == Model%imp_physics_mg) then
       allocate (Interstitial%ncgl       (IM,Model%levs))
       allocate (Interstitial%ncpr       (IM,Model%levs))
       allocate (Interstitial%ncps       (IM,Model%levs))
       allocate (Interstitial%qgl        (IM,Model%levs))
       allocate (Interstitial%qrn        (IM,Model%levs))
       allocate (Interstitial%qsnw       (IM,Model%levs))
    end if
    if (Model%do_shoc) then
       if (.not. associated(Interstitial%qrn))  allocate (Interstitial%qrn  (IM,Model%levs))
       if (.not. associated(Interstitial%qsnw)) allocate (Interstitial%qsnw (IM,Model%levs))
       if (.not. associated(Interstitial%qgl))  allocate (Interstitial%qgl  (IM,Model%levs))
       allocate (Interstitial%ncpi (IM,Model%levs))
       allocate (Interstitial%ncpl (IM,Model%levs))
    end if
    ! Set components that do not change

    Interstitial%im               = IM
    Interstitial%ipr              = min(IM,10)
    Interstitial%ix               = IM
    Interstitial%latidxprnt       = 1
    Interstitial%levi             = Model%levs+1
    Interstitial%nsteps_per_reset = nint(Model%avg_max_length/Model%dtp)
    Interstitial%levh2o           = levh2o
    Interstitial%levozp           = levozp
    Interstitial%lm               = Model%levr
    Interstitial%lmk              = Model%levr+LTP
    Interstitial%lmp              = Model%levr+1+LTP
    Interstitial%h2o_coeff        = h2o_coeff
    Interstitial%oz_coeff         = oz_coeff
    ! h2o_pres and oz_pres do not change during the run, but
    ! need to be set later in GFS_phys_time_vary_init (after
    ! h2o_pres/oz_pres are read in read_h2odata/read_o3data)
    Interstitial%h2o_pres         = clear_val
    Interstitial%oz_pres          = clear_val
    !
    Interstitial%skip_macro       = .false.
    ! The value phys_hydrostatic from dynamics does not match the
    ! hardcoded value for calling GFDL MP in GFS_physics_driver.F90,
    ! which is set to .true.
    Interstitial%phys_hydrostatic = .true.
    !
#ifdef HYBRID
    Interstitial%non_uniform_blocks = non_uniform_blocks
#endif
    !
    ! Reset all other variables
    call Interstitial%rad_reset ()
    call Interstitial%phys_reset (Model)
    !
  end subroutine interstitial_create

  subroutine interstitial_setup_tracers(Interstitial, Model)
    !
    implicit none
    !
    class(GFS_interstitial_type)       :: Interstitial
    type(GFS_control_type), intent(in) :: Model
    integer :: n, tracers

    !first, initialize the values (in case the values don't get initialized within if statements below)
    Interstitial%nncl             = Model%ncld
    Interstitial%nvdiff           = Model%ntrac
    Interstitial%mg3_as_mg2       = .false.
    Interstitial%nn               = Model%ntrac + 1
    Interstitial%ntk              = 0
    Interstitial%ntkev            = 0
    Interstitial%tracers_total    = 0
    Interstitial%otspt(:,:)       = .true.
    Interstitial%nsamftrac        = 0
    Interstitial%ncstrac          = 0

    if (Model%imp_physics == Model%imp_physics_thompson) then
      if (Model%ltaerosol) then
        Interstitial%nvdiff = 8
      else
        Interstitial%nvdiff = 5
      endif
      if (Model%satmedmf) then
        Interstitial%nvdiff = Interstitial%nvdiff + 1
      endif
      Interstitial%nncl = 5
    elseif (Model%imp_physics == Model%imp_physics_wsm6) then
      Interstitial%nvdiff = Model%ntrac -3
      Interstitial%nncl = 5
    elseif (Model%ntclamt > 0) then             ! for GFDL MP don't diffuse cloud amount
      Interstitial%nvdiff = Model%ntrac - 1
    endif

    if (Model%imp_physics == Model%imp_physics_gfdl) then
      Interstitial%nncl = 5
    endif

    if (Model%imp_physics == Model%imp_physics_mg) then
      if (abs(Model%fprcp) == 1) then
        Interstitial%nncl = 4                          ! MG2 with rain and snow
        Interstitial%mg3_as_mg2 = .false.
      elseif (Model%fprcp >= 2) then
        if(Model%ntgl > 0 .and. (Model%mg_do_graupel .or. Model%mg_do_hail)) then
          Interstitial%nncl = 5                        ! MG3 with rain and snow and grapuel/hail
          Interstitial%mg3_as_mg2 = .false.
        else                              ! MG3 code run without graupel/hail i.e. as MG2
          Interstitial%nncl = 4
          Interstitial%mg3_as_mg2 = .true.
        endif
      endif
    endif

    if (Interstitial%nvdiff == Model%ntrac) then
      Interstitial%ntiwx = Model%ntiw
    else
      if (Model%imp_physics == Model%imp_physics_wsm6) then
        Interstitial%ntiwx = 3
      elseif (Model%imp_physics == Model%imp_physics_thompson) then
        if(Model%ltaerosol) then
          Interstitial%ntiwx = 3
        else
          Interstitial%ntiwx = 3
        endif
      elseif (Model%imp_physics == Model%imp_physics_gfdl) then
        Interstitial%ntiwx = 3
      else
        Interstitial%ntiwx = 0
      endif
    end if

    if (Model%imp_physics == Model%imp_physics_zhao_carr) then
      if (Model%cplchm) Interstitial%nvdiff = 3
    end if

    Interstitial%ntkev = Interstitial%nvdiff

    if (Model%ntiw > 0) then
      if (Model%ntclamt > 0) then
        Interstitial%nn = Model%ntrac - 2
      else
        Interstitial%nn = Model%ntrac - 1
      endif
    elseif (Model%ntcw > 0) then
      Interstitial%nn = Model%ntrac
    else
      Interstitial%nn = Model%ntrac + 1
    endif

    if (Model%cscnv .or. Model%satmedmf .or. Model%trans_trac ) then
      Interstitial%otspt(:,:)   = .true.     ! otspt is used only for cscnv
      Interstitial%otspt(1:3,:) = .false.    ! this is for sp.hum, ice and liquid water
      tracers = 2
      do n=2,Model%ntrac
        if ( n /= Model%ntcw  .and. n /= Model%ntiw  .and. n /= Model%ntclamt .and. &
             n /= Model%ntrw  .and. n /= Model%ntsw  .and. n /= Model%ntrnc   .and. &
             n /= Model%ntsnc .and. n /= Model%ntgl  .and. n /= Model%ntgnc) then
          tracers = tracers + 1
          if (Model%ntke  == n ) then
            Interstitial%otspt(tracers+1,1) = .false.
            Interstitial%ntk = tracers
          endif
          if (Model%ntlnc == n .or. Model%ntinc == n .or. Model%ntrnc == n .or. Model%ntsnc == n .or. Model%ntgnc == n)    &
!           if (ntlnc == n .or. ntinc == n .or. ntrnc == n .or. ntsnc == n .or.&
!               ntrw  == n .or. ntsw  == n .or. ntgl  == n)                    &
                  Interstitial%otspt(tracers+1,1) = .false.
        endif
      enddo
      Interstitial%tracers_total = tracers - 2
    endif   ! end if_ras or cfscnv or samf
    if(.not. Model%satmedmf .and. .not. Model%trans_trac) then
       Interstitial%nsamftrac = 0
    else
       Interstitial%nsamftrac = Interstitial%tracers_total
    endif
    Interstitial%ncstrac = Interstitial%tracers_total + 3

  end subroutine interstitial_setup_tracers

  subroutine interstitial_rad_reset (Interstitial)
    !
    implicit none
    !
    class(GFS_interstitial_type) :: Interstitial
    !
    Interstitial%aerodp       = clear_val
    Interstitial%alb1d        = clear_val
    Interstitial%cldsa        = clear_val
    Interstitial%cldtaulw     = clear_val
    Interstitial%cldtausw     = clear_val
    Interstitial%clouds       = clear_val
    Interstitial%de_lgth      = clear_val
    Interstitial%delr         = clear_val
    Interstitial%dzlyr        = clear_val
    Interstitial%faerlw       = clear_val
    Interstitial%faersw       = clear_val
    Interstitial%gasvmr       = clear_val
    Interstitial%idxday       = 0
    Interstitial%kb           = 0
    Interstitial%kd           = 0
    Interstitial%kt           = 0
    Interstitial%mbota        = 0
    Interstitial%mtopa        = 0
    Interstitial%nday         = 0
    Interstitial%olyr         = clear_val
    Interstitial%plvl         = clear_val
    Interstitial%plyr         = clear_val
    Interstitial%qlyr         = clear_val
    Interstitial%raddt        = clear_val
    Interstitial%scmpsw%uvbfc = clear_val
    Interstitial%scmpsw%uvbf0 = clear_val
    Interstitial%scmpsw%nirbm = clear_val
    Interstitial%scmpsw%nirdf = clear_val
    Interstitial%scmpsw%visbm = clear_val
    Interstitial%scmpsw%visdf = clear_val
    Interstitial%sfcalb       = clear_val
    Interstitial%tlvl         = clear_val
    Interstitial%tlyr         = clear_val
    Interstitial%tsfa         = clear_val
    Interstitial%tsfg         = clear_val
    !
  end subroutine interstitial_rad_reset

  subroutine interstitial_phys_reset (Interstitial, Model)
    !
    implicit none
    !
    class(GFS_interstitial_type) :: Interstitial
    type(GFS_control_type), intent(in) :: Model
    !
    Interstitial%adjnirbmd    = clear_val
    Interstitial%adjnirbmu    = clear_val
    Interstitial%adjnirdfd    = clear_val
    Interstitial%adjnirdfu    = clear_val
    Interstitial%adjvisbmd    = clear_val
    Interstitial%adjvisbmu    = clear_val
    Interstitial%adjvisdfu    = clear_val
    Interstitial%adjvisdfd    = clear_val
    Interstitial%bexp1d       = clear_val
    Interstitial%cd           = clear_val
    Interstitial%cdq          = clear_val
    Interstitial%cf_upi       = clear_val
    Interstitial%clcn         = clear_val
    Interstitial%cld1d        = clear_val
    Interstitial%cldf         = clear_val
    Interstitial%clw          = clear_val
    Interstitial%clw(:,:,2)   = -999.9
    Interstitial%clx          = clear_val
    Interstitial%cnv_dqldt    = clear_val
    Interstitial%cnv_fice     = clear_val
    Interstitial%cnv_mfd      = clear_val
    Interstitial%cnv_ndrop    = clear_val
    Interstitial%cnv_nice     = clear_val
    Interstitial%cnvc         = clear_val
    Interstitial%cnvw         = clear_val
    Interstitial%ctei_r       = clear_val
    Interstitial%ctei_rml     = clear_val
    Interstitial%cumabs       = clear_val
    Interstitial%dd_mf        = clear_val
    Interstitial%del          = clear_val
    Interstitial%del_gz       = clear_val
    Interstitial%dkt          = clear_val
    Interstitial%dlength      = clear_val
    Interstitial%dqdt         = clear_val
    Interstitial%dqsfc1       = clear_val
    Interstitial%drain        = clear_val
    Interstitial%dt_mf        = clear_val
    Interstitial%dtdt         = clear_val
    Interstitial%dtdtc        = clear_val
    Interstitial%dtsfc1       = clear_val
    Interstitial%dtzm         = clear_val
    Interstitial%dudt         = clear_val
    Interstitial%dusfcg       = clear_val
    Interstitial%dusfc1       = clear_val
    Interstitial%dvdftra      = clear_val
    Interstitial%dvdt         = clear_val
    Interstitial%dvsfcg       = clear_val
    Interstitial%dvsfc1       = clear_val
    Interstitial%elvmax       = clear_val
    Interstitial%ep1d         = clear_val
    Interstitial%evap         = clear_val
    Interstitial%evbs         = clear_val
    Interstitial%evcw         = clear_val
    Interstitial%fh2          = clear_val
    Interstitial%flag_cice    = .false.
    Interstitial%flag_guess   = .false.
    Interstitial%flag_iter    = .true.
    Interstitial%fm10         = clear_val
    Interstitial%frain        = clear_val
    Interstitial%frland       = clear_val
    Interstitial%fscav        = clear_val
    Interstitial%fswtr        = clear_val
    Interstitial%gabsbdlw     = clear_val
    Interstitial%gamma        = clear_val
    Interstitial%gamq         = clear_val
    Interstitial%gamt         = clear_val
    Interstitial%gflx         = clear_val
    Interstitial%gwdcu        = clear_val
    Interstitial%gwdcv        = clear_val
    Interstitial%hflx         = clear_val
    Interstitial%hprime1      = clear_val
    Interstitial%islmsk       = 0
    Interstitial%kbot         = Model%levs
    Interstitial%kcnv         = 0
    Interstitial%kinver       = Model%levs
    Interstitial%kpbl         = 0
    Interstitial%ktop         = 1
    Interstitial%oa4          = clear_val
    Interstitial%oc           = clear_val
    Interstitial%prcpmp       = clear_val
    Interstitial%prnum        = clear_val
    Interstitial%qicn         = clear_val
    Interstitial%qlcn         = clear_val
    Interstitial%qss          = clear_val
    Interstitial%raincd       = clear_val
    Interstitial%raincs       = clear_val
    Interstitial%rainmcadj    = clear_val
    Interstitial%rainp        = clear_val
    Interstitial%rb           = clear_val
    Interstitial%rhc          = clear_val
    Interstitial%rhcbot       = clear_val
    Interstitial%rhcpbl       = clear_val
    Interstitial%rhctop       = clear_val
    Interstitial%runoff       = clear_val
    Interstitial%save_q       = clear_val
    Interstitial%save_t       = clear_val
    Interstitial%save_u       = clear_val
    Interstitial%save_v       = clear_val
    Interstitial%sbsno        = clear_val
    Interstitial%sigma        = clear_val
    Interstitial%sigmaf       = clear_val
    Interstitial%sigmafrac    = clear_val
    Interstitial%sigmatot     = clear_val
    Interstitial%slopetype    = 0
    Interstitial%snowc        = clear_val
    Interstitial%snohf        = clear_val
    Interstitial%snowmt       = clear_val
    Interstitial%soiltype     = 0
    Interstitial%stress       = clear_val
    Interstitial%theta        = clear_val
    Interstitial%trans        = clear_val
    Interstitial%tseal        = clear_val
    Interstitial%tsurf        = clear_val
    Interstitial%ud_mf        = clear_val
    Interstitial%ulwsfc_cice  = clear_val
    Interstitial%vdftra       = clear_val
    Interstitial%vegf1d       = clear_val
    Interstitial%vegtype      = 0
    Interstitial%w_upi        = clear_val
    Interstitial%wcbmax       = clear_val
    Interstitial%wind         = clear_val
    Interstitial%work1        = clear_val
    Interstitial%work2        = clear_val
    Interstitial%work3        = clear_val
    Interstitial%xcosz        = clear_val
    Interstitial%xlai1d       = clear_val
    Interstitial%xmu          = clear_val
    Interstitial%z01d         = clear_val
    Interstitial%zt1d         = clear_val
    ! Reset fields that are conditional on physics choices
    if (Model%imp_physics == Model%imp_physics_gfdl .or. Model%imp_physics == Model%imp_physics_thompson) then
       Interstitial%graupelmp = clear_val
       Interstitial%icemp     = clear_val
       Interstitial%rainmp    = clear_val
       Interstitial%snowmp    = clear_val
    else if (Model%imp_physics == Model%imp_physics_mg) then
       Interstitial%ncgl      = clear_val
       Interstitial%ncpr      = clear_val
       Interstitial%ncps      = clear_val
       Interstitial%qgl       = clear_val
       Interstitial%qrn       = clear_val
       Interstitial%qsnw      = clear_val
    end if
    if (Model%do_shoc) then
       Interstitial%qrn       = clear_val
       Interstitial%qsnw      = clear_val
       Interstitial%qgl       = clear_val
       Interstitial%ncpi      = clear_val
       Interstitial%ncpl      = clear_val
    end if
    !
  end subroutine interstitial_phys_reset

  subroutine interstitial_print(Interstitial, Model, mpirank, omprank, blkno)
    !
    implicit none
    !
    class(GFS_interstitial_type) :: Interstitial
    type(GFS_control_type), intent(in) :: Model
    integer, intent(in) :: mpirank, omprank, blkno
    !
    ! Print static variables
    write (0,'(a,3i6)') 'Interstitial_print for mpirank, omprank, blkno: ', mpirank, omprank, blkno
    write (0,*) 'Interstitial_print: values that do not change'
    write (0,*) 'Interstitial%h2o_coeff         = ', Interstitial%h2o_coeff
    write (0,*) 'sum(Interstitial%h2o_pres)     = ', sum(Interstitial%h2o_pres)
    write (0,*) 'Interstitial%im                = ', Interstitial%im
    write (0,*) 'Interstitial%ipr               = ', Interstitial%ipr
    write (0,*) 'Interstitial%ix                = ', Interstitial%ix
    write (0,*) 'Interstitial%latidxprnt        = ', Interstitial%latidxprnt
    write (0,*) 'Interstitial%levi              = ', Interstitial%levi
    write (0,*) 'Interstitial%levh2o            = ', Interstitial%levh2o
    write (0,*) 'Interstitial%levozp            = ', Interstitial%levozp
    write (0,*) 'Interstitial%lm                = ', Interstitial%lm
    write (0,*) 'Interstitial%lmk               = ', Interstitial%lmk
    write (0,*) 'Interstitial%lmp               = ', Interstitial%lmp
    write (0,*) 'Interstitial%nsamftrac         = ', Interstitial%nsamftrac
    write (0,*) 'Interstitial%nsteps_per_reset  = ', Interstitial%nsteps_per_reset
    write (0,*) 'Interstitial%ntiwx             = ', Interstitial%ntiwx
    write (0,*) 'Interstitial%nvdiff            = ', Interstitial%nvdiff
    write (0,*) 'Interstitial%oz_coeff          = ', Interstitial%oz_coeff
    write (0,*) 'sum(Interstitial%oz_pres)      = ', sum(Interstitial%oz_pres)
    write (0,*) 'Interstitial%phys_hydrostatic  = ', Interstitial%phys_hydrostatic
    write (0,*) 'Interstitial%skip_macro        = ', Interstitial%skip_macro
#ifdef HYBRID
    write (0,*) 'Interstitial%non_uniform_blocks= ', Interstitial%non_uniform_blocks
#endif
    ! Print all other variables
    write (0,*) 'Interstitial_print: values that change'
    write (0,*) 'sum(Interstitial%adjnirbmd   ) = ', sum(Interstitial%adjnirbmd   )
    write (0,*) 'sum(Interstitial%adjnirbmu   ) = ', sum(Interstitial%adjnirbmu   )
    write (0,*) 'sum(Interstitial%adjnirdfd   ) = ', sum(Interstitial%adjnirdfd   )
    write (0,*) 'sum(Interstitial%adjnirdfu   ) = ', sum(Interstitial%adjnirdfu   )
    write (0,*) 'sum(Interstitial%adjvisbmd   ) = ', sum(Interstitial%adjvisbmd   )
    write (0,*) 'sum(Interstitial%adjvisbmu   ) = ', sum(Interstitial%adjvisbmu   )
    write (0,*) 'sum(Interstitial%adjvisdfu   ) = ', sum(Interstitial%adjvisdfu   )
    write (0,*) 'sum(Interstitial%adjvisdfd   ) = ', sum(Interstitial%adjvisdfd   )
    write (0,*) 'sum(Interstitial%aerodp      ) = ', sum(Interstitial%aerodp      )
    write (0,*) 'sum(Interstitial%alb1d       ) = ', sum(Interstitial%alb1d       )
    write (0,*) 'sum(Interstitial%bexp1d      ) = ', sum(Interstitial%bexp1d      )
    write (0,*) 'sum(Interstitial%cd          ) = ', sum(Interstitial%cd          )
    write (0,*) 'sum(Interstitial%cdq         ) = ', sum(Interstitial%cdq         )
    write (0,*) 'sum(Interstitial%cf_upi      ) = ', sum(Interstitial%cf_upi      )
    write (0,*) 'sum(Interstitial%clcn        ) = ', sum(Interstitial%clcn        )
    write (0,*) 'sum(Interstitial%cldf        ) = ', sum(Interstitial%cldf        )
    write (0,*) 'sum(Interstitial%cldsa       ) = ', sum(Interstitial%cldsa       )
    write (0,*) 'sum(Interstitial%cldtaulw    ) = ', sum(Interstitial%cldtaulw    )
    write (0,*) 'sum(Interstitial%cldtausw    ) = ', sum(Interstitial%cldtausw    )
    write (0,*) 'sum(Interstitial%cld1d       ) = ', sum(Interstitial%cld1d       )
    write (0,*) 'sum(Interstitial%clw         ) = ', sum(Interstitial%clw         )
    write (0,*) 'sum(Interstitial%clx         ) = ', sum(Interstitial%clx         )
    write (0,*) 'sum(Interstitial%clouds      ) = ', sum(Interstitial%clouds      )
    write (0,*) 'sum(Interstitial%cnv_dqldt   ) = ', sum(Interstitial%cnv_dqldt   )
    write (0,*) 'sum(Interstitial%cnv_fice    ) = ', sum(Interstitial%cnv_fice    )
    write (0,*) 'sum(Interstitial%cnv_mfd     ) = ', sum(Interstitial%cnv_mfd     )
    write (0,*) 'sum(Interstitial%cnv_ndrop   ) = ', sum(Interstitial%cnv_ndrop   )
    write (0,*) 'sum(Interstitial%cnv_nice    ) = ', sum(Interstitial%cnv_nice    )
    write (0,*) 'sum(Interstitial%cnvc        ) = ', sum(Interstitial%cnvc        )
    write (0,*) 'sum(Interstitial%cnvw        ) = ', sum(Interstitial%cnvw        )
    write (0,*) 'sum(Interstitial%ctei_r      ) = ', sum(Interstitial%ctei_r      )
    write (0,*) 'sum(Interstitial%ctei_rml    ) = ', sum(Interstitial%ctei_rml    )
    write (0,*) 'sum(Interstitial%cumabs      ) = ', sum(Interstitial%cumabs      )
    write (0,*) 'sum(Interstitial%dd_mf       ) = ', sum(Interstitial%dd_mf       )
    write (0,*) 'sum(Interstitial%de_lgth     ) = ', sum(Interstitial%de_lgth     )
    write (0,*) 'sum(Interstitial%del         ) = ', sum(Interstitial%del         )
    write (0,*) 'sum(Interstitial%del_gz      ) = ', sum(Interstitial%del_gz      )
    write (0,*) 'sum(Interstitial%delr        ) = ', sum(Interstitial%delr        )
    write (0,*) 'sum(Interstitial%dkt         ) = ', sum(Interstitial%dkt         )
    write (0,*) 'sum(Interstitial%dlength     ) = ', sum(Interstitial%dlength     )
    write (0,*) 'sum(Interstitial%dqdt        ) = ', sum(Interstitial%dqdt        )
    write (0,*) 'sum(Interstitial%dqsfc1      ) = ', sum(Interstitial%dqsfc1      )
    write (0,*) 'sum(Interstitial%drain       ) = ', sum(Interstitial%drain       )
    write (0,*) 'sum(Interstitial%dtdt        ) = ', sum(Interstitial%dtdt        )
    write (0,*) 'sum(Interstitial%dtdtc       ) = ', sum(Interstitial%dtdtc       )
    write (0,*) 'sum(Interstitial%dtsfc1      ) = ', sum(Interstitial%dtsfc1      )
    write (0,*) 'sum(Interstitial%dtzm        ) = ', sum(Interstitial%dtzm        )
    write (0,*) 'sum(Interstitial%dt_mf       ) = ', sum(Interstitial%dt_mf       )
    write (0,*) 'sum(Interstitial%dudt        ) = ', sum(Interstitial%dudt        )
    write (0,*) 'sum(Interstitial%dusfcg      ) = ', sum(Interstitial%dusfcg      )
    write (0,*) 'sum(Interstitial%dusfc1      ) = ', sum(Interstitial%dusfc1      )
    write (0,*) 'sum(Interstitial%dvdftra     ) = ', sum(Interstitial%dvdftra     )
    write (0,*) 'sum(Interstitial%dvdt        ) = ', sum(Interstitial%dvdt        )
    write (0,*) 'sum(Interstitial%dvsfcg      ) = ', sum(Interstitial%dvsfcg      )
    write (0,*) 'sum(Interstitial%dvsfc1      ) = ', sum(Interstitial%dvsfc1      )
    write (0,*) 'sum(Interstitial%dzlyr       ) = ', sum(Interstitial%dzlyr       )
    write (0,*) 'sum(Interstitial%elvmax      ) = ', sum(Interstitial%elvmax      )
    write (0,*) 'sum(Interstitial%ep1d        ) = ', sum(Interstitial%ep1d        )
    write (0,*) 'sum(Interstitial%evap        ) = ', sum(Interstitial%evap        )
    write (0,*) 'sum(Interstitial%evbs        ) = ', sum(Interstitial%evbs        )
    write (0,*) 'sum(Interstitial%evcw        ) = ', sum(Interstitial%evcw        )
    write (0,*) 'sum(Interstitial%faerlw      ) = ', sum(Interstitial%faerlw      )
    write (0,*) 'sum(Interstitial%faersw      ) = ', sum(Interstitial%faersw      )
    write (0,*) 'sum(Interstitial%fh2         ) = ', sum(Interstitial%fh2         )
    write (0,*) 'Interstitial%flag_cice(1)      = ', Interstitial%flag_cice(1)
    write (0,*) 'Interstitial%flag_guess(1)     = ', Interstitial%flag_guess(1)
    write (0,*) 'Interstitial%flag_iter(1)      = ', Interstitial%flag_iter(1)
    write (0,*) 'sum(Interstitial%fm10        ) = ', sum(Interstitial%fm10        )
    write (0,*) 'Interstitial%frain             = ', Interstitial%frain
    write (0,*) 'sum(Interstitial%frland      ) = ', sum(Interstitial%frland      )
    write (0,*) 'sum(Interstitial%fscav       ) = ', sum(Interstitial%fscav       )
    write (0,*) 'sum(Interstitial%fswtr       ) = ', sum(Interstitial%fswtr       )
    write (0,*) 'sum(Interstitial%gabsbdlw    ) = ', sum(Interstitial%gabsbdlw    )
    write (0,*) 'sum(Interstitial%gamma       ) = ', sum(Interstitial%gamma       )
    write (0,*) 'sum(Interstitial%gamq        ) = ', sum(Interstitial%gamq        )
    write (0,*) 'sum(Interstitial%gamt        ) = ', sum(Interstitial%gamt        )
    write (0,*) 'sum(Interstitial%gasvmr      ) = ', sum(Interstitial%gasvmr      )
    write (0,*) 'sum(Interstitial%gflx        ) = ', sum(Interstitial%gflx        )
    write (0,*) 'sum(Interstitial%gwdcu       ) = ', sum(Interstitial%gwdcu       )
    write (0,*) 'sum(Interstitial%gwdcv       ) = ', sum(Interstitial%gwdcv       )
    write (0,*) 'sum(Interstitial%hflx        ) = ', sum(Interstitial%hflx        )
    write (0,*) 'sum(Interstitial%hprime1     ) = ', sum(Interstitial%hprime1     )
    write (0,*) 'sum(Interstitial%idxday      ) = ', sum(Interstitial%idxday      )
    write (0,*) 'sum(Interstitial%islmsk      ) = ', sum(Interstitial%islmsk      )
    write (0,*) 'Interstitial%kb                = ', Interstitial%kb
    write (0,*) 'sum(Interstitial%kbot        ) = ', sum(Interstitial%kbot        )
    write (0,*) 'sum(Interstitial%kcnv        ) = ', sum(Interstitial%kcnv        )
    write (0,*) 'Interstitial%kd                = ', Interstitial%kd
    write (0,*) 'sum(Interstitial%kinver      ) = ', sum(Interstitial%kinver      )
    write (0,*) 'sum(Interstitial%kpbl        ) = ', sum(Interstitial%kpbl        )
    write (0,*) 'Interstitial%kt                = ', Interstitial%kt
    write (0,*) 'sum(Interstitial%ktop        ) = ', sum(Interstitial%ktop        )
    write (0,*) 'sum(Interstitial%mbota       ) = ', sum(Interstitial%mbota       )
    write (0,*) 'sum(Interstitial%mtopa       ) = ', sum(Interstitial%mtopa       )
    write (0,*) 'Interstitial%nday              = ', Interstitial%nday
    write (0,*) 'sum(Interstitial%oa4         ) = ', sum(Interstitial%oa4         )
    write (0,*) 'sum(Interstitial%oc          ) = ', sum(Interstitial%oc          )
    write (0,*) 'sum(Interstitial%olyr        ) = ', sum(Interstitial%olyr        )
    write (0,*) 'sum(Interstitial%plvl        ) = ', sum(Interstitial%plvl        )
    write (0,*) 'sum(Interstitial%plyr        ) = ', sum(Interstitial%plyr        )
    write (0,*) 'sum(Interstitial%prcpmp      ) = ', sum(Interstitial%prcpmp      )
    write (0,*) 'sum(Interstitial%prnum       ) = ', sum(Interstitial%prnum       )
    write (0,*) 'sum(Interstitial%qicn        ) = ', sum(Interstitial%qicn        )
    write (0,*) 'sum(Interstitial%qlcn        ) = ', sum(Interstitial%qlcn        )
    write (0,*) 'sum(Interstitial%qlyr        ) = ', sum(Interstitial%qlyr        )
    write (0,*) 'sum(Interstitial%qss         ) = ', sum(Interstitial%qss         )
    write (0,*) 'Interstitial%raddt             = ', Interstitial%raddt
    write (0,*) 'sum(Interstitial%raincd      ) = ', sum(Interstitial%raincd      )
    write (0,*) 'sum(Interstitial%raincs      ) = ', sum(Interstitial%raincs      )
    write (0,*) 'sum(Interstitial%rainmcadj   ) = ', sum(Interstitial%rainmcadj   )
    write (0,*) 'sum(Interstitial%rainp       ) = ', sum(Interstitial%rainp       )
    write (0,*) 'sum(Interstitial%rb          ) = ', sum(Interstitial%rb          )
    write (0,*) 'sum(Interstitial%rhc         ) = ', sum(Interstitial%rhc         )
    write (0,*) 'Interstitial%rhcbot            = ', Interstitial%rhcbot
    write (0,*) 'Interstitial%rhcpbl            = ', Interstitial%rhcpbl
    write (0,*) 'Interstitial%rhctop            = ', Interstitial%rhctop
    write (0,*) 'sum(Interstitial%runoff      ) = ', sum(Interstitial%runoff      )
    write (0,*) 'sum(Interstitial%save_q      ) = ', sum(Interstitial%save_q      )
    write (0,*) 'sum(Interstitial%save_t      ) = ', sum(Interstitial%save_t      )
    write (0,*) 'sum(Interstitial%save_u      ) = ', sum(Interstitial%save_u      )
    write (0,*) 'sum(Interstitial%save_v      ) = ', sum(Interstitial%save_v      )
    write (0,*) 'sum(Interstitial%sbsno       ) = ', sum(Interstitial%sbsno       )
    write (0,*) 'sum(Interstitial%scmpsw%uvbfc) = ', sum(Interstitial%scmpsw%uvbfc)
    write (0,*) 'sum(Interstitial%scmpsw%uvbf0) = ', sum(Interstitial%scmpsw%uvbf0)
    write (0,*) 'sum(Interstitial%scmpsw%nirbm) = ', sum(Interstitial%scmpsw%nirbm)
    write (0,*) 'sum(Interstitial%scmpsw%nirdf) = ', sum(Interstitial%scmpsw%nirdf)
    write (0,*) 'sum(Interstitial%scmpsw%visbm) = ', sum(Interstitial%scmpsw%visbm)
    write (0,*) 'sum(Interstitial%scmpsw%visdf) = ', sum(Interstitial%scmpsw%visdf)
    write (0,*) 'sum(Interstitial%sfcalb      ) = ', sum(Interstitial%sfcalb      )
    write (0,*) 'sum(Interstitial%sigma       ) = ', sum(Interstitial%sigma       )
    write (0,*) 'sum(Interstitial%sigmaf      ) = ', sum(Interstitial%sigmaf      )
    write (0,*) 'sum(Interstitial%sigmafrac   ) = ', sum(Interstitial%sigmafrac   )
    write (0,*) 'sum(Interstitial%sigmatot    ) = ', sum(Interstitial%sigmatot    )
    write (0,*) 'sum(Interstitial%slopetype   ) = ', sum(Interstitial%slopetype   )
    write (0,*) 'sum(Interstitial%snowc       ) = ', sum(Interstitial%snowc       )
    write (0,*) 'sum(Interstitial%snohf       ) = ', sum(Interstitial%snohf       )
    write (0,*) 'sum(Interstitial%snowmt      ) = ', sum(Interstitial%snowmt      )
    write (0,*) 'sum(Interstitial%soiltype    ) = ', sum(Interstitial%soiltype    )
    write (0,*) 'sum(Interstitial%stress      ) = ', sum(Interstitial%stress      )
    write (0,*) 'sum(Interstitial%theta       ) = ', sum(Interstitial%theta       )
    write (0,*) 'sum(Interstitial%tlvl        ) = ', sum(Interstitial%tlvl        )
    write (0,*) 'sum(Interstitial%tlyr        ) = ', sum(Interstitial%tlyr        )
    write (0,*) 'sum(Interstitial%trans       ) = ', sum(Interstitial%trans       )
    write (0,*) 'sum(Interstitial%tseal       ) = ', sum(Interstitial%tseal       )
    write (0,*) 'sum(Interstitial%tsfa        ) = ', sum(Interstitial%tsfa        )
    write (0,*) 'sum(Interstitial%tsfg        ) = ', sum(Interstitial%tsfg        )
    write (0,*) 'sum(Interstitial%tsurf       ) = ', sum(Interstitial%tsurf       )
    write (0,*) 'sum(Interstitial%ud_mf       ) = ', sum(Interstitial%ud_mf       )
    write (0,*) 'sum(Interstitial%ulwsfc_cice ) = ', sum(Interstitial%ulwsfc_cice )
    write (0,*) 'sum(Interstitial%vdftra      ) = ', sum(Interstitial%vdftra      )
    write (0,*) 'sum(Interstitial%vegf1d      ) = ', sum(Interstitial%vegf1d      )
    write (0,*) 'sum(Interstitial%vegtype     ) = ', sum(Interstitial%vegtype     )
    write (0,*) 'sum(Interstitial%w_upi       ) = ', sum(Interstitial%w_upi       )
    write (0,*) 'sum(Interstitial%wcbmax      ) = ', sum(Interstitial%wcbmax      )
    write (0,*) 'sum(Interstitial%wind        ) = ', sum(Interstitial%wind        )
    write (0,*) 'sum(Interstitial%work1       ) = ', sum(Interstitial%work1       )
    write (0,*) 'sum(Interstitial%work2       ) = ', sum(Interstitial%work2       )
    write (0,*) 'sum(Interstitial%work3       ) = ', sum(Interstitial%work3       )
    write (0,*) 'sum(Interstitial%xcosz       ) = ', sum(Interstitial%xcosz       )
    write (0,*) 'sum(Interstitial%xlai1d      ) = ', sum(Interstitial%xlai1d      )
    write (0,*) 'sum(Interstitial%xmu         ) = ', sum(Interstitial%xmu         )
    write (0,*) 'sum(Interstitial%z01d        ) = ', sum(Interstitial%z01d        )
    write (0,*) 'sum(Interstitial%zt1d        ) = ', sum(Interstitial%zt1d        )
    ! Print arrays that are conditional on physics choices
    if (Model%imp_physics == Model%imp_physics_gfdl .or. Model%imp_physics == Model%imp_physics_thompson) then
       write (0,*) 'Interstitial_print: values specific to GFDL/Thompson microphysics'
       write (0,*) 'sum(Interstitial%graupelmp) = ', sum(Interstitial%graupelmp   )
       write (0,*) 'sum(Interstitial%icemp    ) = ', sum(Interstitial%icemp       )
       write (0,*) 'sum(Interstitial%rainmp   ) = ', sum(Interstitial%rainmp      )
       write (0,*) 'sum(Interstitial%snowmp   ) = ', sum(Interstitial%snowmp      )
    else if (Model%imp_physics == Model%imp_physics_mg) then
       write (0,*) 'Interstitial_print: values specific to MG microphysics'
       write (0,*) 'sum(Interstitial%ncgl     ) = ', sum(Interstitial%ncgl        )
       write (0,*) 'sum(Interstitial%ncpr     ) = ', sum(Interstitial%ncpr        )
       write (0,*) 'sum(Interstitial%ncps     ) = ', sum(Interstitial%ncps        )
       write (0,*) 'sum(Interstitial%qgl      ) = ', sum(Interstitial%qgl         )
       write (0,*) 'sum(Interstitial%qrn      ) = ', sum(Interstitial%qrn         )
       write (0,*) 'sum(Interstitial%qsnw     ) = ', sum(Interstitial%qsnw        )
    end if
    if (Model%do_shoc) then
       write (0,*) 'Interstitial_print: values specific to SHOC'
       write (0,*) 'sum(Interstitial%ncgl     ) = ', sum(Interstitial%ncgl        )
       write (0,*) 'sum(Interstitial%qrn      ) = ', sum(Interstitial%qrn         )
       write (0,*) 'sum(Interstitial%qsnw     ) = ', sum(Interstitial%qsnw        )
       write (0,*) 'sum(Interstitial%qgl      ) = ', sum(Interstitial%qgl         )
       write (0,*) 'sum(Interstitial%ncpi     ) = ', sum(Interstitial%ncpi        )
       write (0,*) 'sum(Interstitial%ncpl     ) = ', sum(Interstitial%ncpl        )
    end if
    write (0,*) 'Interstitial_print: end'
    !
  end subroutine interstitial_print
#endif

end module GFS_typedefs

!!  This file contains the RUC land surface scheme driver.

module lsm_ruc

        use machine,           only: kind_phys

        use namelist_soilveg_ruc
        use set_soilveg_ruc_mod,  only: set_soilveg_ruc
        use module_soil_pre
        use module_sf_ruclsm

        implicit none

        private

        public :: lsm_ruc_init, lsm_ruc_run, lsm_ruc_finalize

      contains

!! \section arg_table_lsm_ruc_init Argument Table
!! | local_name     | standard_name                                               | long_name                                               | units      | rank | type      |    kind   | intent | optional |
!! |----------------|-------------------------------------------------------------|---------------------------------------------------------|------------|------|-----------|-----------|--------|----------|
!! | me             | mpi_rank                                                    | current MPI-rank                                        | index      |    0 | integer   |           | in     | F        |
!! | isot           | soil_type_dataset_choice                                    | soil type dataset choice                                | index      |    0 | integer   |           | in     | F        |
!! | ivegsrc        | vegetation_type_dataset_choice                              | land use dataset choice                                 | index      |    0 | integer   |           | in     | F        |
!! | nlunit         | iounit_namelist                                             | fortran unit number for file opens                      | none       |    0 | integer   |           | in     | F        |
!! | errmsg         | ccpp_error_message                                          | error message for error handling in CCPP                | none       |    0 | character | len=*     | out    | F        |
!! | errflg         | ccpp_error_flag                                             | error flag for error handling in CCPP                   | flag       |    0 | integer   |           | out    | F        |
!!
      subroutine lsm_ruc_init (me, isot, ivegsrc, nlunit,  &
     &                         errmsg, errflg)

      implicit none

      integer,              intent(in)  :: me, isot, ivegsrc, nlunit
      character(len=*),     intent(out) :: errmsg
      integer,              intent(out) :: errflg

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0 

      !--- initialize soil vegetation
      call set_soilveg_ruc(me, isot, ivegsrc, nlunit)

      end subroutine lsm_ruc_init

!---
!! \section arg_table_lsm_ruc_finalize Argument Table
!! | local_name     | standard_name                                               | long_name                                  | units      | rank | type      |  kind     | intent | optional |
!! |----------------|-------------------------------------------------------------|--------------------------------------------|------------|------|-----------|-----------|--------|----------|
!! | errmsg         | ccpp_error_message                                          | error message for error handling in CCPP   | none       |    0 | character |  len=*    | out    | F        |
!! | errflg         | ccpp_error_flag                                             | error flag for error handling in CCPP      | flag       |    0 | integer   |           | out    | F        |
!!
      subroutine lsm_ruc_finalize (errmsg, errflg)

      implicit none

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

  end subroutine lsm_ruc_finalize

! ===================================================================== !
!  lsm_ruc_run:                                                         !
!  RUC Surface Model - WRF4.0 version                                   ! 
!  program history log:                                                 !
!    may  2018  -- tanya smirnova                                       !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                       size   !
!     im       - integer, horiz dimention and num of used pts      1    !
!     km       - integer, vertical soil layer dimension            9    !
!     ps       - real, surface pressure (pa)                       im   !
!     u1, v1   - real, u/v component of surface layer wind         im   !
!     t1       - real, surface layer mean temperature (k)          im   !
!     q1       - real, surface layer mean specific humidity        im   !
!     soiltyp  - integer, soil type (integer index)                im   !
!     vegtype  - integer, vegetation type (integer index)          im   !
!     sigmaf   - real, areal fractional cover of green vegetation  im   !
!     sfcemis  - real, sfc lw emissivity ( fraction )              im   !
!     dlwflx   - real, total sky sfc downward lw flux ( w/m**2 )   im   !
!     dswflx   - real, total sky sfc downward sw flux ( w/m**2 )   im   !
!     snet     - real, total sky sfc netsw flx into ground(w/m**2) im   !
!     delt     - real, time interval (second)                      1    !
!     tg3      - real, deep soil temperature (k)                   im   !
!     cm       - real, surface exchange coeff for momentum (m/s)   im   !
!     ch       - real, surface exchange coeff heat & moisture(m/s) im   !
!     prsl1    - real, sfc layer 1 mean pressure (pa)              im   !
!     prslki   - real, dimensionless exner function at layer 1     im   !
!     zf       - real, height of bottom layer (m)                  im   !
!     islmsk  - integer, sea/land/ice mask (=0/1/2)                im   !
!     slopetyp - integer, class of sfc slope (integer index)       im   !
!     shdmin   - real, min fractional coverage of green veg        im   !
!     shdmax   - real, max fractnl cover of green veg (not used)   im   !
!     snoalb   - real, upper bound on max albedo over deep snow    im   !
!     sfalb    - real, mean sfc diffused sw albedo with effect          !
!                      of snow (fractional)                        im   !
!     flag_iter- logical,                                          im   !
!     flag_guess-logical,                                          im   !
!     isot     - integer, sfc soil type data source zobler or statsgo   !
!     ivegsrc  - integer, sfc veg type data source umd or igbp          !
!     smois    - real, total soil moisture content (fractional)   im,km !
!                                                                       !
!  input/outputs:                                                       !
!     weasd    - real, water equivalent accumulated snow depth (mm) im  !
!     snwdph   - real, snow depth (water equiv) over land          im   !
!     tskin    - real, ground surface skin temperature ( k )       im   !
!     tprcp    - real, total precipitation                         im   !
!     srflag   - real, snow/rain flag for precipitation            im   !
!     sr       - real, mixed-phase precipitation fraction          im   !
!     tslb     - real, soil temp (k)                              im,km !
!     sh2o     - real, liquid soil moisture                       im,km !
!     canopy   - real, canopy moisture content (mm)                im   !
!     trans    - real, total plant transpiration (m/s)             im   !
!     tsurf    - real, surface skin temperature (after iteration)  im   !
!                                                                       !
!  outputs:                                                             !
!     sncovr1  - real, snow cover over land (fractional)           im   !
!     qsurf    - real, specific humidity at sfc                    im   !
!     gflux    - real, soil heat flux (w/m**2)                     im   !
!     drain    - real, subsurface runoff (m/s)                     im   !
!     evap     - real, latent heat flux in kg kg-1 m s-1           im   !
!     runof    - real, surface runoff (m/s)                        im   !
!     evbs     - real, direct soil evaporation (m/s)               im   !
!     evcw     - real, canopy water evaporation (m/s)              im   !
!     sbsno    - real, sublimation/deposit from snopack (m/s)      im   !
!     stm      - real, total soil column moisture content (m)      im   !
!     z0rl     - real, surface roughness                           im   !
!     wet1     - real, normalized soil wetness                     im   !
!                                                                       !
!  ====================    end of description    =====================  !

#if 0
!! \section arg_table_lsm_ruc_run Argument Table
!! | local_name      | standard_name                                                                | long_name                                                       | units         | rank | type      |    kind   | intent | optional |
!! |----------------|------------------------------------------------------------------------------|-----------------------------------------------------------------|---------------|------|-----------|-----------|--------|----------|
!! | delt           | time_step_for_dynamics                                                       | physics time step                                               | s             |    0 | real      | kind_phys | in     | F        |
!! | me             | mpi_rank                                                                     | current MPI-rank                                                | index         |    0 | integer   |           | in     | F        |
!! | kdt            | index_of_time_step                                                           | current number of time steps                                    | index         |    0 | integer   |           | in     | F        |
!! | iter           | ccpp_loop_counter                                                            | loop counter for subcycling loops in CCPP                       | index         |    0 | integer   |           | in     | F        |
!! | im             | horizontal_loop_extent                                                       | horizontal loop extent                                          | count         |    0 | integer   |           | in     | F        |
!! | nlev           | vertical_dimension                                                           | number of vertical levels                                       | count         |    0 | integer   |           | in     | F        |
!! | lsm_ruc        | flag_for_ruc_land_surface_scheme                                             | flag for RUC land surface model                                 | flag          |    0 | integer   |           | in     | F        |
!! | lsm            | flag_for_land_surface_scheme                                                 | flag for land surface model                                     | flag          |    0 | integer   |           | in     | F        |
!! | do_mynnsfclay  | do_mynnsfclay                                                                | flag to activate MYNN surface layer                             | flag          |    0 | logical   |           | in     | F        |
!! | lsoil_ruc      | soil_vertical_dimension_for_land_surface_model                               | number of soil layers internal to land surface model            | count         |    0 | integer   |           | in     | F        |
!! | lsoil          | soil_vertical_dimension                                                      | soil vertical layer dimension                                   | count         |    0 | integer   |           | in     | F        |
!! | zs             | depth_of_soil_levels_for_land_surface_model                                  | depth of soil levels for land surface model                     | m             |    1 | real      | kind_phys | inout  | F        |
!! | islmsk         | sea_land_ice_mask                                                            | landmask: sea/land/ice=0/1/2                                    | flag          |    1 | integer   |           | in     | F        |
!! | con_cp         | specific_heat_of_dry_air_at_constant_pressure                                | specific heat !of dry air at constant pressure                  | J kg-1 K-1    |    0 | real      | kind_phys | in     | F        |
!! | con_g          | gravitational_acceleration                                                   | gravitational acceleration                                      | m s-2         |    0 | real      | kind_phys | in     | F        |
!! | con_pi         | pi                                                                           | ratio of a circle's circumference to its diameter               | radians       |    0 | real      | kind_phys | in     | F        |
!! | con_rd         | gas_constant_dry_air                                                         | ideal gas constant for dry air                                  | J kg-1 K-1    |    0 | real      | kind_phys | in     | F        |
!! | con_rv         | gas_constant_water_vapor                                                     | ideal gas constant for water vapor                              | J kg-1 K-1    |    0 | real      | kind_phys | in     | F        |
!! | con_hvap       | latent_heat_of_vaporization_of_water_at_0C                                   | latent heat of vaporization/sublimation (hvap)                  | J kg-1        |    0 | real      | kind_phys | in     | F        |
!! | con_fvirt      | ratio_of_vapor_to_dry_air_gas_constants_minus_one                            | rv/rd - 1 (rv = ideal gas constant for water vapor)             | none          |    0 | real      | kind_phys | in     | F        |
!! | con_tice       | freezing_point_temperature_of_seawater                                       | freezing point temperature of seawater                          | K             |    0 | real      | kind_phys | in     | F        |
!! | land           | flag_nonzero_land_surface_fraction                                           | flag indicating presence of some land surface area fraction     | flag          |    1 | logical   |           | in     | F        |
!! | icy            | flag_nonzero_sea_ice_surface_fraction                                        | flag indicating presence of some sea ice surface area fraction  | flag          |    1 | logical   |           | in     | F        |
!! | rainnc         | lwe_thickness_of_explicit_rainfall_amount_from_previous_timestep             | explicit rainfall from previous timestep                        | m             |    1 | real      | kind_phys | in     | F        |
!! | rainc          | lwe_thickness_of_convective_precipitation_amount_from_previous_timestep      | convective_precipitation_amount from previous timestep          | m             |    1 | real      | kind_phys | in     | F        |
!! | ice            | lwe_thickness_of_ice_amount_from_previous_timestep                           | ice amount from previous timestep                               | m             |    1 | real      | kind_phys | in     | F        |
!! | snow           | lwe_thickness_of_snow_amount_from_previous_timestep                          | snow amount from previous timestep                              | m             |    1 | real      | kind_phys | in     | F        |
!! | graupel        | lwe_thickness_of_graupel_amount_from_previous_timestep                       | graupel amount from previous timestep                           | m             |    1 | real      | kind_phys | in     | F        |
!! | srflag         | flag_for_precipitation_type                                                  | snow/rain flag for precipitation                                | flag          |    1 | real      | kind_phys | in     | F        |
!! | sr             | ratio_of_snowfall_to_rainfall                                                | snow ratio: ratio of snow to total precipitation                | frac          |    1 | real      | kind_phys | in     | F        |
!! | rhosnf         | density_of_frozen_precipitation                                              | density of frozen precipitation                                 | kg m-3        |    1 | real      | kind_phys | out    | F        |
!! | zf             | height_above_ground_at_lowest_model_layer                                    | layer 1 height above ground (not MSL)                           | m             |    1 | real      | kind_phys | in     | F        |
!! | u1             | x_wind_at_lowest_model_layer                                                 | zonal wind at lowest model layer                                | m s-1         |    1 | real      | kind_phys | in     | F        |
!! | v1             | y_wind_at_lowest_model_layer                                                 | meridional wind at lowest model layer                           | m s-1         |    1 | real      | kind_phys | in     | F        |
!! | prsl1          | air_pressure_at_lowest_model_layer                                           | mean pressure at lowest model layer                             | Pa            |    1 | real      | kind_phys | in     | F        |
!! | prslki         | ratio_of_exner_function_between_midlayer_and_interface_at_lowest_model_layer | Exner function ratio bt midlayer and interface at 1st layer     | ratio         |    1 | real      | kind_phys | in     | F        |
!! | ddvel          | surface_wind_enhancement_due_to_convection                                   | surface wind enhancement due to convection                      | m s-1         |    1 | real      | kind_phys | in     | F        |
!! | t1             | air_temperature_at_lowest_model_layer                                        | mean temperature at lowest model layer                          | K             |    1 | real      | kind_phys | in     | F        |
!! | q1             | water_vapor_specific_humidity_at_lowest_model_layer                          | water vapor specific humidity at lowest model layer             | kg kg-1       |    1 | real      | kind_phys | in     | F        |
!! | qc             | cloud_condensed_water_mixing_ratio_at_lowest_model_layer | moist (dry+vapor, no condensates) mixing ratio of cloud water at lowest model layer | kg kg-1       |    1 | real      | kind_phys | in     | F        |
!! | dlwflx         | surface_downwelling_longwave_flux                                            | surface downwelling longwave flux at current time               | W m-2         |    1 | real      | kind_phys | in     | F        |
!! | dswsfc         | surface_downwelling_shortwave_flux                                           | surface downwelling shortwave flux at current time              | W m-2         |    1 | real      | kind_phys | in     | F        |
!! | snet           | surface_net_downwelling_shortwave_flux                                       | surface net downwelling shortwave flux at current time          | W m-2         |    1 | real      | kind_phys | in     | F        |
!! | sfcemis        | surface_longwave_emissivity                                                  | surface lw emissivity in fraction                               | frac          |    1 | real      | kind_phys | inout  | F        |
!! | wspd           | wind_speed_at_lowest_model_layer                                             | wind speed at lowest model level                                | m s-1         |    1 | real      | kind_phys | inout  | F        |
!! | wet1           | normalized_soil_wetness                                                      | normalized soil wetness                                         | frac          |    1 | real      | kind_phys | inout  | F        |
!! | chh_lnd        | surface_drag_mass_flux_for_heat_and_moisture_in_air_over_land                | thermal exchange coefficient over land                          | kg m-2 s-1    |    1 | real      | kind_phys | inout  | F        |
!! | chh_ice        | surface_drag_mass_flux_for_heat_and_moisture_in_air_over_ice                 | thermal exchange coefficient over ice                           | kg m-2 s-1    |    1 | real      | kind_phys | inout  | F        |
!! | cmm_lnd        | surface_drag_wind_speed_for_momentum_in_air_over_land                        | momentum exchange coefficient over land                         | m s-1         |    1 | real      | kind_phys | inout  | F        |
!! | cmm_ice        | surface_drag_wind_speed_for_momentum_in_air_over_ice                         | momentum exchange coefficient over ice                          | m s-1         |    1 | real      | kind_phys | inout  | F        |
!! | canopy         | canopy_water_amount                                                          | canopy water amount                                             | kg m-2        |    1 | real      | kind_phys | inout  | F        |
!! | sigmaf         | vegetation_area_fraction                                                     | areal fractional cover of green vegetation                      | frac          |    1 | real      | kind_phys | in     | F        |
!! | alvwf          | mean_vis_albedo_with_weak_cosz_dependency                                    | mean vis albedo with weak cosz dependency                       | frac          |    1 | real      | kind_phys | in     | F        |
!! | alnwf          | mean_nir_albedo_with_weak_cosz_dependency                                    | mean nir albedo with weak cosz dependency                       | frac          |    1 | real      | kind_phys | in     | F        |
!! | snoalb         | upper_bound_on_max_albedo_over_deep_snow                                     | maximum snow albedo                                             | frac          |    1 | real      | kind_phys | in     | F        |
!! | sfalb_lnd      | surface_diffused_shortwave_albedo_over_land                                  | mean surface diffused sw albedo over land                       | frac          |    1 | real      | kind_phys | inout  | F        |
!! | sfalb_ice      | surface_diffused_shortwave_albedo_over_ice                                   | mean surface diffused sw albedo over ice                        | frac          |    1 | real      | kind_phys | inout  | F        |
!! | tg3            | deep_soil_temperature                                                        | deep soil temperature                                           | K             |    1 | real      | kind_phys | in     | F        |
!! | smc            | volume_fraction_of_soil_moisture                                             | total soil moisture                                             | frac          |    2 | real      | kind_phys | inout  | F        |
!! | slc            | volume_fraction_of_unfrozen_soil_moisture                                    | liquid soil moisture                                            | frac          |    2 | real      | kind_phys | inout  | F        |
!! | stc            | soil_temperature                                                             | soil temperature                                                | K             |    2 | real      | kind_phys | inout  | F        |
!! | smcwlt2        | volume_fraction_of_condensed_water_in_soil_at_wilting_point                  | soil water fraction at wilting point                            | frac          |    1 | real      | kind_phys | inout  | F        |
!! | smcref2        | threshold_volume_fraction_of_condensed_water_in_soil                         | soil moisture threshold                                         | frac          |    1 | real      | kind_phys | inout  | F        | 
!! | vegtype        | vegetation_type_classification                                               | vegetation type at each grid cell                               | index         |    1 | integer   |           | in     | F        |
!! | soiltyp        | soil_type_classification                                                     | soil type at each grid cell                                     | index         |    1 | integer   |           | in     | F        |
!! | isot           | soil_type_dataset_choice                                                     | soil type dataset choice                                        | index         |    0 | integer   |           | in     | F        |
!! | ivegsrc        | vegetation_type_dataset_choice                                               | land use dataset choice                                         | index         |    0 | integer   |           | in     | F        |
!! | lndfrac        | land_area_fraction                                                           | fraction of horizontal grid area occupied by land               | frac          |    1 | real      | kind_phys | in     | F        |
!! | fice           | sea_ice_concentration                                                        | ice fraction over open water                                    | frac          |    1 | real      | kind_phys | in     | F        |
!! | tice           | sea_ice_temperature                                                          | sea ice surface skin temperature                                | K             |    1 | real      | kind_phys | inout  | F        |
!! | keepfr         | flag_for_frozen_soil_physics                                                 | flag for frozen soil physics (RUC)                              | flag          |    2 | real      | kind_phys | inout  | F        |
!! | smois          | volume_fraction_of_soil_moisture_for_land_surface_model                      | volumetric fraction of soil moisture for lsm                    | frac          |    2 | real      | kind_phys | inout  | F        |
!! | sh2o           | volume_fraction_of_unfrozen_soil_moisture_for_land_surface_model             | volume fraction of unfrozen soil moisture for lsm               | frac          |    2 | real      | kind_phys | inout  | F        |
!! | smfrkeep       | volume_fraction_of_frozen_soil_moisture_for_land_surface_model               | volume fraction of frozen soil moisture for lsm                 | frac          |    2 | real      | kind_phys | inout  | F        |
!! | tslb           | soil_temperature_for_land_surface_model                                      | soil temperature for land surface model                         | K             |    2 | real      | kind_phys | inout  | F        |
!! | tslb_lnd       | soil_temperature_over_land                                                   | soil temperature over land                                      | K             |    2 | real      | kind_phys | inout  | F        |
!! | tsice          | ice_temperature_uncoupled                                                    | ice temperature for uncoupled ice model                         | K             |    2 | real      | kind_phys | inout  | F        |
!! | stm            | soil_moisture_content                                                        | soil moisture content                                           | kg m-2        |    1 | real      | kind_phys | inout  | F        |
!! | acsnow         | accumulated_water_equivalent_of_frozen_precip                                | snow water equivalent of run-total frozen precip                | kg m-2        |    1 | real      | kind_phys | inout  | F        |
!! | snowfallac_lnd | total_accumulated_snowfall_over_land                                         | run-total snow accumulation on the land                         | kg m-2        |    1 | real      | kind_phys | inout  | F        |
!! | snowfallac_ice | total_accumulated_snowfall_over_ice                                          | run-total snow accumulation on the ice                          | kg m-2        |    1 | real      | kind_phys | inout  | F        |
!! | sncovr1_lnd    | surface_snow_area_fraction_over_land                                         | surface snow area fraction over land                            | frac          |    1 | real      | kind_phys | inout  | F        |
!! | sncovr1_ice    | surface_snow_area_fraction_over_ice                                          | surface snow area fraction iocer ice                            | frac          |    1 | real      | kind_phys | inout  | F        |
!! | snowc          | surface_snow_area_fraction                                                   | surface snow area fraction                                      | frac          |    1 | real      | kind_phys | inout  | F        |
!! | weasd_lnd      | water_equivalent_accumulated_snow_depth_over_land                            | water equiv of acc snow depth over land                         | mm            |    1 | real      | kind_phys | in     | F        |
!! | weasd_ice      | water_equivalent_accumulated_snow_depth_over_ice                             | water equiv of acc snow depth over ice                          | mm            |    1 | real      | kind_phys | in     | F        |
!! | snwdph_lnd     | surface_snow_thickness_water_equivalent_over_land                            | water equivalent snow depth over land                           | mm            |    1 | real      | kind_phys | in     | F        |
!! | snwdph_ice     | surface_snow_thickness_water_equivalent_over_ice                             | water equivalent snow depth over ice                            | mm            |    1 | real      | kind_phys | in     | F        |
!! | tsnow_lnd      | snow_temperature_bottom_first_layer_over_land                                | snow temperature at the bottom of first snow layer over land    | K             |    1 | real      | kind_phys | inout  | F        |
!! | tsnow_ice      | snow_temperature_bottom_first_layer_over_ice                                 | snow temperature at the bottom of first snow layer over ice     | K             |    1 | real      | kind_phys | inout  | F        |
!! | tskin_lnd      | surface_skin_temperature_over_land_interstitial                              | surface skin temperature over land  (temporary use as interstitial) | K         |    1 | real      | kind_phys | in     | F        |
!! | tskin_ice      | surface_skin_temperature_over_ice_interstitial                               | surface skin temperature over ice   (temporary use as interstitial) | K         |    1 | real      | kind_phys | in     | F        |
!! | tskin_ocn      | surface_skin_temperature_over_ocean_interstitial                             | surface skin temperature over ocean (temporary use as interstitial) | K         |    1 | real      | kind_phys | in     | F        |
!! | tsurf_lnd      | surface_skin_temperature_after_iteration_over_land                           | surface skin temperature after iteration over land              | K             |    1 | real      | kind_phys | in     | F        |
!! | tsurf_ice      | surface_skin_temperature_after_iteration_over_ice                            | surface skin temperature after iteration over ice               | K             |    1 | real      | kind_phys | in     | F        |
!! | z0rl_lnd       | surface_roughness_length_over_land_interstitial                              | surface roughness length over land  (temporary use as interstitial) | cm        |    1 | real      | kind_phys | inout  | F        |
!! | z0rl_ice       | surface_roughness_length_over_ice_interstitial                               | surface roughness length over ice   (temporary use as interstitial) | cm        |    1 | real      | kind_phys | inout  | F        |
!! | cm_lnd         | surface_drag_coefficient_for_momentum_in_air_over_land                       | surface exchange coeff for momentum over land                   | none          |    1 | real      | kind_phys | inout  | F        |
!! | cm_ice         | surface_drag_coefficient_for_momentum_in_air_over_ice                        | surface exchange coeff for momentum over ice                    | none          |    1 | real      | kind_phys | inout  | F        |
!! | ch_lnd         | surface_drag_coefficient_for_heat_and_moisture_in_air_over_land              | surface exchange coeff heat & moisture over land                | none          |    1 | real      | kind_phys | inout  | F        |
!! | ch_ice         | surface_drag_coefficient_for_heat_and_moisture_in_air_over_ice               | surface exchange coeff heat & moisture over ice                 | none          |    1 | real      | kind_phys | inout  | F        |
!! | ch_ocn         | surface_drag_coefficient_for_heat_and_moisture_in_air_over_ocean             | surface exchange coeff heat & moisture over ocean               | none          |    1 | real      | kind_phys | inout  | F        |
!! | qsurf_lnd      | surface_specific_humidity_over_land                                          | surface air saturation specific humidity over land              | kg kg-1       |    1 | real      | kind_phys | inout  | F        |
!! | qsurf_ice      | surface_specific_humidity_over_ice                                           | surface air saturation specific humidity over ice               | kg kg-1       |    1 | real      | kind_phys | inout  | F        |
!! | sfcqc_lnd      | cloud_condensed_water_mixing_ratio_at_surface_over_land                      | moist cloud water mixing ratio at surface over land             | kg kg-1       |    1 | real      | kind_phys | inout  | F        |
!! | sfcqc_ice      | cloud_condensed_water_mixing_ratio_at_surface_over_ice                       | moist cloud water mixing ratio at surface over ice              | kg kg-1       |    1 | real      | kind_phys | inout  | F        |
!! | sfcqv_lnd      | water_vapor_mixing_ratio_at_surface_over_land                                | water vapor mixing ratio at surface over land                   | kg kg-1       |    1 | real      | kind_phys | inout  | F        |
!! | sfcqv_ice      | water_vapor_mixing_ratio_at_surface_over_ice                                 | water vapor mixing ratio at surface over ice                    | kg kg-1       |    1 | real      | kind_phys | inout  | F        |
!! | sfcdew_lnd     | surface_condensation_mass_over_land                                          | surface condensation mass over land                             | kg m-2        |    1 | real      | kind_phys | inout  | F        |
!! | sfcdew_ice     | surface_condensation_mass_over_ice                                           | surface condensation mass over ice                              | kg m-2        |    1 | real      | kind_phys | inout  | F        |
!! | evap_lnd       | kinematic_surface_upward_latent_heat_flux_over_land                          | kinematic surface upward latent heat flux over land             | kg kg-1 m s-1 |    1 | real      | kind_phys | in     | F        |
!! | evap_ice       | kinematic_surface_upward_latent_heat_flux_over_ice                           | kinematic surface upward latent heat flux over ice              | kg kg-1 m s-1 |    1 | real      | kind_phys | in     | F        |
!! | hflx_lnd       | kinematic_surface_upward_sensible_heat_flux_over_land                        | kinematic surface upward sensible heat flux over land           | K m s-1       |    1 | real      | kind_phys | in     | F        |
!! | hflx_ice       | kinematic_surface_upward_sensible_heat_flux_over_ice                         | kinematic surface upward sensible heat flux over ice            | K m s-1       |    1 | real      | kind_phys | in     | F        |
!! | sbsno_lnd      | snow_deposition_sublimation_upward_latent_heat_flux_over_land                | latent heat flux from snow depo/subl over land                  | W m-2         |    1 | real      | kind_phys | out    | F        |
!! | sbsno_ice      | snow_deposition_sublimation_upward_latent_heat_flux_over_ice                 | latent heat flux from snow depo/subl over ice                   | W m-2         |    1 | real      | kind_phys | out    | F        |
!! | evbs           | soil_upward_latent_heat_flux                                                 | soil upward latent heat flux                                    | W m-2         |    1 | real      | kind_phys | out    | F        |
!! | evcw           | canopy_upward_latent_heat_flux                                               | canopy upward latent heat flux                                  | W m-2         |    1 | real      | kind_phys | out    | F        |
!! | trans          | transpiration_flux                                                           | total plant transpiration rate                                  | kg m-2 s-1    |    1 | real      | kind_phys | out    | F        |
!! | runof          | surface_runoff_flux                                                          | surface runoff flux                                             | g m-2 s-1     |    1 | real      | kind_phys | out    | F        |
!! | drain          | subsurface_runoff_flux                                                       | subsurface runoff flux                                          | g m-2 s-1     |    1 | real      | kind_phys | out    | F        |
!! | runoff         | total_runoff                                                                 | total water runoff                                              | kg m-2        |    1 | real      | kind_phys | inout  | F        |
!! | srunoff        | surface_runoff                                                               | surface water runoff (from lsm)                                 | kg m-2        |    1 | real      | kind_phys | inout  | F        |
!! | gflux_lnd      | upward_heat_flux_in_soil_over_land                                           | soil heat flux over land                                        | W m-2         |    1 | real      | kind_phys | out    | F        |
!! | gflux_ice      | upward_heat_flux_in_soil_over_ice                                            | ice heat flux                                                   | W m-2         |    1 | real      | kind_phys | out    | F        |
!! | shdmin         | minimum_vegetation_area_fraction                                             | min fractional coverage of green vegetation                     | frac          |    1 | real      | kind_phys | in     | F        |
!! | shdmax         | maximum_vegetation_area_fraction                                             | max fractional coverage of green vegetation                     | frac          |    1 | real      | kind_phys | in     | F        |
!! | flag_iter      | flag_for_iteration                                                           | flag for iteration                                              | flag          |    1 | logical   |           | in     | F        |
!! | flag_guess     | flag_for_guess_run                                                           | flag for guess run                                              | flag          |    1 | logical   |           | in     | F        |
!! | flag_init      | flag_for_first_time_step                                                     | flag signaling first time step for time integration loop        | flag          |    0 | logical   |           | in     | F        |
!! | flag_restart   | flag_for_restart                                                             | flag for restart (warmstart) or coldstart                       | flag          |    0 | logical   |           | in     | F        |
!! | errmsg         | ccpp_error_message                                                           | error message for error handling in CCPP                        | none          |    0 | character | len=*     | out    | F        |
!! | errflg         | ccpp_error_flag                                                              | error flag for error handling in CCPP                           | flag          |    0 | integer   |           | out    | F        |
!!
#endif

      subroutine lsm_ruc_run                                            &
! --- inputs
     &     ( iter, me, kdt, im, nlev, lsoil_ruc, lsoil, zs,             &
     &       u1, v1, t1, q1, qc, soiltyp, vegtype, sigmaf,              &
     &       sfcemis, dlwflx, dswsfc, snet, delt, tg3, cm_lnd, cm_ice,  &
     &       ch_lnd, ch_ice, ch_ocn, tskin_ocn,                         &
     &       prsl1, prslki, zf, ddvel, shdmin, shdmax, alvwf, alnwf,    &
     &       snoalb, flag_iter, flag_guess, isot, ivegsrc, fice,        &
     &       smc, stc, slc, lsm_ruc, lsm, islmsk, land, icy,            &
     &       smcwlt2, smcref2, wspd, do_mynnsfclay,                     &
! --- constants
     &       con_cp, con_rv, con_rd, con_g, con_pi, con_hvap,           &
     &       con_fvirt, con_tice,                                       &
! --- in/outs
     &       weasd_lnd, snwdph_lnd, tskin_lnd, sfalb_lnd,               &
     &       weasd_ice, snwdph_ice, tskin_ice, sfalb_ice,               &
! --- in
     &       rainnc, rainc, ice, snow, graupel,                         &
! --- in/outs
     &       srflag, sr, snowc,                                         &
     ! for land
     &       smois, tslb, tslb_lnd, sh2o, keepfr, smfrkeep,             & ! on RUC levels
     &       canopy, trans, tsurf_lnd, tsnow_lnd, z0rl_lnd,             &
     &       sfcqc_lnd, sfcdew_lnd, sfcqv_lnd, lndfrac,                 &
     ! for ice
     &       sfcqc_ice, sfcdew_ice, sfcqv_ice,                          &
     &       tice, tsice, tsurf_ice, tsnow_ice, z0rl_ice,               &
! --- outputs
     &       acsnow, rhosnf,                                            & 
     ! for land
     &       sncovr1_lnd, qsurf_lnd, gflux_lnd, evap_lnd, hflx_lnd,     &
     &       runof, runoff, srunoff, drain,                             &
     &       chh_lnd, cmm_lnd, evbs, evcw, sbsno_lnd, stm, wet1,        &
     &       snowfallac_lnd,                                            &
     ! for ice
     &       sncovr1_ice, qsurf_ice, gflux_ice, evap_ice, hflx_ice,     &
     &       chh_ice, cmm_ice, sbsno_ice, snowfallac_ice,               &
     !
     &       flag_init, flag_restart, errmsg, errflg                    &
     &     )

      implicit none

!  ---  constant parameters:
      real(kind=kind_phys), parameter :: rhoh2o  = 1000.0
      real(kind=kind_phys), parameter :: stbolt  = 5.670400e-8
      real(kind=kind_phys), parameter :: cimin   = 0.15 !  --- minimum ice concentration

!  ---  input:
      integer, intent(in) :: me
      integer, intent(in) :: im, nlev, iter, lsoil_ruc, lsoil, kdt, isot, ivegsrc
      integer, intent(in) :: lsm_ruc, lsm

      real (kind=kind_phys), dimension(im,lsoil), intent(in) :: smc,stc,slc

      real (kind=kind_phys), dimension(im), intent(in) :: u1, v1,        &
     &       t1, sigmaf, sfcemis, dlwflx, dswsfc, snet, tg3,             &
     &       cm_lnd, cm_ice,                                             &
     &       ch_lnd, ch_ocn, prsl1, prslki, ddvel, shdmin, shdmax,       &
     &       snoalb, alvwf, alnwf, zf, qc, q1, wspd
      real (kind=kind_phys), dimension(im), intent(inout) :: ch_ice

      integer, dimension(im), intent(in) :: islmsk
      real (kind=kind_phys),  intent(in) :: delt
      real (kind=kind_phys),  intent(in) :: con_tice ! =271.2 K
      real (kind=kind_phys),  intent(in) :: con_cp, con_rv, con_g,       &
                                            con_pi, con_rd,              &
                                            con_hvap, con_fvirt

      logical, dimension(im), intent(in) :: flag_iter, flag_guess, land, icy
      logical,                intent(in) :: do_mynnsfclay

!  ---  in/out:
      integer, dimension(im), intent(inout) :: soiltyp, vegtype
      real (kind=kind_phys), dimension(lsoil_ruc) :: dzs
      real (kind=kind_phys), dimension(lsoil_ruc), intent(inout   ) :: zs
      real (kind=kind_phys), dimension(im), intent(inout) :: srflag, sr, &
     &       canopy, trans, smcwlt2, smcref2,                            & 
     ! for land
     &       weasd_lnd, snwdph_lnd, tskin_lnd, tskin_ocn,                &
     &       tsurf_lnd, z0rl_lnd, tsnow_lnd, lndfrac,                    &
     &       sfcqc_lnd, sfcqv_lnd, sfcdew_lnd, sfalb_lnd,                &
     ! for ice
     &       weasd_ice, snwdph_ice, tskin_ice,                           &
     &       tsurf_ice, z0rl_ice, tsnow_ice,                             &
     &       sfcqc_ice, sfcqv_ice, sfcdew_ice, fice, tice, sfalb_ice

!  ---  in
      real (kind=kind_phys), dimension(im), intent(in) ::                &
     &       rainnc, rainc, ice, snow, graupel
!  ---  in/out:
!  --- on RUC levels
      real (kind=kind_phys), dimension(im,lsoil_ruc), intent(inout) ::   &
     &       smois, tslb, tslb_lnd, sh2o, keepfr, smfrkeep, tsice

!  ---  output:
      real (kind=kind_phys), dimension(im), intent(inout) ::  acsnow,    &
     &       rhosnf, runof, drain, runoff, srunoff, evbs, evcw,          &
     &       stm, wet1, snowc,                                           &
     ! for land
     &       sncovr1_lnd, qsurf_lnd, gflux_lnd, evap_lnd,                &
     &       hflx_lnd, cmm_lnd, chh_lnd, sbsno_lnd,                      &
     &       snowfallac_lnd,                                             &
     ! for ice
     &       sncovr1_ice, qsurf_ice, gflux_ice, evap_ice,                &
     &       hflx_ice, cmm_ice, chh_ice, sbsno_ice,                      &
     &       snowfallac_ice

      logical,          intent(in)  :: flag_init, flag_restart
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

!  ---  locals:
      real (kind=kind_phys), dimension(im) :: rch, rho,                 &
     &       q0, qs1, wind,                                             &
     &       tprcp_old, srflag_old, sr_old, canopy_old,                 &
     &       acsnow_old, wet1_old,                                      &
     ! for land
     &       weasd_lnd_old, snwdph_lnd_old, tskin_lnd_old,              &
     &       tsnow_lnd_old, snowfallac_lnd_old, sfalb_lnd_old,          &
     &       sfcqv_lnd_old, sfcqc_lnd_old, z0rl_lnd_old,                &
     &       sncovr1_lnd_old,                                           &
     ! for ice
     &       weasd_ice_old, snwdph_ice_old, tskin_ice_old,              &
     &       tsnow_ice_old, snowfallac_ice_old, sfalb_ice_old,          &
     &       sfcqv_ice_old, sfcqc_ice_old, z0rl_ice_old,                &
     &       sncovr1_ice_old


      real (kind=kind_phys), dimension(lsoil_ruc) :: et

      real (kind=kind_phys), dimension(im,lsoil_ruc,1) :: smsoil,       &
           slsoil, stsoil, smfrsoil, keepfrsoil, stsice

      real (kind=kind_phys), dimension(im,lsoil_ruc) :: smois_old,      &
     &       tslb_old, tslb_lnd_old, sh2o_old,                          &
     &       keepfr_old, smfrkeep_old, tsice_old

      real (kind=kind_phys),dimension (im,1,1)      ::                  &
     &     conflx2, sfcprs, sfctmp, q2, qcatm, rho2 
      real (kind=kind_phys),dimension (im,1)        ::                  &
     &     albbck_lnd, alb_lnd, chs_lnd, flhc_lnd, flqc_lnd,            &
     &     wet, wet_ice, smmax, cmc, drip,  ec, edir, ett,              &
     &     dew_lnd, lh_lnd, esnow_lnd, etp, qfx_lnd, acceta,            &
     &     ffrozp, lwdn, prcp, xland, xland_ocn, xice, xice_lnd,        &
     &     graupelncv, snowncv, rainncv, raincv,                        &
     &     solnet_lnd, sfcexc,                                          &
     &     runoff1, runoff2, acrunoff,                                  &
     &     sfcems_lnd, hfx_lnd, shdfac, shdmin1d, shdmax1d,             &
     &     sneqv_lnd, snoalb1d_lnd, snowh_lnd, snoh_lnd, tsnav_lnd,     &
     &     snomlt_lnd, sncovr_lnd, soilw, soilm, ssoil_lnd,             &
     &     soilt_lnd, tbot,                                             &
     &     xlai, swdn, z0_lnd, znt_lnd, rhosnfr, infiltr,               &
     &     precipfr, snfallac_lnd, acsn,                                &
     &     qsfc_lnd, qsg_lnd, qvg_lnd, qcg_lnd, soilt1_lnd, chklowq
     ! ice
      real (kind=kind_phys),dimension (im,1)        ::                  &
     &     albbck_ice, alb_ice, chs_ice, flhc_ice, flqc_ice,            &
     &     dew_ice, lh_ice, esnow_ice, qfx_ice,                         &
     &     solnet_ice, sfcems_ice, hfx_ice,                             &
     &     sneqv_ice, snoalb1d_ice, snowh_ice, snoh_ice, tsnav_ice,     &
     &     snomlt_ice, sncovr_ice, ssoil_ice, soilt_ice,                &
     &     z0_ice, znt_ice, snfallac_ice,                               &
     &     qsfc_ice, qsg_ice, qvg_ice, qcg_ice, soilt1_ice


      real (kind=kind_phys) :: xice_threshold
      real (kind=kind_phys) :: focean, qss_ocean, evap_ocean, hflx_ocean

      character(len=256) :: llanduse  ! Land-use dataset.  Valid values are :
                                      ! "USGS" (USGS 24/27 category dataset) and
                                      ! "MODIFIED_IGBP_MODIS_NOAH" (MODIS 20-category dataset)

      integer :: nscat, nlcat
      real (kind=kind_phys), dimension(:,:,:), allocatable :: landusef ! fractional landuse
      real (kind=kind_phys), dimension(:,:,:), allocatable :: soilctop ! fractional soil type

      integer :: nsoil, iswater, isice
      integer, dimension (1:im,1:1) :: stype_ocn, vtype_ocn
      integer, dimension (1:im,1:1) :: stype_lnd, vtype_lnd
      integer, dimension (1:im,1:1) :: stype_ice, vtype_ice
      integer :: ipr

! local
      integer :: ims,ime, its,ite, jms,jme, jts,jte, kms,kme, kts,kte
      integer :: l, k, i, j,  fractional_seaice

      logical :: flag(im), flag_ice_uncoupled(im)
      logical :: rdlai2d, myj, frpcpn
      logical :: debug_print
!
      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      ipr = 10

      debug_print=.false.

      chklowq = 1.

      do i  = 1, im ! i - horizontal loop
        ! - Set flag for ice points for uncoupled model (islmsk(i) == 4 when coupled to CICE)
        flag_ice_uncoupled(i) = (icy(i) .and. (islmsk(i) == 2))
      enddo


      if (isot == 1) then
        nscat = 19 ! stasgo
      else
        nscat = 9  ! zobler
      endif
      allocate(soilctop(im,nscat,1))

      if(ivegsrc == 1) then
        nlcat = 20  ! IGBP - "MODI-RUC"
      else
        nlcat = 13
      endif
      allocate(landusef(im,nlcat,1))

      if(debug_print) then
        print *,'RUC LSM run'
        print *,'noah soil temp',ipr,stc(ipr,:)
        print *,'noah soil mois',ipr,smc(ipr,:)
        print *,'soiltyp=',ipr,soiltyp(ipr)
        print *,'vegtype=',ipr,vegtype(ipr)
        print *,'kdt, iter =',kdt,iter
        print *,'flag_init =',flag_init
        print *,'flag_restart =',flag_restart
      endif

! RUC initialization
      if(flag_init .and. iter==1) then
        !print *,'RUC LSM initialization, kdt=', kdt

        call rucinit          (flag_restart, im, lsoil_ruc, lsoil, nlev, & ! in
                               isot, soiltyp, vegtype, fice,             & ! in
                               land, flag_ice_uncoupled,                 & ! in
                               islmsk, tskin_lnd, tskin_ocn, tg3,        & ! in
                               smc, slc, stc,                            & ! in
                               smcref2, smcwlt2,                         & ! inout
                               lsm_ruc, lsm,                             & ! in
                               zs, sh2o, smfrkeep, tslb, smois, wet1,    & ! out
                               errmsg, errflg)

        do i  = 1, im ! n - horizontal loop
          do k = 1, lsoil_ruc
          ! - at initial time set sea ice T (tsice) and soil temperature (tslb_lnd)
          !   equal to TSLB, initialized from the Noah STC variable
             tslb_lnd(i,k) = tslb(i,k)
             tsice   (i,k) = tslb(i,k)
          enddo

          !tgs - Initialize land and ice surface albedo
          if(land(i)) then
            if (weasd_lnd(i) > 0.) then
              sfalb_lnd(i) = snoalb(i)
            else
              sfalb_lnd(i) = max(0.01, 0.5 * (alvwf(i) + alnwf(i)))
            endif
          elseif(flag_ice_uncoupled(i)) then
            if (weasd_ice(i) > 0.) then
              sfalb_ice(i) = 0.75
            else
              sfalb_ice(i) = 0.55
            endif
          endif
        enddo ! i

      endif ! flag_init=.true.,iter=1
!-- end of initialization

      ims = 1
      its = 1
      ime = 1
      ite = 1
      jms = 1
      jts = 1
      jme = 1
      jte = 1
      kms = 1
      kts = 1
      kme = 1
      kte = 1

      ! mosaic_lu=mosaic_soil=0, set in set_soilveg_ruc.F90
      ! set mosaic_lu=mosaic_soil=1 when fractional land and soil 
      ! categories available
      ! for now set fractions of differnet landuse and soil types 
      ! in the grid cell to zero

      landusef (:,:,:) = 0.0
      soilctop (:,:,:) = 0.0

      !> -- number of soil categories          
      !if(isot == 1) then
      !nscat = 19 ! stasgo
      !else
      !nscat = 9  ! zobler
      !endif
      !> -- set parameters for IGBP land-use data
      if(ivegsrc == 1) then
        llanduse = 'MODI-RUC'  ! IGBP
        iswater = 17
        isice = 15
      endif

      fractional_seaice = 1
      if ( fractional_seaice == 0 ) then
        xice_threshold = 0.5
      else if ( fractional_seaice == 1 ) then
        ! xice_threshold = 0.02
        xice_threshold = 0.15 ! consistent with GFS physics
      endif

      nsoil = lsoil_ruc

      do i  = 1, im ! i - horizontal loop
        ! reassign smcref2 and smcwlt2 to RUC values
        if(.not. land(i)) then
          !water and sea ice
          smcref2 (i) = 1.
          smcwlt2 (i) = 0.
        else
          !land 
          smcref2 (i) = REFSMC(soiltyp(i))
          smcwlt2 (i) = WLTSMC(soiltyp(i))
        endif
      enddo

      do i  = 1, im ! i - horizontal loop
        !> - Set flag for land and ice points.
        !- 10may19 - ice points are turned off.
        flag(i) = land(i) .or. flag_ice_uncoupled(i)
      enddo

      do i  = 1, im ! i - horizontal loop
        if (flag(i) .and. flag_guess(i)) then
          !> - Save land-related prognostic fields for guess run.
          !if(me==0 .and. i==ipr) print *,'before call to RUC guess run', i
          wet1_old(i)            = wet1(i)
          canopy_old(i)          = canopy(i)
          srflag_old(i)          = srflag(i)
          acsnow_old(i)          = acsnow(i)
          sr_old(i)              = sr(i)
          !tprcp_old(i)          = tprcp(i)
          ! for land
          weasd_lnd_old(i)       = weasd_lnd(i)
          snwdph_lnd_old(i)      = snwdph_lnd(i)
          tskin_lnd_old(i)       = tskin_lnd(i)
          tsnow_lnd_old(i)       = tsnow_lnd(i)
          snowfallac_lnd_old(i)  = snowfallac_lnd(i)
          sfalb_lnd_old(i)       = sfalb_lnd(i)
          sfcqv_lnd_old(i)       = sfcqv_lnd(i)
          sfcqc_lnd_old(i)       = sfcqc_lnd(i)
          z0rl_lnd_old(i)        = z0rl_lnd(i)
          sncovr1_lnd_old(i)     = sncovr1_lnd(i)
          ! for ice
          weasd_ice_old(i)       = weasd_ice(i)
          snwdph_ice_old(i)      = snwdph_ice(i)
          tskin_ice_old(i)       = tskin_ice(i)
          tsnow_ice_old(i)       = tsnow_ice(i)
          snowfallac_ice_old(i)  = snowfallac_ice(i)
          sfalb_ice_old(i)       = sfalb_ice(i)
          sfcqv_ice_old(i)       = sfcqv_ice(i)
          sfcqc_ice_old(i)       = sfcqc_ice(i)
          z0rl_ice_old(i)        = z0rl_ice(i)
          sncovr1_ice_old(i)     = sncovr1_ice(i)

          do k = 1, lsoil_ruc
            smois_old(i,k)  = smois(i,k)
            tslb_old(i,k)   = tslb(i,k)
            sh2o_old(i,k)   = sh2o(i,k)
            keepfr_old(i,k) = keepfr(i,k)
            smfrkeep_old(i,k) = smfrkeep(i,k)
            tslb_lnd_old(i,k) = tslb_lnd(i,k)
            ! for ice
            tsice_old(i,k)   = tsice(i,k)
          enddo
        endif
      enddo

!  --- ...  initialization block

      do j  = 1, 1
      do i  = 1, im ! i - horizontal loop
        if (flag_iter(i) .and. flag(i)) then
          !if(me==0 .and. i==ipr) print *,'iter run', iter, i, flag_iter(i),flag_guess(i)
          evap_lnd(i)  = 0.0
          evap_ice(i)  = 0.0
          hflx_lnd (i)  = 0.0
          hflx_ice (i)  = 0.0
          gflux_lnd(i)  = 0.0
          gflux_ice(i)  = 0.0
          drain(i)  = 0.0
          canopy(i) = max(canopy(i), 0.0)
          sfcdew_lnd(i) = 0.0
          sfcdew_ice(i) = 0.0

          evbs (i)  = 0.0
          evcw (i)  = 0.0
          trans(i)  = 0.0
          sbsno_lnd(i)  = 0.0
          sbsno_ice(i)  = 0.0

          !local i,j arrays
          dew_lnd(i,j)      = 0.0
          dew_ice(i,j)      = 0.0
          soilm(i,j)        = 0.0
          smmax(i,j)        = 0.0
          hfx_lnd(i,j)      = 0.0
          hfx_ice(i,j)      = 0.0
          qfx_lnd(i,j)      = 0.0
          qfx_ice(i,j)      = 0.0
          lh_lnd(i,j)       = 0.0
          lh_ice(i,j)       = 0.0
          acsn(i,j)         = 0.0
          sfcexc(i,j)       = 0.0
          acceta(i,j)       = 0.0
          ssoil_lnd(i,j)    = 0.0
          ssoil_ice(i,j)    = 0.0
          snomlt_lnd(i,j)   = 0.0
          snomlt_ice(i,j)   = 0.0
          infiltr(i,j)      = 0.0
          runoff1(i,j)      = 0.0
          runoff2(i,j)      = 0.0
          acrunoff(i,j)     = 0.0
          snfallac_lnd(i,j) = 0.0
          snfallac_ice(i,j) = 0.0
          rhosnfr(i,j)      = 0.0
          precipfr(i,j)     = 0.0

        endif
      enddo ! i=1,im
      enddo

!  --- ...  initialize atm. forcing variables

      do i  = 1, im
        if (flag_iter(i) .and. flag(i)) then
          !if (do_mynnsfclay) then
          ! WARNING - used of wspd computed in MYNN sfc leads to massive cooling.
          !  wind(i) = wspd(i)
          !else
            wind(i) = max(sqrt( u1(i)*u1(i) + v1(i)*v1(i) )               &
                    + max(0.0, min(ddvel(i), 30.0)), 1.0)
          !endif
          q0(i)   = max(q1(i)/(1.-q1(i)), 1.e-8)   !* q1=specific humidity at level 1 (kg/kg)

          rho(i) = prsl1(i) / (con_rd*t1(i)*(1.0+con_fvirt*q0(i)))
          qs1(i) = rslf(prsl1(i),t1(i))  !* qs1=sat. mixing ratio at level 1 (kg/kg)
          q0 (i) = min(qs1(i), q0(i))
        endif
      enddo ! i

!> - Prepare variables to run RUC LSM: 
!!  -   1. configuration information (c):
!!----------------------------------------
!! ffrozp  - fraction of frozen precipitation
!! frpcpn  - .true. if mixed phase precipitation available
!! 1:im - horizontal_loop_extent
!! fice    - fraction of sea-ice in the grid cell
!! delt    - timestep (sec) (dt should not exceed 3600 secs) 
!! conflx2 - height (\f$m\f$) above ground of atmospheric forcing variables
!! lsoil_ruc - number of soil layers (= 6 or 9)
!! zs      - the depth of each soil level (\f$m\f$)

      ! DH* TODO - TEST FOR DIFFERENT PHYSICS AND SET ACCORDINGLY?
      frpcpn = .true.                 ! .true. if mixed phase precipitation available (Thompson)

      do j  = 1, 1    ! 1:1
      do i  = 1, im   ! i - horizontal loop
        if (flag_iter(i) .and. flag(i)) then

        if(.not.frpcpn) then ! no mixed-phase precipitation available
          if (srflag(i) == 1.0) then  ! snow phase
            ffrozp(i,j) = 1.0
          elseif (srflag(i) == 0.0) then  ! rain phase
            ffrozp(i,j) = 0.0
          endif
        else ! mixed-phase precipitation is available
            ffrozp(i,j) = sr(i)
        endif ! frpcpn

        !tgs - for now set rdlai2d to .false., WRF has LAI maps, and RUC LSM
        !      uses rdlai2d = .true.
        rdlai2d = .false.
        !if( .not. rdlai2d) xlai = lai_data(vtype)

        conflx2(i,1,j)  = zf(i) ! first atm. level above ground surface

!  -   2. forcing data (f):
!!---------------------------------------
!! sfcprs  - pressure at height zf above ground (pascals)
!! sfctmp  - air temperature (\f$K\f$) at height zf above ground
!! q2      - pressure at height zf above ground (pascals)
!! qcatm   - cloud water mising ration at height zf above ground (\f$kg !kg^{-1}\f$)
!! rho2    - air density at height zf above ground (pascals)

        sfcprs(i,1,j)  = prsl1(i)
        sfctmp(i,1,j)  = t1(i)
        q2(i,1,j)      = q0(i)
        qcatm(i,1,j)   = max(0., qc(i))
        rho2(i,1,j)    = rho(i)

!!---------------------------------------
!! lwdn    - lw dw radiation flux at surface (\f$W m^{-2}\f$)
!! swdn    - sw dw radiation flux at surface (\f$W m^{-2}\f$)
!! prcp    - time-step total precip (\f$kg m^{-2} \f$)
!! raincv  - time-step convective precip (\f$kg m^{-2} \f$)
!! rainncv - time-step non-convective precip (\f$kg m^{-2} \f$)
!! graupelncv - time-step graupel (\f$kg m^{-2} \f$)
!! snowncv - time-step snow (\f$kg m^{-2} \f$)
!! precipfr - time-step precipitation in solod form (\f$kg m^{-2} \f$)
!! shdfac  - areal fractional coverage of green vegetation (0.0-1.0)
!! shdmin  - minimum areal fractional coverage of green vegetation -> shdmin1d
!! shdmax  - maximum areal fractional coverage of green vegetation -> shdmax1d
!! tbot    - bottom soil temperature (local yearly-mean sfc air temp)

        lwdn(i,j)   = dlwflx(i)         !..downward lw flux at sfc in w/m2
        swdn(i,j)   = dswsfc(i)         !..downward sw flux at sfc in w/m2

        ! all precip input to RUC LSM is in [mm]
        !prcp(i,j)       = rhoh2o * tprcp(i)                   ! tprcp in [m] - convective plus explicit
        !raincv(i,j)     = rhoh2o * rainc(i)                   ! total time-step convective precip
        !rainncv(i,j)    = rhoh2o * max(rain(i)-rainc(i),0.0)  ! total time-step explicit precip 
        !graupelncv(i,j) = rhoh2o * graupel(i)
        !snowncv(i,j)    = rhoh2o * snow(i)
        prcp(i,j)       = rhoh2o * (rainc(i)+rainnc(i))        ! tprcp in [m] - convective plus explicit
        raincv(i,j)     = rhoh2o * rainc(i)                    ! total time-step convective precip
        rainncv(i,j)    = rhoh2o * rainnc(i)                   ! total time-step explicit precip 
        graupelncv(i,j) = rhoh2o * graupel(i)
        snowncv(i,j)    = rhoh2o * snow(i)
        ! ice not used
        ! precipfr(i,j)   = rainncv(i,j) * ffrozp(i,j)
        acsn(i,j)         = acsnow(i)

        ! --- units %
        shdfac(i,j) = sigmaf(i)*100.
        shdmin1d(i,j) = shdmin(i)*100.
        shdmax1d(i,j) = shdmax(i)*100.

        tbot(i,j) = tg3(i)

!  -   3. canopy/soil characteristics (s):
!!------------------------------------
!! vegtyp  - vegetation type (integer index)                   -> vtype
!! vtype_ocn, _lnd, _ice  - vegetation type index for ocean, ice and land
!! portion of a grid cell
!! soiltyp - soil type (integer index)                         -> stype
!! stype_ocn, _lnd, _ice  - soil type index for ocean, ice and land
!! portion of a grid cell
!! sfcems  -  surface emmisivity                                   -> sfcemis
!! 0.5*(alvwf + alnwf) - backround snow-free surface albedo (fraction)         -> albbck
!! snoalb  - upper bound on maximum albedo over deep snow          -> snoalb1d
!! sfalb  - surface albedo including snow effect (unitless fraction) -> alb

        if(ivegsrc == 1) then   ! IGBP - MODIS
            vtype_ocn(i,j) = 17 ! 17 - water (oceans and lakes) in MODIS
            stype_ocn(i,j) = 14
            xland_ocn(i,j) = 2. ! xland = 2 for water
            vtype_lnd(i,j) = vegtype(i)
            stype_lnd(i,j) = soiltyp(i)
            vtype_ice(i,j) = 15 ! MODIS
            if(isot == 0) then
              stype_ice(i,j) = 9  ! ZOBLER
            else
              stype_ice(i,j) = 16 ! STASGO
            endif
        !> - Prepare land/ice/water masks for RUC LSM
        !SLMSK0   - SEA(0),LAND(1),ICE(2) MASK
          !if(islmsk(i) == 0.) then
          if(.not.land(i) .and. .not. flag_ice_uncoupled(i)) then
            xice(i,j)  = 0.
          !elseif(islmsk(i) == 1.) then ! land
          elseif(land(i)) then ! some land
            xland(i,j) = 1.
            xice_lnd(i,j) = 0.
          elseif(flag_ice_uncoupled(i)) then  ! some ice
            xland(i,j) = 1.
            xice(i,j)  = 1. !fice(i)  ! fraction of sea-ice
          endif
        else
          print *,'MODIS landuse is not available'
        endif

! - Call RUC LSM lsmruc(). 
   if(.not.land(i) .and. .not. flag_ice_uncoupled(i)) return

   if (land(i)) then ! at least some land in the grid cell

!  -   4. history (state) variables (h):
!! ------------------------------
!! cmc        - canopy moisture content (\f$mm\f$)
!! soilt = tskin - ground/canopy/snowpack effective skin temperature (\f$K\f$)
!! soilt1 = snowpack temperature at the bottom of the 1st layer (\f$K\f$)
!! tslb(lsoil_ruc) - soil temp (\f$K\f$)                                    -> stsoil
!! smois(lsoil_ruc) - total soil moisture content (volumetric fraction)     -> smsoil
!! sh2o(lsoil_ruc) - unfrozen soil moisture content (volumetric fraction)   -> slsoil
!! smfrsoil(lsoil_ruc) - frozen soil moisture content (volumetric fraction) -> smfrsoil
!! keepfrflag(lsoil_ruc) - flag for frozen soil physics: 0. or 1.
!! wet        - soil moisture availability at surface
!! snowh      - actual snow depth (\f$m\f$)
!! sneqv      - liquid water-equivalent snow depth (\f$m\f$)
!! sncovr     - fraction of snow in the grid cell
!! ch         - surface exchange coefficient for heat (\f$m s^{-1}\f$)      -> chs
!! z0         - surface roughness (\f$m\f$)     -> z0rl(\f$cm\f$)
!! qsfc    - specific humidity at surface (\f$kg kg^{-1}\f$)
!! qvg     - water vapor mixing ratio at surface (\f$kg kg^{-1}\f$)
!! qsg     - saturated water vapor mixing ratio at surface (\f$kg kg^{-1}\f$)
!! qcg     - cloud water mixing ratio at surface (\f$kg kg^{-1}\f$)
!! solnet  - net sw radiation flux (dn-up) (\f$W m^{-2}\f$)

        solnet_lnd(i,j) = dswsfc(i)*(1.-sfalb_lnd(i)) !snet(i) !..net sw rad flx (dn-up) at sfc in w/m2
        qvg_lnd(i,j)    = sfcqv_lnd(i)
        qsfc_lnd(i,j)   = sfcqv_lnd(i)/(1.+sfcqv_lnd(i))
        qsg_lnd(i,j)    = rslf(prsl1(i),tsurf_lnd(i))
        qcg_lnd(i,j)    = sfcqc_lnd(i) 
        sfcems_lnd(i,j) = sfcemis(i)
        snoalb1d_lnd(i,j) = snoalb(i)
        albbck_lnd(i,j)   = max(0.01, 0.5 * (alvwf(i) + alnwf(i))) 
        ! sfalb takes into account snow on the ground
        alb_lnd(i,j)      = sfalb_lnd(i)
        cmc(i,j) = canopy(i)            !  [mm] 
        soilt_lnd(i,j) = tsurf_lnd(i)            ! clu_q2m_iter
        tsnav_lnd(i,j) = 0.5*(soilt_lnd(i,j) + soilt1_lnd(i,j)) - 273.15
        ! sanity check for snow temperature tsnow
        if (tsnow_lnd(i) > 0. .and. tsnow_lnd(i) < 273.15) then
          soilt1_lnd(i,j) = tsnow_lnd(i)
        else
          soilt1_lnd(i,j) = tsurf_lnd(i)
        endif
        do k = 1, lsoil_ruc
          smsoil  (i,k,j) = smois(i,k)
          slsoil  (i,k,j) = sh2o(i,k)
          stsoil  (i,k,j) = tslb_lnd(i,k)
          smfrsoil(i,k,j) = smfrkeep(i,k)
          keepfrsoil(i,k,j) = keepfr(i,k)
        enddo
         ! land
        if (wet1(i) > 0.) then
         wet(i,j) = wet1(i)
        else
         wet(i,j) = max(0.0001,smsoil(i,1,j)/0.3)
        endif

        chs_lnd(i,j)    = ch_lnd(i) * wind(i) ! compute conductance 
        flhc_lnd(i,j)   = chs_lnd(i,j) * rho(i) * con_cp ! * (1. + 0.84*q2(i,1,j))
        flqc_lnd(i,j)   = chs_lnd(i,j) * rho(i) * wet(i,j)
        ! for output
        cmm_lnd(i)      = cm_lnd(i) * wind(i)
        chh_lnd(i)      = chs_lnd(i,j) * rho(i)
        !
        snowh_lnd(i,j) = snwdph_lnd(i) * 0.001         ! convert from mm to m
        sneqv_lnd(i,j) = weasd_lnd(i)                  ! [mm]
        snfallac_lnd(i,j) = snowfallac_lnd(i)
        sncovr_lnd(i,j) = sncovr1_lnd(i)
        !> -- sanity checks on sneqv and snowh
        if (sneqv_lnd(i,j) /= 0.0 .and. snowh_lnd(i,j) == 0.0) then
          snowh_lnd(i,j) = 0.003 * sneqv_lnd(i,j) ! snow density ~300 kg m-3 
        endif

        if (snowh_lnd(i,j) /= 0.0 .and. sneqv_lnd(i,j) == 0.0) then
          sneqv_lnd(i,j) = 300. * snowh_lnd(i,j) ! snow density ~300 kg m-3 
        endif

        if (sneqv_lnd(i,j) > 0. .and. snowh_lnd(i,j) > 0.) then
          if(sneqv_lnd(i,j)/snowh_lnd(i,j) > 950.) then
            sneqv_lnd(i,j) = 300. * snowh_lnd(i,j)
          endif
        endif
        !  ---- ... outside sflx, roughness uses cm as unit
        z0_lnd(i,j)  = z0rl_lnd(i)/100.
        znt_lnd(i,j) = z0rl_lnd(i)/100.

        if(debug_print) then
          !if(me==0 .and. solnet(i,j) =',i,j,solnet(i,j)
            print *,'sfcems(i,j) =',i,j,sfcems_lnd(i,j)
            print *,'chklowq(i,j) =',i,j,chklowq(i,j)
            print *,'chs(i,j) =',i,j,chs_lnd(i,j)
            print *,'flqc(i,j) =',i,j,flqc_lnd(i,j)
            print *,'flhc(i,j) =',i,j,flhc_lnd(i,j)
            print *,'wet(i,j) =',i,j,wet(i,j)
            print *,'cmc(i,j) =',i,j,cmc(i,j)
            print *,'shdfac(i,j) =',i,j,shdfac(i,j)
            print *,'alb(i,j) =',i,j,alb_lnd(i,j)
            print *,'znt(i,j) =',i,j,znt_lnd(i,j)
            print *,'z0(i,j) =',i,j,z0_lnd(i,j)
            print *,'snoalb1d(i,j) =',i,j,snoalb1d_lnd(i,j)
            print *,'alb(i,j) =',i,j,alb_lnd(i,j)
            print *,'landusef(i,:,j) =',i,j,landusef(i,:,j)
            print *,'soilctop(i,:,j) =',i,j,soilctop(i,:,j)
            print *,'nlcat=',nlcat
            print *,'nscat=',nscat
            print *,'qsfc(i,j) =',i,j,qsfc_lnd(i,j)
            print *,'qvg(i,j) =',i,j,qvg_lnd(i,j)
            print *,'qsg(i,j) =',i,j,qsg_lnd(i,j)
            print *,'qcg(i,j) =',i,j,qcg_lnd(i,j)
            print *,'dew(i,j) =',i,j,dew_lnd(i,j)
            print *,'soilt(i,j) =',i,j,soilt_lnd(i,j)
            print *,'tskin(i) =',i,j,tskin_lnd(i)
            print *,'soilt1(i,j) =',i,j,soilt1_lnd(i,j)
            print *,'tsnav(i,j) =',i,j,tsnav_lnd(i,j)
            print *,'tbot(i,j) =',i,j,tbot(i,j)
            print *,'vtype(i,j) =',i,j,vtype_lnd(i,j)
            print *,'stype(i,j) =',i,j,stype_lnd(i,j)
            print *,'xland(i,j) =',i,j,xland(i,j)
            print *,'xice(i,j) =',i,j,xice(i,j)
            print *,'iswater=',iswater
            print *,'isice=',isice
            print *,'xice_threshold=',xice_threshold
            print *,'con_cp=',con_cp
            print *,'con_rv=',con_rv
            print *,'con_rd=',con_rd
            print *,'con_g=',con_g
            print *,'con_pi=',con_pi
            print *,'con_hvap=',con_hvap
            print *,'stbolt=',stbolt
            print *,'smsoil(i,:,j)=',i,j,smsoil(i,:,j)
            print *,'slsoil(i,:,j)=',i,j,slsoil(i,:,j)
            print *,'stsoil(i,:,j)=',i,j,stsoil(i,:,j)
            print *,'smfrsoil(i,:,j)=',i,j,smfrsoil(i,:,j)
            print *,'keepfrsoil(i,:,j)=',i,j,keepfrsoil(i,:,j)
            print *,'soilm(i,j) =',i,j,soilm(i,j)
            print *,'smmax(i,j) =',i,j,smmax(i,j)
            print *,'hfx(i,j) =',i,j,hfx_lnd(i,j)
            print *,'qfx(i,j) =',i,j,qfx_lnd(i,j)
            print *,'lh(i,j) =',i,j,lh_lnd(i,j)
            print *,'infiltr(i,j) =',i,j,infiltr(i,j)
            print *,'runoff1(i,j) =',i,j,runoff1(i,j)
            print *,'runoff2(i,j) =',i,j,runoff2(i,j)
            print *,'acrunoff(i,j) =',i,j,acrunoff(i,j)
            print *,'sfcexc(i,j) =',i,j,sfcexc(i,j)
            print *,'acceta(i,j) =',i,j,acceta(i,j)
            print *,'ssoil(i,j) =',i,j,ssoil_lnd(i,j)
            print *,'snfallac(i,j) =',i,j,snfallac_lnd(i,j)
            print *,'acsn(i,j) =',i,j,acsn(i,j)
            print *,'snomlt(i,j) =',i,j,snomlt_lnd(i,j)
            print *,'shdmin1d(i,j) =',i,j,shdmin1d(i,j)
            print *,'shdmax1d(i,j) =',i,j,shdmax1d(i,j)
            print *,'rdlai2d =',rdlai2d
          !endif
        endif

      call lsmruc(                                                           &
     &          delt, flag_init, flag_restart, kdt, iter, nsoil,             &
     &          graupelncv(i,j), snowncv(i,j), rainncv(i,j), raincv(i,j),    &
     &          zs, prcp(i,j), sneqv_lnd(i,j), snowh_lnd(i,j),               &
     &          sncovr_lnd(i,j),                                             &
     &          ffrozp(i,j), frpcpn,                                         &
     &          rhosnfr(i,j), precipfr(i,j),                                 &
!  ---  inputs:
     &          conflx2(i,1,j), sfcprs(i,1,j), sfctmp(i,1,j), q2(i,1,j),     &
     &          qcatm(i,1,j), rho2(i,1,j),                                   &
     &          lwdn(i,j), solnet_lnd(i,j), sfcems_lnd(i,j), chklowq(i,j),   &
     &          chs_lnd(i,j), flqc_lnd(i,j), flhc_lnd(i,j),                  &
!  ---  input/outputs:
     &          wet(i,j), cmc(i,j), shdfac(i,j), alb_lnd(i,j), znt_lnd(i,j), &
     &          z0_lnd(i,j), snoalb1d_lnd(i,j), albbck_lnd(i,j),             &
     &          landusef(i,:,j), nlcat,                                      &
!  --- mosaic_lu and mosaic_soil are moved to the namelist
!     &          mosaic_lu, mosaic_soil,                                      &
     &          soilctop(i,:,j), nscat,                                      &
     &          qsfc_lnd(i,j), qsg_lnd(i,j), qvg_lnd(i,j), qcg_lnd(i,j),     &
     &          dew_lnd(i,j), soilt1_lnd(i,j),                               &
     &          tsnav_lnd(i,j), tbot(i,j), vtype_lnd(i,j), stype_lnd(i,j),   &
     &          xland(i,j), iswater, isice, xice_lnd(i,j), xice_threshold,   & !  xice=0. for the land portion of grid area
!  ---  constants
     &          con_cp, con_rv, con_rd, con_g, con_pi, con_hvap, stbolt,     &
!  ---  input/outputs:
     &          smsoil(i,:,j), slsoil(i,:,j), soilm(i,j), smmax(i,j),        &
     &          stsoil(i,:,j), soilt_lnd(i,j),                               &
     &          hfx_lnd(i,j), qfx_lnd(i,j), lh_lnd(i,j),                     &
     &          infiltr(i,j), runoff1(i,j), runoff2(i,j), acrunoff(i,j),     &
     &          sfcexc(i,j), acceta(i,j), ssoil_lnd(i,j),                    &
     &          snfallac_lnd(i,j), acsn(i,j), snomlt_lnd(i,j),               &
     &          smfrsoil(i,:,j),keepfrsoil(i,:,j), .false.,                  &
     &          shdmin1d(i,j), shdmax1d(i,j), rdlai2d,                       &
     &          ims,ime, jms,jme, kms,kme,                                   &
     &          its,ite, jts,jte, kts,kte                                    )
        if(debug_print) then
          print *,'after sneqv(i,j) =',i,j,sneqv_lnd(i,j)
          print *,'after snowh(i,j) =',i,j,snowh_lnd(i,j)
          print *,'after sncovr(i,j) =',i,j,sncovr_lnd(i,j)
          print *,'after vtype(i,j) =',i,j,vtype_lnd(i,j)
          print *,'after stype(i,j) =',i,j,stype_lnd(i,j)
          print *,'after wet(i,j) =',i,j,wet(i,j)
          print *,'after cmc(i,j) =',i,j,cmc(i,j)
          print *,'after qsfc(i,j) =',i,j,qsfc_lnd(i,j)
          print *,'after qvg(i,j) =',i,j,qvg_lnd(i,j)
          print *,'after qsg(i,j) =',i,j,qsg_lnd(i,j)
          print *,'after qcg(i,j) =',i,j,qcg_lnd(i,j)
          print *,'after dew(i,j) =',i,j,dew_lnd(i,j)
          print *,'after soilt(i,j) =',i,j,soilt_lnd(i,j)
          print *,'after tskin(i) =',i,j,tskin_lnd(i)
          print *,'after soilt1(i,j) =',i,j,soilt1_lnd(i,j)
          print *,'after tsnav(i,j) =',i,j,tsnav_lnd(i,j)
          print *,'after smsoil(i,:,j)=',i,j,smsoil(i,:,j)
          print *,'after slsoil(i,:,j)=',i,j,slsoil(i,:,j)
          print *,'after stsoil(i,:,j)=',i,j,stsoil(i,:,j)
          print *,'after smfrsoil(i,:,j)=',i,j,smfrsoil(i,:,j)
          print *,'after keepfrsoil(i,:,j)=',i,j,keepfrsoil(i,:,j)
          print *,'after soilm(i,j) =',i,j,soilm(i,j)
          print *,'after smmax(i,j) =',i,j,smmax(i,j)
          print *,'after hfx(i,j) =',i,j,hfx_lnd(i,j)
          print *,'after qfx(i,j) =',i,j,qfx_lnd(i,j)
          print *,'after lh(i,j) =',i,j,lh_lnd(i,j)
          print *,'after infiltr(i,j) =',i,j,infiltr(i,j)
          print *,'after runoff1(i,j) =',i,j,runoff1(i,j)
          print *,'after runoff2(i,j) =',i,j,runoff2(i,j)
          print *,'after ssoil(i,j) =',i,j,ssoil_lnd(i,j)
          print *,'after snfallac(i,j) =',i,j,snfallac_lnd(i,j)
          print *,'after acsn(i,j) =',i,j,acsn(i,j)
          print *,'after snomlt(i,j) =',i,j,snomlt_lnd(i,j)
        endif

!> - RUC LSM: prepare variables for return to parent model and unit conversion.
!>  -   6. output (o):
!! ------------------------------
!! lh     - actual latent heat flux (\f$W m^{-2}\f$: positive, if upward from sfc)
!! hfx    - sensible heat flux (\f$W m^{-2}\f$: positive, if upward from sfc)
!! ssoil   - soil heat flux (\f$W m^{-2}\f$: negative if downward from surface)
!! runoff1 - surface runoff (\f$m s^{-1}\f$), not infiltrating the surface
!! runoff2 - subsurface runoff (\f$m s^{-1}\f$), drainage out bottom
!! snoh    - phase-change heat flux from snowmelt (w m-2)
!
!  --- ...  do not return the following output fields to parent model
!    ec      - canopy water evaporation (m s-1)
!    edir    - direct soil evaporation (m s-1)
!    et(nsoil)-plant transpiration from a particular root layer (m s-1)
!    ett     - total plant transpiration (m s-1)
!    esnow   - sublimation from (or deposition to if <0) snowpack (m s-1)
!    drip    - through-fall of precip and/or dew in excess of canopy
!              water-holding capacity (m)
!    snomlt  - snow melt (m) (water equivalent)
!    xlai    - leaf area index (dimensionless)
!    soilw   - available soil moisture in root zone (unitless fraction
!              between smcwlt and smcmax)
!    soilm   - total soil column moisture content (frozen+unfrozen) (m)
!    nroot   - number of root layers, a function of veg type, determined
!              in subroutine redprm.


        !evbs(i)  = edir(i,j)
        !evcw(i)  = ec(i,j)
        !trans(i) = ett(i,j)
        !sbsno(i) = esnow(i,j)
        !snohf(i) = snoh(i,j)

        ! Interstitial
        evap_lnd(i)   = qfx_lnd(i,j) / rho(i)           ! kinematic
        hflx_lnd(i)   = hfx_lnd(i,j) / (con_cp*rho(i))  ! kinematic
        gflux_lnd(i)  = ssoil_lnd(i,j)
        sfcdew_lnd(i)  = dew_lnd(i,j)
        qsurf_lnd(i)   = qsfc_lnd(i,j)
        stm(i)         = soilm(i,j) * 1.e-3 ! convert to [m]
        tsurf_lnd(i)   = soilt_lnd(i,j)

        runof (i)  = runoff1(i,j)
        drain (i)  = runoff2(i,j)

        wet1(i) = wet(i,j)

        ! tsnow(i)   = soilt1(i,j)
        sfcqv_lnd(i)  = qvg_lnd(i,j)
        sfcqc_lnd(i)  = qcg_lnd(i,j)
        !  --- ...  units [m/s] = [g m-2 s-1] 
        rhosnf(i) = rhosnfr(i,j)
        acsnow(i) = acsn(i,j)     ! kg m-2

        ! --- ... accumulated total runoff and surface runoff
        runoff(i)  = runoff(i)  + (drain(i)+runof(i)) * delt * 0.001 ! kg m-2
        srunoff(i) = srunoff(i) + runof(i) * delt * 0.001            ! kg m-2

        ! --- ... accumulated frozen precipitation (accumulation in lsmruc)
        snowfallac_lnd(i) = snfallac_lnd(i,j) ! kg m-2
        !  --- ...  unit conversion (from m to mm)
        snwdph_lnd(i)  = snowh_lnd(i,j) * 1000.0

        canopy(i)      = cmc(i,j)   ! mm
        weasd_lnd(i)   = sneqv_lnd(i,j) ! mm
        sncovr1_lnd(i) = sncovr_lnd(i,j)
        !  ---- ... outside RUC LSM, roughness uses cm as unit 
        !  (update after snow's effect)
        z0rl_lnd(i) = znt_lnd(i,j)*100.
        sfalb_lnd(i)= alb_lnd(i,j)

        do k = 1, lsoil_ruc
          smois(i,k)  = smsoil(i,k,j)
          sh2o(i,k)   = slsoil(i,k,j)
          tslb_lnd(i,k) = stsoil(i,k,j)
          keepfr(i,k)   = keepfrsoil(i,k,j)
          smfrkeep(i,k) = smfrsoil(i,k,j)
        enddo
     if(debug_print) then
       print *,'LAND -i,j,stype_lnd,stype_ice,vtype_lnd,vtype_ice,',i,j,stype_lnd(i,j),stype_ice(i,j),vtype_lnd(i,j),vtype_ice(i,j)
       print *,'i,j,tsurf_lnd(i),tsurf_ice(i)',i,j,tsurf_lnd(i),tsurf_ice(i)
       print *,'kdt,iter,stsice(i,:,j),stsoil(i,:,j)',kdt,iter,stsice(i,:,j),stsoil(i,:,j)
     endif
   endif ! land

   if (flag_ice_uncoupled(i)) then ! at least some ice in the grid cell

!!\n  \a qsfc    - specific humidity at surface (\f$kg kg^{-1}\f$)
!!\n  \a qvg     - water vapor mixing ratio at surface (\f$kg kg^{-1}\f$)
!!\n  \a qsg     - saturated water vapor mixing ratio at surface (\f$kg kg^{-1}\f$)
!!\n  \a qcg     - cloud water mixing ratio at surface (\f$kg kg^{-1}\f$)
!!\n  \a solnet  - net sw radiation flux (dn-up) (\f$W m^{-2}\f$)
        solnet_ice(i,j) = dswsfc(i)*(1.-sfalb_ice(i))
        qvg_ice(i,j)    = sfcqv_ice(i)
        qsfc_ice(i,j)   = sfcqv_ice(i)/(1.+sfcqv_ice(i))
        qsg_ice(i,j)    = rslf(prsl1(i),tsurf_ice(i))
        qcg_ice(i,j)    = sfcqc_ice(i)
        sfcems_ice(i,j) = 0.98
        snoalb1d_ice(i,j) = 0.75 ! RAP value for max snow alb on ice
        albbck_ice(i,j)   = 0.55 ! RAP value for ice alb
        alb_ice(i,j)      = sfalb_ice(i)
        soilt_ice(i,j) = tsurf_ice(i)            ! clu_q2m_iter
        tsnav_ice(i,j) = 0.5*(soilt_ice(i,j) + soilt1_ice(i,j)) - 273.15
        if (tsnow_ice(i) > 0. .and. tsnow_ice(i) < 273.15) then
          soilt1_ice(i,j) = tsnow_ice(i)
        else
          soilt1_ice(i,j) = tsurf_ice(i)
        endif
        do k = 1, lsoil_ruc
          stsice  (i,k,j) = tsice(i,k)
        enddo
        wet_ice(i,j) = 1.

        chs_ice(i,j)    = ch_ice(i) * wind(i) ! compute conductance 
        flhc_ice(i,j)   = chs_ice(i,j) * rho(i) * con_cp ! * (1. + 0.84*q2(i,1,j))
        flqc_ice(i,j)   = chs_ice(i,j) * rho(i) * wet_ice(i,j) 
        cmm_ice(i)      = cm_ice(i) * wind(i)
        chh_ice(i)      = chs_ice(i,j) * rho(i)

        sncovr_ice(i,j) = sncovr1_ice(i)

        snowh_ice(i,j) = snwdph_ice(i) * 0.001         ! convert from mm to m
        sneqv_ice(i,j) = weasd_ice(i)                  ! [mm]
        snfallac_ice(i,j) = snowfallac_ice(i)

        !> -- sanity checks on sneqv and snowh
        if (sneqv_ice(i,j) /= 0.0 .and. snowh_ice(i,j) == 0.0) then
          snowh_ice(i,j) = 0.003 * sneqv_ice(i,j) ! snow density ~300 kg m-3 
        endif

        if (snowh_ice(i,j) /= 0.0 .and. sneqv_ice(i,j) == 0.0) then
          sneqv_ice(i,j) = 300. * snowh_ice(i,j) ! snow density ~300 kg m-3 
        endif

        if (sneqv_ice(i,j) > 0. .and. snowh_ice(i,j) > 0.) then
          if(sneqv_ice(i,j)/snowh_ice(i,j) > 950.) then
            sneqv_ice(i,j) = 300. * snowh_ice(i,j)
          endif
        endif

        z0_ice(i,j)  = z0rl_ice(i)/100.
        znt_ice(i,j) = z0rl_ice(i)/100.

      call lsmruc(                                                           &
     &          delt, flag_init, flag_restart, kdt, iter, nsoil,             &
     &          graupelncv(i,j), snowncv(i,j), rainncv(i,j), raincv(i,j),    &
     &          zs, prcp(i,j), sneqv_ice(i,j), snowh_ice(i,j),               &
     &          sncovr_ice(i,j),                                             &
     &          ffrozp(i,j), frpcpn,                                         &
     &          rhosnfr(i,j), precipfr(i,j),                                 &
!  ---  inputs:
     &          conflx2(i,1,j), sfcprs(i,1,j), sfctmp(i,1,j), q2(i,1,j),     &
     &          qcatm(i,1,j), rho2(i,1,j),                                   &
     &          lwdn(i,j), solnet_ice(i,j), sfcems_ice(i,j), chklowq(i,j),   &
     &          chs_ice(i,j), flqc_ice(i,j), flhc_ice(i,j),                  &
!  ---  input/outputs:
     &          wet_ice(i,j), cmc(i,j), shdfac(i,j), alb_ice(i,j),           &
     &          znt_ice(i,j), z0_ice(i,j), snoalb1d_ice(i,j),                &
     &          albbck_ice(i,j), landusef(i,:,j), nlcat,                     &
!  --- mosaic_lu and mosaic_soil are moved to the namelist
!     &          mosaic_lu, mosaic_soil,                                      &
     &          soilctop(i,:,j), nscat,                                      &
     &          qsfc_ice(i,j), qsg_ice(i,j), qvg_ice(i,j), qcg_ice(i,j),     &
     &          dew_ice(i,j), soilt1_ice(i,j),                               &
     &          tsnav_ice(i,j), tbot(i,j), vtype_ice(i,j), stype_ice(i,j),   &
     &          xland(i,j), iswater, isice, xice(i,j), xice_threshold,       &
!  ---  constants
     &          con_cp, con_rv, con_rd, con_g, con_pi, con_hvap, stbolt,     &
!  ---  input/outputs:
     &          smsoil(i,:,j), slsoil(i,:,j), soilm(i,j), smmax(i,j),        &
     &          stsice(i,:,j), soilt_ice(i,j),                               &
     &          hfx_ice(i,j), qfx_ice(i,j), lh_ice(i,j),                     &
     &          infiltr(i,j), runoff1(i,j), runoff2(i,j), acrunoff(i,j),     &
     &          sfcexc(i,j), acceta(i,j), ssoil_ice(i,j),                    &
     &          snfallac_ice(i,j), acsn(i,j), snomlt_ice(i,j),                 &
     &          smfrsoil(i,:,j),keepfrsoil(i,:,j), .false.,                  &
     &          shdmin1d(i,j), shdmax1d(i,j), rdlai2d,                       &
     &          ims,ime, jms,jme, kms,kme,                                   &
     &          its,ite, jts,jte, kts,kte                                    )

!> - RUC LSM: prepare variables for return to parent model and unit conversion.
!>  -   6. output (o):
!!\n  ------------------------------
!!\n \a lh     - actual latent heat flux (\f$W m^{-2}\f$: positive, if upward from sfc)
!!\n \a hfx    - sensible heat flux (\f$W m^{-2}\f$: positive, if upward from sfc)
!!\n \a ssoil   - soil heat flux (\f$W m^{-2}\f$: negative if downward from surface)
!!\n \a snoh    - phase-change heat flux from snowmelt (w m-2)
!

        ! Interstitial
        evap_ice(i)   = qfx_ice(i,j) / rho(i)           ! kinematic
        hflx_ice(i)   = hfx_ice(i,j) / (con_cp*rho(i))  ! kinematic
        gflux_ice(i)  = ssoil_ice(i,j)

        sfcdew_ice(i)  = dew_ice(i,j)
        qsurf_ice(i)   = qsfc_ice(i,j)
        tsurf_ice(i)   = soilt_ice(i,j)
        tice(i)        = tsurf_ice(i)

        sfcqv_ice(i)  = qvg_ice(i,j)
        sfcqc_ice(i)  = qcg_ice(i,j)

        !  --- ...  units [m/s] = [g m-2 s-1] 
        rhosnf(i) = rhosnfr(i,j)
        acsnow(i) = acsn(i,j)     ! kg m-2

        snowfallac_ice(i) = snfallac_ice(i,j) ! kg m-2
        !  --- ...  unit conversion (from m to mm)
        snwdph_ice(i)  = snowh_ice(i,j) * 1000.0
        weasd_ice(i)   = sneqv_ice(i,j) ! mm
        sncovr1_ice(i) = sncovr_ice(i,j)
        z0rl_ice(i) = znt_ice(i,j)*100.
        sfalb_ice(i)= alb_ice(i,j)

        do k = 1, lsoil_ruc
          tsice(i,k)  = stsice(i,k,j)
        enddo
     if(debug_print) then
       print *,'ICE - i,j,stype_lnd,stype_ice,vtype_lnd,vtype_ice,',i,j,stype_lnd(i,j),stype_ice(i,j),vtype_lnd(i,j),vtype_ice(i,j)
       print *,'i,j,tsurf_lnd(i),tsurf_ice(i)',i,j,tsurf_lnd(i),tsurf_ice(i)
       print *,'kdt,iter,stsice(i,:,j),stsoil(i,:,j)',kdt,iter,stsice(i,:,j),stsoil(i,:,j)
     endif

   endif ! ice


        endif   ! end if_flag_iter_and_flag
      enddo   ! j
      enddo   ! i

!> - Restore land-related prognostic fields for guess run.
      do j  = 1, 1
      do i  = 1, im
        if (flag(i)) then
          if(debug_print) print *,'end ',i,flag_guess(i),flag_iter(i)
          if (flag_guess(i)) then
            if(debug_print) print *,'guess run'

            weasd_lnd(i)       = weasd_lnd_old(i)
            snwdph_lnd(i)      = snwdph_lnd_old(i)
            tskin_lnd(i)       = tskin_lnd_old(i)
            canopy(i)          = canopy_old(i)
            !tprcp(i)           = tprcp_old(i)
            srflag(i)          = srflag_old(i)
            sr(i)              = sr_old(i)
            tsnow_lnd(i)       = tsnow_lnd_old(i)
            snowfallac_lnd(i)  = snowfallac_lnd_old(i)
            acsnow(i)          = acsnow_old(i)
            sfalb_lnd(i)       = sfalb_lnd_old(i)
            sfcqv_lnd(i)       = sfcqv_lnd_old(i)
            sfcqc_lnd(i)       = sfcqc_lnd_old(i)
            wet1(i)            = wet1_old(i)
            z0rl_lnd(i)        = z0rl_lnd_old(i)
            sncovr1_lnd(i)     = sncovr1_lnd_old(i)
            !ice
            weasd_ice(i)       = weasd_ice_old(i)
            snwdph_ice(i)      = snwdph_ice_old(i)
            tskin_ice(i)       = tskin_ice_old(i)
            tsnow_ice(i)       = tsnow_ice_old(i)
            snowfallac_ice(i)  = snowfallac_ice_old(i)
            sfalb_ice(i)       = sfalb_ice_old(i)
            sfcqv_ice(i)       = sfcqv_ice_old(i)
            sfcqc_ice(i)       = sfcqc_ice_old(i)
            z0rl_ice(i)        = z0rl_ice_old(i)
            sncovr1_ice(i)     = sncovr1_ice_old(i)

            do k = 1, lsoil_ruc
              smois(i,k)    = smois_old(i,k)
              tslb(i,k)     = tslb_old(i,k)
              tslb_lnd(i,k) = tslb_lnd_old(i,k)
              tsice(i,k)    = tsice_old(i,k)
              sh2o(i,k)     = sh2o_old(i,k)
              keepfr(i,k)   = keepfr_old(i,k)
              smfrkeep(i,k) = smfrkeep_old(i,k)
            enddo
          else
            if(debug_print) print *,'iter run', i,j, tskin_ice(i),tsurf_ice(i)
            tskin_lnd(i) = tsurf_lnd(i)
            tskin_ice(i) = tsurf_ice(i)
            tice(i)      = tsurf_ice(i)
          endif
        endif
      enddo  ! i
      enddo  ! j
!
      do i  = 1, im   ! i - horizontal loop
        if ( flag_iter(i) ) then
        ! Compute composite for a fractional sea ice: fice(i) < 1.
        ! This is needed for the 2-way coupling in ithe GFS physics driver 
        ! in the upcoupled case (when sfc_cice is not used).
          if(flag_ice_uncoupled(i) .and. fice(i) < 1.) then
            print *,'Fractional sea ice at i', i, fice(i)
            focean =  1.0 - fice(i)
           ! Check if ice fraction is below the minimum value: 15% in GFS physics.
            if (fice(i) < cimin) then ! cimin - minimal ice fraction
              print *,'warning: ice fraction is low:', fice(i)
              fice(i) = cimin
              focean  = 1.0 - cimin
              !tskin_ocn(i) = con_tice
              print *,'fix ice fraction: reset it to:', fice(i), tskin_ocn(i)
            endif
           ! Compute fluxes over open water
            qss_ocean  = rslf(prsl1(i),tskin_ocn(i))
            evap_ocean = ch_ocn(i) * wind(i) * (qss_ocean - q0(i))
            hflx_ocean = ch_ocn(i) * wind(i) * (tskin_ocn(i) - t1(i)) * prslki(i)

          ! Compute the composite of ice and open water for 2-way coupling in the
          ! uncoupled sea-ice model. Use ice variables for the composite.
            tskin_ice(i) = tice(i)*fice(i) + tskin_ocn(i)*focean
            ch_ice(i)    = ch_ice(i)*fice(i)   + ch_ocn(i)*focean
            hflx_ice(i)  = hflx_ice(i)*fice(i) + hflx_ocean*focean
            evap_ice(i)  = evap_ice(i)*fice(i) + evap_ocean*focean
            qsurf_ice(i) = q1(i) + evap_ice(i) / (ch_ice(i) * wind(i))
          endif ! flag_ice_uncoupled(i) .and. fice(i) < 1.

      ! Compute composite for ice and land temperature in the grid cell
          if(land(i) .and. flag_ice_uncoupled(i)) then
          ! land and ice
            snowc (i) = (sncovr1_lnd(i) * lndfrac(i) + sncovr1_ice(i) *fice(i))/ &
                        ( lndfrac (i) + fice (i) )
            do k = 1, lsoil_ruc
             if(lndfrac(i) > 0.5) then
              tslb(i,k) = tslb_lnd(i,k)
             else
              tslb(i,k) = tsice(i,k)
             endif
            enddo
          elseif (land(i) .and. (.not. flag_ice_uncoupled(i))) then
        ! land and no ice
            snowc (i) = sncovr1_lnd(i)
            do k = 1, lsoil_ruc
              tslb(i,k) = tslb_lnd(i,k)
            enddo
          endif ! land(i) .and. flag_ice_uncoupled(i)

        endif ! flag_iter
      enddo ! i - horizontal loop

      deallocate(soilctop)
      deallocate(landusef)
!
      return
!...................................
      end subroutine lsm_ruc_run
!-----------------------------------
      subroutine rucinit      (restart, im, lsoil_ruc, lsoil, nlev,   & ! in
                               isot, soiltyp, vegtype, fice,          & ! in
                               land, icy,                             & ! in
                               islmsk, tsurf, tsurf_ocn, tg3,         & ! in
                               smc, slc, stc,                         & ! in
                               smcref2, smcwlt2,                      & ! inout
                               lsm_ruc, lsm,                          & ! in
                               zs, sh2o, smfrkeep, tslb, smois, wet1, & ! out
                               errmsg, errflg)

      implicit none

      logical,                                 intent(in   ) :: restart
      integer,                                 intent(in   ) :: lsm
      integer,                                 intent(in   ) :: lsm_ruc
      integer,                                 intent(in   ) :: isot
      integer,                                 intent(in   ) :: im, nlev
      integer,                                 intent(in   ) :: lsoil_ruc
      integer,                                 intent(in   ) :: lsoil
      logical,               dimension(im),    intent(in   ) :: land, icy
      integer,               dimension(im),    intent(in   ) :: islmsk
      real (kind=kind_phys), dimension(im),    intent(in   ) :: tsurf, tsurf_ocn
      real (kind=kind_phys), dimension(im),    intent(inout) :: smcref2
      real (kind=kind_phys), dimension(im),    intent(inout) :: smcwlt2
      real (kind=kind_phys), dimension(im),    intent(in   ) :: tg3
      real (kind=kind_phys), dimension(im,lsoil),  intent(in   ) :: smc !  Noah
      real (kind=kind_phys), dimension(im,lsoil),  intent(in   ) :: stc !  Noah
      real (kind=kind_phys), dimension(im,lsoil),  intent(in   ) :: slc !  Noah

      integer,               dimension(im),    intent(inout) :: soiltyp
      integer,               dimension(im),    intent(inout) :: vegtype
      real (kind=kind_phys), dimension(im),    intent(inout) :: wet1
      real (kind=kind_phys), dimension(im),    intent(inout) :: fice
      real (kind=kind_phys), dimension(im,lsoil_ruc), intent(inout) :: smois! ruc
      real (kind=kind_phys), dimension(im,lsoil_ruc), intent(inout) :: tslb ! ruc
      real (kind=kind_phys), dimension(im,lsoil_ruc), intent(inout) :: sh2o ! ruc
      real (kind=kind_phys), dimension(im,lsoil_ruc), intent(inout) :: smfrkeep ! ruc

      real (kind=kind_phys), dimension(1:lsoil_ruc), intent (out)  :: zs

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

!> local
      logical :: debug_print
      logical :: smadj ! for soil mosture adjustment
      logical :: swi_init ! for initialization in terms of SWI (soil wetness index)

      integer :: flag_soil_layers, flag_soil_levels, flag_sst
      real (kind=kind_phys),    dimension(1:lsoil_ruc) :: factorsm

      integer , dimension( 1:im , 1:1 )      :: ivgtyp
      integer , dimension( 1:im , 1:1)       :: isltyp
      real (kind=kind_phys),    dimension( 1:im , 1:1 )       :: mavail
      real (kind=kind_phys),    dimension( 1:im , 1:1 )       :: xice
      real (kind=kind_phys),    dimension( 1:im , 1:1 )       :: sst
      real (kind=kind_phys),    dimension( 1:im , 1:1 )       :: landmask
      real (kind=kind_phys),    dimension( 1:im , 1:1 )       :: tsk
      real (kind=kind_phys),    dimension( 1:im , 1:1 )       :: tbot
      real (kind=kind_phys),    dimension( 1:im , 1:1 )       :: smtotn
      real (kind=kind_phys),    dimension( 1:im , 1:1 )       :: smtotr
      real (kind=kind_phys),    dimension( 1:im , 1:lsoil_ruc, 1:1 ) :: dumsm 
      real (kind=kind_phys),    dimension( 1:im , 1:lsoil_ruc, 1:1 ) :: dumt
      real (kind=kind_phys),    dimension( 1:im , 1:lsoil_ruc, 1:1 ) :: smfr
      real (kind=kind_phys),    dimension( 1:im , 1:lsoil_ruc, 1:1 ) :: soilm
      real (kind=kind_phys),    dimension( 1:im , 1:lsoil_ruc, 1:1 ) :: soiltemp
      real (kind=kind_phys),    dimension( 1:im , 1:lsoil_ruc, 1:1 ) :: soilh2o

      real (kind=kind_phys) :: st_input(1:im,1:lsoil_ruc*3,1:1)
      real (kind=kind_phys) :: sm_input(1:im,1:lsoil_ruc*3,1:1)

      integer               :: ids,ide, jds,jde, kds,kde, &
                               ims,ime, jms,jme, kms,kme, &
                               its,ite, jts,jte, kts,kte, &
                               i, j, k, l, num_soil_layers, ipr

      real(kind=kind_phys), dimension(1:lsoil_ruc) :: zs2, dzs 
      integer,              dimension(1:lsoil)  :: st_levels_input ! 4 - for Noah lsm
      integer,              dimension(1:lsoil)  :: sm_levels_input ! 4 - for Noah lsm

      ! Initialize the CCPP error handling variables
      errmsg = ''
      errflg = 0

      debug_print = .false.

      if (lsm/=lsm_ruc) then
        write(errmsg,'(a,i0,a,i0)')                                &
              'ERROR in lsm_ruc_init: namelist variable lsm=',     &
              lsm, ' incompatible with RUC LSM, please set to ', lsm_ruc
        errflg = 1
        return
      else if (debug_print) then
        write(0,*) 'Start of RUC LSM initialization'
      endif

      ipr = 10

      ! Set internal dimensions
      ids = 1
      ims = 1
      its = 1
      ide = im
      ime = im
      ite = im
      jds = 1
      jms = 1
      jts = 1
      jde = 1 
      jme = 1
      jte = 1
      kds = 1
      kms = 1
      kts = 1
      kde = nlev
      kme = nlev
      kte = nlev

      ! Initialize the RUC soil levels, needed for cold starts and warm starts
      CALL init_soil_depth_3 ( zs , dzs , lsoil_ruc )

      !! Check if RUC soil data (tslb, ...) is provided or not
      !if (minval(tslb)==maxval(tslb)) then
      ! For restart runs, can assume that RUC soul data is provided
      if (.not.restart) then

        flag_soil_layers = 1  ! =1 for input from the Noah LSM
        flag_soil_levels = 0  ! =1 for input from RUC LSM
        flag_sst = 0

        num_soil_layers =  lsoil ! 4 - for Noah lsm

        ! for Noah input set smadj and swi_init to .true.
        smadj = .true.
        swi_init = .true.
        
        if(lsoil == 4 ) then ! for Noah input
          st_levels_input = (/ 5, 25, 70, 150/)    ! Noah soil levels
          sm_levels_input = (/ 5, 25, 70, 150/)    ! Noah soil levels
        else
          write(errmsg,'(a,i0,a)')                                   &
                'WARNING in lsm_ruc_init: non-Noah input, lsoil=', lsoil
          errflg = 1
          return
        endif

      else

        ! For RUC input data, return here
        return

      endif

      if(debug_print) then
         print *,'Land mask islmsk(ipr) ==', ipr, islmsk(ipr)
         print *,'Noah smc(ipr,:) ==', ipr, smc(ipr,:)
         print *,'Noah stc(ipr,:) ==', ipr, stc(ipr,:)
         print *,'Noah vegtype(ipr) ==', ipr, vegtype(ipr)
         print *,'Noah soiltyp(ipr) ==', ipr, soiltyp(ipr)
         print *,'its,ite,jts,jte ',its,ite,jts,jte 
      endif

      ! Noah lsm input
      if ( flag_soil_layers == 1 ) then

        do j=jts,jte !
        do i=its,ite ! i = horizontal loop

         if(land(i)) then
          tsk(i,j) = tsurf(i)
         else
          tsk(i,j) = tsurf_ocn(i)
         endif

          tbot(i,j)=tg3(i)

          !SLMSK   - SEA(0),LAND(1),ICE(2) MASK
          if (land(i)) then
            ivgtyp(i,j)=vegtype(i)
            isltyp(i,j)=soiltyp(i)
            landmask(i,j)=1.
            xice(i,j)=0.
          elseif (icy(i)) then
            ivgtyp(i,j)=15 ! MODIS
            !> -- number of soil categories
            if(isot == 1) then
              isltyp(i,j) = 16 ! STATSGO
            else
              isltyp(i,j) = 9  ! ZOBLER
            endif
            landmask(i,j)=1.
            xice(i,j)=fice(i)
          else
            ivgtyp(i,j)= 17 ! 17 - water (oceans and lakes) in MODIS
            isltyp(i,j)=14
            xice(i,j)=0.
            landmask(i,j)=0.
          endif

          sst(i,j) = tsurf_ocn(i)

          st_input(i,1,j)=tsk(i,j)
          sm_input(i,1,j)=0.

          !--- initialize smcwlt2 and smcref2 with Noah values
          if(.not. land(i) ) then
            !water and sea ice
            smcref2 (i) = 1.
            smcwlt2 (i) = 0.
          else
            !land 
            smcref2 (i) = REFSMCnoah(soiltyp(i))
            smcwlt2 (i) = WLTSMCnoah(soiltyp(i))
          endif

          do k=1,lsoil
             st_input(i,k+1,j)=stc(i,k)
             ! convert volumetric soil moisture to SWI (soil wetness index)
             if(swi_init) then
               sm_input(i,k+1,j)=min(1.,max(0.,(smc(i,k) - smcwlt2(i))/  &
                                 (smcref2(i) - smcwlt2(i))))
             else
               sm_input(i,k+1,j)=smc(i,k)
             endif
          enddo
          do k=lsoil+2,lsoil_ruc * 3
             st_input(i,k,j)=0.
             sm_input(i,k,j)=0.
          enddo

        enddo ! i - horizontal loop
        enddo ! jme

        if(debug_print) then
          print *,'st_input=',ipr, st_input(ipr,:,1)
          print *,'sm_input=',ipr, sm_input(ipr,:,1)
        endif

        CALL init_soil_3_real ( tsk , tbot , dumsm , dumt ,             &
                                st_input , sm_input , landmask , sst ,  &
                                zs , dzs ,                              &
                                st_levels_input, sm_levels_input,       &
                                lsoil_ruc , num_soil_layers,            &
                                num_soil_layers,                        &
                                lsoil_ruc * 3 , lsoil_ruc * 3 ,         &
                                flag_sst,                               &
                                flag_soil_layers , flag_soil_levels ,   &
                                ids , ide , jds , jde , kds , kde ,     &
                                ims , ime , jms , jme , kms , kme ,     &
                                its , ite , jts , jte , kts , kte )

        do j=jts,jte
        do i=its,ite
          do k=1,lsoil_ruc
           ! convert from SWI to RUC volumetric soil moisture
           if(swi_init) then
             if(land(i)) then
               !land 
               soilm(i,k,j)= dumsm(i,k,j) *                             &
                 (refsmc(isltyp(i,j))-drysmc(isltyp(i,j)))              &
                 + drysmc(isltyp(i,j))
             else
               soilm(i,k,j)= 1.
             endif
           else
             soilm(i,k,j)= dumsm(i,k,j)
           endif
             soiltemp(i,k,j) = dumt(i,k,j)
          enddo
        enddo
        enddo

        if(debug_print) then
          print *,'tsk(i,j),tbot(i,j),sst(i,j),landmask(i,j)' &
                  ,ipr,1,tsk(ipr,1),tbot(ipr,1),sst(ipr,1),landmask(ipr,1)
          print *,'islmsk(ipr)=',ipr,islmsk(ipr)
          print *,'tsurf(ipr)=',ipr,tsurf(ipr)
          print *,'stc(ipr)=',ipr,stc(ipr,:)
          print *,'smc(ipr)=',ipr,smc(ipr,:)
          print *,'soilt(1,:,ipr)',ipr,soiltemp(ipr,:,1)
          print *,'soilm(1,:,ipr)',ipr,soilm(ipr,:,1)
        endif ! debug_print

        ! smadj should be true when the Noah LSM is used to initialize RUC
        if( smadj ) then
        ! With other LSMs as input, or when RUC soil moisture is cycled, it
        ! should be set to .false.

          do j=jts,jte
          do i=its,ite

          IF ( land(i) ) then  ! Land
            ! initialize factor
            do k=1,lsoil_ruc
               factorsm(k)=1.
            enddo
          
            ! RUC soil moisture bucket
            smtotr(i,j)=0.
            do k=1,lsoil_ruc -1
              smtotr(i,j)=smtotr(i,j) + soilm(i,k,j) *dzs(k)
            enddo
            ! Noah soil moisture bucket 
            smtotn(i,j)=smc(i,1)*0.1 + smc(i,2)*0.2 + smc(i,3)*0.7 + smc(i,4)*1.
            
            if(debug_print) then
              if(i==ipr) then
              print *,'from Noah to RUC: RUC bucket and Noah bucket at',    &
                       i,j,smtotr(i,j),smtotn(i,j)
              print *,'before smois=',i,j,soilm(i,:,j)
              endif
            endif
          
            ! RUC soil moisture correction to match Noah soil moisture bucket
            do k=1,lsoil_ruc-1
              soilm(i,k,j) = max(0.02,soilm(i,k,j)*smtotn(i,j)/(0.9*smtotr(i,j)))
            enddo
          
            if( soilm(i,2,j) > soilm(i,1,j) .and. soilm(i,3,j) > soilm(i,2,j)) then
            ! typical for daytime, no recent precip
              factorsm(1) = 0.75
              factorsm(2) = 0.8
              factorsm(3) = 0.85
              factorsm(4) = 0.9
              factorsm(5) = 0.95
            endif
            do k=1,lsoil_ruc
               soilm(i,k,j) = factorsm(k) * soilm(i,k,j)
            enddo
            if(debug_print) then
               if(i==ipr) print *,'after smois=',i,j,soilm(i,:,j)
            endif
               smtotr(i,j) = 0.
            do k=1,lsoil_ruc - 1
               smtotr(i,j)=smtotr(i,j) + soilm(i,k,j) *dzs(k)
            enddo
            if(debug_print) then
                if(i==ipr)print *,'after correction: RUC bucket and Noah bucket at',  &
                         i,j,smtotr(i,j),smtotn(i,j)
            endif
          ENDIF ! land

          enddo
          enddo

        endif ! smadj==.true.

        ! Initialize liquid and frozen soil moisture from total soil moisture
        ! and soil temperature, and also soil moisture availability in the top
        ! layer
        call ruclsminit( debug_print,                                   &
                   lsoil_ruc, isltyp, ivgtyp, xice, mavail,             &
                   soilh2o, smfr, soiltemp, soilm,                      &
                   ims,ime, jms,jme, kms,kme,                           &
                   its,ite, jts,jte, kts,kte                            )

        do j=jts,jte
        do i=its,ite
          wet1(i) = mavail(i,j)
          do k = 1, lsoil_ruc
            smois(i,k) = soilm(i,k,j)
            tslb(i,k)  = soiltemp(i,k,j)
            sh2o(i,k)  = soilh2o(i,k,j)
            smfrkeep(i,k)  = smfr(i,k,j)
          enddo
        enddo
        enddo

      endif ! flag_soil_layers==1

      end subroutine rucinit

end module lsm_ruc

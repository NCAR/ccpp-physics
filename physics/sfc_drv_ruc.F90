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
!     snowc    - real, fractional snow cover                       im   !
!     stm      - real, total soil column moisture content (m)      im   !
!     zorl     - real, surface roughness                           im   !
!     wet1     - real, normalized soil wetness                     im   !
!                                                                       !
!  ====================    end of description    =====================  !

#if 0
!! \section arg_table_lsm_ruc_run Argument Table
!! | local_name      | standard_name                                                                | long_name                                                       | units         | rank | type      |    kind   | intent | optional |
!! |-----------------|------------------------------------------------------------------------------|-----------------------------------------------------------------|---------------|------|-----------|-----------|--------|----------|
!! | delt            | time_step_for_dynamics                                                       | physics time step                                               | s             |    0 | real      | kind_phys | in     | F        |
!! | me              | mpi_rank                                                                     | current MPI-rank                                                | index         |    0 | integer   |           | in     | F        |
!! | kdt             | index_of_time_step                                                           | current number of time steps                                    | index         |    0 | integer   |           | in     | F        |
!! | iter            | ccpp_loop_counter                                                            | loop counter for subcycling loops in CCPP                       | index         |    0 | integer   |           | in     | F        |
!! | im              | horizontal_loop_extent                                                       | horizontal loop extent                                          | count         |    0 | integer   |           | in     | F        |
!! | nlev            | vertical_dimension                                                           | number of vertical levels                                       | count         |    0 | integer   |           | in     | F        |
!! | lsm_ruc         | flag_for_ruc_land_surface_scheme                                             | flag for RUC land surface model                                 | flag          |    0 | integer   |           | in     | F        |
!! | lsm             | flag_for_land_surface_scheme                                                 | flag for land surface model                                     | flag          |    0 | integer   |           | in     | F        |
!! | do_mynnsfclay   | do_mynnsfclay                                                                | flag to activate MYNN surface layer                             | flag          |    0 | logical   |           | in     | F        |
!! | lsoil_ruc       | soil_vertical_dimension_for_land_surface_model                               | number of soil layers internal to land surface model            | count         |    0 | integer   |           | in     | F        |
!! | lsoil           | soil_vertical_dimension                                                      | soil vertical layer dimension                                   | count         |    0 | integer   |           | in     | F        |
!! | zs              | depth_of_soil_levels_for_land_surface_model                                  | depth of soil levels for land surface model                     | m             |    1 | real      | kind_phys | inout  | F        |
!! | islmsk          | sea_land_ice_mask                                                            | landmask: sea/land/ice=0/1/2                                    | flag          |    1 | integer   |           | in     | F        |
!! | con_cp          | specific_heat_of_dry_air_at_constant_pressure                                | specific heat !of dry air at constant pressure                  | J kg-1 K-1    |    0 | real      | kind_phys | in     | F        |
!! | con_g           | gravitational_acceleration                                                   | gravitational acceleration                                      | m s-2         |    0 | real      | kind_phys | in     | F        |
!! | con_pi          | pi                                                                           | ratio of a circle's circumference to its diameter               | radians       |    0 | real      | kind_phys | in     | F        |
!! | con_rd          | gas_constant_dry_air                                                         | ideal gas constant for dry air                                  | J kg-1 K-1    |    0 | real      | kind_phys | in     | F        |
!! | con_rv          | gas_constant_water_vapor                                                     | ideal gas constant for water vapor                              | J kg-1 K-1    |    0 | real      | kind_phys | in     | F        |
!! | con_hvap        | latent_heat_of_vaporization_of_water_at_0C                                   | latent heat of vaporization/sublimation (hvap)                  | J kg-1        |    0 | real      | kind_phys | in     | F        |
!! | con_fvirt       | ratio_of_vapor_to_dry_air_gas_constants_minus_one                            | rv/rd - 1 (rv = ideal gas constant for water vapor)             | none          |    0 | real      | kind_phys | in     | F        |
!! | rainnc          | lwe_thickness_of_explicit_rainfall_amount_from_previous_timestep             | explicit rainfall from previous timestep                        | m             |    1 | real      | kind_phys | in     | F        |
!! | rainc           | lwe_thickness_of_convective_precipitation_amount_from_previous_timestep      | convective_precipitation_amount from previous timestep          | m             |    1 | real      | kind_phys | in     | F        |
!! | ice             | lwe_thickness_of_ice_amount_from_previous_timestep                           | ice amount from previous timestep                               | m             |    1 | real      | kind_phys | in     | F        |
!! | snow            | lwe_thickness_of_snow_amount_from_previous_timestep                          | snow amount from previous timestep                              | m             |    1 | real      | kind_phys | in     | F        |
!! | graupel         | lwe_thickness_of_graupel_amount_from_previous_timestep                       | graupel amount from previous timestep                           | m             |    1 | real      | kind_phys | in     | F        |
!! | srflag          | flag_for_precipitation_type                                                  | snow/rain flag for precipitation                                | flag          |    1 | real      | kind_phys | in     | F        |
!! | sncovr1         | surface_snow_area_fraction_for_diagnostics                                   | surface snow area fraction                                      | frac          |    1 | real      | kind_phys | inout  | F        |
!! | snowc           | surface_snow_area_fraction                                                   | surface snow area fraction                                      | frac          |    1 | real      | kind_phys | inout  | F        |
!! | weasd           | water_equivalent_accumulated_snow_depth                                      | water equiv of acc snow depth over land and sea ice             | mm            |    1 | real      | kind_phys | inout  | F        |
!! | snwdph          | surface_snow_thickness_water_equivalent                                      | water equivalent snow depth over land                           | mm            |    1 | real      | kind_phys | inout  | F        |
!! | sr              | ratio_of_snowfall_to_rainfall                                                | snow ratio: ratio of snow to total precipitation                | frac          |    1 | real      | kind_phys | in     | F        |
!! | rhosnf          | density_of_frozen_precipitation                                              | density of frozen precipitation                                 | kg m-3        |    1 | real      | kind_phys | out    | F        |
!! | zf              | height_above_ground_at_lowest_model_layer                                    | layer 1 height above ground (not MSL)                           | m             |    1 | real      | kind_phys | in     | F        |
!! | u1              | x_wind_at_lowest_model_layer                                                 | zonal wind at lowest model layer                                | m s-1         |    1 | real      | kind_phys | in     | F        |
!! | v1              | y_wind_at_lowest_model_layer                                                 | meridional wind at lowest model layer                           | m s-1         |    1 | real      | kind_phys | in     | F        |
!! | prsl1           | air_pressure_at_lowest_model_layer                                           | mean pressure at lowest model layer                             | Pa            |    1 | real      | kind_phys | in     | F        |
!! | ddvel           | surface_wind_enhancement_due_to_convection                                   | surface wind enhancement due to convection                      | m s-1         |    1 | real      | kind_phys | in     | F        |
!! | t1              | air_temperature_at_lowest_model_layer                                        | mean temperature at lowest model layer                          | K             |    1 | real      | kind_phys | in     | F        |
!! | q1              | water_vapor_specific_humidity_at_lowest_model_layer                          | water vapor specific humidity at lowest model layer             | kg kg-1       |    1 | real      | kind_phys | in     | F        |
!! | qc              | cloud_condensed_water_mixing_ratio_at_lowest_model_layer | moist (dry+vapor, no condensates) mixing ratio of cloud water at lowest model layer | kg kg-1       |    1 | real      | kind_phys | in     | F        |
!! | dlwflx          | surface_downwelling_longwave_flux                                            | surface downwelling longwave flux at current time               | W m-2         |    1 | real      | kind_phys | in     | F        |
!! | dswsfc          | surface_downwelling_shortwave_flux                                           | surface downwelling shortwave flux at current time              | W m-2         |    1 | real      | kind_phys | in     | F        |
!! | snet            | surface_net_downwelling_shortwave_flux                                       | surface net downwelling shortwave flux at current time          | W m-2         |    1 | real      | kind_phys | in     | F        |
!! | sfcemis         | surface_longwave_emissivity                                                  | surface lw emissivity in fraction                               | frac          |    1 | real      | kind_phys | inout  | F        |
!! | wspd            | wind_speed_at_lowest_model_layer                                             | wind speed at lowest model level                                | m s-1         |    1 | real      | kind_phys | inout  | F        |
!! | cm              | surface_drag_coefficient_for_momentum_in_air                                 | surface exchange coeff for momentum                             | none          |    1 | real      | kind_phys | in     | F        |
!! | ch              | surface_drag_coefficient_for_heat_and_moisture_in_air                        | surface exchange coeff heat & moisture                          | none          |    1 | real      | kind_phys | in     | F        |
!! | chh             | surface_drag_mass_flux_for_heat_and_moisture_in_air                          | surf h&m exch coef time surf wind & density                     | kg m-2 s-1    |    1 | real      | kind_phys | inout  | F        |
!! | cmm             | surface_drag_wind_speed_for_momentum_in_air                                  | surf mom exch coef time mean surf wind                          | m s-1         |    1 | real      | kind_phys | inout  | F        |
!! | wet1            | normalized_soil_wetness                                                      | normalized soil wetness                                         | frac          |    1 | real      | kind_phys | inout  | F        |
!! | canopy          | canopy_water_amount                                                          | canopy water amount                                             | kg m-2        |    1 | real      | kind_phys | inout  | F        |
!! | sigmaf          | vegetation_area_fraction                                                     | areal fractional cover of green vegetation                      | frac          |    1 | real      | kind_phys | in     | F        |
!! | sfalb           | surface_diffused_shortwave_albedo                                            | mean surface diffused sw albedo                                 | frac          |    1 | real      | kind_phys | inout  | F        |
!! | alvwf           | mean_vis_albedo_with_weak_cosz_dependency                                    | mean vis albedo with weak cosz dependency                       | frac          |    1 | real      | kind_phys | in     | F        |
!! | alnwf           | mean_nir_albedo_with_weak_cosz_dependency                                    | mean nir albedo with weak cosz dependency                       | frac          |    1 | real      | kind_phys | in     | F        |
!! | snoalb          | upper_bound_on_max_albedo_over_deep_snow                                     | maximum snow albedo                                             | frac          |    1 | real      | kind_phys | in     | F        |
!! | zorl            | surface_roughness_length                                                     | surface roughness length                                        | cm            |    1 | real      | kind_phys | inout  | F        |
!! | qsurf           | surface_specific_humidity                                                    | surface air saturation specific humidity                        | kg kg-1       |    1 | real      | kind_phys | inout  | F        |
!! | sfcqc           | cloud_condensed_water_mixing_ratio_at_surface                                | moist cloud water mixing ratio at surface                       | kg kg-1       |    1 | real      | kind_phys | inout  | F        |
!! | sfcqv           | water_vapor_mixing_ratio_at_surface                                          | water vapor mixing ratio at surface                             | kg kg-1       |    1 | real      | kind_phys | inout  | F        |
!! | sfcdew          | surface_condensation_mass                                                    | surface condensation mass                                       | kg m-2        |    1 | real      | kind_phys | inout  | F        |
!! | tg3             | deep_soil_temperature                                                        | deep soil temperature                                           | K             |    1 | real      | kind_phys | in     | F        |
!! | smc             | volume_fraction_of_soil_moisture                                             | total soil moisture                                             | frac          |    2 | real      | kind_phys | inout  | F        |
!! | slc             | volume_fraction_of_unfrozen_soil_moisture                                    | liquid soil moisture                                            | frac          |    2 | real      | kind_phys | inout  | F        |
!! | stc             | soil_temperature                                                             | soil temperature                                                | K             |    2 | real      | kind_phys | inout  | F        |
!! | smcwlt2         | volume_fraction_of_condensed_water_in_soil_at_wilting_point                  | soil water fraction at wilting point                            | frac          |    1 | real      | kind_phys | inout  | F        |
!! | smcref2         | threshold_volume_fraction_of_condensed_water_in_soil                         | soil moisture threshold                                         | frac          |    1 | real      | kind_phys | inout  | F        | 
!! | vegtype         | vegetation_type_classification                                               | vegetation type at each grid cell                               | index         |    1 | integer   |           | in     | F        |
!! | soiltyp         | soil_type_classification                                                     | soil type at each grid cell                                     | index         |    1 | integer   |           | in     | F        |
!! | isot            | soil_type_dataset_choice                                                     | soil type dataset choice                                        | index         |    0 | integer   |           | in     | F        |
!! | ivegsrc         | vegetation_type_dataset_choice                                               | land use dataset choice                                         | index         |    0 | integer   |           | in     | F        |
!! | fice            | sea_ice_concentration                                                        | ice fraction over open water                                    | frac          |    1 | real      | kind_phys | in     | F        |
!! | keepfr          | flag_for_frozen_soil_physics                                                 | flag for frozen soil physics (RUC)                              | flag          |    2 | real      | kind_phys | inout  | F        |
!! | smois           | volume_fraction_of_soil_moisture_for_land_surface_model                      | volumetric fraction of soil moisture for lsm                    | frac          |    2 | real      | kind_phys | inout  | F        |
!! | sh2o            | volume_fraction_of_unfrozen_soil_moisture_for_land_surface_model             | volume fraction of unfrozen soil moisture for lsm               | frac          |    2 | real      | kind_phys | inout  | F        |
!! | smfrkeep        | volume_fraction_of_frozen_soil_moisture_for_land_surface_model               | volume fraction of frozen soil moisture for lsm                 | frac          |    2 | real      | kind_phys | inout  | F        |
!! | tslb            | soil_temperature_for_land_surface_model                                      | soil temperature for land surface model                         | K             |    2 | real      | kind_phys | inout  | F        |
!! | stm             | soil_moisture_content                                                        | soil moisture content                                           | kg m-2        |    1 | real      | kind_phys | inout  | F        |
!! | tskin           | surface_skin_temperature                                                     | surface skin temperature                                        | K             |    1 | real      | kind_phys | inout  | F        |
!! | tsurf           | surface_skin_temperature_after_iteration                                     | surface skin temperature after iteration                        | K             |    1 | real      | kind_phys | inout  | F        |
!! | tice            | sea_ice_temperature                                                          | sea ice surface skin temperature                                | K             |    1 | real      | kind_phys | inout  | F        |
!! | tsnow           | snow_temperature_bottom_first_layer                                          | snow temperature at the bottom of first snow layer              | K             |    1 | real      | kind_phys | inout  | F        |
!! | snowfallac      | total_accumulated_snowfall                                                   | run-total snow accumulation on the ground                       | kg m-2        |    1 | real      | kind_phys | inout  | F        |
!! | acsnow          | accumulated_water_equivalent_of_frozen_precip                                | snow water equivalent of run-total frozen precip                | kg m-2        |    1 | real      | kind_phys | inout  | F        |
!! | evap            | kinematic_surface_upward_latent_heat_flux                                    | kinematic surface upward evaporation flux                       | kg kg-1 m s-1 |    1 | real      | kind_phys | out    | F        |
!! | hflx            | kinematic_surface_upward_sensible_heat_flux                                  | kinematic surface upward sensible heat flux                     | K m s-1       |    1 | real      | kind_phys | out    | F        |
!! | evbs            | soil_upward_latent_heat_flux                                                 | soil upward latent heat flux                                    | W m-2         |    1 | real      | kind_phys | out    | F        |
!! | evcw            | canopy_upward_latent_heat_flux                                               | canopy upward latent heat flux                                  | W m-2         |    1 | real      | kind_phys | out    | F        |
!! | sbsno           | snow_deposition_sublimation_upward_latent_heat_flux                          | latent heat flux from snow depo/subl                            | W m-2         |    1 | real      | kind_phys | out    | F        |
!! | trans           | transpiration_flux                                                           | total plant transpiration rate                                  | kg m-2 s-1    |    1 | real      | kind_phys | out    | F        |
!! | runof           | surface_runoff_flux                                                          | surface runoff flux                                             | g m-2 s-1     |    1 | real      | kind_phys | out    | F        |
!! | drain           | subsurface_runoff_flux                                                       | subsurface runoff flux                                          | g m-2 s-1     |    1 | real      | kind_phys | out    | F        |
!! | runoff          | total_runoff                                                                 | total water runoff                                              | kg m-2        |    1 | real      | kind_phys | inout  | F        |
!! | srunoff         | surface_runoff                                                               | surface water runoff (from lsm)                                 | kg m-2        |    1 | real      | kind_phys | inout  | F        |
!! | gflux           | upward_heat_flux_in_soil                                                     | soil heat flux                                                  | W m-2         |    1 | real      | kind_phys | out    | F        |
!! | shdmin          | minimum_vegetation_area_fraction                                             | min fractional coverage of green vegetation                     | frac          |    1 | real      | kind_phys | in     | F        |
!! | shdmax          | maximum_vegetation_area_fraction                                             | max fractional coverage of green vegetation                     | frac          |    1 | real      | kind_phys | in     | F        |
!! | flag_iter       | flag_for_iteration                                                           | flag for iteration                                              | flag          |    1 | logical   |           | in     | F        |
!! | flag_guess      | flag_for_guess_run                                                           | flag for guess run                                              | flag          |    1 | logical   |           | in     | F        |
!! | flag_init       | flag_for_first_time_step                                                     | flag signaling first time step for time integration loop        | flag          |    0 | logical   |           | in     | F        |
!! | flag_restart    | flag_for_restart                                                             | flag for restart (warmstart) or coldstart                       | flag          |    0 | logical   |           | in     | F        |
!! | errmsg          | ccpp_error_message                                                           | error message for error handling in CCPP                        | none          |    0 | character | len=*     | out    | F        |
!! | errflg          | ccpp_error_flag                                                              | error flag for error handling in CCPP                           | flag          |    0 | integer   |           | out    | F        |
!!
#endif

      subroutine lsm_ruc_run                                            &
! --- inputs
     &     ( iter, me, kdt, im, nlev, lsoil_ruc, lsoil, zs,             &
     &       u1, v1, t1, q1, qc, soiltyp, vegtype, sigmaf,              &
     &       sfcemis, dlwflx, dswsfc, snet, delt, tg3, cm, ch,          &
     &       prsl1, zf, islmsk, ddvel, shdmin, shdmax, alvwf, alnwf,    &
     &       snoalb, sfalb, flag_iter, flag_guess, isot, ivegsrc, fice, &
     &       smc, stc, slc, lsm_ruc, lsm,                               &
     &       smcwlt2, smcref2, wspd, do_mynnsfclay,                     &
! --- constants
     &       con_cp, con_rv, con_rd, con_g, con_pi, con_hvap, con_fvirt,&
! --- in/outs
     !&       weasd, snwdph, tskin, tprcp, rain, rainc, snow,            &
     &       weasd, snwdph, tskin,                                      &
! --- in
     &       rainnc, rainc, ice, snow, graupel,                         &
! --- in/outs
     &       srflag, sr,                                                &
     !&       graupel, srflag, sr,                                       &
     &       smois, tslb, sh2o, keepfr, smfrkeep,                       & ! on RUC levels
     &       canopy, trans, tsurf, tsnow, zorl,                         &
     &       sfcqc, sfcdew, tice, sfcqv,                                &
! --- outputs
     &       sncovr1, qsurf, gflux, drain, evap, hflx,                  &
     &       rhosnf, runof, runoff, srunoff,                            &
     &       chh, cmm, evbs, evcw, sbsno, snowc, stm, wet1,             &
     &       acsnow, snowfallac,                                        &
     &       flag_init, flag_restart, errmsg, errflg                    &
     &     )

      implicit none

!  ---  constant parameters:
      real(kind=kind_phys), parameter :: rhoh2o  = 1000.0
      real(kind=kind_phys), parameter :: stbolt  = 5.670400e-8

!  ---  input:
      integer, intent(in) :: me
      integer, intent(in) :: im, nlev, iter, lsoil_ruc, lsoil, kdt, isot, ivegsrc
      integer, intent(in) :: lsm_ruc, lsm

      real (kind=kind_phys), dimension(im,lsoil), intent(inout) :: smc,stc,slc

      real (kind=kind_phys), dimension(im), intent(in) :: u1, v1,&
     &       t1, sigmaf, sfcemis, dlwflx, dswsfc, snet, tg3, cm,        &
     &       ch, prsl1, ddvel, shdmin, shdmax,                          &
     &       snoalb, alvwf, alnwf, zf, qc, q1, wspd

      integer, dimension(im), intent(in) :: islmsk
      real (kind=kind_phys),  intent(in) :: delt
      real (kind=kind_phys),  intent(in) :: con_cp, con_rv, con_g,      &
                                            con_pi, con_rd,             &
                                            con_hvap, con_fvirt

      logical, dimension(im), intent(in) :: flag_iter, flag_guess
      logical,                intent(in) :: do_mynnsfclay

!  ---  in/out:
      integer, dimension(im), intent(inout) :: soiltyp, vegtype
      real (kind=kind_phys), dimension(lsoil_ruc) :: dzs
      real (kind=kind_phys), dimension(lsoil_ruc), intent(inout   ) :: zs
      real (kind=kind_phys), dimension(im), intent(inout) :: weasd,     &
!     &       snwdph, tskin, tprcp, rain, rainc, graupel, snow,          &
     &       snwdph, tskin,                                             &
             srflag, sr, canopy, trans, tsurf, zorl, tsnow,             &
             sfcqc, sfcqv, sfcdew, fice, tice, sfalb, smcwlt2, smcref2
!  ---  in
      real (kind=kind_phys), dimension(im), intent(in) ::               &
     &       rainnc, rainc, ice, snow, graupel
!  ---  in/out:
!  --- on RUC levels
      real (kind=kind_phys), dimension(im,lsoil_ruc), intent(inout) ::         &
     &       smois, tslb, sh2o, keepfr, smfrkeep

!  ---  output:
      real (kind=kind_phys), dimension(im), intent(inout) :: sncovr1,   &
     &       qsurf , gflux , evap , runof , drain ,                     &
     &       runoff, srunoff, hflx, cmm, chh,                           &
     &       rhosnf, evbs, evcw, sbsno, snowc, stm, wet1,               &
     &       acsnow, snowfallac

      logical,          intent(in)  :: flag_init, flag_restart
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

!  ---  locals:
      real (kind=kind_phys), dimension(im) :: rch, rho,                 &
     &       q0, qs1, wind, weasd_old, snwdph_old,                      &
     &       tprcp_old, srflag_old, sr_old, tskin_old, canopy_old

      real (kind=kind_phys), dimension(lsoil_ruc) :: et

      real (kind=kind_phys), dimension(im,lsoil_ruc,1) :: smsoil,              &
           slsoil, stsoil, smfrsoil, keepfrsoil

      real (kind=kind_phys), dimension(im,lsoil_ruc) :: smois_old,               &
     &       tslb_old, sh2o_old, keepfr_old, smfrkeep_old

      real (kind=kind_phys),dimension (im,1,1)      ::                  &
     &     conflx2, sfcprs, sfctmp, q2, qcatm, rho2 
      real (kind=kind_phys),dimension (im,1)        ::                  &
     &     albbck, alb, chs, flhc, flqc, wet, smmax, cmc,               &
     &     dew, drip,  ec, edir, ett, lh, esnow, etp, qfx,              &
     &     acceta, ffrozp, lwdn, prcp, xland, xice,                     &
     &     graupelncv, snowncv, rainncv, raincv,                        &
     &     solnet, sfcexc,                                              &
     &     runoff1, runoff2, acrunoff,                                  &
     &     sfcems, hfx, shdfac, shdmin1d, shdmax1d,                     &
     &     sneqv, snoalb1d, snowh, snoh, tsnav,                         &
     &     snomlt, sncovr, soilw, soilm, ssoil, soilt, tbot,            &
     &     xlai, swdn, z0, znt, rhosnfr, infiltr,                       &
     &     precipfr, snfallac, acsn,                                    &
     &     qsfc, qsg, qvg, qcg, soilt1, chklowq

      real (kind=kind_phys) :: xice_threshold

      character(len=256) :: llanduse  ! Land-use dataset.  Valid values are :
                                      ! "USGS" (USGS 24/27 category dataset) and
                                      ! "MODIFIED_IGBP_MODIS_NOAH" (MODIS 20-category dataset)

      integer :: nscat, nlcat
      real (kind=kind_phys), dimension(:,:,:), allocatable :: landusef ! fractional landuse
      real (kind=kind_phys), dimension(:,:,:), allocatable :: soilctop ! fractional soil type

      integer :: nsoil, iswater, isice
      integer, dimension (1:im,1:1) :: stype, vtype
      integer :: ipr

! local
      integer :: ims,ime, its,ite, jms,jme, jts,jte, kms,kme, kts,kte
      integer :: l, k, i, j,  fractional_seaice

      logical :: flag(im)
      logical :: rdlai2d, myj, frpcpn
      logical :: debug_print
!
      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      ipr = 10

      debug_print=.false.

      chklowq = 1.

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
                               islmsk, tskin, tg3,                       & ! in
                               smc, slc, stc,                            & ! in
                               smcref2, smcwlt2,                         & ! inout
                               lsm_ruc, lsm,                             & ! in
                               zs, sh2o, smfrkeep, tslb, smois, wet1,    & ! out
                               errmsg, errflg)

        do i  = 1, im ! n - horizontal loop
          ! overwrite Noah soil fields with initialized RUC soil fields for output
          do k = 1, lsoil
            smc(i,k)   = smois(i,k)
            slc(i,k)   = sh2o(i,k)
            stc(i,k)   = tslb(i,k)
          enddo
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
        xice_threshold = 0.02
      endif

      nsoil = lsoil_ruc

      do i  = 1, im ! i - horizontal loop
        ! reassign smcref2 and smcwlt2 to RUC values
        if(islmsk(i) == 0 .or. islmsk(i) == 2) then
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
        flag(i) = (islmsk(i) == 1 .or. islmsk(i) == 2)
      enddo

      do i  = 1, im ! i - horizontal loop
        if (flag(i) .and. flag_guess(i)) then
          !if(me==0 .and. i==ipr) print *,'before call to RUC guess run', i
          weasd_old(i)  = weasd(i)
          snwdph_old(i) = snwdph(i)
          tskin_old(i)  = tskin(i)
          canopy_old(i) = canopy(i)
          !tprcp_old(i)  = tprcp(i)
          srflag_old(i) = srflag(i)
          sr_old(i)     = sr(i)

          !> - Save land-related prognostic fields for guess run.
          do k = 1, lsoil_ruc
            smois_old(i,k)  = smois(i,k)
            tslb_old(i,k)   = tslb(i,k)
            sh2o_old(i,k)   = sh2o(i,k)
            keepfr_old(i,k) = keepfr(i,k)
            smfrkeep_old(i,k) = smfrkeep(i,k)
          enddo
        endif
      enddo

!  --- ...  initialization block

      do j  = 1, 1
      do i  = 1, im ! i - horizontal loop
        if (flag_iter(i) .and. flag(i)) then
          !if(me==0 .and. i==ipr) print *,'iter run', iter, i, flag_iter(i),flag_guess(i)
          evap (i)  = 0.0
          hflx (i)  = 0.0
          gflux(i)  = 0.0
          drain(i)  = 0.0
          canopy(i) = max(canopy(i), 0.0)
          sfcdew(i) = 0.0

          evbs (i)  = 0.0
          evcw (i)  = 0.0
          trans(i)  = 0.0
          sbsno(i)  = 0.0
          snowc(i)  = 0.0

          !local i,j arrays
          dew(i,j)      = 0.0
          soilm(i,j)    = 0.0
          smmax(i,j)    = 0.0
          hfx(i,j)      = 0.0
          qfx(i,j)      = 0.0
          lh(i,j)       = 0.0
          acsn(i,j)     = 0.0
          sfcexc(i,j)   = 0.0
          acceta(i,j)   = 0.0
          ssoil(i,j)    = 0.0
          snomlt(i,j)   = 0.0
          infiltr(i,j)  = 0.0
          runoff1(i,j)  = 0.0
          runoff2(i,j)  = 0.0
          acrunoff(i,j) = 0.0
          snfallac(i,j) = 0.0
          rhosnfr(i,j)  = 0.0
          precipfr(i,j) = 0.0

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
!!\n  ----------------------------------------
!!\n  \a ffrozp  - fraction of frozen precipitation
!!\n  \a frpcpn  - .true. if mixed phase precipitation available
!!\n  \a 1:im - horizontal_loop_extent
!!\n  \a fice    - fraction of sea-ice in the grid cell
!!\n  \a delt    - timestep (sec) (dt should not exceed 3600 secs) 
!!\n  \a conflx2 - height (\f$m\f$) above ground of atmospheric forcing variables
!!\n  \a lsoil_ruc - number of soil layers (= 6 or 9)
!!\n  \a zs      - the depth of each soil level (\f$m\f$)

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

!>  -   2. forcing data (f):
!!\n  ---------------------------------------
!!\n  \a sfcprs  - pressure at height zf above ground (pascals)
!!\n  \a sfctmp  - air temperature (\f$K\f$) at height zf above ground
!!\n  \a q2      - pressure at height zf above ground (pascals)
!!\n  \a qcatm   - cloud water mising ration at height zf above ground (\f$kg !kg^{-1}\f$)
!!\n  \a rho2    - air density at height zf above ground (pascals)

        sfcprs(i,1,j)  = prsl1(i)
        sfctmp(i,1,j)  = t1(i)
        q2(i,1,j)      = q0(i)
        qcatm(i,1,j)   = max(0., qc(i))
        rho2(i,1,j)    = rho(i)

!!\n  ---------------------------------------
!!\n  \a lwdn    - lw dw radiation flux at surface (\f$W m^{-2}\f$)
!!\n  \a swdn    - sw dw radiation flux at surface (\f$W m^{-2}\f$)
!!\n  \a solnet  - net sw radiation flux (dn-up) (\f$W m^{-2}\f$)
!!\n  \a prcp    - time-step total precip (\f$kg m^{-2} \f$)
!!\n  \a raincv  - time-step convective precip (\f$kg m^{-2} \f$)
!!\n  \a rainncv - time-step non-convective precip (\f$kg m^{-2} \f$)
!!\n  \a graupelncv - time-step graupel (\f$kg m^{-2} \f$)
!!\n  \a snowncv - time-step snow (\f$kg m^{-2} \f$)
!!\n  \a precipfr - time-step precipitation in solod form (\f$kg m^{-2} \f$)
!!\n  \a qsfc    - specific humidity at surface (\f$kg kg^{-1}\f$)
!!\n  \a qvg     - water vapor mixing ratio at surface (\f$kg kg^{-1}\f$)
!!\n  \a qsg     - saturated water vapor mixing ratio at surface (\f$kg kg^{-1}\f$)
!!\n  \a qcg     - cloud water mixing ratio at surface (\f$kg kg^{-1}\f$)

        lwdn(i,j)   = dlwflx(i)         !..downward lw flux at sfc in w/m2
        swdn(i,j)   = dswsfc(i)         !..downward sw flux at sfc in w/m2
        solnet(i,j) = dswsfc(i)*(1.-sfalb(i)) !snet(i) !..net sw rad flx (dn-up) at sfc in w/m2

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

        qvg(i,j)    = sfcqv(i)
        qsfc(i,j)   = sfcqv(i)/(1.+sfcqv(i))
        qsg(i,j)    = rslf(prsl1(i),tsurf(i))
        qcg(i,j)    = sfcqc(i) 

!>  -   3. canopy/soil characteristics (s):
!!\n      ------------------------------------
!!\n \a vegtyp  - vegetation type (integer index)                   -> vtype
!!\n \a soiltyp - soil type (integer index)                         -> stype
!!\n \a shdfac  - areal fractional coverage of green vegetation (0.0-1.0)
!!\n \a shdmin  - minimum areal fractional coverage of green vegetation -> shdmin1d
!!\n \a shdmax  - maximum areal fractional coverage of green vegetation -> shdmax1d
!!\n \a sfcems  -  surface emmisivity                                   -> sfcemis
!!\n \a 0.5*(alvwf + alnwf) - backround snow-free surface albedo (fraction)         -> albbck
!!\n \a snoalb  - upper bound on maximum albedo over deep snow          -> snoalb1d
!!\n \a sfalb  - surface albedo including snow effect (unitless fraction) -> alb
!!\n \a tbot    - bottom soil temperature (local yearly-mean sfc air temp)

        if(ivegsrc == 1) then   ! IGBP - MODIS
        !> - Prepare land/ice/water masks for RUC LSM
        !SLMSK0   - SEA(0),LAND(1),ICE(2) MASK
          if(islmsk(i) == 0.) then
            vtype(i,j) = 17 ! 17 - water (oceans and lakes) in MODIS
            stype(i,j) = 14
            xland(i,j) = 2. ! xland = 2 for water
            xice(i,j)  = 0.
          elseif(islmsk(i) == 1.) then ! land
            vtype(i,j) = vegtype(i)
            stype(i,j) = soiltyp(i)
            xland(i,j) = 1.
            xice(i,j)  = 0.
          elseif(islmsk(i) == 2) then  ! ice
            vtype(i,j) = 15 ! MODIS
            if(isot == 0) then
              stype(i,j) = 9  ! ZOBLER
            else
              stype(i,j) = 16 ! STASGO
            endif
            xland(i,j) = 1.
            xice(i,j) = fice(i)  ! fraction of sea-ice
          endif
        else
          print *,'MODIS landuse is not available'
        endif

        ! --- units %
        shdfac(i,j) = sigmaf(i)*100.
        shdmin1d(i,j) = shdmin(i)*100.
        shdmax1d(i,j) = shdmax(i)*100.

        sfcems(i,j) = sfcemis(i)

        snoalb1d(i,j) = snoalb(i)
        albbck(i,j)   = max(0.01, 0.5 * (alvwf(i) + alnwf(i))) 
        alb(i,j)      = sfalb(i)

        tbot(i,j) = tg3(i)

!>  -   4. history (state) variables (h):
!!\n      ------------------------------
!!\n \a cmc        - canopy moisture content (\f$mm\f$)
!!\n \a soilt = tskin - ground/canopy/snowpack effective skin temperature (\f$K\f$)
!!\n \a soilt1 = snowpack temperature at the bottom of the 1st layer (\f$K\f$)
!!\n \a tslb(lsoil_ruc) - soil temp (\f$K\f$)                                    -> stsoil
!!\n \a smois(lsoil_ruc) - total soil moisture content (volumetric fraction)     -> smsoil
!!\n \a sh2o(lsoil_ruc) - unfrozen soil moisture content (volumetric fraction)   -> slsoil
!!\n \a smfrsoil(lsoil_ruc) - frozen soil moisture content (volumetric fraction) -> smfrsoil
!!\n \a keepfrflag(lsoil_ruc) - flag for frozen soil physics: 0. or 1.
!!\n \a wet        - soil moisture availability at surface
!!\n \a snowh      - actual snow depth (\f$m\f$)
!!\n \a sneqv      - liquid water-equivalent snow depth (\f$m\f$)
!!\n \a sncovr     - fraction of snow in the grid cell
!!\n \a ch         - surface exchange coefficient for heat (\f$m s^{-1}\f$)      -> chs
!!\n \a z0         - surface roughness (\f$m\f$)     -> zorl(\f$cm\f$)

        cmc(i,j) = canopy(i)            !  [mm] 
        soilt(i,j) = tsurf(i)            ! clu_q2m_iter
        ! sanity check for snow temperature tsnow
        if (tsnow(i) > 0. .and. tsnow(i) < 273.15) then
          soilt1(i,j) = tsnow(i)
        else
          soilt1(i,j) = tsurf(i)
        endif

          tsnav(i,j) = 0.5*(soilt(i,j) + soilt1(i,j)) - 273.15

        do k = 1, lsoil_ruc
          smsoil  (i,k,j) = smois(i,k)
          slsoil  (i,k,j) = sh2o(i,k)
          stsoil  (i,k,j) = tslb(i,k)
          smfrsoil(i,k,j) = smfrkeep(i,k)
          keepfrsoil(i,k,j) = keepfr(i,k)
        enddo

        if(stype(i,j) .ne. 14) then
           ! land
           if (wet1(i) > 0.) then
            wet(i,j) = wet1(i)
           else
            wet(i,j) = max(0.0001,smsoil(i,1,j)/0.3)
           endif
        else
          ! water
          wet(i,j) = 1.
        endif

        snowh(i,j) = snwdph(i) * 0.001         ! convert from mm to m
        sneqv(i,j) = weasd(i)                  ! [mm]

        snfallac(i,j) = snowfallac(i)
        acsn(i,j)     = acsnow(i)

        !> -- sanity checks on sneqv and snowh
        if (sneqv(i,j) /= 0.0 .and. snowh(i,j) == 0.0) then
          snowh(i,j) = 0.003 * sneqv(i,j) ! snow density ~300 kg m-3 
        endif

        if (snowh(i,j) /= 0.0 .and. sneqv(i,j) == 0.0) then
          sneqv(i,j) = 300. * snowh(i,j) ! snow density ~300 kg m-3 
        endif

        if (sneqv(i,j) > 0. .and. snowh(i,j) > 0.) then
          if(sneqv(i,j)/snowh(i,j) > 950.) then
            sneqv(i,j) = 300. * snowh(i,j)
          endif
        endif

        !sncovr(i,j) = snowc(i)
        sncovr(i,j) = sncovr1(i)

        chs(i,j)    = ch(i) * wind(i) ! compute conductance 
        flhc(i,j)   = chs(i,j) * rho(i) * con_cp * (1. + 0.84*q2(i,1,j))
        flqc(i,j)   = chs(i,j) * rho(i) * wet(i,j)
        ! for output
        cmm(i)      = cm(i) * wind(i)
        chh(i)      = chs(i,j) * rho(i)
        !

        !  ---- ... outside sflx, roughness uses cm as unit
        z0(i,j)  = zorl(i)/100.
        znt(i,j) = zorl(i)/100.

        if(debug_print) then
          !if(me==0 .and. i==ipr) then
            print *,'before RUC smsoil = ',smsoil(i,:,j), i,j
            print *,'stsoil = ',stsoil(i,:,j), i,j
            print *,'soilt = ',soilt(i,j), i,j
            print *,'wet = ',wet(i,j), i,j
            print *,'soilt1 = ',soilt1(i,j), i,j
            print *,'delt =',delt
            print *,'kdt =',kdt
            print *,'flag_init =',flag_init
            print *,'flag_restart =',flag_restart
            print *,'nsoil =',nsoil
            print *,'frpcpn =',frpcpn
            print *,'zs =',zs
            print *,'graupelncv(i,j) =',i,j,graupelncv(i,j)
            print *,'snowncv(i,j) =',i,j,snowncv(i,j)
            print *,'rainncv(i,j) =',i,j,rainncv(i,j)
            print *,'raincv(i,j) =',i,j,raincv(i,j)
            print *,'prcp(i,j) =',i,j,prcp(i,j)
            print *,'sneqv(i,j) =',i,j,sneqv(i,j)
            print *,'snowh(i,j) =',i,j,snowh(i,j)
            print *,'sncovr(i,j) =',i,j,sncovr(i,j)
            print *,'ffrozp(i,j) =',i,j,ffrozp(i,j)
            print *,'conflx2(i,1,j) =',i,j,conflx2(i,1,j)
            print *,'sfcprs(i,1,j) =',i,j,sfcprs(i,1,j)
            print *,'sfctmp(i,1,j) =',i,j,sfctmp(i,1,j)
            print *,'q2(i,1,j) =',i,j,q2(i,1,j)
            print *,'qcatm(i,1,j) =',i,j,qcatm(i,1,j)
            print *,'rho2(i,1,j) =',i,j,rho2(i,1,j)
            print *,'lwdn(i,j) =',i,j,lwdn(i,j)
            print *,'solnet(i,j) =',i,j,solnet(i,j)
            print *,'sfcems(i,j) =',i,j,sfcems(i,j)
            print *,'chklowq(i,j) =',i,j,chklowq(i,j)
            print *,'chs(i,j) =',i,j,chs(i,j)
            print *,'flqc(i,j) =',i,j,flqc(i,j)
            print *,'flhc(i,j) =',i,j,flhc(i,j)
            print *,'wet(i,j) =',i,j,wet(i,j)
            print *,'cmc(i,j) =',i,j,cmc(i,j)
            print *,'shdfac(i,j) =',i,j,shdfac(i,j)
            print *,'alb(i,j) =',i,j,alb(i,j)
            print *,'znt(i,j) =',i,j,znt(i,j)
            print *,'z0(i,j) =',i,j,z0(i,j)
            print *,'snoalb1d(i,j) =',i,j,snoalb1d(i,j)
            print *,'alb(i,j) =',i,j,alb(i,j)
            print *,'landusef(i,:,j) =',i,j,landusef(i,:,j)
            print *,'soilctop(i,:,j) =',i,j,soilctop(i,:,j)
            print *,'nlcat=',nlcat
            print *,'nscat=',nscat
            print *,'qsfc(i,j) =',i,j,qsfc(i,j)
            print *,'qvg(i,j) =',i,j,qvg(i,j)
            print *,'qsg(i,j) =',i,j,qsg(i,j)
            print *,'qcg(i,j) =',i,j,qcg(i,j)
            print *,'dew(i,j) =',i,j,dew(i,j)
            print *,'soilt(i,j) =',i,j,soilt(i,j)
            print *,'tskin(i) =',i,j,tskin(i)
            print *,'soilt1(i,j) =',i,j,soilt1(i,j)
            print *,'tsnav(i,j) =',i,j,tsnav(i,j)
            print *,'tbot(i,j) =',i,j,tbot(i,j)
            print *,'vtype(i,j) =',i,j,vtype(i,j)
            print *,'stype(i,j) =',i,j,stype(i,j)
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
            print *,'hfx(i,j) =',i,j,hfx(i,j)
            print *,'qfx(i,j) =',i,j,qfx(i,j)
            print *,'lh(i,j) =',i,j,lh(i,j)
            print *,'infiltr(i,j) =',i,j,infiltr(i,j)
            print *,'runoff1(i,j) =',i,j,runoff1(i,j)
            print *,'runoff2(i,j) =',i,j,runoff2(i,j)
            print *,'acrunoff(i,j) =',i,j,acrunoff(i,j)
            print *,'sfcexc(i,j) =',i,j,sfcexc(i,j)
            print *,'acceta(i,j) =',i,j,acceta(i,j)
            print *,'ssoil(i,j) =',i,j,ssoil(i,j)
            print *,'snfallac(i,j) =',i,j,snfallac(i,j)
            print *,'acsn(i,j) =',i,j,acsn(i,j)
            print *,'snomlt(i,j) =',i,j,snomlt(i,j)
            print *,'shdmin1d(i,j) =',i,j,shdmin1d(i,j)
            print *,'shdmax1d(i,j) =',i,j,shdmax1d(i,j)
            print *,'rdlai2d =',rdlai2d
          !endif
        endif

!> - Call RUC LSM lsmruc(). 
      call lsmruc( delt, flag_init, flag_restart, kdt, iter, nsoil,          &
     &          graupelncv(i,j), snowncv(i,j), rainncv(i,j), raincv(i,j),    &
     &          zs, prcp(i,j), sneqv(i,j), snowh(i,j), sncovr(i,j),          &
     &          ffrozp(i,j), frpcpn,                                         &
     &          rhosnfr(i,j), precipfr(i,j),                                 &
!  ---  inputs:
     &          conflx2(i,1,j), sfcprs(i,1,j), sfctmp(i,1,j), q2(i,1,j),     &
     &          qcatm(i,1,j), rho2(i,1,j),                                   &
     &          lwdn(i,j), solnet(i,j), sfcems(i,j), chklowq(i,j),           &
     &          chs(i,j), flqc(i,j), flhc(i,j),                              &
!  ---  input/outputs:
     &          wet(i,j), cmc(i,j), shdfac(i,j), alb(i,j), znt(i,j),         &
     &          z0(i,j), snoalb1d(i,j), albbck(i,j),                         &
!     &          z0, snoalb1d, alb, xlai,                                     &
     &          landusef(i,:,j), nlcat,                                      &
!  --- mosaic_lu and mosaic_soil are moved to the namelist
!     &          mosaic_lu, mosaic_soil,                                      &
     &          soilctop(i,:,j), nscat,                                      &
     &          qsfc(i,j), qsg(i,j), qvg(i,j), qcg(i,j), dew(i,j),           &
     &          soilt1(i,j),                                                 &
     &          tsnav(i,j), tbot(i,j), vtype(i,j), stype(i,j), xland(i,j),   &
     &          iswater, isice, xice(i,j), xice_threshold,                   &
!  ---  constants
     &          con_cp, con_rv, con_rd, con_g, con_pi, con_hvap, stbolt,     &
!  ---  input/outputs:
     &          smsoil(i,:,j), slsoil(i,:,j), soilm(i,j), smmax(i,j),        &
     &          stsoil(i,:,j), soilt(i,j), hfx(i,j), qfx(i,j), lh(i,j),      &
     &          infiltr(i,j), runoff1(i,j), runoff2(i,j), acrunoff(i,j),     &
     &          sfcexc(i,j), acceta(i,j), ssoil(i,j),                        &
     &          snfallac(i,j), acsn(i,j), snomlt(i,j),                       &
     &          smfrsoil(i,:,j),keepfrsoil(i,:,j), .false.,                  &
     &          shdmin1d(i,j), shdmax1d(i,j), rdlai2d,                       &
     &          ims,ime, jms,jme, kms,kme,                                   &
     &          its,ite, jts,jte, kts,kte                                    )

        if(debug_print) then
          print *,'after sneqv(i,j) =',i,j,sneqv(i,j)
          print *,'after snowh(i,j) =',i,j,snowh(i,j)
          print *,'after sncovr(i,j) =',i,j,sncovr(i,j)
          print *,'after vtype(i,j) =',i,j,vtype(i,j)
          print *,'after stype(i,j) =',i,j,stype(i,j)
          print *,'after wet(i,j) =',i,j,wet(i,j)
          print *,'after cmc(i,j) =',i,j,cmc(i,j)
          print *,'after qsfc(i,j) =',i,j,qsfc(i,j)
          print *,'after qvg(i,j) =',i,j,qvg(i,j)
          print *,'after qsg(i,j) =',i,j,qsg(i,j)
          print *,'after qcg(i,j) =',i,j,qcg(i,j)
          print *,'after dew(i,j) =',i,j,dew(i,j)
          print *,'after soilt(i,j) =',i,j,soilt(i,j)
          print *,'after tskin(i) =',i,j,tskin(i)
          print *,'after soilt1(i,j) =',i,j,soilt1(i,j)
          print *,'after tsnav(i,j) =',i,j,tsnav(i,j)
          print *,'after smsoil(i,:,j)=',i,j,smsoil(i,:,j)
          print *,'after slsoil(i,:,j)=',i,j,slsoil(i,:,j)
          print *,'after stsoil(i,:,j)=',i,j,stsoil(i,:,j)
          print *,'after smfrsoil(i,:,j)=',i,j,smfrsoil(i,:,j)
          print *,'after keepfrsoil(i,:,j)=',i,j,keepfrsoil(i,:,j)
          print *,'after soilm(i,j) =',i,j,soilm(i,j)
          print *,'after smmax(i,j) =',i,j,smmax(i,j)
          print *,'after hfx(i,j) =',i,j,hfx(i,j)
          print *,'after qfx(i,j) =',i,j,qfx(i,j)
          print *,'after lh(i,j) =',i,j,lh(i,j)
          print *,'after infiltr(i,j) =',i,j,infiltr(i,j)
          print *,'after runoff1(i,j) =',i,j,runoff1(i,j)
          print *,'after runoff2(i,j) =',i,j,runoff2(i,j)
          print *,'after ssoil(i,j) =',i,j,ssoil(i,j)
          print *,'after snfallac(i,j) =',i,j,snfallac(i,j)
          print *,'after acsn(i,j) =',i,j,acsn(i,j)
          print *,'after snomlt(i,j) =',i,j,snomlt(i,j)
        endif


!> - RUC LSM: prepare variables for return to parent model and unit conversion.
!>  -   6. output (o):
!!\n  ------------------------------
!!\n \a lh     - actual latent heat flux (\f$W m^{-2}\f$: positive, if upward from sfc)
!!\n \a hfx    - sensible heat flux (\f$W m^{-2}\f$: positive, if upward from sfc)
!!\n \a ssoil   - soil heat flux (\f$W m^{-2}\f$: negative if downward from surface)
!!\n \a runoff1 - surface runoff (\f$m s^{-1}\f$), not infiltrating the surface
!!\n \a runoff2 - subsurface runoff (\f$m s^{-1}\f$), drainage out bottom
!!\n \a snoh    - phase-change heat flux from snowmelt (w m-2)
!
        if(debug_print) then
          !if(me==0.and.i==ipr) then
            print *,'after  RUC smsoil = ',smsoil(i,:,j), i, j
            print *,'stsoil = ',stsoil(i,:,j), i,j
            print *,'soilt = ',soilt(i,j), i,j
            print *,'wet = ',wet(i,j), i,j
            print *,'soilt1 = ',soilt1(i,j), i,j
            print *,'rhosnfr = ',rhosnfr(i,j), i,j
          !endif
        endif

        ! Interstitial
        evap(i)   = qfx(i,j) / rho(i)           ! kinematic
        hflx(i)   = hfx(i,j) / (con_cp*rho(i))  ! kinematic
        gflux(i)  = ssoil(i,j)

        !evbs(i)  = edir(i,j)
        !evcw(i)  = ec(i,j)
        !trans(i) = ett(i,j)
        !sbsno(i) = esnow(i,j)
        !snohf(i) = snoh(i,j)

        sfcdew(i) = dew(i,j)
        qsurf(i)  = qsfc(i,j)
        snowc(i)  = sncovr(i,j)
        stm(i)    = soilm(i,j)
        tsurf(i)  = soilt(i,j)
        tice(i)   = tsurf(i)
        !  --- ...  units [m/s] = [g m-2 s-1] 
        runof (i)  = runoff1(i,j)
        drain (i)  = runoff2(i,j)

        wet1(i) = wet(i,j)

        ! State variables
        tsnow(i)   = soilt1(i,j)
        sfcqc(i)  = qcg(i,j)
        sfcqv(i)  = qvg(i,j)
        rhosnf(i) = rhosnfr(i,j)

        ! --- ... accumulated total runoff and surface runoff
        runoff(i)  = runoff(i)  + (drain(i)+runof(i)) * delt * 0.001 ! kg m-2
        srunoff(i) = srunoff(i) + runof(i) * delt * 0.001            ! kg m-2

        !  --- ...  unit conversion (from m to mm)
        snwdph(i)  = snowh(i,j) * 1000.0

        canopy(i)  = cmc(i,j)   ! mm
        weasd(i)   = sneqv(i,j) ! mm
        sncovr1(i) = sncovr(i,j)
        !  ---- ... outside RUC LSM, roughness uses cm as unit 
        !  (update after snow's effect)
        zorl(i) = znt(i,j)*100.
        sfalb(i)= alb(i,j)

        do k = 1, lsoil_ruc
          smois(i,k)  = smsoil(i,k,j)
          sh2o(i,k)   = slsoil(i,k,j)
          tslb(i,k)   = stsoil(i,k,j)
          keepfr(i,k)   = keepfrsoil(i,k,j)
          smfrkeep(i,k) = smfrsoil(i,k,j)
        enddo

        do k = 1, lsoil
          smc(i,k)   = smsoil(i,k,j)
          slc(i,k)   = slsoil(i,k,j)
          stc(i,k)   = stsoil(i,k,j)
        enddo

!  --- ...  do not return the following output fields to parent model
!    ec      - canopy water evaporation (m s-1)
!    edir    - direct soil evaporation (m s-1)
!    et(nsoil)-plant transpiration from a particular root layer (m s-1)
!    ett     - total plant transpiration (m s-1)
!    esnow   - sublimation from (or deposition to if <0) snowpack (m s-1)
!    drip    - through-fall of precip and/or dew in excess of canopy
!              water-holding capacity (m)
!    dew     - dewfall (or frostfall for t<273.15) (m)
!    snomlt  - snow melt (m) (water equivalent)
!    sncovr  - fractional snow cover (unitless fraction, 0-1)
!              for a given soil layer at the end of a time step
!    xlai    - leaf area index (dimensionless)
!    soilw   - available soil moisture in root zone (unitless fraction
!              between smcwlt and smcmax)
!    soilm   - total soil column moisture content (frozen+unfrozen) (m)
!    nroot   - number of root layers, a function of veg type, determined
!              in subroutine redprm.

        endif   ! end if_flag_iter_and_flag_block
      enddo   ! j
      enddo   ! i

!> - Restore land-related prognostic fields for guess run.
      do j  = 1, 1
      do i  = 1, im
        if (flag(i)) then
          if(debug_print) print *,'end ',i,flag_guess(i),flag_iter(i)
          if (flag_guess(i)) then
            if(debug_print) print *,'guess run'
            weasd(i)  = weasd_old(i)
            snwdph(i) = snwdph_old(i)
            tskin(i)  = tskin_old(i)
            canopy(i) = canopy_old(i)
            !tprcp(i)  = tprcp_old(i)
            srflag(i) = srflag_old(i)

            do k = 1, lsoil_ruc
              smois(i,k) = smois_old(i,k)
              tslb(i,k) = tslb_old(i,k)
              sh2o(i,k) = sh2o_old(i,k)
              keepfr(i,k) = keepfr_old(i,k)
              smfrkeep(i,k) = smfrkeep_old(i,k)
            enddo
          else
            if(debug_print) print *,'iter run', i,j, tskin(i),tsurf(i)
            tskin(i) = tsurf(i)
            tice (i) = tsurf(i)
          endif
        endif
      enddo  ! i
      enddo  ! j
!
      deallocate(soilctop)
      deallocate(landusef)
!
      return
!...................................
      end subroutine lsm_ruc_run
!-----------------------------------
      subroutine rucinit      (restart, im, lsoil_ruc, lsoil, nlev,   & ! in
                               isot, soiltyp, vegtype, fice,          & ! in
                               islmsk, tsurf, tg3,                    & ! in
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
      integer,               dimension(im),    intent(in   ) :: islmsk
      real (kind=kind_phys), dimension(im),    intent(in   ) :: tsurf
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

          tsk(i,j) = tsurf(i)
          tbot(i,j)=tg3(i)

          !SLMSK   - SEA(0),LAND(1),ICE(2) MASK
          if(islmsk(i) == 0) then
            ivgtyp(i,j)= 17 ! 17 - water (oceans and lakes) in MODIS
            isltyp(i,j)=14
            xice(i,j)=0.
            landmask(i,j)=0.
          elseif(islmsk(i) == 1) then ! land
            ivgtyp(i,j)=vegtype(i)
            isltyp(i,j)=soiltyp(i)
            landmask(i,j)=1.
            xice(i,j)=0.
          elseif(islmsk(i) == 2) then  ! ice
            ivgtyp(i,j)=15 ! MODIS
            !> -- number of soil categories
            if(isot == 1) then
              isltyp(i,j) = 16 ! STATSGO
            else
              isltyp(i,j) = 9  ! ZOBLER
            endif
            landmask(i,j)=1.
            xice(i,j)=fice(i)
          endif

          sst(i,j) = tsk(i,j)

          st_input(i,1,j)=tsk(i,j)
          sm_input(i,1,j)=0.

          !--- initialize smcwlt2 and smcref2 with Noah values
          if(islmsk(i) == 0 .or. islmsk(i) == 2) then
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
             if(islmsk(i) == 1) then
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

          IF ( islmsk(i) == 1 ) then  ! Land
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

!> \file GFS_surface_generic.f90
!!  Contains code related to all GFS surface schemes.

      module GFS_surface_generic_pre

      contains

      subroutine GFS_surface_generic_pre_init ()
      end subroutine GFS_surface_generic_pre_init

      subroutine GFS_surface_generic_pre_finalize()
      end subroutine GFS_surface_generic_pre_finalize

!> \section arg_table_GFS_surface_generic_pre_run Argument Table
!! | local_name     | standard_name                                                                | long_name                                                        | units      | rank | type      |    kind   | intent | optional |
!! |----------------|------------------------------------------------------------------------------|------------------------------------------------------------------|------------|------|-----------|-----------|--------|----------|
!! | im             | horizontal_loop_extent                                                       | horizontal loop extent                                           | count      |    0 | integer   |           | in     | F        |
!! | levs           | vertical_dimension                                                           | number of vertical levels                                        | count      |    0 | integer   |           | in     | F        |
!! | vfrac          | vegetation_area_fraction                                                     | areal fractional cover of green vegetation                       | frac       |    1 | real      | kind_phys | in     | F        |
!! | islmsk         | sea_land_ice_mask                                                            | landmask: sea/land/ice=0/1/2                                     | flag       |    1 | integer   |           | in     | F        |
!! | isot           | soil_type_dataset_choice                                                     | soil type dataset choice                                         | index      |    0 | integer   |           | in     | F        |
!! | ivegsrc        | vegetation_type_dataset_choice                                               | land use dataset choice                                          | index      |    0 | integer   |           | in     | F        |
!! | stype          | soil_type_classification_real                                                | soil type for lsm                                                | index      |    1 | real      | kind_phys | in     | F        |
!! | vtype          | vegetation_type_classification_real                                          | vegetation type for lsm                                          | index      |    1 | real      | kind_phys | in     | F        |
!! | slope          | surface_slope_classification_real                                            | sfc slope type for lsm                                           | index      |    1 | real      | kind_phys | in     | F        |
!! | prsik_1        | dimensionless_exner_function_at_lowest_model_interface                       | dimensionless Exner function at lowest model interface           | none       |    1 | real      | kind_phys | in     | F        |
!! | prslk_1        | dimensionless_exner_function_at_lowest_model_layer                           | dimensionless Exner function at lowest model layer               | none       |    1 | real      | kind_phys | in     | F        |
!! | semis          | surface_longwave_emissivity                                                  | surface lw emissivity in fraction                                | frac       |    1 | real      | kind_phys | in     | F        |
!! | adjsfcdlw      | surface_downwelling_longwave_flux                                            | surface downwelling longwave flux at current time                | W m-2      |    1 | real      | kind_phys | in     | F        |
!! | tsfc           | surface_skin_temperature                                                     | surface skin temperature                                         | K          |    1 | real      | kind_phys | in     | F        |
!! | phil           | geopotential                                                                 | geopotential at model layer centers                              | m2 s-2     |    2 | real      | kind_phys | in     | F        |
!! | con_g          | gravitational_acceleration                                                   | gravitational acceleration                                       | m s-2      |    0 | real      | kind_phys | in     | F        |
!! | sigmaf         | bounded_vegetation_area_fraction                                             | areal fractional cover of green vegetation bounded on the bottom | frac       |    1 | real      | kind_phys | inout  | F        |
!! | soiltyp        | soil_type_classification                                                     | soil type at each grid cell                                      | index      |    1 | integer   |           | inout  | F        |
!! | vegtype        | vegetation_type_classification                                               | vegetation type at each grid cell                                | index      |    1 | integer   |           | inout  | F        |
!! | slopetyp       | surface_slope_classification                                                 | surface slope type at each grid cell                             | index      |    1 | integer   |           | inout  | F        |
!! | work3          | ratio_of_exner_function_between_midlayer_and_interface_at_lowest_model_layer | Exner function ratio bt midlayer and interface at 1st layer      | ratio      |    1 | real      | kind_phys | inout  | F        |
!! | gabsbdlw       | surface_downwelling_longwave_flux_absorbed_by_ground                         | total sky surface downward longwave flux absorbed by the ground  | W m-2      |    1 | real      | kind_phys | inout  | F        |
!! | tsurf          | surface_skin_temperature_after_iteration                                     | surface skin temperature after iteration                         | K          |    1 | real      | kind_phys | inout  | F        |
!! | zlvl           | height_above_ground_at_lowest_model_layer                                    | layer 1 height above ground (not MSL)                            | m          |    1 | real      | kind_phys | inout  | F        |
!! | errmsg         | ccpp_error_message                                                           | error message for error handling in CCPP                         | none       |    0 | character | len=*     | out    | F        |
!! | errflg         | ccpp_error_flag                                                              | error flag for error handling in CCPP                            | flag       |    0 | integer   |           | out    | F        |
!!
      subroutine GFS_surface_generic_pre_run (im, levs, vfrac, islmsk, isot, ivegsrc, stype, vtype, slope, &
        prsik_1, prslk_1, semis, adjsfcdlw, tsfc, phil, con_g, sigmaf, soiltyp, vegtype, slopetyp, work3, gabsbdlw, &
        tsurf, zlvl, errmsg, errflg)

        use machine,               only: kind_phys

        implicit none

        integer, intent(in) :: im, levs, isot, ivegsrc
        integer, dimension(im), intent(in) :: islmsk
        integer, dimension(im), intent(inout) :: soiltyp, vegtype, slopetyp

        real(kind=kind_phys), intent(in) :: con_g
        real(kind=kind_phys), dimension(im), intent(in) :: vfrac, stype, vtype, slope, prsik_1, prslk_1, &
          semis, adjsfcdlw, tsfc
        real(kind=kind_phys), dimension(im,levs), intent(in) :: phil

        real(kind=kind_phys), dimension(im), intent(inout) :: sigmaf, work3, gabsbdlw, tsurf, zlvl

        character(len=*), intent(out) :: errmsg
        integer,          intent(out) :: errflg

        integer :: i
        real(kind=kind_phys) :: onebg

        onebg  = 1.0/con_g

        ! Initialize CCPP error handling variables
        errmsg = ''
        errflg = 0

        do i=1,im
          sigmaf(i) = max(vfrac(i),0.01 )
          if (islmsk(i) == 2) then
            if (isot == 1) then
              soiltyp(i)  = 16
            else
              soiltyp(i)  = 9
            endif
            if (ivegsrc == 1) then
              vegtype(i)  = 15
            elseif(ivegsrc == 2) then
              vegtype(i)  = 13
            endif
            slopetyp(i) = 9
          else
            soiltyp(i)  = int( stype(i)+0.5 )
            vegtype(i)  = int( vtype(i)+0.5 )
            slopetyp(i) = int( slope(i)+0.5 )    !! clu: slope -> slopetyp
          endif

          work3(i)   = prsik_1(i) / prslk_1(i)
        end do

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
        gabsbdlw(:) = semis(:) * adjsfcdlw(:)

        do i=1,im
          tsurf(i)        = tsfc(i)
          zlvl(i)    = phil(i,1) * onebg
        end do

      end subroutine GFS_surface_generic_pre_run

      end module GFS_surface_generic_pre

      module GFS_surface_generic_post

      contains

      subroutine GFS_surface_generic_post_init ()
      end subroutine GFS_surface_generic_post_init

      subroutine GFS_surface_generic_post_finalize()
      end subroutine GFS_surface_generic_post_finalize
#if 0
!> \section arg_table_GFS_surface_generic_post_run Argument Table
!! | local_name     | standard_name                                                                                                       | long_name                                                                           | units       | rank | type       |    kind   | intent | optional |
!! |----------------|---------------------------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------|-------------|------|------------|-----------|--------|----------|
!! | im             | horizontal_loop_extent                                                                                              | horizontal loop extent                                                              | count       |    0 | integer    |           | in     | F        |
!! | cplflx         | flag_for_flux_coupling                                                                                              | flag controlling cplflx collection (default off)                                    | flag        |    0 | logical    |           | in     | F        |
!! | lssav          | flag_diagnostics                                                                                                    | logical flag for storing diagnostics                                                | flag        |    0 | logical    |           | in     | F        |
!! | islmsk         | sea_land_ice_mask                                                                                                   | landmask: sea/land/ice=0/1/2                                                        | flag        |    1 | integer    |           | in     | F        |
!! | dtf            | time_step_for_dynamics                                                                                              | dynamics timestep                                                                   | s           |    0 | real       | kind_phys | in     | F        |
!! | ep1d           | surface_upward_potential_latent_heat_flux                                                                           | surface upward potential latent heat flux                                           | W m-2       |    1 | real       | kind_phys | in     | F        |
!! | gflx           | upward_heat_flux_in_soil                                                                                            | upward soil heat flux                                                               | W m-2       |    1 | real       | kind_phys | in     | F        |
!! | tgrs_1         | air_temperature_at_lowest_model_layer                                                                               | mean temperature at lowest model layer                                              | K           |    1 | real       | kind_phys | in     | F        |
!! | qgrs_1         | specific_humidity_at_lowest_model_layer                                                                             | specific humidity at lowest model layer                                             | kg kg-1     |    1 | real       | kind_phys | in     | F        |
!! | ugrs_1         | x_wind_at_lowest_model_layer                                                                                        | zonal wind at lowest model layer                                                    | m s-1       |    1 | real       | kind_phys | in     | F        |
!! | vgrs_1         | y_wind_at_lowest_model_layer                                                                                        | meridional wind at lowest model layer                                               | m s-1       |    1 | real       | kind_phys | in     | F        |
!! | adjsfcdlw      | surface_downwelling_longwave_flux                                                                                   | surface downwelling longwave flux at current time                                   | W m-2       |    1 | real       | kind_phys | in     | F        |
!! | adjsfcdsw      | surface_downwelling_shortwave_flux                                                                                  | surface downwelling shortwave flux at current time                                  | W m-2       |    1 | real       | kind_phys | in     | F        |
!! | adjnirbmd      | surface_downwelling_direct_near_infrared_shortwave_flux                                                             | surface downwelling beam near-infrared shortwave flux at current time               | W m-2       |    1 | real       | kind_phys | in     | F        |
!! | adjnirdfd      | surface_downwelling_diffuse_near_infrared_shortwave_flux                                                            | surface downwelling diffuse near-infrared shortwave flux at current time            | W m-2       |    1 | real       | kind_phys | in     | F        |
!! | adjvisbmd      | surface_downwelling_direct_ultraviolet_and_visible_shortwave_flux                                                   | surface downwelling beam ultraviolet plus visible shortwave flux at current time    | W m-2       |    1 | real       | kind_phys | in     | F        |
!! | adjvisdfd      | surface_downwelling_diffuse_ultraviolet_and_visible_shortwave_flux                                                  | surface downwelling diffuse ultraviolet plus visible shortwave flux at current time | W m-2       |    1 | real       | kind_phys | in     | F        |
!! | adjsfculw      | surface_upwelling_longwave_flux                                                                                     | surface upwelling longwave flux at current time                                     | W m-2       |    1 | real       | kind_phys | in     | F        |
!! | adjnirbmu      | surface_upwelling_direct_near_infrared_shortwave_flux                                                               | surface upwelling beam near-infrared shortwave flux at current time                 | W m-2       |    1 | real       | kind_phys | in     | F        |
!! | adjnirdfu      | surface_upwelling_diffuse_near_infrared_shortwave_flux                                                              | surface upwelling diffuse near-infrared shortwave flux at current time              | W m-2       |    1 | real       | kind_phys | in     | F        |
!! | adjvisbmu      | surface_upwelling_direct_ultraviolet_and_visible_shortwave_flux                                                     | surface upwelling beam ultraviolet plus visible shortwave flux at current time      | W m-2       |    1 | real       | kind_phys | in     | F        |
!! | adjvisdfu      | surface_upwelling_diffuse_ultraviolet_and_visible_shortwave_flux                                                    | surface upwelling diffuse ultraviolet plus visible shortwave flux at current time   | W m-2       |    1 | real       | kind_phys | in     | F        |
!! | t2m            | temperature_at_2m                                                                                                   | 2 meter temperature                                                                 | K           |    1 | real       | kind_phys | in     | F        |
!! | q2m            | specific_humidity_at_2m                                                                                             | 2 meter specific humidity                                                           | kg kg-1     |    1 | real       | kind_phys | in     | F        |
!! | u10m           | x_wind_at_10m                                                                                                       | 10 meter u wind speed                                                               | m s-1       |    1 | real       | kind_phys | in     | F        |
!! | v10m           | y_wind_at_10m                                                                                                       | 10 meter v wind speed                                                               | m s-1       |    1 | real       | kind_phys | in     | F        |
!! | tsfc           | surface_skin_temperature                                                                                            | surface skin temperature                                                            | K           |    1 | real       | kind_phys | in     | F        |
!! | pgr            | surface_air_pressure                                                                                                | surface pressure                                                                    | Pa          |    1 | real       | kind_phys | in     | F        |
!! | xcosz          | instantaneous_cosine_of_zenith_angle                                                                                | cosine of zenith angle at current time                                              | none        |    1 | real       | kind_phys | in     | F        |
!! | evbs           | soil_upward_latent_heat_flux                                                                                        | soil upward latent heat flux                                                        | W m-2       |    1 | real       | kind_phys | in     | F        |
!! | evcw           | canopy_upward_latent_heat_flux                                                                                      | canopy upward latent heat flux                                                      | W m-2       |    1 | real       | kind_phys | in     | F        |
!! | trans          | transpiration_flux                                                                                                  | total plant transpiration rate                                                      | kg m-2 s-1  |    1 | real       | kind_phys | in     | F        |
!! | sbsno          | snow_deposition_sublimation_upward_latent_heat_flux                                                                 | latent heat flux from snow depo/subl                                                | W m-2       |    1 | real       | kind_phys | in     | F        |
!! | snowc          | surface_snow_area_fraction                                                                                          | surface snow area fraction                                                          | frac        |    1 | real       | kind_phys | in     | F        |
!! | snohf          | snow_freezing_rain_upward_latent_heat_flux                                                                          | latent heat flux due to snow and frz rain                                           | W m-2       |    1 | real       | kind_phys | in     | F        |
!! | con_eps        | ratio_of_dry_air_to_water_vapor_gas_constants                                                                       | rd/rv                                                                               | none        |    0 | real       | kind_phys | in     | F        |
!! | con_epsm1      | ratio_of_dry_air_to_water_vapor_gas_constants_minus_one                                                             | (rd/rv) - 1                                                                         | none        |    0 | real       | kind_phys | in     | F        |
!! | epi            | instantaneous_surface_potential_evaporation                                                                         | instantaneous sfc potential evaporation                                             | W m-2       |    1 | real       | kind_phys | inout  | F        |
!! | gfluxi         | instantaneous_surface_ground_heat_flux                                                                              | instantaneous sfc ground heat flux                                                  | W m-2       |    1 | real       | kind_phys | inout  | F        |
!! | t1             | air_temperature_at_lowest_model_layer_for_diag                                                                      | layer 1 temperature for diag                                                        | K           |    1 | real       | kind_phys | inout  | F        |
!! | q1             | specific_humidity_at_lowest_model_layer_for_diag                                                                    | layer 1 specific humidity for diag                                                  | kg kg-1     |    1 | real       | kind_phys | inout  | F        |
!! | u1             | x_wind_at_lowest_model_layer_for_diag                                                                               | layer 1 x wind for diag                                                             | m s-1       |    1 | real       | kind_phys | inout  | F        |
!! | v1             | y_wind_at_lowest_model_layer_for_diag                                                                               | layer 1 y wind for diag                                                             | m s-1       |    1 | real       | kind_phys | inout  | F        |
!! | dlwsfci_cpl    | instantaneous_surface_downwelling_longwave_flux_for_coupling                                                        | instantaneous sfc downward lw flux                                                  | W m-2       |    1 | real       | kind_phys | inout  | F        |
!! | dswsfci_cpl    | instantaneous_surface_downwelling_shortwave_flux_for_coupling                                                       | instantaneous sfc downward sw flux                                                  | W m-2       |    1 | real       | kind_phys | inout  | F        |
!! | dlwsfc_cpl     | cumulative_surface_downwelling_longwave_flux_for_coupling_multiplied_by_timestep                                    | cumulative sfc downward lw flux mulitplied by timestep                              | W m-2 s     |    1 | real       | kind_phys | inout  | F        |
!! | dswsfc_cpl     | cumulative_surface_downwelling_shortwave_flux_for_coupling_multiplied_by_timestep                                   | cumulative sfc downward sw flux multiplied by timestep                              | W m-2 s     |    1 | real       | kind_phys | inout  | F        |
!! | dnirbmi_cpl    | instantaneous_surface_downwelling_direct_near_infrared_shortwave_flux_for_coupling                                  | instantaneous sfc nir beam downward sw flux                                         | W m-2       |    1 | real       | kind_phys | inout  | F        |
!! | dnirdfi_cpl    | instantaneous_surface_downwelling_diffuse_near_infrared_shortwave_flux_for_coupling                                 | instantaneous sfc nir diff downward sw flux                                         | W m-2       |    1 | real       | kind_phys | inout  | F        |
!! | dvisbmi_cpl    | instantaneous_surface_downwelling_direct_ultraviolet_and_visible_shortwave_flux_for_coupling                        | instantaneous sfc uv+vis beam downward sw flux                                      | W m-2       |    1 | real       | kind_phys | inout  | F        |
!! | dvisdfi_cpl    | instantaneous_surface_downwelling_diffuse_ultraviolet_and_visible_shortwave_flux_for_coupling                       | instantaneous sfc uv+vis diff downward sw flux                                      | W m-2       |    1 | real       | kind_phys | inout  | F        |
!! | dnirbm_cpl     | cumulative_surface_downwelling_direct_near_infrared_shortwave_flux_for_coupling_multiplied_by_timestep              | cumulative sfc nir beam downward sw flux multiplied by timestep                     | W m-2 s     |    1 | real       | kind_phys | inout  | F        |
!! | dnirdf_cpl     | cumulative_surface_downwelling_diffuse_near_infrared_shortwave_flux_for_coupling_multiplied_by_timestep             | cumulative sfc nir diff downward sw flux multiplied by timestep                     | W m-2 s     |    1 | real       | kind_phys | inout  | F        |
!! | dvisbm_cpl     | cumulative_surface_downwelling_direct_ultraviolet_and_visible_shortwave_flux_for_coupling_multiplied_by_timestep    | cumulative sfc uv+vis beam dnwd sw flux multiplied by timestep                      | W m-2 s     |    1 | real       | kind_phys | inout  | F        |
!! | dvisdf_cpl     | cumulative_surface_downwelling_diffuse_ultraviolet_and_visible_shortwave_flux_for_coupling_multiplied_by_timestep   | cumulative sfc uv+vis diff dnwd sw flux multiplied by timestep                      | W m-2 s     |    1 | real       | kind_phys | inout  | F        |
!! | nlwsfci_cpl    | instantaneous_surface_net_downward_longwave_flux_for_coupling                                                       | instantaneous net sfc downward lw flux                                              | W m-2       |    1 | real       | kind_phys | inout  | F        |
!! | nlwsfc_cpl     | cumulative_surface_net_downward_longwave_flux_for_coupling_multiplied_by_timestep                                   | cumulative net downward lw flux multiplied by timestep                              | W m-2 s     |    1 | real       | kind_phys | inout  | F        |
!! | t2mi_cpl       | instantaneous_temperature_at_2m_for_coupling                                                                        | instantaneous T2m                                                                   | K           |    1 | real       | kind_phys | inout  | F        |
!! | q2mi_cpl       | instantaneous_specific_humidity_at_2m_for_coupling                                                                  | instantaneous Q2m                                                                   | kg kg-1     |    1 | real       | kind_phys | inout  | F        |
!! | u10mi_cpl      | instantaneous_x_wind_at_10m_for_coupling                                                                            | instantaneous U10m                                                                  | m s-1       |    1 | real       | kind_phys | inout  | F        |
!! | v10mi_cpl      | instantaneous_y_wind_at_10m_for_coupling                                                                            | instantaneous V10m                                                                  | m s-1       |    1 | real       | kind_phys | inout  | F        |
!! | tsfci_cpl      | instantaneous_surface_skin_temperature_for_coupling                                                                 | instantaneous sfc temperature                                                       | K           |    1 | real       | kind_phys | inout  | F        |
!! | psurfi_cpl     | instantaneous_surface_air_pressure_for_coupling                                                                     | instantaneous sfc pressure                                                          | Pa          |    1 | real       | kind_phys | inout  | F        |
!! | nnirbmi_cpl    | instantaneous_surface_net_downward_direct_near_infrared_shortwave_flux_for_coupling                                 | instantaneous net nir beam sfc downward sw flux                                     | W m-2       |    1 | real       | kind_phys | inout  | F        |
!! | nnirdfi_cpl    | instantaneous_surface_net_downward_diffuse_near_infrared_shortwave_flux_for_coupling                                | instantaneous net nir diff sfc downward sw flux                                     | W m-2       |    1 | real       | kind_phys | inout  | F        |
!! | nvisbmi_cpl    | instantaneous_surface_net_downward_direct_ultraviolet_and_visible_shortwave_flux_for_coupling                       | instantaneous net uv+vis beam downward sw flux                                      | W m-2       |    1 | real       | kind_phys | inout  | F        |
!! | nvisdfi_cpl    | instantaneous_surface_net_downward_diffuse_ultraviolet_and_visible_shortwave_flux_for_coupling                      | instantaneous net uv+vis diff downward sw flux                                      | W m-2       |    1 | real       | kind_phys | inout  | F        |
!! | nswsfci_cpl    | instantaneous_surface_net_downward_shortwave_flux_for_coupling                                                      | instantaneous net sfc downward sw flux                                              | W m-2       |    1 | real       | kind_phys | inout  | F        |
!! | nswsfc_cpl     | cumulative_surface_net_downward_shortwave_flux_for_coupling_multiplied_by_timestep                                  | cumulative net downward sw flux multiplied by timestep                              | W m-2 s     |    1 | real       | kind_phys | inout  | F        |
!! | nnirbm_cpl     | cumulative_surface_net_downward_direct_near_infrared_shortwave_flux_for_coupling_multiplied_by_timestep             | cumulative net nir beam downward sw flux multiplied by timestep                     | W m-2 s     |    1 | real       | kind_phys | inout  | F        |
!! | nnirdf_cpl     | cumulative_surface_net_downward_diffuse_near_infrared_shortwave_flux_for_coupling_multiplied_by_timestep            | cumulative net nir diff downward sw flux multiplied by timestep                     | W m-2 s     |    1 | real       | kind_phys | inout  | F        |
!! | nvisbm_cpl     | cumulative_surface_net_downward_direct_ultraviolet_and_visible_shortwave_flux_for_coupling_multiplied_by_timestep   | cumulative net uv+vis beam downward sw rad flux multiplied by timestep              | W m-2 s     |    1 | real       | kind_phys | inout  | F        |
!! | nvisdf_cpl     | cumulative_surface_net_downward_diffuse_ultraviolet_and_visible_shortwave_flux_for_coupling_multiplied_by_timestep  | cumulative net uv+vis diff downward sw rad flux multiplied by timestep              | W m-2 s     |    1 | real       | kind_phys | inout  | F        |
!! | gflux          | cumulative_surface_ground_heat_flux_multiplied_by_timestep                                                          | cumulative groud conductive heat flux multiplied by timestep                        | W m-2 s     |    1 | real       | kind_phys | inout  | F        |
!! | evbsa          | cumulative_soil_upward_latent_heat_flux_multiplied_by_timestep                                                      | cumulative soil upward latent heat flux multiplied by timestep                      | W m-2 s     |    1 | real       | kind_phys | inout  | F        |
!! | evcwa          | cumulative_canopy_upward_latent_heat_flu_multiplied_by_timestep                                                     | cumulative canopy upward latent heat flux multiplied by timestep                    | W m-2 s     |    1 | real       | kind_phys | inout  | F        |
!! | transa         | cumulative_transpiration_flux_multiplied_by_timestep                                                                | cumulative total plant transpiration rate multiplied by timestep                    | kg m-2      |    1 | real       | kind_phys | inout  | F        |
!! | sbsnoa         | cumulative_snow_deposition_sublimation_upward_latent_heat_flux_multiplied_by_timestep                               | cumulative latent heat flux from snow depo/subl multiplied by timestep              | W m-2 s     |    1 | real       | kind_phys | inout  | F        |
!! | snowca         | cumulative_surface_snow_area_fraction_multiplied_by_timestep                                                        | cumulative surface snow area fraction multiplied by timestep                        | s           |    1 | real       | kind_phys | inout  | F        |
!! | snohfa         | cumulative_snow_freezing_rain_upward_latent_heat_flux_multiplied_by_timestep                                        | cumulative latent heat flux due to snow and frz rain multiplied by timestep         | W m-2 s     |    1 | real       | kind_phys | inout  | F        |
!! | ep             | cumulative_surface_upward_potential_latent_heat_flux_multiplied_by_timestep                                         | cumulative surface upward potential latent heat flux multiplied by timestep         | W m-2 s     |    1 | real       | kind_phys | inout  | F        |
!! | tmpmin         | minimum_temperature_at_2m                                                                                           | min temperature at 2m height                                                        | K           |    1 | real       | kind_phys | inout  | F        |
!! | tmpmax         | maximum_temperature_at_2m                                                                                           | max temperature at 2m height                                                        | K           |    1 | real       | kind_phys | inout  | F        |
!! | spfhmin        | minimum_specific_humidity_at_2m                                                                                     | minimum specific humidity at 2m height                                              | kg kg-1     |    1 | real       | kind_phys | inout  | F        |
!! | spfhmax        | maximum_specific_humidity_at_2m                                                                                     | maximum specific humidity at 2m height                                              | kg kg-1     |    1 | real       | kind_phys | inout  | F        |
!! | wind10mmax     | maximum_wind_at_10m                                                                                                 | maximum wind speed at 10 m                                                          | m s-1       |    1 | real       | kind_phys | inout  | F        |
!! | u10mmax        | maximum_x_wind_at_10m                                                                                               | maximum x wind at 10 m                                                              | m s-1       |    1 | real       | kind_phys | inout  | F        |
!! | v10mmax        | maximum_y_wind_at_10m                                                                                               | maximum y wind at 10 m                                                              | m s-1       |    1 | real       | kind_phys | inout  | F        |
!! | dpt2m          | dewpoint_temperature_at_2m                                                                                          | 2 meter dewpoint temperature                                                        | K           |    1 | real       | kind_phys | inout  | F        |
!! | errmsg         | ccpp_error_message                                                                                                  | error message for error handling in CCPP                                            | none        |    0 | character  | len=*     | out    | F        |
!! | errflg         | ccpp_error_flag                                                                                                     | error flag for error handling in CCPP                                               | flag        |    0 | integer    |           | out    | F        |
!!
#endif
      subroutine GFS_surface_generic_post_run (im, cplflx, lssav, islmsk, dtf, ep1d, gflx, tgrs_1, qgrs_1, ugrs_1, vgrs_1,          &
        adjsfcdlw, adjsfcdsw, adjnirbmd, adjnirdfd, adjvisbmd, adjvisdfd, adjsfculw, adjnirbmu, adjnirdfu, adjvisbmu, adjvisdfu,    &
        t2m, q2m, u10m, v10m, tsfc, pgr, xcosz, evbs, evcw, trans, sbsno, snowc, snohf, con_eps, con_epsm1,                         &
        epi, gfluxi, t1, q1, u1, v1, dlwsfci_cpl, dswsfci_cpl, dlwsfc_cpl, dswsfc_cpl, dnirbmi_cpl, dnirdfi_cpl, dvisbmi_cpl,       &
        dvisdfi_cpl, dnirbm_cpl, dnirdf_cpl, dvisbm_cpl, dvisdf_cpl, nlwsfci_cpl, nlwsfc_cpl, t2mi_cpl, q2mi_cpl, u10mi_cpl,        &
        v10mi_cpl, tsfci_cpl, psurfi_cpl, nnirbmi_cpl, nnirdfi_cpl, nvisbmi_cpl, nvisdfi_cpl, nswsfci_cpl, nswsfc_cpl, nnirbm_cpl,  &
        nnirdf_cpl, nvisbm_cpl, nvisdf_cpl, gflux, evbsa, evcwa, transa, sbsnoa, snowca, snohfa, ep, tmpmin, tmpmax, spfhmin,       &
        spfhmax, wind10mmax, u10mmax, v10mmax, dpt2m, errmsg, errflg)

        use machine,               only: kind_phys

        implicit none

        integer,                              intent(in) :: im
        logical,                              intent(in) :: cplflx, lssav
        integer, dimension(im),               intent(in) :: islmsk

        real(kind=kind_phys),                 intent(in) :: dtf, con_eps, con_epsm1

        real(kind=kind_phys), dimension(im),  intent(in)  :: ep1d, gflx, tgrs_1, qgrs_1, ugrs_1, vgrs_1, adjsfcdlw, adjsfcdsw, &
          adjnirbmd, adjnirdfd, adjvisbmd, adjvisdfd, adjsfculw, adjnirbmu, adjnirdfu, adjvisbmu, adjvisdfu,                   &
          t2m, q2m, u10m, v10m, tsfc, pgr, xcosz, evbs, evcw, trans, sbsno, snowc, snohf

        real(kind=kind_phys), dimension(im),  intent(inout) :: epi, gfluxi, t1, q1, u1, v1, dlwsfci_cpl, dswsfci_cpl, dlwsfc_cpl, &
          dswsfc_cpl, dnirbmi_cpl, dnirdfi_cpl, dvisbmi_cpl, dvisdfi_cpl, dnirbm_cpl, dnirdf_cpl, dvisbm_cpl, dvisdf_cpl, &
          nlwsfci_cpl, nlwsfc_cpl, t2mi_cpl, q2mi_cpl, u10mi_cpl, v10mi_cpl, tsfci_cpl, psurfi_cpl, nnirbmi_cpl, nnirdfi_cpl, &
          nvisbmi_cpl, nvisdfi_cpl, nswsfci_cpl, nswsfc_cpl, nnirbm_cpl, nnirdf_cpl, nvisbm_cpl, nvisdf_cpl, gflux, evbsa, &
          evcwa, transa, sbsnoa, snowca, snohfa, ep, tmpmin, tmpmax, spfhmin, spfhmax, wind10mmax, u10mmax, v10mmax, dpt2m

        character(len=*),                     intent(out) :: errmsg
        integer,                              intent(out) :: errflg

        real(kind=kind_phys), parameter :: albdf   = 0.06

        integer :: i
        real(kind=kind_phys) :: tem, xcosz_loc, ocalnirdf_cpl, ocalnirbm_cpl, ocalvisdf_cpl, ocalvisbm_cpl

        ! Initialize CCPP error handling variables
        errmsg = ''
        errflg = 0

        do i=1,im
          epi(i)     = ep1d(i)
          gfluxi(i)  = gflx(i)
          t1(i)      = tgrs_1(i)
          q1(i)      = qgrs_1(i)
          u1(i)      = ugrs_1(i)
          v1(i)      = vgrs_1(i)
        enddo

        if (cplflx) then
          do i=1,im
            dlwsfci_cpl (i) = adjsfcdlw(i)
            dswsfci_cpl (i) = adjsfcdsw(i)
            dlwsfc_cpl  (i) = dlwsfc_cpl(i) + adjsfcdlw(i)*dtf
            dswsfc_cpl  (i) = dswsfc_cpl(i) + adjsfcdsw(i)*dtf
            dnirbmi_cpl (i) = adjnirbmd(i)
            dnirdfi_cpl (i) = adjnirdfd(i)
            dvisbmi_cpl (i) = adjvisbmd(i)
            dvisdfi_cpl (i) = adjvisdfd(i)
            dnirbm_cpl  (i) = dnirbm_cpl(i) + adjnirbmd(i)*dtf
            dnirdf_cpl  (i) = dnirdf_cpl(i) + adjnirdfd(i)*dtf
            dvisbm_cpl  (i) = dvisbm_cpl(i) + adjvisbmd(i)*dtf
            dvisdf_cpl  (i) = dvisdf_cpl(i) + adjvisdfd(i)*dtf
            nlwsfci_cpl (i) = adjsfcdlw(i) - adjsfculw(i)
            nlwsfc_cpl  (i) = nlwsfc_cpl(i) + nlwsfci_cpl(i)*dtf
            t2mi_cpl    (i) = t2m(i)
            q2mi_cpl    (i) = q2m(i)
            u10mi_cpl   (i) = u10m(i)
            v10mi_cpl   (i) = v10m(i)
            tsfci_cpl   (i) = tsfc(i)
            psurfi_cpl  (i) = pgr(i)
          enddo

  !  ---  estimate mean albedo for ocean point without ice cover and apply
  !       them to net SW heat fluxes

          do i=1,im
            if (islmsk(i) /= 1) then  ! not a land point
  !  ---  compute open water albedo
              xcosz_loc = max( 0.0, min( 1.0, xcosz(i) ))
              ocalnirdf_cpl = 0.06
              ocalnirbm_cpl = max(albdf, 0.026/(xcosz_loc**1.7+0.065)  &
       &                       + 0.15 * (xcosz_loc-0.1) * (xcosz_loc-0.5) &
       &                       * (xcosz_loc-1.0))
              ocalvisdf_cpl = 0.06
              ocalvisbm_cpl = ocalnirbm_cpl

              nnirbmi_cpl(i) = adjnirbmd(i)-adjnirbmd(i)*ocalnirbm_cpl
              nnirdfi_cpl(i) = adjnirdfd(i)-adjnirdfd(i)*ocalnirdf_cpl
              nvisbmi_cpl(i) = adjvisbmd(i)-adjvisbmd(i)*ocalvisbm_cpl
              nvisdfi_cpl(i) = adjvisdfd(i)-adjvisdfd(i)*ocalvisdf_cpl
            else
              nnirbmi_cpl(i) = adjnirbmd(i) - adjnirbmu(i)
              nnirdfi_cpl(i) = adjnirdfd(i) - adjnirdfu(i)
              nvisbmi_cpl(i) = adjvisbmd(i) - adjvisbmu(i)
              nvisdfi_cpl(i) = adjvisdfd(i) - adjvisdfu(i)
            endif
            nswsfci_cpl(i) = nnirbmi_cpl(i) + nnirdfi_cpl(i)   &
                            + nvisbmi_cpl(i) + nvisdfi_cpl(i)
            nswsfc_cpl(i)  = nswsfc_cpl(i)  + nswsfci_cpl(i)*dtf
            nnirbm_cpl(i)  = nnirbm_cpl(i)  + nnirbmi_cpl(i)*dtf
            nnirdf_cpl(i)  = nnirdf_cpl(i)  + nnirdfi_cpl(i)*dtf
            nvisbm_cpl(i)  = nvisbm_cpl(i)  + nvisbmi_cpl(i)*dtf
            nvisdf_cpl(i)  = nvisdf_cpl(i)  + nvisdfi_cpl(i)*dtf
          enddo
        endif

        if (lssav) then
          do i=1,im
            gflux(i)   = gflux(i)  + gflx(i)  * dtf
            evbsa(i)   = evbsa(i)  + evbs(i)  * dtf
            evcwa(i)   = evcwa(i)  + evcw(i)  * dtf
            transa(i)  = transa(i) + trans(i) * dtf
            sbsnoa(i)  = sbsnoa(i) + sbsno(i) * dtf
            snowca(i)  = snowca(i) + snowc(i) * dtf
            snohfa(i)  = snohfa(i) + snohf(i) * dtf
            ep(i)      = ep(i)     + ep1d(i)  * dtf

            tmpmax(i)  = max(tmpmax(i),t2m(i))
            tmpmin(i)  = min(tmpmin(i),t2m(i))

            spfhmax(i) = max(spfhmax(i),q2m(i))
            spfhmin(i) = min(spfhmin(i),q2m(i))
          enddo

          do i=1, im
  ! find max wind speed then decompose
             tem = sqrt(u10m(i)*u10m(i) + v10m(i)*v10m(i))
             if (tem > wind10mmax(i)) then
                wind10mmax(i) = tem
                u10mmax(i)    = u10m(i)
                v10mmax(i)    = v10m(i)
             endif

  ! Compute dew point, first using vapor pressure
             tem = max(pgr(i) * q2m(i) / ( con_eps - con_epsm1 *q2m(i)), 1.e-8)
             dpt2m(i) = 243.5 / ( ( 17.67 / log(tem/611.2) ) - 1.) + 273.14
          enddo

        endif

      end subroutine GFS_surface_generic_post_run

      end module GFS_surface_generic_post

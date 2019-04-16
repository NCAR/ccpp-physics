!> \file GFS_surface_generic.F90
!!  Contains code related to all GFS surface schemes.

      module GFS_surface_generic_pre

      contains

      subroutine GFS_surface_generic_pre_init ()
      end subroutine GFS_surface_generic_pre_init

      subroutine GFS_surface_generic_pre_finalize()
      end subroutine GFS_surface_generic_pre_finalize

#if 0
!> \section arg_table_GFS_surface_generic_pre_run Argument Table
!! | local_name     | standard_name                                                                | long_name                                                            | units      | rank | type      |    kind   | intent | optional |
!! |----------------|------------------------------------------------------------------------------|----------------------------------------------------------------------|------------|------|-----------|-----------|--------|----------|
!! | im             | horizontal_loop_extent                                                       | horizontal loop extent                                               | count      |    0 | integer   |           | in     | F        |
!! | levs           | vertical_dimension                                                           | number of vertical levels                                            | count      |    0 | integer   |           | in     | F        |
!! | cplflx         | flag_for_flux_coupling                                                       | flag controlling cplflx collection (default off)                     | flag       |    0 | logical   |           | in     | F        |
!! | vfrac          | vegetation_area_fraction                                                     | areal fractional cover of green vegetation                           | frac       |    1 | real      | kind_phys | in     | F        |
!! | islmsk         | sea_land_ice_mask                                                            | landmask: sea/land/ice=0/1/2                                         | flag       |    1 | integer   |           | in     | F        |
!! | lndfrac        | land_area_fraction                                                           | fraction of horizontal grid area occupied by land                    | frac       |    1 | real      | kind_phys | in     | F        |
!! | lakfrac        | lake_area_fraction                                                           | fraction of horizontal grid area occupied by lake                    | frac       |    1 | real      | kind_phys | in     | F        |
!! | ocnfrac        | sea_area_fraction                                                            | fraction of horizontal grid area occupied by ocean                   | frac       |    1 | real      | kind_phys | in     | F        |
!! | idry           | flag_nonzero_land_surface_fraction                                           | flag indicating presence of some land surface area fraction          | flag       |    1 | integer   |           | inout  | F        |
!! | iice           | flag_nonzero_sea_ice_surface_fraction                                        | flag indicating presence of some sea ice surface area fraction       | flag       |    1 | integer   |           | inout  | F        |
!! | ilak           | flag_nonzero_lake_surface_fraction                                           | flag indicating presence of some lake surface area fraction          | flag       |    1 | integer   |           | inout  | F        |
!! | iocn           | flag_nonzero_ocean_surface_fraction                                          | flag indicating presence of some ocean surface area fraction         | flag       |    1 | integer   |           | inout  | F        |
!! | iwet           | flag_nonzero_wet_surface_fraction                                            | flag indicating presence of some ocean or lake surface area fraction | flag       |    1 | integer   |           | inout  | F        |
!! | fice           | sea_ice_concentration                                                        | ice fraction over open water                                         | frac       |    1 | real      | kind_phys | in     | F        |
!! | cimin          | minimum_sea_ice_concentration                                                | minimum sea ice concentration                                        | frac       |    0 | real      | kind_phys | in     | F        |
!! | zorl           | surface_roughness_length                                                     | surface roughness length                                             | cm         |    1 | real      | kind_phys | in     | F        |
!! | zorlo          | surface_roughness_length_over_ocean                                          | surface roughness length over ocean                                  | cm         |    1 | real      | kind_phys | inout  | F        |
!! | zorll          | surface_roughness_length_over_land                                           | surface roughness length over land                                   | cm         |    1 | real      | kind_phys | inout  | F        |
!! | zorl_ocn       | surface_roughness_length_over_ocean_interstitial                             | surface roughness length over ocean (temporary use as interstitial)  | cm         |    1 | real      | kind_phys | inout  | F        |
!! | zorl_lnd       | surface_roughness_length_over_land_interstitial                              | surface roughness length over land  (temporary use as interstitial)  | cm         |    1 | real      | kind_phys | inout  | F        |
!! | zorl_ice       | surface_roughness_length_over_ice_interstitial                               | surface roughness length over ice   (temporary use as interstitial)  | cm         |    1 | real      | kind_phys | inout  | F        |
!! | snowd          | surface_snow_thickness_water_equivalent                                      | water equivalent snow depth                                          | mm         |    1 | real      | kind_phys | in     | F        |
!! | snowd_ocn      | surface_snow_thickness_water_equivalent_over_ocean                           | water equivalent snow depth over ocean                               | mm         |    1 | real      | kind_phys | inout  | F        |
!! | snowd_lnd      | surface_snow_thickness_water_equivalent_over_land                            | water equivalent snow depth over land                                | mm         |    1 | real      | kind_phys | inout  | F        |
!! | snowd_ice      | surface_snow_thickness_water_equivalent_over_ice                             | water equivalent snow depth over ice                                 | mm         |    1 | real      | kind_phys | inout  | F        |
!! | tprcp          | nonnegative_lwe_thickness_of_precipitation_amount_on_dynamics_timestep       | total precipitation amount in each time step                         | m          |    1 | real      | kind_phys | in     | F        |
!! | tprcp_ocn      | nonnegative_lwe_thickness_of_precipitation_amount_on_dynamics_timestep_over_ocean | total precipitation amount in each time step over ocean         | m          |    1 | real      | kind_phys | inout  | F        |
!! | tprcp_lnd      | nonnegative_lwe_thickness_of_precipitation_amount_on_dynamics_timestep_over_land  | total precipitation amount in each time step over land          | m          |    1 | real      | kind_phys | inout  | F        |
!! | tprcp_ice      | nonnegative_lwe_thickness_of_precipitation_amount_on_dynamics_timestep_over_ice   | total precipitation amount in each time step over ice           | m          |    1 | real      | kind_phys | inout  | F        |
!! | uustar         | surface_friction_velocity                                                    | boundary layer parameter                                             | m s-1      |    1 | real      | kind_phys | in     | F        |
!! | uustar_lnd     | surface_friction_velocity_over_land                                          | surface friction velocity over land                                  | m s-1      |    1 | real      | kind_phys | inout  | F        |
!! | uustar_ice     | surface_friction_velocity_over_ice                                           | surface friction velocity over ice                                   | m s-1      |    1 | real      | kind_phys | inout  | F        |
!! | weasd          | water_equivalent_accumulated_snow_depth                                      | water equiv of acc snow depth over land and sea ice                  | mm         |    1 | real      | kind_phys | in     | F        |
!! | weasd_lnd      | water_equivalent_accumulated_snow_depth_over_land                            | water equiv of acc snow depth over land                              | mm         |    1 | real      | kind_phys | inout  | F        |
!! | weasd_ice      | water_equivalent_accumulated_snow_depth_over_ice                             | water equiv of acc snow depth over ice                               | mm         |    1 | real      | kind_phys | inout  | F        |
!! | evap_ocn       | kinematic_surface_upward_latent_heat_flux_over_ocean                         | kinematic surface upward latent heat flux over ocean                 | kg kg-1 m s-1 |    1 | real   | kind_phys | inout  | F        |
!! | hflx_ocn       | kinematic_surface_upward_sensible_heat_flux_over_ocean                       | kinematic surface upward sensible heat flux over ocean               | K m s-1    |    1 | real      | kind_phys | inout  | F        |
!! | stress_ocn     | surface_wind_stress_over_ocean                                               | surface wind stress over ocean                                       | m2 s-2     |    1 | real      | kind_phys | inout  | F        |
!! | ep1d_ice       | surface_upward_potential_latent_heat_flux_over_ice                           | surface upward potential latent heat flux over ice                   | W m-2      |    1 | real      | kind_phys | inout  | F        |
!! | isot           | soil_type_dataset_choice                                                     | soil type dataset choice                                             | index      |    0 | integer   |           | in     | F        |
!! | ivegsrc        | vegetation_type_dataset_choice                                               | land use dataset choice                                              | index      |    0 | integer   |           | in     | F        |
!! | stype          | soil_type_classification_real                                                | soil type for lsm                                                    | index      |    1 | real      | kind_phys | in     | F        |
!! | vtype          | vegetation_type_classification_real                                          | vegetation type for lsm                                              | index      |    1 | real      | kind_phys | in     | F        |
!! | slope          | surface_slope_classification_real                                            | sfc slope type for lsm                                               | index      |    1 | real      | kind_phys | in     | F        |
!! | prsik_1        | dimensionless_exner_function_at_lowest_model_interface                       | dimensionless Exner function at lowest model interface               | none       |    1 | real      | kind_phys | in     | F        |
!! | prslk_1        | dimensionless_exner_function_at_lowest_model_layer                           | dimensionless Exner function at lowest model layer                   | none       |    1 | real      | kind_phys | in     | F        |
!! | semis          | surface_longwave_emissivity                                                  | surface lw emissivity in fraction                                    | frac       |    1 | real      | kind_phys | in     | F        |
!! | adjsfcdlw      | surface_downwelling_longwave_flux                                            | surface downwelling longwave flux at current time                    | W m-2      |    1 | real      | kind_phys | in     | F        |
!! | tsfc           | surface_skin_temperature                                                     | surface skin temperature                                             | K          |    1 | real      | kind_phys | in     | F        |
!! | tsfco          | sea_surface_temperature                                                      | sea surface temperature                                              | K          |    1 | real      | kind_phys | inout  | F        |
!! | tsfcl          | surface_skin_temperature_over_land                                           | surface skin temperature over land                                   | K          |    1 | real      | kind_phys | inout  | F        |
!! | tsfc_ocn       | surface_skin_temperature_over_ocean_interstitial                             | surface skin temperature over ocean (temporary use as interstitial)  | K          |    1 | real      | kind_phys | inout  | F        |
!! | tsfc_lnd       | surface_skin_temperature_over_land_interstitial                              | surface skin temperature over land  (temporary use as interstitial)  | K          |    1 | real      | kind_phys | inout  | F        |
!! | tsfc_ice       | surface_skin_temperature_over_ice_interstitial                               | surface skin temperature over ice   (temporary use as interstitial)  | K          |    1 | real      | kind_phys | inout  | F        |
!! | tisfc          | sea_ice_temperature                                                          | sea ice surface skin temperature                                     | K          |    1 | real      | kind_phys | inout  | F        |
!! | phil           | geopotential                                                                 | geopotential at model layer centers                                  | m2 s-2     |    2 | real      | kind_phys | in     | F        |
!! | con_g          | gravitational_acceleration                                                   | gravitational acceleration                                           | m s-2      |    0 | real      | kind_phys | in     | F        |
!! | sigmaf         | bounded_vegetation_area_fraction                                             | areal fractional cover of green vegetation bounded on the bottom     | frac       |    1 | real      | kind_phys | inout  | F        |
!! | soiltyp        | soil_type_classification                                                     | soil type at each grid cell                                          | index      |    1 | integer   |           | inout  | F        |
!! | vegtype        | vegetation_type_classification                                               | vegetation type at each grid cell                                    | index      |    1 | integer   |           | inout  | F        |
!! | slopetyp       | surface_slope_classification                                                 | surface slope type at each grid cell                                 | index      |    1 | integer   |           | inout  | F        |
!! | work3          | ratio_of_exner_function_between_midlayer_and_interface_at_lowest_model_layer | Exner function ratio bt midlayer and interface at 1st layer          | ratio      |    1 | real      | kind_phys | inout  | F        |
!! | gabsbdlw       | surface_downwelling_longwave_flux_absorbed_by_ground                         | total sky surface downward longwave flux absorbed by the ground      | W m-2      |    1 | real      | kind_phys | inout  | F        |
!! | tsurf          | surface_skin_temperature_after_iteration                                     | surface skin temperature after iteration                             | K          |    1 | real      | kind_phys | inout  | F        |
!! | tsurf_ocn      | surface_skin_temperature_after_iteration_over_ocean                          | surface skin temperature after iteration over ocean                  | K          |    1 | real      | kind_phys | inout  | F        |
!! | tsurf_lnd      | surface_skin_temperature_after_iteration_over_land                           | surface skin temperature after iteration over land                   | K          |    1 | real      | kind_phys | inout  | F        |
!! | tsurf_ice      | surface_skin_temperature_after_iteration_over_ice                            | surface skin temperature after iteration over ice                    | K          |    1 | real      | kind_phys | inout  | F        |
!! | zlvl           | height_above_ground_at_lowest_model_layer                                    | layer 1 height above ground (not MSL)                                | m          |    1 | real      | kind_phys | inout  | F        |
!! | do_sppt        | flag_for_stochastic_surface_physics_perturbations                            | flag for stochastic surface physics perturbations                    | flag       |    0 | logical   |           | in     | F        |
!! | dtdtr          | tendency_of_air_temperature_due_to_radiative_heating_on_physics_time_step    | temp. change due to radiative heating per time step                  | K          |    2 | real      | kind_phys | out    | F        |
!! | drain_cpl      | tendency_of_lwe_thickness_of_precipitation_amount_for_coupling               | change in rain_cpl (coupling_type)                                   | m          |    1 | real      | kind_phys | out    | F        |
!! | dsnow_cpl      | tendency_of_lwe_thickness_of_snow_amount_for_coupling                        | change in show_cpl (coupling_type)                                   | m          |    1 | real      | kind_phys | out    | F        |
!! | rain_cpl       | lwe_thickness_of_precipitation_amount_for_coupling                           | total rain precipitation                                             | m          |    1 | real      | kind_phys | in     | F        |
!! | snow_cpl       | lwe_thickness_of_snow_amount_for_coupling                                    | total snow precipitation                                             | m          |    1 | real      | kind_phys | in     | F        |
!! | do_sfcperts    | flag_for_stochastic_surface_perturbations                                    | flag for stochastic surface perturbations option                     | flag       |    0 | logical   |           | in     | F        |
!! | nsfcpert       | number_of_surface_perturbations                                              | number of surface perturbations                                      | count      |    0 | integer   |           | in     | F        |
!! | sfc_wts        | weights_for_stochastic_surface_physics_perturbation                          | weights for stochastic surface physics perturbation                  | none       |    2 | real      | kind_phys | in     | F        |
!! | pertz0         | magnitude_of_perturbation_of_momentum_roughness_length                       | magnitude of perturbation of momentum roughness length               | frac       |    1 | real      | kind_phys | in     | F        |
!! | pertzt         | magnitude_of_perturbation_of_heat_to_momentum_roughness_length_ratio         | magnitude of perturbation of heat to momentum roughness length r.    | frac       |    1 | real      | kind_phys | in     | F        |
!! | pertshc        | magnitude_of_perturbation_of_soil_type_b_parameter                           | magnitude of perturbation of soil type b parameter                   | frac       |    1 | real      | kind_phys | in     | F        |
!! | pertlai        | magnitude_of_perturbation_of_leaf_area_index                                 | magnitude of perturbation of leaf area index                         | frac       |    1 | real      | kind_phys | in     | F        |
!! | pertvegf       | magnitude_of_perturbation_of_vegetation_fraction                             | magnitude of perturbation of vegetation fraction                     | frac       |    1 | real      | kind_phys | in     | F        |
!! | z01d           | perturbation_of_momentum_roughness_length                                    | perturbation of momentum roughness length                            | frac       |    1 | real      | kind_phys | out    | F        |
!! | zt1d           | perturbation_of_heat_to_momentum_roughness_length_ratio                      | perturbation of heat to momentum roughness length ratio              | frac       |    1 | real      | kind_phys | out    | F        |
!! | bexp1d         | perturbation_of_soil_type_b_parameter                                        | perturbation of soil type "b" parameter                              | frac       |    1 | real      | kind_phys | out    | F        |
!! | xlai1d         | perturbation_of_leaf_area_index                                              | perturbation of leaf area index                                      | frac       |    1 | real      | kind_phys | out    | F        |
!! | vegf1d         | perturbation_of_vegetation_fraction                                          | perturbation of vegetation fraction                                  | frac       |    1 | real      | kind_phys | out    | F        |
!! | errmsg         | ccpp_error_message                                                           | error message for error handling in CCPP                             | none       |    0 | character | len=*     | out    | F        |
!! | errflg         | ccpp_error_flag                                                              | error flag for error handling in CCPP                                | flag       |    0 | integer   |           | out    | F        |
!!
#endif
      subroutine GFS_surface_generic_pre_run (im, levs, cplflx, vfrac, islmsk, lndfrac, lakfrac, ocnfrac,  &
                          idry, iice, ilak, iocn, iwet, fice, cimin, zorl, zorlo, zorll, zorl_ocn,         &
                          zorl_lnd, zorl_ice, snowd, snowd_ocn, snowd_lnd, snowd_ice, tprcp, tprcp_ocn,    &
                          tprcp_lnd, tprcp_ice, uustar, uustar_lnd, uustar_ice, weasd, weasd_lnd,          &
                          weasd_ice, evap_ocn, hflx_ocn, stress_ocn, ep1d_ice, isot, ivegsrc, stype, vtype,&
                          slope, prsik_1, prslk_1, semis, adjsfcdlw, tsfc, tsfco, tsfcl, tsfc_ocn,         &
                          tsfc_lnd, tsfc_ice, tisfc, phil, con_g, sigmaf, soiltyp, vegtype, slopetyp,      &
                          work3, gabsbdlw, tsurf, tsurf_ocn, tsurf_lnd, tsurf_ice, zlvl, do_sppt, dtdtr,   &
                          drain_cpl, dsnow_cpl, rain_cpl, snow_cpl, do_sfcperts, nsfcpert, sfc_wts,        &
                          pertz0, pertzt, pertshc, pertlai, pertvegf, z01d, zt1d, bexp1d, xlai1d, vegf1d,  &
                          errmsg, errflg)

        use machine,               only: kind_phys
        use surface_perturbation,  only: cdfnor

        implicit none

        ! Interface variables
        logical, intent(in) :: cplflx

        integer, intent(in) :: im, levs, isot, ivegsrc
        integer, dimension(im), intent(in) :: islmsk
        integer, dimension(im), intent(inout) :: idry, iice, ilak, iocn, iwet
        integer, dimension(im), intent(inout) :: soiltyp, vegtype, slopetyp

        real(kind=kind_phys), intent(in) :: con_g, cimin
        real(kind=kind_phys), dimension(im), intent(in) :: lndfrac, lakfrac, ocnfrac, fice
        real(kind=kind_phys), dimension(im), intent(in) :: zorl, snowd, tprcp, uustar, weasd
        real(kind=kind_phys), dimension(im), intent(in) :: vfrac, stype, vtype, slope, prsik_1, prslk_1, &
          semis, adjsfcdlw, tsfc
        real(kind=kind_phys), dimension(im,levs), intent(in) :: phil

        real(kind=kind_phys), dimension(im), intent(inout) :: sigmaf, work3, gabsbdlw, tsurf, zlvl
        real(kind=kind_phys), dimension(im), intent(inout) :: zorlo, zorll, tsfco, tsfcl, tisfc
        real(kind=kind_phys), dimension(im), intent(inout) :: snowd_ocn, snowd_lnd, snowd_ice, tprcp_ocn, &
          tprcp_lnd, tprcp_ice, zorl_ocn, zorl_lnd, zorl_ice, tsfc_ocn, tsfc_lnd, tsfc_ice, tsurf_ocn,    &
          tsurf_lnd, tsurf_ice, uustar_lnd, uustar_ice, weasd_lnd, weasd_ice, evap_ocn, hflx_ocn,         &
          stress_ocn, ep1d_ice

        ! Stochastic physics / surface perturbations
        logical, intent(in) :: do_sppt
        real(kind=kind_phys), dimension(im,levs),     intent(out) :: dtdtr
        real(kind=kind_phys), dimension(im),          intent(out) :: drain_cpl
        real(kind=kind_phys), dimension(im),          intent(out) :: dsnow_cpl
        real(kind=kind_phys), dimension(im),          intent(in)  :: rain_cpl
        real(kind=kind_phys), dimension(im),          intent(in)  :: snow_cpl
        logical, intent(in) :: do_sfcperts
        integer, intent(in) :: nsfcpert
        real(kind=kind_phys), dimension(im,nsfcpert), intent(in)  :: sfc_wts
        real(kind=kind_phys), dimension(:),           intent(in)  :: pertz0
        real(kind=kind_phys), dimension(:),           intent(in)  :: pertzt
        real(kind=kind_phys), dimension(:),           intent(in)  :: pertshc
        real(kind=kind_phys), dimension(:),           intent(in)  :: pertlai
        real(kind=kind_phys), dimension(:),           intent(in)  :: pertvegf
        real(kind=kind_phys), dimension(im),          intent(out) :: z01d
        real(kind=kind_phys), dimension(im),          intent(out) :: zt1d
        real(kind=kind_phys), dimension(im),          intent(out) :: bexp1d
        real(kind=kind_phys), dimension(im),          intent(out) :: xlai1d
        real(kind=kind_phys), dimension(im),          intent(out) :: vegf1d

        ! CCPP error handling
        character(len=*), intent(out) :: errmsg
        integer,          intent(out) :: errflg

        ! Local variables
        integer              :: i
        real(kind=kind_phys) :: onebg
        real(kind=kind_phys) :: cdfz

        ! Set constants
        onebg  = 1.0/con_g

        ! Initialize CCPP error handling variables
        errmsg = ''
        errflg = 0

        ! Set initial quantities for stochastic physics deltas
        if (do_sppt) then
          dtdtr     = 0.0
          do i=1,im
            drain_cpl(i) = rain_cpl (i)
            dsnow_cpl(i) = snow_cpl (i)
          enddo
        endif

        ! Scale random patterns for surface perturbations with perturbation size
        ! Turn vegetation fraction pattern into percentile pattern
        if (do_sfcperts) then
          if (pertz0(1) > 0.) then
            z01d(:) = pertz0(1) * sfc_wts(:,1)
  !          if (me == 0) print*,'sfc_wts(:,1) min and max',minval(sfc_wts(:,1)),maxval(sfc_wts(:,1))
  !          if (me == 0) print*,'z01d min and max ',minval(z01d),maxval(z01d)
          endif
          if (pertzt(1) > 0.) then
            zt1d(:) = pertzt(1) * sfc_wts(:,2)
          endif
          if (pertshc(1) > 0.) then
            bexp1d(:) = pertshc(1) * sfc_wts(:,3)
          endif
          if (pertlai(1) > 0.) then
            xlai1d(:) = pertlai(1) * sfc_wts(:,4)
          endif
  ! --- do the albedo percentile calculation in GFS_radiation_driver instead --- !
  !        if (pertalb(1) > 0.) then
  !          do i=1,im
  !            call cdfnor(sfc_wts(i,5),cdfz)
  !            alb1d(i) = cdfz
  !          enddo
  !        endif
          if (pertvegf(1) > 0.) then
            do i=1,im
              call cdfnor(sfc_wts(i,6),cdfz)
              vegf1d(i) = cdfz
            enddo
          endif
        endif

        ! End of stochastic physics / surface perturbation

        do i = 1, im
          if(lndfrac(i)<1.)    iwet(i) = 1
          if(lndfrac(i)>0.)    idry(i) = 1
          if(ocnfrac(i)>0.)    iocn(i) = 1
          if(lakfrac(i)>0.)    ilak(i) = 1
          if(iwet(i) == 1 .and. fice(i) >= cimin) iice(i) = 1
        enddo

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
          tsurf(i)   = tsfc(i)
          zlvl(i)    = phil(i,1) * onebg
        end do

        do i=1,im
          if (.not. cplflx) then
            zorll(i) = zorl(i)
            zorlo(i) = zorl(i)
            tsfcl(i) = tsfc(i)
            tsfco(i) = tsfc(i)
            tisfc(i) = tsfc(i)
          end if
          if(iwet(i) == 1) then
            snowd_ocn(i) = snowd(i)
            tprcp_ocn(i) = tprcp(i)
            zorl_ocn(i) = zorlo(i)
            tsfc_ocn(i) = tsfco(i)
            tsurf_ocn(i)= tsfco(i)
            evap_ocn(i) = 0.
            hflx_ocn(i) = 0.
            stress_ocn(i) = 0.
          endif
  !
          if (idry(i) == 1) then
            uustar_lnd(i) = uustar(i)
            weasd_lnd(i) = weasd(i)
            tprcp_lnd(i) = tprcp(i)
            zorl_lnd(i) = zorll(i)
            tsfc_lnd(i) = tsfcl(i)
            tsurf_lnd(i) = tsfcl(i)
            snowd_lnd(i) = snowd(i)
          end if
  !
          if (iice(i) == 1) then
            uustar_ice(i) = uustar(i)
            weasd_ice(i) = weasd(i)
            tprcp_ice(i) = tprcp(i)
            zorl_ice(i) = zorll(i)
            tsfc_ice(i) = tisfc(i)
            tsurf_ice(i)= tisfc(i)
            snowd_ice(i) = snowd(i)
            ep1d_ice(i) = 0.
          end if
        enddo

      end subroutine GFS_surface_generic_pre_run

      end module GFS_surface_generic_pre

      module GFS_surface_generic_post
        
      use machine,               only: kind_phys
      
      implicit none

      private

      public GFS_surface_generic_post_init, GFS_surface_generic_post_finalize, GFS_surface_generic_post_run

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
!! | cplwav         | flag_for_wave_coupling                                                                                              | flag controlling cplwav collection (default off)                                    | flag        |    0 | logical    |           | in     | F        |
!! | lssav          | flag_diagnostics                                                                                                    | logical flag for storing diagnostics                                                | flag        |    0 | logical    |           | in     | F        |
!! | islmsk         | sea_land_ice_mask                                                                                                   | landmask: sea/land/ice=0/1/2                                                        | flag        |    1 | integer    |           | in     | F        |
!! | idry           | flag_nonzero_land_surface_fraction                                                                                  | flag indicating presence of some land surface area fraction                         | flag        |    1 | integer    |           | in     | F        |
!! | iwet           | flag_nonzero_wet_surface_fraction                                                                                   | flag indicating presence of some ocean or lake surface area fraction                | flag        |    1 | integer    |           | in     | F        |
!! | iice           | flag_nonzero_sea_ice_surface_fraction                                                                               | flag indicating presence of some sea ice surface area fraction                      | flag        |    1 | integer    |           | in     | F        |
!! | ilak           | flag_nonzero_lake_surface_fraction                                                                                  | flag indicating presence of some lake surface area fraction                         | flag        |    1 | integer    |           | in     | F        |
!! | dtf            | time_step_for_dynamics                                                                                              | dynamics timestep                                                                   | s           |    0 | real       | kind_phys | in     | F        |
!! | tgrs_1         | air_temperature_at_lowest_model_layer                                                                               | mean temperature at lowest model layer                                              | K           |    1 | real       | kind_phys | in     | F        |
!! | qgrs_1         | water_vapor_specific_humidity_at_lowest_model_layer                                                                 | specific humidity at lowest model layer                                             | kg kg-1     |    1 | real       | kind_phys | in     | F        |
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
!! | pgr            | surface_air_pressure                                                                                                | surface pressure                                                                    | Pa          |    1 | real       | kind_phys | in     | F        |
!! | xcosz          | instantaneous_cosine_of_zenith_angle                                                                                | cosine of zenith angle at current time                                              | none        |    1 | real       | kind_phys | in     | F        |
!! | evbs           | soil_upward_latent_heat_flux                                                                                        | soil upward latent heat flux                                                        | W m-2       |    1 | real       | kind_phys | in     | F        |
!! | evcw           | canopy_upward_latent_heat_flux                                                                                      | canopy upward latent heat flux                                                      | W m-2       |    1 | real       | kind_phys | in     | F        |
!! | trans          | transpiration_flux                                                                                                  | total plant transpiration rate                                                      | kg m-2 s-1  |    1 | real       | kind_phys | in     | F        |
!! | sbsno          | snow_deposition_sublimation_upward_latent_heat_flux                                                                 | latent heat flux from snow depo/subl                                                | W m-2       |    1 | real       | kind_phys | in     | F        |
!! | snowc          | surface_snow_area_fraction                                                                                          | surface snow area fraction                                                          | frac        |    1 | real       | kind_phys | in     | F        |
!! | snohf          | snow_freezing_rain_upward_latent_heat_flux                                                                          | latent heat flux due to snow and frz rain                                           | W m-2       |    1 | real       | kind_phys | in     | F        |
!! | epi            | instantaneous_surface_potential_evaporation                                                                         | instantaneous sfc potential evaporation                                             | W m-2       |    1 | real       | kind_phys | inout  | F        |
!! | gfluxi         | instantaneous_surface_ground_heat_flux                                                                              | instantaneous sfc ground heat flux                                                  | W m-2       |    1 | real       | kind_phys | inout  | F        |
!! | t1             | air_temperature_at_lowest_model_layer_for_diag                                                                      | layer 1 temperature for diag                                                        | K           |    1 | real       | kind_phys | inout  | F        |
!! | q1             | water_vapor_specific_humidity_at_lowest_model_layer_for_diag                                                        | layer 1 specific humidity for diag                                                  | kg kg-1     |    1 | real       | kind_phys | inout  | F        |
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
!! | runoff         | total_runoff                                                                                                        | total water runoff                                                                  | kg m-2      |    1 | real       | kind_phys | inout  | F        |
!! | srunoff        | surface_runoff                                                                                                      | surface water runoff (from lsm)                                                     | kg m-2      |    1 | real       | kind_phys | inout  | F        |
!! | runof          | surface_runoff_flux                                                                                                 | surface runoff flux                                                                 | g m-2 s-1   |    1 | real       | kind_phys | in     | F        |
!! | drain          | subsurface_runoff_flux                                                                                              | subsurface runoff flux                                                              | g m-2 s-1   |    1 | real       | kind_phys | in     | F        |
!! | lndfrac        | land_area_fraction                                                                                                  | fraction of horizontal grid area occupied by land                                   | frac        |    1 | real       | kind_phys | in     | F        |
!! | lakfrac        | lake_area_fraction                                                                                                  | fraction of horizontal grid area occupied by lake                                   | frac        |    1 | real       | kind_phys | in     | F        |
!! | ocnfrac        | sea_area_fraction                                                                                                   | fraction of horizontal grid area occupied by ocean                                  | frac        |    1 | real       | kind_phys | in     | F        |
!! | cice           | sea_ice_concentration                                                                                               | ice fraction over open water                                                        | frac        |    1 | real       | kind_phys | in     | F        |
!! | zorl           | surface_roughness_length                                                                                            | surface roughness length                                                            | cm          |    1 | real       | kind_phys | inout  | F        |
!! | zorlo          | surface_roughness_length_over_ocean                                                                                 | surface roughness length over ocean                                                 | cm          |    1 | real       | kind_phys | inout  | F        |
!! | zorll          | surface_roughness_length_over_land                                                                                  | surface roughness length over land                                                  | cm          |    1 | real       | kind_phys | inout  | F        |
!! | zorl_ocn       | surface_roughness_length_over_ocean_interstitial                                                                    | surface roughness length over ocean (temporary use as interstitial)                 | cm          |    1 | real       | kind_phys | in     | F        |
!! | zorl_lnd       | surface_roughness_length_over_land_interstitial                                                                     | surface roughness length over land  (temporary use as interstitial)                 | cm          |    1 | real       | kind_phys | in     | F        |
!! | zorl_ice       | surface_roughness_length_over_ice_interstitial                                                                      | surface roughness length over ice   (temporary use as interstitial)                 | cm          |    1 | real       | kind_phys | in     | F        |
!! | cd             | surface_drag_coefficient_for_momentum_in_air                                                                        | surface exchange coeff for momentum                                                 | none        |    1 | real       | kind_phys | inout  | F        |
!! | cd_ocn         | surface_drag_coefficient_for_momentum_in_air_over_ocean                                                             | surface exchange coeff for momentum over ocean                                      | none        |    1 | real       | kind_phys | in     | F        |
!! | cd_lnd         | surface_drag_coefficient_for_momentum_in_air_over_land                                                              | surface exchange coeff for momentum over land                                       | none        |    1 | real       | kind_phys | in     | F        |
!! | cd_ice         | surface_drag_coefficient_for_momentum_in_air_over_ice                                                               | surface exchange coeff for momentum over ice                                        | none        |    1 | real       | kind_phys | in     | F        |
!! | cdq            | surface_drag_coefficient_for_heat_and_moisture_in_air                                                               | surface exchange coeff heat & moisture                                              | none        |    1 | real       | kind_phys | inout  | F        |
!! | cdq_ocn        | surface_drag_coefficient_for_heat_and_moisture_in_air_over_ocean                                                    | surface exchange coeff heat & moisture over ocean                                   | none        |    1 | real       | kind_phys | in     | F        |
!! | cdq_lnd        | surface_drag_coefficient_for_heat_and_moisture_in_air_over_land                                                     | surface exchange coeff heat & moisture over land                                    | none        |    1 | real       | kind_phys | in     | F        |
!! | cdq_ice        | surface_drag_coefficient_for_heat_and_moisture_in_air_over_ice                                                      | surface exchange coeff heat & moisture over ice                                     | none        |    1 | real       | kind_phys | in     | F        |
!! | rb             | bulk_richardson_number_at_lowest_model_level                                                                        | bulk Richardson number at the surface                                               | none        |    1 | real       | kind_phys | inout  | F        |
!! | rb_ocn         | bulk_richardson_number_at_lowest_model_level_over_ocean                                                             | bulk Richardson number at the surface over ocean                                    | none        |    1 | real       | kind_phys | in     | F        |
!! | rb_lnd         | bulk_richardson_number_at_lowest_model_level_over_land                                                              | bulk Richardson number at the surface over land                                     | none        |    1 | real       | kind_phys | in     | F        |
!! | rb_ice         | bulk_richardson_number_at_lowest_model_level_over_ice                                                               | bulk Richardson number at the surface over ice                                      | none        |    1 | real       | kind_phys | in     | F        |
!! | stress         | surface_wind_stress                                                                                                 | surface wind stress                                                                 | m2 s-2      |    1 | real       | kind_phys | inout  | F        |
!! | stress_ocn     | surface_wind_stress_over_ocean                                                                                      | surface wind stress over ocean                                                      | m2 s-2      |    1 | real       | kind_phys | in     | F        |
!! | stress_lnd     | surface_wind_stress_over_land                                                                                       | surface wind stress over land                                                       | m2 s-2      |    1 | real       | kind_phys | in     | F        |
!! | stress_ice     | surface_wind_stress_over_ice                                                                                        | surface wind stress over ice                                                        | m2 s-2      |    1 | real       | kind_phys | in     | F        |
!! | ffmm           | Monin-Obukhov_similarity_function_for_momentum                                                                      | Monin-Obukhov similarity function for momentum                                      | none        |    1 | real       | kind_phys | inout  | F        |
!! | ffmm_ocn       | Monin-Obukhov_similarity_function_for_momentum_over_ocean                                                           | Monin-Obukhov similarity function for momentum over ocean                           | none        |    1 | real       | kind_phys | in     | F        |
!! | ffmm_lnd       | Monin-Obukhov_similarity_function_for_momentum_over_land                                                            | Monin-Obukhov similarity function for momentum over land                            | none        |    1 | real       | kind_phys | in     | F        |
!! | ffmm_ice       | Monin-Obukhov_similarity_function_for_momentum_over_ice                                                             | Monin-Obukhov similarity function for momentum over ice                             | none        |    1 | real       | kind_phys | in     | F        |
!! | ffhh           | Monin-Obukhov_similarity_function_for_heat                                                                          | Monin-Obukhov similarity function for heat                                          | none        |    1 | real       | kind_phys | inout  | F        |
!! | ffhh_ocn       | Monin-Obukhov_similarity_function_for_heat_over_ocean                                                               | Monin-Obukhov similarity function for heat over ocean                               | none        |    1 | real       | kind_phys | in     | F        |
!! | ffhh_lnd       | Monin-Obukhov_similarity_function_for_heat_over_land                                                                | Monin-Obukhov similarity function for heat over land                                | none        |    1 | real       | kind_phys | in     | F        |
!! | ffhh_ice       | Monin-Obukhov_similarity_function_for_heat_over_ice                                                                 | Monin-Obukhov similarity function for heat over ice                                 | none        |    1 | real       | kind_phys | in     | F        |
!! | uustar         | surface_friction_velocity                                                                                           | boundary layer parameter                                                            | m s-1       |    1 | real       | kind_phys | inout  | F        |
!! | uustar_ocn     | surface_friction_velocity_over_ocean                                                                                | surface friction velocity over ocean                                                | m s-1       |    1 | real       | kind_phys | in     | F        |
!! | uustar_lnd     | surface_friction_velocity_over_land                                                                                 | surface friction velocity over land                                                 | m s-1       |    1 | real       | kind_phys | in     | F        |
!! | uustar_ice     | surface_friction_velocity_over_ice                                                                                  | surface friction velocity over ice                                                  | m s-1       |    1 | real       | kind_phys | in     | F        |
!! | fm10           | Monin-Obukhov_similarity_function_for_momentum_at_10m                                                               | Monin-Obukhov similarity parameter for momentum at 10m                              | none        |    1 | real       | kind_phys | inout  | F        |
!! | fm10_ocn       | Monin-Obukhov_similarity_function_for_momentum_at_10m_over_ocean                                                    | Monin-Obukhov similarity parameter for momentum at 10m over ocean                   | none        |    1 | real       | kind_phys | in     | F        |
!! | fm10_lnd       | Monin-Obukhov_similarity_function_for_momentum_at_10m_over_land                                                     | Monin-Obukhov similarity parameter for momentum at 10m over land                    | none        |    1 | real       | kind_phys | in     | F        |
!! | fm10_ice       | Monin-Obukhov_similarity_function_for_momentum_at_10m_over_ice                                                      | Monin-Obukhov similarity parameter for momentum at 10m over ice                     | none        |    1 | real       | kind_phys | in     | F        |
!! | fh2            | Monin-Obukhov_similarity_function_for_heat_at_2m                                                                    | Monin-Obukhov similarity parameter for heat at 2m                                   | none        |    1 | real       | kind_phys | inout  | F        |
!! | fh2_ocn        | Monin-Obukhov_similarity_function_for_heat_at_2m_over_ocean                                                         | Monin-Obukhov similarity parameter for heat at 2m over ocean                        | none        |    1 | real       | kind_phys | in     | F        |
!! | fh2_lnd        | Monin-Obukhov_similarity_function_for_heat_at_2m_over_land                                                          | Monin-Obukhov similarity parameter for heat at 2m over land                         | none        |    1 | real       | kind_phys | in     | F        |
!! | fh2_ice        | Monin-Obukhov_similarity_function_for_heat_at_2m_over_ice                                                           | Monin-Obukhov similarity parameter for heat at 2m over ice                          | none        |    1 | real       | kind_phys | in     | F        |
!! | tsurf          | surface_skin_temperature_after_iteration                                                                            | surface skin temperature after iteration                                            | K           |    1 | real       | kind_phys | inout  | F        |
!! | tsurf_ocn      | surface_skin_temperature_after_iteration_over_ocean                                                                 | surface skin temperature after iteration over ocean                                 | K           |    1 | real       | kind_phys | in     | F        |
!! | tsurf_lnd      | surface_skin_temperature_after_iteration_over_land                                                                  | surface skin temperature after iteration over land                                  | K           |    1 | real       | kind_phys | in     | F        |
!! | tsurf_ice      | surface_skin_temperature_after_iteration_over_ice                                                                   | surface skin temperature after iteration over ice                                   | K           |    1 | real       | kind_phys | in     | F        |
!! | cmm            | surface_drag_wind_speed_for_momentum_in_air                                                                         | momentum exchange coefficient                                                       | m s-1       |    1 | real       | kind_phys | inout  | F        |
!! | cmm_ocn        | surface_drag_wind_speed_for_momentum_in_air_over_ocean                                                              | momentum exchange coefficient over ocean                                            | m s-1       |    1 | real       | kind_phys | in     | F        |
!! | cmm_lnd        | surface_drag_wind_speed_for_momentum_in_air_over_land                                                               | momentum exchange coefficient over land                                             | m s-1       |    1 | real       | kind_phys | in     | F        |
!! | cmm_ice        | surface_drag_wind_speed_for_momentum_in_air_over_ice                                                                | momentum exchange coefficient over ice                                              | m s-1       |    1 | real       | kind_phys | in     | F        |
!! | chh            | surface_drag_mass_flux_for_heat_and_moisture_in_air                                                                 | thermal exchange coefficient                                                        | kg m-2 s-1  |    1 | real       | kind_phys | inout  | F        |
!! | chh_ocn        | surface_drag_mass_flux_for_heat_and_moisture_in_air_over_ocean                                                      | thermal exchange coefficient over ocean                                             | kg m-2 s-1  |    1 | real       | kind_phys | in     | F        |
!! | chh_lnd        | surface_drag_mass_flux_for_heat_and_moisture_in_air_over_land                                                       | thermal exchange coefficient over land                                              | kg m-2 s-1  |    1 | real       | kind_phys | in     | F        |
!! | chh_ice        | surface_drag_mass_flux_for_heat_and_moisture_in_air_over_ice                                                        | thermal exchange coefficient over ice                                               | kg m-2 s-1  |    1 | real       | kind_phys | in     | F        |
!! | gflx           | upward_heat_flux_in_soil                                                                                            | soil heat flux                                                                      | W m-2       |    1 | real       | kind_phys | inout  | F        |
!! | gflx_ocn       | upward_heat_flux_in_soil_over_ocean                                                                                 | soil heat flux over ocean                                                           | W m-2       |    1 | real       | kind_phys | in     | F        |
!! | gflx_lnd       | upward_heat_flux_in_soil_over_land                                                                                  | soil heat flux over land                                                            | W m-2       |    1 | real       | kind_phys | in     | F        |
!! | gflx_ice       | upward_heat_flux_in_soil_over_ice                                                                                   | soil heat flux over ice                                                             | W m-2       |    1 | real       | kind_phys | in     | F        |
!! | ep1d           | surface_upward_potential_latent_heat_flux                                                                           | surface upward potential latent heat flux                                           | W m-2       |    1 | real       | kind_phys | inout  | F        |
!! | ep1d_ocn       | surface_upward_potential_latent_heat_flux_over_ocean                                                                | surface upward potential latent heat flux over ocean                                | W m-2       |    1 | real       | kind_phys | in     | F        |
!! | ep1d_lnd       | surface_upward_potential_latent_heat_flux_over_land                                                                 | surface upward potential latent heat flux over land                                 | W m-2       |    1 | real       | kind_phys | in     | F        |
!! | ep1d_ice       | surface_upward_potential_latent_heat_flux_over_ice                                                                  | surface upward potential latent heat flux over ice                                  | W m-2       |    1 | real       | kind_phys | in     | F        |
!! | weasd          | water_equivalent_accumulated_snow_depth                                                                             | water equiv of acc snow depth over land and sea ice                                 | mm          |    1 | real       | kind_phys | inout  | F        |
!! | weasd_lnd      | water_equivalent_accumulated_snow_depth_over_land                                                                   | water equiv of acc snow depth over land                                             | mm          |    1 | real       | kind_phys | in     | F        |
!! | weasd_ice      | water_equivalent_accumulated_snow_depth_over_ice                                                                    | water equiv of acc snow depth over ice                                              | mm          |    1 | real       | kind_phys | in     | F        |
!! | snowd          | surface_snow_thickness_water_equivalent                                                                             | water equivalent snow depth                                                         | mm          |    1 | real       | kind_phys | inout  | F        |
!! | snowd_ocn      | surface_snow_thickness_water_equivalent_over_ocean                                                                  | water equivalent snow depth over ocean                                              | mm          |    1 | real       | kind_phys | in     | F        |
!! | snowd_lnd      | surface_snow_thickness_water_equivalent_over_land                                                                   | water equivalent snow depth over land                                               | mm          |    1 | real       | kind_phys | in     | F        |
!! | snowd_ice      | surface_snow_thickness_water_equivalent_over_ice                                                                    | water equivalent snow depth over ice                                                | mm          |    1 | real       | kind_phys | in     | F        |
!! | tprcp          | nonnegative_lwe_thickness_of_precipitation_amount_on_dynamics_timestep                                              | total precipitation amount in each time step                                        | m           |    1 | real       | kind_phys | inout  | F        |
!! | tprcp_ocn      | nonnegative_lwe_thickness_of_precipitation_amount_on_dynamics_timestep_over_ocean                                   | total precipitation amount in each time step over ocean                             | m           |    1 | real       | kind_phys | in     | F        |
!! | tprcp_lnd      | nonnegative_lwe_thickness_of_precipitation_amount_on_dynamics_timestep_over_land                                    | total precipitation amount in each time step over land                              | m           |    1 | real       | kind_phys | in     | F        |
!! | tprcp_ice      | nonnegative_lwe_thickness_of_precipitation_amount_on_dynamics_timestep_over_ice                                     | total precipitation amount in each time step over ice                               | m           |    1 | real       | kind_phys | in     | F        |
!! | evap           | kinematic_surface_upward_latent_heat_flux                                                                           | kinematic surface upward latent heat flux                                           | kg kg-1 m s-1 |  1 | real       | kind_phys | inout  | F        |
!! | evap_ocn       | kinematic_surface_upward_latent_heat_flux_over_ocean                                                                | kinematic surface upward latent heat flux over ocean                                | kg kg-1 m s-1 |  1 | real       | kind_phys | in     | F        |
!! | evap_lnd       | kinematic_surface_upward_latent_heat_flux_over_land                                                                 | kinematic surface upward latent heat flux over land                                 | kg kg-1 m s-1 |  1 | real       | kind_phys | in     | F        |
!! | evap_ice       | kinematic_surface_upward_latent_heat_flux_over_ice                                                                  | kinematic surface upward latent heat flux over ice                                  | kg kg-1 m s-1 |  1 | real       | kind_phys | in     | F        |
!! | hflx           | kinematic_surface_upward_sensible_heat_flux                                                                         | kinematic surface upward sensible heat flux                                         | K m s-1     |    1 | real       | kind_phys | inout  | F        |
!! | hflx_ocn       | kinematic_surface_upward_sensible_heat_flux_over_ocean                                                              | kinematic surface upward sensible heat flux over ocean                              | K m s-1     |    1 | real       | kind_phys | in     | F        |
!! | hflx_lnd       | kinematic_surface_upward_sensible_heat_flux_over_land                                                               | kinematic surface upward sensible heat flux over land                               | K m s-1     |    1 | real       | kind_phys | in     | F        |
!! | hflx_ice       | kinematic_surface_upward_sensible_heat_flux_over_ice                                                                | kinematic surface upward sensible heat flux over ice                                | K m s-1     |    1 | real       | kind_phys | in     | F        |
!! | qss            | surface_specific_humidity                                                                                           | surface air saturation specific humidity                                            | kg kg-1     |    1 | real       | kind_phys | inout  | F        |
!! | qss_ocn        | surface_specific_humidity_over_ocean                                                                                | surface air saturation specific humidity over ocean                                 | kg kg-1     |    1 | real       | kind_phys | in     | F        |
!! | qss_lnd        | surface_specific_humidity_over_land                                                                                 | surface air saturation specific humidity over land                                  | kg kg-1     |    1 | real       | kind_phys | in     | F        |
!! | qss_ice        | surface_specific_humidity_over_ice                                                                                  | surface air saturation specific humidity over ice                                   | kg kg-1     |    1 | real       | kind_phys | in     | F        |
!! | tsfc           | surface_skin_temperature                                                                                            | surface skin temperature                                                            | K           |    1 | real       | kind_phys | inout  | F        |
!! | tsfco          | sea_surface_temperature                                                                                             | sea surface temperature                                                             | K           |    1 | real       | kind_phys | inout  | F        |
!! | tsfcl          | surface_skin_temperature_over_land                                                                                  | surface skin temperature over land                                                  | K           |    1 | real       | kind_phys | inout  | F        |
!! | tsfc_ocn       | surface_skin_temperature_over_ocean_interstitial                                                                    | surface skin temperature over ocean (temporary use as interstitial)                 | K           |    1 | real       | kind_phys | in     | F        |
!! | tsfc_lnd       | surface_skin_temperature_over_land_interstitial                                                                     | surface skin temperature over land  (temporary use as interstitial)                 | K           |    1 | real       | kind_phys | in     | F        |
!! | tsfc_ice       | surface_skin_temperature_over_ice_interstitial                                                                      | surface skin temperature over ice   (temporary use as interstitial)                 | K           |    1 | real       | kind_phys | in     | F        |
!! | tisfc          | sea_ice_temperature                                                                                                 | sea ice surface skin temperature                                                    | K           |    1 | real       | kind_phys | inout  | F        |
!! | errmsg         | ccpp_error_message                                                                                                  | error message for error handling in CCPP                                            | none        |    0 | character  | len=*     | out    | F        |
!! | errflg         | ccpp_error_flag                                                                                                     | error flag for error handling in CCPP                                               | flag        |    0 | integer    |           | out    | F        |
!!
#endif
      subroutine GFS_surface_generic_post_run (im, cplflx, cplwav, lssav, idry, iwet, iice, ilak, dtf, tgrs_1, qgrs_1, ugrs_1, vgrs_1, &
        adjsfcdlw, adjsfcdsw, adjnirbmd, adjnirdfd, adjvisbmd, adjvisdfd, adjsfculw, adjnirbmu, adjnirdfu, adjvisbmu, adjvisdfu,    &
        t2m, q2m, u10m, v10m, pgr, xcosz, evbs, evcw, trans, sbsno, snowc, snohf,                                             &
        epi, gfluxi, t1, q1, u1, v1, dlwsfci_cpl, dswsfci_cpl, dlwsfc_cpl, dswsfc_cpl, dnirbmi_cpl, dnirdfi_cpl, dvisbmi_cpl,       &
        dvisdfi_cpl, dnirbm_cpl, dnirdf_cpl, dvisbm_cpl, dvisdf_cpl, nlwsfci_cpl, nlwsfc_cpl, t2mi_cpl, q2mi_cpl, u10mi_cpl,        &
        v10mi_cpl, tsfci_cpl, psurfi_cpl, nnirbmi_cpl, nnirdfi_cpl, nvisbmi_cpl, nvisdfi_cpl, nswsfci_cpl, nswsfc_cpl, nnirbm_cpl,  &
        nnirdf_cpl, nvisbm_cpl, nvisdf_cpl, gflux, evbsa, evcwa, transa, sbsnoa, snowca, snohfa, ep,                                &
        runoff, srunoff, runof, drain, lndfrac, lakfrac, ocnfrac, cice, zorl, zorlo, zorll, zorl_ocn, zorl_lnd, zorl_ice,           &
        cd, cd_ocn, cd_lnd, cd_ice, cdq, cdq_ocn, cdq_lnd, cdq_ice, rb, rb_ocn, rb_lnd, rb_ice, stress, stress_ocn, stress_lnd,     &
        stress_ice, ffmm, ffmm_ocn, ffmm_lnd, ffmm_ice, ffhh, ffhh_ocn, ffhh_lnd, ffhh_ice, uustar, uustar_ocn, uustar_lnd,         &
        uustar_ice, fm10, fm10_ocn, fm10_lnd, fm10_ice, fh2, fh2_ocn, fh2_lnd, fh2_ice, tsurf, tsurf_ocn, tsurf_lnd, tsurf_ice,     &
        cmm, cmm_ocn, cmm_lnd, cmm_ice, chh, chh_ocn, chh_lnd, chh_ice, gflx, gflx_ocn, gflx_lnd, gflx_ice, ep1d, ep1d_ocn,         &
        ep1d_lnd, ep1d_ice, weasd, weasd_lnd, weasd_ice, snowd, snowd_ocn, snowd_lnd, snowd_ice, tprcp, tprcp_ocn, tprcp_lnd,       &
        tprcp_ice, evap, evap_ocn, evap_lnd, evap_ice, hflx, hflx_ocn, hflx_lnd, hflx_ice, qss, qss_ocn, qss_lnd, qss_ice,          &
        tsfc, tsfco, tsfcl, tsfc_ocn, tsfc_lnd, tsfc_ice, tisfc, errmsg, errflg)

        use machine,               only: kind_phys

        implicit none

        integer,                              intent(in) :: im
        logical,                              intent(in) :: cplflx, cplwav, lssav
        integer, dimension(im),               intent(in) :: idry, iwet, iice, ilak

        real(kind=kind_phys),                 intent(in) :: dtf

        real(kind=kind_phys), dimension(im),  intent(in)  :: tgrs_1, qgrs_1, ugrs_1, vgrs_1, adjsfcdlw, adjsfcdsw,    &
          adjnirbmd, adjnirdfd, adjvisbmd, adjvisdfd, adjsfculw, adjnirbmu, adjnirdfu, adjvisbmu, adjvisdfu,                      &
          t2m, q2m, u10m, v10m, pgr, xcosz, evbs, evcw, trans, sbsno, snowc, snohf, lndfrac, lakfrac, ocnfrac, cice,        &
          zorl_ocn, zorl_lnd, zorl_ice, cd_ocn, cd_lnd, cd_ice, cdq_ocn, cdq_lnd, cdq_ice, rb_ocn, rb_lnd, rb_ice, stress_ocn,    &
          stress_lnd, stress_ice, ffmm_ocn, ffmm_lnd, ffmm_ice, ffhh_ocn, ffhh_lnd, ffhh_ice, uustar_ocn, uustar_lnd, uustar_ice, &
          fm10_ocn, fm10_lnd, fm10_ice, fh2_ocn, fh2_lnd, fh2_ice, tsurf_ocn, tsurf_lnd, tsurf_ice, cmm_ocn, cmm_lnd, cmm_ice,    &
          chh_ocn, chh_lnd, chh_ice, gflx_ocn, gflx_lnd, gflx_ice, ep1d_ocn, ep1d_lnd, ep1d_ice, weasd_lnd, weasd_ice, snowd_ocn, &
          snowd_lnd, snowd_ice,tprcp_ocn, tprcp_lnd, tprcp_ice, evap_ocn, evap_lnd, evap_ice, hflx_ocn, hflx_lnd, hflx_ice,       &
          qss_ocn, qss_lnd, qss_ice, tsfc_ocn, tsfc_lnd, tsfc_ice


        real(kind=kind_phys), dimension(im),  intent(inout) :: epi, gfluxi, t1, q1, u1, v1, dlwsfci_cpl, dswsfci_cpl, dlwsfc_cpl, &
          dswsfc_cpl, dnirbmi_cpl, dnirdfi_cpl, dvisbmi_cpl, dvisdfi_cpl, dnirbm_cpl, dnirdf_cpl, dvisbm_cpl, dvisdf_cpl, &
          nlwsfci_cpl, nlwsfc_cpl, t2mi_cpl, q2mi_cpl, u10mi_cpl, v10mi_cpl, tsfci_cpl, psurfi_cpl, nnirbmi_cpl, nnirdfi_cpl, &
          nvisbmi_cpl, nvisdfi_cpl, nswsfci_cpl, nswsfc_cpl, nnirbm_cpl, nnirdf_cpl, nvisbm_cpl, nvisdf_cpl, gflux, evbsa, &
          evcwa, transa, sbsnoa, snowca, snohfa, ep, zorl, zorlo, zorll, cd, cdq, rb, stress, ffmm, ffhh, uustar, fm10, fh2,      &
          tsurf, cmm, chh, gflx, ep1d, weasd, snowd, tprcp, evap, hflx, qss, tsfc, tsfco, tsfcl, tisfc

        real(kind=kind_phys), dimension(im), intent(inout) :: runoff, srunoff
        real(kind=kind_phys), dimension(im), intent(in)    :: drain, runof

        character(len=*), intent(out) :: errmsg
        integer,          intent(out) :: errflg

        real(kind=kind_phys), parameter :: albdf   = 0.06

        integer :: i
        real(kind=kind_phys) :: tem, xcosz_loc, ocalnirdf_cpl, ocalnirbm_cpl, ocalvisdf_cpl, ocalvisbm_cpl

        ! Initialize CCPP error handling variables
        errmsg = ''
        errflg = 0

        ! --- generate ocean/land/ice composites

        do i=1, im
  !
  ! Three-way composites (fields from sfc_diff_ocean, sfc_diff_land, sfc_diff_ice)
          zorl(i)           = cmposit3(ocnfrac(i), lndfrac(i),                &
                                lakfrac(i),cice(i),                           &
                                zorl_ocn(i),  zorl_lnd(i),  zorl_ice(i))
          cd(i)             = cmposit3(ocnfrac(i),lndfrac(i),                 &
                                       lakfrac(i),cice(i),                    &
                                  cd_ocn(i),    cd_lnd(i),    cd_ice(i))
          cdq(i)            = cmposit3(ocnfrac(i),lndfrac(i),                 &
                                       lakfrac(i),cice(i),                    &
                                 cdq_ocn(i),   cdq_lnd(i),   cdq_ice(i))
          rb(i)             = cmposit3(ocnfrac(i),lndfrac(i),                 &
                                       lakfrac(i),cice(i),                    &
                                  rb_ocn(i),    rb_lnd(i),    rb_ice(i))
          stress(i)         = cmposit3(ocnfrac(i),lndfrac(i),                 &
                                       lakfrac(i),cice(i),                    &
                              stress_ocn(i),stress_lnd(i),stress_ice(i))
          ffmm(i)           = cmposit3(ocnfrac(i),lndfrac(i),                 &
                                       lakfrac(i),cice(i),                    &
                                ffmm_ocn(i),  ffmm_lnd(i),  ffmm_ice(i))
          ffhh(i)           = cmposit3(ocnfrac(i),lndfrac(i),                 &
                                       lakfrac(i),cice(i),                    &
                                ffhh_ocn(i),  ffhh_lnd(i),  ffhh_ice(i))
          uustar(i)         = cmposit3(ocnfrac(i),lndfrac(i),                 &
                                       lakfrac(i),cice(i),                    &
                              uustar_ocn(i),uustar_lnd(i),uustar_ice(i))
          fm10(i)           = cmposit3(ocnfrac(i),lndfrac(i),                 &
                                       lakfrac(i),cice(i),                    &
                                fm10_ocn(i),  fm10_lnd(i),  fm10_ice(i))
          fh2(i)            = cmposit3(ocnfrac(i),lndfrac(i),                 &
                                       lakfrac(i),cice(i),                    &
                                 fh2_ocn(i),   fh2_lnd(i),   fh2_ice(i))
          tsurf(i)          = cmposit3(ocnfrac(i),lndfrac(i),                 &
                                       lakfrac(i),cice(i),                    &
                               tsurf_ocn(i), tsurf_lnd(i), tsurf_ice(i))
          cmm(i)            = cmposit3(ocnfrac(i),lndfrac(i),                 &
                                       lakfrac(i),cice(i),                    &
                                 cmm_ocn(i),   cmm_lnd(i),   cmm_ice(i))
          chh(i)            = cmposit3(ocnfrac(i),lndfrac(i),                 &
                                       lakfrac(i),cice(i),                    &
                                 chh_ocn(i),   chh_lnd(i),   chh_ice(i))
          gflx(i)           = cmposit3(ocnfrac(i),lndfrac(i),                 &
                                       lakfrac(i),cice(i),                    &
                                gflx_ocn(i),  gflx_lnd(i),  gflx_ice(i))
          ep1d(i)           = cmposit3(ocnfrac(i),lndfrac(i),                 &
                                       lakfrac(i),cice(i),                    &
                                ep1d_ocn(i),  ep1d_lnd(i),  ep1d_ice(i))
          weasd(i)           = cmposit3(ocnfrac(i),lndfrac(i),                &
                                       lakfrac(i),cice(i),                    &
                                weasd(i), weasd_lnd(i), weasd_ice(i))
          snowd(i)           = cmposit3(ocnfrac(i),lndfrac(i),                &
                                       lakfrac(i),cice(i),                    &
                               snowd_ocn(i), snowd_lnd(i), snowd_ice(i))
          tprcp(i)           = cmposit3(ocnfrac(i),lndfrac(i),                &
                                       lakfrac(i),cice(i),                    &
                               tprcp_ocn(i), tprcp_lnd(i), tprcp_ice(i))

  ! Two-way composites (fields already composited in sfc_sice)
  ! Three-way composites when coupled
          if(cplflx .and. ilak(i) == 0) then  ! Lakes in coupled mode use sice?
            evap(i)           = cmposit3(ocnfrac(i),lndfrac(i),               &
                                         lakfrac(i),cice(i),                  &
                                  evap_ocn(i),  evap_lnd(i),  evap_ice(i))
            hflx(i)           = cmposit3(ocnfrac(i),lndfrac(i),               &
                                         lakfrac(i),cice(i),                  &
                                  hflx_ocn(i),  hflx_lnd(i),  hflx_ice(i))
            qss(i)            = cmposit3(ocnfrac(i),lndfrac(i),               &
                                         lakfrac(i),cice(i),                  &
                                   qss_ocn(i),   qss_lnd(i),   qss_ice(i))
            tsfc(i)           = cmposit3(ocnfrac(i),lndfrac(i),               &
                                         lakfrac(i),cice(i),                  &
                                  tsfc_ocn(i),  tsfc_lnd(i),  tsfc_ice(i))
          else
            evap(i)           = cmposit2(ocnfrac(i),lndfrac(i),               &
                                         lakfrac(i),cice(i),                  &
                                  evap_ocn(i),  evap_lnd(i),  evap_ice(i))
            hflx(i)           = cmposit2(ocnfrac(i),lndfrac(i),               &
                                         lakfrac(i),cice(i),                  &
                                  hflx_ocn(i),  hflx_lnd(i),  hflx_ice(i))
            qss(i)            = cmposit2(ocnfrac(i),lndfrac(i),               &
                                         lakfrac(i),cice(i),                  &
                                   qss_ocn(i),   qss_lnd(i),   qss_ice(i))
            tsfc(i)           = cmposit2(ocnfrac(i),lndfrac(i), &
                                         lakfrac(i),cice(i),            &
                                  tsfc_ocn(i),  tsfc_lnd(i),  tsfc_ice(i))
            if(iice(i) == 1 .and. .not. cplflx) then
              cmm(i)           = cmm_ice(i)
              chh(i)           = chh_ice(i)
              gflx(i)          = gflx_ice(i)
              ep1d(i)          = ep1d_ice(i)
              weasd(i)         = weasd_ice(i)
              snowd(i)         = snowd_ice(i)
            end if
          endif

          zorll(i) = zorl_lnd(i)
          zorlo(i) = zorl_ocn(i)

          if (idry(i)==1) tsfcl(i) = tsfc_lnd(i)
          if (iwet(i)==1) then
            tsfco(i) = tsfc_ocn(i)
            tisfc(i) = tsfc_ice(i)
          end if

        end do

  ! --- compositing done


        do i=1,im
          epi(i)     = ep1d(i)
          gfluxi(i)  = gflx(i)
          t1(i)      = tgrs_1(i)
          q1(i)      = qgrs_1(i)
          u1(i)      = ugrs_1(i)
          v1(i)      = vgrs_1(i)
        enddo

        if (cplflx .or. cplwav) then
          do i=1,im
            u10mi_cpl   (i) = u10m(i)
            v10mi_cpl   (i) = v10m(i)
          enddo
        endif

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
            tsfci_cpl   (i) = tsfc(i)
            psurfi_cpl  (i) = pgr(i)
          enddo

  !  ---  estimate mean albedo for ocean point without ice cover and apply
  !       them to net SW heat fluxes

          do i=1,im
            if(lndfrac(i) < 1.) then ! Not 100% land
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
          enddo
        endif

!  --- ...  total runoff is composed of drainage into water table and
!           runoff at the surface and is accumulated in unit of meters
        if (lssav) then
          tem = dtf * 0.001
          do i=1,im
            runoff(i)  = runoff(i)  + (drain(i)+runof(i)) * tem
            srunoff(i) = srunoff(i) + runof(i) * tem
          enddo
        endif

      end subroutine GFS_surface_generic_post_run

      real function cmposit2(frac_ocn,frac_dry,frac_lak,frac_ice,ocnval,lndval,iceval)
  ! --- 2-way compositing (use with ice/non-ice composited variables)
      implicit none
      real(kind=kind_phys),intent(IN) :: frac_ocn,frac_dry,frac_lak,frac_ice,ocnval,lndval,iceval
      real(kind=kind_phys)            :: frac_wet

      frac_wet=max(frac_lak,frac_ocn)
      if (frac_ice.eq.0.) then
        cmposit2 = frac_dry*lndval + frac_wet*ocnval
      else
        cmposit2 = frac_dry*lndval + frac_wet*iceval
      end if
      return
      end function cmposit2

      real function cmposit3(frac_ocn,frac_dry,frac_lak,frac_ice,ocnval,lndval,iceval)
  ! --- 3-way compositing
      implicit none
      real(kind=kind_phys),intent(IN) :: frac_ocn,frac_dry,frac_lak,frac_ice,ocnval,lndval,iceval

      cmposit3 = frac_dry*lndval + frac_ice*iceval + (1.-frac_dry-frac_ice)*ocnval
      return
      end function cmposit3

      end module GFS_surface_generic_post

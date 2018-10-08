!> \file GFS_PBL_generic.F90
!!  Contains code related to PBL schemes to be used within the GFS physics suite.

      module GFS_PBL_generic_pre

      contains

      subroutine GFS_PBL_generic_pre_init ()
      end subroutine GFS_PBL_generic_pre_init

      subroutine GFS_PBL_generic_pre_finalize()
      end subroutine GFS_PBL_generic_pre_finalize

!> \brief This scheme sets up the vertically diffused tracer array for any PBL scheme based on the microphysics scheme chosen
#if 0
!! \section arg_table_GFS_PBL_generic_pre_run Argument Table
!! | local_name                   | standard_name                                          | long_name                                                                           | units         | rank | type      |    kind   | intent | optional |
!! |------------------------------|--------------------------------------------------------|-------------------------------------------------------------------------------------|---------------|------|-----------|-----------|--------|----------|
!! | im                           | horizontal_loop_extent                                 | horizontal loop extent                                                              | count         |    0 | integer   |           | in     | F        |
!! | levs                         | vertical_dimension                                     | vertical layer dimension                                                            | count         |    0 | integer   |           | in     | F        |
!! | nvdiff                       | number_of_vertical_diffusion_tracers                   | number of tracers to diffuse vertically                                             | count         |    0 | integer   |           | in     | F        |
!! | ntrac                        | number_of_tracers                                      | number of tracers                                                                   | count         |    0 | integer   |           | in     | F        |
!! | imp_physics                  | flag_for_microphysics_scheme                           | choice of microphysics scheme                                                       | flag          |    0 | integer   |           | in     | F        |
!! | imp_physics_gfdl             | flag_for_gfdl_microphysics_scheme                      | choice of GFDL microphysics scheme                                                  | flag          |    0 | integer   |           | in     | F        |
!! | imp_physics_thompson         | flag_for_thompson_microphysics_scheme                  | choice of Thompson microphysics scheme                                              | flag          |    0 | integer   |           | in     | F        |
!! | imp_physics_wsm6             | flag_for_wsm6_microphysics_scheme                      | choice of WSM6 microphysics scheme                                                  | flag          |    0 | integer   |           | in     | F        |
!! | ltaerosol                    | flag_for_aerosol_physics                               | flag for aerosol physics                                                            | flag          |    0 | logical   |           | in     | F        |
!! | qgrs                         | tracer_concentration                                   | model layer mean tracer concentration                                               | kg kg-1       |    3 | real      | kind_phys | in     | F        |
!! | qgrs_water_vapor             | water_vapor_specific_humidity                          | water vapor specific humidity                                                       | kg kg-1       |    2 | real      | kind_phys | in     | F        |
!! | qgrs_liquid_cloud            | cloud_condensed_water_mixing_ratio                     | moist (dry+vapor, no condensates) mixing ratio of cloud water (condensate)          | kg kg-1       |    2 | real      | kind_phys | in     | F        |
!! | qgrs_ice_cloud               | ice_water_mixing_ratio                                 | moist (dry+vapor, no condensates) mixing ratio of ice water                         | kg kg-1       |    2 | real      | kind_phys | in     | F        |
!! | qgrs_ozone                   | ozone_mixing_ratio                                     | ozone mixing ratio                                                                  | kg kg-1       |    2 | real      | kind_phys | in     | F        |
!! | qgrs_cloud_droplet_num_conc  | cloud_droplet_number_concentration                     | number concentration of cloud droplets (liquid)                                     | kg-1          |    2 | real      | kind_phys | in     | F        |
!! | qgrs_cloud_ice_num_conc      | ice_number_concentration                               | number concentration of ice                                                         | kg-1          |    2 | real      | kind_phys | in     | F        |
!! | qgrs_water_aer_num_conc      | water_friendly_aerosol_number_concentration            | number concentration of water-friendly aerosols                                     | kg-1          |    2 | real      | kind_phys | in     | F        |
!! | qgrs_ice_aer_num_conc        | ice_friendly_aerosol_number_concentration              | number concentration of ice-friendly aerosols                                       | kg-1          |    2 | real      | kind_phys | in     | F        |
!! | qgrs_rain                    | rain_water_mixing_ratio                                | moist (dry+vapor, no condensates) mixing ratio of rain water                        | kg kg-1       |    2 | real      | kind_phys | in     | F        |
!! | qgrs_snow                    | snow_water_mixing_ratio                                | moist (dry+vapor, no condensates) mixing ratio of snow water                        | kg kg-1       |    2 | real      | kind_phys | in     | F        |
!! | qgrs_graupel                 | graupel_mixing_ratio                                   | moist (dry+vapor, no condensates) mixing ratio of graupel                           | kg kg-1       |    2 | real      | kind_phys | in     | F        |
!! | vdftra                       | vertically_diffused_tracer_concentration               | tracer concentration diffused by PBL scheme                                         | kg kg-1       |    3 | real      | kind_phys | inout  | F        |
!! | errmsg                       | ccpp_error_message                                     | error message for error handling in CCPP                                            | none          |    0 | character | len=*     | out    | F        |
!! | errflg                       | ccpp_error_flag                                        | error flag for error handling in CCPP                                               | flag          |    0 | integer   |           | out    | F        |
!!
#endif
      subroutine GFS_PBL_generic_pre_run (im, levs, nvdiff, ntrac, imp_physics, imp_physics_gfdl, imp_physics_thompson, &
        imp_physics_wsm6, ltaerosol, qgrs, qgrs_water_vapor, qgrs_liquid_cloud, qgrs_ice_cloud, qgrs_ozone, &
        qgrs_cloud_droplet_num_conc, qgrs_cloud_ice_num_conc, qgrs_water_aer_num_conc, qgrs_ice_aer_num_conc, qgrs_rain, &
        qgrs_snow, qgrs_graupel, vdftra, errmsg, errflg)

      use machine, only : kind_phys

      implicit none

      integer, intent(in) :: im, levs, nvdiff, ntrac, imp_physics, imp_physics_gfdl, imp_physics_thompson, imp_physics_wsm6
      logical, intent(in) :: ltaerosol

      real(kind=kind_phys), dimension(im, levs, ntrac), intent(in) :: qgrs
      real(kind=kind_phys), dimension(im, levs), intent(in) :: qgrs_water_vapor, qgrs_liquid_cloud, qgrs_ice_cloud, &
        qgrs_ozone, qgrs_cloud_droplet_num_conc, qgrs_cloud_ice_num_conc, qgrs_water_aer_num_conc, qgrs_ice_aer_num_conc, &
        qgrs_rain, qgrs_snow, qgrs_graupel
      real(kind=kind_phys), dimension(im, levs, nvdiff), intent(inout) :: vdftra

      character(len=*), intent(out) :: errmsg
      integer, intent(out) :: errflg

      !local variables
      integer :: i, k

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      if(nvdiff == ntrac) then
        vdftra = qgrs
      else
        if (imp_physics == imp_physics_wsm6) then
  ! WSM6
          do k=1,levs
            do i=1,im
              vdftra(i,k,1) = qgrs_water_vapor(i,k)
              vdftra(i,k,2) = qgrs_liquid_cloud(i,k)
              vdftra(i,k,3) = qgrs_ice_cloud(i,k)
              vdftra(i,k,4) = qgrs_ozone(i,k)
            enddo
          enddo
        elseif (imp_physics == imp_physics_thompson) then
  ! Thompson
          if(ltaerosol) then
            do k=1,levs
              do i=1,im
                vdftra(i,k,1) = qgrs_water_vapor(i,k)
                vdftra(i,k,2) = qgrs_liquid_cloud(i,k)
                vdftra(i,k,3) = qgrs_ice_cloud(i,k)
                vdftra(i,k,4) = qgrs_cloud_droplet_num_conc(i,k)
                vdftra(i,k,5) = qgrs_cloud_ice_num_conc(i,k)
                vdftra(i,k,6) = qgrs_ozone(i,k)
                vdftra(i,k,7) = qgrs_water_aer_num_conc(i,k)
                vdftra(i,k,8) = qgrs_ice_aer_num_conc(i,k)
              enddo
            enddo
          else
            do k=1,levs
              do i=1,im
                vdftra(i,k,1) = qgrs_water_vapor(i,k)
                vdftra(i,k,2) = qgrs_liquid_cloud(i,k)
                vdftra(i,k,3) = qgrs_ice_cloud(i,k)
                vdftra(i,k,4) = qgrs_cloud_ice_num_conc(i,k)
                vdftra(i,k,5) = qgrs_ozone(i,k)
              enddo
            enddo
          endif
  !
        elseif (imp_physics == imp_physics_gfdl) then
  ! GFDL MP
          do k=1,levs
            do i=1,im
              vdftra(i,k,1) = qgrs_water_vapor(i,k)
              vdftra(i,k,2) = qgrs_liquid_cloud(i,k)
              vdftra(i,k,3) = qgrs_ice_cloud(i,k)
              vdftra(i,k,4) = qgrs_rain(i,k)
              vdftra(i,k,5) = qgrs_snow(i,k)
              vdftra(i,k,6) = qgrs_graupel(i,k)
              vdftra(i,k,7) = qgrs_ozone(i,k)
            enddo
          enddo
        endif
      endif

    end subroutine GFS_PBL_generic_pre_run

    end module GFS_PBL_generic_pre

    module GFS_PBL_generic_post

    contains

    subroutine GFS_PBL_generic_post_init ()
    end subroutine GFS_PBL_generic_post_init

    subroutine GFS_PBL_generic_post_finalize ()
    end subroutine GFS_PBL_generic_post_finalize


#if 0
!> \section arg_table_GFS_PBL_generic_post_run Argument Table
!! | local_name                   | standard_name                                                                     | long_name                                                                                   | units         | rank | type      |    kind   | intent | optional |
!! |------------------------------|-----------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------|---------------|------|-----------|-----------|--------|----------|
!! | im                           | horizontal_loop_extent                                                            | horizontal loop extent                                                                      | count         |    0 | integer   |           | in     | F        |
!! | levs                         | vertical_dimension                                                                | vertical layer dimension                                                                    | count         |    0 | integer   |           | in     | F        |
!! | nvdiff                       | number_of_vertical_diffusion_tracers                                              | number of tracers to diffuse vertically                                                     | count         |    0 | integer   |           | in     | F        |
!! | ntrac                        | number_of_tracers                                                                 | number of tracers                                                                           | count         |    0 | integer   |           | in     | F        |
!! | ntoz                         | index_for_ozone                                                                   | tracer index for ozone mixing ratio                                                         | index         |    0 | integer   |           | in     | F        |
!! | imp_physics                  | flag_for_microphysics_scheme                                                      | choice of microphysics scheme                                                               | flag          |    0 | integer   |           | in     | F        |
!! | imp_physics_gfdl             | flag_for_gfdl_microphysics_scheme                                                 | choice of GFDL microphysics scheme                                                          | flag          |    0 | integer   |           | in     | F        |
!! | imp_physics_thompson         | flag_for_thompson_microphysics_scheme                                             | choice of Thompson microphysics scheme                                                      | flag          |    0 | integer   |           | in     | F        |
!! | imp_physics_wsm6             | flag_for_wsm6_microphysics_scheme                                                 | choice of WSM6 microphysics scheme                                                          | flag          |    0 | integer   |           | in     | F        |
!! | ltaerosol                    | flag_for_aerosol_physics                                                          | flag for aerosol physics                                                                    | flag          |    0 | logical   |           | in     | F        |
!! | cplflx                       | flag_for_flux_coupling                                                            | flag controlling cplflx collection (default off)                                            | flag          |    0 | logical   |           | in     | F        |
!! | lssav                        | flag_diagnostics                                                                  | logical flag for storing diagnostics                                                        | flag          |    0 | logical   |           | in     | F        |
!! | ldiag3d                      | flag_diagnostics_3D                                                               | flag for 3d diagnostic fields                                                               | flag          |    0 | logical   |           | in     | F        |
!! | lsidea                       | flag_idealized_physics                                                            | flag for idealized physics                                                                  | flag          |    0 | logical   |           | in     | F        |
!! | hybedmf                      | flag_for_hedmf                                                                    | flag for hybrid edmf pbl scheme (moninedmf)                                                 | flag          |    0 | logical   |           | in     | F        |
!! | do_shoc                      | flag_for_shoc                                                                     | flag for SHOC                                                                               | flag          |    0 | logical   |           | in     | F        |
!! | dvdftra                      | tendency_of_vertically_diffused_tracer_concentration                              | updated tendency of the tracers due to vertical diffusion in PBL scheme                     | kg kg-1 s-1   |    3 | real      | kind_phys | in     | F        |
!! | dusfc1                       | instantaneous_surface_x_momentum_flux                                             | surface momentum flux in the x-direction valid for current call                             | Pa            |    1 | real      | kind_phys | in     | F        |
!! | dvsfc1                       | instantaneous_surface_y_momentum_flux                                             | surface momentum flux in the y-direction valid for current call                             | Pa            |    1 | real      | kind_phys | in     | F        |
!! | dtsfc1                       | instantaneous_surface_upward_sensible_heat_flux                                   | surface upward sensible heat flux valid for current call                                    | W m-2         |    1 | real      | kind_phys | in     | F        |
!! | dqsfc1                       | instantaneous_surface_upward_latent_heat_flux                                     | surface upward latent heat flux valid for current call                                      | W m-2         |    1 | real      | kind_phys | in     | F        |
!! | dtf                          | time_step_for_dynamics                                                            | dynamics timestep                                                                           | s             |    0 | real      | kind_phys | in     | F        |
!! | dudt                         | tendency_of_x_wind_due_to_model_physics                                           | updated tendency of the x wind                                                              | m s-2         |    2 | real      | kind_phys | in     | F        |
!! | dvdt                         | tendency_of_y_wind_due_to_model_physics                                           | updated tendency of the y wind                                                              | m s-2         |    2 | real      | kind_phys | in     | F        |
!! | dtdt                         | tendency_of_air_temperature_due_to_model_physics                                  | updated tendency of the temperature                                                         | K s-1         |    2 | real      | kind_phys | in     | F        |
!! | htrsw                        | tendency_of_air_temperature_due_to_shortwave_heating_on_radiation_timestep        | total sky sw heating rate                                                                   | K s-1         |    2 | real      | kind_phys | in     | F        |
!! | htrlw                        | tendency_of_air_temperature_due_to_longwave_heating_on_radiation_timestep         | total sky lw heating rate                                                                   | K s-1         |    2 | real      | kind_phys | in     | F        |
!! | xmu                          | zenith_angle_temporal_adjustment_factor_for_shortwave_fluxes                      | zenith angle temporal adjustment factor for shortwave                                       | none          |    1 | real      | kind_phys | in     | F        |
!! | dqdt                         | tendency_of_tracers_due_to_model_physics                                          | updated tendency of the tracers due to model physics                                        | kg kg-1 s-1   |    3 | real      | kind_phys | inout  | F        |
!! | dqdt_water_vapor             | tendency_of_water_vapor_specific_humidity_due_to_model_physics                    | water vapor specific humidity tendency due to model physics                                 | kg kg-1 s-1   |    2 | real      | kind_phys | inout  | F        |
!! | dqdt_liquid_cloud            | tendency_of_liquid_cloud_water_mixing_ratio_due_to_model_physics                  | cloud condensed water mixing ratio tendency due to model physics                            | kg kg-1 s-1   |    2 | real      | kind_phys | inout  | F        |
!! | dqdt_ice_cloud               | tendency_of_ice_cloud_water_mixing_ratio_due_to_model_physics                     | cloud condensed water mixing ratio tendency due to model physics                            | kg kg-1 s-1   |    2 | real      | kind_phys | inout  | F        |
!! | dqdt_ozone                   | tendency_of_ozone_mixing_ratio_due_to_model_physics                               | ozone mixing ratio tendency due to model physics                                            | kg kg-1 s-1   |    2 | real      | kind_phys | inout  | F        |
!! | dqdt_cloud_droplet_num_conc  | tendency_of_cloud_droplet_number_concentration_due_to_model_physics               | number concentration of cloud droplets (liquid) tendency due to model physics               | kg-1 s-1      |    2 | real      | kind_phys | inout  | F        |
!! | dqdt_ice_num_conc            | tendency_of_ice_number_concentration_due_to_model_physics                         | number concentration of ice tendency due to model physics                                   | kg-1 s-1      |    2 | real      | kind_phys | inout  | F        |
!! | dqdt_water_aer_num_conc      | tendency_of_water_friendly_aerosol_number_concentration_due_to_model_physics      | number concentration of water-friendly aerosols tendency due to model physics               | kg-1 s-1      |    2 | real      | kind_phys | inout  | F        |
!! | dqdt_ice_aer_num_conc        | tendency_of_ice_friendly_aerosol_number_concentration_due_to_model_physics        | number concentration of ice-friendly aerosols tendency due to model physics                 | kg-1 s-1      |    2 | real      | kind_phys | inout  | F        |
!! | dqdt_rain                    | tendency_of_rain_water_mixing_ratio_due_to_model_physics                          | moist (dry+vapor, no condensates) mixing ratio of rain water tendency due to model physics  | kg kg-1 s-1   |    2 | real      | kind_phys | inout  | F        |
!! | dqdt_snow                    | tendency_of_snow_water_mixing_ratio_due_to_model_physics                          | moist (dry+vapor, no condensates) mixing ratio of snow water tendency due to model physics  | kg kg-1 s-1   |    2 | real      | kind_phys | inout  | F        |
!! | dqdt_graupel                 | tendency_of_graupel_mixing_ratio_due_to_model_physics                             | moist (dry+vapor, no condensates) mixing ratio of graupel tendency due to model physics     | kg kg-1 s-1   |    2 | real      | kind_phys | inout  | F        |
!! | dusfc_cpl                    | cumulative_surface_x_momentum_flux_for_coupling_multiplied_by_timestep            | cumulative sfc u momentum flux multiplied by timestep                                       | Pa s          |    1 | real      | kind_phys | inout  | F        |
!! | dvsfc_cpl                    | cumulative_surface_y_momentum_flux_for_coupling_multiplied_by_timestep            | cumulative sfc v momentum flux multiplied by timestep                                       | Pa s          |    1 | real      | kind_phys | inout  | F        |
!! | dtsfc_cpl                    | cumulative_surface_upward_sensible_heat_flux_for_coupling_multiplied_by_timestep  | cumulative sfc sensible heat flux multiplied by timestep                                    | W m-2 s       |    1 | real      | kind_phys | inout  | F        |
!! | dqsfc_cpl                    | cumulative_surface_upward_latent_heat_flux_for_coupling_multiplied_by_timestep    | cumulative sfc latent heat flux multiplied by timestep                                      | W m-2 s       |    1 | real      | kind_phys | inout  | F        |
!! | dusfci_cpl                   | instantaneous_surface_x_momentum_flux_for_coupling                                | instantaneous sfc u momentum flux                                                           | Pa            |    1 | real      | kind_phys | inout  | F        |
!! | dvsfci_cpl                   | instantaneous_surface_y_momentum_flux_for_coupling                                | instantaneous sfc v momentum flux                                                           | Pa            |    1 | real      | kind_phys | inout  | F        |
!! | dtsfci_cpl                   | instantaneous_surface_upward_sensible_heat_flux_for_coupling                      | instantaneous sfc sensible heat flux                                                        | W m-2         |    1 | real      | kind_phys | inout  | F        |
!! | dqsfci_cpl                   | instantaneous_surface_upward_latent_heat_flux_for_coupling                        | instantaneous sfc latent heat flux                                                          | W m-2         |    1 | real      | kind_phys | inout  | F        |
!! | dusfc_diag                   | cumulative_surface_x_momentum_flux_for_diag_multiplied_by_timestep                | cumulative sfc x momentum flux multiplied by timestep                                       | Pa s          |    1 | real      | kind_phys | inout  | F        |
!! | dvsfc_diag                   | cumulative_surface_y_momentum_flux_for_diag_multiplied_by_timestep                | cumulative sfc y momentum flux multiplied by timestep                                       | Pa s          |    1 | real      | kind_phys | inout  | F        |
!! | dtsfc_diag                   | cumulative_surface_upward_sensible_heat_flux_for_diag_multiplied_by_timestep      | cumulative sfc sensible heat flux multiplied by timestep                                    | W m-2 s       |    1 | real      | kind_phys | inout  | F        |
!! | dqsfc_diag                   | cumulative_surface_upward_latent_heat_flux_for_diag_multiplied_by_timestep        | cumulative sfc latent heat flux multiplied by timestep                                      | W m-2 s       |    1 | real      | kind_phys | inout  | F        |
!! | dusfci_diag                  | instantaneous_surface_x_momentum_flux_for_diag                                    | instantaneous sfc x momentum flux multiplied by timestep                                    | Pa            |    1 | real      | kind_phys | inout  | F        |
!! | dvsfci_diag                  | instantaneous_surface_y_momentum_flux_for_diag                                    | instantaneous sfc y momentum flux multiplied by timestep                                    | Pa            |    1 | real      | kind_phys | inout  | F        |
!! | dtsfci_diag                  | instantaneous_surface_upward_sensible_heat_flux_for_diag                          | instantaneous sfc sensible heat flux multiplied by timestep                                 | W m-2         |    1 | real      | kind_phys | inout  | F        |
!! | dqsfci_diag                  | instantaneous_surface_upward_latent_heat_flux_for_diag                            | instantaneous sfc latent heat flux multiplied by timestep                                   | W m-2         |    1 | real      | kind_phys | inout  | F        |
!! | dt3dt                        | cumulative_change_in_temperature_due_to_PBL                                       | cumulative change in temperature due to PBL                                                 | K             |    2 | real      | kind_phys | inout  | F        |
!! | du3dt_PBL                    | cumulative_change_in_x_wind_due_to_PBL                                            | cumulative change in x wind due to PBL                                                      | m s-1         |    2 | real      | kind_phys | inout  | F        |
!! | du3dt_OGWD                   | cumulative_change_in_x_wind_due_to_orographic_gravity_wave_drag                   | cumulative change in x wind due to orographic gravity wave drag                             | m s-1         |    2 | real      | kind_phys | inout  | F        |
!! | dv3dt_PBL                    | cumulative_change_in_y_wind_due_to_PBL                                            | cumulative change in y wind due to PBL                                                      | m s-1         |    2 | real      | kind_phys | inout  | F        |
!! | dv3dt_OGWD                   | cumulative_change_in_y_wind_due_to_orographic_gravity_wave_drag                   | cumulative change in y wind due to orographic gravity wave drag                             | m s-1         |    2 | real      | kind_phys | inout  | F        |
!! | dq3dt                        | cumulative_change_in_water_vapor_specific_humidity_due_to_PBL                     | cumulative change in water vapor specific humidity due to PBL                               | kg kg-1       |    2 | real      | kind_phys | inout  | F        |
!! | dq3dt_ozone                  | cumulative_change_in_ozone_mixing_ratio_due_to_PBL                                | cumulative change in ozone mixing ratio due to PBL                                          | kg kg-1       |    2 | real      | kind_phys | inout  | F        |
!! | errmsg                       | ccpp_error_message                                                                | error message for error handling in CCPP                                                    | none          |    0 | character | len=*     | out    | F        |
!! | errflg                       | ccpp_error_flag                                                                   | error flag for error handling in CCPP                                                       | flag          |    0 | integer   |           | out    | F        |
!!
#endif
      subroutine GFS_PBL_generic_post_run (im, levs, nvdiff, ntrac, ntoz, imp_physics, imp_physics_gfdl, imp_physics_thompson, &
        imp_physics_wsm6, ltaerosol, cplflx, lssav, ldiag3d, lsidea, hybedmf, do_shoc, dvdftra, dusfc1, dvsfc1, dtsfc1, dqsfc1, &
        dtf, dudt, dvdt, dtdt, htrsw, htrlw, xmu,&
        dqdt, dqdt_water_vapor, dqdt_liquid_cloud, dqdt_ice_cloud, dqdt_ozone, dqdt_cloud_droplet_num_conc, dqdt_ice_num_conc,&
        dqdt_water_aer_num_conc, dqdt_ice_aer_num_conc, dqdt_rain, dqdt_snow, dqdt_graupel, dusfc_cpl, dvsfc_cpl, dtsfc_cpl, &
        dqsfc_cpl, dusfci_cpl, dvsfci_cpl, dtsfci_cpl, dqsfci_cpl, dusfc_diag, dvsfc_diag, dtsfc_diag, dqsfc_diag, &
        dusfci_diag, dvsfci_diag, dtsfci_diag, dqsfci_diag, dt3dt, du3dt_PBL, du3dt_OGWD, dv3dt_PBL, dv3dt_OGWD, dq3dt, &
        dq3dt_ozone, errmsg, errflg)

      use machine,               only: kind_phys

      implicit none

      integer, intent(in) :: im, levs, nvdiff, ntrac, ntoz, imp_physics, imp_physics_gfdl, imp_physics_thompson, imp_physics_wsm6
      logical, intent(in) :: ltaerosol, cplflx, lssav, ldiag3d, lsidea, hybedmf, do_shoc

      real(kind=kind_phys), intent(in) :: dtf
      real(kind=kind_phys), dimension(im, levs, nvdiff), intent(in) :: dvdftra
      real(kind=kind_phys), dimension(im), intent(in) :: dusfc1, dvsfc1, dtsfc1, dqsfc1, xmu
      real(kind=kind_phys), dimension(im, levs), intent(in) :: dudt, dvdt, dtdt, htrsw, htrlw

      real(kind=kind_phys), dimension(im, levs, ntrac), intent(inout) :: dqdt
      real(kind=kind_phys), dimension(im, levs), intent(inout) :: dqdt_water_vapor, dqdt_liquid_cloud, dqdt_ice_cloud, dqdt_ozone, &
        dqdt_cloud_droplet_num_conc, dqdt_ice_num_conc, dqdt_water_aer_num_conc, dqdt_ice_aer_num_conc, dqdt_rain,&
        dqdt_snow, dqdt_graupel, dt3dt, du3dt_PBL, du3dt_OGWD, dv3dt_PBL, dv3dt_OGWD, dq3dt, dq3dt_ozone
      real(kind=kind_phys), dimension(im), intent(inout) :: dusfc_cpl, dvsfc_cpl, dtsfc_cpl, dqsfc_cpl, dusfci_cpl, dvsfci_cpl, &
        dtsfci_cpl, dqsfci_cpl, dusfc_diag, dvsfc_diag, dtsfc_diag, dqsfc_diag, dusfci_diag, dvsfci_diag, dtsfci_diag, dqsfci_diag

      character(len=*), intent(out) :: errmsg
      integer, intent(out) :: errflg

      integer :: i, k
      real(kind=kind_phys) :: tem

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0
!GJF: dvdftra is only used if nvdiff != ntrac or (nvdiff == ntrac .and. )
      if (nvdiff == ntrac .and. (hybedmf .or. do_shoc)) then
        dqdt = dvdftra
      elseif (nvdiff /= ntrac) then
        if (imp_physics == imp_physics_wsm6) then
  ! WSM6
          do k=1,levs
            do i=1,im
              dqdt_water_vapor(i,k)  = dvdftra(i,k,1)
              dqdt_liquid_cloud(i,k) = dvdftra(i,k,2)
              dqdt_ice_cloud(i,k)    = dvdftra(i,k,3)
              dqdt_ozone(i,k)        = dvdftra(i,k,4)
            enddo
          enddo
        elseif (imp_physics == imp_physics_thompson) then
  ! Thompson
          if(ltaerosol) then
            do k=1,levs
              do i=1,im
                dqdt_water_vapor(i,k)             = dvdftra(i,k,1)
                dqdt_liquid_cloud(i,k)            = dvdftra(i,k,2)
                dqdt_ice_cloud(i,k)               = dvdftra(i,k,3)
                dqdt_cloud_droplet_num_conc(i,k)  = dvdftra(i,k,4)
                dqdt_ice_num_conc(i,k)            = dvdftra(i,k,5)
                dqdt_ozone(i,k)                   = dvdftra(i,k,6)
                dqdt_water_aer_num_conc(i,k)      = dvdftra(i,k,7)
                dqdt_ice_aer_num_conc(i,k)        = dvdftra(i,k,8)
              enddo
            enddo
          else
            do k=1,levs
              do i=1,im
                dqdt_water_vapor(i,k)   = dvdftra(i,k,1)
                dqdt_liquid_cloud(i,k)  = dvdftra(i,k,2)
                dqdt_ice_cloud(i,k)     = dvdftra(i,k,3)
                dqdt_ice_num_conc(i,k)  = dvdftra(i,k,4)
                dqdt_ozone(i,k)         = dvdftra(i,k,5)
              enddo
            enddo
          endif
        elseif (imp_physics == imp_physics_gfdl) then
  ! GFDL MP
          do k=1,levs
            do i=1,im
              dqdt_water_vapor(i,k)   = dvdftra(i,k,1)
              dqdt_liquid_cloud(i,k)  = dvdftra(i,k,2)
              dqdt_ice_cloud(i,k)     = dvdftra(i,k,3)
              dqdt_rain(i,k)          = dvdftra(i,k,4)
              dqdt_snow(i,k)          = dvdftra(i,k,5)
              dqdt_graupel(i,k)       = dvdftra(i,k,6)
              dqdt_ozone(i,k)         = dvdftra(i,k,7)
            enddo
          enddo
        endif
      endif ! nvdiff == ntrac

!     if (lprnt) then
!       write(0,*) ' dusfc1=',dusfc1(ipr),' kdt=',kdt,' lat=',lat
!       write(0,*)' dtsfc1=',dtsfc1(ipr)
!       write(0,*)' dqsfc1=',dqsfc1(ipr)
!       write(0,*)' dtdtc=',(dtdt(ipr,k),k=1,15)
!       write(0,*)' dqdtc=',(dqdt(ipr,k,1),k=1,15)
!       print *,' dudtm=',dudt(ipr,:)
!     endif

!  --- ...  coupling insertion

      if (cplflx) then
        do i=1,im
          dusfc_cpl (i) = dusfc_cpl(i) + dusfc1(i)*dtf
          dvsfc_cpl (i) = dvsfc_cpl(i) + dvsfc1(i)*dtf
          dtsfc_cpl (i) = dtsfc_cpl(i) + dtsfc1(i)*dtf
          dqsfc_cpl (i) = dqsfc_cpl(i) + dqsfc1(i)*dtf
          dusfci_cpl(i) = dusfc1(i)
          dvsfci_cpl(i) = dvsfc1(i)
          dtsfci_cpl(i) = dtsfc1(i)
          dqsfci_cpl(i) = dqsfc1(i)
        enddo
      endif
!-------------------------------------------------------lssav if loop ----------
      if (lssav) then
        do i=1,im
          dusfc_diag (i) = dusfc_diag(i) + dusfc1(i)*dtf
          dvsfc_diag (i) = dvsfc_diag(i) + dvsfc1(i)*dtf
          dtsfc_diag (i) = dtsfc_diag(i) + dtsfc1(i)*dtf
          dqsfc_diag (i) = dqsfc_diag(i) + dqsfc1(i)*dtf
          dusfci_diag(i) = dusfc1(i)
          dvsfci_diag(i) = dvsfc1(i)
          dtsfci_diag(i) = dtsfc1(i)
          dqsfci_diag(i) = dqsfc1(i)
        enddo
  !       if (lprnt) then
  !         write(0,*)' dusfc=',dusfc(ipr),' dusfc1=',dusfc1(ipr),' dtf=',
  !    &     dtf,' kdt=',kdt,' lat=',lat
  !       endif

        if (ldiag3d) then
          if (lsidea) then
            dt3dt(1:im,:) = dt3dt(1:im,:) + dtdt(1:im,:)*dtf
          else
            do k=1,levs
              do i=1,im
                tem  = dtdt(i,k) - (htrlw(i,k)+htrsw(i,k)*xmu(i))
                dt3dt(i,k) = dt3dt(i,k) + tem*dtf
              enddo
            enddo
          endif
          do k=1,levs
            do i=1,im
              du3dt_PBL(i,k) = du3dt_PBL(i,k) + dudt(i,k) * dtf
              du3dt_OGWD(i,k) = du3dt_OGWD(i,k) - dudt(i,k) * dtf
              dv3dt_PBL(i,k) = dv3dt_PBL(i,k) + dvdt(i,k) * dtf
              dv3dt_OGWD(i,k) = dv3dt_OGWD(i,k) - dvdt(i,k) * dtf
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
          do k=1,levs
            do i=1,im
              tem  = dqdt_water_vapor(i,k) * dtf
              dq3dt(i,k) = dq3dt(i,k) + tem
            enddo
          enddo
          if (ntoz > 0) then
            do k=1,levs
              do i=1,im
                dq3dt_ozone(i,k) = dq3dt_ozone(i,k) + dqdt_ozone(i,k) * dtf
              enddo
            enddo
          endif
        endif

      endif   ! end if_lssav

      end subroutine GFS_PBL_generic_post_run

      end module GFS_PBL_generic_post

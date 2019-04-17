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
!! | ntqv                         | index_for_water_vapor                                  | tracer index for water vapor (specific humidity)                                    | index         |    0 | integer   |           | in     | F        |
!! | ntcw                         | index_for_liquid_cloud_condensate                      | tracer index for cloud condensate (or liquid water)                                 | index         |    0 | integer   |           | in     | F        |
!! | ntiw                         | index_for_ice_cloud_condensate                         | tracer index for  ice water                                                         | index         |    0 | integer   |           | in     | F        |
!! | ntrw                         | index_for_rain_water                                   | tracer index for rain water                                                         | index         |    0 | integer   |           | in     | F        |
!! | ntsw                         | index_for_snow_water                                   | tracer index for snow water                                                         | index         |    0 | integer   |           | in     | F        |
!! | ntlnc                        | index_for_liquid_cloud_number_concentration            | tracer index for liquid number concentration                                        | index         |    0 | integer   |           | in     | F        |
!! | ntinc                        | index_for_ice_cloud_number_concentration               | tracer index for ice    number concentration                                        | index         |    0 | integer   |           | in     | F        |
!! | ntwa                         | index_for_water_friendly_aerosols                      | tracer index for water friendly aerosol                                             | index         |    0 | integer   |           | in     | F        |
!! | ntia                         | index_for_ice_friendly_aerosols                        | tracer index for ice friendly aerosol                                               | index         |    0 | integer   |           | in     | F        |
!! | ntgl                         | index_for_graupel                                      | tracer index for graupel                                                            | index         |    0 | integer   |           | in     | F        |
!! | ntoz                         | index_for_ozone                                        | tracer index for ozone mixing ratio                                                 | index         |    0 | integer   |           | in     | F        |
!! | ntke                         | index_for_turbulent_kinetic_energy                     | tracer index for turbulent kinetic energy                                           | index         |    0 | integer   |           | in     | F        |
!! | ntkev                        | index_for_turbulent_kinetic_energy_vertical_diffusion_tracer | index for turbulent kinetic energy in the vertically diffused tracer array    | index         |    0 | integer   |           | in     | F        |
!! | imp_physics                  | flag_for_microphysics_scheme                           | choice of microphysics scheme                                                       | flag          |    0 | integer   |           | in     | F        |
!! | imp_physics_gfdl             | flag_for_gfdl_microphysics_scheme                      | choice of GFDL microphysics scheme                                                  | flag          |    0 | integer   |           | in     | F        |
!! | imp_physics_thompson         | flag_for_thompson_microphysics_scheme                  | choice of Thompson microphysics scheme                                              | flag          |    0 | integer   |           | in     | F        |
!! | imp_physics_wsm6             | flag_for_wsm6_microphysics_scheme                      | choice of WSM6 microphysics scheme                                                  | flag          |    0 | integer   |           | in     | F        |
!! | imp_physics_zhao_carr        | flag_for_zhao_carr_microphysics_scheme                 | choice of Zhao-Carr microphysics scheme                                             | flag          |    0 | integer   |           | in     | F        |
!! | cplchm                       | flag_for_chemistry_coupling                            | flag controlling cplchm collection (default off)                                    | flag          |    0 | logical   |           | in     | F        |
!! | ltaerosol                    | flag_for_aerosol_physics                               | flag for aerosol physics                                                            | flag          |    0 | logical   |           | in     | F        |
!! | hybedmf                      | flag_for_hedmf                                         | flag for hybrid edmf pbl scheme (moninedmf)                                         | flag          |    0 | logical   |           | in     | F        |
!! | do_shoc                      | flag_for_shoc                                          | flag for SHOC                                                                       | flag          |    0 | logical   |           | in     | F        |
!! | satmedmf                     | flag_for_scale_aware_TKE_moist_EDMF_PBL                | flag for scale-aware TKE moist EDMF PBL scheme                                      | flag          |    0 | logical   |           | in     | F        |
!! | qgrs                         | tracer_concentration                                   | model layer mean tracer concentration                                               | kg kg-1       |    3 | real      | kind_phys | in     | F        |
!! | vdftra                       | vertically_diffused_tracer_concentration               | tracer concentration diffused by PBL scheme                                         | kg kg-1       |    3 | real      | kind_phys | inout  | F        |
!! | errmsg                       | ccpp_error_message                                     | error message for error handling in CCPP                                            | none          |    0 | character | len=*     | out    | F        |
!! | errflg                       | ccpp_error_flag                                        | error flag for error handling in CCPP                                               | flag          |    0 | integer   |           | out    | F        |
!!
#endif
      subroutine GFS_PBL_generic_pre_run (im, levs, nvdiff, ntrac,                       &
        ntqv, ntcw, ntiw, ntrw, ntsw, ntlnc, ntinc, ntwa, ntia, ntgl, ntoz, ntke, ntkev, &
        imp_physics, imp_physics_gfdl, imp_physics_thompson, imp_physics_wsm6,           &
        imp_physics_zhao_carr, cplchm, ltaerosol, hybedmf, do_shoc, satmedmf,            &
        qgrs, vdftra, errmsg, errflg)

      use machine, only : kind_phys

      implicit none

      integer, intent(in) :: im, levs, nvdiff, ntrac
      integer, intent(in) :: ntqv, ntcw, ntiw, ntrw, ntsw, ntlnc, ntinc, ntwa, ntia, ntgl, ntoz, ntke, ntkev
      integer, intent(in) :: imp_physics, imp_physics_gfdl, imp_physics_thompson, imp_physics_wsm6
      integer, intent(in) :: imp_physics_zhao_carr
      logical, intent(in) :: cplchm, ltaerosol, hybedmf, do_shoc, satmedmf

      real(kind=kind_phys), dimension(im, levs, ntrac), intent(in) :: qgrs
      real(kind=kind_phys), dimension(im, levs, nvdiff), intent(inout) :: vdftra

      character(len=*), intent(out) :: errmsg
      integer, intent(out) :: errflg

      !local variables
      integer :: i, k

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

!DH: dvdftra is only used if nvdiff != ntrac or (nvdiff == ntrac .and. )
      if (nvdiff == ntrac .and. (hybedmf .or. do_shoc .or. satmedmf)) then
        vdftra = qgrs
      else
        if (imp_physics == imp_physics_wsm6) then
  ! WSM6
          do k=1,levs
            do i=1,im
              vdftra(i,k,1) = qgrs(i,k,ntqv)
              vdftra(i,k,2) = qgrs(i,k,ntcw)
              vdftra(i,k,3) = qgrs(i,k,ntiw)
              vdftra(i,k,4) = qgrs(i,k,ntoz)
            enddo
          enddo
        elseif (imp_physics == imp_physics_thompson) then
  ! Thompson
          ! DH* Thompson ntrw and ntsw?
          if(ltaerosol) then
            do k=1,levs
              do i=1,im
                vdftra(i,k,1) = qgrs(i,k,ntqv)
                vdftra(i,k,2) = qgrs(i,k,ntcw)
                vdftra(i,k,3) = qgrs(i,k,ntiw)
                vdftra(i,k,4) = qgrs(i,k,ntlnc)
                vdftra(i,k,5) = qgrs(i,k,ntinc)
                vdftra(i,k,6) = qgrs(i,k,ntoz)
                vdftra(i,k,7) = qgrs(i,k,ntwa)
                vdftra(i,k,8) = qgrs(i,k,ntia)
              enddo
            enddo
          else
            do k=1,levs
              do i=1,im
                vdftra(i,k,1) = qgrs(i,k,ntqv)
                vdftra(i,k,2) = qgrs(i,k,ntcw)
                vdftra(i,k,3) = qgrs(i,k,ntiw)
                vdftra(i,k,4) = qgrs(i,k,ntinc)
                vdftra(i,k,5) = qgrs(i,k,ntoz)
              enddo
            enddo
          endif
  !
        elseif (imp_physics == imp_physics_gfdl) then
  ! GFDL MP
          do k=1,levs
            do i=1,im
              vdftra(i,k,1) = qgrs(i,k,ntqv)
              vdftra(i,k,2) = qgrs(i,k,ntcw)
              vdftra(i,k,3) = qgrs(i,k,ntiw)
              vdftra(i,k,4) = qgrs(i,k,ntrw)
              vdftra(i,k,5) = qgrs(i,k,ntsw)
              vdftra(i,k,6) = qgrs(i,k,ntgl)
              vdftra(i,k,7) = qgrs(i,k,ntoz)
            enddo
          enddo
        elseif (imp_physics == imp_physics_zhao_carr) then
! Zhao/Carr/Sundqvist
          if (cplchm) then
            do k=1,levs
              do i=1,im
                vdftra(i,k,1) = qgrs(i,k,ntqv)
                vdftra(i,k,2) = qgrs(i,k,ntcw)
                vdftra(i,k,3) = qgrs(i,k,ntoz)
              enddo
            enddo
          endif
        endif

        if (satmedmf) then
          do k=1,levs
            do i=1,im
              vdftra(i,k,ntkev) = qgrs(i,k,ntke)
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
!! | ntqv                         | index_for_water_vapor                                                             | tracer index for water vapor (specific humidity)                                            | index         |    0 | integer   |           | in     | F        |
!! | ntcw                         | index_for_liquid_cloud_condensate                                                 | tracer index for cloud condensate (or liquid water)                                         | index         |    0 | integer   |           | in     | F        |
!! | ntiw                         | index_for_ice_cloud_condensate                                                    | tracer index for  ice water                                                                 | index         |    0 | integer   |           | in     | F        |
!! | ntrw                         | index_for_rain_water                                                              | tracer index for rain water                                                                 | index         |    0 | integer   |           | in     | F        |
!! | ntsw                         | index_for_snow_water                                                              | tracer index for snow water                                                                 | index         |    0 | integer   |           | in     | F        |
!! | ntlnc                        | index_for_liquid_cloud_number_concentration                                       | tracer index for liquid number concentration                                                | index         |    0 | integer   |           | in     | F        |
!! | ntinc                        | index_for_ice_cloud_number_concentration                                          | tracer index for ice    number concentration                                                | index         |    0 | integer   |           | in     | F        |
!! | ntwa                         | index_for_water_friendly_aerosols                                                 | tracer index for water friendly aerosol                                                     | index         |    0 | integer   |           | in     | F        |
!! | ntia                         | index_for_ice_friendly_aerosols                                                   | tracer index for ice friendly aerosol                                                       | index         |    0 | integer   |           | in     | F        |
!! | ntgl                         | index_for_graupel                                                                 | tracer index for graupel                                                                    | index         |    0 | integer   |           | in     | F        |
!! | ntoz                         | index_for_ozone                                                                   | tracer index for ozone mixing ratio                                                         | index         |    0 | integer   |           | in     | F        |
!! | ntke                         | index_for_turbulent_kinetic_energy                                                | tracer index for turbulent kinetic energy                                                   | index         |    0 | integer   |           | in     | F        |
!! | ntkev                        | index_for_turbulent_kinetic_energy_vertical_diffusion_tracer                      | index for turbulent kinetic energy in the vertically diffused tracer array                  | index         |    0 | integer   |           | in     | F        |
!! | imp_physics                  | flag_for_microphysics_scheme                                                      | choice of microphysics scheme                                                               | flag          |    0 | integer   |           | in     | F        |
!! | imp_physics_gfdl             | flag_for_gfdl_microphysics_scheme                                                 | choice of GFDL microphysics scheme                                                          | flag          |    0 | integer   |           | in     | F        |
!! | imp_physics_thompson         | flag_for_thompson_microphysics_scheme                                             | choice of Thompson microphysics scheme                                                      | flag          |    0 | integer   |           | in     | F        |
!! | imp_physics_wsm6             | flag_for_wsm6_microphysics_scheme                                                 | choice of WSM6 microphysics scheme                                                          | flag          |    0 | integer   |           | in     | F        |
!! | imp_physics_zhao_carr        | flag_for_zhao_carr_microphysics_scheme                                            | choice of Zhao-Carr microphysics scheme                                                     | flag          |    0 | integer   |           | in     | F        |
!! | ltaerosol                    | flag_for_aerosol_physics                                                          | flag for aerosol physics                                                                    | flag          |    0 | logical   |           | in     | F        |
!! | cplflx                       | flag_for_flux_coupling                                                            | flag controlling cplflx collection (default off)                                            | flag          |    0 | logical   |           | in     | F        |
!! | cplchm                       | flag_for_chemistry_coupling                                                       | flag controlling cplchm collection (default off)                                            | flag          |    0 | logical   |           | in     | F        |
!! | lssav                        | flag_diagnostics                                                                  | logical flag for storing diagnostics                                                        | flag          |    0 | logical   |           | in     | F        |
!! | ldiag3d                      | flag_diagnostics_3D                                                               | flag for 3d diagnostic fields                                                               | flag          |    0 | logical   |           | in     | F        |
!! | lsidea                       | flag_idealized_physics                                                            | flag for idealized physics                                                                  | flag          |    0 | logical   |           | in     | F        |
!! | hybedmf                      | flag_for_hedmf                                                                    | flag for hybrid edmf pbl scheme (moninedmf)                                                 | flag          |    0 | logical   |           | in     | F        |
!! | do_shoc                      | flag_for_shoc                                                                     | flag for SHOC                                                                               | flag          |    0 | logical   |           | in     | F        |
!! | satmedmf                     | flag_for_scale_aware_TKE_moist_EDMF_PBL                                           | flag for scale-aware TKE moist EDMF PBL scheme                                              | flag          |    0 | logical   |           | in     | F        |
!! | shinhong                     | flag_for_scale_aware_Shinhong_PBL                                                 | flag for scale-aware Shinhong PBL scheme                                                    | flag          |    0 | logical   |           | in     | F        |
!! | do_ysu                       | flag_for_ysu                                                                      | flag for YSU PBL scheme                                                                     | flag          |    0 | logical   |           | in     | F        |
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
      subroutine GFS_PBL_generic_post_run (im, levs, nvdiff, ntrac,                                                            &
        ntqv, ntcw, ntiw, ntrw, ntsw, ntlnc, ntinc, ntwa, ntia, ntgl, ntoz, ntke, ntkev,                                       &
        imp_physics, imp_physics_gfdl, imp_physics_thompson, imp_physics_wsm6, imp_physics_zhao_carr,                          &
        ltaerosol, cplflx, cplchm, lssav, ldiag3d, lsidea, hybedmf, do_shoc, satmedmf, shinhong, do_ysu,                       &
        dvdftra, dusfc1, dvsfc1, dtsfc1, dqsfc1, dtf, dudt, dvdt, dtdt, htrsw, htrlw, xmu,                                     &
        dqdt, dusfc_cpl, dvsfc_cpl, dtsfc_cpl,                                                                                 &
        dqsfc_cpl, dusfci_cpl, dvsfci_cpl, dtsfci_cpl, dqsfci_cpl, dusfc_diag, dvsfc_diag, dtsfc_diag, dqsfc_diag,             &
        dusfci_diag, dvsfci_diag, dtsfci_diag, dqsfci_diag, dt3dt, du3dt_PBL, du3dt_OGWD, dv3dt_PBL, dv3dt_OGWD, dq3dt,        &
        dq3dt_ozone, errmsg, errflg)

      use machine,               only: kind_phys

      implicit none

      integer, intent(in) :: im, levs, nvdiff, ntrac
      integer, intent(in) :: ntqv, ntcw, ntiw, ntrw, ntsw, ntlnc, ntinc, ntwa, ntia, ntgl, ntoz, ntke, ntkev
      integer, intent(in) :: imp_physics, imp_physics_gfdl, imp_physics_thompson, imp_physics_wsm6
      integer, intent(in) :: imp_physics_zhao_carr
      logical, intent(in) :: ltaerosol, cplflx, cplchm, lssav, ldiag3d, lsidea
      logical, intent(in) :: hybedmf, do_shoc, satmedmf, shinhong, do_ysu

      real(kind=kind_phys), intent(in) :: dtf
      real(kind=kind_phys), dimension(im, levs, nvdiff), intent(in) :: dvdftra
      real(kind=kind_phys), dimension(im), intent(in) :: dusfc1, dvsfc1, dtsfc1, dqsfc1, xmu
      real(kind=kind_phys), dimension(im, levs), intent(in) :: dudt, dvdt, dtdt, htrsw, htrlw

      real(kind=kind_phys), dimension(im, levs, ntrac), intent(inout) :: dqdt

      ! DH* The following arrays may not be allocated, depending on certain flags (cplflx, ...).
      ! Since Intel 15 crashes when passing unallocated arrays to arrays defined with explicit shape,
      ! use assumed-shape arrays. Note that Intel 18 and GNU 6.2.0-8.1.0 tolerate explicit-shape arrays
      ! as long as these do not get used when not allocated.
      real(kind=kind_phys), dimension(:,:), intent(inout) :: dt3dt, du3dt_PBL, du3dt_OGWD, dv3dt_PBL, dv3dt_OGWD, dq3dt, dq3dt_ozone
      real(kind=kind_phys), dimension(:), intent(inout) :: dusfc_cpl, dvsfc_cpl, dtsfc_cpl, dqsfc_cpl, dusfci_cpl, dvsfci_cpl, &
        dtsfci_cpl, dqsfci_cpl, dusfc_diag, dvsfc_diag, dtsfc_diag, dqsfc_diag, dusfci_diag, dvsfci_diag, dtsfci_diag, dqsfci_diag
      ! *DH

      character(len=*), intent(out) :: errmsg
      integer, intent(out) :: errflg

      integer :: i, k
      real(kind=kind_phys) :: tem

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0
!GJF: dvdftra is only used if nvdiff != ntrac or (nvdiff == ntrac .and. )
      if (nvdiff == ntrac .and. (hybedmf .or. do_shoc .or. satmedmf)) then
        dqdt = dvdftra
      elseif (nvdiff /= ntrac .and. .not. shinhong .and. .not. do_ysu) then
        if (imp_physics == imp_physics_wsm6) then
  ! WSM6
          do k=1,levs
            do i=1,im
              dqdt(i,k,ntqv)  = dvdftra(i,k,1)
              dqdt(i,k,ntcw)  = dvdftra(i,k,2)
              dqdt(i,k,ntiw)  = dvdftra(i,k,3)
              dqdt(i,k,ntoz)  = dvdftra(i,k,4)
            enddo
          enddo
        elseif (imp_physics == imp_physics_thompson) then
  ! Thompson
          ! DH* - Thompson ntrw, ntsw?
          if(ltaerosol) then
            do k=1,levs
              do i=1,im
                dqdt(i,k,ntqv)  = dvdftra(i,k,1)
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
                dqdt(i,k,ntqv)  = dvdftra(i,k,1)
                dqdt(i,k,ntcw)  = dvdftra(i,k,2)
                dqdt(i,k,ntiw)  = dvdftra(i,k,3)
                dqdt(i,k,ntinc) = dvdftra(i,k,4)
                dqdt(i,k,ntoz)  = dvdftra(i,k,5)
              enddo
            enddo
          endif
        elseif (imp_physics == imp_physics_gfdl) then
  ! GFDL MP
          do k=1,levs
            do i=1,im
              dqdt(i,k,ntqv) = dvdftra(i,k,1)
              dqdt(i,k,ntcw) = dvdftra(i,k,2)
              dqdt(i,k,ntiw) = dvdftra(i,k,3)
              dqdt(i,k,ntrw) = dvdftra(i,k,4)
              dqdt(i,k,ntsw) = dvdftra(i,k,5)
              dqdt(i,k,ntgl) = dvdftra(i,k,6)
              dqdt(i,k,ntoz) = dvdftra(i,k,7)
            enddo
          enddo
        elseif (imp_physics == imp_physics_zhao_carr) then
          if (cplchm) then
            do k=1,levs
              do i=1,im
                dqdt(i,k,1)    = dvdftra(i,k,1)
                dqdt(i,k,ntcw) = dvdftra(i,k,2)
                dqdt(i,k,ntoz) = dvdftra(i,k,3)
              enddo
            enddo
          endif
        endif

        if (satmedmf) then
          do k=1,levs
            do i=1,im
              dqdt(i,k,ntke)  = dvdftra(i,k,ntkev)
            enddo
          enddo
        endif

      endif ! nvdiff == ntrac

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
!          do k=1,levs
!            do i=1,im
!              tem  = dqdt(i,k,ntqv) * dtf
!              dq3dt(i,k) = dq3dt(i,k) + tem
!            enddo
!          enddo
!          if (ntoz > 0) then
!            do k=1,levs
!              do i=1,im
!                dq3dt_ozone(i,k) = dq3dt_ozone(i,k) + dqdt(i,k,ntoz) * dtf
!              enddo
!            enddo
!          endif
        endif

      endif   ! end if_lssav

      end subroutine GFS_PBL_generic_post_run

      end module GFS_PBL_generic_post

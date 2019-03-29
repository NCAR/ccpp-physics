
!
! This work (Common Community Physics Package), identified by NOAA, NCAR,
! CU/CIRES, is free of known copyright restrictions and is placed in the
! public domain.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
! THE AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
! IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
! CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
!

!>
!! @brief Auto-generated cap module for the m_micro scheme
!!
!
module m_micro_cap

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_f_pointer, c_ptr, c_int32_t
    use            :: ccpp_types,                                      &
                      only: ccpp_t, CCPP_GENERIC_KIND
    use            :: ccpp_fields,                                     &
                      only: ccpp_field_get
    use            :: ccpp_errors,                                     &
                      only: ccpp_error, ccpp_debug
    use            :: m_micro, &
                      only: m_micro_run,m_micro_finalize,m_micro_init
    ! Other modules required, e.g. type definitions
    use machine, only: kind_phys

    implicit none

    private
    public :: m_micro_run_cap,m_micro_finalize_cap,m_micro_init_cap

    contains


    function m_micro_run_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind
        integer, pointer :: im
        integer, pointer :: ix
        integer, pointer :: lm
        logical, pointer :: flipv
        real(kind_phys), pointer :: dt_i
        real(kind_phys), pointer :: prsl_i(:,:)
        real(kind_phys), pointer :: prsi_i(:,:)
        real(kind_phys), pointer :: phil(:,:)
        real(kind_phys), pointer :: phii(:,:)
        real(kind_phys), pointer :: omega_i(:,:)
        real(kind_phys), pointer :: qlls_i(:,:)
        real(kind_phys), pointer :: qlcn_i(:,:)
        real(kind_phys), pointer :: qils_i(:,:)
        real(kind_phys), pointer :: qicn_i(:,:)
        real(kind_phys), pointer :: lwheat_i(:,:)
        real(kind_phys), pointer :: swheat_i(:,:)
        real(kind_phys), pointer :: w_upi(:,:)
        real(kind_phys), pointer :: cf_upi(:,:)
        real(kind_phys), pointer :: frland(:)
        real(kind_phys), pointer :: zpbl(:)
        real(kind_phys), pointer :: cnv_mfd_i(:,:)
        real(kind_phys), pointer :: cnv_dqldt_i(:,:)
        real(kind_phys), pointer :: clcn_i(:,:)
        real(kind_phys), pointer :: u_i(:,:)
        real(kind_phys), pointer :: v_i(:,:)
        real(kind_phys), pointer :: taugwx(:)
        real(kind_phys), pointer :: taugwy(:)
        real(kind_phys), pointer :: tauorox(:)
        real(kind_phys), pointer :: tauoroy(:)
        real(kind_phys), pointer :: cnv_fice_i(:,:)
        real(kind_phys), pointer :: cnv_ndrop_i(:,:)
        real(kind_phys), pointer :: cnv_nice_i(:,:)
        real(kind_phys), pointer :: q_io(:,:)
        real(kind_phys), pointer :: lwm_o(:,:)
        real(kind_phys), pointer :: qi_o(:,:)
        real(kind_phys), pointer :: t_io(:,:)
        real(kind_phys), pointer :: rn_o(:)
        real(kind_phys), pointer :: sr_o(:)
        real(kind_phys), pointer :: ncpl_io(:,:)
        real(kind_phys), pointer :: ncpi_io(:,:)
        integer, pointer :: fprcp
        real(kind_phys), pointer :: rnw_io(:,:)
        real(kind_phys), pointer :: snw_io(:,:)
        real(kind_phys), pointer :: qgl_io(:,:)
        real(kind_phys), pointer :: ncpr_io(:,:)
        real(kind_phys), pointer :: ncps_io(:,:)
        real(kind_phys), pointer :: ncgl_io(:,:)
        real(kind_phys), pointer :: clls_io(:,:)
        integer, pointer :: kcbl(:)
        real(kind_phys), pointer :: cldreffl(:,:)
        real(kind_phys), pointer :: cldreffi(:,:)
        real(kind_phys), pointer :: cldreffr(:,:)
        real(kind_phys), pointer :: cldreffs(:,:)
        real(kind_phys), pointer :: cldreffg(:,:)
        real(kind_phys), pointer :: aerfld_i(:,:,:)
        logical, pointer :: aero_in
        real(kind_phys), pointer :: naai_i(:,:)
        real(kind_phys), pointer :: npccn_i(:,:)
        logical, pointer :: iccn
        logical, pointer :: skip_macro
        logical, pointer :: lprnt
        real(kind_phys), pointer :: alf_fac
        real(kind_phys), pointer :: qc_min(:)
        integer, pointer :: pdfflag
        integer, pointer :: ipr
        integer, pointer :: kdt
        real(kind_phys), pointer :: xlat(:)
        real(kind_phys), pointer :: xlon(:)
        real(kind_phys), pointer :: rhc_i(:,:)

        ierr = 0

        call c_f_pointer(ptr, cdata)


        call ccpp_field_get(cdata, 'horizontal_loop_extent', im, ierr=ierr, kind=ckind, index=366)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve horizontal_loop_extent from CCPP data structure')
            return
        end if
        if (kind(im).ne.ckind) then
            call ccpp_error('Kind mismatch for variable horizontal_loop_extent')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'horizontal_dimension', ix, ierr=ierr, kind=ckind, index=364)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve horizontal_dimension from CCPP data structure')
            return
        end if
        if (kind(ix).ne.ckind) then
            call ccpp_error('Kind mismatch for variable horizontal_dimension')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'vertical_dimension', lm, ierr=ierr, kind=ckind, index=817)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve vertical_dimension from CCPP data structure')
            return
        end if
        if (kind(lm).ne.ckind) then
            call ccpp_error('Kind mismatch for variable vertical_dimension')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'flag_flip', flipv, ierr=ierr, kind=ckind, index=265)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_flip from CCPP data structure')
            return
        end if
        if (kind(flipv).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_flip')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'time_step_for_physics', dt_i, ierr=ierr, kind=ckind, index=793)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve time_step_for_physics from CCPP data structure')
            return
        end if
        if (kind(dt_i).ne.ckind) then
            call ccpp_error('Kind mismatch for variable time_step_for_physics')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'air_pressure', prsl_i, ierr=ierr, dims=cdims, kind=ckind, index=44)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve air_pressure from CCPP data structure')
            return
        end if
        if (kind(prsl_i).ne.ckind) then
            call ccpp_error('Kind mismatch for variable air_pressure')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'air_pressure_at_interface', prsi_i, ierr=ierr, dims=cdims, kind=ckind, index=45)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve air_pressure_at_interface from CCPP data structure')
            return
        end if
        if (kind(prsi_i).ne.ckind) then
            call ccpp_error('Kind mismatch for variable air_pressure_at_interface')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'geopotential', phil, ierr=ierr, dims=cdims, kind=ckind, index=348)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve geopotential from CCPP data structure')
            return
        end if
        if (kind(phil).ne.ckind) then
            call ccpp_error('Kind mismatch for variable geopotential')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'geopotential_at_interface', phii, ierr=ierr, dims=cdims, kind=ckind, index=349)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve geopotential_at_interface from CCPP data structure')
            return
        end if
        if (kind(phii).ne.ckind) then
            call ccpp_error('Kind mismatch for variable geopotential_at_interface')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'omega', omega_i, ierr=ierr, dims=cdims, kind=ckind, index=593)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve omega from CCPP data structure')
            return
        end if
        if (kind(omega_i).ne.ckind) then
            call ccpp_error('Kind mismatch for variable omega')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cloud_condensed_water_mixing_ratio_convective_transport_tracer', qlls_i, ierr=ierr, dims=cdims, kind=ckind, index=93)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cloud_condensed_water_mixing_ratio_convective_transport_tracer from CCPP data structure')
            return
        end if
        if (kind(qlls_i).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cloud_condensed_water_mixing_ratio_convective_transport_tracer')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'mass_fraction_of_convective_cloud_liquid_water', qlcn_i, ierr=ierr, dims=cdims, kind=ckind, index=502)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve mass_fraction_of_convective_cloud_liquid_water from CCPP data structure')
            return
        end if
        if (kind(qlcn_i).ne.ckind) then
            call ccpp_error('Kind mismatch for variable mass_fraction_of_convective_cloud_liquid_water')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'ice_water_mixing_ratio_convective_transport_tracer', qils_i, ierr=ierr, dims=cdims, kind=ckind, index=374)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve ice_water_mixing_ratio_convective_transport_tracer from CCPP data structure')
            return
        end if
        if (kind(qils_i).ne.ckind) then
            call ccpp_error('Kind mismatch for variable ice_water_mixing_ratio_convective_transport_tracer')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'mass_fraction_of_convective_cloud_ice', qicn_i, ierr=ierr, dims=cdims, kind=ckind, index=501)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve mass_fraction_of_convective_cloud_ice from CCPP data structure')
            return
        end if
        if (kind(qicn_i).ne.ckind) then
            call ccpp_error('Kind mismatch for variable mass_fraction_of_convective_cloud_ice')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'tendency_of_air_temperature_due_to_longwave_heating_on_radiation_timestep', lwheat_i, ierr=ierr, dims=cdims, kind=ckind, index=757)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve tendency_of_air_temperature_due_to_longwave_heating_on_radiation_timestep from CCPP data structure')
            return
        end if
        if (kind(lwheat_i).ne.ckind) then
            call ccpp_error('Kind mismatch for variable tendency_of_air_temperature_due_to_longwave_heating_on_radiation_timestep')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'tendency_of_air_temperature_due_to_shortwave_heating_on_radiation_timestep', swheat_i, ierr=ierr, dims=cdims, kind=ckind, index=763)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve tendency_of_air_temperature_due_to_shortwave_heating_on_radiation_timestep from CCPP data structure')
            return
        end if
        if (kind(swheat_i).ne.ckind) then
            call ccpp_error('Kind mismatch for variable tendency_of_air_temperature_due_to_shortwave_heating_on_radiation_timestep')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'vertical_velocity_for_updraft', w_upi, ierr=ierr, dims=cdims, kind=ckind, index=830)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve vertical_velocity_for_updraft from CCPP data structure')
            return
        end if
        if (kind(w_upi).ne.ckind) then
            call ccpp_error('Kind mismatch for variable vertical_velocity_for_updraft')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'convective_cloud_fraction_for_microphysics', cf_upi, ierr=ierr, dims=cdims, kind=ckind, index=125)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve convective_cloud_fraction_for_microphysics from CCPP data structure')
            return
        end if
        if (kind(cf_upi).ne.ckind) then
            call ccpp_error('Kind mismatch for variable convective_cloud_fraction_for_microphysics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'land_area_fraction', frland, ierr=ierr, dims=cdims, kind=ckind, index=456)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve land_area_fraction from CCPP data structure')
            return
        end if
        if (kind(frland).ne.ckind) then
            call ccpp_error('Kind mismatch for variable land_area_fraction')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'atmosphere_boundary_layer_thickness', zpbl, ierr=ierr, dims=cdims, kind=ckind, index=66)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve atmosphere_boundary_layer_thickness from CCPP data structure')
            return
        end if
        if (kind(zpbl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable atmosphere_boundary_layer_thickness')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'detrained_mass_flux', cnv_mfd_i, ierr=ierr, dims=cdims, kind=ckind, index=213)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve detrained_mass_flux from CCPP data structure')
            return
        end if
        if (kind(cnv_mfd_i).ne.ckind) then
            call ccpp_error('Kind mismatch for variable detrained_mass_flux')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'tendency_of_cloud_water_due_to_convective_microphysics', cnv_dqldt_i, ierr=ierr, dims=cdims, kind=ckind, index=766)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve tendency_of_cloud_water_due_to_convective_microphysics from CCPP data structure')
            return
        end if
        if (kind(cnv_dqldt_i).ne.ckind) then
            call ccpp_error('Kind mismatch for variable tendency_of_cloud_water_due_to_convective_microphysics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'convective_cloud_volume_fraction', clcn_i, ierr=ierr, dims=cdims, kind=ckind, index=127)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve convective_cloud_volume_fraction from CCPP data structure')
            return
        end if
        if (kind(clcn_i).ne.ckind) then
            call ccpp_error('Kind mismatch for variable convective_cloud_volume_fraction')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'x_wind_updated_by_physics', u_i, ierr=ierr, dims=cdims, kind=ckind, index=877)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve x_wind_updated_by_physics from CCPP data structure')
            return
        end if
        if (kind(u_i).ne.ckind) then
            call ccpp_error('Kind mismatch for variable x_wind_updated_by_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'y_wind_updated_by_physics', v_i, ierr=ierr, dims=cdims, kind=ckind, index=884)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve y_wind_updated_by_physics from CCPP data structure')
            return
        end if
        if (kind(v_i).ne.ckind) then
            call ccpp_error('Kind mismatch for variable y_wind_updated_by_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cumulative_surface_x_momentum_flux_for_diag_multiplied_by_timestep', taugwx, ierr=ierr, dims=cdims, kind=ckind, index=202)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cumulative_surface_x_momentum_flux_for_diag_multiplied_by_timestep from CCPP data structure')
            return
        end if
        if (kind(taugwx).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cumulative_surface_x_momentum_flux_for_diag_multiplied_by_timestep')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cumulative_surface_y_momentum_flux_for_diag_multiplied_by_timestep', taugwy, ierr=ierr, dims=cdims, kind=ckind, index=204)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cumulative_surface_y_momentum_flux_for_diag_multiplied_by_timestep from CCPP data structure')
            return
        end if
        if (kind(taugwy).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cumulative_surface_y_momentum_flux_for_diag_multiplied_by_timestep')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'instantaneous_surface_x_momentum_flux', tauorox, ierr=ierr, dims=cdims, kind=ckind, index=437)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve instantaneous_surface_x_momentum_flux from CCPP data structure')
            return
        end if
        if (kind(tauorox).ne.ckind) then
            call ccpp_error('Kind mismatch for variable instantaneous_surface_x_momentum_flux')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'instantaneous_surface_y_momentum_flux', tauoroy, ierr=ierr, dims=cdims, kind=ckind, index=440)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve instantaneous_surface_y_momentum_flux from CCPP data structure')
            return
        end if
        if (kind(tauoroy).ne.ckind) then
            call ccpp_error('Kind mismatch for variable instantaneous_surface_y_momentum_flux')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'ice_fraction_in_convective_tower', cnv_fice_i, ierr=ierr, dims=cdims, kind=ckind, index=367)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve ice_fraction_in_convective_tower from CCPP data structure')
            return
        end if
        if (kind(cnv_fice_i).ne.ckind) then
            call ccpp_error('Kind mismatch for variable ice_fraction_in_convective_tower')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'number_concentration_of_cloud_liquid_water_particles_for_detrainment', cnv_ndrop_i, ierr=ierr, dims=cdims, kind=ckind, index=569)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve number_concentration_of_cloud_liquid_water_particles_for_detrainment from CCPP data structure')
            return
        end if
        if (kind(cnv_ndrop_i).ne.ckind) then
            call ccpp_error('Kind mismatch for variable number_concentration_of_cloud_liquid_water_particles_for_detrainment')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'number_concentration_of_ice_crystals_for_detrainment', cnv_nice_i, ierr=ierr, dims=cdims, kind=ckind, index=570)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve number_concentration_of_ice_crystals_for_detrainment from CCPP data structure')
            return
        end if
        if (kind(cnv_nice_i).ne.ckind) then
            call ccpp_error('Kind mismatch for variable number_concentration_of_ice_crystals_for_detrainment')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'water_vapor_specific_humidity_updated_by_physics', q_io, ierr=ierr, dims=cdims, kind=ckind, index=860)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve water_vapor_specific_humidity_updated_by_physics from CCPP data structure')
            return
        end if
        if (kind(q_io).ne.ckind) then
            call ccpp_error('Kind mismatch for variable water_vapor_specific_humidity_updated_by_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cloud_condensed_water_mixing_ratio_updated_by_physics', lwm_o, ierr=ierr, dims=cdims, kind=ckind, index=95)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cloud_condensed_water_mixing_ratio_updated_by_physics from CCPP data structure')
            return
        end if
        if (kind(lwm_o).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cloud_condensed_water_mixing_ratio_updated_by_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'ice_water_mixing_ratio_updated_by_physics', qi_o, ierr=ierr, dims=cdims, kind=ckind, index=376)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve ice_water_mixing_ratio_updated_by_physics from CCPP data structure')
            return
        end if
        if (kind(qi_o).ne.ckind) then
            call ccpp_error('Kind mismatch for variable ice_water_mixing_ratio_updated_by_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'air_temperature_updated_by_physics', t_io, ierr=ierr, dims=cdims, kind=ckind, index=59)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve air_temperature_updated_by_physics from CCPP data structure')
            return
        end if
        if (kind(t_io).ne.ckind) then
            call ccpp_error('Kind mismatch for variable air_temperature_updated_by_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'lwe_thickness_of_explicit_precipitation_amount', rn_o, ierr=ierr, dims=cdims, kind=ckind, index=480)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve lwe_thickness_of_explicit_precipitation_amount from CCPP data structure')
            return
        end if
        if (kind(rn_o).ne.ckind) then
            call ccpp_error('Kind mismatch for variable lwe_thickness_of_explicit_precipitation_amount')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'ratio_of_snowfall_to_rainfall', sr_o, ierr=ierr, dims=cdims, kind=ckind, index=624)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve ratio_of_snowfall_to_rainfall from CCPP data structure')
            return
        end if
        if (kind(sr_o).ne.ckind) then
            call ccpp_error('Kind mismatch for variable ratio_of_snowfall_to_rainfall')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cloud_droplet_number_concentration_updated_by_physics', ncpl_io, ierr=ierr, dims=cdims, kind=ckind, index=98)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cloud_droplet_number_concentration_updated_by_physics from CCPP data structure')
            return
        end if
        if (kind(ncpl_io).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cloud_droplet_number_concentration_updated_by_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'ice_number_concentration_updated_by_physics', ncpi_io, ierr=ierr, dims=cdims, kind=ckind, index=371)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve ice_number_concentration_updated_by_physics from CCPP data structure')
            return
        end if
        if (kind(ncpi_io).ne.ckind) then
            call ccpp_error('Kind mismatch for variable ice_number_concentration_updated_by_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'number_of_frozen_precipitation_species', fprcp, ierr=ierr, kind=ckind, index=577)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve number_of_frozen_precipitation_species from CCPP data structure')
            return
        end if
        if (kind(fprcp).ne.ckind) then
            call ccpp_error('Kind mismatch for variable number_of_frozen_precipitation_species')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'local_rain_water_mixing_ratio', rnw_io, ierr=ierr, dims=cdims, kind=ckind, index=470)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve local_rain_water_mixing_ratio from CCPP data structure')
            return
        end if
        if (kind(rnw_io).ne.ckind) then
            call ccpp_error('Kind mismatch for variable local_rain_water_mixing_ratio')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'local_snow_water_mixing_ratio', snw_io, ierr=ierr, dims=cdims, kind=ckind, index=472)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve local_snow_water_mixing_ratio from CCPP data structure')
            return
        end if
        if (kind(snw_io).ne.ckind) then
            call ccpp_error('Kind mismatch for variable local_snow_water_mixing_ratio')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'local_graupel_mixing_ratio', qgl_io, ierr=ierr, dims=cdims, kind=ckind, index=467)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve local_graupel_mixing_ratio from CCPP data structure')
            return
        end if
        if (kind(qgl_io).ne.ckind) then
            call ccpp_error('Kind mismatch for variable local_graupel_mixing_ratio')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'local_rain_number_concentration', ncpr_io, ierr=ierr, dims=cdims, kind=ckind, index=469)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve local_rain_number_concentration from CCPP data structure')
            return
        end if
        if (kind(ncpr_io).ne.ckind) then
            call ccpp_error('Kind mismatch for variable local_rain_number_concentration')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'local_snow_number_concentration', ncps_io, ierr=ierr, dims=cdims, kind=ckind, index=471)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve local_snow_number_concentration from CCPP data structure')
            return
        end if
        if (kind(ncps_io).ne.ckind) then
            call ccpp_error('Kind mismatch for variable local_snow_number_concentration')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'local_graupel_number_concentration', ncgl_io, ierr=ierr, dims=cdims, kind=ckind, index=468)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve local_graupel_number_concentration from CCPP data structure')
            return
        end if
        if (kind(ncgl_io).ne.ckind) then
            call ccpp_error('Kind mismatch for variable local_graupel_number_concentration')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cloud_fraction_for_MG', clls_io, ierr=ierr, dims=cdims, kind=ckind, index=99)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cloud_fraction_for_MG from CCPP data structure')
            return
        end if
        if (kind(clls_io).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cloud_fraction_for_MG')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'vertical_index_at_cloud_base', kcbl, ierr=ierr, dims=cdims, kind=ckind, index=820)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve vertical_index_at_cloud_base from CCPP data structure')
            return
        end if
        if (kind(kcbl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable vertical_index_at_cloud_base')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'effective_radius_of_stratiform_cloud_liquid_water_particle_in_um', cldreffl, ierr=ierr, dims=cdims, kind=ckind, index=244)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve effective_radius_of_stratiform_cloud_liquid_water_particle_in_um from CCPP data structure')
            return
        end if
        if (kind(cldreffl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable effective_radius_of_stratiform_cloud_liquid_water_particle_in_um')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'effective_radius_of_stratiform_cloud_ice_particle_in_um', cldreffi, ierr=ierr, dims=cdims, kind=ckind, index=243)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve effective_radius_of_stratiform_cloud_ice_particle_in_um from CCPP data structure')
            return
        end if
        if (kind(cldreffi).ne.ckind) then
            call ccpp_error('Kind mismatch for variable effective_radius_of_stratiform_cloud_ice_particle_in_um')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'effective_radius_of_stratiform_cloud_rain_particle_in_um', cldreffr, ierr=ierr, dims=cdims, kind=ckind, index=245)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve effective_radius_of_stratiform_cloud_rain_particle_in_um from CCPP data structure')
            return
        end if
        if (kind(cldreffr).ne.ckind) then
            call ccpp_error('Kind mismatch for variable effective_radius_of_stratiform_cloud_rain_particle_in_um')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'effective_radius_of_stratiform_cloud_snow_particle_in_um', cldreffs, ierr=ierr, dims=cdims, kind=ckind, index=246)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve effective_radius_of_stratiform_cloud_snow_particle_in_um from CCPP data structure')
            return
        end if
        if (kind(cldreffs).ne.ckind) then
            call ccpp_error('Kind mismatch for variable effective_radius_of_stratiform_cloud_snow_particle_in_um')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'effective_radius_of_stratiform_cloud_graupel_particle_in_um', cldreffg, ierr=ierr, dims=cdims, kind=ckind, index=242)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve effective_radius_of_stratiform_cloud_graupel_particle_in_um from CCPP data structure')
            return
        end if
        if (kind(cldreffg).ne.ckind) then
            call ccpp_error('Kind mismatch for variable effective_radius_of_stratiform_cloud_graupel_particle_in_um')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'aerosol_number_concentration_from_gocart_aerosol_climatology', aerfld_i, ierr=ierr, dims=cdims, kind=ckind, index=37)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve aerosol_number_concentration_from_gocart_aerosol_climatology from CCPP data structure')
            return
        end if
        if (kind(aerfld_i).ne.ckind) then
            call ccpp_error('Kind mismatch for variable aerosol_number_concentration_from_gocart_aerosol_climatology')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'flag_for_aerosol_input_MG', aero_in, ierr=ierr, kind=ckind, index=270)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_aerosol_input_MG from CCPP data structure')
            return
        end if
        if (kind(aero_in).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_aerosol_input_MG')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'in_number_concentration', naai_i, ierr=ierr, dims=cdims, kind=ckind, index=377)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve in_number_concentration from CCPP data structure')
            return
        end if
        if (kind(naai_i).ne.ckind) then
            call ccpp_error('Kind mismatch for variable in_number_concentration')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'ccn_number_concentration', npccn_i, ierr=ierr, dims=cdims, kind=ckind, index=81)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve ccn_number_concentration from CCPP data structure')
            return
        end if
        if (kind(npccn_i).ne.ckind) then
            call ccpp_error('Kind mismatch for variable ccn_number_concentration')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'flag_for_in_ccn_forcing_for_morrison_gettelman_microphysics', iccn, ierr=ierr, kind=ckind, index=285)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_in_ccn_forcing_for_morrison_gettelman_microphysics from CCPP data structure')
            return
        end if
        if (kind(iccn).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_in_ccn_forcing_for_morrison_gettelman_microphysics')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'flag_skip_macro', skip_macro, ierr=ierr, kind=ckind, index=334)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_skip_macro from CCPP data structure')
            return
        end if
        if (kind(skip_macro).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_skip_macro')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'flag_print', lprnt, ierr=ierr, kind=ckind, index=332)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_print from CCPP data structure')
            return
        end if
        if (kind(lprnt).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_print')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'mg_tuning_factor_for_alphas', alf_fac, ierr=ierr, kind=ckind, index=542)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve mg_tuning_factor_for_alphas from CCPP data structure')
            return
        end if
        if (kind(alf_fac).ne.ckind) then
            call ccpp_error('Kind mismatch for variable mg_tuning_factor_for_alphas')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'mg_minimum_cloud_condensed_water_and_ice_mixing_ratio', qc_min, ierr=ierr, dims=cdims, kind=ckind, index=539)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve mg_minimum_cloud_condensed_water_and_ice_mixing_ratio from CCPP data structure')
            return
        end if
        if (kind(qc_min).ne.ckind) then
            call ccpp_error('Kind mismatch for variable mg_minimum_cloud_condensed_water_and_ice_mixing_ratio')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'flag_for_pdf_for_morrison_gettelman_microphysics_scheme', pdfflag, ierr=ierr, kind=ckind, index=302)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_pdf_for_morrison_gettelman_microphysics_scheme from CCPP data structure')
            return
        end if
        if (kind(pdfflag).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_pdf_for_morrison_gettelman_microphysics_scheme')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'horizontal_index_of_printed_column', ipr, ierr=ierr, kind=ckind, index=365)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve horizontal_index_of_printed_column from CCPP data structure')
            return
        end if
        if (kind(ipr).ne.ckind) then
            call ccpp_error('Kind mismatch for variable horizontal_index_of_printed_column')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'index_of_time_step', kdt, ierr=ierr, kind=ckind, index=398)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve index_of_time_step from CCPP data structure')
            return
        end if
        if (kind(kdt).ne.ckind) then
            call ccpp_error('Kind mismatch for variable index_of_time_step')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'latitude', xlat, ierr=ierr, dims=cdims, kind=ckind, index=460)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve latitude from CCPP data structure')
            return
        end if
        if (kind(xlat).ne.ckind) then
            call ccpp_error('Kind mismatch for variable latitude')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'longitude', xlon, ierr=ierr, dims=cdims, kind=ckind, index=473)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve longitude from CCPP data structure')
            return
        end if
        if (kind(xlon).ne.ckind) then
            call ccpp_error('Kind mismatch for variable longitude')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'critical_relative_humidity', rhc_i, ierr=ierr, dims=cdims, kind=ckind, index=141)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve critical_relative_humidity from CCPP data structure')
            return
        end if
        if (kind(rhc_i).ne.ckind) then
            call ccpp_error('Kind mismatch for variable critical_relative_humidity')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call m_micro_run(im=im,ix=ix,lm=lm,flipv=flipv,dt_i=dt_i,prsl_i=prsl_i,prsi_i=prsi_i,phil=phil, &
                  phii=phii,omega_i=omega_i,qlls_i=qlls_i,qlcn_i=qlcn_i,qils_i=qils_i,qicn_i=qicn_i, &
                  lwheat_i=lwheat_i,swheat_i=swheat_i,w_upi=w_upi,cf_upi=cf_upi,frland=frland, &
                  zpbl=zpbl,cnv_mfd_i=cnv_mfd_i,cnv_dqldt_i=cnv_dqldt_i,clcn_i=clcn_i,u_i=u_i, &
                  v_i=v_i,taugwx=taugwx,taugwy=taugwy,tauorox=tauorox,tauoroy=tauoroy,cnv_fice_i=cnv_fice_i, &
                  cnv_ndrop_i=cnv_ndrop_i,cnv_nice_i=cnv_nice_i,q_io=q_io,lwm_o=lwm_o,qi_o=qi_o, &
                  t_io=t_io,rn_o=rn_o,sr_o=sr_o,ncpl_io=ncpl_io,ncpi_io=ncpi_io,fprcp=fprcp, &
                  rnw_io=rnw_io,snw_io=snw_io,qgl_io=qgl_io,ncpr_io=ncpr_io,ncps_io=ncps_io, &
                  ncgl_io=ncgl_io,clls_io=clls_io,kcbl=kcbl,cldreffl=cldreffl,cldreffi=cldreffi, &
                  cldreffr=cldreffr,cldreffs=cldreffs,cldreffg=cldreffg,aerfld_i=aerfld_i, &
                  aero_in=aero_in,naai_i=naai_i,npccn_i=npccn_i,iccn=iccn,skip_macro=skip_macro, &
                  lprnt=lprnt,alf_fac=alf_fac,qc_min=qc_min,pdfflag=pdfflag,ipr=ipr,kdt=kdt, &
                  xlat=xlat,xlon=xlon,rhc_i=rhc_i,errmsg=cdata%errmsg,errflg=cdata%errflg)
        ierr=cdata%errflg

    end function m_micro_run_cap

    function m_micro_finalize_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call m_micro_finalize()
        

    end function m_micro_finalize_cap

    function m_micro_init_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind
        integer, pointer :: imp_physics
        integer, pointer :: imp_physics_mg
        integer, pointer :: fprcp
        real(kind_phys), pointer :: gravit
        real(kind_phys), pointer :: rair
        real(kind_phys), pointer :: rh2o
        real(kind_phys), pointer :: cpair
        real(kind_phys), pointer :: tmelt
        real(kind_phys), pointer :: latvap
        real(kind_phys), pointer :: latice
        real(kind_phys), pointer :: mg_dcs
        real(kind_phys), pointer :: mg_qcvar
        real(kind_phys), pointer :: mg_ts_auto_ice(:)
        real(kind_phys), pointer :: mg_rhmini
        logical, pointer :: microp_uniform
        logical, pointer :: do_cldice
        logical, pointer :: hetfrz_classnuc
        character(len=16), pointer :: mg_precip_frac_method
        real(kind_phys), pointer :: mg_berg_eff_factor
        logical, pointer :: sed_supersat
        logical, pointer :: do_sb_physics
        logical, pointer :: mg_do_hail
        logical, pointer :: mg_do_graupel
        logical, pointer :: mg_nccons
        logical, pointer :: mg_nicons
        logical, pointer :: mg_ngcons
        real(kind_phys), pointer :: mg_ncnst
        real(kind_phys), pointer :: mg_ninst
        real(kind_phys), pointer :: mg_ngnst
        logical, pointer :: mg_do_ice_gmao
        logical, pointer :: mg_do_liq_liu

        ierr = 0

        call c_f_pointer(ptr, cdata)


        call ccpp_field_get(cdata, 'flag_for_microphysics_scheme', imp_physics, ierr=ierr, kind=ckind, index=294)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_microphysics_scheme from CCPP data structure')
            return
        end if
        if (kind(imp_physics).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_microphysics_scheme')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'flag_for_morrison_gettelman_microphysics_scheme', imp_physics_mg, ierr=ierr, kind=ckind, index=297)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_morrison_gettelman_microphysics_scheme from CCPP data structure')
            return
        end if
        if (kind(imp_physics_mg).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_morrison_gettelman_microphysics_scheme')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'number_of_frozen_precipitation_species', fprcp, ierr=ierr, kind=ckind, index=577)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve number_of_frozen_precipitation_species from CCPP data structure')
            return
        end if
        if (kind(fprcp).ne.ckind) then
            call ccpp_error('Kind mismatch for variable number_of_frozen_precipitation_species')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'gravitational_acceleration', gravit, ierr=ierr, kind=ckind, index=355)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve gravitational_acceleration from CCPP data structure')
            return
        end if
        if (kind(gravit).ne.ckind) then
            call ccpp_error('Kind mismatch for variable gravitational_acceleration')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'gas_constant_dry_air', rair, ierr=ierr, kind=ckind, index=346)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve gas_constant_dry_air from CCPP data structure')
            return
        end if
        if (kind(rair).ne.ckind) then
            call ccpp_error('Kind mismatch for variable gas_constant_dry_air')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'gas_constant_water_vapor', rh2o, ierr=ierr, kind=ckind, index=347)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve gas_constant_water_vapor from CCPP data structure')
            return
        end if
        if (kind(rh2o).ne.ckind) then
            call ccpp_error('Kind mismatch for variable gas_constant_water_vapor')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'specific_heat_of_dry_air_at_constant_pressure', cpair, ierr=ierr, kind=ckind, index=664)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve specific_heat_of_dry_air_at_constant_pressure from CCPP data structure')
            return
        end if
        if (kind(cpair).ne.ckind) then
            call ccpp_error('Kind mismatch for variable specific_heat_of_dry_air_at_constant_pressure')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'triple_point_temperature_of_water', tmelt, ierr=ierr, kind=ckind, index=805)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve triple_point_temperature_of_water from CCPP data structure')
            return
        end if
        if (kind(tmelt).ne.ckind) then
            call ccpp_error('Kind mismatch for variable triple_point_temperature_of_water')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'latent_heat_of_vaporization_of_water_at_0C', latvap, ierr=ierr, kind=ckind, index=459)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve latent_heat_of_vaporization_of_water_at_0C from CCPP data structure')
            return
        end if
        if (kind(latvap).ne.ckind) then
            call ccpp_error('Kind mismatch for variable latent_heat_of_vaporization_of_water_at_0C')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'latent_heat_of_fusion_of_water_at_0C', latice, ierr=ierr, kind=ckind, index=458)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve latent_heat_of_fusion_of_water_at_0C from CCPP data structure')
            return
        end if
        if (kind(latice).ne.ckind) then
            call ccpp_error('Kind mismatch for variable latent_heat_of_fusion_of_water_at_0C')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'mg_autoconversion_size_threshold_ice_snow', mg_dcs, ierr=ierr, kind=ckind, index=522)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve mg_autoconversion_size_threshold_ice_snow from CCPP data structure')
            return
        end if
        if (kind(mg_dcs).ne.ckind) then
            call ccpp_error('Kind mismatch for variable mg_autoconversion_size_threshold_ice_snow')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'mg_cloud_water_variance', mg_qcvar, ierr=ierr, kind=ckind, index=524)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve mg_cloud_water_variance from CCPP data structure')
            return
        end if
        if (kind(mg_qcvar).ne.ckind) then
            call ccpp_error('Kind mismatch for variable mg_cloud_water_variance')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'mg_time_scale_for_autoconversion_of_ice', mg_ts_auto_ice, ierr=ierr, dims=cdims, kind=ckind, index=541)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve mg_time_scale_for_autoconversion_of_ice from CCPP data structure')
            return
        end if
        if (kind(mg_ts_auto_ice).ne.ckind) then
            call ccpp_error('Kind mismatch for variable mg_time_scale_for_autoconversion_of_ice')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'mg_minimum_rh_for_ice', mg_rhmini, ierr=ierr, kind=ckind, index=540)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve mg_minimum_rh_for_ice from CCPP data structure')
            return
        end if
        if (kind(mg_rhmini).ne.ckind) then
            call ccpp_error('Kind mismatch for variable mg_minimum_rh_for_ice')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'mg_flag_for_uniform_subcolumns', microp_uniform, ierr=ierr, kind=ckind, index=534)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve mg_flag_for_uniform_subcolumns from CCPP data structure')
            return
        end if
        if (kind(microp_uniform).ne.ckind) then
            call ccpp_error('Kind mismatch for variable mg_flag_for_uniform_subcolumns')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'mg_flag_for_cloud_ice_processes', do_cldice, ierr=ierr, kind=ckind, index=527)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve mg_flag_for_cloud_ice_processes from CCPP data structure')
            return
        end if
        if (kind(do_cldice).ne.ckind) then
            call ccpp_error('Kind mismatch for variable mg_flag_for_cloud_ice_processes')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'mg_flag_for_heterogeneous_freezing', hetfrz_classnuc, ierr=ierr, kind=ckind, index=531)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve mg_flag_for_heterogeneous_freezing from CCPP data structure')
            return
        end if
        if (kind(hetfrz_classnuc).ne.ckind) then
            call ccpp_error('Kind mismatch for variable mg_flag_for_heterogeneous_freezing')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'mg_type_of_precip_fraction_method', mg_precip_frac_method, ierr=ierr, kind=ckind, index=543)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve mg_type_of_precip_fraction_method from CCPP data structure')
            return
        end if
        if (kind(mg_precip_frac_method).ne.ckind) then
            call ccpp_error('Kind mismatch for variable mg_type_of_precip_fraction_method')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'mg_bergeron_efficiency_factor', mg_berg_eff_factor, ierr=ierr, kind=ckind, index=523)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve mg_bergeron_efficiency_factor from CCPP data structure')
            return
        end if
        if (kind(mg_berg_eff_factor).ne.ckind) then
            call ccpp_error('Kind mismatch for variable mg_bergeron_efficiency_factor')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'mg_allow_supersat_after_sed', sed_supersat, ierr=ierr, kind=ckind, index=521)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve mg_allow_supersat_after_sed from CCPP data structure')
            return
        end if
        if (kind(sed_supersat).ne.ckind) then
            call ccpp_error('Kind mismatch for variable mg_allow_supersat_after_sed')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'mg_flag_for_sb2001_autoconversion', do_sb_physics, ierr=ierr, kind=ckind, index=533)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve mg_flag_for_sb2001_autoconversion from CCPP data structure')
            return
        end if
        if (kind(do_sb_physics).ne.ckind) then
            call ccpp_error('Kind mismatch for variable mg_flag_for_sb2001_autoconversion')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'mg_flag_for_hail', mg_do_hail, ierr=ierr, kind=ckind, index=530)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve mg_flag_for_hail from CCPP data structure')
            return
        end if
        if (kind(mg_do_hail).ne.ckind) then
            call ccpp_error('Kind mismatch for variable mg_flag_for_hail')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'mg_flag_for_graupel', mg_do_graupel, ierr=ierr, kind=ckind, index=529)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve mg_flag_for_graupel from CCPP data structure')
            return
        end if
        if (kind(mg_do_graupel).ne.ckind) then
            call ccpp_error('Kind mismatch for variable mg_flag_for_graupel')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'mg_flag_drop_concentration_constant', mg_nccons, ierr=ierr, kind=ckind, index=526)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve mg_flag_drop_concentration_constant from CCPP data structure')
            return
        end if
        if (kind(mg_nccons).ne.ckind) then
            call ccpp_error('Kind mismatch for variable mg_flag_drop_concentration_constant')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'mg_flag_ice_concentration_constant', mg_nicons, ierr=ierr, kind=ckind, index=536)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve mg_flag_ice_concentration_constant from CCPP data structure')
            return
        end if
        if (kind(mg_nicons).ne.ckind) then
            call ccpp_error('Kind mismatch for variable mg_flag_ice_concentration_constant')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'mg_flag_graupel_concentration_constant', mg_ngcons, ierr=ierr, kind=ckind, index=535)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve mg_flag_graupel_concentration_constant from CCPP data structure')
            return
        end if
        if (kind(mg_ngcons).ne.ckind) then
            call ccpp_error('Kind mismatch for variable mg_flag_graupel_concentration_constant')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'mg_drop_concentration_constant', mg_ncnst, ierr=ierr, kind=ckind, index=525)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve mg_drop_concentration_constant from CCPP data structure')
            return
        end if
        if (kind(mg_ncnst).ne.ckind) then
            call ccpp_error('Kind mismatch for variable mg_drop_concentration_constant')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'mg_ice_concentration_constant', mg_ninst, ierr=ierr, kind=ckind, index=538)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve mg_ice_concentration_constant from CCPP data structure')
            return
        end if
        if (kind(mg_ninst).ne.ckind) then
            call ccpp_error('Kind mismatch for variable mg_ice_concentration_constant')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'mg_graupel_concentration_constant', mg_ngnst, ierr=ierr, kind=ckind, index=537)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve mg_graupel_concentration_constant from CCPP data structure')
            return
        end if
        if (kind(mg_ngnst).ne.ckind) then
            call ccpp_error('Kind mismatch for variable mg_graupel_concentration_constant')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'mg_flag_for_gmao_ice_formulation', mg_do_ice_gmao, ierr=ierr, kind=ckind, index=528)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve mg_flag_for_gmao_ice_formulation from CCPP data structure')
            return
        end if
        if (kind(mg_do_ice_gmao).ne.ckind) then
            call ccpp_error('Kind mismatch for variable mg_flag_for_gmao_ice_formulation')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'mg_flag_for_liu_liquid_treatment', mg_do_liq_liu, ierr=ierr, kind=ckind, index=532)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve mg_flag_for_liu_liquid_treatment from CCPP data structure')
            return
        end if
        if (kind(mg_do_liq_liu).ne.ckind) then
            call ccpp_error('Kind mismatch for variable mg_flag_for_liu_liquid_treatment')
            ierr = 1
            return
        end if
#endif
        

        call m_micro_init(imp_physics=imp_physics,imp_physics_mg=imp_physics_mg,fprcp=fprcp,gravit=gravit, &
                  rair=rair,rh2o=rh2o,cpair=cpair,tmelt=tmelt,latvap=latvap,latice=latice, &
                  mg_dcs=mg_dcs,mg_qcvar=mg_qcvar,mg_ts_auto_ice=mg_ts_auto_ice,mg_rhmini=mg_rhmini, &
                  microp_uniform=microp_uniform,do_cldice=do_cldice,hetfrz_classnuc=hetfrz_classnuc, &
                  mg_precip_frac_method=mg_precip_frac_method,mg_berg_eff_factor=mg_berg_eff_factor, &
                  sed_supersat=sed_supersat,do_sb_physics=do_sb_physics,mg_do_hail=mg_do_hail, &
                  mg_do_graupel=mg_do_graupel,mg_nccons=mg_nccons,mg_nicons=mg_nicons,mg_ngcons=mg_ngcons, &
                  mg_ncnst=mg_ncnst,mg_ninst=mg_ninst,mg_ngnst=mg_ngnst,mg_do_ice_gmao=mg_do_ice_gmao, &
                  mg_do_liq_liu=mg_do_liq_liu,errmsg=cdata%errmsg,errflg=cdata%errflg)
        ierr=cdata%errflg

    end function m_micro_init_cap
end module m_micro_cap

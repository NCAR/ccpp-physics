
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
!! @brief Auto-generated cap module for the moninshoc scheme
!!
!
module moninshoc_cap

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_f_pointer, c_ptr, c_int32_t
    use            :: ccpp_types,                                      &
                      only: ccpp_t, CCPP_GENERIC_KIND
    use            :: ccpp_fields,                                     &
                      only: ccpp_field_get
    use            :: ccpp_errors,                                     &
                      only: ccpp_error, ccpp_debug
    use            :: moninshoc, &
                      only: moninshoc_run,moninshoc_init,moninshoc_finalize
    ! Other modules required, e.g. type definitions
    use machine, only: kind_phys

    implicit none

    private
    public :: moninshoc_run_cap,moninshoc_init_cap,moninshoc_finalize_cap

    contains


    function moninshoc_run_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind
        integer, pointer :: ix
        integer, pointer :: im
        integer, pointer :: km
        integer, pointer :: ntrac
        integer, pointer :: ntcw
        integer, pointer :: ncnd
        real(kind_phys), pointer :: dv(:,:)
        real(kind_phys), pointer :: du(:,:)
        real(kind_phys), pointer :: tau(:,:)
        real(kind_phys), pointer :: rtg(:,:,:)
        real(kind_phys), pointer :: u1(:,:)
        real(kind_phys), pointer :: v1(:,:)
        real(kind_phys), pointer :: t1(:,:)
        real(kind_phys), pointer :: q1(:,:,:)
        real(kind_phys), pointer :: tkh(:,:)
        real(kind_phys), pointer :: prnum(:,:)
        integer, pointer :: ntke
        real(kind_phys), pointer :: psk(:)
        real(kind_phys), pointer :: rbsoil(:)
        real(kind_phys), pointer :: zorl(:)
        real(kind_phys), pointer :: u10m(:)
        real(kind_phys), pointer :: v10m(:)
        real(kind_phys), pointer :: fm(:)
        real(kind_phys), pointer :: fh(:)
        real(kind_phys), pointer :: tsea(:)
        real(kind_phys), pointer :: heat(:)
        real(kind_phys), pointer :: evap(:)
        real(kind_phys), pointer :: stress(:)
        real(kind_phys), pointer :: spd1(:)
        integer, pointer :: kpbl(:)
        real(kind_phys), pointer :: prsi(:,:)
        real(kind_phys), pointer :: del(:,:)
        real(kind_phys), pointer :: prsl(:,:)
        real(kind_phys), pointer :: prslk(:,:)
        real(kind_phys), pointer :: phii(:,:)
        real(kind_phys), pointer :: phil(:,:)
        real(kind_phys), pointer :: delt
        real(kind_phys), pointer :: dusfc(:)
        real(kind_phys), pointer :: dvsfc(:)
        real(kind_phys), pointer :: dtsfc(:)
        real(kind_phys), pointer :: dqsfc(:)
        real(kind_phys), pointer :: dkt(:,:)
        real(kind_phys), pointer :: hpbl(:)
        integer, pointer :: kinver(:)
        real(kind_phys), pointer :: xkzm_m
        real(kind_phys), pointer :: xkzm_h
        real(kind_phys), pointer :: xkzm_s
        logical, pointer :: lprnt
        integer, pointer :: ipr
        integer, pointer :: me
        real(kind_phys), pointer :: grav
        real(kind_phys), pointer :: rd
        real(kind_phys), pointer :: cp
        real(kind_phys), pointer :: hvap
        real(kind_phys), pointer :: fv

        ierr = 0

        call c_f_pointer(ptr, cdata)


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
        

        call ccpp_field_get(cdata, 'vertical_dimension', km, ierr=ierr, kind=ckind, index=817)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve vertical_dimension from CCPP data structure')
            return
        end if
        if (kind(km).ne.ckind) then
            call ccpp_error('Kind mismatch for variable vertical_dimension')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'number_of_vertical_diffusion_tracers', ntrac, ierr=ierr, kind=ckind, index=590)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve number_of_vertical_diffusion_tracers from CCPP data structure')
            return
        end if
        if (kind(ntrac).ne.ckind) then
            call ccpp_error('Kind mismatch for variable number_of_vertical_diffusion_tracers')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'index_for_liquid_cloud_condensate', ntcw, ierr=ierr, kind=ckind, index=384)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve index_for_liquid_cloud_condensate from CCPP data structure')
            return
        end if
        if (kind(ntcw).ne.ckind) then
            call ccpp_error('Kind mismatch for variable index_for_liquid_cloud_condensate')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'number_of_tracers_for_cloud_condensate', ncnd, ierr=ierr, kind=ckind, index=586)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve number_of_tracers_for_cloud_condensate from CCPP data structure')
            return
        end if
        if (kind(ncnd).ne.ckind) then
            call ccpp_error('Kind mismatch for variable number_of_tracers_for_cloud_condensate')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'tendency_of_y_wind_due_to_model_physics', dv, ierr=ierr, dims=cdims, kind=ckind, index=785)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve tendency_of_y_wind_due_to_model_physics from CCPP data structure')
            return
        end if
        if (kind(dv).ne.ckind) then
            call ccpp_error('Kind mismatch for variable tendency_of_y_wind_due_to_model_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'tendency_of_x_wind_due_to_model_physics', du, ierr=ierr, dims=cdims, kind=ckind, index=782)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve tendency_of_x_wind_due_to_model_physics from CCPP data structure')
            return
        end if
        if (kind(du).ne.ckind) then
            call ccpp_error('Kind mismatch for variable tendency_of_x_wind_due_to_model_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'tendency_of_air_temperature_due_to_model_physics', tau, ierr=ierr, dims=cdims, kind=ckind, index=758)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve tendency_of_air_temperature_due_to_model_physics from CCPP data structure')
            return
        end if
        if (kind(tau).ne.ckind) then
            call ccpp_error('Kind mismatch for variable tendency_of_air_temperature_due_to_model_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'tendency_of_vertically_diffused_tracer_concentration', rtg, ierr=ierr, dims=cdims, kind=ckind, index=777)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve tendency_of_vertically_diffused_tracer_concentration from CCPP data structure')
            return
        end if
        if (kind(rtg).ne.ckind) then
            call ccpp_error('Kind mismatch for variable tendency_of_vertically_diffused_tracer_concentration')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'x_wind', u1, ierr=ierr, dims=cdims, kind=ckind, index=871)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve x_wind from CCPP data structure')
            return
        end if
        if (kind(u1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable x_wind')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'y_wind', v1, ierr=ierr, dims=cdims, kind=ckind, index=878)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve y_wind from CCPP data structure')
            return
        end if
        if (kind(v1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable y_wind')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'air_temperature', t1, ierr=ierr, dims=cdims, kind=ckind, index=50)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve air_temperature from CCPP data structure')
            return
        end if
        if (kind(t1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable air_temperature')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'vertically_diffused_tracer_concentration', q1, ierr=ierr, dims=cdims, kind=ckind, index=831)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve vertically_diffused_tracer_concentration from CCPP data structure')
            return
        end if
        if (kind(q1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable vertically_diffused_tracer_concentration')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'atmosphere_heat_diffusivity_from_shoc', tkh, ierr=ierr, dims=cdims, kind=ckind, index=72)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve atmosphere_heat_diffusivity_from_shoc from CCPP data structure')
            return
        end if
        if (kind(tkh).ne.ckind) then
            call ccpp_error('Kind mismatch for variable atmosphere_heat_diffusivity_from_shoc')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'prandtl_number', prnum, ierr=ierr, dims=cdims, kind=ckind, index=608)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve prandtl_number from CCPP data structure')
            return
        end if
        if (kind(prnum).ne.ckind) then
            call ccpp_error('Kind mismatch for variable prandtl_number')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'index_for_turbulent_kinetic_energy', ntke, ierr=ierr, kind=ckind, index=391)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve index_for_turbulent_kinetic_energy from CCPP data structure')
            return
        end if
        if (kind(ntke).ne.ckind) then
            call ccpp_error('Kind mismatch for variable index_for_turbulent_kinetic_energy')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'dimensionless_exner_function_at_lowest_model_interface', psk, ierr=ierr, dims=cdims, kind=ckind, index=220)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve dimensionless_exner_function_at_lowest_model_interface from CCPP data structure')
            return
        end if
        if (kind(psk).ne.ckind) then
            call ccpp_error('Kind mismatch for variable dimensionless_exner_function_at_lowest_model_interface')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'bulk_richardson_number_at_lowest_model_level', rbsoil, ierr=ierr, dims=cdims, kind=ckind, index=78)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve bulk_richardson_number_at_lowest_model_level from CCPP data structure')
            return
        end if
        if (kind(rbsoil).ne.ckind) then
            call ccpp_error('Kind mismatch for variable bulk_richardson_number_at_lowest_model_level')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_roughness_length', zorl, ierr=ierr, dims=cdims, kind=ckind, index=718)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_roughness_length from CCPP data structure')
            return
        end if
        if (kind(zorl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_roughness_length')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'x_wind_at_10m', u10m, ierr=ierr, dims=cdims, kind=ckind, index=872)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve x_wind_at_10m from CCPP data structure')
            return
        end if
        if (kind(u10m).ne.ckind) then
            call ccpp_error('Kind mismatch for variable x_wind_at_10m')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'y_wind_at_10m', v10m, ierr=ierr, dims=cdims, kind=ckind, index=879)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve y_wind_at_10m from CCPP data structure')
            return
        end if
        if (kind(v10m).ne.ckind) then
            call ccpp_error('Kind mismatch for variable y_wind_at_10m')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'Monin-Obukhov_similarity_function_for_momentum', fm, ierr=ierr, dims=cdims, kind=ckind, index=18)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve Monin-Obukhov_similarity_function_for_momentum from CCPP data structure')
            return
        end if
        if (kind(fm).ne.ckind) then
            call ccpp_error('Kind mismatch for variable Monin-Obukhov_similarity_function_for_momentum')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'Monin-Obukhov_similarity_function_for_heat', fh, ierr=ierr, dims=cdims, kind=ckind, index=16)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve Monin-Obukhov_similarity_function_for_heat from CCPP data structure')
            return
        end if
        if (kind(fh).ne.ckind) then
            call ccpp_error('Kind mismatch for variable Monin-Obukhov_similarity_function_for_heat')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_skin_temperature', tsea, ierr=ierr, dims=cdims, kind=ckind, index=721)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_skin_temperature from CCPP data structure')
            return
        end if
        if (kind(tsea).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_skin_temperature')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'kinematic_surface_upward_sensible_heat_flux', heat, ierr=ierr, dims=cdims, kind=ckind, index=455)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve kinematic_surface_upward_sensible_heat_flux from CCPP data structure')
            return
        end if
        if (kind(heat).ne.ckind) then
            call ccpp_error('Kind mismatch for variable kinematic_surface_upward_sensible_heat_flux')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'kinematic_surface_upward_latent_heat_flux', evap, ierr=ierr, dims=cdims, kind=ckind, index=454)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve kinematic_surface_upward_latent_heat_flux from CCPP data structure')
            return
        end if
        if (kind(evap).ne.ckind) then
            call ccpp_error('Kind mismatch for variable kinematic_surface_upward_latent_heat_flux')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_wind_stress', stress, ierr=ierr, dims=cdims, kind=ckind, index=745)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_wind_stress from CCPP data structure')
            return
        end if
        if (kind(stress).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_wind_stress')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'wind_speed_at_lowest_model_layer', spd1, ierr=ierr, dims=cdims, kind=ckind, index=870)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve wind_speed_at_lowest_model_layer from CCPP data structure')
            return
        end if
        if (kind(spd1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable wind_speed_at_lowest_model_layer')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'vertical_index_at_top_of_atmosphere_boundary_layer', kpbl, ierr=ierr, dims=cdims, kind=ckind, index=822)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve vertical_index_at_top_of_atmosphere_boundary_layer from CCPP data structure')
            return
        end if
        if (kind(kpbl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable vertical_index_at_top_of_atmosphere_boundary_layer')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'air_pressure_at_interface', prsi, ierr=ierr, dims=cdims, kind=ckind, index=45)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve air_pressure_at_interface from CCPP data structure')
            return
        end if
        if (kind(prsi).ne.ckind) then
            call ccpp_error('Kind mismatch for variable air_pressure_at_interface')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'air_pressure_difference_between_midlayers', del, ierr=ierr, dims=cdims, kind=ckind, index=49)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve air_pressure_difference_between_midlayers from CCPP data structure')
            return
        end if
        if (kind(del).ne.ckind) then
            call ccpp_error('Kind mismatch for variable air_pressure_difference_between_midlayers')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'air_pressure', prsl, ierr=ierr, dims=cdims, kind=ckind, index=44)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve air_pressure from CCPP data structure')
            return
        end if
        if (kind(prsl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable air_pressure')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'dimensionless_exner_function_at_model_layers', prslk, ierr=ierr, dims=cdims, kind=ckind, index=222)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve dimensionless_exner_function_at_model_layers from CCPP data structure')
            return
        end if
        if (kind(prslk).ne.ckind) then
            call ccpp_error('Kind mismatch for variable dimensionless_exner_function_at_model_layers')
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
        

        call ccpp_field_get(cdata, 'time_step_for_physics', delt, ierr=ierr, kind=ckind, index=793)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve time_step_for_physics from CCPP data structure')
            return
        end if
        if (kind(delt).ne.ckind) then
            call ccpp_error('Kind mismatch for variable time_step_for_physics')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'instantaneous_surface_x_momentum_flux', dusfc, ierr=ierr, dims=cdims, kind=ckind, index=437)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve instantaneous_surface_x_momentum_flux from CCPP data structure')
            return
        end if
        if (kind(dusfc).ne.ckind) then
            call ccpp_error('Kind mismatch for variable instantaneous_surface_x_momentum_flux')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'instantaneous_surface_y_momentum_flux', dvsfc, ierr=ierr, dims=cdims, kind=ckind, index=440)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve instantaneous_surface_y_momentum_flux from CCPP data structure')
            return
        end if
        if (kind(dvsfc).ne.ckind) then
            call ccpp_error('Kind mismatch for variable instantaneous_surface_y_momentum_flux')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'instantaneous_surface_upward_sensible_heat_flux', dtsfc, ierr=ierr, dims=cdims, kind=ckind, index=434)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve instantaneous_surface_upward_sensible_heat_flux from CCPP data structure')
            return
        end if
        if (kind(dtsfc).ne.ckind) then
            call ccpp_error('Kind mismatch for variable instantaneous_surface_upward_sensible_heat_flux')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'instantaneous_surface_upward_latent_heat_flux', dqsfc, ierr=ierr, dims=cdims, kind=ckind, index=431)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve instantaneous_surface_upward_latent_heat_flux from CCPP data structure')
            return
        end if
        if (kind(dqsfc).ne.ckind) then
            call ccpp_error('Kind mismatch for variable instantaneous_surface_upward_latent_heat_flux')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'atmosphere_heat_diffusivity', dkt, ierr=ierr, dims=cdims, kind=ckind, index=68)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve atmosphere_heat_diffusivity from CCPP data structure')
            return
        end if
        if (kind(dkt).ne.ckind) then
            call ccpp_error('Kind mismatch for variable atmosphere_heat_diffusivity')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'atmosphere_boundary_layer_thickness', hpbl, ierr=ierr, dims=cdims, kind=ckind, index=66)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve atmosphere_boundary_layer_thickness from CCPP data structure')
            return
        end if
        if (kind(hpbl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable atmosphere_boundary_layer_thickness')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'index_of_highest_temperature_inversion', kinver, ierr=ierr, dims=cdims, kind=ckind, index=397)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve index_of_highest_temperature_inversion from CCPP data structure')
            return
        end if
        if (kind(kinver).ne.ckind) then
            call ccpp_error('Kind mismatch for variable index_of_highest_temperature_inversion')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'atmosphere_momentum_diffusivity_background', xkzm_m, ierr=ierr, kind=ckind, index=73)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve atmosphere_momentum_diffusivity_background from CCPP data structure')
            return
        end if
        if (kind(xkzm_m).ne.ckind) then
            call ccpp_error('Kind mismatch for variable atmosphere_momentum_diffusivity_background')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'atmosphere_heat_diffusivity_background', xkzm_h, ierr=ierr, kind=ckind, index=69)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve atmosphere_heat_diffusivity_background from CCPP data structure')
            return
        end if
        if (kind(xkzm_h).ne.ckind) then
            call ccpp_error('Kind mismatch for variable atmosphere_heat_diffusivity_background')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'diffusivity_background_sigma_level', xkzm_s, ierr=ierr, kind=ckind, index=219)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve diffusivity_background_sigma_level from CCPP data structure')
            return
        end if
        if (kind(xkzm_s).ne.ckind) then
            call ccpp_error('Kind mismatch for variable diffusivity_background_sigma_level')
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
        

        call ccpp_field_get(cdata, 'mpi_rank', me, ierr=ierr, kind=ckind, index=558)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve mpi_rank from CCPP data structure')
            return
        end if
        if (kind(me).ne.ckind) then
            call ccpp_error('Kind mismatch for variable mpi_rank')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'gravitational_acceleration', grav, ierr=ierr, kind=ckind, index=355)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve gravitational_acceleration from CCPP data structure')
            return
        end if
        if (kind(grav).ne.ckind) then
            call ccpp_error('Kind mismatch for variable gravitational_acceleration')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'gas_constant_dry_air', rd, ierr=ierr, kind=ckind, index=346)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve gas_constant_dry_air from CCPP data structure')
            return
        end if
        if (kind(rd).ne.ckind) then
            call ccpp_error('Kind mismatch for variable gas_constant_dry_air')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'specific_heat_of_dry_air_at_constant_pressure', cp, ierr=ierr, kind=ckind, index=664)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve specific_heat_of_dry_air_at_constant_pressure from CCPP data structure')
            return
        end if
        if (kind(cp).ne.ckind) then
            call ccpp_error('Kind mismatch for variable specific_heat_of_dry_air_at_constant_pressure')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'latent_heat_of_vaporization_of_water_at_0C', hvap, ierr=ierr, kind=ckind, index=459)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve latent_heat_of_vaporization_of_water_at_0C from CCPP data structure')
            return
        end if
        if (kind(hvap).ne.ckind) then
            call ccpp_error('Kind mismatch for variable latent_heat_of_vaporization_of_water_at_0C')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'ratio_of_vapor_to_dry_air_gas_constants_minus_one', fv, ierr=ierr, kind=ckind, index=625)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve ratio_of_vapor_to_dry_air_gas_constants_minus_one from CCPP data structure')
            return
        end if
        if (kind(fv).ne.ckind) then
            call ccpp_error('Kind mismatch for variable ratio_of_vapor_to_dry_air_gas_constants_minus_one')
            ierr = 1
            return
        end if
#endif
        

        call moninshoc_run(ix=ix,im=im,km=km,ntrac=ntrac,ntcw=ntcw,ncnd=ncnd,dv=dv,du=du,tau=tau,rtg=rtg, &
                  u1=u1,v1=v1,t1=t1,q1=q1,tkh=tkh,prnum=prnum,ntke=ntke,psk=psk,rbsoil=rbsoil, &
                  zorl=zorl,u10m=u10m,v10m=v10m,fm=fm,fh=fh,tsea=tsea,heat=heat,evap=evap, &
                  stress=stress,spd1=spd1,kpbl=kpbl,prsi=prsi,del=del,prsl=prsl,prslk=prslk, &
                  phii=phii,phil=phil,delt=delt,dusfc=dusfc,dvsfc=dvsfc,dtsfc=dtsfc,dqsfc=dqsfc, &
                  dkt=dkt,hpbl=hpbl,kinver=kinver,xkzm_m=xkzm_m,xkzm_h=xkzm_h,xkzm_s=xkzm_s, &
                  lprnt=lprnt,ipr=ipr,me=me,grav=grav,rd=rd,cp=cp,hvap=hvap,fv=fv,errmsg=cdata%errmsg, &
                  errflg=cdata%errflg)
        ierr=cdata%errflg

    end function moninshoc_run_cap

    function moninshoc_init_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call moninshoc_init()
        

    end function moninshoc_init_cap

    function moninshoc_finalize_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call moninshoc_finalize()
        

    end function moninshoc_finalize_cap
end module moninshoc_cap

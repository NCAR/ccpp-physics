
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
!! @brief Auto-generated cap module for the GFS_suite_interstitial_2 scheme
!!
!
module GFS_suite_interstitial_2_cap

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_f_pointer, c_ptr, c_int32_t
    use            :: ccpp_types,                                      &
                      only: ccpp_t, CCPP_GENERIC_KIND
    use            :: ccpp_fields,                                     &
                      only: ccpp_field_get
    use            :: ccpp_errors,                                     &
                      only: ccpp_error, ccpp_debug
    use            :: GFS_suite_interstitial_2, &
                      only: GFS_suite_interstitial_2_init,GFS_suite_interstitial_2_finalize,GFS_suite_interstitial_2_run
    ! Other modules required, e.g. type definitions
    use machine, only: kind_phys

    implicit none

    private
    public :: GFS_suite_interstitial_2_init_cap,GFS_suite_interstitial_2_finalize_cap,GFS_suite_interstitial_2_run_cap

    contains


    function GFS_suite_interstitial_2_init_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call GFS_suite_interstitial_2_init()
        

    end function GFS_suite_interstitial_2_init_cap

    function GFS_suite_interstitial_2_finalize_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call GFS_suite_interstitial_2_finalize()
        

    end function GFS_suite_interstitial_2_finalize_cap

    function GFS_suite_interstitial_2_run_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind
        integer, pointer :: im
        integer, pointer :: levs
        logical, pointer :: lssav
        logical, pointer :: ldiag3d
        logical, pointer :: lsidea
        logical, pointer :: cplflx
        logical, pointer :: flag_cice(:)
        logical, pointer :: shal_cnv
        logical, pointer :: old_monin
        logical, pointer :: mstrat
        logical, pointer :: do_shoc
        integer, pointer :: imfshalcnv
        real(kind_phys), pointer :: dtf
        real(kind_phys), pointer :: xcosz(:)
        real(kind_phys), pointer :: adjsfcdsw(:)
        real(kind_phys), pointer :: adjsfcdlw(:)
        real(kind_phys), pointer :: pgr(:)
        real(kind_phys), pointer :: ulwsfc_cice(:)
        real(kind_phys), pointer :: lwhd(:,:,:)
        real(kind_phys), pointer :: htrsw(:,:)
        real(kind_phys), pointer :: htrlw(:,:)
        real(kind_phys), pointer :: xmu(:)
        real(kind_phys), pointer :: ctei_rm(:)
        real(kind_phys), pointer :: work1(:)
        real(kind_phys), pointer :: work2(:)
        real(kind_phys), pointer :: prsi(:,:)
        real(kind_phys), pointer :: tgrs(:,:)
        real(kind_phys), pointer :: prsl(:,:)
        real(kind_phys), pointer :: qgrs_water_vapor(:,:)
        real(kind_phys), pointer :: qgrs_cloud_water(:,:)
        real(kind_phys), pointer :: cp
        real(kind_phys), pointer :: hvap
        real(kind_phys), pointer :: prslk(:,:)
        real(kind_phys), pointer :: suntim(:)
        real(kind_phys), pointer :: adjsfculw(:)
        real(kind_phys), pointer :: dlwsfc(:)
        real(kind_phys), pointer :: ulwsfc(:)
        real(kind_phys), pointer :: psmean(:)
        real(kind_phys), pointer :: dt3dt_lw(:,:)
        real(kind_phys), pointer :: dt3dt_sw(:,:)
        real(kind_phys), pointer :: dt3dt_pbl(:,:)
        real(kind_phys), pointer :: dt3dt_dcnv(:,:)
        real(kind_phys), pointer :: dt3dt_scnv(:,:)
        real(kind_phys), pointer :: dt3dt_mp(:,:)
        real(kind_phys), pointer :: ctei_rml(:)
        real(kind_phys), pointer :: ctei_r(:)
        integer, pointer :: kinver(:)

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
        

        call ccpp_field_get(cdata, 'vertical_dimension', levs, ierr=ierr, kind=ckind, index=817)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve vertical_dimension from CCPP data structure')
            return
        end if
        if (kind(levs).ne.ckind) then
            call ccpp_error('Kind mismatch for variable vertical_dimension')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'flag_diagnostics', lssav, ierr=ierr, kind=ckind, index=263)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_diagnostics from CCPP data structure')
            return
        end if
        if (kind(lssav).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_diagnostics')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'flag_diagnostics_3D', ldiag3d, ierr=ierr, kind=ckind, index=264)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_diagnostics_3D from CCPP data structure')
            return
        end if
        if (kind(ldiag3d).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_diagnostics_3D')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'flag_idealized_physics', lsidea, ierr=ierr, kind=ckind, index=330)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_idealized_physics from CCPP data structure')
            return
        end if
        if (kind(lsidea).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_idealized_physics')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'flag_for_flux_coupling', cplflx, ierr=ierr, kind=ckind, index=278)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_flux_coupling from CCPP data structure')
            return
        end if
        if (kind(cplflx).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_flux_coupling')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'flag_for_cice', flag_cice, ierr=ierr, dims=cdims, kind=ckind, index=273)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_cice from CCPP data structure')
            return
        end if
        if (kind(flag_cice).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_cice')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'flag_for_shallow_convection', shal_cnv, ierr=ierr, kind=ckind, index=312)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_shallow_convection from CCPP data structure')
            return
        end if
        if (kind(shal_cnv).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_shallow_convection')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'flag_for_old_PBL_scheme', old_monin, ierr=ierr, kind=ckind, index=300)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_old_PBL_scheme from CCPP data structure')
            return
        end if
        if (kind(old_monin).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_old_PBL_scheme')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'flag_for_moorthi_stratus', mstrat, ierr=ierr, kind=ckind, index=296)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_moorthi_stratus from CCPP data structure')
            return
        end if
        if (kind(mstrat).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_moorthi_stratus')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'flag_for_shoc', do_shoc, ierr=ierr, kind=ckind, index=313)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_shoc from CCPP data structure')
            return
        end if
        if (kind(do_shoc).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_shoc')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'flag_for_mass_flux_shallow_convection_scheme', imfshalcnv, ierr=ierr, kind=ckind, index=291)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_mass_flux_shallow_convection_scheme from CCPP data structure')
            return
        end if
        if (kind(imfshalcnv).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_mass_flux_shallow_convection_scheme')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'time_step_for_dynamics', dtf, ierr=ierr, kind=ckind, index=792)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve time_step_for_dynamics from CCPP data structure')
            return
        end if
        if (kind(dtf).ne.ckind) then
            call ccpp_error('Kind mismatch for variable time_step_for_dynamics')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'instantaneous_cosine_of_zenith_angle', xcosz, ierr=ierr, dims=cdims, kind=ckind, index=408)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve instantaneous_cosine_of_zenith_angle from CCPP data structure')
            return
        end if
        if (kind(xcosz).ne.ckind) then
            call ccpp_error('Kind mismatch for variable instantaneous_cosine_of_zenith_angle')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_downwelling_shortwave_flux', adjsfcdsw, ierr=ierr, dims=cdims, kind=ckind, index=700)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_downwelling_shortwave_flux from CCPP data structure')
            return
        end if
        if (kind(adjsfcdsw).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_downwelling_shortwave_flux')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_downwelling_longwave_flux', adjsfcdlw, ierr=ierr, dims=cdims, kind=ckind, index=697)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_downwelling_longwave_flux from CCPP data structure')
            return
        end if
        if (kind(adjsfcdlw).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_downwelling_longwave_flux')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_air_pressure', pgr, ierr=ierr, dims=cdims, kind=ckind, index=677)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_air_pressure from CCPP data structure')
            return
        end if
        if (kind(pgr).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_air_pressure')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_upwelling_longwave_flux_for_cice', ulwsfc_cice, ierr=ierr, dims=cdims, kind=ckind, index=742)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_upwelling_longwave_flux_for_cice from CCPP data structure')
            return
        end if
        if (kind(ulwsfc_cice).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_upwelling_longwave_flux_for_cice')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'tendency_of_air_temperature_due_to_longwave_heating_for_idea', lwhd, ierr=ierr, dims=cdims, kind=ckind, index=755)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve tendency_of_air_temperature_due_to_longwave_heating_for_idea from CCPP data structure')
            return
        end if
        if (kind(lwhd).ne.ckind) then
            call ccpp_error('Kind mismatch for variable tendency_of_air_temperature_due_to_longwave_heating_for_idea')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'tendency_of_air_temperature_due_to_shortwave_heating_on_radiation_timestep', htrsw, ierr=ierr, dims=cdims, kind=ckind, index=763)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve tendency_of_air_temperature_due_to_shortwave_heating_on_radiation_timestep from CCPP data structure')
            return
        end if
        if (kind(htrsw).ne.ckind) then
            call ccpp_error('Kind mismatch for variable tendency_of_air_temperature_due_to_shortwave_heating_on_radiation_timestep')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'tendency_of_air_temperature_due_to_longwave_heating_on_radiation_timestep', htrlw, ierr=ierr, dims=cdims, kind=ckind, index=757)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve tendency_of_air_temperature_due_to_longwave_heating_on_radiation_timestep from CCPP data structure')
            return
        end if
        if (kind(htrlw).ne.ckind) then
            call ccpp_error('Kind mismatch for variable tendency_of_air_temperature_due_to_longwave_heating_on_radiation_timestep')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'zenith_angle_temporal_adjustment_factor_for_shortwave_fluxes', xmu, ierr=ierr, dims=cdims, kind=ckind, index=885)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve zenith_angle_temporal_adjustment_factor_for_shortwave_fluxes from CCPP data structure')
            return
        end if
        if (kind(xmu).ne.ckind) then
            call ccpp_error('Kind mismatch for variable zenith_angle_temporal_adjustment_factor_for_shortwave_fluxes')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'critical_cloud_top_entrainment_instability_criteria', ctei_rm, ierr=ierr, dims=cdims, kind=ckind, index=140)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve critical_cloud_top_entrainment_instability_criteria from CCPP data structure')
            return
        end if
        if (kind(ctei_rm).ne.ckind) then
            call ccpp_error('Kind mismatch for variable critical_cloud_top_entrainment_instability_criteria')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'grid_size_related_coefficient_used_in_scale-sensitive_schemes', work1, ierr=ierr, dims=cdims, kind=ckind, index=357)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve grid_size_related_coefficient_used_in_scale-sensitive_schemes from CCPP data structure')
            return
        end if
        if (kind(work1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable grid_size_related_coefficient_used_in_scale-sensitive_schemes')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'grid_size_related_coefficient_used_in_scale-sensitive_schemes_complement', work2, ierr=ierr, dims=cdims, kind=ckind, index=358)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve grid_size_related_coefficient_used_in_scale-sensitive_schemes_complement from CCPP data structure')
            return
        end if
        if (kind(work2).ne.ckind) then
            call ccpp_error('Kind mismatch for variable grid_size_related_coefficient_used_in_scale-sensitive_schemes_complement')
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
        

        call ccpp_field_get(cdata, 'air_temperature', tgrs, ierr=ierr, dims=cdims, kind=ckind, index=50)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve air_temperature from CCPP data structure')
            return
        end if
        if (kind(tgrs).ne.ckind) then
            call ccpp_error('Kind mismatch for variable air_temperature')
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
        

        call ccpp_field_get(cdata, 'water_vapor_specific_humidity', qgrs_water_vapor, ierr=ierr, dims=cdims, kind=ckind, index=852)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve water_vapor_specific_humidity from CCPP data structure')
            return
        end if
        if (kind(qgrs_water_vapor).ne.ckind) then
            call ccpp_error('Kind mismatch for variable water_vapor_specific_humidity')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cloud_condensed_water_mixing_ratio', qgrs_cloud_water, ierr=ierr, dims=cdims, kind=ckind, index=90)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cloud_condensed_water_mixing_ratio from CCPP data structure')
            return
        end if
        if (kind(qgrs_cloud_water).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cloud_condensed_water_mixing_ratio')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

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
        

        call ccpp_field_get(cdata, 'duration_of_sunshine', suntim, ierr=ierr, dims=cdims, kind=ckind, index=235)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve duration_of_sunshine from CCPP data structure')
            return
        end if
        if (kind(suntim).ne.ckind) then
            call ccpp_error('Kind mismatch for variable duration_of_sunshine')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_upwelling_longwave_flux', adjsfculw, ierr=ierr, dims=cdims, kind=ckind, index=741)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_upwelling_longwave_flux from CCPP data structure')
            return
        end if
        if (kind(adjsfculw).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_upwelling_longwave_flux')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cumulative_surface_downwelling_longwave_flux_multiplied_by_timestep', dlwsfc, ierr=ierr, dims=cdims, kind=ckind, index=184)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cumulative_surface_downwelling_longwave_flux_multiplied_by_timestep from CCPP data structure')
            return
        end if
        if (kind(dlwsfc).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cumulative_surface_downwelling_longwave_flux_multiplied_by_timestep')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cumulative_surface_upwelling_longwave_flux_multiplied_by_timestep', ulwsfc, ierr=ierr, dims=cdims, kind=ckind, index=200)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cumulative_surface_upwelling_longwave_flux_multiplied_by_timestep from CCPP data structure')
            return
        end if
        if (kind(ulwsfc).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cumulative_surface_upwelling_longwave_flux_multiplied_by_timestep')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cumulative_surface_pressure_multiplied_by_timestep', psmean, ierr=ierr, dims=cdims, kind=ckind, index=193)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cumulative_surface_pressure_multiplied_by_timestep from CCPP data structure')
            return
        end if
        if (kind(psmean).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cumulative_surface_pressure_multiplied_by_timestep')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cumulative_change_in_temperature_due_to_longwave_radiation', dt3dt_lw, ierr=ierr, dims=cdims, kind=ckind, index=157)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cumulative_change_in_temperature_due_to_longwave_radiation from CCPP data structure')
            return
        end if
        if (kind(dt3dt_lw).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cumulative_change_in_temperature_due_to_longwave_radiation')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cumulative_change_in_temperature_due_to_shortwave_radiation_and_orographic_gravity_wave_drag', dt3dt_sw, ierr=ierr, dims=cdims, kind=ckind, index=160)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cumulative_change_in_temperature_due_to_shortwave_radiation_and_orographic_gravity_wave_drag from CCPP data structure')
            return
        end if
        if (kind(dt3dt_sw).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cumulative_change_in_temperature_due_to_shortwave_radiation_and_orographic_gravity_wave_drag')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cumulative_change_in_temperature_due_to_PBL', dt3dt_pbl, ierr=ierr, dims=cdims, kind=ckind, index=155)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cumulative_change_in_temperature_due_to_PBL from CCPP data structure')
            return
        end if
        if (kind(dt3dt_pbl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cumulative_change_in_temperature_due_to_PBL')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cumulative_change_in_temperature_due_to_deep_convection', dt3dt_dcnv, ierr=ierr, dims=cdims, kind=ckind, index=156)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cumulative_change_in_temperature_due_to_deep_convection from CCPP data structure')
            return
        end if
        if (kind(dt3dt_dcnv).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cumulative_change_in_temperature_due_to_deep_convection')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cumulative_change_in_temperature_due_to_shal_convection', dt3dt_scnv, ierr=ierr, dims=cdims, kind=ckind, index=159)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cumulative_change_in_temperature_due_to_shal_convection from CCPP data structure')
            return
        end if
        if (kind(dt3dt_scnv).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cumulative_change_in_temperature_due_to_shal_convection')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cumulative_change_in_temperature_due_to_microphysics', dt3dt_mp, ierr=ierr, dims=cdims, kind=ckind, index=158)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cumulative_change_in_temperature_due_to_microphysics from CCPP data structure')
            return
        end if
        if (kind(dt3dt_mp).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cumulative_change_in_temperature_due_to_microphysics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'grid_sensitive_critical_cloud_top_entrainment_instability_criteria', ctei_rml, ierr=ierr, dims=cdims, kind=ckind, index=356)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve grid_sensitive_critical_cloud_top_entrainment_instability_criteria from CCPP data structure')
            return
        end if
        if (kind(ctei_rml).ne.ckind) then
            call ccpp_error('Kind mismatch for variable grid_sensitive_critical_cloud_top_entrainment_instability_criteria')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cloud_top_entrainment_instability_value', ctei_r, ierr=ierr, dims=cdims, kind=ckind, index=110)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cloud_top_entrainment_instability_value from CCPP data structure')
            return
        end if
        if (kind(ctei_r).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cloud_top_entrainment_instability_value')
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
        

        call GFS_suite_interstitial_2_run(im=im,levs=levs,lssav=lssav,ldiag3d=ldiag3d,lsidea=lsidea,cplflx=cplflx, &
                  flag_cice=flag_cice,shal_cnv=shal_cnv,old_monin=old_monin,mstrat=mstrat, &
                  do_shoc=do_shoc,imfshalcnv=imfshalcnv,dtf=dtf,xcosz=xcosz,adjsfcdsw=adjsfcdsw, &
                  adjsfcdlw=adjsfcdlw,pgr=pgr,ulwsfc_cice=ulwsfc_cice,lwhd=lwhd,htrsw=htrsw, &
                  htrlw=htrlw,xmu=xmu,ctei_rm=ctei_rm,work1=work1,work2=work2,prsi=prsi,tgrs=tgrs, &
                  prsl=prsl,qgrs_water_vapor=qgrs_water_vapor,qgrs_cloud_water=qgrs_cloud_water, &
                  cp=cp,hvap=hvap,prslk=prslk,suntim=suntim,adjsfculw=adjsfculw,dlwsfc=dlwsfc, &
                  ulwsfc=ulwsfc,psmean=psmean,dt3dt_lw=dt3dt_lw,dt3dt_sw=dt3dt_sw,dt3dt_pbl=dt3dt_pbl, &
                  dt3dt_dcnv=dt3dt_dcnv,dt3dt_scnv=dt3dt_scnv,dt3dt_mp=dt3dt_mp,ctei_rml=ctei_rml, &
                  ctei_r=ctei_r,kinver=kinver,errmsg=cdata%errmsg,errflg=cdata%errflg)
        ierr=cdata%errflg

    end function GFS_suite_interstitial_2_run_cap
end module GFS_suite_interstitial_2_cap

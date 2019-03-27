
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
!! @brief Auto-generated cap module for the GFS_surface_generic_pre scheme
!!
!
module GFS_surface_generic_pre_cap

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_f_pointer, c_ptr, c_int32_t
    use            :: ccpp_types,                                      &
                      only: ccpp_t, CCPP_GENERIC_KIND
    use            :: ccpp_fields,                                     &
                      only: ccpp_field_get
    use            :: ccpp_errors,                                     &
                      only: ccpp_error, ccpp_debug
    use            :: GFS_surface_generic_pre, &
                      only: GFS_surface_generic_pre_init,GFS_surface_generic_pre_run,GFS_surface_generic_pre_finalize
    ! Other modules required, e.g. type definitions
    use machine, only: kind_phys

    implicit none

    private
    public :: GFS_surface_generic_pre_init_cap,GFS_surface_generic_pre_run_cap,GFS_surface_generic_pre_finalize_cap

    contains


    function GFS_surface_generic_pre_init_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call GFS_surface_generic_pre_init()
        

    end function GFS_surface_generic_pre_init_cap

    function GFS_surface_generic_pre_run_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind
        integer, pointer :: im
        integer, pointer :: levs
        real(kind_phys), pointer :: vfrac(:)
        integer, pointer :: islmsk(:)
        integer, pointer :: isot
        integer, pointer :: ivegsrc
        real(kind_phys), pointer :: stype(:)
        real(kind_phys), pointer :: vtype(:)
        real(kind_phys), pointer :: slope(:)
        real(kind_phys), pointer :: prsik_1(:)
        real(kind_phys), pointer :: prslk_1(:)
        real(kind_phys), pointer :: semis(:)
        real(kind_phys), pointer :: adjsfcdlw(:)
        real(kind_phys), pointer :: tsfc(:)
        real(kind_phys), pointer :: phil(:,:)
        real(kind_phys), pointer :: con_g
        real(kind_phys), pointer :: sigmaf(:)
        integer, pointer :: soiltyp(:)
        integer, pointer :: vegtype(:)
        integer, pointer :: slopetyp(:)
        real(kind_phys), pointer :: work3(:)
        real(kind_phys), pointer :: gabsbdlw(:)
        real(kind_phys), pointer :: tsurf(:)
        real(kind_phys), pointer :: zlvl(:)
        logical, pointer :: do_sppt
        real(kind_phys), pointer :: dtdtr(:,:)
        real(kind_phys), pointer :: drain_cpl(:)
        real(kind_phys), pointer :: dsnow_cpl(:)
        real(kind_phys), pointer :: rain_cpl(:)
        real(kind_phys), pointer :: snow_cpl(:)
        logical, pointer :: do_sfcperts
        integer, pointer :: nsfcpert
        real(kind_phys), pointer :: sfc_wts(:,:)
        real(kind_phys), pointer :: pertz0(:)
        real(kind_phys), pointer :: pertzt(:)
        real(kind_phys), pointer :: pertshc(:)
        real(kind_phys), pointer :: pertlai(:)
        real(kind_phys), pointer :: pertvegf(:)
        real(kind_phys), pointer :: z01d(:)
        real(kind_phys), pointer :: zt1d(:)
        real(kind_phys), pointer :: bexp1d(:)
        real(kind_phys), pointer :: xlai1d(:)
        real(kind_phys), pointer :: vegf1d(:)

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
        

        call ccpp_field_get(cdata, 'vegetation_area_fraction', vfrac, ierr=ierr, dims=cdims, kind=ckind, index=813)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve vegetation_area_fraction from CCPP data structure')
            return
        end if
        if (kind(vfrac).ne.ckind) then
            call ccpp_error('Kind mismatch for variable vegetation_area_fraction')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'sea_land_ice_mask', islmsk, ierr=ierr, dims=cdims, kind=ckind, index=631)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve sea_land_ice_mask from CCPP data structure')
            return
        end if
        if (kind(islmsk).ne.ckind) then
            call ccpp_error('Kind mismatch for variable sea_land_ice_mask')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'soil_type_dataset_choice', isot, ierr=ierr, kind=ckind, index=659)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve soil_type_dataset_choice from CCPP data structure')
            return
        end if
        if (kind(isot).ne.ckind) then
            call ccpp_error('Kind mismatch for variable soil_type_dataset_choice')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'vegetation_type_dataset_choice', ivegsrc, ierr=ierr, kind=ckind, index=816)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve vegetation_type_dataset_choice from CCPP data structure')
            return
        end if
        if (kind(ivegsrc).ne.ckind) then
            call ccpp_error('Kind mismatch for variable vegetation_type_dataset_choice')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'soil_type_classification_real', stype, ierr=ierr, dims=cdims, kind=ckind, index=658)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve soil_type_classification_real from CCPP data structure')
            return
        end if
        if (kind(stype).ne.ckind) then
            call ccpp_error('Kind mismatch for variable soil_type_classification_real')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'vegetation_type_classification_real', vtype, ierr=ierr, dims=cdims, kind=ckind, index=815)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve vegetation_type_classification_real from CCPP data structure')
            return
        end if
        if (kind(vtype).ne.ckind) then
            call ccpp_error('Kind mismatch for variable vegetation_type_classification_real')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_slope_classification_real', slope, ierr=ierr, dims=cdims, kind=ckind, index=725)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_slope_classification_real from CCPP data structure')
            return
        end if
        if (kind(slope).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_slope_classification_real')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'dimensionless_exner_function_at_lowest_model_interface', prsik_1, ierr=ierr, dims=cdims, kind=ckind, index=220)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve dimensionless_exner_function_at_lowest_model_interface from CCPP data structure')
            return
        end if
        if (kind(prsik_1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable dimensionless_exner_function_at_lowest_model_interface')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'dimensionless_exner_function_at_lowest_model_layer', prslk_1, ierr=ierr, dims=cdims, kind=ckind, index=221)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve dimensionless_exner_function_at_lowest_model_layer from CCPP data structure')
            return
        end if
        if (kind(prslk_1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable dimensionless_exner_function_at_lowest_model_layer')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_longwave_emissivity', semis, ierr=ierr, dims=cdims, kind=ckind, index=714)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_longwave_emissivity from CCPP data structure')
            return
        end if
        if (kind(semis).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_longwave_emissivity')
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
        

        call ccpp_field_get(cdata, 'surface_skin_temperature', tsfc, ierr=ierr, dims=cdims, kind=ckind, index=721)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_skin_temperature from CCPP data structure')
            return
        end if
        if (kind(tsfc).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_skin_temperature')
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
        

        call ccpp_field_get(cdata, 'gravitational_acceleration', con_g, ierr=ierr, kind=ckind, index=355)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve gravitational_acceleration from CCPP data structure')
            return
        end if
        if (kind(con_g).ne.ckind) then
            call ccpp_error('Kind mismatch for variable gravitational_acceleration')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'bounded_vegetation_area_fraction', sigmaf, ierr=ierr, dims=cdims, kind=ckind, index=77)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve bounded_vegetation_area_fraction from CCPP data structure')
            return
        end if
        if (kind(sigmaf).ne.ckind) then
            call ccpp_error('Kind mismatch for variable bounded_vegetation_area_fraction')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'soil_type_classification', soiltyp, ierr=ierr, dims=cdims, kind=ckind, index=657)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve soil_type_classification from CCPP data structure')
            return
        end if
        if (kind(soiltyp).ne.ckind) then
            call ccpp_error('Kind mismatch for variable soil_type_classification')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'vegetation_type_classification', vegtype, ierr=ierr, dims=cdims, kind=ckind, index=814)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve vegetation_type_classification from CCPP data structure')
            return
        end if
        if (kind(vegtype).ne.ckind) then
            call ccpp_error('Kind mismatch for variable vegetation_type_classification')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_slope_classification', slopetyp, ierr=ierr, dims=cdims, kind=ckind, index=724)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_slope_classification from CCPP data structure')
            return
        end if
        if (kind(slopetyp).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_slope_classification')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'ratio_of_exner_function_between_midlayer_and_interface_at_lowest_model_layer', work3, ierr=ierr, dims=cdims, kind=ckind, index=623)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve ratio_of_exner_function_between_midlayer_and_interface_at_lowest_model_layer from CCPP data structure')
            return
        end if
        if (kind(work3).ne.ckind) then
            call ccpp_error('Kind mismatch for variable ratio_of_exner_function_between_midlayer_and_interface_at_lowest_model_layer')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_downwelling_longwave_flux_absorbed_by_ground', gabsbdlw, ierr=ierr, dims=cdims, kind=ckind, index=698)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_downwelling_longwave_flux_absorbed_by_ground from CCPP data structure')
            return
        end if
        if (kind(gabsbdlw).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_downwelling_longwave_flux_absorbed_by_ground')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_skin_temperature_after_iteration', tsurf, ierr=ierr, dims=cdims, kind=ckind, index=722)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_skin_temperature_after_iteration from CCPP data structure')
            return
        end if
        if (kind(tsurf).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_skin_temperature_after_iteration')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'height_above_ground_at_lowest_model_layer', zlvl, ierr=ierr, dims=cdims, kind=ckind, index=360)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve height_above_ground_at_lowest_model_layer from CCPP data structure')
            return
        end if
        if (kind(zlvl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable height_above_ground_at_lowest_model_layer')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'flag_for_stochastic_surface_physics_perturbations', do_sppt, ierr=ierr, kind=ckind, index=319)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_stochastic_surface_physics_perturbations from CCPP data structure')
            return
        end if
        if (kind(do_sppt).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_stochastic_surface_physics_perturbations')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'tendency_of_air_temperature_due_to_radiative_heating_on_physics_time_step', dtdtr, ierr=ierr, dims=cdims, kind=ckind, index=760)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve tendency_of_air_temperature_due_to_radiative_heating_on_physics_time_step from CCPP data structure')
            return
        end if
        if (kind(dtdtr).ne.ckind) then
            call ccpp_error('Kind mismatch for variable tendency_of_air_temperature_due_to_radiative_heating_on_physics_time_step')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'tendency_of_lwe_thickness_of_precipitation_amount_for_coupling', drain_cpl, ierr=ierr, dims=cdims, kind=ckind, index=772)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve tendency_of_lwe_thickness_of_precipitation_amount_for_coupling from CCPP data structure')
            return
        end if
        if (kind(drain_cpl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable tendency_of_lwe_thickness_of_precipitation_amount_for_coupling')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'tendency_of_lwe_thickness_of_snow_amount_for_coupling', dsnow_cpl, ierr=ierr, dims=cdims, kind=ckind, index=773)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve tendency_of_lwe_thickness_of_snow_amount_for_coupling from CCPP data structure')
            return
        end if
        if (kind(dsnow_cpl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable tendency_of_lwe_thickness_of_snow_amount_for_coupling')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'lwe_thickness_of_precipitation_amount_for_coupling', rain_cpl, ierr=ierr, dims=cdims, kind=ckind, index=489)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve lwe_thickness_of_precipitation_amount_for_coupling from CCPP data structure')
            return
        end if
        if (kind(rain_cpl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable lwe_thickness_of_precipitation_amount_for_coupling')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'lwe_thickness_of_snow_amount_for_coupling', snow_cpl, ierr=ierr, dims=cdims, kind=ckind, index=493)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve lwe_thickness_of_snow_amount_for_coupling from CCPP data structure')
            return
        end if
        if (kind(snow_cpl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable lwe_thickness_of_snow_amount_for_coupling')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'flag_for_stochastic_surface_perturbations', do_sfcperts, ierr=ierr, kind=ckind, index=318)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_stochastic_surface_perturbations from CCPP data structure')
            return
        end if
        if (kind(do_sfcperts).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_stochastic_surface_perturbations')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'number_of_surface_perturbations', nsfcpert, ierr=ierr, kind=ckind, index=582)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve number_of_surface_perturbations from CCPP data structure')
            return
        end if
        if (kind(nsfcpert).ne.ckind) then
            call ccpp_error('Kind mismatch for variable number_of_surface_perturbations')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'weights_for_stochastic_surface_physics_perturbation', sfc_wts, ierr=ierr, dims=cdims, kind=ckind, index=869)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve weights_for_stochastic_surface_physics_perturbation from CCPP data structure')
            return
        end if
        if (kind(sfc_wts).ne.ckind) then
            call ccpp_error('Kind mismatch for variable weights_for_stochastic_surface_physics_perturbation')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'magnitude_of_perturbation_of_momentum_roughness_length', pertz0, ierr=ierr, dims=cdims, kind=ckind, index=498)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve magnitude_of_perturbation_of_momentum_roughness_length from CCPP data structure')
            return
        end if
        if (kind(pertz0).ne.ckind) then
            call ccpp_error('Kind mismatch for variable magnitude_of_perturbation_of_momentum_roughness_length')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'magnitude_of_perturbation_of_heat_to_momentum_roughness_length_ratio', pertzt, ierr=ierr, dims=cdims, kind=ckind, index=496)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve magnitude_of_perturbation_of_heat_to_momentum_roughness_length_ratio from CCPP data structure')
            return
        end if
        if (kind(pertzt).ne.ckind) then
            call ccpp_error('Kind mismatch for variable magnitude_of_perturbation_of_heat_to_momentum_roughness_length_ratio')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'magnitude_of_perturbation_of_soil_type_b_parameter', pertshc, ierr=ierr, dims=cdims, kind=ckind, index=499)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve magnitude_of_perturbation_of_soil_type_b_parameter from CCPP data structure')
            return
        end if
        if (kind(pertshc).ne.ckind) then
            call ccpp_error('Kind mismatch for variable magnitude_of_perturbation_of_soil_type_b_parameter')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'magnitude_of_perturbation_of_leaf_area_index', pertlai, ierr=ierr, dims=cdims, kind=ckind, index=497)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve magnitude_of_perturbation_of_leaf_area_index from CCPP data structure')
            return
        end if
        if (kind(pertlai).ne.ckind) then
            call ccpp_error('Kind mismatch for variable magnitude_of_perturbation_of_leaf_area_index')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'magnitude_of_perturbation_of_vegetation_fraction', pertvegf, ierr=ierr, dims=cdims, kind=ckind, index=500)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve magnitude_of_perturbation_of_vegetation_fraction from CCPP data structure')
            return
        end if
        if (kind(pertvegf).ne.ckind) then
            call ccpp_error('Kind mismatch for variable magnitude_of_perturbation_of_vegetation_fraction')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'perturbation_of_momentum_roughness_length', z01d, ierr=ierr, dims=cdims, kind=ckind, index=603)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve perturbation_of_momentum_roughness_length from CCPP data structure')
            return
        end if
        if (kind(z01d).ne.ckind) then
            call ccpp_error('Kind mismatch for variable perturbation_of_momentum_roughness_length')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'perturbation_of_heat_to_momentum_roughness_length_ratio', zt1d, ierr=ierr, dims=cdims, kind=ckind, index=601)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve perturbation_of_heat_to_momentum_roughness_length_ratio from CCPP data structure')
            return
        end if
        if (kind(zt1d).ne.ckind) then
            call ccpp_error('Kind mismatch for variable perturbation_of_heat_to_momentum_roughness_length_ratio')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'perturbation_of_soil_type_b_parameter', bexp1d, ierr=ierr, dims=cdims, kind=ckind, index=604)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve perturbation_of_soil_type_b_parameter from CCPP data structure')
            return
        end if
        if (kind(bexp1d).ne.ckind) then
            call ccpp_error('Kind mismatch for variable perturbation_of_soil_type_b_parameter')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'perturbation_of_leaf_area_index', xlai1d, ierr=ierr, dims=cdims, kind=ckind, index=602)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve perturbation_of_leaf_area_index from CCPP data structure')
            return
        end if
        if (kind(xlai1d).ne.ckind) then
            call ccpp_error('Kind mismatch for variable perturbation_of_leaf_area_index')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'perturbation_of_vegetation_fraction', vegf1d, ierr=ierr, dims=cdims, kind=ckind, index=605)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve perturbation_of_vegetation_fraction from CCPP data structure')
            return
        end if
        if (kind(vegf1d).ne.ckind) then
            call ccpp_error('Kind mismatch for variable perturbation_of_vegetation_fraction')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call GFS_surface_generic_pre_run(im=im,levs=levs,vfrac=vfrac,islmsk=islmsk,isot=isot,ivegsrc=ivegsrc,stype=stype, &
                  vtype=vtype,slope=slope,prsik_1=prsik_1,prslk_1=prslk_1,semis=semis,adjsfcdlw=adjsfcdlw, &
                  tsfc=tsfc,phil=phil,con_g=con_g,sigmaf=sigmaf,soiltyp=soiltyp,vegtype=vegtype, &
                  slopetyp=slopetyp,work3=work3,gabsbdlw=gabsbdlw,tsurf=tsurf,zlvl=zlvl,do_sppt=do_sppt, &
                  dtdtr=dtdtr,drain_cpl=drain_cpl,dsnow_cpl=dsnow_cpl,rain_cpl=rain_cpl,snow_cpl=snow_cpl, &
                  do_sfcperts=do_sfcperts,nsfcpert=nsfcpert,sfc_wts=sfc_wts,pertz0=pertz0, &
                  pertzt=pertzt,pertshc=pertshc,pertlai=pertlai,pertvegf=pertvegf,z01d=z01d, &
                  zt1d=zt1d,bexp1d=bexp1d,xlai1d=xlai1d,vegf1d=vegf1d,errmsg=cdata%errmsg, &
                  errflg=cdata%errflg)
        ierr=cdata%errflg

    end function GFS_surface_generic_pre_run_cap

    function GFS_surface_generic_pre_finalize_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call GFS_surface_generic_pre_finalize()
        

    end function GFS_surface_generic_pre_finalize_cap
end module GFS_surface_generic_pre_cap

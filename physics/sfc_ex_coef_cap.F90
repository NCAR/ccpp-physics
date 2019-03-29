
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
!! @brief Auto-generated cap module for the sfc_ex_coef scheme
!!
!
module sfc_ex_coef_cap

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_f_pointer, c_ptr, c_int32_t
    use            :: ccpp_types,                                      &
                      only: ccpp_t, CCPP_GENERIC_KIND
    use            :: ccpp_fields,                                     &
                      only: ccpp_field_get
    use            :: ccpp_errors,                                     &
                      only: ccpp_error, ccpp_debug
    use            :: sfc_ex_coef, &
                      only: sfc_ex_coef_init,sfc_ex_coef_run,sfc_ex_coef_finalize
    ! Other modules required, e.g. type definitions
    use machine, only: kind_phys

    implicit none

    private
    public :: sfc_ex_coef_init_cap,sfc_ex_coef_run_cap,sfc_ex_coef_finalize_cap

    contains


    function sfc_ex_coef_init_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call sfc_ex_coef_init()
        

    end function sfc_ex_coef_init_cap

    function sfc_ex_coef_run_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind
        integer, pointer :: im
        real(kind_phys), pointer :: ps(:)
        real(kind_phys), pointer :: u1(:)
        real(kind_phys), pointer :: v1(:)
        real(kind_phys), pointer :: t1(:)
        real(kind_phys), pointer :: q1(:)
        real(kind_phys), pointer :: z1(:)
        real(kind_phys), pointer :: snwdph(:)
        real(kind_phys), pointer :: tskin(:)
        real(kind_phys), pointer :: z0rl(:)
        real(kind_phys), pointer :: cm(:)
        real(kind_phys), pointer :: ch(:)
        real(kind_phys), pointer :: rb(:)
        real(kind_phys), pointer :: prsl1(:)
        real(kind_phys), pointer :: prslki(:)
        integer, pointer :: islimsk(:)
        real(kind_phys), pointer :: stress(:)
        real(kind_phys), pointer :: fm(:)
        real(kind_phys), pointer :: fh(:)
        real(kind_phys), pointer :: ustar(:)
        real(kind_phys), pointer :: wind(:)
        real(kind_phys), pointer :: ddvel(:)
        real(kind_phys), pointer :: fm10(:)
        real(kind_phys), pointer :: fh2(:)
        real(kind_phys), pointer :: sigmaf(:)
        integer, pointer :: vegtype(:)
        real(kind_phys), pointer :: shdmax(:)
        integer, pointer :: ivegsrc
        real(kind_phys), pointer :: z0pert(:)
        real(kind_phys), pointer :: ztpert(:)
        real(kind_phys), pointer :: tsurf(:)
        logical, pointer :: flag_iter(:)
        logical, pointer :: redrag

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
        

        call ccpp_field_get(cdata, 'surface_air_pressure', ps, ierr=ierr, dims=cdims, kind=ckind, index=677)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_air_pressure from CCPP data structure')
            return
        end if
        if (kind(ps).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_air_pressure')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'x_wind_at_lowest_model_layer', u1, ierr=ierr, dims=cdims, kind=ckind, index=873)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve x_wind_at_lowest_model_layer from CCPP data structure')
            return
        end if
        if (kind(u1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable x_wind_at_lowest_model_layer')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'y_wind_at_lowest_model_layer', v1, ierr=ierr, dims=cdims, kind=ckind, index=880)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve y_wind_at_lowest_model_layer from CCPP data structure')
            return
        end if
        if (kind(v1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable y_wind_at_lowest_model_layer')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'air_temperature_at_lowest_model_layer', t1, ierr=ierr, dims=cdims, kind=ckind, index=53)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve air_temperature_at_lowest_model_layer from CCPP data structure')
            return
        end if
        if (kind(t1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable air_temperature_at_lowest_model_layer')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'water_vapor_specific_humidity_at_lowest_model_layer', q1, ierr=ierr, dims=cdims, kind=ckind, index=854)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve water_vapor_specific_humidity_at_lowest_model_layer from CCPP data structure')
            return
        end if
        if (kind(q1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable water_vapor_specific_humidity_at_lowest_model_layer')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'height_above_ground_at_lowest_model_layer', z1, ierr=ierr, dims=cdims, kind=ckind, index=360)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve height_above_ground_at_lowest_model_layer from CCPP data structure')
            return
        end if
        if (kind(z1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable height_above_ground_at_lowest_model_layer')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_snow_thickness_water_equivalent', snwdph, ierr=ierr, dims=cdims, kind=ckind, index=729)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_snow_thickness_water_equivalent from CCPP data structure')
            return
        end if
        if (kind(snwdph).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_snow_thickness_water_equivalent')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_skin_temperature', tskin, ierr=ierr, dims=cdims, kind=ckind, index=721)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_skin_temperature from CCPP data structure')
            return
        end if
        if (kind(tskin).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_skin_temperature')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_roughness_length', z0rl, ierr=ierr, dims=cdims, kind=ckind, index=718)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_roughness_length from CCPP data structure')
            return
        end if
        if (kind(z0rl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_roughness_length')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_drag_coefficient_for_momentum_in_air', cm, ierr=ierr, dims=cdims, kind=ckind, index=703)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_drag_coefficient_for_momentum_in_air from CCPP data structure')
            return
        end if
        if (kind(cm).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_drag_coefficient_for_momentum_in_air')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_drag_coefficient_for_heat_and_moisture_in_air', ch, ierr=ierr, dims=cdims, kind=ckind, index=702)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_drag_coefficient_for_heat_and_moisture_in_air from CCPP data structure')
            return
        end if
        if (kind(ch).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_drag_coefficient_for_heat_and_moisture_in_air')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'bulk_richardson_number_at_lowest_model_level', rb, ierr=ierr, dims=cdims, kind=ckind, index=78)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve bulk_richardson_number_at_lowest_model_level from CCPP data structure')
            return
        end if
        if (kind(rb).ne.ckind) then
            call ccpp_error('Kind mismatch for variable bulk_richardson_number_at_lowest_model_level')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'air_pressure_at_lowest_model_layer', prsl1, ierr=ierr, dims=cdims, kind=ckind, index=48)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve air_pressure_at_lowest_model_layer from CCPP data structure')
            return
        end if
        if (kind(prsl1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable air_pressure_at_lowest_model_layer')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'ratio_of_exner_function_between_midlayer_and_interface_at_lowest_model_layer', prslki, ierr=ierr, dims=cdims, kind=ckind, index=623)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve ratio_of_exner_function_between_midlayer_and_interface_at_lowest_model_layer from CCPP data structure')
            return
        end if
        if (kind(prslki).ne.ckind) then
            call ccpp_error('Kind mismatch for variable ratio_of_exner_function_between_midlayer_and_interface_at_lowest_model_layer')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'sea_land_ice_mask', islimsk, ierr=ierr, dims=cdims, kind=ckind, index=631)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve sea_land_ice_mask from CCPP data structure')
            return
        end if
        if (kind(islimsk).ne.ckind) then
            call ccpp_error('Kind mismatch for variable sea_land_ice_mask')
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
        

        call ccpp_field_get(cdata, 'surface_friction_velocity', ustar, ierr=ierr, dims=cdims, kind=ckind, index=710)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_friction_velocity from CCPP data structure')
            return
        end if
        if (kind(ustar).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_friction_velocity')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'wind_speed_at_lowest_model_layer', wind, ierr=ierr, dims=cdims, kind=ckind, index=870)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve wind_speed_at_lowest_model_layer from CCPP data structure')
            return
        end if
        if (kind(wind).ne.ckind) then
            call ccpp_error('Kind mismatch for variable wind_speed_at_lowest_model_layer')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_wind_enhancement_due_to_convection', ddvel, ierr=ierr, dims=cdims, kind=ckind, index=744)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_wind_enhancement_due_to_convection from CCPP data structure')
            return
        end if
        if (kind(ddvel).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_wind_enhancement_due_to_convection')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'Monin-Obukhov_similarity_function_for_momentum_at_10m', fm10, ierr=ierr, dims=cdims, kind=ckind, index=19)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve Monin-Obukhov_similarity_function_for_momentum_at_10m from CCPP data structure')
            return
        end if
        if (kind(fm10).ne.ckind) then
            call ccpp_error('Kind mismatch for variable Monin-Obukhov_similarity_function_for_momentum_at_10m')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'Monin-Obukhov_similarity_function_for_heat_at_2m', fh2, ierr=ierr, dims=cdims, kind=ckind, index=17)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve Monin-Obukhov_similarity_function_for_heat_at_2m from CCPP data structure')
            return
        end if
        if (kind(fh2).ne.ckind) then
            call ccpp_error('Kind mismatch for variable Monin-Obukhov_similarity_function_for_heat_at_2m')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

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
        

        call ccpp_field_get(cdata, 'maximum_vegetation_area_fraction', shdmax, ierr=ierr, dims=cdims, kind=ckind, index=510)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve maximum_vegetation_area_fraction from CCPP data structure')
            return
        end if
        if (kind(shdmax).ne.ckind) then
            call ccpp_error('Kind mismatch for variable maximum_vegetation_area_fraction')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

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
        

        call ccpp_field_get(cdata, 'perturbation_of_momentum_roughness_length', z0pert, ierr=ierr, dims=cdims, kind=ckind, index=603)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve perturbation_of_momentum_roughness_length from CCPP data structure')
            return
        end if
        if (kind(z0pert).ne.ckind) then
            call ccpp_error('Kind mismatch for variable perturbation_of_momentum_roughness_length')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'perturbation_of_heat_to_momentum_roughness_length_ratio', ztpert, ierr=ierr, dims=cdims, kind=ckind, index=601)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve perturbation_of_heat_to_momentum_roughness_length_ratio from CCPP data structure')
            return
        end if
        if (kind(ztpert).ne.ckind) then
            call ccpp_error('Kind mismatch for variable perturbation_of_heat_to_momentum_roughness_length_ratio')
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
        

        call ccpp_field_get(cdata, 'flag_for_iteration', flag_iter, ierr=ierr, dims=cdims, kind=ckind, index=287)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_iteration from CCPP data structure')
            return
        end if
        if (kind(flag_iter).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_iteration')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'flag_for_reduced_drag_coefficient_over_sea', redrag, ierr=ierr, kind=ckind, index=308)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_reduced_drag_coefficient_over_sea from CCPP data structure')
            return
        end if
        if (kind(redrag).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_reduced_drag_coefficient_over_sea')
            ierr = 1
            return
        end if
#endif
        

        call sfc_ex_coef_run(im=im,ps=ps,u1=u1,v1=v1,t1=t1,q1=q1,z1=z1,snwdph=snwdph,tskin=tskin,z0rl=z0rl, &
                  cm=cm,ch=ch,rb=rb,prsl1=prsl1,prslki=prslki,islimsk=islimsk,stress=stress, &
                  fm=fm,fh=fh,ustar=ustar,wind=wind,ddvel=ddvel,fm10=fm10,fh2=fh2,sigmaf=sigmaf, &
                  vegtype=vegtype,shdmax=shdmax,ivegsrc=ivegsrc,z0pert=z0pert,ztpert=ztpert, &
                  tsurf=tsurf,flag_iter=flag_iter,redrag=redrag,errmsg=cdata%errmsg,errflg=cdata%errflg)
        ierr=cdata%errflg

    end function sfc_ex_coef_run_cap

    function sfc_ex_coef_finalize_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call sfc_ex_coef_finalize()
        

    end function sfc_ex_coef_finalize_cap
end module sfc_ex_coef_cap

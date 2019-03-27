
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
!! @brief Auto-generated cap module for the mynnsfc_wrapper scheme
!!
!
module mynnsfc_wrapper_cap

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_f_pointer, c_ptr, c_int32_t
    use            :: ccpp_types,                                      &
                      only: ccpp_t, CCPP_GENERIC_KIND
    use            :: ccpp_fields,                                     &
                      only: ccpp_field_get
    use            :: ccpp_errors,                                     &
                      only: ccpp_error, ccpp_debug
    use            :: mynnsfc_wrapper, &
                      only: mynnsfc_wrapper_finalize,mynnsfc_wrapper_init,mynnsfc_wrapper_run
    ! Other modules required, e.g. type definitions
    use machine, only: kind_phys

    implicit none

    private
    public :: mynnsfc_wrapper_finalize_cap,mynnsfc_wrapper_init_cap,mynnsfc_wrapper_run_cap

    contains


    function mynnsfc_wrapper_finalize_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call mynnsfc_wrapper_finalize()
        

    end function mynnsfc_wrapper_finalize_cap

    function mynnsfc_wrapper_init_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call mynnsfc_wrapper_init()
        

    end function mynnsfc_wrapper_init_cap

    function mynnsfc_wrapper_run_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind
        integer, pointer :: ix
        integer, pointer :: im
        integer, pointer :: levs
        logical, pointer :: flag_init
        logical, pointer :: flag_restart
        real(kind_phys), pointer :: delt
        real(kind_phys), pointer :: dx(:)
        real(kind_phys), pointer :: u(:,:)
        real(kind_phys), pointer :: v(:,:)
        real(kind_phys), pointer :: t3d(:,:)
        real(kind_phys), pointer :: qvsh(:,:)
        real(kind_phys), pointer :: qc(:,:)
        real(kind_phys), pointer :: prsl(:,:)
        real(kind_phys), pointer :: phii(:,:)
        real(kind_phys), pointer :: exner(:,:)
        real(kind_phys), pointer :: tsq(:,:)
        real(kind_phys), pointer :: qsq(:,:)
        real(kind_phys), pointer :: cov(:,:)
        real(kind_phys), pointer :: el_pbl(:,:)
        real(kind_phys), pointer :: Sh3D(:,:)
        real(kind_phys), pointer :: QC_BL(:,:)
        real(kind_phys), pointer :: CLDFRA_BL(:,:)
        real(kind_phys), pointer :: ps(:)
        real(kind_phys), pointer :: PBLH(:)
        real(kind_phys), pointer :: slmsk(:)
        real(kind_phys), pointer :: tsk(:)
        real(kind_phys), pointer :: qsfc(:)
        real(kind_phys), pointer :: snowd(:)
        real(kind_phys), pointer :: zorl(:)
        real(kind_phys), pointer :: ust(:)
        real(kind_phys), pointer :: ustm(:)
        real(kind_phys), pointer :: zol(:)
        real(kind_phys), pointer :: mol(:)
        real(kind_phys), pointer :: rmol(:)
        real(kind_phys), pointer :: fm(:)
        real(kind_phys), pointer :: fh(:)
        real(kind_phys), pointer :: fm10(:)
        real(kind_phys), pointer :: fh2(:)
        real(kind_phys), pointer :: wspd(:)
        real(kind_phys), pointer :: br(:)
        real(kind_phys), pointer :: ch(:)
        real(kind_phys), pointer :: hflx(:)
        real(kind_phys), pointer :: QFX(:)
        real(kind_phys), pointer :: lh(:)
        real(kind_phys), pointer :: flhc(:)
        real(kind_phys), pointer :: flqc(:)
        real(kind_phys), pointer :: u10(:)
        real(kind_phys), pointer :: v10(:)
        real(kind_phys), pointer :: th2(:)
        real(kind_phys), pointer :: t2(:)
        real(kind_phys), pointer :: q2(:)
        real(kind_phys), pointer :: wstar(:)
        real(kind_phys), pointer :: chs2(:)
        real(kind_phys), pointer :: cqs2(:)
        real(kind_phys), pointer :: cda(:)
        real(kind_phys), pointer :: cka(:)
        real(kind_phys), pointer :: stress(:)
        integer, pointer :: bl_mynn_cloudpdf
        integer, pointer :: icloud_bl
        logical, pointer :: lprnt

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
        

        call ccpp_field_get(cdata, 'flag_for_first_time_step', flag_init, ierr=ierr, kind=ckind, index=277)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_first_time_step from CCPP data structure')
            return
        end if
        if (kind(flag_init).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_first_time_step')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'flag_for_restart', flag_restart, ierr=ierr, kind=ckind, index=309)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_restart from CCPP data structure')
            return
        end if
        if (kind(flag_restart).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_restart')
            ierr = 1
            return
        end if
#endif
        

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
        

        call ccpp_field_get(cdata, 'cell_size', dx, ierr=ierr, dims=cdims, kind=ckind, index=84)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cell_size from CCPP data structure')
            return
        end if
        if (kind(dx).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cell_size')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'x_wind', u, ierr=ierr, dims=cdims, kind=ckind, index=871)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve x_wind from CCPP data structure')
            return
        end if
        if (kind(u).ne.ckind) then
            call ccpp_error('Kind mismatch for variable x_wind')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'y_wind', v, ierr=ierr, dims=cdims, kind=ckind, index=878)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve y_wind from CCPP data structure')
            return
        end if
        if (kind(v).ne.ckind) then
            call ccpp_error('Kind mismatch for variable y_wind')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'air_temperature', t3d, ierr=ierr, dims=cdims, kind=ckind, index=50)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve air_temperature from CCPP data structure')
            return
        end if
        if (kind(t3d).ne.ckind) then
            call ccpp_error('Kind mismatch for variable air_temperature')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'water_vapor_specific_humidity', qvsh, ierr=ierr, dims=cdims, kind=ckind, index=852)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve water_vapor_specific_humidity from CCPP data structure')
            return
        end if
        if (kind(qvsh).ne.ckind) then
            call ccpp_error('Kind mismatch for variable water_vapor_specific_humidity')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cloud_condensed_water_mixing_ratio', qc, ierr=ierr, dims=cdims, kind=ckind, index=90)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cloud_condensed_water_mixing_ratio from CCPP data structure')
            return
        end if
        if (kind(qc).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cloud_condensed_water_mixing_ratio')
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
        

        call ccpp_field_get(cdata, 'dimensionless_exner_function_at_model_layers', exner, ierr=ierr, dims=cdims, kind=ckind, index=222)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve dimensionless_exner_function_at_model_layers from CCPP data structure')
            return
        end if
        if (kind(exner).ne.ckind) then
            call ccpp_error('Kind mismatch for variable dimensionless_exner_function_at_model_layers')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 't_prime_squared', tsq, ierr=ierr, dims=cdims, kind=ckind, index=749)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve t_prime_squared from CCPP data structure')
            return
        end if
        if (kind(tsq).ne.ckind) then
            call ccpp_error('Kind mismatch for variable t_prime_squared')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'q_prime_squared', qsq, ierr=ierr, dims=cdims, kind=ckind, index=612)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve q_prime_squared from CCPP data structure')
            return
        end if
        if (kind(qsq).ne.ckind) then
            call ccpp_error('Kind mismatch for variable q_prime_squared')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 't_prime_q_prime', cov, ierr=ierr, dims=cdims, kind=ckind, index=748)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve t_prime_q_prime from CCPP data structure')
            return
        end if
        if (kind(cov).ne.ckind) then
            call ccpp_error('Kind mismatch for variable t_prime_q_prime')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'mixing_length', el_pbl, ierr=ierr, dims=cdims, kind=ckind, index=549)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve mixing_length from CCPP data structure')
            return
        end if
        if (kind(el_pbl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable mixing_length')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'stability_function_for_heat', Sh3D, ierr=ierr, dims=cdims, kind=ckind, index=668)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve stability_function_for_heat from CCPP data structure')
            return
        end if
        if (kind(Sh3D).ne.ckind) then
            call ccpp_error('Kind mismatch for variable stability_function_for_heat')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'subgrid_cloud_mixing_ratio_pbl', QC_BL, ierr=ierr, dims=cdims, kind=ckind, index=674)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve subgrid_cloud_mixing_ratio_pbl from CCPP data structure')
            return
        end if
        if (kind(QC_BL).ne.ckind) then
            call ccpp_error('Kind mismatch for variable subgrid_cloud_mixing_ratio_pbl')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'subgrid_cloud_fraction_pbl', CLDFRA_BL, ierr=ierr, dims=cdims, kind=ckind, index=673)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve subgrid_cloud_fraction_pbl from CCPP data structure')
            return
        end if
        if (kind(CLDFRA_BL).ne.ckind) then
            call ccpp_error('Kind mismatch for variable subgrid_cloud_fraction_pbl')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

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
        

        call ccpp_field_get(cdata, 'atmosphere_boundary_layer_thickness', PBLH, ierr=ierr, dims=cdims, kind=ckind, index=66)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve atmosphere_boundary_layer_thickness from CCPP data structure')
            return
        end if
        if (kind(PBLH).ne.ckind) then
            call ccpp_error('Kind mismatch for variable atmosphere_boundary_layer_thickness')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'sea_land_ice_mask_real', slmsk, ierr=ierr, dims=cdims, kind=ckind, index=632)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve sea_land_ice_mask_real from CCPP data structure')
            return
        end if
        if (kind(slmsk).ne.ckind) then
            call ccpp_error('Kind mismatch for variable sea_land_ice_mask_real')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_skin_temperature', tsk, ierr=ierr, dims=cdims, kind=ckind, index=721)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_skin_temperature from CCPP data structure')
            return
        end if
        if (kind(tsk).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_skin_temperature')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_specific_humidity', qsfc, ierr=ierr, dims=cdims, kind=ckind, index=730)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_specific_humidity from CCPP data structure')
            return
        end if
        if (kind(qsfc).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_specific_humidity')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_snow_thickness_water_equivalent', snowd, ierr=ierr, dims=cdims, kind=ckind, index=729)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_snow_thickness_water_equivalent from CCPP data structure')
            return
        end if
        if (kind(snowd).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_snow_thickness_water_equivalent')
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
        

        call ccpp_field_get(cdata, 'surface_friction_velocity', ust, ierr=ierr, dims=cdims, kind=ckind, index=710)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_friction_velocity from CCPP data structure')
            return
        end if
        if (kind(ust).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_friction_velocity')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_friction_velocity_drag', ustm, ierr=ierr, dims=cdims, kind=ckind, index=711)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_friction_velocity_drag from CCPP data structure')
            return
        end if
        if (kind(ustm).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_friction_velocity_drag')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_stability_parameter', zol, ierr=ierr, dims=cdims, kind=ckind, index=731)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_stability_parameter from CCPP data structure')
            return
        end if
        if (kind(zol).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_stability_parameter')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'theta_star', mol, ierr=ierr, dims=cdims, kind=ckind, index=787)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve theta_star from CCPP data structure')
            return
        end if
        if (kind(mol).ne.ckind) then
            call ccpp_error('Kind mismatch for variable theta_star')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'reciprocal_of_obukhov_length', rmol, ierr=ierr, dims=cdims, kind=ckind, index=627)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve reciprocal_of_obukhov_length from CCPP data structure')
            return
        end if
        if (kind(rmol).ne.ckind) then
            call ccpp_error('Kind mismatch for variable reciprocal_of_obukhov_length')
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
        

        call ccpp_field_get(cdata, 'wind_speed_at_lowest_model_layer', wspd, ierr=ierr, dims=cdims, kind=ckind, index=870)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve wind_speed_at_lowest_model_layer from CCPP data structure')
            return
        end if
        if (kind(wspd).ne.ckind) then
            call ccpp_error('Kind mismatch for variable wind_speed_at_lowest_model_layer')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'bulk_richardson_number_at_lowest_model_level', br, ierr=ierr, dims=cdims, kind=ckind, index=78)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve bulk_richardson_number_at_lowest_model_level from CCPP data structure')
            return
        end if
        if (kind(br).ne.ckind) then
            call ccpp_error('Kind mismatch for variable bulk_richardson_number_at_lowest_model_level')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_drag_wind_speed_for_momentum_in_air', ch, ierr=ierr, dims=cdims, kind=ckind, index=705)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_drag_wind_speed_for_momentum_in_air from CCPP data structure')
            return
        end if
        if (kind(ch).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_drag_wind_speed_for_momentum_in_air')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'kinematic_surface_upward_sensible_heat_flux', hflx, ierr=ierr, dims=cdims, kind=ckind, index=455)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve kinematic_surface_upward_sensible_heat_flux from CCPP data structure')
            return
        end if
        if (kind(hflx).ne.ckind) then
            call ccpp_error('Kind mismatch for variable kinematic_surface_upward_sensible_heat_flux')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'kinematic_surface_upward_latent_heat_flux', QFX, ierr=ierr, dims=cdims, kind=ckind, index=454)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve kinematic_surface_upward_latent_heat_flux from CCPP data structure')
            return
        end if
        if (kind(QFX).ne.ckind) then
            call ccpp_error('Kind mismatch for variable kinematic_surface_upward_latent_heat_flux')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_latent_heat', lh, ierr=ierr, dims=cdims, kind=ckind, index=713)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_latent_heat from CCPP data structure')
            return
        end if
        if (kind(lh).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_latent_heat')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_exchange_coefficient_for_heat', flhc, ierr=ierr, dims=cdims, kind=ckind, index=706)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_exchange_coefficient_for_heat from CCPP data structure')
            return
        end if
        if (kind(flhc).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_exchange_coefficient_for_heat')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_exchange_coefficient_for_moisture', flqc, ierr=ierr, dims=cdims, kind=ckind, index=708)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_exchange_coefficient_for_moisture from CCPP data structure')
            return
        end if
        if (kind(flqc).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_exchange_coefficient_for_moisture')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'x_wind_at_10m', u10, ierr=ierr, dims=cdims, kind=ckind, index=872)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve x_wind_at_10m from CCPP data structure')
            return
        end if
        if (kind(u10).ne.ckind) then
            call ccpp_error('Kind mismatch for variable x_wind_at_10m')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'y_wind_at_10m', v10, ierr=ierr, dims=cdims, kind=ckind, index=879)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve y_wind_at_10m from CCPP data structure')
            return
        end if
        if (kind(v10).ne.ckind) then
            call ccpp_error('Kind mismatch for variable y_wind_at_10m')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'potential_temperature_at_2m', th2, ierr=ierr, dims=cdims, kind=ckind, index=607)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve potential_temperature_at_2m from CCPP data structure')
            return
        end if
        if (kind(th2).ne.ckind) then
            call ccpp_error('Kind mismatch for variable potential_temperature_at_2m')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'temperature_at_2m', t2, ierr=ierr, dims=cdims, kind=ckind, index=750)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve temperature_at_2m from CCPP data structure')
            return
        end if
        if (kind(t2).ne.ckind) then
            call ccpp_error('Kind mismatch for variable temperature_at_2m')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'specific_humidity_at_2m', q2, ierr=ierr, dims=cdims, kind=ckind, index=667)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve specific_humidity_at_2m from CCPP data structure')
            return
        end if
        if (kind(q2).ne.ckind) then
            call ccpp_error('Kind mismatch for variable specific_humidity_at_2m')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_wind_enhancement_due_to_convection', wstar, ierr=ierr, dims=cdims, kind=ckind, index=744)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_wind_enhancement_due_to_convection from CCPP data structure')
            return
        end if
        if (kind(wstar).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_wind_enhancement_due_to_convection')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_exchange_coefficient_for_heat_at_2m', chs2, ierr=ierr, dims=cdims, kind=ckind, index=707)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_exchange_coefficient_for_heat_at_2m from CCPP data structure')
            return
        end if
        if (kind(chs2).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_exchange_coefficient_for_heat_at_2m')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_exchange_coefficient_for_moisture_at_2m', cqs2, ierr=ierr, dims=cdims, kind=ckind, index=709)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_exchange_coefficient_for_moisture_at_2m from CCPP data structure')
            return
        end if
        if (kind(cqs2).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_exchange_coefficient_for_moisture_at_2m')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_drag_coefficient_for_momentum_in_air', cda, ierr=ierr, dims=cdims, kind=ckind, index=703)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_drag_coefficient_for_momentum_in_air from CCPP data structure')
            return
        end if
        if (kind(cda).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_drag_coefficient_for_momentum_in_air')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_drag_coefficient_for_heat_and_moisture_in_air', cka, ierr=ierr, dims=cdims, kind=ckind, index=702)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_drag_coefficient_for_heat_and_moisture_in_air from CCPP data structure')
            return
        end if
        if (kind(cka).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_drag_coefficient_for_heat_and_moisture_in_air')
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
        

        call ccpp_field_get(cdata, 'cloudpdf', bl_mynn_cloudpdf, ierr=ierr, kind=ckind, index=112)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cloudpdf from CCPP data structure')
            return
        end if
        if (kind(bl_mynn_cloudpdf).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cloudpdf')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'couple_sgs_clouds_to_radiation_flag', icloud_bl, ierr=ierr, kind=ckind, index=139)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve couple_sgs_clouds_to_radiation_flag from CCPP data structure')
            return
        end if
        if (kind(icloud_bl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable couple_sgs_clouds_to_radiation_flag')
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
        

        call mynnsfc_wrapper_run(ix=ix,im=im,levs=levs,iter=cdata%loop_cnt,flag_init=flag_init,flag_restart=flag_restart, &
                  delt=delt,dx=dx,u=u,v=v,t3d=t3d,qvsh=qvsh,qc=qc,prsl=prsl,phii=phii,exner=exner, &
                  tsq=tsq,qsq=qsq,cov=cov,el_pbl=el_pbl,Sh3D=Sh3D,QC_BL=QC_BL,CLDFRA_BL=CLDFRA_BL, &
                  ps=ps,PBLH=PBLH,slmsk=slmsk,tsk=tsk,qsfc=qsfc,snowd=snowd,zorl=zorl,ust=ust, &
                  ustm=ustm,zol=zol,mol=mol,rmol=rmol,fm=fm,fh=fh,fm10=fm10,fh2=fh2,wspd=wspd, &
                  br=br,ch=ch,hflx=hflx,QFX=QFX,lh=lh,flhc=flhc,flqc=flqc,u10=u10,v10=v10, &
                  th2=th2,t2=t2,q2=q2,wstar=wstar,chs2=chs2,cqs2=cqs2,cda=cda,cka=cka,stress=stress, &
                  bl_mynn_cloudpdf=bl_mynn_cloudpdf,icloud_bl=icloud_bl,lprnt=lprnt,errmsg=cdata%errmsg, &
                  errflg=cdata%errflg)
        ierr=cdata%errflg

    end function mynnsfc_wrapper_run_cap
end module mynnsfc_wrapper_cap

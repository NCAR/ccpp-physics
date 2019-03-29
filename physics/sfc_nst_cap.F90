
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
!! @brief Auto-generated cap module for the sfc_nst scheme
!!
!
module sfc_nst_cap

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_f_pointer, c_ptr, c_int32_t
    use            :: ccpp_types,                                      &
                      only: ccpp_t, CCPP_GENERIC_KIND
    use            :: ccpp_fields,                                     &
                      only: ccpp_field_get
    use            :: ccpp_errors,                                     &
                      only: ccpp_error, ccpp_debug
    use            :: sfc_nst, &
                      only: sfc_nst_init,sfc_nst_finalize,sfc_nst_run
    ! Other modules required, e.g. type definitions
    use machine, only: kind_phys

    implicit none

    private
    public :: sfc_nst_init_cap,sfc_nst_finalize_cap,sfc_nst_run_cap

    contains


    function sfc_nst_init_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call sfc_nst_init()
        

    end function sfc_nst_init_cap

    function sfc_nst_finalize_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call sfc_nst_finalize()
        

    end function sfc_nst_finalize_cap

    function sfc_nst_run_cap(ptr) bind(c) result(ierr)

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
        real(kind_phys), pointer :: tref(:)
        real(kind_phys), pointer :: cm(:)
        real(kind_phys), pointer :: ch(:)
        real(kind_phys), pointer :: prsl1(:)
        real(kind_phys), pointer :: prslki(:)
        integer, pointer :: islimsk(:)
        real(kind_phys), pointer :: xlon(:)
        real(kind_phys), pointer :: sinlat(:)
        real(kind_phys), pointer :: stress(:)
        real(kind_phys), pointer :: sfcemis(:)
        real(kind_phys), pointer :: dlwflx(:)
        real(kind_phys), pointer :: sfcnsw(:)
        real(kind_phys), pointer :: rain(:)
        real(kind_phys), pointer :: timestep
        integer, pointer :: kdt
        real(kind_phys), pointer :: solhr
        real(kind_phys), pointer :: xcosz(:)
        real(kind_phys), pointer :: ddvel(:)
        logical, pointer :: flag_iter(:)
        logical, pointer :: flag_guess(:)
        integer, pointer :: nstf_name1
        integer, pointer :: nstf_name4
        integer, pointer :: nstf_name5
        logical, pointer :: lprnt
        integer, pointer :: ipr
        real(kind_phys), pointer :: tskin(:)
        real(kind_phys), pointer :: tsurf(:)
        real(kind_phys), pointer :: xt(:)
        real(kind_phys), pointer :: xs(:)
        real(kind_phys), pointer :: xu(:)
        real(kind_phys), pointer :: xv(:)
        real(kind_phys), pointer :: xz(:)
        real(kind_phys), pointer :: zm(:)
        real(kind_phys), pointer :: xtts(:)
        real(kind_phys), pointer :: xzts(:)
        real(kind_phys), pointer :: dt_cool(:)
        real(kind_phys), pointer :: z_c(:)
        real(kind_phys), pointer :: c_0(:)
        real(kind_phys), pointer :: c_d(:)
        real(kind_phys), pointer :: w_0(:)
        real(kind_phys), pointer :: w_d(:)
        real(kind_phys), pointer :: d_conv(:)
        real(kind_phys), pointer :: ifd(:)
        real(kind_phys), pointer :: qrain(:)
        real(kind_phys), pointer :: qsurf(:)
        real(kind_phys), pointer :: gflux(:)
        real(kind_phys), pointer :: cmm(:)
        real(kind_phys), pointer :: chh(:)
        real(kind_phys), pointer :: evap(:)
        real(kind_phys), pointer :: hflx(:)
        real(kind_phys), pointer :: ep(:)

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
        

        call ccpp_field_get(cdata, 'sea_surface_reference_temperature', tref, ierr=ierr, dims=cdims, kind=ckind, index=633)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve sea_surface_reference_temperature from CCPP data structure')
            return
        end if
        if (kind(tref).ne.ckind) then
            call ccpp_error('Kind mismatch for variable sea_surface_reference_temperature')
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
        

        call ccpp_field_get(cdata, 'sine_of_latitude', sinlat, ierr=ierr, dims=cdims, kind=ckind, index=645)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve sine_of_latitude from CCPP data structure')
            return
        end if
        if (kind(sinlat).ne.ckind) then
            call ccpp_error('Kind mismatch for variable sine_of_latitude')
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
        

        call ccpp_field_get(cdata, 'surface_longwave_emissivity', sfcemis, ierr=ierr, dims=cdims, kind=ckind, index=714)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_longwave_emissivity from CCPP data structure')
            return
        end if
        if (kind(sfcemis).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_longwave_emissivity')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_downwelling_longwave_flux_absorbed_by_ground', dlwflx, ierr=ierr, dims=cdims, kind=ckind, index=698)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_downwelling_longwave_flux_absorbed_by_ground from CCPP data structure')
            return
        end if
        if (kind(dlwflx).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_downwelling_longwave_flux_absorbed_by_ground')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_net_downwelling_shortwave_flux', sfcnsw, ierr=ierr, dims=cdims, kind=ckind, index=716)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_net_downwelling_shortwave_flux from CCPP data structure')
            return
        end if
        if (kind(sfcnsw).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_net_downwelling_shortwave_flux')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'nonnegative_lwe_thickness_of_precipitation_amount_on_dynamics_timestep', rain, ierr=ierr, dims=cdims, kind=ckind, index=567)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve nonnegative_lwe_thickness_of_precipitation_amount_on_dynamics_timestep from CCPP data structure')
            return
        end if
        if (kind(rain).ne.ckind) then
            call ccpp_error('Kind mismatch for variable nonnegative_lwe_thickness_of_precipitation_amount_on_dynamics_timestep')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'time_step_for_dynamics', timestep, ierr=ierr, kind=ckind, index=792)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve time_step_for_dynamics from CCPP data structure')
            return
        end if
        if (kind(timestep).ne.ckind) then
            call ccpp_error('Kind mismatch for variable time_step_for_dynamics')
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
        

        call ccpp_field_get(cdata, 'forecast_hour', solhr, ierr=ierr, kind=ckind, index=338)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve forecast_hour from CCPP data structure')
            return
        end if
        if (kind(solhr).ne.ckind) then
            call ccpp_error('Kind mismatch for variable forecast_hour')
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
        

        call ccpp_field_get(cdata, 'flag_for_guess_run', flag_guess, ierr=ierr, dims=cdims, kind=ckind, index=281)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_guess_run from CCPP data structure')
            return
        end if
        if (kind(flag_guess).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_guess_run')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'flag_for_nsstm_run', nstf_name1, ierr=ierr, kind=ckind, index=299)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_nsstm_run from CCPP data structure')
            return
        end if
        if (kind(nstf_name1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_nsstm_run')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'vertical_temperature_average_range_lower_bound', nstf_name4, ierr=ierr, kind=ckind, index=828)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve vertical_temperature_average_range_lower_bound from CCPP data structure')
            return
        end if
        if (kind(nstf_name4).ne.ckind) then
            call ccpp_error('Kind mismatch for variable vertical_temperature_average_range_lower_bound')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'vertical_temperature_average_range_upper_bound', nstf_name5, ierr=ierr, kind=ckind, index=829)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve vertical_temperature_average_range_upper_bound from CCPP data structure')
            return
        end if
        if (kind(nstf_name5).ne.ckind) then
            call ccpp_error('Kind mismatch for variable vertical_temperature_average_range_upper_bound')
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
        

        call ccpp_field_get(cdata, 'surface_skin_temperature_for_nsst', tskin, ierr=ierr, dims=cdims, kind=ckind, index=723)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_skin_temperature_for_nsst from CCPP data structure')
            return
        end if
        if (kind(tskin).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_skin_temperature_for_nsst')
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
        

        call ccpp_field_get(cdata, 'diurnal_thermocline_layer_heat_content', xt, ierr=ierr, dims=cdims, kind=ckind, index=224)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve diurnal_thermocline_layer_heat_content from CCPP data structure')
            return
        end if
        if (kind(xt).ne.ckind) then
            call ccpp_error('Kind mismatch for variable diurnal_thermocline_layer_heat_content')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'sea_water_salinity', xs, ierr=ierr, dims=cdims, kind=ckind, index=634)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve sea_water_salinity from CCPP data structure')
            return
        end if
        if (kind(xs).ne.ckind) then
            call ccpp_error('Kind mismatch for variable sea_water_salinity')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'diurnal_thermocline_layer_x_current', xu, ierr=ierr, dims=cdims, kind=ckind, index=226)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve diurnal_thermocline_layer_x_current from CCPP data structure')
            return
        end if
        if (kind(xu).ne.ckind) then
            call ccpp_error('Kind mismatch for variable diurnal_thermocline_layer_x_current')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'diurnal_thermocline_layer_y_current', xv, ierr=ierr, dims=cdims, kind=ckind, index=227)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve diurnal_thermocline_layer_y_current from CCPP data structure')
            return
        end if
        if (kind(xv).ne.ckind) then
            call ccpp_error('Kind mismatch for variable diurnal_thermocline_layer_y_current')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'diurnal_thermocline_layer_thickness', xz, ierr=ierr, dims=cdims, kind=ckind, index=225)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve diurnal_thermocline_layer_thickness from CCPP data structure')
            return
        end if
        if (kind(xz).ne.ckind) then
            call ccpp_error('Kind mismatch for variable diurnal_thermocline_layer_thickness')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'ocean_mixed_layer_thickness', zm, ierr=ierr, dims=cdims, kind=ckind, index=592)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve ocean_mixed_layer_thickness from CCPP data structure')
            return
        end if
        if (kind(zm).ne.ckind) then
            call ccpp_error('Kind mismatch for variable ocean_mixed_layer_thickness')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'sensitivity_of_dtl_heat_content_to_surface_temperature', xtts, ierr=ierr, dims=cdims, kind=ckind, index=638)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve sensitivity_of_dtl_heat_content_to_surface_temperature from CCPP data structure')
            return
        end if
        if (kind(xtts).ne.ckind) then
            call ccpp_error('Kind mismatch for variable sensitivity_of_dtl_heat_content_to_surface_temperature')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'sensitivity_of_dtl_thickness_to_surface_temperature', xzts, ierr=ierr, dims=cdims, kind=ckind, index=639)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve sensitivity_of_dtl_thickness_to_surface_temperature from CCPP data structure')
            return
        end if
        if (kind(xzts).ne.ckind) then
            call ccpp_error('Kind mismatch for variable sensitivity_of_dtl_thickness_to_surface_temperature')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'sub-layer_cooling_amount', dt_cool, ierr=ierr, dims=cdims, kind=ckind, index=671)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve sub-layer_cooling_amount from CCPP data structure')
            return
        end if
        if (kind(dt_cool).ne.ckind) then
            call ccpp_error('Kind mismatch for variable sub-layer_cooling_amount')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'sub-layer_cooling_thickness', z_c, ierr=ierr, dims=cdims, kind=ckind, index=672)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve sub-layer_cooling_thickness from CCPP data structure')
            return
        end if
        if (kind(z_c).ne.ckind) then
            call ccpp_error('Kind mismatch for variable sub-layer_cooling_thickness')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'coefficient_c_0', c_0, ierr=ierr, dims=cdims, kind=ckind, index=113)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve coefficient_c_0 from CCPP data structure')
            return
        end if
        if (kind(c_0).ne.ckind) then
            call ccpp_error('Kind mismatch for variable coefficient_c_0')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'coefficient_c_d', c_d, ierr=ierr, dims=cdims, kind=ckind, index=114)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve coefficient_c_d from CCPP data structure')
            return
        end if
        if (kind(c_d).ne.ckind) then
            call ccpp_error('Kind mismatch for variable coefficient_c_d')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'coefficient_w_0', w_0, ierr=ierr, dims=cdims, kind=ckind, index=118)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve coefficient_w_0 from CCPP data structure')
            return
        end if
        if (kind(w_0).ne.ckind) then
            call ccpp_error('Kind mismatch for variable coefficient_w_0')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'coefficient_w_d', w_d, ierr=ierr, dims=cdims, kind=ckind, index=119)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve coefficient_w_d from CCPP data structure')
            return
        end if
        if (kind(w_d).ne.ckind) then
            call ccpp_error('Kind mismatch for variable coefficient_w_d')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'free_convection_layer_thickness', d_conv, ierr=ierr, dims=cdims, kind=ckind, index=344)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve free_convection_layer_thickness from CCPP data structure')
            return
        end if
        if (kind(d_conv).ne.ckind) then
            call ccpp_error('Kind mismatch for variable free_convection_layer_thickness')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'index_of_dtlm_start', ifd, ierr=ierr, dims=cdims, kind=ckind, index=396)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve index_of_dtlm_start from CCPP data structure')
            return
        end if
        if (kind(ifd).ne.ckind) then
            call ccpp_error('Kind mismatch for variable index_of_dtlm_start')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'sensible_heat_flux_due_to_rainfall', qrain, ierr=ierr, dims=cdims, kind=ckind, index=637)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve sensible_heat_flux_due_to_rainfall from CCPP data structure')
            return
        end if
        if (kind(qrain).ne.ckind) then
            call ccpp_error('Kind mismatch for variable sensible_heat_flux_due_to_rainfall')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_specific_humidity', qsurf, ierr=ierr, dims=cdims, kind=ckind, index=730)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_specific_humidity from CCPP data structure')
            return
        end if
        if (kind(qsurf).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_specific_humidity')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'upward_heat_flux_in_soil', gflux, ierr=ierr, dims=cdims, kind=ckind, index=812)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve upward_heat_flux_in_soil from CCPP data structure')
            return
        end if
        if (kind(gflux).ne.ckind) then
            call ccpp_error('Kind mismatch for variable upward_heat_flux_in_soil')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_drag_wind_speed_for_momentum_in_air', cmm, ierr=ierr, dims=cdims, kind=ckind, index=705)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_drag_wind_speed_for_momentum_in_air from CCPP data structure')
            return
        end if
        if (kind(cmm).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_drag_wind_speed_for_momentum_in_air')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_drag_mass_flux_for_heat_and_moisture_in_air', chh, ierr=ierr, dims=cdims, kind=ckind, index=704)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_drag_mass_flux_for_heat_and_moisture_in_air from CCPP data structure')
            return
        end if
        if (kind(chh).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_drag_mass_flux_for_heat_and_moisture_in_air')
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
        

        call ccpp_field_get(cdata, 'surface_upward_potential_latent_heat_flux', ep, ierr=ierr, dims=cdims, kind=ckind, index=732)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_upward_potential_latent_heat_flux from CCPP data structure')
            return
        end if
        if (kind(ep).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_upward_potential_latent_heat_flux')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call sfc_nst_run(im=im,ps=ps,u1=u1,v1=v1,t1=t1,q1=q1,tref=tref,cm=cm,ch=ch,prsl1=prsl1,prslki=prslki, &
                  islimsk=islimsk,xlon=xlon,sinlat=sinlat,stress=stress,sfcemis=sfcemis,dlwflx=dlwflx, &
                  sfcnsw=sfcnsw,rain=rain,timestep=timestep,kdt=kdt,solhr=solhr,xcosz=xcosz, &
                  ddvel=ddvel,flag_iter=flag_iter,flag_guess=flag_guess,nstf_name1=nstf_name1, &
                  nstf_name4=nstf_name4,nstf_name5=nstf_name5,lprnt=lprnt,ipr=ipr,tskin=tskin, &
                  tsurf=tsurf,xt=xt,xs=xs,xu=xu,xv=xv,xz=xz,zm=zm,xtts=xtts,xzts=xzts,dt_cool=dt_cool, &
                  z_c=z_c,c_0=c_0,c_d=c_d,w_0=w_0,w_d=w_d,d_conv=d_conv,ifd=ifd,qrain=qrain, &
                  qsurf=qsurf,gflux=gflux,cmm=cmm,chh=chh,evap=evap,hflx=hflx,ep=ep,errmsg=cdata%errmsg, &
                  errflg=cdata%errflg)
        ierr=cdata%errflg

    end function sfc_nst_run_cap
end module sfc_nst_cap

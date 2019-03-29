
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
!! @brief Auto-generated cap module for the zhaocarr_precpd scheme
!!
!
module zhaocarr_precpd_cap

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_f_pointer, c_ptr, c_int32_t
    use            :: ccpp_types,                                      &
                      only: ccpp_t, CCPP_GENERIC_KIND
    use            :: ccpp_fields,                                     &
                      only: ccpp_field_get
    use            :: ccpp_errors,                                     &
                      only: ccpp_error, ccpp_debug
    use            :: zhaocarr_precpd, &
                      only: zhaocarr_precpd_init,zhaocarr_precpd_finalize,zhaocarr_precpd_run
    ! Other modules required, e.g. type definitions
    use machine, only: kind_phys

    implicit none

    private
    public :: zhaocarr_precpd_init_cap,zhaocarr_precpd_finalize_cap,zhaocarr_precpd_run_cap

    contains


    function zhaocarr_precpd_init_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call zhaocarr_precpd_init()
        

    end function zhaocarr_precpd_init_cap

    function zhaocarr_precpd_finalize_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call zhaocarr_precpd_finalize()
        

    end function zhaocarr_precpd_finalize_cap

    function zhaocarr_precpd_run_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind
        integer, pointer :: im
        integer, pointer :: ix
        integer, pointer :: km
        real(kind_phys), pointer :: dt
        real(kind_phys), pointer :: del(:,:)
        real(kind_phys), pointer :: prsl(:,:)
        real(kind_phys), pointer :: q(:,:)
        real(kind_phys), pointer :: cwm(:,:)
        real(kind_phys), pointer :: t(:,:)
        real(kind_phys), pointer :: rn(:)
        real(kind_phys), pointer :: sr(:)
        real(kind_phys), pointer :: rainp(:,:)
        real(kind_phys), pointer :: u00k(:,:)
        real(kind_phys), pointer :: psautco(:)
        real(kind_phys), pointer :: prautco(:)
        real(kind_phys), pointer :: evpco
        real(kind_phys), pointer :: wminco(:)
        real(kind_phys), pointer :: wk1(:)
        logical, pointer :: lprnt
        integer, pointer :: jpr

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
        

        call ccpp_field_get(cdata, 'time_step_for_physics', dt, ierr=ierr, kind=ckind, index=793)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve time_step_for_physics from CCPP data structure')
            return
        end if
        if (kind(dt).ne.ckind) then
            call ccpp_error('Kind mismatch for variable time_step_for_physics')
            ierr = 1
            return
        end if
#endif
        

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
        

        call ccpp_field_get(cdata, 'water_vapor_specific_humidity_updated_by_physics', q, ierr=ierr, dims=cdims, kind=ckind, index=860)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve water_vapor_specific_humidity_updated_by_physics from CCPP data structure')
            return
        end if
        if (kind(q).ne.ckind) then
            call ccpp_error('Kind mismatch for variable water_vapor_specific_humidity_updated_by_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cloud_condensed_water_mixing_ratio_updated_by_physics', cwm, ierr=ierr, dims=cdims, kind=ckind, index=95)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cloud_condensed_water_mixing_ratio_updated_by_physics from CCPP data structure')
            return
        end if
        if (kind(cwm).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cloud_condensed_water_mixing_ratio_updated_by_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'air_temperature_updated_by_physics', t, ierr=ierr, dims=cdims, kind=ckind, index=59)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve air_temperature_updated_by_physics from CCPP data structure')
            return
        end if
        if (kind(t).ne.ckind) then
            call ccpp_error('Kind mismatch for variable air_temperature_updated_by_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'lwe_thickness_of_explicit_precipitation_amount', rn, ierr=ierr, dims=cdims, kind=ckind, index=480)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve lwe_thickness_of_explicit_precipitation_amount from CCPP data structure')
            return
        end if
        if (kind(rn).ne.ckind) then
            call ccpp_error('Kind mismatch for variable lwe_thickness_of_explicit_precipitation_amount')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'ratio_of_snowfall_to_rainfall', sr, ierr=ierr, dims=cdims, kind=ckind, index=624)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve ratio_of_snowfall_to_rainfall from CCPP data structure')
            return
        end if
        if (kind(sr).ne.ckind) then
            call ccpp_error('Kind mismatch for variable ratio_of_snowfall_to_rainfall')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'tendency_of_rain_water_mixing_ratio_due_to_microphysics', rainp, ierr=ierr, dims=cdims, kind=ckind, index=775)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve tendency_of_rain_water_mixing_ratio_due_to_microphysics from CCPP data structure')
            return
        end if
        if (kind(rainp).ne.ckind) then
            call ccpp_error('Kind mismatch for variable tendency_of_rain_water_mixing_ratio_due_to_microphysics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'critical_relative_humidity', u00k, ierr=ierr, dims=cdims, kind=ckind, index=141)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve critical_relative_humidity from CCPP data structure')
            return
        end if
        if (kind(u00k).ne.ckind) then
            call ccpp_error('Kind mismatch for variable critical_relative_humidity')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'coefficient_from_cloud_ice_to_snow', psautco, ierr=ierr, dims=cdims, kind=ckind, index=116)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve coefficient_from_cloud_ice_to_snow from CCPP data structure')
            return
        end if
        if (kind(psautco).ne.ckind) then
            call ccpp_error('Kind mismatch for variable coefficient_from_cloud_ice_to_snow')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'coefficient_from_cloud_water_to_rain', prautco, ierr=ierr, dims=cdims, kind=ckind, index=117)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve coefficient_from_cloud_water_to_rain from CCPP data structure')
            return
        end if
        if (kind(prautco).ne.ckind) then
            call ccpp_error('Kind mismatch for variable coefficient_from_cloud_water_to_rain')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'coefficient_for_evaporation_of_rainfall', evpco, ierr=ierr, kind=ckind, index=115)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve coefficient_for_evaporation_of_rainfall from CCPP data structure')
            return
        end if
        if (kind(evpco).ne.ckind) then
            call ccpp_error('Kind mismatch for variable coefficient_for_evaporation_of_rainfall')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'cloud_condensed_water_conversion_threshold', wminco, ierr=ierr, dims=cdims, kind=ckind, index=89)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cloud_condensed_water_conversion_threshold from CCPP data structure')
            return
        end if
        if (kind(wminco).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cloud_condensed_water_conversion_threshold')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'grid_size_related_coefficient_used_in_scale-sensitive_schemes', wk1, ierr=ierr, dims=cdims, kind=ckind, index=357)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve grid_size_related_coefficient_used_in_scale-sensitive_schemes from CCPP data structure')
            return
        end if
        if (kind(wk1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable grid_size_related_coefficient_used_in_scale-sensitive_schemes')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

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
        

        call ccpp_field_get(cdata, 'horizontal_index_of_printed_column', jpr, ierr=ierr, kind=ckind, index=365)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve horizontal_index_of_printed_column from CCPP data structure')
            return
        end if
        if (kind(jpr).ne.ckind) then
            call ccpp_error('Kind mismatch for variable horizontal_index_of_printed_column')
            ierr = 1
            return
        end if
#endif
        

        call zhaocarr_precpd_run(im=im,ix=ix,km=km,dt=dt,del=del,prsl=prsl,q=q,cwm=cwm,t=t,rn=rn,sr=sr,rainp=rainp, &
                  u00k=u00k,psautco=psautco,prautco=prautco,evpco=evpco,wminco=wminco,wk1=wk1, &
                  lprnt=lprnt,jpr=jpr,errmsg=cdata%errmsg,errflg=cdata%errflg)
        ierr=cdata%errflg

    end function zhaocarr_precpd_run_cap
end module zhaocarr_precpd_cap

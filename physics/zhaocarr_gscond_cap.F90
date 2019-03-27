
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
!! @brief Auto-generated cap module for the zhaocarr_gscond scheme
!!
!
module zhaocarr_gscond_cap

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_f_pointer, c_ptr, c_int32_t
    use            :: ccpp_types,                                      &
                      only: ccpp_t, CCPP_GENERIC_KIND
    use            :: ccpp_fields,                                     &
                      only: ccpp_field_get
    use            :: ccpp_errors,                                     &
                      only: ccpp_error, ccpp_debug
    use            :: zhaocarr_gscond, &
                      only: zhaocarr_gscond_finalize,zhaocarr_gscond_run,zhaocarr_gscond_init
    ! Other modules required, e.g. type definitions
    use machine, only: kind_phys

    implicit none

    private
    public :: zhaocarr_gscond_finalize_cap,zhaocarr_gscond_run_cap,zhaocarr_gscond_init_cap

    contains


    function zhaocarr_gscond_finalize_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call zhaocarr_gscond_finalize()
        

    end function zhaocarr_gscond_finalize_cap

    function zhaocarr_gscond_run_cap(ptr) bind(c) result(ierr)

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
        real(kind_phys), pointer :: dtf
        real(kind_phys), pointer :: prsl(:,:)
        real(kind_phys), pointer :: ps(:)
        real(kind_phys), pointer :: q(:,:)
        real(kind_phys), pointer :: clw1(:,:)
        real(kind_phys), pointer :: clw2(:,:)
        real(kind_phys), pointer :: cwm(:,:)
        real(kind_phys), pointer :: t(:,:)
        real(kind_phys), pointer :: tp(:,:)
        real(kind_phys), pointer :: qp(:,:)
        real(kind_phys), pointer :: psp(:)
        real(kind_phys), pointer :: tp1(:,:)
        real(kind_phys), pointer :: qp1(:,:)
        real(kind_phys), pointer :: psp1(:)
        real(kind_phys), pointer :: u(:,:)
        logical, pointer :: lprnt
        integer, pointer :: ipr

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
        

        call ccpp_field_get(cdata, 'ice_water_mixing_ratio_convective_transport_tracer', clw1, ierr=ierr, dims=cdims, kind=ckind, index=374)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve ice_water_mixing_ratio_convective_transport_tracer from CCPP data structure')
            return
        end if
        if (kind(clw1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable ice_water_mixing_ratio_convective_transport_tracer')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cloud_condensed_water_mixing_ratio_convective_transport_tracer', clw2, ierr=ierr, dims=cdims, kind=ckind, index=93)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cloud_condensed_water_mixing_ratio_convective_transport_tracer from CCPP data structure')
            return
        end if
        if (kind(clw2).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cloud_condensed_water_mixing_ratio_convective_transport_tracer')
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
        

        call ccpp_field_get(cdata, 'air_temperature_two_time_steps_back', tp, ierr=ierr, dims=cdims, kind=ckind, index=58)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve air_temperature_two_time_steps_back from CCPP data structure')
            return
        end if
        if (kind(tp).ne.ckind) then
            call ccpp_error('Kind mismatch for variable air_temperature_two_time_steps_back')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'water_vapor_specific_humidity_two_time_steps_back', qp, ierr=ierr, dims=cdims, kind=ckind, index=859)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve water_vapor_specific_humidity_two_time_steps_back from CCPP data structure')
            return
        end if
        if (kind(qp).ne.ckind) then
            call ccpp_error('Kind mismatch for variable water_vapor_specific_humidity_two_time_steps_back')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_air_pressure_two_time_steps_back', psp, ierr=ierr, dims=cdims, kind=ckind, index=680)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_air_pressure_two_time_steps_back from CCPP data structure')
            return
        end if
        if (kind(psp).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_air_pressure_two_time_steps_back')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'air_temperature_at_previous_time_step', tp1, ierr=ierr, dims=cdims, kind=ckind, index=56)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve air_temperature_at_previous_time_step from CCPP data structure')
            return
        end if
        if (kind(tp1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable air_temperature_at_previous_time_step')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'water_vapor_specific_humidity_at_previous_time_step', qp1, ierr=ierr, dims=cdims, kind=ckind, index=857)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve water_vapor_specific_humidity_at_previous_time_step from CCPP data structure')
            return
        end if
        if (kind(qp1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable water_vapor_specific_humidity_at_previous_time_step')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_air_pressure_at_previous_time_step', psp1, ierr=ierr, dims=cdims, kind=ckind, index=678)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_air_pressure_at_previous_time_step from CCPP data structure')
            return
        end if
        if (kind(psp1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_air_pressure_at_previous_time_step')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'critical_relative_humidity', u, ierr=ierr, dims=cdims, kind=ckind, index=141)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve critical_relative_humidity from CCPP data structure')
            return
        end if
        if (kind(u).ne.ckind) then
            call ccpp_error('Kind mismatch for variable critical_relative_humidity')
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
        

        call zhaocarr_gscond_run(im=im,ix=ix,km=km,dt=dt,dtf=dtf,prsl=prsl,ps=ps,q=q,clw1=clw1,clw2=clw2, &
                  cwm=cwm,t=t,tp=tp,qp=qp,psp=psp,tp1=tp1,qp1=qp1,psp1=psp1,u=u,lprnt=lprnt, &
                  ipr=ipr,errmsg=cdata%errmsg,errflg=cdata%errflg)
        ierr=cdata%errflg

    end function zhaocarr_gscond_run_cap

    function zhaocarr_gscond_init_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call zhaocarr_gscond_init()
        

    end function zhaocarr_gscond_init_cap
end module zhaocarr_gscond_cap

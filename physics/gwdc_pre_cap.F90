
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
!! @brief Auto-generated cap module for the gwdc_pre scheme
!!
!
module gwdc_pre_cap

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_f_pointer, c_ptr, c_int32_t
    use            :: ccpp_types,                                      &
                      only: ccpp_t, CCPP_GENERIC_KIND
    use            :: ccpp_fields,                                     &
                      only: ccpp_field_get
    use            :: ccpp_errors,                                     &
                      only: ccpp_error, ccpp_debug
    use            :: gwdc_pre, &
                      only: gwdc_pre_finalize,gwdc_pre_init,gwdc_pre_run
    ! Other modules required, e.g. type definitions
    use machine, only: kind_phys

    implicit none

    private
    public :: gwdc_pre_finalize_cap,gwdc_pre_init_cap,gwdc_pre_run_cap

    contains


    function gwdc_pre_finalize_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call gwdc_pre_finalize()
        

    end function gwdc_pre_finalize_cap

    function gwdc_pre_init_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call gwdc_pre_init()
        

    end function gwdc_pre_init_cap

    function gwdc_pre_run_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind
        integer, pointer :: im
        real(kind_phys), pointer :: cgwf(:)
        real(kind_phys), pointer :: dx(:)
        real(kind_phys), pointer :: work1(:)
        real(kind_phys), pointer :: work2(:)
        real(kind_phys), pointer :: dlength(:)
        real(kind_phys), pointer :: cldf(:)
        integer, pointer :: levs
        integer, pointer :: kbot(:)
        integer, pointer :: ktop(:)
        real(kind_phys), pointer :: dtp
        real(kind_phys), pointer :: gt0(:,:)
        real(kind_phys), pointer :: gt0_init(:,:)
        real(kind_phys), pointer :: del(:,:)
        real(kind_phys), pointer :: cumabs(:)

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
        

        call ccpp_field_get(cdata, 'multiplication_factors_for_convective_gravity_wave_drag', cgwf, ierr=ierr, dims=cdims, kind=ckind, index=561)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve multiplication_factors_for_convective_gravity_wave_drag from CCPP data structure')
            return
        end if
        if (kind(cgwf).ne.ckind) then
            call ccpp_error('Kind mismatch for variable multiplication_factors_for_convective_gravity_wave_drag')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

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
        

        call ccpp_field_get(cdata, 'characteristic_grid_length_scale', dlength, ierr=ierr, dims=cdims, kind=ckind, index=85)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve characteristic_grid_length_scale from CCPP data structure')
            return
        end if
        if (kind(dlength).ne.ckind) then
            call ccpp_error('Kind mismatch for variable characteristic_grid_length_scale')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cloud_area_fraction', cldf, ierr=ierr, dims=cdims, kind=ckind, index=86)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cloud_area_fraction from CCPP data structure')
            return
        end if
        if (kind(cldf).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cloud_area_fraction')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

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
        

        call ccpp_field_get(cdata, 'vertical_index_at_cloud_base', kbot, ierr=ierr, dims=cdims, kind=ckind, index=820)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve vertical_index_at_cloud_base from CCPP data structure')
            return
        end if
        if (kind(kbot).ne.ckind) then
            call ccpp_error('Kind mismatch for variable vertical_index_at_cloud_base')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'vertical_index_at_cloud_top', ktop, ierr=ierr, dims=cdims, kind=ckind, index=821)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve vertical_index_at_cloud_top from CCPP data structure')
            return
        end if
        if (kind(ktop).ne.ckind) then
            call ccpp_error('Kind mismatch for variable vertical_index_at_cloud_top')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'time_step_for_physics', dtp, ierr=ierr, kind=ckind, index=793)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve time_step_for_physics from CCPP data structure')
            return
        end if
        if (kind(dtp).ne.ckind) then
            call ccpp_error('Kind mismatch for variable time_step_for_physics')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'air_temperature_updated_by_physics', gt0, ierr=ierr, dims=cdims, kind=ckind, index=59)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve air_temperature_updated_by_physics from CCPP data structure')
            return
        end if
        if (kind(gt0).ne.ckind) then
            call ccpp_error('Kind mismatch for variable air_temperature_updated_by_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'air_temperature_save', gt0_init, ierr=ierr, dims=cdims, kind=ckind, index=57)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve air_temperature_save from CCPP data structure')
            return
        end if
        if (kind(gt0_init).ne.ckind) then
            call ccpp_error('Kind mismatch for variable air_temperature_save')
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
        

        call ccpp_field_get(cdata, 'maximum_column_heating_rate', cumabs, ierr=ierr, dims=cdims, kind=ckind, index=503)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve maximum_column_heating_rate from CCPP data structure')
            return
        end if
        if (kind(cumabs).ne.ckind) then
            call ccpp_error('Kind mismatch for variable maximum_column_heating_rate')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call gwdc_pre_run(im=im,cgwf=cgwf,dx=dx,work1=work1,work2=work2,dlength=dlength,cldf=cldf, &
                  levs=levs,kbot=kbot,ktop=ktop,dtp=dtp,gt0=gt0,gt0_init=gt0_init,del=del, &
                  cumabs=cumabs,errmsg=cdata%errmsg,errflg=cdata%errflg)
        ierr=cdata%errflg

    end function gwdc_pre_run_cap
end module gwdc_pre_cap


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
!! @brief Auto-generated cap module for the cs_conv_pre scheme
!!
!
module cs_conv_pre_cap

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_f_pointer, c_ptr, c_int32_t
    use            :: ccpp_types,                                      &
                      only: ccpp_t, CCPP_GENERIC_KIND
    use            :: ccpp_fields,                                     &
                      only: ccpp_field_get
    use            :: ccpp_errors,                                     &
                      only: ccpp_error, ccpp_debug
    use            :: cs_conv_pre, &
                      only: cs_conv_pre_run,cs_conv_pre_finalize,cs_conv_pre_init
    ! Other modules required, e.g. type definitions
    use machine, only: kind_phys

    implicit none

    private
    public :: cs_conv_pre_run_cap,cs_conv_pre_finalize_cap,cs_conv_pre_init_cap

    contains


    function cs_conv_pre_run_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind
        integer, pointer :: im
        integer, pointer :: levs
        integer, pointer :: ntrac
        integer, pointer :: ncld
        real(kind_phys), pointer :: q(:,:)
        real(kind_phys), pointer :: clw1(:,:)
        real(kind_phys), pointer :: clw2(:,:)
        real(kind_phys), pointer :: work1(:)
        real(kind_phys), pointer :: work2(:)
        real(kind_phys), pointer :: cs_parm1
        real(kind_phys), pointer :: cs_parm2
        real(kind_phys), pointer :: wcbmax(:)
        real(kind_phys), pointer :: fswtr(:)
        real(kind_phys), pointer :: fscav(:)
        real(kind_phys), pointer :: save_q1(:,:)
        real(kind_phys), pointer :: save_q2(:,:)
        real(kind_phys), pointer :: save_q3(:,:)

        ierr = 0

        call c_f_pointer(ptr, cdata)


        call ccpp_field_get(cdata, 'horizontal_dimension', im, ierr=ierr, kind=ckind, index=364)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve horizontal_dimension from CCPP data structure')
            return
        end if
        if (kind(im).ne.ckind) then
            call ccpp_error('Kind mismatch for variable horizontal_dimension')
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
        

        call ccpp_field_get(cdata, 'number_of_tracers', ntrac, ierr=ierr, kind=ckind, index=584)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve number_of_tracers from CCPP data structure')
            return
        end if
        if (kind(ntrac).ne.ckind) then
            call ccpp_error('Kind mismatch for variable number_of_tracers')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'number_of_hydrometeors', ncld, ierr=ierr, kind=ckind, index=578)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve number_of_hydrometeors from CCPP data structure')
            return
        end if
        if (kind(ncld).ne.ckind) then
            call ccpp_error('Kind mismatch for variable number_of_hydrometeors')
            ierr = 1
            return
        end if
#endif
        

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
        

        call ccpp_field_get(cdata, 'updraft_velocity_tunable_parameter_1_CS', cs_parm1, ierr=ierr, kind=ckind, index=809)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve updraft_velocity_tunable_parameter_1_CS from CCPP data structure')
            return
        end if
        if (kind(cs_parm1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable updraft_velocity_tunable_parameter_1_CS')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'updraft_velocity_tunable_parameter_2_CS', cs_parm2, ierr=ierr, kind=ckind, index=810)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve updraft_velocity_tunable_parameter_2_CS from CCPP data structure')
            return
        end if
        if (kind(cs_parm2).ne.ckind) then
            call ccpp_error('Kind mismatch for variable updraft_velocity_tunable_parameter_2_CS')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'maximum_updraft_velocity_at_cloud_base', wcbmax, ierr=ierr, dims=cdims, kind=ckind, index=509)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve maximum_updraft_velocity_at_cloud_base from CCPP data structure')
            return
        end if
        if (kind(wcbmax).ne.ckind) then
            call ccpp_error('Kind mismatch for variable maximum_updraft_velocity_at_cloud_base')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'fraction_of_cloud_top_water_scavenged', fswtr, ierr=ierr, dims=cdims, kind=ckind, index=340)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve fraction_of_cloud_top_water_scavenged from CCPP data structure')
            return
        end if
        if (kind(fswtr).ne.ckind) then
            call ccpp_error('Kind mismatch for variable fraction_of_cloud_top_water_scavenged')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'fraction_of_tracer_scavenged', fscav, ierr=ierr, dims=cdims, kind=ckind, index=343)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve fraction_of_tracer_scavenged from CCPP data structure')
            return
        end if
        if (kind(fscav).ne.ckind) then
            call ccpp_error('Kind mismatch for variable fraction_of_tracer_scavenged')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'water_vapor_specific_humidity_save', save_q1, ierr=ierr, dims=cdims, kind=ckind, index=858)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve water_vapor_specific_humidity_save from CCPP data structure')
            return
        end if
        if (kind(save_q1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable water_vapor_specific_humidity_save')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cloud_condensed_water_mixing_ratio_save', save_q2, ierr=ierr, dims=cdims, kind=ckind, index=94)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cloud_condensed_water_mixing_ratio_save from CCPP data structure')
            return
        end if
        if (kind(save_q2).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cloud_condensed_water_mixing_ratio_save')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'ice_water_mixing_ratio_save', save_q3, ierr=ierr, dims=cdims, kind=ckind, index=375)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve ice_water_mixing_ratio_save from CCPP data structure')
            return
        end if
        if (kind(save_q3).ne.ckind) then
            call ccpp_error('Kind mismatch for variable ice_water_mixing_ratio_save')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call cs_conv_pre_run(im=im,levs=levs,ntrac=ntrac,ncld=ncld,q=q,clw1=clw1,clw2=clw2,work1=work1, &
                  work2=work2,cs_parm1=cs_parm1,cs_parm2=cs_parm2,wcbmax=wcbmax,fswtr=fswtr, &
                  fscav=fscav,save_q1=save_q1,save_q2=save_q2,save_q3=save_q3,errmsg=cdata%errmsg, &
                  errflg=cdata%errflg)
        ierr=cdata%errflg

    end function cs_conv_pre_run_cap

    function cs_conv_pre_finalize_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call cs_conv_pre_finalize()
        

    end function cs_conv_pre_finalize_cap

    function cs_conv_pre_init_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call cs_conv_pre_init()
        

    end function cs_conv_pre_init_cap
end module cs_conv_pre_cap

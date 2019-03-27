
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
!! @brief Auto-generated cap module for the GFS_MP_generic_pre scheme
!!
!
module GFS_MP_generic_pre_cap

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_f_pointer, c_ptr, c_int32_t
    use            :: ccpp_types,                                      &
                      only: ccpp_t, CCPP_GENERIC_KIND
    use            :: ccpp_fields,                                     &
                      only: ccpp_field_get
    use            :: ccpp_errors,                                     &
                      only: ccpp_error, ccpp_debug
    use            :: GFS_MP_generic_pre, &
                      only: GFS_MP_generic_pre_init,GFS_MP_generic_pre_finalize,GFS_MP_generic_pre_run
    ! Other modules required, e.g. type definitions
    use machine, only: kind_phys

    implicit none

    private
    public :: GFS_MP_generic_pre_init_cap,GFS_MP_generic_pre_finalize_cap,GFS_MP_generic_pre_run_cap

    contains


    function GFS_MP_generic_pre_init_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call GFS_MP_generic_pre_init()
        

    end function GFS_MP_generic_pre_init_cap

    function GFS_MP_generic_pre_finalize_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call GFS_MP_generic_pre_finalize()
        

    end function GFS_MP_generic_pre_finalize_cap

    function GFS_MP_generic_pre_run_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind
        integer, pointer :: im
        integer, pointer :: levs
        logical, pointer :: ldiag3d
        logical, pointer :: do_aw
        integer, pointer :: ntcw
        integer, pointer :: nncl
        integer, pointer :: ntrac
        real(kind_phys), pointer :: gt0(:,:)
        real(kind_phys), pointer :: gq0(:,:,:)
        real(kind_phys), pointer :: save_t(:,:)
        real(kind_phys), pointer :: save_q(:,:,:)

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
        

        call ccpp_field_get(cdata, 'flag_for_Arakawa_Wu_adjustment', do_aw, ierr=ierr, kind=ckind, index=267)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_Arakawa_Wu_adjustment from CCPP data structure')
            return
        end if
        if (kind(do_aw).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_Arakawa_Wu_adjustment')
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
        

        call ccpp_field_get(cdata, 'number_of_tracers_for_cloud_condensate', nncl, ierr=ierr, kind=ckind, index=586)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve number_of_tracers_for_cloud_condensate from CCPP data structure')
            return
        end if
        if (kind(nncl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable number_of_tracers_for_cloud_condensate')
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
        

        call ccpp_field_get(cdata, 'tracer_concentration_updated_by_physics', gq0, ierr=ierr, dims=cdims, kind=ckind, index=803)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve tracer_concentration_updated_by_physics from CCPP data structure')
            return
        end if
        if (kind(gq0).ne.ckind) then
            call ccpp_error('Kind mismatch for variable tracer_concentration_updated_by_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'air_temperature_save', save_t, ierr=ierr, dims=cdims, kind=ckind, index=57)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve air_temperature_save from CCPP data structure')
            return
        end if
        if (kind(save_t).ne.ckind) then
            call ccpp_error('Kind mismatch for variable air_temperature_save')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'tracer_concentration_save', save_q, ierr=ierr, dims=cdims, kind=ckind, index=802)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve tracer_concentration_save from CCPP data structure')
            return
        end if
        if (kind(save_q).ne.ckind) then
            call ccpp_error('Kind mismatch for variable tracer_concentration_save')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call GFS_MP_generic_pre_run(im=im,levs=levs,ldiag3d=ldiag3d,do_aw=do_aw,ntcw=ntcw,nncl=nncl,ntrac=ntrac, &
                  gt0=gt0,gq0=gq0,save_t=save_t,save_q=save_q,errmsg=cdata%errmsg,errflg=cdata%errflg)
        ierr=cdata%errflg

    end function GFS_MP_generic_pre_run_cap
end module GFS_MP_generic_pre_cap

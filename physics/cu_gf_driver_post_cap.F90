
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
!! @brief Auto-generated cap module for the cu_gf_driver_post scheme
!!
!
module cu_gf_driver_post_cap

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_f_pointer, c_ptr, c_int32_t
    use            :: ccpp_types,                                      &
                      only: ccpp_t, CCPP_GENERIC_KIND
    use            :: ccpp_fields,                                     &
                      only: ccpp_field_get
    use            :: ccpp_errors,                                     &
                      only: ccpp_error, ccpp_debug
    use            :: cu_gf_driver_post, &
                      only: cu_gf_driver_post_init,cu_gf_driver_post_finalize,cu_gf_driver_post_run
    ! Other modules required, e.g. type definitions
    use machine, only: kind_phys

    implicit none

    private
    public :: cu_gf_driver_post_init_cap,cu_gf_driver_post_finalize_cap,cu_gf_driver_post_run_cap

    contains


    function cu_gf_driver_post_init_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call cu_gf_driver_post_init()
        

    end function cu_gf_driver_post_init_cap

    function cu_gf_driver_post_finalize_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call cu_gf_driver_post_finalize()
        

    end function cu_gf_driver_post_finalize_cap

    function cu_gf_driver_post_run_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind
        integer, pointer :: im
        real(kind_phys), pointer :: t(:,:)
        real(kind_phys), pointer :: q(:,:)
        real(kind_phys), pointer :: prevst(:,:)
        real(kind_phys), pointer :: prevsq(:,:)
        integer, pointer :: cactiv(:)
        real(kind_phys), pointer :: conv_act(:)

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
        

        call ccpp_field_get(cdata, 'temperature_from_previous_timestep', prevst, ierr=ierr, dims=cdims, kind=ckind, index=752)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve temperature_from_previous_timestep from CCPP data structure')
            return
        end if
        if (kind(prevst).ne.ckind) then
            call ccpp_error('Kind mismatch for variable temperature_from_previous_timestep')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'moisture_from_previous_timestep', prevsq, ierr=ierr, dims=cdims, kind=ckind, index=553)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve moisture_from_previous_timestep from CCPP data structure')
            return
        end if
        if (kind(prevsq).ne.ckind) then
            call ccpp_error('Kind mismatch for variable moisture_from_previous_timestep')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'conv_activity_counter', cactiv, ierr=ierr, dims=cdims, kind=ckind, index=122)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve conv_activity_counter from CCPP data structure')
            return
        end if
        if (kind(cactiv).ne.ckind) then
            call ccpp_error('Kind mismatch for variable conv_activity_counter')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'gf_memory_counter', conv_act, ierr=ierr, dims=cdims, kind=ckind, index=351)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve gf_memory_counter from CCPP data structure')
            return
        end if
        if (kind(conv_act).ne.ckind) then
            call ccpp_error('Kind mismatch for variable gf_memory_counter')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call cu_gf_driver_post_run(im=im,t=t,q=q,prevst=prevst,prevsq=prevsq,cactiv=cactiv,conv_act=conv_act, &
                  errmsg=cdata%errmsg,errflg=cdata%errflg)
        ierr=cdata%errflg

    end function cu_gf_driver_post_run_cap
end module cu_gf_driver_post_cap


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
!! @brief Auto-generated cap module for the GFS_surface_loop_control_part2 scheme
!!
!
module GFS_surface_loop_control_part2_cap

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_f_pointer, c_ptr, c_int32_t
    use            :: ccpp_types,                                      &
                      only: ccpp_t, CCPP_GENERIC_KIND
    use            :: ccpp_fields,                                     &
                      only: ccpp_field_get
    use            :: ccpp_errors,                                     &
                      only: ccpp_error, ccpp_debug
    use            :: GFS_surface_loop_control_part2, &
                      only: GFS_surface_loop_control_part2_finalize,GFS_surface_loop_control_part2_init,GFS_surface_loop_control_part2_run
    ! Other modules required, e.g. type definitions
    use machine, only: kind_phys

    implicit none

    private
    public :: GFS_surface_loop_control_part2_finalize_cap,GFS_surface_loop_control_part2_init_cap,GFS_surface_loop_control_part2_run_cap

    contains


    function GFS_surface_loop_control_part2_finalize_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call GFS_surface_loop_control_part2_finalize()
        

    end function GFS_surface_loop_control_part2_finalize_cap

    function GFS_surface_loop_control_part2_init_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call GFS_surface_loop_control_part2_init()
        

    end function GFS_surface_loop_control_part2_init_cap

    function GFS_surface_loop_control_part2_run_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind
        integer, pointer :: im
        real(kind_phys), pointer :: wind(:)
        logical, pointer :: flag_guess(:)
        logical, pointer :: flag_iter(:)
        integer, pointer :: islmsk(:)
        integer, pointer :: nstf_name1

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
        

        call GFS_surface_loop_control_part2_run(im=im,iter=cdata%loop_cnt,wind=wind,flag_guess=flag_guess,flag_iter=flag_iter, &
                  islmsk=islmsk,nstf_name1=nstf_name1,errmsg=cdata%errmsg,errflg=cdata%errflg)
        ierr=cdata%errflg

    end function GFS_surface_loop_control_part2_run_cap
end module GFS_surface_loop_control_part2_cap

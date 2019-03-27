
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
!! @brief Auto-generated cap module for the cs_conv_post scheme
!!
!
module cs_conv_post_cap

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_f_pointer, c_ptr, c_int32_t
    use            :: ccpp_types,                                      &
                      only: ccpp_t, CCPP_GENERIC_KIND
    use            :: ccpp_fields,                                     &
                      only: ccpp_field_get
    use            :: ccpp_errors,                                     &
                      only: ccpp_error, ccpp_debug
    use            :: cs_conv_post, &
                      only: cs_conv_post_init,cs_conv_post_finalize,cs_conv_post_run
    ! Other modules required, e.g. type definitions
    use machine, only: kind_phys

    implicit none

    private
    public :: cs_conv_post_init_cap,cs_conv_post_finalize_cap,cs_conv_post_run_cap

    contains


    function cs_conv_post_init_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call cs_conv_post_init()
        

    end function cs_conv_post_init_cap

    function cs_conv_post_finalize_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call cs_conv_post_finalize()
        

    end function cs_conv_post_finalize_cap

    function cs_conv_post_run_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind
        integer, pointer :: im
        integer, pointer :: kmax
        logical, pointer :: do_aw
        real(kind_phys), pointer :: sigmatot(:,:)
        real(kind_phys), pointer :: sigmafrac(:,:)

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
        

        call ccpp_field_get(cdata, 'vertical_dimension', kmax, ierr=ierr, kind=ckind, index=817)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve vertical_dimension from CCPP data structure')
            return
        end if
        if (kind(kmax).ne.ckind) then
            call ccpp_error('Kind mismatch for variable vertical_dimension')
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
        

        call ccpp_field_get(cdata, 'convective_updraft_area_fraction_at_model_interfaces', sigmatot, ierr=ierr, dims=cdims, kind=ckind, index=132)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve convective_updraft_area_fraction_at_model_interfaces from CCPP data structure')
            return
        end if
        if (kind(sigmatot).ne.ckind) then
            call ccpp_error('Kind mismatch for variable convective_updraft_area_fraction_at_model_interfaces')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'convective_updraft_area_fraction', sigmafrac, ierr=ierr, dims=cdims, kind=ckind, index=131)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve convective_updraft_area_fraction from CCPP data structure')
            return
        end if
        if (kind(sigmafrac).ne.ckind) then
            call ccpp_error('Kind mismatch for variable convective_updraft_area_fraction')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call cs_conv_post_run(im=im,kmax=kmax,do_aw=do_aw,sigmatot=sigmatot,sigmafrac=sigmafrac,errmsg=cdata%errmsg, &
                  errflg=cdata%errflg)
        ierr=cdata%errflg

    end function cs_conv_post_run_cap
end module cs_conv_post_cap

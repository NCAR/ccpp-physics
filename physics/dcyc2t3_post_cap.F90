
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
!! @brief Auto-generated cap module for the dcyc2t3_post scheme
!!
!
module dcyc2t3_post_cap

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_f_pointer, c_ptr, c_int32_t
    use            :: ccpp_types,                                      &
                      only: ccpp_t, CCPP_GENERIC_KIND
    use            :: ccpp_fields,                                     &
                      only: ccpp_field_get
    use            :: ccpp_errors,                                     &
                      only: ccpp_error, ccpp_debug
    use            :: dcyc2t3_post, &
                      only: dcyc2t3_post_init,dcyc2t3_post_finalize,dcyc2t3_post_run
    ! Other modules required, e.g. type definitions
    use machine, only: kind_phys

    implicit none

    private
    public :: dcyc2t3_post_init_cap,dcyc2t3_post_finalize_cap,dcyc2t3_post_run_cap

    contains


    function dcyc2t3_post_init_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call dcyc2t3_post_init()
        

    end function dcyc2t3_post_init_cap

    function dcyc2t3_post_finalize_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call dcyc2t3_post_finalize()
        

    end function dcyc2t3_post_finalize_cap

    function dcyc2t3_post_run_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind
        integer, pointer :: im
        real(kind_phys), pointer :: adjsfcdsw(:)
        real(kind_phys), pointer :: adjsfcnsw(:)
        real(kind_phys), pointer :: adjsfcusw(:)

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
        

        call ccpp_field_get(cdata, 'surface_downwelling_shortwave_flux', adjsfcdsw, ierr=ierr, dims=cdims, kind=ckind, index=700)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_downwelling_shortwave_flux from CCPP data structure')
            return
        end if
        if (kind(adjsfcdsw).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_downwelling_shortwave_flux')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_net_downwelling_shortwave_flux', adjsfcnsw, ierr=ierr, dims=cdims, kind=ckind, index=716)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_net_downwelling_shortwave_flux from CCPP data structure')
            return
        end if
        if (kind(adjsfcnsw).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_net_downwelling_shortwave_flux')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_upwelling_shortwave_flux', adjsfcusw, ierr=ierr, dims=cdims, kind=ckind, index=743)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_upwelling_shortwave_flux from CCPP data structure')
            return
        end if
        if (kind(adjsfcusw).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_upwelling_shortwave_flux')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call dcyc2t3_post_run(im=im,adjsfcdsw=adjsfcdsw,adjsfcnsw=adjsfcnsw,adjsfcusw=adjsfcusw,errmsg=cdata%errmsg, &
                  errflg=cdata%errflg)
        ierr=cdata%errflg

    end function dcyc2t3_post_run_cap
end module dcyc2t3_post_cap

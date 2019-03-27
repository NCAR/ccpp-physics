
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
!! @brief Auto-generated cap module for the GFS_suite_interstitial_phys_reset scheme
!!
!
module GFS_suite_interstitial_phys_reset_cap

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_f_pointer, c_ptr, c_int32_t
    use            :: ccpp_types,                                      &
                      only: ccpp_t, CCPP_GENERIC_KIND
    use            :: ccpp_fields,                                     &
                      only: ccpp_field_get
    use            :: ccpp_errors,                                     &
                      only: ccpp_error, ccpp_debug
    use            :: GFS_suite_interstitial_phys_reset, &
                      only: GFS_suite_interstitial_phys_reset_finalize,GFS_suite_interstitial_phys_reset_run,GFS_suite_interstitial_phys_reset_init
    ! Other modules required, e.g. type definitions
    use GFS_typedefs, only: GFS_interstitial_type
    use GFS_typedefs, only: GFS_control_type

    implicit none

    private
    public :: GFS_suite_interstitial_phys_reset_finalize_cap,GFS_suite_interstitial_phys_reset_run_cap,GFS_suite_interstitial_phys_reset_init_cap

    contains


    function GFS_suite_interstitial_phys_reset_finalize_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call GFS_suite_interstitial_phys_reset_finalize()
        

    end function GFS_suite_interstitial_phys_reset_finalize_cap

    function GFS_suite_interstitial_phys_reset_run_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind
        type(GFS_interstitial_type), pointer     :: Interstitial
        type(GFS_control_type), pointer     :: Model

        ierr = 0

        call c_f_pointer(ptr, cdata)


        call ccpp_field_get(cdata, 'GFS_interstitial_type_instance', cptr, ierr=ierr, kind=ckind, index=7)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve GFS_interstitial_type_instance from CCPP data structure')
            return
        end if
        if (ckind.ne.CCPP_GENERIC_KIND) then
            call ccpp_error('Kind mismatch for variable GFS_interstitial_type_instance')
            ierr = 1
            return
        end if
#endif
        call c_f_pointer(cptr, Interstitial)

        call ccpp_field_get(cdata, 'GFS_control_type_instance', cptr, ierr=ierr, kind=ckind, index=2)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve GFS_control_type_instance from CCPP data structure')
            return
        end if
        if (ckind.ne.CCPP_GENERIC_KIND) then
            call ccpp_error('Kind mismatch for variable GFS_control_type_instance')
            ierr = 1
            return
        end if
#endif
        call c_f_pointer(cptr, Model)

        call GFS_suite_interstitial_phys_reset_run(Interstitial=Interstitial,Model=Model,errmsg=cdata%errmsg,errflg=cdata%errflg)
        ierr=cdata%errflg

    end function GFS_suite_interstitial_phys_reset_run_cap

    function GFS_suite_interstitial_phys_reset_init_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call GFS_suite_interstitial_phys_reset_init()
        

    end function GFS_suite_interstitial_phys_reset_init_cap
end module GFS_suite_interstitial_phys_reset_cap

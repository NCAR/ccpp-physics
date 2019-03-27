
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
!! @brief Auto-generated cap module for the stochastic_physics_sfc scheme
!!
!
module stochastic_physics_sfc_cap

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_f_pointer, c_ptr, c_int32_t
    use            :: ccpp_types,                                      &
                      only: ccpp_t, CCPP_GENERIC_KIND
    use            :: ccpp_fields,                                     &
                      only: ccpp_field_get
    use            :: ccpp_errors,                                     &
                      only: ccpp_error, ccpp_debug
    use            :: stochastic_physics_sfc, &
                      only: stochastic_physics_sfc_init,stochastic_physics_sfc_finalize,stochastic_physics_sfc_run
    ! Other modules required, e.g. type definitions
    use GFS_typedefs, only: GFS_control_type
    use GFS_typedefs, only: GFS_data_type

    implicit none

    private
    public :: stochastic_physics_sfc_init_cap,stochastic_physics_sfc_finalize_cap,stochastic_physics_sfc_run_cap

    contains


    function stochastic_physics_sfc_init_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind
        type(GFS_control_type), pointer     :: Model
        type(GFS_data_type), pointer     :: Data(:)

        ierr = 0

        call c_f_pointer(ptr, cdata)


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

        call ccpp_field_get(cdata, 'GFS_data_type_instance_all_blocks', cptr, ierr=ierr, dims=cdims, kind=ckind, index=4)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve GFS_data_type_instance_all_blocks from CCPP data structure')
            return
        end if
        if (ckind.ne.CCPP_GENERIC_KIND) then
            call ccpp_error('Kind mismatch for variable GFS_data_type_instance_all_blocks')
            ierr = 1
            return
        end if
#endif
        call c_f_pointer(cptr, Data, cdims)
        deallocate(cdims)
        

        call stochastic_physics_sfc_init(Model=Model,Data=Data,errmsg=cdata%errmsg,errflg=cdata%errflg)
        ierr=cdata%errflg

    end function stochastic_physics_sfc_init_cap

    function stochastic_physics_sfc_finalize_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call stochastic_physics_sfc_finalize()
        

    end function stochastic_physics_sfc_finalize_cap

    function stochastic_physics_sfc_run_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call stochastic_physics_sfc_run()
        

    end function stochastic_physics_sfc_run_cap
end module stochastic_physics_sfc_cap

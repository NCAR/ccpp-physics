
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
!! @brief Auto-generated cap module for the GFS_phys_time_vary scheme
!!
!
module GFS_phys_time_vary_cap

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_f_pointer, c_ptr, c_int32_t
    use            :: ccpp_types,                                      &
                      only: ccpp_t, CCPP_GENERIC_KIND
    use            :: ccpp_fields,                                     &
                      only: ccpp_field_get
    use            :: ccpp_errors,                                     &
                      only: ccpp_error, ccpp_debug
    use            :: GFS_phys_time_vary, &
                      only: GFS_phys_time_vary_init,GFS_phys_time_vary_finalize,GFS_phys_time_vary_run
    ! Other modules required, e.g. type definitions
    use GFS_typedefs, only: GFS_data_type
    use GFS_typedefs, only: GFS_control_type
    use GFS_typedefs, only: GFS_interstitial_type

    implicit none

    private
    public :: GFS_phys_time_vary_init_cap,GFS_phys_time_vary_finalize_cap,GFS_phys_time_vary_run_cap

    contains


    function GFS_phys_time_vary_init_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind
        type(GFS_data_type), pointer     :: Data(:)
        type(GFS_control_type), pointer     :: Model
        type(GFS_interstitial_type), pointer     :: Interstitial(:)
        integer, pointer :: nthrds

        ierr = 0

        call c_f_pointer(ptr, cdata)


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

        call ccpp_field_get(cdata, 'GFS_interstitial_type_instance_all_threads', cptr, ierr=ierr, dims=cdims, kind=ckind, index=8)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve GFS_interstitial_type_instance_all_threads from CCPP data structure')
            return
        end if
        if (ckind.ne.CCPP_GENERIC_KIND) then
            call ccpp_error('Kind mismatch for variable GFS_interstitial_type_instance_all_threads')
            ierr = 1
            return
        end if
#endif
        call c_f_pointer(cptr, Interstitial, cdims)
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'omp_threads', nthrds, ierr=ierr, kind=ckind, index=594)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve omp_threads from CCPP data structure')
            return
        end if
        if (kind(nthrds).ne.ckind) then
            call ccpp_error('Kind mismatch for variable omp_threads')
            ierr = 1
            return
        end if
#endif
        

        call GFS_phys_time_vary_init(Data=Data,Model=Model,Interstitial=Interstitial,nthrds=nthrds,errmsg=cdata%errmsg, &
                  errflg=cdata%errflg)
        ierr=cdata%errflg

    end function GFS_phys_time_vary_init_cap

    function GFS_phys_time_vary_finalize_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call GFS_phys_time_vary_finalize(errmsg=cdata%errmsg,errflg=cdata%errflg)
        ierr=cdata%errflg

    end function GFS_phys_time_vary_finalize_cap

    function GFS_phys_time_vary_run_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind
        type(GFS_data_type), pointer     :: Data(:)
        type(GFS_control_type), pointer     :: Model
        integer, pointer :: nthrds

        ierr = 0

        call c_f_pointer(ptr, cdata)


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

        call ccpp_field_get(cdata, 'omp_threads', nthrds, ierr=ierr, kind=ckind, index=594)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve omp_threads from CCPP data structure')
            return
        end if
        if (kind(nthrds).ne.ckind) then
            call ccpp_error('Kind mismatch for variable omp_threads')
            ierr = 1
            return
        end if
#endif
        

        call GFS_phys_time_vary_run(Data=Data,Model=Model,nthrds=nthrds,errmsg=cdata%errmsg,errflg=cdata%errflg)
        ierr=cdata%errflg

    end function GFS_phys_time_vary_run_cap
end module GFS_phys_time_vary_cap

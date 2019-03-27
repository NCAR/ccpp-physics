
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
!! @brief Auto-generated cap module for the memcheck scheme
!!
!
module memcheck_cap

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_f_pointer, c_ptr, c_int32_t
    use            :: ccpp_types,                                      &
                      only: ccpp_t, CCPP_GENERIC_KIND
    use            :: ccpp_fields,                                     &
                      only: ccpp_field_get
    use            :: ccpp_errors,                                     &
                      only: ccpp_error, ccpp_debug
    use            :: memcheck, &
                      only: memcheck_init,memcheck_finalize,memcheck_run
    ! Other modules required, e.g. type definitions
    

    implicit none

    private
    public :: memcheck_init_cap,memcheck_finalize_cap,memcheck_run_cap

    contains


    function memcheck_init_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call memcheck_init()
        

    end function memcheck_init_cap

    function memcheck_finalize_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call memcheck_finalize()
        

    end function memcheck_finalize_cap

    function memcheck_run_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind
        integer, pointer :: mpicomm
        integer, pointer :: mpirank
        integer, pointer :: mpisize
        integer, pointer :: mpiroot
        integer, pointer :: ompthreads

        ierr = 0

        call c_f_pointer(ptr, cdata)


        call ccpp_field_get(cdata, 'mpi_comm', mpicomm, ierr=ierr, kind=ckind, index=557)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve mpi_comm from CCPP data structure')
            return
        end if
        if (kind(mpicomm).ne.ckind) then
            call ccpp_error('Kind mismatch for variable mpi_comm')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'mpi_rank', mpirank, ierr=ierr, kind=ckind, index=558)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve mpi_rank from CCPP data structure')
            return
        end if
        if (kind(mpirank).ne.ckind) then
            call ccpp_error('Kind mismatch for variable mpi_rank')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'mpi_size', mpisize, ierr=ierr, kind=ckind, index=560)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve mpi_size from CCPP data structure')
            return
        end if
        if (kind(mpisize).ne.ckind) then
            call ccpp_error('Kind mismatch for variable mpi_size')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'mpi_root', mpiroot, ierr=ierr, kind=ckind, index=559)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve mpi_root from CCPP data structure')
            return
        end if
        if (kind(mpiroot).ne.ckind) then
            call ccpp_error('Kind mismatch for variable mpi_root')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'omp_threads', ompthreads, ierr=ierr, kind=ckind, index=594)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve omp_threads from CCPP data structure')
            return
        end if
        if (kind(ompthreads).ne.ckind) then
            call ccpp_error('Kind mismatch for variable omp_threads')
            ierr = 1
            return
        end if
#endif
        

        call memcheck_run(mpicomm=mpicomm,mpirank=mpirank,mpisize=mpisize,mpiroot=mpiroot,ompthreads=ompthreads, &
                  errmsg=cdata%errmsg,errflg=cdata%errflg)
        ierr=cdata%errflg

    end function memcheck_run_cap
end module memcheck_cap

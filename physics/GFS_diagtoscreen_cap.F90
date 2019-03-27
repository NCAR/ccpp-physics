
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
!! @brief Auto-generated cap module for the GFS_diagtoscreen scheme
!!
!
module GFS_diagtoscreen_cap

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_f_pointer, c_ptr, c_int32_t
    use            :: ccpp_types,                                      &
                      only: ccpp_t, CCPP_GENERIC_KIND
    use            :: ccpp_fields,                                     &
                      only: ccpp_field_get
    use            :: ccpp_errors,                                     &
                      only: ccpp_error, ccpp_debug
    use            :: GFS_diagtoscreen, &
                      only: GFS_diagtoscreen_init,GFS_diagtoscreen_finalize,GFS_diagtoscreen_run
    ! Other modules required, e.g. type definitions
    use GFS_typedefs, only: GFS_control_type
    use GFS_typedefs, only: GFS_statein_type
    use GFS_typedefs, only: GFS_stateout_type
    use GFS_typedefs, only: GFS_sfcprop_type
    use GFS_typedefs, only: GFS_coupling_type
    use GFS_typedefs, only: GFS_grid_type
    use GFS_typedefs, only: GFS_tbd_type
    use GFS_typedefs, only: GFS_cldprop_type
    use GFS_typedefs, only: GFS_radtend_type
    use GFS_typedefs, only: GFS_diag_type
    use GFS_typedefs, only: GFS_interstitial_type

    implicit none

    private
    public :: GFS_diagtoscreen_init_cap,GFS_diagtoscreen_finalize_cap,GFS_diagtoscreen_run_cap

    contains


    function GFS_diagtoscreen_init_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call GFS_diagtoscreen_init()
        

    end function GFS_diagtoscreen_init_cap

    function GFS_diagtoscreen_finalize_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call GFS_diagtoscreen_finalize()
        

    end function GFS_diagtoscreen_finalize_cap

    function GFS_diagtoscreen_run_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind
        type(GFS_control_type), pointer     :: Model
        type(GFS_statein_type), pointer     :: Statein
        type(GFS_stateout_type), pointer     :: Stateout
        type(GFS_sfcprop_type), pointer     :: Sfcprop
        type(GFS_coupling_type), pointer     :: Coupling
        type(GFS_grid_type), pointer     :: Grid
        type(GFS_tbd_type), pointer     :: Tbd
        type(GFS_cldprop_type), pointer     :: Cldprop
        type(GFS_radtend_type), pointer     :: Radtend
        type(GFS_diag_type), pointer     :: Diag
        type(GFS_interstitial_type), pointer     :: Interstitial
        integer, pointer :: nthreads
        integer, pointer :: blkno

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

        call ccpp_field_get(cdata, 'GFS_statein_type_instance', cptr, ierr=ierr, kind=ckind, index=13)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve GFS_statein_type_instance from CCPP data structure')
            return
        end if
        if (ckind.ne.CCPP_GENERIC_KIND) then
            call ccpp_error('Kind mismatch for variable GFS_statein_type_instance')
            ierr = 1
            return
        end if
#endif
        call c_f_pointer(cptr, Statein)

        call ccpp_field_get(cdata, 'GFS_stateout_type_instance', cptr, ierr=ierr, kind=ckind, index=14)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve GFS_stateout_type_instance from CCPP data structure')
            return
        end if
        if (ckind.ne.CCPP_GENERIC_KIND) then
            call ccpp_error('Kind mismatch for variable GFS_stateout_type_instance')
            ierr = 1
            return
        end if
#endif
        call c_f_pointer(cptr, Stateout)

        call ccpp_field_get(cdata, 'GFS_sfcprop_type_instance', cptr, ierr=ierr, kind=ckind, index=11)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve GFS_sfcprop_type_instance from CCPP data structure')
            return
        end if
        if (ckind.ne.CCPP_GENERIC_KIND) then
            call ccpp_error('Kind mismatch for variable GFS_sfcprop_type_instance')
            ierr = 1
            return
        end if
#endif
        call c_f_pointer(cptr, Sfcprop)

        call ccpp_field_get(cdata, 'GFS_coupling_type_instance', cptr, ierr=ierr, kind=ckind, index=3)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve GFS_coupling_type_instance from CCPP data structure')
            return
        end if
        if (ckind.ne.CCPP_GENERIC_KIND) then
            call ccpp_error('Kind mismatch for variable GFS_coupling_type_instance')
            ierr = 1
            return
        end if
#endif
        call c_f_pointer(cptr, Coupling)

        call ccpp_field_get(cdata, 'GFS_grid_type_instance', cptr, ierr=ierr, kind=ckind, index=6)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve GFS_grid_type_instance from CCPP data structure')
            return
        end if
        if (ckind.ne.CCPP_GENERIC_KIND) then
            call ccpp_error('Kind mismatch for variable GFS_grid_type_instance')
            ierr = 1
            return
        end if
#endif
        call c_f_pointer(cptr, Grid)

        call ccpp_field_get(cdata, 'GFS_tbd_type_instance', cptr, ierr=ierr, kind=ckind, index=15)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve GFS_tbd_type_instance from CCPP data structure')
            return
        end if
        if (ckind.ne.CCPP_GENERIC_KIND) then
            call ccpp_error('Kind mismatch for variable GFS_tbd_type_instance')
            ierr = 1
            return
        end if
#endif
        call c_f_pointer(cptr, Tbd)

        call ccpp_field_get(cdata, 'GFS_cldprop_type_instance', cptr, ierr=ierr, kind=ckind, index=1)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve GFS_cldprop_type_instance from CCPP data structure')
            return
        end if
        if (ckind.ne.CCPP_GENERIC_KIND) then
            call ccpp_error('Kind mismatch for variable GFS_cldprop_type_instance')
            ierr = 1
            return
        end if
#endif
        call c_f_pointer(cptr, Cldprop)

        call ccpp_field_get(cdata, 'GFS_radtend_type_instance', cptr, ierr=ierr, kind=ckind, index=9)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve GFS_radtend_type_instance from CCPP data structure')
            return
        end if
        if (ckind.ne.CCPP_GENERIC_KIND) then
            call ccpp_error('Kind mismatch for variable GFS_radtend_type_instance')
            ierr = 1
            return
        end if
#endif
        call c_f_pointer(cptr, Radtend)

        call ccpp_field_get(cdata, 'GFS_diag_type_instance', cptr, ierr=ierr, kind=ckind, index=5)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve GFS_diag_type_instance from CCPP data structure')
            return
        end if
        if (ckind.ne.CCPP_GENERIC_KIND) then
            call ccpp_error('Kind mismatch for variable GFS_diag_type_instance')
            ierr = 1
            return
        end if
#endif
        call c_f_pointer(cptr, Diag)

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

        call ccpp_field_get(cdata, 'omp_threads', nthreads, ierr=ierr, kind=ckind, index=594)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve omp_threads from CCPP data structure')
            return
        end if
        if (kind(nthreads).ne.ckind) then
            call ccpp_error('Kind mismatch for variable omp_threads')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'ccpp_block_number', blkno, ierr=ierr, kind=ckind, index=82)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve ccpp_block_number from CCPP data structure')
            return
        end if
        if (kind(blkno).ne.ckind) then
            call ccpp_error('Kind mismatch for variable ccpp_block_number')
            ierr = 1
            return
        end if
#endif
        

        call GFS_diagtoscreen_run(Model=Model,Statein=Statein,Stateout=Stateout,Sfcprop=Sfcprop,Coupling=Coupling, &
                  Grid=Grid,Tbd=Tbd,Cldprop=Cldprop,Radtend=Radtend,Diag=Diag,Interstitial=Interstitial, &
                  nthreads=nthreads,blkno=blkno,errmsg=cdata%errmsg,errflg=cdata%errflg)
        ierr=cdata%errflg

    end function GFS_diagtoscreen_run_cap
end module GFS_diagtoscreen_cap


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
!! @brief Auto-generated cap module for the sfc_nst_pre scheme
!!
!
module sfc_nst_pre_cap

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_f_pointer, c_ptr, c_int32_t
    use            :: ccpp_types,                                      &
                      only: ccpp_t, CCPP_GENERIC_KIND
    use            :: ccpp_fields,                                     &
                      only: ccpp_field_get
    use            :: ccpp_errors,                                     &
                      only: ccpp_error, ccpp_debug
    use            :: sfc_nst_pre, &
                      only: sfc_nst_pre_init,sfc_nst_pre_finalize,sfc_nst_pre_run
    ! Other modules required, e.g. type definitions
    use machine, only: kind_phys

    implicit none

    private
    public :: sfc_nst_pre_init_cap,sfc_nst_pre_finalize_cap,sfc_nst_pre_run_cap

    contains


    function sfc_nst_pre_init_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call sfc_nst_pre_init()
        

    end function sfc_nst_pre_init_cap

    function sfc_nst_pre_finalize_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call sfc_nst_pre_finalize()
        

    end function sfc_nst_pre_finalize_cap

    function sfc_nst_pre_run_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind
        integer, pointer :: im
        integer, pointer :: islimsk(:)
        real(kind_phys), pointer :: oro(:)
        real(kind_phys), pointer :: oro_uf(:)
        real(kind_phys), pointer :: tsfc(:)
        real(kind_phys), pointer :: tsurf(:)
        real(kind_phys), pointer :: tskin(:)

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
        

        call ccpp_field_get(cdata, 'sea_land_ice_mask', islimsk, ierr=ierr, dims=cdims, kind=ckind, index=631)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve sea_land_ice_mask from CCPP data structure')
            return
        end if
        if (kind(islimsk).ne.ckind) then
            call ccpp_error('Kind mismatch for variable sea_land_ice_mask')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'orography', oro, ierr=ierr, dims=cdims, kind=ckind, index=595)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve orography from CCPP data structure')
            return
        end if
        if (kind(oro).ne.ckind) then
            call ccpp_error('Kind mismatch for variable orography')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'orography_unfiltered', oro_uf, ierr=ierr, dims=cdims, kind=ckind, index=596)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve orography_unfiltered from CCPP data structure')
            return
        end if
        if (kind(oro_uf).ne.ckind) then
            call ccpp_error('Kind mismatch for variable orography_unfiltered')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_skin_temperature', tsfc, ierr=ierr, dims=cdims, kind=ckind, index=721)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_skin_temperature from CCPP data structure')
            return
        end if
        if (kind(tsfc).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_skin_temperature')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_skin_temperature_after_iteration', tsurf, ierr=ierr, dims=cdims, kind=ckind, index=722)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_skin_temperature_after_iteration from CCPP data structure')
            return
        end if
        if (kind(tsurf).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_skin_temperature_after_iteration')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_skin_temperature_for_nsst', tskin, ierr=ierr, dims=cdims, kind=ckind, index=723)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_skin_temperature_for_nsst from CCPP data structure')
            return
        end if
        if (kind(tskin).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_skin_temperature_for_nsst')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call sfc_nst_pre_run(im=im,islimsk=islimsk,oro=oro,oro_uf=oro_uf,tsfc=tsfc,tsurf=tsurf,tskin=tskin, &
                  errmsg=cdata%errmsg,errflg=cdata%errflg)
        ierr=cdata%errflg

    end function sfc_nst_pre_run_cap
end module sfc_nst_pre_cap

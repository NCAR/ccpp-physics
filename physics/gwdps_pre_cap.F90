
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
!! @brief Auto-generated cap module for the gwdps_pre scheme
!!
!
module gwdps_pre_cap

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_f_pointer, c_ptr, c_int32_t
    use            :: ccpp_types,                                      &
                      only: ccpp_t, CCPP_GENERIC_KIND
    use            :: ccpp_fields,                                     &
                      only: ccpp_field_get
    use            :: ccpp_errors,                                     &
                      only: ccpp_error, ccpp_debug
    use            :: gwdps_pre, &
                      only: gwdps_pre_run,gwdps_pre_init,gwdps_pre_finalize
    ! Other modules required, e.g. type definitions
    use machine, only: kind_phys

    implicit none

    private
    public :: gwdps_pre_run_cap,gwdps_pre_init_cap,gwdps_pre_finalize_cap

    contains


    function gwdps_pre_run_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind
        integer, pointer :: im
        integer, pointer :: nmtvr
        real(kind_phys), pointer :: mntvar(:,:)
        real(kind_phys), pointer :: hprime(:)
        real(kind_phys), pointer :: oc(:)
        real(kind_phys), pointer :: oa4(:,:)
        real(kind_phys), pointer :: clx(:,:)
        real(kind_phys), pointer :: theta(:)
        real(kind_phys), pointer :: sigma(:)
        real(kind_phys), pointer :: gamma(:)
        real(kind_phys), pointer :: elvmax(:)

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
        

        call ccpp_field_get(cdata, 'number_of_statistical_measures_of_subgrid_orography', nmtvr, ierr=ierr, kind=ckind, index=581)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve number_of_statistical_measures_of_subgrid_orography from CCPP data structure')
            return
        end if
        if (kind(nmtvr).ne.ckind) then
            call ccpp_error('Kind mismatch for variable number_of_statistical_measures_of_subgrid_orography')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'statistical_measures_of_subgrid_orography', mntvar, ierr=ierr, dims=cdims, kind=ckind, index=670)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve statistical_measures_of_subgrid_orography from CCPP data structure')
            return
        end if
        if (kind(mntvar).ne.ckind) then
            call ccpp_error('Kind mismatch for variable statistical_measures_of_subgrid_orography')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'standard_deviation_of_subgrid_orography', hprime, ierr=ierr, dims=cdims, kind=ckind, index=669)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve standard_deviation_of_subgrid_orography from CCPP data structure')
            return
        end if
        if (kind(hprime).ne.ckind) then
            call ccpp_error('Kind mismatch for variable standard_deviation_of_subgrid_orography')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'convexity_of_subgrid_orography', oc, ierr=ierr, dims=cdims, kind=ckind, index=133)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve convexity_of_subgrid_orography from CCPP data structure')
            return
        end if
        if (kind(oc).ne.ckind) then
            call ccpp_error('Kind mismatch for variable convexity_of_subgrid_orography')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'asymmetry_of_subgrid_orography', oa4, ierr=ierr, dims=cdims, kind=ckind, index=65)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve asymmetry_of_subgrid_orography from CCPP data structure')
            return
        end if
        if (kind(oa4).ne.ckind) then
            call ccpp_error('Kind mismatch for variable asymmetry_of_subgrid_orography')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'fraction_of_grid_box_with_subgrid_orography_higher_than_critical_height', clx, ierr=ierr, dims=cdims, kind=ckind, index=342)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve fraction_of_grid_box_with_subgrid_orography_higher_than_critical_height from CCPP data structure')
            return
        end if
        if (kind(clx).ne.ckind) then
            call ccpp_error('Kind mismatch for variable fraction_of_grid_box_with_subgrid_orography_higher_than_critical_height')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'angle_from_east_of_maximum_subgrid_orographic_variations', theta, ierr=ierr, dims=cdims, kind=ckind, index=60)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve angle_from_east_of_maximum_subgrid_orographic_variations from CCPP data structure')
            return
        end if
        if (kind(theta).ne.ckind) then
            call ccpp_error('Kind mismatch for variable angle_from_east_of_maximum_subgrid_orographic_variations')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'slope_of_subgrid_orography', sigma, ierr=ierr, dims=cdims, kind=ckind, index=647)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve slope_of_subgrid_orography from CCPP data structure')
            return
        end if
        if (kind(sigma).ne.ckind) then
            call ccpp_error('Kind mismatch for variable slope_of_subgrid_orography')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'anisotropy_of_subgrid_orography', gamma, ierr=ierr, dims=cdims, kind=ckind, index=61)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve anisotropy_of_subgrid_orography from CCPP data structure')
            return
        end if
        if (kind(gamma).ne.ckind) then
            call ccpp_error('Kind mismatch for variable anisotropy_of_subgrid_orography')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'maximum_subgrid_orography', elvmax, ierr=ierr, dims=cdims, kind=ckind, index=507)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve maximum_subgrid_orography from CCPP data structure')
            return
        end if
        if (kind(elvmax).ne.ckind) then
            call ccpp_error('Kind mismatch for variable maximum_subgrid_orography')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call gwdps_pre_run(im=im,nmtvr=nmtvr,mntvar=mntvar,hprime=hprime,oc=oc,oa4=oa4,clx=clx,theta=theta, &
                  sigma=sigma,gamma=gamma,elvmax=elvmax,errmsg=cdata%errmsg,errflg=cdata%errflg)
        ierr=cdata%errflg

    end function gwdps_pre_run_cap

    function gwdps_pre_init_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call gwdps_pre_init()
        

    end function gwdps_pre_init_cap

    function gwdps_pre_finalize_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call gwdps_pre_finalize()
        

    end function gwdps_pre_finalize_cap
end module gwdps_pre_cap

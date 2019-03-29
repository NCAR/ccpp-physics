
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
!! @brief Auto-generated cap module for the cnvc90 scheme
!!
!
module cnvc90_cap

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_f_pointer, c_ptr, c_int32_t
    use            :: ccpp_types,                                      &
                      only: ccpp_t, CCPP_GENERIC_KIND
    use            :: ccpp_fields,                                     &
                      only: ccpp_field_get
    use            :: ccpp_errors,                                     &
                      only: ccpp_error, ccpp_debug
    use            :: cnvc90, &
                      only: cnvc90_finalize,cnvc90_init,cnvc90_run
    ! Other modules required, e.g. type definitions
    use machine, only: kind_phys

    implicit none

    private
    public :: cnvc90_finalize_cap,cnvc90_init_cap,cnvc90_run_cap

    contains


    function cnvc90_finalize_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call cnvc90_finalize()
        

    end function cnvc90_finalize_cap

    function cnvc90_init_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call cnvc90_init()
        

    end function cnvc90_init_cap

    function cnvc90_run_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind
        real(kind_phys), pointer :: clstp
        integer, pointer :: im
        integer, pointer :: ix
        real(kind_phys), pointer :: rn(:)
        integer, pointer :: kbot(:)
        integer, pointer :: ktop(:)
        integer, pointer :: km
        real(kind_phys), pointer :: prsi(:,:)
        real(kind_phys), pointer :: acv(:)
        real(kind_phys), pointer :: acvb(:)
        real(kind_phys), pointer :: acvt(:)
        real(kind_phys), pointer :: cv(:)
        real(kind_phys), pointer :: cvb(:)
        real(kind_phys), pointer :: cvt(:)

        ierr = 0

        call c_f_pointer(ptr, cdata)


        call ccpp_field_get(cdata, 'convective_cloud_switch', clstp, ierr=ierr, kind=ckind, index=126)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve convective_cloud_switch from CCPP data structure')
            return
        end if
        if (kind(clstp).ne.ckind) then
            call ccpp_error('Kind mismatch for variable convective_cloud_switch')
            ierr = 1
            return
        end if
#endif
        

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
        

        call ccpp_field_get(cdata, 'horizontal_dimension', ix, ierr=ierr, kind=ckind, index=364)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve horizontal_dimension from CCPP data structure')
            return
        end if
        if (kind(ix).ne.ckind) then
            call ccpp_error('Kind mismatch for variable horizontal_dimension')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'lwe_thickness_of_convective_precipitation_amount_on_dynamics_timestep', rn, ierr=ierr, dims=cdims, kind=ckind, index=478)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve lwe_thickness_of_convective_precipitation_amount_on_dynamics_timestep from CCPP data structure')
            return
        end if
        if (kind(rn).ne.ckind) then
            call ccpp_error('Kind mismatch for variable lwe_thickness_of_convective_precipitation_amount_on_dynamics_timestep')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'vertical_index_at_cloud_base', kbot, ierr=ierr, dims=cdims, kind=ckind, index=820)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve vertical_index_at_cloud_base from CCPP data structure')
            return
        end if
        if (kind(kbot).ne.ckind) then
            call ccpp_error('Kind mismatch for variable vertical_index_at_cloud_base')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'vertical_index_at_cloud_top', ktop, ierr=ierr, dims=cdims, kind=ckind, index=821)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve vertical_index_at_cloud_top from CCPP data structure')
            return
        end if
        if (kind(ktop).ne.ckind) then
            call ccpp_error('Kind mismatch for variable vertical_index_at_cloud_top')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'vertical_dimension', km, ierr=ierr, kind=ckind, index=817)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve vertical_dimension from CCPP data structure')
            return
        end if
        if (kind(km).ne.ckind) then
            call ccpp_error('Kind mismatch for variable vertical_dimension')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'air_pressure_at_interface', prsi, ierr=ierr, dims=cdims, kind=ckind, index=45)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve air_pressure_at_interface from CCPP data structure')
            return
        end if
        if (kind(prsi).ne.ckind) then
            call ccpp_error('Kind mismatch for variable air_pressure_at_interface')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'accumulated_lwe_thickness_of_convective_precipitation_amount_cnvc90', acv, ierr=ierr, dims=cdims, kind=ckind, index=21)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve accumulated_lwe_thickness_of_convective_precipitation_amount_cnvc90 from CCPP data structure')
            return
        end if
        if (kind(acv).ne.ckind) then
            call ccpp_error('Kind mismatch for variable accumulated_lwe_thickness_of_convective_precipitation_amount_cnvc90')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'smallest_cloud_base_vertical_index_encountered_thus_far', acvb, ierr=ierr, dims=cdims, kind=ckind, index=648)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve smallest_cloud_base_vertical_index_encountered_thus_far from CCPP data structure')
            return
        end if
        if (kind(acvb).ne.ckind) then
            call ccpp_error('Kind mismatch for variable smallest_cloud_base_vertical_index_encountered_thus_far')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'largest_cloud_top_vertical_index_encountered_thus_far', acvt, ierr=ierr, dims=cdims, kind=ckind, index=457)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve largest_cloud_top_vertical_index_encountered_thus_far from CCPP data structure')
            return
        end if
        if (kind(acvt).ne.ckind) then
            call ccpp_error('Kind mismatch for variable largest_cloud_top_vertical_index_encountered_thus_far')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'fraction_of_convective_cloud', cv, ierr=ierr, dims=cdims, kind=ckind, index=341)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve fraction_of_convective_cloud from CCPP data structure')
            return
        end if
        if (kind(cv).ne.ckind) then
            call ccpp_error('Kind mismatch for variable fraction_of_convective_cloud')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'pressure_at_bottom_of_convective_cloud', cvb, ierr=ierr, dims=cdims, kind=ckind, index=609)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve pressure_at_bottom_of_convective_cloud from CCPP data structure')
            return
        end if
        if (kind(cvb).ne.ckind) then
            call ccpp_error('Kind mismatch for variable pressure_at_bottom_of_convective_cloud')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'pressure_at_top_of_convective_cloud', cvt, ierr=ierr, dims=cdims, kind=ckind, index=610)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve pressure_at_top_of_convective_cloud from CCPP data structure')
            return
        end if
        if (kind(cvt).ne.ckind) then
            call ccpp_error('Kind mismatch for variable pressure_at_top_of_convective_cloud')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call cnvc90_run(clstp=clstp,im=im,ix=ix,rn=rn,kbot=kbot,ktop=ktop,km=km,prsi=prsi,acv=acv, &
                  acvb=acvb,acvt=acvt,cv=cv,cvb=cvb,cvt=cvt,errmsg=cdata%errmsg,errflg=cdata%errflg)
        ierr=cdata%errflg

    end function cnvc90_run_cap
end module cnvc90_cap

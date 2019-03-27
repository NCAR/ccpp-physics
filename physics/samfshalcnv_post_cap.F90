
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
!! @brief Auto-generated cap module for the samfshalcnv_post scheme
!!
!
module samfshalcnv_post_cap

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_f_pointer, c_ptr, c_int32_t
    use            :: ccpp_types,                                      &
                      only: ccpp_t, CCPP_GENERIC_KIND
    use            :: ccpp_fields,                                     &
                      only: ccpp_field_get
    use            :: ccpp_errors,                                     &
                      only: ccpp_error, ccpp_debug
    use            :: samfshalcnv_post, &
                      only: samfshalcnv_post_run,samfshalcnv_post_finalize,samfshalcnv_post_init
    ! Other modules required, e.g. type definitions
    use machine, only: kind_phys

    implicit none

    private
    public :: samfshalcnv_post_run_cap,samfshalcnv_post_finalize_cap,samfshalcnv_post_init_cap

    contains


    function samfshalcnv_post_run_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind
        integer, pointer :: im
        integer, pointer :: levs
        logical, pointer :: lssav
        logical, pointer :: shcnvcw
        real(kind_phys), pointer :: frain
        real(kind_phys), pointer :: rain1(:)
        integer, pointer :: npdf3d
        integer, pointer :: num_p3d
        integer, pointer :: ncnvcld3d
        real(kind_phys), pointer :: cnvc(:,:)
        real(kind_phys), pointer :: cnvw(:,:)
        real(kind_phys), pointer :: rainc(:)
        real(kind_phys), pointer :: cnvprcp(:)
        real(kind_phys), pointer :: cnvprcpb(:)
        real(kind_phys), pointer :: cnvw_phy_f3d(:,:)
        real(kind_phys), pointer :: cnvc_phy_f3d(:,:)

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
        

        call ccpp_field_get(cdata, 'vertical_dimension', levs, ierr=ierr, kind=ckind, index=817)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve vertical_dimension from CCPP data structure')
            return
        end if
        if (kind(levs).ne.ckind) then
            call ccpp_error('Kind mismatch for variable vertical_dimension')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'flag_diagnostics', lssav, ierr=ierr, kind=ckind, index=263)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_diagnostics from CCPP data structure')
            return
        end if
        if (kind(lssav).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_diagnostics')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'flag_shallow_convective_cloud', shcnvcw, ierr=ierr, kind=ckind, index=333)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_shallow_convective_cloud from CCPP data structure')
            return
        end if
        if (kind(shcnvcw).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_shallow_convective_cloud')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'dynamics_to_physics_timestep_ratio', frain, ierr=ierr, kind=ckind, index=236)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve dynamics_to_physics_timestep_ratio from CCPP data structure')
            return
        end if
        if (kind(frain).ne.ckind) then
            call ccpp_error('Kind mismatch for variable dynamics_to_physics_timestep_ratio')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'lwe_thickness_of_shallow_convective_precipitation_amount', rain1, ierr=ierr, dims=cdims, kind=ckind, index=491)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve lwe_thickness_of_shallow_convective_precipitation_amount from CCPP data structure')
            return
        end if
        if (kind(rain1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable lwe_thickness_of_shallow_convective_precipitation_amount')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'number_of_3d_arrays_associated_with_pdf-based_clouds', npdf3d, ierr=ierr, kind=ckind, index=571)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve number_of_3d_arrays_associated_with_pdf-based_clouds from CCPP data structure')
            return
        end if
        if (kind(npdf3d).ne.ckind) then
            call ccpp_error('Kind mismatch for variable number_of_3d_arrays_associated_with_pdf-based_clouds')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'array_dimension_of_3d_arrays_for_microphysics', num_p3d, ierr=ierr, kind=ckind, index=63)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve array_dimension_of_3d_arrays_for_microphysics from CCPP data structure')
            return
        end if
        if (kind(num_p3d).ne.ckind) then
            call ccpp_error('Kind mismatch for variable array_dimension_of_3d_arrays_for_microphysics')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'number_of_convective_3d_cloud_fields', ncnvcld3d, ierr=ierr, kind=ckind, index=575)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve number_of_convective_3d_cloud_fields from CCPP data structure')
            return
        end if
        if (kind(ncnvcld3d).ne.ckind) then
            call ccpp_error('Kind mismatch for variable number_of_convective_3d_cloud_fields')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'convective_cloud_cover', cnvc, ierr=ierr, dims=cdims, kind=ckind, index=123)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve convective_cloud_cover from CCPP data structure')
            return
        end if
        if (kind(cnvc).ne.ckind) then
            call ccpp_error('Kind mismatch for variable convective_cloud_cover')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'convective_cloud_water_mixing_ratio', cnvw, ierr=ierr, dims=cdims, kind=ckind, index=128)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve convective_cloud_water_mixing_ratio from CCPP data structure')
            return
        end if
        if (kind(cnvw).ne.ckind) then
            call ccpp_error('Kind mismatch for variable convective_cloud_water_mixing_ratio')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'lwe_thickness_of_convective_precipitation_amount_on_dynamics_timestep', rainc, ierr=ierr, dims=cdims, kind=ckind, index=478)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve lwe_thickness_of_convective_precipitation_amount_on_dynamics_timestep from CCPP data structure')
            return
        end if
        if (kind(rainc).ne.ckind) then
            call ccpp_error('Kind mismatch for variable lwe_thickness_of_convective_precipitation_amount_on_dynamics_timestep')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cumulative_lwe_thickness_of_convective_precipitation_amount', cnvprcp, ierr=ierr, dims=cdims, kind=ckind, index=174)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cumulative_lwe_thickness_of_convective_precipitation_amount from CCPP data structure')
            return
        end if
        if (kind(cnvprcp).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cumulative_lwe_thickness_of_convective_precipitation_amount')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cumulative_lwe_thickness_of_convective_precipitation_amount_in_bucket', cnvprcpb, ierr=ierr, dims=cdims, kind=ckind, index=175)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cumulative_lwe_thickness_of_convective_precipitation_amount_in_bucket from CCPP data structure')
            return
        end if
        if (kind(cnvprcpb).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cumulative_lwe_thickness_of_convective_precipitation_amount_in_bucket')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'convective_cloud_water_mixing_ratio_in_phy_f3d', cnvw_phy_f3d, ierr=ierr, dims=cdims, kind=ckind, index=129)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve convective_cloud_water_mixing_ratio_in_phy_f3d from CCPP data structure')
            return
        end if
        if (kind(cnvw_phy_f3d).ne.ckind) then
            call ccpp_error('Kind mismatch for variable convective_cloud_water_mixing_ratio_in_phy_f3d')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'convective_cloud_cover_in_phy_f3d', cnvc_phy_f3d, ierr=ierr, dims=cdims, kind=ckind, index=124)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve convective_cloud_cover_in_phy_f3d from CCPP data structure')
            return
        end if
        if (kind(cnvc_phy_f3d).ne.ckind) then
            call ccpp_error('Kind mismatch for variable convective_cloud_cover_in_phy_f3d')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call samfshalcnv_post_run(im=im,levs=levs,lssav=lssav,shcnvcw=shcnvcw,frain=frain,rain1=rain1,npdf3d=npdf3d, &
                  num_p3d=num_p3d,ncnvcld3d=ncnvcld3d,cnvc=cnvc,cnvw=cnvw,rainc=rainc,cnvprcp=cnvprcp, &
                  cnvprcpb=cnvprcpb,cnvw_phy_f3d=cnvw_phy_f3d,cnvc_phy_f3d=cnvc_phy_f3d,errmsg=cdata%errmsg, &
                  errflg=cdata%errflg)
        ierr=cdata%errflg

    end function samfshalcnv_post_run_cap

    function samfshalcnv_post_finalize_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call samfshalcnv_post_finalize()
        

    end function samfshalcnv_post_finalize_cap

    function samfshalcnv_post_init_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call samfshalcnv_post_init()
        

    end function samfshalcnv_post_init_cap
end module samfshalcnv_post_cap

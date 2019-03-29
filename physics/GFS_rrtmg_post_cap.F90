
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
!! @brief Auto-generated cap module for the GFS_rrtmg_post scheme
!!
!
module GFS_rrtmg_post_cap

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_f_pointer, c_ptr, c_int32_t
    use            :: ccpp_types,                                      &
                      only: ccpp_t, CCPP_GENERIC_KIND
    use            :: ccpp_fields,                                     &
                      only: ccpp_field_get
    use            :: ccpp_errors,                                     &
                      only: ccpp_error, ccpp_debug
    use            :: GFS_rrtmg_post, &
                      only: GFS_rrtmg_post_init,GFS_rrtmg_post_run,GFS_rrtmg_post_finalize
    ! Other modules required, e.g. type definitions
    use GFS_typedefs, only: GFS_control_type
    use GFS_typedefs, only: GFS_grid_type
    use GFS_typedefs, only: GFS_diag_type
    use GFS_typedefs, only: GFS_radtend_type
    use GFS_typedefs, only: GFS_statein_type
    use GFS_typedefs, only: GFS_coupling_type
    use GFS_typedefs, only: cmpfsw_type
    use machine, only: kind_phys

    implicit none

    private
    public :: GFS_rrtmg_post_init_cap,GFS_rrtmg_post_run_cap,GFS_rrtmg_post_finalize_cap

    contains


    function GFS_rrtmg_post_init_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call GFS_rrtmg_post_init()
        

    end function GFS_rrtmg_post_init_cap

    function GFS_rrtmg_post_run_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind
        type(GFS_control_type), pointer     :: Model
        type(GFS_grid_type), pointer     :: Grid
        type(GFS_diag_type), pointer     :: Diag
        type(GFS_radtend_type), pointer     :: Radtend
        type(GFS_statein_type), pointer     :: Statein
        type(GFS_coupling_type), pointer     :: Coupling
        type(cmpfsw_type), pointer     :: scmpsw(:)
        integer, pointer :: im
        integer, pointer :: lm
        integer, pointer :: ltp
        integer, pointer :: kt
        integer, pointer :: kb
        integer, pointer :: kd
        real(kind_phys), pointer :: raddt
        real(kind_phys), pointer :: aerodp(:,:)
        real(kind_phys), pointer :: cldsa(:,:)
        integer, pointer :: mtopa(:,:)
        integer, pointer :: mbota(:,:)
        real(kind_phys), pointer :: clouds1(:,:)
        real(kind_phys), pointer :: cldtaulw(:,:)
        real(kind_phys), pointer :: cldtausw(:,:)

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

        call ccpp_field_get(cdata, 'components_of_surface_downward_shortwave_fluxes', cptr, ierr=ierr, dims=cdims, kind=ckind, index=121)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve components_of_surface_downward_shortwave_fluxes from CCPP data structure')
            return
        end if
        if (ckind.ne.CCPP_GENERIC_KIND) then
            call ccpp_error('Kind mismatch for variable components_of_surface_downward_shortwave_fluxes')
            ierr = 1
            return
        end if
#endif
        call c_f_pointer(cptr, scmpsw, cdims)
        deallocate(cdims)
        

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
        

        call ccpp_field_get(cdata, 'vertical_layer_dimension_for_radiation', lm, ierr=ierr, kind=ckind, index=826)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve vertical_layer_dimension_for_radiation from CCPP data structure')
            return
        end if
        if (kind(lm).ne.ckind) then
            call ccpp_error('Kind mismatch for variable vertical_layer_dimension_for_radiation')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'extra_top_layer', ltp, ierr=ierr, kind=ckind, index=257)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve extra_top_layer from CCPP data structure')
            return
        end if
        if (kind(ltp).ne.ckind) then
            call ccpp_error('Kind mismatch for variable extra_top_layer')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'vertical_index_difference_between_layer_and_upper_bound', kt, ierr=ierr, kind=ckind, index=825)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve vertical_index_difference_between_layer_and_upper_bound from CCPP data structure')
            return
        end if
        if (kind(kt).ne.ckind) then
            call ccpp_error('Kind mismatch for variable vertical_index_difference_between_layer_and_upper_bound')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'vertical_index_difference_between_layer_and_lower_bound', kb, ierr=ierr, kind=ckind, index=824)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve vertical_index_difference_between_layer_and_lower_bound from CCPP data structure')
            return
        end if
        if (kind(kb).ne.ckind) then
            call ccpp_error('Kind mismatch for variable vertical_index_difference_between_layer_and_lower_bound')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'vertical_index_difference_between_inout_and_local', kd, ierr=ierr, kind=ckind, index=823)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve vertical_index_difference_between_inout_and_local from CCPP data structure')
            return
        end if
        if (kind(kd).ne.ckind) then
            call ccpp_error('Kind mismatch for variable vertical_index_difference_between_inout_and_local')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'time_step_for_radiation', raddt, ierr=ierr, kind=ckind, index=794)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve time_step_for_radiation from CCPP data structure')
            return
        end if
        if (kind(raddt).ne.ckind) then
            call ccpp_error('Kind mismatch for variable time_step_for_radiation')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'atmosphere_optical_thickness_due_to_ambient_aerosol_particles', aerodp, ierr=ierr, dims=cdims, kind=ckind, index=75)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve atmosphere_optical_thickness_due_to_ambient_aerosol_particles from CCPP data structure')
            return
        end if
        if (kind(aerodp).ne.ckind) then
            call ccpp_error('Kind mismatch for variable atmosphere_optical_thickness_due_to_ambient_aerosol_particles')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cloud_area_fraction_for_radiation', cldsa, ierr=ierr, dims=cdims, kind=ckind, index=87)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cloud_area_fraction_for_radiation from CCPP data structure')
            return
        end if
        if (kind(cldsa).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cloud_area_fraction_for_radiation')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'model_layer_number_at_cloud_top', mtopa, ierr=ierr, dims=cdims, kind=ckind, index=552)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve model_layer_number_at_cloud_top from CCPP data structure')
            return
        end if
        if (kind(mtopa).ne.ckind) then
            call ccpp_error('Kind mismatch for variable model_layer_number_at_cloud_top')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'model_layer_number_at_cloud_base', mbota, ierr=ierr, dims=cdims, kind=ckind, index=551)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve model_layer_number_at_cloud_base from CCPP data structure')
            return
        end if
        if (kind(mbota).ne.ckind) then
            call ccpp_error('Kind mismatch for variable model_layer_number_at_cloud_base')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'total_cloud_fraction', clouds1, ierr=ierr, dims=cdims, kind=ckind, index=799)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve total_cloud_fraction from CCPP data structure')
            return
        end if
        if (kind(clouds1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable total_cloud_fraction')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cloud_optical_depth_layers_at_10mu_band', cldtaulw, ierr=ierr, dims=cdims, kind=ckind, index=104)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cloud_optical_depth_layers_at_10mu_band from CCPP data structure')
            return
        end if
        if (kind(cldtaulw).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cloud_optical_depth_layers_at_10mu_band')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cloud_optical_depth_layers_at_0.55mu_band', cldtausw, ierr=ierr, dims=cdims, kind=ckind, index=103)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cloud_optical_depth_layers_at_0.55mu_band from CCPP data structure')
            return
        end if
        if (kind(cldtausw).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cloud_optical_depth_layers_at_0.55mu_band')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call GFS_rrtmg_post_run(Model=Model,Grid=Grid,Diag=Diag,Radtend=Radtend,Statein=Statein,Coupling=Coupling, &
                  scmpsw=scmpsw,im=im,lm=lm,ltp=ltp,kt=kt,kb=kb,kd=kd,raddt=raddt,aerodp=aerodp, &
                  cldsa=cldsa,mtopa=mtopa,mbota=mbota,clouds1=clouds1,cldtaulw=cldtaulw,cldtausw=cldtausw, &
                  errmsg=cdata%errmsg,errflg=cdata%errflg)
        ierr=cdata%errflg

    end function GFS_rrtmg_post_run_cap

    function GFS_rrtmg_post_finalize_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call GFS_rrtmg_post_finalize()
        

    end function GFS_rrtmg_post_finalize_cap
end module GFS_rrtmg_post_cap

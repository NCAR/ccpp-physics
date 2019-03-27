
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
!! @brief Auto-generated cap module for the GFS_DCNV_generic_post scheme
!!
!
module GFS_DCNV_generic_post_cap

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_f_pointer, c_ptr, c_int32_t
    use            :: ccpp_types,                                      &
                      only: ccpp_t, CCPP_GENERIC_KIND
    use            :: ccpp_fields,                                     &
                      only: ccpp_field_get
    use            :: ccpp_errors,                                     &
                      only: ccpp_error, ccpp_debug
    use            :: GFS_DCNV_generic_post, &
                      only: GFS_DCNV_generic_post_init,GFS_DCNV_generic_post_run,GFS_DCNV_generic_post_finalize
    ! Other modules required, e.g. type definitions
    use machine, only: kind_phys

    implicit none

    private
    public :: GFS_DCNV_generic_post_init_cap,GFS_DCNV_generic_post_run_cap,GFS_DCNV_generic_post_finalize_cap

    contains


    function GFS_DCNV_generic_post_init_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call GFS_DCNV_generic_post_init()
        

    end function GFS_DCNV_generic_post_init_cap

    function GFS_DCNV_generic_post_run_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind
        integer, pointer :: im
        integer, pointer :: levs
        logical, pointer :: lssav
        logical, pointer :: ldiag3d
        logical, pointer :: lgocart
        logical, pointer :: ras
        logical, pointer :: cscnv
        real(kind_phys), pointer :: frain
        real(kind_phys), pointer :: rain1(:)
        real(kind_phys), pointer :: dtf
        real(kind_phys), pointer :: cld1d(:)
        real(kind_phys), pointer :: save_u(:,:)
        real(kind_phys), pointer :: save_v(:,:)
        real(kind_phys), pointer :: save_t(:,:)
        real(kind_phys), pointer :: save_qv(:,:)
        real(kind_phys), pointer :: gu0(:,:)
        real(kind_phys), pointer :: gv0(:,:)
        real(kind_phys), pointer :: gt0(:,:)
        real(kind_phys), pointer :: gq0_water_vapor(:,:)
        real(kind_phys), pointer :: ud_mf(:,:)
        real(kind_phys), pointer :: dd_mf(:,:)
        real(kind_phys), pointer :: dt_mf(:,:)
        real(kind_phys), pointer :: con_g
        real(kind_phys), pointer :: clw_ice(:,:)
        real(kind_phys), pointer :: clw_liquid(:,:)
        integer, pointer :: npdf3d
        integer, pointer :: num_p3d
        integer, pointer :: ncnvcld3d
        real(kind_phys), pointer :: rainc(:)
        real(kind_phys), pointer :: cldwrk(:)
        real(kind_phys), pointer :: cnvprcp(:)
        real(kind_phys), pointer :: cnvprcpb(:)
        real(kind_phys), pointer :: dt3dt(:,:)
        real(kind_phys), pointer :: dq3dt(:,:)
        real(kind_phys), pointer :: du3dt(:,:)
        real(kind_phys), pointer :: dv3dt(:,:)
        real(kind_phys), pointer :: upd_mf(:,:)
        real(kind_phys), pointer :: dwn_mf(:,:)
        real(kind_phys), pointer :: det_mf(:,:)
        real(kind_phys), pointer :: dqdti(:,:)
        real(kind_phys), pointer :: cnvqci(:,:)
        real(kind_phys), pointer :: upd_mfi(:,:)
        real(kind_phys), pointer :: dwn_mfi(:,:)
        real(kind_phys), pointer :: det_mfi(:,:)
        real(kind_phys), pointer :: cnvw(:,:)
        real(kind_phys), pointer :: cnvc(:,:)
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
        

        call ccpp_field_get(cdata, 'flag_diagnostics_3D', ldiag3d, ierr=ierr, kind=ckind, index=264)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_diagnostics_3D from CCPP data structure')
            return
        end if
        if (kind(ldiag3d).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_diagnostics_3D')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'flag_gocart', lgocart, ierr=ierr, kind=ckind, index=329)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_gocart from CCPP data structure')
            return
        end if
        if (kind(lgocart).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_gocart')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'flag_for_ras_deep_convection', ras, ierr=ierr, kind=ckind, index=307)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_ras_deep_convection from CCPP data structure')
            return
        end if
        if (kind(ras).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_ras_deep_convection')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'flag_for_Chikira_Sugiyama_deep_convection', cscnv, ierr=ierr, kind=ckind, index=269)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_Chikira_Sugiyama_deep_convection from CCPP data structure')
            return
        end if
        if (kind(cscnv).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_Chikira_Sugiyama_deep_convection')
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
        

        call ccpp_field_get(cdata, 'lwe_thickness_of_deep_convective_precipitation_amount', rain1, ierr=ierr, dims=cdims, kind=ckind, index=479)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve lwe_thickness_of_deep_convective_precipitation_amount from CCPP data structure')
            return
        end if
        if (kind(rain1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable lwe_thickness_of_deep_convective_precipitation_amount')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'time_step_for_dynamics', dtf, ierr=ierr, kind=ckind, index=792)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve time_step_for_dynamics from CCPP data structure')
            return
        end if
        if (kind(dtf).ne.ckind) then
            call ccpp_error('Kind mismatch for variable time_step_for_dynamics')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'cloud_work_function', cld1d, ierr=ierr, dims=cdims, kind=ckind, index=111)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cloud_work_function from CCPP data structure')
            return
        end if
        if (kind(cld1d).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cloud_work_function')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'x_wind_save', save_u, ierr=ierr, dims=cdims, kind=ckind, index=876)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve x_wind_save from CCPP data structure')
            return
        end if
        if (kind(save_u).ne.ckind) then
            call ccpp_error('Kind mismatch for variable x_wind_save')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'y_wind_save', save_v, ierr=ierr, dims=cdims, kind=ckind, index=883)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve y_wind_save from CCPP data structure')
            return
        end if
        if (kind(save_v).ne.ckind) then
            call ccpp_error('Kind mismatch for variable y_wind_save')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'air_temperature_save', save_t, ierr=ierr, dims=cdims, kind=ckind, index=57)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve air_temperature_save from CCPP data structure')
            return
        end if
        if (kind(save_t).ne.ckind) then
            call ccpp_error('Kind mismatch for variable air_temperature_save')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'water_vapor_specific_humidity_save', save_qv, ierr=ierr, dims=cdims, kind=ckind, index=858)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve water_vapor_specific_humidity_save from CCPP data structure')
            return
        end if
        if (kind(save_qv).ne.ckind) then
            call ccpp_error('Kind mismatch for variable water_vapor_specific_humidity_save')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'x_wind_updated_by_physics', gu0, ierr=ierr, dims=cdims, kind=ckind, index=877)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve x_wind_updated_by_physics from CCPP data structure')
            return
        end if
        if (kind(gu0).ne.ckind) then
            call ccpp_error('Kind mismatch for variable x_wind_updated_by_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'y_wind_updated_by_physics', gv0, ierr=ierr, dims=cdims, kind=ckind, index=884)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve y_wind_updated_by_physics from CCPP data structure')
            return
        end if
        if (kind(gv0).ne.ckind) then
            call ccpp_error('Kind mismatch for variable y_wind_updated_by_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'air_temperature_updated_by_physics', gt0, ierr=ierr, dims=cdims, kind=ckind, index=59)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve air_temperature_updated_by_physics from CCPP data structure')
            return
        end if
        if (kind(gt0).ne.ckind) then
            call ccpp_error('Kind mismatch for variable air_temperature_updated_by_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'water_vapor_specific_humidity_updated_by_physics', gq0_water_vapor, ierr=ierr, dims=cdims, kind=ckind, index=860)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve water_vapor_specific_humidity_updated_by_physics from CCPP data structure')
            return
        end if
        if (kind(gq0_water_vapor).ne.ckind) then
            call ccpp_error('Kind mismatch for variable water_vapor_specific_humidity_updated_by_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'instantaneous_atmosphere_updraft_convective_mass_flux', ud_mf, ierr=ierr, dims=cdims, kind=ckind, index=403)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve instantaneous_atmosphere_updraft_convective_mass_flux from CCPP data structure')
            return
        end if
        if (kind(ud_mf).ne.ckind) then
            call ccpp_error('Kind mismatch for variable instantaneous_atmosphere_updraft_convective_mass_flux')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'instantaneous_atmosphere_downdraft_convective_mass_flux', dd_mf, ierr=ierr, dims=cdims, kind=ckind, index=401)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve instantaneous_atmosphere_downdraft_convective_mass_flux from CCPP data structure')
            return
        end if
        if (kind(dd_mf).ne.ckind) then
            call ccpp_error('Kind mismatch for variable instantaneous_atmosphere_downdraft_convective_mass_flux')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'instantaneous_atmosphere_detrainment_convective_mass_flux', dt_mf, ierr=ierr, dims=cdims, kind=ckind, index=399)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve instantaneous_atmosphere_detrainment_convective_mass_flux from CCPP data structure')
            return
        end if
        if (kind(dt_mf).ne.ckind) then
            call ccpp_error('Kind mismatch for variable instantaneous_atmosphere_detrainment_convective_mass_flux')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'gravitational_acceleration', con_g, ierr=ierr, kind=ckind, index=355)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve gravitational_acceleration from CCPP data structure')
            return
        end if
        if (kind(con_g).ne.ckind) then
            call ccpp_error('Kind mismatch for variable gravitational_acceleration')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'ice_water_mixing_ratio_convective_transport_tracer', clw_ice, ierr=ierr, dims=cdims, kind=ckind, index=374)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve ice_water_mixing_ratio_convective_transport_tracer from CCPP data structure')
            return
        end if
        if (kind(clw_ice).ne.ckind) then
            call ccpp_error('Kind mismatch for variable ice_water_mixing_ratio_convective_transport_tracer')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cloud_condensed_water_mixing_ratio_convective_transport_tracer', clw_liquid, ierr=ierr, dims=cdims, kind=ckind, index=93)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cloud_condensed_water_mixing_ratio_convective_transport_tracer from CCPP data structure')
            return
        end if
        if (kind(clw_liquid).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cloud_condensed_water_mixing_ratio_convective_transport_tracer')
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
        

        call ccpp_field_get(cdata, 'cumulative_cloud_work_function', cldwrk, ierr=ierr, dims=cdims, kind=ckind, index=173)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cumulative_cloud_work_function from CCPP data structure')
            return
        end if
        if (kind(cldwrk).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cumulative_cloud_work_function')
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
        

        call ccpp_field_get(cdata, 'cumulative_change_in_temperature_due_to_deep_convection', dt3dt, ierr=ierr, dims=cdims, kind=ckind, index=156)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cumulative_change_in_temperature_due_to_deep_convection from CCPP data structure')
            return
        end if
        if (kind(dt3dt).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cumulative_change_in_temperature_due_to_deep_convection')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cumulative_change_in_water_vapor_specific_humidity_due_to_deep_convection', dq3dt, ierr=ierr, dims=cdims, kind=ckind, index=162)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cumulative_change_in_water_vapor_specific_humidity_due_to_deep_convection from CCPP data structure')
            return
        end if
        if (kind(dq3dt).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cumulative_change_in_water_vapor_specific_humidity_due_to_deep_convection')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cumulative_change_in_x_wind_due_to_deep_convection', du3dt, ierr=ierr, dims=cdims, kind=ckind, index=167)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cumulative_change_in_x_wind_due_to_deep_convection from CCPP data structure')
            return
        end if
        if (kind(du3dt).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cumulative_change_in_x_wind_due_to_deep_convection')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cumulative_change_in_y_wind_due_to_deep_convection', dv3dt, ierr=ierr, dims=cdims, kind=ckind, index=171)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cumulative_change_in_y_wind_due_to_deep_convection from CCPP data structure')
            return
        end if
        if (kind(dv3dt).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cumulative_change_in_y_wind_due_to_deep_convection')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cumulative_atmosphere_updraft_convective_mass_flux', upd_mf, ierr=ierr, dims=cdims, kind=ckind, index=148)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cumulative_atmosphere_updraft_convective_mass_flux from CCPP data structure')
            return
        end if
        if (kind(upd_mf).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cumulative_atmosphere_updraft_convective_mass_flux')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cumulative_atmosphere_downdraft_convective_mass_flux', dwn_mf, ierr=ierr, dims=cdims, kind=ckind, index=147)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cumulative_atmosphere_downdraft_convective_mass_flux from CCPP data structure')
            return
        end if
        if (kind(dwn_mf).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cumulative_atmosphere_downdraft_convective_mass_flux')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cumulative_atmosphere_detrainment_convective_mass_flux', det_mf, ierr=ierr, dims=cdims, kind=ckind, index=146)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cumulative_atmosphere_detrainment_convective_mass_flux from CCPP data structure')
            return
        end if
        if (kind(det_mf).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cumulative_atmosphere_detrainment_convective_mass_flux')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'instantaneous_water_vapor_specific_humidity_tendency_due_to_convection', dqdti, ierr=ierr, dims=cdims, kind=ckind, index=444)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve instantaneous_water_vapor_specific_humidity_tendency_due_to_convection from CCPP data structure')
            return
        end if
        if (kind(dqdti).ne.ckind) then
            call ccpp_error('Kind mismatch for variable instantaneous_water_vapor_specific_humidity_tendency_due_to_convection')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'instantaneous_deep_convective_cloud_condensate_mixing_ratio_on_dynamics_time_step', cnvqci, ierr=ierr, dims=cdims, kind=ckind, index=409)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve instantaneous_deep_convective_cloud_condensate_mixing_ratio_on_dynamics_time_step from CCPP data structure')
            return
        end if
        if (kind(cnvqci).ne.ckind) then
            call ccpp_error('Kind mismatch for variable instantaneous_deep_convective_cloud_condensate_mixing_ratio_on_dynamics_time_step')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'instantaneous_atmosphere_updraft_convective_mass_flux_on_dynamics_timestep', upd_mfi, ierr=ierr, dims=cdims, kind=ckind, index=404)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve instantaneous_atmosphere_updraft_convective_mass_flux_on_dynamics_timestep from CCPP data structure')
            return
        end if
        if (kind(upd_mfi).ne.ckind) then
            call ccpp_error('Kind mismatch for variable instantaneous_atmosphere_updraft_convective_mass_flux_on_dynamics_timestep')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'instantaneous_atmosphere_downdraft_convective_mass_flux_on_dynamics_timestep', dwn_mfi, ierr=ierr, dims=cdims, kind=ckind, index=402)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve instantaneous_atmosphere_downdraft_convective_mass_flux_on_dynamics_timestep from CCPP data structure')
            return
        end if
        if (kind(dwn_mfi).ne.ckind) then
            call ccpp_error('Kind mismatch for variable instantaneous_atmosphere_downdraft_convective_mass_flux_on_dynamics_timestep')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'instantaneous_atmosphere_detrainment_convective_mass_flux_on_dynamics_timestep', det_mfi, ierr=ierr, dims=cdims, kind=ckind, index=400)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve instantaneous_atmosphere_detrainment_convective_mass_flux_on_dynamics_timestep from CCPP data structure')
            return
        end if
        if (kind(det_mfi).ne.ckind) then
            call ccpp_error('Kind mismatch for variable instantaneous_atmosphere_detrainment_convective_mass_flux_on_dynamics_timestep')
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
        

        call GFS_DCNV_generic_post_run(im=im,levs=levs,lssav=lssav,ldiag3d=ldiag3d,lgocart=lgocart,ras=ras,cscnv=cscnv, &
                  frain=frain,rain1=rain1,dtf=dtf,cld1d=cld1d,save_u=save_u,save_v=save_v, &
                  save_t=save_t,save_qv=save_qv,gu0=gu0,gv0=gv0,gt0=gt0,gq0_water_vapor=gq0_water_vapor, &
                  ud_mf=ud_mf,dd_mf=dd_mf,dt_mf=dt_mf,con_g=con_g,clw_ice=clw_ice,clw_liquid=clw_liquid, &
                  npdf3d=npdf3d,num_p3d=num_p3d,ncnvcld3d=ncnvcld3d,rainc=rainc,cldwrk=cldwrk, &
                  cnvprcp=cnvprcp,cnvprcpb=cnvprcpb,dt3dt=dt3dt,dq3dt=dq3dt,du3dt=du3dt,dv3dt=dv3dt, &
                  upd_mf=upd_mf,dwn_mf=dwn_mf,det_mf=det_mf,dqdti=dqdti,cnvqci=cnvqci,upd_mfi=upd_mfi, &
                  dwn_mfi=dwn_mfi,det_mfi=det_mfi,cnvw=cnvw,cnvc=cnvc,cnvw_phy_f3d=cnvw_phy_f3d, &
                  cnvc_phy_f3d=cnvc_phy_f3d,errmsg=cdata%errmsg,errflg=cdata%errflg)
        ierr=cdata%errflg

    end function GFS_DCNV_generic_post_run_cap

    function GFS_DCNV_generic_post_finalize_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call GFS_DCNV_generic_post_finalize()
        

    end function GFS_DCNV_generic_post_finalize_cap
end module GFS_DCNV_generic_post_cap

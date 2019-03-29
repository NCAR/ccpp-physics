
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
!! @brief Auto-generated cap module for the mynnedmf_wrapper scheme
!!
!
module mynnedmf_wrapper_cap

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_f_pointer, c_ptr, c_int32_t
    use            :: ccpp_types,                                      &
                      only: ccpp_t, CCPP_GENERIC_KIND
    use            :: ccpp_fields,                                     &
                      only: ccpp_field_get
    use            :: ccpp_errors,                                     &
                      only: ccpp_error, ccpp_debug
    use            :: mynnedmf_wrapper, &
                      only: mynnedmf_wrapper_init,mynnedmf_wrapper_finalize,mynnedmf_wrapper_run
    ! Other modules required, e.g. type definitions
    use machine, only: kind_phys

    implicit none

    private
    public :: mynnedmf_wrapper_init_cap,mynnedmf_wrapper_finalize_cap,mynnedmf_wrapper_run_cap

    contains


    function mynnedmf_wrapper_init_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call mynnedmf_wrapper_init()
        

    end function mynnedmf_wrapper_init_cap

    function mynnedmf_wrapper_finalize_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call mynnedmf_wrapper_finalize()
        

    end function mynnedmf_wrapper_finalize_cap

    function mynnedmf_wrapper_run_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind
        integer, pointer :: ix
        integer, pointer :: im
        integer, pointer :: levs
        logical, pointer :: flag_init
        logical, pointer :: flag_restart
        real(kind_phys), pointer :: delt
        real(kind_phys), pointer :: dx(:)
        real(kind_phys), pointer :: zorl(:)
        real(kind_phys), pointer :: phii(:,:)
        real(kind_phys), pointer :: U(:,:)
        real(kind_phys), pointer :: V(:,:)
        real(kind_phys), pointer :: omega(:,:)
        real(kind_phys), pointer :: T3D(:,:)
        real(kind_phys), pointer :: qgrs_water_vapor(:,:)
        real(kind_phys), pointer :: qgrs_liquid_cloud(:,:)
        real(kind_phys), pointer :: qgrs_ice_cloud(:,:)
        real(kind_phys), pointer :: qgrs_cloud_droplet_num_conc(:,:)
        real(kind_phys), pointer :: qgrs_cloud_ice_num_conc(:,:)
        real(kind_phys), pointer :: qgrs_ozone(:,:)
        real(kind_phys), pointer :: qgrs_water_aer_num_conc(:,:)
        real(kind_phys), pointer :: qgrs_ice_aer_num_conc(:,:)
        real(kind_phys), pointer :: prsl(:,:)
        real(kind_phys), pointer :: exner(:,:)
        real(kind_phys), pointer :: slmsk(:)
        real(kind_phys), pointer :: tsurf(:)
        real(kind_phys), pointer :: qsfc(:)
        real(kind_phys), pointer :: ps(:)
        real(kind_phys), pointer :: ust(:)
        real(kind_phys), pointer :: ch(:)
        real(kind_phys), pointer :: hflx(:)
        real(kind_phys), pointer :: qflx(:)
        real(kind_phys), pointer :: wspd(:)
        real(kind_phys), pointer :: rb(:)
        real(kind_phys), pointer :: dtsfc1(:)
        real(kind_phys), pointer :: dqsfc1(:)
        real(kind_phys), pointer :: dtsfci_diag(:)
        real(kind_phys), pointer :: dqsfci_diag(:)
        real(kind_phys), pointer :: dtsfc_diag(:)
        real(kind_phys), pointer :: dqsfc_diag(:)
        real(kind_phys), pointer :: recmol(:)
        real(kind_phys), pointer :: qke(:,:)
        real(kind_phys), pointer :: qke_adv(:,:)
        real(kind_phys), pointer :: tsq(:,:)
        real(kind_phys), pointer :: qsq(:,:)
        real(kind_phys), pointer :: cov(:,:)
        real(kind_phys), pointer :: el_pbl(:,:)
        real(kind_phys), pointer :: Sh3D(:,:)
        real(kind_phys), pointer :: exch_h(:,:)
        real(kind_phys), pointer :: exch_m(:,:)
        real(kind_phys), pointer :: PBLH(:)
        integer, pointer :: kpbl(:)
        real(kind_phys), pointer :: QC_BL(:,:)
        real(kind_phys), pointer :: CLDFRA_BL(:,:)
        real(kind_phys), pointer :: edmf_a(:,:)
        real(kind_phys), pointer :: edmf_w(:,:)
        real(kind_phys), pointer :: edmf_qt(:,:)
        real(kind_phys), pointer :: edmf_thl(:,:)
        real(kind_phys), pointer :: edmf_ent(:,:)
        real(kind_phys), pointer :: edmf_qc(:,:)
        integer, pointer :: nupdraft(:)
        real(kind_phys), pointer :: maxMF(:)
        integer, pointer :: ktop_shallow(:)
        real(kind_phys), pointer :: RTHRATEN(:,:)
        real(kind_phys), pointer :: dudt(:,:)
        real(kind_phys), pointer :: dvdt(:,:)
        real(kind_phys), pointer :: dtdt(:,:)
        real(kind_phys), pointer :: dqdt_water_vapor(:,:)
        real(kind_phys), pointer :: dqdt_liquid_cloud(:,:)
        real(kind_phys), pointer :: dqdt_ice_cloud(:,:)
        real(kind_phys), pointer :: dqdt_ozone(:,:)
        real(kind_phys), pointer :: dqdt_cloud_droplet_num_conc(:,:)
        real(kind_phys), pointer :: dqdt_ice_num_conc(:,:)
        real(kind_phys), pointer :: dqdt_water_aer_num_conc(:,:)
        real(kind_phys), pointer :: dqdt_ice_aer_num_conc(:,:)
        integer, pointer :: grav_settling
        integer, pointer :: bl_mynn_tkebudget
        logical, pointer :: bl_mynn_tkeadvect
        integer, pointer :: bl_mynn_cloudpdf
        integer, pointer :: bl_mynn_mixlength
        integer, pointer :: bl_mynn_edmf
        integer, pointer :: bl_mynn_edmf_mom
        integer, pointer :: bl_mynn_edmf_tke
        integer, pointer :: bl_mynn_edmf_part
        integer, pointer :: bl_mynn_cloudmix
        integer, pointer :: bl_mynn_mixqt
        integer, pointer :: icloud_bl
        logical, pointer :: do_mynnsfclay
        integer, pointer :: imp_physics
        integer, pointer :: imp_physics_gfdl
        integer, pointer :: imp_physics_thompson
        integer, pointer :: imp_physics_wsm6
        logical, pointer :: ltaerosol
        logical, pointer :: lprnt

        ierr = 0

        call c_f_pointer(ptr, cdata)


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
        

        call ccpp_field_get(cdata, 'flag_for_first_time_step', flag_init, ierr=ierr, kind=ckind, index=277)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_first_time_step from CCPP data structure')
            return
        end if
        if (kind(flag_init).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_first_time_step')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'flag_for_restart', flag_restart, ierr=ierr, kind=ckind, index=309)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_restart from CCPP data structure')
            return
        end if
        if (kind(flag_restart).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_restart')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'time_step_for_physics', delt, ierr=ierr, kind=ckind, index=793)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve time_step_for_physics from CCPP data structure')
            return
        end if
        if (kind(delt).ne.ckind) then
            call ccpp_error('Kind mismatch for variable time_step_for_physics')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'cell_size', dx, ierr=ierr, dims=cdims, kind=ckind, index=84)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cell_size from CCPP data structure')
            return
        end if
        if (kind(dx).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cell_size')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_roughness_length', zorl, ierr=ierr, dims=cdims, kind=ckind, index=718)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_roughness_length from CCPP data structure')
            return
        end if
        if (kind(zorl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_roughness_length')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'geopotential_at_interface', phii, ierr=ierr, dims=cdims, kind=ckind, index=349)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve geopotential_at_interface from CCPP data structure')
            return
        end if
        if (kind(phii).ne.ckind) then
            call ccpp_error('Kind mismatch for variable geopotential_at_interface')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'x_wind', U, ierr=ierr, dims=cdims, kind=ckind, index=871)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve x_wind from CCPP data structure')
            return
        end if
        if (kind(U).ne.ckind) then
            call ccpp_error('Kind mismatch for variable x_wind')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'y_wind', V, ierr=ierr, dims=cdims, kind=ckind, index=878)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve y_wind from CCPP data structure')
            return
        end if
        if (kind(V).ne.ckind) then
            call ccpp_error('Kind mismatch for variable y_wind')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'omega', omega, ierr=ierr, dims=cdims, kind=ckind, index=593)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve omega from CCPP data structure')
            return
        end if
        if (kind(omega).ne.ckind) then
            call ccpp_error('Kind mismatch for variable omega')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'air_temperature', T3D, ierr=ierr, dims=cdims, kind=ckind, index=50)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve air_temperature from CCPP data structure')
            return
        end if
        if (kind(T3D).ne.ckind) then
            call ccpp_error('Kind mismatch for variable air_temperature')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'water_vapor_specific_humidity', qgrs_water_vapor, ierr=ierr, dims=cdims, kind=ckind, index=852)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve water_vapor_specific_humidity from CCPP data structure')
            return
        end if
        if (kind(qgrs_water_vapor).ne.ckind) then
            call ccpp_error('Kind mismatch for variable water_vapor_specific_humidity')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cloud_condensed_water_mixing_ratio', qgrs_liquid_cloud, ierr=ierr, dims=cdims, kind=ckind, index=90)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cloud_condensed_water_mixing_ratio from CCPP data structure')
            return
        end if
        if (kind(qgrs_liquid_cloud).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cloud_condensed_water_mixing_ratio')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'ice_water_mixing_ratio', qgrs_ice_cloud, ierr=ierr, dims=cdims, kind=ckind, index=373)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve ice_water_mixing_ratio from CCPP data structure')
            return
        end if
        if (kind(qgrs_ice_cloud).ne.ckind) then
            call ccpp_error('Kind mismatch for variable ice_water_mixing_ratio')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cloud_droplet_number_concentration', qgrs_cloud_droplet_num_conc, ierr=ierr, dims=cdims, kind=ckind, index=97)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cloud_droplet_number_concentration from CCPP data structure')
            return
        end if
        if (kind(qgrs_cloud_droplet_num_conc).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cloud_droplet_number_concentration')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'ice_number_concentration', qgrs_cloud_ice_num_conc, ierr=ierr, dims=cdims, kind=ckind, index=370)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve ice_number_concentration from CCPP data structure')
            return
        end if
        if (kind(qgrs_cloud_ice_num_conc).ne.ckind) then
            call ccpp_error('Kind mismatch for variable ice_number_concentration')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'ozone_mixing_ratio', qgrs_ozone, ierr=ierr, dims=cdims, kind=ckind, index=600)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve ozone_mixing_ratio from CCPP data structure')
            return
        end if
        if (kind(qgrs_ozone).ne.ckind) then
            call ccpp_error('Kind mismatch for variable ozone_mixing_ratio')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'water_friendly_aerosol_number_concentration', qgrs_water_aer_num_conc, ierr=ierr, dims=cdims, kind=ckind, index=849)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve water_friendly_aerosol_number_concentration from CCPP data structure')
            return
        end if
        if (kind(qgrs_water_aer_num_conc).ne.ckind) then
            call ccpp_error('Kind mismatch for variable water_friendly_aerosol_number_concentration')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'ice_friendly_aerosol_number_concentration', qgrs_ice_aer_num_conc, ierr=ierr, dims=cdims, kind=ckind, index=368)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve ice_friendly_aerosol_number_concentration from CCPP data structure')
            return
        end if
        if (kind(qgrs_ice_aer_num_conc).ne.ckind) then
            call ccpp_error('Kind mismatch for variable ice_friendly_aerosol_number_concentration')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'air_pressure', prsl, ierr=ierr, dims=cdims, kind=ckind, index=44)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve air_pressure from CCPP data structure')
            return
        end if
        if (kind(prsl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable air_pressure')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'dimensionless_exner_function_at_model_layers', exner, ierr=ierr, dims=cdims, kind=ckind, index=222)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve dimensionless_exner_function_at_model_layers from CCPP data structure')
            return
        end if
        if (kind(exner).ne.ckind) then
            call ccpp_error('Kind mismatch for variable dimensionless_exner_function_at_model_layers')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'sea_land_ice_mask_real', slmsk, ierr=ierr, dims=cdims, kind=ckind, index=632)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve sea_land_ice_mask_real from CCPP data structure')
            return
        end if
        if (kind(slmsk).ne.ckind) then
            call ccpp_error('Kind mismatch for variable sea_land_ice_mask_real')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_skin_temperature', tsurf, ierr=ierr, dims=cdims, kind=ckind, index=721)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_skin_temperature from CCPP data structure')
            return
        end if
        if (kind(tsurf).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_skin_temperature')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_specific_humidity', qsfc, ierr=ierr, dims=cdims, kind=ckind, index=730)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_specific_humidity from CCPP data structure')
            return
        end if
        if (kind(qsfc).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_specific_humidity')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_air_pressure', ps, ierr=ierr, dims=cdims, kind=ckind, index=677)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_air_pressure from CCPP data structure')
            return
        end if
        if (kind(ps).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_air_pressure')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_friction_velocity', ust, ierr=ierr, dims=cdims, kind=ckind, index=710)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_friction_velocity from CCPP data structure')
            return
        end if
        if (kind(ust).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_friction_velocity')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_drag_wind_speed_for_momentum_in_air', ch, ierr=ierr, dims=cdims, kind=ckind, index=705)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_drag_wind_speed_for_momentum_in_air from CCPP data structure')
            return
        end if
        if (kind(ch).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_drag_wind_speed_for_momentum_in_air')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'kinematic_surface_upward_sensible_heat_flux', hflx, ierr=ierr, dims=cdims, kind=ckind, index=455)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve kinematic_surface_upward_sensible_heat_flux from CCPP data structure')
            return
        end if
        if (kind(hflx).ne.ckind) then
            call ccpp_error('Kind mismatch for variable kinematic_surface_upward_sensible_heat_flux')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'kinematic_surface_upward_latent_heat_flux', qflx, ierr=ierr, dims=cdims, kind=ckind, index=454)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve kinematic_surface_upward_latent_heat_flux from CCPP data structure')
            return
        end if
        if (kind(qflx).ne.ckind) then
            call ccpp_error('Kind mismatch for variable kinematic_surface_upward_latent_heat_flux')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'wind_speed_at_lowest_model_layer', wspd, ierr=ierr, dims=cdims, kind=ckind, index=870)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve wind_speed_at_lowest_model_layer from CCPP data structure')
            return
        end if
        if (kind(wspd).ne.ckind) then
            call ccpp_error('Kind mismatch for variable wind_speed_at_lowest_model_layer')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'bulk_richardson_number_at_lowest_model_level', rb, ierr=ierr, dims=cdims, kind=ckind, index=78)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve bulk_richardson_number_at_lowest_model_level from CCPP data structure')
            return
        end if
        if (kind(rb).ne.ckind) then
            call ccpp_error('Kind mismatch for variable bulk_richardson_number_at_lowest_model_level')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'instantaneous_surface_upward_sensible_heat_flux', dtsfc1, ierr=ierr, dims=cdims, kind=ckind, index=434)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve instantaneous_surface_upward_sensible_heat_flux from CCPP data structure')
            return
        end if
        if (kind(dtsfc1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable instantaneous_surface_upward_sensible_heat_flux')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'instantaneous_surface_upward_latent_heat_flux', dqsfc1, ierr=ierr, dims=cdims, kind=ckind, index=431)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve instantaneous_surface_upward_latent_heat_flux from CCPP data structure')
            return
        end if
        if (kind(dqsfc1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable instantaneous_surface_upward_latent_heat_flux')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'instantaneous_surface_upward_sensible_heat_flux_for_diag', dtsfci_diag, ierr=ierr, dims=cdims, kind=ckind, index=436)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve instantaneous_surface_upward_sensible_heat_flux_for_diag from CCPP data structure')
            return
        end if
        if (kind(dtsfci_diag).ne.ckind) then
            call ccpp_error('Kind mismatch for variable instantaneous_surface_upward_sensible_heat_flux_for_diag')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'instantaneous_surface_upward_latent_heat_flux_for_diag', dqsfci_diag, ierr=ierr, dims=cdims, kind=ckind, index=433)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve instantaneous_surface_upward_latent_heat_flux_for_diag from CCPP data structure')
            return
        end if
        if (kind(dqsfci_diag).ne.ckind) then
            call ccpp_error('Kind mismatch for variable instantaneous_surface_upward_latent_heat_flux_for_diag')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cumulative_surface_upward_sensible_heat_flux_for_diag_multiplied_by_timestep', dtsfc_diag, ierr=ierr, dims=cdims, kind=ckind, index=199)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cumulative_surface_upward_sensible_heat_flux_for_diag_multiplied_by_timestep from CCPP data structure')
            return
        end if
        if (kind(dtsfc_diag).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cumulative_surface_upward_sensible_heat_flux_for_diag_multiplied_by_timestep')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cumulative_surface_upward_latent_heat_flux_for_diag_multiplied_by_timestep', dqsfc_diag, ierr=ierr, dims=cdims, kind=ckind, index=196)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cumulative_surface_upward_latent_heat_flux_for_diag_multiplied_by_timestep from CCPP data structure')
            return
        end if
        if (kind(dqsfc_diag).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cumulative_surface_upward_latent_heat_flux_for_diag_multiplied_by_timestep')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'reciprocal_of_obukhov_length', recmol, ierr=ierr, dims=cdims, kind=ckind, index=627)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve reciprocal_of_obukhov_length from CCPP data structure')
            return
        end if
        if (kind(recmol).ne.ckind) then
            call ccpp_error('Kind mismatch for variable reciprocal_of_obukhov_length')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'tke_at_mass_points', qke, ierr=ierr, dims=cdims, kind=ckind, index=796)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve tke_at_mass_points from CCPP data structure')
            return
        end if
        if (kind(qke).ne.ckind) then
            call ccpp_error('Kind mismatch for variable tke_at_mass_points')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'turbulent_kinetic_energy', qke_adv, ierr=ierr, dims=cdims, kind=ckind, index=807)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve turbulent_kinetic_energy from CCPP data structure')
            return
        end if
        if (kind(qke_adv).ne.ckind) then
            call ccpp_error('Kind mismatch for variable turbulent_kinetic_energy')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 't_prime_squared', tsq, ierr=ierr, dims=cdims, kind=ckind, index=749)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve t_prime_squared from CCPP data structure')
            return
        end if
        if (kind(tsq).ne.ckind) then
            call ccpp_error('Kind mismatch for variable t_prime_squared')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'q_prime_squared', qsq, ierr=ierr, dims=cdims, kind=ckind, index=612)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve q_prime_squared from CCPP data structure')
            return
        end if
        if (kind(qsq).ne.ckind) then
            call ccpp_error('Kind mismatch for variable q_prime_squared')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 't_prime_q_prime', cov, ierr=ierr, dims=cdims, kind=ckind, index=748)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve t_prime_q_prime from CCPP data structure')
            return
        end if
        if (kind(cov).ne.ckind) then
            call ccpp_error('Kind mismatch for variable t_prime_q_prime')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'mixing_length', el_pbl, ierr=ierr, dims=cdims, kind=ckind, index=549)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve mixing_length from CCPP data structure')
            return
        end if
        if (kind(el_pbl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable mixing_length')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'stability_function_for_heat', Sh3D, ierr=ierr, dims=cdims, kind=ckind, index=668)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve stability_function_for_heat from CCPP data structure')
            return
        end if
        if (kind(Sh3D).ne.ckind) then
            call ccpp_error('Kind mismatch for variable stability_function_for_heat')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'atmosphere_heat_diffusivity_for_mynnpbl', exch_h, ierr=ierr, dims=cdims, kind=ckind, index=71)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve atmosphere_heat_diffusivity_for_mynnpbl from CCPP data structure')
            return
        end if
        if (kind(exch_h).ne.ckind) then
            call ccpp_error('Kind mismatch for variable atmosphere_heat_diffusivity_for_mynnpbl')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'atmosphere_momentum_diffusivity_for_mynnpbl', exch_m, ierr=ierr, dims=cdims, kind=ckind, index=74)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve atmosphere_momentum_diffusivity_for_mynnpbl from CCPP data structure')
            return
        end if
        if (kind(exch_m).ne.ckind) then
            call ccpp_error('Kind mismatch for variable atmosphere_momentum_diffusivity_for_mynnpbl')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'atmosphere_boundary_layer_thickness', PBLH, ierr=ierr, dims=cdims, kind=ckind, index=66)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve atmosphere_boundary_layer_thickness from CCPP data structure')
            return
        end if
        if (kind(PBLH).ne.ckind) then
            call ccpp_error('Kind mismatch for variable atmosphere_boundary_layer_thickness')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'vertical_index_at_top_of_atmosphere_boundary_layer', kpbl, ierr=ierr, dims=cdims, kind=ckind, index=822)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve vertical_index_at_top_of_atmosphere_boundary_layer from CCPP data structure')
            return
        end if
        if (kind(kpbl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable vertical_index_at_top_of_atmosphere_boundary_layer')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'subgrid_cloud_mixing_ratio_pbl', QC_BL, ierr=ierr, dims=cdims, kind=ckind, index=674)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve subgrid_cloud_mixing_ratio_pbl from CCPP data structure')
            return
        end if
        if (kind(QC_BL).ne.ckind) then
            call ccpp_error('Kind mismatch for variable subgrid_cloud_mixing_ratio_pbl')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'subgrid_cloud_fraction_pbl', CLDFRA_BL, ierr=ierr, dims=cdims, kind=ckind, index=673)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve subgrid_cloud_fraction_pbl from CCPP data structure')
            return
        end if
        if (kind(CLDFRA_BL).ne.ckind) then
            call ccpp_error('Kind mismatch for variable subgrid_cloud_fraction_pbl')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'emdf_updraft_area', edmf_a, ierr=ierr, dims=cdims, kind=ckind, index=247)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve emdf_updraft_area from CCPP data structure')
            return
        end if
        if (kind(edmf_a).ne.ckind) then
            call ccpp_error('Kind mismatch for variable emdf_updraft_area')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'emdf_updraft_vertical_velocity', edmf_w, ierr=ierr, dims=cdims, kind=ckind, index=252)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve emdf_updraft_vertical_velocity from CCPP data structure')
            return
        end if
        if (kind(edmf_w).ne.ckind) then
            call ccpp_error('Kind mismatch for variable emdf_updraft_vertical_velocity')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'emdf_updraft_total_water', edmf_qt, ierr=ierr, dims=cdims, kind=ckind, index=251)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve emdf_updraft_total_water from CCPP data structure')
            return
        end if
        if (kind(edmf_qt).ne.ckind) then
            call ccpp_error('Kind mismatch for variable emdf_updraft_total_water')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'emdf_updraft_theta_l', edmf_thl, ierr=ierr, dims=cdims, kind=ckind, index=250)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve emdf_updraft_theta_l from CCPP data structure')
            return
        end if
        if (kind(edmf_thl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable emdf_updraft_theta_l')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'emdf_updraft_entrainment_rate', edmf_ent, ierr=ierr, dims=cdims, kind=ckind, index=249)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve emdf_updraft_entrainment_rate from CCPP data structure')
            return
        end if
        if (kind(edmf_ent).ne.ckind) then
            call ccpp_error('Kind mismatch for variable emdf_updraft_entrainment_rate')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'emdf_updraft_cloud_water', edmf_qc, ierr=ierr, dims=cdims, kind=ckind, index=248)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve emdf_updraft_cloud_water from CCPP data structure')
            return
        end if
        if (kind(edmf_qc).ne.ckind) then
            call ccpp_error('Kind mismatch for variable emdf_updraft_cloud_water')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'number_of_plumes', nupdraft, ierr=ierr, dims=cdims, kind=ckind, index=580)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve number_of_plumes from CCPP data structure')
            return
        end if
        if (kind(nupdraft).ne.ckind) then
            call ccpp_error('Kind mismatch for variable number_of_plumes')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'maximum_mass_flux', maxMF, ierr=ierr, dims=cdims, kind=ckind, index=505)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve maximum_mass_flux from CCPP data structure')
            return
        end if
        if (kind(maxMF).ne.ckind) then
            call ccpp_error('Kind mismatch for variable maximum_mass_flux')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'k_level_of_highest_reaching_plume', ktop_shallow, ierr=ierr, dims=cdims, kind=ckind, index=452)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve k_level_of_highest_reaching_plume from CCPP data structure')
            return
        end if
        if (kind(ktop_shallow).ne.ckind) then
            call ccpp_error('Kind mismatch for variable k_level_of_highest_reaching_plume')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'tendency_of_air_temperature_due_to_longwave_heating_on_radiation_time_step', RTHRATEN, ierr=ierr, dims=cdims, kind=ckind, index=756)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve tendency_of_air_temperature_due_to_longwave_heating_on_radiation_time_step from CCPP data structure')
            return
        end if
        if (kind(RTHRATEN).ne.ckind) then
            call ccpp_error('Kind mismatch for variable tendency_of_air_temperature_due_to_longwave_heating_on_radiation_time_step')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'tendency_of_x_wind_due_to_model_physics', dudt, ierr=ierr, dims=cdims, kind=ckind, index=782)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve tendency_of_x_wind_due_to_model_physics from CCPP data structure')
            return
        end if
        if (kind(dudt).ne.ckind) then
            call ccpp_error('Kind mismatch for variable tendency_of_x_wind_due_to_model_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'tendency_of_y_wind_due_to_model_physics', dvdt, ierr=ierr, dims=cdims, kind=ckind, index=785)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve tendency_of_y_wind_due_to_model_physics from CCPP data structure')
            return
        end if
        if (kind(dvdt).ne.ckind) then
            call ccpp_error('Kind mismatch for variable tendency_of_y_wind_due_to_model_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'tendency_of_air_temperature_due_to_model_physics', dtdt, ierr=ierr, dims=cdims, kind=ckind, index=758)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve tendency_of_air_temperature_due_to_model_physics from CCPP data structure')
            return
        end if
        if (kind(dtdt).ne.ckind) then
            call ccpp_error('Kind mismatch for variable tendency_of_air_temperature_due_to_model_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'tendency_of_water_vapor_specific_humidity_due_to_model_physics', dqdt_water_vapor, ierr=ierr, dims=cdims, kind=ckind, index=780)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve tendency_of_water_vapor_specific_humidity_due_to_model_physics from CCPP data structure')
            return
        end if
        if (kind(dqdt_water_vapor).ne.ckind) then
            call ccpp_error('Kind mismatch for variable tendency_of_water_vapor_specific_humidity_due_to_model_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'tendency_of_liquid_cloud_water_mixing_ratio_due_to_model_physics', dqdt_liquid_cloud, ierr=ierr, dims=cdims, kind=ckind, index=771)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve tendency_of_liquid_cloud_water_mixing_ratio_due_to_model_physics from CCPP data structure')
            return
        end if
        if (kind(dqdt_liquid_cloud).ne.ckind) then
            call ccpp_error('Kind mismatch for variable tendency_of_liquid_cloud_water_mixing_ratio_due_to_model_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'tendency_of_ice_cloud_water_mixing_ratio_due_to_model_physics', dqdt_ice_cloud, ierr=ierr, dims=cdims, kind=ckind, index=767)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve tendency_of_ice_cloud_water_mixing_ratio_due_to_model_physics from CCPP data structure')
            return
        end if
        if (kind(dqdt_ice_cloud).ne.ckind) then
            call ccpp_error('Kind mismatch for variable tendency_of_ice_cloud_water_mixing_ratio_due_to_model_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'tendency_of_ozone_mixing_ratio_due_to_model_physics', dqdt_ozone, ierr=ierr, dims=cdims, kind=ckind, index=774)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve tendency_of_ozone_mixing_ratio_due_to_model_physics from CCPP data structure')
            return
        end if
        if (kind(dqdt_ozone).ne.ckind) then
            call ccpp_error('Kind mismatch for variable tendency_of_ozone_mixing_ratio_due_to_model_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'tendency_of_cloud_droplet_number_concentration_due_to_model_physics', dqdt_cloud_droplet_num_conc, ierr=ierr, dims=cdims, kind=ckind, index=765)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve tendency_of_cloud_droplet_number_concentration_due_to_model_physics from CCPP data structure')
            return
        end if
        if (kind(dqdt_cloud_droplet_num_conc).ne.ckind) then
            call ccpp_error('Kind mismatch for variable tendency_of_cloud_droplet_number_concentration_due_to_model_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'tendency_of_ice_number_concentration_due_to_model_physics', dqdt_ice_num_conc, ierr=ierr, dims=cdims, kind=ckind, index=770)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve tendency_of_ice_number_concentration_due_to_model_physics from CCPP data structure')
            return
        end if
        if (kind(dqdt_ice_num_conc).ne.ckind) then
            call ccpp_error('Kind mismatch for variable tendency_of_ice_number_concentration_due_to_model_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'tendency_of_water_friendly_aerosol_number_concentration_due_to_model_physics', dqdt_water_aer_num_conc, ierr=ierr, dims=cdims, kind=ckind, index=778)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve tendency_of_water_friendly_aerosol_number_concentration_due_to_model_physics from CCPP data structure')
            return
        end if
        if (kind(dqdt_water_aer_num_conc).ne.ckind) then
            call ccpp_error('Kind mismatch for variable tendency_of_water_friendly_aerosol_number_concentration_due_to_model_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'tendency_of_ice_friendly_aerosol_number_concentration_due_to_model_physics', dqdt_ice_aer_num_conc, ierr=ierr, dims=cdims, kind=ckind, index=768)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve tendency_of_ice_friendly_aerosol_number_concentration_due_to_model_physics from CCPP data structure')
            return
        end if
        if (kind(dqdt_ice_aer_num_conc).ne.ckind) then
            call ccpp_error('Kind mismatch for variable tendency_of_ice_friendly_aerosol_number_concentration_due_to_model_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'grav_settling', grav_settling, ierr=ierr, kind=ckind, index=354)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve grav_settling from CCPP data structure')
            return
        end if
        if (kind(grav_settling).ne.ckind) then
            call ccpp_error('Kind mismatch for variable grav_settling')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'tke_budget', bl_mynn_tkebudget, ierr=ierr, kind=ckind, index=797)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve tke_budget from CCPP data structure')
            return
        end if
        if (kind(bl_mynn_tkebudget).ne.ckind) then
            call ccpp_error('Kind mismatch for variable tke_budget')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'tke_advect', bl_mynn_tkeadvect, ierr=ierr, kind=ckind, index=795)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve tke_advect from CCPP data structure')
            return
        end if
        if (kind(bl_mynn_tkeadvect).ne.ckind) then
            call ccpp_error('Kind mismatch for variable tke_advect')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'cloudpdf', bl_mynn_cloudpdf, ierr=ierr, kind=ckind, index=112)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cloudpdf from CCPP data structure')
            return
        end if
        if (kind(bl_mynn_cloudpdf).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cloudpdf')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'mixing_length_flag', bl_mynn_mixlength, ierr=ierr, kind=ckind, index=550)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve mixing_length_flag from CCPP data structure')
            return
        end if
        if (kind(bl_mynn_mixlength).ne.ckind) then
            call ccpp_error('Kind mismatch for variable mixing_length_flag')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'edmf_flag', bl_mynn_edmf, ierr=ierr, kind=ckind, index=238)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve edmf_flag from CCPP data structure')
            return
        end if
        if (kind(bl_mynn_edmf).ne.ckind) then
            call ccpp_error('Kind mismatch for variable edmf_flag')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'edmf_momentum_transport_flag', bl_mynn_edmf_mom, ierr=ierr, kind=ckind, index=239)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve edmf_momentum_transport_flag from CCPP data structure')
            return
        end if
        if (kind(bl_mynn_edmf_mom).ne.ckind) then
            call ccpp_error('Kind mismatch for variable edmf_momentum_transport_flag')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'edmf_tke_transport_flag', bl_mynn_edmf_tke, ierr=ierr, kind=ckind, index=241)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve edmf_tke_transport_flag from CCPP data structure')
            return
        end if
        if (kind(bl_mynn_edmf_tke).ne.ckind) then
            call ccpp_error('Kind mismatch for variable edmf_tke_transport_flag')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'edmf_partition_flag', bl_mynn_edmf_part, ierr=ierr, kind=ckind, index=240)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve edmf_partition_flag from CCPP data structure')
            return
        end if
        if (kind(bl_mynn_edmf_part).ne.ckind) then
            call ccpp_error('Kind mismatch for variable edmf_partition_flag')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'cloud_specie_mix_flag', bl_mynn_cloudmix, ierr=ierr, kind=ckind, index=109)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cloud_specie_mix_flag from CCPP data structure')
            return
        end if
        if (kind(bl_mynn_cloudmix).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cloud_specie_mix_flag')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'mix_total_water_flag', bl_mynn_mixqt, ierr=ierr, kind=ckind, index=548)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve mix_total_water_flag from CCPP data structure')
            return
        end if
        if (kind(bl_mynn_mixqt).ne.ckind) then
            call ccpp_error('Kind mismatch for variable mix_total_water_flag')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'couple_sgs_clouds_to_radiation_flag', icloud_bl, ierr=ierr, kind=ckind, index=139)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve couple_sgs_clouds_to_radiation_flag from CCPP data structure')
            return
        end if
        if (kind(icloud_bl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable couple_sgs_clouds_to_radiation_flag')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'do_mynnsfclay', do_mynnsfclay, ierr=ierr, kind=ckind, index=228)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve do_mynnsfclay from CCPP data structure')
            return
        end if
        if (kind(do_mynnsfclay).ne.ckind) then
            call ccpp_error('Kind mismatch for variable do_mynnsfclay')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'flag_for_microphysics_scheme', imp_physics, ierr=ierr, kind=ckind, index=294)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_microphysics_scheme from CCPP data structure')
            return
        end if
        if (kind(imp_physics).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_microphysics_scheme')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'flag_for_gfdl_microphysics_scheme', imp_physics_gfdl, ierr=ierr, kind=ckind, index=280)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_gfdl_microphysics_scheme from CCPP data structure')
            return
        end if
        if (kind(imp_physics_gfdl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_gfdl_microphysics_scheme')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'flag_for_thompson_microphysics_scheme', imp_physics_thompson, ierr=ierr, kind=ckind, index=322)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_thompson_microphysics_scheme from CCPP data structure')
            return
        end if
        if (kind(imp_physics_thompson).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_thompson_microphysics_scheme')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'flag_for_wsm6_microphysics_scheme', imp_physics_wsm6, ierr=ierr, kind=ckind, index=326)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_wsm6_microphysics_scheme from CCPP data structure')
            return
        end if
        if (kind(imp_physics_wsm6).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_wsm6_microphysics_scheme')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'flag_for_aerosol_physics', ltaerosol, ierr=ierr, kind=ckind, index=271)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_aerosol_physics from CCPP data structure')
            return
        end if
        if (kind(ltaerosol).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_aerosol_physics')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'flag_print', lprnt, ierr=ierr, kind=ckind, index=332)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_print from CCPP data structure')
            return
        end if
        if (kind(lprnt).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_print')
            ierr = 1
            return
        end if
#endif
        

        call mynnedmf_wrapper_run(ix=ix,im=im,levs=levs,flag_init=flag_init,flag_restart=flag_restart,delt=delt, &
                  dx=dx,zorl=zorl,phii=phii,U=U,V=V,omega=omega,T3D=T3D,qgrs_water_vapor=qgrs_water_vapor, &
                  qgrs_liquid_cloud=qgrs_liquid_cloud,qgrs_ice_cloud=qgrs_ice_cloud,qgrs_cloud_droplet_num_conc=qgrs_cloud_droplet_num_conc, &
                  qgrs_cloud_ice_num_conc=qgrs_cloud_ice_num_conc,qgrs_ozone=qgrs_ozone,qgrs_water_aer_num_conc=qgrs_water_aer_num_conc, &
                  qgrs_ice_aer_num_conc=qgrs_ice_aer_num_conc,prsl=prsl,exner=exner,slmsk=slmsk, &
                  tsurf=tsurf,qsfc=qsfc,ps=ps,ust=ust,ch=ch,hflx=hflx,qflx=qflx,wspd=wspd, &
                  rb=rb,dtsfc1=dtsfc1,dqsfc1=dqsfc1,dtsfci_diag=dtsfci_diag,dqsfci_diag=dqsfci_diag, &
                  dtsfc_diag=dtsfc_diag,dqsfc_diag=dqsfc_diag,recmol=recmol,qke=qke,qke_adv=qke_adv, &
                  tsq=tsq,qsq=qsq,cov=cov,el_pbl=el_pbl,Sh3D=Sh3D,exch_h=exch_h,exch_m=exch_m, &
                  PBLH=PBLH,kpbl=kpbl,QC_BL=QC_BL,CLDFRA_BL=CLDFRA_BL,edmf_a=edmf_a,edmf_w=edmf_w, &
                  edmf_qt=edmf_qt,edmf_thl=edmf_thl,edmf_ent=edmf_ent,edmf_qc=edmf_qc,nupdraft=nupdraft, &
                  maxMF=maxMF,ktop_shallow=ktop_shallow,RTHRATEN=RTHRATEN,dudt=dudt,dvdt=dvdt, &
                  dtdt=dtdt,dqdt_water_vapor=dqdt_water_vapor,dqdt_liquid_cloud=dqdt_liquid_cloud, &
                  dqdt_ice_cloud=dqdt_ice_cloud,dqdt_ozone=dqdt_ozone,dqdt_cloud_droplet_num_conc=dqdt_cloud_droplet_num_conc, &
                  dqdt_ice_num_conc=dqdt_ice_num_conc,dqdt_water_aer_num_conc=dqdt_water_aer_num_conc, &
                  dqdt_ice_aer_num_conc=dqdt_ice_aer_num_conc,grav_settling=grav_settling, &
                  bl_mynn_tkebudget=bl_mynn_tkebudget,bl_mynn_tkeadvect=bl_mynn_tkeadvect, &
                  bl_mynn_cloudpdf=bl_mynn_cloudpdf,bl_mynn_mixlength=bl_mynn_mixlength,bl_mynn_edmf=bl_mynn_edmf, &
                  bl_mynn_edmf_mom=bl_mynn_edmf_mom,bl_mynn_edmf_tke=bl_mynn_edmf_tke,bl_mynn_edmf_part=bl_mynn_edmf_part, &
                  bl_mynn_cloudmix=bl_mynn_cloudmix,bl_mynn_mixqt=bl_mynn_mixqt,icloud_bl=icloud_bl, &
                  do_mynnsfclay=do_mynnsfclay,imp_physics=imp_physics,imp_physics_gfdl=imp_physics_gfdl, &
                  imp_physics_thompson=imp_physics_thompson,imp_physics_wsm6=imp_physics_wsm6, &
                  ltaerosol=ltaerosol,lprnt=lprnt,errmsg=cdata%errmsg,errflg=cdata%errflg)
        ierr=cdata%errflg

    end function mynnedmf_wrapper_run_cap
end module mynnedmf_wrapper_cap

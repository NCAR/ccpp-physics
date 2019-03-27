
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
!! @brief Auto-generated cap module for the GFS_PBL_generic_post scheme
!!
!
module GFS_PBL_generic_post_cap

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_f_pointer, c_ptr, c_int32_t
    use            :: ccpp_types,                                      &
                      only: ccpp_t, CCPP_GENERIC_KIND
    use            :: ccpp_fields,                                     &
                      only: ccpp_field_get
    use            :: ccpp_errors,                                     &
                      only: ccpp_error, ccpp_debug
    use            :: GFS_PBL_generic_post, &
                      only: GFS_PBL_generic_post_init,GFS_PBL_generic_post_finalize,GFS_PBL_generic_post_run
    ! Other modules required, e.g. type definitions
    use machine, only: kind_phys

    implicit none

    private
    public :: GFS_PBL_generic_post_init_cap,GFS_PBL_generic_post_finalize_cap,GFS_PBL_generic_post_run_cap

    contains


    function GFS_PBL_generic_post_init_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call GFS_PBL_generic_post_init()
        

    end function GFS_PBL_generic_post_init_cap

    function GFS_PBL_generic_post_finalize_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call GFS_PBL_generic_post_finalize()
        

    end function GFS_PBL_generic_post_finalize_cap

    function GFS_PBL_generic_post_run_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind
        integer, pointer :: im
        integer, pointer :: levs
        integer, pointer :: nvdiff
        integer, pointer :: ntrac
        integer, pointer :: ntqv
        integer, pointer :: ntcw
        integer, pointer :: ntiw
        integer, pointer :: ntrw
        integer, pointer :: ntsw
        integer, pointer :: ntlnc
        integer, pointer :: ntinc
        integer, pointer :: ntwa
        integer, pointer :: ntia
        integer, pointer :: ntgl
        integer, pointer :: ntoz
        integer, pointer :: ntke
        integer, pointer :: ntkev
        integer, pointer :: imp_physics
        integer, pointer :: imp_physics_gfdl
        integer, pointer :: imp_physics_thompson
        integer, pointer :: imp_physics_wsm6
        logical, pointer :: ltaerosol
        logical, pointer :: cplflx
        logical, pointer :: lssav
        logical, pointer :: ldiag3d
        logical, pointer :: lsidea
        logical, pointer :: hybedmf
        logical, pointer :: do_shoc
        logical, pointer :: satmedmf
        real(kind_phys), pointer :: dvdftra(:,:,:)
        real(kind_phys), pointer :: dusfc1(:)
        real(kind_phys), pointer :: dvsfc1(:)
        real(kind_phys), pointer :: dtsfc1(:)
        real(kind_phys), pointer :: dqsfc1(:)
        real(kind_phys), pointer :: dtf
        real(kind_phys), pointer :: dudt(:,:)
        real(kind_phys), pointer :: dvdt(:,:)
        real(kind_phys), pointer :: dtdt(:,:)
        real(kind_phys), pointer :: htrsw(:,:)
        real(kind_phys), pointer :: htrlw(:,:)
        real(kind_phys), pointer :: xmu(:)
        real(kind_phys), pointer :: dqdt(:,:,:)
        real(kind_phys), pointer :: dusfc_cpl(:)
        real(kind_phys), pointer :: dvsfc_cpl(:)
        real(kind_phys), pointer :: dtsfc_cpl(:)
        real(kind_phys), pointer :: dqsfc_cpl(:)
        real(kind_phys), pointer :: dusfci_cpl(:)
        real(kind_phys), pointer :: dvsfci_cpl(:)
        real(kind_phys), pointer :: dtsfci_cpl(:)
        real(kind_phys), pointer :: dqsfci_cpl(:)
        real(kind_phys), pointer :: dusfc_diag(:)
        real(kind_phys), pointer :: dvsfc_diag(:)
        real(kind_phys), pointer :: dtsfc_diag(:)
        real(kind_phys), pointer :: dqsfc_diag(:)
        real(kind_phys), pointer :: dusfci_diag(:)
        real(kind_phys), pointer :: dvsfci_diag(:)
        real(kind_phys), pointer :: dtsfci_diag(:)
        real(kind_phys), pointer :: dqsfci_diag(:)
        real(kind_phys), pointer :: dt3dt(:,:)
        real(kind_phys), pointer :: du3dt_PBL(:,:)
        real(kind_phys), pointer :: du3dt_OGWD(:,:)
        real(kind_phys), pointer :: dv3dt_PBL(:,:)
        real(kind_phys), pointer :: dv3dt_OGWD(:,:)
        real(kind_phys), pointer :: dq3dt(:,:)
        real(kind_phys), pointer :: dq3dt_ozone(:,:)

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
        

        call ccpp_field_get(cdata, 'number_of_vertical_diffusion_tracers', nvdiff, ierr=ierr, kind=ckind, index=590)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve number_of_vertical_diffusion_tracers from CCPP data structure')
            return
        end if
        if (kind(nvdiff).ne.ckind) then
            call ccpp_error('Kind mismatch for variable number_of_vertical_diffusion_tracers')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'number_of_tracers', ntrac, ierr=ierr, kind=ckind, index=584)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve number_of_tracers from CCPP data structure')
            return
        end if
        if (kind(ntrac).ne.ckind) then
            call ccpp_error('Kind mismatch for variable number_of_tracers')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'index_for_water_vapor', ntqv, ierr=ierr, kind=ckind, index=395)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve index_for_water_vapor from CCPP data structure')
            return
        end if
        if (kind(ntqv).ne.ckind) then
            call ccpp_error('Kind mismatch for variable index_for_water_vapor')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'index_for_liquid_cloud_condensate', ntcw, ierr=ierr, kind=ckind, index=384)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve index_for_liquid_cloud_condensate from CCPP data structure')
            return
        end if
        if (kind(ntcw).ne.ckind) then
            call ccpp_error('Kind mismatch for variable index_for_liquid_cloud_condensate')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'index_for_ice_cloud_condensate', ntiw, ierr=ierr, kind=ckind, index=381)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve index_for_ice_cloud_condensate from CCPP data structure')
            return
        end if
        if (kind(ntiw).ne.ckind) then
            call ccpp_error('Kind mismatch for variable index_for_ice_cloud_condensate')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'index_for_rain_water', ntrw, ierr=ierr, kind=ckind, index=388)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve index_for_rain_water from CCPP data structure')
            return
        end if
        if (kind(ntrw).ne.ckind) then
            call ccpp_error('Kind mismatch for variable index_for_rain_water')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'index_for_snow_water', ntsw, ierr=ierr, kind=ckind, index=390)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve index_for_snow_water from CCPP data structure')
            return
        end if
        if (kind(ntsw).ne.ckind) then
            call ccpp_error('Kind mismatch for variable index_for_snow_water')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'index_for_liquid_cloud_number_concentration', ntlnc, ierr=ierr, kind=ckind, index=385)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve index_for_liquid_cloud_number_concentration from CCPP data structure')
            return
        end if
        if (kind(ntlnc).ne.ckind) then
            call ccpp_error('Kind mismatch for variable index_for_liquid_cloud_number_concentration')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'index_for_ice_cloud_number_concentration', ntinc, ierr=ierr, kind=ckind, index=382)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve index_for_ice_cloud_number_concentration from CCPP data structure')
            return
        end if
        if (kind(ntinc).ne.ckind) then
            call ccpp_error('Kind mismatch for variable index_for_ice_cloud_number_concentration')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'index_for_water_friendly_aerosols', ntwa, ierr=ierr, kind=ckind, index=394)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve index_for_water_friendly_aerosols from CCPP data structure')
            return
        end if
        if (kind(ntwa).ne.ckind) then
            call ccpp_error('Kind mismatch for variable index_for_water_friendly_aerosols')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'index_for_ice_friendly_aerosols', ntia, ierr=ierr, kind=ckind, index=383)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve index_for_ice_friendly_aerosols from CCPP data structure')
            return
        end if
        if (kind(ntia).ne.ckind) then
            call ccpp_error('Kind mismatch for variable index_for_ice_friendly_aerosols')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'index_for_graupel', ntgl, ierr=ierr, kind=ckind, index=379)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve index_for_graupel from CCPP data structure')
            return
        end if
        if (kind(ntgl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable index_for_graupel')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'index_for_ozone', ntoz, ierr=ierr, kind=ckind, index=386)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve index_for_ozone from CCPP data structure')
            return
        end if
        if (kind(ntoz).ne.ckind) then
            call ccpp_error('Kind mismatch for variable index_for_ozone')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'index_for_turbulent_kinetic_energy', ntke, ierr=ierr, kind=ckind, index=391)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve index_for_turbulent_kinetic_energy from CCPP data structure')
            return
        end if
        if (kind(ntke).ne.ckind) then
            call ccpp_error('Kind mismatch for variable index_for_turbulent_kinetic_energy')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'index_for_turbulent_kinetic_energy_vertical_diffusion_tracer', ntkev, ierr=ierr, kind=ckind, index=393)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve index_for_turbulent_kinetic_energy_vertical_diffusion_tracer from CCPP data structure')
            return
        end if
        if (kind(ntkev).ne.ckind) then
            call ccpp_error('Kind mismatch for variable index_for_turbulent_kinetic_energy_vertical_diffusion_tracer')
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
        

        call ccpp_field_get(cdata, 'flag_for_flux_coupling', cplflx, ierr=ierr, kind=ckind, index=278)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_flux_coupling from CCPP data structure')
            return
        end if
        if (kind(cplflx).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_flux_coupling')
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
        

        call ccpp_field_get(cdata, 'flag_idealized_physics', lsidea, ierr=ierr, kind=ckind, index=330)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_idealized_physics from CCPP data structure')
            return
        end if
        if (kind(lsidea).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_idealized_physics')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'flag_for_hedmf', hybedmf, ierr=ierr, kind=ckind, index=282)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_hedmf from CCPP data structure')
            return
        end if
        if (kind(hybedmf).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_hedmf')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'flag_for_shoc', do_shoc, ierr=ierr, kind=ckind, index=313)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_shoc from CCPP data structure')
            return
        end if
        if (kind(do_shoc).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_shoc')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'flag_for_scale_aware_TKE_moist_EDMF_PBL', satmedmf, ierr=ierr, kind=ckind, index=311)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_scale_aware_TKE_moist_EDMF_PBL from CCPP data structure')
            return
        end if
        if (kind(satmedmf).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_scale_aware_TKE_moist_EDMF_PBL')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'tendency_of_vertically_diffused_tracer_concentration', dvdftra, ierr=ierr, dims=cdims, kind=ckind, index=777)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve tendency_of_vertically_diffused_tracer_concentration from CCPP data structure')
            return
        end if
        if (kind(dvdftra).ne.ckind) then
            call ccpp_error('Kind mismatch for variable tendency_of_vertically_diffused_tracer_concentration')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'instantaneous_surface_x_momentum_flux', dusfc1, ierr=ierr, dims=cdims, kind=ckind, index=437)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve instantaneous_surface_x_momentum_flux from CCPP data structure')
            return
        end if
        if (kind(dusfc1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable instantaneous_surface_x_momentum_flux')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'instantaneous_surface_y_momentum_flux', dvsfc1, ierr=ierr, dims=cdims, kind=ckind, index=440)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve instantaneous_surface_y_momentum_flux from CCPP data structure')
            return
        end if
        if (kind(dvsfc1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable instantaneous_surface_y_momentum_flux')
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
        

        call ccpp_field_get(cdata, 'tendency_of_air_temperature_due_to_shortwave_heating_on_radiation_timestep', htrsw, ierr=ierr, dims=cdims, kind=ckind, index=763)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve tendency_of_air_temperature_due_to_shortwave_heating_on_radiation_timestep from CCPP data structure')
            return
        end if
        if (kind(htrsw).ne.ckind) then
            call ccpp_error('Kind mismatch for variable tendency_of_air_temperature_due_to_shortwave_heating_on_radiation_timestep')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'tendency_of_air_temperature_due_to_longwave_heating_on_radiation_timestep', htrlw, ierr=ierr, dims=cdims, kind=ckind, index=757)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve tendency_of_air_temperature_due_to_longwave_heating_on_radiation_timestep from CCPP data structure')
            return
        end if
        if (kind(htrlw).ne.ckind) then
            call ccpp_error('Kind mismatch for variable tendency_of_air_temperature_due_to_longwave_heating_on_radiation_timestep')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'zenith_angle_temporal_adjustment_factor_for_shortwave_fluxes', xmu, ierr=ierr, dims=cdims, kind=ckind, index=885)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve zenith_angle_temporal_adjustment_factor_for_shortwave_fluxes from CCPP data structure')
            return
        end if
        if (kind(xmu).ne.ckind) then
            call ccpp_error('Kind mismatch for variable zenith_angle_temporal_adjustment_factor_for_shortwave_fluxes')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'tendency_of_tracers_due_to_model_physics', dqdt, ierr=ierr, dims=cdims, kind=ckind, index=776)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve tendency_of_tracers_due_to_model_physics from CCPP data structure')
            return
        end if
        if (kind(dqdt).ne.ckind) then
            call ccpp_error('Kind mismatch for variable tendency_of_tracers_due_to_model_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cumulative_surface_x_momentum_flux_for_coupling_multiplied_by_timestep', dusfc_cpl, ierr=ierr, dims=cdims, kind=ckind, index=201)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cumulative_surface_x_momentum_flux_for_coupling_multiplied_by_timestep from CCPP data structure')
            return
        end if
        if (kind(dusfc_cpl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cumulative_surface_x_momentum_flux_for_coupling_multiplied_by_timestep')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cumulative_surface_y_momentum_flux_for_coupling_multiplied_by_timestep', dvsfc_cpl, ierr=ierr, dims=cdims, kind=ckind, index=203)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cumulative_surface_y_momentum_flux_for_coupling_multiplied_by_timestep from CCPP data structure')
            return
        end if
        if (kind(dvsfc_cpl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cumulative_surface_y_momentum_flux_for_coupling_multiplied_by_timestep')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cumulative_surface_upward_sensible_heat_flux_for_coupling_multiplied_by_timestep', dtsfc_cpl, ierr=ierr, dims=cdims, kind=ckind, index=198)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cumulative_surface_upward_sensible_heat_flux_for_coupling_multiplied_by_timestep from CCPP data structure')
            return
        end if
        if (kind(dtsfc_cpl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cumulative_surface_upward_sensible_heat_flux_for_coupling_multiplied_by_timestep')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cumulative_surface_upward_latent_heat_flux_for_coupling_multiplied_by_timestep', dqsfc_cpl, ierr=ierr, dims=cdims, kind=ckind, index=195)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cumulative_surface_upward_latent_heat_flux_for_coupling_multiplied_by_timestep from CCPP data structure')
            return
        end if
        if (kind(dqsfc_cpl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cumulative_surface_upward_latent_heat_flux_for_coupling_multiplied_by_timestep')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'instantaneous_surface_x_momentum_flux_for_coupling', dusfci_cpl, ierr=ierr, dims=cdims, kind=ckind, index=438)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve instantaneous_surface_x_momentum_flux_for_coupling from CCPP data structure')
            return
        end if
        if (kind(dusfci_cpl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable instantaneous_surface_x_momentum_flux_for_coupling')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'instantaneous_surface_y_momentum_flux_for_coupling', dvsfci_cpl, ierr=ierr, dims=cdims, kind=ckind, index=441)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve instantaneous_surface_y_momentum_flux_for_coupling from CCPP data structure')
            return
        end if
        if (kind(dvsfci_cpl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable instantaneous_surface_y_momentum_flux_for_coupling')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'instantaneous_surface_upward_sensible_heat_flux_for_coupling', dtsfci_cpl, ierr=ierr, dims=cdims, kind=ckind, index=435)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve instantaneous_surface_upward_sensible_heat_flux_for_coupling from CCPP data structure')
            return
        end if
        if (kind(dtsfci_cpl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable instantaneous_surface_upward_sensible_heat_flux_for_coupling')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'instantaneous_surface_upward_latent_heat_flux_for_coupling', dqsfci_cpl, ierr=ierr, dims=cdims, kind=ckind, index=432)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve instantaneous_surface_upward_latent_heat_flux_for_coupling from CCPP data structure')
            return
        end if
        if (kind(dqsfci_cpl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable instantaneous_surface_upward_latent_heat_flux_for_coupling')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cumulative_surface_x_momentum_flux_for_diag_multiplied_by_timestep', dusfc_diag, ierr=ierr, dims=cdims, kind=ckind, index=202)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cumulative_surface_x_momentum_flux_for_diag_multiplied_by_timestep from CCPP data structure')
            return
        end if
        if (kind(dusfc_diag).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cumulative_surface_x_momentum_flux_for_diag_multiplied_by_timestep')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cumulative_surface_y_momentum_flux_for_diag_multiplied_by_timestep', dvsfc_diag, ierr=ierr, dims=cdims, kind=ckind, index=204)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cumulative_surface_y_momentum_flux_for_diag_multiplied_by_timestep from CCPP data structure')
            return
        end if
        if (kind(dvsfc_diag).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cumulative_surface_y_momentum_flux_for_diag_multiplied_by_timestep')
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
        

        call ccpp_field_get(cdata, 'instantaneous_surface_x_momentum_flux_for_diag', dusfci_diag, ierr=ierr, dims=cdims, kind=ckind, index=439)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve instantaneous_surface_x_momentum_flux_for_diag from CCPP data structure')
            return
        end if
        if (kind(dusfci_diag).ne.ckind) then
            call ccpp_error('Kind mismatch for variable instantaneous_surface_x_momentum_flux_for_diag')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'instantaneous_surface_y_momentum_flux_for_diag', dvsfci_diag, ierr=ierr, dims=cdims, kind=ckind, index=442)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve instantaneous_surface_y_momentum_flux_for_diag from CCPP data structure')
            return
        end if
        if (kind(dvsfci_diag).ne.ckind) then
            call ccpp_error('Kind mismatch for variable instantaneous_surface_y_momentum_flux_for_diag')
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
        

        call ccpp_field_get(cdata, 'cumulative_change_in_temperature_due_to_PBL', dt3dt, ierr=ierr, dims=cdims, kind=ckind, index=155)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cumulative_change_in_temperature_due_to_PBL from CCPP data structure')
            return
        end if
        if (kind(dt3dt).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cumulative_change_in_temperature_due_to_PBL')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cumulative_change_in_x_wind_due_to_PBL', du3dt_PBL, ierr=ierr, dims=cdims, kind=ckind, index=165)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cumulative_change_in_x_wind_due_to_PBL from CCPP data structure')
            return
        end if
        if (kind(du3dt_PBL).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cumulative_change_in_x_wind_due_to_PBL')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cumulative_change_in_x_wind_due_to_orographic_gravity_wave_drag', du3dt_OGWD, ierr=ierr, dims=cdims, kind=ckind, index=168)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cumulative_change_in_x_wind_due_to_orographic_gravity_wave_drag from CCPP data structure')
            return
        end if
        if (kind(du3dt_OGWD).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cumulative_change_in_x_wind_due_to_orographic_gravity_wave_drag')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cumulative_change_in_y_wind_due_to_PBL', dv3dt_PBL, ierr=ierr, dims=cdims, kind=ckind, index=169)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cumulative_change_in_y_wind_due_to_PBL from CCPP data structure')
            return
        end if
        if (kind(dv3dt_PBL).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cumulative_change_in_y_wind_due_to_PBL')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cumulative_change_in_y_wind_due_to_orographic_gravity_wave_drag', dv3dt_OGWD, ierr=ierr, dims=cdims, kind=ckind, index=172)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cumulative_change_in_y_wind_due_to_orographic_gravity_wave_drag from CCPP data structure')
            return
        end if
        if (kind(dv3dt_OGWD).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cumulative_change_in_y_wind_due_to_orographic_gravity_wave_drag')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cumulative_change_in_water_vapor_specific_humidity_due_to_PBL', dq3dt, ierr=ierr, dims=cdims, kind=ckind, index=161)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cumulative_change_in_water_vapor_specific_humidity_due_to_PBL from CCPP data structure')
            return
        end if
        if (kind(dq3dt).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cumulative_change_in_water_vapor_specific_humidity_due_to_PBL')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cumulative_change_in_ozone_mixing_ratio_due_to_PBL', dq3dt_ozone, ierr=ierr, dims=cdims, kind=ckind, index=154)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cumulative_change_in_ozone_mixing_ratio_due_to_PBL from CCPP data structure')
            return
        end if
        if (kind(dq3dt_ozone).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cumulative_change_in_ozone_mixing_ratio_due_to_PBL')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call GFS_PBL_generic_post_run(im=im,levs=levs,nvdiff=nvdiff,ntrac=ntrac,ntqv=ntqv,ntcw=ntcw,ntiw=ntiw, &
                  ntrw=ntrw,ntsw=ntsw,ntlnc=ntlnc,ntinc=ntinc,ntwa=ntwa,ntia=ntia,ntgl=ntgl, &
                  ntoz=ntoz,ntke=ntke,ntkev=ntkev,imp_physics=imp_physics,imp_physics_gfdl=imp_physics_gfdl, &
                  imp_physics_thompson=imp_physics_thompson,imp_physics_wsm6=imp_physics_wsm6, &
                  ltaerosol=ltaerosol,cplflx=cplflx,lssav=lssav,ldiag3d=ldiag3d,lsidea=lsidea, &
                  hybedmf=hybedmf,do_shoc=do_shoc,satmedmf=satmedmf,dvdftra=dvdftra,dusfc1=dusfc1, &
                  dvsfc1=dvsfc1,dtsfc1=dtsfc1,dqsfc1=dqsfc1,dtf=dtf,dudt=dudt,dvdt=dvdt,dtdt=dtdt, &
                  htrsw=htrsw,htrlw=htrlw,xmu=xmu,dqdt=dqdt,dusfc_cpl=dusfc_cpl,dvsfc_cpl=dvsfc_cpl, &
                  dtsfc_cpl=dtsfc_cpl,dqsfc_cpl=dqsfc_cpl,dusfci_cpl=dusfci_cpl,dvsfci_cpl=dvsfci_cpl, &
                  dtsfci_cpl=dtsfci_cpl,dqsfci_cpl=dqsfci_cpl,dusfc_diag=dusfc_diag,dvsfc_diag=dvsfc_diag, &
                  dtsfc_diag=dtsfc_diag,dqsfc_diag=dqsfc_diag,dusfci_diag=dusfci_diag,dvsfci_diag=dvsfci_diag, &
                  dtsfci_diag=dtsfci_diag,dqsfci_diag=dqsfci_diag,dt3dt=dt3dt,du3dt_PBL=du3dt_PBL, &
                  du3dt_OGWD=du3dt_OGWD,dv3dt_PBL=dv3dt_PBL,dv3dt_OGWD=dv3dt_OGWD,dq3dt=dq3dt, &
                  dq3dt_ozone=dq3dt_ozone,errmsg=cdata%errmsg,errflg=cdata%errflg)
        ierr=cdata%errflg

    end function GFS_PBL_generic_post_run_cap
end module GFS_PBL_generic_post_cap

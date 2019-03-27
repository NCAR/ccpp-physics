
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
!! @brief Auto-generated cap module for the fv_sat_adj scheme
!!
!
module fv_sat_adj_cap

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_f_pointer, c_ptr, c_int32_t
    use            :: ccpp_types,                                      &
                      only: ccpp_t, CCPP_GENERIC_KIND
    use            :: ccpp_fields,                                     &
                      only: ccpp_field_get
    use            :: ccpp_errors,                                     &
                      only: ccpp_error, ccpp_debug
    use            :: fv_sat_adj, &
                      only: fv_sat_adj_init,fv_sat_adj_finalize,fv_sat_adj_run
    ! Other modules required, e.g. type definitions
    use machine, only: kind_dyn
    use machine, only: kind_grid

    implicit none

    private
    public :: fv_sat_adj_init_cap,fv_sat_adj_finalize_cap,fv_sat_adj_run_cap

    contains


    function fv_sat_adj_init_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind
        logical, pointer :: do_sat_adj
        integer, pointer :: kmp

        ierr = 0

        call c_f_pointer(ptr, cdata)


        call ccpp_field_get(cdata, 'flag_for_saturation_adjustment_for_microphysics_in_dynamics', do_sat_adj, ierr=ierr, kind=ckind, index=20)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_saturation_adjustment_for_microphysics_in_dynamics from CCPP data structure')
            return
        end if
        if (kind(do_sat_adj).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_saturation_adjustment_for_microphysics_in_dynamics')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'top_layer_index_for_fast_physics', kmp, ierr=ierr, kind=ckind, index=37)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve top_layer_index_for_fast_physics from CCPP data structure')
            return
        end if
        if (kind(kmp).ne.ckind) then
            call ccpp_error('Kind mismatch for variable top_layer_index_for_fast_physics')
            ierr = 1
            return
        end if
#endif
        

        call fv_sat_adj_init(do_sat_adj=do_sat_adj,kmp=kmp,errmsg=cdata%errmsg,errflg=cdata%errflg)
        ierr=cdata%errflg

    end function fv_sat_adj_init_cap

    function fv_sat_adj_finalize_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call fv_sat_adj_finalize(errmsg=cdata%errmsg,errflg=cdata%errflg)
        ierr=cdata%errflg

    end function fv_sat_adj_finalize_cap

    function fv_sat_adj_run_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind
        real(kind_dyn), pointer :: mdt
        real(kind_dyn), pointer :: zvir
        integer, pointer :: is
        integer, pointer :: ie
        integer, pointer :: isd
        integer, pointer :: ied
        integer, pointer :: kmp
        integer, pointer :: km
        integer, pointer :: kmdelz
        integer, pointer :: js
        integer, pointer :: je
        integer, pointer :: jsd
        integer, pointer :: jed
        integer, pointer :: ng
        logical, pointer :: hydrostatic
        logical, pointer :: fast_mp_consv
        real(kind_dyn), pointer :: te0_2d(:,:)
        real(kind_dyn), pointer :: te0(:,:,:)
        real(kind_dyn), pointer :: qv(:,:,:)
        real(kind_dyn), pointer :: ql(:,:,:)
        real(kind_dyn), pointer :: qi(:,:,:)
        real(kind_dyn), pointer :: qr(:,:,:)
        real(kind_dyn), pointer :: qs(:,:,:)
        real(kind_dyn), pointer :: qg(:,:,:)
        real(kind_dyn), pointer :: hs(:,:)
        real(kind_dyn), pointer :: peln(:,:,:)
        real(kind_dyn), pointer :: delz(:,:,:)
        real(kind_dyn), pointer :: delp(:,:,:)
        real(kind_dyn), pointer :: pt(:,:,:)
        real(kind_dyn), pointer :: pkz(:,:,:)
        real(kind_dyn), pointer :: q_con(:,:,:)
        real(kind_dyn), pointer :: akap
        real(kind_dyn), pointer :: cappa(:,:,:)
        real(kind_grid), pointer :: area(:,:)
        real(kind_dyn), pointer :: dtdt(:,:,:)
        logical, pointer :: out_dt
        logical, pointer :: last_step
        logical, pointer :: do_qa
        real(kind_dyn), pointer :: qa(:,:,:)
        integer, pointer :: nthreads

        ierr = 0

        call c_f_pointer(ptr, cdata)


        call ccpp_field_get(cdata, 'time_step_for_remapping_for_fast_physics', mdt, ierr=ierr, kind=ckind, index=36)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve time_step_for_remapping_for_fast_physics from CCPP data structure')
            return
        end if
        if (kind(mdt).ne.ckind) then
            call ccpp_error('Kind mismatch for variable time_step_for_remapping_for_fast_physics')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'ratio_of_vapor_to_dry_air_gas_constants_minus_one_default_kind', zvir, ierr=ierr, kind=ckind, index=28)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve ratio_of_vapor_to_dry_air_gas_constants_minus_one_default_kind from CCPP data structure')
            return
        end if
        if (kind(zvir).ne.ckind) then
            call ccpp_error('Kind mismatch for variable ratio_of_vapor_to_dry_air_gas_constants_minus_one_default_kind')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'starting_x_direction_index', is, ierr=ierr, kind=ckind, index=29)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve starting_x_direction_index from CCPP data structure')
            return
        end if
        if (kind(is).ne.ckind) then
            call ccpp_error('Kind mismatch for variable starting_x_direction_index')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'ending_x_direction_index', ie, ierr=ierr, kind=ckind, index=12)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve ending_x_direction_index from CCPP data structure')
            return
        end if
        if (kind(ie).ne.ckind) then
            call ccpp_error('Kind mismatch for variable ending_x_direction_index')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'starting_x_direction_index_domain', isd, ierr=ierr, kind=ckind, index=30)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve starting_x_direction_index_domain from CCPP data structure')
            return
        end if
        if (kind(isd).ne.ckind) then
            call ccpp_error('Kind mismatch for variable starting_x_direction_index_domain')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'ending_x_direction_index_domain', ied, ierr=ierr, kind=ckind, index=13)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve ending_x_direction_index_domain from CCPP data structure')
            return
        end if
        if (kind(ied).ne.ckind) then
            call ccpp_error('Kind mismatch for variable ending_x_direction_index_domain')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'top_layer_index_for_fast_physics', kmp, ierr=ierr, kind=ckind, index=37)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve top_layer_index_for_fast_physics from CCPP data structure')
            return
        end if
        if (kind(kmp).ne.ckind) then
            call ccpp_error('Kind mismatch for variable top_layer_index_for_fast_physics')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'vertical_dimension_for_fast_physics', km, ierr=ierr, kind=ckind, index=38)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve vertical_dimension_for_fast_physics from CCPP data structure')
            return
        end if
        if (kind(km).ne.ckind) then
            call ccpp_error('Kind mismatch for variable vertical_dimension_for_fast_physics')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'vertical_dimension_for_thickness_at_Lagrangian_surface', kmdelz, ierr=ierr, kind=ckind, index=39)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve vertical_dimension_for_thickness_at_Lagrangian_surface from CCPP data structure')
            return
        end if
        if (kind(kmdelz).ne.ckind) then
            call ccpp_error('Kind mismatch for variable vertical_dimension_for_thickness_at_Lagrangian_surface')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'starting_y_direction_index', js, ierr=ierr, kind=ckind, index=31)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve starting_y_direction_index from CCPP data structure')
            return
        end if
        if (kind(js).ne.ckind) then
            call ccpp_error('Kind mismatch for variable starting_y_direction_index')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'ending_y_direction_index', je, ierr=ierr, kind=ckind, index=14)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve ending_y_direction_index from CCPP data structure')
            return
        end if
        if (kind(je).ne.ckind) then
            call ccpp_error('Kind mismatch for variable ending_y_direction_index')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'starting_y_direction_index_domain', jsd, ierr=ierr, kind=ckind, index=32)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve starting_y_direction_index_domain from CCPP data structure')
            return
        end if
        if (kind(jsd).ne.ckind) then
            call ccpp_error('Kind mismatch for variable starting_y_direction_index_domain')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'ending_y_direction_index_domain', jed, ierr=ierr, kind=ckind, index=15)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve ending_y_direction_index_domain from CCPP data structure')
            return
        end if
        if (kind(jed).ne.ckind) then
            call ccpp_error('Kind mismatch for variable ending_y_direction_index_domain')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'number_of_ghost_zones', ng, ierr=ierr, kind=ckind, index=25)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve number_of_ghost_zones from CCPP data structure')
            return
        end if
        if (kind(ng).ne.ckind) then
            call ccpp_error('Kind mismatch for variable number_of_ghost_zones')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'flag_for_hydrostatic_solver', hydrostatic, ierr=ierr, kind=ckind, index=18)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_hydrostatic_solver from CCPP data structure')
            return
        end if
        if (kind(hydrostatic).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_hydrostatic_solver')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'flag_for_fast_microphysics_energy_conservation', fast_mp_consv, ierr=ierr, kind=ckind, index=17)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_fast_microphysics_energy_conservation from CCPP data structure')
            return
        end if
        if (kind(fast_mp_consv).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_fast_microphysics_energy_conservation')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'atmosphere_energy_content_in_column', te0_2d, ierr=ierr, dims=cdims, kind=ckind, index=2)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve atmosphere_energy_content_in_column from CCPP data structure')
            return
        end if
        if (kind(te0_2d).ne.ckind) then
            call ccpp_error('Kind mismatch for variable atmosphere_energy_content_in_column')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'atmosphere_energy_content_at_Lagrangian_surface', te0, ierr=ierr, dims=cdims, kind=ckind, index=1)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve atmosphere_energy_content_at_Lagrangian_surface from CCPP data structure')
            return
        end if
        if (kind(te0).ne.ckind) then
            call ccpp_error('Kind mismatch for variable atmosphere_energy_content_at_Lagrangian_surface')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'water_vapor_specific_humidity_at_Lagrangian_surface', qv, ierr=ierr, dims=cdims, kind=ckind, index=41)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve water_vapor_specific_humidity_at_Lagrangian_surface from CCPP data structure')
            return
        end if
        if (kind(qv).ne.ckind) then
            call ccpp_error('Kind mismatch for variable water_vapor_specific_humidity_at_Lagrangian_surface')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cloud_liquid_water_specific_humidity_at_Lagrangian_surface', ql, ierr=ierr, dims=cdims, kind=ckind, index=9)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cloud_liquid_water_specific_humidity_at_Lagrangian_surface from CCPP data structure')
            return
        end if
        if (kind(ql).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cloud_liquid_water_specific_humidity_at_Lagrangian_surface')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cloud_ice_specific_humidity_at_Lagrangian_surface', qi, ierr=ierr, dims=cdims, kind=ckind, index=8)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cloud_ice_specific_humidity_at_Lagrangian_surface from CCPP data structure')
            return
        end if
        if (kind(qi).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cloud_ice_specific_humidity_at_Lagrangian_surface')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cloud_rain_specific_humidity_at_Lagrangian_surface', qr, ierr=ierr, dims=cdims, kind=ckind, index=10)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cloud_rain_specific_humidity_at_Lagrangian_surface from CCPP data structure')
            return
        end if
        if (kind(qr).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cloud_rain_specific_humidity_at_Lagrangian_surface')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cloud_snow_specific_humidity_at_Lagrangian_surface', qs, ierr=ierr, dims=cdims, kind=ckind, index=11)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cloud_snow_specific_humidity_at_Lagrangian_surface from CCPP data structure')
            return
        end if
        if (kind(qs).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cloud_snow_specific_humidity_at_Lagrangian_surface')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cloud_graupel_specific_humidity_at_Lagrangian_surface', qg, ierr=ierr, dims=cdims, kind=ckind, index=7)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cloud_graupel_specific_humidity_at_Lagrangian_surface from CCPP data structure')
            return
        end if
        if (kind(qg).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cloud_graupel_specific_humidity_at_Lagrangian_surface')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_geopotential_at_Lagrangian_surface', hs, ierr=ierr, dims=cdims, kind=ckind, index=33)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_geopotential_at_Lagrangian_surface from CCPP data structure')
            return
        end if
        if (kind(hs).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_geopotential_at_Lagrangian_surface')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'log_pressure_at_Lagrangian_surface', peln, ierr=ierr, dims=cdims, kind=ckind, index=24)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve log_pressure_at_Lagrangian_surface from CCPP data structure')
            return
        end if
        if (kind(peln).ne.ckind) then
            call ccpp_error('Kind mismatch for variable log_pressure_at_Lagrangian_surface')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'thickness_at_Lagrangian_surface', delz, ierr=ierr, dims=cdims, kind=ckind, index=35)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve thickness_at_Lagrangian_surface from CCPP data structure')
            return
        end if
        if (kind(delz).ne.ckind) then
            call ccpp_error('Kind mismatch for variable thickness_at_Lagrangian_surface')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'pressure_thickness_at_Lagrangian_surface', delp, ierr=ierr, dims=cdims, kind=ckind, index=27)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve pressure_thickness_at_Lagrangian_surface from CCPP data structure')
            return
        end if
        if (kind(delp).ne.ckind) then
            call ccpp_error('Kind mismatch for variable pressure_thickness_at_Lagrangian_surface')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'virtual_temperature_at_Lagrangian_surface', pt, ierr=ierr, dims=cdims, kind=ckind, index=40)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve virtual_temperature_at_Lagrangian_surface from CCPP data structure')
            return
        end if
        if (kind(pt).ne.ckind) then
            call ccpp_error('Kind mismatch for variable virtual_temperature_at_Lagrangian_surface')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'finite-volume_mean_edge_pressure_raised_to_the_power_of_kappa', pkz, ierr=ierr, dims=cdims, kind=ckind, index=16)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve finite-volume_mean_edge_pressure_raised_to_the_power_of_kappa from CCPP data structure')
            return
        end if
        if (kind(pkz).ne.ckind) then
            call ccpp_error('Kind mismatch for variable finite-volume_mean_edge_pressure_raised_to_the_power_of_kappa')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cloud_condensed_water_specific_humidity_at_Lagrangian_surface', q_con, ierr=ierr, dims=cdims, kind=ckind, index=5)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cloud_condensed_water_specific_humidity_at_Lagrangian_surface from CCPP data structure')
            return
        end if
        if (kind(q_con).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cloud_condensed_water_specific_humidity_at_Lagrangian_surface')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'kappa_dry_for_fast_physics', akap, ierr=ierr, kind=ckind, index=23)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve kappa_dry_for_fast_physics from CCPP data structure')
            return
        end if
        if (kind(akap).ne.ckind) then
            call ccpp_error('Kind mismatch for variable kappa_dry_for_fast_physics')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'cappa_moist_gas_constant_at_Lagrangian_surface', cappa, ierr=ierr, dims=cdims, kind=ckind, index=3)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cappa_moist_gas_constant_at_Lagrangian_surface from CCPP data structure')
            return
        end if
        if (kind(cappa).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cappa_moist_gas_constant_at_Lagrangian_surface')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cell_area_for_fast_physics', area, ierr=ierr, dims=cdims, kind=ckind, index=4)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cell_area_for_fast_physics from CCPP data structure')
            return
        end if
        if (kind(area).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cell_area_for_fast_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'tendency_of_air_temperature_at_Lagrangian_surface', dtdt, ierr=ierr, dims=cdims, kind=ckind, index=34)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve tendency_of_air_temperature_at_Lagrangian_surface from CCPP data structure')
            return
        end if
        if (kind(dtdt).ne.ckind) then
            call ccpp_error('Kind mismatch for variable tendency_of_air_temperature_at_Lagrangian_surface')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'flag_for_tendency_of_air_temperature_at_Lagrangian_surface', out_dt, ierr=ierr, kind=ckind, index=21)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_tendency_of_air_temperature_at_Lagrangian_surface from CCPP data structure')
            return
        end if
        if (kind(out_dt).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_tendency_of_air_temperature_at_Lagrangian_surface')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'flag_for_the_last_step_of_k_split_remapping', last_step, ierr=ierr, kind=ckind, index=22)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_the_last_step_of_k_split_remapping from CCPP data structure')
            return
        end if
        if (kind(last_step).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_the_last_step_of_k_split_remapping')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'flag_for_inline_cloud_fraction_calculation', do_qa, ierr=ierr, kind=ckind, index=19)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_inline_cloud_fraction_calculation from CCPP data structure')
            return
        end if
        if (kind(do_qa).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_inline_cloud_fraction_calculation')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'cloud_fraction_at_Lagrangian_surface', qa, ierr=ierr, dims=cdims, kind=ckind, index=6)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cloud_fraction_at_Lagrangian_surface from CCPP data structure')
            return
        end if
        if (kind(qa).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cloud_fraction_at_Lagrangian_surface')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'omp_threads', nthreads, ierr=ierr, kind=ckind, index=26)
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
        

        call fv_sat_adj_run(mdt=mdt,zvir=zvir,is=is,ie=ie,isd=isd,ied=ied,kmp=kmp,km=km,kmdelz=kmdelz, &
                  js=js,je=je,jsd=jsd,jed=jed,ng=ng,hydrostatic=hydrostatic,fast_mp_consv=fast_mp_consv, &
                  te0_2d=te0_2d,te0=te0,qv=qv,ql=ql,qi=qi,qr=qr,qs=qs,qg=qg,hs=hs,peln=peln, &
                  delz=delz,delp=delp,pt=pt,pkz=pkz,q_con=q_con,akap=akap,cappa=cappa,area=area, &
                  dtdt=dtdt,out_dt=out_dt,last_step=last_step,do_qa=do_qa,qa=qa,nthreads=nthreads, &
                  errmsg=cdata%errmsg,errflg=cdata%errflg)
        ierr=cdata%errflg

    end function fv_sat_adj_run_cap
end module fv_sat_adj_cap

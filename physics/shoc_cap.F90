
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
!! @brief Auto-generated cap module for the shoc scheme
!!
!
module shoc_cap

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_f_pointer, c_ptr, c_int32_t
    use            :: ccpp_types,                                      &
                      only: ccpp_t, CCPP_GENERIC_KIND
    use            :: ccpp_fields,                                     &
                      only: ccpp_field_get
    use            :: ccpp_errors,                                     &
                      only: ccpp_error, ccpp_debug
    use            :: shoc, &
                      only: shoc_finalize,shoc_init,shoc_run
    ! Other modules required, e.g. type definitions
    use machine, only: kind_phys

    implicit none

    private
    public :: shoc_finalize_cap,shoc_init_cap,shoc_run_cap

    contains


    function shoc_finalize_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call shoc_finalize()
        

    end function shoc_finalize_cap

    function shoc_init_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call shoc_init()
        

    end function shoc_init_cap

    function shoc_run_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind
        integer, pointer :: ix
        integer, pointer :: nx
        integer, pointer :: nzm
        logical, pointer :: do_shoc
        logical, pointer :: shocaftcnv
        logical, pointer :: mg3_as_mg2
        integer, pointer :: imp_physics
        integer, pointer :: imp_physics_gfdl
        integer, pointer :: imp_physics_zhao_carr
        integer, pointer :: imp_physics_zhao_carr_pdf
        integer, pointer :: imp_physics_mg
        integer, pointer :: fprcp
        real(kind_phys), pointer :: tcr
        real(kind_phys), pointer :: tcrf
        real(kind_phys), pointer :: con_cp
        real(kind_phys), pointer :: con_g
        real(kind_phys), pointer :: con_hvap
        real(kind_phys), pointer :: con_hfus
        real(kind_phys), pointer :: con_rv
        real(kind_phys), pointer :: con_rd
        real(kind_phys), pointer :: con_pi
        real(kind_phys), pointer :: con_fvirt
        real(kind_phys), pointer :: gq0_cloud_ice(:,:)
        real(kind_phys), pointer :: gq0_rain(:,:)
        real(kind_phys), pointer :: gq0_snow(:,:)
        real(kind_phys), pointer :: gq0_graupel(:,:)
        real(kind_phys), pointer :: dtp
        integer, pointer :: me
        real(kind_phys), pointer :: prsl(:,:)
        real(kind_phys), pointer :: phii(:,:)
        real(kind_phys), pointer :: phil(:,:)
        real(kind_phys), pointer :: u(:,:)
        real(kind_phys), pointer :: v(:,:)
        real(kind_phys), pointer :: omega(:,:)
        real(kind_phys), pointer :: rhc(:,:)
        real(kind_phys), pointer :: supice
        real(kind_phys), pointer :: pcrit
        real(kind_phys), pointer :: cefac
        real(kind_phys), pointer :: cesfac
        real(kind_phys), pointer :: tkef1
        real(kind_phys), pointer :: dis_opt
        real(kind_phys), pointer :: hflx(:)
        real(kind_phys), pointer :: evap(:)
        real(kind_phys), pointer :: prnum(:,:)
        logical, pointer :: skip_macro
        real(kind_phys), pointer :: clw_ice(:,:)
        real(kind_phys), pointer :: clw_liquid(:,:)
        real(kind_phys), pointer :: gq0_cloud_liquid(:,:)
        real(kind_phys), pointer :: ncpl(:,:)
        real(kind_phys), pointer :: ncpi(:,:)
        real(kind_phys), pointer :: gt0(:,:)
        real(kind_phys), pointer :: gq0_water_vapor(:,:)
        real(kind_phys), pointer :: cld_sgs(:,:)
        real(kind_phys), pointer :: tke(:,:)
        real(kind_phys), pointer :: tkh(:,:)
        real(kind_phys), pointer :: wthv_sec(:,:)

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
        

        call ccpp_field_get(cdata, 'horizontal_loop_extent', nx, ierr=ierr, kind=ckind, index=366)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve horizontal_loop_extent from CCPP data structure')
            return
        end if
        if (kind(nx).ne.ckind) then
            call ccpp_error('Kind mismatch for variable horizontal_loop_extent')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'vertical_dimension', nzm, ierr=ierr, kind=ckind, index=817)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve vertical_dimension from CCPP data structure')
            return
        end if
        if (kind(nzm).ne.ckind) then
            call ccpp_error('Kind mismatch for variable vertical_dimension')
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
        

        call ccpp_field_get(cdata, 'flag_for_shoc_after_convection', shocaftcnv, ierr=ierr, kind=ckind, index=314)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_shoc_after_convection from CCPP data structure')
            return
        end if
        if (kind(shocaftcnv).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_shoc_after_convection')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'flag_mg3_as_mg2', mg3_as_mg2, ierr=ierr, kind=ckind, index=331)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_mg3_as_mg2 from CCPP data structure')
            return
        end if
        if (kind(mg3_as_mg2).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_mg3_as_mg2')
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
        

        call ccpp_field_get(cdata, 'flag_for_zhao_carr_microphysics_scheme', imp_physics_zhao_carr, ierr=ierr, kind=ckind, index=327)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_zhao_carr_microphysics_scheme from CCPP data structure')
            return
        end if
        if (kind(imp_physics_zhao_carr).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_zhao_carr_microphysics_scheme')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'flag_for_zhao_carr_pdf_microphysics_scheme', imp_physics_zhao_carr_pdf, ierr=ierr, kind=ckind, index=328)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_zhao_carr_pdf_microphysics_scheme from CCPP data structure')
            return
        end if
        if (kind(imp_physics_zhao_carr_pdf).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_zhao_carr_pdf_microphysics_scheme')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'flag_for_morrison_gettelman_microphysics_scheme', imp_physics_mg, ierr=ierr, kind=ckind, index=297)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_morrison_gettelman_microphysics_scheme from CCPP data structure')
            return
        end if
        if (kind(imp_physics_mg).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_morrison_gettelman_microphysics_scheme')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'number_of_frozen_precipitation_species', fprcp, ierr=ierr, kind=ckind, index=577)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve number_of_frozen_precipitation_species from CCPP data structure')
            return
        end if
        if (kind(fprcp).ne.ckind) then
            call ccpp_error('Kind mismatch for variable number_of_frozen_precipitation_species')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'cloud_phase_transition_threshold_temperature', tcr, ierr=ierr, kind=ckind, index=106)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cloud_phase_transition_threshold_temperature from CCPP data structure')
            return
        end if
        if (kind(tcr).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cloud_phase_transition_threshold_temperature')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'cloud_phase_transition_denominator', tcrf, ierr=ierr, kind=ckind, index=105)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cloud_phase_transition_denominator from CCPP data structure')
            return
        end if
        if (kind(tcrf).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cloud_phase_transition_denominator')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'specific_heat_of_dry_air_at_constant_pressure', con_cp, ierr=ierr, kind=ckind, index=664)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve specific_heat_of_dry_air_at_constant_pressure from CCPP data structure')
            return
        end if
        if (kind(con_cp).ne.ckind) then
            call ccpp_error('Kind mismatch for variable specific_heat_of_dry_air_at_constant_pressure')
            ierr = 1
            return
        end if
#endif
        

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
        

        call ccpp_field_get(cdata, 'latent_heat_of_vaporization_of_water_at_0C', con_hvap, ierr=ierr, kind=ckind, index=459)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve latent_heat_of_vaporization_of_water_at_0C from CCPP data structure')
            return
        end if
        if (kind(con_hvap).ne.ckind) then
            call ccpp_error('Kind mismatch for variable latent_heat_of_vaporization_of_water_at_0C')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'latent_heat_of_fusion_of_water_at_0C', con_hfus, ierr=ierr, kind=ckind, index=458)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve latent_heat_of_fusion_of_water_at_0C from CCPP data structure')
            return
        end if
        if (kind(con_hfus).ne.ckind) then
            call ccpp_error('Kind mismatch for variable latent_heat_of_fusion_of_water_at_0C')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'gas_constant_water_vapor', con_rv, ierr=ierr, kind=ckind, index=347)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve gas_constant_water_vapor from CCPP data structure')
            return
        end if
        if (kind(con_rv).ne.ckind) then
            call ccpp_error('Kind mismatch for variable gas_constant_water_vapor')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'gas_constant_dry_air', con_rd, ierr=ierr, kind=ckind, index=346)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve gas_constant_dry_air from CCPP data structure')
            return
        end if
        if (kind(con_rd).ne.ckind) then
            call ccpp_error('Kind mismatch for variable gas_constant_dry_air')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'pi', con_pi, ierr=ierr, kind=ckind, index=606)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve pi from CCPP data structure')
            return
        end if
        if (kind(con_pi).ne.ckind) then
            call ccpp_error('Kind mismatch for variable pi')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'ratio_of_vapor_to_dry_air_gas_constants_minus_one', con_fvirt, ierr=ierr, kind=ckind, index=625)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve ratio_of_vapor_to_dry_air_gas_constants_minus_one from CCPP data structure')
            return
        end if
        if (kind(con_fvirt).ne.ckind) then
            call ccpp_error('Kind mismatch for variable ratio_of_vapor_to_dry_air_gas_constants_minus_one')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'ice_water_mixing_ratio_updated_by_physics', gq0_cloud_ice, ierr=ierr, dims=cdims, kind=ckind, index=376)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve ice_water_mixing_ratio_updated_by_physics from CCPP data structure')
            return
        end if
        if (kind(gq0_cloud_ice).ne.ckind) then
            call ccpp_error('Kind mismatch for variable ice_water_mixing_ratio_updated_by_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'rain_water_mixing_ratio_updated_by_physics', gq0_rain, ierr=ierr, dims=cdims, kind=ckind, index=619)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve rain_water_mixing_ratio_updated_by_physics from CCPP data structure')
            return
        end if
        if (kind(gq0_rain).ne.ckind) then
            call ccpp_error('Kind mismatch for variable rain_water_mixing_ratio_updated_by_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'snow_water_mixing_ratio_updated_by_physics', gq0_snow, ierr=ierr, dims=cdims, kind=ckind, index=653)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve snow_water_mixing_ratio_updated_by_physics from CCPP data structure')
            return
        end if
        if (kind(gq0_snow).ne.ckind) then
            call ccpp_error('Kind mismatch for variable snow_water_mixing_ratio_updated_by_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'graupel_mixing_ratio_updated_by_physics', gq0_graupel, ierr=ierr, dims=cdims, kind=ckind, index=352)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve graupel_mixing_ratio_updated_by_physics from CCPP data structure')
            return
        end if
        if (kind(gq0_graupel).ne.ckind) then
            call ccpp_error('Kind mismatch for variable graupel_mixing_ratio_updated_by_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'time_step_for_physics', dtp, ierr=ierr, kind=ckind, index=793)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve time_step_for_physics from CCPP data structure')
            return
        end if
        if (kind(dtp).ne.ckind) then
            call ccpp_error('Kind mismatch for variable time_step_for_physics')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'mpi_rank', me, ierr=ierr, kind=ckind, index=558)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve mpi_rank from CCPP data structure')
            return
        end if
        if (kind(me).ne.ckind) then
            call ccpp_error('Kind mismatch for variable mpi_rank')
            ierr = 1
            return
        end if
#endif
        

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
        

        call ccpp_field_get(cdata, 'geopotential', phil, ierr=ierr, dims=cdims, kind=ckind, index=348)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve geopotential from CCPP data structure')
            return
        end if
        if (kind(phil).ne.ckind) then
            call ccpp_error('Kind mismatch for variable geopotential')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'x_wind_updated_by_physics', u, ierr=ierr, dims=cdims, kind=ckind, index=877)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve x_wind_updated_by_physics from CCPP data structure')
            return
        end if
        if (kind(u).ne.ckind) then
            call ccpp_error('Kind mismatch for variable x_wind_updated_by_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'y_wind_updated_by_physics', v, ierr=ierr, dims=cdims, kind=ckind, index=884)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve y_wind_updated_by_physics from CCPP data structure')
            return
        end if
        if (kind(v).ne.ckind) then
            call ccpp_error('Kind mismatch for variable y_wind_updated_by_physics')
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
        

        call ccpp_field_get(cdata, 'critical_relative_humidity', rhc, ierr=ierr, dims=cdims, kind=ckind, index=141)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve critical_relative_humidity from CCPP data structure')
            return
        end if
        if (kind(rhc).ne.ckind) then
            call ccpp_error('Kind mismatch for variable critical_relative_humidity')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'ice_supersaturation_threshold', supice, ierr=ierr, kind=ckind, index=372)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve ice_supersaturation_threshold from CCPP data structure')
            return
        end if
        if (kind(supice).ne.ckind) then
            call ccpp_error('Kind mismatch for variable ice_supersaturation_threshold')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'shoc_tke_dissipatation_pressure_threshold', pcrit, ierr=ierr, kind=ckind, index=642)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve shoc_tke_dissipatation_pressure_threshold from CCPP data structure')
            return
        end if
        if (kind(pcrit).ne.ckind) then
            call ccpp_error('Kind mismatch for variable shoc_tke_dissipatation_pressure_threshold')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'shoc_tke_dissipation_tunable_parameter', cefac, ierr=ierr, kind=ckind, index=643)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve shoc_tke_dissipation_tunable_parameter from CCPP data structure')
            return
        end if
        if (kind(cefac).ne.ckind) then
            call ccpp_error('Kind mismatch for variable shoc_tke_dissipation_tunable_parameter')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'shoc_tke_dissipation_tunable_parameter_near_surface', cesfac, ierr=ierr, kind=ckind, index=644)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve shoc_tke_dissipation_tunable_parameter_near_surface from CCPP data structure')
            return
        end if
        if (kind(cesfac).ne.ckind) then
            call ccpp_error('Kind mismatch for variable shoc_tke_dissipation_tunable_parameter_near_surface')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'shoc_implicit_TKE_integration_uncentering_term', tkef1, ierr=ierr, kind=ckind, index=641)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve shoc_implicit_TKE_integration_uncentering_term from CCPP data structure')
            return
        end if
        if (kind(tkef1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable shoc_implicit_TKE_integration_uncentering_term')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'shoc_flag_for_optional_surface_TKE_dissipation', dis_opt, ierr=ierr, kind=ckind, index=640)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve shoc_flag_for_optional_surface_TKE_dissipation from CCPP data structure')
            return
        end if
        if (kind(dis_opt).ne.ckind) then
            call ccpp_error('Kind mismatch for variable shoc_flag_for_optional_surface_TKE_dissipation')
            ierr = 1
            return
        end if
#endif
        

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
        

        call ccpp_field_get(cdata, 'kinematic_surface_upward_latent_heat_flux', evap, ierr=ierr, dims=cdims, kind=ckind, index=454)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve kinematic_surface_upward_latent_heat_flux from CCPP data structure')
            return
        end if
        if (kind(evap).ne.ckind) then
            call ccpp_error('Kind mismatch for variable kinematic_surface_upward_latent_heat_flux')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'prandtl_number', prnum, ierr=ierr, dims=cdims, kind=ckind, index=608)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve prandtl_number from CCPP data structure')
            return
        end if
        if (kind(prnum).ne.ckind) then
            call ccpp_error('Kind mismatch for variable prandtl_number')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'flag_skip_macro', skip_macro, ierr=ierr, kind=ckind, index=334)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_skip_macro from CCPP data structure')
            return
        end if
        if (kind(skip_macro).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_skip_macro')
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
        

        call ccpp_field_get(cdata, 'cloud_condensed_water_mixing_ratio_updated_by_physics', gq0_cloud_liquid, ierr=ierr, dims=cdims, kind=ckind, index=95)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cloud_condensed_water_mixing_ratio_updated_by_physics from CCPP data structure')
            return
        end if
        if (kind(gq0_cloud_liquid).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cloud_condensed_water_mixing_ratio_updated_by_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cloud_droplet_number_concentration_updated_by_physics', ncpl, ierr=ierr, dims=cdims, kind=ckind, index=98)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cloud_droplet_number_concentration_updated_by_physics from CCPP data structure')
            return
        end if
        if (kind(ncpl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cloud_droplet_number_concentration_updated_by_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'ice_number_concentration_updated_by_physics', ncpi, ierr=ierr, dims=cdims, kind=ckind, index=371)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve ice_number_concentration_updated_by_physics from CCPP data structure')
            return
        end if
        if (kind(ncpi).ne.ckind) then
            call ccpp_error('Kind mismatch for variable ice_number_concentration_updated_by_physics')
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
        

        call ccpp_field_get(cdata, 'subgrid_scale_cloud_fraction_from_shoc', cld_sgs, ierr=ierr, dims=cdims, kind=ckind, index=675)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve subgrid_scale_cloud_fraction_from_shoc from CCPP data structure')
            return
        end if
        if (kind(cld_sgs).ne.ckind) then
            call ccpp_error('Kind mismatch for variable subgrid_scale_cloud_fraction_from_shoc')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'turbulent_kinetic_energy_convective_transport_tracer', tke, ierr=ierr, dims=cdims, kind=ckind, index=808)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve turbulent_kinetic_energy_convective_transport_tracer from CCPP data structure')
            return
        end if
        if (kind(tke).ne.ckind) then
            call ccpp_error('Kind mismatch for variable turbulent_kinetic_energy_convective_transport_tracer')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'atmosphere_heat_diffusivity_from_shoc', tkh, ierr=ierr, dims=cdims, kind=ckind, index=72)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve atmosphere_heat_diffusivity_from_shoc from CCPP data structure')
            return
        end if
        if (kind(tkh).ne.ckind) then
            call ccpp_error('Kind mismatch for variable atmosphere_heat_diffusivity_from_shoc')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'kinematic_buoyancy_flux_from_shoc', wthv_sec, ierr=ierr, dims=cdims, kind=ckind, index=453)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve kinematic_buoyancy_flux_from_shoc from CCPP data structure')
            return
        end if
        if (kind(wthv_sec).ne.ckind) then
            call ccpp_error('Kind mismatch for variable kinematic_buoyancy_flux_from_shoc')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call shoc_run(ix=ix,nx=nx,nzm=nzm,do_shoc=do_shoc,shocaftcnv=shocaftcnv,mg3_as_mg2=mg3_as_mg2, &
                  imp_physics=imp_physics,imp_physics_gfdl=imp_physics_gfdl,imp_physics_zhao_carr=imp_physics_zhao_carr, &
                  imp_physics_zhao_carr_pdf=imp_physics_zhao_carr_pdf,imp_physics_mg=imp_physics_mg, &
                  fprcp=fprcp,tcr=tcr,tcrf=tcrf,con_cp=con_cp,con_g=con_g,con_hvap=con_hvap, &
                  con_hfus=con_hfus,con_rv=con_rv,con_rd=con_rd,con_pi=con_pi,con_fvirt=con_fvirt, &
                  gq0_cloud_ice=gq0_cloud_ice,gq0_rain=gq0_rain,gq0_snow=gq0_snow,gq0_graupel=gq0_graupel, &
                  dtp=dtp,me=me,prsl=prsl,phii=phii,phil=phil,u=u,v=v,omega=omega,rhc=rhc, &
                  supice=supice,pcrit=pcrit,cefac=cefac,cesfac=cesfac,tkef1=tkef1,dis_opt=dis_opt, &
                  hflx=hflx,evap=evap,prnum=prnum,skip_macro=skip_macro,clw_ice=clw_ice,clw_liquid=clw_liquid, &
                  gq0_cloud_liquid=gq0_cloud_liquid,ncpl=ncpl,ncpi=ncpi,gt0=gt0,gq0_water_vapor=gq0_water_vapor, &
                  cld_sgs=cld_sgs,tke=tke,tkh=tkh,wthv_sec=wthv_sec,errmsg=cdata%errmsg,errflg=cdata%errflg)
        ierr=cdata%errflg

    end function shoc_run_cap
end module shoc_cap


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
!! @brief Auto-generated cap module for the GFS_MP_generic_post scheme
!!
!
module GFS_MP_generic_post_cap

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_f_pointer, c_ptr, c_int32_t
    use            :: ccpp_types,                                      &
                      only: ccpp_t, CCPP_GENERIC_KIND
    use            :: ccpp_fields,                                     &
                      only: ccpp_field_get
    use            :: ccpp_errors,                                     &
                      only: ccpp_error, ccpp_debug
    use            :: GFS_MP_generic_post, &
                      only: GFS_MP_generic_post_init,GFS_MP_generic_post_run,GFS_MP_generic_post_finalize
    ! Other modules required, e.g. type definitions
    use machine, only: kind_phys

    implicit none

    private
    public :: GFS_MP_generic_post_init_cap,GFS_MP_generic_post_run_cap,GFS_MP_generic_post_finalize_cap

    contains


    function GFS_MP_generic_post_init_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call GFS_MP_generic_post_init()
        

    end function GFS_MP_generic_post_init_cap

    function GFS_MP_generic_post_run_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind
        integer, pointer :: im
        integer, pointer :: ix
        integer, pointer :: levs
        integer, pointer :: kdt
        integer, pointer :: nrcm
        integer, pointer :: ncld
        integer, pointer :: nncl
        integer, pointer :: ntcw
        integer, pointer :: ntrac
        integer, pointer :: imp_physics
        integer, pointer :: imp_physics_gfdl
        integer, pointer :: imp_physics_thompson
        logical, pointer :: cal_pre
        logical, pointer :: lssav
        logical, pointer :: ldiag3d
        logical, pointer :: cplflx
        logical, pointer :: cplchm
        real(kind_phys), pointer :: con_g
        real(kind_phys), pointer :: dtf
        real(kind_phys), pointer :: frain
        real(kind_phys), pointer :: rainc(:)
        real(kind_phys), pointer :: rain1(:)
        real(kind_phys), pointer :: rann(:,:)
        real(kind_phys), pointer :: xlat(:)
        real(kind_phys), pointer :: xlon(:)
        real(kind_phys), pointer :: gt0(:,:)
        real(kind_phys), pointer :: gq0(:,:,:)
        real(kind_phys), pointer :: prsl(:,:)
        real(kind_phys), pointer :: prsi(:,:)
        real(kind_phys), pointer :: phii(:,:)
        real(kind_phys), pointer :: tsfc(:)
        real(kind_phys), pointer :: ice(:)
        real(kind_phys), pointer :: snow(:)
        real(kind_phys), pointer :: graupel(:)
        real(kind_phys), pointer :: save_t(:,:)
        real(kind_phys), pointer :: save_qv(:,:)
        real(kind_phys), pointer :: ice0(:)
        real(kind_phys), pointer :: snow0(:)
        real(kind_phys), pointer :: graupel0(:)
        real(kind_phys), pointer :: del(:,:)
        real(kind_phys), pointer :: rain(:)
        real(kind_phys), pointer :: domr_diag(:)
        real(kind_phys), pointer :: domzr_diag(:)
        real(kind_phys), pointer :: domip_diag(:)
        real(kind_phys), pointer :: doms_diag(:)
        real(kind_phys), pointer :: tprcp(:)
        real(kind_phys), pointer :: srflag(:)
        real(kind_phys), pointer :: totprcp(:)
        real(kind_phys), pointer :: totice(:)
        real(kind_phys), pointer :: totsnw(:)
        real(kind_phys), pointer :: totgrp(:)
        real(kind_phys), pointer :: totprcpb(:)
        real(kind_phys), pointer :: toticeb(:)
        real(kind_phys), pointer :: totsnwb(:)
        real(kind_phys), pointer :: totgrpb(:)
        real(kind_phys), pointer :: dt3dt(:,:)
        real(kind_phys), pointer :: dq3dt(:,:)
        real(kind_phys), pointer :: rain_cpl(:)
        real(kind_phys), pointer :: rainc_cpl(:)
        real(kind_phys), pointer :: snow_cpl(:)
        real(kind_phys), pointer :: pwat(:)
        logical, pointer :: do_sppt
        real(kind_phys), pointer :: dtdtr(:,:)
        real(kind_phys), pointer :: dtdtc(:,:)
        real(kind_phys), pointer :: drain_cpl(:)
        real(kind_phys), pointer :: dsnow_cpl(:)
        integer, pointer :: lsm
        integer, pointer :: lsm_ruc
        real(kind_phys), pointer :: raincprv(:)
        real(kind_phys), pointer :: rainncprv(:)
        real(kind_phys), pointer :: iceprv(:)
        real(kind_phys), pointer :: snowprv(:)
        real(kind_phys), pointer :: graupelprv(:)

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
        

        call ccpp_field_get(cdata, 'index_of_time_step', kdt, ierr=ierr, kind=ckind, index=398)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve index_of_time_step from CCPP data structure')
            return
        end if
        if (kind(kdt).ne.ckind) then
            call ccpp_error('Kind mismatch for variable index_of_time_step')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'array_dimension_of_random_number', nrcm, ierr=ierr, kind=ckind, index=64)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve array_dimension_of_random_number from CCPP data structure')
            return
        end if
        if (kind(nrcm).ne.ckind) then
            call ccpp_error('Kind mismatch for variable array_dimension_of_random_number')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'number_of_hydrometeors', ncld, ierr=ierr, kind=ckind, index=578)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve number_of_hydrometeors from CCPP data structure')
            return
        end if
        if (kind(ncld).ne.ckind) then
            call ccpp_error('Kind mismatch for variable number_of_hydrometeors')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'number_of_tracers_for_cloud_condensate', nncl, ierr=ierr, kind=ckind, index=586)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve number_of_tracers_for_cloud_condensate from CCPP data structure')
            return
        end if
        if (kind(nncl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable number_of_tracers_for_cloud_condensate')
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
        

        call ccpp_field_get(cdata, 'flag_for_precipitation_type_algorithm', cal_pre, ierr=ierr, kind=ckind, index=305)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_precipitation_type_algorithm from CCPP data structure')
            return
        end if
        if (kind(cal_pre).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_precipitation_type_algorithm')
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
        

        call ccpp_field_get(cdata, 'flag_for_chemistry_coupling', cplchm, ierr=ierr, kind=ckind, index=272)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_chemistry_coupling from CCPP data structure')
            return
        end if
        if (kind(cplchm).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_chemistry_coupling')
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
        

        call ccpp_field_get(cdata, 'lwe_thickness_of_explicit_precipitation_amount', rain1, ierr=ierr, dims=cdims, kind=ckind, index=480)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve lwe_thickness_of_explicit_precipitation_amount from CCPP data structure')
            return
        end if
        if (kind(rain1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable lwe_thickness_of_explicit_precipitation_amount')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'random_number_array', rann, ierr=ierr, dims=cdims, kind=ckind, index=620)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve random_number_array from CCPP data structure')
            return
        end if
        if (kind(rann).ne.ckind) then
            call ccpp_error('Kind mismatch for variable random_number_array')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'latitude', xlat, ierr=ierr, dims=cdims, kind=ckind, index=460)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve latitude from CCPP data structure')
            return
        end if
        if (kind(xlat).ne.ckind) then
            call ccpp_error('Kind mismatch for variable latitude')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'longitude', xlon, ierr=ierr, dims=cdims, kind=ckind, index=473)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve longitude from CCPP data structure')
            return
        end if
        if (kind(xlon).ne.ckind) then
            call ccpp_error('Kind mismatch for variable longitude')
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
        

        call ccpp_field_get(cdata, 'tracer_concentration_updated_by_physics', gq0, ierr=ierr, dims=cdims, kind=ckind, index=803)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve tracer_concentration_updated_by_physics from CCPP data structure')
            return
        end if
        if (kind(gq0).ne.ckind) then
            call ccpp_error('Kind mismatch for variable tracer_concentration_updated_by_physics')
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
        

        call ccpp_field_get(cdata, 'lwe_thickness_of_ice_amount_on_dynamics_timestep', ice, ierr=ierr, dims=cdims, kind=ckind, index=488)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve lwe_thickness_of_ice_amount_on_dynamics_timestep from CCPP data structure')
            return
        end if
        if (kind(ice).ne.ckind) then
            call ccpp_error('Kind mismatch for variable lwe_thickness_of_ice_amount_on_dynamics_timestep')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'lwe_thickness_of_snow_amount_on_dynamics_timestep', snow, ierr=ierr, dims=cdims, kind=ckind, index=495)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve lwe_thickness_of_snow_amount_on_dynamics_timestep from CCPP data structure')
            return
        end if
        if (kind(snow).ne.ckind) then
            call ccpp_error('Kind mismatch for variable lwe_thickness_of_snow_amount_on_dynamics_timestep')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'lwe_thickness_of_graupel_amount_on_dynamics_timestep', graupel, ierr=ierr, dims=cdims, kind=ckind, index=485)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve lwe_thickness_of_graupel_amount_on_dynamics_timestep from CCPP data structure')
            return
        end if
        if (kind(graupel).ne.ckind) then
            call ccpp_error('Kind mismatch for variable lwe_thickness_of_graupel_amount_on_dynamics_timestep')
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
        

        call ccpp_field_get(cdata, 'lwe_thickness_of_ice_amount', ice0, ierr=ierr, dims=cdims, kind=ckind, index=486)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve lwe_thickness_of_ice_amount from CCPP data structure')
            return
        end if
        if (kind(ice0).ne.ckind) then
            call ccpp_error('Kind mismatch for variable lwe_thickness_of_ice_amount')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'lwe_thickness_of_snow_amount', snow0, ierr=ierr, dims=cdims, kind=ckind, index=492)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve lwe_thickness_of_snow_amount from CCPP data structure')
            return
        end if
        if (kind(snow0).ne.ckind) then
            call ccpp_error('Kind mismatch for variable lwe_thickness_of_snow_amount')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'lwe_thickness_of_graupel_amount', graupel0, ierr=ierr, dims=cdims, kind=ckind, index=483)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve lwe_thickness_of_graupel_amount from CCPP data structure')
            return
        end if
        if (kind(graupel0).ne.ckind) then
            call ccpp_error('Kind mismatch for variable lwe_thickness_of_graupel_amount')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'air_pressure_difference_between_midlayers', del, ierr=ierr, dims=cdims, kind=ckind, index=49)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve air_pressure_difference_between_midlayers from CCPP data structure')
            return
        end if
        if (kind(del).ne.ckind) then
            call ccpp_error('Kind mismatch for variable air_pressure_difference_between_midlayers')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'lwe_thickness_of_precipitation_amount_on_dynamics_timestep', rain, ierr=ierr, dims=cdims, kind=ckind, index=490)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve lwe_thickness_of_precipitation_amount_on_dynamics_timestep from CCPP data structure')
            return
        end if
        if (kind(rain).ne.ckind) then
            call ccpp_error('Kind mismatch for variable lwe_thickness_of_precipitation_amount_on_dynamics_timestep')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'dominant_rain_type', domr_diag, ierr=ierr, dims=cdims, kind=ckind, index=230)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve dominant_rain_type from CCPP data structure')
            return
        end if
        if (kind(domr_diag).ne.ckind) then
            call ccpp_error('Kind mismatch for variable dominant_rain_type')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'dominant_freezing_rain_type', domzr_diag, ierr=ierr, dims=cdims, kind=ckind, index=229)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve dominant_freezing_rain_type from CCPP data structure')
            return
        end if
        if (kind(domzr_diag).ne.ckind) then
            call ccpp_error('Kind mismatch for variable dominant_freezing_rain_type')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'dominant_sleet_type', domip_diag, ierr=ierr, dims=cdims, kind=ckind, index=231)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve dominant_sleet_type from CCPP data structure')
            return
        end if
        if (kind(domip_diag).ne.ckind) then
            call ccpp_error('Kind mismatch for variable dominant_sleet_type')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'dominant_snow_type', doms_diag, ierr=ierr, dims=cdims, kind=ckind, index=232)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve dominant_snow_type from CCPP data structure')
            return
        end if
        if (kind(doms_diag).ne.ckind) then
            call ccpp_error('Kind mismatch for variable dominant_snow_type')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'nonnegative_lwe_thickness_of_precipitation_amount_on_dynamics_timestep', tprcp, ierr=ierr, dims=cdims, kind=ckind, index=567)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve nonnegative_lwe_thickness_of_precipitation_amount_on_dynamics_timestep from CCPP data structure')
            return
        end if
        if (kind(tprcp).ne.ckind) then
            call ccpp_error('Kind mismatch for variable nonnegative_lwe_thickness_of_precipitation_amount_on_dynamics_timestep')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'flag_for_precipitation_type', srflag, ierr=ierr, dims=cdims, kind=ckind, index=304)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_precipitation_type from CCPP data structure')
            return
        end if
        if (kind(srflag).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_precipitation_type')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'accumulated_lwe_thickness_of_precipitation_amount', totprcp, ierr=ierr, dims=cdims, kind=ckind, index=26)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve accumulated_lwe_thickness_of_precipitation_amount from CCPP data structure')
            return
        end if
        if (kind(totprcp).ne.ckind) then
            call ccpp_error('Kind mismatch for variable accumulated_lwe_thickness_of_precipitation_amount')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'accumulated_lwe_thickness_of_ice_amount', totice, ierr=ierr, dims=cdims, kind=ckind, index=24)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve accumulated_lwe_thickness_of_ice_amount from CCPP data structure')
            return
        end if
        if (kind(totice).ne.ckind) then
            call ccpp_error('Kind mismatch for variable accumulated_lwe_thickness_of_ice_amount')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'accumulated_lwe_thickness_of_snow_amount', totsnw, ierr=ierr, dims=cdims, kind=ckind, index=28)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve accumulated_lwe_thickness_of_snow_amount from CCPP data structure')
            return
        end if
        if (kind(totsnw).ne.ckind) then
            call ccpp_error('Kind mismatch for variable accumulated_lwe_thickness_of_snow_amount')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'accumulated_lwe_thickness_of_graupel_amount', totgrp, ierr=ierr, dims=cdims, kind=ckind, index=22)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve accumulated_lwe_thickness_of_graupel_amount from CCPP data structure')
            return
        end if
        if (kind(totgrp).ne.ckind) then
            call ccpp_error('Kind mismatch for variable accumulated_lwe_thickness_of_graupel_amount')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'accumulated_lwe_thickness_of_precipitation_amount_in_bucket', totprcpb, ierr=ierr, dims=cdims, kind=ckind, index=27)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve accumulated_lwe_thickness_of_precipitation_amount_in_bucket from CCPP data structure')
            return
        end if
        if (kind(totprcpb).ne.ckind) then
            call ccpp_error('Kind mismatch for variable accumulated_lwe_thickness_of_precipitation_amount_in_bucket')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'accumulated_lwe_thickness_of_ice_amount_in_bucket', toticeb, ierr=ierr, dims=cdims, kind=ckind, index=25)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve accumulated_lwe_thickness_of_ice_amount_in_bucket from CCPP data structure')
            return
        end if
        if (kind(toticeb).ne.ckind) then
            call ccpp_error('Kind mismatch for variable accumulated_lwe_thickness_of_ice_amount_in_bucket')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'accumulated_lwe_thickness_of_snow_amount_in_bucket', totsnwb, ierr=ierr, dims=cdims, kind=ckind, index=29)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve accumulated_lwe_thickness_of_snow_amount_in_bucket from CCPP data structure')
            return
        end if
        if (kind(totsnwb).ne.ckind) then
            call ccpp_error('Kind mismatch for variable accumulated_lwe_thickness_of_snow_amount_in_bucket')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'accumulated_lwe_thickness_of_graupel_amount_in_bucket', totgrpb, ierr=ierr, dims=cdims, kind=ckind, index=23)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve accumulated_lwe_thickness_of_graupel_amount_in_bucket from CCPP data structure')
            return
        end if
        if (kind(totgrpb).ne.ckind) then
            call ccpp_error('Kind mismatch for variable accumulated_lwe_thickness_of_graupel_amount_in_bucket')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cumulative_change_in_temperature_due_to_microphysics', dt3dt, ierr=ierr, dims=cdims, kind=ckind, index=158)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cumulative_change_in_temperature_due_to_microphysics from CCPP data structure')
            return
        end if
        if (kind(dt3dt).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cumulative_change_in_temperature_due_to_microphysics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cumulative_change_in_water_vapor_specific_humidity_due_to_microphysics', dq3dt, ierr=ierr, dims=cdims, kind=ckind, index=163)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cumulative_change_in_water_vapor_specific_humidity_due_to_microphysics from CCPP data structure')
            return
        end if
        if (kind(dq3dt).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cumulative_change_in_water_vapor_specific_humidity_due_to_microphysics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'lwe_thickness_of_precipitation_amount_for_coupling', rain_cpl, ierr=ierr, dims=cdims, kind=ckind, index=489)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve lwe_thickness_of_precipitation_amount_for_coupling from CCPP data structure')
            return
        end if
        if (kind(rain_cpl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable lwe_thickness_of_precipitation_amount_for_coupling')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'lwe_thickness_of_convective_precipitation_amount_for_coupling', rainc_cpl, ierr=ierr, dims=cdims, kind=ckind, index=476)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve lwe_thickness_of_convective_precipitation_amount_for_coupling from CCPP data structure')
            return
        end if
        if (kind(rainc_cpl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable lwe_thickness_of_convective_precipitation_amount_for_coupling')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'lwe_thickness_of_snow_amount_for_coupling', snow_cpl, ierr=ierr, dims=cdims, kind=ckind, index=493)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve lwe_thickness_of_snow_amount_for_coupling from CCPP data structure')
            return
        end if
        if (kind(snow_cpl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable lwe_thickness_of_snow_amount_for_coupling')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'column_precipitable_water', pwat, ierr=ierr, dims=cdims, kind=ckind, index=120)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve column_precipitable_water from CCPP data structure')
            return
        end if
        if (kind(pwat).ne.ckind) then
            call ccpp_error('Kind mismatch for variable column_precipitable_water')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'flag_for_stochastic_surface_physics_perturbations', do_sppt, ierr=ierr, kind=ckind, index=319)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_stochastic_surface_physics_perturbations from CCPP data structure')
            return
        end if
        if (kind(do_sppt).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_stochastic_surface_physics_perturbations')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'tendency_of_air_temperature_due_to_radiative_heating_on_physics_time_step', dtdtr, ierr=ierr, dims=cdims, kind=ckind, index=760)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve tendency_of_air_temperature_due_to_radiative_heating_on_physics_time_step from CCPP data structure')
            return
        end if
        if (kind(dtdtr).ne.ckind) then
            call ccpp_error('Kind mismatch for variable tendency_of_air_temperature_due_to_radiative_heating_on_physics_time_step')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'tendency_of_air_temperature_due_to_radiative_heating_assuming_clear_sky', dtdtc, ierr=ierr, dims=cdims, kind=ckind, index=759)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve tendency_of_air_temperature_due_to_radiative_heating_assuming_clear_sky from CCPP data structure')
            return
        end if
        if (kind(dtdtc).ne.ckind) then
            call ccpp_error('Kind mismatch for variable tendency_of_air_temperature_due_to_radiative_heating_assuming_clear_sky')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'tendency_of_lwe_thickness_of_precipitation_amount_for_coupling', drain_cpl, ierr=ierr, dims=cdims, kind=ckind, index=772)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve tendency_of_lwe_thickness_of_precipitation_amount_for_coupling from CCPP data structure')
            return
        end if
        if (kind(drain_cpl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable tendency_of_lwe_thickness_of_precipitation_amount_for_coupling')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'tendency_of_lwe_thickness_of_snow_amount_for_coupling', dsnow_cpl, ierr=ierr, dims=cdims, kind=ckind, index=773)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve tendency_of_lwe_thickness_of_snow_amount_for_coupling from CCPP data structure')
            return
        end if
        if (kind(dsnow_cpl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable tendency_of_lwe_thickness_of_snow_amount_for_coupling')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'flag_for_land_surface_scheme', lsm, ierr=ierr, kind=ckind, index=288)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_land_surface_scheme from CCPP data structure')
            return
        end if
        if (kind(lsm).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_land_surface_scheme')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'flag_for_ruc_land_surface_scheme', lsm_ruc, ierr=ierr, kind=ckind, index=310)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_ruc_land_surface_scheme from CCPP data structure')
            return
        end if
        if (kind(lsm_ruc).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_ruc_land_surface_scheme')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'lwe_thickness_of_explicit_rainfall_amount_from_previous_timestep', raincprv, ierr=ierr, dims=cdims, kind=ckind, index=482)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve lwe_thickness_of_explicit_rainfall_amount_from_previous_timestep from CCPP data structure')
            return
        end if
        if (kind(raincprv).ne.ckind) then
            call ccpp_error('Kind mismatch for variable lwe_thickness_of_explicit_rainfall_amount_from_previous_timestep')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'lwe_thickness_of_convective_precipitation_amount_from_previous_timestep', rainncprv, ierr=ierr, dims=cdims, kind=ckind, index=477)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve lwe_thickness_of_convective_precipitation_amount_from_previous_timestep from CCPP data structure')
            return
        end if
        if (kind(rainncprv).ne.ckind) then
            call ccpp_error('Kind mismatch for variable lwe_thickness_of_convective_precipitation_amount_from_previous_timestep')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'lwe_thickness_of_ice_amount_from_previous_timestep', iceprv, ierr=ierr, dims=cdims, kind=ckind, index=487)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve lwe_thickness_of_ice_amount_from_previous_timestep from CCPP data structure')
            return
        end if
        if (kind(iceprv).ne.ckind) then
            call ccpp_error('Kind mismatch for variable lwe_thickness_of_ice_amount_from_previous_timestep')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'lwe_thickness_of_snow_amount_from_previous_timestep', snowprv, ierr=ierr, dims=cdims, kind=ckind, index=494)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve lwe_thickness_of_snow_amount_from_previous_timestep from CCPP data structure')
            return
        end if
        if (kind(snowprv).ne.ckind) then
            call ccpp_error('Kind mismatch for variable lwe_thickness_of_snow_amount_from_previous_timestep')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'lwe_thickness_of_graupel_amount_from_previous_timestep', graupelprv, ierr=ierr, dims=cdims, kind=ckind, index=484)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve lwe_thickness_of_graupel_amount_from_previous_timestep from CCPP data structure')
            return
        end if
        if (kind(graupelprv).ne.ckind) then
            call ccpp_error('Kind mismatch for variable lwe_thickness_of_graupel_amount_from_previous_timestep')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call GFS_MP_generic_post_run(im=im,ix=ix,levs=levs,kdt=kdt,nrcm=nrcm,ncld=ncld,nncl=nncl,ntcw=ntcw,ntrac=ntrac, &
                  imp_physics=imp_physics,imp_physics_gfdl=imp_physics_gfdl,imp_physics_thompson=imp_physics_thompson, &
                  cal_pre=cal_pre,lssav=lssav,ldiag3d=ldiag3d,cplflx=cplflx,cplchm=cplchm, &
                  con_g=con_g,dtf=dtf,frain=frain,rainc=rainc,rain1=rain1,rann=rann,xlat=xlat, &
                  xlon=xlon,gt0=gt0,gq0=gq0,prsl=prsl,prsi=prsi,phii=phii,tsfc=tsfc,ice=ice, &
                  snow=snow,graupel=graupel,save_t=save_t,save_qv=save_qv,ice0=ice0,snow0=snow0, &
                  graupel0=graupel0,del=del,rain=rain,domr_diag=domr_diag,domzr_diag=domzr_diag, &
                  domip_diag=domip_diag,doms_diag=doms_diag,tprcp=tprcp,srflag=srflag,totprcp=totprcp, &
                  totice=totice,totsnw=totsnw,totgrp=totgrp,totprcpb=totprcpb,toticeb=toticeb, &
                  totsnwb=totsnwb,totgrpb=totgrpb,dt3dt=dt3dt,dq3dt=dq3dt,rain_cpl=rain_cpl, &
                  rainc_cpl=rainc_cpl,snow_cpl=snow_cpl,pwat=pwat,do_sppt=do_sppt,dtdtr=dtdtr, &
                  dtdtc=dtdtc,drain_cpl=drain_cpl,dsnow_cpl=dsnow_cpl,lsm=lsm,lsm_ruc=lsm_ruc, &
                  raincprv=raincprv,rainncprv=rainncprv,iceprv=iceprv,snowprv=snowprv,graupelprv=graupelprv, &
                  errmsg=cdata%errmsg,errflg=cdata%errflg)
        ierr=cdata%errflg

    end function GFS_MP_generic_post_run_cap

    function GFS_MP_generic_post_finalize_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call GFS_MP_generic_post_finalize()
        

    end function GFS_MP_generic_post_finalize_cap
end module GFS_MP_generic_post_cap

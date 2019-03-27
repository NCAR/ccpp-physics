
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
!! @brief Auto-generated cap module for the lsm_ruc scheme
!!
!
module lsm_ruc_cap

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_f_pointer, c_ptr, c_int32_t
    use            :: ccpp_types,                                      &
                      only: ccpp_t, CCPP_GENERIC_KIND
    use            :: ccpp_fields,                                     &
                      only: ccpp_field_get
    use            :: ccpp_errors,                                     &
                      only: ccpp_error, ccpp_debug
    use            :: lsm_ruc, &
                      only: lsm_ruc_finalize,lsm_ruc_run,lsm_ruc_init
    ! Other modules required, e.g. type definitions
    use machine, only: kind_phys

    implicit none

    private
    public :: lsm_ruc_finalize_cap,lsm_ruc_run_cap,lsm_ruc_init_cap

    contains


    function lsm_ruc_finalize_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call lsm_ruc_finalize(errmsg=cdata%errmsg,errflg=cdata%errflg)
        ierr=cdata%errflg

    end function lsm_ruc_finalize_cap

    function lsm_ruc_run_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind
        real(kind_phys), pointer :: delt
        integer, pointer :: me
        integer, pointer :: kdt
        integer, pointer :: im
        integer, pointer :: nlev
        integer, pointer :: lsm_ruc
        integer, pointer :: lsm
        logical, pointer :: do_mynnsfclay
        integer, pointer :: lsoil_ruc
        integer, pointer :: lsoil
        real(kind_phys), pointer :: zs(:)
        integer, pointer :: islmsk(:)
        real(kind_phys), pointer :: con_cp
        real(kind_phys), pointer :: con_g
        real(kind_phys), pointer :: con_pi
        real(kind_phys), pointer :: con_rd
        real(kind_phys), pointer :: con_rv
        real(kind_phys), pointer :: con_hvap
        real(kind_phys), pointer :: con_fvirt
        real(kind_phys), pointer :: rainnc(:)
        real(kind_phys), pointer :: rainc(:)
        real(kind_phys), pointer :: ice(:)
        real(kind_phys), pointer :: snow(:)
        real(kind_phys), pointer :: graupel(:)
        real(kind_phys), pointer :: srflag(:)
        real(kind_phys), pointer :: sncovr1(:)
        real(kind_phys), pointer :: snowc(:)
        real(kind_phys), pointer :: weasd(:)
        real(kind_phys), pointer :: snwdph(:)
        real(kind_phys), pointer :: sr(:)
        real(kind_phys), pointer :: rhosnf(:)
        real(kind_phys), pointer :: zf(:)
        real(kind_phys), pointer :: u1(:)
        real(kind_phys), pointer :: v1(:)
        real(kind_phys), pointer :: prsl1(:)
        real(kind_phys), pointer :: ddvel(:)
        real(kind_phys), pointer :: t1(:)
        real(kind_phys), pointer :: q1(:)
        real(kind_phys), pointer :: qc(:)
        real(kind_phys), pointer :: dlwflx(:)
        real(kind_phys), pointer :: dswsfc(:)
        real(kind_phys), pointer :: snet(:)
        real(kind_phys), pointer :: sfcemis(:)
        real(kind_phys), pointer :: wspd(:)
        real(kind_phys), pointer :: cm(:)
        real(kind_phys), pointer :: ch(:)
        real(kind_phys), pointer :: chh(:)
        real(kind_phys), pointer :: cmm(:)
        real(kind_phys), pointer :: wet1(:)
        real(kind_phys), pointer :: canopy(:)
        real(kind_phys), pointer :: sigmaf(:)
        real(kind_phys), pointer :: sfalb(:)
        real(kind_phys), pointer :: alvwf(:)
        real(kind_phys), pointer :: alnwf(:)
        real(kind_phys), pointer :: snoalb(:)
        real(kind_phys), pointer :: zorl(:)
        real(kind_phys), pointer :: qsurf(:)
        real(kind_phys), pointer :: sfcqc(:)
        real(kind_phys), pointer :: sfcqv(:)
        real(kind_phys), pointer :: sfcdew(:)
        real(kind_phys), pointer :: tg3(:)
        real(kind_phys), pointer :: smc(:,:)
        real(kind_phys), pointer :: slc(:,:)
        real(kind_phys), pointer :: stc(:,:)
        real(kind_phys), pointer :: smcwlt2(:)
        real(kind_phys), pointer :: smcref2(:)
        integer, pointer :: vegtype(:)
        integer, pointer :: soiltyp(:)
        integer, pointer :: isot
        integer, pointer :: ivegsrc
        real(kind_phys), pointer :: fice(:)
        real(kind_phys), pointer :: keepfr(:,:)
        real(kind_phys), pointer :: smois(:,:)
        real(kind_phys), pointer :: sh2o(:,:)
        real(kind_phys), pointer :: smfrkeep(:,:)
        real(kind_phys), pointer :: tslb(:,:)
        real(kind_phys), pointer :: stm(:)
        real(kind_phys), pointer :: tskin(:)
        real(kind_phys), pointer :: tsurf(:)
        real(kind_phys), pointer :: tice(:)
        real(kind_phys), pointer :: tsnow(:)
        real(kind_phys), pointer :: snowfallac(:)
        real(kind_phys), pointer :: acsnow(:)
        real(kind_phys), pointer :: evap(:)
        real(kind_phys), pointer :: hflx(:)
        real(kind_phys), pointer :: evbs(:)
        real(kind_phys), pointer :: evcw(:)
        real(kind_phys), pointer :: sbsno(:)
        real(kind_phys), pointer :: trans(:)
        real(kind_phys), pointer :: runof(:)
        real(kind_phys), pointer :: drain(:)
        real(kind_phys), pointer :: runoff(:)
        real(kind_phys), pointer :: srunoff(:)
        real(kind_phys), pointer :: gflux(:)
        real(kind_phys), pointer :: shdmin(:)
        real(kind_phys), pointer :: shdmax(:)
        logical, pointer :: flag_iter(:)
        logical, pointer :: flag_guess(:)
        logical, pointer :: flag_init
        logical, pointer :: flag_restart

        ierr = 0

        call c_f_pointer(ptr, cdata)


        call ccpp_field_get(cdata, 'time_step_for_dynamics', delt, ierr=ierr, kind=ckind, index=792)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve time_step_for_dynamics from CCPP data structure')
            return
        end if
        if (kind(delt).ne.ckind) then
            call ccpp_error('Kind mismatch for variable time_step_for_dynamics')
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
        

        call ccpp_field_get(cdata, 'vertical_dimension', nlev, ierr=ierr, kind=ckind, index=817)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve vertical_dimension from CCPP data structure')
            return
        end if
        if (kind(nlev).ne.ckind) then
            call ccpp_error('Kind mismatch for variable vertical_dimension')
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
        

        call ccpp_field_get(cdata, 'soil_vertical_dimension_for_land_surface_model', lsoil_ruc, ierr=ierr, kind=ckind, index=662)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve soil_vertical_dimension_for_land_surface_model from CCPP data structure')
            return
        end if
        if (kind(lsoil_ruc).ne.ckind) then
            call ccpp_error('Kind mismatch for variable soil_vertical_dimension_for_land_surface_model')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'soil_vertical_dimension', lsoil, ierr=ierr, kind=ckind, index=661)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve soil_vertical_dimension from CCPP data structure')
            return
        end if
        if (kind(lsoil).ne.ckind) then
            call ccpp_error('Kind mismatch for variable soil_vertical_dimension')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'depth_of_soil_levels_for_land_surface_model', zs, ierr=ierr, dims=cdims, kind=ckind, index=212)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve depth_of_soil_levels_for_land_surface_model from CCPP data structure')
            return
        end if
        if (kind(zs).ne.ckind) then
            call ccpp_error('Kind mismatch for variable depth_of_soil_levels_for_land_surface_model')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'sea_land_ice_mask', islmsk, ierr=ierr, dims=cdims, kind=ckind, index=631)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve sea_land_ice_mask from CCPP data structure')
            return
        end if
        if (kind(islmsk).ne.ckind) then
            call ccpp_error('Kind mismatch for variable sea_land_ice_mask')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

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
        

        call ccpp_field_get(cdata, 'lwe_thickness_of_explicit_rainfall_amount_from_previous_timestep', rainnc, ierr=ierr, dims=cdims, kind=ckind, index=482)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve lwe_thickness_of_explicit_rainfall_amount_from_previous_timestep from CCPP data structure')
            return
        end if
        if (kind(rainnc).ne.ckind) then
            call ccpp_error('Kind mismatch for variable lwe_thickness_of_explicit_rainfall_amount_from_previous_timestep')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'lwe_thickness_of_convective_precipitation_amount_from_previous_timestep', rainc, ierr=ierr, dims=cdims, kind=ckind, index=477)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve lwe_thickness_of_convective_precipitation_amount_from_previous_timestep from CCPP data structure')
            return
        end if
        if (kind(rainc).ne.ckind) then
            call ccpp_error('Kind mismatch for variable lwe_thickness_of_convective_precipitation_amount_from_previous_timestep')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'lwe_thickness_of_ice_amount_from_previous_timestep', ice, ierr=ierr, dims=cdims, kind=ckind, index=487)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve lwe_thickness_of_ice_amount_from_previous_timestep from CCPP data structure')
            return
        end if
        if (kind(ice).ne.ckind) then
            call ccpp_error('Kind mismatch for variable lwe_thickness_of_ice_amount_from_previous_timestep')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'lwe_thickness_of_snow_amount_from_previous_timestep', snow, ierr=ierr, dims=cdims, kind=ckind, index=494)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve lwe_thickness_of_snow_amount_from_previous_timestep from CCPP data structure')
            return
        end if
        if (kind(snow).ne.ckind) then
            call ccpp_error('Kind mismatch for variable lwe_thickness_of_snow_amount_from_previous_timestep')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'lwe_thickness_of_graupel_amount_from_previous_timestep', graupel, ierr=ierr, dims=cdims, kind=ckind, index=484)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve lwe_thickness_of_graupel_amount_from_previous_timestep from CCPP data structure')
            return
        end if
        if (kind(graupel).ne.ckind) then
            call ccpp_error('Kind mismatch for variable lwe_thickness_of_graupel_amount_from_previous_timestep')
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
        

        call ccpp_field_get(cdata, 'surface_snow_area_fraction_for_diagnostics', sncovr1, ierr=ierr, dims=cdims, kind=ckind, index=727)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_snow_area_fraction_for_diagnostics from CCPP data structure')
            return
        end if
        if (kind(sncovr1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_snow_area_fraction_for_diagnostics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_snow_area_fraction', snowc, ierr=ierr, dims=cdims, kind=ckind, index=726)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_snow_area_fraction from CCPP data structure')
            return
        end if
        if (kind(snowc).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_snow_area_fraction')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'water_equivalent_accumulated_snow_depth', weasd, ierr=ierr, dims=cdims, kind=ckind, index=848)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve water_equivalent_accumulated_snow_depth from CCPP data structure')
            return
        end if
        if (kind(weasd).ne.ckind) then
            call ccpp_error('Kind mismatch for variable water_equivalent_accumulated_snow_depth')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_snow_thickness_water_equivalent', snwdph, ierr=ierr, dims=cdims, kind=ckind, index=729)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_snow_thickness_water_equivalent from CCPP data structure')
            return
        end if
        if (kind(snwdph).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_snow_thickness_water_equivalent')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'ratio_of_snowfall_to_rainfall', sr, ierr=ierr, dims=cdims, kind=ckind, index=624)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve ratio_of_snowfall_to_rainfall from CCPP data structure')
            return
        end if
        if (kind(sr).ne.ckind) then
            call ccpp_error('Kind mismatch for variable ratio_of_snowfall_to_rainfall')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'density_of_frozen_precipitation', rhosnf, ierr=ierr, dims=cdims, kind=ckind, index=211)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve density_of_frozen_precipitation from CCPP data structure')
            return
        end if
        if (kind(rhosnf).ne.ckind) then
            call ccpp_error('Kind mismatch for variable density_of_frozen_precipitation')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'height_above_ground_at_lowest_model_layer', zf, ierr=ierr, dims=cdims, kind=ckind, index=360)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve height_above_ground_at_lowest_model_layer from CCPP data structure')
            return
        end if
        if (kind(zf).ne.ckind) then
            call ccpp_error('Kind mismatch for variable height_above_ground_at_lowest_model_layer')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'x_wind_at_lowest_model_layer', u1, ierr=ierr, dims=cdims, kind=ckind, index=873)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve x_wind_at_lowest_model_layer from CCPP data structure')
            return
        end if
        if (kind(u1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable x_wind_at_lowest_model_layer')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'y_wind_at_lowest_model_layer', v1, ierr=ierr, dims=cdims, kind=ckind, index=880)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve y_wind_at_lowest_model_layer from CCPP data structure')
            return
        end if
        if (kind(v1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable y_wind_at_lowest_model_layer')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'air_pressure_at_lowest_model_layer', prsl1, ierr=ierr, dims=cdims, kind=ckind, index=48)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve air_pressure_at_lowest_model_layer from CCPP data structure')
            return
        end if
        if (kind(prsl1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable air_pressure_at_lowest_model_layer')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_wind_enhancement_due_to_convection', ddvel, ierr=ierr, dims=cdims, kind=ckind, index=744)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_wind_enhancement_due_to_convection from CCPP data structure')
            return
        end if
        if (kind(ddvel).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_wind_enhancement_due_to_convection')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'air_temperature_at_lowest_model_layer', t1, ierr=ierr, dims=cdims, kind=ckind, index=53)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve air_temperature_at_lowest_model_layer from CCPP data structure')
            return
        end if
        if (kind(t1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable air_temperature_at_lowest_model_layer')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'water_vapor_specific_humidity_at_lowest_model_layer', q1, ierr=ierr, dims=cdims, kind=ckind, index=854)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve water_vapor_specific_humidity_at_lowest_model_layer from CCPP data structure')
            return
        end if
        if (kind(q1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable water_vapor_specific_humidity_at_lowest_model_layer')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cloud_condensed_water_mixing_ratio_at_lowest_model_layer', qc, ierr=ierr, dims=cdims, kind=ckind, index=91)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cloud_condensed_water_mixing_ratio_at_lowest_model_layer from CCPP data structure')
            return
        end if
        if (kind(qc).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cloud_condensed_water_mixing_ratio_at_lowest_model_layer')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_downwelling_longwave_flux', dlwflx, ierr=ierr, dims=cdims, kind=ckind, index=697)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_downwelling_longwave_flux from CCPP data structure')
            return
        end if
        if (kind(dlwflx).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_downwelling_longwave_flux')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_downwelling_shortwave_flux', dswsfc, ierr=ierr, dims=cdims, kind=ckind, index=700)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_downwelling_shortwave_flux from CCPP data structure')
            return
        end if
        if (kind(dswsfc).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_downwelling_shortwave_flux')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_net_downwelling_shortwave_flux', snet, ierr=ierr, dims=cdims, kind=ckind, index=716)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_net_downwelling_shortwave_flux from CCPP data structure')
            return
        end if
        if (kind(snet).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_net_downwelling_shortwave_flux')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_longwave_emissivity', sfcemis, ierr=ierr, dims=cdims, kind=ckind, index=714)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_longwave_emissivity from CCPP data structure')
            return
        end if
        if (kind(sfcemis).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_longwave_emissivity')
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
        

        call ccpp_field_get(cdata, 'surface_drag_coefficient_for_momentum_in_air', cm, ierr=ierr, dims=cdims, kind=ckind, index=703)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_drag_coefficient_for_momentum_in_air from CCPP data structure')
            return
        end if
        if (kind(cm).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_drag_coefficient_for_momentum_in_air')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_drag_coefficient_for_heat_and_moisture_in_air', ch, ierr=ierr, dims=cdims, kind=ckind, index=702)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_drag_coefficient_for_heat_and_moisture_in_air from CCPP data structure')
            return
        end if
        if (kind(ch).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_drag_coefficient_for_heat_and_moisture_in_air')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_drag_mass_flux_for_heat_and_moisture_in_air', chh, ierr=ierr, dims=cdims, kind=ckind, index=704)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_drag_mass_flux_for_heat_and_moisture_in_air from CCPP data structure')
            return
        end if
        if (kind(chh).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_drag_mass_flux_for_heat_and_moisture_in_air')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_drag_wind_speed_for_momentum_in_air', cmm, ierr=ierr, dims=cdims, kind=ckind, index=705)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_drag_wind_speed_for_momentum_in_air from CCPP data structure')
            return
        end if
        if (kind(cmm).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_drag_wind_speed_for_momentum_in_air')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'normalized_soil_wetness', wet1, ierr=ierr, dims=cdims, kind=ckind, index=568)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve normalized_soil_wetness from CCPP data structure')
            return
        end if
        if (kind(wet1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable normalized_soil_wetness')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'canopy_water_amount', canopy, ierr=ierr, dims=cdims, kind=ckind, index=80)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve canopy_water_amount from CCPP data structure')
            return
        end if
        if (kind(canopy).ne.ckind) then
            call ccpp_error('Kind mismatch for variable canopy_water_amount')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'vegetation_area_fraction', sigmaf, ierr=ierr, dims=cdims, kind=ckind, index=813)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve vegetation_area_fraction from CCPP data structure')
            return
        end if
        if (kind(sigmaf).ne.ckind) then
            call ccpp_error('Kind mismatch for variable vegetation_area_fraction')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_diffused_shortwave_albedo', sfalb, ierr=ierr, dims=cdims, kind=ckind, index=688)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_diffused_shortwave_albedo from CCPP data structure')
            return
        end if
        if (kind(sfalb).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_diffused_shortwave_albedo')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'mean_vis_albedo_with_weak_cosz_dependency', alvwf, ierr=ierr, dims=cdims, kind=ckind, index=520)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve mean_vis_albedo_with_weak_cosz_dependency from CCPP data structure')
            return
        end if
        if (kind(alvwf).ne.ckind) then
            call ccpp_error('Kind mismatch for variable mean_vis_albedo_with_weak_cosz_dependency')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'mean_nir_albedo_with_weak_cosz_dependency', alnwf, ierr=ierr, dims=cdims, kind=ckind, index=519)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve mean_nir_albedo_with_weak_cosz_dependency from CCPP data structure')
            return
        end if
        if (kind(alnwf).ne.ckind) then
            call ccpp_error('Kind mismatch for variable mean_nir_albedo_with_weak_cosz_dependency')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'upper_bound_on_max_albedo_over_deep_snow', snoalb, ierr=ierr, dims=cdims, kind=ckind, index=811)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve upper_bound_on_max_albedo_over_deep_snow from CCPP data structure')
            return
        end if
        if (kind(snoalb).ne.ckind) then
            call ccpp_error('Kind mismatch for variable upper_bound_on_max_albedo_over_deep_snow')
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
        

        call ccpp_field_get(cdata, 'surface_specific_humidity', qsurf, ierr=ierr, dims=cdims, kind=ckind, index=730)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_specific_humidity from CCPP data structure')
            return
        end if
        if (kind(qsurf).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_specific_humidity')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cloud_condensed_water_mixing_ratio_at_surface', sfcqc, ierr=ierr, dims=cdims, kind=ckind, index=92)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cloud_condensed_water_mixing_ratio_at_surface from CCPP data structure')
            return
        end if
        if (kind(sfcqc).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cloud_condensed_water_mixing_ratio_at_surface')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'water_vapor_mixing_ratio_at_surface', sfcqv, ierr=ierr, dims=cdims, kind=ckind, index=851)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve water_vapor_mixing_ratio_at_surface from CCPP data structure')
            return
        end if
        if (kind(sfcqv).ne.ckind) then
            call ccpp_error('Kind mismatch for variable water_vapor_mixing_ratio_at_surface')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_condensation_mass', sfcdew, ierr=ierr, dims=cdims, kind=ckind, index=687)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_condensation_mass from CCPP data structure')
            return
        end if
        if (kind(sfcdew).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_condensation_mass')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'deep_soil_temperature', tg3, ierr=ierr, dims=cdims, kind=ckind, index=210)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve deep_soil_temperature from CCPP data structure')
            return
        end if
        if (kind(tg3).ne.ckind) then
            call ccpp_error('Kind mismatch for variable deep_soil_temperature')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'volume_fraction_of_soil_moisture', smc, ierr=ierr, dims=cdims, kind=ckind, index=834)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve volume_fraction_of_soil_moisture from CCPP data structure')
            return
        end if
        if (kind(smc).ne.ckind) then
            call ccpp_error('Kind mismatch for variable volume_fraction_of_soil_moisture')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'volume_fraction_of_unfrozen_soil_moisture', slc, ierr=ierr, dims=cdims, kind=ckind, index=836)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve volume_fraction_of_unfrozen_soil_moisture from CCPP data structure')
            return
        end if
        if (kind(slc).ne.ckind) then
            call ccpp_error('Kind mismatch for variable volume_fraction_of_unfrozen_soil_moisture')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'soil_temperature', stc, ierr=ierr, dims=cdims, kind=ckind, index=655)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve soil_temperature from CCPP data structure')
            return
        end if
        if (kind(stc).ne.ckind) then
            call ccpp_error('Kind mismatch for variable soil_temperature')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'volume_fraction_of_condensed_water_in_soil_at_wilting_point', smcwlt2, ierr=ierr, dims=cdims, kind=ckind, index=832)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve volume_fraction_of_condensed_water_in_soil_at_wilting_point from CCPP data structure')
            return
        end if
        if (kind(smcwlt2).ne.ckind) then
            call ccpp_error('Kind mismatch for variable volume_fraction_of_condensed_water_in_soil_at_wilting_point')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'threshold_volume_fraction_of_condensed_water_in_soil', smcref2, ierr=ierr, dims=cdims, kind=ckind, index=788)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve threshold_volume_fraction_of_condensed_water_in_soil from CCPP data structure')
            return
        end if
        if (kind(smcref2).ne.ckind) then
            call ccpp_error('Kind mismatch for variable threshold_volume_fraction_of_condensed_water_in_soil')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'vegetation_type_classification', vegtype, ierr=ierr, dims=cdims, kind=ckind, index=814)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve vegetation_type_classification from CCPP data structure')
            return
        end if
        if (kind(vegtype).ne.ckind) then
            call ccpp_error('Kind mismatch for variable vegetation_type_classification')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'soil_type_classification', soiltyp, ierr=ierr, dims=cdims, kind=ckind, index=657)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve soil_type_classification from CCPP data structure')
            return
        end if
        if (kind(soiltyp).ne.ckind) then
            call ccpp_error('Kind mismatch for variable soil_type_classification')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'soil_type_dataset_choice', isot, ierr=ierr, kind=ckind, index=659)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve soil_type_dataset_choice from CCPP data structure')
            return
        end if
        if (kind(isot).ne.ckind) then
            call ccpp_error('Kind mismatch for variable soil_type_dataset_choice')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'vegetation_type_dataset_choice', ivegsrc, ierr=ierr, kind=ckind, index=816)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve vegetation_type_dataset_choice from CCPP data structure')
            return
        end if
        if (kind(ivegsrc).ne.ckind) then
            call ccpp_error('Kind mismatch for variable vegetation_type_dataset_choice')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'sea_ice_concentration', fice, ierr=ierr, dims=cdims, kind=ckind, index=628)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve sea_ice_concentration from CCPP data structure')
            return
        end if
        if (kind(fice).ne.ckind) then
            call ccpp_error('Kind mismatch for variable sea_ice_concentration')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'flag_for_frozen_soil_physics', keepfr, ierr=ierr, dims=cdims, kind=ckind, index=279)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_frozen_soil_physics from CCPP data structure')
            return
        end if
        if (kind(keepfr).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_frozen_soil_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'volume_fraction_of_soil_moisture_for_land_surface_model', smois, ierr=ierr, dims=cdims, kind=ckind, index=835)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve volume_fraction_of_soil_moisture_for_land_surface_model from CCPP data structure')
            return
        end if
        if (kind(smois).ne.ckind) then
            call ccpp_error('Kind mismatch for variable volume_fraction_of_soil_moisture_for_land_surface_model')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'volume_fraction_of_unfrozen_soil_moisture_for_land_surface_model', sh2o, ierr=ierr, dims=cdims, kind=ckind, index=837)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve volume_fraction_of_unfrozen_soil_moisture_for_land_surface_model from CCPP data structure')
            return
        end if
        if (kind(sh2o).ne.ckind) then
            call ccpp_error('Kind mismatch for variable volume_fraction_of_unfrozen_soil_moisture_for_land_surface_model')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'volume_fraction_of_frozen_soil_moisture_for_land_surface_model', smfrkeep, ierr=ierr, dims=cdims, kind=ckind, index=833)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve volume_fraction_of_frozen_soil_moisture_for_land_surface_model from CCPP data structure')
            return
        end if
        if (kind(smfrkeep).ne.ckind) then
            call ccpp_error('Kind mismatch for variable volume_fraction_of_frozen_soil_moisture_for_land_surface_model')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'soil_temperature_for_land_surface_model', tslb, ierr=ierr, dims=cdims, kind=ckind, index=656)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve soil_temperature_for_land_surface_model from CCPP data structure')
            return
        end if
        if (kind(tslb).ne.ckind) then
            call ccpp_error('Kind mismatch for variable soil_temperature_for_land_surface_model')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'soil_moisture_content', stm, ierr=ierr, dims=cdims, kind=ckind, index=654)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve soil_moisture_content from CCPP data structure')
            return
        end if
        if (kind(stm).ne.ckind) then
            call ccpp_error('Kind mismatch for variable soil_moisture_content')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_skin_temperature', tskin, ierr=ierr, dims=cdims, kind=ckind, index=721)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_skin_temperature from CCPP data structure')
            return
        end if
        if (kind(tskin).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_skin_temperature')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_skin_temperature_after_iteration', tsurf, ierr=ierr, dims=cdims, kind=ckind, index=722)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_skin_temperature_after_iteration from CCPP data structure')
            return
        end if
        if (kind(tsurf).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_skin_temperature_after_iteration')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'sea_ice_temperature', tice, ierr=ierr, dims=cdims, kind=ckind, index=629)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve sea_ice_temperature from CCPP data structure')
            return
        end if
        if (kind(tice).ne.ckind) then
            call ccpp_error('Kind mismatch for variable sea_ice_temperature')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'snow_temperature_bottom_first_layer', tsnow, ierr=ierr, dims=cdims, kind=ckind, index=652)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve snow_temperature_bottom_first_layer from CCPP data structure')
            return
        end if
        if (kind(tsnow).ne.ckind) then
            call ccpp_error('Kind mismatch for variable snow_temperature_bottom_first_layer')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'total_accumulated_snowfall', snowfallac, ierr=ierr, dims=cdims, kind=ckind, index=798)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve total_accumulated_snowfall from CCPP data structure')
            return
        end if
        if (kind(snowfallac).ne.ckind) then
            call ccpp_error('Kind mismatch for variable total_accumulated_snowfall')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'accumulated_water_equivalent_of_frozen_precip', acsnow, ierr=ierr, dims=cdims, kind=ckind, index=30)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve accumulated_water_equivalent_of_frozen_precip from CCPP data structure')
            return
        end if
        if (kind(acsnow).ne.ckind) then
            call ccpp_error('Kind mismatch for variable accumulated_water_equivalent_of_frozen_precip')
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
        

        call ccpp_field_get(cdata, 'soil_upward_latent_heat_flux', evbs, ierr=ierr, dims=cdims, kind=ckind, index=660)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve soil_upward_latent_heat_flux from CCPP data structure')
            return
        end if
        if (kind(evbs).ne.ckind) then
            call ccpp_error('Kind mismatch for variable soil_upward_latent_heat_flux')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'canopy_upward_latent_heat_flux', evcw, ierr=ierr, dims=cdims, kind=ckind, index=79)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve canopy_upward_latent_heat_flux from CCPP data structure')
            return
        end if
        if (kind(evcw).ne.ckind) then
            call ccpp_error('Kind mismatch for variable canopy_upward_latent_heat_flux')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'snow_deposition_sublimation_upward_latent_heat_flux', sbsno, ierr=ierr, dims=cdims, kind=ckind, index=649)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve snow_deposition_sublimation_upward_latent_heat_flux from CCPP data structure')
            return
        end if
        if (kind(sbsno).ne.ckind) then
            call ccpp_error('Kind mismatch for variable snow_deposition_sublimation_upward_latent_heat_flux')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'transpiration_flux', trans, ierr=ierr, dims=cdims, kind=ckind, index=804)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve transpiration_flux from CCPP data structure')
            return
        end if
        if (kind(trans).ne.ckind) then
            call ccpp_error('Kind mismatch for variable transpiration_flux')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_runoff_flux', runof, ierr=ierr, dims=cdims, kind=ckind, index=720)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_runoff_flux from CCPP data structure')
            return
        end if
        if (kind(runof).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_runoff_flux')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'subsurface_runoff_flux', drain, ierr=ierr, dims=cdims, kind=ckind, index=676)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve subsurface_runoff_flux from CCPP data structure')
            return
        end if
        if (kind(drain).ne.ckind) then
            call ccpp_error('Kind mismatch for variable subsurface_runoff_flux')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'total_runoff', runoff, ierr=ierr, dims=cdims, kind=ckind, index=800)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve total_runoff from CCPP data structure')
            return
        end if
        if (kind(runoff).ne.ckind) then
            call ccpp_error('Kind mismatch for variable total_runoff')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_runoff', srunoff, ierr=ierr, dims=cdims, kind=ckind, index=719)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_runoff from CCPP data structure')
            return
        end if
        if (kind(srunoff).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_runoff')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'upward_heat_flux_in_soil', gflux, ierr=ierr, dims=cdims, kind=ckind, index=812)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve upward_heat_flux_in_soil from CCPP data structure')
            return
        end if
        if (kind(gflux).ne.ckind) then
            call ccpp_error('Kind mismatch for variable upward_heat_flux_in_soil')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'minimum_vegetation_area_fraction', shdmin, ierr=ierr, dims=cdims, kind=ckind, index=547)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve minimum_vegetation_area_fraction from CCPP data structure')
            return
        end if
        if (kind(shdmin).ne.ckind) then
            call ccpp_error('Kind mismatch for variable minimum_vegetation_area_fraction')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'maximum_vegetation_area_fraction', shdmax, ierr=ierr, dims=cdims, kind=ckind, index=510)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve maximum_vegetation_area_fraction from CCPP data structure')
            return
        end if
        if (kind(shdmax).ne.ckind) then
            call ccpp_error('Kind mismatch for variable maximum_vegetation_area_fraction')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'flag_for_iteration', flag_iter, ierr=ierr, dims=cdims, kind=ckind, index=287)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_iteration from CCPP data structure')
            return
        end if
        if (kind(flag_iter).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_iteration')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'flag_for_guess_run', flag_guess, ierr=ierr, dims=cdims, kind=ckind, index=281)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_guess_run from CCPP data structure')
            return
        end if
        if (kind(flag_guess).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_guess_run')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

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
        

        call lsm_ruc_run(delt=delt,me=me,kdt=kdt,iter=cdata%loop_cnt,im=im,nlev=nlev,lsm_ruc=lsm_ruc, &
                  lsm=lsm,do_mynnsfclay=do_mynnsfclay,lsoil_ruc=lsoil_ruc,lsoil=lsoil,zs=zs, &
                  islmsk=islmsk,con_cp=con_cp,con_g=con_g,con_pi=con_pi,con_rd=con_rd,con_rv=con_rv, &
                  con_hvap=con_hvap,con_fvirt=con_fvirt,rainnc=rainnc,rainc=rainc,ice=ice, &
                  snow=snow,graupel=graupel,srflag=srflag,sncovr1=sncovr1,snowc=snowc,weasd=weasd, &
                  snwdph=snwdph,sr=sr,rhosnf=rhosnf,zf=zf,u1=u1,v1=v1,prsl1=prsl1,ddvel=ddvel, &
                  t1=t1,q1=q1,qc=qc,dlwflx=dlwflx,dswsfc=dswsfc,snet=snet,sfcemis=sfcemis, &
                  wspd=wspd,cm=cm,ch=ch,chh=chh,cmm=cmm,wet1=wet1,canopy=canopy,sigmaf=sigmaf, &
                  sfalb=sfalb,alvwf=alvwf,alnwf=alnwf,snoalb=snoalb,zorl=zorl,qsurf=qsurf, &
                  sfcqc=sfcqc,sfcqv=sfcqv,sfcdew=sfcdew,tg3=tg3,smc=smc,slc=slc,stc=stc,smcwlt2=smcwlt2, &
                  smcref2=smcref2,vegtype=vegtype,soiltyp=soiltyp,isot=isot,ivegsrc=ivegsrc, &
                  fice=fice,keepfr=keepfr,smois=smois,sh2o=sh2o,smfrkeep=smfrkeep,tslb=tslb, &
                  stm=stm,tskin=tskin,tsurf=tsurf,tice=tice,tsnow=tsnow,snowfallac=snowfallac, &
                  acsnow=acsnow,evap=evap,hflx=hflx,evbs=evbs,evcw=evcw,sbsno=sbsno,trans=trans, &
                  runof=runof,drain=drain,runoff=runoff,srunoff=srunoff,gflux=gflux,shdmin=shdmin, &
                  shdmax=shdmax,flag_iter=flag_iter,flag_guess=flag_guess,flag_init=flag_init, &
                  flag_restart=flag_restart,errmsg=cdata%errmsg,errflg=cdata%errflg)
        ierr=cdata%errflg

    end function lsm_ruc_run_cap

    function lsm_ruc_init_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind
        integer, pointer :: me
        integer, pointer :: isot
        integer, pointer :: ivegsrc
        integer, pointer :: nlunit

        ierr = 0

        call c_f_pointer(ptr, cdata)


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
        

        call ccpp_field_get(cdata, 'soil_type_dataset_choice', isot, ierr=ierr, kind=ckind, index=659)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve soil_type_dataset_choice from CCPP data structure')
            return
        end if
        if (kind(isot).ne.ckind) then
            call ccpp_error('Kind mismatch for variable soil_type_dataset_choice')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'vegetation_type_dataset_choice', ivegsrc, ierr=ierr, kind=ckind, index=816)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve vegetation_type_dataset_choice from CCPP data structure')
            return
        end if
        if (kind(ivegsrc).ne.ckind) then
            call ccpp_error('Kind mismatch for variable vegetation_type_dataset_choice')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'iounit_namelist', nlunit, ierr=ierr, kind=ckind, index=451)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve iounit_namelist from CCPP data structure')
            return
        end if
        if (kind(nlunit).ne.ckind) then
            call ccpp_error('Kind mismatch for variable iounit_namelist')
            ierr = 1
            return
        end if
#endif
        

        call lsm_ruc_init(me=me,isot=isot,ivegsrc=ivegsrc,nlunit=nlunit,errmsg=cdata%errmsg,errflg=cdata%errflg)
        ierr=cdata%errflg

    end function lsm_ruc_init_cap
end module lsm_ruc_cap

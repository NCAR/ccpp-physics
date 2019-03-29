
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
!! @brief Auto-generated cap module for the lsm_noah scheme
!!
!
module lsm_noah_cap

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_f_pointer, c_ptr, c_int32_t
    use            :: ccpp_types,                                      &
                      only: ccpp_t, CCPP_GENERIC_KIND
    use            :: ccpp_fields,                                     &
                      only: ccpp_field_get
    use            :: ccpp_errors,                                     &
                      only: ccpp_error, ccpp_debug
    use            :: lsm_noah, &
                      only: lsm_noah_init,lsm_noah_finalize,lsm_noah_run
    ! Other modules required, e.g. type definitions
    use machine, only: kind_phys

    implicit none

    private
    public :: lsm_noah_init_cap,lsm_noah_finalize_cap,lsm_noah_run_cap

    contains


    function lsm_noah_init_cap(ptr) bind(c) result(ierr)

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
        

        call lsm_noah_init(me=me,isot=isot,ivegsrc=ivegsrc,nlunit=nlunit,errmsg=cdata%errmsg,errflg=cdata%errflg)
        ierr=cdata%errflg

    end function lsm_noah_init_cap

    function lsm_noah_finalize_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call lsm_noah_finalize(errmsg=cdata%errmsg,errflg=cdata%errflg)
        ierr=cdata%errflg

    end function lsm_noah_finalize_cap

    function lsm_noah_run_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind
        integer, pointer :: im
        integer, pointer :: km
        real(kind_phys), pointer :: ps(:)
        real(kind_phys), pointer :: u1(:)
        real(kind_phys), pointer :: v1(:)
        real(kind_phys), pointer :: t1(:)
        real(kind_phys), pointer :: q1(:)
        integer, pointer :: soiltyp(:)
        integer, pointer :: vegtype(:)
        real(kind_phys), pointer :: sigmaf(:)
        real(kind_phys), pointer :: sfcemis(:)
        real(kind_phys), pointer :: dlwflx(:)
        real(kind_phys), pointer :: dswsfc(:)
        real(kind_phys), pointer :: snet(:)
        real(kind_phys), pointer :: delt
        real(kind_phys), pointer :: tg3(:)
        real(kind_phys), pointer :: cm(:)
        real(kind_phys), pointer :: ch(:)
        real(kind_phys), pointer :: prsl1(:)
        real(kind_phys), pointer :: prslki(:)
        real(kind_phys), pointer :: zf(:)
        integer, pointer :: islimsk(:)
        real(kind_phys), pointer :: ddvel(:)
        integer, pointer :: slopetyp(:)
        real(kind_phys), pointer :: shdmin(:)
        real(kind_phys), pointer :: shdmax(:)
        real(kind_phys), pointer :: snoalb(:)
        real(kind_phys), pointer :: sfalb(:)
        logical, pointer :: flag_iter(:)
        logical, pointer :: flag_guess(:)
        integer, pointer :: isot
        integer, pointer :: ivegsrc
        real(kind_phys), pointer :: bexppert(:)
        real(kind_phys), pointer :: xlaipert(:)
        real(kind_phys), pointer :: vegfpert(:)
        real(kind_phys), pointer :: pertvegf(:)
        real(kind_phys), pointer :: weasd(:)
        real(kind_phys), pointer :: snwdph(:)
        real(kind_phys), pointer :: tskin(:)
        real(kind_phys), pointer :: tprcp(:)
        real(kind_phys), pointer :: srflag(:)
        real(kind_phys), pointer :: smc(:,:)
        real(kind_phys), pointer :: stc(:,:)
        real(kind_phys), pointer :: slc(:,:)
        real(kind_phys), pointer :: canopy(:)
        real(kind_phys), pointer :: trans(:)
        real(kind_phys), pointer :: tsurf(:)
        real(kind_phys), pointer :: zorl(:)
        real(kind_phys), pointer :: sncovr1(:)
        real(kind_phys), pointer :: qsurf(:)
        real(kind_phys), pointer :: gflux(:)
        real(kind_phys), pointer :: drain(:)
        real(kind_phys), pointer :: evap(:)
        real(kind_phys), pointer :: hflx(:)
        real(kind_phys), pointer :: ep(:)
        real(kind_phys), pointer :: runoff(:)
        real(kind_phys), pointer :: cmm(:)
        real(kind_phys), pointer :: chh(:)
        real(kind_phys), pointer :: evbs(:)
        real(kind_phys), pointer :: evcw(:)
        real(kind_phys), pointer :: sbsno(:)
        real(kind_phys), pointer :: snowc(:)
        real(kind_phys), pointer :: stm(:)
        real(kind_phys), pointer :: snohf(:)
        real(kind_phys), pointer :: smcwlt2(:)
        real(kind_phys), pointer :: smcref2(:)
        real(kind_phys), pointer :: wet1(:)

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
        

        call ccpp_field_get(cdata, 'soil_vertical_dimension', km, ierr=ierr, kind=ckind, index=661)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve soil_vertical_dimension from CCPP data structure')
            return
        end if
        if (kind(km).ne.ckind) then
            call ccpp_error('Kind mismatch for variable soil_vertical_dimension')
            ierr = 1
            return
        end if
#endif
        

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
        

        call ccpp_field_get(cdata, 'bounded_vegetation_area_fraction', sigmaf, ierr=ierr, dims=cdims, kind=ckind, index=77)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve bounded_vegetation_area_fraction from CCPP data structure')
            return
        end if
        if (kind(sigmaf).ne.ckind) then
            call ccpp_error('Kind mismatch for variable bounded_vegetation_area_fraction')
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
        

        call ccpp_field_get(cdata, 'surface_downwelling_longwave_flux_absorbed_by_ground', dlwflx, ierr=ierr, dims=cdims, kind=ckind, index=698)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_downwelling_longwave_flux_absorbed_by_ground from CCPP data structure')
            return
        end if
        if (kind(dlwflx).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_downwelling_longwave_flux_absorbed_by_ground')
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
        

        call ccpp_field_get(cdata, 'ratio_of_exner_function_between_midlayer_and_interface_at_lowest_model_layer', prslki, ierr=ierr, dims=cdims, kind=ckind, index=623)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve ratio_of_exner_function_between_midlayer_and_interface_at_lowest_model_layer from CCPP data structure')
            return
        end if
        if (kind(prslki).ne.ckind) then
            call ccpp_error('Kind mismatch for variable ratio_of_exner_function_between_midlayer_and_interface_at_lowest_model_layer')
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
        

        call ccpp_field_get(cdata, 'sea_land_ice_mask', islimsk, ierr=ierr, dims=cdims, kind=ckind, index=631)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve sea_land_ice_mask from CCPP data structure')
            return
        end if
        if (kind(islimsk).ne.ckind) then
            call ccpp_error('Kind mismatch for variable sea_land_ice_mask')
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
        

        call ccpp_field_get(cdata, 'surface_slope_classification', slopetyp, ierr=ierr, dims=cdims, kind=ckind, index=724)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_slope_classification from CCPP data structure')
            return
        end if
        if (kind(slopetyp).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_slope_classification')
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
        

        call ccpp_field_get(cdata, 'perturbation_of_soil_type_b_parameter', bexppert, ierr=ierr, dims=cdims, kind=ckind, index=604)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve perturbation_of_soil_type_b_parameter from CCPP data structure')
            return
        end if
        if (kind(bexppert).ne.ckind) then
            call ccpp_error('Kind mismatch for variable perturbation_of_soil_type_b_parameter')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'perturbation_of_leaf_area_index', xlaipert, ierr=ierr, dims=cdims, kind=ckind, index=602)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve perturbation_of_leaf_area_index from CCPP data structure')
            return
        end if
        if (kind(xlaipert).ne.ckind) then
            call ccpp_error('Kind mismatch for variable perturbation_of_leaf_area_index')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'perturbation_of_vegetation_fraction', vegfpert, ierr=ierr, dims=cdims, kind=ckind, index=605)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve perturbation_of_vegetation_fraction from CCPP data structure')
            return
        end if
        if (kind(vegfpert).ne.ckind) then
            call ccpp_error('Kind mismatch for variable perturbation_of_vegetation_fraction')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'magnitude_of_perturbation_of_vegetation_fraction', pertvegf, ierr=ierr, dims=cdims, kind=ckind, index=500)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve magnitude_of_perturbation_of_vegetation_fraction from CCPP data structure')
            return
        end if
        if (kind(pertvegf).ne.ckind) then
            call ccpp_error('Kind mismatch for variable magnitude_of_perturbation_of_vegetation_fraction')
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
        

        call ccpp_field_get(cdata, 'surface_upward_potential_latent_heat_flux', ep, ierr=ierr, dims=cdims, kind=ckind, index=732)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_upward_potential_latent_heat_flux from CCPP data structure')
            return
        end if
        if (kind(ep).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_upward_potential_latent_heat_flux')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_runoff_flux', runoff, ierr=ierr, dims=cdims, kind=ckind, index=720)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_runoff_flux from CCPP data structure')
            return
        end if
        if (kind(runoff).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_runoff_flux')
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
        

        call ccpp_field_get(cdata, 'snow_freezing_rain_upward_latent_heat_flux', snohf, ierr=ierr, dims=cdims, kind=ckind, index=650)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve snow_freezing_rain_upward_latent_heat_flux from CCPP data structure')
            return
        end if
        if (kind(snohf).ne.ckind) then
            call ccpp_error('Kind mismatch for variable snow_freezing_rain_upward_latent_heat_flux')
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
        

        call lsm_noah_run(im=im,km=km,ps=ps,u1=u1,v1=v1,t1=t1,q1=q1,soiltyp=soiltyp,vegtype=vegtype, &
                  sigmaf=sigmaf,sfcemis=sfcemis,dlwflx=dlwflx,dswsfc=dswsfc,snet=snet,delt=delt, &
                  tg3=tg3,cm=cm,ch=ch,prsl1=prsl1,prslki=prslki,zf=zf,islimsk=islimsk,ddvel=ddvel, &
                  slopetyp=slopetyp,shdmin=shdmin,shdmax=shdmax,snoalb=snoalb,sfalb=sfalb, &
                  flag_iter=flag_iter,flag_guess=flag_guess,isot=isot,ivegsrc=ivegsrc,bexppert=bexppert, &
                  xlaipert=xlaipert,vegfpert=vegfpert,pertvegf=pertvegf,weasd=weasd,snwdph=snwdph, &
                  tskin=tskin,tprcp=tprcp,srflag=srflag,smc=smc,stc=stc,slc=slc,canopy=canopy, &
                  trans=trans,tsurf=tsurf,zorl=zorl,sncovr1=sncovr1,qsurf=qsurf,gflux=gflux, &
                  drain=drain,evap=evap,hflx=hflx,ep=ep,runoff=runoff,cmm=cmm,chh=chh,evbs=evbs, &
                  evcw=evcw,sbsno=sbsno,snowc=snowc,stm=stm,snohf=snohf,smcwlt2=smcwlt2,smcref2=smcref2, &
                  wet1=wet1,errmsg=cdata%errmsg,errflg=cdata%errflg)
        ierr=cdata%errflg

    end function lsm_noah_run_cap
end module lsm_noah_cap

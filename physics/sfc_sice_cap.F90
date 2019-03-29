
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
!! @brief Auto-generated cap module for the sfc_sice scheme
!!
!
module sfc_sice_cap

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_f_pointer, c_ptr, c_int32_t
    use            :: ccpp_types,                                      &
                      only: ccpp_t, CCPP_GENERIC_KIND
    use            :: ccpp_fields,                                     &
                      only: ccpp_field_get
    use            :: ccpp_errors,                                     &
                      only: ccpp_error, ccpp_debug
    use            :: sfc_sice, &
                      only: sfc_sice_init,sfc_sice_finalize,sfc_sice_run
    ! Other modules required, e.g. type definitions
    use machine, only: kind_phys

    implicit none

    private
    public :: sfc_sice_init_cap,sfc_sice_finalize_cap,sfc_sice_run_cap

    contains


    function sfc_sice_init_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call sfc_sice_init()
        

    end function sfc_sice_init_cap

    function sfc_sice_finalize_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call sfc_sice_finalize()
        

    end function sfc_sice_finalize_cap

    function sfc_sice_run_cap(ptr) bind(c) result(ierr)

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
        real(kind_phys), pointer :: delt
        real(kind_phys), pointer :: sfcemis(:)
        real(kind_phys), pointer :: dlwflx(:)
        real(kind_phys), pointer :: sfcnsw(:)
        real(kind_phys), pointer :: sfcdsw(:)
        real(kind_phys), pointer :: srflag(:)
        real(kind_phys), pointer :: cm(:)
        real(kind_phys), pointer :: ch(:)
        real(kind_phys), pointer :: prsl1(:)
        real(kind_phys), pointer :: prslki(:)
        integer, pointer :: islimsk(:)
        real(kind_phys), pointer :: ddvel(:)
        logical, pointer :: flag_iter(:)
        logical, pointer :: mom4ice
        integer, pointer :: lsm
        logical, pointer :: lprnt
        integer, pointer :: ipr
        real(kind_phys), pointer :: hice(:)
        real(kind_phys), pointer :: fice(:)
        real(kind_phys), pointer :: tice(:)
        real(kind_phys), pointer :: weasd(:)
        real(kind_phys), pointer :: tskin(:)
        real(kind_phys), pointer :: tprcp(:)
        real(kind_phys), pointer :: stc(:,:)
        real(kind_phys), pointer :: ep(:)
        real(kind_phys), pointer :: snwdph(:)
        real(kind_phys), pointer :: qsurf(:)
        real(kind_phys), pointer :: snowmt(:)
        real(kind_phys), pointer :: gflux(:)
        real(kind_phys), pointer :: cmm(:)
        real(kind_phys), pointer :: chh(:)
        real(kind_phys), pointer :: evap(:)
        real(kind_phys), pointer :: hflx(:)

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
        

        call ccpp_field_get(cdata, 'surface_net_downwelling_shortwave_flux', sfcnsw, ierr=ierr, dims=cdims, kind=ckind, index=716)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_net_downwelling_shortwave_flux from CCPP data structure')
            return
        end if
        if (kind(sfcnsw).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_net_downwelling_shortwave_flux')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_downwelling_shortwave_flux', sfcdsw, ierr=ierr, dims=cdims, kind=ckind, index=700)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_downwelling_shortwave_flux from CCPP data structure')
            return
        end if
        if (kind(sfcdsw).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_downwelling_shortwave_flux')
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
        

        call ccpp_field_get(cdata, 'flag_for_mom4_coupling', mom4ice, ierr=ierr, kind=ckind, index=295)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_mom4_coupling from CCPP data structure')
            return
        end if
        if (kind(mom4ice).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_mom4_coupling')
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
        

        call ccpp_field_get(cdata, 'horizontal_index_of_printed_column', ipr, ierr=ierr, kind=ckind, index=365)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve horizontal_index_of_printed_column from CCPP data structure')
            return
        end if
        if (kind(ipr).ne.ckind) then
            call ccpp_error('Kind mismatch for variable horizontal_index_of_printed_column')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'sea_ice_thickness', hice, ierr=ierr, dims=cdims, kind=ckind, index=630)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve sea_ice_thickness from CCPP data structure')
            return
        end if
        if (kind(hice).ne.ckind) then
            call ccpp_error('Kind mismatch for variable sea_ice_thickness')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

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
        

        call ccpp_field_get(cdata, 'surface_snow_melt', snowmt, ierr=ierr, dims=cdims, kind=ckind, index=728)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_snow_melt from CCPP data structure')
            return
        end if
        if (kind(snowmt).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_snow_melt')
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
        

        call sfc_sice_run(im=im,km=km,ps=ps,u1=u1,v1=v1,t1=t1,q1=q1,delt=delt,sfcemis=sfcemis,dlwflx=dlwflx, &
                  sfcnsw=sfcnsw,sfcdsw=sfcdsw,srflag=srflag,cm=cm,ch=ch,prsl1=prsl1,prslki=prslki, &
                  islimsk=islimsk,ddvel=ddvel,flag_iter=flag_iter,mom4ice=mom4ice,lsm=lsm, &
                  lprnt=lprnt,ipr=ipr,hice=hice,fice=fice,tice=tice,weasd=weasd,tskin=tskin, &
                  tprcp=tprcp,stc=stc,ep=ep,snwdph=snwdph,qsurf=qsurf,snowmt=snowmt,gflux=gflux, &
                  cmm=cmm,chh=chh,evap=evap,hflx=hflx,errmsg=cdata%errmsg,errflg=cdata%errflg)
        ierr=cdata%errflg

    end function sfc_sice_run_cap
end module sfc_sice_cap

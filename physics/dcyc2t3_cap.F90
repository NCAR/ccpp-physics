
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
!! @brief Auto-generated cap module for the dcyc2t3 scheme
!!
!
module dcyc2t3_cap

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_f_pointer, c_ptr, c_int32_t
    use            :: ccpp_types,                                      &
                      only: ccpp_t, CCPP_GENERIC_KIND
    use            :: ccpp_fields,                                     &
                      only: ccpp_field_get
    use            :: ccpp_errors,                                     &
                      only: ccpp_error, ccpp_debug
    use            :: dcyc2t3, &
                      only: dcyc2t3_finalize,dcyc2t3_run,dcyc2t3_init
    ! Other modules required, e.g. type definitions
    use machine, only: kind_phys

    implicit none

    private
    public :: dcyc2t3_finalize_cap,dcyc2t3_run_cap,dcyc2t3_init_cap

    contains


    function dcyc2t3_finalize_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call dcyc2t3_finalize()
        

    end function dcyc2t3_finalize_cap

    function dcyc2t3_run_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind
        real(kind_phys), pointer :: solhr
        real(kind_phys), pointer :: slag
        real(kind_phys), pointer :: sdec
        real(kind_phys), pointer :: cdec
        real(kind_phys), pointer :: sinlat(:)
        real(kind_phys), pointer :: coslat(:)
        real(kind_phys), pointer :: xlon(:)
        real(kind_phys), pointer :: coszen(:)
        real(kind_phys), pointer :: tsea(:)
        real(kind_phys), pointer :: tf(:)
        real(kind_phys), pointer :: tsflw(:)
        real(kind_phys), pointer :: sfcemis(:)
        real(kind_phys), pointer :: sfcdsw(:)
        real(kind_phys), pointer :: sfcnsw(:)
        real(kind_phys), pointer :: sfcdlw(:)
        real(kind_phys), pointer :: swh(:,:)
        real(kind_phys), pointer :: swhc(:,:)
        real(kind_phys), pointer :: hlw(:,:)
        real(kind_phys), pointer :: hlwc(:,:)
        real(kind_phys), pointer :: sfcnirbmu(:)
        real(kind_phys), pointer :: sfcnirdfu(:)
        real(kind_phys), pointer :: sfcvisbmu(:)
        real(kind_phys), pointer :: sfcvisdfu(:)
        real(kind_phys), pointer :: sfcnirbmd(:)
        real(kind_phys), pointer :: sfcnirdfd(:)
        real(kind_phys), pointer :: sfcvisbmd(:)
        real(kind_phys), pointer :: sfcvisdfd(:)
        integer, pointer :: ix
        integer, pointer :: im
        integer, pointer :: levs
        real(kind_phys), pointer :: deltim
        real(kind_phys), pointer :: dtdt(:,:)
        real(kind_phys), pointer :: dtdtc(:,:)
        real(kind_phys), pointer :: adjsfcdsw(:)
        real(kind_phys), pointer :: adjsfcnsw(:)
        real(kind_phys), pointer :: adjsfcdlw(:)
        real(kind_phys), pointer :: adjsfculw(:)
        real(kind_phys), pointer :: xmu(:)
        real(kind_phys), pointer :: xcosz(:)
        real(kind_phys), pointer :: adjnirbmu(:)
        real(kind_phys), pointer :: adjnirdfu(:)
        real(kind_phys), pointer :: adjvisbmu(:)
        real(kind_phys), pointer :: adjvisdfu(:)
        real(kind_phys), pointer :: adjnirbmd(:)
        real(kind_phys), pointer :: adjnirdfd(:)
        real(kind_phys), pointer :: adjvisbmd(:)
        real(kind_phys), pointer :: adjvisdfd(:)

        ierr = 0

        call c_f_pointer(ptr, cdata)


        call ccpp_field_get(cdata, 'forecast_hour', solhr, ierr=ierr, kind=ckind, index=338)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve forecast_hour from CCPP data structure')
            return
        end if
        if (kind(solhr).ne.ckind) then
            call ccpp_error('Kind mismatch for variable forecast_hour')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'equation_of_time', slag, ierr=ierr, kind=ckind, index=256)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve equation_of_time from CCPP data structure')
            return
        end if
        if (kind(slag).ne.ckind) then
            call ccpp_error('Kind mismatch for variable equation_of_time')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'sine_of_solar_declination_angle', sdec, ierr=ierr, kind=ckind, index=646)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve sine_of_solar_declination_angle from CCPP data structure')
            return
        end if
        if (kind(sdec).ne.ckind) then
            call ccpp_error('Kind mismatch for variable sine_of_solar_declination_angle')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'cosine_of_solar_declination_angle', cdec, ierr=ierr, kind=ckind, index=135)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cosine_of_solar_declination_angle from CCPP data structure')
            return
        end if
        if (kind(cdec).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cosine_of_solar_declination_angle')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'sine_of_latitude', sinlat, ierr=ierr, dims=cdims, kind=ckind, index=645)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve sine_of_latitude from CCPP data structure')
            return
        end if
        if (kind(sinlat).ne.ckind) then
            call ccpp_error('Kind mismatch for variable sine_of_latitude')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cosine_of_latitude', coslat, ierr=ierr, dims=cdims, kind=ckind, index=134)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cosine_of_latitude from CCPP data structure')
            return
        end if
        if (kind(coslat).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cosine_of_latitude')
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
        

        call ccpp_field_get(cdata, 'cosine_of_zenith_angle', coszen, ierr=ierr, dims=cdims, kind=ckind, index=136)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cosine_of_zenith_angle from CCPP data structure')
            return
        end if
        if (kind(coszen).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cosine_of_zenith_angle')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_skin_temperature', tsea, ierr=ierr, dims=cdims, kind=ckind, index=721)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_skin_temperature from CCPP data structure')
            return
        end if
        if (kind(tsea).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_skin_temperature')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'air_temperature_at_lowest_model_layer', tf, ierr=ierr, dims=cdims, kind=ckind, index=53)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve air_temperature_at_lowest_model_layer from CCPP data structure')
            return
        end if
        if (kind(tf).ne.ckind) then
            call ccpp_error('Kind mismatch for variable air_temperature_at_lowest_model_layer')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_midlayer_air_temperature_in_longwave_radiation', tsflw, ierr=ierr, dims=cdims, kind=ckind, index=715)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_midlayer_air_temperature_in_longwave_radiation from CCPP data structure')
            return
        end if
        if (kind(tsflw).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_midlayer_air_temperature_in_longwave_radiation')
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
        

        call ccpp_field_get(cdata, 'surface_downwelling_shortwave_flux_on_radiation_time_step', sfcdsw, ierr=ierr, dims=cdims, kind=ckind, index=701)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_downwelling_shortwave_flux_on_radiation_time_step from CCPP data structure')
            return
        end if
        if (kind(sfcdsw).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_downwelling_shortwave_flux_on_radiation_time_step')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_net_downwelling_shortwave_flux_on_radiation_time_step', sfcnsw, ierr=ierr, dims=cdims, kind=ckind, index=717)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_net_downwelling_shortwave_flux_on_radiation_time_step from CCPP data structure')
            return
        end if
        if (kind(sfcnsw).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_net_downwelling_shortwave_flux_on_radiation_time_step')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_downwelling_longwave_flux_on_radiation_time_step', sfcdlw, ierr=ierr, dims=cdims, kind=ckind, index=699)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_downwelling_longwave_flux_on_radiation_time_step from CCPP data structure')
            return
        end if
        if (kind(sfcdlw).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_downwelling_longwave_flux_on_radiation_time_step')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'tendency_of_air_temperature_due_to_shortwave_heating_on_radiation_time_step', swh, ierr=ierr, dims=cdims, kind=ckind, index=762)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve tendency_of_air_temperature_due_to_shortwave_heating_on_radiation_time_step from CCPP data structure')
            return
        end if
        if (kind(swh).ne.ckind) then
            call ccpp_error('Kind mismatch for variable tendency_of_air_temperature_due_to_shortwave_heating_on_radiation_time_step')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'tendency_of_air_temperature_due_to_shortwave_heating_assuming_clear_sky_on_radiation_time_step', swhc, ierr=ierr, dims=cdims, kind=ckind, index=761)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve tendency_of_air_temperature_due_to_shortwave_heating_assuming_clear_sky_on_radiation_time_step from CCPP data structure')
            return
        end if
        if (kind(swhc).ne.ckind) then
            call ccpp_error('Kind mismatch for variable tendency_of_air_temperature_due_to_shortwave_heating_assuming_clear_sky_on_radiation_time_step')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'tendency_of_air_temperature_due_to_longwave_heating_on_radiation_time_step', hlw, ierr=ierr, dims=cdims, kind=ckind, index=756)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve tendency_of_air_temperature_due_to_longwave_heating_on_radiation_time_step from CCPP data structure')
            return
        end if
        if (kind(hlw).ne.ckind) then
            call ccpp_error('Kind mismatch for variable tendency_of_air_temperature_due_to_longwave_heating_on_radiation_time_step')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'tendency_of_air_temperature_due_to_longwave_heating_assuming_clear_sky_on_radiation_time_step', hlwc, ierr=ierr, dims=cdims, kind=ckind, index=754)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve tendency_of_air_temperature_due_to_longwave_heating_assuming_clear_sky_on_radiation_time_step from CCPP data structure')
            return
        end if
        if (kind(hlwc).ne.ckind) then
            call ccpp_error('Kind mismatch for variable tendency_of_air_temperature_due_to_longwave_heating_assuming_clear_sky_on_radiation_time_step')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_upwelling_direct_near_infrared_shortwave_flux_on_radiation_time_step', sfcnirbmu, ierr=ierr, dims=cdims, kind=ckind, index=738)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_upwelling_direct_near_infrared_shortwave_flux_on_radiation_time_step from CCPP data structure')
            return
        end if
        if (kind(sfcnirbmu).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_upwelling_direct_near_infrared_shortwave_flux_on_radiation_time_step')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_upwelling_diffuse_near_infrared_shortwave_flux_on_radiation_time_step', sfcnirdfu, ierr=ierr, dims=cdims, kind=ckind, index=734)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_upwelling_diffuse_near_infrared_shortwave_flux_on_radiation_time_step from CCPP data structure')
            return
        end if
        if (kind(sfcnirdfu).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_upwelling_diffuse_near_infrared_shortwave_flux_on_radiation_time_step')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_upwelling_direct_ultraviolet_and_visible_shortwave_flux_on_radiation_time_step', sfcvisbmu, ierr=ierr, dims=cdims, kind=ckind, index=740)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_upwelling_direct_ultraviolet_and_visible_shortwave_flux_on_radiation_time_step from CCPP data structure')
            return
        end if
        if (kind(sfcvisbmu).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_upwelling_direct_ultraviolet_and_visible_shortwave_flux_on_radiation_time_step')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_upwelling_diffuse_ultraviolet_and_visible_shortwave_flux_on_radiation_time_step', sfcvisdfu, ierr=ierr, dims=cdims, kind=ckind, index=736)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_upwelling_diffuse_ultraviolet_and_visible_shortwave_flux_on_radiation_time_step from CCPP data structure')
            return
        end if
        if (kind(sfcvisdfu).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_upwelling_diffuse_ultraviolet_and_visible_shortwave_flux_on_radiation_time_step')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_downwelling_direct_near_infrared_shortwave_flux_on_radiation_time_step', sfcnirbmd, ierr=ierr, dims=cdims, kind=ckind, index=694)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_downwelling_direct_near_infrared_shortwave_flux_on_radiation_time_step from CCPP data structure')
            return
        end if
        if (kind(sfcnirbmd).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_downwelling_direct_near_infrared_shortwave_flux_on_radiation_time_step')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_downwelling_diffuse_near_infrared_shortwave_flux_on_radiation_time_step', sfcnirdfd, ierr=ierr, dims=cdims, kind=ckind, index=690)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_downwelling_diffuse_near_infrared_shortwave_flux_on_radiation_time_step from CCPP data structure')
            return
        end if
        if (kind(sfcnirdfd).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_downwelling_diffuse_near_infrared_shortwave_flux_on_radiation_time_step')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_downwelling_direct_ultraviolet_and_visible_shortwave_flux_on_radiation_time_step', sfcvisbmd, ierr=ierr, dims=cdims, kind=ckind, index=696)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_downwelling_direct_ultraviolet_and_visible_shortwave_flux_on_radiation_time_step from CCPP data structure')
            return
        end if
        if (kind(sfcvisbmd).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_downwelling_direct_ultraviolet_and_visible_shortwave_flux_on_radiation_time_step')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_downwelling_diffuse_ultraviolet_and_visible_shortwave_flux_on_radiation_time_step', sfcvisdfd, ierr=ierr, dims=cdims, kind=ckind, index=692)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_downwelling_diffuse_ultraviolet_and_visible_shortwave_flux_on_radiation_time_step from CCPP data structure')
            return
        end if
        if (kind(sfcvisdfd).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_downwelling_diffuse_ultraviolet_and_visible_shortwave_flux_on_radiation_time_step')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

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
        

        call ccpp_field_get(cdata, 'time_step_for_dynamics', deltim, ierr=ierr, kind=ckind, index=792)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve time_step_for_dynamics from CCPP data structure')
            return
        end if
        if (kind(deltim).ne.ckind) then
            call ccpp_error('Kind mismatch for variable time_step_for_dynamics')
            ierr = 1
            return
        end if
#endif
        

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
        

        call ccpp_field_get(cdata, 'surface_downwelling_shortwave_flux', adjsfcdsw, ierr=ierr, dims=cdims, kind=ckind, index=700)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_downwelling_shortwave_flux from CCPP data structure')
            return
        end if
        if (kind(adjsfcdsw).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_downwelling_shortwave_flux')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_net_downwelling_shortwave_flux', adjsfcnsw, ierr=ierr, dims=cdims, kind=ckind, index=716)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_net_downwelling_shortwave_flux from CCPP data structure')
            return
        end if
        if (kind(adjsfcnsw).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_net_downwelling_shortwave_flux')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_downwelling_longwave_flux', adjsfcdlw, ierr=ierr, dims=cdims, kind=ckind, index=697)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_downwelling_longwave_flux from CCPP data structure')
            return
        end if
        if (kind(adjsfcdlw).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_downwelling_longwave_flux')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_upwelling_longwave_flux', adjsfculw, ierr=ierr, dims=cdims, kind=ckind, index=741)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_upwelling_longwave_flux from CCPP data structure')
            return
        end if
        if (kind(adjsfculw).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_upwelling_longwave_flux')
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
        

        call ccpp_field_get(cdata, 'instantaneous_cosine_of_zenith_angle', xcosz, ierr=ierr, dims=cdims, kind=ckind, index=408)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve instantaneous_cosine_of_zenith_angle from CCPP data structure')
            return
        end if
        if (kind(xcosz).ne.ckind) then
            call ccpp_error('Kind mismatch for variable instantaneous_cosine_of_zenith_angle')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_upwelling_direct_near_infrared_shortwave_flux', adjnirbmu, ierr=ierr, dims=cdims, kind=ckind, index=737)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_upwelling_direct_near_infrared_shortwave_flux from CCPP data structure')
            return
        end if
        if (kind(adjnirbmu).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_upwelling_direct_near_infrared_shortwave_flux')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_upwelling_diffuse_near_infrared_shortwave_flux', adjnirdfu, ierr=ierr, dims=cdims, kind=ckind, index=733)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_upwelling_diffuse_near_infrared_shortwave_flux from CCPP data structure')
            return
        end if
        if (kind(adjnirdfu).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_upwelling_diffuse_near_infrared_shortwave_flux')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_upwelling_direct_ultraviolet_and_visible_shortwave_flux', adjvisbmu, ierr=ierr, dims=cdims, kind=ckind, index=739)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_upwelling_direct_ultraviolet_and_visible_shortwave_flux from CCPP data structure')
            return
        end if
        if (kind(adjvisbmu).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_upwelling_direct_ultraviolet_and_visible_shortwave_flux')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_upwelling_diffuse_ultraviolet_and_visible_shortwave_flux', adjvisdfu, ierr=ierr, dims=cdims, kind=ckind, index=735)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_upwelling_diffuse_ultraviolet_and_visible_shortwave_flux from CCPP data structure')
            return
        end if
        if (kind(adjvisdfu).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_upwelling_diffuse_ultraviolet_and_visible_shortwave_flux')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_downwelling_direct_near_infrared_shortwave_flux', adjnirbmd, ierr=ierr, dims=cdims, kind=ckind, index=693)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_downwelling_direct_near_infrared_shortwave_flux from CCPP data structure')
            return
        end if
        if (kind(adjnirbmd).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_downwelling_direct_near_infrared_shortwave_flux')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_downwelling_diffuse_near_infrared_shortwave_flux', adjnirdfd, ierr=ierr, dims=cdims, kind=ckind, index=689)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_downwelling_diffuse_near_infrared_shortwave_flux from CCPP data structure')
            return
        end if
        if (kind(adjnirdfd).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_downwelling_diffuse_near_infrared_shortwave_flux')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_downwelling_direct_ultraviolet_and_visible_shortwave_flux', adjvisbmd, ierr=ierr, dims=cdims, kind=ckind, index=695)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_downwelling_direct_ultraviolet_and_visible_shortwave_flux from CCPP data structure')
            return
        end if
        if (kind(adjvisbmd).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_downwelling_direct_ultraviolet_and_visible_shortwave_flux')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_downwelling_diffuse_ultraviolet_and_visible_shortwave_flux', adjvisdfd, ierr=ierr, dims=cdims, kind=ckind, index=691)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_downwelling_diffuse_ultraviolet_and_visible_shortwave_flux from CCPP data structure')
            return
        end if
        if (kind(adjvisdfd).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_downwelling_diffuse_ultraviolet_and_visible_shortwave_flux')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call dcyc2t3_run(solhr=solhr,slag=slag,sdec=sdec,cdec=cdec,sinlat=sinlat,coslat=coslat,xlon=xlon, &
                  coszen=coszen,tsea=tsea,tf=tf,tsflw=tsflw,sfcemis=sfcemis,sfcdsw=sfcdsw, &
                  sfcnsw=sfcnsw,sfcdlw=sfcdlw,swh=swh,swhc=swhc,hlw=hlw,hlwc=hlwc,sfcnirbmu=sfcnirbmu, &
                  sfcnirdfu=sfcnirdfu,sfcvisbmu=sfcvisbmu,sfcvisdfu=sfcvisdfu,sfcnirbmd=sfcnirbmd, &
                  sfcnirdfd=sfcnirdfd,sfcvisbmd=sfcvisbmd,sfcvisdfd=sfcvisdfd,ix=ix,im=im, &
                  levs=levs,deltim=deltim,dtdt=dtdt,dtdtc=dtdtc,adjsfcdsw=adjsfcdsw,adjsfcnsw=adjsfcnsw, &
                  adjsfcdlw=adjsfcdlw,adjsfculw=adjsfculw,xmu=xmu,xcosz=xcosz,adjnirbmu=adjnirbmu, &
                  adjnirdfu=adjnirdfu,adjvisbmu=adjvisbmu,adjvisdfu=adjvisdfu,adjnirbmd=adjnirbmd, &
                  adjnirdfd=adjnirdfd,adjvisbmd=adjvisbmd,adjvisdfd=adjvisdfd,errmsg=cdata%errmsg, &
                  errflg=cdata%errflg)
        ierr=cdata%errflg

    end function dcyc2t3_run_cap

    function dcyc2t3_init_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call dcyc2t3_init()
        

    end function dcyc2t3_init_cap
end module dcyc2t3_cap

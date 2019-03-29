
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
!! @brief Auto-generated cap module for the GFS_surface_generic_post scheme
!!
!
module GFS_surface_generic_post_cap

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_f_pointer, c_ptr, c_int32_t
    use            :: ccpp_types,                                      &
                      only: ccpp_t, CCPP_GENERIC_KIND
    use            :: ccpp_fields,                                     &
                      only: ccpp_field_get
    use            :: ccpp_errors,                                     &
                      only: ccpp_error, ccpp_debug
    use            :: GFS_surface_generic_post, &
                      only: GFS_surface_generic_post_finalize,GFS_surface_generic_post_init,GFS_surface_generic_post_run
    ! Other modules required, e.g. type definitions
    use machine, only: kind_phys

    implicit none

    private
    public :: GFS_surface_generic_post_finalize_cap,GFS_surface_generic_post_init_cap,GFS_surface_generic_post_run_cap

    contains


    function GFS_surface_generic_post_finalize_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call GFS_surface_generic_post_finalize()
        

    end function GFS_surface_generic_post_finalize_cap

    function GFS_surface_generic_post_init_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call GFS_surface_generic_post_init()
        

    end function GFS_surface_generic_post_init_cap

    function GFS_surface_generic_post_run_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind
        integer, pointer :: im
        logical, pointer :: cplflx
        logical, pointer :: lssav
        integer, pointer :: islmsk(:)
        real(kind_phys), pointer :: dtf
        real(kind_phys), pointer :: ep1d(:)
        real(kind_phys), pointer :: gflx(:)
        real(kind_phys), pointer :: tgrs_1(:)
        real(kind_phys), pointer :: qgrs_1(:)
        real(kind_phys), pointer :: ugrs_1(:)
        real(kind_phys), pointer :: vgrs_1(:)
        real(kind_phys), pointer :: adjsfcdlw(:)
        real(kind_phys), pointer :: adjsfcdsw(:)
        real(kind_phys), pointer :: adjnirbmd(:)
        real(kind_phys), pointer :: adjnirdfd(:)
        real(kind_phys), pointer :: adjvisbmd(:)
        real(kind_phys), pointer :: adjvisdfd(:)
        real(kind_phys), pointer :: adjsfculw(:)
        real(kind_phys), pointer :: adjnirbmu(:)
        real(kind_phys), pointer :: adjnirdfu(:)
        real(kind_phys), pointer :: adjvisbmu(:)
        real(kind_phys), pointer :: adjvisdfu(:)
        real(kind_phys), pointer :: t2m(:)
        real(kind_phys), pointer :: q2m(:)
        real(kind_phys), pointer :: u10m(:)
        real(kind_phys), pointer :: v10m(:)
        real(kind_phys), pointer :: tsfc(:)
        real(kind_phys), pointer :: pgr(:)
        real(kind_phys), pointer :: xcosz(:)
        real(kind_phys), pointer :: evbs(:)
        real(kind_phys), pointer :: evcw(:)
        real(kind_phys), pointer :: trans(:)
        real(kind_phys), pointer :: sbsno(:)
        real(kind_phys), pointer :: snowc(:)
        real(kind_phys), pointer :: snohf(:)
        real(kind_phys), pointer :: epi(:)
        real(kind_phys), pointer :: gfluxi(:)
        real(kind_phys), pointer :: t1(:)
        real(kind_phys), pointer :: q1(:)
        real(kind_phys), pointer :: u1(:)
        real(kind_phys), pointer :: v1(:)
        real(kind_phys), pointer :: dlwsfci_cpl(:)
        real(kind_phys), pointer :: dswsfci_cpl(:)
        real(kind_phys), pointer :: dlwsfc_cpl(:)
        real(kind_phys), pointer :: dswsfc_cpl(:)
        real(kind_phys), pointer :: dnirbmi_cpl(:)
        real(kind_phys), pointer :: dnirdfi_cpl(:)
        real(kind_phys), pointer :: dvisbmi_cpl(:)
        real(kind_phys), pointer :: dvisdfi_cpl(:)
        real(kind_phys), pointer :: dnirbm_cpl(:)
        real(kind_phys), pointer :: dnirdf_cpl(:)
        real(kind_phys), pointer :: dvisbm_cpl(:)
        real(kind_phys), pointer :: dvisdf_cpl(:)
        real(kind_phys), pointer :: nlwsfci_cpl(:)
        real(kind_phys), pointer :: nlwsfc_cpl(:)
        real(kind_phys), pointer :: t2mi_cpl(:)
        real(kind_phys), pointer :: q2mi_cpl(:)
        real(kind_phys), pointer :: u10mi_cpl(:)
        real(kind_phys), pointer :: v10mi_cpl(:)
        real(kind_phys), pointer :: tsfci_cpl(:)
        real(kind_phys), pointer :: psurfi_cpl(:)
        real(kind_phys), pointer :: nnirbmi_cpl(:)
        real(kind_phys), pointer :: nnirdfi_cpl(:)
        real(kind_phys), pointer :: nvisbmi_cpl(:)
        real(kind_phys), pointer :: nvisdfi_cpl(:)
        real(kind_phys), pointer :: nswsfci_cpl(:)
        real(kind_phys), pointer :: nswsfc_cpl(:)
        real(kind_phys), pointer :: nnirbm_cpl(:)
        real(kind_phys), pointer :: nnirdf_cpl(:)
        real(kind_phys), pointer :: nvisbm_cpl(:)
        real(kind_phys), pointer :: nvisdf_cpl(:)
        real(kind_phys), pointer :: gflux(:)
        real(kind_phys), pointer :: evbsa(:)
        real(kind_phys), pointer :: evcwa(:)
        real(kind_phys), pointer :: transa(:)
        real(kind_phys), pointer :: sbsnoa(:)
        real(kind_phys), pointer :: snowca(:)
        real(kind_phys), pointer :: snohfa(:)
        real(kind_phys), pointer :: ep(:)
        real(kind_phys), pointer :: runoff(:)
        real(kind_phys), pointer :: srunoff(:)
        real(kind_phys), pointer :: runof(:)
        real(kind_phys), pointer :: drain(:)

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
        

        call ccpp_field_get(cdata, 'surface_upward_potential_latent_heat_flux', ep1d, ierr=ierr, dims=cdims, kind=ckind, index=732)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_upward_potential_latent_heat_flux from CCPP data structure')
            return
        end if
        if (kind(ep1d).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_upward_potential_latent_heat_flux')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'upward_heat_flux_in_soil', gflx, ierr=ierr, dims=cdims, kind=ckind, index=812)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve upward_heat_flux_in_soil from CCPP data structure')
            return
        end if
        if (kind(gflx).ne.ckind) then
            call ccpp_error('Kind mismatch for variable upward_heat_flux_in_soil')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'air_temperature_at_lowest_model_layer', tgrs_1, ierr=ierr, dims=cdims, kind=ckind, index=53)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve air_temperature_at_lowest_model_layer from CCPP data structure')
            return
        end if
        if (kind(tgrs_1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable air_temperature_at_lowest_model_layer')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'water_vapor_specific_humidity_at_lowest_model_layer', qgrs_1, ierr=ierr, dims=cdims, kind=ckind, index=854)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve water_vapor_specific_humidity_at_lowest_model_layer from CCPP data structure')
            return
        end if
        if (kind(qgrs_1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable water_vapor_specific_humidity_at_lowest_model_layer')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'x_wind_at_lowest_model_layer', ugrs_1, ierr=ierr, dims=cdims, kind=ckind, index=873)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve x_wind_at_lowest_model_layer from CCPP data structure')
            return
        end if
        if (kind(ugrs_1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable x_wind_at_lowest_model_layer')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'y_wind_at_lowest_model_layer', vgrs_1, ierr=ierr, dims=cdims, kind=ckind, index=880)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve y_wind_at_lowest_model_layer from CCPP data structure')
            return
        end if
        if (kind(vgrs_1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable y_wind_at_lowest_model_layer')
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
        

        call ccpp_field_get(cdata, 'temperature_at_2m', t2m, ierr=ierr, dims=cdims, kind=ckind, index=750)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve temperature_at_2m from CCPP data structure')
            return
        end if
        if (kind(t2m).ne.ckind) then
            call ccpp_error('Kind mismatch for variable temperature_at_2m')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'specific_humidity_at_2m', q2m, ierr=ierr, dims=cdims, kind=ckind, index=667)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve specific_humidity_at_2m from CCPP data structure')
            return
        end if
        if (kind(q2m).ne.ckind) then
            call ccpp_error('Kind mismatch for variable specific_humidity_at_2m')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'x_wind_at_10m', u10m, ierr=ierr, dims=cdims, kind=ckind, index=872)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve x_wind_at_10m from CCPP data structure')
            return
        end if
        if (kind(u10m).ne.ckind) then
            call ccpp_error('Kind mismatch for variable x_wind_at_10m')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'y_wind_at_10m', v10m, ierr=ierr, dims=cdims, kind=ckind, index=879)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve y_wind_at_10m from CCPP data structure')
            return
        end if
        if (kind(v10m).ne.ckind) then
            call ccpp_error('Kind mismatch for variable y_wind_at_10m')
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
        

        call ccpp_field_get(cdata, 'surface_air_pressure', pgr, ierr=ierr, dims=cdims, kind=ckind, index=677)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_air_pressure from CCPP data structure')
            return
        end if
        if (kind(pgr).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_air_pressure')
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
        

        call ccpp_field_get(cdata, 'instantaneous_surface_potential_evaporation', epi, ierr=ierr, dims=cdims, kind=ckind, index=429)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve instantaneous_surface_potential_evaporation from CCPP data structure')
            return
        end if
        if (kind(epi).ne.ckind) then
            call ccpp_error('Kind mismatch for variable instantaneous_surface_potential_evaporation')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'instantaneous_surface_ground_heat_flux', gfluxi, ierr=ierr, dims=cdims, kind=ckind, index=422)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve instantaneous_surface_ground_heat_flux from CCPP data structure')
            return
        end if
        if (kind(gfluxi).ne.ckind) then
            call ccpp_error('Kind mismatch for variable instantaneous_surface_ground_heat_flux')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'air_temperature_at_lowest_model_layer_for_diag', t1, ierr=ierr, dims=cdims, kind=ckind, index=54)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve air_temperature_at_lowest_model_layer_for_diag from CCPP data structure')
            return
        end if
        if (kind(t1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable air_temperature_at_lowest_model_layer_for_diag')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'water_vapor_specific_humidity_at_lowest_model_layer_for_diag', q1, ierr=ierr, dims=cdims, kind=ckind, index=855)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve water_vapor_specific_humidity_at_lowest_model_layer_for_diag from CCPP data structure')
            return
        end if
        if (kind(q1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable water_vapor_specific_humidity_at_lowest_model_layer_for_diag')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'x_wind_at_lowest_model_layer_for_diag', u1, ierr=ierr, dims=cdims, kind=ckind, index=874)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve x_wind_at_lowest_model_layer_for_diag from CCPP data structure')
            return
        end if
        if (kind(u1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable x_wind_at_lowest_model_layer_for_diag')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'y_wind_at_lowest_model_layer_for_diag', v1, ierr=ierr, dims=cdims, kind=ckind, index=881)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve y_wind_at_lowest_model_layer_for_diag from CCPP data structure')
            return
        end if
        if (kind(v1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable y_wind_at_lowest_model_layer_for_diag')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'instantaneous_surface_downwelling_longwave_flux_for_coupling', dlwsfci_cpl, ierr=ierr, dims=cdims, kind=ckind, index=420)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve instantaneous_surface_downwelling_longwave_flux_for_coupling from CCPP data structure')
            return
        end if
        if (kind(dlwsfci_cpl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable instantaneous_surface_downwelling_longwave_flux_for_coupling')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'instantaneous_surface_downwelling_shortwave_flux_for_coupling', dswsfci_cpl, ierr=ierr, dims=cdims, kind=ckind, index=421)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve instantaneous_surface_downwelling_shortwave_flux_for_coupling from CCPP data structure')
            return
        end if
        if (kind(dswsfci_cpl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable instantaneous_surface_downwelling_shortwave_flux_for_coupling')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cumulative_surface_downwelling_longwave_flux_for_coupling_multiplied_by_timestep', dlwsfc_cpl, ierr=ierr, dims=cdims, kind=ckind, index=183)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cumulative_surface_downwelling_longwave_flux_for_coupling_multiplied_by_timestep from CCPP data structure')
            return
        end if
        if (kind(dlwsfc_cpl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cumulative_surface_downwelling_longwave_flux_for_coupling_multiplied_by_timestep')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cumulative_surface_downwelling_shortwave_flux_for_coupling_multiplied_by_timestep', dswsfc_cpl, ierr=ierr, dims=cdims, kind=ckind, index=185)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cumulative_surface_downwelling_shortwave_flux_for_coupling_multiplied_by_timestep from CCPP data structure')
            return
        end if
        if (kind(dswsfc_cpl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cumulative_surface_downwelling_shortwave_flux_for_coupling_multiplied_by_timestep')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'instantaneous_surface_downwelling_direct_near_infrared_shortwave_flux_for_coupling', dnirbmi_cpl, ierr=ierr, dims=cdims, kind=ckind, index=418)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve instantaneous_surface_downwelling_direct_near_infrared_shortwave_flux_for_coupling from CCPP data structure')
            return
        end if
        if (kind(dnirbmi_cpl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable instantaneous_surface_downwelling_direct_near_infrared_shortwave_flux_for_coupling')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'instantaneous_surface_downwelling_diffuse_near_infrared_shortwave_flux_for_coupling', dnirdfi_cpl, ierr=ierr, dims=cdims, kind=ckind, index=416)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve instantaneous_surface_downwelling_diffuse_near_infrared_shortwave_flux_for_coupling from CCPP data structure')
            return
        end if
        if (kind(dnirdfi_cpl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable instantaneous_surface_downwelling_diffuse_near_infrared_shortwave_flux_for_coupling')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'instantaneous_surface_downwelling_direct_ultraviolet_and_visible_shortwave_flux_for_coupling', dvisbmi_cpl, ierr=ierr, dims=cdims, kind=ckind, index=419)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve instantaneous_surface_downwelling_direct_ultraviolet_and_visible_shortwave_flux_for_coupling from CCPP data structure')
            return
        end if
        if (kind(dvisbmi_cpl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable instantaneous_surface_downwelling_direct_ultraviolet_and_visible_shortwave_flux_for_coupling')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'instantaneous_surface_downwelling_diffuse_ultraviolet_and_visible_shortwave_flux_for_coupling', dvisdfi_cpl, ierr=ierr, dims=cdims, kind=ckind, index=417)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve instantaneous_surface_downwelling_diffuse_ultraviolet_and_visible_shortwave_flux_for_coupling from CCPP data structure')
            return
        end if
        if (kind(dvisdfi_cpl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable instantaneous_surface_downwelling_diffuse_ultraviolet_and_visible_shortwave_flux_for_coupling')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cumulative_surface_downwelling_direct_near_infrared_shortwave_flux_for_coupling_multiplied_by_timestep', dnirbm_cpl, ierr=ierr, dims=cdims, kind=ckind, index=181)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cumulative_surface_downwelling_direct_near_infrared_shortwave_flux_for_coupling_multiplied_by_timestep from CCPP data structure')
            return
        end if
        if (kind(dnirbm_cpl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cumulative_surface_downwelling_direct_near_infrared_shortwave_flux_for_coupling_multiplied_by_timestep')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cumulative_surface_downwelling_diffuse_near_infrared_shortwave_flux_for_coupling_multiplied_by_timestep', dnirdf_cpl, ierr=ierr, dims=cdims, kind=ckind, index=179)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cumulative_surface_downwelling_diffuse_near_infrared_shortwave_flux_for_coupling_multiplied_by_timestep from CCPP data structure')
            return
        end if
        if (kind(dnirdf_cpl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cumulative_surface_downwelling_diffuse_near_infrared_shortwave_flux_for_coupling_multiplied_by_timestep')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cumulative_surface_downwelling_direct_ultraviolet_and_visible_shortwave_flux_for_coupling_multiplied_by_timestep', dvisbm_cpl, ierr=ierr, dims=cdims, kind=ckind, index=182)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cumulative_surface_downwelling_direct_ultraviolet_and_visible_shortwave_flux_for_coupling_multiplied_by_timestep from CCPP data structure')
            return
        end if
        if (kind(dvisbm_cpl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cumulative_surface_downwelling_direct_ultraviolet_and_visible_shortwave_flux_for_coupling_multiplied_by_timestep')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cumulative_surface_downwelling_diffuse_ultraviolet_and_visible_shortwave_flux_for_coupling_multiplied_by_timestep', dvisdf_cpl, ierr=ierr, dims=cdims, kind=ckind, index=180)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cumulative_surface_downwelling_diffuse_ultraviolet_and_visible_shortwave_flux_for_coupling_multiplied_by_timestep from CCPP data structure')
            return
        end if
        if (kind(dvisdf_cpl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cumulative_surface_downwelling_diffuse_ultraviolet_and_visible_shortwave_flux_for_coupling_multiplied_by_timestep')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'instantaneous_surface_net_downward_longwave_flux_for_coupling', nlwsfci_cpl, ierr=ierr, dims=cdims, kind=ckind, index=427)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve instantaneous_surface_net_downward_longwave_flux_for_coupling from CCPP data structure')
            return
        end if
        if (kind(nlwsfci_cpl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable instantaneous_surface_net_downward_longwave_flux_for_coupling')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cumulative_surface_net_downward_longwave_flux_for_coupling_multiplied_by_timestep', nlwsfc_cpl, ierr=ierr, dims=cdims, kind=ckind, index=191)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cumulative_surface_net_downward_longwave_flux_for_coupling_multiplied_by_timestep from CCPP data structure')
            return
        end if
        if (kind(nlwsfc_cpl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cumulative_surface_net_downward_longwave_flux_for_coupling_multiplied_by_timestep')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'instantaneous_temperature_at_2m_for_coupling', t2mi_cpl, ierr=ierr, dims=cdims, kind=ckind, index=443)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve instantaneous_temperature_at_2m_for_coupling from CCPP data structure')
            return
        end if
        if (kind(t2mi_cpl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable instantaneous_temperature_at_2m_for_coupling')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'instantaneous_specific_humidity_at_2m_for_coupling', q2mi_cpl, ierr=ierr, dims=cdims, kind=ckind, index=414)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve instantaneous_specific_humidity_at_2m_for_coupling from CCPP data structure')
            return
        end if
        if (kind(q2mi_cpl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable instantaneous_specific_humidity_at_2m_for_coupling')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'instantaneous_x_wind_at_10m_for_coupling', u10mi_cpl, ierr=ierr, dims=cdims, kind=ckind, index=446)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve instantaneous_x_wind_at_10m_for_coupling from CCPP data structure')
            return
        end if
        if (kind(u10mi_cpl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable instantaneous_x_wind_at_10m_for_coupling')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'instantaneous_y_wind_at_10m_for_coupling', v10mi_cpl, ierr=ierr, dims=cdims, kind=ckind, index=448)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve instantaneous_y_wind_at_10m_for_coupling from CCPP data structure')
            return
        end if
        if (kind(v10mi_cpl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable instantaneous_y_wind_at_10m_for_coupling')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'instantaneous_surface_skin_temperature_for_coupling', tsfci_cpl, ierr=ierr, dims=cdims, kind=ckind, index=430)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve instantaneous_surface_skin_temperature_for_coupling from CCPP data structure')
            return
        end if
        if (kind(tsfci_cpl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable instantaneous_surface_skin_temperature_for_coupling')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'instantaneous_surface_air_pressure_for_coupling', psurfi_cpl, ierr=ierr, dims=cdims, kind=ckind, index=415)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve instantaneous_surface_air_pressure_for_coupling from CCPP data structure')
            return
        end if
        if (kind(psurfi_cpl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable instantaneous_surface_air_pressure_for_coupling')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'instantaneous_surface_net_downward_direct_near_infrared_shortwave_flux_for_coupling', nnirbmi_cpl, ierr=ierr, dims=cdims, kind=ckind, index=425)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve instantaneous_surface_net_downward_direct_near_infrared_shortwave_flux_for_coupling from CCPP data structure')
            return
        end if
        if (kind(nnirbmi_cpl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable instantaneous_surface_net_downward_direct_near_infrared_shortwave_flux_for_coupling')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'instantaneous_surface_net_downward_diffuse_near_infrared_shortwave_flux_for_coupling', nnirdfi_cpl, ierr=ierr, dims=cdims, kind=ckind, index=423)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve instantaneous_surface_net_downward_diffuse_near_infrared_shortwave_flux_for_coupling from CCPP data structure')
            return
        end if
        if (kind(nnirdfi_cpl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable instantaneous_surface_net_downward_diffuse_near_infrared_shortwave_flux_for_coupling')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'instantaneous_surface_net_downward_direct_ultraviolet_and_visible_shortwave_flux_for_coupling', nvisbmi_cpl, ierr=ierr, dims=cdims, kind=ckind, index=426)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve instantaneous_surface_net_downward_direct_ultraviolet_and_visible_shortwave_flux_for_coupling from CCPP data structure')
            return
        end if
        if (kind(nvisbmi_cpl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable instantaneous_surface_net_downward_direct_ultraviolet_and_visible_shortwave_flux_for_coupling')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'instantaneous_surface_net_downward_diffuse_ultraviolet_and_visible_shortwave_flux_for_coupling', nvisdfi_cpl, ierr=ierr, dims=cdims, kind=ckind, index=424)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve instantaneous_surface_net_downward_diffuse_ultraviolet_and_visible_shortwave_flux_for_coupling from CCPP data structure')
            return
        end if
        if (kind(nvisdfi_cpl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable instantaneous_surface_net_downward_diffuse_ultraviolet_and_visible_shortwave_flux_for_coupling')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'instantaneous_surface_net_downward_shortwave_flux_for_coupling', nswsfci_cpl, ierr=ierr, dims=cdims, kind=ckind, index=428)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve instantaneous_surface_net_downward_shortwave_flux_for_coupling from CCPP data structure')
            return
        end if
        if (kind(nswsfci_cpl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable instantaneous_surface_net_downward_shortwave_flux_for_coupling')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cumulative_surface_net_downward_shortwave_flux_for_coupling_multiplied_by_timestep', nswsfc_cpl, ierr=ierr, dims=cdims, kind=ckind, index=192)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cumulative_surface_net_downward_shortwave_flux_for_coupling_multiplied_by_timestep from CCPP data structure')
            return
        end if
        if (kind(nswsfc_cpl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cumulative_surface_net_downward_shortwave_flux_for_coupling_multiplied_by_timestep')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cumulative_surface_net_downward_direct_near_infrared_shortwave_flux_for_coupling_multiplied_by_timestep', nnirbm_cpl, ierr=ierr, dims=cdims, kind=ckind, index=189)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cumulative_surface_net_downward_direct_near_infrared_shortwave_flux_for_coupling_multiplied_by_timestep from CCPP data structure')
            return
        end if
        if (kind(nnirbm_cpl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cumulative_surface_net_downward_direct_near_infrared_shortwave_flux_for_coupling_multiplied_by_timestep')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cumulative_surface_net_downward_diffuse_near_infrared_shortwave_flux_for_coupling_multiplied_by_timestep', nnirdf_cpl, ierr=ierr, dims=cdims, kind=ckind, index=187)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cumulative_surface_net_downward_diffuse_near_infrared_shortwave_flux_for_coupling_multiplied_by_timestep from CCPP data structure')
            return
        end if
        if (kind(nnirdf_cpl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cumulative_surface_net_downward_diffuse_near_infrared_shortwave_flux_for_coupling_multiplied_by_timestep')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cumulative_surface_net_downward_direct_ultraviolet_and_visible_shortwave_flux_for_coupling_multiplied_by_timestep', nvisbm_cpl, ierr=ierr, dims=cdims, kind=ckind, index=190)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cumulative_surface_net_downward_direct_ultraviolet_and_visible_shortwave_flux_for_coupling_multiplied_by_timestep from CCPP data structure')
            return
        end if
        if (kind(nvisbm_cpl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cumulative_surface_net_downward_direct_ultraviolet_and_visible_shortwave_flux_for_coupling_multiplied_by_timestep')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cumulative_surface_net_downward_diffuse_ultraviolet_and_visible_shortwave_flux_for_coupling_multiplied_by_timestep', nvisdf_cpl, ierr=ierr, dims=cdims, kind=ckind, index=188)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cumulative_surface_net_downward_diffuse_ultraviolet_and_visible_shortwave_flux_for_coupling_multiplied_by_timestep from CCPP data structure')
            return
        end if
        if (kind(nvisdf_cpl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cumulative_surface_net_downward_diffuse_ultraviolet_and_visible_shortwave_flux_for_coupling_multiplied_by_timestep')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cumulative_surface_ground_heat_flux_multiplied_by_timestep', gflux, ierr=ierr, dims=cdims, kind=ckind, index=186)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cumulative_surface_ground_heat_flux_multiplied_by_timestep from CCPP data structure')
            return
        end if
        if (kind(gflux).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cumulative_surface_ground_heat_flux_multiplied_by_timestep')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cumulative_soil_upward_latent_heat_flux_multiplied_by_timestep', evbsa, ierr=ierr, dims=cdims, kind=ckind, index=178)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cumulative_soil_upward_latent_heat_flux_multiplied_by_timestep from CCPP data structure')
            return
        end if
        if (kind(evbsa).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cumulative_soil_upward_latent_heat_flux_multiplied_by_timestep')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cumulative_canopy_upward_latent_heat_flu_multiplied_by_timestep', evcwa, ierr=ierr, dims=cdims, kind=ckind, index=149)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cumulative_canopy_upward_latent_heat_flu_multiplied_by_timestep from CCPP data structure')
            return
        end if
        if (kind(evcwa).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cumulative_canopy_upward_latent_heat_flu_multiplied_by_timestep')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cumulative_transpiration_flux_multiplied_by_timestep', transa, ierr=ierr, dims=cdims, kind=ckind, index=205)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cumulative_transpiration_flux_multiplied_by_timestep from CCPP data structure')
            return
        end if
        if (kind(transa).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cumulative_transpiration_flux_multiplied_by_timestep')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cumulative_snow_deposition_sublimation_upward_latent_heat_flux_multiplied_by_timestep', sbsnoa, ierr=ierr, dims=cdims, kind=ckind, index=176)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cumulative_snow_deposition_sublimation_upward_latent_heat_flux_multiplied_by_timestep from CCPP data structure')
            return
        end if
        if (kind(sbsnoa).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cumulative_snow_deposition_sublimation_upward_latent_heat_flux_multiplied_by_timestep')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cumulative_surface_snow_area_fraction_multiplied_by_timestep', snowca, ierr=ierr, dims=cdims, kind=ckind, index=194)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cumulative_surface_snow_area_fraction_multiplied_by_timestep from CCPP data structure')
            return
        end if
        if (kind(snowca).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cumulative_surface_snow_area_fraction_multiplied_by_timestep')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cumulative_snow_freezing_rain_upward_latent_heat_flux_multiplied_by_timestep', snohfa, ierr=ierr, dims=cdims, kind=ckind, index=177)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cumulative_snow_freezing_rain_upward_latent_heat_flux_multiplied_by_timestep from CCPP data structure')
            return
        end if
        if (kind(snohfa).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cumulative_snow_freezing_rain_upward_latent_heat_flux_multiplied_by_timestep')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cumulative_surface_upward_potential_latent_heat_flux_multiplied_by_timestep', ep, ierr=ierr, dims=cdims, kind=ckind, index=197)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cumulative_surface_upward_potential_latent_heat_flux_multiplied_by_timestep from CCPP data structure')
            return
        end if
        if (kind(ep).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cumulative_surface_upward_potential_latent_heat_flux_multiplied_by_timestep')
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
        

        call GFS_surface_generic_post_run(im=im,cplflx=cplflx,lssav=lssav,islmsk=islmsk,dtf=dtf,ep1d=ep1d,gflx=gflx, &
                  tgrs_1=tgrs_1,qgrs_1=qgrs_1,ugrs_1=ugrs_1,vgrs_1=vgrs_1,adjsfcdlw=adjsfcdlw, &
                  adjsfcdsw=adjsfcdsw,adjnirbmd=adjnirbmd,adjnirdfd=adjnirdfd,adjvisbmd=adjvisbmd, &
                  adjvisdfd=adjvisdfd,adjsfculw=adjsfculw,adjnirbmu=adjnirbmu,adjnirdfu=adjnirdfu, &
                  adjvisbmu=adjvisbmu,adjvisdfu=adjvisdfu,t2m=t2m,q2m=q2m,u10m=u10m,v10m=v10m, &
                  tsfc=tsfc,pgr=pgr,xcosz=xcosz,evbs=evbs,evcw=evcw,trans=trans,sbsno=sbsno, &
                  snowc=snowc,snohf=snohf,epi=epi,gfluxi=gfluxi,t1=t1,q1=q1,u1=u1,v1=v1,dlwsfci_cpl=dlwsfci_cpl, &
                  dswsfci_cpl=dswsfci_cpl,dlwsfc_cpl=dlwsfc_cpl,dswsfc_cpl=dswsfc_cpl,dnirbmi_cpl=dnirbmi_cpl, &
                  dnirdfi_cpl=dnirdfi_cpl,dvisbmi_cpl=dvisbmi_cpl,dvisdfi_cpl=dvisdfi_cpl, &
                  dnirbm_cpl=dnirbm_cpl,dnirdf_cpl=dnirdf_cpl,dvisbm_cpl=dvisbm_cpl,dvisdf_cpl=dvisdf_cpl, &
                  nlwsfci_cpl=nlwsfci_cpl,nlwsfc_cpl=nlwsfc_cpl,t2mi_cpl=t2mi_cpl,q2mi_cpl=q2mi_cpl, &
                  u10mi_cpl=u10mi_cpl,v10mi_cpl=v10mi_cpl,tsfci_cpl=tsfci_cpl,psurfi_cpl=psurfi_cpl, &
                  nnirbmi_cpl=nnirbmi_cpl,nnirdfi_cpl=nnirdfi_cpl,nvisbmi_cpl=nvisbmi_cpl, &
                  nvisdfi_cpl=nvisdfi_cpl,nswsfci_cpl=nswsfci_cpl,nswsfc_cpl=nswsfc_cpl,nnirbm_cpl=nnirbm_cpl, &
                  nnirdf_cpl=nnirdf_cpl,nvisbm_cpl=nvisbm_cpl,nvisdf_cpl=nvisdf_cpl,gflux=gflux, &
                  evbsa=evbsa,evcwa=evcwa,transa=transa,sbsnoa=sbsnoa,snowca=snowca,snohfa=snohfa, &
                  ep=ep,runoff=runoff,srunoff=srunoff,runof=runof,drain=drain,errmsg=cdata%errmsg, &
                  errflg=cdata%errflg)
        ierr=cdata%errflg

    end function GFS_surface_generic_post_run_cap
end module GFS_surface_generic_post_cap

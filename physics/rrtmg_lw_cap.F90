
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
!! @brief Auto-generated cap module for the rrtmg_lw scheme
!!
!
module rrtmg_lw_cap

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_f_pointer, c_ptr, c_int32_t
    use            :: ccpp_types,                                      &
                      only: ccpp_t, CCPP_GENERIC_KIND
    use            :: ccpp_fields,                                     &
                      only: ccpp_field_get
    use            :: ccpp_errors,                                     &
                      only: ccpp_error, ccpp_debug
    use            :: rrtmg_lw, &
                      only: rrtmg_lw_init,rrtmg_lw_finalize,rrtmg_lw_run
    ! Other modules required, e.g. type definitions
    use machine, only: kind_phys
    use GFS_typedefs, only: topflw_type
    use GFS_typedefs, only: sfcflw_type

    implicit none

    private
    public :: rrtmg_lw_init_cap,rrtmg_lw_finalize_cap,rrtmg_lw_run_cap

    contains


    function rrtmg_lw_init_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call rrtmg_lw_init()
        

    end function rrtmg_lw_init_cap

    function rrtmg_lw_finalize_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call rrtmg_lw_finalize()
        

    end function rrtmg_lw_finalize_cap

    function rrtmg_lw_run_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind
        real(kind_phys), pointer :: plyr(:,:)
        real(kind_phys), pointer :: plvl(:,:)
        real(kind_phys), pointer :: tlyr(:,:)
        real(kind_phys), pointer :: tlvl(:,:)
        real(kind_phys), pointer :: qlyr(:,:)
        real(kind_phys), pointer :: olyr(:,:)
        real(kind_phys), pointer :: gasvmr_co2(:,:)
        real(kind_phys), pointer :: gasvmr_n2o(:,:)
        real(kind_phys), pointer :: gasvmr_ch4(:,:)
        real(kind_phys), pointer :: gasvmr_o2(:,:)
        real(kind_phys), pointer :: gasvmr_co(:,:)
        real(kind_phys), pointer :: gasvmr_cfc11(:,:)
        real(kind_phys), pointer :: gasvmr_cfc12(:,:)
        real(kind_phys), pointer :: gasvmr_cfc22(:,:)
        real(kind_phys), pointer :: gasvmr_ccl4(:,:)
        integer, pointer :: icseed(:)
        real(kind_phys), pointer :: aeraod(:,:,:)
        real(kind_phys), pointer :: aerssa(:,:,:)
        real(kind_phys), pointer :: sfemis(:)
        real(kind_phys), pointer :: sfgtmp(:)
        real(kind_phys), pointer :: dzlyr(:,:)
        real(kind_phys), pointer :: delpin(:,:)
        real(kind_phys), pointer :: de_lgth(:)
        integer, pointer :: npts
        integer, pointer :: nlay
        integer, pointer :: nlp1
        logical, pointer :: lprnt
        real(kind_phys), pointer :: cld_cf(:,:)
        logical, pointer :: lslwr
        real(kind_phys), pointer :: hlwc(:,:)
        type(topflw_type), pointer     :: topflx(:)
        type(sfcflw_type), pointer     :: sfcflx(:)
        real(kind_phys), pointer :: cldtau(:,:)
        real(kind_phys), pointer :: hlw0(:,:)
        real(kind_phys), pointer :: cld_lwp(:,:)
        real(kind_phys), pointer :: cld_ref_liq(:,:)
        real(kind_phys), pointer :: cld_iwp(:,:)
        real(kind_phys), pointer :: cld_ref_ice(:,:)
        real(kind_phys), pointer :: cld_rwp(:,:)
        real(kind_phys), pointer :: cld_ref_rain(:,:)
        real(kind_phys), pointer :: cld_swp(:,:)
        real(kind_phys), pointer :: cld_ref_snow(:,:)

        ierr = 0

        call c_f_pointer(ptr, cdata)


        call ccpp_field_get(cdata, 'air_pressure_at_layer_for_radiation_in_hPa', plyr, ierr=ierr, dims=cdims, kind=ckind, index=47)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve air_pressure_at_layer_for_radiation_in_hPa from CCPP data structure')
            return
        end if
        if (kind(plyr).ne.ckind) then
            call ccpp_error('Kind mismatch for variable air_pressure_at_layer_for_radiation_in_hPa')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'air_pressure_at_interface_for_radiation_in_hPa', plvl, ierr=ierr, dims=cdims, kind=ckind, index=46)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve air_pressure_at_interface_for_radiation_in_hPa from CCPP data structure')
            return
        end if
        if (kind(plvl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable air_pressure_at_interface_for_radiation_in_hPa')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'air_temperature_at_layer_for_radiation', tlyr, ierr=ierr, dims=cdims, kind=ckind, index=52)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve air_temperature_at_layer_for_radiation from CCPP data structure')
            return
        end if
        if (kind(tlyr).ne.ckind) then
            call ccpp_error('Kind mismatch for variable air_temperature_at_layer_for_radiation')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'air_temperature_at_interface_for_radiation', tlvl, ierr=ierr, dims=cdims, kind=ckind, index=51)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve air_temperature_at_interface_for_radiation from CCPP data structure')
            return
        end if
        if (kind(tlvl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable air_temperature_at_interface_for_radiation')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'water_vapor_specific_humidity_at_layer_for_radiation', qlyr, ierr=ierr, dims=cdims, kind=ckind, index=853)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve water_vapor_specific_humidity_at_layer_for_radiation from CCPP data structure')
            return
        end if
        if (kind(qlyr).ne.ckind) then
            call ccpp_error('Kind mismatch for variable water_vapor_specific_humidity_at_layer_for_radiation')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'ozone_concentration_at_layer_for_radiation', olyr, ierr=ierr, dims=cdims, kind=ckind, index=597)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve ozone_concentration_at_layer_for_radiation from CCPP data structure')
            return
        end if
        if (kind(olyr).ne.ckind) then
            call ccpp_error('Kind mismatch for variable ozone_concentration_at_layer_for_radiation')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'volume_mixing_ratio_co2', gasvmr_co2, ierr=ierr, dims=cdims, kind=ckind, index=845)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve volume_mixing_ratio_co2 from CCPP data structure')
            return
        end if
        if (kind(gasvmr_co2).ne.ckind) then
            call ccpp_error('Kind mismatch for variable volume_mixing_ratio_co2')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'volume_mixing_ratio_n2o', gasvmr_n2o, ierr=ierr, dims=cdims, kind=ckind, index=846)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve volume_mixing_ratio_n2o from CCPP data structure')
            return
        end if
        if (kind(gasvmr_n2o).ne.ckind) then
            call ccpp_error('Kind mismatch for variable volume_mixing_ratio_n2o')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'volume_mixing_ratio_ch4', gasvmr_ch4, ierr=ierr, dims=cdims, kind=ckind, index=843)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve volume_mixing_ratio_ch4 from CCPP data structure')
            return
        end if
        if (kind(gasvmr_ch4).ne.ckind) then
            call ccpp_error('Kind mismatch for variable volume_mixing_ratio_ch4')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'volume_mixing_ratio_o2', gasvmr_o2, ierr=ierr, dims=cdims, kind=ckind, index=847)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve volume_mixing_ratio_o2 from CCPP data structure')
            return
        end if
        if (kind(gasvmr_o2).ne.ckind) then
            call ccpp_error('Kind mismatch for variable volume_mixing_ratio_o2')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'volume_mixing_ratio_co', gasvmr_co, ierr=ierr, dims=cdims, kind=ckind, index=844)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve volume_mixing_ratio_co from CCPP data structure')
            return
        end if
        if (kind(gasvmr_co).ne.ckind) then
            call ccpp_error('Kind mismatch for variable volume_mixing_ratio_co')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'volume_mixing_ratio_cfc11', gasvmr_cfc11, ierr=ierr, dims=cdims, kind=ckind, index=839)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve volume_mixing_ratio_cfc11 from CCPP data structure')
            return
        end if
        if (kind(gasvmr_cfc11).ne.ckind) then
            call ccpp_error('Kind mismatch for variable volume_mixing_ratio_cfc11')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'volume_mixing_ratio_cfc12', gasvmr_cfc12, ierr=ierr, dims=cdims, kind=ckind, index=841)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve volume_mixing_ratio_cfc12 from CCPP data structure')
            return
        end if
        if (kind(gasvmr_cfc12).ne.ckind) then
            call ccpp_error('Kind mismatch for variable volume_mixing_ratio_cfc12')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'volume_mixing_ratio_cfc22', gasvmr_cfc22, ierr=ierr, dims=cdims, kind=ckind, index=842)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve volume_mixing_ratio_cfc22 from CCPP data structure')
            return
        end if
        if (kind(gasvmr_cfc22).ne.ckind) then
            call ccpp_error('Kind mismatch for variable volume_mixing_ratio_cfc22')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'volume_mixing_ratio_ccl4', gasvmr_ccl4, ierr=ierr, dims=cdims, kind=ckind, index=838)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve volume_mixing_ratio_ccl4 from CCPP data structure')
            return
        end if
        if (kind(gasvmr_ccl4).ne.ckind) then
            call ccpp_error('Kind mismatch for variable volume_mixing_ratio_ccl4')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'seed_random_numbers_lw', icseed, ierr=ierr, dims=cdims, kind=ckind, index=635)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve seed_random_numbers_lw from CCPP data structure')
            return
        end if
        if (kind(icseed).ne.ckind) then
            call ccpp_error('Kind mismatch for variable seed_random_numbers_lw')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'aerosol_optical_depth_for_longwave_bands_01-16', aeraod, ierr=ierr, dims=cdims, kind=ckind, index=38)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve aerosol_optical_depth_for_longwave_bands_01-16 from CCPP data structure')
            return
        end if
        if (kind(aeraod).ne.ckind) then
            call ccpp_error('Kind mismatch for variable aerosol_optical_depth_for_longwave_bands_01-16')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'aerosol_single_scattering_albedo_for_longwave_bands_01-16', aerssa, ierr=ierr, dims=cdims, kind=ckind, index=42)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve aerosol_single_scattering_albedo_for_longwave_bands_01-16 from CCPP data structure')
            return
        end if
        if (kind(aerssa).ne.ckind) then
            call ccpp_error('Kind mismatch for variable aerosol_single_scattering_albedo_for_longwave_bands_01-16')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_longwave_emissivity', sfemis, ierr=ierr, dims=cdims, kind=ckind, index=714)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_longwave_emissivity from CCPP data structure')
            return
        end if
        if (kind(sfemis).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_longwave_emissivity')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_ground_temperature_for_radiation', sfgtmp, ierr=ierr, dims=cdims, kind=ckind, index=712)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_ground_temperature_for_radiation from CCPP data structure')
            return
        end if
        if (kind(sfgtmp).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_ground_temperature_for_radiation')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'layer_thickness_for_radiation', dzlyr, ierr=ierr, dims=cdims, kind=ckind, index=464)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve layer_thickness_for_radiation from CCPP data structure')
            return
        end if
        if (kind(dzlyr).ne.ckind) then
            call ccpp_error('Kind mismatch for variable layer_thickness_for_radiation')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'layer_pressure_thickness_for_radiation', delpin, ierr=ierr, dims=cdims, kind=ckind, index=463)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve layer_pressure_thickness_for_radiation from CCPP data structure')
            return
        end if
        if (kind(delpin).ne.ckind) then
            call ccpp_error('Kind mismatch for variable layer_pressure_thickness_for_radiation')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cloud_decorrelation_length', de_lgth, ierr=ierr, dims=cdims, kind=ckind, index=96)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cloud_decorrelation_length from CCPP data structure')
            return
        end if
        if (kind(de_lgth).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cloud_decorrelation_length')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'horizontal_loop_extent', npts, ierr=ierr, kind=ckind, index=366)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve horizontal_loop_extent from CCPP data structure')
            return
        end if
        if (kind(npts).ne.ckind) then
            call ccpp_error('Kind mismatch for variable horizontal_loop_extent')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'adjusted_vertical_layer_dimension_for_radiation', nlay, ierr=ierr, kind=ckind, index=31)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve adjusted_vertical_layer_dimension_for_radiation from CCPP data structure')
            return
        end if
        if (kind(nlay).ne.ckind) then
            call ccpp_error('Kind mismatch for variable adjusted_vertical_layer_dimension_for_radiation')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'adjusted_vertical_level_dimension_for_radiation', nlp1, ierr=ierr, kind=ckind, index=32)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve adjusted_vertical_level_dimension_for_radiation from CCPP data structure')
            return
        end if
        if (kind(nlp1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable adjusted_vertical_level_dimension_for_radiation')
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
        

        call ccpp_field_get(cdata, 'total_cloud_fraction', cld_cf, ierr=ierr, dims=cdims, kind=ckind, index=799)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve total_cloud_fraction from CCPP data structure')
            return
        end if
        if (kind(cld_cf).ne.ckind) then
            call ccpp_error('Kind mismatch for variable total_cloud_fraction')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'flag_to_calc_lw', lslwr, ierr=ierr, kind=ckind, index=335)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_to_calc_lw from CCPP data structure')
            return
        end if
        if (kind(lslwr).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_to_calc_lw')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'tendency_of_air_temperature_due_to_longwave_heating_on_radiation_time_step', hlwc, ierr=ierr, dims=cdims, kind=ckind, index=756)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve tendency_of_air_temperature_due_to_longwave_heating_on_radiation_time_step from CCPP data structure')
            return
        end if
        if (kind(hlwc).ne.ckind) then
            call ccpp_error('Kind mismatch for variable tendency_of_air_temperature_due_to_longwave_heating_on_radiation_time_step')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'lw_fluxes_top_atmosphere', cptr, ierr=ierr, dims=cdims, kind=ckind, index=475)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve lw_fluxes_top_atmosphere from CCPP data structure')
            return
        end if
        if (ckind.ne.CCPP_GENERIC_KIND) then
            call ccpp_error('Kind mismatch for variable lw_fluxes_top_atmosphere')
            ierr = 1
            return
        end if
#endif
        call c_f_pointer(cptr, topflx, cdims)
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'lw_fluxes_sfc', cptr, ierr=ierr, dims=cdims, kind=ckind, index=474)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve lw_fluxes_sfc from CCPP data structure')
            return
        end if
        if (ckind.ne.CCPP_GENERIC_KIND) then
            call ccpp_error('Kind mismatch for variable lw_fluxes_sfc')
            ierr = 1
            return
        end if
#endif
        call c_f_pointer(cptr, sfcflx, cdims)
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cloud_optical_depth_layers_at_10mu_band', cldtau, ierr=ierr, dims=cdims, kind=ckind, index=104)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cloud_optical_depth_layers_at_10mu_band from CCPP data structure')
            return
        end if
        if (kind(cldtau).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cloud_optical_depth_layers_at_10mu_band')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'tendency_of_air_temperature_due_to_longwave_heating_assuming_clear_sky_on_radiation_time_step', hlw0, ierr=ierr, dims=cdims, kind=ckind, index=754)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve tendency_of_air_temperature_due_to_longwave_heating_assuming_clear_sky_on_radiation_time_step from CCPP data structure')
            return
        end if
        if (kind(hlw0).ne.ckind) then
            call ccpp_error('Kind mismatch for variable tendency_of_air_temperature_due_to_longwave_heating_assuming_clear_sky_on_radiation_time_step')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cloud_liquid_water_path', cld_lwp, ierr=ierr, dims=cdims, kind=ckind, index=102)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cloud_liquid_water_path from CCPP data structure')
            return
        end if
        if (kind(cld_lwp).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cloud_liquid_water_path')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'mean_effective_radius_for_liquid_cloud', cld_ref_liq, ierr=ierr, dims=cdims, kind=ckind, index=516)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve mean_effective_radius_for_liquid_cloud from CCPP data structure')
            return
        end if
        if (kind(cld_ref_liq).ne.ckind) then
            call ccpp_error('Kind mismatch for variable mean_effective_radius_for_liquid_cloud')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cloud_ice_water_path', cld_iwp, ierr=ierr, dims=cdims, kind=ckind, index=101)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cloud_ice_water_path from CCPP data structure')
            return
        end if
        if (kind(cld_iwp).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cloud_ice_water_path')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'mean_effective_radius_for_ice_cloud', cld_ref_ice, ierr=ierr, dims=cdims, kind=ckind, index=515)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve mean_effective_radius_for_ice_cloud from CCPP data structure')
            return
        end if
        if (kind(cld_ref_ice).ne.ckind) then
            call ccpp_error('Kind mismatch for variable mean_effective_radius_for_ice_cloud')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cloud_rain_water_path', cld_rwp, ierr=ierr, dims=cdims, kind=ckind, index=107)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cloud_rain_water_path from CCPP data structure')
            return
        end if
        if (kind(cld_rwp).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cloud_rain_water_path')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'mean_effective_radius_for_rain_drop', cld_ref_rain, ierr=ierr, dims=cdims, kind=ckind, index=517)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve mean_effective_radius_for_rain_drop from CCPP data structure')
            return
        end if
        if (kind(cld_ref_rain).ne.ckind) then
            call ccpp_error('Kind mismatch for variable mean_effective_radius_for_rain_drop')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cloud_snow_water_path', cld_swp, ierr=ierr, dims=cdims, kind=ckind, index=108)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cloud_snow_water_path from CCPP data structure')
            return
        end if
        if (kind(cld_swp).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cloud_snow_water_path')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'mean_effective_radius_for_snow_flake', cld_ref_snow, ierr=ierr, dims=cdims, kind=ckind, index=518)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve mean_effective_radius_for_snow_flake from CCPP data structure')
            return
        end if
        if (kind(cld_ref_snow).ne.ckind) then
            call ccpp_error('Kind mismatch for variable mean_effective_radius_for_snow_flake')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call rrtmg_lw_run(plyr=plyr,plvl=plvl,tlyr=tlyr,tlvl=tlvl,qlyr=qlyr,olyr=olyr,gasvmr_co2=gasvmr_co2, &
                  gasvmr_n2o=gasvmr_n2o,gasvmr_ch4=gasvmr_ch4,gasvmr_o2=gasvmr_o2,gasvmr_co=gasvmr_co, &
                  gasvmr_cfc11=gasvmr_cfc11,gasvmr_cfc12=gasvmr_cfc12,gasvmr_cfc22=gasvmr_cfc22, &
                  gasvmr_ccl4=gasvmr_ccl4,icseed=icseed,aeraod=aeraod,aerssa=aerssa,sfemis=sfemis, &
                  sfgtmp=sfgtmp,dzlyr=dzlyr,delpin=delpin,de_lgth=de_lgth,npts=npts,nlay=nlay, &
                  nlp1=nlp1,lprnt=lprnt,cld_cf=cld_cf,lslwr=lslwr,hlwc=hlwc,topflx=topflx, &
                  sfcflx=sfcflx,cldtau=cldtau,hlw0=hlw0,cld_lwp=cld_lwp,cld_ref_liq=cld_ref_liq, &
                  cld_iwp=cld_iwp,cld_ref_ice=cld_ref_ice,cld_rwp=cld_rwp,cld_ref_rain=cld_ref_rain, &
                  cld_swp=cld_swp,cld_ref_snow=cld_ref_snow,errmsg=cdata%errmsg,errflg=cdata%errflg)
        ierr=cdata%errflg

    end function rrtmg_lw_run_cap
end module rrtmg_lw_cap

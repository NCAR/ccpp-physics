
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
!! @brief Auto-generated cap module for the GFS_rrtmg_pre scheme
!!
!
module GFS_rrtmg_pre_cap

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_f_pointer, c_ptr, c_int32_t
    use            :: ccpp_types,                                      &
                      only: ccpp_t, CCPP_GENERIC_KIND
    use            :: ccpp_fields,                                     &
                      only: ccpp_field_get
    use            :: ccpp_errors,                                     &
                      only: ccpp_error, ccpp_debug
    use            :: GFS_rrtmg_pre, &
                      only: GFS_rrtmg_pre_finalize,GFS_rrtmg_pre_init,GFS_rrtmg_pre_run
    ! Other modules required, e.g. type definitions
    use GFS_typedefs, only: GFS_control_type
    use GFS_typedefs, only: GFS_grid_type
    use GFS_typedefs, only: GFS_sfcprop_type
    use GFS_typedefs, only: GFS_statein_type
    use GFS_typedefs, only: GFS_tbd_type
    use GFS_typedefs, only: GFS_cldprop_type
    use GFS_typedefs, only: GFS_coupling_type
    use GFS_typedefs, only: GFS_radtend_type
    use machine, only: kind_phys

    implicit none

    private
    public :: GFS_rrtmg_pre_finalize_cap,GFS_rrtmg_pre_init_cap,GFS_rrtmg_pre_run_cap

    contains


    function GFS_rrtmg_pre_finalize_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call GFS_rrtmg_pre_finalize()
        

    end function GFS_rrtmg_pre_finalize_cap

    function GFS_rrtmg_pre_init_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call GFS_rrtmg_pre_init()
        

    end function GFS_rrtmg_pre_init_cap

    function GFS_rrtmg_pre_run_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind
        type(GFS_control_type), pointer     :: Model
        type(GFS_grid_type), pointer     :: Grid
        type(GFS_sfcprop_type), pointer     :: Sfcprop
        type(GFS_statein_type), pointer     :: Statein
        type(GFS_tbd_type), pointer     :: Tbd
        type(GFS_cldprop_type), pointer     :: Cldprop
        type(GFS_coupling_type), pointer     :: Coupling
        type(GFS_radtend_type), pointer     :: Radtend
        integer, pointer :: lm
        integer, pointer :: im
        integer, pointer :: lmk
        integer, pointer :: lmp
        integer, pointer :: kd
        integer, pointer :: kt
        integer, pointer :: kb
        real(kind_phys), pointer :: raddt
        real(kind_phys), pointer :: delp(:,:)
        real(kind_phys), pointer :: dz(:,:)
        real(kind_phys), pointer :: plvl(:,:)
        real(kind_phys), pointer :: plyr(:,:)
        real(kind_phys), pointer :: tlvl(:,:)
        real(kind_phys), pointer :: tlyr(:,:)
        real(kind_phys), pointer :: tsfg(:)
        real(kind_phys), pointer :: tsfa(:)
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
        real(kind_phys), pointer :: gasvmr_cfc113(:,:)
        real(kind_phys), pointer :: faersw1(:,:,:)
        real(kind_phys), pointer :: faersw2(:,:,:)
        real(kind_phys), pointer :: faersw3(:,:,:)
        real(kind_phys), pointer :: faerlw1(:,:,:)
        real(kind_phys), pointer :: faerlw2(:,:,:)
        real(kind_phys), pointer :: faerlw3(:,:,:)
        real(kind_phys), pointer :: aerodp(:,:)
        real(kind_phys), pointer :: clouds1(:,:)
        real(kind_phys), pointer :: clouds2(:,:)
        real(kind_phys), pointer :: clouds3(:,:)
        real(kind_phys), pointer :: clouds4(:,:)
        real(kind_phys), pointer :: clouds5(:,:)
        real(kind_phys), pointer :: clouds6(:,:)
        real(kind_phys), pointer :: clouds7(:,:)
        real(kind_phys), pointer :: clouds8(:,:)
        real(kind_phys), pointer :: clouds9(:,:)
        real(kind_phys), pointer :: cldsa(:,:)
        integer, pointer :: mtopa(:,:)
        integer, pointer :: mbota(:,:)
        real(kind_phys), pointer :: de_lgth(:)
        real(kind_phys), pointer :: alb1d(:)

        ierr = 0

        call c_f_pointer(ptr, cdata)


        call ccpp_field_get(cdata, 'GFS_control_type_instance', cptr, ierr=ierr, kind=ckind, index=2)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve GFS_control_type_instance from CCPP data structure')
            return
        end if
        if (ckind.ne.CCPP_GENERIC_KIND) then
            call ccpp_error('Kind mismatch for variable GFS_control_type_instance')
            ierr = 1
            return
        end if
#endif
        call c_f_pointer(cptr, Model)

        call ccpp_field_get(cdata, 'GFS_grid_type_instance', cptr, ierr=ierr, kind=ckind, index=6)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve GFS_grid_type_instance from CCPP data structure')
            return
        end if
        if (ckind.ne.CCPP_GENERIC_KIND) then
            call ccpp_error('Kind mismatch for variable GFS_grid_type_instance')
            ierr = 1
            return
        end if
#endif
        call c_f_pointer(cptr, Grid)

        call ccpp_field_get(cdata, 'GFS_sfcprop_type_instance', cptr, ierr=ierr, kind=ckind, index=11)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve GFS_sfcprop_type_instance from CCPP data structure')
            return
        end if
        if (ckind.ne.CCPP_GENERIC_KIND) then
            call ccpp_error('Kind mismatch for variable GFS_sfcprop_type_instance')
            ierr = 1
            return
        end if
#endif
        call c_f_pointer(cptr, Sfcprop)

        call ccpp_field_get(cdata, 'GFS_statein_type_instance', cptr, ierr=ierr, kind=ckind, index=13)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve GFS_statein_type_instance from CCPP data structure')
            return
        end if
        if (ckind.ne.CCPP_GENERIC_KIND) then
            call ccpp_error('Kind mismatch for variable GFS_statein_type_instance')
            ierr = 1
            return
        end if
#endif
        call c_f_pointer(cptr, Statein)

        call ccpp_field_get(cdata, 'GFS_tbd_type_instance', cptr, ierr=ierr, kind=ckind, index=15)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve GFS_tbd_type_instance from CCPP data structure')
            return
        end if
        if (ckind.ne.CCPP_GENERIC_KIND) then
            call ccpp_error('Kind mismatch for variable GFS_tbd_type_instance')
            ierr = 1
            return
        end if
#endif
        call c_f_pointer(cptr, Tbd)

        call ccpp_field_get(cdata, 'GFS_cldprop_type_instance', cptr, ierr=ierr, kind=ckind, index=1)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve GFS_cldprop_type_instance from CCPP data structure')
            return
        end if
        if (ckind.ne.CCPP_GENERIC_KIND) then
            call ccpp_error('Kind mismatch for variable GFS_cldprop_type_instance')
            ierr = 1
            return
        end if
#endif
        call c_f_pointer(cptr, Cldprop)

        call ccpp_field_get(cdata, 'GFS_coupling_type_instance', cptr, ierr=ierr, kind=ckind, index=3)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve GFS_coupling_type_instance from CCPP data structure')
            return
        end if
        if (ckind.ne.CCPP_GENERIC_KIND) then
            call ccpp_error('Kind mismatch for variable GFS_coupling_type_instance')
            ierr = 1
            return
        end if
#endif
        call c_f_pointer(cptr, Coupling)

        call ccpp_field_get(cdata, 'GFS_radtend_type_instance', cptr, ierr=ierr, kind=ckind, index=9)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve GFS_radtend_type_instance from CCPP data structure')
            return
        end if
        if (ckind.ne.CCPP_GENERIC_KIND) then
            call ccpp_error('Kind mismatch for variable GFS_radtend_type_instance')
            ierr = 1
            return
        end if
#endif
        call c_f_pointer(cptr, Radtend)

        call ccpp_field_get(cdata, 'vertical_layer_dimension_for_radiation', lm, ierr=ierr, kind=ckind, index=826)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve vertical_layer_dimension_for_radiation from CCPP data structure')
            return
        end if
        if (kind(lm).ne.ckind) then
            call ccpp_error('Kind mismatch for variable vertical_layer_dimension_for_radiation')
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
        

        call ccpp_field_get(cdata, 'adjusted_vertical_layer_dimension_for_radiation', lmk, ierr=ierr, kind=ckind, index=31)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve adjusted_vertical_layer_dimension_for_radiation from CCPP data structure')
            return
        end if
        if (kind(lmk).ne.ckind) then
            call ccpp_error('Kind mismatch for variable adjusted_vertical_layer_dimension_for_radiation')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'adjusted_vertical_level_dimension_for_radiation', lmp, ierr=ierr, kind=ckind, index=32)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve adjusted_vertical_level_dimension_for_radiation from CCPP data structure')
            return
        end if
        if (kind(lmp).ne.ckind) then
            call ccpp_error('Kind mismatch for variable adjusted_vertical_level_dimension_for_radiation')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'vertical_index_difference_between_inout_and_local', kd, ierr=ierr, kind=ckind, index=823)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve vertical_index_difference_between_inout_and_local from CCPP data structure')
            return
        end if
        if (kind(kd).ne.ckind) then
            call ccpp_error('Kind mismatch for variable vertical_index_difference_between_inout_and_local')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'vertical_index_difference_between_layer_and_upper_bound', kt, ierr=ierr, kind=ckind, index=825)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve vertical_index_difference_between_layer_and_upper_bound from CCPP data structure')
            return
        end if
        if (kind(kt).ne.ckind) then
            call ccpp_error('Kind mismatch for variable vertical_index_difference_between_layer_and_upper_bound')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'vertical_index_difference_between_layer_and_lower_bound', kb, ierr=ierr, kind=ckind, index=824)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve vertical_index_difference_between_layer_and_lower_bound from CCPP data structure')
            return
        end if
        if (kind(kb).ne.ckind) then
            call ccpp_error('Kind mismatch for variable vertical_index_difference_between_layer_and_lower_bound')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'time_step_for_radiation', raddt, ierr=ierr, kind=ckind, index=794)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve time_step_for_radiation from CCPP data structure')
            return
        end if
        if (kind(raddt).ne.ckind) then
            call ccpp_error('Kind mismatch for variable time_step_for_radiation')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'layer_pressure_thickness_for_radiation', delp, ierr=ierr, dims=cdims, kind=ckind, index=463)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve layer_pressure_thickness_for_radiation from CCPP data structure')
            return
        end if
        if (kind(delp).ne.ckind) then
            call ccpp_error('Kind mismatch for variable layer_pressure_thickness_for_radiation')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'layer_thickness_for_radiation', dz, ierr=ierr, dims=cdims, kind=ckind, index=464)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve layer_thickness_for_radiation from CCPP data structure')
            return
        end if
        if (kind(dz).ne.ckind) then
            call ccpp_error('Kind mismatch for variable layer_thickness_for_radiation')
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
        

        call ccpp_field_get(cdata, 'surface_ground_temperature_for_radiation', tsfg, ierr=ierr, dims=cdims, kind=ckind, index=712)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_ground_temperature_for_radiation from CCPP data structure')
            return
        end if
        if (kind(tsfg).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_ground_temperature_for_radiation')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_air_temperature_for_radiation', tsfa, ierr=ierr, dims=cdims, kind=ckind, index=681)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_air_temperature_for_radiation from CCPP data structure')
            return
        end if
        if (kind(tsfa).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_air_temperature_for_radiation')
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
        

        call ccpp_field_get(cdata, 'volume_mixing_ratio_cfc113', gasvmr_cfc113, ierr=ierr, dims=cdims, kind=ckind, index=840)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve volume_mixing_ratio_cfc113 from CCPP data structure')
            return
        end if
        if (kind(gasvmr_cfc113).ne.ckind) then
            call ccpp_error('Kind mismatch for variable volume_mixing_ratio_cfc113')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'aerosol_optical_depth_for_shortwave_bands_01-16', faersw1, ierr=ierr, dims=cdims, kind=ckind, index=39)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve aerosol_optical_depth_for_shortwave_bands_01-16 from CCPP data structure')
            return
        end if
        if (kind(faersw1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable aerosol_optical_depth_for_shortwave_bands_01-16')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'aerosol_single_scattering_albedo_for_shortwave_bands_01-16', faersw2, ierr=ierr, dims=cdims, kind=ckind, index=43)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve aerosol_single_scattering_albedo_for_shortwave_bands_01-16 from CCPP data structure')
            return
        end if
        if (kind(faersw2).ne.ckind) then
            call ccpp_error('Kind mismatch for variable aerosol_single_scattering_albedo_for_shortwave_bands_01-16')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'aerosol_asymmetry_parameter_for_shortwave_bands_01-16', faersw3, ierr=ierr, dims=cdims, kind=ckind, index=34)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve aerosol_asymmetry_parameter_for_shortwave_bands_01-16 from CCPP data structure')
            return
        end if
        if (kind(faersw3).ne.ckind) then
            call ccpp_error('Kind mismatch for variable aerosol_asymmetry_parameter_for_shortwave_bands_01-16')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'aerosol_optical_depth_for_longwave_bands_01-16', faerlw1, ierr=ierr, dims=cdims, kind=ckind, index=38)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve aerosol_optical_depth_for_longwave_bands_01-16 from CCPP data structure')
            return
        end if
        if (kind(faerlw1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable aerosol_optical_depth_for_longwave_bands_01-16')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'aerosol_single_scattering_albedo_for_longwave_bands_01-16', faerlw2, ierr=ierr, dims=cdims, kind=ckind, index=42)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve aerosol_single_scattering_albedo_for_longwave_bands_01-16 from CCPP data structure')
            return
        end if
        if (kind(faerlw2).ne.ckind) then
            call ccpp_error('Kind mismatch for variable aerosol_single_scattering_albedo_for_longwave_bands_01-16')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'aerosol_asymmetry_parameter_for_longwave_bands_01-16', faerlw3, ierr=ierr, dims=cdims, kind=ckind, index=33)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve aerosol_asymmetry_parameter_for_longwave_bands_01-16 from CCPP data structure')
            return
        end if
        if (kind(faerlw3).ne.ckind) then
            call ccpp_error('Kind mismatch for variable aerosol_asymmetry_parameter_for_longwave_bands_01-16')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'atmosphere_optical_thickness_due_to_ambient_aerosol_particles', aerodp, ierr=ierr, dims=cdims, kind=ckind, index=75)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve atmosphere_optical_thickness_due_to_ambient_aerosol_particles from CCPP data structure')
            return
        end if
        if (kind(aerodp).ne.ckind) then
            call ccpp_error('Kind mismatch for variable atmosphere_optical_thickness_due_to_ambient_aerosol_particles')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'total_cloud_fraction', clouds1, ierr=ierr, dims=cdims, kind=ckind, index=799)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve total_cloud_fraction from CCPP data structure')
            return
        end if
        if (kind(clouds1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable total_cloud_fraction')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cloud_liquid_water_path', clouds2, ierr=ierr, dims=cdims, kind=ckind, index=102)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cloud_liquid_water_path from CCPP data structure')
            return
        end if
        if (kind(clouds2).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cloud_liquid_water_path')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'mean_effective_radius_for_liquid_cloud', clouds3, ierr=ierr, dims=cdims, kind=ckind, index=516)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve mean_effective_radius_for_liquid_cloud from CCPP data structure')
            return
        end if
        if (kind(clouds3).ne.ckind) then
            call ccpp_error('Kind mismatch for variable mean_effective_radius_for_liquid_cloud')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cloud_ice_water_path', clouds4, ierr=ierr, dims=cdims, kind=ckind, index=101)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cloud_ice_water_path from CCPP data structure')
            return
        end if
        if (kind(clouds4).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cloud_ice_water_path')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'mean_effective_radius_for_ice_cloud', clouds5, ierr=ierr, dims=cdims, kind=ckind, index=515)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve mean_effective_radius_for_ice_cloud from CCPP data structure')
            return
        end if
        if (kind(clouds5).ne.ckind) then
            call ccpp_error('Kind mismatch for variable mean_effective_radius_for_ice_cloud')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cloud_rain_water_path', clouds6, ierr=ierr, dims=cdims, kind=ckind, index=107)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cloud_rain_water_path from CCPP data structure')
            return
        end if
        if (kind(clouds6).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cloud_rain_water_path')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'mean_effective_radius_for_rain_drop', clouds7, ierr=ierr, dims=cdims, kind=ckind, index=517)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve mean_effective_radius_for_rain_drop from CCPP data structure')
            return
        end if
        if (kind(clouds7).ne.ckind) then
            call ccpp_error('Kind mismatch for variable mean_effective_radius_for_rain_drop')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cloud_snow_water_path', clouds8, ierr=ierr, dims=cdims, kind=ckind, index=108)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cloud_snow_water_path from CCPP data structure')
            return
        end if
        if (kind(clouds8).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cloud_snow_water_path')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'mean_effective_radius_for_snow_flake', clouds9, ierr=ierr, dims=cdims, kind=ckind, index=518)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve mean_effective_radius_for_snow_flake from CCPP data structure')
            return
        end if
        if (kind(clouds9).ne.ckind) then
            call ccpp_error('Kind mismatch for variable mean_effective_radius_for_snow_flake')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cloud_area_fraction_for_radiation', cldsa, ierr=ierr, dims=cdims, kind=ckind, index=87)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cloud_area_fraction_for_radiation from CCPP data structure')
            return
        end if
        if (kind(cldsa).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cloud_area_fraction_for_radiation')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'model_layer_number_at_cloud_top', mtopa, ierr=ierr, dims=cdims, kind=ckind, index=552)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve model_layer_number_at_cloud_top from CCPP data structure')
            return
        end if
        if (kind(mtopa).ne.ckind) then
            call ccpp_error('Kind mismatch for variable model_layer_number_at_cloud_top')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'model_layer_number_at_cloud_base', mbota, ierr=ierr, dims=cdims, kind=ckind, index=551)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve model_layer_number_at_cloud_base from CCPP data structure')
            return
        end if
        if (kind(mbota).ne.ckind) then
            call ccpp_error('Kind mismatch for variable model_layer_number_at_cloud_base')
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
        

        call ccpp_field_get(cdata, 'surface_albedo_perturbation', alb1d, ierr=ierr, dims=cdims, kind=ckind, index=686)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_albedo_perturbation from CCPP data structure')
            return
        end if
        if (kind(alb1d).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_albedo_perturbation')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call GFS_rrtmg_pre_run(Model=Model,Grid=Grid,Sfcprop=Sfcprop,Statein=Statein,Tbd=Tbd,Cldprop=Cldprop, &
                  Coupling=Coupling,Radtend=Radtend,lm=lm,im=im,lmk=lmk,lmp=lmp,kd=kd,kt=kt, &
                  kb=kb,raddt=raddt,delp=delp,dz=dz,plvl=plvl,plyr=plyr,tlvl=tlvl,tlyr=tlyr, &
                  tsfg=tsfg,tsfa=tsfa,qlyr=qlyr,olyr=olyr,gasvmr_co2=gasvmr_co2,gasvmr_n2o=gasvmr_n2o, &
                  gasvmr_ch4=gasvmr_ch4,gasvmr_o2=gasvmr_o2,gasvmr_co=gasvmr_co,gasvmr_cfc11=gasvmr_cfc11, &
                  gasvmr_cfc12=gasvmr_cfc12,gasvmr_cfc22=gasvmr_cfc22,gasvmr_ccl4=gasvmr_ccl4, &
                  gasvmr_cfc113=gasvmr_cfc113,faersw1=faersw1,faersw2=faersw2,faersw3=faersw3, &
                  faerlw1=faerlw1,faerlw2=faerlw2,faerlw3=faerlw3,aerodp=aerodp,clouds1=clouds1, &
                  clouds2=clouds2,clouds3=clouds3,clouds4=clouds4,clouds5=clouds5,clouds6=clouds6, &
                  clouds7=clouds7,clouds8=clouds8,clouds9=clouds9,cldsa=cldsa,mtopa=mtopa, &
                  mbota=mbota,de_lgth=de_lgth,alb1d=alb1d,errmsg=cdata%errmsg,errflg=cdata%errflg)
        ierr=cdata%errflg

    end function GFS_rrtmg_pre_run_cap
end module GFS_rrtmg_pre_cap

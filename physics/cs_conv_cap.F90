
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
!! @brief Auto-generated cap module for the cs_conv scheme
!!
!
module cs_conv_cap

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_f_pointer, c_ptr, c_int32_t
    use            :: ccpp_types,                                      &
                      only: ccpp_t, CCPP_GENERIC_KIND
    use            :: ccpp_fields,                                     &
                      only: ccpp_field_get
    use            :: ccpp_errors,                                     &
                      only: ccpp_error, ccpp_debug
    use            :: cs_conv, &
                      only: cs_conv_init,cs_conv_finalize,cs_conv_run
    ! Other modules required, e.g. type definitions
    use machine, only: kind_phys

    implicit none

    private
    public :: cs_conv_init_cap,cs_conv_finalize_cap,cs_conv_run_cap

    contains


    function cs_conv_init_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call cs_conv_init()
        

    end function cs_conv_init_cap

    function cs_conv_finalize_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call cs_conv_finalize()
        

    end function cs_conv_finalize_cap

    function cs_conv_run_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind
        integer, pointer :: im
        integer, pointer :: ijsdim
        integer, pointer :: kmax
        integer, pointer :: ntracp1
        integer, pointer :: nn
        integer, pointer :: ntr
        integer, pointer :: nctp
        logical, pointer :: otspt(:,:)
        integer, pointer :: lat
        integer, pointer :: kdt
        real(kind_phys), pointer :: t(:,:)
        real(kind_phys), pointer :: q(:,:)
        real(kind_phys), pointer :: rain1(:)
        real(kind_phys), pointer :: clw(:,:,:)
        real(kind_phys), pointer :: zm(:,:)
        real(kind_phys), pointer :: zi(:,:)
        real(kind_phys), pointer :: pap(:,:)
        real(kind_phys), pointer :: paph(:,:)
        real(kind_phys), pointer :: delta
        real(kind_phys), pointer :: delti
        real(kind_phys), pointer :: ud_mf(:,:)
        real(kind_phys), pointer :: dd_mf(:,:)
        real(kind_phys), pointer :: dt_mf(:,:)
        real(kind_phys), pointer :: u(:,:)
        real(kind_phys), pointer :: v(:,:)
        real(kind_phys), pointer :: fscav(:)
        real(kind_phys), pointer :: fswtr(:)
        real(kind_phys), pointer :: cbmfx(:,:)
        integer, pointer :: mype
        real(kind_phys), pointer :: wcbmaxm(:)
        real(kind_phys), pointer :: precz0in
        real(kind_phys), pointer :: preczhin
        real(kind_phys), pointer :: clmdin
        real(kind_phys), pointer :: sigma(:,:)
        logical, pointer :: do_aw
        logical, pointer :: do_awdd
        logical, pointer :: flx_form
        logical, pointer :: lprnt
        integer, pointer :: ipr
        integer, pointer :: kcnv(:)
        real(kind_phys), pointer :: qlcn(:,:)
        real(kind_phys), pointer :: qicn(:,:)
        real(kind_phys), pointer :: w_upi(:,:)
        real(kind_phys), pointer :: cf_upi(:,:)
        real(kind_phys), pointer :: cnv_mfd(:,:)
        real(kind_phys), pointer :: cnv_dqldt(:,:)
        real(kind_phys), pointer :: clcn(:,:)
        real(kind_phys), pointer :: cnv_fice(:,:)
        real(kind_phys), pointer :: cnv_ndrop(:,:)
        real(kind_phys), pointer :: cnv_nice(:,:)
        integer, pointer :: mp_phys

        ierr = 0

        call c_f_pointer(ptr, cdata)


        call ccpp_field_get(cdata, 'horizontal_dimension', im, ierr=ierr, kind=ckind, index=364)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve horizontal_dimension from CCPP data structure')
            return
        end if
        if (kind(im).ne.ckind) then
            call ccpp_error('Kind mismatch for variable horizontal_dimension')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'horizontal_loop_extent', ijsdim, ierr=ierr, kind=ckind, index=366)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve horizontal_loop_extent from CCPP data structure')
            return
        end if
        if (kind(ijsdim).ne.ckind) then
            call ccpp_error('Kind mismatch for variable horizontal_loop_extent')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'vertical_dimension', kmax, ierr=ierr, kind=ckind, index=817)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve vertical_dimension from CCPP data structure')
            return
        end if
        if (kind(kmax).ne.ckind) then
            call ccpp_error('Kind mismatch for variable vertical_dimension')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'number_of_tracers_plus_one', ntracp1, ierr=ierr, kind=ckind, index=589)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve number_of_tracers_plus_one from CCPP data structure')
            return
        end if
        if (kind(ntracp1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable number_of_tracers_plus_one')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'number_of_tracers_for_convective_transport', nn, ierr=ierr, kind=ckind, index=587)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve number_of_tracers_for_convective_transport from CCPP data structure')
            return
        end if
        if (kind(nn).ne.ckind) then
            call ccpp_error('Kind mismatch for variable number_of_tracers_for_convective_transport')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'number_of_tracers_for_CS', ntr, ierr=ierr, kind=ckind, index=585)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve number_of_tracers_for_CS from CCPP data structure')
            return
        end if
        if (kind(ntr).ne.ckind) then
            call ccpp_error('Kind mismatch for variable number_of_tracers_for_CS')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'number_of_cloud_types_CS', nctp, ierr=ierr, kind=ckind, index=572)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve number_of_cloud_types_CS from CCPP data structure')
            return
        end if
        if (kind(nctp).ne.ckind) then
            call ccpp_error('Kind mismatch for variable number_of_cloud_types_CS')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'flag_convective_tracer_transport', otspt, ierr=ierr, dims=cdims, kind=ckind, index=261)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_convective_tracer_transport from CCPP data structure')
            return
        end if
        if (kind(otspt).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_convective_tracer_transport')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'latitude_index_in_debug_printouts', lat, ierr=ierr, kind=ckind, index=462)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve latitude_index_in_debug_printouts from CCPP data structure')
            return
        end if
        if (kind(lat).ne.ckind) then
            call ccpp_error('Kind mismatch for variable latitude_index_in_debug_printouts')
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
        

        call ccpp_field_get(cdata, 'air_temperature_updated_by_physics', t, ierr=ierr, dims=cdims, kind=ckind, index=59)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve air_temperature_updated_by_physics from CCPP data structure')
            return
        end if
        if (kind(t).ne.ckind) then
            call ccpp_error('Kind mismatch for variable air_temperature_updated_by_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'water_vapor_specific_humidity_updated_by_physics', q, ierr=ierr, dims=cdims, kind=ckind, index=860)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve water_vapor_specific_humidity_updated_by_physics from CCPP data structure')
            return
        end if
        if (kind(q).ne.ckind) then
            call ccpp_error('Kind mismatch for variable water_vapor_specific_humidity_updated_by_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'lwe_thickness_of_deep_convective_precipitation_amount', rain1, ierr=ierr, dims=cdims, kind=ckind, index=479)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve lwe_thickness_of_deep_convective_precipitation_amount from CCPP data structure')
            return
        end if
        if (kind(rain1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable lwe_thickness_of_deep_convective_precipitation_amount')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'convective_transportable_tracers', clw, ierr=ierr, dims=cdims, kind=ckind, index=130)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve convective_transportable_tracers from CCPP data structure')
            return
        end if
        if (kind(clw).ne.ckind) then
            call ccpp_error('Kind mismatch for variable convective_transportable_tracers')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'geopotential', zm, ierr=ierr, dims=cdims, kind=ckind, index=348)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve geopotential from CCPP data structure')
            return
        end if
        if (kind(zm).ne.ckind) then
            call ccpp_error('Kind mismatch for variable geopotential')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'geopotential_at_interface', zi, ierr=ierr, dims=cdims, kind=ckind, index=349)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve geopotential_at_interface from CCPP data structure')
            return
        end if
        if (kind(zi).ne.ckind) then
            call ccpp_error('Kind mismatch for variable geopotential_at_interface')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'air_pressure', pap, ierr=ierr, dims=cdims, kind=ckind, index=44)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve air_pressure from CCPP data structure')
            return
        end if
        if (kind(pap).ne.ckind) then
            call ccpp_error('Kind mismatch for variable air_pressure')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'air_pressure_at_interface', paph, ierr=ierr, dims=cdims, kind=ckind, index=45)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve air_pressure_at_interface from CCPP data structure')
            return
        end if
        if (kind(paph).ne.ckind) then
            call ccpp_error('Kind mismatch for variable air_pressure_at_interface')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'time_step_for_physics', delta, ierr=ierr, kind=ckind, index=793)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve time_step_for_physics from CCPP data structure')
            return
        end if
        if (kind(delta).ne.ckind) then
            call ccpp_error('Kind mismatch for variable time_step_for_physics')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'time_step_for_dynamics', delti, ierr=ierr, kind=ckind, index=792)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve time_step_for_dynamics from CCPP data structure')
            return
        end if
        if (kind(delti).ne.ckind) then
            call ccpp_error('Kind mismatch for variable time_step_for_dynamics')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'instantaneous_atmosphere_updraft_convective_mass_flux', ud_mf, ierr=ierr, dims=cdims, kind=ckind, index=403)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve instantaneous_atmosphere_updraft_convective_mass_flux from CCPP data structure')
            return
        end if
        if (kind(ud_mf).ne.ckind) then
            call ccpp_error('Kind mismatch for variable instantaneous_atmosphere_updraft_convective_mass_flux')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'instantaneous_atmosphere_downdraft_convective_mass_flux', dd_mf, ierr=ierr, dims=cdims, kind=ckind, index=401)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve instantaneous_atmosphere_downdraft_convective_mass_flux from CCPP data structure')
            return
        end if
        if (kind(dd_mf).ne.ckind) then
            call ccpp_error('Kind mismatch for variable instantaneous_atmosphere_downdraft_convective_mass_flux')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'instantaneous_atmosphere_detrainment_convective_mass_flux', dt_mf, ierr=ierr, dims=cdims, kind=ckind, index=399)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve instantaneous_atmosphere_detrainment_convective_mass_flux from CCPP data structure')
            return
        end if
        if (kind(dt_mf).ne.ckind) then
            call ccpp_error('Kind mismatch for variable instantaneous_atmosphere_detrainment_convective_mass_flux')
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
        

        call ccpp_field_get(cdata, 'fraction_of_tracer_scavenged', fscav, ierr=ierr, dims=cdims, kind=ckind, index=343)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve fraction_of_tracer_scavenged from CCPP data structure')
            return
        end if
        if (kind(fscav).ne.ckind) then
            call ccpp_error('Kind mismatch for variable fraction_of_tracer_scavenged')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'fraction_of_cloud_top_water_scavenged', fswtr, ierr=ierr, dims=cdims, kind=ckind, index=340)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve fraction_of_cloud_top_water_scavenged from CCPP data structure')
            return
        end if
        if (kind(fswtr).ne.ckind) then
            call ccpp_error('Kind mismatch for variable fraction_of_cloud_top_water_scavenged')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cloud_base_mass_flux', cbmfx, ierr=ierr, dims=cdims, kind=ckind, index=88)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cloud_base_mass_flux from CCPP data structure')
            return
        end if
        if (kind(cbmfx).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cloud_base_mass_flux')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'mpi_rank', mype, ierr=ierr, kind=ckind, index=558)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve mpi_rank from CCPP data structure')
            return
        end if
        if (kind(mype).ne.ckind) then
            call ccpp_error('Kind mismatch for variable mpi_rank')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'maximum_updraft_velocity_at_cloud_base', wcbmaxm, ierr=ierr, dims=cdims, kind=ckind, index=509)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve maximum_updraft_velocity_at_cloud_base from CCPP data structure')
            return
        end if
        if (kind(wcbmaxm).ne.ckind) then
            call ccpp_error('Kind mismatch for variable maximum_updraft_velocity_at_cloud_base')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'detrainment_and_precipitation_tunable_parameter_3_CS', precz0in, ierr=ierr, kind=ckind, index=214)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve detrainment_and_precipitation_tunable_parameter_3_CS from CCPP data structure')
            return
        end if
        if (kind(precz0in).ne.ckind) then
            call ccpp_error('Kind mismatch for variable detrainment_and_precipitation_tunable_parameter_3_CS')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'detrainment_and_precipitation_tunable_parameter_4_CS', preczhin, ierr=ierr, kind=ckind, index=215)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve detrainment_and_precipitation_tunable_parameter_4_CS from CCPP data structure')
            return
        end if
        if (kind(preczhin).ne.ckind) then
            call ccpp_error('Kind mismatch for variable detrainment_and_precipitation_tunable_parameter_4_CS')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'entrainment_efficiency_tunable_parameter_9_CS', clmdin, ierr=ierr, kind=ckind, index=253)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve entrainment_efficiency_tunable_parameter_9_CS from CCPP data structure')
            return
        end if
        if (kind(clmdin).ne.ckind) then
            call ccpp_error('Kind mismatch for variable entrainment_efficiency_tunable_parameter_9_CS')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'convective_updraft_area_fraction_at_model_interfaces', sigma, ierr=ierr, dims=cdims, kind=ckind, index=132)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve convective_updraft_area_fraction_at_model_interfaces from CCPP data structure')
            return
        end if
        if (kind(sigma).ne.ckind) then
            call ccpp_error('Kind mismatch for variable convective_updraft_area_fraction_at_model_interfaces')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'flag_for_Arakawa_Wu_adjustment', do_aw, ierr=ierr, kind=ckind, index=267)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_Arakawa_Wu_adjustment from CCPP data structure')
            return
        end if
        if (kind(do_aw).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_Arakawa_Wu_adjustment')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'flag_arakawa_wu_downdraft', do_awdd, ierr=ierr, kind=ckind, index=259)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_arakawa_wu_downdraft from CCPP data structure')
            return
        end if
        if (kind(do_awdd).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_arakawa_wu_downdraft')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'flag_flux_form_CS', flx_form, ierr=ierr, kind=ckind, index=266)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_flux_form_CS from CCPP data structure')
            return
        end if
        if (kind(flx_form).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_flux_form_CS')
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
        

        call ccpp_field_get(cdata, 'flag_deep_convection', kcnv, ierr=ierr, dims=cdims, kind=ckind, index=262)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_deep_convection from CCPP data structure')
            return
        end if
        if (kind(kcnv).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_deep_convection')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'mass_fraction_of_convective_cloud_liquid_water', qlcn, ierr=ierr, dims=cdims, kind=ckind, index=502)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve mass_fraction_of_convective_cloud_liquid_water from CCPP data structure')
            return
        end if
        if (kind(qlcn).ne.ckind) then
            call ccpp_error('Kind mismatch for variable mass_fraction_of_convective_cloud_liquid_water')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'mass_fraction_of_convective_cloud_ice', qicn, ierr=ierr, dims=cdims, kind=ckind, index=501)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve mass_fraction_of_convective_cloud_ice from CCPP data structure')
            return
        end if
        if (kind(qicn).ne.ckind) then
            call ccpp_error('Kind mismatch for variable mass_fraction_of_convective_cloud_ice')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'vertical_velocity_for_updraft', w_upi, ierr=ierr, dims=cdims, kind=ckind, index=830)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve vertical_velocity_for_updraft from CCPP data structure')
            return
        end if
        if (kind(w_upi).ne.ckind) then
            call ccpp_error('Kind mismatch for variable vertical_velocity_for_updraft')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'convective_cloud_fraction_for_microphysics', cf_upi, ierr=ierr, dims=cdims, kind=ckind, index=125)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve convective_cloud_fraction_for_microphysics from CCPP data structure')
            return
        end if
        if (kind(cf_upi).ne.ckind) then
            call ccpp_error('Kind mismatch for variable convective_cloud_fraction_for_microphysics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'detrained_mass_flux', cnv_mfd, ierr=ierr, dims=cdims, kind=ckind, index=213)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve detrained_mass_flux from CCPP data structure')
            return
        end if
        if (kind(cnv_mfd).ne.ckind) then
            call ccpp_error('Kind mismatch for variable detrained_mass_flux')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'tendency_of_cloud_water_due_to_convective_microphysics', cnv_dqldt, ierr=ierr, dims=cdims, kind=ckind, index=766)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve tendency_of_cloud_water_due_to_convective_microphysics from CCPP data structure')
            return
        end if
        if (kind(cnv_dqldt).ne.ckind) then
            call ccpp_error('Kind mismatch for variable tendency_of_cloud_water_due_to_convective_microphysics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'convective_cloud_volume_fraction', clcn, ierr=ierr, dims=cdims, kind=ckind, index=127)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve convective_cloud_volume_fraction from CCPP data structure')
            return
        end if
        if (kind(clcn).ne.ckind) then
            call ccpp_error('Kind mismatch for variable convective_cloud_volume_fraction')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'ice_fraction_in_convective_tower', cnv_fice, ierr=ierr, dims=cdims, kind=ckind, index=367)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve ice_fraction_in_convective_tower from CCPP data structure')
            return
        end if
        if (kind(cnv_fice).ne.ckind) then
            call ccpp_error('Kind mismatch for variable ice_fraction_in_convective_tower')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'number_concentration_of_cloud_liquid_water_particles_for_detrainment', cnv_ndrop, ierr=ierr, dims=cdims, kind=ckind, index=569)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve number_concentration_of_cloud_liquid_water_particles_for_detrainment from CCPP data structure')
            return
        end if
        if (kind(cnv_ndrop).ne.ckind) then
            call ccpp_error('Kind mismatch for variable number_concentration_of_cloud_liquid_water_particles_for_detrainment')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'number_concentration_of_ice_crystals_for_detrainment', cnv_nice, ierr=ierr, dims=cdims, kind=ckind, index=570)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve number_concentration_of_ice_crystals_for_detrainment from CCPP data structure')
            return
        end if
        if (kind(cnv_nice).ne.ckind) then
            call ccpp_error('Kind mismatch for variable number_concentration_of_ice_crystals_for_detrainment')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'flag_for_microphysics_scheme', mp_phys, ierr=ierr, kind=ckind, index=294)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_microphysics_scheme from CCPP data structure')
            return
        end if
        if (kind(mp_phys).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_microphysics_scheme')
            ierr = 1
            return
        end if
#endif
        

        call cs_conv_run(im=im,ijsdim=ijsdim,kmax=kmax,ntracp1=ntracp1,nn=nn,ntr=ntr,nctp=nctp,otspt=otspt, &
                  lat=lat,kdt=kdt,t=t,q=q,rain1=rain1,clw=clw,zm=zm,zi=zi,pap=pap,paph=paph, &
                  delta=delta,delti=delti,ud_mf=ud_mf,dd_mf=dd_mf,dt_mf=dt_mf,u=u,v=v,fscav=fscav, &
                  fswtr=fswtr,cbmfx=cbmfx,mype=mype,wcbmaxm=wcbmaxm,precz0in=precz0in,preczhin=preczhin, &
                  clmdin=clmdin,sigma=sigma,do_aw=do_aw,do_awdd=do_awdd,flx_form=flx_form, &
                  lprnt=lprnt,ipr=ipr,kcnv=kcnv,qlcn=qlcn,qicn=qicn,w_upi=w_upi,cf_upi=cf_upi, &
                  cnv_mfd=cnv_mfd,cnv_dqldt=cnv_dqldt,clcn=clcn,cnv_fice=cnv_fice,cnv_ndrop=cnv_ndrop, &
                  cnv_nice=cnv_nice,mp_phys=mp_phys,errmsg=cdata%errmsg,errflg=cdata%errflg)
        ierr=cdata%errflg

    end function cs_conv_run_cap
end module cs_conv_cap

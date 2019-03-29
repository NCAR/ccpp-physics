
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
!! @brief Auto-generated cap module for the GFS_stochastics scheme
!!
!
module GFS_stochastics_cap

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_f_pointer, c_ptr, c_int32_t
    use            :: ccpp_types,                                      &
                      only: ccpp_t, CCPP_GENERIC_KIND
    use            :: ccpp_fields,                                     &
                      only: ccpp_field_get
    use            :: ccpp_errors,                                     &
                      only: ccpp_error, ccpp_debug
    use            :: GFS_stochastics, &
                      only: GFS_stochastics_run,GFS_stochastics_finalize,GFS_stochastics_init
    ! Other modules required, e.g. type definitions
    use machine, only: kind_phys

    implicit none

    private
    public :: GFS_stochastics_run_cap,GFS_stochastics_finalize_cap,GFS_stochastics_init_cap

    contains


    function GFS_stochastics_run_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind
        integer, pointer :: im
        integer, pointer :: km
        logical, pointer :: do_sppt
        logical, pointer :: use_zmtnblck
        logical, pointer :: do_shum
        logical, pointer :: do_skeb
        real(kind_phys), pointer :: zmtnblck(:)
        real(kind_phys), pointer :: sppt_wts(:,:)
        real(kind_phys), pointer :: skebu_wts(:,:)
        real(kind_phys), pointer :: skebv_wts(:,:)
        real(kind_phys), pointer :: shum_wts(:,:)
        real(kind_phys), pointer :: sppt_wts_inv(:,:)
        real(kind_phys), pointer :: skebu_wts_inv(:,:)
        real(kind_phys), pointer :: skebv_wts_inv(:,:)
        real(kind_phys), pointer :: shum_wts_inv(:,:)
        real(kind_phys), pointer :: diss_est(:,:)
        real(kind_phys), pointer :: ugrs(:,:)
        real(kind_phys), pointer :: vgrs(:,:)
        real(kind_phys), pointer :: tgrs(:,:)
        real(kind_phys), pointer :: qgrs(:,:)
        real(kind_phys), pointer :: gu0(:,:)
        real(kind_phys), pointer :: gv0(:,:)
        real(kind_phys), pointer :: gt0(:,:)
        real(kind_phys), pointer :: gq0(:,:)
        real(kind_phys), pointer :: dtdtr(:,:)
        real(kind_phys), pointer :: rain(:)
        real(kind_phys), pointer :: rainc(:)
        real(kind_phys), pointer :: tprcp(:)
        real(kind_phys), pointer :: totprcp(:)
        real(kind_phys), pointer :: cnvprcp(:)
        real(kind_phys), pointer :: totprcpb(:)
        real(kind_phys), pointer :: cnvprcpb(:)
        logical, pointer :: cplflx
        real(kind_phys), pointer :: rain_cpl(:)
        real(kind_phys), pointer :: snow_cpl(:)
        real(kind_phys), pointer :: drain_cpl(:)
        real(kind_phys), pointer :: dsnow_cpl(:)

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
        

        call ccpp_field_get(cdata, 'vertical_dimension', km, ierr=ierr, kind=ckind, index=817)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve vertical_dimension from CCPP data structure')
            return
        end if
        if (kind(km).ne.ckind) then
            call ccpp_error('Kind mismatch for variable vertical_dimension')
            ierr = 1
            return
        end if
#endif
        

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
        

        call ccpp_field_get(cdata, 'flag_for_mountain_blocking', use_zmtnblck, ierr=ierr, kind=ckind, index=298)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_mountain_blocking from CCPP data structure')
            return
        end if
        if (kind(use_zmtnblck).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_mountain_blocking')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'flag_for_stochastic_shum_option', do_shum, ierr=ierr, kind=ckind, index=316)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_stochastic_shum_option from CCPP data structure')
            return
        end if
        if (kind(do_shum).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_stochastic_shum_option')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'flag_for_stochastic_skeb_option', do_skeb, ierr=ierr, kind=ckind, index=317)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_stochastic_skeb_option from CCPP data structure')
            return
        end if
        if (kind(do_skeb).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_stochastic_skeb_option')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'level_of_dividing_streamline', zmtnblck, ierr=ierr, dims=cdims, kind=ckind, index=465)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve level_of_dividing_streamline from CCPP data structure')
            return
        end if
        if (kind(zmtnblck).ne.ckind) then
            call ccpp_error('Kind mismatch for variable level_of_dividing_streamline')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'weights_for_stochastic_sppt_perturbation', sppt_wts, ierr=ierr, dims=cdims, kind=ckind, index=867)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve weights_for_stochastic_sppt_perturbation from CCPP data structure')
            return
        end if
        if (kind(sppt_wts).ne.ckind) then
            call ccpp_error('Kind mismatch for variable weights_for_stochastic_sppt_perturbation')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'weights_for_stochastic_skeb_perturbation_of_x_wind', skebu_wts, ierr=ierr, dims=cdims, kind=ckind, index=863)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve weights_for_stochastic_skeb_perturbation_of_x_wind from CCPP data structure')
            return
        end if
        if (kind(skebu_wts).ne.ckind) then
            call ccpp_error('Kind mismatch for variable weights_for_stochastic_skeb_perturbation_of_x_wind')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'weights_for_stochastic_skeb_perturbation_of_y_wind', skebv_wts, ierr=ierr, dims=cdims, kind=ckind, index=865)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve weights_for_stochastic_skeb_perturbation_of_y_wind from CCPP data structure')
            return
        end if
        if (kind(skebv_wts).ne.ckind) then
            call ccpp_error('Kind mismatch for variable weights_for_stochastic_skeb_perturbation_of_y_wind')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'weights_for_stochastic_shum_perturbation', shum_wts, ierr=ierr, dims=cdims, kind=ckind, index=861)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve weights_for_stochastic_shum_perturbation from CCPP data structure')
            return
        end if
        if (kind(shum_wts).ne.ckind) then
            call ccpp_error('Kind mismatch for variable weights_for_stochastic_shum_perturbation')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'weights_for_stochastic_sppt_perturbation_flipped', sppt_wts_inv, ierr=ierr, dims=cdims, kind=ckind, index=868)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve weights_for_stochastic_sppt_perturbation_flipped from CCPP data structure')
            return
        end if
        if (kind(sppt_wts_inv).ne.ckind) then
            call ccpp_error('Kind mismatch for variable weights_for_stochastic_sppt_perturbation_flipped')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'weights_for_stochastic_skeb_perturbation_of_x_wind_flipped', skebu_wts_inv, ierr=ierr, dims=cdims, kind=ckind, index=864)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve weights_for_stochastic_skeb_perturbation_of_x_wind_flipped from CCPP data structure')
            return
        end if
        if (kind(skebu_wts_inv).ne.ckind) then
            call ccpp_error('Kind mismatch for variable weights_for_stochastic_skeb_perturbation_of_x_wind_flipped')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'weights_for_stochastic_skeb_perturbation_of_y_wind_flipped', skebv_wts_inv, ierr=ierr, dims=cdims, kind=ckind, index=866)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve weights_for_stochastic_skeb_perturbation_of_y_wind_flipped from CCPP data structure')
            return
        end if
        if (kind(skebv_wts_inv).ne.ckind) then
            call ccpp_error('Kind mismatch for variable weights_for_stochastic_skeb_perturbation_of_y_wind_flipped')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'weights_for_stochastic_shum_perturbation_flipped', shum_wts_inv, ierr=ierr, dims=cdims, kind=ckind, index=862)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve weights_for_stochastic_shum_perturbation_flipped from CCPP data structure')
            return
        end if
        if (kind(shum_wts_inv).ne.ckind) then
            call ccpp_error('Kind mismatch for variable weights_for_stochastic_shum_perturbation_flipped')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'dissipation_estimate_of_air_temperature_at_model_layers', diss_est, ierr=ierr, dims=cdims, kind=ckind, index=223)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve dissipation_estimate_of_air_temperature_at_model_layers from CCPP data structure')
            return
        end if
        if (kind(diss_est).ne.ckind) then
            call ccpp_error('Kind mismatch for variable dissipation_estimate_of_air_temperature_at_model_layers')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'x_wind', ugrs, ierr=ierr, dims=cdims, kind=ckind, index=871)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve x_wind from CCPP data structure')
            return
        end if
        if (kind(ugrs).ne.ckind) then
            call ccpp_error('Kind mismatch for variable x_wind')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'y_wind', vgrs, ierr=ierr, dims=cdims, kind=ckind, index=878)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve y_wind from CCPP data structure')
            return
        end if
        if (kind(vgrs).ne.ckind) then
            call ccpp_error('Kind mismatch for variable y_wind')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'air_temperature', tgrs, ierr=ierr, dims=cdims, kind=ckind, index=50)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve air_temperature from CCPP data structure')
            return
        end if
        if (kind(tgrs).ne.ckind) then
            call ccpp_error('Kind mismatch for variable air_temperature')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'water_vapor_specific_humidity', qgrs, ierr=ierr, dims=cdims, kind=ckind, index=852)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve water_vapor_specific_humidity from CCPP data structure')
            return
        end if
        if (kind(qgrs).ne.ckind) then
            call ccpp_error('Kind mismatch for variable water_vapor_specific_humidity')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'x_wind_updated_by_physics', gu0, ierr=ierr, dims=cdims, kind=ckind, index=877)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve x_wind_updated_by_physics from CCPP data structure')
            return
        end if
        if (kind(gu0).ne.ckind) then
            call ccpp_error('Kind mismatch for variable x_wind_updated_by_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'y_wind_updated_by_physics', gv0, ierr=ierr, dims=cdims, kind=ckind, index=884)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve y_wind_updated_by_physics from CCPP data structure')
            return
        end if
        if (kind(gv0).ne.ckind) then
            call ccpp_error('Kind mismatch for variable y_wind_updated_by_physics')
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
        

        call ccpp_field_get(cdata, 'water_vapor_specific_humidity_updated_by_physics', gq0, ierr=ierr, dims=cdims, kind=ckind, index=860)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve water_vapor_specific_humidity_updated_by_physics from CCPP data structure')
            return
        end if
        if (kind(gq0).ne.ckind) then
            call ccpp_error('Kind mismatch for variable water_vapor_specific_humidity_updated_by_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

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
        

        call ccpp_field_get(cdata, 'cumulative_lwe_thickness_of_convective_precipitation_amount', cnvprcp, ierr=ierr, dims=cdims, kind=ckind, index=174)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cumulative_lwe_thickness_of_convective_precipitation_amount from CCPP data structure')
            return
        end if
        if (kind(cnvprcp).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cumulative_lwe_thickness_of_convective_precipitation_amount')
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
        

        call ccpp_field_get(cdata, 'cumulative_lwe_thickness_of_convective_precipitation_amount_in_bucket', cnvprcpb, ierr=ierr, dims=cdims, kind=ckind, index=175)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cumulative_lwe_thickness_of_convective_precipitation_amount_in_bucket from CCPP data structure')
            return
        end if
        if (kind(cnvprcpb).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cumulative_lwe_thickness_of_convective_precipitation_amount_in_bucket')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

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
        

        call GFS_stochastics_run(im=im,km=km,do_sppt=do_sppt,use_zmtnblck=use_zmtnblck,do_shum=do_shum,do_skeb=do_skeb, &
                  zmtnblck=zmtnblck,sppt_wts=sppt_wts,skebu_wts=skebu_wts,skebv_wts=skebv_wts, &
                  shum_wts=shum_wts,sppt_wts_inv=sppt_wts_inv,skebu_wts_inv=skebu_wts_inv, &
                  skebv_wts_inv=skebv_wts_inv,shum_wts_inv=shum_wts_inv,diss_est=diss_est, &
                  ugrs=ugrs,vgrs=vgrs,tgrs=tgrs,qgrs=qgrs,gu0=gu0,gv0=gv0,gt0=gt0,gq0=gq0, &
                  dtdtr=dtdtr,rain=rain,rainc=rainc,tprcp=tprcp,totprcp=totprcp,cnvprcp=cnvprcp, &
                  totprcpb=totprcpb,cnvprcpb=cnvprcpb,cplflx=cplflx,rain_cpl=rain_cpl,snow_cpl=snow_cpl, &
                  drain_cpl=drain_cpl,dsnow_cpl=dsnow_cpl,errmsg=cdata%errmsg,errflg=cdata%errflg)
        ierr=cdata%errflg

    end function GFS_stochastics_run_cap

    function GFS_stochastics_finalize_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call GFS_stochastics_finalize()
        

    end function GFS_stochastics_finalize_cap

    function GFS_stochastics_init_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call GFS_stochastics_init()
        

    end function GFS_stochastics_init_cap
end module GFS_stochastics_cap

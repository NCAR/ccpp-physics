
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
!! @brief Auto-generated cap module for the cires_ugwp scheme
!!
!
module cires_ugwp_cap

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_f_pointer, c_ptr, c_int32_t
    use            :: ccpp_types,                                      &
                      only: ccpp_t, CCPP_GENERIC_KIND
    use            :: ccpp_fields,                                     &
                      only: ccpp_field_get
    use            :: ccpp_errors,                                     &
                      only: ccpp_error, ccpp_debug
    use            :: cires_ugwp, &
                      only: cires_ugwp_run,cires_ugwp_init,cires_ugwp_finalize
    ! Other modules required, e.g. type definitions
    use GFS_typedefs, only: GFS_control_type
    use machine, only: kind_phys
    use GFS_typedefs, only: GFS_statein_type
    use GFS_typedefs, only: GFS_sfcprop_type

    implicit none

    private
    public :: cires_ugwp_run_cap,cires_ugwp_init_cap,cires_ugwp_finalize_cap

    contains


    function cires_ugwp_run_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind
        type(GFS_control_type), pointer     :: Model
        integer, pointer :: me
        integer, pointer :: master
        integer, pointer :: im
        integer, pointer :: levs
        integer, pointer :: ntrac
        integer, pointer :: nvoro
        real(kind_phys), pointer :: dtp
        integer, pointer :: kdt
        integer, pointer :: lonr
        integer, pointer :: nvoro
        logical, pointer :: do_tofd
        real(kind_phys), pointer :: cdmbgwd(:)
        real(kind_phys), pointer :: xlat(:)
        real(kind_phys), pointer :: xlat_d(:)
        real(kind_phys), pointer :: sinlat(:)
        real(kind_phys), pointer :: coslat(:)
        type(GFS_statein_type), pointer     :: Statein
        type(GFS_sfcprop_type), pointer     :: Sfcprop
        real(kind_phys), pointer :: del(:,:)
        real(kind_phys), pointer :: orostat(:,:)
        integer, pointer :: kpbl(:)
        real(kind_phys), pointer :: dusfcg(:)
        real(kind_phys), pointer :: dvsfcg(:)
        real(kind_phys), pointer :: gw_dudt(:,:)
        real(kind_phys), pointer :: gw_dvdt(:,:)
        real(kind_phys), pointer :: gw_dtdt(:,:)
        real(kind_phys), pointer :: gw_kdis(:,:)
        real(kind_phys), pointer :: tau_tofd(:,:)
        real(kind_phys), pointer :: tau_mtb(:,:)
        real(kind_phys), pointer :: tau_ogw(:,:)
        real(kind_phys), pointer :: tau_ngw(:,:)
        real(kind_phys), pointer :: zmtb(:)
        real(kind_phys), pointer :: zlwb(:)
        real(kind_phys), pointer :: zogw(:)
        real(kind_phys), pointer :: du3dt_mtb(:,:)
        real(kind_phys), pointer :: du3dt_ogw(:,:)
        real(kind_phys), pointer :: du3dt_tms(:,:)
        real(kind_phys), pointer :: rdxzb(:)

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
        

        call ccpp_field_get(cdata, 'mpi_root', master, ierr=ierr, kind=ckind, index=559)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve mpi_root from CCPP data structure')
            return
        end if
        if (kind(master).ne.ckind) then
            call ccpp_error('Kind mismatch for variable mpi_root')
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
        

        call ccpp_field_get(cdata, 'number_of_statistical_measures_of_subgrid_orography', nvoro, ierr=ierr, kind=ckind, index=581)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve number_of_statistical_measures_of_subgrid_orography from CCPP data structure')
            return
        end if
        if (kind(nvoro).ne.ckind) then
            call ccpp_error('Kind mismatch for variable number_of_statistical_measures_of_subgrid_orography')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'time_step_for_physics', dtp, ierr=ierr, kind=ckind, index=793)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve time_step_for_physics from CCPP data structure')
            return
        end if
        if (kind(dtp).ne.ckind) then
            call ccpp_error('Kind mismatch for variable time_step_for_physics')
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
        

        call ccpp_field_get(cdata, 'number_of_equatorial_longitude_points', lonr, ierr=ierr, kind=ckind, index=576)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve number_of_equatorial_longitude_points from CCPP data structure')
            return
        end if
        if (kind(lonr).ne.ckind) then
            call ccpp_error('Kind mismatch for variable number_of_equatorial_longitude_points')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'number_of_statistical_measures_of_subgrid_orography', nvoro, ierr=ierr, kind=ckind, index=581)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve number_of_statistical_measures_of_subgrid_orography from CCPP data structure')
            return
        end if
        if (kind(nvoro).ne.ckind) then
            call ccpp_error('Kind mismatch for variable number_of_statistical_measures_of_subgrid_orography')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'turb_oro_form_drag_flag', do_tofd, ierr=ierr, kind=ckind, index=806)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve turb_oro_form_drag_flag from CCPP data structure')
            return
        end if
        if (kind(do_tofd).ne.ckind) then
            call ccpp_error('Kind mismatch for variable turb_oro_form_drag_flag')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'multiplication_factors_for_mountain_blocking_and_orographic_gravity_wave_drag', cdmbgwd, ierr=ierr, dims=cdims, kind=ckind, index=562)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve multiplication_factors_for_mountain_blocking_and_orographic_gravity_wave_drag from CCPP data structure')
            return
        end if
        if (kind(cdmbgwd).ne.ckind) then
            call ccpp_error('Kind mismatch for variable multiplication_factors_for_mountain_blocking_and_orographic_gravity_wave_drag')
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
        

        call ccpp_field_get(cdata, 'latitude_degree', xlat_d, ierr=ierr, dims=cdims, kind=ckind, index=461)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve latitude_degree from CCPP data structure')
            return
        end if
        if (kind(xlat_d).ne.ckind) then
            call ccpp_error('Kind mismatch for variable latitude_degree')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

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
        

        call ccpp_field_get(cdata, 'GFS_statein_type', cptr, ierr=ierr, kind=ckind, index=12)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve GFS_statein_type from CCPP data structure')
            return
        end if
        if (ckind.ne.CCPP_GENERIC_KIND) then
            call ccpp_error('Kind mismatch for variable GFS_statein_type')
            ierr = 1
            return
        end if
#endif
        call c_f_pointer(cptr, Statein)

        call ccpp_field_get(cdata, 'GFS_sfcprop_type', cptr, ierr=ierr, kind=ckind, index=10)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve GFS_sfcprop_type from CCPP data structure')
            return
        end if
        if (ckind.ne.CCPP_GENERIC_KIND) then
            call ccpp_error('Kind mismatch for variable GFS_sfcprop_type')
            ierr = 1
            return
        end if
#endif
        call c_f_pointer(cptr, Sfcprop)

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
        

        call ccpp_field_get(cdata, 'statistical_measures_of_subgrid_orography', orostat, ierr=ierr, dims=cdims, kind=ckind, index=670)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve statistical_measures_of_subgrid_orography from CCPP data structure')
            return
        end if
        if (kind(orostat).ne.ckind) then
            call ccpp_error('Kind mismatch for variable statistical_measures_of_subgrid_orography')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'vertical_index_at_top_of_atmosphere_boundary_layer', kpbl, ierr=ierr, dims=cdims, kind=ckind, index=822)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve vertical_index_at_top_of_atmosphere_boundary_layer from CCPP data structure')
            return
        end if
        if (kind(kpbl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable vertical_index_at_top_of_atmosphere_boundary_layer')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'instantaneous_x_stress_due_to_gravity_wave_drag', dusfcg, ierr=ierr, dims=cdims, kind=ckind, index=445)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve instantaneous_x_stress_due_to_gravity_wave_drag from CCPP data structure')
            return
        end if
        if (kind(dusfcg).ne.ckind) then
            call ccpp_error('Kind mismatch for variable instantaneous_x_stress_due_to_gravity_wave_drag')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'instantaneous_y_stress_due_to_gravity_wave_drag', dvsfcg, ierr=ierr, dims=cdims, kind=ckind, index=447)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve instantaneous_y_stress_due_to_gravity_wave_drag from CCPP data structure')
            return
        end if
        if (kind(dvsfcg).ne.ckind) then
            call ccpp_error('Kind mismatch for variable instantaneous_y_stress_due_to_gravity_wave_drag')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'tendency_of_x_wind_due_to_ugwp', gw_dudt, ierr=ierr, dims=cdims, kind=ckind, index=783)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve tendency_of_x_wind_due_to_ugwp from CCPP data structure')
            return
        end if
        if (kind(gw_dudt).ne.ckind) then
            call ccpp_error('Kind mismatch for variable tendency_of_x_wind_due_to_ugwp')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'tendency_of_y_wind_due_to_ugwp', gw_dvdt, ierr=ierr, dims=cdims, kind=ckind, index=786)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve tendency_of_y_wind_due_to_ugwp from CCPP data structure')
            return
        end if
        if (kind(gw_dvdt).ne.ckind) then
            call ccpp_error('Kind mismatch for variable tendency_of_y_wind_due_to_ugwp')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'tendency_of_air_temperature_due_to_ugwp', gw_dtdt, ierr=ierr, dims=cdims, kind=ckind, index=764)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve tendency_of_air_temperature_due_to_ugwp from CCPP data structure')
            return
        end if
        if (kind(gw_dtdt).ne.ckind) then
            call ccpp_error('Kind mismatch for variable tendency_of_air_temperature_due_to_ugwp')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'eddy_mixing_due_to_ugwp', gw_kdis, ierr=ierr, dims=cdims, kind=ckind, index=237)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve eddy_mixing_due_to_ugwp from CCPP data structure')
            return
        end if
        if (kind(gw_kdis).ne.ckind) then
            call ccpp_error('Kind mismatch for variable eddy_mixing_due_to_ugwp')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'instantaneous_momentum_flux_due_to_turbulent_orographic_form_drag', tau_tofd, ierr=ierr, dims=cdims, kind=ckind, index=413)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve instantaneous_momentum_flux_due_to_turbulent_orographic_form_drag from CCPP data structure')
            return
        end if
        if (kind(tau_tofd).ne.ckind) then
            call ccpp_error('Kind mismatch for variable instantaneous_momentum_flux_due_to_turbulent_orographic_form_drag')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'instantaneous_momentum_flux_due_to_mountain_blocking_drag', tau_mtb, ierr=ierr, dims=cdims, kind=ckind, index=410)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve instantaneous_momentum_flux_due_to_mountain_blocking_drag from CCPP data structure')
            return
        end if
        if (kind(tau_mtb).ne.ckind) then
            call ccpp_error('Kind mismatch for variable instantaneous_momentum_flux_due_to_mountain_blocking_drag')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'instantaneous_momentum_flux_due_to_orographic_gravity_wave_drag', tau_ogw, ierr=ierr, dims=cdims, kind=ckind, index=412)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve instantaneous_momentum_flux_due_to_orographic_gravity_wave_drag from CCPP data structure')
            return
        end if
        if (kind(tau_ogw).ne.ckind) then
            call ccpp_error('Kind mismatch for variable instantaneous_momentum_flux_due_to_orographic_gravity_wave_drag')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'instantaneous_momentum_flux_due_to_nonstationary_gravity_wave', tau_ngw, ierr=ierr, dims=cdims, kind=ckind, index=411)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve instantaneous_momentum_flux_due_to_nonstationary_gravity_wave from CCPP data structure')
            return
        end if
        if (kind(tau_ngw).ne.ckind) then
            call ccpp_error('Kind mismatch for variable instantaneous_momentum_flux_due_to_nonstationary_gravity_wave')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'height_of_mountain_blocking', zmtb, ierr=ierr, dims=cdims, kind=ckind, index=363)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve height_of_mountain_blocking from CCPP data structure')
            return
        end if
        if (kind(zmtb).ne.ckind) then
            call ccpp_error('Kind mismatch for variable height_of_mountain_blocking')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'height_of_low_level_wave_breaking', zlwb, ierr=ierr, dims=cdims, kind=ckind, index=362)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve height_of_low_level_wave_breaking from CCPP data structure')
            return
        end if
        if (kind(zlwb).ne.ckind) then
            call ccpp_error('Kind mismatch for variable height_of_low_level_wave_breaking')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'height_of_launch_level_of_orographic_gravity_wave', zogw, ierr=ierr, dims=cdims, kind=ckind, index=361)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve height_of_launch_level_of_orographic_gravity_wave from CCPP data structure')
            return
        end if
        if (kind(zogw).ne.ckind) then
            call ccpp_error('Kind mismatch for variable height_of_launch_level_of_orographic_gravity_wave')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'instantaneous_change_in_x_wind_due_to_mountain_blocking_drag', du3dt_mtb, ierr=ierr, dims=cdims, kind=ckind, index=405)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve instantaneous_change_in_x_wind_due_to_mountain_blocking_drag from CCPP data structure')
            return
        end if
        if (kind(du3dt_mtb).ne.ckind) then
            call ccpp_error('Kind mismatch for variable instantaneous_change_in_x_wind_due_to_mountain_blocking_drag')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'instantaneous_change_in_x_wind_due_to_orographic_gravity_wave_drag', du3dt_ogw, ierr=ierr, dims=cdims, kind=ckind, index=406)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve instantaneous_change_in_x_wind_due_to_orographic_gravity_wave_drag from CCPP data structure')
            return
        end if
        if (kind(du3dt_ogw).ne.ckind) then
            call ccpp_error('Kind mismatch for variable instantaneous_change_in_x_wind_due_to_orographic_gravity_wave_drag')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'instantaneous_change_in_x_wind_due_to_turbulent_orographic_form_drag', du3dt_tms, ierr=ierr, dims=cdims, kind=ckind, index=407)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve instantaneous_change_in_x_wind_due_to_turbulent_orographic_form_drag from CCPP data structure')
            return
        end if
        if (kind(du3dt_tms).ne.ckind) then
            call ccpp_error('Kind mismatch for variable instantaneous_change_in_x_wind_due_to_turbulent_orographic_form_drag')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'level_of_dividing_streamline', rdxzb, ierr=ierr, dims=cdims, kind=ckind, index=465)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve level_of_dividing_streamline from CCPP data structure')
            return
        end if
        if (kind(rdxzb).ne.ckind) then
            call ccpp_error('Kind mismatch for variable level_of_dividing_streamline')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call cires_ugwp_run(Model=Model,me=me,master=master,im=im,levs=levs,ntrac=ntrac,nvoro=nvoro, &
                  dtp=dtp,kdt=kdt,lonr=lonr,nvoro=nvoro,do_tofd=do_tofd,cdmbgwd=cdmbgwd,xlat=xlat, &
                  xlat_d=xlat_d,sinlat=sinlat,coslat=coslat,Statein=Statein,Sfcprop=Sfcprop, &
                  del=del,orostat=orostat,kpbl=kpbl,dusfcg=dusfcg,dvsfcg=dvsfcg,gw_dudt=gw_dudt, &
                  gw_dvdt=gw_dvdt,gw_dtdt=gw_dtdt,gw_kdis=gw_kdis,tau_tofd=tau_tofd,tau_mtb=tau_mtb, &
                  tau_ogw=tau_ogw,tau_ngw=tau_ngw,zmtb=zmtb,zlwb=zlwb,zogw=zogw,du3dt_mtb=du3dt_mtb, &
                  du3dt_ogw=du3dt_ogw,du3dt_tms=du3dt_tms,rdxzb=rdxzb,errmsg=cdata%errmsg, &
                  errflg=cdata%errflg)
        ierr=cdata%errflg

    end function cires_ugwp_run_cap

    function cires_ugwp_init_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind
        integer, pointer :: me
        integer, pointer :: master
        integer, pointer :: nlunit
        integer, pointer :: logunit
        character(len=64), pointer :: fn_nml2
        integer, pointer :: lonr
        integer, pointer :: latr
        integer, pointer :: levs
        real(kind_phys), pointer :: ak(:)
        real(kind_phys), pointer :: bk(:)
        real(kind_phys), pointer :: dtp
        real(kind_phys), pointer :: cdmvgwd(:)
        real(kind_phys), pointer :: cgwf(:)

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
        

        call ccpp_field_get(cdata, 'mpi_root', master, ierr=ierr, kind=ckind, index=559)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve mpi_root from CCPP data structure')
            return
        end if
        if (kind(master).ne.ckind) then
            call ccpp_error('Kind mismatch for variable mpi_root')
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
        

        call ccpp_field_get(cdata, 'iounit_log', logunit, ierr=ierr, kind=ckind, index=450)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve iounit_log from CCPP data structure')
            return
        end if
        if (kind(logunit).ne.ckind) then
            call ccpp_error('Kind mismatch for variable iounit_log')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'namelist_filename', fn_nml2, ierr=ierr, kind=ckind, index=563)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve namelist_filename from CCPP data structure')
            return
        end if
        if (kind(fn_nml2).ne.ckind) then
            call ccpp_error('Kind mismatch for variable namelist_filename')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'number_of_equatorial_longitude_points', lonr, ierr=ierr, kind=ckind, index=576)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve number_of_equatorial_longitude_points from CCPP data structure')
            return
        end if
        if (kind(lonr).ne.ckind) then
            call ccpp_error('Kind mismatch for variable number_of_equatorial_longitude_points')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'number_of_laitude_points', latr, ierr=ierr, kind=ckind, index=579)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve number_of_laitude_points from CCPP data structure')
            return
        end if
        if (kind(latr).ne.ckind) then
            call ccpp_error('Kind mismatch for variable number_of_laitude_points')
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
        

        call ccpp_field_get(cdata, 'a_parameter_of_the_hybrid_coordinate', ak, ierr=ierr, dims=cdims, kind=ckind, index=20)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve a_parameter_of_the_hybrid_coordinate from CCPP data structure')
            return
        end if
        if (kind(ak).ne.ckind) then
            call ccpp_error('Kind mismatch for variable a_parameter_of_the_hybrid_coordinate')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'b_parameter_of_the_hybrid_coordinate', bk, ierr=ierr, dims=cdims, kind=ckind, index=76)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve b_parameter_of_the_hybrid_coordinate from CCPP data structure')
            return
        end if
        if (kind(bk).ne.ckind) then
            call ccpp_error('Kind mismatch for variable b_parameter_of_the_hybrid_coordinate')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'time_step_for_physics', dtp, ierr=ierr, kind=ckind, index=793)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve time_step_for_physics from CCPP data structure')
            return
        end if
        if (kind(dtp).ne.ckind) then
            call ccpp_error('Kind mismatch for variable time_step_for_physics')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'multiplication_factors_for_mountain_blocking_and_orographic_gravity_wave_drag', cdmvgwd, ierr=ierr, dims=cdims, kind=ckind, index=562)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve multiplication_factors_for_mountain_blocking_and_orographic_gravity_wave_drag from CCPP data structure')
            return
        end if
        if (kind(cdmvgwd).ne.ckind) then
            call ccpp_error('Kind mismatch for variable multiplication_factors_for_mountain_blocking_and_orographic_gravity_wave_drag')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'multiplication_factors_for_convective_gravity_wave_drag', cgwf, ierr=ierr, dims=cdims, kind=ckind, index=561)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve multiplication_factors_for_convective_gravity_wave_drag from CCPP data structure')
            return
        end if
        if (kind(cgwf).ne.ckind) then
            call ccpp_error('Kind mismatch for variable multiplication_factors_for_convective_gravity_wave_drag')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call cires_ugwp_init(me=me,master=master,nlunit=nlunit,logunit=logunit,fn_nml2=fn_nml2,lonr=lonr, &
                  latr=latr,levs=levs,ak=ak,bk=bk,dtp=dtp,cdmvgwd=cdmvgwd,cgwf=cgwf,errmsg=cdata%errmsg, &
                  errflg=cdata%errflg)
        ierr=cdata%errflg

    end function cires_ugwp_init_cap

    function cires_ugwp_finalize_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call cires_ugwp_finalize(errmsg=cdata%errmsg,errflg=cdata%errflg)
        ierr=cdata%errflg

    end function cires_ugwp_finalize_cap
end module cires_ugwp_cap


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
!! @brief Auto-generated cap module for the gwdps scheme
!!
!
module gwdps_cap

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_f_pointer, c_ptr, c_int32_t
    use            :: ccpp_types,                                      &
                      only: ccpp_t, CCPP_GENERIC_KIND
    use            :: ccpp_fields,                                     &
                      only: ccpp_field_get
    use            :: ccpp_errors,                                     &
                      only: ccpp_error, ccpp_debug
    use            :: gwdps, &
                      only: gwdps_init,gwdps_run,gwdps_finalize
    ! Other modules required, e.g. type definitions
    use machine, only: kind_phys

    implicit none

    private
    public :: gwdps_init_cap,gwdps_run_cap,gwdps_finalize_cap

    contains


    function gwdps_init_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call gwdps_init()
        

    end function gwdps_init_cap

    function gwdps_run_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind
        integer, pointer :: im
        integer, pointer :: ix
        integer, pointer :: km
        real(kind_phys), pointer :: A(:,:)
        real(kind_phys), pointer :: B(:,:)
        real(kind_phys), pointer :: C(:,:)
        real(kind_phys), pointer :: u1(:,:)
        real(kind_phys), pointer :: v1(:,:)
        real(kind_phys), pointer :: t1(:,:)
        real(kind_phys), pointer :: q1(:,:)
        integer, pointer :: kpbl(:)
        real(kind_phys), pointer :: prsi(:,:)
        real(kind_phys), pointer :: del(:,:)
        real(kind_phys), pointer :: prsl(:,:)
        real(kind_phys), pointer :: prslk(:,:)
        real(kind_phys), pointer :: phii(:,:)
        real(kind_phys), pointer :: phil(:,:)
        real(kind_phys), pointer :: deltim
        integer, pointer :: kdt
        real(kind_phys), pointer :: hprime(:)
        real(kind_phys), pointer :: oc(:)
        real(kind_phys), pointer :: oa4(:,:)
        real(kind_phys), pointer :: clx4(:,:)
        real(kind_phys), pointer :: theta(:)
        real(kind_phys), pointer :: sigma(:)
        real(kind_phys), pointer :: gamma(:)
        real(kind_phys), pointer :: elvmax(:)
        real(kind_phys), pointer :: dusfc(:)
        real(kind_phys), pointer :: dvsfc(:)
        real(kind_phys), pointer :: g
        real(kind_phys), pointer :: cp
        real(kind_phys), pointer :: rd
        real(kind_phys), pointer :: rv
        integer, pointer :: imx
        integer, pointer :: nmtvr
        real(kind_phys), pointer :: cdmbgwd(:)
        integer, pointer :: me
        logical, pointer :: lprnt
        integer, pointer :: ipr
        real(kind_phys), pointer :: rdxzb(:)

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
        

        call ccpp_field_get(cdata, 'tendency_of_y_wind_due_to_model_physics', A, ierr=ierr, dims=cdims, kind=ckind, index=785)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve tendency_of_y_wind_due_to_model_physics from CCPP data structure')
            return
        end if
        if (kind(A).ne.ckind) then
            call ccpp_error('Kind mismatch for variable tendency_of_y_wind_due_to_model_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'tendency_of_x_wind_due_to_model_physics', B, ierr=ierr, dims=cdims, kind=ckind, index=782)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve tendency_of_x_wind_due_to_model_physics from CCPP data structure')
            return
        end if
        if (kind(B).ne.ckind) then
            call ccpp_error('Kind mismatch for variable tendency_of_x_wind_due_to_model_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'tendency_of_air_temperature_due_to_model_physics', C, ierr=ierr, dims=cdims, kind=ckind, index=758)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve tendency_of_air_temperature_due_to_model_physics from CCPP data structure')
            return
        end if
        if (kind(C).ne.ckind) then
            call ccpp_error('Kind mismatch for variable tendency_of_air_temperature_due_to_model_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'x_wind', u1, ierr=ierr, dims=cdims, kind=ckind, index=871)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve x_wind from CCPP data structure')
            return
        end if
        if (kind(u1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable x_wind')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'y_wind', v1, ierr=ierr, dims=cdims, kind=ckind, index=878)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve y_wind from CCPP data structure')
            return
        end if
        if (kind(v1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable y_wind')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'air_temperature', t1, ierr=ierr, dims=cdims, kind=ckind, index=50)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve air_temperature from CCPP data structure')
            return
        end if
        if (kind(t1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable air_temperature')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'water_vapor_specific_humidity', q1, ierr=ierr, dims=cdims, kind=ckind, index=852)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve water_vapor_specific_humidity from CCPP data structure')
            return
        end if
        if (kind(q1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable water_vapor_specific_humidity')
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
        

        call ccpp_field_get(cdata, 'air_pressure_at_interface', prsi, ierr=ierr, dims=cdims, kind=ckind, index=45)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve air_pressure_at_interface from CCPP data structure')
            return
        end if
        if (kind(prsi).ne.ckind) then
            call ccpp_error('Kind mismatch for variable air_pressure_at_interface')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

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
        

        call ccpp_field_get(cdata, 'air_pressure', prsl, ierr=ierr, dims=cdims, kind=ckind, index=44)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve air_pressure from CCPP data structure')
            return
        end if
        if (kind(prsl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable air_pressure')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'dimensionless_exner_function_at_model_layers', prslk, ierr=ierr, dims=cdims, kind=ckind, index=222)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve dimensionless_exner_function_at_model_layers from CCPP data structure')
            return
        end if
        if (kind(prslk).ne.ckind) then
            call ccpp_error('Kind mismatch for variable dimensionless_exner_function_at_model_layers')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'geopotential_at_interface', phii, ierr=ierr, dims=cdims, kind=ckind, index=349)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve geopotential_at_interface from CCPP data structure')
            return
        end if
        if (kind(phii).ne.ckind) then
            call ccpp_error('Kind mismatch for variable geopotential_at_interface')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'geopotential', phil, ierr=ierr, dims=cdims, kind=ckind, index=348)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve geopotential from CCPP data structure')
            return
        end if
        if (kind(phil).ne.ckind) then
            call ccpp_error('Kind mismatch for variable geopotential')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'time_step_for_physics', deltim, ierr=ierr, kind=ckind, index=793)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve time_step_for_physics from CCPP data structure')
            return
        end if
        if (kind(deltim).ne.ckind) then
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
        

        call ccpp_field_get(cdata, 'standard_deviation_of_subgrid_orography', hprime, ierr=ierr, dims=cdims, kind=ckind, index=669)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve standard_deviation_of_subgrid_orography from CCPP data structure')
            return
        end if
        if (kind(hprime).ne.ckind) then
            call ccpp_error('Kind mismatch for variable standard_deviation_of_subgrid_orography')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'convexity_of_subgrid_orography', oc, ierr=ierr, dims=cdims, kind=ckind, index=133)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve convexity_of_subgrid_orography from CCPP data structure')
            return
        end if
        if (kind(oc).ne.ckind) then
            call ccpp_error('Kind mismatch for variable convexity_of_subgrid_orography')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'asymmetry_of_subgrid_orography', oa4, ierr=ierr, dims=cdims, kind=ckind, index=65)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve asymmetry_of_subgrid_orography from CCPP data structure')
            return
        end if
        if (kind(oa4).ne.ckind) then
            call ccpp_error('Kind mismatch for variable asymmetry_of_subgrid_orography')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'fraction_of_grid_box_with_subgrid_orography_higher_than_critical_height', clx4, ierr=ierr, dims=cdims, kind=ckind, index=342)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve fraction_of_grid_box_with_subgrid_orography_higher_than_critical_height from CCPP data structure')
            return
        end if
        if (kind(clx4).ne.ckind) then
            call ccpp_error('Kind mismatch for variable fraction_of_grid_box_with_subgrid_orography_higher_than_critical_height')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'angle_from_east_of_maximum_subgrid_orographic_variations', theta, ierr=ierr, dims=cdims, kind=ckind, index=60)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve angle_from_east_of_maximum_subgrid_orographic_variations from CCPP data structure')
            return
        end if
        if (kind(theta).ne.ckind) then
            call ccpp_error('Kind mismatch for variable angle_from_east_of_maximum_subgrid_orographic_variations')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'slope_of_subgrid_orography', sigma, ierr=ierr, dims=cdims, kind=ckind, index=647)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve slope_of_subgrid_orography from CCPP data structure')
            return
        end if
        if (kind(sigma).ne.ckind) then
            call ccpp_error('Kind mismatch for variable slope_of_subgrid_orography')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'anisotropy_of_subgrid_orography', gamma, ierr=ierr, dims=cdims, kind=ckind, index=61)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve anisotropy_of_subgrid_orography from CCPP data structure')
            return
        end if
        if (kind(gamma).ne.ckind) then
            call ccpp_error('Kind mismatch for variable anisotropy_of_subgrid_orography')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'maximum_subgrid_orography', elvmax, ierr=ierr, dims=cdims, kind=ckind, index=507)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve maximum_subgrid_orography from CCPP data structure')
            return
        end if
        if (kind(elvmax).ne.ckind) then
            call ccpp_error('Kind mismatch for variable maximum_subgrid_orography')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'instantaneous_x_stress_due_to_gravity_wave_drag', dusfc, ierr=ierr, dims=cdims, kind=ckind, index=445)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve instantaneous_x_stress_due_to_gravity_wave_drag from CCPP data structure')
            return
        end if
        if (kind(dusfc).ne.ckind) then
            call ccpp_error('Kind mismatch for variable instantaneous_x_stress_due_to_gravity_wave_drag')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'instantaneous_y_stress_due_to_gravity_wave_drag', dvsfc, ierr=ierr, dims=cdims, kind=ckind, index=447)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve instantaneous_y_stress_due_to_gravity_wave_drag from CCPP data structure')
            return
        end if
        if (kind(dvsfc).ne.ckind) then
            call ccpp_error('Kind mismatch for variable instantaneous_y_stress_due_to_gravity_wave_drag')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'gravitational_acceleration', g, ierr=ierr, kind=ckind, index=355)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve gravitational_acceleration from CCPP data structure')
            return
        end if
        if (kind(g).ne.ckind) then
            call ccpp_error('Kind mismatch for variable gravitational_acceleration')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'specific_heat_of_dry_air_at_constant_pressure', cp, ierr=ierr, kind=ckind, index=664)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve specific_heat_of_dry_air_at_constant_pressure from CCPP data structure')
            return
        end if
        if (kind(cp).ne.ckind) then
            call ccpp_error('Kind mismatch for variable specific_heat_of_dry_air_at_constant_pressure')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'gas_constant_dry_air', rd, ierr=ierr, kind=ckind, index=346)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve gas_constant_dry_air from CCPP data structure')
            return
        end if
        if (kind(rd).ne.ckind) then
            call ccpp_error('Kind mismatch for variable gas_constant_dry_air')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'gas_constant_water_vapor', rv, ierr=ierr, kind=ckind, index=347)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve gas_constant_water_vapor from CCPP data structure')
            return
        end if
        if (kind(rv).ne.ckind) then
            call ccpp_error('Kind mismatch for variable gas_constant_water_vapor')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'number_of_equatorial_longitude_points', imx, ierr=ierr, kind=ckind, index=576)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve number_of_equatorial_longitude_points from CCPP data structure')
            return
        end if
        if (kind(imx).ne.ckind) then
            call ccpp_error('Kind mismatch for variable number_of_equatorial_longitude_points')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'number_of_statistical_measures_of_subgrid_orography', nmtvr, ierr=ierr, kind=ckind, index=581)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve number_of_statistical_measures_of_subgrid_orography from CCPP data structure')
            return
        end if
        if (kind(nmtvr).ne.ckind) then
            call ccpp_error('Kind mismatch for variable number_of_statistical_measures_of_subgrid_orography')
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
        

        call gwdps_run(im=im,ix=ix,km=km,A=A,B=B,C=C,u1=u1,v1=v1,t1=t1,q1=q1,kpbl=kpbl,prsi=prsi, &
                  del=del,prsl=prsl,prslk=prslk,phii=phii,phil=phil,deltim=deltim,kdt=kdt, &
                  hprime=hprime,oc=oc,oa4=oa4,clx4=clx4,theta=theta,sigma=sigma,gamma=gamma, &
                  elvmax=elvmax,dusfc=dusfc,dvsfc=dvsfc,g=g,cp=cp,rd=rd,rv=rv,imx=imx,nmtvr=nmtvr, &
                  cdmbgwd=cdmbgwd,me=me,lprnt=lprnt,ipr=ipr,rdxzb=rdxzb,errmsg=cdata%errmsg, &
                  errflg=cdata%errflg)
        ierr=cdata%errflg

    end function gwdps_run_cap

    function gwdps_finalize_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call gwdps_finalize()
        

    end function gwdps_finalize_cap
end module gwdps_cap

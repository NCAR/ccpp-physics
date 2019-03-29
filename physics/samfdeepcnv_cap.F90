
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
!! @brief Auto-generated cap module for the samfdeepcnv scheme
!!
!
module samfdeepcnv_cap

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_f_pointer, c_ptr, c_int32_t
    use            :: ccpp_types,                                      &
                      only: ccpp_t, CCPP_GENERIC_KIND
    use            :: ccpp_fields,                                     &
                      only: ccpp_field_get
    use            :: ccpp_errors,                                     &
                      only: ccpp_error, ccpp_debug
    use            :: samfdeepcnv, &
                      only: samfdeepcnv_run,samfdeepcnv_init,samfdeepcnv_finalize
    ! Other modules required, e.g. type definitions
    use machine, only: kind_phys

    implicit none

    private
    public :: samfdeepcnv_run_cap,samfdeepcnv_init_cap,samfdeepcnv_finalize_cap

    contains


    function samfdeepcnv_run_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind
        integer, pointer :: im
        integer, pointer :: ix
        integer, pointer :: km
        real(kind_phys), pointer :: cliq
        real(kind_phys), pointer :: cp
        real(kind_phys), pointer :: cvap
        real(kind_phys), pointer :: eps
        real(kind_phys), pointer :: epsm1
        real(kind_phys), pointer :: fv
        real(kind_phys), pointer :: grav
        real(kind_phys), pointer :: hvap
        real(kind_phys), pointer :: rd
        real(kind_phys), pointer :: rv
        real(kind_phys), pointer :: t0c
        real(kind_phys), pointer :: delt
        integer, pointer :: ntk
        integer, pointer :: ntr
        real(kind_phys), pointer :: delp(:,:)
        real(kind_phys), pointer :: prslp(:,:)
        real(kind_phys), pointer :: psp(:)
        real(kind_phys), pointer :: phil(:,:)
        real(kind_phys), pointer :: qtr(:,:,:)
        real(kind_phys), pointer :: q1(:,:)
        real(kind_phys), pointer :: t1(:,:)
        real(kind_phys), pointer :: u1(:,:)
        real(kind_phys), pointer :: v1(:,:)
        real(kind_phys), pointer :: cldwrk(:)
        real(kind_phys), pointer :: rn(:)
        integer, pointer :: kbot(:)
        integer, pointer :: ktop(:)
        integer, pointer :: kcnv(:)
        integer, pointer :: islimsk(:)
        real(kind_phys), pointer :: garea(:)
        real(kind_phys), pointer :: dot(:,:)
        integer, pointer :: ncloud
        real(kind_phys), pointer :: ud_mf(:,:)
        real(kind_phys), pointer :: dd_mf(:,:)
        real(kind_phys), pointer :: dt_mf(:,:)
        real(kind_phys), pointer :: cnvw(:,:)
        real(kind_phys), pointer :: cnvc(:,:)
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
        integer, pointer :: mp_phys_mg
        real(kind_phys), pointer :: clam
        real(kind_phys), pointer :: c0s
        real(kind_phys), pointer :: c1
        real(kind_phys), pointer :: betal
        real(kind_phys), pointer :: betas
        real(kind_phys), pointer :: evfact
        real(kind_phys), pointer :: evfactl
        real(kind_phys), pointer :: pgcon
        real(kind_phys), pointer :: asolfac

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
        

        call ccpp_field_get(cdata, 'specific_heat_of_liquid_water_at_constant_pressure', cliq, ierr=ierr, kind=ckind, index=665)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve specific_heat_of_liquid_water_at_constant_pressure from CCPP data structure')
            return
        end if
        if (kind(cliq).ne.ckind) then
            call ccpp_error('Kind mismatch for variable specific_heat_of_liquid_water_at_constant_pressure')
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
        

        call ccpp_field_get(cdata, 'specific_heat_of_water_vapor_at_constant_pressure', cvap, ierr=ierr, kind=ckind, index=666)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve specific_heat_of_water_vapor_at_constant_pressure from CCPP data structure')
            return
        end if
        if (kind(cvap).ne.ckind) then
            call ccpp_error('Kind mismatch for variable specific_heat_of_water_vapor_at_constant_pressure')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'ratio_of_dry_air_to_water_vapor_gas_constants', eps, ierr=ierr, kind=ckind, index=621)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve ratio_of_dry_air_to_water_vapor_gas_constants from CCPP data structure')
            return
        end if
        if (kind(eps).ne.ckind) then
            call ccpp_error('Kind mismatch for variable ratio_of_dry_air_to_water_vapor_gas_constants')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'ratio_of_dry_air_to_water_vapor_gas_constants_minus_one', epsm1, ierr=ierr, kind=ckind, index=622)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve ratio_of_dry_air_to_water_vapor_gas_constants_minus_one from CCPP data structure')
            return
        end if
        if (kind(epsm1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable ratio_of_dry_air_to_water_vapor_gas_constants_minus_one')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'ratio_of_vapor_to_dry_air_gas_constants_minus_one', fv, ierr=ierr, kind=ckind, index=625)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve ratio_of_vapor_to_dry_air_gas_constants_minus_one from CCPP data structure')
            return
        end if
        if (kind(fv).ne.ckind) then
            call ccpp_error('Kind mismatch for variable ratio_of_vapor_to_dry_air_gas_constants_minus_one')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'gravitational_acceleration', grav, ierr=ierr, kind=ckind, index=355)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve gravitational_acceleration from CCPP data structure')
            return
        end if
        if (kind(grav).ne.ckind) then
            call ccpp_error('Kind mismatch for variable gravitational_acceleration')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'latent_heat_of_vaporization_of_water_at_0C', hvap, ierr=ierr, kind=ckind, index=459)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve latent_heat_of_vaporization_of_water_at_0C from CCPP data structure')
            return
        end if
        if (kind(hvap).ne.ckind) then
            call ccpp_error('Kind mismatch for variable latent_heat_of_vaporization_of_water_at_0C')
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
        

        call ccpp_field_get(cdata, 'temperature_at_zero_celsius', t0c, ierr=ierr, kind=ckind, index=751)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve temperature_at_zero_celsius from CCPP data structure')
            return
        end if
        if (kind(t0c).ne.ckind) then
            call ccpp_error('Kind mismatch for variable temperature_at_zero_celsius')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'time_step_for_physics', delt, ierr=ierr, kind=ckind, index=793)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve time_step_for_physics from CCPP data structure')
            return
        end if
        if (kind(delt).ne.ckind) then
            call ccpp_error('Kind mismatch for variable time_step_for_physics')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'index_for_turbulent_kinetic_energy_convective_transport_tracer', ntk, ierr=ierr, kind=ckind, index=392)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve index_for_turbulent_kinetic_energy_convective_transport_tracer from CCPP data structure')
            return
        end if
        if (kind(ntk).ne.ckind) then
            call ccpp_error('Kind mismatch for variable index_for_turbulent_kinetic_energy_convective_transport_tracer')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'number_of_tracers_for_samf', ntr, ierr=ierr, kind=ckind, index=588)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve number_of_tracers_for_samf from CCPP data structure')
            return
        end if
        if (kind(ntr).ne.ckind) then
            call ccpp_error('Kind mismatch for variable number_of_tracers_for_samf')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'air_pressure_difference_between_midlayers', delp, ierr=ierr, dims=cdims, kind=ckind, index=49)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve air_pressure_difference_between_midlayers from CCPP data structure')
            return
        end if
        if (kind(delp).ne.ckind) then
            call ccpp_error('Kind mismatch for variable air_pressure_difference_between_midlayers')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'air_pressure', prslp, ierr=ierr, dims=cdims, kind=ckind, index=44)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve air_pressure from CCPP data structure')
            return
        end if
        if (kind(prslp).ne.ckind) then
            call ccpp_error('Kind mismatch for variable air_pressure')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_air_pressure', psp, ierr=ierr, dims=cdims, kind=ckind, index=677)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_air_pressure from CCPP data structure')
            return
        end if
        if (kind(psp).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_air_pressure')
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
        

        call ccpp_field_get(cdata, 'convective_transportable_tracers', qtr, ierr=ierr, dims=cdims, kind=ckind, index=130)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve convective_transportable_tracers from CCPP data structure')
            return
        end if
        if (kind(qtr).ne.ckind) then
            call ccpp_error('Kind mismatch for variable convective_transportable_tracers')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'water_vapor_specific_humidity_updated_by_physics', q1, ierr=ierr, dims=cdims, kind=ckind, index=860)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve water_vapor_specific_humidity_updated_by_physics from CCPP data structure')
            return
        end if
        if (kind(q1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable water_vapor_specific_humidity_updated_by_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'air_temperature_updated_by_physics', t1, ierr=ierr, dims=cdims, kind=ckind, index=59)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve air_temperature_updated_by_physics from CCPP data structure')
            return
        end if
        if (kind(t1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable air_temperature_updated_by_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'x_wind_updated_by_physics', u1, ierr=ierr, dims=cdims, kind=ckind, index=877)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve x_wind_updated_by_physics from CCPP data structure')
            return
        end if
        if (kind(u1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable x_wind_updated_by_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'y_wind_updated_by_physics', v1, ierr=ierr, dims=cdims, kind=ckind, index=884)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve y_wind_updated_by_physics from CCPP data structure')
            return
        end if
        if (kind(v1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable y_wind_updated_by_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cloud_work_function', cldwrk, ierr=ierr, dims=cdims, kind=ckind, index=111)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cloud_work_function from CCPP data structure')
            return
        end if
        if (kind(cldwrk).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cloud_work_function')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'lwe_thickness_of_deep_convective_precipitation_amount', rn, ierr=ierr, dims=cdims, kind=ckind, index=479)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve lwe_thickness_of_deep_convective_precipitation_amount from CCPP data structure')
            return
        end if
        if (kind(rn).ne.ckind) then
            call ccpp_error('Kind mismatch for variable lwe_thickness_of_deep_convective_precipitation_amount')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'vertical_index_at_cloud_base', kbot, ierr=ierr, dims=cdims, kind=ckind, index=820)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve vertical_index_at_cloud_base from CCPP data structure')
            return
        end if
        if (kind(kbot).ne.ckind) then
            call ccpp_error('Kind mismatch for variable vertical_index_at_cloud_base')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'vertical_index_at_cloud_top', ktop, ierr=ierr, dims=cdims, kind=ckind, index=821)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve vertical_index_at_cloud_top from CCPP data structure')
            return
        end if
        if (kind(ktop).ne.ckind) then
            call ccpp_error('Kind mismatch for variable vertical_index_at_cloud_top')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

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
        

        call ccpp_field_get(cdata, 'cell_area', garea, ierr=ierr, dims=cdims, kind=ckind, index=83)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cell_area from CCPP data structure')
            return
        end if
        if (kind(garea).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cell_area')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'omega', dot, ierr=ierr, dims=cdims, kind=ckind, index=593)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve omega from CCPP data structure')
            return
        end if
        if (kind(dot).ne.ckind) then
            call ccpp_error('Kind mismatch for variable omega')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'number_of_hydrometeors', ncloud, ierr=ierr, kind=ckind, index=578)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve number_of_hydrometeors from CCPP data structure')
            return
        end if
        if (kind(ncloud).ne.ckind) then
            call ccpp_error('Kind mismatch for variable number_of_hydrometeors')
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
        

        call ccpp_field_get(cdata, 'convective_cloud_water_mixing_ratio', cnvw, ierr=ierr, dims=cdims, kind=ckind, index=128)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve convective_cloud_water_mixing_ratio from CCPP data structure')
            return
        end if
        if (kind(cnvw).ne.ckind) then
            call ccpp_error('Kind mismatch for variable convective_cloud_water_mixing_ratio')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'convective_cloud_cover', cnvc, ierr=ierr, dims=cdims, kind=ckind, index=123)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve convective_cloud_cover from CCPP data structure')
            return
        end if
        if (kind(cnvc).ne.ckind) then
            call ccpp_error('Kind mismatch for variable convective_cloud_cover')
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
        

        call ccpp_field_get(cdata, 'flag_for_morrison_gettelman_microphysics_scheme', mp_phys_mg, ierr=ierr, kind=ckind, index=297)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_morrison_gettelman_microphysics_scheme from CCPP data structure')
            return
        end if
        if (kind(mp_phys_mg).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_morrison_gettelman_microphysics_scheme')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'entrainment_rate_coefficient_deep_convection', clam, ierr=ierr, kind=ckind, index=254)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve entrainment_rate_coefficient_deep_convection from CCPP data structure')
            return
        end if
        if (kind(clam).ne.ckind) then
            call ccpp_error('Kind mismatch for variable entrainment_rate_coefficient_deep_convection')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'rain_conversion_parameter_deep_convection', c0s, ierr=ierr, kind=ckind, index=614)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve rain_conversion_parameter_deep_convection from CCPP data structure')
            return
        end if
        if (kind(c0s).ne.ckind) then
            call ccpp_error('Kind mismatch for variable rain_conversion_parameter_deep_convection')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'detrainment_conversion_parameter_deep_convection', c1, ierr=ierr, kind=ckind, index=216)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve detrainment_conversion_parameter_deep_convection from CCPP data structure')
            return
        end if
        if (kind(c1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable detrainment_conversion_parameter_deep_convection')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'downdraft_fraction_reaching_surface_over_land_deep_convection', betal, ierr=ierr, kind=ckind, index=233)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve downdraft_fraction_reaching_surface_over_land_deep_convection from CCPP data structure')
            return
        end if
        if (kind(betal).ne.ckind) then
            call ccpp_error('Kind mismatch for variable downdraft_fraction_reaching_surface_over_land_deep_convection')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'downdraft_fraction_reaching_surface_over_ocean_deep_convection', betas, ierr=ierr, kind=ckind, index=234)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve downdraft_fraction_reaching_surface_over_ocean_deep_convection from CCPP data structure')
            return
        end if
        if (kind(betas).ne.ckind) then
            call ccpp_error('Kind mismatch for variable downdraft_fraction_reaching_surface_over_ocean_deep_convection')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'rain_evaporation_coefficient_deep_convection', evfact, ierr=ierr, kind=ckind, index=616)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve rain_evaporation_coefficient_deep_convection from CCPP data structure')
            return
        end if
        if (kind(evfact).ne.ckind) then
            call ccpp_error('Kind mismatch for variable rain_evaporation_coefficient_deep_convection')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'rain_evaporation_coefficient_over_land_deep_convection', evfactl, ierr=ierr, kind=ckind, index=617)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve rain_evaporation_coefficient_over_land_deep_convection from CCPP data structure')
            return
        end if
        if (kind(evfactl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable rain_evaporation_coefficient_over_land_deep_convection')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'momentum_transport_reduction_factor_pgf_deep_convection', pgcon, ierr=ierr, kind=ckind, index=555)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve momentum_transport_reduction_factor_pgf_deep_convection from CCPP data structure')
            return
        end if
        if (kind(pgcon).ne.ckind) then
            call ccpp_error('Kind mismatch for variable momentum_transport_reduction_factor_pgf_deep_convection')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'aerosol_aware_parameter_deep_convection', asolfac, ierr=ierr, kind=ckind, index=35)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve aerosol_aware_parameter_deep_convection from CCPP data structure')
            return
        end if
        if (kind(asolfac).ne.ckind) then
            call ccpp_error('Kind mismatch for variable aerosol_aware_parameter_deep_convection')
            ierr = 1
            return
        end if
#endif
        

        call samfdeepcnv_run(im=im,ix=ix,km=km,cliq=cliq,cp=cp,cvap=cvap,eps=eps,epsm1=epsm1,fv=fv,grav=grav, &
                  hvap=hvap,rd=rd,rv=rv,t0c=t0c,delt=delt,ntk=ntk,ntr=ntr,delp=delp,prslp=prslp, &
                  psp=psp,phil=phil,qtr=qtr,q1=q1,t1=t1,u1=u1,v1=v1,cldwrk=cldwrk,rn=rn,kbot=kbot, &
                  ktop=ktop,kcnv=kcnv,islimsk=islimsk,garea=garea,dot=dot,ncloud=ncloud,ud_mf=ud_mf, &
                  dd_mf=dd_mf,dt_mf=dt_mf,cnvw=cnvw,cnvc=cnvc,qlcn=qlcn,qicn=qicn,w_upi=w_upi, &
                  cf_upi=cf_upi,cnv_mfd=cnv_mfd,cnv_dqldt=cnv_dqldt,clcn=clcn,cnv_fice=cnv_fice, &
                  cnv_ndrop=cnv_ndrop,cnv_nice=cnv_nice,mp_phys=mp_phys,mp_phys_mg=mp_phys_mg, &
                  clam=clam,c0s=c0s,c1=c1,betal=betal,betas=betas,evfact=evfact,evfactl=evfactl, &
                  pgcon=pgcon,asolfac=asolfac,errmsg=cdata%errmsg,errflg=cdata%errflg)
        ierr=cdata%errflg

    end function samfdeepcnv_run_cap

    function samfdeepcnv_init_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call samfdeepcnv_init()
        

    end function samfdeepcnv_init_cap

    function samfdeepcnv_finalize_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call samfdeepcnv_finalize()
        

    end function samfdeepcnv_finalize_cap
end module samfdeepcnv_cap

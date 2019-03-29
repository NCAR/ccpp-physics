
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
!! @brief Auto-generated cap module for the cu_gf_driver scheme
!!
!
module cu_gf_driver_cap

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_f_pointer, c_ptr, c_int32_t
    use            :: ccpp_types,                                      &
                      only: ccpp_t, CCPP_GENERIC_KIND
    use            :: ccpp_fields,                                     &
                      only: ccpp_field_get
    use            :: ccpp_errors,                                     &
                      only: ccpp_error, ccpp_debug
    use            :: cu_gf_driver, &
                      only: cu_gf_driver_init,cu_gf_driver_run,cu_gf_driver_finalize
    ! Other modules required, e.g. type definitions
    use machine, only: kind_phys

    implicit none

    private
    public :: cu_gf_driver_init_cap,cu_gf_driver_run_cap,cu_gf_driver_finalize_cap

    contains


    function cu_gf_driver_init_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind
        integer, pointer :: mpirank
        integer, pointer :: mpiroot

        ierr = 0

        call c_f_pointer(ptr, cdata)


        call ccpp_field_get(cdata, 'mpi_rank', mpirank, ierr=ierr, kind=ckind, index=558)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve mpi_rank from CCPP data structure')
            return
        end if
        if (kind(mpirank).ne.ckind) then
            call ccpp_error('Kind mismatch for variable mpi_rank')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'mpi_root', mpiroot, ierr=ierr, kind=ckind, index=559)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve mpi_root from CCPP data structure')
            return
        end if
        if (kind(mpiroot).ne.ckind) then
            call ccpp_error('Kind mismatch for variable mpi_root')
            ierr = 1
            return
        end if
#endif
        

        call cu_gf_driver_init(mpirank=mpirank,mpiroot=mpiroot,errmsg=cdata%errmsg,errflg=cdata%errflg)
        ierr=cdata%errflg

    end function cu_gf_driver_init_cap

    function cu_gf_driver_run_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind
        integer, pointer :: tottracer
        integer, pointer :: ntrac
        real(kind_phys), pointer :: garea(:)
        integer, pointer :: im
        integer, pointer :: ix
        integer, pointer :: km
        real(kind_phys), pointer :: dt
        integer, pointer :: cactiv(:)
        real(kind_phys), pointer :: forcet(:,:)
        real(kind_phys), pointer :: forceqv_spechum(:,:)
        real(kind_phys), pointer :: phil(:,:)
        real(kind_phys), pointer :: raincv(:)
        real(kind_phys), pointer :: qv_spechum(:,:)
        real(kind_phys), pointer :: t(:,:)
        real(kind_phys), pointer :: cld1d(:)
        real(kind_phys), pointer :: us(:,:)
        real(kind_phys), pointer :: vs(:,:)
        real(kind_phys), pointer :: t2di(:,:)
        real(kind_phys), pointer :: w(:,:)
        real(kind_phys), pointer :: qv2di_spechum(:,:)
        real(kind_phys), pointer :: p2di(:,:)
        real(kind_phys), pointer :: psuri(:)
        integer, pointer :: hbot(:)
        integer, pointer :: htop(:)
        integer, pointer :: kcnv(:)
        integer, pointer :: xland(:)
        real(kind_phys), pointer :: hfx2(:)
        real(kind_phys), pointer :: qfx2(:)
        real(kind_phys), pointer :: clw(:,:,:)
        real(kind_phys), pointer :: pbl(:)
        real(kind_phys), pointer :: ud_mf(:,:)
        real(kind_phys), pointer :: dd_mf(:,:)
        real(kind_phys), pointer :: dt_mf(:,:)
        real(kind_phys), pointer :: cnvw_moist(:,:)
        real(kind_phys), pointer :: cnvc(:,:)
        integer, pointer :: imfshalcnv

        ierr = 0

        call c_f_pointer(ptr, cdata)


        call ccpp_field_get(cdata, 'number_of_total_tracers', tottracer, ierr=ierr, kind=ckind, index=583)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve number_of_total_tracers from CCPP data structure')
            return
        end if
        if (kind(tottracer).ne.ckind) then
            call ccpp_error('Kind mismatch for variable number_of_total_tracers')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'number_of_vertical_diffusion_tracers', ntrac, ierr=ierr, kind=ckind, index=590)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve number_of_vertical_diffusion_tracers from CCPP data structure')
            return
        end if
        if (kind(ntrac).ne.ckind) then
            call ccpp_error('Kind mismatch for variable number_of_vertical_diffusion_tracers')
            ierr = 1
            return
        end if
#endif
        

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
        

        call ccpp_field_get(cdata, 'time_step_for_physics', dt, ierr=ierr, kind=ckind, index=793)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve time_step_for_physics from CCPP data structure')
            return
        end if
        if (kind(dt).ne.ckind) then
            call ccpp_error('Kind mismatch for variable time_step_for_physics')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'conv_activity_counter', cactiv, ierr=ierr, dims=cdims, kind=ckind, index=122)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve conv_activity_counter from CCPP data structure')
            return
        end if
        if (kind(cactiv).ne.ckind) then
            call ccpp_error('Kind mismatch for variable conv_activity_counter')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'temperature_tendency_due_to_dynamics', forcet, ierr=ierr, dims=cdims, kind=ckind, index=753)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve temperature_tendency_due_to_dynamics from CCPP data structure')
            return
        end if
        if (kind(forcet).ne.ckind) then
            call ccpp_error('Kind mismatch for variable temperature_tendency_due_to_dynamics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'moisture_tendency_due_to_dynamics', forceqv_spechum, ierr=ierr, dims=cdims, kind=ckind, index=554)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve moisture_tendency_due_to_dynamics from CCPP data structure')
            return
        end if
        if (kind(forceqv_spechum).ne.ckind) then
            call ccpp_error('Kind mismatch for variable moisture_tendency_due_to_dynamics')
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
        

        call ccpp_field_get(cdata, 'lwe_thickness_of_deep_convective_precipitation_amount', raincv, ierr=ierr, dims=cdims, kind=ckind, index=479)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve lwe_thickness_of_deep_convective_precipitation_amount from CCPP data structure')
            return
        end if
        if (kind(raincv).ne.ckind) then
            call ccpp_error('Kind mismatch for variable lwe_thickness_of_deep_convective_precipitation_amount')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'water_vapor_specific_humidity_updated_by_physics', qv_spechum, ierr=ierr, dims=cdims, kind=ckind, index=860)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve water_vapor_specific_humidity_updated_by_physics from CCPP data structure')
            return
        end if
        if (kind(qv_spechum).ne.ckind) then
            call ccpp_error('Kind mismatch for variable water_vapor_specific_humidity_updated_by_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

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
        

        call ccpp_field_get(cdata, 'cloud_work_function', cld1d, ierr=ierr, dims=cdims, kind=ckind, index=111)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cloud_work_function from CCPP data structure')
            return
        end if
        if (kind(cld1d).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cloud_work_function')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'x_wind_updated_by_physics', us, ierr=ierr, dims=cdims, kind=ckind, index=877)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve x_wind_updated_by_physics from CCPP data structure')
            return
        end if
        if (kind(us).ne.ckind) then
            call ccpp_error('Kind mismatch for variable x_wind_updated_by_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'y_wind_updated_by_physics', vs, ierr=ierr, dims=cdims, kind=ckind, index=884)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve y_wind_updated_by_physics from CCPP data structure')
            return
        end if
        if (kind(vs).ne.ckind) then
            call ccpp_error('Kind mismatch for variable y_wind_updated_by_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'air_temperature', t2di, ierr=ierr, dims=cdims, kind=ckind, index=50)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve air_temperature from CCPP data structure')
            return
        end if
        if (kind(t2di).ne.ckind) then
            call ccpp_error('Kind mismatch for variable air_temperature')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'omega', w, ierr=ierr, dims=cdims, kind=ckind, index=593)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve omega from CCPP data structure')
            return
        end if
        if (kind(w).ne.ckind) then
            call ccpp_error('Kind mismatch for variable omega')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'water_vapor_specific_humidity', qv2di_spechum, ierr=ierr, dims=cdims, kind=ckind, index=852)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve water_vapor_specific_humidity from CCPP data structure')
            return
        end if
        if (kind(qv2di_spechum).ne.ckind) then
            call ccpp_error('Kind mismatch for variable water_vapor_specific_humidity')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'air_pressure', p2di, ierr=ierr, dims=cdims, kind=ckind, index=44)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve air_pressure from CCPP data structure')
            return
        end if
        if (kind(p2di).ne.ckind) then
            call ccpp_error('Kind mismatch for variable air_pressure')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_air_pressure', psuri, ierr=ierr, dims=cdims, kind=ckind, index=677)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_air_pressure from CCPP data structure')
            return
        end if
        if (kind(psuri).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_air_pressure')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'vertical_index_at_cloud_base', hbot, ierr=ierr, dims=cdims, kind=ckind, index=820)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve vertical_index_at_cloud_base from CCPP data structure')
            return
        end if
        if (kind(hbot).ne.ckind) then
            call ccpp_error('Kind mismatch for variable vertical_index_at_cloud_base')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'vertical_index_at_cloud_top', htop, ierr=ierr, dims=cdims, kind=ckind, index=821)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve vertical_index_at_cloud_top from CCPP data structure')
            return
        end if
        if (kind(htop).ne.ckind) then
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
        

        call ccpp_field_get(cdata, 'sea_land_ice_mask', xland, ierr=ierr, dims=cdims, kind=ckind, index=631)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve sea_land_ice_mask from CCPP data structure')
            return
        end if
        if (kind(xland).ne.ckind) then
            call ccpp_error('Kind mismatch for variable sea_land_ice_mask')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'kinematic_surface_upward_sensible_heat_flux', hfx2, ierr=ierr, dims=cdims, kind=ckind, index=455)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve kinematic_surface_upward_sensible_heat_flux from CCPP data structure')
            return
        end if
        if (kind(hfx2).ne.ckind) then
            call ccpp_error('Kind mismatch for variable kinematic_surface_upward_sensible_heat_flux')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'kinematic_surface_upward_latent_heat_flux', qfx2, ierr=ierr, dims=cdims, kind=ckind, index=454)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve kinematic_surface_upward_latent_heat_flux from CCPP data structure')
            return
        end if
        if (kind(qfx2).ne.ckind) then
            call ccpp_error('Kind mismatch for variable kinematic_surface_upward_latent_heat_flux')
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
        

        call ccpp_field_get(cdata, 'atmosphere_boundary_layer_thickness', pbl, ierr=ierr, dims=cdims, kind=ckind, index=66)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve atmosphere_boundary_layer_thickness from CCPP data structure')
            return
        end if
        if (kind(pbl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable atmosphere_boundary_layer_thickness')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

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
        

        call ccpp_field_get(cdata, 'convective_cloud_water_mixing_ratio', cnvw_moist, ierr=ierr, dims=cdims, kind=ckind, index=128)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve convective_cloud_water_mixing_ratio from CCPP data structure')
            return
        end if
        if (kind(cnvw_moist).ne.ckind) then
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
        

        call ccpp_field_get(cdata, 'flag_for_mass_flux_shallow_convection_scheme', imfshalcnv, ierr=ierr, kind=ckind, index=291)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_mass_flux_shallow_convection_scheme from CCPP data structure')
            return
        end if
        if (kind(imfshalcnv).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_mass_flux_shallow_convection_scheme')
            ierr = 1
            return
        end if
#endif
        

        call cu_gf_driver_run(tottracer=tottracer,ntrac=ntrac,garea=garea,im=im,ix=ix,km=km,dt=dt,cactiv=cactiv, &
                  forcet=forcet,forceqv_spechum=forceqv_spechum,phil=phil,raincv=raincv,qv_spechum=qv_spechum, &
                  t=t,cld1d=cld1d,us=us,vs=vs,t2di=t2di,w=w,qv2di_spechum=qv2di_spechum,p2di=p2di, &
                  psuri=psuri,hbot=hbot,htop=htop,kcnv=kcnv,xland=xland,hfx2=hfx2,qfx2=qfx2, &
                  clw=clw,pbl=pbl,ud_mf=ud_mf,dd_mf=dd_mf,dt_mf=dt_mf,cnvw_moist=cnvw_moist, &
                  cnvc=cnvc,imfshalcnv=imfshalcnv,errmsg=cdata%errmsg,errflg=cdata%errflg)
        ierr=cdata%errflg

    end function cu_gf_driver_run_cap

    function cu_gf_driver_finalize_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call cu_gf_driver_finalize()
        

    end function cu_gf_driver_finalize_cap
end module cu_gf_driver_cap

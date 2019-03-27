
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
!! @brief Auto-generated cap module for the mp_thompson_hrrr scheme
!!
!
module mp_thompson_hrrr_cap

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_f_pointer, c_ptr, c_int32_t
    use            :: ccpp_types,                                      &
                      only: ccpp_t, CCPP_GENERIC_KIND
    use            :: ccpp_fields,                                     &
                      only: ccpp_field_get
    use            :: ccpp_errors,                                     &
                      only: ccpp_error, ccpp_debug
    use            :: mp_thompson_hrrr, &
                      only: mp_thompson_hrrr_finalize,mp_thompson_hrrr_init,mp_thompson_hrrr_run
    ! Other modules required, e.g. type definitions
    use machine, only: kind_phys

    implicit none

    private
    public :: mp_thompson_hrrr_finalize_cap,mp_thompson_hrrr_init_cap,mp_thompson_hrrr_run_cap

    contains


    function mp_thompson_hrrr_finalize_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call mp_thompson_hrrr_finalize(errmsg=cdata%errmsg,errflg=cdata%errflg)
        ierr=cdata%errflg

    end function mp_thompson_hrrr_finalize_cap

    function mp_thompson_hrrr_init_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind
        integer, pointer :: ncol
        integer, pointer :: nlev
        logical, pointer :: is_aerosol_aware
        real(kind_phys), pointer :: nwfa2d(:)
        real(kind_phys), pointer :: nifa2d(:)
        real(kind_phys), pointer :: nwfa(:,:)
        real(kind_phys), pointer :: nifa(:,:)
        integer, pointer :: mpicomm
        integer, pointer :: mpirank
        integer, pointer :: mpiroot
        integer, pointer :: threads
        integer, pointer :: imp_physics
        integer, pointer :: imp_physics_thompson

        ierr = 0

        call c_f_pointer(ptr, cdata)


        call ccpp_field_get(cdata, 'horizontal_loop_extent', ncol, ierr=ierr, kind=ckind, index=366)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve horizontal_loop_extent from CCPP data structure')
            return
        end if
        if (kind(ncol).ne.ckind) then
            call ccpp_error('Kind mismatch for variable horizontal_loop_extent')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'vertical_dimension', nlev, ierr=ierr, kind=ckind, index=817)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve vertical_dimension from CCPP data structure')
            return
        end if
        if (kind(nlev).ne.ckind) then
            call ccpp_error('Kind mismatch for variable vertical_dimension')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'flag_for_aerosol_physics', is_aerosol_aware, ierr=ierr, kind=ckind, index=271)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_aerosol_physics from CCPP data structure')
            return
        end if
        if (kind(is_aerosol_aware).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_aerosol_physics')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'tendency_of_water_friendly_aerosols_at_surface', nwfa2d, ierr=ierr, dims=cdims, kind=ckind, index=779)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve tendency_of_water_friendly_aerosols_at_surface from CCPP data structure')
            return
        end if
        if (kind(nwfa2d).ne.ckind) then
            call ccpp_error('Kind mismatch for variable tendency_of_water_friendly_aerosols_at_surface')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'tendency_of_ice_friendly_aerosols_at_surface', nifa2d, ierr=ierr, dims=cdims, kind=ckind, index=769)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve tendency_of_ice_friendly_aerosols_at_surface from CCPP data structure')
            return
        end if
        if (kind(nifa2d).ne.ckind) then
            call ccpp_error('Kind mismatch for variable tendency_of_ice_friendly_aerosols_at_surface')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'water_friendly_aerosol_number_concentration', nwfa, ierr=ierr, dims=cdims, kind=ckind, index=849)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve water_friendly_aerosol_number_concentration from CCPP data structure')
            return
        end if
        if (kind(nwfa).ne.ckind) then
            call ccpp_error('Kind mismatch for variable water_friendly_aerosol_number_concentration')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'ice_friendly_aerosol_number_concentration', nifa, ierr=ierr, dims=cdims, kind=ckind, index=368)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve ice_friendly_aerosol_number_concentration from CCPP data structure')
            return
        end if
        if (kind(nifa).ne.ckind) then
            call ccpp_error('Kind mismatch for variable ice_friendly_aerosol_number_concentration')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'mpi_comm', mpicomm, ierr=ierr, kind=ckind, index=557)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve mpi_comm from CCPP data structure')
            return
        end if
        if (kind(mpicomm).ne.ckind) then
            call ccpp_error('Kind mismatch for variable mpi_comm')
            ierr = 1
            return
        end if
#endif
        

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
        

        call ccpp_field_get(cdata, 'omp_threads', threads, ierr=ierr, kind=ckind, index=594)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve omp_threads from CCPP data structure')
            return
        end if
        if (kind(threads).ne.ckind) then
            call ccpp_error('Kind mismatch for variable omp_threads')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'flag_for_microphysics_scheme', imp_physics, ierr=ierr, kind=ckind, index=294)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_microphysics_scheme from CCPP data structure')
            return
        end if
        if (kind(imp_physics).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_microphysics_scheme')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'flag_for_thompson_microphysics_scheme', imp_physics_thompson, ierr=ierr, kind=ckind, index=322)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_thompson_microphysics_scheme from CCPP data structure')
            return
        end if
        if (kind(imp_physics_thompson).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_thompson_microphysics_scheme')
            ierr = 1
            return
        end if
#endif
        

        call mp_thompson_hrrr_init(ncol=ncol,nlev=nlev,is_aerosol_aware=is_aerosol_aware,nwfa2d=nwfa2d,nifa2d=nifa2d, &
                  nwfa=nwfa,nifa=nifa,mpicomm=mpicomm,mpirank=mpirank,mpiroot=mpiroot,threads=threads, &
                  imp_physics=imp_physics,imp_physics_thompson=imp_physics_thompson,errmsg=cdata%errmsg, &
                  errflg=cdata%errflg)
        ierr=cdata%errflg

    end function mp_thompson_hrrr_init_cap

    function mp_thompson_hrrr_run_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind
        integer, pointer :: ncol
        integer, pointer :: nlev
        real(kind_phys), pointer :: con_g
        real(kind_phys), pointer :: con_rd
        real(kind_phys), pointer :: spechum(:,:)
        real(kind_phys), pointer :: qc(:,:)
        real(kind_phys), pointer :: qr(:,:)
        real(kind_phys), pointer :: qi(:,:)
        real(kind_phys), pointer :: qs(:,:)
        real(kind_phys), pointer :: qg(:,:)
        real(kind_phys), pointer :: ni(:,:)
        real(kind_phys), pointer :: nr(:,:)
        logical, pointer :: is_aerosol_aware
        real(kind_phys), pointer :: nc(:,:)
        real(kind_phys), pointer :: nwfa(:,:)
        real(kind_phys), pointer :: nifa(:,:)
        real(kind_phys), pointer :: nwfa2d(:)
        real(kind_phys), pointer :: nifa2d(:)
        real(kind_phys), pointer :: tgrs(:,:)
        real(kind_phys), pointer :: prsl(:,:)
        real(kind_phys), pointer :: phii(:,:)
        real(kind_phys), pointer :: omega(:,:)
        real(kind_phys), pointer :: dtp
        real(kind_phys), pointer :: prcp(:)
        real(kind_phys), pointer :: rain(:)
        real(kind_phys), pointer :: graupel(:)
        real(kind_phys), pointer :: ice(:)
        real(kind_phys), pointer :: snow(:)
        real(kind_phys), pointer :: sr(:)
        real(kind_phys), pointer :: refl_10cm(:,:)
        logical, pointer :: do_radar_ref
        real(kind_phys), pointer :: re_cloud(:,:)
        real(kind_phys), pointer :: re_ice(:,:)
        real(kind_phys), pointer :: re_snow(:,:)
        integer, pointer :: mpicomm
        integer, pointer :: mpirank
        integer, pointer :: mpiroot

        ierr = 0

        call c_f_pointer(ptr, cdata)


        call ccpp_field_get(cdata, 'horizontal_loop_extent', ncol, ierr=ierr, kind=ckind, index=366)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve horizontal_loop_extent from CCPP data structure')
            return
        end if
        if (kind(ncol).ne.ckind) then
            call ccpp_error('Kind mismatch for variable horizontal_loop_extent')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'vertical_dimension', nlev, ierr=ierr, kind=ckind, index=817)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve vertical_dimension from CCPP data structure')
            return
        end if
        if (kind(nlev).ne.ckind) then
            call ccpp_error('Kind mismatch for variable vertical_dimension')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'gravitational_acceleration', con_g, ierr=ierr, kind=ckind, index=355)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve gravitational_acceleration from CCPP data structure')
            return
        end if
        if (kind(con_g).ne.ckind) then
            call ccpp_error('Kind mismatch for variable gravitational_acceleration')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'gas_constant_dry_air', con_rd, ierr=ierr, kind=ckind, index=346)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve gas_constant_dry_air from CCPP data structure')
            return
        end if
        if (kind(con_rd).ne.ckind) then
            call ccpp_error('Kind mismatch for variable gas_constant_dry_air')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'water_vapor_specific_humidity_updated_by_physics', spechum, ierr=ierr, dims=cdims, kind=ckind, index=860)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve water_vapor_specific_humidity_updated_by_physics from CCPP data structure')
            return
        end if
        if (kind(spechum).ne.ckind) then
            call ccpp_error('Kind mismatch for variable water_vapor_specific_humidity_updated_by_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cloud_condensed_water_mixing_ratio_updated_by_physics', qc, ierr=ierr, dims=cdims, kind=ckind, index=95)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cloud_condensed_water_mixing_ratio_updated_by_physics from CCPP data structure')
            return
        end if
        if (kind(qc).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cloud_condensed_water_mixing_ratio_updated_by_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'rain_water_mixing_ratio_updated_by_physics', qr, ierr=ierr, dims=cdims, kind=ckind, index=619)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve rain_water_mixing_ratio_updated_by_physics from CCPP data structure')
            return
        end if
        if (kind(qr).ne.ckind) then
            call ccpp_error('Kind mismatch for variable rain_water_mixing_ratio_updated_by_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'ice_water_mixing_ratio_updated_by_physics', qi, ierr=ierr, dims=cdims, kind=ckind, index=376)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve ice_water_mixing_ratio_updated_by_physics from CCPP data structure')
            return
        end if
        if (kind(qi).ne.ckind) then
            call ccpp_error('Kind mismatch for variable ice_water_mixing_ratio_updated_by_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'snow_water_mixing_ratio_updated_by_physics', qs, ierr=ierr, dims=cdims, kind=ckind, index=653)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve snow_water_mixing_ratio_updated_by_physics from CCPP data structure')
            return
        end if
        if (kind(qs).ne.ckind) then
            call ccpp_error('Kind mismatch for variable snow_water_mixing_ratio_updated_by_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'graupel_mixing_ratio_updated_by_physics', qg, ierr=ierr, dims=cdims, kind=ckind, index=352)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve graupel_mixing_ratio_updated_by_physics from CCPP data structure')
            return
        end if
        if (kind(qg).ne.ckind) then
            call ccpp_error('Kind mismatch for variable graupel_mixing_ratio_updated_by_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'ice_number_concentration_updated_by_physics', ni, ierr=ierr, dims=cdims, kind=ckind, index=371)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve ice_number_concentration_updated_by_physics from CCPP data structure')
            return
        end if
        if (kind(ni).ne.ckind) then
            call ccpp_error('Kind mismatch for variable ice_number_concentration_updated_by_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'rain_number_concentration_updated_by_physics', nr, ierr=ierr, dims=cdims, kind=ckind, index=618)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve rain_number_concentration_updated_by_physics from CCPP data structure')
            return
        end if
        if (kind(nr).ne.ckind) then
            call ccpp_error('Kind mismatch for variable rain_number_concentration_updated_by_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'flag_for_aerosol_physics', is_aerosol_aware, ierr=ierr, kind=ckind, index=271)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_aerosol_physics from CCPP data structure')
            return
        end if
        if (kind(is_aerosol_aware).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_aerosol_physics')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'cloud_droplet_number_concentration_updated_by_physics', nc, ierr=ierr, dims=cdims, kind=ckind, index=98)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cloud_droplet_number_concentration_updated_by_physics from CCPP data structure')
            return
        end if
        if (kind(nc).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cloud_droplet_number_concentration_updated_by_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'water_friendly_aerosol_number_concentration_updated_by_physics', nwfa, ierr=ierr, dims=cdims, kind=ckind, index=850)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve water_friendly_aerosol_number_concentration_updated_by_physics from CCPP data structure')
            return
        end if
        if (kind(nwfa).ne.ckind) then
            call ccpp_error('Kind mismatch for variable water_friendly_aerosol_number_concentration_updated_by_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'ice_friendly_aerosol_number_concentration_updated_by_physics', nifa, ierr=ierr, dims=cdims, kind=ckind, index=369)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve ice_friendly_aerosol_number_concentration_updated_by_physics from CCPP data structure')
            return
        end if
        if (kind(nifa).ne.ckind) then
            call ccpp_error('Kind mismatch for variable ice_friendly_aerosol_number_concentration_updated_by_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'tendency_of_water_friendly_aerosols_at_surface', nwfa2d, ierr=ierr, dims=cdims, kind=ckind, index=779)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve tendency_of_water_friendly_aerosols_at_surface from CCPP data structure')
            return
        end if
        if (kind(nwfa2d).ne.ckind) then
            call ccpp_error('Kind mismatch for variable tendency_of_water_friendly_aerosols_at_surface')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'tendency_of_ice_friendly_aerosols_at_surface', nifa2d, ierr=ierr, dims=cdims, kind=ckind, index=769)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve tendency_of_ice_friendly_aerosols_at_surface from CCPP data structure')
            return
        end if
        if (kind(nifa2d).ne.ckind) then
            call ccpp_error('Kind mismatch for variable tendency_of_ice_friendly_aerosols_at_surface')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'air_temperature_updated_by_physics', tgrs, ierr=ierr, dims=cdims, kind=ckind, index=59)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve air_temperature_updated_by_physics from CCPP data structure')
            return
        end if
        if (kind(tgrs).ne.ckind) then
            call ccpp_error('Kind mismatch for variable air_temperature_updated_by_physics')
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
        

        call ccpp_field_get(cdata, 'omega', omega, ierr=ierr, dims=cdims, kind=ckind, index=593)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve omega from CCPP data structure')
            return
        end if
        if (kind(omega).ne.ckind) then
            call ccpp_error('Kind mismatch for variable omega')
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
        

        call ccpp_field_get(cdata, 'lwe_thickness_of_explicit_precipitation_amount', prcp, ierr=ierr, dims=cdims, kind=ckind, index=480)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve lwe_thickness_of_explicit_precipitation_amount from CCPP data structure')
            return
        end if
        if (kind(prcp).ne.ckind) then
            call ccpp_error('Kind mismatch for variable lwe_thickness_of_explicit_precipitation_amount')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'lwe_thickness_of_explicit_rain_amount', rain, ierr=ierr, dims=cdims, kind=ckind, index=481)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve lwe_thickness_of_explicit_rain_amount from CCPP data structure')
            return
        end if
        if (kind(rain).ne.ckind) then
            call ccpp_error('Kind mismatch for variable lwe_thickness_of_explicit_rain_amount')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'lwe_thickness_of_graupel_amount', graupel, ierr=ierr, dims=cdims, kind=ckind, index=483)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve lwe_thickness_of_graupel_amount from CCPP data structure')
            return
        end if
        if (kind(graupel).ne.ckind) then
            call ccpp_error('Kind mismatch for variable lwe_thickness_of_graupel_amount')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'lwe_thickness_of_ice_amount', ice, ierr=ierr, dims=cdims, kind=ckind, index=486)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve lwe_thickness_of_ice_amount from CCPP data structure')
            return
        end if
        if (kind(ice).ne.ckind) then
            call ccpp_error('Kind mismatch for variable lwe_thickness_of_ice_amount')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'lwe_thickness_of_snow_amount', snow, ierr=ierr, dims=cdims, kind=ckind, index=492)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve lwe_thickness_of_snow_amount from CCPP data structure')
            return
        end if
        if (kind(snow).ne.ckind) then
            call ccpp_error('Kind mismatch for variable lwe_thickness_of_snow_amount')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'ratio_of_snowfall_to_rainfall', sr, ierr=ierr, dims=cdims, kind=ckind, index=624)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve ratio_of_snowfall_to_rainfall from CCPP data structure')
            return
        end if
        if (kind(sr).ne.ckind) then
            call ccpp_error('Kind mismatch for variable ratio_of_snowfall_to_rainfall')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'radar_reflectivity_10cm', refl_10cm, ierr=ierr, dims=cdims, kind=ckind, index=613)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve radar_reflectivity_10cm from CCPP data structure')
            return
        end if
        if (kind(refl_10cm).ne.ckind) then
            call ccpp_error('Kind mismatch for variable radar_reflectivity_10cm')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'flag_for_radar_reflectivity', do_radar_ref, ierr=ierr, kind=ckind, index=306)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_radar_reflectivity from CCPP data structure')
            return
        end if
        if (kind(do_radar_ref).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_radar_reflectivity')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'effective_radius_of_stratiform_cloud_liquid_water_particle_in_um', re_cloud, ierr=ierr, dims=cdims, kind=ckind, index=244)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve effective_radius_of_stratiform_cloud_liquid_water_particle_in_um from CCPP data structure')
            return
        end if
        if (kind(re_cloud).ne.ckind) then
            call ccpp_error('Kind mismatch for variable effective_radius_of_stratiform_cloud_liquid_water_particle_in_um')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'effective_radius_of_stratiform_cloud_ice_particle_in_um', re_ice, ierr=ierr, dims=cdims, kind=ckind, index=243)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve effective_radius_of_stratiform_cloud_ice_particle_in_um from CCPP data structure')
            return
        end if
        if (kind(re_ice).ne.ckind) then
            call ccpp_error('Kind mismatch for variable effective_radius_of_stratiform_cloud_ice_particle_in_um')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'effective_radius_of_stratiform_cloud_snow_particle_in_um', re_snow, ierr=ierr, dims=cdims, kind=ckind, index=246)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve effective_radius_of_stratiform_cloud_snow_particle_in_um from CCPP data structure')
            return
        end if
        if (kind(re_snow).ne.ckind) then
            call ccpp_error('Kind mismatch for variable effective_radius_of_stratiform_cloud_snow_particle_in_um')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'mpi_comm', mpicomm, ierr=ierr, kind=ckind, index=557)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve mpi_comm from CCPP data structure')
            return
        end if
        if (kind(mpicomm).ne.ckind) then
            call ccpp_error('Kind mismatch for variable mpi_comm')
            ierr = 1
            return
        end if
#endif
        

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
        

        call mp_thompson_hrrr_run(ncol=ncol,nlev=nlev,con_g=con_g,con_rd=con_rd,spechum=spechum,qc=qc,qr=qr, &
                  qi=qi,qs=qs,qg=qg,ni=ni,nr=nr,is_aerosol_aware=is_aerosol_aware,nc=nc,nwfa=nwfa, &
                  nifa=nifa,nwfa2d=nwfa2d,nifa2d=nifa2d,tgrs=tgrs,prsl=prsl,phii=phii,omega=omega, &
                  dtp=dtp,prcp=prcp,rain=rain,graupel=graupel,ice=ice,snow=snow,sr=sr,refl_10cm=refl_10cm, &
                  do_radar_ref=do_radar_ref,re_cloud=re_cloud,re_ice=re_ice,re_snow=re_snow, &
                  mpicomm=mpicomm,mpirank=mpirank,mpiroot=mpiroot,errmsg=cdata%errmsg,errflg=cdata%errflg)
        ierr=cdata%errflg

    end function mp_thompson_hrrr_run_cap
end module mp_thompson_hrrr_cap


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
!! @brief Auto-generated cap module for the mp_thompson_hrrr_pre scheme
!!
!
module mp_thompson_hrrr_pre_cap

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_f_pointer, c_ptr, c_int32_t
    use            :: ccpp_types,                                      &
                      only: ccpp_t, CCPP_GENERIC_KIND
    use            :: ccpp_fields,                                     &
                      only: ccpp_field_get
    use            :: ccpp_errors,                                     &
                      only: ccpp_error, ccpp_debug
    use            :: mp_thompson_hrrr_pre, &
                      only: mp_thompson_hrrr_pre_init,mp_thompson_hrrr_pre_run,mp_thompson_hrrr_pre_finalize
    ! Other modules required, e.g. type definitions
    use machine, only: kind_phys

    implicit none

    private
    public :: mp_thompson_hrrr_pre_init_cap,mp_thompson_hrrr_pre_run_cap,mp_thompson_hrrr_pre_finalize_cap

    contains


    function mp_thompson_hrrr_pre_init_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call mp_thompson_hrrr_pre_init()
        

    end function mp_thompson_hrrr_pre_init_cap

    function mp_thompson_hrrr_pre_run_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind
        integer, pointer :: ncol
        integer, pointer :: nlev
        integer, pointer :: kdt
        real(kind_phys), pointer :: con_g
        real(kind_phys), pointer :: con_rd
        logical, pointer :: is_aerosol_aware
        real(kind_phys), pointer :: nwfa(:,:)
        real(kind_phys), pointer :: nifa(:,:)
        real(kind_phys), pointer :: nwfa2d(:)
        real(kind_phys), pointer :: nifa2d(:)
        real(kind_phys), pointer :: tgrs(:,:)
        real(kind_phys), pointer :: tgrs_save(:,:)
        real(kind_phys), pointer :: prsl(:,:)
        real(kind_phys), pointer :: phil(:,:)
        real(kind_phys), pointer :: area(:)
        integer, pointer :: mpicomm
        integer, pointer :: mpirank
        integer, pointer :: mpiroot
        integer, pointer :: blkno

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
        

        call ccpp_field_get(cdata, 'air_temperature_save', tgrs_save, ierr=ierr, dims=cdims, kind=ckind, index=57)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve air_temperature_save from CCPP data structure')
            return
        end if
        if (kind(tgrs_save).ne.ckind) then
            call ccpp_error('Kind mismatch for variable air_temperature_save')
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
        

        call ccpp_field_get(cdata, 'cell_area', area, ierr=ierr, dims=cdims, kind=ckind, index=83)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cell_area from CCPP data structure')
            return
        end if
        if (kind(area).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cell_area')
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
        

        call ccpp_field_get(cdata, 'ccpp_block_number', blkno, ierr=ierr, kind=ckind, index=82)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve ccpp_block_number from CCPP data structure')
            return
        end if
        if (kind(blkno).ne.ckind) then
            call ccpp_error('Kind mismatch for variable ccpp_block_number')
            ierr = 1
            return
        end if
#endif
        

        call mp_thompson_hrrr_pre_run(ncol=ncol,nlev=nlev,kdt=kdt,con_g=con_g,con_rd=con_rd,is_aerosol_aware=is_aerosol_aware, &
                  nwfa=nwfa,nifa=nifa,nwfa2d=nwfa2d,nifa2d=nifa2d,tgrs=tgrs,tgrs_save=tgrs_save, &
                  prsl=prsl,phil=phil,area=area,mpicomm=mpicomm,mpirank=mpirank,mpiroot=mpiroot, &
                  blkno=blkno,errmsg=cdata%errmsg,errflg=cdata%errflg)
        ierr=cdata%errflg

    end function mp_thompson_hrrr_pre_run_cap

    function mp_thompson_hrrr_pre_finalize_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call mp_thompson_hrrr_pre_finalize()
        

    end function mp_thompson_hrrr_pre_finalize_cap
end module mp_thompson_hrrr_pre_cap

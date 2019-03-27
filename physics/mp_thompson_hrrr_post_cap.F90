
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
!! @brief Auto-generated cap module for the mp_thompson_hrrr_post scheme
!!
!
module mp_thompson_hrrr_post_cap

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_f_pointer, c_ptr, c_int32_t
    use            :: ccpp_types,                                      &
                      only: ccpp_t, CCPP_GENERIC_KIND
    use            :: ccpp_fields,                                     &
                      only: ccpp_field_get
    use            :: ccpp_errors,                                     &
                      only: ccpp_error, ccpp_debug
    use            :: mp_thompson_hrrr_post, &
                      only: mp_thompson_hrrr_post_init,mp_thompson_hrrr_post_run,mp_thompson_hrrr_post_finalize
    ! Other modules required, e.g. type definitions
    use machine, only: kind_phys

    implicit none

    private
    public :: mp_thompson_hrrr_post_init_cap,mp_thompson_hrrr_post_run_cap,mp_thompson_hrrr_post_finalize_cap

    contains


    function mp_thompson_hrrr_post_init_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind
        integer, pointer :: ncol
        real(kind_phys), pointer :: area(:)
        real(kind_phys), pointer :: ttendlim

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
        

        call ccpp_field_get(cdata, 'limit_for_temperature_tendency_for_microphysics', ttendlim, ierr=ierr, kind=ckind, index=466)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve limit_for_temperature_tendency_for_microphysics from CCPP data structure')
            return
        end if
        if (kind(ttendlim).ne.ckind) then
            call ccpp_error('Kind mismatch for variable limit_for_temperature_tendency_for_microphysics')
            ierr = 1
            return
        end if
#endif
        

        call mp_thompson_hrrr_post_init(ncol=ncol,area=area,ttendlim=ttendlim,errmsg=cdata%errmsg,errflg=cdata%errflg)
        ierr=cdata%errflg

    end function mp_thompson_hrrr_post_init_cap

    function mp_thompson_hrrr_post_run_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind
        integer, pointer :: ncol
        integer, pointer :: nlev
        real(kind_phys), pointer :: tgrs_save(:,:)
        real(kind_phys), pointer :: tgrs(:,:)
        real(kind_phys), pointer :: prslk(:,:)
        real(kind_phys), pointer :: dtp
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
        

        call mp_thompson_hrrr_post_run(ncol=ncol,nlev=nlev,tgrs_save=tgrs_save,tgrs=tgrs,prslk=prslk,dtp=dtp,mpicomm=mpicomm, &
                  mpirank=mpirank,mpiroot=mpiroot,errmsg=cdata%errmsg,errflg=cdata%errflg)
        ierr=cdata%errflg

    end function mp_thompson_hrrr_post_run_cap

    function mp_thompson_hrrr_post_finalize_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call mp_thompson_hrrr_post_finalize(errmsg=cdata%errmsg,errflg=cdata%errflg)
        ierr=cdata%errflg

    end function mp_thompson_hrrr_post_finalize_cap
end module mp_thompson_hrrr_post_cap

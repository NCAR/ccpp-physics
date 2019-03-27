
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
!! @brief Auto-generated cap module for the ozphys_2015 scheme
!!
!
module ozphys_2015_cap

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_f_pointer, c_ptr, c_int32_t
    use            :: ccpp_types,                                      &
                      only: ccpp_t, CCPP_GENERIC_KIND
    use            :: ccpp_fields,                                     &
                      only: ccpp_field_get
    use            :: ccpp_errors,                                     &
                      only: ccpp_error, ccpp_debug
    use            :: ozphys_2015, &
                      only: ozphys_2015_init,ozphys_2015_run,ozphys_2015_finalize
    ! Other modules required, e.g. type definitions
    use machine, only: kind_phys

    implicit none

    private
    public :: ozphys_2015_init_cap,ozphys_2015_run_cap,ozphys_2015_finalize_cap

    contains


    function ozphys_2015_init_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call ozphys_2015_init()
        

    end function ozphys_2015_init_cap

    function ozphys_2015_run_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind
        integer, pointer :: ix
        integer, pointer :: im
        integer, pointer :: levs
        integer, pointer :: ko3
        real(kind_phys), pointer :: dt
        real(kind_phys), pointer :: oz(:,:)
        real(kind_phys), pointer :: tin(:,:)
        real(kind_phys), pointer :: po3(:)
        real(kind_phys), pointer :: prsl(:,:)
        real(kind_phys), pointer :: prdout(:,:,:)
        integer, pointer :: pl_coeff
        real(kind_phys), pointer :: delp(:,:)
        logical, pointer :: ldiag3d
        real(kind_phys), pointer :: ozp1(:,:)
        real(kind_phys), pointer :: ozp2(:,:)
        real(kind_phys), pointer :: ozp3(:,:)
        real(kind_phys), pointer :: ozp4(:,:)
        real(kind_phys), pointer :: con_g
        integer, pointer :: me

        ierr = 0

        call c_f_pointer(ptr, cdata)


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
        

        call ccpp_field_get(cdata, 'vertical_dimension_of_ozone_forcing_data', ko3, ierr=ierr, kind=ckind, index=819)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve vertical_dimension_of_ozone_forcing_data from CCPP data structure')
            return
        end if
        if (kind(ko3).ne.ckind) then
            call ccpp_error('Kind mismatch for variable vertical_dimension_of_ozone_forcing_data')
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
        

        call ccpp_field_get(cdata, 'ozone_concentration_updated_by_physics', oz, ierr=ierr, dims=cdims, kind=ckind, index=598)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve ozone_concentration_updated_by_physics from CCPP data structure')
            return
        end if
        if (kind(oz).ne.ckind) then
            call ccpp_error('Kind mismatch for variable ozone_concentration_updated_by_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'air_temperature_updated_by_physics', tin, ierr=ierr, dims=cdims, kind=ckind, index=59)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve air_temperature_updated_by_physics from CCPP data structure')
            return
        end if
        if (kind(tin).ne.ckind) then
            call ccpp_error('Kind mismatch for variable air_temperature_updated_by_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'natural_log_of_ozone_forcing_data_pressure_levels', po3, ierr=ierr, dims=cdims, kind=ckind, index=566)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve natural_log_of_ozone_forcing_data_pressure_levels from CCPP data structure')
            return
        end if
        if (kind(po3).ne.ckind) then
            call ccpp_error('Kind mismatch for variable natural_log_of_ozone_forcing_data_pressure_levels')
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
        

        call ccpp_field_get(cdata, 'ozone_forcing', prdout, ierr=ierr, dims=cdims, kind=ckind, index=599)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve ozone_forcing from CCPP data structure')
            return
        end if
        if (kind(prdout).ne.ckind) then
            call ccpp_error('Kind mismatch for variable ozone_forcing')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'number_of_coefficients_in_ozone_forcing_data', pl_coeff, ierr=ierr, kind=ckind, index=574)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve number_of_coefficients_in_ozone_forcing_data from CCPP data structure')
            return
        end if
        if (kind(pl_coeff).ne.ckind) then
            call ccpp_error('Kind mismatch for variable number_of_coefficients_in_ozone_forcing_data')
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
        

        call ccpp_field_get(cdata, 'flag_diagnostics_3D', ldiag3d, ierr=ierr, kind=ckind, index=264)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_diagnostics_3D from CCPP data structure')
            return
        end if
        if (kind(ldiag3d).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_diagnostics_3D')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'cumulative_change_in_ozone_concentration_due_to_production_and_loss_rate', ozp1, ierr=ierr, dims=cdims, kind=ckind, index=152)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cumulative_change_in_ozone_concentration_due_to_production_and_loss_rate from CCPP data structure')
            return
        end if
        if (kind(ozp1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cumulative_change_in_ozone_concentration_due_to_production_and_loss_rate')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cumulative_change_in_ozone_concentration_due_to_ozone_mixing_ratio', ozp2, ierr=ierr, dims=cdims, kind=ckind, index=151)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cumulative_change_in_ozone_concentration_due_to_ozone_mixing_ratio from CCPP data structure')
            return
        end if
        if (kind(ozp2).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cumulative_change_in_ozone_concentration_due_to_ozone_mixing_ratio')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cumulative_change_in_ozone_concentration_due_to_temperature', ozp3, ierr=ierr, dims=cdims, kind=ckind, index=153)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cumulative_change_in_ozone_concentration_due_to_temperature from CCPP data structure')
            return
        end if
        if (kind(ozp3).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cumulative_change_in_ozone_concentration_due_to_temperature')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cumulative_change_in_ozone_concentration_due_to_overhead_ozone_column', ozp4, ierr=ierr, dims=cdims, kind=ckind, index=150)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cumulative_change_in_ozone_concentration_due_to_overhead_ozone_column from CCPP data structure')
            return
        end if
        if (kind(ozp4).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cumulative_change_in_ozone_concentration_due_to_overhead_ozone_column')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

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
        

        call ozphys_2015_run(ix=ix,im=im,levs=levs,ko3=ko3,dt=dt,oz=oz,tin=tin,po3=po3,prsl=prsl,prdout=prdout, &
                  pl_coeff=pl_coeff,delp=delp,ldiag3d=ldiag3d,ozp1=ozp1,ozp2=ozp2,ozp3=ozp3, &
                  ozp4=ozp4,con_g=con_g,me=me,errmsg=cdata%errmsg,errflg=cdata%errflg)
        ierr=cdata%errflg

    end function ozphys_2015_run_cap

    function ozphys_2015_finalize_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call ozphys_2015_finalize()
        

    end function ozphys_2015_finalize_cap
end module ozphys_2015_cap

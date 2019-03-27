
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
!! @brief Auto-generated cap module for the rayleigh_damp scheme
!!
!
module rayleigh_damp_cap

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_f_pointer, c_ptr, c_int32_t
    use            :: ccpp_types,                                      &
                      only: ccpp_t, CCPP_GENERIC_KIND
    use            :: ccpp_fields,                                     &
                      only: ccpp_field_get
    use            :: ccpp_errors,                                     &
                      only: ccpp_error, ccpp_debug
    use            :: rayleigh_damp, &
                      only: rayleigh_damp_init,rayleigh_damp_finalize,rayleigh_damp_run
    ! Other modules required, e.g. type definitions
    use machine, only: kind_phys

    implicit none

    private
    public :: rayleigh_damp_init_cap,rayleigh_damp_finalize_cap,rayleigh_damp_run_cap

    contains


    function rayleigh_damp_init_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call rayleigh_damp_init()
        

    end function rayleigh_damp_init_cap

    function rayleigh_damp_finalize_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call rayleigh_damp_finalize()
        

    end function rayleigh_damp_finalize_cap

    function rayleigh_damp_run_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind
        logical, pointer :: lsidea
        integer, pointer :: im
        integer, pointer :: ix
        integer, pointer :: km
        real(kind_phys), pointer :: A(:,:)
        real(kind_phys), pointer :: B(:,:)
        real(kind_phys), pointer :: C(:,:)
        real(kind_phys), pointer :: u1(:,:)
        real(kind_phys), pointer :: v1(:,:)
        real(kind_phys), pointer :: dt
        real(kind_phys), pointer :: cp
        integer, pointer :: levr
        real(kind_phys), pointer :: pgr(:)
        real(kind_phys), pointer :: prsl(:,:)
        real(kind_phys), pointer :: prslrd0
        real(kind_phys), pointer :: ral_ts

        ierr = 0

        call c_f_pointer(ptr, cdata)


        call ccpp_field_get(cdata, 'flag_idealized_physics', lsidea, ierr=ierr, kind=ckind, index=330)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_idealized_physics from CCPP data structure')
            return
        end if
        if (kind(lsidea).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_idealized_physics')
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
        

        call ccpp_field_get(cdata, 'number_of_vertical_layers_for_radiation_calculations', levr, ierr=ierr, kind=ckind, index=591)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve number_of_vertical_layers_for_radiation_calculations from CCPP data structure')
            return
        end if
        if (kind(levr).ne.ckind) then
            call ccpp_error('Kind mismatch for variable number_of_vertical_layers_for_radiation_calculations')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'surface_air_pressure', pgr, ierr=ierr, dims=cdims, kind=ckind, index=677)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_air_pressure from CCPP data structure')
            return
        end if
        if (kind(pgr).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_air_pressure')
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
        

        call ccpp_field_get(cdata, 'pressure_cutoff_for_rayleigh_damping', prslrd0, ierr=ierr, kind=ckind, index=611)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve pressure_cutoff_for_rayleigh_damping from CCPP data structure')
            return
        end if
        if (kind(prslrd0).ne.ckind) then
            call ccpp_error('Kind mismatch for variable pressure_cutoff_for_rayleigh_damping')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'time_scale_for_rayleigh_damping', ral_ts, ierr=ierr, kind=ckind, index=791)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve time_scale_for_rayleigh_damping from CCPP data structure')
            return
        end if
        if (kind(ral_ts).ne.ckind) then
            call ccpp_error('Kind mismatch for variable time_scale_for_rayleigh_damping')
            ierr = 1
            return
        end if
#endif
        

        call rayleigh_damp_run(lsidea=lsidea,im=im,ix=ix,km=km,A=A,B=B,C=C,u1=u1,v1=v1,dt=dt,cp=cp,levr=levr, &
                  pgr=pgr,prsl=prsl,prslrd0=prslrd0,ral_ts=ral_ts,errmsg=cdata%errmsg,errflg=cdata%errflg)
        ierr=cdata%errflg

    end function rayleigh_damp_run_cap
end module rayleigh_damp_cap

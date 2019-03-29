
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
!! @brief Auto-generated cap module for the gwdps_post scheme
!!
!
module gwdps_post_cap

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_f_pointer, c_ptr, c_int32_t
    use            :: ccpp_types,                                      &
                      only: ccpp_t, CCPP_GENERIC_KIND
    use            :: ccpp_fields,                                     &
                      only: ccpp_field_get
    use            :: ccpp_errors,                                     &
                      only: ccpp_error, ccpp_debug
    use            :: gwdps_post, &
                      only: gwdps_post_init,gwdps_post_finalize,gwdps_post_run
    ! Other modules required, e.g. type definitions
    use machine, only: kind_phys

    implicit none

    private
    public :: gwdps_post_init_cap,gwdps_post_finalize_cap,gwdps_post_run_cap

    contains


    function gwdps_post_init_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call gwdps_post_init()
        

    end function gwdps_post_init_cap

    function gwdps_post_finalize_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call gwdps_post_finalize()
        

    end function gwdps_post_finalize_cap

    function gwdps_post_run_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind
        logical, pointer :: lssav
        logical, pointer :: ldiag3d
        real(kind_phys), pointer :: dtf
        real(kind_phys), pointer :: dusfcg(:)
        real(kind_phys), pointer :: dvsfcg(:)
        real(kind_phys), pointer :: dudt(:,:)
        real(kind_phys), pointer :: dvdt(:,:)
        real(kind_phys), pointer :: dtdt(:,:)
        real(kind_phys), pointer :: dugwd(:)
        real(kind_phys), pointer :: dvgwd(:)
        real(kind_phys), pointer :: du3dt(:,:)
        real(kind_phys), pointer :: dv3dt(:,:)
        real(kind_phys), pointer :: dt3dt(:,:)

        ierr = 0

        call c_f_pointer(ptr, cdata)


        call ccpp_field_get(cdata, 'flag_diagnostics', lssav, ierr=ierr, kind=ckind, index=263)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_diagnostics from CCPP data structure')
            return
        end if
        if (kind(lssav).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_diagnostics')
            ierr = 1
            return
        end if
#endif
        

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
        

        call ccpp_field_get(cdata, 'time_step_for_dynamics', dtf, ierr=ierr, kind=ckind, index=792)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve time_step_for_dynamics from CCPP data structure')
            return
        end if
        if (kind(dtf).ne.ckind) then
            call ccpp_error('Kind mismatch for variable time_step_for_dynamics')
            ierr = 1
            return
        end if
#endif
        

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
        

        call ccpp_field_get(cdata, 'tendency_of_x_wind_due_to_model_physics', dudt, ierr=ierr, dims=cdims, kind=ckind, index=782)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve tendency_of_x_wind_due_to_model_physics from CCPP data structure')
            return
        end if
        if (kind(dudt).ne.ckind) then
            call ccpp_error('Kind mismatch for variable tendency_of_x_wind_due_to_model_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'tendency_of_y_wind_due_to_model_physics', dvdt, ierr=ierr, dims=cdims, kind=ckind, index=785)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve tendency_of_y_wind_due_to_model_physics from CCPP data structure')
            return
        end if
        if (kind(dvdt).ne.ckind) then
            call ccpp_error('Kind mismatch for variable tendency_of_y_wind_due_to_model_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'tendency_of_air_temperature_due_to_model_physics', dtdt, ierr=ierr, dims=cdims, kind=ckind, index=758)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve tendency_of_air_temperature_due_to_model_physics from CCPP data structure')
            return
        end if
        if (kind(dtdt).ne.ckind) then
            call ccpp_error('Kind mismatch for variable tendency_of_air_temperature_due_to_model_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'time_integral_of_x_stress_due_to_gravity_wave_drag', dugwd, ierr=ierr, dims=cdims, kind=ckind, index=789)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve time_integral_of_x_stress_due_to_gravity_wave_drag from CCPP data structure')
            return
        end if
        if (kind(dugwd).ne.ckind) then
            call ccpp_error('Kind mismatch for variable time_integral_of_x_stress_due_to_gravity_wave_drag')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'time_integral_of_y_stress_due_to_gravity_wave_drag', dvgwd, ierr=ierr, dims=cdims, kind=ckind, index=790)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve time_integral_of_y_stress_due_to_gravity_wave_drag from CCPP data structure')
            return
        end if
        if (kind(dvgwd).ne.ckind) then
            call ccpp_error('Kind mismatch for variable time_integral_of_y_stress_due_to_gravity_wave_drag')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cumulative_change_in_x_wind_due_to_orographic_gravity_wave_drag', du3dt, ierr=ierr, dims=cdims, kind=ckind, index=168)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cumulative_change_in_x_wind_due_to_orographic_gravity_wave_drag from CCPP data structure')
            return
        end if
        if (kind(du3dt).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cumulative_change_in_x_wind_due_to_orographic_gravity_wave_drag')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cumulative_change_in_y_wind_due_to_orographic_gravity_wave_drag', dv3dt, ierr=ierr, dims=cdims, kind=ckind, index=172)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cumulative_change_in_y_wind_due_to_orographic_gravity_wave_drag from CCPP data structure')
            return
        end if
        if (kind(dv3dt).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cumulative_change_in_y_wind_due_to_orographic_gravity_wave_drag')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cumulative_change_in_temperature_due_to_shortwave_radiation_and_orographic_gravity_wave_drag', dt3dt, ierr=ierr, dims=cdims, kind=ckind, index=160)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cumulative_change_in_temperature_due_to_shortwave_radiation_and_orographic_gravity_wave_drag from CCPP data structure')
            return
        end if
        if (kind(dt3dt).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cumulative_change_in_temperature_due_to_shortwave_radiation_and_orographic_gravity_wave_drag')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call gwdps_post_run(lssav=lssav,ldiag3d=ldiag3d,dtf=dtf,dusfcg=dusfcg,dvsfcg=dvsfcg,dudt=dudt, &
                  dvdt=dvdt,dtdt=dtdt,dugwd=dugwd,dvgwd=dvgwd,du3dt=du3dt,dv3dt=dv3dt,dt3dt=dt3dt, &
                  errmsg=cdata%errmsg,errflg=cdata%errflg)
        ierr=cdata%errflg

    end function gwdps_post_run_cap
end module gwdps_post_cap

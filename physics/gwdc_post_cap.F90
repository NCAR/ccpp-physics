
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
!! @brief Auto-generated cap module for the gwdc_post scheme
!!
!
module gwdc_post_cap

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_f_pointer, c_ptr, c_int32_t
    use            :: ccpp_types,                                      &
                      only: ccpp_t, CCPP_GENERIC_KIND
    use            :: ccpp_fields,                                     &
                      only: ccpp_field_get
    use            :: ccpp_errors,                                     &
                      only: ccpp_error, ccpp_debug
    use            :: gwdc_post, &
                      only: gwdc_post_finalize,gwdc_post_run,gwdc_post_init
    ! Other modules required, e.g. type definitions
    use machine, only: kind_phys

    implicit none

    private
    public :: gwdc_post_finalize_cap,gwdc_post_run_cap,gwdc_post_init_cap

    contains


    function gwdc_post_finalize_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call gwdc_post_finalize()
        

    end function gwdc_post_finalize_cap

    function gwdc_post_run_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind
        integer, pointer :: im
        integer, pointer :: levs
        logical, pointer :: lssav
        logical, pointer :: ldiag3d
        real(kind_phys), pointer :: dtf
        real(kind_phys), pointer :: dtp
        real(kind_phys), pointer :: con_cp
        real(kind_phys), pointer :: tauctx(:)
        real(kind_phys), pointer :: taucty(:)
        real(kind_phys), pointer :: gwdcu(:,:)
        real(kind_phys), pointer :: gwdcv(:,:)
        real(kind_phys), pointer :: dugwd(:)
        real(kind_phys), pointer :: dvgwd(:)
        real(kind_phys), pointer :: du3dt(:,:)
        real(kind_phys), pointer :: dv3dt(:,:)
        real(kind_phys), pointer :: gu0(:,:)
        real(kind_phys), pointer :: gv0(:,:)
        real(kind_phys), pointer :: gt0(:,:)

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
        

        call ccpp_field_get(cdata, 'specific_heat_of_dry_air_at_constant_pressure', con_cp, ierr=ierr, kind=ckind, index=664)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve specific_heat_of_dry_air_at_constant_pressure from CCPP data structure')
            return
        end if
        if (kind(con_cp).ne.ckind) then
            call ccpp_error('Kind mismatch for variable specific_heat_of_dry_air_at_constant_pressure')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'instantaneous_x_stress_due_to_gravity_wave_drag', tauctx, ierr=ierr, dims=cdims, kind=ckind, index=445)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve instantaneous_x_stress_due_to_gravity_wave_drag from CCPP data structure')
            return
        end if
        if (kind(tauctx).ne.ckind) then
            call ccpp_error('Kind mismatch for variable instantaneous_x_stress_due_to_gravity_wave_drag')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'instantaneous_y_stress_due_to_gravity_wave_drag', taucty, ierr=ierr, dims=cdims, kind=ckind, index=447)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve instantaneous_y_stress_due_to_gravity_wave_drag from CCPP data structure')
            return
        end if
        if (kind(taucty).ne.ckind) then
            call ccpp_error('Kind mismatch for variable instantaneous_y_stress_due_to_gravity_wave_drag')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'tendency_of_x_wind_due_to_convective_gravity_wave_drag', gwdcu, ierr=ierr, dims=cdims, kind=ckind, index=781)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve tendency_of_x_wind_due_to_convective_gravity_wave_drag from CCPP data structure')
            return
        end if
        if (kind(gwdcu).ne.ckind) then
            call ccpp_error('Kind mismatch for variable tendency_of_x_wind_due_to_convective_gravity_wave_drag')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'tendency_of_y_wind_due_to_convective_gravity_wave_drag', gwdcv, ierr=ierr, dims=cdims, kind=ckind, index=784)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve tendency_of_y_wind_due_to_convective_gravity_wave_drag from CCPP data structure')
            return
        end if
        if (kind(gwdcv).ne.ckind) then
            call ccpp_error('Kind mismatch for variable tendency_of_y_wind_due_to_convective_gravity_wave_drag')
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
        

        call ccpp_field_get(cdata, 'cumulative_change_in_x_wind_due_to_convective_gravity_wave_drag', du3dt, ierr=ierr, dims=cdims, kind=ckind, index=166)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cumulative_change_in_x_wind_due_to_convective_gravity_wave_drag from CCPP data structure')
            return
        end if
        if (kind(du3dt).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cumulative_change_in_x_wind_due_to_convective_gravity_wave_drag')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cumulative_change_in_y_wind_due_to_convective_gravity_wave_drag', dv3dt, ierr=ierr, dims=cdims, kind=ckind, index=170)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cumulative_change_in_y_wind_due_to_convective_gravity_wave_drag from CCPP data structure')
            return
        end if
        if (kind(dv3dt).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cumulative_change_in_y_wind_due_to_convective_gravity_wave_drag')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'x_wind_updated_by_physics', gu0, ierr=ierr, dims=cdims, kind=ckind, index=877)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve x_wind_updated_by_physics from CCPP data structure')
            return
        end if
        if (kind(gu0).ne.ckind) then
            call ccpp_error('Kind mismatch for variable x_wind_updated_by_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'y_wind_updated_by_physics', gv0, ierr=ierr, dims=cdims, kind=ckind, index=884)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve y_wind_updated_by_physics from CCPP data structure')
            return
        end if
        if (kind(gv0).ne.ckind) then
            call ccpp_error('Kind mismatch for variable y_wind_updated_by_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'air_temperature_updated_by_physics', gt0, ierr=ierr, dims=cdims, kind=ckind, index=59)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve air_temperature_updated_by_physics from CCPP data structure')
            return
        end if
        if (kind(gt0).ne.ckind) then
            call ccpp_error('Kind mismatch for variable air_temperature_updated_by_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call gwdc_post_run(im=im,levs=levs,lssav=lssav,ldiag3d=ldiag3d,dtf=dtf,dtp=dtp,con_cp=con_cp, &
                  tauctx=tauctx,taucty=taucty,gwdcu=gwdcu,gwdcv=gwdcv,dugwd=dugwd,dvgwd=dvgwd, &
                  du3dt=du3dt,dv3dt=dv3dt,gu0=gu0,gv0=gv0,gt0=gt0,errmsg=cdata%errmsg,errflg=cdata%errflg)
        ierr=cdata%errflg

    end function gwdc_post_run_cap

    function gwdc_post_init_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call gwdc_post_init()
        

    end function gwdc_post_init_cap
end module gwdc_post_cap

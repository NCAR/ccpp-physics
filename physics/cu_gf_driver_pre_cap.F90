
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
!! @brief Auto-generated cap module for the cu_gf_driver_pre scheme
!!
!
module cu_gf_driver_pre_cap

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_f_pointer, c_ptr, c_int32_t
    use            :: ccpp_types,                                      &
                      only: ccpp_t, CCPP_GENERIC_KIND
    use            :: ccpp_fields,                                     &
                      only: ccpp_field_get
    use            :: ccpp_errors,                                     &
                      only: ccpp_error, ccpp_debug
    use            :: cu_gf_driver_pre, &
                      only: cu_gf_driver_pre_init,cu_gf_driver_pre_finalize,cu_gf_driver_pre_run
    ! Other modules required, e.g. type definitions
    use machine, only: kind_phys

    implicit none

    private
    public :: cu_gf_driver_pre_init_cap,cu_gf_driver_pre_finalize_cap,cu_gf_driver_pre_run_cap

    contains


    function cu_gf_driver_pre_init_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call cu_gf_driver_pre_init()
        

    end function cu_gf_driver_pre_init_cap

    function cu_gf_driver_pre_finalize_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call cu_gf_driver_pre_finalize()
        

    end function cu_gf_driver_pre_finalize_cap

    function cu_gf_driver_pre_run_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind
        logical, pointer :: flag_init
        logical, pointer :: flag_restart
        integer, pointer :: kdt
        real(kind_phys), pointer :: fhour
        real(kind_phys), pointer :: dtp
        real(kind_phys), pointer :: t(:,:)
        real(kind_phys), pointer :: q(:,:)
        real(kind_phys), pointer :: prevst(:,:)
        real(kind_phys), pointer :: prevsq(:,:)
        real(kind_phys), pointer :: forcet(:,:)
        real(kind_phys), pointer :: forceq(:,:)
        integer, pointer :: cactiv(:)
        real(kind_phys), pointer :: conv_act(:)

        ierr = 0

        call c_f_pointer(ptr, cdata)


        call ccpp_field_get(cdata, 'flag_for_first_time_step', flag_init, ierr=ierr, kind=ckind, index=277)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_first_time_step from CCPP data structure')
            return
        end if
        if (kind(flag_init).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_first_time_step')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'flag_for_restart', flag_restart, ierr=ierr, kind=ckind, index=309)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_restart from CCPP data structure')
            return
        end if
        if (kind(flag_restart).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_restart')
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
        

        call ccpp_field_get(cdata, 'forecast_time', fhour, ierr=ierr, kind=ckind, index=339)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve forecast_time from CCPP data structure')
            return
        end if
        if (kind(fhour).ne.ckind) then
            call ccpp_error('Kind mismatch for variable forecast_time')
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
        

        call ccpp_field_get(cdata, 'air_temperature', t, ierr=ierr, dims=cdims, kind=ckind, index=50)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve air_temperature from CCPP data structure')
            return
        end if
        if (kind(t).ne.ckind) then
            call ccpp_error('Kind mismatch for variable air_temperature')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'water_vapor_specific_humidity', q, ierr=ierr, dims=cdims, kind=ckind, index=852)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve water_vapor_specific_humidity from CCPP data structure')
            return
        end if
        if (kind(q).ne.ckind) then
            call ccpp_error('Kind mismatch for variable water_vapor_specific_humidity')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'temperature_from_previous_timestep', prevst, ierr=ierr, dims=cdims, kind=ckind, index=752)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve temperature_from_previous_timestep from CCPP data structure')
            return
        end if
        if (kind(prevst).ne.ckind) then
            call ccpp_error('Kind mismatch for variable temperature_from_previous_timestep')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'moisture_from_previous_timestep', prevsq, ierr=ierr, dims=cdims, kind=ckind, index=553)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve moisture_from_previous_timestep from CCPP data structure')
            return
        end if
        if (kind(prevsq).ne.ckind) then
            call ccpp_error('Kind mismatch for variable moisture_from_previous_timestep')
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
        

        call ccpp_field_get(cdata, 'moisture_tendency_due_to_dynamics', forceq, ierr=ierr, dims=cdims, kind=ckind, index=554)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve moisture_tendency_due_to_dynamics from CCPP data structure')
            return
        end if
        if (kind(forceq).ne.ckind) then
            call ccpp_error('Kind mismatch for variable moisture_tendency_due_to_dynamics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

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
        

        call ccpp_field_get(cdata, 'gf_memory_counter', conv_act, ierr=ierr, dims=cdims, kind=ckind, index=351)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve gf_memory_counter from CCPP data structure')
            return
        end if
        if (kind(conv_act).ne.ckind) then
            call ccpp_error('Kind mismatch for variable gf_memory_counter')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call cu_gf_driver_pre_run(flag_init=flag_init,flag_restart=flag_restart,kdt=kdt,fhour=fhour,dtp=dtp, &
                  t=t,q=q,prevst=prevst,prevsq=prevsq,forcet=forcet,forceq=forceq,cactiv=cactiv, &
                  conv_act=conv_act,errmsg=cdata%errmsg,errflg=cdata%errflg)
        ierr=cdata%errflg

    end function cu_gf_driver_pre_run_cap
end module cu_gf_driver_pre_cap


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
!! @brief Auto-generated cap module for the gwdc scheme
!!
!
module gwdc_cap

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_f_pointer, c_ptr, c_int32_t
    use            :: ccpp_types,                                      &
                      only: ccpp_t, CCPP_GENERIC_KIND
    use            :: ccpp_fields,                                     &
                      only: ccpp_field_get
    use            :: ccpp_errors,                                     &
                      only: ccpp_error, ccpp_debug
    use            :: gwdc, &
                      only: gwdc_init,gwdc_finalize,gwdc_run
    ! Other modules required, e.g. type definitions
    use machine, only: kind_phys

    implicit none

    private
    public :: gwdc_init_cap,gwdc_finalize_cap,gwdc_run_cap

    contains


    function gwdc_init_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call gwdc_init()
        

    end function gwdc_init_cap

    function gwdc_finalize_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call gwdc_finalize()
        

    end function gwdc_finalize_cap

    function gwdc_run_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind
        integer, pointer :: im
        integer, pointer :: ix
        integer, pointer :: km
        integer, pointer :: lat
        real(kind_phys), pointer :: u1(:,:)
        real(kind_phys), pointer :: v1(:,:)
        real(kind_phys), pointer :: t1(:,:)
        real(kind_phys), pointer :: q1(:,:)
        real(kind_phys), pointer :: deltim
        real(kind_phys), pointer :: pmid1(:,:)
        real(kind_phys), pointer :: pint1(:,:)
        real(kind_phys), pointer :: dpmid1(:,:)
        real(kind_phys), pointer :: qmax(:)
        integer, pointer :: ktop(:)
        integer, pointer :: kbot(:)
        integer, pointer :: kcnv(:)
        real(kind_phys), pointer :: cldf(:)
        real(kind_phys), pointer :: grav
        real(kind_phys), pointer :: cp
        real(kind_phys), pointer :: rd
        real(kind_phys), pointer :: fv
        real(kind_phys), pointer :: pi
        real(kind_phys), pointer :: dlength(:)
        logical, pointer :: lprnt
        integer, pointer :: ipr
        real(kind_phys), pointer :: fhour
        real(kind_phys), pointer :: utgwc(:,:)
        real(kind_phys), pointer :: vtgwc(:,:)
        real(kind_phys), pointer :: tauctx(:)
        real(kind_phys), pointer :: taucty(:)

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
        

        call ccpp_field_get(cdata, 'latitude_index_in_debug_printouts', lat, ierr=ierr, kind=ckind, index=462)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve latitude_index_in_debug_printouts from CCPP data structure')
            return
        end if
        if (kind(lat).ne.ckind) then
            call ccpp_error('Kind mismatch for variable latitude_index_in_debug_printouts')
            ierr = 1
            return
        end if
#endif
        

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
        

        call ccpp_field_get(cdata, 'air_temperature', t1, ierr=ierr, dims=cdims, kind=ckind, index=50)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve air_temperature from CCPP data structure')
            return
        end if
        if (kind(t1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable air_temperature')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'water_vapor_specific_humidity', q1, ierr=ierr, dims=cdims, kind=ckind, index=852)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve water_vapor_specific_humidity from CCPP data structure')
            return
        end if
        if (kind(q1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable water_vapor_specific_humidity')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'time_step_for_physics', deltim, ierr=ierr, kind=ckind, index=793)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve time_step_for_physics from CCPP data structure')
            return
        end if
        if (kind(deltim).ne.ckind) then
            call ccpp_error('Kind mismatch for variable time_step_for_physics')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'air_pressure', pmid1, ierr=ierr, dims=cdims, kind=ckind, index=44)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve air_pressure from CCPP data structure')
            return
        end if
        if (kind(pmid1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable air_pressure')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'air_pressure_at_interface', pint1, ierr=ierr, dims=cdims, kind=ckind, index=45)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve air_pressure_at_interface from CCPP data structure')
            return
        end if
        if (kind(pint1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable air_pressure_at_interface')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'air_pressure_difference_between_midlayers', dpmid1, ierr=ierr, dims=cdims, kind=ckind, index=49)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve air_pressure_difference_between_midlayers from CCPP data structure')
            return
        end if
        if (kind(dpmid1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable air_pressure_difference_between_midlayers')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'maximum_column_heating_rate', qmax, ierr=ierr, dims=cdims, kind=ckind, index=503)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve maximum_column_heating_rate from CCPP data structure')
            return
        end if
        if (kind(qmax).ne.ckind) then
            call ccpp_error('Kind mismatch for variable maximum_column_heating_rate')
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
        

        call ccpp_field_get(cdata, 'cloud_area_fraction', cldf, ierr=ierr, dims=cdims, kind=ckind, index=86)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cloud_area_fraction from CCPP data structure')
            return
        end if
        if (kind(cldf).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cloud_area_fraction')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

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
        

        call ccpp_field_get(cdata, 'pi', pi, ierr=ierr, kind=ckind, index=606)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve pi from CCPP data structure')
            return
        end if
        if (kind(pi).ne.ckind) then
            call ccpp_error('Kind mismatch for variable pi')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'characteristic_grid_length_scale', dlength, ierr=ierr, dims=cdims, kind=ckind, index=85)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve characteristic_grid_length_scale from CCPP data structure')
            return
        end if
        if (kind(dlength).ne.ckind) then
            call ccpp_error('Kind mismatch for variable characteristic_grid_length_scale')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'flag_print', lprnt, ierr=ierr, kind=ckind, index=332)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_print from CCPP data structure')
            return
        end if
        if (kind(lprnt).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_print')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'horizontal_index_of_printed_column', ipr, ierr=ierr, kind=ckind, index=365)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve horizontal_index_of_printed_column from CCPP data structure')
            return
        end if
        if (kind(ipr).ne.ckind) then
            call ccpp_error('Kind mismatch for variable horizontal_index_of_printed_column')
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
        

        call ccpp_field_get(cdata, 'tendency_of_x_wind_due_to_convective_gravity_wave_drag', utgwc, ierr=ierr, dims=cdims, kind=ckind, index=781)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve tendency_of_x_wind_due_to_convective_gravity_wave_drag from CCPP data structure')
            return
        end if
        if (kind(utgwc).ne.ckind) then
            call ccpp_error('Kind mismatch for variable tendency_of_x_wind_due_to_convective_gravity_wave_drag')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'tendency_of_y_wind_due_to_convective_gravity_wave_drag', vtgwc, ierr=ierr, dims=cdims, kind=ckind, index=784)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve tendency_of_y_wind_due_to_convective_gravity_wave_drag from CCPP data structure')
            return
        end if
        if (kind(vtgwc).ne.ckind) then
            call ccpp_error('Kind mismatch for variable tendency_of_y_wind_due_to_convective_gravity_wave_drag')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

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
        

        call gwdc_run(im=im,ix=ix,km=km,lat=lat,u1=u1,v1=v1,t1=t1,q1=q1,deltim=deltim,pmid1=pmid1, &
                  pint1=pint1,dpmid1=dpmid1,qmax=qmax,ktop=ktop,kbot=kbot,kcnv=kcnv,cldf=cldf, &
                  grav=grav,cp=cp,rd=rd,fv=fv,pi=pi,dlength=dlength,lprnt=lprnt,ipr=ipr,fhour=fhour, &
                  utgwc=utgwc,vtgwc=vtgwc,tauctx=tauctx,taucty=taucty,errmsg=cdata%errmsg, &
                  errflg=cdata%errflg)
        ierr=cdata%errflg

    end function gwdc_run_cap
end module gwdc_cap

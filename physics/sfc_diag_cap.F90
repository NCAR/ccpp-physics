
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
!! @brief Auto-generated cap module for the sfc_diag scheme
!!
!
module sfc_diag_cap

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_f_pointer, c_ptr, c_int32_t
    use            :: ccpp_types,                                      &
                      only: ccpp_t, CCPP_GENERIC_KIND
    use            :: ccpp_fields,                                     &
                      only: ccpp_field_get
    use            :: ccpp_errors,                                     &
                      only: ccpp_error, ccpp_debug
    use            :: sfc_diag, &
                      only: sfc_diag_init,sfc_diag_finalize,sfc_diag_run
    ! Other modules required, e.g. type definitions
    use machine, only: kind_phys

    implicit none

    private
    public :: sfc_diag_init_cap,sfc_diag_finalize_cap,sfc_diag_run_cap

    contains


    function sfc_diag_init_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call sfc_diag_init()
        

    end function sfc_diag_init_cap

    function sfc_diag_finalize_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call sfc_diag_finalize()
        

    end function sfc_diag_finalize_cap

    function sfc_diag_run_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind
        integer, pointer :: im
        real(kind_phys), pointer :: grav
        real(kind_phys), pointer :: cp
        real(kind_phys), pointer :: eps
        real(kind_phys), pointer :: epsm1
        real(kind_phys), pointer :: ps(:)
        real(kind_phys), pointer :: u1(:)
        real(kind_phys), pointer :: v1(:)
        real(kind_phys), pointer :: t1(:)
        real(kind_phys), pointer :: q1(:)
        real(kind_phys), pointer :: tskin(:)
        real(kind_phys), pointer :: qsurf(:)
        real(kind_phys), pointer :: f10m(:)
        real(kind_phys), pointer :: u10m(:)
        real(kind_phys), pointer :: v10m(:)
        real(kind_phys), pointer :: t2m(:)
        real(kind_phys), pointer :: q2m(:)
        real(kind_phys), pointer :: prslki(:)
        real(kind_phys), pointer :: evap(:)
        real(kind_phys), pointer :: fm(:)
        real(kind_phys), pointer :: fh(:)
        real(kind_phys), pointer :: fm10(:)
        real(kind_phys), pointer :: fh2(:)

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
        

        call ccpp_field_get(cdata, 'ratio_of_dry_air_to_water_vapor_gas_constants', eps, ierr=ierr, kind=ckind, index=621)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve ratio_of_dry_air_to_water_vapor_gas_constants from CCPP data structure')
            return
        end if
        if (kind(eps).ne.ckind) then
            call ccpp_error('Kind mismatch for variable ratio_of_dry_air_to_water_vapor_gas_constants')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'ratio_of_dry_air_to_water_vapor_gas_constants_minus_one', epsm1, ierr=ierr, kind=ckind, index=622)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve ratio_of_dry_air_to_water_vapor_gas_constants_minus_one from CCPP data structure')
            return
        end if
        if (kind(epsm1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable ratio_of_dry_air_to_water_vapor_gas_constants_minus_one')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'surface_air_pressure', ps, ierr=ierr, dims=cdims, kind=ckind, index=677)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_air_pressure from CCPP data structure')
            return
        end if
        if (kind(ps).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_air_pressure')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'x_wind_at_lowest_model_layer_updated_by_physics', u1, ierr=ierr, dims=cdims, kind=ckind, index=875)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve x_wind_at_lowest_model_layer_updated_by_physics from CCPP data structure')
            return
        end if
        if (kind(u1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable x_wind_at_lowest_model_layer_updated_by_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'y_wind_at_lowest_model_layer_updated_by_physics', v1, ierr=ierr, dims=cdims, kind=ckind, index=882)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve y_wind_at_lowest_model_layer_updated_by_physics from CCPP data structure')
            return
        end if
        if (kind(v1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable y_wind_at_lowest_model_layer_updated_by_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'air_temperature_at_lowest_model_layer_updated_by_physics', t1, ierr=ierr, dims=cdims, kind=ckind, index=55)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve air_temperature_at_lowest_model_layer_updated_by_physics from CCPP data structure')
            return
        end if
        if (kind(t1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable air_temperature_at_lowest_model_layer_updated_by_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'water_vapor_specific_humidity_at_lowest_model_layer_updated_by_physics', q1, ierr=ierr, dims=cdims, kind=ckind, index=856)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve water_vapor_specific_humidity_at_lowest_model_layer_updated_by_physics from CCPP data structure')
            return
        end if
        if (kind(q1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable water_vapor_specific_humidity_at_lowest_model_layer_updated_by_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_skin_temperature', tskin, ierr=ierr, dims=cdims, kind=ckind, index=721)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_skin_temperature from CCPP data structure')
            return
        end if
        if (kind(tskin).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_skin_temperature')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_specific_humidity', qsurf, ierr=ierr, dims=cdims, kind=ckind, index=730)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_specific_humidity from CCPP data structure')
            return
        end if
        if (kind(qsurf).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_specific_humidity')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'ratio_of_wind_at_lowest_model_layer_and_wind_at_10m', f10m, ierr=ierr, dims=cdims, kind=ckind, index=626)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve ratio_of_wind_at_lowest_model_layer_and_wind_at_10m from CCPP data structure')
            return
        end if
        if (kind(f10m).ne.ckind) then
            call ccpp_error('Kind mismatch for variable ratio_of_wind_at_lowest_model_layer_and_wind_at_10m')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'x_wind_at_10m', u10m, ierr=ierr, dims=cdims, kind=ckind, index=872)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve x_wind_at_10m from CCPP data structure')
            return
        end if
        if (kind(u10m).ne.ckind) then
            call ccpp_error('Kind mismatch for variable x_wind_at_10m')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'y_wind_at_10m', v10m, ierr=ierr, dims=cdims, kind=ckind, index=879)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve y_wind_at_10m from CCPP data structure')
            return
        end if
        if (kind(v10m).ne.ckind) then
            call ccpp_error('Kind mismatch for variable y_wind_at_10m')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'temperature_at_2m', t2m, ierr=ierr, dims=cdims, kind=ckind, index=750)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve temperature_at_2m from CCPP data structure')
            return
        end if
        if (kind(t2m).ne.ckind) then
            call ccpp_error('Kind mismatch for variable temperature_at_2m')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'specific_humidity_at_2m', q2m, ierr=ierr, dims=cdims, kind=ckind, index=667)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve specific_humidity_at_2m from CCPP data structure')
            return
        end if
        if (kind(q2m).ne.ckind) then
            call ccpp_error('Kind mismatch for variable specific_humidity_at_2m')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'ratio_of_exner_function_between_midlayer_and_interface_at_lowest_model_layer', prslki, ierr=ierr, dims=cdims, kind=ckind, index=623)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve ratio_of_exner_function_between_midlayer_and_interface_at_lowest_model_layer from CCPP data structure')
            return
        end if
        if (kind(prslki).ne.ckind) then
            call ccpp_error('Kind mismatch for variable ratio_of_exner_function_between_midlayer_and_interface_at_lowest_model_layer')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'kinematic_surface_upward_latent_heat_flux', evap, ierr=ierr, dims=cdims, kind=ckind, index=454)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve kinematic_surface_upward_latent_heat_flux from CCPP data structure')
            return
        end if
        if (kind(evap).ne.ckind) then
            call ccpp_error('Kind mismatch for variable kinematic_surface_upward_latent_heat_flux')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'Monin-Obukhov_similarity_function_for_momentum', fm, ierr=ierr, dims=cdims, kind=ckind, index=18)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve Monin-Obukhov_similarity_function_for_momentum from CCPP data structure')
            return
        end if
        if (kind(fm).ne.ckind) then
            call ccpp_error('Kind mismatch for variable Monin-Obukhov_similarity_function_for_momentum')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'Monin-Obukhov_similarity_function_for_heat', fh, ierr=ierr, dims=cdims, kind=ckind, index=16)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve Monin-Obukhov_similarity_function_for_heat from CCPP data structure')
            return
        end if
        if (kind(fh).ne.ckind) then
            call ccpp_error('Kind mismatch for variable Monin-Obukhov_similarity_function_for_heat')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'Monin-Obukhov_similarity_function_for_momentum_at_10m', fm10, ierr=ierr, dims=cdims, kind=ckind, index=19)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve Monin-Obukhov_similarity_function_for_momentum_at_10m from CCPP data structure')
            return
        end if
        if (kind(fm10).ne.ckind) then
            call ccpp_error('Kind mismatch for variable Monin-Obukhov_similarity_function_for_momentum_at_10m')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'Monin-Obukhov_similarity_function_for_heat_at_2m', fh2, ierr=ierr, dims=cdims, kind=ckind, index=17)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve Monin-Obukhov_similarity_function_for_heat_at_2m from CCPP data structure')
            return
        end if
        if (kind(fh2).ne.ckind) then
            call ccpp_error('Kind mismatch for variable Monin-Obukhov_similarity_function_for_heat_at_2m')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call sfc_diag_run(im=im,grav=grav,cp=cp,eps=eps,epsm1=epsm1,ps=ps,u1=u1,v1=v1,t1=t1,q1=q1, &
                  tskin=tskin,qsurf=qsurf,f10m=f10m,u10m=u10m,v10m=v10m,t2m=t2m,q2m=q2m,prslki=prslki, &
                  evap=evap,fm=fm,fh=fh,fm10=fm10,fh2=fh2,errmsg=cdata%errmsg,errflg=cdata%errflg)
        ierr=cdata%errflg

    end function sfc_diag_run_cap
end module sfc_diag_cap

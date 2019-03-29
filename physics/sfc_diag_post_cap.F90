
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
!! @brief Auto-generated cap module for the sfc_diag_post scheme
!!
!
module sfc_diag_post_cap

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_f_pointer, c_ptr, c_int32_t
    use            :: ccpp_types,                                      &
                      only: ccpp_t, CCPP_GENERIC_KIND
    use            :: ccpp_fields,                                     &
                      only: ccpp_field_get
    use            :: ccpp_errors,                                     &
                      only: ccpp_error, ccpp_debug
    use            :: sfc_diag_post, &
                      only: sfc_diag_post_init,sfc_diag_post_run,sfc_diag_post_finalize
    ! Other modules required, e.g. type definitions
    use machine, only: kind_phys

    implicit none

    private
    public :: sfc_diag_post_init_cap,sfc_diag_post_run_cap,sfc_diag_post_finalize_cap

    contains


    function sfc_diag_post_init_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call sfc_diag_post_init()
        

    end function sfc_diag_post_init_cap

    function sfc_diag_post_run_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind
        integer, pointer :: im
        logical, pointer :: lssav
        real(kind_phys), pointer :: dtf
        real(kind_phys), pointer :: con_eps
        real(kind_phys), pointer :: con_epsm1
        real(kind_phys), pointer :: pgr(:)
        real(kind_phys), pointer :: t2m(:)
        real(kind_phys), pointer :: q2m(:)
        real(kind_phys), pointer :: u10m(:)
        real(kind_phys), pointer :: v10m(:)
        real(kind_phys), pointer :: tmpmin(:)
        real(kind_phys), pointer :: tmpmax(:)
        real(kind_phys), pointer :: spfhmin(:)
        real(kind_phys), pointer :: spfhmax(:)
        real(kind_phys), pointer :: wind10mmax(:)
        real(kind_phys), pointer :: u10mmax(:)
        real(kind_phys), pointer :: v10mmax(:)
        real(kind_phys), pointer :: dpt2m(:)

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
        

        call ccpp_field_get(cdata, 'ratio_of_dry_air_to_water_vapor_gas_constants', con_eps, ierr=ierr, kind=ckind, index=621)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve ratio_of_dry_air_to_water_vapor_gas_constants from CCPP data structure')
            return
        end if
        if (kind(con_eps).ne.ckind) then
            call ccpp_error('Kind mismatch for variable ratio_of_dry_air_to_water_vapor_gas_constants')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'ratio_of_dry_air_to_water_vapor_gas_constants_minus_one', con_epsm1, ierr=ierr, kind=ckind, index=622)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve ratio_of_dry_air_to_water_vapor_gas_constants_minus_one from CCPP data structure')
            return
        end if
        if (kind(con_epsm1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable ratio_of_dry_air_to_water_vapor_gas_constants_minus_one')
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
        

        call ccpp_field_get(cdata, 'minimum_temperature_at_2m', tmpmin, ierr=ierr, dims=cdims, kind=ckind, index=546)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve minimum_temperature_at_2m from CCPP data structure')
            return
        end if
        if (kind(tmpmin).ne.ckind) then
            call ccpp_error('Kind mismatch for variable minimum_temperature_at_2m')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'maximum_temperature_at_2m', tmpmax, ierr=ierr, dims=cdims, kind=ckind, index=508)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve maximum_temperature_at_2m from CCPP data structure')
            return
        end if
        if (kind(tmpmax).ne.ckind) then
            call ccpp_error('Kind mismatch for variable maximum_temperature_at_2m')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'minimum_specific_humidity_at_2m', spfhmin, ierr=ierr, dims=cdims, kind=ckind, index=545)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve minimum_specific_humidity_at_2m from CCPP data structure')
            return
        end if
        if (kind(spfhmin).ne.ckind) then
            call ccpp_error('Kind mismatch for variable minimum_specific_humidity_at_2m')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'maximum_specific_humidity_at_2m', spfhmax, ierr=ierr, dims=cdims, kind=ckind, index=506)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve maximum_specific_humidity_at_2m from CCPP data structure')
            return
        end if
        if (kind(spfhmax).ne.ckind) then
            call ccpp_error('Kind mismatch for variable maximum_specific_humidity_at_2m')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'maximum_wind_at_10m', wind10mmax, ierr=ierr, dims=cdims, kind=ckind, index=511)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve maximum_wind_at_10m from CCPP data structure')
            return
        end if
        if (kind(wind10mmax).ne.ckind) then
            call ccpp_error('Kind mismatch for variable maximum_wind_at_10m')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'maximum_x_wind_at_10m', u10mmax, ierr=ierr, dims=cdims, kind=ckind, index=512)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve maximum_x_wind_at_10m from CCPP data structure')
            return
        end if
        if (kind(u10mmax).ne.ckind) then
            call ccpp_error('Kind mismatch for variable maximum_x_wind_at_10m')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'maximum_y_wind_at_10m', v10mmax, ierr=ierr, dims=cdims, kind=ckind, index=513)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve maximum_y_wind_at_10m from CCPP data structure')
            return
        end if
        if (kind(v10mmax).ne.ckind) then
            call ccpp_error('Kind mismatch for variable maximum_y_wind_at_10m')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'dewpoint_temperature_at_2m', dpt2m, ierr=ierr, dims=cdims, kind=ckind, index=218)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve dewpoint_temperature_at_2m from CCPP data structure')
            return
        end if
        if (kind(dpt2m).ne.ckind) then
            call ccpp_error('Kind mismatch for variable dewpoint_temperature_at_2m')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call sfc_diag_post_run(im=im,lssav=lssav,dtf=dtf,con_eps=con_eps,con_epsm1=con_epsm1,pgr=pgr,t2m=t2m, &
                  q2m=q2m,u10m=u10m,v10m=v10m,tmpmin=tmpmin,tmpmax=tmpmax,spfhmin=spfhmin, &
                  spfhmax=spfhmax,wind10mmax=wind10mmax,u10mmax=u10mmax,v10mmax=v10mmax,dpt2m=dpt2m, &
                  errmsg=cdata%errmsg,errflg=cdata%errflg)
        ierr=cdata%errflg

    end function sfc_diag_post_run_cap

    function sfc_diag_post_finalize_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call sfc_diag_post_finalize()
        

    end function sfc_diag_post_finalize_cap
end module sfc_diag_post_cap


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
!! @brief Auto-generated cap module for the rrtmg_sw_pre scheme
!!
!
module rrtmg_sw_pre_cap

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_f_pointer, c_ptr, c_int32_t
    use            :: ccpp_types,                                      &
                      only: ccpp_t, CCPP_GENERIC_KIND
    use            :: ccpp_fields,                                     &
                      only: ccpp_field_get
    use            :: ccpp_errors,                                     &
                      only: ccpp_error, ccpp_debug
    use            :: rrtmg_sw_pre, &
                      only: rrtmg_sw_pre_init,rrtmg_sw_pre_finalize,rrtmg_sw_pre_run
    ! Other modules required, e.g. type definitions
    use GFS_typedefs, only: GFS_control_type
    use GFS_typedefs, only: GFS_grid_type
    use GFS_typedefs, only: GFS_sfcprop_type
    use GFS_typedefs, only: GFS_radtend_type
    use machine, only: kind_phys

    implicit none

    private
    public :: rrtmg_sw_pre_init_cap,rrtmg_sw_pre_finalize_cap,rrtmg_sw_pre_run_cap

    contains


    function rrtmg_sw_pre_init_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call rrtmg_sw_pre_init()
        

    end function rrtmg_sw_pre_init_cap

    function rrtmg_sw_pre_finalize_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call rrtmg_sw_pre_finalize()
        

    end function rrtmg_sw_pre_finalize_cap

    function rrtmg_sw_pre_run_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind
        type(GFS_control_type), pointer     :: Model
        type(GFS_grid_type), pointer     :: Grid
        type(GFS_sfcprop_type), pointer     :: Sfcprop
        type(GFS_radtend_type), pointer     :: Radtend
        integer, pointer :: im
        integer, pointer :: nday
        integer, pointer :: idxday(:)
        real(kind_phys), pointer :: tsfg(:)
        real(kind_phys), pointer :: tsfa(:)
        real(kind_phys), pointer :: sfcalb1(:)
        real(kind_phys), pointer :: sfcalb2(:)
        real(kind_phys), pointer :: sfcalb3(:)
        real(kind_phys), pointer :: sfcalb4(:)
        real(kind_phys), pointer :: alb1d(:)

        ierr = 0

        call c_f_pointer(ptr, cdata)


        call ccpp_field_get(cdata, 'GFS_control_type_instance', cptr, ierr=ierr, kind=ckind, index=2)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve GFS_control_type_instance from CCPP data structure')
            return
        end if
        if (ckind.ne.CCPP_GENERIC_KIND) then
            call ccpp_error('Kind mismatch for variable GFS_control_type_instance')
            ierr = 1
            return
        end if
#endif
        call c_f_pointer(cptr, Model)

        call ccpp_field_get(cdata, 'GFS_grid_type_instance', cptr, ierr=ierr, kind=ckind, index=6)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve GFS_grid_type_instance from CCPP data structure')
            return
        end if
        if (ckind.ne.CCPP_GENERIC_KIND) then
            call ccpp_error('Kind mismatch for variable GFS_grid_type_instance')
            ierr = 1
            return
        end if
#endif
        call c_f_pointer(cptr, Grid)

        call ccpp_field_get(cdata, 'GFS_sfcprop_type_instance', cptr, ierr=ierr, kind=ckind, index=11)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve GFS_sfcprop_type_instance from CCPP data structure')
            return
        end if
        if (ckind.ne.CCPP_GENERIC_KIND) then
            call ccpp_error('Kind mismatch for variable GFS_sfcprop_type_instance')
            ierr = 1
            return
        end if
#endif
        call c_f_pointer(cptr, Sfcprop)

        call ccpp_field_get(cdata, 'GFS_radtend_type_instance', cptr, ierr=ierr, kind=ckind, index=9)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve GFS_radtend_type_instance from CCPP data structure')
            return
        end if
        if (ckind.ne.CCPP_GENERIC_KIND) then
            call ccpp_error('Kind mismatch for variable GFS_radtend_type_instance')
            ierr = 1
            return
        end if
#endif
        call c_f_pointer(cptr, Radtend)

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
        

        call ccpp_field_get(cdata, 'daytime_points_dimension', nday, ierr=ierr, kind=ckind, index=209)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve daytime_points_dimension from CCPP data structure')
            return
        end if
        if (kind(nday).ne.ckind) then
            call ccpp_error('Kind mismatch for variable daytime_points_dimension')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'daytime_points', idxday, ierr=ierr, dims=cdims, kind=ckind, index=208)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve daytime_points from CCPP data structure')
            return
        end if
        if (kind(idxday).ne.ckind) then
            call ccpp_error('Kind mismatch for variable daytime_points')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_ground_temperature_for_radiation', tsfg, ierr=ierr, dims=cdims, kind=ckind, index=712)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_ground_temperature_for_radiation from CCPP data structure')
            return
        end if
        if (kind(tsfg).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_ground_temperature_for_radiation')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_air_temperature_for_radiation', tsfa, ierr=ierr, dims=cdims, kind=ckind, index=681)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_air_temperature_for_radiation from CCPP data structure')
            return
        end if
        if (kind(tsfa).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_air_temperature_for_radiation')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_albedo_due_to_near_IR_direct', sfcalb1, ierr=ierr, dims=cdims, kind=ckind, index=685)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_albedo_due_to_near_IR_direct from CCPP data structure')
            return
        end if
        if (kind(sfcalb1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_albedo_due_to_near_IR_direct')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_albedo_due_to_near_IR_diffused', sfcalb2, ierr=ierr, dims=cdims, kind=ckind, index=684)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_albedo_due_to_near_IR_diffused from CCPP data structure')
            return
        end if
        if (kind(sfcalb2).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_albedo_due_to_near_IR_diffused')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_albedo_due_to_UV_and_VIS_direct', sfcalb3, ierr=ierr, dims=cdims, kind=ckind, index=683)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_albedo_due_to_UV_and_VIS_direct from CCPP data structure')
            return
        end if
        if (kind(sfcalb3).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_albedo_due_to_UV_and_VIS_direct')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_albedo_due_to_UV_and_VIS_diffused', sfcalb4, ierr=ierr, dims=cdims, kind=ckind, index=682)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_albedo_due_to_UV_and_VIS_diffused from CCPP data structure')
            return
        end if
        if (kind(sfcalb4).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_albedo_due_to_UV_and_VIS_diffused')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_albedo_perturbation', alb1d, ierr=ierr, dims=cdims, kind=ckind, index=686)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_albedo_perturbation from CCPP data structure')
            return
        end if
        if (kind(alb1d).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_albedo_perturbation')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call rrtmg_sw_pre_run(Model=Model,Grid=Grid,Sfcprop=Sfcprop,Radtend=Radtend,im=im,nday=nday,idxday=idxday, &
                  tsfg=tsfg,tsfa=tsfa,sfcalb1=sfcalb1,sfcalb2=sfcalb2,sfcalb3=sfcalb3,sfcalb4=sfcalb4, &
                  alb1d=alb1d,errmsg=cdata%errmsg,errflg=cdata%errflg)
        ierr=cdata%errflg

    end function rrtmg_sw_pre_run_cap
end module rrtmg_sw_pre_cap


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
!! @brief Auto-generated cap module for the rrtmg_lw_post scheme
!!
!
module rrtmg_lw_post_cap

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_f_pointer, c_ptr, c_int32_t
    use            :: ccpp_types,                                      &
                      only: ccpp_t, CCPP_GENERIC_KIND
    use            :: ccpp_fields,                                     &
                      only: ccpp_field_get
    use            :: ccpp_errors,                                     &
                      only: ccpp_error, ccpp_debug
    use            :: rrtmg_lw_post, &
                      only: rrtmg_lw_post_init,rrtmg_lw_post_run,rrtmg_lw_post_finalize
    ! Other modules required, e.g. type definitions
    use GFS_typedefs, only: GFS_control_type
    use GFS_typedefs, only: GFS_grid_type
    use GFS_typedefs, only: GFS_radtend_type
    use GFS_typedefs, only: GFS_coupling_type
    use machine, only: kind_phys

    implicit none

    private
    public :: rrtmg_lw_post_init_cap,rrtmg_lw_post_run_cap,rrtmg_lw_post_finalize_cap

    contains


    function rrtmg_lw_post_init_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call rrtmg_lw_post_init()
        

    end function rrtmg_lw_post_init_cap

    function rrtmg_lw_post_run_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind
        type(GFS_control_type), pointer     :: Model
        type(GFS_grid_type), pointer     :: Grid
        type(GFS_radtend_type), pointer     :: Radtend
        type(GFS_coupling_type), pointer     :: Coupling
        integer, pointer :: im
        integer, pointer :: ltp
        integer, pointer :: lm
        integer, pointer :: kd
        real(kind_phys), pointer :: tsfa(:)
        real(kind_phys), pointer :: htlwc(:,:)
        real(kind_phys), pointer :: htlw0(:,:)

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

        call ccpp_field_get(cdata, 'GFS_coupling_type_instance', cptr, ierr=ierr, kind=ckind, index=3)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve GFS_coupling_type_instance from CCPP data structure')
            return
        end if
        if (ckind.ne.CCPP_GENERIC_KIND) then
            call ccpp_error('Kind mismatch for variable GFS_coupling_type_instance')
            ierr = 1
            return
        end if
#endif
        call c_f_pointer(cptr, Coupling)

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
        

        call ccpp_field_get(cdata, 'extra_top_layer', ltp, ierr=ierr, kind=ckind, index=257)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve extra_top_layer from CCPP data structure')
            return
        end if
        if (kind(ltp).ne.ckind) then
            call ccpp_error('Kind mismatch for variable extra_top_layer')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'vertical_layer_dimension_for_radiation', lm, ierr=ierr, kind=ckind, index=826)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve vertical_layer_dimension_for_radiation from CCPP data structure')
            return
        end if
        if (kind(lm).ne.ckind) then
            call ccpp_error('Kind mismatch for variable vertical_layer_dimension_for_radiation')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'vertical_index_difference_between_inout_and_local', kd, ierr=ierr, kind=ckind, index=823)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve vertical_index_difference_between_inout_and_local from CCPP data structure')
            return
        end if
        if (kind(kd).ne.ckind) then
            call ccpp_error('Kind mismatch for variable vertical_index_difference_between_inout_and_local')
            ierr = 1
            return
        end if
#endif
        

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
        

        call ccpp_field_get(cdata, 'tendency_of_air_temperature_due_to_longwave_heating_on_radiation_time_step', htlwc, ierr=ierr, dims=cdims, kind=ckind, index=756)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve tendency_of_air_temperature_due_to_longwave_heating_on_radiation_time_step from CCPP data structure')
            return
        end if
        if (kind(htlwc).ne.ckind) then
            call ccpp_error('Kind mismatch for variable tendency_of_air_temperature_due_to_longwave_heating_on_radiation_time_step')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'tendency_of_air_temperature_due_to_longwave_heating_assuming_clear_sky_on_radiation_time_step', htlw0, ierr=ierr, dims=cdims, kind=ckind, index=754)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve tendency_of_air_temperature_due_to_longwave_heating_assuming_clear_sky_on_radiation_time_step from CCPP data structure')
            return
        end if
        if (kind(htlw0).ne.ckind) then
            call ccpp_error('Kind mismatch for variable tendency_of_air_temperature_due_to_longwave_heating_assuming_clear_sky_on_radiation_time_step')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call rrtmg_lw_post_run(Model=Model,Grid=Grid,Radtend=Radtend,Coupling=Coupling,im=im,ltp=ltp,lm=lm, &
                  kd=kd,tsfa=tsfa,htlwc=htlwc,htlw0=htlw0,errmsg=cdata%errmsg,errflg=cdata%errflg)
        ierr=cdata%errflg

    end function rrtmg_lw_post_run_cap

    function rrtmg_lw_post_finalize_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call rrtmg_lw_post_finalize()
        

    end function rrtmg_lw_post_finalize_cap
end module rrtmg_lw_post_cap

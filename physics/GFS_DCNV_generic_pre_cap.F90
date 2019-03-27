
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
!! @brief Auto-generated cap module for the GFS_DCNV_generic_pre scheme
!!
!
module GFS_DCNV_generic_pre_cap

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_f_pointer, c_ptr, c_int32_t
    use            :: ccpp_types,                                      &
                      only: ccpp_t, CCPP_GENERIC_KIND
    use            :: ccpp_fields,                                     &
                      only: ccpp_field_get
    use            :: ccpp_errors,                                     &
                      only: ccpp_error, ccpp_debug
    use            :: GFS_DCNV_generic_pre, &
                      only: GFS_DCNV_generic_pre_run,GFS_DCNV_generic_pre_init,GFS_DCNV_generic_pre_finalize
    ! Other modules required, e.g. type definitions
    use machine, only: kind_phys

    implicit none

    private
    public :: GFS_DCNV_generic_pre_run_cap,GFS_DCNV_generic_pre_init_cap,GFS_DCNV_generic_pre_finalize_cap

    contains


    function GFS_DCNV_generic_pre_run_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind
        integer, pointer :: im
        integer, pointer :: levs
        logical, pointer :: ldiag3d
        logical, pointer :: cnvgwd
        logical, pointer :: lgocart
        real(kind_phys), pointer :: gu0(:,:)
        real(kind_phys), pointer :: gv0(:,:)
        real(kind_phys), pointer :: gt0(:,:)
        real(kind_phys), pointer :: gq0_water_vapor(:,:)
        real(kind_phys), pointer :: save_u(:,:)
        real(kind_phys), pointer :: save_v(:,:)
        real(kind_phys), pointer :: save_t(:,:)
        real(kind_phys), pointer :: save_qv(:,:)

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
        

        call ccpp_field_get(cdata, 'flag_convective_gravity_wave_drag', cnvgwd, ierr=ierr, kind=ckind, index=260)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_convective_gravity_wave_drag from CCPP data structure')
            return
        end if
        if (kind(cnvgwd).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_convective_gravity_wave_drag')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'flag_gocart', lgocart, ierr=ierr, kind=ckind, index=329)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_gocart from CCPP data structure')
            return
        end if
        if (kind(lgocart).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_gocart')
            ierr = 1
            return
        end if
#endif
        

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
        

        call ccpp_field_get(cdata, 'water_vapor_specific_humidity_updated_by_physics', gq0_water_vapor, ierr=ierr, dims=cdims, kind=ckind, index=860)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve water_vapor_specific_humidity_updated_by_physics from CCPP data structure')
            return
        end if
        if (kind(gq0_water_vapor).ne.ckind) then
            call ccpp_error('Kind mismatch for variable water_vapor_specific_humidity_updated_by_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'x_wind_save', save_u, ierr=ierr, dims=cdims, kind=ckind, index=876)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve x_wind_save from CCPP data structure')
            return
        end if
        if (kind(save_u).ne.ckind) then
            call ccpp_error('Kind mismatch for variable x_wind_save')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'y_wind_save', save_v, ierr=ierr, dims=cdims, kind=ckind, index=883)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve y_wind_save from CCPP data structure')
            return
        end if
        if (kind(save_v).ne.ckind) then
            call ccpp_error('Kind mismatch for variable y_wind_save')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'air_temperature_save', save_t, ierr=ierr, dims=cdims, kind=ckind, index=57)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve air_temperature_save from CCPP data structure')
            return
        end if
        if (kind(save_t).ne.ckind) then
            call ccpp_error('Kind mismatch for variable air_temperature_save')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'water_vapor_specific_humidity_save', save_qv, ierr=ierr, dims=cdims, kind=ckind, index=858)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve water_vapor_specific_humidity_save from CCPP data structure')
            return
        end if
        if (kind(save_qv).ne.ckind) then
            call ccpp_error('Kind mismatch for variable water_vapor_specific_humidity_save')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call GFS_DCNV_generic_pre_run(im=im,levs=levs,ldiag3d=ldiag3d,cnvgwd=cnvgwd,lgocart=lgocart,gu0=gu0,gv0=gv0, &
                  gt0=gt0,gq0_water_vapor=gq0_water_vapor,save_u=save_u,save_v=save_v,save_t=save_t, &
                  save_qv=save_qv,errmsg=cdata%errmsg,errflg=cdata%errflg)
        ierr=cdata%errflg

    end function GFS_DCNV_generic_pre_run_cap

    function GFS_DCNV_generic_pre_init_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call GFS_DCNV_generic_pre_init()
        

    end function GFS_DCNV_generic_pre_init_cap

    function GFS_DCNV_generic_pre_finalize_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call GFS_DCNV_generic_pre_finalize()
        

    end function GFS_DCNV_generic_pre_finalize_cap
end module GFS_DCNV_generic_pre_cap

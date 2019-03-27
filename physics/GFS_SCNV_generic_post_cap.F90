
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
!! @brief Auto-generated cap module for the GFS_SCNV_generic_post scheme
!!
!
module GFS_SCNV_generic_post_cap

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_f_pointer, c_ptr, c_int32_t
    use            :: ccpp_types,                                      &
                      only: ccpp_t, CCPP_GENERIC_KIND
    use            :: ccpp_fields,                                     &
                      only: ccpp_field_get
    use            :: ccpp_errors,                                     &
                      only: ccpp_error, ccpp_debug
    use            :: GFS_SCNV_generic_post, &
                      only: GFS_SCNV_generic_post_finalize,GFS_SCNV_generic_post_run,GFS_SCNV_generic_post_init
    ! Other modules required, e.g. type definitions
    use machine, only: kind_phys

    implicit none

    private
    public :: GFS_SCNV_generic_post_finalize_cap,GFS_SCNV_generic_post_run_cap,GFS_SCNV_generic_post_init_cap

    contains


    function GFS_SCNV_generic_post_finalize_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call GFS_SCNV_generic_post_finalize()
        

    end function GFS_SCNV_generic_post_finalize_cap

    function GFS_SCNV_generic_post_run_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind
        integer, pointer :: im
        integer, pointer :: levs
        integer, pointer :: nn
        logical, pointer :: lssav
        logical, pointer :: ldiag3d
        logical, pointer :: lgocart
        real(kind_phys), pointer :: frain
        real(kind_phys), pointer :: gt0(:,:)
        real(kind_phys), pointer :: gq0_water_vapor(:,:)
        real(kind_phys), pointer :: save_t(:,:)
        real(kind_phys), pointer :: save_qv(:,:)
        real(kind_phys), pointer :: dqdti(:,:)
        real(kind_phys), pointer :: dt3dt(:,:)
        real(kind_phys), pointer :: dq3dt(:,:)
        real(kind_phys), pointer :: clw(:,:,:)

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
        

        call ccpp_field_get(cdata, 'number_of_tracers_for_convective_transport', nn, ierr=ierr, kind=ckind, index=587)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve number_of_tracers_for_convective_transport from CCPP data structure')
            return
        end if
        if (kind(nn).ne.ckind) then
            call ccpp_error('Kind mismatch for variable number_of_tracers_for_convective_transport')
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
        

        call ccpp_field_get(cdata, 'dynamics_to_physics_timestep_ratio', frain, ierr=ierr, kind=ckind, index=236)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve dynamics_to_physics_timestep_ratio from CCPP data structure')
            return
        end if
        if (kind(frain).ne.ckind) then
            call ccpp_error('Kind mismatch for variable dynamics_to_physics_timestep_ratio')
            ierr = 1
            return
        end if
#endif
        

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
        

        call ccpp_field_get(cdata, 'instantaneous_water_vapor_specific_humidity_tendency_due_to_convection', dqdti, ierr=ierr, dims=cdims, kind=ckind, index=444)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve instantaneous_water_vapor_specific_humidity_tendency_due_to_convection from CCPP data structure')
            return
        end if
        if (kind(dqdti).ne.ckind) then
            call ccpp_error('Kind mismatch for variable instantaneous_water_vapor_specific_humidity_tendency_due_to_convection')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cumulative_change_in_temperature_due_to_shal_convection', dt3dt, ierr=ierr, dims=cdims, kind=ckind, index=159)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cumulative_change_in_temperature_due_to_shal_convection from CCPP data structure')
            return
        end if
        if (kind(dt3dt).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cumulative_change_in_temperature_due_to_shal_convection')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cumulative_change_in_water_vapor_specific_humidity_due_to_shal_convection', dq3dt, ierr=ierr, dims=cdims, kind=ckind, index=164)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cumulative_change_in_water_vapor_specific_humidity_due_to_shal_convection from CCPP data structure')
            return
        end if
        if (kind(dq3dt).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cumulative_change_in_water_vapor_specific_humidity_due_to_shal_convection')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'convective_transportable_tracers', clw, ierr=ierr, dims=cdims, kind=ckind, index=130)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve convective_transportable_tracers from CCPP data structure')
            return
        end if
        if (kind(clw).ne.ckind) then
            call ccpp_error('Kind mismatch for variable convective_transportable_tracers')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call GFS_SCNV_generic_post_run(im=im,levs=levs,nn=nn,lssav=lssav,ldiag3d=ldiag3d,lgocart=lgocart,frain=frain, &
                  gt0=gt0,gq0_water_vapor=gq0_water_vapor,save_t=save_t,save_qv=save_qv,dqdti=dqdti, &
                  dt3dt=dt3dt,dq3dt=dq3dt,clw=clw,errmsg=cdata%errmsg,errflg=cdata%errflg)
        ierr=cdata%errflg

    end function GFS_SCNV_generic_post_run_cap

    function GFS_SCNV_generic_post_init_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call GFS_SCNV_generic_post_init()
        

    end function GFS_SCNV_generic_post_init_cap
end module GFS_SCNV_generic_post_cap

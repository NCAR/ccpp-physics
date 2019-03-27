
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
!! @brief Auto-generated cap module for the m_micro_post scheme
!!
!
module m_micro_post_cap

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_f_pointer, c_ptr, c_int32_t
    use            :: ccpp_types,                                      &
                      only: ccpp_t, CCPP_GENERIC_KIND
    use            :: ccpp_fields,                                     &
                      only: ccpp_field_get
    use            :: ccpp_errors,                                     &
                      only: ccpp_error, ccpp_debug
    use            :: m_micro_post, &
                      only: m_micro_post_init,m_micro_post_run,m_micro_post_finalize
    ! Other modules required, e.g. type definitions
    use machine, only: kind_phys

    implicit none

    private
    public :: m_micro_post_init_cap,m_micro_post_run_cap,m_micro_post_finalize_cap

    contains


    function m_micro_post_init_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call m_micro_post_init()
        

    end function m_micro_post_init_cap

    function m_micro_post_run_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind
        integer, pointer :: im
        integer, pointer :: levs
        integer, pointer :: fprcp
        logical, pointer :: mg3_as_mg2
        real(kind_phys), pointer :: ncpr(:,:)
        real(kind_phys), pointer :: ncps(:,:)
        real(kind_phys), pointer :: ncgl(:,:)
        real(kind_phys), pointer :: qrn(:,:)
        real(kind_phys), pointer :: qsnw(:,:)
        real(kind_phys), pointer :: qgl(:,:)
        real(kind_phys), pointer :: gq0_rain(:,:)
        real(kind_phys), pointer :: gq0_snow(:,:)
        real(kind_phys), pointer :: gq0_graupel(:,:)
        real(kind_phys), pointer :: gq0_rain_nc(:,:)
        real(kind_phys), pointer :: gq0_snow_nc(:,:)
        real(kind_phys), pointer :: gq0_graupel_nc(:,:)

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
        

        call ccpp_field_get(cdata, 'number_of_frozen_precipitation_species', fprcp, ierr=ierr, kind=ckind, index=577)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve number_of_frozen_precipitation_species from CCPP data structure')
            return
        end if
        if (kind(fprcp).ne.ckind) then
            call ccpp_error('Kind mismatch for variable number_of_frozen_precipitation_species')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'flag_mg3_as_mg2', mg3_as_mg2, ierr=ierr, kind=ckind, index=331)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_mg3_as_mg2 from CCPP data structure')
            return
        end if
        if (kind(mg3_as_mg2).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_mg3_as_mg2')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'local_rain_number_concentration', ncpr, ierr=ierr, dims=cdims, kind=ckind, index=469)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve local_rain_number_concentration from CCPP data structure')
            return
        end if
        if (kind(ncpr).ne.ckind) then
            call ccpp_error('Kind mismatch for variable local_rain_number_concentration')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'local_snow_number_concentration', ncps, ierr=ierr, dims=cdims, kind=ckind, index=471)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve local_snow_number_concentration from CCPP data structure')
            return
        end if
        if (kind(ncps).ne.ckind) then
            call ccpp_error('Kind mismatch for variable local_snow_number_concentration')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'local_graupel_number_concentration', ncgl, ierr=ierr, dims=cdims, kind=ckind, index=468)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve local_graupel_number_concentration from CCPP data structure')
            return
        end if
        if (kind(ncgl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable local_graupel_number_concentration')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'local_rain_water_mixing_ratio', qrn, ierr=ierr, dims=cdims, kind=ckind, index=470)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve local_rain_water_mixing_ratio from CCPP data structure')
            return
        end if
        if (kind(qrn).ne.ckind) then
            call ccpp_error('Kind mismatch for variable local_rain_water_mixing_ratio')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'local_snow_water_mixing_ratio', qsnw, ierr=ierr, dims=cdims, kind=ckind, index=472)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve local_snow_water_mixing_ratio from CCPP data structure')
            return
        end if
        if (kind(qsnw).ne.ckind) then
            call ccpp_error('Kind mismatch for variable local_snow_water_mixing_ratio')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'local_graupel_mixing_ratio', qgl, ierr=ierr, dims=cdims, kind=ckind, index=467)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve local_graupel_mixing_ratio from CCPP data structure')
            return
        end if
        if (kind(qgl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable local_graupel_mixing_ratio')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'rain_water_mixing_ratio_updated_by_physics', gq0_rain, ierr=ierr, dims=cdims, kind=ckind, index=619)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve rain_water_mixing_ratio_updated_by_physics from CCPP data structure')
            return
        end if
        if (kind(gq0_rain).ne.ckind) then
            call ccpp_error('Kind mismatch for variable rain_water_mixing_ratio_updated_by_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'snow_water_mixing_ratio_updated_by_physics', gq0_snow, ierr=ierr, dims=cdims, kind=ckind, index=653)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve snow_water_mixing_ratio_updated_by_physics from CCPP data structure')
            return
        end if
        if (kind(gq0_snow).ne.ckind) then
            call ccpp_error('Kind mismatch for variable snow_water_mixing_ratio_updated_by_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'graupel_mixing_ratio_updated_by_physics', gq0_graupel, ierr=ierr, dims=cdims, kind=ckind, index=352)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve graupel_mixing_ratio_updated_by_physics from CCPP data structure')
            return
        end if
        if (kind(gq0_graupel).ne.ckind) then
            call ccpp_error('Kind mismatch for variable graupel_mixing_ratio_updated_by_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'rain_number_concentration_updated_by_physics', gq0_rain_nc, ierr=ierr, dims=cdims, kind=ckind, index=618)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve rain_number_concentration_updated_by_physics from CCPP data structure')
            return
        end if
        if (kind(gq0_rain_nc).ne.ckind) then
            call ccpp_error('Kind mismatch for variable rain_number_concentration_updated_by_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'snow_number_concentration_updated_by_physics', gq0_snow_nc, ierr=ierr, dims=cdims, kind=ckind, index=651)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve snow_number_concentration_updated_by_physics from CCPP data structure')
            return
        end if
        if (kind(gq0_snow_nc).ne.ckind) then
            call ccpp_error('Kind mismatch for variable snow_number_concentration_updated_by_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'graupel_number_concentration_updated_by_physics', gq0_graupel_nc, ierr=ierr, dims=cdims, kind=ckind, index=353)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve graupel_number_concentration_updated_by_physics from CCPP data structure')
            return
        end if
        if (kind(gq0_graupel_nc).ne.ckind) then
            call ccpp_error('Kind mismatch for variable graupel_number_concentration_updated_by_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call m_micro_post_run(im=im,levs=levs,fprcp=fprcp,mg3_as_mg2=mg3_as_mg2,ncpr=ncpr,ncps=ncps,ncgl=ncgl, &
                  qrn=qrn,qsnw=qsnw,qgl=qgl,gq0_rain=gq0_rain,gq0_snow=gq0_snow,gq0_graupel=gq0_graupel, &
                  gq0_rain_nc=gq0_rain_nc,gq0_snow_nc=gq0_snow_nc,gq0_graupel_nc=gq0_graupel_nc, &
                  errmsg=cdata%errmsg,errflg=cdata%errflg)
        ierr=cdata%errflg

    end function m_micro_post_run_cap

    function m_micro_post_finalize_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call m_micro_post_finalize()
        

    end function m_micro_post_finalize_cap
end module m_micro_post_cap

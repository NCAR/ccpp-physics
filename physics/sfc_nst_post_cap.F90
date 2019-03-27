
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
!! @brief Auto-generated cap module for the sfc_nst_post scheme
!!
!
module sfc_nst_post_cap

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_f_pointer, c_ptr, c_int32_t
    use            :: ccpp_types,                                      &
                      only: ccpp_t, CCPP_GENERIC_KIND
    use            :: ccpp_fields,                                     &
                      only: ccpp_field_get
    use            :: ccpp_errors,                                     &
                      only: ccpp_error, ccpp_debug
    use            :: sfc_nst_post, &
                      only: sfc_nst_post_init,sfc_nst_post_run,sfc_nst_post_finalize
    ! Other modules required, e.g. type definitions
    use machine, only: kind_phys

    implicit none

    private
    public :: sfc_nst_post_init_cap,sfc_nst_post_run_cap,sfc_nst_post_finalize_cap

    contains


    function sfc_nst_post_init_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call sfc_nst_post_init()
        

    end function sfc_nst_post_init_cap

    function sfc_nst_post_run_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind
        integer, pointer :: im
        integer, pointer :: islimsk(:)
        real(kind_phys), pointer :: oro(:)
        real(kind_phys), pointer :: oro_uf(:)
        integer, pointer :: nstf_name1
        integer, pointer :: nstf_name4
        integer, pointer :: nstf_name5
        real(kind_phys), pointer :: xt(:)
        real(kind_phys), pointer :: xz(:)
        real(kind_phys), pointer :: dt_cool(:)
        real(kind_phys), pointer :: z_c(:)
        real(kind_phys), pointer :: rslimsk(:)
        real(kind_phys), pointer :: tref(:)
        real(kind_phys), pointer :: xlon(:)
        real(kind_phys), pointer :: tsurf(:)
        real(kind_phys), pointer :: dtzm(:)
        real(kind_phys), pointer :: tsfc(:)

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
        

        call ccpp_field_get(cdata, 'sea_land_ice_mask', islimsk, ierr=ierr, dims=cdims, kind=ckind, index=631)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve sea_land_ice_mask from CCPP data structure')
            return
        end if
        if (kind(islimsk).ne.ckind) then
            call ccpp_error('Kind mismatch for variable sea_land_ice_mask')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'orography', oro, ierr=ierr, dims=cdims, kind=ckind, index=595)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve orography from CCPP data structure')
            return
        end if
        if (kind(oro).ne.ckind) then
            call ccpp_error('Kind mismatch for variable orography')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'orography_unfiltered', oro_uf, ierr=ierr, dims=cdims, kind=ckind, index=596)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve orography_unfiltered from CCPP data structure')
            return
        end if
        if (kind(oro_uf).ne.ckind) then
            call ccpp_error('Kind mismatch for variable orography_unfiltered')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'flag_for_nsstm_run', nstf_name1, ierr=ierr, kind=ckind, index=299)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_nsstm_run from CCPP data structure')
            return
        end if
        if (kind(nstf_name1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_nsstm_run')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'vertical_temperature_average_range_lower_bound', nstf_name4, ierr=ierr, kind=ckind, index=828)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve vertical_temperature_average_range_lower_bound from CCPP data structure')
            return
        end if
        if (kind(nstf_name4).ne.ckind) then
            call ccpp_error('Kind mismatch for variable vertical_temperature_average_range_lower_bound')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'vertical_temperature_average_range_upper_bound', nstf_name5, ierr=ierr, kind=ckind, index=829)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve vertical_temperature_average_range_upper_bound from CCPP data structure')
            return
        end if
        if (kind(nstf_name5).ne.ckind) then
            call ccpp_error('Kind mismatch for variable vertical_temperature_average_range_upper_bound')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'diurnal_thermocline_layer_heat_content', xt, ierr=ierr, dims=cdims, kind=ckind, index=224)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve diurnal_thermocline_layer_heat_content from CCPP data structure')
            return
        end if
        if (kind(xt).ne.ckind) then
            call ccpp_error('Kind mismatch for variable diurnal_thermocline_layer_heat_content')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'diurnal_thermocline_layer_thickness', xz, ierr=ierr, dims=cdims, kind=ckind, index=225)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve diurnal_thermocline_layer_thickness from CCPP data structure')
            return
        end if
        if (kind(xz).ne.ckind) then
            call ccpp_error('Kind mismatch for variable diurnal_thermocline_layer_thickness')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'sub-layer_cooling_amount', dt_cool, ierr=ierr, dims=cdims, kind=ckind, index=671)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve sub-layer_cooling_amount from CCPP data structure')
            return
        end if
        if (kind(dt_cool).ne.ckind) then
            call ccpp_error('Kind mismatch for variable sub-layer_cooling_amount')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'sub-layer_cooling_thickness', z_c, ierr=ierr, dims=cdims, kind=ckind, index=672)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve sub-layer_cooling_thickness from CCPP data structure')
            return
        end if
        if (kind(z_c).ne.ckind) then
            call ccpp_error('Kind mismatch for variable sub-layer_cooling_thickness')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'sea_land_ice_mask_real', rslimsk, ierr=ierr, dims=cdims, kind=ckind, index=632)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve sea_land_ice_mask_real from CCPP data structure')
            return
        end if
        if (kind(rslimsk).ne.ckind) then
            call ccpp_error('Kind mismatch for variable sea_land_ice_mask_real')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'sea_surface_reference_temperature', tref, ierr=ierr, dims=cdims, kind=ckind, index=633)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve sea_surface_reference_temperature from CCPP data structure')
            return
        end if
        if (kind(tref).ne.ckind) then
            call ccpp_error('Kind mismatch for variable sea_surface_reference_temperature')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'longitude', xlon, ierr=ierr, dims=cdims, kind=ckind, index=473)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve longitude from CCPP data structure')
            return
        end if
        if (kind(xlon).ne.ckind) then
            call ccpp_error('Kind mismatch for variable longitude')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_skin_temperature_after_iteration', tsurf, ierr=ierr, dims=cdims, kind=ckind, index=722)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_skin_temperature_after_iteration from CCPP data structure')
            return
        end if
        if (kind(tsurf).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_skin_temperature_after_iteration')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'mean_change_over_depth_in_sea_water_temperature', dtzm, ierr=ierr, dims=cdims, kind=ckind, index=514)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve mean_change_over_depth_in_sea_water_temperature from CCPP data structure')
            return
        end if
        if (kind(dtzm).ne.ckind) then
            call ccpp_error('Kind mismatch for variable mean_change_over_depth_in_sea_water_temperature')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_skin_temperature', tsfc, ierr=ierr, dims=cdims, kind=ckind, index=721)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_skin_temperature from CCPP data structure')
            return
        end if
        if (kind(tsfc).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_skin_temperature')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call sfc_nst_post_run(im=im,islimsk=islimsk,oro=oro,oro_uf=oro_uf,nstf_name1=nstf_name1,nstf_name4=nstf_name4, &
                  nstf_name5=nstf_name5,xt=xt,xz=xz,dt_cool=dt_cool,z_c=z_c,rslimsk=rslimsk, &
                  tref=tref,xlon=xlon,tsurf=tsurf,dtzm=dtzm,tsfc=tsfc,errmsg=cdata%errmsg, &
                  errflg=cdata%errflg)
        ierr=cdata%errflg

    end function sfc_nst_post_run_cap

    function sfc_nst_post_finalize_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call sfc_nst_post_finalize()
        

    end function sfc_nst_post_finalize_cap
end module sfc_nst_post_cap

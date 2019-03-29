
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
!! @brief Auto-generated cap module for the mynnrad_pre scheme
!!
!
module mynnrad_pre_cap

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_f_pointer, c_ptr, c_int32_t
    use            :: ccpp_types,                                      &
                      only: ccpp_t, CCPP_GENERIC_KIND
    use            :: ccpp_fields,                                     &
                      only: ccpp_field_get
    use            :: ccpp_errors,                                     &
                      only: ccpp_error, ccpp_debug
    use            :: mynnrad_pre, &
                      only: mynnrad_pre_init,mynnrad_pre_run,mynnrad_pre_finalize
    ! Other modules required, e.g. type definitions
    use machine, only: kind_phys

    implicit none

    private
    public :: mynnrad_pre_init_cap,mynnrad_pre_run_cap,mynnrad_pre_finalize_cap

    contains


    function mynnrad_pre_init_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call mynnrad_pre_init()
        

    end function mynnrad_pre_init_cap

    function mynnrad_pre_run_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind
        integer, pointer :: ix
        integer, pointer :: im
        integer, pointer :: levs
        real(kind_phys), pointer :: qc(:,:)
        real(kind_phys), pointer :: qi(:,:)
        real(kind_phys), pointer :: T3D(:,:)
        real(kind_phys), pointer :: qc_save(:,:)
        real(kind_phys), pointer :: qi_save(:,:)
        real(kind_phys), pointer :: QC_BL(:,:)
        real(kind_phys), pointer :: CLDFRA_BL(:,:)
        real(kind_phys), pointer :: delp(:,:)
        real(kind_phys), pointer :: clouds1(:,:)
        real(kind_phys), pointer :: clouds2(:,:)
        real(kind_phys), pointer :: clouds3(:,:)
        real(kind_phys), pointer :: clouds4(:,:)
        real(kind_phys), pointer :: clouds5(:,:)
        real(kind_phys), pointer :: slmsk(:)

        ierr = 0

        call c_f_pointer(ptr, cdata)


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
        

        call ccpp_field_get(cdata, 'cloud_condensed_water_mixing_ratio', qc, ierr=ierr, dims=cdims, kind=ckind, index=90)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cloud_condensed_water_mixing_ratio from CCPP data structure')
            return
        end if
        if (kind(qc).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cloud_condensed_water_mixing_ratio')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'ice_water_mixing_ratio', qi, ierr=ierr, dims=cdims, kind=ckind, index=373)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve ice_water_mixing_ratio from CCPP data structure')
            return
        end if
        if (kind(qi).ne.ckind) then
            call ccpp_error('Kind mismatch for variable ice_water_mixing_ratio')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'air_temperature', T3D, ierr=ierr, dims=cdims, kind=ckind, index=50)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve air_temperature from CCPP data structure')
            return
        end if
        if (kind(T3D).ne.ckind) then
            call ccpp_error('Kind mismatch for variable air_temperature')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cloud_condensed_water_mixing_ratio_save', qc_save, ierr=ierr, dims=cdims, kind=ckind, index=94)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cloud_condensed_water_mixing_ratio_save from CCPP data structure')
            return
        end if
        if (kind(qc_save).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cloud_condensed_water_mixing_ratio_save')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'ice_water_mixing_ratio_save', qi_save, ierr=ierr, dims=cdims, kind=ckind, index=375)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve ice_water_mixing_ratio_save from CCPP data structure')
            return
        end if
        if (kind(qi_save).ne.ckind) then
            call ccpp_error('Kind mismatch for variable ice_water_mixing_ratio_save')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'subgrid_cloud_mixing_ratio_pbl', QC_BL, ierr=ierr, dims=cdims, kind=ckind, index=674)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve subgrid_cloud_mixing_ratio_pbl from CCPP data structure')
            return
        end if
        if (kind(QC_BL).ne.ckind) then
            call ccpp_error('Kind mismatch for variable subgrid_cloud_mixing_ratio_pbl')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'subgrid_cloud_fraction_pbl', CLDFRA_BL, ierr=ierr, dims=cdims, kind=ckind, index=673)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve subgrid_cloud_fraction_pbl from CCPP data structure')
            return
        end if
        if (kind(CLDFRA_BL).ne.ckind) then
            call ccpp_error('Kind mismatch for variable subgrid_cloud_fraction_pbl')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'layer_pressure_thickness_for_radiation', delp, ierr=ierr, dims=cdims, kind=ckind, index=463)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve layer_pressure_thickness_for_radiation from CCPP data structure')
            return
        end if
        if (kind(delp).ne.ckind) then
            call ccpp_error('Kind mismatch for variable layer_pressure_thickness_for_radiation')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'total_cloud_fraction', clouds1, ierr=ierr, dims=cdims, kind=ckind, index=799)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve total_cloud_fraction from CCPP data structure')
            return
        end if
        if (kind(clouds1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable total_cloud_fraction')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cloud_liquid_water_path', clouds2, ierr=ierr, dims=cdims, kind=ckind, index=102)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cloud_liquid_water_path from CCPP data structure')
            return
        end if
        if (kind(clouds2).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cloud_liquid_water_path')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'mean_effective_radius_for_liquid_cloud', clouds3, ierr=ierr, dims=cdims, kind=ckind, index=516)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve mean_effective_radius_for_liquid_cloud from CCPP data structure')
            return
        end if
        if (kind(clouds3).ne.ckind) then
            call ccpp_error('Kind mismatch for variable mean_effective_radius_for_liquid_cloud')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cloud_ice_water_path', clouds4, ierr=ierr, dims=cdims, kind=ckind, index=101)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cloud_ice_water_path from CCPP data structure')
            return
        end if
        if (kind(clouds4).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cloud_ice_water_path')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'mean_effective_radius_for_ice_cloud', clouds5, ierr=ierr, dims=cdims, kind=ckind, index=515)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve mean_effective_radius_for_ice_cloud from CCPP data structure')
            return
        end if
        if (kind(clouds5).ne.ckind) then
            call ccpp_error('Kind mismatch for variable mean_effective_radius_for_ice_cloud')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'sea_land_ice_mask_real', slmsk, ierr=ierr, dims=cdims, kind=ckind, index=632)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve sea_land_ice_mask_real from CCPP data structure')
            return
        end if
        if (kind(slmsk).ne.ckind) then
            call ccpp_error('Kind mismatch for variable sea_land_ice_mask_real')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call mynnrad_pre_run(ix=ix,im=im,levs=levs,qc=qc,qi=qi,T3D=T3D,qc_save=qc_save,qi_save=qi_save, &
                  QC_BL=QC_BL,CLDFRA_BL=CLDFRA_BL,delp=delp,clouds1=clouds1,clouds2=clouds2, &
                  clouds3=clouds3,clouds4=clouds4,clouds5=clouds5,slmsk=slmsk,errmsg=cdata%errmsg, &
                  errflg=cdata%errflg)
        ierr=cdata%errflg

    end function mynnrad_pre_run_cap

    function mynnrad_pre_finalize_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call mynnrad_pre_finalize()
        

    end function mynnrad_pre_finalize_cap
end module mynnrad_pre_cap

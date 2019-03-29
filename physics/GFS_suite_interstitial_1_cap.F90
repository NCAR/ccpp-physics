
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
!! @brief Auto-generated cap module for the GFS_suite_interstitial_1 scheme
!!
!
module GFS_suite_interstitial_1_cap

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_f_pointer, c_ptr, c_int32_t
    use            :: ccpp_types,                                      &
                      only: ccpp_t, CCPP_GENERIC_KIND
    use            :: ccpp_fields,                                     &
                      only: ccpp_field_get
    use            :: ccpp_errors,                                     &
                      only: ccpp_error, ccpp_debug
    use            :: GFS_suite_interstitial_1, &
                      only: GFS_suite_interstitial_1_init,GFS_suite_interstitial_1_finalize,GFS_suite_interstitial_1_run
    ! Other modules required, e.g. type definitions
    use machine, only: kind_phys

    implicit none

    private
    public :: GFS_suite_interstitial_1_init_cap,GFS_suite_interstitial_1_finalize_cap,GFS_suite_interstitial_1_run_cap

    contains


    function GFS_suite_interstitial_1_init_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call GFS_suite_interstitial_1_init()
        

    end function GFS_suite_interstitial_1_init_cap

    function GFS_suite_interstitial_1_finalize_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call GFS_suite_interstitial_1_finalize()
        

    end function GFS_suite_interstitial_1_finalize_cap

    function GFS_suite_interstitial_1_run_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind
        integer, pointer :: im
        integer, pointer :: levs
        integer, pointer :: ntrac
        real(kind_phys), pointer :: crtrh(:)
        real(kind_phys), pointer :: dtf
        real(kind_phys), pointer :: dtp
        real(kind_phys), pointer :: slmsk(:)
        real(kind_phys), pointer :: area(:)
        real(kind_phys), pointer :: dxmin
        real(kind_phys), pointer :: dxinv
        real(kind_phys), pointer :: pgr(:)
        real(kind_phys), pointer :: rhbbot
        real(kind_phys), pointer :: rhpbl
        real(kind_phys), pointer :: rhbtop
        real(kind_phys), pointer :: frain
        integer, pointer :: islmsk(:)
        real(kind_phys), pointer :: frland(:)
        real(kind_phys), pointer :: work1(:)
        real(kind_phys), pointer :: work2(:)
        real(kind_phys), pointer :: psurf(:)
        real(kind_phys), pointer :: dudt(:,:)
        real(kind_phys), pointer :: dvdt(:,:)
        real(kind_phys), pointer :: dtdt(:,:)
        real(kind_phys), pointer :: dtdtc(:,:)
        real(kind_phys), pointer :: dqdt(:,:,:)

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
        

        call ccpp_field_get(cdata, 'number_of_tracers', ntrac, ierr=ierr, kind=ckind, index=584)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve number_of_tracers from CCPP data structure')
            return
        end if
        if (kind(ntrac).ne.ckind) then
            call ccpp_error('Kind mismatch for variable number_of_tracers')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'critical_relative_humidity_at_sfc_pbltop_toa', crtrh, ierr=ierr, dims=cdims, kind=ckind, index=143)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve critical_relative_humidity_at_sfc_pbltop_toa from CCPP data structure')
            return
        end if
        if (kind(crtrh).ne.ckind) then
            call ccpp_error('Kind mismatch for variable critical_relative_humidity_at_sfc_pbltop_toa')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

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
        

        call ccpp_field_get(cdata, 'time_step_for_physics', dtp, ierr=ierr, kind=ckind, index=793)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve time_step_for_physics from CCPP data structure')
            return
        end if
        if (kind(dtp).ne.ckind) then
            call ccpp_error('Kind mismatch for variable time_step_for_physics')
            ierr = 1
            return
        end if
#endif
        

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
        

        call ccpp_field_get(cdata, 'cell_area', area, ierr=ierr, dims=cdims, kind=ckind, index=83)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cell_area from CCPP data structure')
            return
        end if
        if (kind(area).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cell_area')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'minimum_scaling_factor_for_critical_relative_humidity', dxmin, ierr=ierr, kind=ckind, index=544)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve minimum_scaling_factor_for_critical_relative_humidity from CCPP data structure')
            return
        end if
        if (kind(dxmin).ne.ckind) then
            call ccpp_error('Kind mismatch for variable minimum_scaling_factor_for_critical_relative_humidity')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'inverse_scaling_factor_for_critical_relative_humidity', dxinv, ierr=ierr, kind=ckind, index=449)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve inverse_scaling_factor_for_critical_relative_humidity from CCPP data structure')
            return
        end if
        if (kind(dxinv).ne.ckind) then
            call ccpp_error('Kind mismatch for variable inverse_scaling_factor_for_critical_relative_humidity')
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
        

        call ccpp_field_get(cdata, 'critical_relative_humidity_at_surface', rhbbot, ierr=ierr, kind=ckind, index=144)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve critical_relative_humidity_at_surface from CCPP data structure')
            return
        end if
        if (kind(rhbbot).ne.ckind) then
            call ccpp_error('Kind mismatch for variable critical_relative_humidity_at_surface')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'critical_relative_humidity_at_PBL_top', rhpbl, ierr=ierr, kind=ckind, index=142)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve critical_relative_humidity_at_PBL_top from CCPP data structure')
            return
        end if
        if (kind(rhpbl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable critical_relative_humidity_at_PBL_top')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'critical_relative_humidity_at_top_of_atmosphere', rhbtop, ierr=ierr, kind=ckind, index=145)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve critical_relative_humidity_at_top_of_atmosphere from CCPP data structure')
            return
        end if
        if (kind(rhbtop).ne.ckind) then
            call ccpp_error('Kind mismatch for variable critical_relative_humidity_at_top_of_atmosphere')
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
        

        call ccpp_field_get(cdata, 'sea_land_ice_mask', islmsk, ierr=ierr, dims=cdims, kind=ckind, index=631)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve sea_land_ice_mask from CCPP data structure')
            return
        end if
        if (kind(islmsk).ne.ckind) then
            call ccpp_error('Kind mismatch for variable sea_land_ice_mask')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'land_area_fraction', frland, ierr=ierr, dims=cdims, kind=ckind, index=456)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve land_area_fraction from CCPP data structure')
            return
        end if
        if (kind(frland).ne.ckind) then
            call ccpp_error('Kind mismatch for variable land_area_fraction')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'grid_size_related_coefficient_used_in_scale-sensitive_schemes', work1, ierr=ierr, dims=cdims, kind=ckind, index=357)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve grid_size_related_coefficient_used_in_scale-sensitive_schemes from CCPP data structure')
            return
        end if
        if (kind(work1).ne.ckind) then
            call ccpp_error('Kind mismatch for variable grid_size_related_coefficient_used_in_scale-sensitive_schemes')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'grid_size_related_coefficient_used_in_scale-sensitive_schemes_complement', work2, ierr=ierr, dims=cdims, kind=ckind, index=358)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve grid_size_related_coefficient_used_in_scale-sensitive_schemes_complement from CCPP data structure')
            return
        end if
        if (kind(work2).ne.ckind) then
            call ccpp_error('Kind mismatch for variable grid_size_related_coefficient_used_in_scale-sensitive_schemes_complement')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'surface_air_pressure_diag', psurf, ierr=ierr, dims=cdims, kind=ckind, index=679)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve surface_air_pressure_diag from CCPP data structure')
            return
        end if
        if (kind(psurf).ne.ckind) then
            call ccpp_error('Kind mismatch for variable surface_air_pressure_diag')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'tendency_of_x_wind_due_to_model_physics', dudt, ierr=ierr, dims=cdims, kind=ckind, index=782)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve tendency_of_x_wind_due_to_model_physics from CCPP data structure')
            return
        end if
        if (kind(dudt).ne.ckind) then
            call ccpp_error('Kind mismatch for variable tendency_of_x_wind_due_to_model_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'tendency_of_y_wind_due_to_model_physics', dvdt, ierr=ierr, dims=cdims, kind=ckind, index=785)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve tendency_of_y_wind_due_to_model_physics from CCPP data structure')
            return
        end if
        if (kind(dvdt).ne.ckind) then
            call ccpp_error('Kind mismatch for variable tendency_of_y_wind_due_to_model_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'tendency_of_air_temperature_due_to_model_physics', dtdt, ierr=ierr, dims=cdims, kind=ckind, index=758)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve tendency_of_air_temperature_due_to_model_physics from CCPP data structure')
            return
        end if
        if (kind(dtdt).ne.ckind) then
            call ccpp_error('Kind mismatch for variable tendency_of_air_temperature_due_to_model_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'tendency_of_air_temperature_due_to_radiative_heating_assuming_clear_sky', dtdtc, ierr=ierr, dims=cdims, kind=ckind, index=759)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve tendency_of_air_temperature_due_to_radiative_heating_assuming_clear_sky from CCPP data structure')
            return
        end if
        if (kind(dtdtc).ne.ckind) then
            call ccpp_error('Kind mismatch for variable tendency_of_air_temperature_due_to_radiative_heating_assuming_clear_sky')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'tendency_of_tracers_due_to_model_physics', dqdt, ierr=ierr, dims=cdims, kind=ckind, index=776)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve tendency_of_tracers_due_to_model_physics from CCPP data structure')
            return
        end if
        if (kind(dqdt).ne.ckind) then
            call ccpp_error('Kind mismatch for variable tendency_of_tracers_due_to_model_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call GFS_suite_interstitial_1_run(im=im,levs=levs,ntrac=ntrac,crtrh=crtrh,dtf=dtf,dtp=dtp,slmsk=slmsk,area=area, &
                  dxmin=dxmin,dxinv=dxinv,pgr=pgr,rhbbot=rhbbot,rhpbl=rhpbl,rhbtop=rhbtop, &
                  frain=frain,islmsk=islmsk,frland=frland,work1=work1,work2=work2,psurf=psurf, &
                  dudt=dudt,dvdt=dvdt,dtdt=dtdt,dtdtc=dtdtc,dqdt=dqdt,errmsg=cdata%errmsg, &
                  errflg=cdata%errflg)
        ierr=cdata%errflg

    end function GFS_suite_interstitial_1_run_cap
end module GFS_suite_interstitial_1_cap


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
!! @brief Auto-generated cap module for the GFS_suite_interstitial_4 scheme
!!
!
module GFS_suite_interstitial_4_cap

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_f_pointer, c_ptr, c_int32_t
    use            :: ccpp_types,                                      &
                      only: ccpp_t, CCPP_GENERIC_KIND
    use            :: ccpp_fields,                                     &
                      only: ccpp_field_get
    use            :: ccpp_errors,                                     &
                      only: ccpp_error, ccpp_debug
    use            :: GFS_suite_interstitial_4, &
                      only: GFS_suite_interstitial_4_finalize,GFS_suite_interstitial_4_init,GFS_suite_interstitial_4_run
    ! Other modules required, e.g. type definitions
    use machine, only: kind_phys

    implicit none

    private
    public :: GFS_suite_interstitial_4_finalize_cap,GFS_suite_interstitial_4_init_cap,GFS_suite_interstitial_4_run_cap

    contains


    function GFS_suite_interstitial_4_finalize_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call GFS_suite_interstitial_4_finalize()
        

    end function GFS_suite_interstitial_4_finalize_cap

    function GFS_suite_interstitial_4_init_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call GFS_suite_interstitial_4_init()
        

    end function GFS_suite_interstitial_4_init_cap

    function GFS_suite_interstitial_4_run_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind
        integer, pointer :: im
        integer, pointer :: levs
        logical, pointer :: ltaerosol
        logical, pointer :: lgocart
        integer, pointer :: tracers_total
        integer, pointer :: ntrac
        integer, pointer :: ntcw
        integer, pointer :: ntiw
        integer, pointer :: ntclamt
        integer, pointer :: ntrw
        integer, pointer :: ntsw
        integer, pointer :: ntrnc
        integer, pointer :: ntsnc
        integer, pointer :: ntgl
        integer, pointer :: ntgnc
        integer, pointer :: ntlnc
        integer, pointer :: ntinc
        integer, pointer :: nn
        integer, pointer :: imp_physics
        integer, pointer :: imp_physics_gfdl
        integer, pointer :: imp_physics_thompson
        integer, pointer :: imp_physics_zhao_carr
        integer, pointer :: imp_physics_zhao_carr_pdf
        real(kind_phys), pointer :: dtf
        real(kind_phys), pointer :: save_qc(:,:)
        real(kind_phys), pointer :: save_qi(:,:)
        real(kind_phys), pointer :: con_pi
        real(kind_phys), pointer :: gq0(:,:,:)
        real(kind_phys), pointer :: clw(:,:,:)
        real(kind_phys), pointer :: dqdti(:,:)

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
        

        call ccpp_field_get(cdata, 'flag_for_aerosol_physics', ltaerosol, ierr=ierr, kind=ckind, index=271)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_aerosol_physics from CCPP data structure')
            return
        end if
        if (kind(ltaerosol).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_aerosol_physics')
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
        

        call ccpp_field_get(cdata, 'number_of_total_tracers', tracers_total, ierr=ierr, kind=ckind, index=583)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve number_of_total_tracers from CCPP data structure')
            return
        end if
        if (kind(tracers_total).ne.ckind) then
            call ccpp_error('Kind mismatch for variable number_of_total_tracers')
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
        

        call ccpp_field_get(cdata, 'index_for_liquid_cloud_condensate', ntcw, ierr=ierr, kind=ckind, index=384)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve index_for_liquid_cloud_condensate from CCPP data structure')
            return
        end if
        if (kind(ntcw).ne.ckind) then
            call ccpp_error('Kind mismatch for variable index_for_liquid_cloud_condensate')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'index_for_ice_cloud_condensate', ntiw, ierr=ierr, kind=ckind, index=381)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve index_for_ice_cloud_condensate from CCPP data structure')
            return
        end if
        if (kind(ntiw).ne.ckind) then
            call ccpp_error('Kind mismatch for variable index_for_ice_cloud_condensate')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'index_for_cloud_amount', ntclamt, ierr=ierr, kind=ckind, index=378)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve index_for_cloud_amount from CCPP data structure')
            return
        end if
        if (kind(ntclamt).ne.ckind) then
            call ccpp_error('Kind mismatch for variable index_for_cloud_amount')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'index_for_rain_water', ntrw, ierr=ierr, kind=ckind, index=388)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve index_for_rain_water from CCPP data structure')
            return
        end if
        if (kind(ntrw).ne.ckind) then
            call ccpp_error('Kind mismatch for variable index_for_rain_water')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'index_for_snow_water', ntsw, ierr=ierr, kind=ckind, index=390)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve index_for_snow_water from CCPP data structure')
            return
        end if
        if (kind(ntsw).ne.ckind) then
            call ccpp_error('Kind mismatch for variable index_for_snow_water')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'index_for_rain_number_concentration', ntrnc, ierr=ierr, kind=ckind, index=387)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve index_for_rain_number_concentration from CCPP data structure')
            return
        end if
        if (kind(ntrnc).ne.ckind) then
            call ccpp_error('Kind mismatch for variable index_for_rain_number_concentration')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'index_for_snow_number_concentration', ntsnc, ierr=ierr, kind=ckind, index=389)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve index_for_snow_number_concentration from CCPP data structure')
            return
        end if
        if (kind(ntsnc).ne.ckind) then
            call ccpp_error('Kind mismatch for variable index_for_snow_number_concentration')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'index_for_graupel', ntgl, ierr=ierr, kind=ckind, index=379)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve index_for_graupel from CCPP data structure')
            return
        end if
        if (kind(ntgl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable index_for_graupel')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'index_for_graupel_number_concentration', ntgnc, ierr=ierr, kind=ckind, index=380)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve index_for_graupel_number_concentration from CCPP data structure')
            return
        end if
        if (kind(ntgnc).ne.ckind) then
            call ccpp_error('Kind mismatch for variable index_for_graupel_number_concentration')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'index_for_liquid_cloud_number_concentration', ntlnc, ierr=ierr, kind=ckind, index=385)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve index_for_liquid_cloud_number_concentration from CCPP data structure')
            return
        end if
        if (kind(ntlnc).ne.ckind) then
            call ccpp_error('Kind mismatch for variable index_for_liquid_cloud_number_concentration')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'index_for_ice_cloud_number_concentration', ntinc, ierr=ierr, kind=ckind, index=382)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve index_for_ice_cloud_number_concentration from CCPP data structure')
            return
        end if
        if (kind(ntinc).ne.ckind) then
            call ccpp_error('Kind mismatch for variable index_for_ice_cloud_number_concentration')
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
        

        call ccpp_field_get(cdata, 'flag_for_microphysics_scheme', imp_physics, ierr=ierr, kind=ckind, index=294)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_microphysics_scheme from CCPP data structure')
            return
        end if
        if (kind(imp_physics).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_microphysics_scheme')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'flag_for_gfdl_microphysics_scheme', imp_physics_gfdl, ierr=ierr, kind=ckind, index=280)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_gfdl_microphysics_scheme from CCPP data structure')
            return
        end if
        if (kind(imp_physics_gfdl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_gfdl_microphysics_scheme')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'flag_for_thompson_microphysics_scheme', imp_physics_thompson, ierr=ierr, kind=ckind, index=322)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_thompson_microphysics_scheme from CCPP data structure')
            return
        end if
        if (kind(imp_physics_thompson).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_thompson_microphysics_scheme')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'flag_for_zhao_carr_microphysics_scheme', imp_physics_zhao_carr, ierr=ierr, kind=ckind, index=327)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_zhao_carr_microphysics_scheme from CCPP data structure')
            return
        end if
        if (kind(imp_physics_zhao_carr).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_zhao_carr_microphysics_scheme')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'flag_for_zhao_carr_pdf_microphysics_scheme', imp_physics_zhao_carr_pdf, ierr=ierr, kind=ckind, index=328)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_zhao_carr_pdf_microphysics_scheme from CCPP data structure')
            return
        end if
        if (kind(imp_physics_zhao_carr_pdf).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_zhao_carr_pdf_microphysics_scheme')
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
        

        call ccpp_field_get(cdata, 'cloud_condensed_water_mixing_ratio_save', save_qc, ierr=ierr, dims=cdims, kind=ckind, index=94)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cloud_condensed_water_mixing_ratio_save from CCPP data structure')
            return
        end if
        if (kind(save_qc).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cloud_condensed_water_mixing_ratio_save')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'ice_water_mixing_ratio_save', save_qi, ierr=ierr, dims=cdims, kind=ckind, index=375)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve ice_water_mixing_ratio_save from CCPP data structure')
            return
        end if
        if (kind(save_qi).ne.ckind) then
            call ccpp_error('Kind mismatch for variable ice_water_mixing_ratio_save')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'pi', con_pi, ierr=ierr, kind=ckind, index=606)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve pi from CCPP data structure')
            return
        end if
        if (kind(con_pi).ne.ckind) then
            call ccpp_error('Kind mismatch for variable pi')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'tracer_concentration_updated_by_physics', gq0, ierr=ierr, dims=cdims, kind=ckind, index=803)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve tracer_concentration_updated_by_physics from CCPP data structure')
            return
        end if
        if (kind(gq0).ne.ckind) then
            call ccpp_error('Kind mismatch for variable tracer_concentration_updated_by_physics')
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
        

        call GFS_suite_interstitial_4_run(im=im,levs=levs,ltaerosol=ltaerosol,lgocart=lgocart,tracers_total=tracers_total, &
                  ntrac=ntrac,ntcw=ntcw,ntiw=ntiw,ntclamt=ntclamt,ntrw=ntrw,ntsw=ntsw,ntrnc=ntrnc, &
                  ntsnc=ntsnc,ntgl=ntgl,ntgnc=ntgnc,ntlnc=ntlnc,ntinc=ntinc,nn=nn,imp_physics=imp_physics, &
                  imp_physics_gfdl=imp_physics_gfdl,imp_physics_thompson=imp_physics_thompson, &
                  imp_physics_zhao_carr=imp_physics_zhao_carr,imp_physics_zhao_carr_pdf=imp_physics_zhao_carr_pdf, &
                  dtf=dtf,save_qc=save_qc,save_qi=save_qi,con_pi=con_pi,gq0=gq0,clw=clw,dqdti=dqdti, &
                  errmsg=cdata%errmsg,errflg=cdata%errflg)
        ierr=cdata%errflg

    end function GFS_suite_interstitial_4_run_cap
end module GFS_suite_interstitial_4_cap

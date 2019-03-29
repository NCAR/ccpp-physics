
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
!! @brief Auto-generated cap module for the GFS_suite_interstitial_3 scheme
!!
!
module GFS_suite_interstitial_3_cap

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_f_pointer, c_ptr, c_int32_t
    use            :: ccpp_types,                                      &
                      only: ccpp_t, CCPP_GENERIC_KIND
    use            :: ccpp_fields,                                     &
                      only: ccpp_field_get
    use            :: ccpp_errors,                                     &
                      only: ccpp_error, ccpp_debug
    use            :: GFS_suite_interstitial_3, &
                      only: GFS_suite_interstitial_3_run,GFS_suite_interstitial_3_init,GFS_suite_interstitial_3_finalize
    ! Other modules required, e.g. type definitions
    use machine, only: kind_phys

    implicit none

    private
    public :: GFS_suite_interstitial_3_run_cap,GFS_suite_interstitial_3_init_cap,GFS_suite_interstitial_3_finalize_cap

    contains


    function GFS_suite_interstitial_3_run_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind
        integer, pointer :: im
        integer, pointer :: levs
        integer, pointer :: nn
        logical, pointer :: cscnv
        logical, pointer :: satmedmf
        logical, pointer :: trans_trac
        logical, pointer :: do_shoc
        logical, pointer :: ltaerosol
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
        real(kind_phys), pointer :: xlat(:)
        real(kind_phys), pointer :: gq0(:,:,:)
        integer, pointer :: imp_physics
        integer, pointer :: imp_physics_mg
        integer, pointer :: imp_physics_zhao_carr
        integer, pointer :: imp_physics_zhao_carr_pdf
        integer, pointer :: imp_physics_gfdl
        integer, pointer :: imp_physics_thompson
        integer, pointer :: imp_physics_wsm6
        real(kind_phys), pointer :: prsi(:,:)
        real(kind_phys), pointer :: prsl(:,:)
        real(kind_phys), pointer :: prslk(:,:)
        real(kind_phys), pointer :: rhcbot
        real(kind_phys), pointer :: rhcpbl
        real(kind_phys), pointer :: rhctop
        real(kind_phys), pointer :: rhcmax
        integer, pointer :: islmsk(:)
        real(kind_phys), pointer :: work1(:)
        real(kind_phys), pointer :: work2(:)
        integer, pointer :: kpbl(:)
        real(kind_phys), pointer :: clw(:,:,:)
        real(kind_phys), pointer :: rhc(:,:)
        real(kind_phys), pointer :: save_qc(:,:)
        real(kind_phys), pointer :: save_qi(:,:)

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
        

        call ccpp_field_get(cdata, 'flag_for_Chikira_Sugiyama_deep_convection', cscnv, ierr=ierr, kind=ckind, index=269)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_Chikira_Sugiyama_deep_convection from CCPP data structure')
            return
        end if
        if (kind(cscnv).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_Chikira_Sugiyama_deep_convection')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'flag_for_scale_aware_TKE_moist_EDMF_PBL', satmedmf, ierr=ierr, kind=ckind, index=311)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_scale_aware_TKE_moist_EDMF_PBL from CCPP data structure')
            return
        end if
        if (kind(satmedmf).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_scale_aware_TKE_moist_EDMF_PBL')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'flag_for_convective_transport_of_tracers', trans_trac, ierr=ierr, kind=ckind, index=275)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_convective_transport_of_tracers from CCPP data structure')
            return
        end if
        if (kind(trans_trac).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_convective_transport_of_tracers')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'flag_for_shoc', do_shoc, ierr=ierr, kind=ckind, index=313)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_shoc from CCPP data structure')
            return
        end if
        if (kind(do_shoc).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_shoc')
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
        

        call ccpp_field_get(cdata, 'latitude', xlat, ierr=ierr, dims=cdims, kind=ckind, index=460)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve latitude from CCPP data structure')
            return
        end if
        if (kind(xlat).ne.ckind) then
            call ccpp_error('Kind mismatch for variable latitude')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

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
        

        call ccpp_field_get(cdata, 'flag_for_morrison_gettelman_microphysics_scheme', imp_physics_mg, ierr=ierr, kind=ckind, index=297)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_morrison_gettelman_microphysics_scheme from CCPP data structure')
            return
        end if
        if (kind(imp_physics_mg).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_morrison_gettelman_microphysics_scheme')
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
        

        call ccpp_field_get(cdata, 'flag_for_wsm6_microphysics_scheme', imp_physics_wsm6, ierr=ierr, kind=ckind, index=326)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_wsm6_microphysics_scheme from CCPP data structure')
            return
        end if
        if (kind(imp_physics_wsm6).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_wsm6_microphysics_scheme')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'air_pressure_at_interface', prsi, ierr=ierr, dims=cdims, kind=ckind, index=45)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve air_pressure_at_interface from CCPP data structure')
            return
        end if
        if (kind(prsi).ne.ckind) then
            call ccpp_error('Kind mismatch for variable air_pressure_at_interface')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'air_pressure', prsl, ierr=ierr, dims=cdims, kind=ckind, index=44)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve air_pressure from CCPP data structure')
            return
        end if
        if (kind(prsl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable air_pressure')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'dimensionless_exner_function_at_model_layers', prslk, ierr=ierr, dims=cdims, kind=ckind, index=222)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve dimensionless_exner_function_at_model_layers from CCPP data structure')
            return
        end if
        if (kind(prslk).ne.ckind) then
            call ccpp_error('Kind mismatch for variable dimensionless_exner_function_at_model_layers')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'critical_relative_humidity_at_surface', rhcbot, ierr=ierr, kind=ckind, index=144)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve critical_relative_humidity_at_surface from CCPP data structure')
            return
        end if
        if (kind(rhcbot).ne.ckind) then
            call ccpp_error('Kind mismatch for variable critical_relative_humidity_at_surface')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'critical_relative_humidity_at_PBL_top', rhcpbl, ierr=ierr, kind=ckind, index=142)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve critical_relative_humidity_at_PBL_top from CCPP data structure')
            return
        end if
        if (kind(rhcpbl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable critical_relative_humidity_at_PBL_top')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'critical_relative_humidity_at_top_of_atmosphere', rhctop, ierr=ierr, kind=ckind, index=145)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve critical_relative_humidity_at_top_of_atmosphere from CCPP data structure')
            return
        end if
        if (kind(rhctop).ne.ckind) then
            call ccpp_error('Kind mismatch for variable critical_relative_humidity_at_top_of_atmosphere')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'maximum_critical_relative_humidity', rhcmax, ierr=ierr, kind=ckind, index=504)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve maximum_critical_relative_humidity from CCPP data structure')
            return
        end if
        if (kind(rhcmax).ne.ckind) then
            call ccpp_error('Kind mismatch for variable maximum_critical_relative_humidity')
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
        

        call ccpp_field_get(cdata, 'vertical_index_at_top_of_atmosphere_boundary_layer', kpbl, ierr=ierr, dims=cdims, kind=ckind, index=822)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve vertical_index_at_top_of_atmosphere_boundary_layer from CCPP data structure')
            return
        end if
        if (kind(kpbl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable vertical_index_at_top_of_atmosphere_boundary_layer')
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
        

        call ccpp_field_get(cdata, 'critical_relative_humidity', rhc, ierr=ierr, dims=cdims, kind=ckind, index=141)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve critical_relative_humidity from CCPP data structure')
            return
        end if
        if (kind(rhc).ne.ckind) then
            call ccpp_error('Kind mismatch for variable critical_relative_humidity')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

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
        

        call GFS_suite_interstitial_3_run(im=im,levs=levs,nn=nn,cscnv=cscnv,satmedmf=satmedmf,trans_trac=trans_trac, &
                  do_shoc=do_shoc,ltaerosol=ltaerosol,ntrac=ntrac,ntcw=ntcw,ntiw=ntiw,ntclamt=ntclamt, &
                  ntrw=ntrw,ntsw=ntsw,ntrnc=ntrnc,ntsnc=ntsnc,ntgl=ntgl,ntgnc=ntgnc,xlat=xlat, &
                  gq0=gq0,imp_physics=imp_physics,imp_physics_mg=imp_physics_mg,imp_physics_zhao_carr=imp_physics_zhao_carr, &
                  imp_physics_zhao_carr_pdf=imp_physics_zhao_carr_pdf,imp_physics_gfdl=imp_physics_gfdl, &
                  imp_physics_thompson=imp_physics_thompson,imp_physics_wsm6=imp_physics_wsm6, &
                  prsi=prsi,prsl=prsl,prslk=prslk,rhcbot=rhcbot,rhcpbl=rhcpbl,rhctop=rhctop, &
                  rhcmax=rhcmax,islmsk=islmsk,work1=work1,work2=work2,kpbl=kpbl,clw=clw,rhc=rhc, &
                  save_qc=save_qc,save_qi=save_qi,errmsg=cdata%errmsg,errflg=cdata%errflg)
        ierr=cdata%errflg

    end function GFS_suite_interstitial_3_run_cap

    function GFS_suite_interstitial_3_init_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call GFS_suite_interstitial_3_init()
        

    end function GFS_suite_interstitial_3_init_cap

    function GFS_suite_interstitial_3_finalize_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call GFS_suite_interstitial_3_finalize()
        

    end function GFS_suite_interstitial_3_finalize_cap
end module GFS_suite_interstitial_3_cap


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
!! @brief Auto-generated cap module for the GFS_PBL_generic_pre scheme
!!
!
module GFS_PBL_generic_pre_cap

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_f_pointer, c_ptr, c_int32_t
    use            :: ccpp_types,                                      &
                      only: ccpp_t, CCPP_GENERIC_KIND
    use            :: ccpp_fields,                                     &
                      only: ccpp_field_get
    use            :: ccpp_errors,                                     &
                      only: ccpp_error, ccpp_debug
    use            :: GFS_PBL_generic_pre, &
                      only: GFS_PBL_generic_pre_run,GFS_PBL_generic_pre_finalize,GFS_PBL_generic_pre_init
    ! Other modules required, e.g. type definitions
    use machine, only: kind_phys

    implicit none

    private
    public :: GFS_PBL_generic_pre_run_cap,GFS_PBL_generic_pre_finalize_cap,GFS_PBL_generic_pre_init_cap

    contains


    function GFS_PBL_generic_pre_run_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind
        integer, pointer :: im
        integer, pointer :: levs
        integer, pointer :: nvdiff
        integer, pointer :: ntrac
        integer, pointer :: ntqv
        integer, pointer :: ntcw
        integer, pointer :: ntiw
        integer, pointer :: ntrw
        integer, pointer :: ntsw
        integer, pointer :: ntlnc
        integer, pointer :: ntinc
        integer, pointer :: ntwa
        integer, pointer :: ntia
        integer, pointer :: ntgl
        integer, pointer :: ntoz
        integer, pointer :: ntke
        integer, pointer :: ntkev
        integer, pointer :: imp_physics
        integer, pointer :: imp_physics_gfdl
        integer, pointer :: imp_physics_thompson
        integer, pointer :: imp_physics_wsm6
        logical, pointer :: ltaerosol
        logical, pointer :: satmedmf
        real(kind_phys), pointer :: qgrs(:,:,:)
        real(kind_phys), pointer :: vdftra(:,:,:)

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
        

        call ccpp_field_get(cdata, 'number_of_vertical_diffusion_tracers', nvdiff, ierr=ierr, kind=ckind, index=590)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve number_of_vertical_diffusion_tracers from CCPP data structure')
            return
        end if
        if (kind(nvdiff).ne.ckind) then
            call ccpp_error('Kind mismatch for variable number_of_vertical_diffusion_tracers')
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
        

        call ccpp_field_get(cdata, 'index_for_water_vapor', ntqv, ierr=ierr, kind=ckind, index=395)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve index_for_water_vapor from CCPP data structure')
            return
        end if
        if (kind(ntqv).ne.ckind) then
            call ccpp_error('Kind mismatch for variable index_for_water_vapor')
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
        

        call ccpp_field_get(cdata, 'index_for_water_friendly_aerosols', ntwa, ierr=ierr, kind=ckind, index=394)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve index_for_water_friendly_aerosols from CCPP data structure')
            return
        end if
        if (kind(ntwa).ne.ckind) then
            call ccpp_error('Kind mismatch for variable index_for_water_friendly_aerosols')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'index_for_ice_friendly_aerosols', ntia, ierr=ierr, kind=ckind, index=383)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve index_for_ice_friendly_aerosols from CCPP data structure')
            return
        end if
        if (kind(ntia).ne.ckind) then
            call ccpp_error('Kind mismatch for variable index_for_ice_friendly_aerosols')
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
        

        call ccpp_field_get(cdata, 'index_for_ozone', ntoz, ierr=ierr, kind=ckind, index=386)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve index_for_ozone from CCPP data structure')
            return
        end if
        if (kind(ntoz).ne.ckind) then
            call ccpp_error('Kind mismatch for variable index_for_ozone')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'index_for_turbulent_kinetic_energy', ntke, ierr=ierr, kind=ckind, index=391)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve index_for_turbulent_kinetic_energy from CCPP data structure')
            return
        end if
        if (kind(ntke).ne.ckind) then
            call ccpp_error('Kind mismatch for variable index_for_turbulent_kinetic_energy')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'index_for_turbulent_kinetic_energy_vertical_diffusion_tracer', ntkev, ierr=ierr, kind=ckind, index=393)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve index_for_turbulent_kinetic_energy_vertical_diffusion_tracer from CCPP data structure')
            return
        end if
        if (kind(ntkev).ne.ckind) then
            call ccpp_error('Kind mismatch for variable index_for_turbulent_kinetic_energy_vertical_diffusion_tracer')
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
        

        call ccpp_field_get(cdata, 'tracer_concentration', qgrs, ierr=ierr, dims=cdims, kind=ckind, index=801)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve tracer_concentration from CCPP data structure')
            return
        end if
        if (kind(qgrs).ne.ckind) then
            call ccpp_error('Kind mismatch for variable tracer_concentration')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'vertically_diffused_tracer_concentration', vdftra, ierr=ierr, dims=cdims, kind=ckind, index=831)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve vertically_diffused_tracer_concentration from CCPP data structure')
            return
        end if
        if (kind(vdftra).ne.ckind) then
            call ccpp_error('Kind mismatch for variable vertically_diffused_tracer_concentration')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call GFS_PBL_generic_pre_run(im=im,levs=levs,nvdiff=nvdiff,ntrac=ntrac,ntqv=ntqv,ntcw=ntcw,ntiw=ntiw, &
                  ntrw=ntrw,ntsw=ntsw,ntlnc=ntlnc,ntinc=ntinc,ntwa=ntwa,ntia=ntia,ntgl=ntgl, &
                  ntoz=ntoz,ntke=ntke,ntkev=ntkev,imp_physics=imp_physics,imp_physics_gfdl=imp_physics_gfdl, &
                  imp_physics_thompson=imp_physics_thompson,imp_physics_wsm6=imp_physics_wsm6, &
                  ltaerosol=ltaerosol,satmedmf=satmedmf,qgrs=qgrs,vdftra=vdftra,errmsg=cdata%errmsg, &
                  errflg=cdata%errflg)
        ierr=cdata%errflg

    end function GFS_PBL_generic_pre_run_cap

    function GFS_PBL_generic_pre_finalize_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call GFS_PBL_generic_pre_finalize()
        

    end function GFS_PBL_generic_pre_finalize_cap

    function GFS_PBL_generic_pre_init_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call GFS_PBL_generic_pre_init()
        

    end function GFS_PBL_generic_pre_init_cap
end module GFS_PBL_generic_pre_cap

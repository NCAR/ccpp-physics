
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
!! @brief Auto-generated cap module for the GFS_rrtmg_setup scheme
!!
!
module GFS_rrtmg_setup_cap

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_f_pointer, c_ptr, c_int32_t
    use            :: ccpp_types,                                      &
                      only: ccpp_t, CCPP_GENERIC_KIND
    use            :: ccpp_fields,                                     &
                      only: ccpp_field_get
    use            :: ccpp_errors,                                     &
                      only: ccpp_error, ccpp_debug
    use            :: GFS_rrtmg_setup, &
                      only: GFS_rrtmg_setup_run,GFS_rrtmg_setup_finalize,GFS_rrtmg_setup_init
    ! Other modules required, e.g. type definitions
    use machine, only: kind_phys

    implicit none

    private
    public :: GFS_rrtmg_setup_run_cap,GFS_rrtmg_setup_finalize_cap,GFS_rrtmg_setup_init_cap

    contains


    function GFS_rrtmg_setup_run_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind
        integer, pointer :: idate(:)
        integer, pointer :: jdate(:)
        real(kind_phys), pointer :: deltsw
        real(kind_phys), pointer :: deltim
        logical, pointer :: lsswr
        integer, pointer :: me
        real(kind_phys), pointer :: slag
        real(kind_phys), pointer :: sdec
        real(kind_phys), pointer :: cdec
        real(kind_phys), pointer :: solcon

        ierr = 0

        call c_f_pointer(ptr, cdata)


        call ccpp_field_get(cdata, 'date_and_time_at_model_initialization', idate, ierr=ierr, dims=cdims, kind=ckind, index=206)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve date_and_time_at_model_initialization from CCPP data structure')
            return
        end if
        if (kind(idate).ne.ckind) then
            call ccpp_error('Kind mismatch for variable date_and_time_at_model_initialization')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'forecast_date_and_time', jdate, ierr=ierr, dims=cdims, kind=ckind, index=337)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve forecast_date_and_time from CCPP data structure')
            return
        end if
        if (kind(jdate).ne.ckind) then
            call ccpp_error('Kind mismatch for variable forecast_date_and_time')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'frequency_for_shortwave_radiation', deltsw, ierr=ierr, kind=ckind, index=345)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve frequency_for_shortwave_radiation from CCPP data structure')
            return
        end if
        if (kind(deltsw).ne.ckind) then
            call ccpp_error('Kind mismatch for variable frequency_for_shortwave_radiation')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'time_step_for_dynamics', deltim, ierr=ierr, kind=ckind, index=792)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve time_step_for_dynamics from CCPP data structure')
            return
        end if
        if (kind(deltim).ne.ckind) then
            call ccpp_error('Kind mismatch for variable time_step_for_dynamics')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'flag_to_calc_sw', lsswr, ierr=ierr, kind=ckind, index=336)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_to_calc_sw from CCPP data structure')
            return
        end if
        if (kind(lsswr).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_to_calc_sw')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'mpi_rank', me, ierr=ierr, kind=ckind, index=558)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve mpi_rank from CCPP data structure')
            return
        end if
        if (kind(me).ne.ckind) then
            call ccpp_error('Kind mismatch for variable mpi_rank')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'equation_of_time', slag, ierr=ierr, kind=ckind, index=256)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve equation_of_time from CCPP data structure')
            return
        end if
        if (kind(slag).ne.ckind) then
            call ccpp_error('Kind mismatch for variable equation_of_time')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'sine_of_solar_declination_angle', sdec, ierr=ierr, kind=ckind, index=646)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve sine_of_solar_declination_angle from CCPP data structure')
            return
        end if
        if (kind(sdec).ne.ckind) then
            call ccpp_error('Kind mismatch for variable sine_of_solar_declination_angle')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'cosine_of_solar_declination_angle', cdec, ierr=ierr, kind=ckind, index=135)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cosine_of_solar_declination_angle from CCPP data structure')
            return
        end if
        if (kind(cdec).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cosine_of_solar_declination_angle')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'solar_constant', solcon, ierr=ierr, kind=ckind, index=663)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve solar_constant from CCPP data structure')
            return
        end if
        if (kind(solcon).ne.ckind) then
            call ccpp_error('Kind mismatch for variable solar_constant')
            ierr = 1
            return
        end if
#endif
        

        call GFS_rrtmg_setup_run(idate=idate,jdate=jdate,deltsw=deltsw,deltim=deltim,lsswr=lsswr,me=me,slag=slag, &
                  sdec=sdec,cdec=cdec,solcon=solcon,errmsg=cdata%errmsg,errflg=cdata%errflg)
        ierr=cdata%errflg

    end function GFS_rrtmg_setup_run_cap

    function GFS_rrtmg_setup_finalize_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call GFS_rrtmg_setup_finalize(errmsg=cdata%errmsg,errflg=cdata%errflg)
        ierr=cdata%errflg

    end function GFS_rrtmg_setup_finalize_cap

    function GFS_rrtmg_setup_init_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind
        real(kind_phys), pointer :: si(:)
        integer, pointer :: levr
        integer, pointer :: ictm
        integer, pointer :: isol
        integer, pointer :: ico2
        integer, pointer :: iaer
        integer, pointer :: ialb
        integer, pointer :: iems
        integer, pointer :: ntcw
        integer, pointer :: num_p2d
        integer, pointer :: num_p3d
        integer, pointer :: npdf3d
        integer, pointer :: ntoz
        integer, pointer :: iovr_sw
        integer, pointer :: iovr_lw
        integer, pointer :: isubc_sw
        integer, pointer :: isubc_lw
        integer, pointer :: icliq_sw
        logical, pointer :: crick_proof
        logical, pointer :: ccnorm
        integer, pointer :: imp_physics
        logical, pointer :: norad_precip
        integer, pointer :: idate(:)
        integer, pointer :: iflip
        integer, pointer :: im
        real(kind_phys), pointer :: faerlw(:,:,:,:)
        real(kind_phys), pointer :: faersw(:,:,:,:)
        real(kind_phys), pointer :: aerodp(:,:)
        integer, pointer :: me

        ierr = 0

        call c_f_pointer(ptr, cdata)


        call ccpp_field_get(cdata, 'vertical_sigma_coordinate_for_radiation_initialization', si, ierr=ierr, dims=cdims, kind=ckind, index=827)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve vertical_sigma_coordinate_for_radiation_initialization from CCPP data structure')
            return
        end if
        if (kind(si).ne.ckind) then
            call ccpp_error('Kind mismatch for variable vertical_sigma_coordinate_for_radiation_initialization')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'number_of_vertical_layers_for_radiation_calculations', levr, ierr=ierr, kind=ckind, index=591)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve number_of_vertical_layers_for_radiation_calculations from CCPP data structure')
            return
        end if
        if (kind(levr).ne.ckind) then
            call ccpp_error('Kind mismatch for variable number_of_vertical_layers_for_radiation_calculations')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'flag_for_initial_time-date_control', ictm, ierr=ierr, kind=ckind, index=286)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_initial_time-date_control from CCPP data structure')
            return
        end if
        if (kind(ictm).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_initial_time-date_control')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'flag_for_solar_constant', isol, ierr=ierr, kind=ckind, index=315)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_solar_constant from CCPP data structure')
            return
        end if
        if (kind(isol).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_solar_constant')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'flag_for_using_prescribed_global_mean_co2_value', ico2, ierr=ierr, kind=ckind, index=324)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_using_prescribed_global_mean_co2_value from CCPP data structure')
            return
        end if
        if (kind(ico2).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_using_prescribed_global_mean_co2_value')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'flag_for_default_aerosol_effect_in_shortwave_radiation', iaer, ierr=ierr, kind=ckind, index=276)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_default_aerosol_effect_in_shortwave_radiation from CCPP data structure')
            return
        end if
        if (kind(iaer).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_default_aerosol_effect_in_shortwave_radiation')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'flag_for_using_climatology_albedo', ialb, ierr=ierr, kind=ckind, index=323)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_using_climatology_albedo from CCPP data structure')
            return
        end if
        if (kind(ialb).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_using_climatology_albedo')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'flag_for_surface_emissivity_control', iems, ierr=ierr, kind=ckind, index=320)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_surface_emissivity_control from CCPP data structure')
            return
        end if
        if (kind(iems).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_surface_emissivity_control')
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
        

        call ccpp_field_get(cdata, 'array_dimension_of_2d_arrays_for_microphysics', num_p2d, ierr=ierr, kind=ckind, index=62)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve array_dimension_of_2d_arrays_for_microphysics from CCPP data structure')
            return
        end if
        if (kind(num_p2d).ne.ckind) then
            call ccpp_error('Kind mismatch for variable array_dimension_of_2d_arrays_for_microphysics')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'array_dimension_of_3d_arrays_for_microphysics', num_p3d, ierr=ierr, kind=ckind, index=63)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve array_dimension_of_3d_arrays_for_microphysics from CCPP data structure')
            return
        end if
        if (kind(num_p3d).ne.ckind) then
            call ccpp_error('Kind mismatch for variable array_dimension_of_3d_arrays_for_microphysics')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'number_of_3d_arrays_associated_with_pdf-based_clouds', npdf3d, ierr=ierr, kind=ckind, index=571)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve number_of_3d_arrays_associated_with_pdf-based_clouds from CCPP data structure')
            return
        end if
        if (kind(npdf3d).ne.ckind) then
            call ccpp_error('Kind mismatch for variable number_of_3d_arrays_associated_with_pdf-based_clouds')
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
        

        call ccpp_field_get(cdata, 'flag_for_max-random_overlap_clouds_for_shortwave_radiation', iovr_sw, ierr=ierr, kind=ckind, index=293)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_max-random_overlap_clouds_for_shortwave_radiation from CCPP data structure')
            return
        end if
        if (kind(iovr_sw).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_max-random_overlap_clouds_for_shortwave_radiation')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'flag_for_max-random_overlap_clouds_for_longwave_radiation', iovr_lw, ierr=ierr, kind=ckind, index=292)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_max-random_overlap_clouds_for_longwave_radiation from CCPP data structure')
            return
        end if
        if (kind(iovr_lw).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_max-random_overlap_clouds_for_longwave_radiation')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'flag_for_sw_clouds_without_sub-grid_approximation', isubc_sw, ierr=ierr, kind=ckind, index=321)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_sw_clouds_without_sub-grid_approximation from CCPP data structure')
            return
        end if
        if (kind(isubc_sw).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_sw_clouds_without_sub-grid_approximation')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'flag_for_lw_clouds_without_sub-grid_approximation', isubc_lw, ierr=ierr, kind=ckind, index=289)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_lw_clouds_without_sub-grid_approximation from CCPP data structure')
            return
        end if
        if (kind(isubc_lw).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_lw_clouds_without_sub-grid_approximation')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'flag_for_optical_property_for_liquid_clouds_for_shortwave_radiation', icliq_sw, ierr=ierr, kind=ckind, index=301)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_optical_property_for_liquid_clouds_for_shortwave_radiation from CCPP data structure')
            return
        end if
        if (kind(icliq_sw).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_optical_property_for_liquid_clouds_for_shortwave_radiation')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'flag_for_CRICK-proof_cloud_water', crick_proof, ierr=ierr, kind=ckind, index=268)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_CRICK-proof_cloud_water from CCPP data structure')
            return
        end if
        if (kind(crick_proof).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_CRICK-proof_cloud_water')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'flag_for_cloud_condensate_normalized_by_cloud_cover', ccnorm, ierr=ierr, kind=ckind, index=274)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_cloud_condensate_normalized_by_cloud_cover from CCPP data structure')
            return
        end if
        if (kind(ccnorm).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_cloud_condensate_normalized_by_cloud_cover')
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
        

        call ccpp_field_get(cdata, 'flag_for_precipitation_effect_on_radiation', norad_precip, ierr=ierr, kind=ckind, index=303)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_precipitation_effect_on_radiation from CCPP data structure')
            return
        end if
        if (kind(norad_precip).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_precipitation_effect_on_radiation')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'date_and_time_at_model_initialization_reordered', idate, ierr=ierr, dims=cdims, kind=ckind, index=207)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve date_and_time_at_model_initialization_reordered from CCPP data structure')
            return
        end if
        if (kind(idate).ne.ckind) then
            call ccpp_error('Kind mismatch for variable date_and_time_at_model_initialization_reordered')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'flag_for_vertical_index_direction_control', iflip, ierr=ierr, kind=ckind, index=325)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_vertical_index_direction_control from CCPP data structure')
            return
        end if
        if (kind(iflip).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_vertical_index_direction_control')
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
        

        call ccpp_field_get(cdata, 'aerosol_optical_properties_for_longwave_bands_01-16', faerlw, ierr=ierr, dims=cdims, kind=ckind, index=40)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve aerosol_optical_properties_for_longwave_bands_01-16 from CCPP data structure')
            return
        end if
        if (kind(faerlw).ne.ckind) then
            call ccpp_error('Kind mismatch for variable aerosol_optical_properties_for_longwave_bands_01-16')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'aerosol_optical_properties_for_shortwave_bands_01-16', faersw, ierr=ierr, dims=cdims, kind=ckind, index=41)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve aerosol_optical_properties_for_shortwave_bands_01-16 from CCPP data structure')
            return
        end if
        if (kind(faersw).ne.ckind) then
            call ccpp_error('Kind mismatch for variable aerosol_optical_properties_for_shortwave_bands_01-16')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'atmosphere_optical_thickness_due_to_ambient_aerosol_particles', aerodp, ierr=ierr, dims=cdims, kind=ckind, index=75)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve atmosphere_optical_thickness_due_to_ambient_aerosol_particles from CCPP data structure')
            return
        end if
        if (kind(aerodp).ne.ckind) then
            call ccpp_error('Kind mismatch for variable atmosphere_optical_thickness_due_to_ambient_aerosol_particles')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'mpi_rank', me, ierr=ierr, kind=ckind, index=558)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve mpi_rank from CCPP data structure')
            return
        end if
        if (kind(me).ne.ckind) then
            call ccpp_error('Kind mismatch for variable mpi_rank')
            ierr = 1
            return
        end if
#endif
        

        call GFS_rrtmg_setup_init(si=si,levr=levr,ictm=ictm,isol=isol,ico2=ico2,iaer=iaer,ialb=ialb,iems=iems, &
                  ntcw=ntcw,num_p2d=num_p2d,num_p3d=num_p3d,npdf3d=npdf3d,ntoz=ntoz,iovr_sw=iovr_sw, &
                  iovr_lw=iovr_lw,isubc_sw=isubc_sw,isubc_lw=isubc_lw,icliq_sw=icliq_sw,crick_proof=crick_proof, &
                  ccnorm=ccnorm,imp_physics=imp_physics,norad_precip=norad_precip,idate=idate, &
                  iflip=iflip,im=im,faerlw=faerlw,faersw=faersw,aerodp=aerodp,me=me,errmsg=cdata%errmsg, &
                  errflg=cdata%errflg)
        ierr=cdata%errflg

    end function GFS_rrtmg_setup_init_cap
end module GFS_rrtmg_setup_cap

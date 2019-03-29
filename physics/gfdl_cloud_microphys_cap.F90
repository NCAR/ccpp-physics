
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
!! @brief Auto-generated cap module for the gfdl_cloud_microphys scheme
!!
!
module gfdl_cloud_microphys_cap

    use, intrinsic :: iso_c_binding,                                   &
                      only: c_f_pointer, c_ptr, c_int32_t
    use            :: ccpp_types,                                      &
                      only: ccpp_t, CCPP_GENERIC_KIND
    use            :: ccpp_fields,                                     &
                      only: ccpp_field_get
    use            :: ccpp_errors,                                     &
                      only: ccpp_error, ccpp_debug
    use            :: gfdl_cloud_microphys, &
                      only: gfdl_cloud_microphys_init,gfdl_cloud_microphys_run,gfdl_cloud_microphys_finalize
    ! Other modules required, e.g. type definitions
    use machine, only: kind_phys

    implicit none

    private
    public :: gfdl_cloud_microphys_init_cap,gfdl_cloud_microphys_run_cap,gfdl_cloud_microphys_finalize_cap

    contains


    function gfdl_cloud_microphys_init_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind
        integer, pointer :: me
        integer, pointer :: master
        integer, pointer :: nlunit
        character(len=256), pointer :: input_nml_file(:)
        integer, pointer :: logunit
        character(len=64), pointer :: fn_nml
        integer, pointer :: imp_physics
        integer, pointer :: imp_physics_gfdl
        logical, pointer :: do_shoc

        ierr = 0

        call c_f_pointer(ptr, cdata)


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
        

        call ccpp_field_get(cdata, 'mpi_root', master, ierr=ierr, kind=ckind, index=559)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve mpi_root from CCPP data structure')
            return
        end if
        if (kind(master).ne.ckind) then
            call ccpp_error('Kind mismatch for variable mpi_root')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'iounit_namelist', nlunit, ierr=ierr, kind=ckind, index=451)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve iounit_namelist from CCPP data structure')
            return
        end if
        if (kind(nlunit).ne.ckind) then
            call ccpp_error('Kind mismatch for variable iounit_namelist')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'namelist_filename_for_internal_file_reads', input_nml_file, ierr=ierr, dims=cdims, kind=ckind, index=564)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve namelist_filename_for_internal_file_reads from CCPP data structure')
            return
        end if
        if (kind(input_nml_file).ne.ckind) then
            call ccpp_error('Kind mismatch for variable namelist_filename_for_internal_file_reads')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'iounit_log', logunit, ierr=ierr, kind=ckind, index=450)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve iounit_log from CCPP data structure')
            return
        end if
        if (kind(logunit).ne.ckind) then
            call ccpp_error('Kind mismatch for variable iounit_log')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'namelist_filename', fn_nml, ierr=ierr, kind=ckind, index=563)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve namelist_filename from CCPP data structure')
            return
        end if
        if (kind(fn_nml).ne.ckind) then
            call ccpp_error('Kind mismatch for variable namelist_filename')
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
        

        call gfdl_cloud_microphys_init(me=me,master=master,nlunit=nlunit,input_nml_file=input_nml_file,logunit=logunit, &
                  fn_nml=fn_nml,imp_physics=imp_physics,imp_physics_gfdl=imp_physics_gfdl, &
                  do_shoc=do_shoc,errmsg=cdata%errmsg,errflg=cdata%errflg)
        ierr=cdata%errflg

    end function gfdl_cloud_microphys_init_cap

    function gfdl_cloud_microphys_run_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind
        integer, pointer :: levs
        integer, pointer :: im
        real(kind_phys), pointer :: con_g
        real(kind_phys), pointer :: con_fvirt
        real(kind_phys), pointer :: con_rd
        real(kind_phys), pointer :: frland(:)
        real(kind_phys), pointer :: garea(:)
        real(kind_phys), pointer :: gq0(:,:)
        real(kind_phys), pointer :: gq0_ntcw(:,:)
        real(kind_phys), pointer :: gq0_ntrw(:,:)
        real(kind_phys), pointer :: gq0_ntiw(:,:)
        real(kind_phys), pointer :: gq0_ntsw(:,:)
        real(kind_phys), pointer :: gq0_ntgl(:,:)
        real(kind_phys), pointer :: gq0_ntclamt(:,:)
        real(kind_phys), pointer :: gt0(:,:)
        real(kind_phys), pointer :: gu0(:,:)
        real(kind_phys), pointer :: gv0(:,:)
        real(kind_phys), pointer :: vvl(:,:)
        real(kind_phys), pointer :: prsl(:,:)
        real(kind_phys), pointer :: phii(:,:)
        real(kind_phys), pointer :: del(:,:)
        real(kind_phys), pointer :: rain0(:)
        real(kind_phys), pointer :: ice0(:)
        real(kind_phys), pointer :: snow0(:)
        real(kind_phys), pointer :: graupel0(:)
        real(kind_phys), pointer :: prcp0(:)
        real(kind_phys), pointer :: sr(:)
        real(kind_phys), pointer :: dtp
        logical, pointer :: hydrostatic
        logical, pointer :: phys_hydrostatic

        ierr = 0

        call c_f_pointer(ptr, cdata)


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
        

        call ccpp_field_get(cdata, 'gravitational_acceleration', con_g, ierr=ierr, kind=ckind, index=355)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve gravitational_acceleration from CCPP data structure')
            return
        end if
        if (kind(con_g).ne.ckind) then
            call ccpp_error('Kind mismatch for variable gravitational_acceleration')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'ratio_of_vapor_to_dry_air_gas_constants_minus_one', con_fvirt, ierr=ierr, kind=ckind, index=625)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve ratio_of_vapor_to_dry_air_gas_constants_minus_one from CCPP data structure')
            return
        end if
        if (kind(con_fvirt).ne.ckind) then
            call ccpp_error('Kind mismatch for variable ratio_of_vapor_to_dry_air_gas_constants_minus_one')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'gas_constant_dry_air', con_rd, ierr=ierr, kind=ckind, index=346)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve gas_constant_dry_air from CCPP data structure')
            return
        end if
        if (kind(con_rd).ne.ckind) then
            call ccpp_error('Kind mismatch for variable gas_constant_dry_air')
            ierr = 1
            return
        end if
#endif
        

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
        

        call ccpp_field_get(cdata, 'cell_area', garea, ierr=ierr, dims=cdims, kind=ckind, index=83)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cell_area from CCPP data structure')
            return
        end if
        if (kind(garea).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cell_area')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'water_vapor_specific_humidity_updated_by_physics', gq0, ierr=ierr, dims=cdims, kind=ckind, index=860)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve water_vapor_specific_humidity_updated_by_physics from CCPP data structure')
            return
        end if
        if (kind(gq0).ne.ckind) then
            call ccpp_error('Kind mismatch for variable water_vapor_specific_humidity_updated_by_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cloud_condensed_water_mixing_ratio_updated_by_physics', gq0_ntcw, ierr=ierr, dims=cdims, kind=ckind, index=95)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cloud_condensed_water_mixing_ratio_updated_by_physics from CCPP data structure')
            return
        end if
        if (kind(gq0_ntcw).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cloud_condensed_water_mixing_ratio_updated_by_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'rain_water_mixing_ratio_updated_by_physics', gq0_ntrw, ierr=ierr, dims=cdims, kind=ckind, index=619)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve rain_water_mixing_ratio_updated_by_physics from CCPP data structure')
            return
        end if
        if (kind(gq0_ntrw).ne.ckind) then
            call ccpp_error('Kind mismatch for variable rain_water_mixing_ratio_updated_by_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'ice_water_mixing_ratio_updated_by_physics', gq0_ntiw, ierr=ierr, dims=cdims, kind=ckind, index=376)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve ice_water_mixing_ratio_updated_by_physics from CCPP data structure')
            return
        end if
        if (kind(gq0_ntiw).ne.ckind) then
            call ccpp_error('Kind mismatch for variable ice_water_mixing_ratio_updated_by_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'snow_water_mixing_ratio_updated_by_physics', gq0_ntsw, ierr=ierr, dims=cdims, kind=ckind, index=653)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve snow_water_mixing_ratio_updated_by_physics from CCPP data structure')
            return
        end if
        if (kind(gq0_ntsw).ne.ckind) then
            call ccpp_error('Kind mismatch for variable snow_water_mixing_ratio_updated_by_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'graupel_mixing_ratio_updated_by_physics', gq0_ntgl, ierr=ierr, dims=cdims, kind=ckind, index=352)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve graupel_mixing_ratio_updated_by_physics from CCPP data structure')
            return
        end if
        if (kind(gq0_ntgl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable graupel_mixing_ratio_updated_by_physics')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'cloud_fraction_updated_by_physics', gq0_ntclamt, ierr=ierr, dims=cdims, kind=ckind, index=100)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve cloud_fraction_updated_by_physics from CCPP data structure')
            return
        end if
        if (kind(gq0_ntclamt).ne.ckind) then
            call ccpp_error('Kind mismatch for variable cloud_fraction_updated_by_physics')
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
        

        call ccpp_field_get(cdata, 'omega', vvl, ierr=ierr, dims=cdims, kind=ckind, index=593)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve omega from CCPP data structure')
            return
        end if
        if (kind(vvl).ne.ckind) then
            call ccpp_error('Kind mismatch for variable omega')
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
        

        call ccpp_field_get(cdata, 'geopotential_at_interface', phii, ierr=ierr, dims=cdims, kind=ckind, index=349)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve geopotential_at_interface from CCPP data structure')
            return
        end if
        if (kind(phii).ne.ckind) then
            call ccpp_error('Kind mismatch for variable geopotential_at_interface')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'air_pressure_difference_between_midlayers', del, ierr=ierr, dims=cdims, kind=ckind, index=49)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve air_pressure_difference_between_midlayers from CCPP data structure')
            return
        end if
        if (kind(del).ne.ckind) then
            call ccpp_error('Kind mismatch for variable air_pressure_difference_between_midlayers')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'lwe_thickness_of_explicit_rain_amount', rain0, ierr=ierr, dims=cdims, kind=ckind, index=481)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve lwe_thickness_of_explicit_rain_amount from CCPP data structure')
            return
        end if
        if (kind(rain0).ne.ckind) then
            call ccpp_error('Kind mismatch for variable lwe_thickness_of_explicit_rain_amount')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'lwe_thickness_of_ice_amount', ice0, ierr=ierr, dims=cdims, kind=ckind, index=486)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve lwe_thickness_of_ice_amount from CCPP data structure')
            return
        end if
        if (kind(ice0).ne.ckind) then
            call ccpp_error('Kind mismatch for variable lwe_thickness_of_ice_amount')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'lwe_thickness_of_snow_amount', snow0, ierr=ierr, dims=cdims, kind=ckind, index=492)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve lwe_thickness_of_snow_amount from CCPP data structure')
            return
        end if
        if (kind(snow0).ne.ckind) then
            call ccpp_error('Kind mismatch for variable lwe_thickness_of_snow_amount')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'lwe_thickness_of_graupel_amount', graupel0, ierr=ierr, dims=cdims, kind=ckind, index=483)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve lwe_thickness_of_graupel_amount from CCPP data structure')
            return
        end if
        if (kind(graupel0).ne.ckind) then
            call ccpp_error('Kind mismatch for variable lwe_thickness_of_graupel_amount')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'lwe_thickness_of_explicit_precipitation_amount', prcp0, ierr=ierr, dims=cdims, kind=ckind, index=480)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve lwe_thickness_of_explicit_precipitation_amount from CCPP data structure')
            return
        end if
        if (kind(prcp0).ne.ckind) then
            call ccpp_error('Kind mismatch for variable lwe_thickness_of_explicit_precipitation_amount')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

        call ccpp_field_get(cdata, 'ratio_of_snowfall_to_rainfall', sr, ierr=ierr, dims=cdims, kind=ckind, index=624)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve ratio_of_snowfall_to_rainfall from CCPP data structure')
            return
        end if
        if (kind(sr).ne.ckind) then
            call ccpp_error('Kind mismatch for variable ratio_of_snowfall_to_rainfall')
            ierr = 1
            return
        end if
#endif
        deallocate(cdims)
        

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
        

        call ccpp_field_get(cdata, 'flag_for_hydrostatic_solver', hydrostatic, ierr=ierr, kind=ckind, index=284)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_hydrostatic_solver from CCPP data structure')
            return
        end if
        if (kind(hydrostatic).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_hydrostatic_solver')
            ierr = 1
            return
        end if
#endif
        

        call ccpp_field_get(cdata, 'flag_for_hydrostatic_heating_from_physics', phys_hydrostatic, ierr=ierr, kind=ckind, index=283)
#ifdef DEBUG
        if (ierr /= 0) then
            call ccpp_error('Unable to retrieve flag_for_hydrostatic_heating_from_physics from CCPP data structure')
            return
        end if
        if (kind(phys_hydrostatic).ne.ckind) then
            call ccpp_error('Kind mismatch for variable flag_for_hydrostatic_heating_from_physics')
            ierr = 1
            return
        end if
#endif
        

        call gfdl_cloud_microphys_run(levs=levs,im=im,con_g=con_g,con_fvirt=con_fvirt,con_rd=con_rd,frland=frland, &
                  garea=garea,gq0=gq0,gq0_ntcw=gq0_ntcw,gq0_ntrw=gq0_ntrw,gq0_ntiw=gq0_ntiw, &
                  gq0_ntsw=gq0_ntsw,gq0_ntgl=gq0_ntgl,gq0_ntclamt=gq0_ntclamt,gt0=gt0,gu0=gu0, &
                  gv0=gv0,vvl=vvl,prsl=prsl,phii=phii,del=del,rain0=rain0,ice0=ice0,snow0=snow0, &
                  graupel0=graupel0,prcp0=prcp0,sr=sr,dtp=dtp,hydrostatic=hydrostatic,phys_hydrostatic=phys_hydrostatic, &
                  errmsg=cdata%errmsg,errflg=cdata%errflg)
        ierr=cdata%errflg

    end function gfdl_cloud_microphys_run_cap

    function gfdl_cloud_microphys_finalize_cap(ptr) bind(c) result(ierr)

        integer(c_int32_t)         :: ierr
        type(c_ptr), intent(inout) :: ptr

        type(ccpp_t), pointer           :: cdata
        type(c_ptr)                     :: cptr
        integer, allocatable            :: cdims(:)
        integer                         :: ckind


        ierr = 0

        call c_f_pointer(ptr, cdata)



        call gfdl_cloud_microphys_finalize(errmsg=cdata%errmsg,errflg=cdata%errflg)
        ierr=cdata%errflg

    end function gfdl_cloud_microphys_finalize_cap
end module gfdl_cloud_microphys_cap

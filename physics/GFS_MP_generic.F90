!> \file GFS_MP_generic.F90
!! This file contains the subroutines that calculate diagnotics variables
!! before/after calling any microphysics scheme:

!> This module contains the CCPP-compliant MP generic pre interstitial codes.
      module GFS_MP_generic_pre
      contains


!! \section arg_table_GFS_MP_generic_pre_init Argument Table
!!
      subroutine GFS_MP_generic_pre_init
      end subroutine GFS_MP_generic_pre_init


!> \section arg_table_GFS_MP_generic_pre_run Argument Table
!! | local_name     | standard_name                                          | long_name                                                                 | units       | rank |  type     |   kind    | intent | optional |
!! |----------------|--------------------------------------------------------|---------------------------------------------------------------------------|-------------|------|-----------|-----------|--------|----------|
!! | im             | horizontal_loop_extent                                 | horizontal loop extent                                                    | count       |    0 | integer   |           | in     | F        |
!! | levs           | vertical_dimension                                     | vertical layer dimension                                                  | count       |    0 | integer   |           | in     | F        |
!! | ldiag3d        | flag_diagnostics_3D                                    | logical flag for 3D diagnostics                                           | flag        |    0 | logical   |           | in     | F        |
!! | do_aw          | flag_for_Arakawa_Wu_adjustment                         | flag for Arakawa Wu scale-aware adjustment                                | flag        |    0 | logical   |           | in     | F        |
!! | ntcw           | index_for_liquid_cloud_condensate                      | tracer index for cloud condensate (or liquid water)                       | index       |    0 | integer   |           | in     | F        |
!! | nncl           | number_of_tracers_for_cloud_condensate                 | number of tracers for cloud condensate                                    | count       |    0 | integer   |           | in     | F        |
!! | ntrac          | number_of_tracers                                      | number of tracers                                                         | count       |    0 | integer   |           | in     | F        |
!! | gt0            | air_temperature_updated_by_physics                     | temperature updated by physics                                            | K           |    2 | real      | kind_phys | in     | F        |
!! | gq0            | tracer_concentration_updated_by_physics                | tracer concentration updated by physics                                   | kg kg-1     |    3 | real      | kind_phys | in     | F        |
!! | save_t         | air_temperature_save                                   | air temperature before entering a physics scheme                          | K           |    2 | real      | kind_phys | inout  | F        |
!! | save_q         | tracer_concentration_save                              | tracer concentration before entering a physics scheme                     | kg kg-1     |    3 | real      | kind_phys | inout  | F        |
!! | errmsg         | ccpp_error_message                                     | error message for error handling in CCPP                                  | none        |    0 | character | len=*     | out    | F        |
!! | errflg         | ccpp_error_flag                                        | error flag for error handling in CCPP                                     | flag        |    0 | integer   |           | out    | F        |
!!
      subroutine GFS_MP_generic_pre_run(im, levs, ldiag3d, do_aw, ntcw, nncl, ntrac, gt0, gq0, save_t, save_q, errmsg, errflg)
!
      use machine,               only: kind_phys

      implicit none
      integer,                                          intent(in) :: im, levs, ntcw, nncl, ntrac
      logical,                                          intent(in) :: ldiag3d, do_aw
      real(kind=kind_phys), dimension(im, levs),        intent(in) :: gt0
      real(kind=kind_phys), dimension(im, levs, ntrac), intent(in) :: gq0

      real(kind=kind_phys), dimension(im, levs),        intent(inout) :: save_t
      real(kind=kind_phys), dimension(im, levs, ntrac), intent(inout) :: save_q

      character(len=*), intent(out) :: errmsg
      integer, intent(out) :: errflg

      integer :: i, k, n

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      if (ldiag3d .or. do_aw) then
        do k=1,levs
          do i=1,im
            save_t(i,k) = gt0(i,k)
            save_q(1:im,:,1) = gq0(1:im,:,1)
          enddo
        enddo
        do n=ntcw,ntcw+nncl-1
          save_q(1:im,:,n) = gq0(1:im,:,n)
        enddo
      endif

      end subroutine GFS_MP_generic_pre_run

!> \section arg_table_GFS_MP_generic_pre_finalize Argument Table
!!
      subroutine GFS_MP_generic_pre_finalize
      end subroutine GFS_MP_generic_pre_finalize

      end module GFS_MP_generic_pre

!> This module contains the subroutine that calculates 
!! precipitation type and its post, which provides precipitation forcing
!! to LSM.
      module GFS_MP_generic_post
      contains

!! \section arg_table_GFS_MP_generic_post_init Argument Table
!!
      subroutine GFS_MP_generic_post_init
      end subroutine GFS_MP_generic_post_init

!>\defgroup gfs_calpreciptype GFS Precipitation Type Diagnostics Module
!! \brief If dominant precip type is requested (i.e., Zhao-Carr MP scheme), 4 more algorithms in calpreciptype()
!! will be called.  the tallies are then summed in calwxt_dominant(). For GFDL cloud MP scheme, determine convective 
!! rain/snow by surface temperature;  and determine explicit rain/snow by rain/snow coming out directly from MP.
!! 
!! \section arg_table_GFS_MP_generic_post_run Argument Table
!! | local_name       | standard_name                                                           | long_name                                                               | units       | rank |  type      |   kind    | intent | optional |
!! |------------------|-------------------------------------------------------------------------|-------------------------------------------------------------------------|-------------|------|------------|-----------|--------|----------|
!! | im               | horizontal_loop_extent                                                  | horizontal loop extent                                                  | count       |    0 | integer    |           | in     | F        |
!! | ix               | horizontal_dimension                                                    | horizontal dimension                                                    | count       |    0 | integer    |           | in     | F        |
!! | levs             | vertical_dimension                                                      | vertical layer dimension                                                | count       |    0 | integer    |           | in     | F        |
!! | kdt              | index_of_time_step                                                      | current time step index                                                 | index       |    0 | integer    |           | in     | F        |
!! | nrcm             | array_dimension_of_random_number                                        | second dimension of random number stream for RAS                        | count       |    0 | integer    |           | in     | F        |
!! | ncld             | number_of_hydrometeors                                                  | choice of cloud scheme / number of hydrometeors                         | count       |    0 | integer    |           | in     | F        |
!! | nncl             | number_of_tracers_for_cloud_condensate                                  | number of tracers for cloud condensate                                  | count       |    0 | integer    |           | in     | F        |
!! | ntcw             | index_for_liquid_cloud_condensate                                       | tracer index for cloud condensate (or liquid water)                     | index       |    0 | integer    |           | in     | F        |
!! | ntrac            | number_of_tracers                                                       | number of tracers                                                       | count       |    0 | integer    |           | in     | F        |
!! | imp_physics      | flag_for_microphysics_scheme                                            | choice of microphysics scheme                                           | flag        |    0 | integer    |           | in     | F        |
!! | imp_physics_gfdl | flag_for_gfdl_microphysics_scheme                                       | choice of GFDL microphysics scheme                                      | flag        |    0 | integer    |           | in     | F        |
!! | imp_physics_thompson | flag_for_thompson_microphysics_scheme                               | choice of Thompson microphysics scheme                                  | flag        |    0 | integer    |           | in     | F        |
!! | imp_physics_mg   | flag_for_morrison_gettelman_microphysics_scheme                         | choice of Morrison-Gettelman microphysics scheme                        | flag        |    0 | integer    |           | in     | F        |
!! | cal_pre          | flag_for_precipitation_type_algorithm                                   | flag controls precip type algorithm                                     | flag        |    0 | logical    |           | in     | F        |
!! | lssav            | flag_diagnostics                                                        | logical flag for storing diagnostics                                    | flag        |    0 | logical    |           | in     | F        |
!! | ldiag3d          | flag_diagnostics_3D                                                     | flag for 3d diagnostic fields                                           | flag        |    0 | logical    |           | in     | F        |
!! | cplflx           | flag_for_flux_coupling                                                  | flag controlling cplflx collection (default off)                        | flag        |    0 | logical    |           | in     | F        |
!! | cplchm           | flag_for_chemistry_coupling                                             | flag controlling cplchm collection (default off)                        | flag        |    0 | logical    |           | in     | F        |
!! | con_g            | gravitational_acceleration                                              | gravitational acceleration                                              | m s-2       |    0 | real       | kind_phys | in     | F        |
!! | dtf              | time_step_for_dynamics                                                  | dynamics timestep                                                       | s           |    0 | real       | kind_phys | in     | F        |
!! | frain            | dynamics_to_physics_timestep_ratio                                      | ratio of dynamics timestep to physics timestep                          | none        |    0 | real       | kind_phys | in     | F        |
!! | rainc            | lwe_thickness_of_convective_precipitation_amount_on_dynamics_timestep   | convective rain at this time step                                       | m           |    1 | real       | kind_phys | in     | F        |
!! | rain1            | lwe_thickness_of_explicit_precipitation_amount                          | explicit rainfall amount on physics timestep                            | m           |    1 | real       | kind_phys | in     | F        |
!! | rann             | random_number_array                                                     | random number array (0-1)                                               | none        |    2 | real       | kind_phys | in     | F        |
!! | xlat             | latitude                                                                | latitude                                                                | radians     |    1 | real       | kind_phys | in     | F        |
!! | xlon             | longitude                                                               | longitude                                                               | radians     |    1 | real       | kind_phys | in     | F        |
!! | gt0              | air_temperature_updated_by_physics                                      | temperature updated by physics                                          | K           |    2 | real       | kind_phys | in     | F        |
!! | gq0              | tracer_concentration_updated_by_physics                                 | tracer concentration updated by physics                                 | kg kg-1     |    3 | real       | kind_phys | in     | F        |
!! | prsl             | air_pressure                                                            | layer mean pressure                                                     | Pa          |    2 | real       | kind_phys | in     | F        |
!! | prsi             | air_pressure_at_interface                                               | pressure at layer interface                                             | Pa          |    2 | real       | kind_phys | in     | F        |
!! | phii             | geopotential_at_interface                                               | geopotential at model layer interfaces                                  | m2 s-2      |    2 | real       | kind_phys | in     | F        |
!! | tsfc             | surface_skin_temperature                                                | surface skin temperature                                                | K           |    1 | real       | kind_phys | in     | F        |
!! | ice              | lwe_thickness_of_ice_amount_on_dynamics_timestep                        | ice fall at this time step                                              | m           |    1 | real       | kind_phys | inout  | F        |
!! | snow             | lwe_thickness_of_snow_amount_on_dynamics_timestep                       | snow fall at this time step                                             | m           |    1 | real       | kind_phys | inout  | F        |
!! | graupel          | lwe_thickness_of_graupel_amount_on_dynamics_timestep                    | graupel fall at this time step                                          | m           |    1 | real       | kind_phys | inout  | F        |
!! | save_t           | air_temperature_save                                                    | air temperature before entering a physics scheme                        | K           |    2 | real       | kind_phys | in     | F        |
!! | save_qv          | water_vapor_specific_humidity_save                                      | water vapor specific humidity before entering a physics scheme          | kg kg-1     |    2 | real       | kind_phys | in     | F        |
!! | rain0            | lwe_thickness_of_explicit_rain_amount                                   | explicit rain on physics timestep                                       | m           |    1 | real       | kind_phys | in     | F        |
!! | ice0             | lwe_thickness_of_ice_amount                                             | ice fall on physics timestep                                            | m           |    1 | real       | kind_phys | in     | F        |
!! | snow0            | lwe_thickness_of_snow_amount                                            | snow fall on physics timestep                                           | m           |    1 | real       | kind_phys | in     | F        |
!! | graupel0         | lwe_thickness_of_graupel_amount                                         | graupel fall on physics timestep                                        | m           |    1 | real       | kind_phys | in     | F        |
!! | del              | air_pressure_difference_between_midlayers                               | air pressure difference between midlayers                               | Pa          |    2 | real       | kind_phys | in     | F        |
!! | rain             | lwe_thickness_of_precipitation_amount_on_dynamics_timestep              | total rain at this time step                                            | m           |    1 | real       | kind_phys | inout  | F        |
!! | domr_diag        | dominant_rain_type                                                      | dominant rain type                                                      | none        |    1 | real       | kind_phys | inout  | F        |
!! | domzr_diag       | dominant_freezing_rain_type                                             | dominant freezing rain type                                             | none        |    1 | real       | kind_phys | inout  | F        |
!! | domip_diag       | dominant_sleet_type                                                     | dominant sleet type                                                     | none        |    1 | real       | kind_phys | inout  | F        |
!! | doms_diag        | dominant_snow_type                                                      | dominant snow type                                                      | none        |    1 | real       | kind_phys | inout  | F        |
!! | tprcp            | nonnegative_lwe_thickness_of_precipitation_amount_on_dynamics_timestep  | total precipitation amount in each time step                            | m           |    1 | real       | kind_phys | inout  | F        |
!! | srflag           | flag_for_precipitation_type                                             | snow/rain flag for precipitation                                        | flag        |    1 | real       | kind_phys | inout  | F        |
!! | sr               | ratio_of_snowfall_to_rainfall                                           | snow ratio: ratio of snow to total precipitation                        | frac        |    1 | real       | kind_phys | in     | F        |
!! | cnvprcp          | cumulative_lwe_thickness_of_convective_precipitation_amount             | cumulative convective precipitation                                     | m           |    1 | real       | kind_phys | inout  | F        |
!! | totprcp          | accumulated_lwe_thickness_of_precipitation_amount                       | accumulated total precipitation                                         | m           |    1 | real       | kind_phys | inout  | F        |
!! | totice           | accumulated_lwe_thickness_of_ice_amount                                 | accumulated ice precipitation                                           | kg m-2      |    1 | real       | kind_phys | inout  | F        |
!! | totsnw           | accumulated_lwe_thickness_of_snow_amount                                | accumulated snow precipitation                                          | kg m-2      |    1 | real       | kind_phys | inout  | F        |
!! | totgrp           | accumulated_lwe_thickness_of_graupel_amount                             | accumulated graupel precipitation                                       | kg m-2      |    1 | real       | kind_phys | inout  | F        |
!! | cnvprcpb         | cumulative_lwe_thickness_of_convective_precipitation_amount_in_bucket   | cumulative convective precipitation in bucket                           | m           |    1 | real       | kind_phys | inout  | F        |
!! | totprcpb         | accumulated_lwe_thickness_of_precipitation_amount_in_bucket             | accumulated total precipitation in bucket                               | m           |    1 | real       | kind_phys | inout  | F        |
!! | toticeb          | accumulated_lwe_thickness_of_ice_amount_in_bucket                       | accumulated ice precipitation in bucket                                 | kg m-2      |    1 | real       | kind_phys | inout  | F        |
!! | totsnwb          | accumulated_lwe_thickness_of_snow_amount_in_bucket                      | accumulated snow precipitation in bucket                                | kg m-2      |    1 | real       | kind_phys | inout  | F        |
!! | totgrpb          | accumulated_lwe_thickness_of_graupel_amount_in_bucket                   | accumulated graupel precipitation in bucket                             | kg m-2      |    1 | real       | kind_phys | inout  | F        |
!! | dt3dt            | cumulative_change_in_temperature_due_to_microphysics                    | cumulative change in temperature due to microphysics                    | K           |    2 | real       | kind_phys | inout  | F        |
!! | dq3dt            | cumulative_change_in_water_vapor_specific_humidity_due_to_microphysics  | cumulative change in water vapor specific humidity due to microphysics  | kg kg-1     |    2 | real       | kind_phys | inout  | F        |
!! | rain_cpl         | lwe_thickness_of_precipitation_amount_for_coupling                      | total rain precipitation                                                | m           |    1 | real       | kind_phys | inout  | F        |
!! | rainc_cpl        | lwe_thickness_of_convective_precipitation_amount_for_coupling           | total convective precipitation                                          | m           |    1 | real       | kind_phys | inout  | F        |
!! | snow_cpl         | lwe_thickness_of_snow_amount_for_coupling                               | total snow precipitation                                                | m           |    1 | real       | kind_phys | inout  | F        |
!! | pwat             | column_precipitable_water                                               | precipitable water                                                      | kg m-2      |    1 | real       | kind_phys | inout  | F        |
!! | do_sppt          | flag_for_stochastic_surface_physics_perturbations                       | flag for stochastic surface physics perturbations                       | flag        |    0 | logical    |           | in     | F        |
!! | dtdtr            | tendency_of_air_temperature_due_to_radiative_heating_on_physics_time_step| temp. change due to radiative heating per time step                    | K           |    2 | real       | kind_phys | inout  | F        |
!! | dtdtc            | tendency_of_air_temperature_due_to_radiative_heating_assuming_clear_sky | clear sky radiative (shortwave + longwave) heating rate at current time | K s-1       |    2 | real       | kind_phys | in     | F        |
!! | drain_cpl        | tendency_of_lwe_thickness_of_precipitation_amount_for_coupling          | change in rain_cpl (coupling_type)                                      | m           |    1 | real       | kind_phys | inout  | F        |
!! | dsnow_cpl        | tendency_of_lwe_thickness_of_snow_amount_for_coupling                   | change in show_cpl (coupling_type)                                      | m           |    1 | real       | kind_phys | inout  | F        |
!! | lsm              | flag_for_land_surface_scheme                                            | flag for land surface model                                             | flag        |    0 | integer    |           | in     | F        |
!! | lsm_ruc          | flag_for_ruc_land_surface_scheme                                        | flag for RUC land surface model                                         | flag        |    0 | integer    |           | in     | F        |
!! | raincprv         | lwe_thickness_of_convective_precipitation_amount_from_previous_timestep | convective_precipitation_amount from previous timestep                  | m           |    1 | real       | kind_phys | inout  | F        |
!! | rainncprv        | lwe_thickness_of_explicit_rainfall_amount_from_previous_timestep        | explicit rainfall from previous timestep                                | m           |    1 | real       | kind_phys | inout  | F        |
!! | iceprv           | lwe_thickness_of_ice_amount_from_previous_timestep                      | ice amount from previous timestep                                       | m           |    1 | real       | kind_phys | inout  | F        |
!! | snowprv          | lwe_thickness_of_snow_amount_from_previous_timestep                     | snow amount from previous timestep                                      | m           |    1 | real       | kind_phys | inout  | F        |
!! | graupelprv       | lwe_thickness_of_graupel_amount_from_previous_timestep                  | graupel amount from previous timestep                                   | m           |    1 | real       | kind_phys | inout  | F        |
!! | dtp              | time_step_for_physics                                                   | physics timestep                                                        | s           |    0 | real       | kind_phys | in     | F        |
!! | errmsg           | ccpp_error_message                                                      | error message for error handling in CCPP                                | none        |    0 | character  | len=*     | out    | F        |
!! | errflg           | ccpp_error_flag                                                         | error flag for error handling in CCPP                                   | flag        |    0 | integer    |           | out    | F        |
!!
!> \section gfs_mp_gen GFS MP Generic Post General Algorithm
!> @{
      subroutine GFS_MP_generic_post_run(im, ix, levs, kdt, nrcm, ncld, nncl, ntcw, ntrac, imp_physics, imp_physics_gfdl, &
        imp_physics_thompson, imp_physics_mg, cal_pre, lssav, ldiag3d, cplflx, cplchm, con_g, dtf, frain, rainc, rain1,   &
        rann, xlat, xlon, gt0, gq0, prsl, prsi, phii, tsfc, ice, snow, graupel, save_t, save_qv, rain0, ice0, snow0,      &
        graupel0, del, rain, domr_diag, domzr_diag, domip_diag, doms_diag, tprcp, srflag, sr, cnvprcp, totprcp, totice,   &
        totsnw, totgrp, cnvprcpb, totprcpb, toticeb, totsnwb, totgrpb, dt3dt, dq3dt, rain_cpl, rainc_cpl, snow_cpl, pwat, &
        do_sppt, dtdtr, dtdtc, drain_cpl, dsnow_cpl, lsm, lsm_ruc, raincprv, rainncprv, iceprv, snowprv, graupelprv,      &
        dtp, errmsg, errflg)
!
      use machine, only: kind_phys

      implicit none

      integer, intent(in) :: im, ix, levs, kdt, nrcm, ncld, nncl, ntcw, ntrac
      integer, intent(in) :: imp_physics, imp_physics_gfdl, imp_physics_thompson, imp_physics_mg
      logical, intent(in) :: cal_pre, lssav, ldiag3d, cplflx, cplchm

      real(kind=kind_phys),                           intent(in)    :: dtf, frain, con_g
      real(kind=kind_phys), dimension(im),            intent(in)    :: rainc, rain1, xlat, xlon, tsfc
      real(kind=kind_phys), dimension(im),            intent(inout) :: ice, snow, graupel
      real(kind=kind_phys), dimension(im),            intent(in)    :: rain0, ice0, snow0, graupel0
      real(kind=kind_phys), dimension(ix,nrcm),       intent(in)    :: rann
      real(kind=kind_phys), dimension(im,levs),       intent(in)    :: gt0, prsl, save_t, save_qv, del
      real(kind=kind_phys), dimension(im,levs+1),     intent(in)    :: prsi, phii
      real(kind=kind_phys), dimension(im,levs,ntrac), intent(in)    :: gq0

      real(kind=kind_phys), dimension(im),      intent(in   ) :: sr
      real(kind=kind_phys), dimension(im),      intent(inout) :: rain, domr_diag, domzr_diag, domip_diag, doms_diag, tprcp,  &
                                                                 srflag, cnvprcp, totprcp, totice, totsnw, totgrp, cnvprcpb, &
                                                                 totprcpb, toticeb, totsnwb, totgrpb, rain_cpl, rainc_cpl,   &
                                                                 snow_cpl, pwat
      ! These arrays are only allocated if ldiag3d is .true.
      real(kind=kind_phys), dimension(:,:),     intent(inout) :: dt3dt, dq3dt

      ! Stochastic physics / surface perturbations
      logical, intent(in) :: do_sppt
      real(kind=kind_phys), dimension(im,levs), intent(inout) :: dtdtr
      real(kind=kind_phys), dimension(im,levs), intent(in)    :: dtdtc
      real(kind=kind_phys), dimension(im),      intent(inout) :: drain_cpl
      real(kind=kind_phys), dimension(im),      intent(inout) :: dsnow_cpl

      ! Rainfall variables previous time step (update for RUC LSM)
      integer, intent(in) :: lsm, lsm_ruc
      real(kind=kind_phys), dimension(im),      intent(inout) :: raincprv
      real(kind=kind_phys), dimension(im),      intent(inout) :: rainncprv
      real(kind=kind_phys), dimension(im),      intent(inout) :: iceprv
      real(kind=kind_phys), dimension(im),      intent(inout) :: snowprv
      real(kind=kind_phys), dimension(im),      intent(inout) :: graupelprv

      real(kind=kind_phys),                     intent(in)    :: dtp

      ! CCPP error handling
      character(len=*), intent(out) :: errmsg
      integer, intent(out) :: errflg

      ! DH* TODO: CLEANUP, all of these should be coming in through the argument list
      real(kind=kind_phys), parameter :: con_p001= 0.001d0
      real(kind=kind_phys), parameter :: con_day = 86400.0d0
      real(kind=kind_phys), parameter :: rainmin = 1.0d-13
      real(kind=kind_phys), parameter :: p850    = 85000.0d0
      ! *DH

      integer :: i, k, ic

      real(kind=kind_phys), parameter :: zero = 0.0d0, one = 1.0d0
      real(kind=kind_phys) :: crain, csnow, onebg, tem, total_precip
      real(kind=kind_phys), dimension(im) :: domr, domzr, domip, doms, t850, work1

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      onebg = one/con_g

      do i = 1, im
          rain(i) = rainc(i) + frain * rain1(i) ! time-step convective plus explicit
      enddo

!> - If requested (e.g. Zhao-Carr MP scheme), call calpreciptype() to calculate dominant 
!! precipitation type.
      ! DH* TODO - Fix wrong code in non-CCPP build (GFS_physics_driver)
      ! and use commented lines here (keep wrong version for bit-for-bit):
      ! put ice, snow, graupel on dynamics timestep. The way the code in
      ! GFS_physics_driver is written, Diag%{graupel,ice,snow} are on the
      ! physics timestep, while Diag%{rain,rainc} and all totprecip etc
      ! are on the dynamics timestep. Confusing, but works if frain=1. *DH
      if (imp_physics == imp_physics_gfdl) then
        tprcp   = max(0., rain)               ! clu: rain -> tprcp
        !graupel = frain*graupel0
        !ice     = frain*ice0
        !snow    = frain*snow0
        graupel = graupel0
        ice     = ice0
        snow    = snow0
      ! Do it right from the beginning for Thompson
      else if (imp_physics == imp_physics_thompson) then
        tprcp   = max (0.,rainc + frain * rain1) ! time-step convective and explicit precip
        graupel = frain*graupel0              ! time-step graupel
        ice     = frain*ice0                  ! time-step ice
        snow    = frain*snow0                 ! time-step snow
      end if

      if (lsm==lsm_ruc) then
        if (imp_physics == imp_physics_gfdl .or. imp_physics == imp_physics_thompson) then
            raincprv(:)   = rainc(:)
            rainncprv(:)  = frain * rain1(:)
            iceprv(:)     = ice(:)
            snowprv(:)    = snow(:)
            graupelprv(:) = graupel(:)
        end if
      end if

      if (cal_pre) then       ! hchuang: add dominant precipitation type algorithm
!
        call calpreciptype (kdt, nrcm, im, ix, levs, levs+1,           &
                            rann, xlat, xlon, gt0,    &
                            gq0(:,:,1), prsl, prsi,        &
                            rain, phii, tsfc,           &  !input
                            domr, domzr, domip, doms)                           ! output
!
!        if (lprnt) print*,'debug calpreciptype: DOMR,DOMZR,DOMIP,DOMS '
!     &,DOMR(ipr),DOMZR(ipr),DOMIP(ipr),DOMS(ipr)
!        do i=1,im
!         if (abs(xlon(i)*57.29578-114.0) .lt. 0.2  .and.
!     &    abs(xlat(i)*57.29578-40.0) .lt. 0.2)
!     &    print*,'debug calpreciptype: DOMR,DOMZR,DOMIP,DOMS ',
!     &    DOMR(i),DOMZR(i),DOMIP(i),DOMS(i)
!       end do
!       HCHUANG: use new precipitation type to decide snow flag for LSM snow accumulation

        if (imp_physics /= imp_physics_gfdl .and. imp_physics /= imp_physics_thompson) then
          do i=1,im
            tprcp(i)  = max(0.0, rain(i) )
            if(doms(i) > 0.0 .or. domip(i) > 0.0) then
              srflag(i) = 1.
            else
              srflag(i) = 0.
            end if
          enddo
        endif
        if (lssav) then
          do i=1,im
              domr_diag(i)  = domr_diag(i)  + domr(i)  * dtf
              domzr_diag(i) = domzr_diag(i) + domzr(i) * dtf
              domip_diag(i) = domip_diag(i) + domip(i) * dtf
              doms_diag(i)  = doms_diag(i)  + doms(i)  * dtf
          enddo
        endif

      endif

      if (lssav) then
!        if (Model%me == 0) print *,'in phys drive, kdt=',Model%kdt, &
!          'totprcpb=', Diag%totprcpb(1),'totprcp=',Diag%totprcp(1), &
!          'rain=',Diag%rain(1)
        do i=1,im
          cnvprcp (i) = cnvprcp (i) + rainc(i)
          totprcp (i) = totprcp (i) + rain(i)
          totice  (i) = totice  (i) + ice(i)
          totsnw  (i) = totsnw  (i) + snow(i)
          totgrp  (i) = totgrp  (i) + graupel(i)

          cnvprcpb(i) = cnvprcpb(i) + rainc(i)
          totprcpb(i) = totprcpb(i) + rain(i)
          toticeb (i) = toticeb (i) + ice(i)
          totsnwb (i) = totsnwb (i) + snow(i)
          totgrpb (i) = totgrpb (i) + graupel(i)
        enddo

        if (ldiag3d) then
          do k=1,levs
            do i=1,im
              dt3dt(i,k) = dt3dt(i,k) + (gt0(i,k)-save_t(i,k)) * frain
!              dq3dt(i,k) = dq3dt(i,k) + (gq0(i,k,1)-save_qv(i,k)) * frain
            enddo
          enddo
        endif
      endif

      t850(1:im) = gt0(1:im,1)

      do k = 1, levs-1
        do i = 1, im
          if (prsl(i,k) > p850 .and. prsl(i,k+1) <= p850) then
            t850(i) = gt0(i,k) - (prsl(i,k)-p850) / &
                      (prsl(i,k)-prsl(i,k+1)) *      &
                      (gt0(i,k)-gt0(i,k+1))
          endif
        enddo
      enddo

      ! Conversion factor mm per physics timestep to m per day
      tem = dtp * con_p001 / con_day

!> - For GFDL and Thompson MP scheme, determine convective snow by surface temperature;
!! and determine explicit rain/snow by snow/ice/graupel coming out directly from MP
!! and convective rainfall from the cumulus scheme if the surface temperature is below
!! \f$0^oC\f$.
      if (imp_physics == imp_physics_gfdl .or. imp_physics == imp_physics_thompson) then
! determine convective rain/snow by surface temperature
! determine large-scale rain/snow by rain/snow coming out directly from MP
        do i = 1, im
          !tprcp(i)  = max(0.0, rain(i) )! clu: rain -> tprcp ! DH now lines 245-250
          srflag(i) = 0.                     ! clu: default srflag as 'rain' (i.e. 0)
          if (tsfc(i) >= 273.15) then
            crain = rainc(i)
            csnow = 0.0
          else
            crain = 0.0
            csnow = rainc(i)
          endif
!          if (snow0(i,1)+ice0(i,1)+graupel0(i,1)+csnow > rain0(i,1)+crain) then
!          if (snow0(i)+ice0(i)+graupel0(i)+csnow > 0.0) then
!            Sfcprop%srflag(i) = 1.                   ! clu: set srflag to 'snow' (i.e. 1)
!          endif
! compute fractional srflag
          total_precip = snow0(i)+ice0(i)+graupel0(i)+rain0(i)+rainc(i)
          if (total_precip > rainmin) then
            srflag(i) = (snow0(i)+ice0(i)+graupel0(i)+csnow)/total_precip
          endif
        enddo
      elseif( .not. cal_pre) then
        if (imp_physics == imp_physics_mg) then              ! MG microphysics
          do i=1,im
            if (rain(i)*tem > rainmin) then
              srflag(i) = max(zero, min(one, (rain(i)-rainc(i))*sr(i)/rain(i)))
            else
              srflag(i) = 0.0
            endif
          enddo
        else
          do i = 1, im
            tprcp(i)  = max(0.0, rain(i) )! clu: rain -> tprcp
            srflag(i) = 0.0                    ! clu: default srflag as 'rain' (i.e. 0)
            if (t850(i) <= 273.16) then
              srflag(i) = 1.0                  ! clu: set srflag to 'snow' (i.e. 1)
            endif
          enddo
        endif
      endif

      if (cplflx .or. cplchm) then
        do i = 1, im
          rain_cpl(i) = rain_cpl(i) + rain(i) * (one-srflag(i))
          snow_cpl(i) = snow_cpl(i) + rain(i) * srflag(i)
        enddo
      endif

      if (cplchm) then
        do i = 1, im
             rainc_cpl(i) = rainc_cpl(i) + rainc(i)
        enddo
      endif

      pwat(:) = 0.0
      do k = 1, levs
        do i=1, im
          work1(i) = 0.0
        enddo
        if (ncld > 0) then
          do ic = ntcw, ntcw+nncl-1
            do i=1,im
              work1(i) = work1(i) + gq0(i,k,ic)
            enddo
          enddo
        endif
        do i=1,im
          pwat(i) = pwat(i) + del(i,k)*(gq0(i,k,1)+work1(i))
        enddo
!     if (lprnt .and. i == ipr) write(0,*)' gq0=',
!    &gq0(i,k,1),' qgrs=',qgrs(i,k,1),' work2=',work2(i),' k=',k
      enddo
      do i=1,im
        pwat(i) = pwat(i) * onebg
      enddo

      ! Stochastic physics / surface perturbations
      if (do_sppt) then
!--- radiation heating rate
        dtdtr(1:im,:) = dtdtr(1:im,:) + dtdtc(1:im,:)*dtf
        do i = 1, im
          if (t850(i) > 273.16) then
!--- change in change in rain precip
             drain_cpl(i) = rain(i) - drain_cpl(i)
          else
!--- change in change in snow precip
             dsnow_cpl(i) = rain(i) - dsnow_cpl(i)
          endif
        enddo
      endif

    end subroutine GFS_MP_generic_post_run
!> @}

!> \section arg_table_GFS_MP_generic_post_finalize Argument Table
!!
      subroutine GFS_MP_generic_post_finalize
      end subroutine GFS_MP_generic_post_finalize
      end module GFS_MP_generic_post

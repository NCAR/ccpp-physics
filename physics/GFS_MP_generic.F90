!> \file GFS_MP_generic.F90
!! This file contains the subroutines that calculate diagnotics variables
!! before/after calling any microphysics scheme:

      module GFS_MP_generic_pre
      contains

!> \defgroup GFS_MP_generic_pre GFS MP generic pre
!! @{
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
!! | gq0_water_vapor| water_vapor_specific_humidity_updated_by_physics       | water vapor specific humidity updated by physics                          | kg kg-1     |    2 | real      | kind_phys | in     | F        |
!! | gq0            | tracer_concentration_updated_by_physics                | tracer concentration updated by physics                                   | kg kg-1     |    3 | real      | kind_phys | in     | F        |
!! | save_t         | air_temperature_save                                   | air temperature before entering a physics scheme                          | K           |    2 | real      | kind_phys | inout  | F        |
!! | save_qv        | water_vapor_specific_humidity_save                     | water vapor specific humidity before entering a physics scheme            | kg kg-1     |    2 | real      | kind_phys | inout  | F        |
!! | save_q         | tracer_concentration_save                              | tracer concentration before entering a physics scheme                     | kg kg-1     |    3 | real      | kind_phys | inout  | F        |
!! | errmsg         | ccpp_error_message                                     | error message for error handling in CCPP                                  | none        |    0 | character | len=*     | out    | F        |
!! | errflg         | ccpp_error_flag                                        | error flag for error handling in CCPP                                     | flag        |    0 | integer   |           | out    | F        |
!!
      subroutine GFS_MP_generic_pre_run(im, levs, ldiag3d, do_aw, ntcw, nncl, ntrac, gt0, gq0_water_vapor, gq0, &
        save_t, save_qv, save_q, errmsg, errflg)
!
      use machine,               only: kind_phys

      implicit none

      integer,                                          intent(in) :: im, levs, ntcw, nncl, ntrac
      logical,                                          intent(in) :: ldiag3d, do_aw
      real(kind=kind_phys), dimension(im, levs),        intent(in) :: gt0, gq0_water_vapor
      real(kind=kind_phys), dimension(im, levs, ntrac), intent(in) :: gq0

      real(kind=kind_phys), dimension(im, levs),        intent(inout) :: save_t, save_qv
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
            save_t(i,k)   = gt0(i,k)
            save_qv(i,k)  = gq0_water_vapor(i,k)
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
!! @}
      end module GFS_MP_generic_pre

      module GFS_MP_generic_post
      contains

!> \defgroup GFS_MP_generic_post GFS MP generic post
!! @{
!! \section arg_table_GFS_MP_generic_post_init Argument Table
!!
      subroutine GFS_MP_generic_post_init
      end subroutine GFS_MP_generic_post_init


!> \section arg_table_GFS_MP_generic_post_run Argument Table
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
!! | cal_pre          | flag_for_precipitation_type_algorithm                                   | flag controls precip type algorithm                                     | flag        |    0 | logical    |           | in     | F        |
!! | lssav            | flag_diagnostics                                                        | logical flag for storing diagnostics                                    | flag        |    0 | logical    |           | in     | F        |
!! | ldiag3d          | flag_diagnostics_3D                                                     | flag for 3d diagnostic fields                                           | flag        |    0 | logical    |           | in     | F        |
!! | cplflx           | flag_for_flux_coupling                                                  | flag controlling cplflx collection (default off)                        | flag        |    0 | logical    |           | in     | F        |
!! | cplchm           | flag_for_chemistry_coupling                                             | flag controlling cplchm collection (default off)                        | flag        |    0 | logical    |           | in     | F        |
!! | con_g            | gravitational_acceleration                                              | gravitational acceleration                                              | m s-2       |    0 | real       | kind_phys | in     | F        |
!! | dtf              | time_step_for_dynamics                                                  | dynamics timestep                                                       | s           |    0 | real       | kind_phys | in     | F        |
!! | frain            | dynamics_to_physics_timestep_ratio                                      | ratio of dynamics timestep to physics timestep                          | none        |    0 | real       | kind_phys | in     | F        |
!! | rainc            | lwe_thickness_of_convective_precipitation_amount_on_dynamics_timestep   | convective rain at this time step                                       | m           |    1 | real       | kind_phys | in     | F        |
!! | rain1            | lwe_thickness_of_stratiform_precipitation_amount                        | stratiform rainfall amount on physics timestep                          | m           |    1 | real       | kind_phys | in     | F        |
!! | rann             | random_number_array                                                     | random number array (0-1)                                               | none        |    2 | real       | kind_phys | in     | F        |
!! | xlat             | latitude                                                                | latitude                                                                | radians     |    1 | real       | kind_phys | in     | F        |
!! | xlon             | longitude                                                               | longitude                                                               | radians     |    1 | real       | kind_phys | in     | F        |
!! | gt0              | air_temperature_updated_by_physics                                      | temperature updated by physics                                          | K           |    2 | real       | kind_phys | in     | F        |
!! | gq0              | tracer_concentration_updated_by_physics                                 | tracer concentration updated by physics                                 | kg kg-1     |    3 | real       | kind_phys | in     | F        |
!! | gq0_water_vapor  | water_vapor_specific_humidity_updated_by_physics                        | water vapor specific humidity updated by physics                        | kg kg-1     |    2 | real       | kind_phys | in     | F        |
!! | prsl             | air_pressure                                                            | layer mean pressure                                                     | Pa          |    2 | real       | kind_phys | in     | F        |
!! | prsi             | air_pressure_at_interface                                               | pressure at layer interface                                             | Pa          |    2 | real       | kind_phys | in     | F        |
!! | phii             | geopotential_at_interface                                               | geopotential at model layer interfaces                                  | m2 s-2      |    2 | real       | kind_phys | in     | F        |
!! | tsfc             | surface_skin_temperature                                                | surface skin temperature                                                | K           |    1 | real       | kind_phys | in     | F        |
!! | ice              | lwe_thickness_of_ice_amount_on_dynamics_timestep                        | ice fall at this time step                                              | m           |    1 | real       | kind_phys | in     | F        |
!! | snow             | lwe_thickness_of_snow_amount_on_dynamics_timestep                       | snow fall at this time step                                             | m           |    1 | real       | kind_phys | in     | F        |
!! | graupel          | lwe_thickness_of_graupel_amount_on_dynamics_timestep                    | graupel fall at this time step                                          | m           |    1 | real       | kind_phys | in     | F        |
!! | save_t           | air_temperature_save                                                    | air temperature before entering a physics scheme                        | K           |    2 | real       | kind_phys | in     | F        |
!! | save_qv          | water_vapor_specific_humidity_save                                      | water vapor specific humidity before entering a physics scheme          | kg kg-1     |    2 | real       | kind_phys | in     | F        |
!! | ice0             | lwe_thickness_of_ice_amount_per_day                                     | ice fall over 24h period                                                | mm          |    1 | real       | kind_phys | in     | F        |
!! | snow0            | lwe_thickness_of_snow_amount_per_day                                    | snow fall over 24h period                                               | mm          |    1 | real       | kind_phys | in     | F        |
!! | graupel0         | lwe_thickness_of_graupel_amount_per_day                                 | graupel fall over 24h period                                            | mm          |    1 | real       | kind_phys | in     | F        |
!! | del              | air_pressure_difference_between_midlayers                               | air pressure difference between midlayers                               | Pa          |    2 | real       | kind_phys | in     | F        |
!! | rain             | lwe_thickness_of_precipitation_amount_on_dynamics_timestep              | total rain at this time step                                            | m           |    1 | real       | kind_phys | inout  | F        |
!! | domr_diag        | dominant_rain_type                                                      | dominant rain type                                                      | none        |    1 | real       | kind_phys | inout  | F        |
!! | domzr_diag       | dominant_freezing_rain_type                                             | dominant freezing rain type                                             | none        |    1 | real       | kind_phys | inout  | F        |
!! | domip_diag       | dominant_sleet_type                                                     | dominant sleet type                                                     | none        |    1 | real       | kind_phys | inout  | F        |
!! | doms_diag        | dominant_snow_type                                                      | dominant snow type                                                      | none        |    1 | real       | kind_phys | inout  | F        |
!! | tprcp            | nonnegative_lwe_thickness_of_precipitation_amount_on_dynamics_timestep  | total precipitation amount in each time step                            | m           |    1 | real       | kind_phys | inout  | F        |
!! | srflag           | flag_for_precipitation_type                                             | snow/rain flag for precipitation                                        | flag        |    1 | real       | kind_phys | inout  | F        |
!! | totprcp          | accumulated_lwe_thickness_of_precipitation_amount                       | accumulated total precipitation                                         | m           |    1 | real       | kind_phys | inout  | F        |
!! | totice           | accumulated_lwe_thickness_of_ice_amount                                 | accumulated ice precipitation                                           | kg m-2      |    1 | real       | kind_phys | inout  | F        |
!! | totsnw           | accumulated_lwe_thickness_of_snow_amount                                | accumulated snow precipitation                                          | kg m-2      |    1 | real       | kind_phys | inout  | F        |
!! | totgrp           | accumulated_lwe_thickness_of_graupel_amount                             | accumulated graupel precipitation                                       | kg m-2      |    1 | real       | kind_phys | inout  | F        |
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
!! | errmsg           | ccpp_error_message                                                      | error message for error handling in CCPP                                | none        |    0 | character  | len=*     | out    | F        |
!! | errflg           | ccpp_error_flag                                                         | error flag for error handling in CCPP                                   | flag        |    0 | integer    |           | out    | F        |
!!
      subroutine GFS_MP_generic_post_run(im, ix, levs, kdt, nrcm, ncld, nncl, ntcw, ntrac, imp_physics, imp_physics_gfdl,       &
        cal_pre, lssav, ldiag3d, cplflx, cplchm, con_g, dtf, frain, rainc, rain1, rann, xlat, xlon, gt0, gq0, gq0_water_vapor,  &
        prsl, prsi, phii, tsfc, ice, snow, graupel, save_t, save_qv, ice0, snow0, graupel0, del,                                &
        rain, domr_diag, domzr_diag, domip_diag, doms_diag, tprcp, srflag, totprcp, totice, totsnw,                             &
        totgrp, totprcpb, toticeb, totsnwb, totgrpb, dt3dt, dq3dt, rain_cpl, rainc_cpl, snow_cpl, pwat, errmsg, errflg)
!
      use machine,               only: kind_phys

      implicit none

      integer, intent(in) :: im, ix, levs, kdt, nrcm, ncld, nncl, ntcw, ntrac, imp_physics, imp_physics_gfdl
      logical, intent(in) :: cal_pre, lssav, ldiag3d, cplflx, cplchm

      real(kind=kind_phys),                           intent(in) :: dtf, frain, con_g
      real(kind=kind_phys), dimension(im),            intent(in) :: rainc, rain1, xlat, xlon, tsfc, ice, snow, graupel,         &
                                                                    ice0, snow0, graupel0
      real(kind=kind_phys), dimension(ix,nrcm),       intent(in) :: rann
      real(kind=kind_phys), dimension(im,levs),       intent(in) :: gt0, gq0_water_vapor, prsl, save_t, save_qv, del
      real(kind=kind_phys), dimension(im,levs+1),     intent(in) :: prsi, phii
      real(kind=kind_phys), dimension(im,levs,ntrac), intent(in) :: gq0


      real(kind=kind_phys), dimension(im),      intent(inout) :: rain, domr_diag, domzr_diag, domip_diag, doms_diag, tprcp,     &
        srflag, totprcp, totice, totsnw, totgrp, totprcpb, toticeb, totsnwb, totgrpb, rain_cpl, rainc_cpl, snow_cpl, pwat
      real(kind=kind_phys), dimension(im,levs), intent(inout) :: dt3dt, dq3dt

      character(len=*), intent(out) :: errmsg
      integer, intent(out) :: errflg

      real(kind=kind_phys), parameter :: p850    = 85000.0

      integer :: i, k, ic
      real(kind=kind_phys) :: crain, csnow, onebg
      real(kind=kind_phys), dimension(im) :: domr, domzr, domip, doms, t850, work1

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      onebg = 1.0d0/con_g

      do i = 1, im
          rain(i) = rainc(i) + frain * rain1(i)
      enddo

      if (cal_pre) then       ! hchuang: add dominant precipitation type algorithm
!
        call calpreciptype (kdt, nrcm, im, ix, levs, levs+1,           &
                            rann, xlat, xlon, gt0,    &
                            gq0_water_vapor, prsl, prsi,        &
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

        if (imp_physics /= imp_physics_gfdl) then
          do i=1,im
            tprcp(i)  = max(0.0, rain(i) )
            if(doms(i) > 0.0 .or. domip(i) > 0.0) then
              srflag(i) = 1.
            else
              srflag(i) = 0.
            end if
          enddo
        endif

      endif

      if (lssav) then
!        if (Model%me == 0) print *,'in phys drive, kdt=',Model%kdt, &
!          'totprcpb=', Diag%totprcpb(1),'totprcp=',Diag%totprcp(1), &
!          'rain=',Diag%rain(1)
        do i=1,im
          totprcp (i) = totprcp (i) + rain(i)
          totice  (i) = totice  (i) + ice(i)
          totsnw  (i) = totsnw  (i) + snow(i)
          totgrp  (i) = totgrp  (i) + graupel(i)
          totprcpb(i) = totprcpb(i) + rain(i)
          toticeb (i) = toticeb (i) + ice(i)
          totsnwb (i) = totsnwb (i) + snow(i)
          totgrpb (i) = totgrpb (i) + graupel(i)
!
          if (cal_pre) then
            domr_diag(i)  = domr_diag(i)  + domr(i)  * dtf
            domzr_diag(i) = domzr_diag(i) + domzr(i) * dtf
            domip_diag(i) = domip_diag(i) + domip(i) * dtf
            doms_diag(i)  = doms_diag(i)  + doms(i)  * dtf
          endif
        enddo

        if (ldiag3d) then
          do k=1,levs
            do i=1,im
              dt3dt(i,k) = dt3dt(i,k) + (gt0(i,k)-save_t(i,k)) * frain
              dq3dt(i,k) = dq3dt(i,k) + (gq0_water_vapor(i,k)-save_qv(i,k)) * frain
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

      if (imp_physics == imp_physics_gfdl) then
! determine convective rain/snow by surface temperature
! determine large-scale rain/snow by rain/snow coming out directly from MP
        do i = 1, im
          tprcp(i)  = max(0.0, rain(i) )! clu: rain -> tprcp
          srflag(i) = 0.                     ! clu: default srflag as 'rain' (i.e. 0)
          if (tsfc(i) .ge. 273.15) then
            crain = rainc(i)
            csnow = 0.0
          else
            crain = 0.0
            csnow = rainc(i)
          endif
!         if ((snow0(i,1)+ice0(i,1)+graupel0(i,1)+csnow) > (rain0(i,1)+crain)) then
          if ((snow0(i)+ice0(i)+graupel0(i)+csnow) > 0.0) then
            srflag(i) = 1.                   ! clu: set srflag to 'snow' (i.e. 1)
          endif
        enddo
      elseif( .not. cal_pre) then
        do i = 1, im
          tprcp(i)  = max(0.0, rain(i) )! clu: rain -> tprcp
          srflag(i) = 0.                     ! clu: default srflag as 'rain' (i.e. 0)
          if (t850(i) <= 273.16) then
            srflag(i) = 1.                   ! clu: set srflag to 'snow' (i.e. 1)
          endif
        enddo
      endif

      if (cplflx) then
        do i = 1, im
          if (t850(i) > 273.16) then
             rain_cpl(i) = rain_cpl(i) + rain(i)
          else
             snow_cpl(i) = snow_cpl(i) + rain(i)
          endif
        enddo
      endif

      if ((cplchm).and.(.not. cplflx)) then
        do i = 1, im
             rain_cpl(i) = rain_cpl(i) + rain(i)
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
          pwat(i) = pwat(i) + del(i,k)*(gq0_water_vapor(i,k)+work1(i))
        enddo
!     if (lprnt .and. i == ipr) write(0,*)' gq0=',
!    &gq0(i,k,1),' qgrs=',qgrs(i,k,1),' work2=',work2(i),' k=',k
      enddo
      do i=1,im
        pwat(i) = pwat(i) * onebg
      enddo

    end subroutine GFS_MP_generic_post_run

!> \section arg_table_GFS_MP_generic_post_finalize Argument Table
!!
      subroutine GFS_MP_generic_post_finalize
      end subroutine GFS_MP_generic_post_finalize
!! @}
      end module GFS_MP_generic_post

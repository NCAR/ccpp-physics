!> \file GFS_DCNV_generic.F90
!!  Contains code related to deep convective schemes to be used within the GFS physics suite.

      module GFS_DCNV_generic_pre

      contains

      subroutine GFS_DCNV_generic_pre_init ()
      end subroutine GFS_DCNV_generic_pre_init

      subroutine GFS_DCNV_generic_pre_finalize()
      end subroutine GFS_DCNV_generic_pre_finalize

!> \brief Interstitial scheme called prior to any deep convective scheme to save state variables for calculating tendencies after the deep convective scheme is executed
!! \section arg_table_GFS_DCNV_generic_pre_run Argument Table
!! | local_name      | standard_name                                          | long_name                                                                 | units   | rank | type      |    kind   | intent | optional |
!! |-----------------|--------------------------------------------------------|---------------------------------------------------------------------------|---------|------|-----------|-----------|--------|----------|
!! | im              | horizontal_loop_extent                                 | horizontal loop extent                                                    | count   |    0 | integer   |           | in     | F        |
!! | levs            | vertical_dimension                                     | vertical layer dimension                                                  | count   |    0 | integer   |           | in     | F        |
!! | ldiag3d         | flag_diagnostics_3D                                    | flag for 3d diagnostic fields                                             | flag    |    0 | logical   |           | in     | F        |
!! | cnvgwd          | flag_convective_gravity_wave_drag                      | flag for conv gravity wave drag                                           | flag    |    0 | logical   |           | in     | F        |
!! | lgocart         | flag_gocart                                            | flag for 3d diagnostic fields for gocart 1                                | flag    |    0 | logical   |           | in     | F        |
!! | gu0             | x_wind_updated_by_physics                              | zonal wind updated by physics                                             | m s-1   |    2 | real      | kind_phys | in     | F        |
!! | gv0             | y_wind_updated_by_physics                              | meridional wind updated by physics                                        | m s-1   |    2 | real      | kind_phys | in     | F        |
!! | gt0             | air_temperature_updated_by_physics                     | temperature updated by physics                                            | K       |    2 | real      | kind_phys | in     | F        |
!! | gq0_water_vapor | water_vapor_specific_humidity_updated_by_physics       | water vapor specific humidity updated by physics                          | kg kg-1 |    2 | real      | kind_phys | in     | F        |
!! | save_u          | x_wind_save                                            | x-wind before entering a physics scheme                                   | m s-1   |    2 | real      | kind_phys | inout  | F        |
!! | save_v          | y_wind_save                                            | y-wind before entering a physics scheme                                   | m s-1   |    2 | real      | kind_phys | inout  | F        |
!! | save_t          | air_temperature_save                                   | air temperature before entering a physics scheme                          | K       |    2 | real      | kind_phys | inout  | F        |
!! | save_qv         | water_vapor_specific_humidity_save                     | water vapor specific humidity before entering a physics scheme            | kg kg-1 |    2 | real      | kind_phys | inout  | F        |
!! | errmsg          | ccpp_error_message                                     | error message for error handling in CCPP                                  | none    |    0 | character | len=*     | out    | F        |
!! | errflg          | ccpp_error_flag                                        | error flag for error handling in CCPP                                     | flag    |    0 | integer   |           | out    | F        |
!!
    subroutine GFS_DCNV_generic_pre_run (im, levs, ldiag3d, cnvgwd, lgocart, gu0, gv0, gt0, gq0_water_vapor, &
      save_u, save_v, save_t, save_qv, errmsg, errflg)

      use machine,               only: kind_phys

      implicit none

      integer, intent(in) :: im, levs
      logical, intent(in) :: ldiag3d, cnvgwd, lgocart
      real(kind=kind_phys), dimension(im,levs), intent(in) :: gu0
      real(kind=kind_phys), dimension(im,levs), intent(in) :: gv0
      real(kind=kind_phys), dimension(im,levs), intent(in) :: gt0
      real(kind=kind_phys), dimension(im,levs), intent(in) :: gq0_water_vapor
      real(kind=kind_phys), dimension(im,levs), intent(inout) :: save_u
      real(kind=kind_phys), dimension(im,levs), intent(inout) :: save_v
      real(kind=kind_phys), dimension(im,levs), intent(inout) :: save_t
      real(kind=kind_phys), dimension(im,levs), intent(inout) :: save_qv
      character(len=*), intent(out) :: errmsg
      integer, intent(out) :: errflg

      integer :: i, k

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      if (ldiag3d) then
        do k=1,levs
          do i=1,im
            save_t(i,k) = gt0(i,k)
            save_u(i,k) = gu0(i,k)
            save_v(i,k) = gv0(i,k)
          enddo
        enddo
      elseif (cnvgwd) then
        save_t(1:im,:) = gt0(1:im,:)
      endif   ! end if_ldiag3d/cnvgwd

      if (ldiag3d .or. lgocart) then
        save_qv(1:im,:) = gq0_water_vapor(1:im,:)
      endif   ! end if_ldiag3d/lgocart

    end subroutine GFS_DCNV_generic_pre_run

    end module GFS_DCNV_generic_pre

    module GFS_DCNV_generic_post

    contains

    subroutine GFS_DCNV_generic_post_init ()
    end subroutine GFS_DCNV_generic_post_init

    subroutine GFS_DCNV_generic_post_finalize ()
    end subroutine GFS_DCNV_generic_post_finalize

!> \section arg_table_GFS_DCNV_generic_post_run Argument Table
!! | local_name      | standard_name                                                                               | long_name                                                            | units         | rank | type              |    kind   | intent | optional |
!! |-----------------|---------------------------------------------------------------------------------------------|----------------------------------------------------------------------|---------------|------|-------------------|-----------|--------|----------|
!! | im              | horizontal_loop_extent                                                                      | horizontal loop extent                                               | count         |    0 | integer           |           | in     | F        |
!! | levs            | vertical_dimension                                                                          | vertical layer dimension                                             | count         |    0 | integer           |           | in     | F        |
!! | lssav           | flag_diagnostics                                                                            | logical flag for storing diagnostics                                 | flag          |    0 | logical           |           | in     | F        |
!! | ldiag3d         | flag_diagnostics_3D                                                                         | flag for 3d diagnostic fields                                        | flag          |    0 | logical           |           | in     | F        |
!! | lgocart         | flag_gocart                                                                                 | flag for 3d diagnostic fields for gocart 1                           | flag          |    0 | logical           |           | in     | F        |
!! | frain           | dynamics_to_physics_timestep_ratio                                                          | ratio of dynamics timestep to physics timestep                       | none          |    0 | real              | kind_phys | in     | F        |
!! | rain1           | lwe_thickness_of_deep_convective_precipitation_amount                                       | deep convective rainfall amount on physics timestep                  | m             |    1 | real              | kind_phys | in     | F        |
!! | dtf             | time_step_for_dynamics                                                                      | dynamics timestep                                                    | s             |    0 | real              | kind_phys | in     | F        |
!! | cld1d           | cloud_work_function                                                                         | cloud work function                                                  | m2 s-2        |    1 | real              | kind_phys | in     | F        |
!! | save_u          | x_wind_save                                                                                 | x-wind before entering a physics scheme                              | m s-1         |    2 | real              | kind_phys | in     | F        |
!! | save_v          | y_wind_save                                                                                 | y-wind before entering a physics scheme                              | m s-1         |    2 | real              | kind_phys | in     | F        |
!! | save_t          | air_temperature_save                                                                        | air temperature before entering a physics scheme                     | K             |    2 | real              | kind_phys | in     | F        |
!! | save_qv         | water_vapor_specific_humidity_save                                                          | water vapor specific humidity before entering a physics scheme       | kg kg-1       |    2 | real              | kind_phys | in     | F        |
!! | gu0             | x_wind_updated_by_physics                                                                   | zonal wind updated by physics                                        | m s-1         |    2 | real              | kind_phys | in     | F        |
!! | gv0             | y_wind_updated_by_physics                                                                   | meridional wind updated by physics                                   | m s-1         |    2 | real              | kind_phys | in     | F        |
!! | gt0             | air_temperature_updated_by_physics                                                          | temperature updated by physics                                       | K             |    2 | real              | kind_phys | in     | F        |
!! | gq0_water_vapor | water_vapor_specific_humidity_updated_by_physics                                            | water vapor specific humidity updated by physics                     | kg kg-1       |    2 | real              | kind_phys | in     | F        |
!! | ud_mf           | instantaneous_atmosphere_updraft_convective_mass_flux                                       | (updraft mass flux) * delt                                           | kg m-2        |    2 | real              | kind_phys | in     | F        |
!! | dd_mf           | instantaneous_atmosphere_downdraft_convective_mass_flux                                     | (downdraft mass flux) * delt                                         | kg m-2        |    2 | real              | kind_phys | in     | F        |
!! | dt_mf           | instantaneous_atmosphere_detrainment_convective_mass_flux                                   | (detrainment mass flux) * delt                                       | kg m-2        |    2 | real              | kind_phys | in     | F        |
!! | con_g           | gravitational_acceleration                                                                  | gravitational acceleration                                           | m s-2         |    0 | real              | kind_phys | in     | F        |
!! | clw_ice         | cloud_ice_mixing_ratio                                                                      | moist cloud ice mixing ratio                                         | kg kg-1       |    2 | real              | kind_phys | in     | F        |
!! | clw_liquid      | cloud_liquid_water_mixing_ratio                                                             | moist cloud water mixing ratio                                       | kg kg-1       |    2 | real              | kind_phys | in     | F        |
!! | npdf3d          | number_of_3d_arrays_associated_with_pdf-based_clouds                                        | number of 3d arrays associated with pdf based clouds/mp              | count         |    0 | integer           |           | in     | F        |
!! | num_p3d         | array_dimension_of_3d_arrays_for_microphysics                                               | number of 3D arrays needed for microphysics                          | count         |    0 | integer           |           | in     | F        |
!! | ncnvcld3d       | number_of_convective_3d_cloud_fields                                                        | number of convective 3d clouds fields                                | count         |    0 | integer           |           | in     | F        |
!! | rainc           | lwe_thickness_of_convective_precipitation_amount_on_dynamics_timestep                       | convective rain at this time step                                    | m             |    1 | real              | kind_phys | inout  | F        |
!! | cldwrk          | cumulative_cloud_work_function                                                              | cumulative cloud work function (valid only with sas)                 | m2 s-1        |    1 | real              | kind_phys | inout  | F        |
!! | cnvprcp         | cumulative_lwe_thickness_of_convective_precipitation_amount                                 | cumulative convective precipitation                                  | m             |    1 | real              | kind_phys | inout  | F        |
!! | cnvprcpb        | cumulative_lwe_thickness_of_convective_precipitation_amount_in_bucket                       | cumulative convective precipitation in bucket                        | m             |    1 | real              | kind_phys | inout  | F        |
!! | dt3dt           | cumulative_change_in_temperature_due_to_deep_convection                                     | cumulative change in temperature due to deep conv.                   | K             |    2 | real              | kind_phys | inout  | F        |
!! | dq3dt           | cumulative_change_in_water_vapor_specific_humidity_due_to_deep_convection                   | cumulative change in water vapor specific humidity due to deep conv. | kg kg-1       |    2 | real              | kind_phys | inout  | F        |
!! | du3dt           | cumulative_change_in_x_wind_due_to_deep_convection                                          | cumulative change in x wind due to deep convection                   | m s-1         |    2 | real              | kind_phys | inout  | F        |
!! | dv3dt           | cumulative_change_in_y_wind_due_to_deep_convection                                          | cumulative change in y wind due to deep convection                   | m s-1         |    2 | real              | kind_phys | inout  | F        |
!! | upd_mf          | cumulative_atmosphere_updraft_convective_mass_flux                                          | cumulative updraft mass flux                                         | Pa            |    2 | real              | kind_phys | inout  | F        |
!! | dwn_mf          | cumulative_atmosphere_downdraft_convective_mass_flux                                        | cumulative downdraft mass flux                                       | Pa            |    2 | real              | kind_phys | inout  | F        |
!! | det_mf          | cumulative_atmosphere_detrainment_convective_mass_flux                                      | cumulative detrainment mass flux                                     | Pa            |    2 | real              | kind_phys | inout  | F        |
!! | dqdti           | instantaneous_water_vapor_specific_humidity_tendency_due_to_convection                      | instantaneous moisture tendency due to convection                    | kg kg-1 s-1   |    2 | real              | kind_phys | inout  | F        |
!! | cnvqci          | instantaneous_deep_convective_cloud_condensate_mixing_ratio_on_dynamics_time_step           | instantaneous total convective condensate mixing ratio               | kg kg-1       |    2 | real              | kind_phys | inout  | F        |
!! | upd_mfi         | instantaneous_atmosphere_updraft_convective_mass_flux_on_dynamics_timestep                  | (updraft mass flux) * delt                                           | kg m-2        |    2 | real              | kind_phys | inout  | F        |
!! | dwn_mfi         | instantaneous_atmosphere_downdraft_convective_mass_flux_on_dynamics_timestep                | (downdraft mass flux) * delt                                         | kg m-2        |    2 | real              | kind_phys | inout  | F        |
!! | det_mfi         | instantaneous_atmosphere_detrainment_convective_mass_flux_on_dynamics_timestep              | (detrainment mass flux) * delt                                       | kg m-2        |    2 | real              | kind_phys | inout  | F        |
!! | cnvw            | convective_cloud_water_mixing_ratio                                                         | moist convective cloud water mixing ratio                            | kg kg-1       |    2 | real              | kind_phys | inout  | F        |
!! | cnvc            | convective_cloud_cover                                                                      | convective cloud cover                                               | frac          |    2 | real              | kind_phys | inout  | F        |
!! | cnvw_phy_f3d    | convective_cloud_water_mixing_ratio_in_phy_f3d                                              | convective cloud water mixing ratio in the phy_f3d array             | kg kg-1       |    2 | real              | kind_phys | inout  | F        |
!! | cnvc_phy_f3d    | convective_cloud_cover_in_phy_f3d                                                           | convective cloud cover in the phy_f3d array                          | frac          |    2 | real              | kind_phys | inout  | F        |
!! | errmsg          | ccpp_error_message                                                                          | error message for error handling in CCPP                             | none          |    0 | character         | len=*     | out    | F        |
!! | errflg          | ccpp_error_flag                                                                             | error flag for error handling in CCPP                                | flag          |    0 | integer           |           | out    | F        |
!!
    subroutine GFS_DCNV_generic_post_run (im, levs, lssav, ldiag3d, lgocart, frain, rain1, dtf, cld1d, &
      save_u, save_v, save_t, save_qv, gu0, gv0, gt0, gq0_water_vapor, ud_mf, dd_mf, dt_mf, con_g, &
      clw_ice, clw_liquid, npdf3d, num_p3d, ncnvcld3d, &
      rainc, cldwrk, cnvprcp, cnvprcpb, dt3dt, dq3dt, du3dt, dv3dt, upd_mf, dwn_mf, det_mf, dqdti, &
      cnvqci, upd_mfi, dwn_mfi, det_mfi, cnvw, cnvc, cnvw_phy_f3d, cnvc_phy_f3d, errmsg, errflg)

      use machine,               only: kind_phys

      implicit none

      integer, intent(in) :: im, levs
      logical, intent(in) :: lssav, ldiag3d, lgocart

      real(kind=kind_phys), intent(in) :: frain, dtf
      real(kind=kind_phys), dimension(im), intent(in) :: rain1, cld1d
      real(kind=kind_phys), dimension(im,levs), intent(in) :: save_u, save_v, save_t, save_qv
      real(kind=kind_phys), dimension(im,levs), intent(in) :: gu0, gv0, gt0, gq0_water_vapor
      real(kind=kind_phys), dimension(im,levs), intent(in) :: ud_mf, dd_mf, dt_mf
      real(kind=kind_phys), intent(in) :: con_g
      real(kind=kind_phys), dimension(im,levs), intent(in) :: clw_ice, clw_liquid
      integer, intent(in) :: npdf3d, num_p3d, ncnvcld3d

      real(kind=kind_phys), dimension(im), intent(inout) :: rainc, cldwrk, cnvprcp, cnvprcpb
      ! dt3dt, dq3dt, du3dt, dv3dt upd_mf, dwn_mf, det_mf only allocated if ldiag3d == .true.
      real(kind=kind_phys), dimension(:,:), intent(inout) :: dt3dt, dq3dt, du3dt, dv3dt
      real(kind=kind_phys), dimension(:,:), intent(inout) :: upd_mf, dwn_mf, det_mf
      ! dqdti, cnvqci, upd_mfi, dwn_mfi, det_mfi only allocated if ldiag3d == .true. or lgocart == .true.
      real(kind=kind_phys), dimension(:,:), intent(inout) :: dqdti, cnvqci, upd_mfi, dwn_mfi, det_mfi
      real(kind=kind_phys), dimension(im,levs), intent(inout) :: cnvw, cnvc, cnvw_phy_f3d, cnvc_phy_f3d

      character(len=*), intent(out) :: errmsg
      integer, intent(out) :: errflg

      integer :: i, k

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      do i=1,im
        rainc(i) = frain * rain1(i)
      enddo
!
      if (lssav) then
        do i=1,im
          cldwrk (i)  = cldwrk (i)  + cld1d(i) * dtf
          cnvprcp(i)  = cnvprcp(i)  + rainc(i)
          cnvprcpb(i) = cnvprcpb(i) + rainc(i)
        enddo

        if (ldiag3d) then
          do k=1,levs
            do i=1,im
              dt3dt(i,k) = dt3dt(i,k) + (gt0(i,k)-save_t(i,k)) * frain
              dq3dt(i,k) = dq3dt(i,k) + (gq0_water_vapor(i,k)-save_qv(i,k)) * frain
              du3dt(i,k) = du3dt(i,k) + (gu0(i,k)-save_u(i,k)) * frain
              dv3dt(i,k) = dv3dt(i,k) + (gv0(i,k)-save_v(i,k)) * frain

              upd_mf(i,k)  = upd_mf(i,k)  + ud_mf(i,k) * (con_g*frain)
              dwn_mf(i,k)  = dwn_mf(i,k)  + dd_mf(i,k) * (con_g*frain)
              det_mf(i,k)  = det_mf(i,k)  + dt_mf(i,k) * (con_g*frain)
            enddo
          enddo
        endif ! if (ldiag3d)

endif   ! end if_lssav

      !update dqdt_v to include moisture tendency due to deep convection
      if (lgocart) then
        do k=1,levs
          do i=1,im
            dqdti  (i,k) = (gq0_water_vapor(i,k)  - save_qv(i,k)) * frain
            upd_mfi(i,k) = upd_mfi(i,k) + ud_mf(i,k)   * frain
            dwn_mfi(i,k) = dwn_mfi(i,k) + dd_mf(i,k)   * frain
            det_mfi(i,k) = det_mfi(i,k) + dt_mf(i,k)   * frain
            cnvqci (i,k) = cnvqci (i,k) + (clw_ice(i,k)+clw_liquid(i,k))*frain
          enddo
        enddo
      endif ! if (lgocart)
!
      if ((npdf3d == 3) .and. (num_p3d == 4)) then
        do k=1,levs
          do i=1,im
            cnvw_phy_f3d(i,k) = cnvw(i,k)
            cnvc_phy_f3d(i,k) = cnvc(i,k)
            cnvw(i,k)             = 0.0
            cnvc(i,k)             = 0.0
          enddo
        enddo
      elseif ((npdf3d == 0) .and. (ncnvcld3d == 1)) then
        do k=1,levs
          do i=1,im
            cnvw_phy_f3d(i,k) = cnvw(i,k)
            cnvw(i,k)             = 0.0
          enddo
        enddo
      endif


    end subroutine GFS_DCNV_generic_post_run

    end module GFS_DCNV_generic_post

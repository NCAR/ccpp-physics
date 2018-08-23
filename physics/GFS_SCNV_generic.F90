!> \file GFS_SCNV_generic.F90
!!  Contains code related to shallow convective schemes to be used within the GFS physics suite.

      module GFS_SCNV_generic_pre

      contains

      subroutine GFS_SCNV_generic_pre_init ()
      end subroutine GFS_SCNV_generic_pre_init

      subroutine GFS_SCNV_generic_pre_finalize()
      end subroutine GFS_SCNV_generic_pre_finalize

!> \section arg_table_GFS_SCNV_generic_pre_run Argument Table
!! | local_name      | standard_name                                          | long_name                                                             | units         | rank | type                  |    kind   | intent | optional |
!! |-----------------|--------------------------------------------------------|-----------------------------------------------------------------------|---------------|------|-----------------------|-----------|--------|----------|
!! | im              | horizontal_loop_extent                                 | horizontal loop extent                                                | count         |    0 | integer               |           | in     | F        |
!! | levs            | vertical_dimension                                     | vertical layer dimension                                              | count         |    0 | integer               |           | in     | F        |
!! | ldiag3d         | flag_diagnostics_3D                                    | flag for 3d diagnostic fields                                         | flag          |    0 | logical               |           | in     | F        |
!! | lgocart         | flag_gocart                                            | flag for 3d diagnostic fields for gocart 1                            | flag          |    0 | logical               |           | in     | F        |
!! | gt0             | air_temperature_updated_by_physics                     | temperature updated by physics                                        | K             |    2 | real                  | kind_phys | in     | F        |
!! | gq0_water_vapor | water_vapor_specific_humidity_updated_by_physics       | water vapor specific humidity updated by physics                      | kg kg-1       |    2 | real                  | kind_phys | in     | F        |
!! | save_t          | air_temperature_save                                   | air temperature before entering a physics scheme                      | K             |    2 | real                  | kind_phys | inout  | F        |
!! | save_qv         | water_vapor_specific_humidity_save                     | water vapor specific humidity before entering a physics scheme        | kg kg-1       |    2 | real                  | kind_phys | inout  | F        |
!! | errmsg          | ccpp_error_message                                     | error message for error handling in CCPP                              | none          |    0 | character             | len=*     | out    | F        |
!! | errflg          | ccpp_error_flag                                        | error flag for error handling in CCPP                                 | flag          |    0 | integer               |           | out    | F        |
!!
      subroutine GFS_SCNV_generic_pre_run (im, levs, ldiag3d, lgocart, gt0, gq0_water_vapor, &
        save_t, save_qv, errmsg, errflg)

        use machine,               only: kind_phys

        implicit none

        integer, intent(in) :: im, levs
        logical, intent(in) :: ldiag3d, lgocart
        real(kind=kind_phys), dimension(im,levs), intent(in) :: gt0, gq0_water_vapor

        real(kind=kind_phys), dimension(im,levs), intent(inout) :: save_t, save_qv
        character(len=*),                 intent(out) :: errmsg
        integer,                          intent(out) :: errflg

        integer :: i, k

        ! Initialize CCPP error handling variables
        errmsg = ''
        errflg = 0

        if (ldiag3d) then
          do k=1,levs
            do i=1,im
              save_t(i,k)   = gt0(i,k)
            enddo
          enddo
        endif
        if (ldiag3d .or. lgocart) then
          do k=1,levs
            do i=1,im
              save_qv(i,k) = gq0_water_vapor(i,k)
            enddo
          enddo
        endif

    end subroutine GFS_SCNV_generic_pre_run

    end module GFS_SCNV_generic_pre

    module GFS_SCNV_generic_post

    contains

    subroutine GFS_SCNV_generic_post_init ()
    end subroutine GFS_SCNV_generic_post_init

    subroutine GFS_SCNV_generic_post_finalize ()
    end subroutine GFS_SCNV_generic_post_finalize

!> \section arg_table_GFS_SCNV_generic_post_run Argument Table
!! | local_name      | standard_name                                                                               | long_name                                                            | units         | rank | type        |    kind   | intent | optional |
!! |-----------------|---------------------------------------------------------------------------------------------|----------------------------------------------------------------------|---------------|------|-------------|-----------|--------|----------|
!! | im              | horizontal_loop_extent                                                                      | horizontal loop extent                                               | count         |    0 | integer     |           | in     | F        |
!! | levs            | vertical_dimension                                                                          | vertical layer dimension                                             | count         |    0 | integer     |           | in     | F        |
!! | nn              | number_of_tracers_for_convective_transport                                                  | number of tracers for convective transport                           | count         |    0 | integer     |           | in     | F        |
!! | lssav           | flag_diagnostics                                                                            | logical flag for storing diagnostics                                 | flag          |    0 | logical     |           | in     | F        |
!! | ldiag3d         | flag_diagnostics_3D                                                                         | flag for 3d diagnostic fields                                        | flag          |    0 | logical     |           | in     | F        |
!! | lgocart         | flag_gocart                                                                                 | flag for 3d diagnostic fields for gocart 1                           | flag          |    0 | logical     |           | in     | F        |
!! | frain           | dynamics_to_physics_timestep_ratio                                                          | ratio of dynamics timestep to physics timestep                       | none          |    0 | real        | kind_phys | in     | F        |
!! | gt0             | air_temperature_updated_by_physics                                                          | temperature updated by physics                                       | K             |    2 | real        | kind_phys | in     | F        |
!! | gq0_water_vapor | water_vapor_specific_humidity_updated_by_physics                                            | water vapor specific humidity updated by physics                     | kg kg-1       |    2 | real        | kind_phys | in     | F        |
!! | save_t          | air_temperature_save                                                                        | air temperature before entering a physics scheme                     | K             |    2 | real        | kind_phys | in     | F        |
!! | save_qv         | water_vapor_specific_humidity_save                                                          | water vapor specific humidity before entering a physics scheme       | kg kg-1       |    2 | real        | kind_phys | in     | F        |
!! | dqdti           | instantaneous_water_vapor_specific_humidity_tendency_due_to_convection                      | instantaneous moisture tendency due to convection                    | kg kg-1 s-1   |    2 | real        | kind_phys | inout  | F        |
!! | dt3dt           | cumulative_change_in_temperature_due_to_shal_convection                                     | cumulative change in temperature due to shal conv.                   | K             |    2 | real        | kind_phys | inout  | F        |
!! | dq3dt           | cumulative_change_in_water_vapor_specific_humidity_due_to_shal_convection                   | cumulative change in water vapor specific humidity due to shal conv. | kg kg-1       |    2 | real        | kind_phys | inout  | F        |
!! | clw             | convective_transportable_tracers                                                            | array to contain cloud water and other convective trans. tracers     | kg kg-1       |    3 | real        | kind_phys | inout  | F        |
!! | errmsg          | ccpp_error_message                                                                          | error message for error handling in CCPP                             | none          |    0 | character   | len=*     | out    | F        |
!! | errflg          | ccpp_error_flag                                                                             | error flag for error handling in CCPP                                | flag          |    0 | integer     |           | out    | F        |
!!
      subroutine GFS_SCNV_generic_post_run (im, levs, nn, lssav, ldiag3d, lgocart, frain, gt0, gq0_water_vapor, &
        save_t, save_qv, dqdti, dt3dt, dq3dt, clw, errmsg, errflg)

      use machine,               only: kind_phys

      implicit none

      integer, intent(in) :: im, levs, nn
      logical, intent(in) :: lssav, ldiag3d, lgocart
      real(kind=kind_phys),                     intent(in) :: frain
      real(kind=kind_phys), dimension(im,levs), intent(in) :: gt0, gq0_water_vapor
      real(kind=kind_phys), dimension(im,levs), intent(in) :: save_t, save_qv

      ! dqdti only allocated if ldiag3d == .true. or lgocart == .true.
      real(kind=kind_phys), dimension(:,:), intent(inout) :: dqdti
      ! dt3dt, dq3dt, only allocated if ldiag3d == .true.
      real(kind=kind_phys), dimension(:,:), intent(inout) :: dt3dt, dq3dt
      real(kind=kind_phys), dimension(im,levs,nn), intent(inout) :: clw

      character(len=*),              intent(out) :: errmsg
      integer,                       intent(out) :: errflg

      integer :: i, k
      real(kind=kind_phys) :: tem

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      if (lssav) then
!          update dqdt_v to include moisture tendency due to shallow convection
        if (lgocart) then
          do k=1,levs
            do i=1,im
              tem  = (gq0_water_vapor(i,k)-save_qv(i,k)) * frain
              dqdti(i,k) = dqdti(i,k)  + tem
            enddo
          enddo
        endif
        if (ldiag3d) then
          do k=1,levs
            do i=1,im
              dt3dt(i,k) = dt3dt(i,k) + (gt0(i,k)  - save_t(i,k))   * frain
              dq3dt(i,k) = dq3dt(i,k) + (gq0_water_vapor(i,k) - save_qv(i,k)) * frain
            enddo
          enddo
        endif
      endif   ! end if_lssav
!
      do k=1,levs
        do i=1,im
          if (clw(i,k,2) <= -999.0) clw(i,k,2) = 0.0
        enddo
      enddo

      end subroutine GFS_SCNV_generic_post_run

      end module GFS_SCNV_generic_post

!> \file GFS_DCNV_generic.f90
!!  Contains code related to deep convective schemes to be used within the GFS physics suite.

      module GFS_DCNV_generic_pre

      contains

      subroutine GFS_DCNV_generic_pre_init ()
      end subroutine GFS_DCNV_generic_pre_init

      subroutine GFS_DCNV_generic_pre_finalize()
      end subroutine GFS_DCNV_generic_pre_finalize

!> \section arg_table_GFS_DCNV_generic_pre_run Argument Table
!! | local var name | longname                                               | description                                                           | units         | rank | type                          |    kind   | intent | optional |
!! |----------------|--------------------------------------------------------|-----------------------------------------------------------------------|---------------|------|-------------------------------|-----------|--------|----------|
!! | Model          | FV3-GFS_Control_type                                   | Fortran DDT containing FV3-GFS model control parameters               | DDT           |    0 | GFS_typedefs%GFS_control_type |           | in     | F        |
!! | Stateout       | FV3-GFS_Stateout_type                                  | Fortran DDT containing FV3-GFS prognostic state to return to dycore   | DDT           |    0 | GFS_typedefs%GFS_stateout_type|           | in     | F        |
!! | Grid           | FV3-GFS_Grid_type                                      | Fortran DDT containing FV3-GFS grid and interpolation related data    | DDT           |    0 | GFS_typedefs%GFS_grid_type    |           | in     | F        |
!! | initial_u      | x_wind_initial                                         | x-wind before entering a physics scheme                               | m s-1         |    2 | real                          | kind_phys | inout  | F        |
!! | initial_v      | y_wind_initial                                         | y-wind before entering a physics scheme                               | m s-1         |    2 | real                          | kind_phys | inout  | F        |
!! | initial_t      | air_temperature_initial                                | air temperature before entering a physics scheme                      | K             |    2 | real                          | kind_phys | inout  | F        |
!! | initial_qv     | water_vapor_specific_humidity_initial                  | water vapor specific humidity before entering a physics scheme        | kg kg-1       |    2 | real                          | kind_phys | inout  | F        |
!!
      subroutine GFS_DCNV_generic_pre_run (Model, Stateout, Grid, initial_u, initial_v, initial_t, initial_qv)
        use machine,               only: kind_phys
        use GFS_typedefs,          only: GFS_control_type, GFS_stateout_type, GFS_grid_type

        type(GFS_control_type),           intent(in) :: Model
        type(GFS_stateout_type),          intent(in) :: Stateout
        type(GFS_grid_type),              intent(in) :: Grid
        real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levs), intent(inout) :: initial_u, initial_v, initial_t, initial_qv

        if (Model%ldiag3d) then
          initial_t(:,:) = Stateout%gt0(:,:)
          initial_u(:,:) = Stateout%gu0(:,:)
          initial_v(:,:) = Stateout%gv0(:,:)
        elseif (Model%cnvgwd) then
          initial_t(:,:) = Stateout%gt0(:,:)
        endif   ! end if_ldiag3d/cnvgwd

        if (Model%ldiag3d .or. Model%lgocart) then
          initial_qv(:,:) = Stateout%gq0(:,:,1)
        endif   ! end if_ldiag3d/lgocart

    end subroutine GFS_DCNV_generic_pre_run

    end module

    module GFS_DCNV_generic_post

    contains

    subroutine GFS_DCNV_generic_post_init ()
    end subroutine GFS_DCNV_generic_post_init

    subroutine GFS_DCNV_generic_post_finalize ()
    end subroutine GFS_DCNV_generic_post_finalize

!> \section arg_table_GFS_PBL_generic_post_run Argument Table
!! | local var name | longname                                                  | description                                                           | units         | rank | type                          |    kind   | intent | optional |
!! |----------------|-----------------------------------------------------------|-----------------------------------------------------------------------|---------------|------|-------------------------------|-----------|--------|----------|
!! | Grid           | FV3-GFS_Grid_type                                         | Fortran DDT containing FV3-GFS grid and interpolation related data    | DDT           |    0 | GFS_typedefs%GFS_grid_type    |           | in     | F        |
!! | Model          | FV3-GFS_Control_type                                      | Fortran DDT containing FV3-GFS model control parameters               | DDT           |    0 | GFS_typedefs%GFS_control_type |           | in     | F        |
!! | frain          | dynamics_to_physics_timestep_ratio                        | ratio of dynamics timestep to physics timestep                        | none          |    0 | real                          | kind_phys | in     | F        |
!! | rain1          | instantaneous_rainfall_amount                             | instantaneous rainfall amount                                         | m             |    1 | real                          | kind_phys | in     | F        |
!! | cld1d          | cloud_work_function                                       | cloud work function                                                   | m2 s-2        |    1 | real                          | kind_phys | in     | F        |
!! | initial_u      | x_wind_initial                                            | x-wind before entering a physics scheme                               | m s-1         |    2 | real                          | kind_phys | in     | F        |
!! | initial_v      | y_wind_initial                                            | y-wind before entering a physics scheme                               | m s-1         |    2 | real                          | kind_phys | in     | F        |
!! | initial_t      | air_temperature_initial                                   | air temperature before entering a physics scheme                      | K             |    2 | real                          | kind_phys | in     | F        |
!! | initial_qv     | water_vapor_specific_humidity_initial                     | water vapor specific humidity before entering a physics scheme        | kg kg-1       |    2 | real                          | kind_phys | in     | F        |
!! | ud_mf          | instantaneous_atmosphere_updraft_convective_mass_flux     | (updraft mass flux) * delt                                            | kg m-2        |    2 | real                          | kind_phys | in     | F        |
!! | dd_mf          | instantaneous_atmosphere_downdraft_convective_mass_flux   | (downdraft mass flux) * delt                                          | kg m-2        |    2 | real                          | kind_phys | in     | F        |
!! | dt_mf          | instantaneous_atmosphere_detrainment_convective_mass_flux | (detrainment mass flux) * delt                                        | kg m-2        |    2 | real                          | kind_phys | in     | F        |
!! | Diag           | FV3-GFS_Diag_type                                         | Fortran DDT containing FV3-GFS fields targeted for diagnostic output  | DDT           |    0 | GFS_typedefs%GFS_diag_type    |           | inout  | F        |
!! | Tbd            | FV3-GFS_Tbd_type                                          | Fortran DDT containing FV3-GFS miscellaneous data                     | DDT           |    0 | GFS_typedefs%GFS_tbd_type     |           | inout  | F        |
!!
      subroutine GFS_DCNV_generic_post_run (Grid, Model, Stateout, frain, rain1, cld1d, initial_u, initial_v, initial_t, initial_qv, &
        ud_mf, dd_mf, dt_mf, cnvw, cnvc, Diag, Tbd)

      use machine,               only: kind_phys
      use GFS_typedefs,          only: GFS_grid_type, GFS_control_type, GFS_stateout_type, GFS_diag_type, GFS_tbd_type
      use physcons,              only: con_g

      type(GFS_grid_type),            intent(in) :: Grid
      type(GFS_control_type),         intent(in) :: Model
      type(GFS_stateout_type),        intent(in) :: Stateout
      type(GFS_diag_type),         intent(inout) :: Diag
      type(GFS_tbd_type),          intent(inout) :: Tbd

      real(kind=kind_phys), intent(in) :: frain
      real(kind=kind_phys), dimension(size(Grid%xlon,1)), intent(in) :: rain1, cld1d
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levs), intent(in) :: initial_u, initial_v, initial_t, initial_qv
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levs), intent(in) :: ud_mf, dd_mf, dt_mf
      real(kind=kind_phys), dimension(size(Grid%xlon,1),Model%levs), intent(in) :: cnvw, cnvc

      integer :: i, num2, num3

        do i = 1, size(Grid%xlon,1)
          Diag%rainc(:) = frain * rain1(:)
        enddo
  !
        if (Model%lssav) then
          Diag%cldwrk (:) = Diag%cldwrk (:) + cld1d(:) * Model%dtf
          Diag%cnvprcp(:) = Diag%cnvprcp(:) + Diag%rainc(:)

          if (Model%ldiag3d) then
            Diag%dt3dt(:,:,4) = Diag%dt3dt(:,:,4) + (Stateout%gt0(:,:)-initial_t(:,:)) * frain
            Diag%dq3dt(:,:,2) = Diag%dq3dt(:,:,2) + (Stateout%gq0(:,:,1)-initial_qv(:,:)) * frain
            Diag%du3dt(:,:,3) = Diag%du3dt(:,:,3) + (Stateout%gu0(:,:)-initial_u(:,:)) * frain
            Diag%dv3dt(:,:,3) = Diag%dv3dt(:,:,3) + (Stateout%gv0(:,:)-initial_v(:,:)) * frain

            Diag%upd_mf(:,:)  = Diag%upd_mf(:,:)  + ud_mf(:,:) * (con_g*frain)
            Diag%dwn_mf(:,:)  = Diag%dwn_mf(:,:)  + dd_mf(:,:) * (con_g*frain)
            Diag%det_mf(:,:)  = Diag%det_mf(:,:)  + dt_mf(:,:) * (con_g*frain)
          endif ! if (ldiag3d)

        endif   ! end if_lssav

        if ((Model%npdf3d == 3) .and. (Model%num_p3d == 4)) then
          num2 = Model%num_p3d + 2
          num3 = num2 + 1
          Tbd%phy_f3d(:,:,num2) = cnvw(:,:)
          Tbd%phy_f3d(:,:,num3) = cnvc(:,:)
        elseif ((Model%npdf3d == 0) .and. (Model%ncnvcld3d == 1)) then
          num2 = Model%num_p3d + 1
          Tbd%phy_f3d(:,:,num2) = cnvw(:,:)
        endif

      end subroutine GFS_DCNV_generic_post_run

      end module

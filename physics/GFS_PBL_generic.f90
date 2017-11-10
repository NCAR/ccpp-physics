!> \file GFS_PBL_generic.f90
!!  Contains code related to PBL schemes to be used within the GFS physics suite.

      module GFS_PBL_generic_pre

      contains

      subroutine GFS_PBL_generic_pre_init ()
      end subroutine GFS_PBL_generic_pre_init

      subroutine GFS_PBL_generic_pre_finalize()
      subroutine GFS_PBL_generic_pre_finalize

!> \section arg_table_GFS_PBL_generic_pre_run Argument Table
!! | local var name | longname                                               | description                                        | units         | rank | type    |    kind   | intent | optional |
!! |----------------|--------------------------------------------------------|----------------------------------------------------|---------------|------|---------|-----------|--------|----------|
!! | im             | horizontal_loop_extent                                 | horizontal loop extent, start at 1                 | index         |    0 | integer |           | in     | F        |
!! | levs           | vertical_dimension                                     | vertical layer dimension                           | index         |    0 | integer |           | in     | F        |
!! | kinver         | index_of_highest_temperature_inversion                 | index of highest temperature inversion             | index         |    1 | integer |           | in     | F        |
!!
      subroutine GFS_PBL_generic_pre_run (im, levs, kinver)

      integer , intent(in) :: im, levs
      integer, dimension(im), intent(inout) :: kinver

      kinver(:) = levs

    end subroutine GFS_PBL_generic_pre_run

    end module

    module GFS_PBL_generic_post

    contains

    subroutine GFS_PBL_generic_post_init ()
    end subroutine GFS_PBL_generic_post_init

    subroutine GFS_PBL_generic_post_finalize ()
    end subroutine GFS_PBL_generic_post_finalize



!> \section arg_table_GFS_PBL_generic_post_run Argument Table
!! | local var name | longname                                               | description                                                           | units         | rank | type                          |    kind   | intent | optional |
!! |----------------|--------------------------------------------------------|-----------------------------------------------------------------------|---------------|------|-------------------------------|-----------|--------|----------|
!! | Grid           | FV3-GFS_Grid_type                                      | Fortran DDT containing FV3-GFS grid and interpolation related data    | DDT           |    0 | GFS_typedefs%GFS_grid_type    |           | in     | F        |
!! | Model          | FV3-GFS_Control_type                                   | Fortran DDT containing FV3-GFS model control parameters               | DDT           |    0 | GFS_typedefs%GFS_control_type |           | in     | F        |
!! | Radtend        | FV3-GFS_Radtend_type                                   | Fortran DDT containing FV3-GFS radiation tendencies needed in physics | DDT           |    0 | GFS_typedefs%GFS_radtend_type |           | in     | F        |
!! | dusfc1         | instantaneous_surface_x_momentum_flux                  | surface momentum flux in the x-direction valid for current call       | Pa            |    1 | real                          | kind_phys | in     | F        |
!! | dvsfc1         | instantaneous_surface_y_momentum_flux                  | surface momentum flux in the y-direction valid for current call       | Pa            |    1 | real                          | kind_phys | in     | F        |
!! | dtsfc1         | surface_upward_sensible_heat_flux                      | surface upward sensible heat flux valid for current call              | W m-2         |    1 | real                          | kind_phys | in     | F        |
!! | dqsfc1         | surface_upward_latent_heat_flux                        | surface upward latent heat flux valid for current call                | W m-2         |    1 | real                          | kind_phys | in     | F        |
!! | dudt           | tendency_of_x_wind_due_to_model_physics                | updated tendency of the x wind                                        | m s-2         |    2 | real                          | kind_phys | in     | F        |
!! | dvdt           | tendency_of_y_wind_due_to_model_physics                | updated tendency of the y wind                                        | m s-2         |    2 | real                          | kind_phys | in     | F        |
!! | dtdt           | tendency_of_air_temperature_due_to_model_physics       | updated tendency of the temperature                                   | K s-1         |    2 | real                          | kind_phys | in     | F        |
!! | dqdt           | tendency_of_tracers_due_to_model_physics               | updated tendency of the tracers                                       | kg kg-1 s-1   |    3 | real                          | kind_phys | in     | F        |
!! | xmu            | time_step_zenith_angle_adjust_factor_for_sw            | time step zenith angle adjust factor for shortwave                    | none          |    2 | real                          | kind_phys | in     | F        |
!! | Diag           | FV3-GFS_diag_type                                      | Fortran DDT containing FV3-GFS fields targeted for diagnostic output  | DDT           |    0 | GFS_typedefs%GFS_diag_type    |           | in     | F        |
!!
      subroutine GFS_PBL_generic_post_run (Grid, Model, Radtend, dusfc1, dvsfc1, dtsfc1, dqsfc1, &
        dudt, dvdt, dtdt, dqdt, xmu, Diag)

      use machine,               only: kind_phys
      use GFS_typedefs,          only: GFS_diag_type, GFS_radtend_type, GFS_model_type, GFS_grid_type

      type(GFS_grid_type),            intent(in) :: Grid
      type(GFS_radtend_type),         intent(in) :: Radtend
      type(GFS_control_type),           intent(in) :: Model
      type(GFS_diag_type),         intent(inout) :: Diag

      real(kind=kind_phys), dimension(size(Grid%xlon,1)), intent(in) :: dusfc1, dvsfc1, dtsfc1, dqsfc1, xmu
      real(kind=kind_phys), dimension(size(Grid%xlon,1), Model%levs), intent(in) :: dudt, dvdt, dtdt
      real(kind=kind_phys), dimension(size(Grid%xlon,1), Model%levs, Model%ntrac), intent(in) :: dqdt

      integer :: i, k
      real(kind=kind_phys) :: tem

      if (Model%lssav) then
          Diag%dusfc (:) = Diag%dusfc(:) + dusfc1(:)*Model%dtf
          Diag%dvsfc (:) = Diag%dvsfc(:) + dvsfc1(:)*Model%dtf
          Diag%dtsfc (:) = Diag%dtsfc(:) + dtsfc1(:)*Model%dtf
          Diag%dqsfc (:) = Diag%dqsfc(:) + dqsfc1(:)*Model%dtf
          Diag%dusfci(:) = dusfc1(:)
          Diag%dvsfci(:) = dvsfc1(:)
          Diag%dtsfci(:) = dtsfc1(:)
          Diag%dqsfci(:) = dqsfc1(:)
  !       if (lprnt) then
  !         write(0,*)' dusfc=',dusfc(ipr),' dusfc1=',dusfc1(ipr),' dtf=',
  !    &     dtf,' kdt=',kdt,' lat=',lat
  !       endif

          if (Model%ldiag3d) then
            do k = 1, Model%levs
              do i = 1, size(Grid%xlon,1)
                tem          = dtdt(i,k) - (Radtend%htrlw(i,k)+Radtend%htrsw(i,k)*xmu(i))
                Diag%dt3dt(i,k,3) = Diag%dt3dt(i,k,3) + tem*Model%dtf
              enddo
            enddo
            Diag%du3dt(:,:,1) = Diag%du3dt(:,:,1) + dudt(:,:) * Model%dtf
            Diag%du3dt(:,:,2) = Diag%du3dt(:,:,2) - dudt(:,:) * Model%dtf
            Diag%dv3dt(:,:,1) = Diag%dv3dt(:,:,1) + dvdt(:,:) * Model%dtf
            Diag%dv3dt(:,:,2) = Diag%dv3dt(:,:,2) - dvdt(:,:) * Model%dtf
  ! update dqdt_v to include moisture tendency due to vertical diffusion
  !         if (lgocart) then
  !           do k = 1, levs
  !             do i = 1, im
  !               dqdt_v(i,k)  = dqdt(i,k,1) * dtf
  !             enddo
  !           enddo
  !         endif
            do k = 1, Model%levs
              do i = 1, size(Grid%xlon,1)
                tem  = dqdt(i,k,1) * Model%dtf
                Diag%dq3dt(i,k,1) = Diag%dq3dt(i,k,1) + tem
              enddo
            enddo
            if (Model%ntoz > 0) then
              Diag%dq3dt(:,:,5) = Diag%dq3dt(:,:,5) + dqdt(i,k,Model%ntoz) * Model%dtf
            endif
          endif

        endif   ! end if_lssav
      end subroutine GFS_PBL_generic_post_run

      end module

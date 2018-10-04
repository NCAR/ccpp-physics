!> \file cu_gf_driver_post.F90
!!  Contains code related to GF convective schemes to be used within the GFS physics suite.

      module cu_gf_driver_post

      contains

      subroutine cu_gf_driver_post_init ()
      end subroutine cu_gf_driver_post_init

      subroutine cu_gf_driver_post_finalize()
      end subroutine cu_gf_driver_post_finalize

!> \section arg_table_cu_gf_driver_post_run Argument Table
!! | local_name     | standard_name                                          | long_name                                                                | units         | rank | type              |    kind   | intent | optional |
!! |----------------|--------------------------------------------------------|--------------------------------------------------------------------------|---------------|------|-------------------|-----------|--------|----------|
!! | Model          | FV3-GFS_Control_type                                   | Fortran DDT containing FV3-GFS model control parameters                  | DDT           |    0 | GFS_control_type  |           | in     | F        |
!! | Stateout       | FV3-GFS_Stateout_type                                  |Fortran DDT containing FV3-GFS prognostic state to return to dycore       | DDT           |    0 | GFS_stateout_type |           | inout  | F        |
!! | Grid           | FV3-GFS_Grid_type                                      | Fortran DDT containing FV3-GFS grid and interpolation related data       | DDT           |    0 | GFS_grid_type     |           | in     | F        |
!! | Tbd            | FV3-GFS_Tbd_type                                       | Fortran DDT containing FV3-GFS miscellaneous data                        | DDT           |    0 | GFS_tbd_type      |           | inout  | F        |
!! | errmsg         | error_message                                          | error message for error handling in CCPP                                 | none          |    0 | character         | len=*     | out    | F        |
!! | errflg         | error_flag                                             | error flag for error handling in CCPP                                    | flag          |    0 | integer           |           | out    | F        |
!!
    subroutine cu_gf_driver_post_run (Model, Stateout, Grid, Tbd, errmsg, errflg)

      use machine,               only: kind_phys
      use GFS_typedefs,          only: GFS_control_type, GFS_stateout_type, GFS_grid_type, GFS_tbd_type

      implicit none

      type(GFS_control_type),           intent(in) :: Model
      type(GFS_stateout_type),       intent(inout) :: Stateout
      type(GFS_grid_type),              intent(in) :: Grid
      type(GFS_tbd_type),            intent(inout) :: Tbd

      character(len=*), intent(out) :: errmsg
      integer, intent(out) :: errflg

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      Tbd%prevst(:,:) = Stateout%gt0(:,:)
      Tbd%prevsq(:,:) = Stateout%gq0(:,:,1)

    end subroutine cu_gf_driver_post_run

    end module cu_gf_driver_post

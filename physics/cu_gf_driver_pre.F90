!> \file cu_gf_driver_pre.F90
!!  Contains code related to GF convective schemes to be used within the GFS physics suite.

      module cu_gf_driver_pre

      contains

      subroutine cu_gf_driver_pre_init ()
      end subroutine cu_gf_driver_pre_init

      subroutine cu_gf_driver_pre_finalize()
      end subroutine cu_gf_driver_pre_finalize

!> \section arg_table_cu_gf_driver_pre_run Argument Table
!! | local_name     | standard_name                                          | long_name                                                                | units         | rank | type              |    kind   | intent | optional |
!! |----------------|--------------------------------------------------------|--------------------------------------------------------------------------|---------------|------|-------------------|-----------|--------|----------|
!! | Model          | FV3-GFS_Control_type                                   | Fortran DDT containing FV3-GFS model control parameters                  | DDT           |    0 | GFS_control_type  |           | in     | F        |
!! | Statein        | FV3-GFS_Statein_type                                   | Fortran DDT containing FV3-GFS prognostic state in from dycore           | DDT           |    0 | GFS_statein_type  |           | in     | F        |
!! | Grid           | FV3-GFS_Grid_type                                      | Fortran DDT containing FV3-GFS grid and interpolation related data       | DDT           |    0 | GFS_grid_type     |           | in     | F        |
!! | Tbd            | FV3-GFS_Tbd_type                                       | Fortran DDT containing FV3-GFS miscellaneous data                        | DDT           |    0 | GFS_tbd_type      |           | inout  | F        |
!! | errmsg         | error_message                                          | error message for error handling in CCPP                                 | none          |    0 | character         | len=*     | out    | F        |
!! | errflg         | error_flag                                             | error flag for error handling in CCPP                                    | flag          |    0 | integer           |           | out    | F        |
!!
    subroutine cu_gf_driver_pre_run (Model, Statein, Grid, Tbd, errmsg, errflg)

      use machine,               only: kind_phys
      use GFS_typedefs,          only: GFS_control_type, GFS_statein_type, GFS_grid_type, GFS_tbd_type

      implicit none

      type(GFS_control_type),           intent(in) :: Model
      type(GFS_statein_type),           intent(in) :: Statein
      type(GFS_grid_type),              intent(in) :: Grid
      type(GFS_tbd_type),            intent(inout) :: Tbd

      !--- hli for GF convective scheme
      real(kind=kind_phys) :: dtdyn
      !--
      character(len=*), intent(out) :: errmsg
      integer, intent(out) :: errflg

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      if(Model%kdt.gt.1) then
       dtdyn=3600.0*(Model%fhour)/Model%kdt   ! G. Firl's suggestion
       if(Model%dtp > dtdyn) then
        Tbd%forcet(:,:)=(Statein%tgrs(:,:)  - Tbd%prevst(:,:))/Model%dtp
        Tbd%forceq(:,:)=(Statein%qgrs(:,:,1)- Tbd%prevsq(:,:))/Model%dtp
       else
        Tbd%forcet(:,:)=(Statein%tgrs(:,:)  - Tbd%prevst(:,:))/dtdyn
        Tbd%forceq(:,:)=(Statein%qgrs(:,:,1)- Tbd%prevsq(:,:))/dtdyn
       endif
      else
       Tbd%forcet(:,:)=0.0
       Tbd%forceq(:,:)=0.0
      endif

    end subroutine cu_gf_driver_pre_run

    end module cu_gf_driver_pre

!> \file GFS_suite_stateout_reset.f90
!!  Contains code to set the values of the physics-updated state to the before-physics state prior to actually being modified by physics.

  module GFS_suite_stateout_reset

  contains

!> \section arg_table_GFS_suite_stateout_reset_run Argument Table
!! \htmlinclude GFS_suite_stateout_reset_run.html
!!
    subroutine GFS_suite_stateout_reset_run (im, levs, ntrac,        &
                                             tgrs, ugrs, vgrs, qgrs, &
                                             gt0 , gu0 , gv0 , gq0 , &
                                             errmsg, errflg)

      use machine,               only: kind_phys

      implicit none

      ! interface variables
      integer,              intent(in )                   :: im
      integer,              intent(in )                   :: levs
      integer,              intent(in )                   :: ntrac
      real(kind=kind_phys), intent(in ), dimension(:,:)   :: tgrs, ugrs, vgrs
      real(kind=kind_phys), intent(in ), dimension(:,:,:) :: qgrs
      real(kind=kind_phys), intent(out), dimension(:,:)   :: gt0, gu0, gv0
      real(kind=kind_phys), intent(out), dimension(:,:,:) :: gq0

      character(len=*),     intent(out)                   :: errmsg
      integer,              intent(out)                   :: errflg

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      gt0(:,:)   = tgrs(:,:)
      gu0(:,:)   = ugrs(:,:)
      gv0(:,:)   = vgrs(:,:)
      gq0(:,:,:) = qgrs(:,:,:)

    end subroutine GFS_suite_stateout_reset_run

  end module GFS_suite_stateout_reset
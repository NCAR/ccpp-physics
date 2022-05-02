!> \file GFS_suite_interstitial_rad_reset.f90
!!  Contains code to reset radiation-related interstitial variables

    module GFS_suite_interstitial_rad_reset

    contains

!> \section arg_table_GFS_suite_interstitial_rad_reset_run Argument Table
!! \htmlinclude GFS_suite_interstitial_rad_reset_run.html
!!
    subroutine GFS_suite_interstitial_rad_reset_run (Interstitial, Model, errmsg, errflg)

      use machine,       only: kind_phys
      use GFS_typedefs,  only: GFS_control_type
      use CCPP_typedefs, only: GFS_interstitial_type

      implicit none

      ! interface variables
      type(GFS_interstitial_type), intent(inout) :: Interstitial
      type(GFS_control_type),      intent(in)    :: Model
      character(len=*),            intent(out)   :: errmsg
      integer,                     intent(out)   :: errflg

      errmsg = ''
      errflg = 0

      call Interstitial%rad_reset(Model)

    end subroutine GFS_suite_interstitial_rad_reset_run

    end module GFS_suite_interstitial_rad_reset
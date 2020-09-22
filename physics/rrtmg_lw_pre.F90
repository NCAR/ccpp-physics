!>\file rrtmg_lw_pre.f90
!! This file contains a call to module_radiation_surface::setemis() to
!! setup surface emissivity for LW radiation.
      module rrtmg_lw_pre
      contains

!>\defgroup rrtmg_lw_pre GFS RRTMG scheme pre
!! @{
!> \section arg_table_rrtmg_lw_pre_init Argument Table
!!
      subroutine rrtmg_lw_pre_init ()
      end subroutine rrtmg_lw_pre_init 

!> \section arg_table_rrtmg_lw_pre_run Argument Table
!! \htmlinclude rrtmg_lw_pre_run.html
!!
      subroutine rrtmg_lw_pre_run (Model, Grid, Sfcprop, Radtend, im, tsfg, tsfa, errmsg, errflg)
    
      use machine,                   only: kind_phys

      use GFS_typedefs,              only: GFS_control_type,           &
                                           GFS_grid_type,              &
                                           GFS_radtend_type,           &
                                           GFS_sfcprop_type         
      use module_radiation_surface,  only: setemis

      implicit none
      type(GFS_control_type),         intent(in)    :: Model
      type(GFS_radtend_type),         intent(inout) :: Radtend
      type(GFS_sfcprop_type),         intent(in)    :: Sfcprop
      type(GFS_grid_type),            intent(in)    :: Grid
      integer, intent(in)                           :: im
      real(kind=kind_phys), dimension(size(Grid%xlon,1)), intent(in) :: tsfa, tsfg
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      if (Model%lslwr) then
!>  - Call module_radiation_surface::setemis(),to setup surface
!! emissivity for LW radiation.
        call setemis (Grid%xlon, Grid%xlat, Sfcprop%slmsk,        &        !  ---  inputs
                     Sfcprop%snowd, Sfcprop%sncovr, Sfcprop%zorl, &
                     tsfg, tsfa, Sfcprop%hprime(:,1), IM,         &
                      Radtend%semis)                                       !  ---  outputs
      endif

      end subroutine rrtmg_lw_pre_run

!> \section arg_table_rrtmg_lw_pre_finalize Argument Table
!!
       subroutine rrtmg_lw_pre_finalize ()
       end subroutine rrtmg_lw_pre_finalize
!! @}
       end module rrtmg_lw_pre

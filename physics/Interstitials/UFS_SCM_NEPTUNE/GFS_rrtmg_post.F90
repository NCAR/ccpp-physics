!>\file GFS_rrtmg_post.F90
!! This file contains the calculation of time averaged output quantities (including total-sky and
!! clear-sky SW and LW fluxes at TOA and surface; conventional
!! 3-domain cloud amount, cloud top and base pressure, and cloud top
!! temperature; aerosols AOD, etc.), store computed results in
!! corresponding slots of array fluxr with appropriate time weights. 

       module GFS_rrtmg_post
       contains

!>\defgroup GFS_rrtmg_post_mod GFS RRTMG Scheme Post
!! This module calculate time averaged output quantities (including total-sky and
!! clear-sky SW and LW fluxes at TOA and surface; conventional
!! 3-domain cloud amount, cloud top and base pressure, and cloud top
!! temperature; aerosols AOD, etc.), store computed results in
!! corresponding slots of array fluxr with appropriate time weights.
!> @{
!> \section arg_table_GFS_rrtmg_post_run Argument Table
!! \htmlinclude GFS_rrtmg_post_run.html
!!
      subroutine GFS_rrtmg_post_run(im, lsswr, topfsw, total_albedo, errmsg, errflg)

      use machine,                             only: kind_phys
      use module_radsw_parameters,             only: topfsw_type

      implicit none

      ! Interface variables
      integer,              intent(in) :: im 
      logical,              intent(in) :: lsswr
      real(kind=kind_phys), dimension(im),        intent(inout) :: total_albedo
      
      type(topfsw_type), dimension(im), intent(in) :: topfsw
      
      character(len=*), intent(out) :: errmsg
      integer, intent(out) :: errflg

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0
      
      if (.not. lsswr) return

!  ---  The total sky (with clouds) shortwave albedo
      total_albedo = 0.0
      if (lsswr) then
        where(topfsw(:)%dnfxc>0) total_albedo(:) = topfsw(:)%upfxc/topfsw(:)%dnfxc
      endif
!
      end subroutine GFS_rrtmg_post_run
!> @}
      end module GFS_rrtmg_post

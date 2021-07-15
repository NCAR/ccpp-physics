!> \file GFS_GWD_generic.F90
!! This file contains the CCPP-compliant orographic gravity wave
!! drag pre interstitial codes.

module GFS_GWD_generic_pre

contains

!! \section arg_table_GFS_GWD_generic_pre_init Argument Table
!! \htmlinclude GFS_GWD_generic_pre_init.html
!!
      subroutine GFS_GWD_generic_pre_init()
      end subroutine GFS_GWD_generic_pre_init

!! \section arg_table_GFS_GWD_generic_pre_run Argument Table
!! \htmlinclude GFS_GWD_generic_pre_run.html
!!
!!  \section general General Algorithm
!!  \section detailed Detailed Algorithm
!!  @{
      subroutine GFS_GWD_generic_pre_run(                               &
     &           im, levs, nmtvr, mntvar,                               &
     &           oc, oa4, clx, theta,                                   &
     &           varss, ocss, oa4ss, clxss,                             &
     &           sigma, gamma, elvmax, lssav, ldiag3d,                  &
     &           dtend, dtidx, index_of_temperature, index_of_x_wind,   &
     &           index_of_y_wind, index_of_process_orographic_gwd,      &
     &           dudt, dvdt, dtdt, dtf,                                 &
     &           flag_for_gwd_generic_tend, errmsg, errflg)

      use machine, only : kind_phys
      implicit none

      integer, intent(in) :: im, levs, nmtvr
      real(kind=kind_phys), intent(in) :: mntvar(:,:)

      real(kind=kind_phys), intent(out) ::                              &
     &  oc(:), oa4(:,:), clx(:,:),                                      &
     &  varss(:), ocss(:), oa4ss(:,:), clxss(:,:),                      &
     &  theta(:), sigma(:), gamma(:), elvmax(:)

      logical, intent(in) :: lssav, ldiag3d, flag_for_gwd_generic_tend
      real(kind=kind_phys), intent(in) :: dtdt(:,:), dudt(:,:), dvdt(:,:)
      ! dtend only allocated only if ldiag3d is .true.
      real(kind=kind_phys), intent(inout) :: dtend(:,:,:)
      integer, intent(in) :: dtidx(:,:), index_of_temperature,          &
     &  index_of_x_wind, index_of_y_wind, index_of_process_orographic_gwd
      real(kind=kind_phys), intent(in) :: dtf

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      integer :: i, k, idtend

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      if (nmtvr == 14) then  ! current operational - as of 2014
        oc(:)     = mntvar(:,2)
        oa4(:,1)  = mntvar(:,3)
        oa4(:,2)  = mntvar(:,4)
        oa4(:,3)  = mntvar(:,5)
        oa4(:,4)  = mntvar(:,6)
        clx(:,1)  = mntvar(:,7)
        clx(:,2)  = mntvar(:,8)
        clx(:,3)  = mntvar(:,9)
        clx(:,4)  = mntvar(:,10)
        theta(:)  = mntvar(:,11)
        gamma(:)  = mntvar(:,12)
        sigma(:)  = mntvar(:,13)
        elvmax(:) = mntvar(:,14)
      elseif (nmtvr == 10) then
        oc(:)     = mntvar(:,2)
        oa4(:,1)  = mntvar(:,3)
        oa4(:,2)  = mntvar(:,4)
        oa4(:,3)  = mntvar(:,5)
        oa4(:,4)  = mntvar(:,6)
        clx(:,1)  = mntvar(:,7)
        clx(:,2)  = mntvar(:,8)
        clx(:,3)  = mntvar(:,9)
        clx(:,4)  = mntvar(:,10)
      elseif (nmtvr == 6) then
        oc(:)     = mntvar(:,2)
        oa4(:,1)  = mntvar(:,3)
        oa4(:,2)  = mntvar(:,4)
        oa4(:,3)  = mntvar(:,5)
        oa4(:,4)  = mntvar(:,6)
        clx(:,1)  = 0.0
        clx(:,2)  = 0.0
        clx(:,3)  = 0.0
        clx(:,4)  = 0.0
      elseif (nmtvr == 24) then   ! GSD_drag_suite and unified_ugwp
        oc(:)       = mntvar(:,2)
        oa4(:,1)    = mntvar(:,3)
        oa4(:,2)    = mntvar(:,4)
        oa4(:,3)    = mntvar(:,5)
        oa4(:,4)    = mntvar(:,6)
        clx(:,1)    = mntvar(:,7)
        clx(:,2)    = mntvar(:,8)
        clx(:,3)    = mntvar(:,9)
        clx(:,4)    = mntvar(:,10)
        theta(:)    = mntvar(:,11)
        gamma(:)    = mntvar(:,12)
        sigma(:)    = mntvar(:,13)
        elvmax(:)   = mntvar(:,14)
        varss(:)    = mntvar(:,15)
        ocss(:)     = mntvar(:,16)
        oa4ss(:,1)  = mntvar(:,17)
        oa4ss(:,2)  = mntvar(:,18)
        oa4ss(:,3)  = mntvar(:,19)
        oa4ss(:,4)  = mntvar(:,20)
        clxss(:,1)  = mntvar(:,21)
        clxss(:,2)  = mntvar(:,22)
        clxss(:,3)  = mntvar(:,23)
        clxss(:,4)  = mntvar(:,24)
      else
        oc     = 0
        oa4    = 0
        clx    = 0
        theta  = 0
        gamma  = 0
        sigma  = 0
        elvmax = 0
      endif   ! end if_nmtvr

      if (lssav .and. ldiag3d .and. flag_for_gwd_generic_tend) then
        idtend = dtidx(index_of_temperature, index_of_process_orographic_gwd)
        if(idtend>=1) then
          dtend(:,:,idtend) = dtend(:,:,idtend) - dtdt*dtf
        endif

        idtend = dtidx(index_of_x_wind, index_of_process_orographic_gwd)
        if(idtend>=1) then
          dtend(:,:,idtend) = dtend(:,:,idtend) - dudt*dtf
        endif

        idtend = dtidx(index_of_y_wind, index_of_process_orographic_gwd)
        if(idtend>=1) then
          dtend(:,:,idtend) = dtend(:,:,idtend) - dvdt*dtf
        endif
      endif

      end subroutine GFS_GWD_generic_pre_run
!> @}

!! \section arg_table_GFS_GWD_generic_pre_finalize Argument Table
!! \htmlinclude GFS_GWD_generic_pre_finalize.html
!!
      subroutine GFS_GWD_generic_pre_finalize()
      end subroutine GFS_GWD_generic_pre_finalize

end module GFS_GWD_generic_pre

!> This module contains the CCPP-compliant orographic gravity wave drag post
!! interstitial codes.
module GFS_GWD_generic_post

contains


      subroutine GFS_GWD_generic_post_init()
      end subroutine GFS_GWD_generic_post_init

!! \section arg_table_GFS_GWD_generic_post_run Argument Table
!! \htmlinclude GFS_GWD_generic_post_run.html
!!
!!  \section general General Algorithm
!!  \section detailed Detailed Algorithm
!!  @{
      subroutine GFS_GWD_generic_post_run(lssav, ldiag3d, dtf, dusfcg, dvsfcg, dudt, dvdt, dtdt,          &
      &  dugwd, dvgwd, flag_for_gwd_generic_tend, dtend, dtidx, index_of_temperature, index_of_x_wind,  &
      &  index_of_y_wind, index_of_process_orographic_gwd, errmsg, errflg)

      use machine, only : kind_phys
      implicit none
      
      logical, intent(in) :: lssav, ldiag3d, flag_for_gwd_generic_tend
      
      real(kind=kind_phys), intent(in) :: dusfcg(:), dvsfcg(:)
      real(kind=kind_phys), intent(in) :: dudt(:,:), dvdt(:,:), dtdt(:,:)
      real(kind=kind_phys), intent(in) :: dtf
      
      real(kind=kind_phys), intent(inout) :: dugwd(:), dvgwd(:)

      ! dtend only allocated only if ldiag3d is .true.
      real(kind=kind_phys), intent(inout) :: dtend(:,:,:)
      integer, intent(in) :: dtidx(:,:), index_of_temperature,          &
     &  index_of_x_wind, index_of_y_wind, index_of_process_orographic_gwd
      
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      integer :: idtend

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      if (lssav) then
        dugwd(:) = dugwd(:) + dusfcg(:)*dtf
        dvgwd(:) = dvgwd(:) + dvsfcg(:)*dtf

        if (ldiag3d .and. flag_for_gwd_generic_tend) then
          idtend = dtidx(index_of_temperature, index_of_process_orographic_gwd)
          if(idtend>=1) then
            dtend(:,:,idtend) = dtend(:,:,idtend) + dtdt*dtf
          endif

          idtend = dtidx(index_of_x_wind, index_of_process_orographic_gwd)
          if(idtend>=1) then
            dtend(:,:,idtend) = dtend(:,:,idtend) + dudt*dtf
          endif

          idtend = dtidx(index_of_y_wind, index_of_process_orographic_gwd)
          if(idtend>=1) then
            dtend(:,:,idtend) = dtend(:,:,idtend) + dvdt*dtf
          endif
        endif
      endif

    end subroutine GFS_GWD_generic_post_run
!> @}
    
!! \section arg_table_GFS_GWD_generic_post_finalize Argument Table
!! \htmlinclude GFS_GWD_generic_post_finalize.html
!!
    subroutine GFS_GWD_generic_post_finalize()
    end subroutine GFS_GWD_generic_post_finalize

end module GFS_GWD_generic_post

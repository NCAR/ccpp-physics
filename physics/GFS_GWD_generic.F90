!> \file GFS_GWD_generic.f
!! This file contains the CCPP-compliant orographic gravity wave
!! drag pre interstitial codes.

module GFS_GWD_generic_pre

contains

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
     &           sigma, gamma, elvmax, lssav, ldiag3d,                  &
     &           dudt, dvdt, dtdt, du3dt, dv3dt, dt3dt, dtf,            &
     &           flag_for_gwd_generic_tend, errmsg, errflg)

      use machine, only : kind_phys
      implicit none

      integer, intent(in) :: im, levs, nmtvr
      real(kind=kind_phys), intent(in) :: mntvar(im,nmtvr)

      real(kind=kind_phys), intent(out) ::                              &
     &  oc(im), oa4(im,4), clx(im,4),                                   &
     &  theta(im), sigma(im), gamma(im), elvmax(im)

      logical, intent(in) :: lssav, ldiag3d, flag_for_gwd_generic_tend
      real(kind=kind_phys), intent(in) :: dtdt(im,levs), dudt(im,levs), dvdt(im,levs)
      ! dt3dt only allocated only if ldiag3d is .true.
      real(kind=kind_phys), intent(inout) :: dt3dt(:,:), du3dt(:,:), dv3dt(:,:)
      real(kind=kind_phys), intent(in) :: dtf

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      integer :: i, k

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
      else
        oc     = 0
        oa4    = 0
        clx    = 0
        theta  = 0
        gamma  = 0
        sigma  = 0
        elvmax = 0
      endif   ! end if_nmtvr

      if (lssav) then
        if (ldiag3d .and. flag_for_gwd_generic_tend) then
          do k=1,levs
            do i=1,im
              dt3dt(i,k) = dt3dt(i,k) - dtdt(i,k)*dtf
              du3dt(i,k) = du3dt(i,k) - dudt(i,k)*dtf
              dv3dt(i,k) = dv3dt(i,k) - dvdt(i,k)*dtf
            enddo
          enddo
        endif
      endif

      end subroutine GFS_GWD_generic_pre_run
!> @}

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
      &  dugwd, dvgwd, du3dt, dv3dt, dt3dt, flag_for_gwd_generic_tend, errmsg, errflg)

      use machine, only : kind_phys
      implicit none
      
      logical, intent(in) :: lssav, ldiag3d, flag_for_gwd_generic_tend
      
      real(kind=kind_phys), intent(in) :: dusfcg(:), dvsfcg(:)
      real(kind=kind_phys), intent(in) :: dudt(:,:), dvdt(:,:), dtdt(:,:)
      real(kind=kind_phys), intent(in) :: dtf
      
      real(kind=kind_phys), intent(inout) :: dugwd(:), dvgwd(:)
      real(kind=kind_phys), intent(inout) :: du3dt(:,:), dv3dt(:,:), dt3dt(:,:)
      
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      if (lssav) then
        dugwd(:) = dugwd(:) + dusfcg(:)*dtf
        dvgwd(:) = dvgwd(:) + dvsfcg(:)*dtf

        if (ldiag3d .and. flag_for_gwd_generic_tend) then
          du3dt(:,:) = du3dt(:,:) + dudt(:,:) * dtf
          dv3dt(:,:) = dv3dt(:,:) + dvdt(:,:) * dtf
          dt3dt(:,:) = dt3dt(:,:) + dtdt(:,:) * dtf
        endif
      endif

    end subroutine GFS_GWD_generic_post_run
!> @}
    
    subroutine GFS_GWD_generic_post_finalize()
    end subroutine GFS_GWD_generic_post_finalize

end module GFS_GWD_generic_post

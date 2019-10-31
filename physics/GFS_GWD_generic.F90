!> \file GFS_GWD_generic.f
!! This file contains the CCPP-compliant orographic gravity wave
!! drag pre interstitial codes.

module GFS_GWD_generic_pre

contains

!> \section arg_table_GFS_GWD_generic_pre_init Argument Table
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
     &           sigma, gamma, elvmax, lssav, ldiag3d,                  &
     &           dtdt, dt3dt, dtf, errmsg, errflg)

      use machine, only : kind_phys
      implicit none

      integer, intent(in) :: im, levs, nmtvr
      real(kind=kind_phys), intent(in) :: mntvar(im,nmtvr)

      real(kind=kind_phys), intent(out) ::                              &
     &  oc(im), oa4(im,4), clx(im,4),                                   &
     &  theta(im), sigma(im), gamma(im), elvmax(im)

      logical, intent(in) :: lssav, ldiag3d
      real(kind=kind_phys), intent(in) :: dtdt(im,levs)
      ! dt3dt only allocated only if ldiag3d is .true.
      real(kind=kind_phys), intent(inout) :: dt3dt(:,:)
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
        if (ldiag3d) then
          do k=1,levs
            do i=1,im
              dt3dt(i,k) = dt3dt(i,k) - dtdt(i,k)*dtf
            enddo
          enddo
        endif
      endif

      end subroutine GFS_GWD_generic_pre_run
!> @}

! \ingroup GFS_ogwd
! \brief Brief description of the subroutine
!
!> \section arg_table_GFS_GWD_generic_pre_finalize Argument Table
!!
      subroutine GFS_GWD_generic_pre_finalize()
      end subroutine GFS_GWD_generic_pre_finalize

end module GFS_GWD_generic_pre

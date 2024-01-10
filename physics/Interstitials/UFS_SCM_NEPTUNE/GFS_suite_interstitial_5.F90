!> \file GFS_suite_interstitial_5.F90
!!  Contains code to update cloud liquid and ice in the convective transportable tracer array before RAS convection.

  module GFS_suite_interstitial_5

  contains

!> \section arg_table_GFS_suite_interstitial_5_run Argument Table
!! \htmlinclude GFS_suite_interstitial_5_run.html
!!
    subroutine GFS_suite_interstitial_5_run (im, levs, ntrac, ntcw, ntiw, nn, gq0, clw, errmsg, errflg)

      use machine, only: kind_phys

      implicit none

      ! interface variables
      integer,              intent(in )                   :: im, levs, ntrac, ntcw, ntiw, nn

      real(kind=kind_phys), intent(in ), dimension(:,:,:) :: gq0

      real(kind=kind_phys), intent(out), dimension(:,:,:) :: clw

      character(len=*),     intent(out)                   :: errmsg
      integer,              intent(out)                   :: errflg

      ! local variables
      integer :: i,k

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      do k=1,levs
        do i=1,im
          clw(i,k,1) = gq0(i,k,ntiw)                    ! ice
          clw(i,k,2) = gq0(i,k,ntcw)                    ! water
        enddo
      enddo

    end subroutine GFS_suite_interstitial_5_run

  end module GFS_suite_interstitial_5

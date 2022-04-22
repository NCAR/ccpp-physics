!> \file GFS_suite_interstitial_1.f90
!!  Contains code to calculate scale-aware variables used in cs_conv, gwdc, and precpd and to reset tendencies used in the
!!  process-split section of GFS-based physics suites.

    module GFS_suite_interstitial_1

    contains

!> \section arg_table_GFS_suite_interstitial_1_run Argument Table
!! \htmlinclude GFS_suite_interstitial_1_run.html
!!
    subroutine GFS_suite_interstitial_1_run (im, levs, ntrac, dtf, dtp, slmsk, area, dxmin, dxinv, pgr, &
      islmsk, work1, work2, psurf, dudt, dvdt, dtdt, dqdt, errmsg, errflg)

      use machine, only: kind_phys

      implicit none

      ! interface variables
      integer,              intent(in )                   :: im, levs, ntrac
      real(kind=kind_phys), intent(in )                   :: dtf, dtp, dxmin, dxinv
      real(kind=kind_phys), intent(in ), dimension(:)     :: slmsk, area, pgr

      integer,              intent(out), dimension(:)     :: islmsk
      real(kind=kind_phys), intent(out), dimension(:)     :: work1, work2, psurf
      real(kind=kind_phys), intent(out), dimension(:,:)   :: dudt, dvdt, dtdt
      real(kind=kind_phys), intent(out), dimension(:,:,:) :: dqdt
      
      character(len=*),     intent(out)                   :: errmsg
      integer,              intent(out)                   :: errflg

      ! local variables
      real(kind=kind_phys), parameter   :: zero = 0.0_kind_phys, one = 1.0_kind_phys
      integer :: i, k, n

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      do i = 1, im
        islmsk(i)   = nint(slmsk(i))

        work1(i) = (log(area(i)) - dxmin) * dxinv
        work1(i) = max(zero, min(one, work1(i)))
        work2(i) = one - work1(i)
        psurf(i) = pgr(i)
      end do

      do k=1,levs
        do i=1,im
          dudt(i,k)  = zero
          dvdt(i,k)  = zero
          dtdt(i,k)  = zero
        enddo
      enddo
      do n=1,ntrac
        do k=1,levs
          do i=1,im
            dqdt(i,k,n) = zero
          enddo
        enddo
      enddo

    end subroutine GFS_suite_interstitial_1_run

  end module GFS_suite_interstitial_1
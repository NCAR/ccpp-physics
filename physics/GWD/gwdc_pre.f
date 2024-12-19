!> \file gwdc_pre.f This file is preparation for the original code for parameterization of
!! stationary convection forced gravity wave drag based on
!! Chun and Baik (1998) \cite chun_and_baik_1998.

      module gwdc_pre
      contains

!! \section arg_table_gwdc_pre_run Argument Table
!! \htmlinclude gwdc_pre_run.html
!!
      subroutine gwdc_pre_run (                                         &
     &  im, cgwf, dx, work1, work2, dlength, cldf,                      &
     &  levs, kbot, ktop, dtp, gt0, gt0_init, del, cumabs,              &
     &  errmsg, errflg )

      use machine, only : kind_phys
      implicit none

      integer, intent(in) :: im, levs
      integer, intent(in) :: kbot(:), ktop(:)
      real(kind=kind_phys), intent(in) :: dtp
      real(kind=kind_phys), intent(in) :: cgwf(:)
      real(kind=kind_phys), intent(in) :: dx(:), work1(:), work2(:)
      real(kind=kind_phys), intent(in) ::                               &
     &  gt0(:,:), gt0_init(:,:), del(:,:)

      real(kind=kind_phys), intent(out) ::                              &
     &  dlength(:), cldf(:), cumabs(:)

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      integer :: i, k
      real(kind=kind_phys) :: tem1, tem2
      real(kind=kind_phys) :: work3(im)

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      do i = 1, im
        tem1       = dx(i)
        tem2       = tem1
        dlength(i) = sqrt( tem1*tem1+tem2*tem2 )
        cldf(i)    = cgwf(1)*work1(i) + cgwf(2)*work2(i)
      enddo

!  --- ...  calculate maximum convective heating rate
!           cuhr = temperature change due to deep convection

      cumabs(:) = 0.0
      work3(:)  = 0.0
      do k = 1, levs
        do i = 1, im
          if (k >= kbot(i) .and. k <= ktop(i)) then
            cumabs(i)                                                   &
     &        = cumabs(i) + (gt0(i,k) - gt0_init(i,k)) * del(i,k)
            work3(i)  = work3(i)  + del(i,k)
          endif
        enddo
      enddo
      do i=1,im
        if (work3(i) > 0.0) cumabs(i) = cumabs(i) / (dtp*work3(i))
      enddo

      end subroutine gwdc_pre_run

      end module gwdc_pre
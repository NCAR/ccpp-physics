!>\file tridi.f
!! These subroutines are originally internal subroutines in moninedmf.f

!>\ingroup HEDMF
!!\brief Routine to solve the tridiagonal system to calculate
!!temperature and moisture at \f$ t + \Delta t \f$; part of two-part
!!process to calculate time tendencies due to vertical diffusion.
      subroutine tridi1(l,n,cl,cm,cu,r1,au,a1)
      !
      use machine     , only : kind_phys
      implicit none
      integer, parameter :: one = 1.0_kind_phys
      integer             k,n,l,i
      real(kind=kind_phys) fk
      !
      real(kind=kind_phys) cl(l,2:n),cm(l,n),cu(l,n-1),r1(l,n),         &
     &                     au(l,n-1),a1(l,n)
      !
      do i=1,l
        fk      = one / cm(i,1)
        au(i,1) = fk*cu(i,1)
        a1(i,1) = fk*r1(i,1)
      enddo
      do k=2,n-1
        do i=1,l
          fk      = one / (cm(i,k)-cl(i,k)*au(i,k-1))
          au(i,k) = fk*cu(i,k)
          a1(i,k) = fk*(r1(i,k)-cl(i,k)*a1(i,k-1))
        enddo
      enddo
      do i=1,l
        fk      = one / (cm(i,n)-cl(i,n)*au(i,n-1))
        a1(i,n) = fk*(r1(i,n)-cl(i,n)*a1(i,n-1))
      enddo
      do k=n-1,1,-1
        do i=1,l
          a1(i,k) = a1(i,k)-au(i,k)*a1(i,k+1)
        enddo
      enddo
      !
      return
      end subroutine tridi1

!-----------------------------------------------------------------------
!>\ingroup satmedmf
!>\ingroup satmedmfvdifq
!> This subroutine ..
      subroutine tridi2(l,n,cl,cm,cu,r1,r2,au,a1,a2)
!
      use machine     , only : kind_phys
      implicit none
      integer, parameter :: one = 1.0_kind_phys
      integer             k,n,l,i
      real(kind=kind_phys) fk
!
      real(kind=kind_phys) cl(l,2:n),cm(l,n),cu(l,n-1),r1(l,n),r2(l,n), &
     &          au(l,n-1),a1(l,n),a2(l,n)
!----------------------------------------------------------------------
      do i=1,l
        fk      = one / cm(i,1)
        au(i,1) = fk*cu(i,1)
        a1(i,1) = fk*r1(i,1)
        a2(i,1) = fk*r2(i,1)
      enddo
      do k=2,n-1
        do i=1,l
          fk      = one / (cm(i,k)-cl(i,k)*au(i,k-1))
          au(i,k) = fk*cu(i,k)
          a1(i,k) = fk*(r1(i,k)-cl(i,k)*a1(i,k-1))
          a2(i,k) = fk*(r2(i,k)-cl(i,k)*a2(i,k-1))
        enddo
      enddo
      do i=1,l
        fk      = one / (cm(i,n)-cl(i,n)*au(i,n-1))
        a1(i,n) = fk*(r1(i,n)-cl(i,n)*a1(i,n-1))
        a2(i,n) = fk*(r2(i,n)-cl(i,n)*a2(i,n-1))
      enddo
      do k=n-1,1,-1
        do i=1,l
          a1(i,k) = a1(i,k)-au(i,k)*a1(i,k+1)
          a2(i,k) = a2(i,k)-au(i,k)*a2(i,k+1)
        enddo
      enddo
!-----------------------------------------------------------------------
      return
      end subroutine tridi2

!-----------------------------------------------------------------------
!>\ingroup satmedmf
!>\ingroup satmedmfvdifq
!>  Routine to solve the tridiagonal system to calculate u- and
!!  v-momentum at \f$ t + \Delta t \f$; part of two-part process to
!!  calculate time tendencies due to vertical diffusion.
      subroutine tridin(l,n,nt,cl,cm,cu,r1,r2,au,a1,a2)
!
      use machine     , only : kind_phys
      implicit none
      integer, parameter :: one = 1.0_kind_phys
      integer             is,k,kk,n,nt,l,i
      real(kind=kind_phys) fk(l)
!
      real(kind=kind_phys) cl(l,2:n), cm(l,n), cu(l,n-1),               &
     &                     r1(l,n),   r2(l,n*nt),                       &
     &                     au(l,n-1), a1(l,n), a2(l,n*nt),              &
     &                     fkk(l,2:n-1)
!-----------------------------------------------------------------------
      do i=1,l
        fk(i)   = one / cm(i,1)
        au(i,1) = fk(i)*cu(i,1)
        a1(i,1) = fk(i)*r1(i,1)
      enddo
      do k = 1, nt
        is = (k-1) * n
        do i = 1, l
          a2(i,1+is) = fk(i) * r2(i,1+is)
        enddo
      enddo
      do k=2,n-1
        do i=1,l
          fkk(i,k) = one / (cm(i,k)-cl(i,k)*au(i,k-1))
          au(i,k)  = fkk(i,k)*cu(i,k)
          a1(i,k)  = fkk(i,k)*(r1(i,k)-cl(i,k)*a1(i,k-1))
        enddo
      enddo
      do kk = 1, nt
        is = (kk-1) * n
        do k=2,n-1
          do i=1,l
            a2(i,k+is) = fkk(i,k)*(r2(i,k+is)-cl(i,k)*a2(i,k+is-1))
          enddo
        enddo
      enddo
      do i=1,l
        fk(i)   = one / (cm(i,n)-cl(i,n)*au(i,n-1))
        a1(i,n) = fk(i)*(r1(i,n)-cl(i,n)*a1(i,n-1))
      enddo
      do k = 1, nt
        is = (k-1) * n
        do i = 1, l
          a2(i,n+is) = fk(i)*(r2(i,n+is)-cl(i,n)*a2(i,n+is-1))
        enddo
      enddo
      do k=n-1,1,-1
        do i=1,l
          a1(i,k) = a1(i,k) - au(i,k)*a1(i,k+1)
        enddo
      enddo
      do kk = 1, nt
        is = (kk-1) * n
        do k=n-1,1,-1
          do i=1,l
            a2(i,k+is) = a2(i,k+is) - au(i,k)*a2(i,k+is+1)
          enddo
        enddo
      enddo
!-----------------------------------------------------------------------
      return
      end subroutine tridin

!-----------------------------------------------------------------------
!>\ingroup satmedmf
!>\ingroup satmedmfvdifq
!! This subroutine solves tridiagonal problem for TKE.
      subroutine tridit(l,n,nt,cl,cm,cu,rt,au,at)
!-----------------------------------------------------------------------
!!
      use machine     , only : kind_phys
      implicit none
      integer, parameter :: one = 1.0_kind_phys
      integer             is,k,kk,n,nt,l,i
      real(kind=kind_phys) fk(l)
!!
      real(kind=kind_phys) cl(l,2:n), cm(l,n), cu(l,n-1),               &
     &                     rt(l,n*nt),                                  &
     &                     au(l,n-1), at(l,n*nt),                       &
     &                     fkk(l,2:n-1)
!-----------------------------------------------------------------------
      do i=1,l
        fk(i)   = one / cm(i,1)
        au(i,1) = fk(i)*cu(i,1)
      enddo
      do k = 1, nt
        is = (k-1) * n
        do i = 1, l
          at(i,1+is) = fk(i) * rt(i,1+is)
        enddo
      enddo
      do k=2,n-1
        do i=1,l
          fkk(i,k) = one / (cm(i,k)-cl(i,k)*au(i,k-1))
          au(i,k)  = fkk(i,k)*cu(i,k)
        enddo
      enddo
      do kk = 1, nt
        is = (kk-1) * n
        do k=2,n-1
          do i=1,l
            at(i,k+is) = fkk(i,k)*(rt(i,k+is)-cl(i,k)*at(i,k+is-1))
          enddo
        enddo
      enddo
      do i=1,l
        fk(i)   = one / (cm(i,n)-cl(i,n)*au(i,n-1))
      enddo
      do k = 1, nt
        is = (k-1) * n
        do i = 1, l
          at(i,n+is) = fk(i)*(rt(i,n+is)-cl(i,n)*at(i,n+is-1))
        enddo
      enddo
      do kk = 1, nt
        is = (kk-1) * n
        do k=n-1,1,-1
          do i=1,l
            at(i,k+is) = at(i,k+is) - au(i,k)*at(i,k+is+1)
          enddo
        enddo
      enddo
!-----------------------------------------------------------------------
      return
      end subroutine tridit
!> @}

!>  \file sfc_diag.f
!!  This file contains the land surface diagnose calculation scheme.

!> \defgroup Sfc_diag Land Surface Diagnose Calculation
!! @{

      module sfc_diag
      contains
  
      subroutine sfc_diag_init
      end subroutine sfc_diag_init
      
      subroutine sfc_diag_finalize
      end subroutine sfc_diag_finalize
      
!> \brief Brief description of the subroutine
!!
!! \section arg_table_sfc_diag_run Arguments
!! \htmlinclude sfc_diag_run.html
!!
!!  \section general General Algorithm
!!  \section detailed Detailed Algorithm
!!  @{
      subroutine sfc_diag_run                                           &
     &                   (im,grav,cp,eps,epsm1,ps,u1,v1,t1,q1,prslki,   &
     &                    evap,fm,fh,fm10,fh2,tskin,qsurf,thsfc_loc,    &
     &                    f10m,u10m,v10m,t2m,q2m,errmsg,errflg          &
     &                   )
!
      use machine , only : kind_phys
      use funcphys, only : fpvs
      implicit none
!
      integer, intent(in) :: im
      logical, intent(in) :: thsfc_loc  ! Flag for reference pot. temp.
      real(kind=kind_phys), intent(in) :: grav,cp,eps,epsm1
      real(kind=kind_phys), dimension(:), intent(in) ::                 &
     &                       ps, u1, v1, t1, q1, tskin,                 &
     &                       qsurf, prslki, evap, fm, fh, fm10, fh2
      real(kind=kind_phys), dimension(:), intent(out) ::                &
     &                       f10m, u10m, v10m, t2m, q2m
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
!
!     locals
!
      real(kind=kind_phys), parameter :: qmin=1.0e-8
      integer :: k,i
!
      real(kind=kind_phys) :: fhi, qss, wrk
!     real(kind=kind_phys) sig2k, fhi, qss
!
!     real, parameter :: g=grav
!
      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0
!
!     estimate sigma ** k at 2 m
!
!     sig2k = 1. - 4. * g * 2. / (cp * 280.)
!
!  initialize variables. all units are supposedly m.k.s. unless specified
!  ps is in pascals
!
!!
      do i = 1, im
        f10m(i) = fm10(i) / fm(i)
!       f10m(i) = min(f10m(i),1.)
        u10m(i) = f10m(i) * u1(i)
        v10m(i) = f10m(i) * v1(i)
        fhi     = fh2(i) / fh(i)
!       t2m(i)  = tskin(i)*(1. - fhi) + t1(i) * prslki(i) * fhi
!       sig2k   = 1. - (grav+grav) / (cp * t2m(i))
!       t2m(i)  = t2m(i) * sig2k
        wrk     = 1.0 - fhi


        if(thsfc_loc) then ! Use local potential temperature
          t2m(i)  = tskin(i)*wrk + t1(i)*prslki(i)*fhi - (grav+grav)/cp
        else ! Use potential temperature referenced to 1000 hPa
          t2m(i)  = tskin(i)*wrk + t1(i)*fhi - (grav+grav)/cp
        endif

        if(evap(i) >= 0.) then !  for evaporation>0, use inferred qsurf to deduce q2m
          q2m(i) = qsurf(i)*wrk + max(qmin,q1(i))*fhi
        else                   !  for dew formation, use saturated q at tskin
          qss    = fpvs(tskin(i))
          qss    = eps * qss / (ps(i) + epsm1 * qss)
          q2m(i) = qss*wrk + max(qmin,q1(i))*fhi
        endif
        qss    = fpvs(t2m(i))
        qss    = eps * qss / (ps(i) + epsm1 * qss)
        q2m(i) = min(q2m(i),qss)
      enddo

      return
      end subroutine sfc_diag_run
!> @}

      end module sfc_diag
!> @}

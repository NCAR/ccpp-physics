!>\file sfc_data.F
!! This file contains an data surface scheme.

!> This module contains the CCPP-compliant CDEPS data scheme 
!! scheme when the model is using data provided by CDEPS.
      module sfc_data
      implicit none
      private
      public :: sfc_data_run

      contains

!>\defgroup gfs_data_main Simple Wrapper for CDEPS inline 
!! This subroutine pass CDEPS inline provided data to other schemes
!! by overwriting the specific part
!>@{
!! \section arg_table_sfc_data_run Argument Table
!! \htmlinclude sfc_data_run.html
!!
!!>\section gen_sfc_data CDEPS Inline data scheme
      subroutine sfc_data_run                                           &
!...................................
!  ---  inputs:
     &     ( im, use_data,                                              &
     &       tsfco_dat, mask_dat, tice_dat, hice_dat, fice_dat,         &
     &       hvap, tgice, cp, eps, epsm1, rvrdm1, rd, ps, t1, q1,       &
     &       cm, ch, prsl1, prslki, prsik1, prslk1,                     &
     &       wind, thsfc_loc,                                           &
!  ---  input/outputs:
     &       hice, fice, tice, tsfc_wat, qss_i, qss_w, cmm, chh,        &
     &       evapi, evapw, hflxi, hflxw,                                &
     &       errmsg, errflg,                                            &
     &       naux2d, aux2d)

! ===================================================================== !
!  description:                                                         !
!                                                                       !
!  usage:                                                               !
!                                                                       !
!    call sfc_data                                                      !
!       inputs:                                                         !
!        ( im, use_data, tsfco_dat                                      !
!       outputs:                                                        !
!          tsfc_wat )                                                   !
!                                                                       !
!                                                                       !
!  subprograms/functions called: fpvs                                   !
!                                                                       !
!                                                                       !
!  program history log:                                                 !
!    feb  2024  -- u. turuncoglu added initially to interact with       !
!                  fv3 cdeps inline                                     !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                       size   !
!     im       - integer, horizontal dimension                     1    !
!     use_data - logical, =T if data provided by cdeps inline      1    !
!     tsfco_dat- real, sea surface temperature ( k )               im   !
!     mask_dat - data model mask                                   1    !
!                                                                       !
!  outputs:                                                             !
!     tsfc_wat - real, sea surface temperature ( k )               im   !
!                                                                       !
! ===================================================================== !
!
      use machine , only : kind_phys
      use funcphys, only : fpvs
!
      implicit none

!  ---  constants:
      real(kind=kind_phys), parameter :: one  = 1.0_kind_phys
      real(kind=kind_phys), parameter :: qmin = 1.0e-8_kind_phys
      real(kind=kind_phys), parameter :: zero  = 0.0_kind_phys

!  ---  inputs:
      integer, intent(in) :: im
      logical, intent(in) :: use_data 
      logical, intent(in) :: thsfc_loc
      real (kind=kind_phys), dimension(:), intent(in) :: mask_dat,      &
     & hice_dat, fice_dat, tice_dat, tsfco_dat
      real (kind=kind_phys), dimension(:), intent(in) :: q1, t1, prslki,&
     & prslk1, prsik1, prsl1, ps, cm, ch, wind
      real (kind=kind_phys), intent(in) :: cp, eps, epsm1, hvap, rd,    &
     & rvrdm1, tgice

!  ---  input/outputs:
      real (kind=kind_phys), dimension(:), intent(inout) :: hice, fice, &
     & tice, tsfc_wat
      integer, intent(in) :: naux2d
      real(kind_phys), intent(inout) :: aux2d(:,:)

!  ---  outputs:
      real (kind=kind_phys), dimension(:), intent(inout) :: cmm, chh,   &
     & evapi, evapw, hflxi, hflxw, qss_i, qss_w
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

!  ---  locals:
      integer :: i
      real (kind=kind_phys) :: cpinv, elocp, hvapi,                     &
     & q0, qs1, qssi, qssw, tem
      real (kind=kind_phys), dimension(im) :: rch, rho, sneti,          &
     & theta1

      ! calculate some constants
      cpinv = one/cp
      elocp = hvap/cp
      hvapi = one/hvap

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      ! Check coupling from component land to atmosphere
      if (.not. use_data) return

      do i = 1, im
         if (mask_dat(i) > 0.0) then
            ! overwrite internal variables
            tice(i) = tice_dat(i)
            hice(i) = hice_dat(i)
            fice(i) = fice_dat(i) 
            tsfc_wat(i) = tsfco_dat(i)

            aux2d(i,1) = cmm(i)
            aux2d(i,2) = chh(i)
            aux2d(i,3) = evapi(i)
            aux2d(i,4) = evapw(i)
            aux2d(i,5) = hflxi(i)
            aux2d(i,6) = hflxw(i)

            ! sfc_sice calculates fluxes only for islmsk == 2
            if (fice(i) > zero) then
               q0 = max(q1(i), qmin)
               if (thsfc_loc) then
                  theta1(i) = t1(i)*prslki(i)
               else
                  theta1(i) = t1(i)/prslk1(i)
               end if
               rho(i) = prsl1(i)/(rd*t1(i)*(one+rvrdm1*q0))
               qs1 = fpvs(t1(i))
               qs1 = max(eps*qs1/(prsl1(i)+epsm1*qs1), qmin)
               q0 = min(qs1, q0)
               qssi = fpvs(tice(i))
               qssi = eps*qssi/(ps(i)+epsm1*qssi)
               qssw = fpvs(tgice)
               qssw = eps*qssw/(ps(i)+epsm1*qssw)
               cmm(i) = cm(i)*wind(i)
               chh(i) = rho(i)*ch(i)*wind(i)
               rch(i) = chh(i)*cp
               evapi(i) = elocp*rch(i)*(qssi-q0)
               evapw(i) = elocp*rch(i)*(qssw-q0)
               if (thsfc_loc) then
                  hflxi(i) = rch(i)*(tice(i)-theta1(i))
                  hflxw(i) = rch(i)*(tgice-theta1(i))
               else
                  tem = one/prsik1(i)
                  hflxi(i) = rch(i)*(tice(i)*tem-theta1(i))
                  hflxw(i) = rch(i)*(tgice*tem-theta1(i))
               end if 
               tsfc_wat(i) = tgice
               qss_i(i) = q1(i)+evapi(i)/(elocp*rch(i))
               qss_w(i) = q1(i)+evapw(i)/(elocp*rch(i))
               !weasd(i) = snowd(i)*1000.0_kind_phys
               !snwdph(i) = weasd(i)*dsi
               tem = one/rho(i)
               hflxi(i) = hflxi(i)*tem*cpinv
               hflxw(i) = hflxw(i)*tem*cpinv
               evapi(i) = evapi(i)*tem*hvapi
               evapw(i) = evapw(i)*tem*hvapi 
            end if

            aux2d(i,7) = cmm(i)
            aux2d(i,8) = chh(i)
            aux2d(i,9) = evapi(i)
            aux2d(i,10) = evapw(i)
            aux2d(i,11) = hflxi(i)
            aux2d(i,12) = hflxw(i)

         end if
      end do
!
      return
!...................................
      end subroutine sfc_data_run
!-----------------------------------
!>@}
      end module sfc_data

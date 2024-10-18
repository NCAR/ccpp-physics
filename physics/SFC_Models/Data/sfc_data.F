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
     &       hice, fice, tice, tsfc_wat, qss_i, qss_w,                  &
     &       cmm_i, cmm_w, chh_i, chh_w,                                &
     &       evap_i, evap_w, hflx_i, hflx_w, gflux,                     &
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
      real (kind=kind_phys), dimension(:), intent(in), optional ::      &
     & mask_dat, fice_dat, hice_dat, tsfco_dat, tice_dat
      real (kind=kind_phys), dimension(:), intent(in) :: q1, t1, prslki,&
     & prslk1, prsik1, prsl1, ps, cm, ch, wind
      real (kind=kind_phys), intent(in) :: cp, eps, epsm1, hvap, rd,    &
     & rvrdm1, tgice

!  ---  input/outputs:
      real (kind=kind_phys), dimension(:), intent(inout) :: hice, fice, &
     & tice, tsfc_wat
      integer, intent(in), optional :: naux2d
      real(kind_phys), intent(inout), optional :: aux2d(:,:)

!  ---  outputs:
      real (kind=kind_phys), dimension(:), intent(inout) ::             &
     & cmm_i, cmm_w, chh_i, chh_w,                                      &
     & evap_i, evap_w, hflx_i, hflx_w, qss_i, qss_w, gflux
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

         aux2d(i,1) = mask_dat(i) 
         aux2d(i,2) = tice_dat(i)
         aux2d(i,3) = hice_dat(i)
         aux2d(i,4) = fice_dat(i)
         aux2d(i,5) = tsfco_dat(i)

         if (mask_dat(i) > 0.0) then
            ! overwrite internal variables
            tice(i) = tice_dat(i)
            hice(i) = hice_dat(i)
            fice(i) = fice_dat(i) 
            tsfc_wat(i) = tsfco_dat(i)

            aux2d(i,6) = cmm_i(i)
            aux2d(i,7) = chh_i(i)
            aux2d(i,8) = evap_i(i)
            aux2d(i,9) = hflx_i(i)
            aux2d(i,10) = qss_i(i)
            aux2d(i,11) = cmm_w(i)
            aux2d(i,12) = chh_w(i)
            aux2d(i,13) = evap_w(i)
            aux2d(i,14) = hflx_w(i)
            aux2d(i,15) = qss_w(i)

            ! sfc_sice calculates fluxes only for islmsk == 2
            if (fice(i) > zero) then ! calculate fluxes over sea-ice
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
               cmm_i(i) = cm(i)*wind(i)
               chh_i(i) = rho(i)*ch(i)*wind(i)
               rch(i) = chh_i(i)*cp
               evap_i(i) = elocp*rch(i)*(qssi-q0)
               if (thsfc_loc) then
                  hflx_i(i) = rch(i)*(tice(i)-theta1(i))
               else
                  tem = one/prsik1(i)
                  hflx_i(i) = rch(i)*(tice(i)*tem-theta1(i))
               end if 
               tsfc_wat(i) = tgice
               qss_i(i) = q1(i)+evap_i(i)/(elocp*rch(i))
               tem = one/rho(i)
               hflx_i(i) = hflx_i(i)*tem*cpinv
               evap_i(i) = evap_i(i)*tem*hvapi
            else ! calculate fluxes over sea
               q0 = max(q1(i), qmin)
               rho(i) = prsl1(i)/(rd*t1(i)*(one+rvrdm1*q0))
               qssw = fpvs(tsfc_wat(i))
               qssw = eps*qssw/(ps(i)+epsm1*qssw)
               cmm_w(i) = cm(i)*wind(i)
               chh_w(i) = rho(i)*ch(i)*wind(i)
               rch(i) = chh_w(i)*cp
               tem = one/rho(i)
               hflx_w(i) = rch(i)*(tsfc_wat(i)-t1(i)*prslki(i)) 
               evap_w(i) = elocp*rch(i)*(qssw-q0)
               hflx_w(i) = hflx_w(i)*tem*cpinv
               evap_w(i) = evap_w(i)*tem*hvapi
               qss_w(i) = qssw
               gflux(i) = zero
            end if

            aux2d(i,16) = cmm_i(i)
            aux2d(i,17) = chh_i(i)
            aux2d(i,18) = evap_i(i)
            aux2d(i,19) = hflx_i(i)
            aux2d(i,20) = qss_i(i)
            aux2d(i,21) = cmm_w(i)
            aux2d(i,22) = chh_w(i)
            aux2d(i,23) = evap_w(i)
            aux2d(i,24) = hflx_w(i)
            aux2d(i,25) = qss_w(i)

         end if
      end do
!
      return
!...................................
      end subroutine sfc_data_run
!-----------------------------------
!>@}
      end module sfc_data

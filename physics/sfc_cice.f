!>  \file sfc_cice.f
!!  This file contains the sfc_sice for coupling to CICE

!> This module contains the CCPP-compliant GFS sea ice post
!! interstitial codes, which returns updated ice thickness and 
!! concentration to global arrays where there is no ice, and 
!! set temperature to surface skin temperature.

!> This module contains the CCPP-compliant GFS sea ice scheme.
      module sfc_cice

      contains

      subroutine sfc_cice_init
      end subroutine sfc_cice_init
!
      subroutine sfc_cice_finalize
      end subroutine sfc_cice_finalize


!> \defgroup sfc_sice for coupling to CICE
!! @{
!!  \section diagram Calling Hierarchy Diagram
!!  \section intraphysics Intraphysics Communication
!!
!> \brief Brief description of the subroutine
!!
!! \section arg_table_sfc_cice_run Arguments
!! \htmlinclude sfc_cice_run.html
!!

!!
!!  \section general General Algorithm
!!  \section detailed Detailed Algorithm
!!  @{


!!      use physcons, only : hvap => con_hvap,  cp => con_cp,           &
!!    &                     rvrdm1 => con_fvirt, rd => con_rd
!
!-----------------------------------
      subroutine sfc_cice_run                                           &
!  ---  inputs:
     &     ( im, cplflx, hvap, cp, rvrdm1, rd,                          &
     &       t1, q1, cm, ch, prsl1,                                     &
     &       wind, flag_cice, flag_iter, dqsfc, dtsfc,                  &
     &       dusfc, dvsfc, snowd,                                       &
!  ---  outputs:
     &       qsurf, cmm, chh, evap, hflx, stress, weasd, snwdph, ep,    &
     &       errmsg, errflg
     &     )

! ===================================================================== !
!  description:                                                         !
!  Sep 2015  --  Xingren Wu created from sfc_sice for coupling to CICE  !
!                                                                       !
!  usage:                                                               !
!                                                                       !
!    call sfc_cice                                                      !
!       inputs:                                                         !
!          ( im, cplflx, hvap, cp, rvrdm1, rd,                          !
!            t1, q1, cm, ch, prsl1,                                     !
!            wind, flag_cice, flag_iter, dqsfc, dtsfc,                  !
!            dusfc, dvsfc, snowd,                                       !
!       outputs:                                                        !
!            qsurf, cmm, chh, evap, hflx, stress, weasd, snwdph, ep)    !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:
!     im, - integer, horiz dimension
!!    u1, v1   - real, u/v component of surface layer wind
!     t1       - real, surface layer mean temperature ( k )
!     q1       - real, surface layer mean specific humidity
!     cm       - real, surface exchange coeff for momentum (m/s)
!     ch       - real, surface exchange coeff heat & moisture(m/s)
!     prsl1    - real, surface layer mean pressure
!     wind     - real, wind speed (m/s)
!     flag_iter- logical
!     dqsfc    - real, latent heat flux
!     dtsfc    - real, sensible heat flux
!     dusfc    - real, zonal momentum stress
!     dvsfc    - real, meridional momentum stress
!     snowd    - real, snow depth from cice
!  outputs:
!     qsurf    - real, specific humidity at sfc
!     cmm      - real, ?
!     chh      - real, ?
!     evap     - real, evaperation from latent heat
!     hflx     - real, sensible heat
!     stress   - real, surface stress
!     weasd    - real, water equivalent accumulated snow depth (mm)
!     snwdph   - real, water equivalent snow depth (mm)
!     ep       - real, potential evaporation
!  ====================    end of description    =====================  !
!
!
      use machine , only : kind_phys
      implicit none

      real(kind=kind_phys), parameter   :: one = 1.0_kind_phys
      real(kind=kind_phys), parameter   :: dsi = one/0.33_kind_phys
      real (kind=kind_phys), intent(in) :: hvap, cp, rvrdm1, rd

!  ---  inputs:
      integer, intent(in) :: im
      logical, intent(in) :: cplflx

!     real (kind=kind_phys), dimension(im), intent(in) :: u1, v1,       &
      real (kind=kind_phys), dimension(im), intent(in) ::               &
     &       t1, q1, cm, ch, prsl1, wind, dqsfc, dtsfc, dusfc, dvsfc
     &,      snowd

      logical,                intent(in) :: flag_cice(im), flag_iter(im)

!  ---  outputs:
      real (kind=kind_phys), dimension(im), intent(inout) :: qsurf,     &
     &                                  cmm, chh, evap, hflx, stress
     &,                                 weasd, snwdph, ep
!
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

!  ---  locals:

      real (kind=kind_phys) :: rho, tem

      real(kind=kind_phys) :: cpinv, hvapi, elocp

      integer :: i

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0
!
      if (.not. cplflx) return
!
      cpinv = one / cp
      hvapi = one / hvap
      elocp = hvap/cp
!
      do i = 1, im
        if (flag_cice(i) .and. flag_iter(i)) then

          rho    = prsl1(i)                                             &
     &      / (rd * t1(i) * (one + rvrdm1*max(q1(i), 1.0e-8_kind_phys)))

          cmm(i) = wind(i) * cm(i)
          chh(i) = wind(i) * ch(i) * rho

          qsurf(i)  = q1(i) + dqsfc(i) / (hvap*chh(i))
          tem       = one / rho
          hflx(i)   = dtsfc(i) * tem * cpinv
          evap(i)   = dqsfc(i) * tem * hvapi
          stress(i) = sqrt(dusfc(i)*dusfc(i) + dvsfc(i)*dvsfc(i)) * tem

          snwdph(i) = snowd(i)  * 1000.0_kind_phys
          weasd(i)  = snwdph(i) * 0.33_kind_phys

!         weasd(i)  = snowd(i) * 1000.0_kind_phys
!         snwdph(i) = weasd(i) * dsi           ! snow depth in mm

          ep(i)     = evap(i)
        endif
      enddo
 
      return
!-----------------------------------
      end subroutine sfc_cice_run
!-----------------------------------

!> @}
      end module sfc_cice

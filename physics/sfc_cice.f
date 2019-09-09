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


!!      use physcons, only : hvap => con_hvap,  cp => con_cp,             &
!!    &                     rvrdm1 => con_fvirt, rd => con_rd
!
!-----------------------------------
      subroutine sfc_cice_run                                           &
     &     ( im, cplflx, cplchm, hvap, cp, rvrdm1, rd,                  & ! ---  inputs:
     &       u1, v1, t1, q1, cm, ch, prsl1, prslki,                     &
     &       flag_cice, ddvel, flag_iter, dqsfc, dtsfc,                 &
     &       qsurf, cmm, chh, evap, hflx,                               & ! ---  outputs:
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
!          ( im, u1, v1, t1, q1, cm, ch, prsl1, prslki,                 !
!            islimsk, ddvel, flag_iter, dqsfc, dtsfc,                   !
!       outputs:                                                        !
!            qsurf, cmm, chh, evap, hflx)                               !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:
!     im, - integer, horiz dimension
!     u1, v1   - real, u/v component of surface layer wind
!     t1       - real, surface layer mean temperature ( k )
!     q1       - real, surface layer mean specific humidity
!     cm       - real, surface exchange coeff for momentum (m/s)
!     ch       - real, surface exchange coeff heat & moisture(m/s)
!     prsl1    - real, surface layer mean pressure
!     prslki   - real, ?
!     islimsk  - integer, sea/land/ice mask
!     ddvel    - real, ?
!     flag_iter- logical
!     dqsfc    - real, latent heat flux
!     dtsfc    - real, sensible heat flux
!  outputs:
!     qsurf    - real, specific humidity at sfc
!     cmm      - real, ?
!     chh      - real, ?
!     evap     - real, evaperation from latent heat
!     hflx     - real, sensible heat
!  ====================    end of description    =====================  !
!
!
      use machine , only : kind_phys
      implicit none


      real (kind=kind_phys), intent(in) :: hvap, cp, rvrdm1, rd

!  ---  inputs:
      integer, intent(in) :: im
      logical, intent(in) :: cplflx
      logical, intent(in) :: cplchm

      real (kind=kind_phys), dimension(im), intent(in) :: u1, v1,       &
     &       t1, q1, cm, ch, prsl1, prslki, ddvel, dqsfc, dtsfc

      logical, dimension(im), intent(in) :: flag_cice

      logical, intent(in) :: flag_iter(im)

!  ---  outputs:
      real (kind=kind_phys), dimension(im), intent(out) :: qsurf,       &
     &       cmm, chh, evap, hflx
!
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

!  ---  locals:
      real (kind=kind_phys), dimension(im) :: q0, rch, rho, tv1, wind

      real (kind=kind_phys) :: tem

      real(kind=kind_phys) :: cpinv, hvapi, elocp

      integer :: i

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0
!
      if ((.not. cplflx) .and. (.not.cplchm)) then
         return
      endif
!
      cpinv = 1.0/cp
      hvapi = 1.0/hvap
      elocp = hvap/cp
!
      do i = 1, im
        if (flag_cice(i) .and. flag_iter(i)) then

          wind(i)   = sqrt(u1(i)*u1(i) + v1(i)*v1(i))                   &
     &              + max(0.0, min(ddvel(i), 30.0))
          wind(i)   = max(wind(i), 1.0)

          q0(i)     = max(q1(i), 1.0e-8)
          tv1(i)    = t1(i) * (1.0 + rvrdm1*q0(i))
          rho(i)    = prsl1(i) / (rd*tv1(i))

          cmm(i) = cm(i)  * wind(i)
          chh(i) = rho(i) * ch(i) * wind(i)
          rch(i) = chh(i) * cp

          qsurf(i) = q1(i) + dqsfc(i) / (elocp*rch(i))
          tem     = 1.0 / rho(i)
          hflx(i) = dtsfc(i) * tem * cpinv
          evap(i) = dqsfc(i) * tem * hvapi
        endif
      enddo
 
      return
!-----------------------------------
      end subroutine sfc_cice_run
!-----------------------------------

!> @}
      end module sfc_cice

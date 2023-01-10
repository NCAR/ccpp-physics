!>  \file sfc_land.f
!!  This file contains the code for coupling to land component

!> This module contains the CCPP-compliant GFS land post
!! interstitial codes, which returns updated surface
!! properties such as latent heat and sensible heat 
!! provided by the component version of land model

!> This module contains the CCPP-compliant GFS land scheme.
      module sfc_land

      contains

!> \defgroup sfc_land for coupling to land
!! @{
!!  \section diagram Calling Hierarchy Diagram
!!  \section intraphysics Intraphysics Communication
!!
!> \brief Brief description of the subroutine
!!
!! \section arg_table_sfc_land_run Arguments
!! \htmlinclude sfc_land_run.html
!!

!!
!!  \section general General Algorithm
!!  \section detailed Detailed Algorithm
!!  @{


!!      use physcons, only : hvap => con_hvap,  cp => con_cp,           &
!!    &                     rvrdm1 => con_fvirt, rd => con_rd
!
!-----------------------------------
      subroutine sfc_land_run                                           &
!  ---  inputs:
     &     ( im, cpllnd, cpllnd2atm, flag_iter, dry,                    &
     &       sncovr1_lnd, qsurf_lnd, evap_lnd, hflx_lnd,                &
     &       ep_lnd, t2mmp_lnd, q2mp_lnd,                               &
!  ---  outputs:
     &       sncovr1, qsurf, evap, hflx, ep, t2mmp, q2mp,               &
     &       errmsg, errflg, naux2d, aux2d
     &     )

! ===================================================================== !
!  description:                                                         !
!  Dec 2022  --  Ufuk Turuncoglu created for coupling to land           !
!                                                                       !
!  usage:                                                               !
!                                                                       !
!    call sfc_land                                                      !
!       inputs:                                                         !
!          ( im, cpllnd, cpllnd2atm, flag_iter, dry,                    !
!            sncovr1_lnd, qsurf_lnd, evap_lnd, hflx_lnd,                !
!            ep_lnd, t2mmp_lnd, q2mp_lnd,                               !
!       outputs:                                                        !
!            sncovr1, qsurf, evap, hflx, ep, t2mmp, q2mp,               !
!            errmsg, errflg)                                            !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:
!     im          - integer, horiz dimension
!     cpllnd      - logical, flag for land coupling
!     cpllnd2atm  - logical, flag for land coupling (lnd->atm)
!     flag_iter   - logical, flag for iteration
!     dry         - logical, eq T if a point with any land 
!     sncovr1_lnd - real   , surface snow area fraction 
!     qsurf_lnd   - real   , specific humidity at sfc      
!     evap_lnd    - real   , evaporation from latent heat
!     hflx_lnd    - real   , sensible heat
!     ep_lnd      - real   , surface upward potential latent heat flux 
!     t2mmp_lnd   - real   , 2m temperature 
!     q2mp_lnd    - real   , 2m specific humidity 
!  outputs:
!     sncovr1     - real   , snow cover over land
!     qsurf       - real   , specific humidity at sfc
!     evap        - real   , evaporation from latent heat
!     hflx        - real   , sensible heat
!     ep          - real   , potential evaporation 
!     t2mmp       - real   , temperature at 2m 
!     q2mp        - real   , specific humidity at 2m 
!  ====================    end of description    =====================  !
!
!
      use machine , only : kind_phys
      implicit none

!  ---  inputs:
      integer, intent(in) :: im
      logical, intent(in) :: cpllnd, cpllnd2atm
      logical, dimension(:), intent(in) :: flag_iter
      logical, dimension(:), intent(in) :: dry

      real (kind=kind_phys), dimension(:), intent(in) ::                &
     &       sncovr1_lnd, qsurf_lnd, evap_lnd, hflx_lnd, ep_lnd,        &
     &       t2mmp_lnd, q2mp_lnd

!  ---  outputs:
      real (kind=kind_phys), dimension(:), intent(out) ::               &
     &       sncovr1, qsurf, evap, hflx, ep, t2mmp, q2mp
!
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      integer, intent(in) :: naux2d
      real(kind_phys), intent(out) :: aux2d(:,:)

!  ---  locals:

      integer :: i

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0
!
      if (.not. cpllnd2atm) return
!
      do i = 1, im
        !if (flag_iter(i) .and. dry(i)) then
        !if (dry(i)) then
           sncovr1(i) = sncovr1_lnd(i)
           qsurf(i)   = qsurf_lnd(i)
           hflx(i)    = hflx_lnd(i)
           evap(i)    = evap_lnd(i)
           ep(i)      = ep_lnd(i)
           t2mmp(i)   = t2mmp_lnd(i)
           q2mp(i)    = q2mp_lnd(i)
        !end if
      enddo

      aux2d(:,1) = dry(:) !sncovr1(:)
      aux2d(:,2) = qsurf(:)
      aux2d(:,3) = hflx(:)
      aux2d(:,4) = evap(:)
      aux2d(:,5) = ep(:)
      aux2d(:,6) = qsurf_lnd(:) !t2mmp(:)
      aux2d(:,7) = q2mp(:)
 
      return
!-----------------------------------
      end subroutine sfc_land_run
!-----------------------------------

!> @}
      end module sfc_land

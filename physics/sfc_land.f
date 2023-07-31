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

!
!-----------------------------------
      subroutine sfc_land_run                                           &
!  ---  inputs:
     &     ( im, cpllnd, cpllnd2atm, flag_iter, dry,                    &
     &       sncovr1_lnd, qsurf_lnd, evap_lnd, hflx_lnd,                &
     &       ep_lnd, t2mmp_lnd, q2mp_lnd, gflux_lnd,                    &
     &       runoff_lnd, drain_lnd, cmm_lnd, chh_lnd, zvfun_lnd,        &
!  ---  outputs:
     &       sncovr1, qsurf, evap, hflx, ep, t2mmp, q2mp,               &
     &       gflux, runoff, drain, cmm, chh, zvfun,                     &
     &       errmsg, errflg
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
!            ep_lnd, t2mmp_lnd, q2mp_lnd, gflux_lnd,                    !
!            runoff_lnd, drain_lnd, cmm_lnd, chh_lnd,                   !
!            zvfun_lnd,                                                 !
!       outputs:                                                        !
!            sncovr1, qsurf, evap, hflx, ep, t2mmp, q2mp,               !
!            gflux, runoff, drain, cmm, chh, zvfun,                     !
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
!     gflux_lnd   - real   , soil heat flux over land
!     runoff_lnd  - real   , surface runoff
!     drain_lnd   - real   , subsurface runoff  
!     cmm_lnd     - real   , surface drag wind speed for momentum
!     chh_lnd     - real   , surface drag mass flux for heat and moisture
!     zvfun_lnd   - real   , function of surface roughness length and green vegetation fraction
!  outputs:
!     sncovr1     - real   , snow cover over land
!     qsurf       - real   , specific humidity at sfc
!     evap        - real   , evaporation from latent heat
!     hflx        - real   , sensible heat
!     ep          - real   , potential evaporation 
!     t2mmp       - real   , temperature at 2m 
!     q2mp        - real   , specific humidity at 2m
!     gflux       - real   , soil heat flux over land
!     runoff      - real   , surface runoff
!     drain       - real   , subsurface runoff
!     cmm         - real   , surface drag wind speed for momentum
!     chh         - real   , surface drag mass flux for heat and moisture
!     zvfun       - real   , function of surface roughness length and green vegetation fraction
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
     &       t2mmp_lnd, q2mp_lnd, gflux_lnd, runoff_lnd, drain_lnd,     &
     &       cmm_lnd, chh_lnd, zvfun_lnd

!  ---  outputs:
      real (kind=kind_phys), dimension(:), intent(inout) ::             &
     &       sncovr1, qsurf, evap, hflx, ep, t2mmp, q2mp, gflux,        &
     &       runoff, drain, cmm, chh, zvfun
!
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

!  ---  locals:

      integer :: i

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0
!
      if (.not. cpllnd2atm) return
!
      do i = 1, im
         sncovr1(i) = sncovr1_lnd(i)
         qsurf(i)   = qsurf_lnd(i)
         hflx(i)    = hflx_lnd(i)
         evap(i)    = evap_lnd(i)
         ep(i)      = ep_lnd(i)
         t2mmp(i)   = t2mmp_lnd(i)
         q2mp(i)    = q2mp_lnd(i)
         gflux(i)   = gflux_lnd(i)
         drain(i)   = drain_lnd(i)
         runoff(i)  = runoff_lnd(i)
         cmm(i)     = cmm_lnd(i)
         chh(i)     = chh_lnd(i)
         zvfun(i)   = zvfun_lnd(i)
      enddo
 
      return
!-----------------------------------
      end subroutine sfc_land_run
!-----------------------------------

!> @}
      end module sfc_land

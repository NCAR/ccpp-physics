!> \file sfc_land.F90
!! This file contains the code for coupling to land component

!> This module contains the CCPP-compliant GFS land post
!! interstitial codes, which returns updated surface
!! properties such as latent heat and sensible heat 
!! provided by the component version of land model

!> This module contains the CCPP-compliant GFS land scheme.
   module sfc_land

   use machine, only : kind_phys
   use funcphys, only : fpvs

   contains

!> \defgroup sfc_land for coupling to land
!! @{
!! \section diagram Calling Hierarchy Diagram
!! \section intraphysics Intraphysics Communication
!!
!> \brief Brief description of the subroutine
!!
!! \section arg_table_sfc_land_run Arguments
!! \htmlinclude sfc_land_run.html
!!

!!
!! \section general General Algorithm
!! \section detailed Detailed Algorithm
!! @{
   subroutine sfc_land_run(im, flag_init, flag_restart,              &
     cpllnd, cpllnd2atm, flag_iter, dry,                             &
     t1, q1, prsl1, prslki, ps, tskin, wind, cm, ch, rd, eps, epsm1, &
     rvrdm1, hvap, cp, sncovr1_lnd, qsurf_lnd,                       &
     evap_lnd, hflx_lnd, ep_lnd, t2mmp_lnd, q2mp_lnd, gflux_lnd,     &
     runoff_lnd, drain_lnd, cmm_lnd, chh_lnd, zvfun_lnd,             &
     sncovr1, qsurf, evap, hflx, ep, t2mmp, q2mp,                    &
     gflux, runoff, drain, cmm, chh, zvfun,                          &
     errmsg, errflg, naux2d, aux2d)

   implicit none

   ! Inputs
   integer             , intent(in) :: im
   logical             , intent(in) :: flag_init
   logical             , intent(in) :: flag_restart
   logical             , intent(in) :: cpllnd
   logical             , intent(in) :: cpllnd2atm
   logical             , intent(in) :: flag_iter(:)
   logical             , intent(in) :: dry(:)
   real(kind=kind_phys), intent(in) :: t1(:)
   real(kind=kind_phys), intent(in) :: q1(:)
   real(kind=kind_phys), intent(in) :: prsl1(:)
   real(kind=kind_phys), intent(in) :: prslki(:)
   real(kind=kind_phys), intent(in) :: ps(:)
   real(kind=kind_phys), intent(in) :: tskin(:)
   real(kind=kind_phys), intent(in) :: wind(:)
   real(kind=kind_phys), intent(in) :: cm(:)
   real(kind=kind_phys), intent(in) :: ch(:)
   real(kind=kind_phys), intent(in) :: rd
   real(kind=kind_phys), intent(in) :: eps
   real(kind=kind_phys), intent(in) :: epsm1
   real(kind=kind_phys), intent(in) :: rvrdm1
   real(kind=kind_phys), intent(in) :: hvap
   real(kind=kind_phys), intent(in) :: cp
   real(kind=kind_phys), intent(in), optional :: sncovr1_lnd(:)
   real(kind=kind_phys), intent(in), optional :: qsurf_lnd(:)
   real(kind=kind_phys), intent(in), optional :: evap_lnd(:)
   real(kind=kind_phys), intent(in), optional :: hflx_lnd(:)
   real(kind=kind_phys), intent(in), optional :: ep_lnd(:)
   real(kind=kind_phys), intent(in), optional :: t2mmp_lnd(:)
   real(kind=kind_phys), intent(in), optional :: q2mp_lnd(:)
   real(kind=kind_phys), intent(in), optional :: gflux_lnd(:)
   real(kind=kind_phys), intent(in), optional :: runoff_lnd(:)
   real(kind=kind_phys), intent(in), optional :: drain_lnd(:)
   real(kind=kind_phys), intent(in), optional :: cmm_lnd(:)
   real(kind=kind_phys), intent(in), optional :: chh_lnd(:)
   real(kind=kind_phys), intent(in), optional :: zvfun_lnd(:)
   ! Inputs/Outputs
   real(kind=kind_phys), intent(inout) :: sncovr1(:)
   real(kind=kind_phys), intent(inout) :: qsurf(:)
   real(kind=kind_phys), intent(inout) :: evap(:)
   real(kind=kind_phys), intent(inout) :: hflx(:)
   real(kind=kind_phys), intent(inout) :: ep(:)
   real(kind=kind_phys), intent(inout), optional :: t2mmp(:)
   real(kind=kind_phys), intent(inout), optional :: q2mp(:)
   real(kind=kind_phys), intent(inout) :: gflux(:)
   real(kind=kind_phys), intent(inout) :: runoff(:)
   real(kind=kind_phys), intent(inout) :: drain(:)
   real(kind=kind_phys), intent(inout) :: cmm(:)
   real(kind=kind_phys), intent(inout) :: chh(:)
   real(kind=kind_phys), intent(inout) :: zvfun(:)
   ! Outputs
   character(len=*)    , intent(out)   :: errmsg
   integer             , intent(out)   :: errflg

   ! Constant parameters
   real(kind=kind_phys), parameter :: &
  &    one  = 1.0_kind_phys, &
  &    zero = 0.0_kind_phys, &
  &    qmin = 1.0e-8_kind_phys

   ! Locals
   integer :: i
   real(kind=kind_phys) :: qss, rch, tem, cpinv, hvapi, elocp
   real(kind=kind_phys), dimension(im) :: rho, q0

   ! Initialize CCPP error handling variables
   errmsg = ''
   errflg = 0

   cpinv = one/cp
   hvapi = one/hvap
   elocp = hvap/cp

   ! Check coupling from component land to atmosphere
   if (.not. cpllnd2atm) return

   ! Check if it is cold or warm run
   if (flag_init .and. .not. flag_restart) then 
      ! Calculate fluxes internally
      do i = 1, im
         if (dry(i)) then
            q0(i) = max(q1(i), qmin)
            rho(i) = prsl1(i)/(rd*t1(i)*(one+rvrdm1*q0(i)))
            qss = fpvs(tskin(i))
            qss = eps*qss/(ps(i)+epsm1*qss)
            rch = rho(i)*cp*ch(i)*wind(i)
            tem = ch(i)*wind(i)
            sncovr1(i) = zero
            qsurf(i) = qss
            hflx(i) = rch*(tskin(i)-t1(i)*prslki(i))
            hflx(i) = hflx(i)*(1.0/rho(i))*cpinv
            evap(i) = elocp*rch*(qss-q0(i))
            ep(i) = evap(i)
            evap(i) = evap(i)*(1.0/rho(i))*hvapi
            t2mmp(i) = tskin(i)
            q2mp(i) = qsurf(i)
            gflux(i) = zero
            drain(i) = zero
            runoff(i) = zero
            cmm(i) = cm(i)*wind(i)
            chh(i) = rho(i)*tem
            zvfun(i) = one
         end if
      enddo
   else
      ! Use fluxes from land component model 
      do i = 1, im
         if (dry(i)) then
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
         end if
      enddo
   endif
 
   end subroutine sfc_land_run

!> @}
   end module sfc_land

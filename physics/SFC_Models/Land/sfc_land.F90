!> \file sfc_land.F90
!! This file contains the code for coupling to land component

!> This module contains the CCPP-compliant GFS land post
!! interstitial codes, which returns updated surface
!! properties such as latent heat and sensible heat 
!! provided by the component version of land model

!> This module contains the CCPP-compliant GFS land scheme.
   module sfc_land

   use machine, only : kind_phys

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
   subroutine sfc_land_run(im, cpllnd, cpllnd2atm, flag_iter, dry,   &
     sncovr1_lnd, qsurf_lnd, evap_lnd, hflx_lnd,                     &
     ep_lnd, t2mmp_lnd, q2mp_lnd, gflux_lnd,                         &
     runoff_lnd, drain_lnd, cmm_lnd, chh_lnd, zvfun_lnd,             &
     sncovr1, qsurf, evap, hflx, ep, t2mmp, q2mp,                    &
     gflux, runoff, drain, cmm, chh, zvfun,                          &
     errmsg, errflg)

   implicit none

   ! Inputs
   integer             , intent(in)    :: im
   logical             , intent(in)    :: cpllnd
   logical             , intent(in)    :: cpllnd2atm
   logical             , intent(in)    :: flag_iter(:)
   logical             , intent(in)    :: dry(:)
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

   ! Locals
   integer :: i

   ! Initialize CCPP error handling variables
   errmsg = ''
   errflg = 0

   ! Check coupling from component land to atmosphere
   if (.not. cpllnd2atm) return

   ! Fill variables
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
 
   end subroutine sfc_land_run

!> @}
   end module sfc_land

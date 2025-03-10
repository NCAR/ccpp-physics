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

!> \brief Brief description of the subroutine
!! \section arg_table_sfc_land_run Arguments
!! \htmlinclude sfc_land_run.html
!!

!!
!! \section general General Algorithm
!! \section detailed Detailed Algorithm
!! @{
   subroutine sfc_land_run(im, flag_init, flag_restart,              &
     cpllnd, cpllnd2atm, flag_iter, dry,                             &
     t1, q1, prsl1, prslki, ps, tskin, wind, cm, ch,                 &
     dlwflx, dswsfc, sfalb, sfcemis,                                 &
     rd, eps, epsm1, rvrdm1, hvap, cp, con_sbc,                      &
     sncovr1_lnd, qsurf_lnd,                                         &
     evap_lnd, hflx_lnd, ep_lnd, t2mmp_lnd, q2mp_lnd, gflux_lnd,     &
     runoff_lnd, drain_lnd, cmm_lnd, chh_lnd, zvfun_lnd, slc,        &
     sncovr1, qsurf, evap, hflx, ep, t2mmp, q2mp,                    &
     gflux, runoff, drain, cmm, chh, zvfun,                          &
     errmsg, errflg)

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
   real(kind=kind_phys), intent(in) :: dlwflx(:)
   real(kind=kind_phys), intent(in) :: dswsfc(:)
   real(kind=kind_phys), intent(in) :: sfalb(:)
   real(kind=kind_phys), intent(in) :: sfcemis(:)
   real(kind=kind_phys), intent(in) :: rd
   real(kind=kind_phys), intent(in) :: eps
   real(kind=kind_phys), intent(in) :: epsm1
   real(kind=kind_phys), intent(in) :: rvrdm1
   real(kind=kind_phys), intent(in) :: hvap
   real(kind=kind_phys), intent(in) :: cp
   real(kind=kind_phys), intent(in) :: con_sbc
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
   real(kind=kind_phys), intent(in), optional :: slc(:,:)
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
  &    qmin = 1.0e-8_kind_phys, &
  & slc_min = 0.05_kind_phys, &   ! estimate dry limit for soil moisture
  & slc_max = 0.50_kind_phys      ! estimate saturated limit for soil moisture

   ! Locals
   integer :: i
   real(kind=kind_phys) :: qss, rch, tem, cpinv, hvapi, elocp
   real(kind=kind_phys) :: available_energy, soil_stress_factor
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
            soil_stress_factor = (slc(i,1)-slc_min)/(slc_max-slc_min)
            soil_stress_factor = min(max(zero,soil_stress_factor),one)
            available_energy = dswsfc(i)*(one-sfalb(i))+dlwflx(i)*sfcemis(i) - &
                               sfcemis(i)*con_sbc*tskin(i)**4
            available_energy = min(max(-200.0,available_energy),1000.0)   ! set some arbitrary limits
            q0(i) = max(q1(i), qmin)
            rho(i) = prsl1(i)/(rd*t1(i)*(one+rvrdm1*q0(i)))
            qss = fpvs(tskin(i))
            qss = eps*qss/(ps(i)+epsm1*qss)
            rch = rho(i)*cp*ch(i)*wind(i)
            tem = ch(i)*wind(i)
            sncovr1(i) = zero
            qsurf(i) = qss
            hflx(i) = rch*(tskin(i)-t1(i)*prslki(i))   ! first guess hflx [W/m2]
            evap(i) = elocp*rch*(qss-q0(i))            ! first guess evap [W/m2]
            evap(i) = evap(i)*soil_stress_factor     ! reduce evap for soil moisture stress
            hflx(i) = min(max(-100.0,hflx(i)),500.0)   ! set some arbitrary limits
            evap(i) = min(max(-100.0,evap(i)),500.0)   ! set some arbitrary limits
            if(evap(i) + hflx(i) /= zero) then
              hflx(i) = available_energy * hflx(i) / (abs(evap(i)) + abs(hflx(i)))
              evap(i) = available_energy * evap(i) / (abs(evap(i)) + abs(hflx(i)))
            else
              hflx(i) = zero
              evap(i) = zero
            end if
            hflx(i) = min(max(-100.0,hflx(i)),500.0)   ! set some arbitrary limits
            evap(i) = min(max(-100.0,evap(i)),500.0)   ! set some arbitrary limits
            hflx(i) = hflx(i)*(1.0/rho(i))*cpinv       ! convert to expected units
            ep(i) = evap(i)
            evap(i) = evap(i)*(1.0/rho(i))*hvapi       ! convert to expected units
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

   end module sfc_land

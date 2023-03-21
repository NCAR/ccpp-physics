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
      subroutine sfc_diag_run (im,xlat_d,xlon_d,                        &
     &                    lsm,lsm_ruc,grav,cp,eps,epsm1,rocp,con_karman,&
     &                    wet,shflx,cdq,wind,                           &
     &                    zf,ps,u1,v1,t1,q1,prslki,evap,fm,fh,fm10,fh2, &
     &                    ust,tskin,qsurf,thsfc_loc,diag_flux,diag_log, &
     &                    f10m,u10m,v10m,t2m,q2m,dpt2m,errmsg,errflg    &
     &                   )
!
      use machine , only : kind_phys
      use funcphys, only : fpvs
      implicit none
!
      integer, intent(in) :: im, lsm, lsm_ruc
      logical, intent(in) :: thsfc_loc  ! Flag for reference pot. temp.
      logical, intent(in) :: diag_flux  ! Flag for flux method in 2-m diagnostics
      logical, intent(in) :: diag_log   ! Flag for 2-m log diagnostics under stable conditions
      real(kind=kind_phys), intent(in) :: grav,cp,eps,epsm1,rocp
      real(kind=kind_phys), intent(in) :: con_karman
      real(kind=kind_phys), dimension(:), intent( in) ::                &
     &                      zf, ps, u1, v1, t1, q1, ust, tskin, wet,    &
     &                      qsurf, prslki, evap, fm, fh, fm10, fh2,     &
     &                      shflx, cdq, wind, xlat_d, xlon_d
      real(kind=kind_phys), dimension(:), intent(out) ::                &
     &                       f10m, u10m, v10m, t2m, q2m, dpt2m
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
!
!     locals
!
      logical :: debug_print
      real(kind=kind_phys), parameter :: qmin=1.0e-8
      real(kind=kind_phys) :: q1c, qv, tem, qv1, th2m, x2m, rho
      real(kind=kind_phys) :: dT, dQ, qsfcmr, qsfcprox, ff, fac, dz1
      real(kind=kind_phys) :: t2_alt, q2_alt
      real(kind=kind_phys) :: thcon, cqs, chs, chs2, cqs2
      real(kind=kind_phys) :: testptlat, testptlon
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

      !--
      testptlat = 35.3 !41.02 !42.05 !39.0 !74.12 !29.5 
      testptlon = 273.0 !284.50 !286.75 !280.6 !164.0 !283.0 
      !--
      debug_print = .false.
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

        thcon    = (1.e5/ps(i))**rocp    
       !-- make sure 1st level q is not higher than saturated value
        qss    = fpvs(t1(i))
        qss    = eps * qss / (ps(i) + epsm1 * qss)
        q1c    = min(q1(i),qss) ! lev 1 spec. humidity

        qv1    = q1c / (1. - q1c) ! lev 1 mixing ratio
        qsfcmr   = qsurf(i)/(1. - qsurf(i)) ! surface mixing ratio
        chs = cdq(i) * wind(i)
        cqs = chs
        chs2 = ust(i)*con_karman/fh2(i)
        cqs2 = chs2
        qsfcprox = max(qmin,qv1 + evap(i)/cqs) ! surface mix. ratio computed from the flux

        if(.not. diag_flux) then
        !-- original method
          if(lsm /= lsm_ruc) then
            if(thsfc_loc) then ! Use local potential temperature
              t2m(i)=tskin(i)*wrk + t1(i)*prslki(i)*fhi - (grav+grav)/cp
            else ! Use potential temperature referenced to 1000 hPa
              t2m(i) = tskin(i)*wrk + t1(i)*fhi - (grav+grav)/cp
            endif
            if(evap(i) >= 0.) then !  for evaporation>0, use inferred qsurf to deduce q2m
              q2m(i) = qsurf(i)*wrk + max(qmin,q1c)*fhi
            else                   !  for dew formation, use saturated q at tskin
              qss    = fpvs(tskin(i))
              qss    = eps * qss/(ps(i) + epsm1 * qss)
              q2m(i) = qss*wrk + max(qmin,q1c)*fhi
            endif
          else
            t2m(i)  = tskin(i)*wrk + t1(i)*fhi - (grav+grav)/cp
            q2m(i) = qsurf(i)*wrk + max(qmin,q1c)*fhi
          endif ! RUC lsm

        else
        !-- flux method
            th2m = tskin(i)*thcon - shflx(i)/chs2
            t2m(i) = th2m/thcon

            x2m = max(qmin,qsfcprox - evap(i)/cqs2) ! mix. ratio
            q2m(i) = x2m/(1. + x2m) ! spec. humidity
        endif ! flux method

        if(diag_log) then
        !-- Alternative logarithmic diagnostics:
          dT = t1(i) - tskin(i)
          dQ = qv1 - qsfcmr
          dz1= zf(i) ! level of atm. forcing
          IF (dT > 0.) THEN
            ff  = MIN(MAX(1.-dT/10.,0.01), 1.0)
            !for now, set zt = 0.05
            fac = LOG((2.  + .05)/(0.05 + ff))/                         &
     &          LOG((dz1 + .05)/(0.05 + ff))
            T2_alt = tskin(i) + fac * dT
          ELSE
          !no alternatives (yet) for unstable conditions
            T2_alt = t2m(i)
          ENDIF

          IF (dQ > 0.) THEN
            ff  = MIN(MAX(1.-dQ/0.003,0.01), 1.0)
            !-- for now, set zt = 0.05
            fac = LOG((2.  + .05)/(0.05 + ff))/                         &
     &            LOG((dz1 + .05)/(0.05 + ff))
            Q2_alt = qsfcmr + fac * dQ ! mix. ratio
            Q2_alt = Q2_alt/(1. + Q2_alt) ! spec. humidity
          ELSE
            !no alternatives (yet) for unstable conditions
            Q2_alt = q2m(i)
          ENDIF
          !-- Note: use of alternative diagnostics will make 
          !   it cooler and drier with stable stratification
          t2m(i) = T2_alt  
          q2m(i) = Q2_alt
        endif ! log method for stable regime

        !-- check that T2m values lie in the range between tskin and t1
        x2m     = max(min(tskin(i),t1(i)) , t2m(i))
        t2m(i)  = min(max(tskin(i),t1(i)) , x2m)
        !-- check that Q2m values lie in the range between qsurf and q1
        x2m    = max(min(qsurf(i),q1c) , q2m(i))
        q2m(i) = min(max(qsurf(i),q1c) , x2m)

        !-- make sure q2m is not oversaturated
        qss    = fpvs(t2m(i))
        qss    = eps * qss/(ps(i) + epsm1 * qss)
        q2m(i) = min(q2m(i),qss)

        if(diag_flux) then
        !-- from WRF, applied in HRRR - Jimy Dudhia
        !   Limit Q2m diagnostic to no more than 5 percent higher than lowest level value
        !   This prevents unrealistic values when QFX is not mostly surface
        !   flux because calculation is based on surface flux only.
        !   Problems occurred in transition periods and weak winds and vegetation source
          q2m(i) = min(q2m(i),1.05*q1c) ! works if qsurf > q1c, evaporation
        endif


        !-- Compute dew point, using vapor pressure
        qv  = max(qmin,(q2m(i)/(1.-q2m(i))))
        tem = max(ps(i) * qv/( eps - epsm1 *qv), 1.e-8)
        dpt2m(i) = 243.5/( ( 17.67 / log(tem/611.2) ) - 1.) + 273.14
        dpt2m(i) = min(dpt2m(i),t2m(i))
       

        if (debug_print) then
        !-- diagnostics for a test point with known lat/lon
          if (abs(xlat_d(i)-testptlat).lt.0.2 .and.                     &
     &      abs(xlon_d(i)-testptlon).lt.0.2)then
            print 100,'(ruc_lsm_diag)  i=',i,                           &
     &      '  lat,lon=',xlat_d(i),xlon_d(i),'zf ',zf(i),               &
     &      'tskin ',tskin(i),'t2m ',t2m(i),'t1',t1(i),'shflx',shflx(i),&
     &      'qsurf ',qsurf(i),'qsfcprox ',qsfcprox,'q2m ',q2m(i),       &
     &      'q1 ',q1(i),'evap ',evap(i),'dpt2m ',dpt2m(i),              &
     &      'chs2 ',chs2,'cqs2 ',cqs2,'cqs ',cqs,'cdq',cdq(i)
          endif
        endif
 100    format (";;; ",a,i4,a,2f14.7/(4(a10,'='es11.4)))

      enddo

      return
      end subroutine sfc_diag_run
!> @}

      end module sfc_diag
!> @}

!>\file wv_saturation.F
!!This file contains common block and statement functions for saturation vapor pressure
!! look-up procedure, J. J. Hack, February 1990
!
! $Id: wv_saturation.F90,v 1.1.12.2.10.1 2014-04-14 16:04:56 dbarahon Exp $
!
!>\ingroup mg2mg3
!>\defgroup wv_saturation_mod Morrison-Gettelman MP wv_saturation Module
!! This module contain some utility functions for saturation vapor pressure.
      module wv_saturation
#ifdef GEOS5
       use MAPL_ConstantsMod, kp => MAPL_R8
#endif
#ifdef NEMS_GSM
       use funcphys,  only : fpvsl, fpvsi, fpvs    ! saturation vapor pressure for water & ice
       use machine,   only : kp => kind_phys
#endif

!++jtb (comm out)



!--jtb

       implicit none
       private
       save
!
! Public interfaces
!
       public gestbl
       public estblf
       public aqsat
       public aqsatd
       public vqsatd
       public fqsatd
       public qsat_water
       public vqsat_water
       public qsat_ice
       public vqsat_ice
       public vqsatd_water
       public aqsat_water
       public vqsatd2_water
       public vqsatd2_water_single
       public vqsatd2_ice_single
       public vqsatd2
       public vqsatd2_single
       public polysvp
!
! Data used by cldwat
!
       public hlatv, tmin, hlatf, rgasv, pcf, cp, epsqs, ttrice
!
! Data
!
       integer plenest
       parameter (plenest=250)
!
! Table of saturation vapor pressure values es from tmin degrees
! to tmax+1 degrees k in one degree increments.  ttrice defines the
! transition region where es is a combination of ice & water values
!
       real(kp) estbl(plenest)
       real(kp) tmin
       real(kp) tmax
       real(kp) ttrice
       real(kp) pcf(6)
       real(kp) epsqs
       real(kp) rgasv
       real(kp) hlatf
       real(kp) hlatv
       real(kp) cp
       real(kp) tmelt
       logical icephs

       integer, parameter :: iulog=6

       contains

       real(kp) function estblf( td )
!
! Saturation vapor pressure table lookup
!
       real(kp), intent(in) :: td
!
       real(kp) :: e
       real(kp) :: ai
       integer :: i
!
       e = max(min(td,tmax),tmin)
       i = int(e-tmin)+1
       ai = aint(e-tmin)
       estblf = (tmin+ai-e+1._kp)* estbl(i)-(tmin+ai-e)* estbl(i+1)
       end function estblf

!>\ingroup wv_saturation_mod
!! This subroutine builds saturation vapor pressure table for later
!! procedure.
!!
!! Method:
!! uses Goff and Gratch (1946) relationships to generate the table
!! according to a set of free paramters defined below. Auxiliary 
!! routines are also included for making rapid estimates (well with 1%)
!! of both es and d(es)/dt for the particular table configuration.
!>\author J. Hack
      subroutine gestbl(tmn ,tmx ,trice ,ip ,epsil , latvap ,latice ,   &
     &                  rh2o ,cpair ,tmeltx )
!------------------------------Arguments--------------------------------
!
! Input arguments
!
       real(kp), intent(in) :: tmn
       real(kp), intent(in) :: tmx
       real(kp), intent(in) :: epsil
       real(kp), intent(in) :: trice
       real(kp), intent(in) :: latvap
       real(kp), intent(in) :: latice
       real(kp), intent(in) :: rh2o
       real(kp), intent(in) :: cpair
       real(kp), intent(in) :: tmeltx
!
!---------------------------Local variables-----------------------------
!
       real(kp) t
       integer n
       integer lentbl
       integer itype
!            1 -> ice phase, no transition
!           -x -> ice phase, x degree transition
       logical ip
!
!-----------------------------------------------------------------------
!
! Set es table parameters
!
       tmin   = tmn
       tmax   = tmx
       ttrice = trice
       icephs = ip
!
! Set physical constants required for es calculation
!
       epsqs = epsil
       hlatv = latvap
       hlatf = latice
       rgasv = rh2o
       cp    = cpair
       tmelt = tmeltx
!
       lentbl = INT(tmax-tmin+2.000001_kp)
       if (lentbl .gt. plenest) then



       write(*,*) "AHHH wv_sat"
       STOP
       end if
!
! Begin building es table.
! Check whether ice phase requested.
! If so, set appropriate transition range for temperature
!
       if (icephs) then
         if (ttrice /= 0.0_kp) then
           itype = -ttrice
         else
           itype = 1
         end if
       else
         itype = 0
       end if
!
       t = tmin - 1.0_kp
       do n=1,lentbl
         t = t + 1.0_kp
         call gffgch(t,estbl(n),tmelt,itype)
       end do
!
       do n=lentbl+1,plenest
         estbl(n) = -99999.0_kp
       end do
!
! Table complete -- Set coefficients for polynomial approximation of
! difference between saturation vapor press over water and saturation
! pressure over ice for -ttrice < t < 0 (degrees C). NOTE: polynomial
! is valid in the range -40 < t < 0 (degrees C).
!
!                  --- Degree 5 approximation ---
!
       pcf(1) = 5.04469588506e-01_kp
       pcf(2) = -5.47288442819e+00_kp
       pcf(3) = -3.67471858735e-01_kp
       pcf(4) = -8.95963532403e-03_kp
       pcf(5) = -7.78053686625e-05_kp
!
!                  --- Degree 6 approximation ---
!
!-----pcf(1) =  7.63285250063e-02
!-----pcf(2) = -5.86048427932e+00
!-----pcf(3) = -4.38660831780e-01
!-----pcf(4) = -1.37898276415e-02
!-----pcf(5) = -2.14444472424e-04
!-----pcf(6) = -1.36639103771e-06
!

!++jtb (comment out)
!   if (masterproc) then
!      !!write(iulog,*)' *** SATURATION VAPOR PRESSURE TABLE COMPLETED ***'
!   end if
!--jtb

       return
!
9000   format('GESTBL: FATAL ERROR *********************************    &
     &',/, ' TMAX AND TMIN REQUIRE A LARGER DIMENSION ON THE LENGTH',   &
     & ' OF THE SATURATION VAPOR PRESSURE TABLE ESTBL(PLENEST)',/,      &
     & ' TMAX, TMIN, AND PLENEST => ', 2f7.2, i3)
!
      end subroutine gestbl

!>\ingroup wv_saturation_mod
!! This subroutine is the utility procedure to look up and returen
!! saturation vapor pressure from precomputed table, calculate and
!! return saturation specific humidity (\f$g g^{-1}\f$), for input 
!! arrays of temperature and pressure (dimensioned ii,kk).
!> \author J. Hack
      subroutine aqsat(t ,p ,es ,qs ,ii , ilen ,kk ,kstart ,kend )
!------------------------------Arguments--------------------------------
!
! Input arguments
!
       integer,  intent(in) :: ii
       integer,  intent(in) :: kk
       integer,  intent(in) :: ilen
       integer,  intent(in) :: kstart
       integer,  intent(in) :: kend
       real(kp), intent(in) :: t(ii,kk)
       real(kp), intent(in) :: p(ii,kk)
!
! Output arguments
!
       real(kp), intent(out) :: es(ii,kk)
       real(kp), intent(out) :: qs(ii,kk)
!
!---------------------------Local workspace-----------------------------
!
       real(kp) omeps
       integer i, k
!
!-----------------------------------------------------------------------
!
       omeps = 1.0_kp - epsqs
       do k=kstart,kend
         do i=1,ilen
           es(i,k) = min(estblf(t(i,k)),p(i,k))
!
! Saturation specific humidity
!
           qs(i,k) = min(1.0_kp, epsqs*es(i,k)/(p(i,k)-omeps*es(i,k)))
!
! The following check is to avoid the generation of negative values
! that can occur in the upper stratosphere and mesosphere
!
!      if (qs(i,k) < 0.0_kp) then
!      qs(i,k) = 1.0_kp
!      es(i,k) = p(i,k)
!      end if

         end do
       end do
!
       return
      end subroutine aqsat

!++xl
!>\ingroup wv_saturation_mod
!! This subroutine includes the utility procedure to look up and return
!! saturation vapor pressure from precomputed table, calculate and return
!! saturation specific humidity (g/g), for input arrays of temperature and 
!! pressure (dimensiond ii, kk). This routine is useful for evaluating only
!! a selected region in the vertical. 
!>\author J. Hack
      subroutine aqsat_water(t, p, es, qs, ii, ilen, kk, kstart,kend)
!------------------------------Arguments--------------------------------
!
! Input arguments
!
       integer,  intent(in) :: ii
       integer,  intent(in) :: kk
       integer,  intent(in) :: ilen
       integer,  intent(in) :: kstart
       integer,  intent(in) :: kend
       real(kp), intent(in) :: t(ii,kk)
       real(kp), intent(in) :: p(ii,kk)
!
! Output arguments
!
       real(kp), intent(out) :: es(ii,kk)
       real(kp), intent(out) :: qs(ii,kk)
!
!---------------------------Local workspace-----------------------------
!
       real(kp) omeps
       integer i, k
!
!-----------------------------------------------------------------------
!
       omeps = 1.0_kp - epsqs
       do k=kstart,kend
         do i=1,ilen
!          es(i,k) = estblf(t(i,k))
#ifdef GEOS5
           es(i,k) = min(polysvp(t(i,k),0), p(i,k))
#endif
#ifdef NEMS_GSM
           es(i,k) = min(fpvsl(t(i,k)), p(i,k))
#endif
!
! Saturation specific humidity
!
           qs(i,k) = min(1.0_kp, epsqs*es(i,k)/(p(i,k)-omeps*es(i,k)))
!
! The following check is to avoid the generation of negative values
! that can occur in the upper stratosphere and mesosphere
!
!          if (qs(i,k) < 0.0_kp) then
!            qs(i,k) = 1.0_kp
!            es(i,k) = p(i,k)
!          end if
         end do
       end do
!
       return
      end subroutine aqsat_water
!--xl

!>\ingroup wv_saturation_mod
!! This subroutine include the utility procedure to look up and returen saturation
!! vapor pressure from precomputed table, calculate and return saturation 
!! specific humidity (g/g)
!!
!! Method:
!! Differs from aqsat by also calculating and returning
!! gamma (l/cp)*(d(qsat)/dT)
!! Input arrays temperature and pressure (dimensioned ii,kk).
!!
!! \author J. Hack
      subroutine aqsatd(t, p, es, qs, gam, ii, ilen, kk, kstart, kend)
!------------------------------Arguments--------------------------------
!
! Input arguments
!
       integer,  intent(in) :: ii
       integer,  intent(in) :: ilen
       integer,  intent(in) :: kk
       integer,  intent(in) :: kstart
       integer,  intent(in) :: kend

       real(kp), intent(in) :: t(ii,kk)
       real(kp), intent(in) :: p(ii,kk)

!
! Output arguments
!
       real(kp), intent(out) :: es(ii,kk)
       real(kp), intent(out) :: qs(ii,kk)
       real(kp), intent(out) :: gam(ii,kk)
!
!---------------------------Local workspace-----------------------------
!
       logical lflg
       integer i
       integer k
       real(kp) omeps
       real(kp) trinv
       real(kp) tc
       real(kp) weight
       real(kp) hltalt
       real(kp) hlatsb
       real(kp) hlatvp
       real(kp) tterm
       real(kp) desdt
!
!-----------------------------------------------------------------------
!
       omeps = 1.0_kp - epsqs
       do k=kstart,kend
         do i=1,ilen
           es(i,k) = min(p(i,k), estblf(t(i,k)))
!
! Saturation specific humidity
!
           qs(i,k) = min(1.0_kp, epsqs*es(i,k)/(p(i,k)-omeps*es(i,k)))
!
! The following check is to avoid the generation of negative qs
! values which can occur in the upper stratosphere and mesosphere
!
!
!          if (qs(i,k) < 0.0_kp) then
!            qs(i,k) = 1.0_kp
!            es(i,k) = p(i,k)
!          end if
         end do
       end do
!
! "generalized" analytic expression for t derivative of es
! accurate to within 1 percent for 173.16 < t < 373.16
!
       trinv = 0.0_kp
       if ((.not. icephs) .or. (ttrice == 0.0_kp)) go to 10
       trinv = 1.0_kp/ttrice
!
       do k=kstart,kend
         do i=1,ilen
!
! Weighting of hlat accounts for transition from water to ice
! polynomial expression approximates difference between es over
! water and es over ice from 0 to -ttrice (C) (min of ttrice is
! -40): required for accurate estimate of es derivative in transition
! range from ice to water also accounting for change of hlatv with t
! above freezing where constant slope is given by -2369 j/(kg c) =cpv - cw
!
           tc     = t(i,k) - tmelt
           lflg   = (tc >= -ttrice .and. tc < 0.0_kp)
           weight = min(-tc*trinv,1.0_kp)
           hlatsb = hlatv + weight*hlatf
           hlatvp = hlatv - 2369.0_kp*tc
           if (t(i,k) < tmelt) then
             hltalt = hlatsb
           else
             hltalt = hlatvp
           end if
           if (lflg) then
             tterm = pcf(1) + tc*(pcf(2) + tc*(pcf(3) + tc*(pcf(4)      &
     &                      + tc*pcf(5))))
           else
             tterm = 0.0_kp
           end if
             desdt = hltalt*es(i,k)/(rgasv*t(i,k)*t(i,k)) + tterm*trinv
             gam(i,k) = hltalt*qs(i,k)*p(i,k)*desdt/(cp*es(i,k)*(p(i,k) &
     &                  - omeps*es(i,k)))
           if (qs(i,k) == 1.0_kp) gam(i,k) = 0.0_kp
         end do
       end do
!
       go to 20
!
! No icephs or water to ice transition
!
10     do k=kstart,kend
         do i=1,ilen
!
! Account for change of hlatv with t above freezing where
! constant slope is given by -2369 j/(kg c) = cpv - cw
!
           hlatvp = hlatv - 2369.0_kp*(t(i,k)-tmelt)
           if (icephs) then
             hlatsb = hlatv + hlatf
           else
             hlatsb = hlatv
           end if
           if (t(i,k) < tmelt) then
             hltalt = hlatsb
           else
             hltalt = hlatvp
           end if
           desdt    = hltalt*es(i,k)/(rgasv*t(i,k)*t(i,k))
           gam(i,k) = hltalt*qs(i,k)*p(i,k)*desdt/(cp*es(i,k)*(p(i,k)   &
     &                  - omeps*es(i,k)))
           if (qs(i,k) == 1.0_kp) gam(i,k) = 0.0_kp
         end do
       end do
!
20     return
      end subroutine aqsatd

!>\ingroup wv_saturation_mod
!! This subroutine is the utility procedure to look up and return 
!! saturation vapor pressure from precomputed table, calculated the 
!! return saturation specific humidity (g/g). and calculate and 
!! return gamma (1/cp)*(d(qsat)/dT). The same function as qsatd,
!! but operates on vectors of temperature and pressure
!>\author J. Hack
      subroutine vqsatd(t ,p ,es ,qs ,gam , len )
!------------------------------Arguments--------------------------------
!
! Input arguments
!
       integer,  intent(in) :: len
       real(kp), intent(in) :: t(len)
       real(kp), intent(in) :: p(len)
!
! Output arguments
!
       real(kp), intent(out) :: es(len)
       real(kp), intent(out) :: qs(len)
       real(kp), intent(out) :: gam(len)
!
!--------------------------Local Variables------------------------------
!
       logical lflg
!
       integer i
!
       real(kp) omeps
       real(kp) trinv
       real(kp) tc
       real(kp) weight
       real(kp) hltalt
!
       real(kp) hlatsb
       real(kp) hlatvp
       real(kp) tterm
       real(kp) desdt
!
!-----------------------------------------------------------------------
!
       omeps = 1.0_kp - epsqs
       do i=1,len
         es(i) = min(estblf(t(i)), p(i))
!
! Saturation specific humidity
!
         qs(i) = epsqs*es(i)/(p(i) - omeps*es(i))
!
! The following check is to avoid the generation of negative
! values that can occur in the upper stratosphere and mesosphere
!
         qs(i) = min(1.0_kp,qs(i))
!
!        if (qs(i) < 0.0_kp) then
!          qs(i) = 1.0_kp
!          es(i) = p(i)
!        end if

       end do
!
! "generalized" analytic expression for t derivative of es
! accurate to within 1 percent for 173.16 < t < 373.16
!
       trinv = 0.0_kp
       if ((.not. icephs) .or. (ttrice.eq.0.0_kp)) go to 10
       trinv = 1.0_kp/ttrice
       do i=1,len
!
! Weighting of hlat accounts for transition from water to ice
! polynomial expression approximates difference between es over
! water and es over ice from 0 to -ttrice (C) (min of ttrice is
! -40): required for accurate estimate of es derivative in transition
! range from ice to water also accounting for change of hlatv with t
! above freezing where const slope is given by -2369 j/(kg c) = cpv - cw
!
       tc = t(i) - tmelt
       lflg = (tc >= -ttrice .and. tc < 0.0_kp)
       weight = min(-tc*trinv,1.0_kp)
       hlatsb = hlatv + weight*hlatf
       hlatvp = hlatv - 2369.0_kp*tc
       if (t(i) < tmelt) then
       hltalt = hlatsb
       else
       hltalt = hlatvp
       end if
       if (lflg) then
       tterm = pcf(1) + tc*(pcf(2) + tc*(pcf(3) + tc*(pcf(4)            &
     &                + tc*pcf(5))))
       else
       tterm = 0.0_kp
       end if
       desdt = hltalt*es(i)/(rgasv*t(i)*t(i)) + tterm*trinv
       gam(i) = hltalt*qs(i)*p(i)*desdt/(cp*es(i)*(p(i) - omeps*es(i)))
       if (qs(i) == 1.0_kp) gam(i) = 0.0_kp
       end do
       return
!
! No icephs or water to ice transition
!
10      do i=1,len
!
! Account for change of hlatv with t above freezing where
! constant slope is given by -2369 j/(kg c) = cpv - cw
!
       hlatvp = hlatv - 2369.0_kp*(t(i)-tmelt)
       if (icephs) then
       hlatsb = hlatv + hlatf
       else
       hlatsb = hlatv
       end if
       if (t(i) < tmelt) then
       hltalt = hlatsb
       else
       hltalt = hlatvp
       end if
       desdt = hltalt*es(i)/(rgasv*t(i)*t(i))
       gam(i) = hltalt*qs(i)*p(i)*desdt/(cp*es(i)*(p(i) - omeps*es(i)))
       if (qs(i) == 1.0_kp) gam(i) = 0.0_kp
       end do
!
       return
!
      end subroutine vqsatd

!++xl
!>\ingroup wv_saturation_mod
!!
      subroutine vqsatd_water(t, p, es, qs, gam, len)

!------------------------------Arguments--------------------------------
!
! Input arguments
!
       integer,  intent(in) :: len
       real(kp), intent(in) :: t(len)
       real(kp), intent(in) :: p(len)

!
! Output arguments
!
       real(kp), intent(out) :: es(len)
       real(kp), intent(out) :: qs(len)
       real(kp), intent(out) :: gam(len)

!
!--------------------------Local Variables------------------------------
!
!
       integer i
!
       real(kp) omeps
       real(kp) hltalt
!
       real(kp) hlatsb
       real(kp) hlatvp
       real(kp) desdt
!
!-----------------------------------------------------------------------
!
       omeps = 1.0_kp - epsqs
       do i=1,len
#ifdef NEMS_GSM
         es(i) = min(fpvsl(t(i)), p(i))
#else
         es(i) = min(polysvp(t(i),0), p(i))
#endif
!
! Saturation specific humidity
!
         qs(i) = min(1.0_kp, epsqs*es(i) / (p(i)-omeps*es(i)))
!
! The following check is to avoid the generation of negative
! values that can occur in the upper stratosphere and mesosphere
!
!        qs(i) = min(1.0_kp,qs(i))
!
!        if (qs(i) < 0.0_kp) then
!          qs(i) = 1.0_kp
!          es(i) = p(i)
!        end if

       end do
!
! No icephs or water to ice transition
!
       do i=1,len
!
! Account for change of hlatv with t above freezing where
! constant slope is given by -2369 j/(kg c) = cpv - cw
!
       hlatvp = hlatv - 2369.0_kp*(t(i)-tmelt)
       hlatsb = hlatv
       if (t(i) < tmelt) then
         hltalt = hlatsb
       else
         hltalt = hlatvp
       end if
       desdt = hltalt*es(i)/(rgasv*t(i)*t(i))
       gam(i) = hltalt*qs(i)*p(i)*desdt/(cp*es(i)*(p(i) - omeps*es(i)))
       if (qs(i) == 1.0_kp) gam(i) = 0.0_kp
       end do
!
       return
!
      end subroutine vqsatd_water

!>\ingroup wv_saturation_mod
      function polysvp (T,typ)
!> This function computes saturation vapor pressure by using
!! function from Goff and Gatch (1946)
!!  Polysvp returned in units of pa.
!!  T is input in units of K.
!!  type refers to saturation with respect to liquid (0) or ice (1)

!!DONIFF Changed to Murphy and Koop (2005) (03/04/14)


       real(kp) dum

       real(kp) t,polysvp

       integer typ


      if (.true.) then
!ice
       if (typ == 1) then
         polysvp = MurphyKoop_svp_ice(t)
       end if
       if (typ == 0) then
         polysvp = MurphyKoop_svp_water(t)
       end if

      else

! ice
       if (typ.eq.1) then



       polysvp = 10._kp**(-9.09718_kp*(273.16_kp/t-1._kp)-3.56654_kp*   &
     & log10(273.16_kp/t)+0.876793_kp*(1._kp-t/273.16_kp)+              &
     & log10(6.1071_kp))*100._kp

       end if



       if (typ.eq.0) then
       polysvp = 10._kp**(-7.90298_kp*(373.16_kp/t-1._kp)+ 5.02808_kp*   &
     &log10(373.16_kp/t)- 1.3816e-7_kp*(10._kp**(11.344_kp*(1._kp-t/     &
     &373.16_kp))-1._kp)+ 8.1328e-3_kp*(10._kp**(-3.49149_kp*(373.16_kp/ &
     &t-1._kp))-1._kp)+ log10(1013.246_kp))*100._kp
       end if

       end if

       end function polysvp
!--xl



      integer function fqsatd(t ,p ,es ,qs ,gam , len )





       integer, intent(in) :: len
       real(kp), intent(in) :: t(len)
       real(kp), intent(in) :: p(len)

       real(kp), intent(out) :: es(len)
       real(kp), intent(out) :: qs(len)
       real(kp), intent(out) :: gam(len)

       call vqsatd(t ,p ,es ,qs ,gam , len )
       fqsatd = 1
       return
      end function fqsatd

      real(kp) function qsat_water(t,p)

       real(kp) t
       real(kp) p
       real(kp) es
       real(kp) ps, ts, e1, e2, f1, f2, f3, f4, f5, f





       ps = 1013.246_kp
       ts = 373.16_kp
       e1 = 11.344_kp*(1.0_kp - t/ts)
       e2 = -3.49149_kp*(ts/t - 1.0_kp)
       f1 = -7.90298_kp*(ts/t - 1.0_kp)
       f2 = 5.02808_kp*log10(ts/t)
       f3 = -1.3816_kp*(10.0_kp**e1 - 1.0_kp)/10000000.0_kp
       f4 = 8.1328_kp*(10.0_kp**e2 - 1.0_kp)/1000.0_kp
       f5 = log10(ps)
       f = f1 + f2 + f3 + f4 + f5
       es = (10.0_kp**f)*100.0_kp

       qsat_water = epsqs*es/(p-(1.-epsqs)*es)
       if(qsat_water < 0.) qsat_water = 1.

      end function qsat_water

!>\ingroup wv_saturation_mod
!! This subroutine 
      subroutine vqsat_water(t,p,qsat_water,len)

       integer, intent(in) :: len
       real(kp) t(len)
       real(kp) p(len)
       real(kp) qsat_water(len)
       real(kp) es
       real(kp), parameter :: t0inv = 1._kp/273._kp
       real(kp) coef
       integer :: i

       coef = hlatv/rgasv
       do i=1,len
       es = 611._kp*exp(coef*(t0inv-1./t(i)))
       qsat_water(i) = epsqs*es/(p(i)-(1.-epsqs)*es)
       if(qsat_water(i) < 0.) qsat_water(i) = 1.
       enddo

       return

      end subroutine vqsat_water

!>\ingroup wv_saturation_mod
      real(kp) function qsat_ice(t,p)

       real(kp) t
       real(kp) p
       real(kp) es
       real(kp), parameter :: t0inv = 1._kp/273._kp
       es = 611.*exp((hlatv+hlatf)/rgasv*(t0inv-1./t))
       qsat_ice = epsqs*es/(p-(1.-epsqs)*es)
       if(qsat_ice < 0.) qsat_ice = 1.

      end function qsat_ice

!>\ingroup wv_saturation_mod
      subroutine vqsat_ice(t,p,qsat_ice,len)

       integer,intent(in) :: len
       real(kp) t(len)
       real(kp) p(len)
       real(kp) qsat_ice(len)
       real(kp) es
       real(kp), parameter :: t0inv = 1._kp/273._kp
       real(kp) coef
       integer :: i

       coef = (hlatv+hlatf)/rgasv
       do i=1,len
       es = 611.*exp(coef*(t0inv-1./t(i)))
       qsat_ice(i) = epsqs*es/(p(i)-(1.-epsqs)*es)
       if(qsat_ice(i) < 0.) qsat_ice(i) = 1.
       enddo

       return

      end subroutine vqsat_ice

! Sungsu
! Below two subroutines (vqsatd2_water,vqsatd2_water_single) are by Sungsu
! Replace 'gam -> dqsdt'
! Sungsu

!>\ingroup wv_saturation_mod
      subroutine vqsatd2_water(t ,p ,es ,qs ,dqsdt , len )

!------------------------------Arguments--------------------------------
!
! Input arguments
!
       integer, intent(in) :: len
       real(kp), intent(in) :: t(len)
       real(kp), intent(in) :: p(len)

!
! Output arguments
!
       real(kp), intent(out) :: es(len)
       real(kp), intent(out) :: qs(len)


       real(kp), intent(out) :: dqsdt(len)


!
!--------------------------Local Variables------------------------------
!
!
       integer i
!
       real(kp) omeps
       real(kp) hltalt
!
       real(kp) hlatsb
       real(kp) hlatvp
       real(kp) desdt


       real(kp) gam(len)


!
!-----------------------------------------------------------------------
!
       omeps = 1.0_kp - epsqs
       do i=1,len
#ifdef GEOS5
       es(i) = min(polysvp(t(i),0), p(i))
#endif
#ifdef NEMS_GSM
       es(i) = min(fpvsl(t(i)), p(i))
#endif
!
! Saturation specific humidity
!
       qs(i) = epsqs*es(i)/(p(i) - omeps*es(i))
!
! The following check is to avoid the generation of negative
! values that can occur in the upper stratosphere and mesosphere
!
       qs(i) = min(1.0_kp,qs(i))
!
!      if (qs(i) < 0.0_kp) then
!        qs(i) = 1.0_kp
!        es(i) = p(i)
!      end if

       end do
!
! No icephs or water to ice transition
!
       do i=1,len
!
! Account for change of hlatv with t above freezing where
! constant slope is given by -2369 j/(kg c) = cpv - cw
!
       hlatvp = hlatv - 2369.0_kp*(t(i)-tmelt)
       hlatsb = hlatv
       if (t(i) < tmelt) then
         hltalt = hlatsb
       else
         hltalt = hlatvp
       end if
       desdt  = hltalt*es(i)/(rgasv*t(i)*t(i))
       gam(i) = hltalt*qs(i)*p(i)*desdt/(cp*es(i)*(p(i) - omeps*es(i)))
       if (qs(i) == 1.0_kp) gam(i) = 0.0_kp

       dqsdt(i) = (cp/hltalt)*gam(i)

       end do
!
       return
!
      end subroutine vqsatd2_water

!>\ingroup wv_saturation_mod
      subroutine vqsatd2_water_single(t ,p ,es ,qs ,dqsdt)

!------------------------------Arguments--------------------------------
!
! Input arguments
!

       real(kp), intent(in) :: t, p

!
! Output arguments
!
       real(kp), intent(out) :: es, qs, dqsdt
!
!--------------------------Local Variables------------------------------
!
!      integer i
!
       real(kp) omeps, hltalt, hlatsb, hlatvp, desdt, gam
!
!-----------------------------------------------------------------------
!
       omeps = 1.0_kp - epsqs
!  do i=1,len
#ifdef GEOS5
       es = min(p, polysvp(t,0))
#endif
#ifdef NEMS_GSM
       es = min(p, fpvsl(t))
#endif
!
! Saturation specific humidity
!
       qs = min(1.0_kp, epsqs*es/(p-omeps*es))
!
! The following check is to avoid the generation of negative
! values that can occur in the upper stratosphere and mesosphere
!
!      if (qs < 0.0_kp) then
!      qs = 1.0_kp
!      es = p
!      end if
!  end do
!
! No icephs or water to ice transition
!
!  do i=1,len
!
! Account for change of hlatv with t above freezing where
! constant slope is given by -2369 j/(kg c) = cpv - cw
!
       hlatvp = hlatv - 2369.0_kp*(t-tmelt)
       hlatsb = hlatv
       if (t < tmelt) then
         hltalt = hlatsb
       else
         hltalt = hlatvp
       end if
       desdt = hltalt*es/(rgasv*t*t)
       gam   = hltalt*qs*p*desdt/(cp*es*(p - omeps*es))
       if (qs >= 1.0_kp) gam = 0.0_kp

       dqsdt = (cp/hltalt)*gam

!  end do
!
       return
!
      end subroutine vqsatd2_water_single

!>\ingroup wv_saturation_mod
      subroutine vqsatd2(t ,p ,es ,qs ,dqsdt , len )
!-----------------------------------------------------------------------
! Sungsu : This is directly copied from 'vqsatd' but 'dqsdt' is output
!          instead of gam for use in Sungsu's equilibrium stratiform
!          macrophysics scheme.
!
! Purpose:
! Utility procedure to look up and return saturation vapor pressure from
! precomputed table, calculate and return saturation specific humidity
! (g/g), and calculate and return gamma (l/cp)*(d(qsat)/dT).  The same
! function as qsatd, but operates on vectors of temperature and pressure
!
! Method:
!
! Author: J. Hack
!
!------------------------------Arguments--------------------------------
!
! Input arguments
!
       integer, intent(in) :: len
       real(kp), intent(in) :: t(len)
       real(kp), intent(in) :: p(len)
!
! Output arguments
!
       real(kp), intent(out) :: es(len)
       real(kp), intent(out) :: qs(len)


       real(kp), intent(out) :: dqsdt(len)


!
!--------------------------Local Variables------------------------------
!
       logical lflg
!
       integer i
!
       real(kp) omeps
       real(kp) trinv
       real(kp) tc
       real(kp) weight
       real(kp) hltalt
!
       real(kp) hlatsb
       real(kp) hlatvp
       real(kp) tterm
       real(kp) desdt


       real(kp) gam(len)

!
!-----------------------------------------------------------------------
!
       omeps = 1.0_kp - epsqs
       do i=1,len
#ifdef GEOS5
         es(i) = min(p(i), estblf(t(i)))
#endif
#ifdef NEMS_GSM
         es(i) = min(p(i), fpvsi(t(i)))
#endif
!
! Saturation specific humidity
!
         qs(i) = epsqs*es(i)/(p(i) - omeps*es(i))
!
! The following check is to avoid the generation of negative
! values that can occur in the upper stratosphere and mesosphere
!
         qs(i) = min(1.0_kp,qs(i))
!
!        if (qs(i) < 0.0_kp) then
!        qs(i) = 1.0_kp
!        es(i) = p(i)
!        end if
       end do
!
! "generalized" analytic expression for t derivative of es
! accurate to within 1 percent for 173.16 < t < 373.16
!
       trinv = 0.0_kp
       if ((.not. icephs) .or. (ttrice == 0.0_kp)) go to 10
       trinv = 1.0_kp/ttrice
       do i=1,len
!
! Weighting of hlat accounts for transition from water to ice
! polynomial expression approximates difference between es over
! water and es over ice from 0 to -ttrice (C) (min of ttrice is
! -40): required for accurate estimate of es derivative in transition
! range from ice to water also accounting for change of hlatv with t
! above freezing where const slope is given by -2369 j/(kg c) = cpv - cw
!
       tc = t(i) - tmelt
       lflg = (tc >= -ttrice .and. tc < 0.0_kp)
       weight = min(-tc*trinv,1.0_kp)
       hlatsb = hlatv + weight*hlatf
       hlatvp = hlatv - 2369.0_kp*tc
       if (t(i) < tmelt) then
         hltalt = hlatsb
       else
         hltalt = hlatvp
       end if
       if (lflg) then
         tterm = pcf(1) + tc*(pcf(2) + tc*(pcf(3) + tc*(pcf(4)          &
     &                  + tc*pcf(5))))
       else
         tterm = 0.0_kp
       end if
       desdt  = hltalt*es(i)/(rgasv*t(i)*t(i)) + tterm*trinv
       gam(i) = hltalt*qs(i)*p(i)*desdt/(cp*es(i)*(p(i) - omeps*es(i)))
       if (qs(i) == 1.0_kp) gam(i) = 0.0_kp

       dqsdt(i) = (cp/hltalt)*gam(i)

       end do
       return
!
! No icephs or water to ice transition
!
10     do i=1,len
!
! Account for change of hlatv with t above freezing where
! constant slope is given by -2369 j/(kg c) = cpv - cw
!
       hlatvp = hlatv - 2369.0_kp*(t(i)-tmelt)
       if (icephs) then
       hlatsb = hlatv + hlatf
       else
       hlatsb = hlatv
       end if
       if (t(i) < tmelt) then
       hltalt = hlatsb
       else
       hltalt = hlatvp
       end if
       desdt = hltalt*es(i)/(rgasv*t(i)*t(i))
       gam(i) = hltalt*qs(i)*p(i)*desdt/(cp*es(i)*(p(i) - omeps*es(i)))
       if (qs(i) == 1.0_kp) gam(i) = 0.0_kp

       dqsdt(i) = (cp/hltalt)*gam(i)

       end do
!
       return
!
      end subroutine vqsatd2


! Below routine is by Sungsu

!>\ingroup wv_saturation_mod
      subroutine vqsatd2_single(t ,p ,es ,qs ,dqsdt)
!-----------------------------------------------------------------------
! Sungsu : This is directly copied from 'vqsatd' but 'dqsdt' is output
!          instead of gam for use in Sungsu's equilibrium stratiform
!          macrophysics scheme.
!
! Purpose:
! Utility procedure to look up and return saturation vapor pressure from
! precomputed table, calculate and return saturation specific humidity
! (g/g), and calculate and return gamma (l/cp)*(d(qsat)/dT).  The same
! function as qsatd, but operates on vectors of temperature and pressure
!
! Method:
!
! Author: J. Hack
!
!------------------------------Arguments--------------------------------
!
! Input arguments
!
       real(kp), intent(in) :: t, p
!
! Output arguments
!
       real(kp), intent(out) :: es, qs, dqsdt
!
!--------------------------Local Variables------------------------------
!
       logical lflg
!
!  integer i      ! index for vector calculations
!
       real(kp) omeps
       real(kp) trinv
       real(kp) tc
       real(kp) weight
       real(kp) hltalt
!
       real(kp) hlatsb
       real(kp) hlatvp
       real(kp) tterm
       real(kp) desdt

       real(kp) gam

!
!-----------------------------------------------------------------------
!
       omeps = 1.0_kp - epsqs

!  do i=1,len

#ifdef GEOS5
       es = estblf(t)
#endif
#ifdef NEMS_GSM
       es = min(fpvs(t), p)
#endif
!
! Saturation specific humidity
!
       qs = epsqs*es/(p - omeps*es)
!
! The following check is to avoid the generation of negative
! values that can occur in the upper stratosphere and mesosphere
!
       qs = min(1.0_kp,qs)
!
!      if (qs < 0.0_kp) then
!        qs = 1.0_kp
!        es = p
!      end if

!  end do
!
! "generalized" analytic expression for t derivative of es
! accurate to within 1 percent for 173.16 < t < 373.16
!
       trinv = 0.0_kp
       if ((.not. icephs) .or. (ttrice == 0.0_kp)) go to 10
       trinv = 1.0_kp/ttrice

!  do i=1,len
!
! Weighting of hlat accounts for transition from water to ice
! polynomial expression approximates difference between es over
! water and es over ice from 0 to -ttrice (C) (min of ttrice is
! -40): required for accurate estimate of es derivative in transition
! range from ice to water also accounting for change of hlatv with t
! above freezing where const slope is given by -2369 j/(kg c) = cpv - cw
!
       tc = t - tmelt
       lflg = (tc >= -ttrice .and. tc < 0.0_kp)
       weight = min(-tc*trinv,1.0_kp)
       hlatsb = hlatv + weight*hlatf
       hlatvp = hlatv - 2369.0_kp*tc
       if (t < tmelt) then
         hltalt = hlatsb
       else
         hltalt = hlatvp
       end if
       if (lflg) then
         tterm = pcf(1) + tc*(pcf(2) + tc*(pcf(3) + tc*(pcf(4)          &
     &                  + tc*pcf(5))))
       else
         tterm = 0.0_kp
       end if
       desdt = hltalt*es/(rgasv*t*t) + tterm*trinv
       gam = hltalt*qs*p*desdt/(cp*es*(p - omeps*es))
       if (qs == 1.0_kp) gam = 0.0_kp

       dqsdt = (cp/hltalt)*gam

!  end do
       return
!
! No icephs or water to ice transition
!

10     continue

!10 do i=1,len
!
! Account for change of hlatv with t above freezing where
! constant slope is given by -2369 j/(kg c) = cpv - cw
!
       hlatvp = hlatv - 2369.0_kp*(t-tmelt)
       if (icephs) then
         hlatsb = hlatv + hlatf
       else
         hlatsb = hlatv
       end if
       if (t < tmelt) then
         hltalt = hlatsb
       else
         hltalt = hlatvp
       end if
       desdt = hltalt*es/(rgasv*t*t)
       gam = hltalt*qs*p*desdt/(cp*es*(p - omeps*es))
       if (qs == 1.0_kp) gam = 0.0_kp

       dqsdt = (cp/hltalt)*gam


!  end do
!
       return
!
      end subroutine vqsatd2_single

!----------------------------------------------------------------------

!----------------------------------------------------------------------

!>\ingroup wv_saturation_mod
      subroutine gffgch(t ,es ,tmelt ,itype )
!-----------------------------------------------------------------------
!
! Purpose:
! Computes saturation vapor pressure over water and/or over ice using
! Goff & Gratch (1946) relationships.
! <Say what the routine does>
!
! Method:
! T (temperature), and itype are input parameters, while es (saturation
! vapor pressure) is an output parameter.  The input parameter itype
! serves two purposes: a value of zero indicates that saturation vapor
! pressures over water are to be returned (regardless of temperature),
! while a value of one indicates that saturation vapor pressures over
! ice should be returned when t is less than freezing degrees.  If itype
! is negative, its absolute value is interpreted to define a temperature
! transition region below freezing in which the returned
! saturation vapor pressure is a weighted average of the respective ice
! and water value.  That is, in the temperature range 0 => -itype
! degrees c, the saturation vapor pressures are assumed to be a weighted
! average of the vapor pressure over supercooled water and ice (all
! water at 0 c; all ice at -itype c).  Maximum transition range => 40 c
!
! Author: J. Hack
!
!-----------------------------------------------------------------------





       implicit none
!------------------------------Arguments--------------------------------
!
! Input arguments
!
       real(kp), intent(in) :: t ,tmelt
!
! Output arguments
!
       integer, intent(inout) :: itype

       real(kp), intent(out) :: es
!
!---------------------------Local variables-----------------------------
!
       real(kp) e1
       real(kp) e2
       real(kp) eswtr
       real(kp) f
       real(kp) f1
       real(kp) f2
       real(kp) f3
       real(kp) f4
       real(kp) f5
       real(kp) ps
       real(kp) t0
       real(kp) term1
       real(kp) term2
       real(kp) term3
       real(kp) tr
       real(kp) ts
       real(kp) weight
       integer itypo
!
!-----------------------------------------------------------------------
!
! Check on whether there is to be a transition region for es
!
       if (itype < 0) then
       tr = abs(real(itype,kp))
       itypo = itype
       itype = 1
       else
       tr = 0.0_kp
       itypo = itype
       end if
       if (tr > 40.0_kp) then
       write(iulog,900) tr

       end if
!
       if(t < (tmelt - tr) .and. itype == 1) go to 10
!
! Water
!
       ps = 1013.246_kp
       ts = 373.16_kp
       e1 = 11.344_kp*(1.0_kp - t/ts)
       e2 = -3.49149_kp*(ts/t - 1.0_kp)
       f1 = -7.90298_kp*(ts/t - 1.0_kp)
       f2 = 5.02808_kp*log10(ts/t)
       f3 = -1.3816_kp*(10.0_kp**e1 - 1.0_kp)/10000000.0_kp
       f4 = 8.1328_kp*(10.0_kp**e2 - 1.0_kp)/1000.0_kp
       f5 = log10(ps)
       f = f1 + f2 + f3 + f4 + f5
       es = (10.0_kp**f)*100.0_kp
       eswtr = es
!
       if(t >= tmelt .or. itype == 0) go to 20
!
! Ice
!
10     continue
       t0 = tmelt
       term1 = 2.01889049_kp/(t0/t)
       term2 = 3.56654_kp*log(t0/t)
       term3 = 20.947031_kp*(t0/t)
       es = 575.185606e10_kp*exp(-(term1 + term2 + term3))
!
       if (t < (tmelt - tr)) go to 20
!
! Weighted transition between water and ice
!
       weight = min((tmelt - t)/tr,1.0_kp)
       es = weight*es + (1.0_kp - weight)*eswtr
!
20     continue
       itype = itypo
       return
!
900    format('GFFGCH: FATAL ERROR ******************************',/,    &
     & 'TRANSITION RANGE FOR WATER TO ICE SATURATION VAPOR', ' PRESSURE, &
     & TR, EXCEEDS MAXIMUM ALLOWABLE VALUE OF', ' 40.0 DEGREES C',/,     &
     & ' TR = ',f7.2)
!
      end subroutine gffgch

!>\ingroup wv_saturation_mod
!!DONIF USe Murphy and Koop (2005) (Written by Andrew Gettelman)
       function MurphyKoop_svp_water(tx) result(es)
       real(kp), intent(in) :: tx
       real(kp) :: es
       real(kp):: t

       t=min(tx, 332.0_kp)
       t=max(123.0_kp, tx)

       es = exp(54.842763_kp - (6763.22_kp / t) - (4.210_kp * log(t)) + &
     & (0.000367_kp * t) + (tanh(0.0415_kp * (t - 218.8_kp)) *          &
     & (53.878_kp - (1331.22_kp / t) - (9.44523_kp * log(t)) +          &
     & 0.014025_kp * t)))

      end function MurphyKoop_svp_water

       function MurphyKoop_svp_ice(tx) result(es)
       real(kp), intent(in) :: tx
       real(kp) :: t
       real(kp) :: es

       t=max(100.0_kp, tx)
       t=min(274.0_kp, tx)


       es = exp(9.550426_kp - (5723.265_kp / t) + (3.53068_kp *         &
     & log(t)) - (0.00728332_kp * t))

      end function MurphyKoop_svp_ice
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!>\ingroup wv_saturation_mod
      subroutine vqsatd2_ice_single(t ,p ,es ,qs ,dqsdt)

!------------------------------Arguments--------------------------------
!
! Input arguments
!
       real(kp), intent(in) :: t, p
!
! Output arguments
!
       real(kp), intent(out) :: es, qs, dqsdt

!
!--------------------------Local Variables------------------------------
!
!      integer i
!
       real(kp) omeps, hltalt, hlatsb, hlatvp, desdt, gam
!
!-----------------------------------------------------------------------
!
       omeps = 1.0_kp - epsqs
!  do i=1,len
#ifdef GEOS5
       es = min(polysvp(t,1),p)
#endif
#ifdef NEMS_GSM
       es = min(fpvsi(t),p)
#endif
!
! Saturation specific humidity
!
       qs = min(1.0_kp, epsqs*es/(p-omeps*es))
!
! The following check is to avoid the generation of negative
! values that can occur in the upper stratosphere and mesosphere
!
!      if (qs < 0.0_kp) then
!        qs = 1.0_kp
!        es = p
!      end if
!  end do
!
! No icephs or water to ice transition
!
!  do i=1,len
!
! Account for change of hlatv with t above freezing where
! constant slope is given by -2369 j/(kg c) = cpv - cw
!
       hltalt = hlatv + hlatf
       desdt  = hltalt*es/(rgasv*t*t)
       if (qs < 1.0_kp) then
         gam = hltalt*qs*p*desdt/(cp*es*(p - omeps*es))
       else
         gam = 0.0_kp
       endif

       dqsdt = (cp/hltalt)*gam

!  end do
!
       return
!
      end subroutine vqsatd2_ice_single

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      end module wv_saturation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

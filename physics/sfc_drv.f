!>  \file sfc_drv.f
!!  This file contains the Noah land surface scheme driver.

!> This module contains the CCPP-compliant Noah land surface scheme driver.
      module lsm_noah

      use set_soilveg_mod,  only: set_soilveg

      implicit none

      private

      public :: lsm_noah_init, lsm_noah_run, lsm_noah_finalize

      contains

!>\ingroup Noah_LSM
!! This subroutine contains the CCPP-compliant lsm_noah_init to initialize soil vegetation.
!! \section arg_table_lsm_noah_init Argument Table
!! \htmlinclude lsm_noah_init.html
!!
      subroutine lsm_noah_init(me, isot, ivegsrc, nlunit,
     &                         errmsg, errflg)

      implicit none

      integer,              intent(in)  :: me, isot, ivegsrc, nlunit
      character(len=*),     intent(out) :: errmsg
      integer,              intent(out) :: errflg

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      !--- initialize soil vegetation
      call set_soilveg(me, isot, ivegsrc, nlunit)

      end subroutine lsm_noah_init


!! \section arg_table_lsm_noah_finalize Argument Table
!! \htmlinclude lsm_noah_finalize.html
!!
      subroutine lsm_noah_finalize(errmsg, errflg)

      implicit none

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      end subroutine lsm_noah_finalize


! ===================================================================== !
!  description:                                                         !
!                                                                       !
!  usage:                                                               !
!                                                                       !
!      call sfc_drv                                                     !
!  ---  inputs:                                                         !
!          ( im, km, ps, t1, q1, soiltyp, vegtype, sigmaf,              !
!            sfcemis, dlwflx, dswsfc, snet, delt, tg3, cm, ch,          !
!            prsl1, prslki, zf, land, wind,  slopetyp,                  !
!            shdmin, shdmax, snoalb, sfalb, flag_iter, flag_guess,      !
!            lheatstrg, isot, ivegsrc,                                  !
!  ---  in/outs:                                                        !
!            weasd, snwdph, tskin, tprcp, srflag, smc, stc, slc,        !
!            canopy, trans, tsurf, zorl,                                !
!  ---  outputs:                                                        !
!            sncovr1, qsurf, gflux, drain, evap, hflx, ep, runoff,      !
!            cmm, chh, evbs, evcw, sbsno, snowc, stm, snohf,            !
!            smcwlt2, smcref2, wet1 )                                   !
!                                                                       !
!                                                                       !
!  subprogram called:  sflx                                             !
!                                                                       !
!  program history log:                                                 !
!         xxxx  --             created                                  !
!         200x  -- sarah lu    modified                                 !
!    oct  2006  -- h. wei      modified                                 !
!    apr  2009  -- y.-t. hou   modified to include surface emissivity   !
!                     effect on lw radiation. replaced the comfussing   !
!                     slrad (net sw + dlw) with sfc net sw snet=dsw-usw !
!    sep  2009  -- s. moorthi modification to remove rcl and unit change!
!    nov  2011  -- sarah lu    corrected wet1 calculation
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                       size   !
!     im       - integer, horiz dimention and num of used pts      1    !
!     km       - integer, vertical soil layer dimension            1    !
!     ps       - real, surface pressure (pa)                       im   !
!     t1       - real, surface layer mean temperature (k)          im   !
!     q1       - real, surface layer mean specific humidity        im   !
!     soiltyp  - integer, soil type (integer index)                im   !
!     vegtype  - integer, vegetation type (integer index)          im   !
!     sigmaf   - real, areal fractional cover of green vegetation  im   !
!     sfcemis  - real, sfc lw emissivity ( fraction )              im   !
!     dlwflx   - real, total sky sfc downward lw flux ( w/m**2 )   im   !
!     dswflx   - real, total sky sfc downward sw flux ( w/m**2 )   im   !
!     snet     - real, total sky sfc netsw flx into ground(w/m**2) im   !
!     delt     - real, time interval (second)                      1    !
!     tg3      - real, deep soil temperature (k)                   im   !
!     cm       - real, surface exchange coeff for momentum (m/s)   im   !
!     ch       - real, surface exchange coeff heat & moisture(m/s) im   !
!     prsl1    - real, sfc layer 1 mean pressure (pa)              im   !
!     prslki   - real, dimensionless exner function at layer 1     im   !
!     zf       - real, height of bottom layer (m)                  im   !
!     land     - logical, = T if a point with any land             im   !
!     wind     - real, wind speed (m/s)                            im   !
!     slopetyp - integer, class of sfc slope (integer index)       im   !
!     shdmin   - real, min fractional coverage of green veg        im   !
!     shdmax   - real, max fractnl cover of green veg (not used)   im   !
!     snoalb   - real, upper bound on max albedo over deep snow    im   !
!     sfalb    - real, mean sfc diffused sw albedo (fractional)    im   !
!     flag_iter- logical,                                          im   !
!     flag_guess-logical,                                          im   !
!     lheatstrg- logical, flag for canopy heat storage             1    !
!                         parameterization                              !
!     isot     - integer, sfc soil type data source zobler or statsgo   !
!     ivegsrc  - integer, sfc veg type data source umd or igbp          !
!                                                                       !
!  input/outputs:                                                       !
!     weasd    - real, water equivalent accumulated snow depth (mm) im  !
!     snwdph   - real, snow depth (water equiv) over land          im   !
!     tskin    - real, ground surface skin temperature ( k )       im   !
!     tprcp    - real, total precipitation                         im   !
!     srflag   - real, snow/rain flag for precipitation            im   !
!     smc      - real, total soil moisture content (fractional)   im,km !
!     stc      - real, soil temp (k)                              im,km !
!     slc      - real, liquid soil moisture                       im,km !
!     canopy   - real, canopy moisture content (m)                 im   !
!     trans    - real, total plant transpiration (m/s)             im   !
!     tsurf    - real, surface skin temperature (after iteration)  im   !
!     zorl     - real, surface roughness                           im   !
!     sncovr1  - real, snow cover over land (fractional)           im   !
!     qsurf    - real, specific humidity at sfc                    im   !
!     gflux    - real, soil heat flux (w/m**2)                     im   !
!     drain    - real, subsurface runoff (mm/s)                    im   !
!     evap     - real, evaperation from latent heat flux           im   !
!     hflx     - real, sensible heat flux                          im   !
!     ep       - real, potential evaporation                       im   !
!     runoff   - real, surface runoff (m/s)                        im   !
!     cmm      - real,                                             im   !
!     chh      - real,                                             im   !
!     evbs     - real, direct soil evaporation (m/s)               im   !
!     evcw     - real, canopy water evaporation (m/s)              im   !
!     sbsno    - real, sublimation/deposit from snopack (m/s)      im   !
!     snowc    - real, fractional snow cover                       im   !
!     stm      - real, total soil column moisture content (m)      im   !
!     snohf    - real, snow/freezing-rain latent heat flux (w/m**2)im   !
!     smcwlt2  - real, dry soil moisture threshold                 im   !
!     smcref2  - real, soil moisture threshold                     im   !
!     wet1     - real, normalized soil wetness                     im   !
!                                                                       !
!  ====================    end of description    =====================  !

!>\defgroup Noah_LSM GFS Noah LSM Model
!! \brief This is Noah LSM driver module, with the functionality of 
!! preparing variables to run Noah LSM gfssflx(), calling Noah LSM and post-processing
!! variables for return to the parent model suite including unit conversion, as well
!! as diagnotics calculation.
!! \section arg_table_lsm_noah_run Argument Table
!! \htmlinclude lsm_noah_run.html
!!
!> \section general_noah_drv GFS sfc_drv General Algorithm
!>  @{
      subroutine lsm_noah_run                                           &
     &     ( im, km, grav, cp, hvap, rd, eps, epsm1, rvrdm1, ps,        & !  ---  inputs:
     &       t1, q1, soiltyp, vegtype, sigmaf,                          &
     &       sfcemis, dlwflx, dswsfc, snet, delt, tg3, cm, ch,          &
     &       prsl1, prslki, zf, land, wind, slopetyp,                   &
     &       shdmin, shdmax, snoalb, sfalb, flag_iter, flag_guess,      &
     &       lheatstrg, isot, ivegsrc,                                  &
     &       bexppert, xlaipert, vegfpert,pertvegf,                     &  ! sfc perts, mgehne
!  ---  in/outs:
     &       weasd, snwdph, tskin, tprcp, srflag, smc, stc, slc,        &
     &       canopy, trans, tsurf, zorl,                                &
!  ---  outputs:
     &       sncovr1, qsurf, gflux, drain, evap, hflx, ep, runoff,      &
     &       cmm, chh, evbs, evcw, sbsno, snowc, stm, snohf,            &
     &       smcwlt2, smcref2, wet1, errmsg, errflg                     &
     &     )
!
      use machine , only : kind_phys
      use funcphys, only : fpvs

      use surface_perturbation, only : ppfbet

      implicit none

!  ---  constant parameters:
      real(kind=kind_phys), parameter :: zero    = 0.0_kind_phys
      real(kind=kind_phys), parameter :: one     = 1.0_kind_phys
      real(kind=kind_phys), parameter :: rhoh2o  = 1000.0_kind_phys
      real(kind=kind_phys), parameter :: a2      = 17.2693882_kind_phys
      real(kind=kind_phys), parameter :: a3      = 273.16_kind_phys
      real(kind=kind_phys), parameter :: a4      = 35.86_kind_phys
      real(kind=kind_phys), parameter :: a23m4   = a2*(a3-a4)
      real(kind=kind_phys), parameter :: qmin    = 1.0e-8_kind_phys

      real(kind=kind_phys), save         :: zsoil_noah(4)
      data zsoil_noah / -0.1_kind_phys, -0.4_kind_phys,                 &
     &                  -1.0_kind_phys, -2.0_kind_phys /

!  ---  input:
      integer, intent(in) :: im, km, isot, ivegsrc
      real (kind=kind_phys), intent(in) :: grav, cp, hvap, rd, eps,     &
     &       epsm1, rvrdm1
      real (kind=kind_phys), intent(in) :: pertvegf

      integer, dimension(im), intent(in) :: soiltyp, vegtype, slopetyp

      real (kind=kind_phys), dimension(im), intent(in) :: ps,           &
     &       t1, q1, sigmaf, sfcemis, dlwflx, dswsfc, snet, tg3, cm,    &
     &       ch, prsl1, prslki, wind, shdmin, shdmax,                   &
     &       snoalb, sfalb, zf,                                         &
     &       bexppert, xlaipert, vegfpert

      real (kind=kind_phys),  intent(in) :: delt

      logical, dimension(im), intent(in) :: flag_iter, flag_guess, land

      logical, intent(in) :: lheatstrg

!  ---  in/out:
      real (kind=kind_phys), dimension(im), intent(inout) :: weasd,     &
     &       snwdph, tskin, tprcp, srflag, canopy, trans, tsurf, zorl

      real (kind=kind_phys), dimension(im,km), intent(inout) ::         &
     &       smc, stc, slc

!  ---  output:
      real (kind=kind_phys), dimension(im), intent(inout) :: sncovr1,   &
     &       qsurf, gflux, drain, evap, hflx, ep, runoff, cmm, chh,     &
     &       evbs, evcw, sbsno, snowc, stm, snohf, smcwlt2, smcref2,    &
     &       wet1

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

!  ---  locals:
      real (kind=kind_phys), dimension(im) :: rch, rho,                 &
     &       q0, qs1, theta1,       weasd_old, snwdph_old,              &
     &       tprcp_old, srflag_old, tskin_old, canopy_old

      real (kind=kind_phys), dimension(km) :: et, sldpth, stsoil,       &
     &       smsoil, slsoil

      real (kind=kind_phys), dimension(im,km) :: zsoil, smc_old,        &
     &       stc_old, slc_old

      real (kind=kind_phys) :: alb, albedo, beta, chx, cmx, cmc,        &
     &       dew, drip, dqsdt2, ec, edir, ett, eta, esnow, etp,         &
     &       flx1, flx2, flx3, ffrozp, lwdn, pc, prcp, ptu, q2,         &
     &       q2sat, solnet, rc, rcs, rct, rcq, rcsoil, rsmin,           &
     &       runoff1, runoff2, runoff3, sfcspd, sfcprs, sfctmp,         &
     &       sfcems, sheat, shdfac, shdmin1d, shdmax1d, smcwlt,         &
     &       smcdry, smcref, smcmax, sneqv, snoalb1d, snowh,            &
     &       snomlt, sncovr, soilw, soilm, ssoil, tsea, th2, tbot,      &
     &       xlai, zlvl, swdn, tem, z0, bexpp, xlaip, vegfp,            &
     &       mv, sv, alphav, betav, vegftmp, cpinv, hvapi, elocp

      integer :: couple, ice, nsoil, nroot, slope, stype, vtype
      integer :: i, k, iflag
!
!===> ...  begin here
!
      cpinv = one/cp
      hvapi = one/hvap
      elocp = hvap/cp

!> - Initialize CCPP error handling variables

      errmsg = ''
      errflg = 0

!> - Save land-related prognostic fields for guess run.

      do i = 1, im
        if (land(i) .and. flag_guess(i)) then
          weasd_old(i)  = weasd(i)
          snwdph_old(i) = snwdph(i)
          tskin_old(i)  = tskin(i)
          canopy_old(i) = canopy(i)
          tprcp_old(i)  = tprcp(i)
          srflag_old(i) = srflag(i)

          do k = 1, km
            smc_old(i,k) = smc(i,k)
            stc_old(i,k) = stc(i,k)
            slc_old(i,k) = slc(i,k)
          enddo
        endif   ! land & flag_guess
      enddo

!  --- ...  initialization block

      do i = 1, im
        if (flag_iter(i) .and. land(i)) then
          ep(i)     = zero
          evap (i)  = zero
          hflx (i)  = zero
          gflux(i)  = zero
          drain(i)  = zero
          canopy(i) = max(canopy(i), zero)

          evbs (i)  = zero
          evcw (i)  = zero
          trans(i)  = zero
          sbsno(i)  = zero
          snowc(i)  = zero
          snohf(i)  = zero

!> - initialize variables wind, q, and rh at level 1.

          q0(i)   = max(q1(i), qmin)   !* q1=specific humidity at level 1 (kg/kg)
          theta1(i) = t1(i) * prslki(i) !* adiabatic temp at level 1 (k)

          rho(i) = prsl1(i) / (rd*t1(i)*(one+rvrdm1*q0(i)))
          qs1(i) = fpvs( t1(i) )        !* qs1=sat. humidity at level 1 (kg/kg)
          qs1(i) = max(eps*qs1(i) / (prsl1(i)+epsm1*qs1(i)), qmin)
          q0 (i) = min(qs1(i), q0(i))

          do k = 1, km
            zsoil(i,k) = zsoil_noah(k)
          enddo

!> - Prepare variables to run Noah LSM:
!!  -   1. configuration information (c):
! couple   couple-uncouple flag (=1: coupled, =0: uncoupled)
! ffrozp   flag for snow-rain detection (1.=snow, 0.=rain)
! ice      sea-ice flag (=1: sea-ice, =0: land)
! dt       timestep (sec) (dt should not exceed 3600 secs) = delt
! zlvl     height (\f$m\f$) above ground of atmospheric forcing variables
! nsoil    number of soil layers (at least 2)
! sldpth   the thickness of each soil layer (\f$m\f$)

          couple = 1                      ! run noah lsm in 'couple' mode
! use srflag directly to allow fractional rain/snow
!          if     (srflag(i) == 1.0) then  ! snow phase
!            ffrozp = 1.0
!          elseif (srflag(i) == 0.0) then  ! rain phase
!            ffrozp = 0.0
!          endif
          ffrozp = srflag(i)
          ice = 0

          zlvl = zf(i)

          nsoil = km
          sldpth(1) = - zsoil(i,1)
          do k = 2, km
            sldpth(k) = zsoil(i,k-1) - zsoil(i,k)
          enddo

!>  -   2. forcing data (f):
! lwdn     lw dw radiation flux (\f$W m^{-2}\f$)
! solnet  - net sw radiation flux (dn-up) (\f$W m^{-2}\f$)
! sfcprs  - pressure at height zlvl above ground (pascals)
! prcp    - precip rate (\f$kg m^{-2} s^{-1}\f$)
! sfctmp  - air temperature (\f$K\f$) at height zlvl above ground
! th2     - air potential temperature (\f$K\f$) at height zlvl above ground
! q2      - mixing ratio at height zlvl above ground (\f$kg kg^{-1}\f$)

          lwdn   = dlwflx(i)         !..downward lw flux at sfc in w/m2
          swdn   = dswsfc(i)         !..downward sw flux at sfc in w/m2
          solnet = snet(i)           !..net sw rad flx (dn-up) at sfc in w/m2
          sfcems = sfcemis(i)

          sfcprs = prsl1(i)
          prcp   = rhoh2o * tprcp(i) / delt
          sfctmp = t1(i)
          th2    = theta1(i)
          q2     = q0(i)

!>  -   3. other forcing (input) data (i):
! sfcspd  - wind speed (\f$m s^{-1}\f$) at height zlvl above ground
! q2sat   - sat mixing ratio at height zlvl above ground (\f$kg kg^{-1}\f$)
! dqsdt2  - slope of sat specific humidity curve at t=sfctmp (\f$kg kg^{-1} k^{-1}\f$)

          sfcspd = wind(i)
          q2sat =  qs1(i)
          dqsdt2 = q2sat * a23m4/(sfctmp-a4)**2

!>  -   4. canopy/soil characteristics (s):
! vegtyp  - vegetation type (integer index)                   -> vtype
! soiltyp - soil type (integer index)                         -> stype
! slopetyp- class of sfc slope (integer index)                -> slope
! shdfac  - areal fractional coverage of green vegetation (0.0-1.0)
! shdmin  - minimum areal fractional coverage of green vegetation -> shdmin1d
! ptu     - photo thermal unit (plant phenology for annuals/crops)
! alb     - backround snow-free surface albedo (fraction)
! snoalb  - upper bound on maximum albedo over deep snow          -> snoalb1d
! tbot    - bottom soil temperature (local yearly-mean sfc air temp)

          vtype = vegtype(i)
          stype = soiltyp(i)
          slope = slopetyp(i)
          shdfac= sigmaf(i)

!>  - Call surface_perturbation::ppfbet() to perturb vegetation fraction that goes into gsflx().
!  perturb vegetation fraction that goes into sflx, use the same
!  perturbation strategy as for albedo (percentile matching)
!! Following Gehne et al. (2018) \cite gehne_et_al_2018, a perturbation of vegetation
!! fraction is added to account for the uncertainty. A percentile matching technique
!! is applied to guarantee the perturbed vegetation fraction is bounded between 0 and
!! 1. The standard deviation of the perturbations is 0.25 for vegetation fraction of
!! 0.5 and the perturbations go to zero as vegetation fraction  approaches its upper
!! or lower bound.
        vegfp  = vegfpert(i)                    ! sfc-perts, mgehne
        if (pertvegf>zero) then
                ! compute beta distribution parameters for vegetation fraction
                mv = shdfac
                sv = pertvegf*mv*(one-mv)
                alphav = mv*mv*(one-mv)/(sv*sv)-mv
                betav  = alphav*(one-mv)/mv
                ! compute beta distribution value corresponding
                ! to the given percentile albPpert to use as new albedo
                call ppfbet(vegfp,alphav,betav,iflag,vegftmp)
                shdfac = vegftmp
        endif
! *** sfc-perts, mgehne

          shdmin1d = shdmin(i)
          shdmax1d = shdmax(i)
          snoalb1d = snoalb(i)

          ptu  = zero
          alb  = sfalb(i)
          tbot = tg3(i)

!>  -   5. history (state) variables (h):
! cmc        - canopy moisture content (\f$m\f$)
! t1         - ground/canopy/snowpack effective skin temperature (\f$K\f$)   -> tsea
! stc(nsoil) - soil temp (\f$K\f$)                                         -> stsoil
! smc(nsoil) - total soil moisture content (volumetric fraction)     -> smsoil
! sh2o(nsoil)- unfrozen soil moisture content (volumetric fraction)  -> slsoil
! snowh      - actual snow depth (\f$m\f$)
! sneqv      - liquid water-equivalent snow depth (\f$m\f$)
! albedo     - surface albedo including snow effect (unitless fraction)
! ch         - surface exchange coefficient for heat and moisture (\f$m s^{-1}\f$) -> chx
! cm         - surface exchange coefficient for momentum (\f$m s^{-1}\f$)          -> cmx
! z0         - surface roughness (\f$m\f$)     -> zorl(\f$cm\f$)

          cmc = canopy(i) * 0.001_kind_phys      ! convert from mm to m
          tsea = tsurf(i)                        ! clu_q2m_iter

          do k = 1, km
            stsoil(k) = stc(i,k)
            smsoil(k) = smc(i,k)
            slsoil(k) = slc(i,k)
          enddo

          snowh = snwdph(i) * 0.001_kind_phys    ! convert from mm to m
          sneqv = weasd(i)  * 0.001_kind_phys    ! convert from mm to m
          if (sneqv /= zero .and. snowh == zero) then
            snowh = 10.0_kind_phys * sneqv
          endif

          chx    = ch(i)  * wind(i)              ! compute conductance
          cmx    = cm(i)  * wind(i)
          chh(i) = chx * rho(i)
          cmm(i) = cmx

!  ---- ... outside sflx, roughness uses cm as unit
          z0 = zorl(i) * 0.01_kind_phys
!  ---- mgehne, sfc-perts
!  - Apply perturbation of soil type b parameter and leave area index.
          bexpp  = bexppert(i)                   ! sfc perts, mgehne
          xlaip  = xlaipert(i)                   ! sfc perts, mgehne

!> - Call Noah LSM gfssflx().

          call gfssflx                                                  & ! ccppdox: these is sflx in mpbl
!  ---  inputs:
     &     ( nsoil, couple, ice, ffrozp, delt, zlvl, sldpth,            &
     &       swdn, solnet, lwdn, sfcems, sfcprs, sfctmp,                &
     &       sfcspd, prcp, q2, q2sat, dqsdt2, th2, ivegsrc,             &
     &       vtype, stype, slope, shdmin1d, alb, snoalb1d,              &
     &       bexpp, xlaip,                                              & ! sfc-perts, mgehne
     &       lheatstrg,                                                 &
!  ---  input/outputs:
     &       tbot, cmc, tsea, stsoil, smsoil, slsoil, sneqv, chx, cmx,  &
     &       z0,                                                        &
!  ---  outputs:
     &       nroot, shdfac, snowh, albedo, eta, sheat, ec,              &
     &       edir, et, ett, esnow, drip, dew, beta, etp, ssoil,         &
     &       flx1, flx2, flx3, runoff1, runoff2, runoff3,               &
     &       snomlt, sncovr, rc, pc, rsmin, xlai, rcs, rct, rcq,        &
     &       rcsoil, soilw, soilm, smcwlt, smcdry, smcref, smcmax)

!> - Noah LSM: prepare variables for return to parent model and unit conversion.
!>  -   6. output (o):
!!\n  eta     - actual latent heat flux (\f$W m^{-2}\f$: positive, if upward from sfc)
!!\n  sheat   - sensible heat flux (\f$W m^{-2}\f$: positive, if upward from sfc)
!!\n  beta    - ratio of actual/potential evap (dimensionless)
!!\n  etp     - potential evaporation (\f$W m^{-2}\f$)
!!\n  ssoil   - soil heat flux (\f$W m^{-2}\f$: negative if downward from surface)
!!\n  runoff1 - surface runoff (\f$m s^{-1}\f$), not infiltrating the surface
!!\n  runoff2 - subsurface runoff (\f$m s^{-1}\f$), drainage out bottom

          evap(i)  = eta
          hflx(i)  = sheat
          gflux(i) = ssoil

          evbs(i)  = edir
          evcw(i)  = ec
          trans(i) = ett
          sbsno(i) = esnow
          snowc(i) = sncovr
          stm(i)   = soilm * 1000.0_kind_phys ! unit conversion (from m to kg m-2)
          snohf(i) = flx1 + flx2 + flx3

          smcwlt2(i) = smcwlt
          smcref2(i) = smcref

          ep(i)      = etp
          tsurf(i)   = tsea

          do k = 1, km
            stc(i,k) = stsoil(k)
            smc(i,k) = smsoil(k)
            slc(i,k) = slsoil(k)
          enddo
          wet1(i) = smsoil(1) / smcmax !Sarah Lu added 09/09/2010 (for GOCART)

!  --- ...  unit conversion (from m s-1 to mm s-1 and kg m-2 s-1)
          runoff(i)  = runoff1 * 1000.0_kind_phys
          drain (i)  = runoff2 * 1000.0_kind_phys

!  --- ...  unit conversion (from m to mm)
          canopy(i)  = cmc   * 1000.0_kind_phys
          snwdph(i)  = snowh * 1000.0_kind_phys
          weasd(i)   = sneqv * 1000.0_kind_phys
          sncovr1(i) = sncovr
!  ---- ... outside sflx, roughness uses cm as unit (update after snow's
!  effect)
          zorl(i) = z0*100.0_kind_phys

!>  - Do not return the following output fields to parent model:
!!\n  ec      - canopy water evaporation (m s-1)
!!\n  edir    - direct soil evaporation (m s-1)
!!\n  et(nsoil)-plant transpiration from a particular root layer (m s-1)
!!\n  ett     - total plant transpiration (m s-1)
!!\n  esnow   - sublimation from (or deposition to if <0) snowpack (m s-1)
!!\n  drip    - through-fall of precip and/or dew in excess of canopy
!!              water-holding capacity (m)
!!\n  dew     - dewfall (or frostfall for t<273.15) (m)
!!\n  beta    - ratio of actual/potential evap (dimensionless)
!!\n  flx1    - precip-snow sfc (w m-2)
!!\n  flx2    - freezing rain latent heat flux (w m-2)
!!\n  flx3    - phase-change heat flux from snowmelt (w m-2)
!!\n  snomlt  - snow melt (m) (water equivalent)
!!\n  sncovr  - fractional snow cover (unitless fraction, 0-1)
!!\n  runoff3 - numerical trunctation in excess of porosity (smcmax)
!!              for a given soil layer at the end of a time step
!!\n  rc      - canopy resistance (s m-1)
!!\n  pc      - plant coefficient (unitless fraction, 0-1) where pc*etp
!!              = actual transp
!!\n  xlai    - leaf area index (dimensionless)
!!\n  rsmin   - minimum canopy resistance (s m-1)
!!\n  rcs     - incoming solar rc factor (dimensionless)
!!\n  rct     - air temperature rc factor (dimensionless)
!!\n  rcq     - atmos vapor pressure deficit rc factor (dimensionless)
!!\n  rcsoil  - soil moisture rc factor (dimensionless)
!!\n  soilw   - available soil moisture in root zone (unitless fraction
!!              between smcwlt and smcmax)
!!\n  soilm   - total soil column moisture content (frozen+unfrozen) (m)
!!\n  smcwlt  - wilting point (volumetric)
!!\n  smcdry  - dry soil moisture threshold where direct evap frm top
!!              layer ends (volumetric)
!!\n  smcref  - soil moisture threshold where transpiration begins to
!!              stress (volumetric)
!!\n  smcmax  - porosity, i.e. saturated value of soil moisture
!!              (volumetric)
!!\n  nroot   - number of root layers, a function of veg type, determined
!!              in subroutine redprm.

!       endif   ! end if flag_iter and flag
!     enddo   ! end do_i_loop

!> - Compute specific humidity at surface (\a qsurf).

          rch(i)   = rho(i) * cp * ch(i) * wind(i)
          qsurf(i) = q1(i)  + evap(i) / (elocp * rch(i))

!> - Compute surface upward sensible heat flux (\a hflx) and evaporation
!! flux (\a evap).
          tem     = one / rho(i)
          hflx(i) = hflx(i) * tem * cpinv
          evap(i) = evap(i) * tem * hvapi

        endif   ! flag_iter & land
      enddo

!> - Restore land-related prognostic fields for guess run.

      do i = 1, im
        if (land(i)) then
          if (flag_guess(i)) then
            weasd(i)  = weasd_old(i)
            snwdph(i) = snwdph_old(i)
            tskin(i)  = tskin_old(i)
            canopy(i) = canopy_old(i)
            tprcp(i)  = tprcp_old(i)
            srflag(i) = srflag_old(i)

            do k = 1, km
              smc(i,k) = smc_old(i,k)
              stc(i,k) = stc_old(i,k)
              slc(i,k) = slc_old(i,k)
            enddo
          else    ! flag_guess = F
            tskin(i) = tsurf(i)
          endif   ! flag_guess
        endif     ! land
      enddo
!
      return
!...................................
      end subroutine lsm_noah_run
!-----------------------------
!> @}

      end module lsm_noah

!>  \file sfc_noahmp_drv.f
!!  This file contains the NoahMP land surface scheme driver.

!>\defgroup NoahMP_LSM NoahMP LSM Model
!! \brief This is the NoahMP LSM driver module, with the functionality of 
!! preparing variables to run the NoahMP LSM subroutine noahmp_sflx(), calling NoahMP LSM and post-processing
!! variables for return to the parent model suite including unit conversion, as well
!! as diagnotics calculation.

!> This module contains the CCPP-compliant NoahMP land surface model driver.
      module noahmpdrv

      implicit none

      private

      public :: noahmpdrv_init, noahmpdrv_run, noahmpdrv_finalize

      contains

!> \ingroup NoahMP_LSM
!! \brief This subroutine is called during the CCPP initialization phase and calls set_soilveg() to 
!! initialize soil and vegetation parameters for the chosen soil and vegetation data sources.
!! \section arg_table_noahmpdrv_init Argument Table
!! \htmlinclude noahmpdrv_init.html
!!
      subroutine noahmpdrv_init(me, isot, ivegsrc, nlunit, pores, resid,
     &                          errmsg, errflg)
        
        use machine,          only: kind_phys
        use set_soilveg_mod,  only: set_soilveg
        use namelist_soilveg
        
        implicit none
      
        integer,              intent(in)  :: me, isot, ivegsrc, nlunit

        real (kind=kind_phys), dimension(:), intent(out) :: pores, resid

        character(len=*),     intent(out) :: errmsg
        integer,              intent(out) :: errflg

        ! Initialize CCPP error handling variables
        errmsg = ''
        errflg = 0

        if (ivegsrc /= 1) then
          errmsg = 'The NOAHMP LSM expects that the ivegsrc physics '//
     &             'namelist parameter is 1. Exiting...'
          errflg = 1
          return
        end if
        if (isot /= 1) then
          errmsg = 'The NOAHMP LSM expects that the isot physics '//
     &             'namelist parameter is 1. Exiting...'
          errflg = 1
          return
        end if

        !--- initialize soil vegetation
        call set_soilveg(me, isot, ivegsrc, nlunit)

        pores (:) = maxsmc (:)
        resid (:) = drysmc (:)
      
      end subroutine noahmpdrv_init

      subroutine noahmpdrv_finalize
      end subroutine noahmpdrv_finalize

!> \ingroup NoahMP_LSM
!! \brief This subroutine is the main CCPP entry point for the NoahMP LSM.
!! \section arg_table_noahmpdrv_run Argument Table
!! \htmlinclude noahmpdrv_run.html
!!
!! \section general_noahmpdrv NoahMP Driver General Algorithm
!!  @{
!!    - Initialize CCPP error handling variables.
!!    - Set a flag to only continue with each grid cell if the fraction of land is non-zero.
!!    - This driver may be called as part of an iterative loop. If called as the first "guess" run, 
!!        save land-related prognostic fields to restore.
!!    - Initialize output variables to zero and prepare variables for input into the NoahMP LSM.
!!    - Call transfer_mp_parameters() to fill a derived datatype for input into the NoahMP LSM.
!!    - Call noahmp_options() to set module-level scheme options for the NoahMP LSM.
!!    - If the vegetation type is ice for the grid cell, call noahmp_options_glacier() to set 
!!        module-level scheme options for NoahMP Glacier and call noahmp_glacier().
!!    - For other vegetation types, call noahmp_sflx(), the entry point of the NoahMP LSM.
!!    - Set output variables from the output of noahmp_glacier() and/or noahmp_sflx().
!!    - Call penman() to calculate potential evaporation.
!!    - Calculate the surface specific humidity and convert surface sensible and latent heat fluxes in W m-2 from their kinematic values.
!!    - If a "guess" run, restore the land-related prognostic fields.
!                                                                       !
!-----------------------------------
      subroutine noahmpdrv_run                                          &
!...................................
!  ---  inputs:
     &     ( im, km, itime, ps, u1, v1, t1, q1, soiltyp, vegtype,       &
     &       sigmaf, sfcemis, dlwflx, dswsfc, snet, delt, tg3, cm, ch,  &
     &       prsl1, prslki, zf, dry, wind, slopetyp,                    &
     &       shdmin, shdmax, snoalb, sfalb, flag_iter, flag_guess,      &
     &       idveg, iopt_crs, iopt_btr, iopt_run, iopt_sfc, iopt_frz,   &
     &       iopt_inf, iopt_rad, iopt_alb, iopt_snf, iopt_tbot,         &
     &       iopt_stc, xlatin, xcoszin, iyrlen, julian,                 &
     &       rainn_mp, rainc_mp, snow_mp, graupel_mp, ice_mp,           &
     &       con_hvap, con_cp, con_jcal, rhoh2o, con_eps, con_epsm1,    &
     &       con_fvirt, con_rd, con_hfus,                               &

!  ---  in/outs:
     &       weasd, snwdph, tskin, tprcp, srflag, smc, stc, slc,        &
     &       canopy, trans, tsurf, zorl,                                &

! --- Noah MP specific

     &       snowxy, tvxy, tgxy, canicexy, canliqxy, eahxy, tahxy, cmxy,&
     &       chxy, fwetxy, sneqvoxy, alboldxy, qsnowxy, wslakexy, zwtxy,&
     &       waxy, wtxy, tsnoxy, zsnsoxy, snicexy, snliqxy, lfmassxy,   &
     &       rtmassxy, stmassxy, woodxy, stblcpxy, fastcpxy, xlaixy,    &
     &       xsaixy, taussxy, smoiseq, smcwtdxy, deeprechxy, rechxy,    &

!  ---  outputs:
     &       sncovr1, qsurf, gflux, drain, evap, hflx, ep, runoff,      &
     &       cmm, chh, evbs, evcw, sbsno, snowc, stm, snohf,            &
     &       smcwlt2, smcref2, wet1, t2mmp, q2mp, errmsg, errflg)     
!
!
      use machine ,   only : kind_phys
!     use date_def,   only : idate
      use funcphys,   only : fpvs

      use module_sf_noahmplsm
      use module_sf_noahmp_glacier
      use noahmp_tables, only : isice_table, co2_table, o2_table,       &
     &                       isurban_table,smcref_table,smcdry_table,   &
     &                       smcmax_table,co2_table,o2_table,           &
     &                       saim_table,laim_table

      implicit none
      
      real(kind=kind_phys), parameter :: a2      = 17.2693882
      real(kind=kind_phys), parameter :: a3      = 273.16
      real(kind=kind_phys), parameter :: a4      = 35.86
      real(kind=kind_phys), parameter :: a23m4   = a2*(a3-a4)
      
      real, parameter                 :: undefined       =  -1.e36

      real                            :: dz8w    =  undefined
      real                            :: dx      =  undefined
      real                            :: qc      =  undefined
      real                            :: foln    =  1.0   ! foliage
      integer                         :: nsoil   = 4   ! hardwired to Noah
      integer                         :: nsnow   = 3   ! max. snow layers
      integer                         :: ist     = 1   ! soil type, 1 soil; 2  lake;  14 is water
      integer                         :: isc     = 4   ! middle day soil color: soil 1-9 lightest

      real(kind=kind_phys), save  :: zsoil(4),sldpth(4)
      data zsoil / -0.1, -0.4, -1.0, -2.0 /
      data sldpth /0.1, 0.3, 0.6, 1.0 /
!     data dzs /0.1, 0.3, 0.6, 1.0 /

!
!  ---  input:
!

      integer, intent(in) :: im, km, itime

      integer, dimension(im), intent(in) :: soiltyp, vegtype, slopetyp

      real (kind=kind_phys), dimension(im), intent(in) :: ps, u1, v1,   &
     &       t1, q1, sigmaf, sfcemis, dlwflx, dswsfc, snet, tg3, cm,    &
     &       ch, prsl1, prslki, wind, shdmin, shdmax,                   &
     &       snoalb, sfalb, zf,                                         &
     &       rainn_mp,rainc_mp,snow_mp,graupel_mp,ice_mp

      logical, dimension(im), intent(in) :: dry

      real (kind=kind_phys),dimension(im),intent(in) :: xlatin,xcoszin

      integer,  intent(in) ::  idveg, iopt_crs,iopt_btr,iopt_run,       &
     &                         iopt_sfc,iopt_frz,iopt_inf,iopt_rad,     &
     &                         iopt_alb,iopt_snf,iopt_tbot,iopt_stc

      real (kind=kind_phys),  intent(in) :: julian
      integer,  intent(in)               :: iyrlen


      real (kind=kind_phys),  intent(in) :: delt
      logical, dimension(im), intent(in) :: flag_iter, flag_guess

      real (kind=kind_phys),  intent(in) :: con_hvap, con_cp, con_jcal, &
     &                      rhoh2o, con_eps, con_epsm1, con_fvirt,      &
     &                      con_rd, con_hfus

!  ---  in/out:
      real (kind=kind_phys), dimension(im), intent(inout) :: weasd,     &
     &       snwdph, tskin, tprcp, srflag, canopy, trans, tsurf, zorl

      real (kind=kind_phys), dimension(im,km), intent(inout) ::         &
     &       smc, stc, slc

      real (kind=kind_phys), dimension(im), intent(inout) :: snowxy,    &
     &                      tvxy,tgxy,canicexy,canliqxy,eahxy,tahxy,    & 
     &                      cmxy,chxy,fwetxy,sneqvoxy,alboldxy,qsnowxy, &
     &                      wslakexy,zwtxy,waxy,wtxy,lfmassxy,rtmassxy, &
     &                      stmassxy,woodxy,stblcpxy,fastcpxy,xlaixy,   &
     &                      xsaixy,taussxy,smcwtdxy,deeprechxy,rechxy  

      real (kind=kind_phys),dimension(im,-2:0),intent(inout) :: tsnoxy 
      real (kind=kind_phys),dimension(im,-2:0),intent(inout) :: snicexy
      real (kind=kind_phys),dimension(im,-2:0),intent(inout) :: snliqxy
      real (kind=kind_phys),dimension(im,1:4), intent(inout) :: smoiseq
      real (kind=kind_phys),dimension(im,-2:4),intent(inout) :: zsnsoxy

      integer, dimension(im)                   :: jsnowxy
      real (kind=kind_phys),dimension(im)      :: snodep
      real (kind=kind_phys),dimension(im,-2:4) :: tsnsoxy
      
!  ---  output:

      real (kind=kind_phys), dimension(im), intent(out) :: sncovr1,     &
     &       qsurf, gflux, drain, evap, hflx, ep, runoff, cmm, chh,     &
     &    evbs, evcw, sbsno, snowc, stm, snohf, smcwlt2, smcref2, wet1
      real (kind=kind_phys), dimension(:), intent(out) :: t2mmp, q2mp

! error messages
      character(len=*), intent(out)    :: errmsg
      integer,          intent(out)    :: errflg

!  ---  locals:
      real (kind=kind_phys), dimension(im) :: rch, rho,                 &
     &       q0, qs1, theta1, tv1,  weasd_old, snwdph_old,              &
     &       tprcp_old, srflag_old, tskin_old, canopy_old

      real (kind=kind_phys), dimension(km) :: et,stsoil,smsoil, slsoil

      real (kind=kind_phys),dimension(im,km) :: smc_old,stc_old,slc_old

      real (kind=kind_phys), dimension(im) :: snow_old, tv_old,tg_old,  &
     &       canice_old,canliq_old,eah_old,tah_old,fwet_old,sneqvo_old, &
     &       albold_old,qsnow_old,wslake_old,zwt_old,wa_old,wt_old,     &
     &       lfmass_old,rtmass_old,stmass_old,wood_old,stblcp_old,      &
     &       fastcp_old,xlai_old,xsai_old,tauss_old,smcwtd_old,         &
     &       deeprech_old,rech_old 

      real(kind=kind_phys),dimension(im,1:4)  :: smoiseq_old
      real(kind=kind_phys),dimension(im,-2:0) :: tsno_old  
      real(kind=kind_phys),dimension(im,-2:0) :: snice_old
      real(kind=kind_phys),dimension(im,-2:0) :: snliq_old 
      real(kind=kind_phys),dimension(im,-2:4) :: zsnso_old
      real(kind=kind_phys),dimension(im,-2:4) :: tsnso_old

 
      real (kind=kind_phys) :: alb, albedo, beta, chx, cmx, cmc,        &
     &       dew, drip, dqsdt2, ec, edir, ett, eta, esnow, etp,         &
     &       flx1, flx2, flx3, ffrozp, lwdn, pc, prcp, ptu, q2,         &
     &       q2sat, solnet, rc, rcs, rct, rcq, rcsoil, rsmin,           &
     &       runoff1, runoff2, runoff3, sfcspd, sfcprs, sfctmp,         &
     &       sfcems, sheat, shdfac, shdmin1d, shdmax1d, smcwlt,         &
     &       smcdry, smcref, smcmax, sneqv, snoalb1d, snowh,            &
     &       snomlt, sncovr, soilw, soilm, ssoil, tsea, th2,            &
     &       xlai, zlvl, swdn, tem, psfc,fdown,t2v,tbot

      real (kind=kind_phys) :: pconv,pnonc,pshcv,psnow,pgrpl,phail
      real (kind=kind_phys) :: lat,cosz,uu,vv,swe
      integer               :: isnowx

      real (kind=kind_phys) :: tvx,tgx,canicex,canliqx,eahx,            &
     &       tahx,fwetx,sneqvox,alboldx,qsnowx,wslakex,zwtx,            &
     &       wax,wtx,lfmassx, rtmassx,stmassx, woodx,stblcpx,           &
     &       fastcpx,xlaix,xsaix,taussx,smcwtdx,deeprechx,rechx,        &
     &       qsfc1d

      real (kind=kind_phys), dimension(-2:0) :: tsnox, snicex, snliqx
      real (kind=kind_phys), dimension(-2:0) :: ficeold
      real (kind=kind_phys), dimension( km ) :: smoiseqx
      real (kind=kind_phys), dimension(-2:4) :: zsnsox
      real (kind=kind_phys), dimension(-2:4) :: tsnsox

      real (kind=kind_phys) :: z0wrf,fsa,fsr,fira,fsh,fcev,fgev,        &
     &                         fctr,ecan,etran,trad,tgb,tgv,t2mv,       &
     &                         t2mb,q2v,q2b,runsrf,runsub,apar,         &
     &                         psn,sav,sag,fsno,nee,gpp,npp,fveg,       &
     &                         qsnbot,ponding,ponding1,ponding2,        &
     &                         rssun,rssha,bgap,wgap,chv,chb,emissi,    &
     &                         shg,shc,shb,evg,evb,ghv,ghb,irg,irc,     &
     &                         irb,tr,evc,chleaf,chuc,chv2,chb2,        &
     &                         fpice,pahv,pahg,pahb,pah,co2pp,o2pp,ch2b

      integer :: i, k, ice, stype, vtype ,slope,nroot,couple
      logical :: flag(im)
      logical :: snowng,frzgra
      
      !  ---  local derived constants:

      real(kind=kind_phys) :: cpinv, hvapi, convrad, elocp
      
      type(noahmp_parameters) :: parameters

!
!===> ...  begin here
!     
      cpinv   = 1.0/con_cp
      hvapi   = 1.0/con_hvap
      convrad = con_jcal*1.e4/60.0
      elocp   = con_hvap/con_cp

! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

!  --- ...  set flag for land points

      do i = 1, im
        flag(i) = dry(i)
      enddo

!  --- ...  save land-related prognostic fields for guess run

      do i = 1, im
        if (flag(i) .and. flag_guess(i)) then
          weasd_old(i)  = weasd(i)
          snwdph_old(i) = snwdph(i)
          tskin_old(i)  = tskin(i)
          canopy_old(i) = canopy(i)
          tprcp_old(i)  = tprcp(i)
          srflag_old(i) = srflag(i)
!
!
          snow_old(i)   = snowxy(i)
          tv_old(i)     = tvxy(i)
          tg_old(i)     = tgxy(i)
          canice_old(i) = canicexy(i)
          canliq_old(i) = canliqxy(i)
          eah_old(i)    = eahxy(i)
          tah_old(i)    = tahxy(i)
          fwet_old(i)   = fwetxy(i)
          sneqvo_old(i) = sneqvoxy(i) 
          albold_old(i) = alboldxy(i)
          qsnow_old(i)  = qsnowxy(i)
          wslake_old(i) = wslakexy(i)
          zwt_old(i)    = zwtxy(i)
          wa_old(i)     = waxy(i)
          wt_old(i)     = wtxy(i)
          lfmass_old(i) = lfmassxy(i)
          rtmass_old(i) = rtmassxy(i)
          stmass_old(i) = stmassxy(i)
          wood_old(i)   = woodxy(i)
          stblcp_old(i) = stblcpxy(i)
          fastcp_old(i) = fastcpxy(i)
          xlai_old(i)   = xlaixy(i)
          xsai_old(i)   = xsaixy(i)
          tauss_old(i)  = taussxy(i)
          smcwtd_old(i) = smcwtdxy(i)
          rech_old(i)   = rechxy(i)

          deeprech_old(i) = deeprechxy(i)
!
          do k = 1, km
           smc_old(i,k) = smc(i,k)
           stc_old(i,k) = stc(i,k)
           slc_old(i,k) = slc(i,k)
          enddo

!
          do k = 1, km
            smoiseq_old(i,k) = smoiseq(i,k)
          enddo

          do k = -2,0
           tsno_old(i,k)  = tsnoxy(i,k)
           snice_old(i,k) = snicexy(i,k)
           snliq_old(i,k) = snliqxy(i,k)
          enddo

          do k = -2,4
            zsnso_old (i,k) = zsnsoxy(i,k)
          enddo

        endif
      enddo

!
!   call to init MP options
!
!    &_________________________________________________________________ &

!  --- ...  initialization block

      do i = 1, im
        if (flag_iter(i) .and. flag(i)) then
          ep(i)     = 0.0
          evap (i)  = 0.0
          hflx (i)  = 0.0
          gflux(i)  = 0.0
          drain(i)  = 0.0
          canopy(i) = max(canopy(i), 0.0)

          evbs (i)  = 0.0
          evcw (i)  = 0.0
          trans(i)  = 0.0
          sbsno(i)  = 0.0
          snowc(i)  = 0.0
          snohf(i)  = 0.0
        endif
      enddo

!  --- ...  initialize variables 

      do i = 1, im
        if (flag_iter(i) .and. flag(i)) then
          q0(i)   = max(q1(i), 1.e-8)   !* q1=specific humidity at level 1 (kg/kg)
          theta1(i) = t1(i) * prslki(i) !* adiabatic temp at level 1 (k)

          tv1(i) = t1(i) * (1.0 + con_fvirt*q0(i))
          rho(i) = prsl1(i) / (con_rd * tv1(i))
          qs1(i) = fpvs( t1(i) )        !* qs1=sat. humidity at level 1 (kg/kg)
          qs1(i) = con_eps*qs1(i) / (prsl1(i) + con_epsm1*qs1(i))
          qs1(i) = max(qs1(i), 1.e-8)
          q0 (i) = min(qs1(i), q0(i))

          if (vegtype(i) == isice_table ) then
              if (weasd(i) < 0.1) then
                 weasd(i)  = 0.1
               endif
          endif                                       

         endif                                       
        enddo

!  --- ...  noah: prepare variables to run noah lsm
!   1. configuration information (c):
!      ------------------------------
!    couple  - couple-uncouple flag (=1: coupled, =0: uncoupled) 
!    ffrozp  - fraction for snow-rain (1.=snow, 0.=rain,  0-1 mixed))                
!    ice     - sea-ice flag (=1: sea-ice, =0: land)
!    dt      - timestep (sec) (dt should not exceed 3600 secs) = delt
!    zlvl    - height (m) above ground of atmospheric forcing variables
!    nsoil   - number of soil layers (at least 2)
!    sldpth  - the thickness of each soil layer (m)

      do i = 1, im

        if (flag_iter(i) .and. flag(i)) then


          couple = 1

          ice    = 0
          nsoil  = km
          snowng = .false.          
          frzgra = .false.


!         if     (srflag(i) == 1.0) then  ! snow phase
!           ffrozp = 1.0
!         elseif (srflag(i) == 0.0) then  ! rain phase
!           ffrozp = 0.0
!         endif
! use srflag directly to allow fractional rain/snow
          ffrozp = srflag(i)

          zlvl = zf(i)

!   2. forcing data (f):
!      -----------------
!    lwdn    - lw dw radiation flux (w/m2)
!    solnet  - net sw radiation flux (dn-up) (w/m2)
!    sfcprs  - pressure at height zlvl above ground (pascals)
!    prcp    - precip rate (kg m-2 s-1)
!    sfctmp  - air temperature (k) at height zlvl above ground
!    th2     - air potential temperature (k) at height zlvl above ground
!    q2      - mixing ratio at height zlvl above ground (kg kg-1)

          lat    = xlatin(i)                 ! in radian
          cosz   = xcoszin(i)

          lwdn   = dlwflx(i)         !..downward lw flux at sfc in w/m2
          swdn   = dswsfc(i)         !..downward sw flux at sfc in w/m2
          solnet = snet(i)           !..net sw rad flx (dn-up) at sfc in w/m2
          sfcems = sfcemis(i)

          sfctmp = t1(i)  
          sfcprs = prsl1(i) 
          psfc   = ps(i) 
          prcp   = rhoh2o * tprcp(i) / delt

        if (prcp > 0.0) then
          if (ffrozp > 0.0) then                   ! rain/snow flag, one condition is enough?
            snowng = .true.
            qsnowxy(i) = ffrozp * prcp/10.0                 !still use rho water?
          else
            if (sfctmp <= 275.15) frzgra = .true.  
          endif
        endif

          th2    = theta1(i)
          q2     = q0(i)

!   3. other forcing (input) data (i):
!      ------------------------------
!    sfcspd  - wind speed (m s-1) at height zlvl above ground
!    q2sat   - sat mixing ratio at height zlvl above ground (kg kg-1)
!    dqsdt2  - slope of sat specific humidity curve at t=sfctmp (kg kg-1 k-1)

          uu     = u1(i)
          vv     = v1(i)

          sfcspd = wind(i)
          q2sat =  qs1(i)
          dqsdt2 = q2sat * a23m4/(sfctmp-a4)**2

!   4. canopy/soil characteristics (s):
!      --------------------------------
!    vegtyp  - vegetation type (integer index)                       -> vtype
!    soiltyp - soil type (integer index)                             -> stype
!    slopetyp- class of sfc slope (integer index)                    -> slope
!    shdfac  - areal fractional coverage of green vegetation (0.0-1.0)
!    shdmin  - minimum areal fractional coverage of green vegetation -> shdmin1d
!    ptu     - photo thermal unit (plant phenology for annuals/crops)
!    alb     - backround snow-free surface albedo (fraction)
!    snoalb  - upper bound on maximum albedo over deep snow          -> snoalb1d
!    tbot    - bottom soil temperature (local yearly-mean sfc air temp)

          vtype = vegtype(i)
          stype = soiltyp(i)
          slope = slopetyp(i)
          shdfac= sigmaf(i)

          shdmin1d = shdmin(i)   
          shdmax1d = shdmax(i)     
          snoalb1d = snoalb(i)

          alb  = sfalb(i)

          tbot = tg3(i)
          ptu  = 0.0


          cmc = canopy(i)/1000.              ! convert from mm to m
          tsea = tsurf(i)                    ! clu_q2m_iter

          snowh = snwdph(i) * 0.001         ! convert from mm to m
          sneqv = weasd(i)  * 0.001         ! convert from mm to m



!   5. history (state) variables (h):
!      ------------------------------
!    cmc     - canopy moisture content (m)
!    t1      - ground/canopy/snowpack) effective skin temperature (k)   -> tsea
!    stc(nsoil) - soil temp (k)                                         -> stsoil
!    smc(nsoil) - total soil moisture content (volumetric fraction)     -> smsoil
!    sh2o(nsoil)- unfrozen soil moisture content (volumetric fraction)  -> slsoil
!    snowh   - actual snow depth (m)
!    sneqv   - liquid water-equivalent snow depth (m)
!    albedo  - surface albedo including snow effect (unitless fraction)
!    ch      - surface exchange coefficient for heat and moisture (m s-1) -> chx
!    cm      - surface exchange coefficient for momentum (m s-1)          -> cmx

          isnowx   = nint(snowxy(i))
          tvx      = tvxy(i)
          tgx      = tgxy(i)
          canliqx  = canliqxy(i)    !in mm
          canicex  = canicexy(i)

          eahxy(i) = (ps(i)*q2)/(0.622+q2) ! use q0 to reinit;
          eahx     = eahxy(i)
          tahx     = tahxy(i)

          co2pp    = co2_table * sfcprs        
          o2pp     = o2_table  * sfcprs
          fwetx    = fwetxy(i)

          sneqvox  = sneqvoxy(i)
          alboldx  = alboldxy(i) 

          qsnowx   = qsnowxy(i)
          wslakex  = wslakexy(i)

          zwtx     = zwtxy(i)
          wax      = waxy(i)
          wtx      = waxy(i)     

          do k = -2,0
           tsnsoxy(i,k) = tsnoxy(i,k)
          enddo

          do k = 1,4
           tsnsoxy(i,k) = stc(i,k)
          enddo

         do k = -2,0
              snicex(k)  = snicexy(i,k)       ! in k/m3; mm
              snliqx(k)  = snliqxy(i,k)       ! in k/m3; mm
              tsnox (k)  = tsnoxy(i,k)

              ficeold(k) = 0.0                ! derived

              if (snicex(k) > 0.0 ) then
               ficeold(k) = snicex(k) /(snicex(k)+snliqx(k)) 
             
              endif
         enddo

          do k = -2, km
            zsnsox(k)  = zsnsoxy(i,k)
            tsnsox(k)  = tsnsoxy(i,k)
          enddo

          lfmassx  = lfmassxy(i)
          rtmassx  = rtmassxy(i)
          stmassx  = stmassxy(i)

          woodx    = woodxy(i)
          stblcpx  = stblcpxy(i)
          fastcpx  = fastcpxy(i)

          xsaix    = xsaixy(i)
          xlaix    = xlaixy(i)

          taussx   = taussxy(i)  

          qsfc1d   = undefined              ! derive later, it is an in/out?
          swe      = weasd(i)

          do k = 1, km
            smoiseqx(k)  = smoiseq(i,k)
          enddo

          smcwtdx        = smcwtdxy(i)
          rechx          = rechxy(i)
          deeprechx      = deeprechxy(i)
!--
!   the optional details for precip
!--

!         pconv = 0.                      !     convective - may introduce later
!         pnonc = (1 - ffrozp) * prcp                    !     large scale total in mm/s;
!         pshcv = 0.
!         psnow = ffrozp * prcp /10.0     !      snow = qsnowx?
!         pgrpl = 0.
!         phail = 0.
          pnonc = rainn_mp(i)
          pconv = rainc_mp(i)
          pshcv = 0.
          psnow = snow_mp(i)
          pgrpl = graupel_mp(i)
          phail = ice_mp(i)
!
!-- old
!
          do k = 1, km
!           stsoil(k) = stc(i,k)
            smsoil(k) = smc(i,k)
            slsoil(k) = slc(i,k)
          enddo

          snowh = snwdph(i) * 0.001         ! convert from mm to m

          if (swe /= 0.0 .and. snowh == 0.0) then
            snowh = 10.0 * swe /1000.0
          endif

          chx    = chxy(i) ! maybe chxy  
          cmx    = cmxy(i)

          chh(i) = ch(i)  * wind(i) * rho(i)
          cmm(i) = cm(i)  * wind(i)



       call transfer_mp_parameters(vtype,stype,slope,isc,parameters)

       call noahmp_options(idveg ,iopt_crs,iopt_btr,iopt_run,iopt_sfc,  &
     & iopt_frz,iopt_inf,iopt_rad,iopt_alb,iopt_snf,iopt_tbot,iopt_stc)

        if ( vtype == isice_table )  then

          ice = -1
          tbot = min(tbot,263.15)

         call noahmp_options_glacier                                    &
     &   (idveg  ,iopt_crs  ,iopt_btr, iopt_run ,iopt_sfc ,iopt_frz,    &
     &   iopt_inf ,iopt_rad ,iopt_alb ,iopt_snf ,iopt_tbot, iopt_stc )

       call noahmp_glacier (                                            &
     &             i       ,1       ,cosz    ,nsnow   ,nsoil   ,delt  , & ! in : time/space/model-related
     &             sfctmp  ,sfcprs  ,uu      ,vv      ,q2      ,swdn  , & ! in : forcing
     &             prcp    ,lwdn    ,tbot    ,zlvl    ,ficeold ,zsoil , & ! in : forcing
     &             qsnowx  ,sneqvox ,alboldx ,cmx     ,chx     ,isnowx, & ! in/out :sneqvox + alboldx -LST 
     &             swe  ,smsoil  ,zsnsox     ,snowh  ,snicex ,snliqx ,  & ! in/out : sneqvx + snowhx are avgd
     &             tgx     ,tsnsox  ,slsoil  ,taussx  ,qsfc1d         , & ! in/out : 
     &             fsa     ,fsr     ,fira    ,fsh     ,fgev  ,ssoil   , & ! out : 
     &             trad    ,edir    ,runsrf  ,runsub  ,sag   ,albedo  , & ! out : albedo is surface albedo
     &             qsnbot  ,ponding ,ponding1,ponding2,t2mb  ,q2b     , & ! out :
#ifdef CCPP
     &             emissi  ,fpice   ,ch2b    ,esnow, errmsg, errflg )
#else
     &             emissi  ,fpice   ,ch2b    ,esnow )
#endif

#ifdef CCPP
       if (errflg /= 0) return
#endif
!
! in/out and outs
!

         fsno    = 1.0

         tvx     = undefined  
         canicex = undefined
         canliqx = undefined
         eahx    = undefined
         tahx    = undefined

         fwetx   = undefined
         wslakex = undefined
         zwtx    = undefined
         wax     = undefined
         wtx     = undefined

         lfmassx = undefined
         rtmassx = undefined
         stmassx = undefined
         woodx   = undefined
         stblcpx = undefined
         fastcpx = undefined
         xlaix   = undefined
         xsaix   = undefined

         smcwtdx = 0.0
         rechx   = 0.0
         deeprechx       = 0.0

         do k = 1,4
         smoiseqx(k)   = smsoil(k)
         enddo

        fctr    = undefined
        fcev    = undefined

        z0wrf  = 0.002
  
        eta    = fgev
        t2mmp(i) = t2mb 
        q2mp(i) = q2b 
!
! Non-glacial case
!
        else
                 ice = 0 

!        write(*,*)'tsnsox(1)=',tsnsox,'tgx=',tgx
       call noahmp_sflx (parameters                                    ,&
     &        i       , 1       , lat     , iyrlen  , julian  , cosz   ,& ! in : time/space-related
     &        delt    , dx      , dz8w    , nsoil   , zsoil   , nsnow  ,& ! in : model configuration 
     &        shdfac  , shdmax1d, vtype   , ice     , ist              ,& ! in : vegetation/soil 
     &        smoiseqx                                                 ,& ! in 
     &        sfctmp  , sfcprs  , psfc    , uu      , vv      , q2     ,& ! in : forcing
     &        qc      , swdn    , lwdn                                 ,& ! in : forcing
     &        pconv   , pnonc   , pshcv   , psnow   , pgrpl   , phail  ,& ! in : forcing
     &        tbot    , co2pp   , o2pp    , foln    , ficeold , zlvl   ,& ! in : forcing
     &        alboldx , sneqvox                                        ,& ! in/out : 
     &        tsnsox  , slsoil  , smsoil  , tahx    , eahx    , fwetx  ,& ! in/out : 
     &        canliqx , canicex , tvx     , tgx     , qsfc1d  , qsnowx ,& ! in/out : 
     &        isnowx  , zsnsox  , snowh   , swe     , snicex  , snliqx ,& ! in/out : 
     &        zwtx    , wax     , wtx     , wslakex , lfmassx , rtmassx,& ! in/out : 
     &        stmassx , woodx   , stblcpx , fastcpx , xlaix   ,xsaix   ,& ! in/out : 
     &        cmx     , chx     , taussx                               ,& ! in/out : 
     &        smcwtdx ,deeprechx, rechx                                ,& ! in/out :
     &        z0wrf                                                    ,& ! out
     &        fsa     , fsr     , fira    , fsh     , ssoil   , fcev   ,& ! out : 
     &        fgev    , fctr    , ecan    , etran   , edir    , trad   ,& ! out :
     &        tgb     , tgv     , t2mv    , t2mb    , q2v     , q2b    ,& ! out :
     &        runsrf  , runsub  , apar    , psn     , sav     , sag    ,& ! out :
     &        fsno    , nee     , gpp     , npp     , fveg    , albedo ,& ! out :
     &        qsnbot  , ponding , ponding1, ponding2, rssun   , rssha  ,& ! out :
     &        bgap    , wgap    , chv     , chb     , emissi           ,& ! out :
     &        shg     , shc     , shb     , evg     , evb     , ghv    ,&! out :
     &        ghb     , irg     , irc     , irb     , tr      , evc    ,& ! out :
     &        chleaf  , chuc    , chv2    , chb2    , fpice   , pahv   ,& ! out
#ifdef CCPP
     &        pahg    , pahb    , pah     , esnow, errmsg, errflg   )
#else
     &        pahg    , pahb    , pah     , esnow   )
#endif

#ifdef CCPP
       if (errflg /= 0) return
#endif

       eta  = fcev + fgev + fctr     ! the flux w/m2

       t2mmp(i) = t2mv*fveg+t2mb*(1-fveg) 
       q2mp(i) = q2v*fveg+q2b*(1-fveg)

      endif          ! glacial split ends

!
! mp in/out
!
             snowxy   (i)                = float(isnowx)
             tvxy     (i)                = tvx
             tgxy     (i)                = tgx
             canliqxy (i)                = canliqx
             canicexy (i)                = canicex
             eahxy    (i)                = eahx
             tahxy    (i)                = tahx

             cmxy     (i)                = cmx
             chxy     (i)                = chx

             fwetxy   (i)                = fwetx
             sneqvoxy (i)                = sneqvox
             alboldxy (i)                = alboldx
             qsnowxy  (i)                = qsnowx

             wslakexy (i)                = wslakex
             zwtxy    (i)                = zwtx
             waxy     (i)                = wax
             wtxy     (i)                = wtx

             do k = -2,0
             tsnoxy   (i,k) = tsnsox(k)
             snicexy  (i,k) = snicex (k)
             snliqxy  (i,k) = snliqx (k)
             enddo

             do k = -2,4
             zsnsoxy  (i,k) = zsnsox(k)
             enddo

             lfmassxy (i)                = lfmassx
             rtmassxy (i)                = rtmassx
             stmassxy (i)                = stmassx
             woodxy   (i)                = woodx
             stblcpxy (i)                = stblcpx
             fastcpxy (i)                = fastcpx

             xlaixy   (i)                = xlaix
             xsaixy   (i)                = xsaix

             taussxy  (i)                = taussx

             rechxy   (i)                = rechx
             deeprechxy(i)               = deeprechx
             smcwtdxy(i)                 = smcwtdx
             smoiseq(i,1:4)              = smoiseqx(1:4)

!
! generic in/outs
!
          do k = 1, km
            stc(i,k) = tsnsox(k) 
            smc(i,k) = smsoil(k)
            slc(i,k) = slsoil(k)
          enddo

          canopy(i)  = canicex + canliqx
          weasd(i)   = swe
          snwdph(i)  = snowh * 1000.0

!         write(*,*) 'swe,snowh,can'
!         write (*,*) swe,snowh*1000.0,canopy(i)
!
          smcmax = smcmax_table(stype) 
          smcref = smcref_table(stype)
          smcwlt = smcdry_table(stype)
!
! outs
!
          wet1(i)    = smsoil(1) / smcmax 
          smcwlt2(i) = smcwlt
          smcref2(i) = smcref

          runoff(i)  = runsrf
          drain(i)   = runsub

          zorl(i)    = z0wrf * 100.0

          sncovr1(i) = fsno
          snowc  (i) = fsno

          sbsno(i)   = esnow
          gflux(i)   = -1.0*ssoil
          hflx(i)    = fsh
          evbs(i)    = fgev
          evcw(i)    = fcev
          trans(i)   = fctr
          evap(i)    = eta

!         write(*,*) 'vtype, stype are',vtype,stype
!         write(*,*) 'fsh,gflx,eta',fsh,ssoil,eta
!         write(*,*) 'esnow,runsrf,runsub',esnow,runsrf,runsub
!         write(*,*) 'evbs,evcw,trans',fgev,fcev,fctr
!         write(*,*) 'snowc',fsno

          tsurf(i)   = trad

          stm(i) = (0.1*smsoil(1)+0.3*smsoil(2)+0.6*smsoil(3)+           &
     &              1.0*smsoil(4))*1000.0  ! unit conversion from m to kg m-2
!
          snohf (i) = qsnbot * con_hfus  ! only part of it but is diagnostic
!         write(*,*) 'snohf',snohf(i)

          fdown     = fsa + lwdn
          t2v       = sfctmp * (1.0 + 0.61*q2)
!         ssoil     = -1.0 *ssoil

       call penman (sfctmp,sfcprs,chx,t2v,th2,prcp,fdown,ssoil,         &
     &   q2,q2sat,etp,snowng,frzgra,ffrozp,dqsdt2,emissi,fsno)

          ep(i) = etp

        endif   ! end if_flag_iter_and_flag_block
      enddo   ! end do_i_loop

!   --- ...  compute qsurf (specific humidity at sfc)

      do i = 1, im
        if (flag_iter(i) .and. flag(i)) then
          rch(i)   = rho(i) * con_cp * ch(i) * wind(i)
          qsurf(i) = q1(i)  + evap(i) / (elocp * rch(i))
        endif
      enddo

      do i = 1, im
        if (flag_iter(i) .and. flag(i)) then
          tem     = 1.0 / rho(i)
          hflx(i) = hflx(i) * tem * cpinv
          evap(i) = evap(i) * tem * hvapi
        endif
      enddo

!  --- ...  restore land-related prognostic fields for guess run

      do i = 1, im
        if (flag(i)) then
          if (flag_guess(i)) then
            weasd(i)  = weasd_old(i)
            snwdph(i) = snwdph_old(i)
            tskin(i)  = tskin_old(i)
            canopy(i) = canopy_old(i)
            tprcp(i)  = tprcp_old(i)
            srflag(i) = srflag_old(i)


            snowxy(i) =  snow_old(i)
            tvxy(i)   =  tv_old(i)
            tgxy(i)   =  tg_old(i)

            canicexy(i)   = canice_old(i)
            canliqxy(i)   = canliq_old(i)
            eahxy(i)      = eah_old(i)
            tahxy(i)      = tah_old(i)
            fwetxy(i)     = fwet_old(i)
            sneqvoxy(i)   = sneqvo_old(i)
            alboldxy(i)   = albold_old(i)
            qsnowxy(i)    = qsnow_old(i)
            wslakexy(i)   = wslake_old(i)
            zwtxy(i)      = zwt_old(i)
            waxy(i)       = wa_old(i)
            wtxy(i)       = wt_old(i)
            lfmassxy(i)   = lfmass_old(i)
            rtmassxy(i)   = rtmass_old(i)
            stmassxy(i)   = stmass_old(i)
            woodxy(i)     = wood_old(i)
            stblcpxy(i)   = stblcp_old(i)
            fastcpxy(i)   = fastcp_old(i)
            xlaixy(i)     = xlai_old(i)
            xsaixy(i)     = xsai_old(i)
            taussxy(i)    = tauss_old(i)
            smcwtdxy(i)   = smcwtd_old(i)
            deeprechxy(i) = deeprech_old(i)
            rechxy(i)     = rech_old(i)

            do k = 1, km
              smc(i,k) = smc_old(i,k)
              stc(i,k) = stc_old(i,k)
              slc(i,k) = slc_old(i,k)
            enddo
!
            do k = 1, km
              smoiseq(i,k) = smoiseq_old(i,k)
            enddo

            do k = -2,0
              tsnoxy(i,k)  = tsno_old(i,k)
              snicexy(i,k) = snice_old(i,k)
              snliqxy(i,k) = snliq_old(i,k)
            enddo

            do k = -2,4
              zsnsoxy(i,k) =  zsnso_old(i,k)
            enddo
       else
            tskin(i) = tsurf(i)    
          endif
        endif
      enddo
!
      return
!...................................
      end subroutine noahmpdrv_run
!> @}
!-----------------------------------

!> \ingroup NoahMP_LSM
!! \brief This subroutine fills in a derived data type of type noahmp_parameters with data
!! from the module \ref noahmp_tables.
      subroutine transfer_mp_parameters (vegtype,soiltype,slopetype,    &
     &                                          soilcolor,parameters)
     
        use noahmp_tables
        use module_sf_noahmplsm
      
        implicit none
      
        integer, intent(in)    :: vegtype
        integer, intent(in)    :: soiltype
        integer, intent(in)    :: slopetype
        integer, intent(in)    :: soilcolor
          
        type (noahmp_parameters), intent(out) :: parameters
          
        real    :: refdk
        real    :: refkdt
        real    :: frzk
        real    :: frzfact
      
        parameters%iswater   =  iswater_table
        parameters%isbarren  =  isbarren_table
        parameters%isice     =  isice_table
        parameters%eblforest =  eblforest_table
      
!-----------------------------------------------------------------------&
        parameters%urban_flag = .false.
        if( vegtype == isurban_table .or. vegtype == 31                 &
     &         .or.vegtype  == 32 .or. vegtype == 33) then
           parameters%urban_flag = .true.
        endif
      
!------------------------------------------------------------------------------------------!
! transfer veg parameters
!------------------------------------------------------------------------------------------!
      
        parameters%ch2op  =  ch2op_table(vegtype)       !maximum intercepted h2o per unit lai+sai (mm)
        parameters%dleaf  =  dleaf_table(vegtype)       !characteristic leaf dimension (m)
        parameters%z0mvt  =  z0mvt_table(vegtype)       !momentum roughness length (m)
        parameters%hvt    =    hvt_table(vegtype)       !top of canopy (m)
        parameters%hvb    =    hvb_table(vegtype)       !bottom of canopy (m)
        parameters%den    =    den_table(vegtype)       !tree density (no. of trunks per m2)
        parameters%rc     =     rc_table(vegtype)       !tree crown radius (m)
        parameters%mfsno  =  mfsno_table(vegtype)       !snowmelt m parameter ()
        parameters%saim   =   saim_table(vegtype,:)     !monthly stem area index, one-sided
        parameters%laim   =   laim_table(vegtype,:)     !monthly leaf area index, one-sided
        parameters%sla    =    sla_table(vegtype)       !single-side leaf area per kg [m2/kg]
        parameters%dilefc = dilefc_table(vegtype)       !coeficient for leaf stress death [1/s]
        parameters%dilefw = dilefw_table(vegtype)       !coeficient for leaf stress death [1/s]
        parameters%fragr  =  fragr_table(vegtype)       !fraction of growth respiration  !original was 0.3 
        parameters%ltovrc = ltovrc_table(vegtype)       !leaf turnover [1/s]
      
        parameters%c3psn  =  c3psn_table(vegtype)       !photosynthetic pathway: 0. = c4, 1. = c3
        parameters%kc25   =   kc25_table(vegtype)       !co2 michaelis-menten constant at 25c (pa)
        parameters%akc    =    akc_table(vegtype)       !q10 for kc25
        parameters%ko25   =   ko25_table(vegtype)       !o2 michaelis-menten constant at 25c (pa)
        parameters%ako    =    ako_table(vegtype)       !q10 for ko25
        parameters%vcmx25 = vcmx25_table(vegtype)       !maximum rate of carboxylation at 25c (umol co2/m**2/s)
        parameters%avcmx  =  avcmx_table(vegtype)       !q10 for vcmx25
        parameters%bp     =     bp_table(vegtype)       !minimum leaf conductance (umol/m**2/s)
        parameters%mp     =     mp_table(vegtype)       !slope of conductance-to-photosynthesis relationship
        parameters%qe25   =   qe25_table(vegtype)       !quantum efficiency at 25c (umol co2 / umol photon)
        parameters%aqe    =    aqe_table(vegtype)       !q10 for qe25
        parameters%rmf25  =  rmf25_table(vegtype)       !leaf maintenance respiration at 25c (umol co2/m**2/s)
        parameters%rms25  =  rms25_table(vegtype)       !stem maintenance respiration at 25c (umol co2/kg bio/s)
        parameters%rmr25  =  rmr25_table(vegtype)       !root maintenance respiration at 25c (umol co2/kg bio/s)
        parameters%arm    =    arm_table(vegtype)       !q10 for maintenance respiration
        parameters%folnmx = folnmx_table(vegtype)       !foliage nitrogen concentration when f(n)=1 (%)
        parameters%tmin   =   tmin_table(vegtype)       !minimum temperature for photosynthesis (k)
      
        parameters%xl     =     xl_table(vegtype)       !leaf/stem orientation index
        parameters%rhol   =   rhol_table(vegtype,:)     !leaf reflectance: 1=vis, 2=nir
        parameters%rhos   =   rhos_table(vegtype,:)     !stem reflectance: 1=vis, 2=nir
        parameters%taul   =   taul_table(vegtype,:)     !leaf transmittance: 1=vis, 2=nir
        parameters%taus   =   taus_table(vegtype,:)     !stem transmittance: 1=vis, 2=nir
      
        parameters%mrp    =    mrp_table(vegtype)       !microbial respiration parameter (umol co2 /kg c/ s)
        parameters%cwpvt  =  cwpvt_table(vegtype)       !empirical canopy wind parameter
      
        parameters%wrrat  =  wrrat_table(vegtype)       !wood to non-wood ratio
        parameters%wdpool = wdpool_table(vegtype)       !wood pool (switch 1 or 0) depending on woody or not [-]
        parameters%tdlef  =  tdlef_table(vegtype)       !characteristic t for leaf freezing [k]
      
        parameters%nroot  =  nroot_table(vegtype)       !number of soil layers with root present
        parameters%rgl    =    rgl_table(vegtype)       !parameter used in radiation stress function
        parameters%rsmin  =     rs_table(vegtype)       !minimum stomatal resistance [s m-1]
        parameters%hs     =     hs_table(vegtype)       !parameter used in vapor pressure deficit function
        parameters%topt   =   topt_table(vegtype)       !optimum transpiration air temperature [k]
        parameters%rsmax  =  rsmax_table(vegtype)       !maximal stomatal resistance [s m-1]
      
!------------------------------------------------------------------------------------------!
! transfer rad parameters
!------------------------------------------------------------------------------------------!
      
         parameters%albsat    = albsat_table(soilcolor,:)
         parameters%albdry    = albdry_table(soilcolor,:)
         parameters%albice    = albice_table
         parameters%alblak    = alblak_table               
         parameters%omegas    = omegas_table
         parameters%betads    = betads_table
         parameters%betais    = betais_table
         parameters%eg        = eg_table
      
!------------------------------------------------------------------------------------------!
! transfer global parameters
!------------------------------------------------------------------------------------------!
      
         parameters%co2       =    co2_table
         parameters%o2        =     o2_table
         parameters%timean    = timean_table
         parameters%fsatmx    = fsatmx_table
         parameters%z0sno     =  z0sno_table
         parameters%ssi       =    ssi_table
         parameters%swemx     =  swemx_table
      
! ----------------------------------------------------------------------
!  transfer soil parameters
! ----------------------------------------------------------------------
      
          parameters%bexp   = bexp_table   (soiltype)
          parameters%dksat  = dksat_table  (soiltype)
          parameters%dwsat  = dwsat_table  (soiltype)
          parameters%f1     = f1_table     (soiltype)
          parameters%psisat = psisat_table (soiltype)
          parameters%quartz = quartz_table (soiltype)
          parameters%smcdry = smcdry_table (soiltype)
          parameters%smcmax = smcmax_table (soiltype)
          parameters%smcref = smcref_table (soiltype)
          parameters%smcwlt = smcwlt_table (soiltype)
      
! ----------------------------------------------------------------------
! transfer genparm parameters
! ----------------------------------------------------------------------
          parameters%csoil  = csoil_table
          parameters%zbot   = zbot_table
          parameters%czil   = czil_table
      
          frzk   = frzk_table
          refdk  = refdk_table
          refkdt = refkdt_table
          parameters%kdt    = refkdt * parameters%dksat / refdk
          parameters%slope  = slope_table(slopetype)
      
          if(parameters%urban_flag)then  ! hardcoding some urban parameters for soil
             parameters%smcmax = 0.45 
             parameters%smcref = 0.42 
             parameters%smcwlt = 0.40 
             parameters%smcdry = 0.40 
             parameters%csoil  = 3.e6
          endif
      
      ! adjust frzk parameter to actual soil type: frzk * frzfact
      
!-----------------------------------------------------------------------&
          if(soiltype /= 14) then
            frzfact = (parameters%smcmax / parameters%smcref)           &
     &      * (0.412 / 0.468)
            parameters%frzx = frzk * frzfact
          end if
      
       end subroutine transfer_mp_parameters

!-----------------------------------------------------------------------&

!> \ingroup NoahMP_LSM
!! brief Calculate potential evaporation for the current point. Various
!! partial sums/products are also calculated and passed back to the
!! calling routine for later use.
      subroutine penman (sfctmp,sfcprs,ch,t2v,th2,prcp,fdown,ssoil,     &
     &                   q2,q2sat,etp,snowng,frzgra,ffrozp,             &
     &                   dqsdt2,emissi_in,sncovr)
 
! etp is calcuated right after ssoil

! ----------------------------------------------------------------------
! subroutine penman
! ----------------------------------------------------------------------
      implicit none
      logical, intent(in)     :: snowng, frzgra
      real, intent(in)        :: ch, dqsdt2,fdown,prcp,ffrozp,          &
     &                           q2, q2sat,ssoil, sfcprs, sfctmp,       &
     &                           t2v, th2,emissi_in,sncovr
      real, intent(out)       :: etp
      real                    :: epsca,flx2,rch,rr,t24
      real                    :: a, delta, fnet,rad,rho,emissi,elcp1,lvs

      real, parameter :: elcp = 2.4888e+3, lsubc = 2.501000e+6,cp = 1004.6
      real, parameter :: lsubs = 2.83e+6, rd = 287.05, cph2o = 4.1855e+3
      real, parameter :: cpice = 2.106e+3, lsubf   = 3.335e5  
      real, parameter :: sigma = 5.6704e-8

! ----------------------------------------------------------------------
! executable code begins here:
! ----------------------------------------------------------------------
! ----------------------------------------------------------------------
! prepare partial quantities for penman equation.
! ----------------------------------------------------------------------
        emissi=emissi_in
!       elcp1  = (1.0-sncovr)*elcp  + sncovr*elcp*lsubs/lsubc
        lvs    = (1.0-sncovr)*lsubc + sncovr*lsubs

      flx2 = 0.0
      delta = elcp * dqsdt2
!     delta = elcp1 * dqsdt2
      t24 = sfctmp * sfctmp * sfctmp * sfctmp
       rr = t24 * 6.48e-8 / (sfcprs * ch) + 1.0
!     rr = emissi*t24 * 6.48e-8 / (sfcprs * ch) + 1.0
      rho = sfcprs / (rd * t2v)

! ----------------------------------------------------------------------
! adjust the partial sums / products with the latent heat
! effects caused by falling precipitation.
! ----------------------------------------------------------------------
      rch = rho * cp * ch
      if (.not. snowng) then
         if (prcp >  0.0) rr = rr + cph2o * prcp / rch
      else
! ---- ...  fractional snowfall/rainfall
        rr = rr + (cpice*ffrozp+cph2o*(1.-ffrozp))                      &
     &       *prcp/rch
      end if

! ----------------------------------------------------------------------
! include the latent heat effects of frzng rain converting to ice on
! impact in the calculation of flx2 and fnet.
! ----------------------------------------------------------------------
!      fnet = fdown - sigma * t24- ssoil
      fnet = fdown -  emissi*sigma * t24- ssoil
      if (frzgra) then
         flx2 = - lsubf * prcp
         fnet = fnet - flx2
! ----------------------------------------------------------------------
! finish penman equation calculations.
! ----------------------------------------------------------------------
      end if
      rad = fnet / rch + th2- sfctmp
       a = elcp * (q2sat - q2)
!     a = elcp1 * (q2sat - q2)
      epsca = (a * rr + rad * delta) / (delta + rr)
       etp = epsca * rch / lsubc
!     etp = epsca * rch / lvs

! ----------------------------------------------------------------------
      end subroutine penman

      end module noahmpdrv

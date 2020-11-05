!>  \file sfc_sice.f
!!  This file contains the GFS three level thermodynamic sea ice model.

!> This module contains the CCPP-compliant GFS sea ice scheme.
      module sfc_sice

      contains

      subroutine sfc_sice_init()
      end subroutine sfc_sice_init
!
      subroutine sfc_sice_finalize()
      end subroutine sfc_sice_finalize

!>\defgroup gfs_sice_main GFS Three-layer Thermodynomics Sea-Ice Scheme Module
!!  \brief  This is three-layer thermodynomics sea-ice model based on Winton (2000) \cite winton_2000.
!! \section arg_table_sfc_sice_run Argument Table
!! \htmlinclude sfc_sice_run.html
!!
!>  \section general_sice_run GFS Sea Ice Driver General Algorithm
!!The model has four prognostic variables: the snow layer thickness \f$h_s\f$, the ice layer thickness
!! \f$h_i\f$, the upper and lower ice layer temperatures located at the midpoints of the layers
!! \f$h_i/4\f$ and \f$3h_i/4\f$ below the ice surface, respectively \f$T_1\f$ and \f$T_2\f$. The temperature of
!! the bottom of the ice is fixed at \f$T_f\f$, the freezing temperature of seawater. The temperature of
!! the top of the ice or snow, \f$T_s\f$, is determined from the surface energy balance.
!! The model consists of a zero-heat-capacity snow layer overlying two equally thick sea ice layers (Figure 1).
!! The upper ice layer has a variable heat capacity to represent brine pockets.
!! \image html GFS_sice_wonton2000_fig1.png "Fig.1  Schematic representation of the three-layer model" width=5cm
!!  The ice model main program ice3lay() performs two functions:
!!  - \b Calculation \b of \b ice \b temperature
!!\n The surface temperature is determined from the diagnostic balance between
!! the upward conduction of heat through snow and/or ice and upward flux of heat
!! from the surface.
!!  - \b Calculation \b of \b ice \b and \b snow \b changes
!!\n In addition to calculating ice temperature changes, the ice model must
!! also readjust the sizes of the snow and ice layers 1) to accommodate
!! mass fluxes at the upper and lower surfaces, 2) to convert snow below
!! the water line to ice, and 3) to equalize the thickness of the two
!! ice layers.
!>  \section detailed_sice_run GFS Sea Ice Driver Detailed Algorithm
!>  @{
      subroutine sfc_sice_run                                           &
     &     ( im, kice, sbc, hvap, tgice, cp, eps, epsm1, rvrdm1, grav,  & !  ---  inputs:
     &       t0c, rd, ps, t1, q1, delt,                                 &
     &       sfcemis, dlwflx, sfcnsw, sfcdsw, srflag,                   &
     &       cm, ch, prsl1, prslki, prsik1, prslk1, wind,               &
     &       flag_iter, lprnt, ipr,                                     &
     &       hice, fice, tice, weasd, tskin, tprcp, tiice, ep,          & !  ---  input/outputs:
     &       snwdph, qsurf, snowmt, gflux, cmm, chh, evap, hflx,        & !  
     &       frac_grid, icy, islmsk_cice,                               &
     &       min_lakeice, min_seaice, oceanfrac,                        &
     &       errmsg, errflg
     &     )

! ===================================================================== !
!  description:                                                         !
!                                                                       !
!  usage:                                                               !
!                                                                       !
!    call sfc_sice                                                      !
!       inputs:                                                         !
!          ( im, kice, ps, t1, q1, delt,                                !
!            sfcemis, dlwflx, sfcnsw, sfcdsw, srflag,                   !
!            cm, ch, prsl1, prslki, prsik1, prslk1, wind,               !
!            flag_iter,                                                 !
!       input/outputs:                                                  !
!            hice, fice, tice, weasd, tskin, tprcp, tiice, ep,          !
!       outputs:                                                        !
!            snwdph, qsurf, snowmt, gflux, cmm, chh, evap, hflx )       !
!                                                                       !
!  subprogram called:  ice3lay.                                         !
!                                                                       !
!>  program history log:                                                 
!!-         2005  --  xingren wu created  from original progtm and added  
!!                     two-layer ice model                               
!!-         200x  -- sarah lu    added flag_iter           
!!-    oct  2006  -- h. wei      added cmm and chh to output     
!!-         2007  -- x. wu modified for mom4 coupling (i.e. cpldice)
!!                                    (not used anymore)
!!-         2007  -- s. moorthi micellaneous changes   
!!-    may  2009  -- y.-t. hou   modified to include surface emissivity  
!!                     effect on lw radiation. replaced the confusing  
!!                     slrad with sfc net sw sfcnsw (dn-up). reformatted
!!                     the code and add program documentation block. 
!!-    sep  2009 -- s. moorthi removed rcl, changed pressure units and 
!!                     further optimized    
!!-    jan  2015 -- x. wu change "cimin = 0.15" for both  
!!                              uncoupled and coupled case 
!                                                                       !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                       size   !
!     im, kice - integer, horiz dimension and num of ice layers    1    !
!     ps       - real, surface pressure                            im   !
!     t1       - real, surface layer mean temperature ( k )        im   !
!     q1       - real, surface layer mean specific humidity        im   !
!     delt     - real, time interval (second)                      1    !
!     sfcemis  - real, sfc lw emissivity ( fraction )              im   !
!     dlwflx   - real, total sky sfc downward lw flux ( w/m**2 )   im   !
!     sfcnsw   - real, total sky sfc netsw flx into ground(w/m**2) im   !
!     sfcdsw   - real, total sky sfc downward sw flux ( w/m**2 )   im   !
!     srflag   - real, snow/rain fraction for precipitation        im   !
!     cm       - real, surface exchange coeff for momentum (m/s)   im   !
!     ch       - real, surface exchange coeff heat & moisture(m/s) im   !
!     prsl1    - real, surface layer mean pressure                 im   !
!     prslki   - real,                                             im   !
!     prsik1   - real,                                             im   !
!     prslk1   - real,                                             im   !
!     islimsk  - integer, sea/land/ice mask (=0/1/2)               im   !
!     wind     - real,                                             im   !
!     flag_iter- logical,                                          im   !
!                                                                       !
!  input/outputs:                                                       !
!     hice     - real, sea-ice thickness                           im   !
!     fice     - real, sea-ice concentration                       im   !
!     tice     - real, sea-ice surface temperature                 im   !
!     weasd    - real, water equivalent accumulated snow depth (mm)im   !
!     tskin    - real, ground surface skin temperature ( k )       im   !
!     tprcp    - real, total precipitation                         im   !
!     tiice    - real, temperature of ice internal (k)          im,kice !
!     ep       - real, potential evaporation                       im   !
!                                                                       !
!  outputs:                                                             !
!     snwdph   - real, water equivalent snow depth (mm)            im   !
!     qsurf    - real, specific humidity at sfc                    im   !
!     snowmt   - real, snow melt (m)                               im   !
!     gflux    - real, soil heat flux (w/m**2)                     im   !
!     cmm      - real, surface exchange coeff for momentum (m/s)   im   !
!     chh      - real, surface exchange coeff heat&moisture (m/s)  im   !
!     evap     - real, evaperation from latent heat flux           im   !
!     hflx     - real, sensible heat flux                          im   !
!                                                                       !
! ===================================================================== !
!
      use machine, only : kind_phys
      use funcphys, only : fpvs
!
      implicit none
!
! - Define constant parameters
      integer,              parameter :: kmi   = 2                  !< 2-layer of ice
      real(kind=kind_phys), parameter :: zero  = 0.0_kind_phys, one = 1.0_kind_phys
      real(kind=kind_phys), parameter :: himax = 8.0_kind_phys      !< maximum ice thickness allowed
      real(kind=kind_phys), parameter :: himin = 0.1_kind_phys      !< minimum ice thickness required
      real(kind=kind_phys), parameter :: hsmax = 2.0_kind_phys      !< maximum snow depth allowed
      real(kind=kind_phys), parameter :: timin = 173.0_kind_phys    !< minimum temperature allowed for snow/ice
      real(kind=kind_phys), parameter :: albfw = 0.06_kind_phys     !< albedo for lead
      real(kind=kind_phys), parameter :: dsi   = one/0.33_kind_phys
      real(kind=kind_phys), parameter :: qmin  = 1.0e-8_kind_phys

!  ---  inputs:
      integer, intent(in) :: im, kice, ipr
      logical, intent(in) :: lprnt
      logical, intent(in) :: frac_grid

      real (kind=kind_phys), intent(in) :: sbc, hvap, tgice, cp, eps,   &
     &       epsm1, grav, rvrdm1, t0c, rd

      real (kind=kind_phys), dimension(im), intent(in) :: ps,           &
     &       t1, q1, sfcemis, dlwflx, sfcnsw, sfcdsw, srflag, cm, ch,   &
     &       prsl1, prslki, prsik1, prslk1, wind, oceanfrac

!     integer, dimension(im), intent(in) :: islimsk
      integer, dimension(im), intent(in) :: islmsk_cice
      real (kind=kind_phys), intent(in)  :: delt, min_seaice,           &
     &                                            min_lakeice

      logical, dimension(im), intent(in) :: flag_iter, icy

!  ---  input/outputs:
      real (kind=kind_phys), dimension(im), intent(inout) :: hice,      &
     &       fice, tice, weasd, tskin, tprcp, ep

      real (kind=kind_phys), dimension(im,kice), intent(inout) :: tiice

!  ---  outputs:
      real (kind=kind_phys), dimension(im), intent(inout) :: snwdph,    &
     &       qsurf, snowmt, gflux, cmm, chh, evap, hflx

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

!  ---  locals:
      real (kind=kind_phys), dimension(im) :: ffw, evapi, evapw,        &
     &       sneti, hfd, hfi,                                           &
!    &       hflxi, hflxw, sneti, snetw, qssi, qssw, hfd, hfi, hfw,     &
     &       focn, snof,                                   rch, rho,    &
     &       snowd, theta1

      real (kind=kind_phys) :: t12, t14, tem, stsice(im,kice)
     &,                        hflxi, hflxw, q0, qs1, qssi, qssw
      real (kind=kind_phys) :: cpinv, hvapi, elocp, snetw, cimin

      integer :: i, k
      integer, dimension(im) :: islmsk_local

      logical :: flag(im)
!
!===> ...  begin here
!
      cpinv = one/cp
      hvapi = one/hvap
      elocp = hvap/cp

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0


      islmsk_local = islmsk_cice
      if (frac_grid) then
        do i=1,im
          if (icy(i) .and. islmsk_local(i) < 2) then
            if (oceanfrac(i) > zero) then
              tem = min_seaice
            else
              tem = min_lakeice
            endif
            if (fice(i) > tem) then
              islmsk_local(i) = 2
              tice(i) =min( tice(i), tgice)
            endif
          endif
        enddo
      endif

!
!> - Set flag for sea-ice.

      do i = 1, im
        flag(i) = (islmsk_local(i) == 2) .and. flag_iter(i)
        if (flag_iter(i) .and. islmsk_local(i) < 2) then
          hice(i) = zero
          fice(i) = zero
        endif
      enddo

      do i = 1, im
        if (flag(i)) then
          if (srflag(i) > zero) then
            ep(i)    = ep(i)*(one-srflag(i))
            weasd(i) = weasd(i) + 1000.0_kind_phys*tprcp(i)*srflag(i)
            tprcp(i) = tprcp(i)*(one-srflag(i))
          endif
        endif
      enddo
!  --- ...  update sea ice temperature

      do k = 1, kice
        do i = 1, im
          if (flag(i)) then
            stsice(i,k) = tiice(i,k)
          endif
        enddo
      enddo

!  --- ...  initialize variables. all units are supposedly m.k.s. unless specifie
!           psurf is in pascals, wind is wind speed, theta1 is adiabatic surface
!           temp from level 1, rho is density, qs1 is sat. hum. at level1 and qss
!           is sat. hum. at surface
!           convert slrad to the civilized unit from langley minute-1 k-4

      do i = 1, im
        if (flag(i)) then
          if (oceanfrac(i) > zero) then
            cimin = min_seaice
          else
            cimin = min_lakeice
          endif
!         psurf(i) = 1000.0 * ps(i)
!         ps1(i)   = 1000.0 * prsl1(i)

!         dlwflx has been given a negative sign for downward longwave
!         sfcnsw is the net shortwave flux (direction: dn-up)

          q0        = max(q1(i), qmin)
!         tsurf(i)  = tskin(i)
#ifdef GSD_SURFACE_FLUXES_BUGFIX
          theta1(i) = t1(i) / prslk1(i) ! potential temperature in middle of lowest atm. layer
#else
          theta1(i) = t1(i) * prslki(i)
#endif
          rho(i)    = prsl1(i) / (rd*t1(i)*(one+rvrdm1*q0))
          qs1       = fpvs(t1(i))
          qs1       = max(eps*qs1 / (prsl1(i) + epsm1*qs1), qmin)
          q0        = min(qs1, q0)

          if (fice(i) < cimin) then
            print *,'warning: ice fraction is low:', fice(i)
            fice(i) = cimin
            tice(i) = tgice
            tskin(i)= tgice
            print *,'fix ice fraction: reset it to:', fice(i)
          endif
          ffw(i)    = one - fice(i)

          qssi = fpvs(tice(i))
          qssi = eps*qssi / (ps(i) + epsm1*qssi)
          qssw = fpvs(tgice)
          qssw = eps*qssw / (ps(i) + epsm1*qssw)

!> - Convert snow depth in water equivalent from mm to m unit.

          snowd(i) = weasd(i) * 0.001_kind_phys
!         flagsnw(i) = .false.

!  --- ...  when snow depth is less than 1 mm, a patchy snow is assumed and
!           soil is allowed to interact with the atmosphere.
!           we should eventually move to a linear combination of soil and
!           snow under the condition of patchy snow.

!  --- ...  rcp = rho cp ch v

          cmm(i) = cm(i)  * wind(i)
          chh(i) = rho(i) * ch(i) * wind(i)
          rch(i) = chh(i) * cp

!> - Calculate sensible and latent heat flux over open water & sea ice.

          evapi(i) = elocp * rch(i) * (qssi - q0)
          evapw(i) = elocp * rch(i) * (qssw - q0)
!         evap(i)  = fice(i)*evapi(i) + ffw(i)*evapw(i)

          snetw    = sfcdsw(i) * (one - albfw)
          snetw    = min(3.0_kind_phys*sfcnsw(i)                        &
     &                  / (one+2.0_kind_phys*ffw(i)), snetw)
!> - Calculate net solar incoming at top \a sneti.
          sneti(i) = (sfcnsw(i) - ffw(i)*snetw) / fice(i)

          t12 = tice(i) * tice(i)
          t14 = t12 * t12

!> - Calculate net non-solar and upir heat flux @ ice surface \a hfi.

#ifdef GSD_SURFACE_FLUXES_BUGFIX
          hfi(i) = -dlwflx(i) + sfcemis(i)*sbc*t14 + evapi(i)           &
     &           + rch(i)*(tice(i)/prsik1(i) - theta1(i))
#else
          hfi(i) = -dlwflx(i) + sfcemis(i)*sbc*t14 + evapi(i)           &
     &           + rch(i)*(tice(i) - theta1(i))
#endif
!> - Calculate heat flux derivative at surface \a hfd.
          hfd(i) = 4.0_kind_phys*sfcemis(i)*sbc*tice(i)*t12             &
     &           + (one + elocp*eps*hvap*qs1/(rd*t12)) * rch(i)

          t12 = tgice * tgice
          t14 = t12 * t12

!  --- ...  hfw = net heat flux @ water surface (within ice)

!         hfw(i) = -dlwflx(i) + sfcemis(i)*sbc*t14 + evapw(i)           &
!    &           + rch(i)*(tgice - theta1(i)) - snetw

!> - Assigin heat flux from ocean \a focn and snowfall rate as constants, which
!! should be from ocean model and other physics.
          focn(i) = 2.0_kind_phys   ! heat flux from ocean - should be from ocn model
          snof(i) = zero    ! snowfall rate - snow accumulates in gbphys

!> - Initialize snow depth \a snowd.
          hice(i) = max( min( hice(i), himax ), himin )
          snowd(i) = min( snowd(i), hsmax )

          if (snowd(i) > (2.0_kind_phys*hice(i))) then
            print *, 'warning: too much snow :',snowd(i)
            snowd(i) = hice(i) + hice(i)
            print *,'fix: decrease snow depth to:',snowd(i)
          endif
        endif
      enddo

!> - Call the three-layer thermodynamics sea ice model ice3lay().
      call ice3lay
!  ---  inputs:                                                         !
     &     ( im, kice, fice, flag, hfi, hfd, sneti, focn, delt,          !
     &       lprnt, ipr,
!  ---  outputs:                                                        !
     &       snowd, hice, stsice, tice, snof, snowmt, gflux )           !

      do i = 1, im
        if (flag(i)) then
          if (tice(i) < timin) then
            print *,'warning: snow/ice temperature is too low:',tice(i)
     &,' i=',i
            tice(i) = timin
            print *,'fix snow/ice temperature: reset it to:',tice(i)
          endif

          if (stsice(i,1) < timin) then
            print *,'warning: layer 1 ice temp is too low:',stsice(i,1)
     &,' i=',i
            stsice(i,1) = timin
            print *,'fix layer 1 ice temp: reset it to:',stsice(i,1)
          endif

          if (stsice(i,2) < timin) then
            print *,'warning: layer 2 ice temp is too low:',stsice(i,2)
            stsice(i,2) = timin
            print *,'fix layer 2 ice temp: reset it to:',stsice(i,2)
          endif

          tskin(i) = tice(i)*fice(i) + tgice*ffw(i)
        endif
      enddo

      do k = 1, kice
        do i = 1, im
          if (flag(i)) then
            tiice(i,k) = min(stsice(i,k), t0c)
          endif
        enddo
      enddo

      do i = 1, im
        if (flag(i)) then
!  --- ...  calculate sensible heat flux (& evap over sea ice)

#ifdef GSD_SURFACE_FLUXES_BUGFIX
          hflxi    = rch(i) * (tice(i)/prsik1(i) - theta1(i))
          hflxw    = rch(i) * (tgice / prsik1(i) - theta1(i))
#else
          hflxi    = rch(i) * (tice(i) - theta1(i))
          hflxw    = rch(i) * (tgice - theta1(i))
#endif
          hflx(i)  = fice(i)*hflxi    + ffw(i)*hflxw
          evap(i)  = fice(i)*evapi(i) + ffw(i)*evapw(i)
!
!  --- ...  the rest of the output

          qsurf(i) = q1(i) + evap(i) / (elocp*rch(i))

!  --- ...  convert snow depth back to mm of water equivalent

          weasd(i)  = snowd(i) * 1000.0_kind_phys
          snwdph(i) = weasd(i) * dsi             ! snow depth in mm

          tem     = one / rho(i)
          hflx(i) = hflx(i) * tem * cpinv
          evap(i) = evap(i) * tem * hvapi
        endif
      enddo
!
      return
!! @}

! =================
      contains
! =================


!-----------------------------------
!> This subroutine is the entity of three-layer sea ice vertical thermodynamics
!! based on Winton(2000) \cite winton_2000 .
!!\ingroup gfs_sice_main
!\param[in] im    integer, horizontal dimension
!\param[in] kmi   integer, number of ice layers (2)
!\param[in] fice  real, sea-ice concentration
!\param[in] flag  logical, ice mask flag
!\param[in] hfi   real, net non-solar and heat flux at surface (\f$W/m^2\f$)
!\param[in] hfd   real, heat flux derivative at surface
!\param[in] sneti real, net solar incoming at top (\f$W/m^2\f$)
!\param[in] focn  real, heat flux from ocean (\f$W/m^2\f$)
!\param[in] delt  real, time step(\f$sec\f$)
!\param[in,out] snowd  real, snow depth
!\param[in,out] hice real, sea-ice thickness
!\param[in,out] stsice real, temperature at mid-point of ice levels (\f$^oC\f$)
!\param[in,out] tice real, surface temperature (\f$^oC\f$)
!\param[in,out] snof real, snowfall rate (\f$ms^{-1}\f$)
!\param[out] snowmt real, snow melt during delt (\f$m\f$)
!\param[out] gflux real, conductive heat flux (\f$W/m^2\f$)
!>\section gen_ice3lay Three-layer Thermodynamics Sea Ice Model General Algorithm
!> @{
      subroutine ice3lay
!...................................
!  ---  inputs:
     &     ( im, kmi, fice, flag, hfi, hfd, sneti, focn, delt,          &
     &       lprnt, ipr,
!  ---  input/outputs:
     &       snowd, hice, stsice, tice, snof,                           &
!  ---  outputs:
     &       snowmt, gflux                                              &
     &     )

!**************************************************************************
!                                                                         *
!            three-layer sea ice vertical thermodynamics                  *
!                                                                         *
! based on:  m. winton, "a reformulated three-layer sea ice model",       *
! journal of atmospheric and oceanic technology, 2000                     *
!                                                                         *
!                                                                         *
!        -> +---------+ <- tice - diagnostic surface temperature ( <= 0c )*
!       /   |         |                                                   *
!   snowd   |  snow   | <- 0-heat capacity snow layer                     *
!       \   |         |                                                   *
!        => +---------+                                                   *
!       /   |         |                                                   *
!      /    |         | <- t1 - upper 1/2 ice temperature; this layer has *
!     /     |         |         a variable (t/s dependent) heat capacity  *
!   hice    |...ice...|                                                   *
!     \     |         |                                                   *
!      \    |         | <- t2 - lower 1/2 ice temp. (fixed heat capacity) *
!       \   |         |                                                   *
!        -> +---------+ <- base of ice fixed at seawater freezing temp.   *
!                                                                         *
!  =====================  defination of variables  =====================  !
!                                                                         !
!  inputs:                                                         size   !
!     im, kmi  - integer, horiz dimension and num of ice layers      1    !
!     fice     - real, sea-ice concentration                         im   !
!     flag     - logical, ice mask flag                              1    !
!     hfi      - real, net non-solar and heat flux @ surface(w/m^2)  im   !
!     hfd      - real, heat flux derivatice @ sfc (w/m^2/deg-c)      im   !
!     sneti    - real, net solar incoming at top  (w/m^2)            im   !
!     focn     - real, heat flux from ocean    (w/m^2)               im   !
!     delt     - real, timestep                (sec)                 1    !
!                                                                         !
!  input/outputs:                                                         !
!     snowd    - real, surface pressure                              im   !
!     hice     - real, sea-ice thickness                             im   !
!     stsice   - real, temp @ midpt of ice levels  (deg c)          im,kmi!
!     tice     - real, surface temperature     (deg c)               im   !
!     snof     - real, snowfall rate           (m/sec)               im   !
!                                                                         !
!  outputs:                                                               !
!     snowmt   - real, snow melt during delt   (m)                   im   !
!     gflux    - real, conductive heat flux    (w/m^2)               im   !
!                                                                         !
!  locals:                                                                !
!     hdi      - real, ice-water interface     (m)                        !
!     hsni     - real, snow-ice                (m)                        !
!                                                                         !
! ======================================================================= !
!

!  ---  constant parameters: (properties of ice, snow, and seawater)
      real (kind=kind_phys), parameter :: ds   = 330.0_kind_phys    !< snow (ov sea ice) density (kg/m^3)
      real (kind=kind_phys), parameter :: dw   =1000.0_kind_phys    !< fresh water density  (kg/m^3)
      real (kind=kind_phys), parameter :: dsdw = ds/dw
      real (kind=kind_phys), parameter :: dwds = dw/ds
      real (kind=kind_phys), parameter :: ks   = 0.31_kind_phys     !< conductivity of snow   (w/mk)
      real (kind=kind_phys), parameter :: i0   = 0.3_kind_phys      !< ice surface penetrating solar fraction
      real (kind=kind_phys), parameter :: ki   = 2.03_kind_phys     !< conductivity of ice  (w/mk)
      real (kind=kind_phys), parameter :: di   = 917.0_kind_phys    !< density of ice   (kg/m^3)
      real (kind=kind_phys), parameter :: didw = di/dw
      real (kind=kind_phys), parameter :: dsdi = ds/di
      real (kind=kind_phys), parameter :: ci   = 2054.0_kind_phys   !< heat capacity of fresh ice (j/kg/k)
      real (kind=kind_phys), parameter :: li   = 3.34e5_kind_phys   !< latent heat of fusion (j/kg-ice)
      real (kind=kind_phys), parameter :: si   = 1.0_kind_phys      !< salinity of sea ice
      real (kind=kind_phys), parameter :: mu   = 0.054_kind_phys    !< relates freezing temp to salinity
      real (kind=kind_phys), parameter :: tfi  = -mu*si             !< sea ice freezing temp = -mu*salinity
      real (kind=kind_phys), parameter :: tfw  = -1.8_kind_phys     !< tfw - seawater freezing temp (c)
      real (kind=kind_phys), parameter :: tfi0 = tfi-0.0001_kind_phys
      real (kind=kind_phys), parameter :: dici = di*ci
      real (kind=kind_phys), parameter :: dili = di*li
      real (kind=kind_phys), parameter :: dsli = ds*li
      real (kind=kind_phys), parameter :: ki4  = ki*4.0_kind_phys
      real (kind=kind_phys), parameter :: zero = 0.0_kind_phys, one  = 1.0_kind_phys

!  ---  inputs:
      integer, intent(in) :: im, kmi, ipr
      logical             :: lprnt

      real (kind=kind_phys), dimension(im), intent(in) :: fice, hfi,    &
     &       hfd, sneti, focn

      real (kind=kind_phys), intent(in) :: delt

      logical, dimension(im), intent(in) :: flag

!  ---  input/outputs:
      real (kind=kind_phys), dimension(im), intent(inout) :: snowd,     &
     &       hice, tice, snof

      real (kind=kind_phys), dimension(im,kmi), intent(inout) :: stsice

!  ---  outputs:
      real (kind=kind_phys), dimension(im), intent(out) :: snowmt,      &
     &       gflux

!  ---  locals:

      real (kind=kind_phys) :: dt2, dt4, dt6, h1, h2, dh, wrk, wrk1,    &
     &                         dt2i, hdi, hsni, ai, bi, a1, b1, a10, b10&
     &,                        c1, ip, k12, k32, tsf, f1, tmelt, bmelt

      integer :: i
!
!===> ...  begin here
!
      dt2  =  delt + delt
      dt4  =  dt2  + dt2
      dt6  =  dt2  + dt4
      dt2i = one / dt2

      do i = 1, im
        if (flag(i)) then
          snowd(i) = snowd(i) * dwds
          hdi      = (dsdw*snowd(i) + didw*hice(i))

          if (hice(i) < hdi) then
            snowd(i) = snowd(i) + hice(i) - hdi
            hsni     = (hdi - hice(i)) * dsdi
            hice (i) = hice(i) + hsni
          endif

          snof(i)     = snof(i) * dwds
          tice(i)     = tice(i) - t0c                  ! convert from K to C
          stsice(i,1) = min(stsice(i,1)-t0c, tfi0)     ! degc
          stsice(i,2) = min(stsice(i,2)-t0c, tfi0)     ! degc

          ip = i0 * sneti(i)         ! ip +v (in winton ip=-i0*sneti as sol -v)
          if (snowd(i) > zero) then
            tsf = zero
            ip  = zero
          else
            tsf = tfi
            ip  = i0 * sneti(i)      ! ip +v here (in winton ip=-i0*sneti)
          endif
          tice(i) = min(tice(i), tsf)

!> - Ice temperature calculation.

          bi   = hfd(i)
          ai   = hfi(i) - sneti(i) + ip - tice(i)*bi  ! +v sol input here
!>  - Calculate the effective conductive coupling of the snow-ice layer
!! between the surface and the upper layer ice temperature \f$h_i/4\f$
!! beneath the snow-ice interface (see \a eq.(5) in Winton (2000) \cite winton_2000).
          k12  = ki4*ks / (ks*hice(i) + ki4*snowd(i))

!>  - Calculate the conductive coupling between the two ice temperature
!! points (see \a eq.(10) in Winton (2000) \cite winton_2000).
          k32  = (ki+ki) / hice(i)

          wrk    = one / (dt6*k32 + dici*hice(i))
          a10    = dici*hice(i)*dt2i + k32*(dt4*k32 + dici*hice(i))*wrk
          b10    = -di*hice(i) * (ci*stsice(i,1) + li*tfi/stsice(i,1))  &
     &           * dt2i - ip                                            &
     &           - k32*(dt4*k32*tfw + dici*hice(i)*stsice(i,2)) * wrk

          wrk1  = k12 / (k12 + bi)
          a1    = a10 + bi * wrk1
          b1    = b10 + ai * wrk1
          c1    = dili * tfi * dt2i * hice(i)

!>  - Calculate the new upper ice temperature following \a eq.(21)
!! in Winton (2000) \cite winton_2000.
          stsice(i,1) = -(sqrt(b1*b1-4.0_kind_phys*a1*c1) + b1)/(a1+a1)
          tice(i) = (k12*stsice(i,1) - ai) / (k12 + bi)

!>  - If the surface temperature is greater than the freezing temperature
!! of snow (when there is snow over) or sea ice (when there is none), the
!! surface temperature is fixed at the melting temperature of snow or sea
!! ice, respectively, and the upper ice temperature is recomputed from
!! \a eq.(21) using the coefficients given by \a eqs. (19),(20), and (18). An energy flux
!! \a eq.(22) is applied toward surface melting thereby balancing the surface
!! energy budget.
          if (tice(i) > tsf) then
            a1 = a10 + k12
            b1 = b10 - k12*tsf
            stsice(i,1) = -(sqrt(b1*b1-4.0_kind_phys*a1*c1) + b1)       &
     &                  / (a1+a1)
            tice(i) = tsf
            tmelt   = (k12*(stsice(i,1)-tsf) - (ai+bi*tsf)) * delt
          else
            tmelt    =zero
            snowd(i) = snowd(i) + snof(i)*delt
          endif
!>  - Calculate the new lower ice temperature following \a eq.(15)
!! in Winton (2000) \cite winton_2000.
          stsice(i,2) = (dt2*k32*(stsice(i,1) + tfw + tfw)              &
     &                +  dici*hice(i)*stsice(i,2)) * wrk

!>  - Calculate the energy for bottom melting (or freezing, if negative)
!! following \a eq.(23), which serves to balance the difference between
!! the oceanic heat flux to the ice bottom and the conductive flux of
!! heat upward from the bottom.
          bmelt = (focn(i) + ki4*(stsice(i,2) - tfw)/hice(i)) * delt

!> - Calculation of ice and snow mass changes.

          h1 = 0.5_kind_phys * hice(i)
          h2 = 0.5_kind_phys * hice(i)

!>  - Calculate the top layer thickness.

          if (tmelt <= snowd(i)*dsli) then
            snowmt(i) = tmelt / dsli
            snowd (i) = snowd(i) - snowmt(i)
          else
            snowmt(i) = snowd(i)
            h1 = h1 - (tmelt - snowd(i)*dsli)                           &
     &         / (di * (ci - li/stsice(i,1)) * (tfi - stsice(i,1)))
            snowd(i) = zero
          endif

!  --- ...  and bottom
!>  - When the energy for bottem melting \f$M_b\f$ is negative (i.e., freezing
!! is happening),calculate the bottom layer thickness \f$h_2\f$ and the new
!! lower layer temperature (see \a eqs.(24)-(26)).
          if (bmelt < zero) then
            dh = -bmelt / (dili + dici*(tfi - tfw))
            stsice(i,2) = (h2*stsice(i,2) + dh*tfw) / (h2 + dh)
            h2 = h2 + dh
          else
            h2 = h2 - bmelt / (dili + dici*(tfi - stsice(i,2)))
          endif

!>  - If ice remains, even up 2 layers, else, pass negative energy back in snow.
!! Calculate the new upper layer temperature (see \a eq.(38)).

          hice(i) = h1 + h2

          if (hice(i) > zero) then
            if (h1 > 0.5_kind_phys*hice(i)) then
              f1 = one - (h2+h2) / hice(i)
              stsice(i,2) = f1 * (stsice(i,1) + li*tfi/(ci*stsice(i,1)))&
     &                    + (one - f1)*stsice(i,2)

              if (stsice(i,2) > tfi) then
                hice(i) = hice(i) - h2*ci*(stsice(i,2) - tfi)/ (li*delt)
                stsice(i,2) = tfi
              endif
            else
              f1 = (h1+h1) / hice(i)
              stsice(i,1) = f1 * (stsice(i,1) + li*tfi/(ci*stsice(i,1)))&
     &                    + (one - f1)*stsice(i,2)
              stsice(i,1) = (stsice(i,1) - sqrt(stsice(i,1)*stsice(i,1) &
     &                    - 4.0_kind_phys*tfi*li/ci)) * 0.5_kind_phys
            endif

            k12      = ki4*ks / (ks*hice(i) + ki4*snowd(i))
            gflux(i) = k12 * (stsice(i,1) - tice(i))
          else
            snowd(i) = snowd(i) + (h1*(ci*(stsice(i,1) - tfi)           &
     &               - li*(one - tfi/stsice(i,1)))                      &
     &               + h2*(ci*(stsice(i,2) - tfi) - li)) / li

            hice(i)     = max(zero, snowd(i)*dsdi)
            snowd(i)    = zero
            stsice(i,1) = tfw
            stsice(i,2) = tfw
            gflux(i)    = zero
          endif   ! end if_hice_block

          gflux(i)    = fice(i) * gflux(i)
          snowmt(i)   = snowmt(i) * dsdw
          snowd(i)    = snowd(i) * dsdw
          tice(i)     = tice(i) + t0c
          stsice(i,1) = stsice(i,1) + t0c
          stsice(i,2) = stsice(i,2) + t0c
        endif   ! end if_flag_block
      enddo   ! end do_i_loop

      return
!...................................
      end subroutine ice3lay
!> @}
!-----------------------------------

! =========================== !
!     end contain programs    !
! =========================== !

!...................................
      end subroutine sfc_sice_run
!-----------------------------------
!> @}
      end module sfc_sice

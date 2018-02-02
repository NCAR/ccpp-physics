!>  \file sfc_sice.f
!!  This file contains the GFS three level thermodynamic sea ice model.

!> \defgroup GFS_Ice GFS Three-layer Thermodynamics Sea Ice
!!  \brief  This is three-layer thermodynomics sea-ice model based on Winton (2000) \cite winton_2000.
!!
!! Sea ice is a thin skin of frozen water covering the polar oceans.  The sea ice strongly interacts with both the atmosphere above and the ocean underneath in the high
!! latitudes. In a coupled weather/climate system, changes in sea ice extent, thickness and concentration
!! regionally or globally would influence oceanic and atmospheric conditions, which in turn affect the 
!! sea ice distribution. The physical and dynamical processes affecting the weather and climate are 
!! considered as follows.
!!
!! The high albedo of the sea ice reflects more solar radiation back to the space. The feedbacks are
!! considered as positive. The broader the sea ice cover, the higher the surface albedo, which result
!! in less amount of solar radiation absorbed at the Earth's surface. A cooler surface would favor more 
!! sea ice to form. The process would be reversed in less sea ice situation.
!!
!! The sea ice restricts the heat/water exchange between the air and ocean. The presence of extensive
!! areas of sea ice would suppress the heat loss in winter and the heat gain in summer by the ocean.
!! Even a thin ice cover influences the turbulent heat transfer significantly between ocean and
!! atmosphere. The surface fluxes of sensible and latent heat can be greater by up to two orders of magnitude
!! at the open water surface of a lead or polynya than that through (snow covered) pack ice.
!!
!! The sea ice modifies air/sea momentum transfer, ocean fresh water balance and ocean circulation.
!! The freezing and melting of the ocean surface and the associated fluxes of salt and heat produce major
!! changes in the density structure of the polar water. Formation of sea ice injects salt into the ocean
!! makes the water heavier and more convectively unstable, conversely when melting occurs, stable and fresh
!! layers can prevent deep covective activity.
!!
!! A sea ice model, in general, may contain subcomponents treating 1)dynamics (ice motion),
!! 2)ice transport, 3) multiple ice thickness categories (including leads), 4) surface albedo,
!! and 5) vertical thermodynamics. GFS sea ice scheme is concerned with a scheme for the 
!! last of these processes.
!!
!! A three-layer thermodynamic sea ice model has been coupled to NCEP GFS. It predicts sea ice/snow thickness,
!! the surface temperature and ice temperature structure. In each model grid box, the heat and moisture
!! fluxes and albedo are treated separately for the ice and the open water.
!!
!!\image html GFS_sice_wonton2000_fig1.png "Fig.1  Schematic representation of the three-layer model" width=5cm
!! The model has four prognostic variables: the snow layer thickness \f$h_s\f$, the ice layer thickness 
!! \f$h_i\f$, the upper and lower ice layer temperatures located at the midpoints of the layers
!! \f$h_i/4\f$ and \f$3h_i/4\f$ below the ice surface, respectively \f$T_1\f$ and \f$T_2\f$. The temperature of 
!! the bottom of the ice is fixed at \f$T_f\f$, the freezing temperature of seawater. The temperature of 
!! the top of the ice or snow, \f$T_s\f$, is determined from the surface energy balance.
!!
!! The model consists of a zero-heat-capacity snow layer overlying two equally thick sea ice layers (Fig.1).
!! The upper ice layer has a variable heat capacity to represent brine pockets. 
!!
!!  \section intraphysics Intraphysics Communication
!!\image html schematic_sice.png "Fig.2  NCEP Sea Ice Model System Diagram" width=10cm

      module sfc_sice

      contains

! \brief This subroutine is empty since there are no procedures that need to be done to initialize the GFS SICE code.
! \section arg_table_sice_init  Argument Table
!
      subroutine sfc_sice_init
      end subroutine sfc_sice_init
!!
! \brief This subroutine is empty since there are no procedures that need to be done to finalize the GFS SICE code.
! \section arg_table_sice_finalize  Argument Table
!
      subroutine sfc_sice_finalize
      end subroutine sfc_sice_finalize


!>\defgroup gfs_sice_main GFS sfc_sice Main
!! @{
!! \ingroup GFS_Ice
!! \section arg_table_sfc_sice_run Arguments
!! | local var name | longname                                                                     | description                                                     | units         | rank | type    |    kind   | intent | optional |
!! |----------------|------------------------------------------------------------------------------|-----------------------------------------------------------------|---------------|------|---------|-----------|--------|----------|
!! | im             | horizontal_loop_extent                                                       | horizontal loop extent                                          | count         |    0 | integer |           | in     | F        |
!! | km             | soil_vertical_dimension                                                      | vertical loop extent for soil levels, start at 1                | count         |    0 | integer |           | in     | F        |
!! | ps             | surface_air_pressure                                                         | surface pressure                                                | Pa            |    1 | real    | kind_phys | in     | F        |
!! | u1             | x_wind_at_lowest_model_layer                                                 | u component of surface layer wind                               | m s-1         |    1 | real    | kind_phys | in     | F        |
!! | v1             | y_wind_at_lowest_model_layer                                                 | v component of surface layer wind                               | m s-1         |    1 | real    | kind_phys | in     | F        |
!! | t1             | air_temperature_at_lowest_model_layer                                        | surface layer mean temperature                                  | K             |    1 | real    | kind_phys | in     | F        |
!! | q1             | specific_humidity_at_lowest_model_layer                                      | surface layer mean specific humidity                            | kg kg-1       |    1 | real    | kind_phys | in     | F        |
!! | delt           | time_step_for_dynamics                                                       | time step                                                       | s             |    0 | real    | kind_phys | in     | F        |
!! | sfcemis        | surface_longwave_emissivity                                                  | sfc lw emissivity                                               | frac          |    1 | real    | kind_phys | in     | F        |
!! | dlwflx         | surface_downwelling_longwave_flux_absorbed_by_ground                         | total sky surface downward longwave flux absorbed by the ground | W m-2         |    1 | real    | kind_phys | in     | F        |
!! | sfcnsw         | surface_net_downwelling_shortwave_flux                                       | total sky sfc netsw flx into ground                             | W m-2         |    1 | real    | kind_phys | in     | F        |
!! | sfcdsw         | surface_downwelling_shortwave_flux                                           | total sky sfc downward sw flux                                  | W m-2         |    1 | real    | kind_phys | in     | F        |
!! | srflag         | flag_for_precipitation_type                                                  | snow/rain flag for precipitation                                | flag          |    1 | real    | kind_phys | in     | F        |
!! | cm             | surface_drag_coefficient_for_momentum_in_air                                 | surface exchange coeff for momentum                             | none          |    1 | real    | kind_phys | in     | F        |
!! | ch             | surface_drag_coefficient_for_heat_and_moisture_in_air                        | surface exchange coeff heat & moisture                          | none          |    1 | real    | kind_phys | in     | F        |
!! | prsl1          | air_pressure_at_lowest_model_layer                                           | surface layer mean pressure                                     | Pa            |    1 | real    | kind_phys | in     | F        |
!! | prslki         | ratio_of_exner_function_between_midlayer_and_interface_at_lowest_model_layer | Exner function ratio bt midlayer and interface at 1st layer     | ratio         |    1 | real    | kind_phys | in     | F        |
!! | islimsk        | sea_land_ice_mask                                                            | sea/land/ice mask (=0/1/2)                                      | flag          |    1 | integer |           | in     | F        |
!! | ddvel          | surface_wind_enhancement_due_to_convection                                   | wind enhancement due to convection                              | m s-1         |    1 | real    | kind_phys | in     | F        |
!! | flag_iter      | flag_for_iteration                                                           | flag for iteration                                              | flag          |    1 | logical |           | in     | F        |
!! | mom4ice        | flag_for_mom4_coupling                                                       | flag for Mom4 coupling                                          | flag          |    0 | logical |           | in     | F        |
!! | lsm            | flag_for_land_surface_scheme                                                 | flag for land sfc scheme =0: osu; =1: noah                      | flag          |    0 | integer |           | in     | F        |
!! | lprnt          | flag_print                                                                   | switch for printing sample column to stdout                     | flag          |    0 | logical |           | in     | F        |
!! | ipr            | horizontal_index_of_printed_column                                           | horizontal index of printed column                              | index         |    0 | integer |           | in     | F        |
!! | hice           | sea_ice_thickness                                                            | sea-ice thickness                                               | m             |    1 | real    | kind_phys | inout  | F        |
!! | fice           | sea_ice_concentration                                                        | sea-ice concentration [0,1]                                     | frac          |    1 | real    | kind_phys | inout  | F        |
!! | tice           | sea_ice_temperature                                                          | sea-ice surface temperature                                     | K             |    1 | real    | kind_phys | inout  | F        |
!! | weasd          | water_equivalent_accumulated_snow_depth                                      | water equivalent accumulated snow depth                         | mm            |    1 | real    | kind_phys | inout  | F        |
!! | tskin          | surface_skin_temperature                                                     | ground surface skin temperature                                 | K             |    1 | real    | kind_phys | inout  | F        |
!! | tprcp          | nonnegative_lwe_thickness_of_precipitation_amount_on_dynamics_timestep       | nonnegative precipitation amount in one dynamics time step      | m             |    1 | real    | kind_phys | inout  | F        |
!! | stc            | soil_temperature                                                             | soil temp                                                       | K             |    2 | real    | kind_phys | inout  | F        |
!! | ep             | surface_upward_potential_latent_heat_flux                                    | potential evaporation                                           | W m-2         |    1 | real    | kind_phys | inout  | F        |
!! | snwdph         | surface_snow_thickness_water_equivalent                                      | water equivalent snow depth                                     | mm            |    1 | real    | kind_phys |   out  | F        |
!! | qsurf          | surface_specific_humidity                                                    | sfc air saturation specific humidity                            | kg kg-1       |    1 | real    | kind_phys |   out  | F        |
!! | snowmt         | surface_snow_melt                                                            | snow melt during timestep                                       | m             |    1 | real    | kind_phys |   out  | F        |
!! | gflux          | upward_heat_flux_in_soil                                                     | soil heat flux                                                  | W m-2         |    1 | real    | kind_phys |   out  | F        |
!! | cmm            | surface_drag_wind_speed_for_momentum_in_air                                  | surf mom exch coef time mean surf wind                          | m s-1         |    1 | real    | kind_phys |   out  | F        |
!! | chh            | surface_drag_mass_flux_for_heat_and_moisture_in_air                          | surf h&m exch coef time surf wind & density                     | kg m-2 s-1    |    1 | real    | kind_phys |   out  | F        |
!! | evap           | kinematic_surface_upward_latent_heat_flux                                    | evaporative latent heat flux                                    | kg kg-1 m s-1 |    1 | real    | kind_phys |   out  | F        |
!! | hflx           | kinematic_surface_upward_sensible_heat_flux                                  | kinematic sensible heat flux                                    | K m s-1       |    1 | real    | kind_phys |   out  | F        |
!!
!!  \section general General Algorithm
!!  \section detailed Detailed Algorithm
!>  @{
      subroutine sfc_sice_run                                           &
     &     ( im, km, ps, u1, v1, t1, q1, delt,                          &
     &       sfcemis, dlwflx, sfcnsw, sfcdsw, srflag,                   &
     &       cm, ch, prsl1, prslki, islimsk, ddvel,                     &
     &       flag_iter, mom4ice, lsm, lprnt, ipr,                       & ! -- inputs from here and above
     &       hice, fice, tice, weasd, tskin, tprcp, stc, ep,            & ! -- in/outs
     &       snwdph, qsurf, snowmt, gflux, cmm, chh, evap, hflx         & ! -- outputs
     &     )

! ===================================================================== !
!  description:                                                         !
!                                                                       !
!  usage:                                                               !
!                                                                       !
!    call sfc_sice                                                      !
!       inputs:                                                         !
!          ( im, km, ps, u1, v1, t1, q1, delt,                          !
!            sfcemis, dlwflx, sfcnsw, sfcdsw, srflag,                   !
!            cm, ch, prsl1, prslki, islimsk, ddvel,                     !
!            flag_iter, mom4ice, lsm,                                   !
!       input/outputs:                                                  !
!            hice, fice, tice, weasd, tskin, tprcp, stc, ep,            !
!       outputs:                                                        !
!            snwdph, qsurf, snowmt, gflux, cmm, chh, evap, hflx )       !
!                                                                       !
!  subprogram called:  ice3lay.                                         !
!                                                                       !
u  program history log:                                                 !
!         2005  --  xingren wu created  from original progtm and added  !
!                     two-layer ice model                               !
!         200x  -- sarah lu    added flag_iter                          !
!    oct  2006  -- h. wei      added cmm and chh to output              !
!         2007  -- x. wu modified for mom4 coupling (i.e. mom4ice)      !
!         2007  -- s. moorthi micellaneous changes                      !
!    may  2009  -- y.-t. hou   modified to include surface emissivity   !
!                     effect on lw radiation. replaced the confusing    !
!                     slrad with sfc net sw sfcnsw (dn-up). reformatted !
!                     the code and add program documentation block.     !
!    sep  2009 -- s. moorthi removed rcl, changed pressure units and    !
!                     further optimized                                 !
!    jan  2015 -- x. wu change "cimin = 0.15" for both                  !
!                              uncoupled and coupled case               !
!                                                                       !
!                                                                       !
!  ====================  defination of variables  ====================  !
!                                                                       !
!  inputs:                                                       size   !
!     im, km   - integer, horiz dimension and num of soil layers   1    !
!     ps       - real, surface pressure                            im   !
!     u1, v1   - real, u/v component of surface layer wind         im   !
!     t1       - real, surface layer mean temperature ( k )        im   !
!     q1       - real, surface layer mean specific humidity        im   !
!     delt     - real, time interval (second)                      1    !
!     sfcemis  - real, sfc lw emissivity ( fraction )              im   !
!     dlwflx   - real, total sky sfc downward lw flux ( w/m**2 )   im   !
!     sfcnsw   - real, total sky sfc netsw flx into ground(w/m**2) im   !
!     sfcdsw   - real, total sky sfc downward sw flux ( w/m**2 )   im   !
!     srflag   - real, snow/rain flag for precipitation            im   !
!     cm       - real, surface exchange coeff for momentum (m/s)   im   !
!     ch       - real, surface exchange coeff heat & moisture(m/s) im   !
!     prsl1    - real, surface layer mean pressure                 im   !
!     prslki   - real,                                             im   !
!     islimsk  - integer, sea/land/ice mask (=0/1/2)               im   !
!     ddvel    - real,                                             im   !
!     flag_iter- logical,                                          im   !
!     mom4ice  - logical,                                          im   !
!     lsm      - integer, flag for land surface model scheme       1    !
!                =0: use osu scheme; =1: use noah scheme                !
!                                                                       !
!  input/outputs:                                                       !
!     hice     - real, sea-ice thickness                           im   !
!     fice     - real, sea-ice concentration                       im   !
!     tice     - real, sea-ice surface temperature                 im   !
!     weasd    - real, water equivalent accumulated snow depth (mm)im   !
!     tskin    - real, ground surface skin temperature ( k )       im   !
!     tprcp    - real, total precipitation                         im   !
!     stc      - real, soil temp (k)                              im,km !
!     ep       - real, potential evaporation                       im   !
!                                                                       !
!  outputs:                                                             !
!     snwdph   - real, water equivalent snow depth (mm)            im   !
!     qsurf    - real, specific humidity at sfc                    im   !
!     snowmt   - real, snow melt (m)                               im   !
!     gflux    - real, soil heat flux (w/m**2)                     im   !
!     cmm      - real,                                             im   !
!     chh      - real,                                             im   !
!     evap     - real, evaperation from latent heat flux           im   !
!     hflx     - real, sensible heat flux                          im   !
!                                                                       !
! ===================================================================== !
!
      use machine , only : kind_phys
      use funcphys, only : fpvs
      use physcons, only : sbc => con_sbc, hvap => con_hvap,            &
     &                     tgice => con_tice, cp => con_cp,             &
     &                     eps => con_eps, epsm1 => con_epsm1,          &
     &                     grav => con_g, rvrdm1 => con_fvirt,          &
     &                     t0c => con_t0c, rd => con_rd
!
      implicit none
!
!  ---  constant parameters:
      integer,              parameter :: kmi   = 2        ! 2-layer of ice
      real(kind=kind_phys), parameter :: cpinv = 1.0/cp
      real(kind=kind_phys), parameter :: hvapi = 1.0/hvap
      real(kind=kind_phys), parameter :: elocp = hvap/cp
      real(kind=kind_phys), parameter :: himax = 8.0      ! maximum ice thickness allowed
      real(kind=kind_phys), parameter :: himin = 0.1      ! minimum ice thickness required
      real(kind=kind_phys), parameter :: hsmax = 2.0      ! maximum snow depth allowed
      real(kind=kind_phys), parameter :: timin = 173.0    ! minimum temperature allowed for snow/ice
      real(kind=kind_phys), parameter :: albfw = 0.06     ! albedo for lead
      real(kind=kind_phys), parameter :: dsi   = 1.0/0.33

!  ---  inputs:
      integer, intent(in) :: im, km, lsm, ipr
      logical, intent(in) :: lprnt

      real (kind=kind_phys), dimension(im), intent(in) :: ps, u1, v1,   &
     &       t1, q1, sfcemis, dlwflx, sfcnsw, sfcdsw, srflag, cm, ch,   &
     &       prsl1, prslki, ddvel

      integer, dimension(im), intent(in) :: islimsk
      real (kind=kind_phys), intent(in)  :: delt

      logical, intent(in) :: flag_iter(im), mom4ice

!  ---  input/outputs:
      real (kind=kind_phys), dimension(im), intent(inout) :: hice,      &
     &       fice, tice, weasd, tskin, tprcp, ep

      real (kind=kind_phys), dimension(im,km), intent(inout) :: stc

!  ---  outputs:
      real (kind=kind_phys), dimension(im), intent(out) :: snwdph,      &
     &       qsurf, snowmt, gflux, cmm, chh, evap, hflx

!  ---  locals:
      real (kind=kind_phys), dimension(im) :: ffw, evapi, evapw,        &
     &       sneti, snetw, hfd, hfi,                                    &
!    &       hflxi, hflxw, sneti, snetw, qssi, qssw, hfd, hfi, hfw,     &
     &       focn, snof, hi_save, hs_save,                 rch, rho,    &
     &       snowd, theta1

      real (kind=kind_phys) :: t12, t14, tem, stsice(im,kmi)
     &,                        hflxi, hflxw, q0, qs1, wind, qssi, qssw
      real (kind=kind_phys), parameter :: cimin=0.15 !  --- minimum ice concentration

      integer :: i, k

      logical :: flag(im)
!
!===> ...  begin here
!
!  --- ...  set flag for sea-ice

      do i = 1, im
        flag(i) = (islimsk(i) >= 2) .and. flag_iter(i)
        if (flag_iter(i) .and. islimsk(i) < 2) then
          hice(i) = 0.0
          fice(i) = 0.0
        endif
      enddo

!  --- ...  update sea ice temperature

      do k = 1, kmi
        do i = 1, im
          if (flag(i)) then
            stsice(i,k) = stc(i,k)
          endif
        enddo
      enddo
!
      if (mom4ice) then
        do i = 1, im
          if (flag(i)) then
            hi_save(i) = hice(i)
            hs_save(i) = weasd(i) * 0.001
          endif
        enddo
      elseif (lsm > 0) then           !  --- ...  snow-rain detection
        do i = 1, im
          if (flag(i)) then
            if (srflag(i) == 1.0) then
              ep(i) = 0.0
              weasd(i) = weasd(i) + 1.e3*tprcp(i)
              tprcp(i)  = 0.0
            endif
          endif
        enddo
      endif

!  --- ...  initialize variables. all units are supposedly m.k.s. unless specifie
!           psurf is in pascals, wind is wind speed, theta1 is adiabatic surface
!           temp from level 1, rho is density, qs1 is sat. hum. at level1 and qss
!           is sat. hum. at surface
!           convert slrad to the civilized unit from langley minute-1 k-4

      do i = 1, im
        if (flag(i)) then
!         psurf(i) = 1000.0 * ps(i)
!         ps1(i)   = 1000.0 * prsl1(i)

!         dlwflx has been given a negative sign for downward longwave
!         sfcnsw is the net shortwave flux (direction: dn-up)

          wind      = max(sqrt(u1(i)*u1(i) + v1(i)*v1(i))               &
     &                  + max(0.0, min(ddvel(i), 30.0)), 1.0)

          q0        = max(q1(i), 1.0e-8)
!         tsurf(i)  = tskin(i)
          theta1(i) = t1(i) * prslki(i)
          rho(i)    = prsl1(i) / (rd*t1(i)*(1.0+rvrdm1*q0))
          qs1       = fpvs(t1(i))
          qs1       = max(eps*qs1 / (prsl1(i) + epsm1*qs1), 1.e-8)
          q0        = min(qs1, q0)

          ffw(i)    = 1.0 - fice(i)
          if (fice(i) < cimin) then
            print *,'warning: ice fraction is low:', fice(i)
            fice(i) = cimin
            ffw (i) = 1.0 - fice(i)
            tice(i) = tgice
            tskin(i)= tgice
            print *,'fix ice fraction: reset it to:', fice(i)
          endif

          qssi = fpvs(tice(i))
          qssi = eps*qssi / (ps(i) + epsm1*qssi)
          qssw = fpvs(tgice)
          qssw = eps*qssw / (ps(i) + epsm1*qssw)

!  --- ...  snow depth in water equivalent is converted from mm to m unit

          if (mom4ice) then
            snowd(i) = weasd(i) * 0.001 / fice(i)
          else
            snowd(i) = weasd(i) * 0.001
          endif
!         flagsnw(i) = .false.

!  --- ...  when snow depth is less than 1 mm, a patchy snow is assumed and
!           soil is allowed to interact with the atmosphere.
!           we should eventually move to a linear combination of soil and
!           snow under the condition of patchy snow.

!  --- ...  rcp = rho cp ch v

          cmm(i) = cm(i)  * wind
          chh(i) = rho(i) * ch(i) * wind
          rch(i) = chh(i) * cp

!  --- ...  sensible and latent heat flux over open water & sea ice

          evapi(i) = elocp * rch(i) * (qssi - q0)
          evapw(i) = elocp * rch(i) * (qssw - q0)
!         evap(i)  = fice(i)*evapi(i) + ffw(i)*evapw(i)

!     if (lprnt) write(0,*)' tice=',tice(ipr)

          snetw(i) = sfcdsw(i) * (1.0 - albfw)
          snetw(i) = min(3.0*sfcnsw(i)/(1.0+2.0*ffw(i)), snetw(i))
          sneti(i) = (sfcnsw(i) - ffw(i)*snetw(i)) / fice(i)

          t12 = tice(i) * tice(i)
          t14 = t12 * t12

!  --- ...  hfi = net non-solar and upir heat flux @ ice surface

          hfi(i) = -dlwflx(i) + sfcemis(i)*sbc*t14 + evapi(i)           &
     &           + rch(i)*(tice(i) - theta1(i))
          hfd(i) = 4.0*sfcemis(i)*sbc*tice(i)*t12                       &
     &           + (1.0 + elocp*eps*hvap*qs1/(rd*t12)) * rch(i)

          t12 = tgice * tgice
          t14 = t12 * t12

!  --- ...  hfw = net heat flux @ water surface (within ice)

!         hfw(i) = -dlwflx(i) + sfcemis(i)*sbc*t14 + evapw(i)           &
!    &           + rch(i)*(tgice - theta1(i)) - snetw(i)

          focn(i) = 2.0     ! heat flux from ocean - should be from ocn model
          snof(i) = 0.0     ! snowfall rate - snow accumulates in gbphys

          hice(i) = max( min( hice(i), himax ), himin )
          snowd(i) = min( snowd(i), hsmax )

          if (snowd(i) > (2.0*hice(i))) then
            print *, 'warning: too much snow :',snowd(i)
            snowd(i) = hice(i) + hice(i)
            print *,'fix: decrease snow depth to:',snowd(i)
          endif
        endif
      enddo

!     if (lprnt) write(0,*)' tice2=',tice(ipr)
      call ice3lay
!  ---  inputs:                                                         !
!    &     ( im, kmi, fice, flag, hfi, hfd, sneti, focn, delt,          !
!  ---  outputs:                                                        !
!    &       snowd, hice, stsice, tice, snof, snowmt, gflux )           !

!     if (lprnt) write(0,*)' tice3=',tice(ipr)
      if (mom4ice) then
        do i = 1, im
          if (flag(i)) then
            hice(i)  = hi_save(i)
            snowd(i) = hs_save(i)
          endif
        enddo
      endif

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

      do k = 1, kmi
        do i = 1, im
          if (flag(i)) then
            stc(i,k) = min(stsice(i,k), t0c)
          endif
        enddo
      enddo

      do i = 1, im
        if (flag(i)) then
!  --- ...  calculate sensible heat flux (& evap over sea ice)

          hflxi    = rch(i) * (tice(i) - theta1(i))
          hflxw    = rch(i) * (tgice - theta1(i))
          hflx(i)  = fice(i)*hflxi    + ffw(i)*hflxw
          evap(i)  = fice(i)*evapi(i) + ffw(i)*evapw(i)
!
!  --- ...  the rest of the output

          qsurf(i) = q1(i) + evap(i) / (elocp*rch(i))

!  --- ...  convert snow depth back to mm of water equivalent

          weasd(i)  = snowd(i) * 1000.0
          snwdph(i) = weasd(i) * dsi             ! snow depth in mm

          tem     = 1.0 / rho(i)
          hflx(i) = hflx(i) * tem * cpinv
          evap(i) = evap(i) * tem * hvapi
        endif
      enddo
!
      return

! =================
      contains
! =================

!> @}

!-----------------------------------
!> This subroutine is the entity of three-layer sea ice vertical thermodynamics 
!! based on Winton(2000) \cite winton_2000 .
      subroutine ice3lay
!...................................
!  ---  inputs:
!    &     ( im, kmi, fice, flag, hfi, hfd, sneti, focn, delt,          &
!  ---  input/outputs:
!    &       snowd, hice, stsice, tice, snof,                           &
!  ---  outputs:
!    &       snowmt, gflux                                              &
!    &     )

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
      real (kind=kind_phys), parameter :: ds   = 330.0    ! snow (ov sea ice) density (kg/m^3)
      real (kind=kind_phys), parameter :: dw   =1000.0    ! fresh water density  (kg/m^3)
      real (kind=kind_phys), parameter :: dsdw = ds/dw
      real (kind=kind_phys), parameter :: dwds = dw/ds
      real (kind=kind_phys), parameter :: t0c  =273.15    ! freezing temp of fresh ice (k)
      real (kind=kind_phys), parameter :: ks   = 0.31     ! conductivity of snow   (w/mk)
      real (kind=kind_phys), parameter :: i0   = 0.3      ! ice surface penetrating solar fraction
      real (kind=kind_phys), parameter :: ki   = 2.03     ! conductivity of ice  (w/mk)
      real (kind=kind_phys), parameter :: di   = 917.0    ! density of ice   (kg/m^3)
      real (kind=kind_phys), parameter :: didw = di/dw
      real (kind=kind_phys), parameter :: dsdi = ds/di
      real (kind=kind_phys), parameter :: ci   = 2054.0   ! heat capacity of fresh ice (j/kg/k)
      real (kind=kind_phys), parameter :: li   = 3.34e5   ! latent heat of fusion (j/kg-ice)
      real (kind=kind_phys), parameter :: si   = 1.0      ! salinity of sea ice
      real (kind=kind_phys), parameter :: mu   = 0.054    ! relates freezing temp to salinity
      real (kind=kind_phys), parameter :: tfi  = -mu*si   ! sea ice freezing temp = -mu*salinity
      real (kind=kind_phys), parameter :: tfw  = -1.8     ! tfw - seawater freezing temp (c)
      real (kind=kind_phys), parameter :: tfi0 = tfi-0.0001
      real (kind=kind_phys), parameter :: dici = di*ci
      real (kind=kind_phys), parameter :: dili = di*li
      real (kind=kind_phys), parameter :: dsli = ds*li
      real (kind=kind_phys), parameter :: ki4  = ki*4.0

!  ---  inputs:
!     integer, intent(in) :: im, kmi

!     real (kind=kind_phys), dimension(im), intent(in) :: fice, hfi,    &
!    &       hfd, sneti, focn

!     real (kind=kind_phys), intent(in) :: delt

!     logical, dimension(im), intent(in) :: flag

!  ---  input/outputs:
!     real (kind=kind_phys), dimension(im), intent(inout) :: snowd,     &
!    &       hice, tice, snof

!     real (kind=kind_phys), dimension(im,kmi), intent(inout) :: stsice

!  ---  outputs:
!     real (kind=kind_phys), dimension(im), intent(out) :: snowmt,      &
!    &       gflux

!  ---  locals:

      real (kind=kind_phys) :: dt2, dt4, dt6, h1, h2, dh, wrk, wrk1,    &
     &                         dt2i, hdi, hsni, ai, bi, a1, b1, a10, b10&
     &,                        c1, ip, k12, k32, tsf, f1, tmelt, bmelt

      integer :: i
!
!===> ...  begin here
!
      dt2  = 2.0 * delt
      dt4  = 4.0 * delt
      dt6  = 6.0 * delt
      dt2i = 1.0 / dt2

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
          if (snowd(i) > 0.0) then
            tsf = 0.0
            ip  = 0.0
          else
            tsf = tfi
            ip  = i0 * sneti(i)      ! ip +v here (in winton ip=-i0*sneti)
          endif
          tice(i) = min(tice(i), tsf)

!  --- ...  compute ice temperature

          bi   = hfd(i)
          ai   = hfi(i) - sneti(i) + ip - tice(i)*bi  ! +v sol input here
          k12  = ki4*ks / (ks*hice(i) + ki4*snowd(i))
          k32  = (ki+ki) / hice(i)

          wrk    = 1.0 / (dt6*k32 + dici*hice(i))
          a10    = dici*hice(i)*dt2i + k32*(dt4*k32 + dici*hice(i))*wrk
          b10    = -di*hice(i) * (ci*stsice(i,1) + li*tfi/stsice(i,1))  &
     &           * dt2i - ip                                            &
     &           - k32*(dt4*k32*tfw + dici*hice(i)*stsice(i,2)) * wrk

          wrk1  = k12 / (k12 + bi)
          a1    = a10 + bi * wrk1
          b1    = b10 + ai * wrk1
          c1    = dili * tfi * dt2i * hice(i)

          stsice(i,1) = -(sqrt(b1*b1 - 4.0*a1*c1) + b1)/(a1+a1)
          tice(i) = (k12*stsice(i,1) - ai) / (k12 + bi)

          if (tice(i) > tsf) then
            a1 = a10 + k12
            b1 = b10 - k12*tsf
            stsice(i,1) = -(sqrt(b1*b1 - 4.0*a1*c1) + b1)/(a1+a1)
            tice(i) = tsf
            tmelt   = (k12*(stsice(i,1)-tsf) - (ai+bi*tsf)) * delt
          else
            tmelt    = 0.0
            snowd(i) = snowd(i) + snof(i)*delt
          endif

          stsice(i,2) = (dt2*k32*(stsice(i,1) + tfw + tfw)              &
     &                +  dici*hice(i)*stsice(i,2)) * wrk

          bmelt = (focn(i) + ki4*(stsice(i,2) - tfw)/hice(i)) * delt

!  --- ...  resize the ice ...

          h1 = 0.5 * hice(i)
          h2 = 0.5 * hice(i)

!  --- ...  top ...

          if (tmelt <= snowd(i)*dsli) then
            snowmt(i) = tmelt / dsli
            snowd (i) = snowd(i) - snowmt(i)
          else
            snowmt(i) = snowd(i)
            h1 = h1 - (tmelt - snowd(i)*dsli)                           &
     &         / (di * (ci - li/stsice(i,1)) * (tfi - stsice(i,1)))
            snowd(i) = 0.0
          endif

!  --- ...  and bottom

          if (bmelt < 0.0) then
            dh = -bmelt / (dili + dici*(tfi - tfw))
            stsice(i,2) = (h2*stsice(i,2) + dh*tfw) / (h2 + dh)
            h2 = h2 + dh
          else
            h2 = h2 - bmelt / (dili + dici*(tfi - stsice(i,2)))
          endif

!  --- ...  if ice remains, even up 2 layers, else, pass negative energy back in snow

          hice(i) = h1 + h2

          if (hice(i) > 0.0) then
            if (h1 > 0.5*hice(i)) then
              f1 = 1.0 - (h2+h2) / hice(i)
              stsice(i,2) = f1 * (stsice(i,1) + li*tfi/(ci*stsice(i,1)))&
     &                    + (1.0 - f1)*stsice(i,2)

              if (stsice(i,2) > tfi) then
                hice(i) = hice(i) - h2*ci*(stsice(i,2) - tfi)/ (li*delt)
                stsice(i,2) = tfi
              endif
            else
              f1 = (h1+h1) / hice(i)
              stsice(i,1) = f1 * (stsice(i,1) + li*tfi/(ci*stsice(i,1)))&
     &                    + (1.0 - f1)*stsice(i,2)
              stsice(i,1) = (stsice(i,1) - sqrt(stsice(i,1)*stsice(i,1) &
     &                    - 4.0*tfi*li/ci)) * 0.5
            endif

            k12      = ki4*ks / (ks*hice(i) + ki4*snowd(i))
            gflux(i) = k12 * (stsice(i,1) - tice(i))
          else
            snowd(i) = snowd(i) + (h1*(ci*(stsice(i,1) - tfi)           &
     &               - li*(1.0 - tfi/stsice(i,1)))                      &
     &               + h2*(ci*(stsice(i,2) - tfi) - li)) / li

            hice(i)     = max(0.0, snowd(i)*dsdi)
            snowd(i)    = 0.0
            stsice(i,1) = tfw
            stsice(i,2) = tfw
            gflux(i)    = 0.0
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
!-----------------------------------

! =========================== !
!     end contain programs    !
! =========================== !

!...................................
      end subroutine sfc_sice_run
!-----------------------------------
      end module sfc_sice


      module sfc_sice_pre

      contains

!!
! \brief This subroutine is empty since there are no procedures needed
! \section arg_table_sfc_sice_pre_init  Argument Table
!!
      subroutine sfc_sice_pre_init
      end subroutine sfc_sice_pre_init

!!
! \brief This subroutine is empty since there are no procedures needed
! \section arg_table_sfc_sice_pre_finalize  Argument Table
!!
      subroutine sfc_sice_pre_finalize
      end subroutine sfc_sice_pre_finalize


!! \section arg_table_sfc_sice_pre_run Argument Table
!! | local var name | longname                                                                     | description                                                 | units         | rank | type    |    kind   | intent | optional |
!! |----------------|------------------------------------------------------------------------------|-------------------------------------------------------------|---------------|------|---------|-----------|--------|----------|
!! | im             | horizontal_loop_extent                                                       | horizontal loop extent                                      | count         |    0 | integer |           | in     | F        |
!! | fice           | sea_ice_concentration                                                        | sea-ice concentration [0,1]                                 | frac          |    1 | real    | kind_phys | in     | F        |
!! | hice           | sea_ice_thickness                                                            | sea-ice thickness                                           | m             |    1 | real    | kind_phys | in     | F        |
!! | tisfc          | sea_ice_temperature                                                          | sea-ice surface temperature                                 | K             |    1 | real    | kind_phys | in     | F        |
!! | prsik          | exner_function_at_lowest_model_interface                                     | Exner function at lowest model interface                    | none          |    1 | real    | kind_phys | in     | F        |
!! | prslk          | dimensionless_exner_function_at_lowest_model_layer                           | dimensionless Exner function at lowest model layer          |  none         |    1 | real    | kind_phys | in     | F        |
!! | cice           | sea_ice_concentration                                                        | sea-ice concentration [0,1]                                 | frac          |    1 | real    | kind_phys |   out  | F        |
!! | zice           | sea_ice_thickness                                                            | sea-ice thickness                                           | m             |    1 | real    | kind_phys |   out  | F        |
!! | tice           | sea_ice_temperature                                                          | sea-ice surface temperature                                 | K             |    1 | real    | kind_phys |   out  | F        |
!! | work3          | ratio_of_exner_function_between_midlayer_and_interface_at_lowest_model_layer | Exner function ratio bt midlayer and interface at 1st layer | ratio         |    1 | real    | kind_phys |   out  | F        |
!!
!! @{
      subroutine sfc_sice_pre_run(im, fice, hice, tisfc , prsik, prslk,       &
     &                            cice, zice, tice, work3)

      use machine, only : kind_phys

      implicit none

! --- inputs
      integer :: im
      real(kind=kind_phys), dimension(im), intent(in) :: fice, hice,    &
     &     tisfc, prsik, prslk

! --- input/output
      real(kind=kind_phys), dimension(im), intent(out) :: cice, zice,   &
     &     tice, work3

! --- locals
      integer :: i

      do i = 1, im
! transfer ice thickness & concentration from global to local variables
        zice(i) = hice(i)
        cice(i) = fice(i)
        tice(i) = tisfc(i)
        work3(i)= prsik(i) / prslk(i)
      enddo

      return
      end subroutine sfc_sice_pre_run

      end module sfc_sice_pre

      module sfc_sice_post

      contains

!!
! \brief This subroutine is empty since there are no procedures needed
! \section arg_table_sfc_sice_post_init  Argument Table
!!
      subroutine sfc_sice_post_init
      end subroutine sfc_sice_post_init

!!
! \brief This subroutine is empty since there are no procedures needed
! \section arg_table_sfc_sice_post_finalize  Argument Table
!!
      subroutine sfc_sice_post_finalize
      end subroutine sfc_sice_post_finalize


!!
!! \section arg_table_sfc_sice_post_run Argument Table
!! | local var name | longname                                              | description                                 | units         | rank | type    |    kind   | intent | optional |
!! |----------------|-------------------------------------------------------|---------------------------------------------|---------------|------|---------|-----------|--------|----------|
!! | im             | horizontal_loop_extent                                | horizontal loop extent                      | count         |    0 | integer |           | in     | F        |
!! | islmsk         | sea_land_ice_mask                                     | sea/land/ice mask (=0/1/2)                  | flag          |    1 | integer |           | in     | F        |
!! | cice           | sea_ice_concentration                                 | sea-ice concentration [0,1]                 | frac          |    1 | real    | kind_phys | in     | F        |
!! | zice           | sea_ice_thickness                                     | sea-ice thickness                           | m             |    1 | real    | kind_phys | in     | F        |
!! | tice           | sea_ice_temperature                                   | sea-ice surface temperature                 | K             |    1 | real    | kind_phys | in     | F        |
!! | tsfc           | surface_skin_temperature                              | surface skin temperature                    | K             |    1 | real    | kind_phys | in     | F        |
!! | fice           | sea_ice_concentration                                 | sea-ice concentration [0,1]                 | frac          |    1 | real    | kind_phys |   out  | F        |
!! | hice           | sea_ice_thickness                                     | sea-ice thickness                           | m             |    1 | real    | kind_phys |   out  | F        |
!! | tisfc          | sea_ice_temperature                                   | sea-ice surface temperature                 | K             |    1 | real    | kind_phys |   out  | F        |
!!
!! @{
      subroutine sfc_sice_post_run(im, islmsk, cice, zice, tice, tsfc,        &
     &                             fice, hice, tisfc)

      use machine, only : kind_phys

      implicit none

! --- input
      integer :: im
      integer, dimension(im) :: islmsk
      real(kind=kind_phys), dimension(im), intent(in) :: cice, zice,    &
     &     tice, tsfc

! --- outputs
      real(kind=kind_phys), dimension(im), intent(out) :: fice, hice,   &
     &     tisfc

! --- locals
      integer :: i

!--- return updated ice thickness & concentration to global arrays
!    where there is no ice, set temperature to surface skin temperature.
      do i = 1, im

        if (islmsk(i) == 2) then
           hice(i) = zice(i)
           fice(i) = cice(i)
           tisfc(i) = tice(i)
        else
           hice(i) = 0.0
           fice(i) = 0.0
           tisfc(i) = tsfc(i)
        endif
      enddo

      end subroutine sfc_sice_post_run
!! @}

      end module  sfc_sice_post
!> @}

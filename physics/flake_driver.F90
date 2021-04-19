!>  \file flake_driver.F90
!!  This file contains the flake scheme driver.

!> This module contains the CCPP-compliant flake scheme driver.
     module flake_driver

      implicit none

      private

      public :: flake_driver_init, flake_driver_run, flake_driver_finalize

      contains

!> \section arg_table_flake_driver_init Argument Table
!! \htmlinclude flake_driver_init.html
!!
      subroutine flake_driver_init (errmsg, errflg)

      implicit none
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg


      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      end subroutine flake_driver_init

!> \section arg_table_flake_driver_finalize Argument Table
!! \htmlinclude flake_driver_finalize.html
!!
      subroutine flake_driver_finalize (errmsg, errflg)

      implicit none

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      end subroutine flake_driver_finalize

!> \section arg_table_flake_driver_run Argument Table
!! \htmlinclude flake_driver_run.html
!!
      SUBROUTINE flake_driver_run (                          &
! ---- Inputs
            im, ps, t1, q1, wind,                            &
            dlwflx, dswsfc, weasd, lakedepth,                &
            use_flake, xlat, delt, zlvl, elev,               &
            wet, flag_iter, yearlen, julian, imon,           &
! ---- in/outs
            snwdph, hice, tsurf, fice, T_sfc, hflx, evap,    &
            ustar, qsfc, ch, cm, chh, cmm,                   &
            errmsg, errflg                     )     

!==============================================================================
!
! Declarations
!      use module_flake_ini, only:flake_init
      use module_FLake
!      use flake_albedo_ref
!      use data_parameters
!      use flake_derivedtypes
!      use flake_paramoptic_ref 
!      use flake_parameters
      use machine , only : kind_phys
!      use funcphys, only : fpvs
!      use physcons, only : grav   => con_g,    cp   => con_cp,          &
!     &                     hvap   => con_hvap, rd   => con_rd,          &
!     &                     eps    => con_eps, epsm1 => con_epsm1,       &
!     &                     rvrdm1 => con_fvirt

!==============================================================================
IMPLICIT NONE

      integer, intent(in) :: im, imon,yearlen
!      integer, dimension(im), intent(in) :: islmsk

      real (kind=kind_phys), dimension(im), intent(in) :: ps, wind,     &
     &           t1, q1, dlwflx, dswsfc, zlvl, elev

      real (kind=kind_phys),  intent(in) :: delt

      real (kind=kind_phys), dimension(im), intent(in) ::               &
     &           xlat, weasd, lakedepth

      real (kind=kind_phys),dimension(im),intent(inout) ::              &
     &           snwdph, hice, tsurf, t_sfc, hflx, evap, fice, ustar, qsfc,   &
     &           ch, cm, chh, cmm           

      real (kind=kind_phys),  intent(in) :: julian

      logical, dimension(im), intent(in) :: flag_iter, wet, use_flake

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

! --- locals

      real (kind=kind_phys) , parameter :: lake_pct_min = 0.1

      real (kind=kind_phys), dimension(im) :: &
      T_snow     , & ! Temperature at the air-snow interface [K]
      T_ice      , & ! Temperature at the snow-ice or air-ice interface [K]
      T_mnw      , & ! Mean temperature of the water column [K]
      T_wML      , & ! Mixed-layer temperature [K]
      T_bot      , & ! Temperature at the water-bottom sediment interface [K]
      T_B1       , & ! Temperature at the upper layer of the sediments [K]
      C_T        , & ! Shape factor (thermocline)
      fetch      , & ! Typical wind fetch [m]
      h_ML       , & ! Thickness of the mixed-layer [m]
      H_B1       , & ! Thickness of the upper layer of bottom sediments [m]
      w_albedo   , & !
      w_extinc 

!  Input (procedure arguments)

REAL (KIND = kind_phys) ::   &

  dMsnowdt_in                       , & ! The rate of snow accumulation [kg m^{-2} s^{-1}]
  I_atm_in                          , & ! Solar radiation flux at the surface [W m^{-2}]
  Q_atm_lw_in                       , & ! Long-wave radiation flux from the atmosphere [W m^{-2}]
  height_u_in                       , & ! Height above the lake surface where the wind speed is measured [m]
  height_tq_in                      , & ! Height where temperature and humidity are measured [m]
  U_a_in                            , & ! Wind speed at z=height_u_in [m s^{-1}]
  T_a_in                            , & ! Air temperature at z=height_tq_in [K]
  q_a_in                            , & ! Air specific humidity at z=height_tq_in
  P_a_in                                ! Surface air pressure [N m^{-2} = kg m^{-1} s^{-2}]

REAL (KIND = kind_phys) ::   &
  depth_w                           , & ! The lake depth [m]
  fetch_in                          , & ! Typical wind fetch [m]
  depth_bs_in                       , & ! Depth of the thermally active layer of the bottom sediments [m]
  T_bs_in                           , & ! Temperature at the outer edge of 
                                        ! the thermally active layer of the bottom sediments [K]
  par_Coriolis                      , & ! The Coriolis parameter [s^{-1}]
  del_time                              ! The model time step [s]

REAL (KIND = kind_phys) ::   &
  T_snow_in                        , & ! Temperature at the air-snow interface [K] 
  T_ice_in                         , & ! Temperature at the snow-ice or air-ice interface [K]
  T_mnw_in                         , & ! Mean temperature of the water column [K]
  T_wML_in                         , & ! Mixed-layer temperature [K]
  T_bot_in                         , & ! Temperature at the water-bottom sediment interface [K]
  T_B1_in                          , & ! Temperature at the bottom of the upper layer of the sediments [K]
  C_T_in                           , & ! Shape factor (thermocline)
  h_snow_in                        , & ! Snow thickness [m]
  h_ice_in                         , & ! Ice thickness [m]
  h_ML_in                          , & ! Thickness of the mixed-layer [m]
  H_B1_in                          , & ! Thickness of the upper layer of bottom sediments [m]
  T_sfc_in                         , & ! Surface temperature at the previous time step [K]  
  ch_in                            , &
  cm_in                            , &
  albedo_water                     , &
  water_extinc

REAL (KIND = kind_phys) ::   &
  T_snow_out                        , & ! Temperature at the air-snow interface [K] 
  T_ice_out                         , & ! Temperature at the snow-ice or air-ice interface [K]
  T_mnw_out                         , & ! Mean temperature of the water column [K]
  T_wML_out                         , & ! Mixed-layer temperature [K]
  T_bot_out                         , & ! Temperature at the water-bottom sediment interface [K]
  T_B1_out                          , & ! Temperature at the bottom of the upper layer of the sediments [K]
  C_T_out                           , & ! Shape factor (thermocline)
  h_snow_out                        , & ! Snow thickness [m]
  h_ice_out                         , & ! Ice thickness [m]
  h_ML_out                          , & ! Thickness of the mixed-layer [m]
  H_B1_out                          , & ! Thickness of the upper layer of bottom sediments [m]
  T_sfc_out                         , & ! surface temperature [K]
  T_sfc_n                           , & ! Updated surface temperature [K]  
  u_star                            , &
  q_sfc                             , &
  chh_out                           , &
  cmm_out

REAL (KIND = kind_phys) ::   &
  Q_momentum             , & ! Momentum flux [N m^{-2}]
  Q_SHT_flx              , & ! Sensible heat flux [W m^{-2}]
  Q_LHT_flx              , & ! Latent heat flux [W m^{-2}]
  Q_watvap                   ! Flux of water vapour [kg m^{-2} s^{-1}]

REAL (KIND = kind_phys) ::   &
  lake_depth_max, T_bot_2_in, T_bot_2_out, dxlat,tb,tr,tt,temp,Kbar, DelK


REAL (KIND = kind_phys) :: x, y !temperarory variables used for Tbot and Tsfc
                                !initilizations 

INTEGER :: i,ipr,iter

LOGICAL :: lflk_botsed_use
logical :: flag(im)
CHARACTER(LEN=*), PARAMETER  :: FMT2 = "(1x,8(F12.4,1x))"

!==============================================================================
!  Start calculations
!------------------------------------------------------------------------------
!       FLake_write need to assign original value to make the model somooth

       lake_depth_max = 60.0
       ipr     = min(im,10)

!  --- ...  set flag for lake points

      do i = 1, im
        flag(i) = (wet(i) .and. flag_iter(i))
      enddo

        Kbar=3.5
        DelK=3.0

      do i = 1, im
        if (flag(i)) then
          if( use_flake(i) ) then
           T_ice(i)    = 273.15
           T_snow(i)   = 273.15
           fetch(i)    = 2.0E+03
           C_T(i)      = 0.50

           dxlat        = 57.29578*abs(xlat(i))
           tt = 29.275+0.0813*dxlat-0.0052*dxlat*dxlat-0.0038*elev(i)+273.15
           tb = 29.075-0.7566*dxlat+0.0051*dxlat*dxlat-0.0038*elev(i)+273.15
!           if(fice(i).le.0.0) then
!                h_ice(i) = 0.0
!                h_snow(i)= 0.0
!           endif
           if(snwdph(i).gt.0.0 .or. hice(i).gt.0.0) then
               if(tsurf(i).lt.T_ice(i)) then
                 T_sfc(i) = T_ice(i)
               else
                 T_sfc(i) = tsurf(i)
               endif
           else
!               if(tsurf(i).lt.tt) then
!                  T_sfc(i)    = tt
!               else
!                  T_sfc(i) = tsurf(i)
!               endif
               T_sfc(i) = 0.1*tt + 0.9* tsurf(i)
           endif
!
!  Add empirical climatology of lake Tsfc and Tbot to the current Tsfc and Tbot
! to make sure Tsfc and Tbot are warmer than Tair in Winter or colder than Tair
! in Summer

           x = 0.03279*julian
           if(xlat(i) .ge. 0.0) then
              y = ((((0.0034*x-0.1241)*x+1.6231)*x-8.8666)*x+17.206)*x-4.2929
              T_sfc(i) = T_sfc(i) + 0.3*y
              tb = tb  + 0.05*y
           else
              y = ((((0.0034*x-0.1241)*x+1.6231)*x-8.8666)*x+17.206)*x-4.2929
              T_sfc(i) = T_sfc(i) - 0.3*y                                                                            
              tb = tb - 0.05*y
           endif
           T_bot(i)    = tb 
           T_B1(i)     = tb 

!           if(lakedepth(i).lt.10.0) then
!               T_bot(i) = T_sfc(i)
!               T_B1(i)  = T_bot(i)
!           endif

           T_mnw(i)    = C_T(i)*T_sfc(i)+(1-C_T(i))*T_bot(i)
           T_wML(i)    = C_T(i)*T_sfc(i)+(1-C_T(i))*T_bot(i)
           h_ML(i)     = C_T(i)* min ( lakedepth(i), lake_depth_max )
           H_B1(i)     = min ( lakedepth(i),4.0)
           hflx(i)     = 0.0
           evap(i)     = 0.0

! compute albedo as a function of julian day and latitute
           temp = 2*3.14159265*(julian-1)/float(yearlen)
           temp = 0.006918-0.399912*cos(temp)+0.070257*sin(temp)-  &
              0.006758*cos(2.0*temp)+0.000907*sin(2.0*temp) -  &
              0.002697*cos(3.0*temp)+0.00148*sin(3.0*temp)
           w_albedo(I) = 0.06/cos((xlat(i)-temp)/1.2)
!           w_albedo(I) = 0.06
! compute water extinction coefficient as a function of julian day
           if(julian.lt.90 .or. julian .gt. 333) then
             w_extinc(i) = Kbar-Kbar/DelK
           else
             w_extinc(i) = Kbar+Kbar/DelK*sin(2*3.14159265*(julian-151)/244)
           endif
!           w_extinc(i) = 3.0

!     write(65,1002) julian,xlat(i),w_albedo(I),w_extinc(i),lakedepth(i),elev(i),tb,tt,tsurf(i),T_sfc(i)
!     print 1002 julian,xlat(i),w_albedo(I),w_extinc(i),lakedepth(i),elev(i),tb,tt,tsurf(i),T_sfc(i)
!     print*,'inside flake driver'
!     print*,  julian,xlat(i),w_albedo(I),w_extinc(i),lakedepth(i),elev(i),tb,tt,tsurf(i),T_sfc(i)

        endif  !lake 
        endif  !flag
      enddo
 1001 format ( 'At icount=', i5, '  x = ', f5.2,5x, 'y = ', &
               1p, e12.3)
! 1002 format ( ' julian= ',F6.2,1x,5(F8.4,1x),3(f11.4,1x))
 1002 format (I4,1x,3(f8.4,1x),6(f11.4,1x))


!  
!  call lake interface
       do i=1,im
          if (flag(i)) then
            if( use_flake(i) ) then
              dMsnowdt_in = weasd(i)/delt
              I_atm_in    = dswsfc(i)
              Q_atm_lw_in = dlwflx(i)
              height_u_in = zlvl(i)
              height_tq_in = zlvl(i)
              U_a_in       = wind(i)
              T_a_in      = t1(i)
              q_a_in      = q1(i)
              P_a_in      = ps(i)
              ch_in       = ch(i)
              cm_in       = cm(i)
              albedo_water= w_albedo(i) 
              water_extinc= w_extinc(i) 

              depth_w     = min ( lakedepth(i), lake_depth_max )
              depth_bs_in = max ( 4.0, min ( depth_w * 0.2, 10.0 ) )
              fetch_in    = fetch(i)
              T_bs_in     = T_bot(i)
              par_Coriolis =  2 * 7.2921 / 100000. * sin ( xlat(i) )
              del_time    = delt

      do iter=1,10                !interation loop
            T_snow_in   = T_snow(i)
            T_ice_in    = T_ice(i)
            T_mnw_in    = T_mnw(i)
            T_wML_in    = T_wML(i)
            T_bot_in    = T_bot(i)
            T_B1_in     = T_B1(i)
            C_T_in      = C_T(i)
            h_snow_in   = snwdph(i)
            h_ice_in    = hice(i)
            h_ML_in     = h_ML(i)
            H_B1_in     = H_B1(i)
            T_sfc_in    = T_sfc(i)

            T_bot_2_in  = T_bot(i)
            Q_SHT_flx   = hflx(i)
            Q_watvap    = evap(i)

!------------------------------------------------------------------------------
!  Set the rate of snow accumulation
!------------------------------------------------------------------------------

       CALL flake_interface(dMsnowdt_in, I_atm_in, Q_atm_lw_in, height_u_in,       &
                  height_tq_in, U_a_in, T_a_in, q_a_in, P_a_in,                    &

                  depth_w, fetch_in, depth_bs_in, T_bs_in, par_Coriolis, del_time, &
                  T_snow_in, T_ice_in, T_mnw_in, T_wML_in, T_bot_in, T_B1_in,      &
                  C_T_in, h_snow_in, h_ice_in, h_ML_in, H_B1_in, T_sfc_in,         &
                  ch_in, cm_in, albedo_water, water_extinc,                        &
!
                  T_snow_out, T_ice_out, T_mnw_out, T_wML_out, T_bot_out,          &
                  T_B1_out, C_T_out, h_snow_out, h_ice_out, h_ML_out,              &
                  H_B1_out, T_sfc_out, Q_SHT_flx, Q_watvap,                        &
!
                  T_bot_2_in, T_bot_2_out,u_star, q_sfc,chh_out,cmm_out ) 

!------------------------------------------------------------------------------
! Update output and values for previous time step
!
              T_snow(i) = T_snow_out
              T_ice(i)  = T_ice_out
              T_mnw(i)  = T_mnw_out
              T_wML(i)  = T_wML_out
              T_sfc(i)  = T_sfc_out
              Tsurf(i)  = T_sfc_out
              T_bot(i)  = T_bot_out
              T_B1(i)   = T_B1_out
              C_T(i)    = C_T_out
              h_ML(i)   = h_ML_out
              H_B1(i)   = H_B1_out
              ustar(i)  = u_star
              qsfc(i)   = q_sfc
              chh(i)    = chh_out
              cmm(i)    = cmm_out
              snwdph(i) = h_snow_out
              hice(i)  = h_ice_out
              evap(i)   = Q_watvap
              hflx(i)   = Q_SHT_flx

              if(hice(i) .gt. 0.0 .or. snwdph(i) .gt. 0.0) then
                fice(i) = 1.0
              else
                fice(i) = 0.0
              endif
          enddo   !iter loop
         endif    !endif of lake
       endif      !endif of flag

       ENDDO

 125   format(1x,i2,1x,i2,1x,i2,1x,6(1x,f14.8))
 126   format(1x,i2,1x,i2,1x,6(1x,f14.8))
 127   format(1x,i2,2(1x,f16.9))
!------------------------------------------------------------------------------
!  End calculations
!==============================================================================

END SUBROUTINE flake_driver_run

!---------------------------------
      end module flake_driver

!>  \file flake_driver.F90
!!  This file contains the flake scheme driver.

!> This module contains the CCPP-compliant flake scheme driver.
     module flake_driver

      implicit none

      private

      public :: flake_driver_run

      contains

!> \section arg_table_flake_driver_run Argument Table
!! \htmlinclude flake_driver_run.html
!!
      SUBROUTINE flake_driver_run (                          &
! ---- Inputs
            im, ps, t1, q1, wind, min_lakeice,               &
            dlwflx, dswsfc, lakedepth,                       &
            use_lake_model, snow, xlat, delt, zlvl, elev,    &
            wet, yearlen, julian, imon,                      &
            flag_iter, first_time_step, flag_restart,        &
            weasd,                                           &
! ---- in/outs
            snwdph, hice, tsurf, t_sfc, fice, hflx, evap,    &
            lflx, gflx, ustar, qsfc, ch, cm, chh, cmm,       &
            h_ML, t_wML, t_mnw, H_B, T_B, t_bot1,            &
            t_bot2, c_t, T_snow, T_ice, tsurf_ice,           &
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

      implicit none
      integer, intent(in) :: im, imon,yearlen
!     integer, dimension(im), intent(in) :: islmsk

      real (kind=kind_phys), dimension(:), intent(in) :: ps, wind,     &
     &           t1, q1, dlwflx, dswsfc, zlvl, elev

      real (kind=kind_phys),  intent(in) :: delt, min_lakeice

      real (kind=kind_phys), dimension(:), intent(in) ::               &
     &           xlat, lakedepth, snow

      real (kind=kind_phys), dimension(:), intent(in) :: weasd

      real (kind=kind_phys),dimension(:),intent(inout) ::                     &
     &           snwdph, hice, tsurf, t_sfc, hflx, evap, fice, ustar, qsfc,    &
     &           ch, cm, chh, cmm, h_ML, t_wML, t_mnw, H_B, T_B,              &
     &           t_bot1, t_bot2, c_t, T_snow, T_ice, tsurf_ice, lflx, gflx

      real (kind=kind_phys),  intent(in) :: julian

      logical, dimension(:), intent(in) :: flag_iter, wet
      integer, dimension(:), intent(in) :: use_lake_model
      logical,               intent(in) :: flag_restart, first_time_step

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

! --- locals
      real (kind=kind_phys), parameter :: lake_pct_min = 0.1

      real (kind=kind_phys), dimension(im) :: &
!      T_snow     , & ! Temperature at the air-snow interface [K]
!      T_ice      , & ! Temperature at the snow-ice or air-ice interface [K]
!      T_mnw      , & ! Mean temperature of the water column [K]
!      T_wML      , & ! Mixed-layer temperature [K]
!      T_bot      , & ! Temperature at the water-bottom sediment interface [K]
!      T_B       , & ! Temperature at the upper layer of the sediments [K]
!      C_T        , & ! Shape factor (thermocline)
      fetch      , & ! Typical wind fetch [m]
!      h_ML       , & ! Thickness of the mixed-layer [m]
!      H_B1       , & ! Thickness of the upper layer of bottom sediments [m]
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
        T_B_in                          , & ! Temperature at the bottom of the upper layer of the sediments [K]
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
        T_B_out                          , & ! Temperature at the bottom of the upper layer of the sediments [K]
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

      REAL (KIND = kind_phys) :: &
        Q_momentum             , & ! Momentum flux [N m^{-2}]
        Q_SHT_flx              , & ! Sensible heat flux [W m^{-2}]
        Q_LHT_flx              , & ! Latent heat flux [W m^{-2}]
        Q_watvap               , & ! Flux of water vapour [kg m^{-2} s^{-1}]
        Q_gflx                 , & ! Flux from ice to water [W m^{-2}]
        Q_lflx                     ! latent fluxes [W m^{-2}]    

      REAL (KIND = kind_phys) ::   &
        lake_depth_max, T_bot_2_in, T_bot_2_out, dlat,tb,tr,tt,temp,temp2

      real (kind=kind_phys), parameter :: pi=4.0_kind_phys*atan(1.0_kind_phys)
      real (kind=kind_phys), parameter :: degrad=180.0_kind_phys/pi
      real (kind=kind_phys), parameter :: Kbar = 3.5_kind_phys, DelK = 3.0_kind_phys, &
                                          KbaroDelK = Kbar / DelK

      REAL (KIND = kind_phys) :: x, y, w !temperarory variables used for Tbot and Tsfc
                                      !initilizations

      INTEGER :: i,ipr,iter

      LOGICAL :: lflk_botsed_use, do_flake
      logical :: flag(im)
!     CHARACTER(LEN=*), PARAMETER  :: FMT2 = "(1x,8(F12.4,1x))"

!==============================================================================
!  Start calculations
!------------------------------------------------------------------------------
!     FLake_write need to assign original value to make the model somooth
      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

!  --- ...  set flag for lake points

      do_flake = .false.
      do i = 1, im
            flag(i) = flag_iter(i) .and. use_lake_model(i) .gt. 0
            do_flake  = flag(i) .or. do_flake
      enddo
      if (.not. do_flake) return

      lake_depth_max = 60.0
      ipr = min(im,10)

      x = 0.03279*julian
      y = ((((0.0034*x-0.1241)*x+1.6231)*x-8.8666)*x+17.206)*x-4.2929

      temp = (pi+pi)*(julian-1)/float(yearlen)
      temp = 0.006918-0.399912*cos(temp)+0.070257*sin(temp)   &
           - 0.006758*cos(2.0*temp)+0.000907*sin(2.0*temp)    &
           - 0.002697*cos(3.0*temp)+0.00148*sin(3.0*temp)

      temp2 = sin((pi+pi)*(julian-151)/244)

      do i = 1, im
        if (flag(i) .and. lakedepth(i) >1.0) then
           if(.not.flag_restart .and. first_time_step) then
              T_ice(i)    = 273.15
              T_snow(i)   = 273.15
              C_T(i)      = 0.50
              dlat       = abs(xlat(i))
              if(dlat .lt. 1.40) then
                 tt = (((21.181*dlat-51.376)*dlat+20.808)*dlat-3.8408)*dlat+29.554
                 tt = tt -0.0038*elev(i)+273.15
                 tb = (((-29.794*dlat+96.91)*dlat-86.129)*dlat-7.1921)*dlat+28.176
                 tb = tb -0.0038*elev(i)+273.15
                 w =  (((2.5467*dlat-7.4683)*dlat+5.2465)*dlat+0.4360)*dlat+0.0643
              else
                 tt = 4.0+273.15-0.0038*elev(i)
                 tb = 0.05+273.15-0.0038*elev(i)
                 w = 0.207312
              endif
              if(tsurf(i) > 400.00) then
                  write(0,*) tsurf(i)
                  write(0,*) 'Surface temperature initial is bad'
                  tsurf(i) = tt
                  write(0,*) tsurf(i)
              endif
              T_sfc(i) = 0.05*tt + 0.95* tsurf(i)

!  Add empirical climatology of lake Tsfc and Tbot to the current Tsfc and Tbot
! to make sure Tsfc and Tbot are warmer than Tair in Winter or colder than Tair
! in Summer

              if (xlat(i) >= 0.0) then
                 T_sfc(i) = T_sfc(i) + 0.05*y*w
                 tb = tb  + 0.005*y*w
              else
                 T_sfc(i) = T_sfc(i) - 0.5*y*w
                 tb = tb - 0.005*y*w
              endif

              t_bot1(i) = tb
              t_bot2(i) = tb
              T_B(i)  = tb

             T_mnw(i) = C_T(i)*T_sfc(i) + (1-C_T(i))*t_bot1(i)
             T_wML(i) = C_T(i)*T_sfc(i) + (1-C_T(i))*t_bot1(i)
             h_ML(i)  = C_T(i)* min ( lakedepth(i), lake_depth_max )
             H_B(i)  = min ( lakedepth(i),4.0)
             hflx(i)  = 0.0
             lflx(i)  = 0.0
             evap(i)  = 0.0
             chh        = ch(i) * wind(i) * 1.225 !(kg/m3)
             cmm        = cm(i) * wind(i)
           endif  !end of .not.flag_restart

          fetch(i)    = 2.0E+03
! compute albedo as a function of julian day and latitude
!          write(0,*) ' xlat= ',xlat(i), temp
          w_albedo(I) = 0.06/cos((xlat(i)-temp)/1.2)
!         w_albedo(I) = 0.06
! compute water extinction coefficient as a function of julian day
          if (julian < 90 .or. julian > 333) then
            w_extinc(i) = Kbar - KbaroDelK
          else
            w_extinc(i) = Kbar + KbaroDelK*temp2
          endif
!         w_extinc(i) = 3.0

!     write(0,1002) julian,xlat(i),w_albedo(I),w_extinc(i),elev(i),tsurf(i),T_sfc(i),t_bot1(i)
!     write(0,1003) use_lake_model(i),i,lakedepth(i), snwdph(i), hice(i), fice(i)        
!     write(0,1004) ps(i), wind(i), t1(i), q1(i), dlwflx(i), dswsfc(i), zlvl(i)

        endif  !flag
      enddo
 1002 format ( 'julian=',F6.2,1x,F8.3,1x,2(E7.2,1x),E7.2,1x,3(E7.2,1x))
 1003 format ( 'use_lake_model=',I2,1x,I3,1x,F6.4,1x,F9.4,1x,2(F8.4,1x),F7.4)
 1004 format ( 'pressure',F12.2,1x,F6.2,1x,F7.2,1x,F7.4,1x,2(F8.2,1x),F8.4)
!
!  call lake interface
       do i=1,im
         if (flag(i) .and. lakedepth(i) > 1.0) then
!         write(0,*) 'flag(i)= ', i, flag(i)
!           if(weasd(i) < 0.0 .or. hice(i) < 0.0) weasd(i) =0.0
           if(snwdph(i) < 0.0) snwdph(i) =0.0
!           dMsnowdt_in  = 10.0*0.001*weasd(i)/delt
!           dMsnowdt_in  = snow(i)/delt
           dMsnowdt_in  = snow(i)*0.001
           if(dMsnowdt_in < 0.0) dMsnowdt_in=0.0
           I_atm_in     = dswsfc(i)
           Q_atm_lw_in  = dlwflx(i)
           height_u_in  = zlvl(i)
           height_tq_in = zlvl(i)
           U_a_in       = wind(i)
           T_a_in       = t1(i)
           q_a_in       = q1(i)
           P_a_in       = ps(i)
           ch_in        = ch(i)
           cm_in        = cm(i)
           albedo_water = w_albedo(i)
           water_extinc = w_extinc(i)

           depth_w      = min ( lakedepth(i), lake_depth_max )
           depth_bs_in  = max ( 4.0, min ( depth_w * 0.2, 10.0 ) )
           fetch_in     = fetch(i)
           T_bs_in      = T_bot1(i)
           par_Coriolis =  2 * 7.2921 / 100000. * sin ( xlat(i) )
           del_time     = delt

!           if(lakedepth(i).lt.10) then
!                T_sfc(i) = t1(i)
!                T_bs_in = T_sfc(i)
!                T_B(i)  = T_bs_in
!           endif

           do iter=1,5                !interation loop
             T_snow_in   = T_snow(i)
             T_ice_in    = T_ice(i)
             T_mnw_in    = T_mnw(i)
             T_wML_in    = T_wML(i)
             T_bot_in    = t_bot1(i)
             T_B_in      = T_B(i)
             C_T_in      = C_T(i)
             h_snow_in   = snwdph(i)
             h_ice_in    = hice(i)
             h_ML_in     = h_ML(i)
             H_B1_in     = H_B(i)
             T_sfc_in    = T_sfc(i)
             tsurf_ice(i)= T_ice(i)

             T_bot_2_in  = t_bot2(i)
             Q_SHT_flx   = hflx(i)
             Q_watvap    = evap(i)
             Q_gflx      = 0.0
             Q_lflx      = 0.0

!------------------------------------------------------------------------------
!  Set the rate of snow accumulation
!------------------------------------------------------------------------------

             CALL flake_interface(dMsnowdt_in, I_atm_in, Q_atm_lw_in, height_u_in, &
                  height_tq_in, U_a_in, T_a_in, q_a_in, P_a_in,                    &

                  depth_w, fetch_in, depth_bs_in, T_bs_in, par_Coriolis, del_time, &
                  T_snow_in, T_ice_in, T_mnw_in, T_wML_in, T_bot_in, T_B_in,       &
                  C_T_in, h_snow_in, h_ice_in, h_ML_in, H_B1_in, T_sfc_in,         &
                  ch_in, cm_in, albedo_water, water_extinc,                        &
!
                  T_snow_out, T_ice_out, T_mnw_out, T_wML_out, T_bot_out,          &
                  T_B_out, C_T_out, h_snow_out, h_ice_out, h_ML_out,               &
                  H_B1_out, T_sfc_out, Q_SHT_flx, Q_watvap, Q_gflx, Q_lflx,        &
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
             tsurf_ice(i) =  T_ice(i)
             t_bot1(i)  = T_bot_out
             t_bot2(i) = T_bot_2_out
             T_B(i)   = T_B_out
             C_T(i)    = C_T_out
             h_ML(i)   = h_ML_out
             H_B(i)   = H_B1_out
             ustar(i)  = u_star
             qsfc(i)   = q_sfc
             chh(i)    = chh_out
             cmm(i)    = cmm_out
             snwdph(i) = h_snow_out
             hice(i)   = h_ice_out
             evap(i)   = Q_watvap
             hflx(i)   = Q_SHT_flx
             gflx(i)   = Q_gflx
             lflx(i)   = Q_lflx
!             if(lflx(i) > 2500.00 .or. Tsurf(i) > 350.00) then
!                 write(0,125) i,lflx(i), Tsurf(i),ps(i), wind(i),     &
!     &           t1(i), q1(i), dlwflx(i), dswsfc(i),hflx(i)
!             endif
!             fice(i)   = fice(i)+0.01*(h_ice_out-h_ice_in)
!             if(fice(i) .lt. min_lakeice ) then
!                fice(i) = 0.0
!             elseif(fice(i) .gt. 1.0) then
!                fice(i) = 1.0
!             endif
           enddo   !iter loop
!           endif    !endif use_lake_model

         endif     !endif of flag

       enddo

125   format(1x,i3,1x,9(1x,f10.3))
!126   format(1x,i2,1x,i2,1x,6(1x,f14.8))
!127   format(1x,i2,2(1x,f16.9))
!------------------------------------------------------------------------------
!  End calculations
!==============================================================================

   END SUBROUTINE flake_driver_run

end module flake_driver

module flake_driver_post
   use machine, only: kind_phys
   implicit none
   private
   public flake_driver_post_init, flake_driver_post_finalize, flake_driver_post_run 

contains
   subroutine flake_driver_post_init()
   end subroutine flake_driver_post_init

   subroutine flake_driver_post_finalize()
   end subroutine flake_driver_post_finalize

!> \section arg_table_flake_driver_post Argument Table
!! \htmlinclude flake_driver_post.html
!!
subroutine flake_driver_post_run (im, use_lake_model, h_ML, T_wML,  &
                           Tsurf, lakedepth, xz, zm, tref, tsfco,   &
                           errmsg, errflg)

!use machine , only : kind_phys
!==============================================================================

      implicit none
      integer, intent(in) :: im
!     integer, dimension(im), intent(in) :: islmsk

      real (kind=kind_phys), dimension(:), intent(in) ::               &
     &           lakedepth, tsurf, h_ML, t_wML

      real (kind=kind_phys),dimension(:),intent(inout) ::              &
     &           xz, zm, tref, tsfco   

      integer, dimension(:), intent(in) :: use_lake_model

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      integer :: i
      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      do I=1, im
         if(use_lake_model(i).eq.2) then
         write(0,*)'flake-post-use-lake-model= ',use_lake_model(i)
            xz(i) = lakedepth(i)
            zm(i) = h_ML(i)
            tref(i) = tsurf(i)
            tsfco(i) = t_wML(i)
         endif
      enddo


end subroutine flake_driver_post_run  

!---------------------------------
end module flake_driver_post

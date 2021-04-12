!>  \file sfc_nst.f
!!  This file contains the GFS NSST model.

!> This module contains the CCPP-compliant GFS near-surface sea temperature scheme.
      module sfc_nst

      contains

! \brief This subroutine is empty since there are no procedures that need to be done to initialize the GFS NSST code.
!! This subroutine is empty since there are no procedures that need to be done to initialize the GFS NSST code.
!!
      subroutine sfc_nst_init
      end subroutine sfc_nst_init

! \brief This subroutine is empty since there are no procedures that need to be done to finalize the GFS NSST code.
!! This subroutine is empty since there are no procedures that need to be done to finalize the GFS NSST code.
!!
      subroutine sfc_nst_finalize
      end subroutine sfc_nst_finalize

!>\defgroup gfs_nst_main GFS Near-Surface Sea Temperature Scheme Module
!> \brief This subroutine calls the Thermal Skin-layer and Diurnal Thermocline models to update the NSST profile.
!! \section arg_table_sfc_nst_run Argument Table
!! \htmlinclude sfc_nst_run.html
!!
!! \section NSST_general_algorithm GFS Near-Surface Sea Temperature Scheme General Algorithm
!> @{
      subroutine sfc_nst_run                                            &
     &     ( im, hvap, cp, hfus, jcal, eps, epsm1, rvrdm1, rd, rhw0,    &  ! --- inputs:
     &       pi, tgice, sbc, ps, u1, v1, t1, q1, tref, cm, ch,                 &
     &       prsl1, prslki, prsik1, prslk1, wet, xlon, sinlat,          &
     &       stress,                                                    &
     &       sfcemis, dlwflx, sfcnsw, rain, timestep, kdt, solhr,xcosz, &
     &       wind, flag_iter, flag_guess, nstf_name1, nstf_name4,       &
     &       nstf_name5, lprnt, ipr,                                    &
     &       tskin, tsurf, xt, xs, xu, xv, xz, zm, xtts, xzts, dt_cool, &  ! --- input/output:
     &       z_c,   c_0,   c_d,   w_0, w_d, d_conv, ifd, qrain,         &
     &       qsurf, gflux, cmm, chh, evap, hflx, ep, errmsg, errflg     &  ! --- outputs:
     &      )
!
! ===================================================================== !
!  description:                                                         !
!                                                                       !
!                                                                       !
!  usage:                                                               !
!                                                                       !
!    call sfc_nst                                                       !
!       inputs:                                                         !
!          ( im, ps, u1, v1, t1, q1, tref, cm, ch,                      !
!            prsl1, prslki, wet, xlon, sinlat, stress,                  !
!            sfcemis, dlwflx, sfcnsw, rain, timestep, kdt,solhr,xcosz,  !
!            wind,  flag_iter, flag_guess, nstf_name1, nstf_name4,      !
!            nstf_name5, lprnt, ipr,                                    !
!       input/outputs:                                                  !
!            tskin, tsurf, xt, xs, xu, xv, xz, zm, xtts, xzts, dt_cool, !
!            z_c, c_0,   c_d,   w_0, w_d, d_conv, ifd, qrain,           !
!  --   outputs:
!            qsurf, gflux, cmm, chh, evap, hflx, ep                     !
!           )
!                                                                       !
!                                                                       !
!  subprogram/functions called: w3movdat, iw3jdn, fpvs, density,        !
!       rhocoef, cool_skin, warm_layer, jacobi_temp.                    !
!                                                                       !
!  program history log:                                                 !
!         2007  -- xu li       createad original code                   !
!         2008  -- s. moorthi  adapted to the parallel version          !
!    may  2009  -- y.-t. hou   modified to include input lw surface     !
!                     emissivity from radiation. also replaced the      !
!                     often comfusing combined sw and lw suface         !
!                     flux with separate sfc net sw flux (defined       !
!                     as dn-up) and lw flux. added a program doc block. !
!    sep  2009 --  s. moorthi removed rcl and additional reformatting   !
!                     and optimization + made pa as input pressure unit.!
!         2009  -- xu li       recreatead the code                      !
!    feb  2010  -- s. moorthi added some changes made to the previous   !
!                  version                                              !
!    Jul  2016  -- X. Li, modify the diurnal warming event reset        !
!                                                                       !
!                                                                       !
!  ====================  definition of variables  ====================  !
!                                                                       !
!  inputs:                                                       size   !
!     im       - integer, horiz dimension                          1    !
!     ps       - real, surface pressure (pa)                       im   !
!     u1, v1   - real, u/v component of surface layer wind (m/s)   im   !
!     t1       - real, surface layer mean temperature ( k )        im   !
!     q1       - real, surface layer mean specific humidity        im   !
!     tref     - real, reference/foundation temperature ( k )      im   !
!     cm       - real, surface exchange coeff for momentum (m/s)   im   !
!     ch       - real, surface exchange coeff heat & moisture(m/s) im   !
!     prsl1    - real, surface layer mean pressure (pa)            im   !
!     prslki   - real,                                             im   !
!     prsik1   - real,                                             im   !
!     prslk1   - real,                                             im   !
!     wet      - logical, =T if any ocn/lake water (F otherwise)   im   !
!     icy      - logical, =T if any ice                            im   !
!     xlon     - real, longitude         (radians)                 im   !
!     sinlat   - real, sin of latitude                             im   !
!     stress   - real, wind stress       (n/m**2)                  im   !
!     sfcemis  - real, sfc lw emissivity (fraction)                im   !
!     dlwflx   - real, total sky sfc downward lw flux (w/m**2)     im   !
!     sfcnsw   - real, total sky sfc netsw flx into ocean (w/m**2) im   !
!     rain     - real, rainfall rate     (kg/m**2/s)               im   !
!     timestep - real, timestep interval (second)                  1    !
!     kdt      - integer, time step counter                        1    !
!     solhr    - real, fcst hour at the end of prev time step      1    !
!     xcosz    - real, consine of solar zenith angle               1    !
!     wind     - real, wind speed (m/s)                            im   !
!     flag_iter- logical, execution or not                         im   !
!                when iter = 1, flag_iter = .true. for all grids   im   !
!                when iter = 2, flag_iter = .true. when wind < 2   im   !
!                for both land and ocean (when nstf_name1 > 0)     im   !
!     flag_guess-logical, .true.=  guess step to get CD et al      im   !
!                when iter = 1, flag_guess = .true. when wind < 2  im   !
!                when iter = 2, flag_guess = .false. for all grids im   !
!     nstf_name - integers , NSST related flag parameters          1    !
!                nstf_name1 : 0 = NSSTM off                        1    !
!                             1 = NSSTM on but uncoupled           1    !
!                             2 = NSSTM on and coupled             1    !
!                nstf_name4 : zsea1 in mm                          1    !
!                nstf_name5 : zsea2 in mm                          1    !
!     lprnt    - logical, control flag for check print out         1    !
!     ipr      - integer, grid index for check print out           1    !
!                                                                       !
!  input/outputs:
! li added for oceanic components
!     tskin    - real, ocean surface skin temperature ( k )        im   !
!     tsurf    - real, the same as tskin ( k ) but for guess run   im   !
!     xt       - real, heat content in dtl                         im   !
!     xs       - real, salinity  content in dtl                    im   !
!     xu       - real, u-current content in dtl                    im   !
!     xv       - real, v-current content in dtl                    im   !
!     xz       - real, dtl thickness                               im   !
!     zm       - real, mxl thickness                               im   !
!     xtts     - real, d(xt)/d(ts)                                 im   !
!     xzts     - real, d(xz)/d(ts)                                 im   !
!     dt_cool  - real, sub-layer cooling amount                    im   !
!     d_conv   - real, thickness of free convection layer (fcl)    im   !
!     z_c      - sub-layer cooling thickness                       im   !
!     c_0      - coefficient1 to calculate d(tz)/d(ts)             im   !
!     c_d      - coefficient2 to calculate d(tz)/d(ts)             im   !
!     w_0      - coefficient3 to calculate d(tz)/d(ts)             im   !
!     w_d      - coefficient4 to calculate d(tz)/d(ts)             im   !
!     ifd      - real, index to start dtlm run or not              im   !
!     qrain    - real, sensible heat flux due to rainfall (watts)  im   !

!  outputs:                                                             !

!     qsurf    - real, surface air saturation specific humidity    im   !
!     gflux    - real, soil heat flux (w/m**2)                     im   !
!     cmm      - real,                                             im   !
!     chh      - real,                                             im   !
!     evap     - real, evaperation from latent heat flux           im   !
!     hflx     - real, sensible heat flux                          im   !
!     ep       - real, potential evaporation                       im   !
!                                                                       !
! ===================================================================== !
      use machine , only : kind_phys
      use funcphys, only : fpvs
      use date_def, only : idate
      use module_nst_water_prop, only: get_dtzm_point
      use module_nst_parameters, only : t0k,cp_w,omg_m,omg_sh,          &
     &    sigma_r,solar_time_6am,ri_c,z_w_max,delz,wd_max,              &
     &    rad2deg,const_rot,tau_min,tw_max,sst_max
      use module_nst_water_prop, only: solar_time_from_julian,          &
     &                                 density,rhocoef,compjd,grv       &
     &,                                sw_ps_9b
      use nst_module, only : cool_skin,dtm_1p,cal_w,cal_ttop,           &
     &                       convdepth,dtm_1p_fca,dtm_1p_tla,           &
     &                       dtm_1p_mwa,dtm_1p_mda,dtm_1p_mta,          &
     &                       dtl_reset
!
      implicit none

      integer, parameter :: kp = kind_phys
!
!  ---  constant parameters:
      real (kind=kind_phys), parameter :: f24   = 24.0_kp     ! hours/day
      real (kind=kind_phys), parameter :: f1440 = 1440.0_kp   ! minutes/day
      real (kind=kind_phys), parameter :: czmin = 0.0001_kp   ! cos(89.994)
      real (kind=kind_phys), parameter :: zero  = 0.0_kp, one = 1.0_kp


!  ---  inputs:
      integer, intent(in) :: im, kdt, ipr, nstf_name1, nstf_name4,      &
     &       nstf_name5
      real (kind=kind_phys), intent(in) :: hvap, cp, hfus, jcal, eps,   &
     &       epsm1, rvrdm1, rd, rhw0, sbc, pi, tgice
      real (kind=kind_phys), dimension(im), intent(in) :: ps, u1, v1,   &
     &       t1, q1, tref, cm, ch, prsl1, prslki, prsik1, prslk1,       &
     &       xlon,xcosz,                                                &
     &       sinlat, stress, sfcemis, dlwflx, sfcnsw, rain, wind
      real (kind=kind_phys), intent(in) :: timestep
      real (kind=kind_phys), intent(in) :: solhr

      logical, dimension(im), intent(in) :: flag_iter, flag_guess, wet
!    &,      icy
      logical,                intent(in) :: lprnt

!  ---  input/outputs:
! control variables of dtl system (5+2) and sl (2) and coefficients for d(tz)/d(ts) calculation
      real (kind=kind_phys), dimension(im), intent(inout) :: tskin,     &
     &      tsurf, xt, xs, xu, xv, xz, zm, xtts, xzts, dt_cool,         &
     &      z_c, c_0, c_d, w_0, w_d, d_conv, ifd, qrain

!  ---  outputs:
      real (kind=kind_phys), dimension(im), intent(inout) ::            &
     &       qsurf, gflux, cmm, chh, evap, hflx, ep

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

!
!     locals
!
      integer :: k,i
!
      real (kind=kind_phys), dimension(im) ::  q0, qss, rch,
     &                     rho_a, theta1, tv1, wndmag

      real(kind=kind_phys) elocp,tem,cpinv,hvapi
!
!    nstm related prognostic fields
!
      logical flag(im)
      real (kind=kind_phys), dimension(im) ::
     &   xt_old, xs_old, xu_old, xv_old, xz_old,zm_old,xtts_old,
     &   xzts_old, ifd_old, tref_old, tskin_old, dt_cool_old,z_c_old

      real(kind=kind_phys) ulwflx(im), nswsfc(im)
!     real(kind=kind_phys) rig(im),
!    &                     ulwflx(im),dlwflx(im),
!    &                     slrad(im),nswsfc(im)
      real(kind=kind_phys) alpha,beta,rho_w,f_nsol,sss,sep,
     &                     cosa,sina,taux,tauy,grav,dz,t0,ttop0,ttop

      real(kind=kind_phys) le,fc,dwat,dtmp,wetc,alfac,ustar_a,rich
      real(kind=kind_phys) rnl_ts,hs_ts,hl_ts,rf_ts,q_ts
      real(kind=kind_phys) fw,q_warm
      real(kind=kind_phys) t12,alon,tsea,sstc,dta,dtz
      real(kind=kind_phys) zsea1,zsea2,soltim

!  external functions called: iw3jdn
      integer :: iw3jdn
!======================================================================================================
cc
      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      cpinv = one/cp
      hvapi = one/hvap
      elocp = hvap/cp

      sss = 34.0_kp             ! temporarily, when sea surface salinity data is not ready
!
! flag for open water and where the iteration is on
!
      do i = 1, im
!       flag(i) = wet(i) .and. .not.icy(i) .and. flag_iter(i)
        flag(i) = wet(i) .and. flag_iter(i)
      enddo
!
!  save nst-related prognostic fields for guess run
!
      do i=1, im
!       if(wet(i) .and. .not.icy(i) .and. flag_guess(i)) then
        if(wet(i) .and. flag_guess(i)) then
          xt_old(i)      = xt(i)
          xs_old(i)      = xs(i)
          xu_old(i)      = xu(i)
          xv_old(i)      = xv(i)
          xz_old(i)      = xz(i)
          zm_old(i)      = zm(i)
          xtts_old(i)    = xtts(i)
          xzts_old(i)    = xzts(i)
          ifd_old(i)     = ifd(i)
          tskin_old(i)   = tskin(i)
          dt_cool_old(i) = dt_cool(i)
          z_c_old(i)     = z_c(i)
        endif
      enddo


!  --- ...  initialize variables. all units are m.k.s. unless specified.
!           ps is in pascals, wind is wind speed, theta1 is surface air
!           estimated from level 1 temperature, rho_a is air density and
!           qss is saturation specific humidity at the water surface
!!
      do i = 1, im
        if ( flag(i) ) then

          nswsfc(i) = sfcnsw(i) ! net solar radiation at the air-sea surface (positive downward)
          wndmag(i) = sqrt(u1(i)*u1(i) + v1(i)*v1(i))

          q0(i)     = max(q1(i), 1.0e-8_kp)
#ifdef GSD_SURFACE_FLUXES_BUGFIX
          theta1(i) = t1(i) / prslk1(i) ! potential temperature at the middle of lowest model layer
#else
          theta1(i) = t1(i) * prslki(i)
#endif
          tv1(i)    = t1(i) * (one + rvrdm1*q0(i))
          rho_a(i)  = prsl1(i) / (rd*tv1(i))
          qss(i)    = fpvs(tsurf(i))                          ! pa
          qss(i)    = eps*qss(i) / (ps(i) + epsm1*qss(i))     ! pa
!
          evap(i)    = zero
          hflx(i)    = zero
          gflux(i)   = zero
          ep(i)      = zero

!  --- ...  rcp = rho cp ch v

          rch(i)     = rho_a(i) * cp * ch(i) * wind(i)
          cmm(i)     = cm (i)   * wind(i)
          chh(i)     = rho_a(i) * ch(i) * wind(i)

!> - Calculate latent and sensible heat flux over open water with tskin.
!           at previous time step
          evap(i)    = elocp * rch(i) * (qss(i) - q0(i))
          qsurf(i)   = qss(i)
#ifdef GSD_SURFACE_FLUXES_BUGFIX
          hflx(i)    = rch(i) * (tsurf(i)/prsik1(i) - theta1(i))
#else
          hflx(i)    = rch(i) * (tsurf(i) - theta1(i))
#endif

!     if (lprnt .and. i == ipr) print *,' tskin=',tskin(i),' theta1=',
!    & theta1(i),' hflx=',hflx(i),' t1=',t1(i),'prslki=',prslki(i)
!    &,' tsurf=',tsurf(i)
        endif
      enddo

! run nst model: dtm + slm
!
      zsea1 = 0.001_kp*real(nstf_name4)
      zsea2 = 0.001_kp*real(nstf_name5)

!> - Call module_nst_water_prop::density() to compute sea water density.
!> - Call module_nst_water_prop::rhocoef() to compute thermal expansion
!! coefficient (\a alpha) and saline contraction coefficient (\a beta).
      do i = 1, im
        if ( flag(i) ) then
          tsea      = tsurf(i)
          t12       = tsea*tsea
          ulwflx(i) = sfcemis(i) * sbc * t12 * t12
          alon      = xlon(i)*rad2deg
          grav      = grv(sinlat(i))
          soltim  = mod(alon/15.0_kp + solhr, 24.0_kp)*3600.0_kp
          call density(tsea,sss,rho_w)                     ! sea water density
          call rhocoef(tsea,sss,rho_w,alpha,beta)          ! alpha & beta
!
!> - Calculate sensible heat flux (\a qrain) due to rainfall.
!
          le       = (2.501_kp-0.00237_kp*tsea)*1e6_kp
          dwat     = 2.11e-5_kp*(t1(i)/t0k)**1.94_kp               ! water vapor diffusivity
          dtmp     = (one+3.309e-3_kp*(t1(i)-t0k)-1.44e-6_kp*(t1(i)-t0k)
     &             * (t1(i)-t0k))*0.02411_kp/(rho_a(i)*cp)         ! heat diffusivity
          wetc     = 622.0_kp*le*qss(i)/(rd*t1(i)*t1(i))
          alfac    = one / (one + (wetc*le*dwat)/(cp*dtmp))        ! wet bulb factor
          tem      = (1.0e3_kp * rain(i) / rho_w) * alfac * cp_w
          qrain(i) =  tem * (tsea-t1(i)+1.0e3_kp*(qss(i)-q0(i))*le/cp)

!> - Calculate input non solar heat flux as upward = positive to models here

          f_nsol   = hflx(i) + evap(i) + ulwflx(i) - dlwflx(i)
     &             + omg_sh*qrain(i)

!     if (lprnt .and. i == ipr) print *,' f_nsol=',f_nsol,' hflx=',
!    &hflx(i),' evap=',evap(i),' ulwflx=',ulwflx(i),' dlwflx=',dlwflx(i)
!    &,' omg_sh=',omg_sh,' qrain=',qrain(i)

          sep      = sss*(evap(i)/le-rain(i))/rho_w
          ustar_a  = sqrt(stress(i)/rho_a(i))          ! air friction velocity
!
!  sensitivities of heat flux components to ts
!
          rnl_ts = 4.0_kp*sfcemis(i)*sbc*tsea*tsea*tsea     ! d(rnl)/d(ts)
          hs_ts  = rch(i)
          hl_ts  = rch(i)*elocp*eps*hvap*qss(i)/(rd*t12)
          rf_ts  = tem * (one+rch(i)*hl_ts)
          q_ts   = rnl_ts + hs_ts + hl_ts + omg_sh*rf_ts
!
!> - Call cool_skin(), which is the sub-layer cooling parameterization
!! (Fairfall et al. (1996) \cite fairall_et_al_1996).
! & calculate c_0, c_d
!
          call cool_skin(ustar_a,f_nsol,nswsfc(i),evap(i),sss,alpha,beta
     &,                  rho_w,rho_a(i),tsea,q_ts,hl_ts,grav,le
     &,                  dt_cool(i),z_c(i),c_0(i),c_d(i))

          tem  = one / wndmag(i)
          cosa = u1(i)*tem
          sina = v1(i)*tem
          taux = max(stress(i),tau_min)*cosa
          tauy = max(stress(i),tau_min)*sina
          fc   = const_rot*sinlat(i)
!
!  Run DTM-1p system.
!
          if ( (soltim > solar_time_6am .and. ifd(i) == zero) ) then
          else
            ifd(i) = one
!
!     calculate fcl thickness with current forcing and previous time's profile
!
!     if (lprnt .and. i == ipr) print *,' beg xz=',xz(i)

!> - Call convdepth() to calculate depth for convective adjustments.
            if ( f_nsol > zero .and. xt(i) > zero ) then
              call convdepth(kdt,timestep,nswsfc(i),f_nsol,sss,sep,rho_w
     &,                      alpha,beta,xt(i),xs(i),xz(i),d_conv(i))
            else
              d_conv(i) = zero
            endif

!     if (lprnt .and. i == ipr) print *,' beg xz1=',xz(i)
!
!    determine rich: wind speed dependent (right now)
!
!           if ( wind(i) < 1.0 ) then
!             rich = 0.25 + 0.03*wind(i)
!           elseif ( wind(i) >= 1.0 .and. wind(i) < 1.5 ) then
!             rich = 0.25 + 0.1*wind(i)
!           elseif ( wind(i) >= 1.5 .and. wind(i) < 6.0 ) then
!             rich = 0.25 + 0.6*wind(i)
!           elseif ( wind(i) >= 6.0 ) then
!             rich = 0.25 + min(0.8*wind(i),0.50)
!           endif

            rich = ri_c

!> - Call the diurnal thermocline layer model dtm_1p().
            call dtm_1p(kdt,timestep,rich,taux,tauy,nswsfc(i),
     &                  f_nsol,sss,sep,q_ts,hl_ts,rho_w,alpha,beta,alon,
     &                  sinlat(i),soltim,grav,le,d_conv(i),
     &                  xt(i),xs(i),xu(i),xv(i),xz(i),xzts(i),xtts(i))

!     if (lprnt .and. i == ipr) print *,' beg xz2=',xz(i)

!  apply mda
            if ( xt(i) > zero ) then
!>  - If \a dtl heat content \a xt > 0.0, call dtm_1p_mda() to apply
!!  minimum depth adjustment (mda).
              call dtm_1p_mda(xt(i),xtts(i),xz(i),xzts(i))
              if ( xz(i) >= z_w_max ) then
!>   - If \a dtl thickness >= module_nst_parameters::z_w_max, call dtl_reset()
!! to reset xt/xs/x/xv to zero, and xz to module_nst_parameters::z_w_max.
                call dtl_reset(xt(i),xs(i),xu(i),xv(i),xz(i),xtts(i),
     &                                                       xzts(i))

!     if (lprnt .and. i == ipr) print *,' beg xz3=',xz(i),' z_w_max='
!    &,z_w_max
              endif

!  apply fca
              if ( d_conv(i) > zero ) then
!>  - If thickness of free convection layer > 0.0, call dtm_1p_fca()
!! to apply free convection adjustment.
!>   - If \a dtl thickness >= module_nst_parameters::z_w_max(), call dtl_reset()
!! to reset xt/xs/x/xv to zero, and xz to module_nst_parameters::z_w_max().
                call dtm_1p_fca(d_conv(i),xt(i),xtts(i),xz(i),xzts(i))
                if ( xz(i) >= z_w_max ) then
                  call dtl_reset
     &              (xt(i),xs(i),xu(i),xv(i),xz(i),xzts(i),xtts(i))
                endif
              endif

!     if (lprnt .and. i == ipr) print *,' beg xz4=',xz(i)

!  apply tla
              dz = min(xz(i),max(d_conv(i),delz))
!
!>  - Call sw_ps_9b() to compute the fraction of the solar radiation
!! absorbed by the depth \a delz (Paulson and Simpson (1981) \cite paulson_and_simpson_1981).
!! And calculate the total heat absorbed in warm layer.
              call sw_ps_9b(delz,fw)
              q_warm = fw*nswsfc(i)-f_nsol    !total heat absorbed in warm layer

!>  - Call cal_ttop() to calculate the diurnal warming amount at the top layer with
!! thickness of \a dz.
              if ( q_warm > zero ) then
                call cal_ttop(kdt,timestep,q_warm,rho_w,dz,
     &                        xt(i),xz(i),ttop0)

!     if (lprnt .and. i == ipr) print *,' d_conv=',d_conv(i),' delz=',
!    &delz,' kdt=',kdt,' timestep=',timestep,' nswsfc=',nswsfc(i),
!    &' f_nsol=',f_nsol,' rho_w=',rho_w,' dz=',dz,' xt=',xt(i),
!    &' xz=',xz(i),' qrain=',qrain(i)

                ttop = ((xt(i)+xt(i))/xz(i))*(one-dz/((xz(i)+xz(i))))

!     if (lprnt .and. i == ipr) print *,' beg xz4a=',xz(i)
!    &,' ttop=',ttop,' ttop0=',ttop0,' xt=',xt(i),' dz=',dz
!    &,' xznew=',(xt(i)+sqrt(xt(i)*(xt(i)-dz*ttop0)))/ttop0

!>  - Call dtm_1p_tla() to apply top layer adjustment.
                if ( ttop > ttop0 ) then
                  call dtm_1p_tla(dz,ttop0,xt(i),xtts(i),xz(i),xzts(i))

!     if (lprnt .and. i == ipr) print *,' beg xz4b=',xz(i),'z_w_max=',
!    &z_w_max
                  if ( xz(i) >= z_w_max ) then
                    call dtl_reset
     &                   (xt(i),xs(i),xu(i),xv(i),xz(i),xzts(i),xtts(i))
                  endif
                endif
              endif           ! if ( q_warm > 0.0 ) then

!     if (lprnt .and. i == ipr) print *,' beg xz5=',xz(i)

!  apply mwa
!>  - Call dt_1p_mwa() to apply maximum warming adjustment.
              t0 = (xt(i)+xt(i))/xz(i)
              if ( t0 > tw_max ) then
                call dtm_1p_mwa(xt(i),xtts(i),xz(i),xzts(i))
                if ( xz(i) >= z_w_max ) then
                  call dtl_reset
     &                 (xt(i),xs(i),xu(i),xv(i),xz(i),xzts(i),xtts(i))
                endif
              endif

!     if (lprnt .and. i == ipr) print *,' beg xz6=',xz(i)

!  apply mta
!>  - Call dtm_1p_mta() to apply maximum temperature adjustment.
       sstc = tref(i) + (xt(i)+xt(i))/xz(i) - dt_cool(i)

              if ( sstc > sst_max ) then
                dta = sstc - sst_max
                call  dtm_1p_mta(dta,xt(i),xtts(i),xz(i),xzts(i))
!               write(*,'(a,f3.0,7f8.3)') 'mta, sstc,dta :',islimsk(i),
!    &          sstc,dta,tref(i),xt(i),xz(i),2.0*xt(i)/xz(i),dt_cool(i)
               if ( xz(i) >= z_w_max ) then
                  call dtl_reset
     &                 (xt(i),xs(i),xu(i),xv(i),xz(i),xzts(i),xtts(i))
                endif
              endif
!
            endif             ! if ( xt(i) > 0.0 ) then
!           reset dtl at midnight and when solar zenith angle > 89.994 degree
            if ( abs(soltim) < 2.0_kp*timestep ) then
              call dtl_reset
     &           (xt(i),xs(i),xu(i),xv(i),xz(i),xzts(i),xtts(i))
            endif

          endif             ! if (solar_time > solar_time_6am .and. ifd(i) == 0.0 ) then: too late to start the first day

!     if (lprnt .and. i == ipr) print *,' beg xz7=',xz(i)

!     update tsurf  (when flag(i) .eqv. .true. )
!>  - Call get_dtzm_point() to computes \a dtz and \a tsurf.
          call get_dtzm_point(xt(i),xz(i),dt_cool(i),z_c(i),
     &                        zsea1,zsea2,dtz)
          tsurf(i) = max(tgice, tref(i) + dtz )

!     if (lprnt .and. i == ipr) print *,' tsurf=',tsurf(i),' tref=',
!    &tref(i),' xz=',xz(i),' dt_cool=',dt_cool(i)

!>  - Call cal_w() to calculate \a w_0 and \a w_d.
          if ( xt(i) > zero ) then
            call cal_w(kdt,xz(i),xt(i),xzts(i),xtts(i),w_0(i),w_d(i))
          else
            w_0(i) = zero
            w_d(i) = zero
          endif

!         if ( xt(i) > 0.0 ) then
!           rig(i) = grav*xz(i)*xz(i)*(alpha*xt(i)-beta*xs(i))
!    &             /(2.0*(xu(i)*xu(i)+xv(i)*xv(i)))
!         else
!           rig(i) = 0.25
!         endif

!         qrain(i) = rig(i)
          zm(i) = wind(i)

        endif
      enddo

! restore nst-related prognostic fields for guess run
      do i=1, im
!       if (wet(i) .and. .not.icy(i)) then
        if (wet(i)) then
          if (flag_guess(i)) then    ! when it is guess of
            xt(i)      = xt_old(i)
            xs(i)      = xs_old(i)
            xu(i)      = xu_old(i)
            xv(i)      = xv_old(i)
            xz(i)      = xz_old(i)
            zm(i)      = zm_old(i)
            xtts(i)    = xtts_old(i)
            xzts(i)    = xzts_old(i)
            ifd(i)     = ifd_old(i)
            tskin(i)   = tskin_old(i)
            dt_cool(i) = dt_cool_old(i)
            z_c(i)     = z_c_old(i)
          else
!
!         update tskin when coupled and not guess run
!         (all other NSST variables have been updated in this case)
!
            if ( nstf_name1 > 1 ) then
              tskin(i) = tsurf(i)
            endif               ! if nstf_name1 > 1 then
          endif                 ! if flag_guess(i) then
        endif                   ! if wet(i) .and. .not.icy(i) then
      enddo

!     if (lprnt .and. i == ipr) print *,' beg xz8=',xz(i)

      if ( nstf_name1 > 1 ) then
!> - Calculate latent and sensible heat flux over open water with updated tskin
!!      for the grids of open water and the iteration is on.
        do i = 1, im
          if ( flag(i) ) then
            qss(i)   = fpvs( tskin(i) )
            qss(i)   = eps*qss(i) / (ps(i) + epsm1*qss(i))
            qsurf(i) = qss(i)
            evap(i)  = elocp*rch(i) * (qss(i) - q0(i))
#ifdef GSD_SURFACE_FLUXES_BUGFIX
            hflx(i)  = rch(i) * (tskin(i)/prsik1(i) - theta1(i))
#else
            hflx(i)  = rch(i) * (tskin(i) - theta1(i))
#endif
          endif
        enddo
      endif                   ! if ( nstf_name1 > 1 ) then

!
      do i=1,im
        if ( flag(i) ) then
          tem     = one / rho_a(i)
          hflx(i) = hflx(i) * tem * cpinv
          evap(i) = evap(i) * tem * hvapi
        endif
      enddo
!
!     if (lprnt) print *,' tskin=',tskin(ipr)

      return
      end subroutine sfc_nst_run
!> @}
      end module sfc_nst

!> This module contains the CCPP-compliant GFS near-surface sea temperature pre
!! interstitial codes.
      module sfc_nst_pre

      contains

! \defgroup GFS_NSST_PRE GFS Near-Surface Sea Temperature Pre
!!
!! The NSST scheme is one of the three schemes used to represent the
!! surface in the GFS physics suite. The other two are the Noah land
!! surface model and the sice simplified ice model.
!!
      subroutine sfc_nst_pre_init
      end subroutine sfc_nst_pre_init

      subroutine sfc_nst_pre_finalize
      end subroutine sfc_nst_pre_finalize

!! \section arg_table_sfc_nst_pre_run Argument Table
!! \htmlinclude sfc_nst_pre_run.html
!!
!> \section NSST_general_pre_algorithm General Algorithm
!! @{
      subroutine sfc_nst_pre_run
     &    (im,wet,tgice,tsfco,tsfc_wat,tsurf_wat,tseal,xt,xz,dt_cool,
     &     z_c, tref, cplflx, oceanfrac, nthreads, errmsg, errflg)

      use machine , only : kind_phys
      use module_nst_water_prop, only: get_dtzm_2d

      implicit none

      integer, parameter :: kp = kind_phys

!  ---  inputs:
      integer, intent(in) :: im, nthreads
      logical, dimension(im), intent(in) :: wet
      real (kind=kind_phys), intent(in) :: tgice
      real (kind=kind_phys), dimension(im), intent(in) ::
     &  tsfco, tsfc_wat, xt, xz, dt_cool, z_c, oceanfrac
      logical, intent(in) :: cplflx

!  ---  input/outputs:
      real (kind=kind_phys), dimension(im), intent(inout) ::
     &    tsurf_wat, tseal, tref

!  ---  outputs:
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

!  ---  locals
      integer :: i
      real(kind=kind_phys), parameter :: zero = 0.0_kp,
     &                                   one  = 1.0_kp,
     &                                   half = 0.5_kp,
     &                                   omz1 = 2.0_kp
      real(kind=kind_phys) :: tem1, tem2, dnsst
      real(kind=kind_phys), dimension(im) :: dtzm,z_c_0

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      do i=1,im
        if (wet(i)) then
!          tem         = (oro(i)-oro_uf(i)) * rlapse
          ! DH* 20190927 simplyfing this code because tem is zero
          !tem          = zero
          !tseal(i)     = tsfc_wat(i)  + tem
          tseal(i)      = tsfc_wat(i)
          !tsurf_wat(i) = tsurf_wat(i) + tem
          ! *DH
        endif
      enddo

!
!   update tsfc & tref with T1 from OGCM & NSST Profile if coupled
!
      if (cplflx) then
        z_c_0 = 0.0
        call get_dtzm_2d (xt,  xz, dt_cool,                             &
     &                    z_c_0, wet, zero, omz1, im, 1, nthreads, dtzm)
        do i=1,im
          if (wet(i) .and. oceanfrac(i) > zero) then
!           dnsst   = tsfc_wat(i) - tref(i)                 !  retrive/get difference of Ts and Tf
            tref(i) = max(tgice, tsfco(i) - dtzm(i))        !  update Tf with T1 and NSST T-Profile
!           tsfc_wat(i) = max(271.2,tref(i) + dnsst)        !  get Ts updated due to Tf update
!           tseal(i)    = tsfc_wat(i)
            if (abs(xz(i)) > zero) then
              tem2 = one / xz(i)
            else
              tem2 = zero
            endif
            tseal(i)     = tref(i) + (xt(i)+xt(i)) * tem2 - dt_cool(i)
            tsurf_wat(i) = tseal(i)
          endif
        enddo
      endif

      return
      end subroutine sfc_nst_pre_run
!! @}
      end module sfc_nst_pre

!> This module contains the CCPP-compliant GFS near-surface sea temperature post
!! interstitial codes.
      module sfc_nst_post

      contains

! \defgroup GFS_NSST_POST GFS Near-Surface Sea Temperature Post
!! \brief Brief description of the parameterization
!!
      subroutine sfc_nst_post_init
      end subroutine sfc_nst_post_init

! \brief Brief description of the subroutine
!!
      subroutine sfc_nst_post_finalize
      end subroutine sfc_nst_post_finalize

!> \brief Brief description of the subroutine
!!
!! \section arg_table_sfc_nst_post_run Argument Table
!! \htmlinclude sfc_nst_post_run.html
!!
! \section NSST_general_post_algorithm General Algorithm
!
! \section NSST_detailed_post_algorithm Detailed Algorithm
! @{
      subroutine sfc_nst_post_run                                       &
     &     ( im, kdt, rlapse, tgice, wet, icy, oro, oro_uf, nstf_name1, &
     &       nstf_name4, nstf_name5, xt, xz, dt_cool, z_c, tref, xlon,  &
     &       tsurf_wat, tsfc_wat, nthreads, dtzm, errmsg, errflg        &
     &     )

      use machine , only : kind_phys
      use module_nst_water_prop, only: get_dtzm_2d

      implicit none

      integer, parameter :: kp = kind_phys

!  ---  inputs:
      integer, intent(in) :: im, kdt, nthreads
      logical, dimension(im), intent(in) :: wet, icy
      real (kind=kind_phys), intent(in) :: rlapse, tgice
      real (kind=kind_phys), dimension(im), intent(in) :: oro, oro_uf
      integer, intent(in) :: nstf_name1, nstf_name4, nstf_name5
      real (kind=kind_phys), dimension(im), intent(in) :: xt, xz,       &
     &      dt_cool, z_c, tref, xlon

!  ---  input/outputs:
      real (kind=kind_phys), dimension(im), intent(inout) :: tsurf_wat, &
     &      tsfc_wat

!  ---  outputs:
      real (kind=kind_phys), dimension(size(xlon,1)), intent(out) ::    &
     &      dtzm

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

!  ---  locals
      integer :: i
      real(kind=kind_phys) :: zsea1, zsea2

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

!     if (lprnt) print *,' tseaz2=',tseal(ipr),' tref=',tref(ipr),
!    &     ' dt_cool=',dt_cool(ipr),' dt_warm=',2.0*xt(ipr)/xz(ipr),
!    &     ' kdt=',kdt

!      do i = 1, im
!        if (wet(i) .and. .not. icy(i)) then
!          tsurf_wat(i) = tsurf_wat(i) - (oro(i)-oro_uf(i)) * rlapse
!        endif
!      enddo

!  --- ...  run nsst model  ... ---

      if (nstf_name1 > 1) then
        zsea1 = 0.001_kp*real(nstf_name4)
        zsea2 = 0.001_kp*real(nstf_name5)
        call get_dtzm_2d (xt, xz, dt_cool, z_c, wet, zsea1, zsea2,      &
     &                    im, 1, nthreads, dtzm)
        do i = 1, im
!         if (wet(i) .and. .not.icy(i)) then
!         if (wet(i) .and. (frac_grid .or. .not. icy(i))) then
          if (wet(i)) then
            tsfc_wat(i) = max(tgice, tref(i) + dtzm(i))
!           tsfc_wat(i) = max(271.2, tref(i) + dtzm(i)) -  &
!                           (oro(i)-oro_uf(i))*rlapse
          endif
        enddo
      endif

!     if (lprnt) print *,' tseaz2=',tsea(ipr),' tref=',tref(ipr),   &
!    &    ' dt_cool=',dt_cool(ipr),' dt_warm=',dt_warm(ipr),' kdt=',kdt

      return
      end subroutine sfc_nst_post_run

      end module sfc_nst_post

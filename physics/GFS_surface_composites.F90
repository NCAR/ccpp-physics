!> \file GFS_surface_composites.F90
!!  Contains code related to generating composites for all GFS surface schemes.

module GFS_surface_composites_pre

   use machine, only: kind_phys
   use physparam, only : iemsflg

   implicit none

   private

   public GFS_surface_composites_pre_init, GFS_surface_composites_pre_finalize, GFS_surface_composites_pre_run

   real(kind=kind_phys), parameter :: zero = 0.0_kind_phys, one = 1.0_kind_phys, epsln = 1.0e-10_kind_phys

   real(kind=kind_phys), parameter :: huge      = 9.9692099683868690E36 ! NetCDF float FillValue

contains

   subroutine GFS_surface_composites_pre_init ()
   end subroutine GFS_surface_composites_pre_init

   subroutine GFS_surface_composites_pre_finalize()
   end subroutine GFS_surface_composites_pre_finalize

!> \section arg_table_GFS_surface_composites_pre_run Argument Table
!! \htmlinclude GFS_surface_composites_pre_run.html
!!
   subroutine GFS_surface_composites_pre_run (im, flag_init, flag_restart, lkm, lsm, lsm_noahmp, lsm_ruc, frac_grid,      &
                                 flag_cice, cplflx, cplwav2atm, landfrac, lakefrac, lakedepth, oceanfrac, frland,         &
                                 dry, icy, use_flake, ocean, wet, hice, cice, zorlo, zorll, zorli,                        &
                                 snowd, snowd_wat, snowd_lnd, snowd_ice, tprcp, tprcp_wat,                                &
                                 tprcp_lnd, tprcp_ice, uustar, uustar_wat, uustar_lnd, uustar_ice,                        &
                                 weasd, weasd_wat, weasd_lnd, weasd_ice, ep1d_ice, tsfc, tsfco, tsfcl, tsfc_wat,          &
                                 tsfc_lnd, tsfc_ice, tisfc, tice, tsurf_wat, tsurf_lnd, tsurf_ice,                        &
                                 gflx_ice, tgice, islmsk, islmsk_cice, slmsk, semis_rad, semis_wat, semis_lnd, semis_ice, &
                                 emis_lnd, emis_ice, qss, qss_wat, qss_lnd, qss_ice, hflx, hflx_wat, hflx_lnd, hflx_ice,  &
                                 min_lakeice, min_seaice, errmsg, errflg)

      implicit none

      ! Interface variables
      integer,                             intent(in   ) :: im, lkm
      integer,                             intent(in   ) :: lsm, lsm_noahmp, lsm_ruc
      logical,                             intent(in   ) :: flag_init, flag_restart, frac_grid, cplflx, cplwav2atm
      logical, dimension(:),              intent(inout)  :: flag_cice
      logical,              dimension(:), intent(inout)  :: dry, icy, use_flake, ocean, wet
      real(kind=kind_phys), dimension(:), intent(in   )  :: landfrac, lakefrac, lakedepth, oceanfrac
      real(kind=kind_phys), dimension(:), intent(inout)  :: cice, hice
      real(kind=kind_phys), dimension(:), intent(  out)  :: frland
      real(kind=kind_phys), dimension(:), intent(in   )  :: snowd, tprcp, uustar, weasd, qss, hflx

      real(kind=kind_phys), dimension(:), intent(inout)  :: tsfc, tsfco, tsfcl, tisfc
      real(kind=kind_phys), dimension(:), intent(inout)  :: snowd_wat, snowd_lnd, snowd_ice, tprcp_wat, &
                    tprcp_lnd, tprcp_ice, tsfc_wat, tsfc_lnd, tsfc_ice, tsurf_wat,tsurf_lnd, tsurf_ice, &
                    uustar_wat, uustar_lnd, uustar_ice, weasd_wat, weasd_lnd, weasd_ice,                &
                    qss_wat, qss_lnd, qss_ice, hflx_wat, hflx_lnd, hflx_ice, ep1d_ice, gflx_ice
      real(kind=kind_phys), dimension(:), intent(  out)  :: tice
      real(kind=kind_phys),                intent(in   ) :: tgice
      integer,              dimension(:), intent(inout)  :: islmsk, islmsk_cice
      real(kind=kind_phys), dimension(:), intent(in   )  :: semis_rad
      real(kind=kind_phys), dimension(:), intent(inout)  :: semis_wat, semis_lnd, semis_ice, slmsk
      real(kind=kind_phys), dimension(:), intent(inout)  :: emis_lnd, emis_ice
      real(kind=kind_phys),                intent(in   ) :: min_lakeice, min_seaice
      !
      real(kind=kind_phys), dimension(:), intent(inout)  :: zorlo, zorll, zorli
      !
      real(kind=kind_phys), parameter :: timin = 173.0_kind_phys  ! minimum temperature allowed for snow/ice

      ! CCPP error handling
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      ! Local variables
      integer :: i

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      if (frac_grid) then  ! cice is ice fraction wrt water area
        do i=1,im
          frland(i) = landfrac(i)
          if (frland(i) > zero) dry(i) = .true.
          if (frland(i) < one) then
            if (oceanfrac(i) > zero) then
              if (cice(i) >= min_seaice) then
                icy(i)  = .true.
                tisfc(i) = max(timin, min(tisfc(i), tgice))
                if (cplflx)  then
                  islmsk_cice(i) = 4
                  flag_cice(i)   = .true.
                else
                  islmsk_cice(i) = 2
                endif
                islmsk(i) = 2
              else
                cice(i)        = zero
                hice(i)        = zero
                flag_cice(i)   = .false.
                islmsk_cice(i) = 0
                islmsk(i)      = 0
              endif
              if (cice(i) < one) then
                wet(i) = .true. ! some open ocean
                if (.not. cplflx .and. icy(i)) tsfco(i) = max(tisfc(i), tgice)
              endif
            else
              if (cice(i) >= min_lakeice) then
                icy(i) = .true.
                islmsk(i) = 2
                tisfc(i) = max(timin, min(tisfc(i), tgice))
              else
                cice(i)   = zero
                hice(i)   = zero
                islmsk(i) = 0
              endif
              islmsk_cice(i) = islmsk(i)
              if (cice(i) < one) then
                wet(i) = .true. ! some open lake
                if (icy(i)) tsfco(i) = max(tisfc(i), tgice)
              endif
            endif
          else            ! all land
            cice(i) = zero
            hice(i) = zero
            islmsk_cice(i) = 1
            islmsk(i)      = 1
          endif
        enddo  

      else

        do i = 1, IM
          if (islmsk(i) == 1) then
!           tsfcl(i) = tsfc(i)
            dry(i)    = .true.
            frland(i) = one
            cice(i)   = zero
            hice(i)   = zero
          else
            frland(i) = zero
            if (oceanfrac(i) > zero) then
              if (cice(i) >= min_seaice) then
                icy(i)   = .true.
                tisfc(i) = max(timin, min(tisfc(i), tgice))
              else
                cice(i)        = zero
                hice(i)        = zero
                flag_cice(i)   = .false.
                islmsk(i)      = 0
                islmsk_cice(i) = 0
              endif
              if (cice(i) < one) then
                wet(i) = .true. ! some open ocean
                if (.not. cplflx .and. icy(i)) tsfco(i) = max(tisfc(i), tgice)
              endif
            else
              if (cice(i) >= min_lakeice) then
                icy(i) = .true.
                tisfc(i) = max(timin, min(tisfc(i), tgice))
              else
                cice(i)   = zero
                hice(i)   = zero
                flag_cice(i) = .false.
                islmsk(i) = 0
              endif
              islmsk_cice(i) = islmsk(i)
              if (cice(i) < one) then
                wet(i) = .true. ! some open lake
                if (icy(i)) tsfco(i) = max(tisfc(i), tgice)
              endif
            endif
          endif
        enddo
      endif

      do i=1,im
        tprcp_wat(i) = tprcp(i)
        tprcp_lnd(i) = tprcp(i)
        tprcp_ice(i) = tprcp(i)
        if (wet(i)) then                   ! Water
          uustar_wat(i) = uustar(i)
            tsfc_wat(i) = tsfco(i)
           tsurf_wat(i) = tsfco(i)
!          weasd_wat(i) = weasd(i)
!          snowd_wat(i) = snowd(i)
           weasd_wat(i) = zero
           snowd_wat(i) = zero
           !-- reference emiss value for surface emissivity in setemis
           !   1-open water, 2-grass/shrub land, 3-bare soil, tundra,
           !   4-sandy desert, 5-rocky desert, 6-forest, 7-ice, 8-snow
           !data  emsref / 0.97, 0.95, 0.94, 0.90, 0.93, 0.96, 0.96, 0.99 /
          if(iemsflg == 2) then
           semis_wat(i) = 0.97_kind_phys ! consistent with setemis
          else
           semis_wat(i) = 0.984_kind_phys
          endif
             qss_wat(i) = qss(i)
            hflx_wat(i) = hflx(i)
        ! DH*
        else
          zorlo(i) = huge
        ! *DH
        endif
        if (dry(i)) then                   ! Land
          uustar_lnd(i) = uustar(i)
           weasd_lnd(i) = weasd(i)
            tsfc_lnd(i) = tsfcl(i)
           tsurf_lnd(i) = tsfcl(i)
           snowd_lnd(i) = snowd(i)
          if (iemsflg == 2 .and. .not. flag_init) then
           !-- use land emissivity from the LSM
           semis_lnd(i) = emis_lnd(i)
          else
           semis_lnd(i) = semis_rad(i)
          endif
             qss_lnd(i) = qss(i)
            hflx_lnd(i) = hflx(i)
        ! DH*
        else
          zorll(i) = huge
        ! *DH
        end if
        if (icy(i)) then                   ! Ice
          uustar_ice(i) = uustar(i)
           weasd_ice(i) = weasd(i)
            tsfc_ice(i) = tisfc(i)
           tsurf_ice(i) = tisfc(i)
           snowd_ice(i) = snowd(i)
            ep1d_ice(i) = zero
            gflx_ice(i) = zero
          if (iemsflg == 2 .and. (.not.flag_init .or. flag_restart) .and. lsm == lsm_ruc) then
           !-- use emis_ice from RUC LSM with snow effect
           semis_ice(i) = emis_ice(i)
          else
           semis_ice(i) = 0.95_kind_phys
          endif
             qss_ice(i) = qss(i)
            hflx_ice(i) = hflx(i)
        ! DH*
        else
          zorli(i) = huge
        ! *DH
        end if
        if (nint(slmsk(i)) /= 1) slmsk(i)  = islmsk(i)
      enddo

! to prepare to separate lake from ocean under water category
      do i = 1, im
        if(wet(i) .and. lkm == 1) then
           if(lakefrac(i) >= 0.15 .and. lakedepth(i) > one) then
              use_flake(i) = .true.
           else
              use_flake(i) = .false.
           endif
        else
           use_flake(i) = .false.
        endif
      enddo

     ! Assign sea ice temperature to interstitial variable
      do i = 1, im
        tice(i) = tisfc(i)
      enddo

   end subroutine GFS_surface_composites_pre_run

end module GFS_surface_composites_pre


module GFS_surface_composites_inter

   use machine, only: kind_phys

   implicit none

   private

   public GFS_surface_composites_inter_init, GFS_surface_composites_inter_finalize, GFS_surface_composites_inter_run

contains

   subroutine GFS_surface_composites_inter_init ()
   end subroutine GFS_surface_composites_inter_init

   subroutine GFS_surface_composites_inter_finalize()
   end subroutine GFS_surface_composites_inter_finalize

!> \section arg_table_GFS_surface_composites_inter_run Argument Table
!! \htmlinclude GFS_surface_composites_inter_run.html
!!
   subroutine GFS_surface_composites_inter_run (im, dry, icy, wet, semis_wat, semis_lnd, semis_ice, adjsfcdlw, &
                                                gabsbdlw_lnd, gabsbdlw_ice, gabsbdlw_wat,                      &
                                                adjsfcusw, adjsfcdsw, adjsfcnsw, errmsg, errflg)

      implicit none

      ! Interface variables
      integer,                            intent(in   ) :: im
      logical,              dimension(:), intent(in   ) :: dry, icy, wet
      real(kind=kind_phys), dimension(:), intent(in   ) :: semis_wat, semis_lnd, semis_ice, adjsfcdlw, &
                                                           adjsfcdsw, adjsfcnsw
      real(kind=kind_phys), dimension(:), intent(inout) :: gabsbdlw_lnd, gabsbdlw_ice, gabsbdlw_wat
      real(kind=kind_phys), dimension(:), intent(out)   :: adjsfcusw

      ! CCPP error handling
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      ! Local variables
      integer :: i

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      !  ---  convert lw fluxes for land/ocean/sea-ice models - requires dcyc2t3 to set adjsfcdlw
      !  note: for sw: adjsfcdsw and adjsfcnsw are zenith angle adjusted downward/net fluxes.
      !        for lw: adjsfcdlw is (sfc temp adjusted) downward fluxe with no emiss effect.
      !                adjsfculw is (sfc temp adjusted) upward fluxe including emiss effect.
      !        one needs to be aware that that the absorbed downward lw flux (used by land/ocean
      !        models as downward flux) is not the same as adjsfcdlw but a value reduced by
      !        the factor of emissivity.  however, the net effects are the same when seeing
      !        it either above the surface interface or below.
      !
      !   - flux above the interface used by atmosphere model:
      !        down: adjsfcdlw;    up: adjsfculw = sfcemis*sigma*T**4 + (1-sfcemis)*adjsfcdlw
      !        net = up - down = sfcemis * (sigma*T**4 - adjsfcdlw)
      !   - flux below the interface used by lnd/oc/ice models:
      !        down: sfcemis*adjsfcdlw;  up: sfcemis*sigma*T**4
      !        net = up - down = sfcemis * (sigma*T**4 - adjsfcdlw)
      ! surface upwelling shortwave flux at current time is in adjsfcusw

      !  --- ...  define the downward lw flux absorbed by ground
      do i=1,im
        if (dry(i)) gabsbdlw_lnd(i) = semis_lnd(i) * adjsfcdlw(i)
        if (icy(i)) gabsbdlw_ice(i) = semis_ice(i) * adjsfcdlw(i)
        if (wet(i)) gabsbdlw_wat(i) = semis_wat(i) * adjsfcdlw(i)
        adjsfcusw(i) = adjsfcdsw(i) - adjsfcnsw(i)
      enddo

   end subroutine GFS_surface_composites_inter_run

end module GFS_surface_composites_inter


module GFS_surface_composites_post

   use machine, only: kind_phys

   ! For consistent calculations of composite surface properties
   use sfc_diff, only: stability

   implicit none

   private

   public GFS_surface_composites_post_init, GFS_surface_composites_post_finalize, GFS_surface_composites_post_run

   real(kind=kind_phys), parameter :: zero = 0.0_kind_phys, one = 1.0_kind_phys, &
                                      half = 0.5_kind_phys, qmin = 1.0e-8_kind_phys

contains

   subroutine GFS_surface_composites_post_init ()
   end subroutine GFS_surface_composites_post_init

   subroutine GFS_surface_composites_post_finalize()
   end subroutine GFS_surface_composites_post_finalize

!> \section arg_table_GFS_surface_composites_post_run Argument Table
!! \htmlinclude GFS_surface_composites_post_run.html
!!
   subroutine GFS_surface_composites_post_run (                                                                                   &
      im, kice, km, rd, rvrdm1, cplflx, cplwav2atm, frac_grid, flag_cice, thsfc_loc, islmsk, dry, wet, icy, wind, t1, q1, prsl1,  &
      landfrac, lakefrac, oceanfrac, zorl, zorlo, zorll, zorli,                                                                   &
      cd, cd_wat, cd_lnd, cd_ice, cdq, cdq_wat, cdq_lnd, cdq_ice, rb, rb_wat, rb_lnd, rb_ice, stress, stress_wat, stress_lnd,     &
      stress_ice, ffmm, ffmm_wat, ffmm_lnd, ffmm_ice, ffhh, ffhh_wat, ffhh_lnd, ffhh_ice, uustar, uustar_wat, uustar_lnd,         &
      uustar_ice, fm10, fm10_wat, fm10_lnd, fm10_ice, fh2, fh2_wat, fh2_lnd, fh2_ice, tsurf_wat, tsurf_lnd, tsurf_ice,            &
      cmm, cmm_wat, cmm_lnd, cmm_ice, chh, chh_wat, chh_lnd, chh_ice, gflx, gflx_wat, gflx_lnd, gflx_ice, ep1d, ep1d_wat,         &
      ep1d_lnd, ep1d_ice, weasd, weasd_wat, weasd_lnd, weasd_ice, snowd, snowd_wat, snowd_lnd, snowd_ice, tprcp, tprcp_wat,       &
      tprcp_lnd, tprcp_ice, evap, evap_wat, evap_lnd, evap_ice, hflx, hflx_wat, hflx_lnd, hflx_ice, qss, qss_wat, qss_lnd,        &
      qss_ice, tsfc, tsfco, tsfcl, tsfc_wat, tsfc_lnd, tsfc_ice, tisfc, tice, hice, cice, min_seaice, tiice, stc,                 &
      grav, prsik1, prslk1, prslki, z1, ztmax_wat, ztmax_lnd, ztmax_ice, errmsg, errflg)

      implicit none

      integer,                              intent(in) :: im, kice, km
      logical,                              intent(in) :: cplflx, frac_grid, cplwav2atm
      logical, dimension(:),                intent(in) :: flag_cice, dry, wet, icy
      integer, dimension(:),                intent(in) :: islmsk
      real(kind=kind_phys), dimension(:),   intent(in) :: wind, t1, q1, prsl1, landfrac, lakefrac, oceanfrac,                   &
        cd_wat, cd_lnd, cd_ice, cdq_wat, cdq_lnd, cdq_ice, rb_wat, rb_lnd, rb_ice, stress_wat,                                  &
        stress_lnd, stress_ice, ffmm_wat, ffmm_lnd, ffmm_ice, ffhh_wat, ffhh_lnd, ffhh_ice, uustar_wat, uustar_lnd, uustar_ice, &
        fm10_wat, fm10_lnd, fm10_ice, fh2_wat, fh2_lnd, fh2_ice, tsurf_wat, tsurf_lnd, tsurf_ice, cmm_wat, cmm_lnd, cmm_ice,    &
        chh_wat, chh_lnd, chh_ice, gflx_wat, gflx_lnd, gflx_ice, ep1d_wat, ep1d_lnd, ep1d_ice, weasd_wat, weasd_lnd, weasd_ice, &
        snowd_wat, snowd_lnd, snowd_ice,tprcp_wat, tprcp_lnd, tprcp_ice, evap_wat, evap_lnd, evap_ice, hflx_wat, hflx_lnd,      &
        hflx_ice, qss_wat, qss_lnd, qss_ice, tsfc_wat, tsfc_lnd, tsfc_ice, zorlo, zorll, zorli

      real(kind=kind_phys), dimension(:),   intent(inout) :: zorl, cd, cdq, rb, stress, ffmm, ffhh, uustar, fm10,               &
        fh2, cmm, chh, gflx, ep1d, weasd, snowd, tprcp, evap, hflx, qss, tsfc, tsfco, tsfcl, tisfc

      real(kind=kind_phys), dimension(:),   intent(in   ) :: tice ! interstitial sea ice temperature
      real(kind=kind_phys), dimension(:),   intent(inout) :: hice, cice
      real(kind=kind_phys),                 intent(in   ) :: min_seaice
      real(kind=kind_phys),                 intent(in   ) :: rd, rvrdm1

      real(kind=kind_phys), dimension(:,:), intent(in   ) :: tiice
      real(kind=kind_phys), dimension(:,:), intent(inout) :: stc

      ! Additional data needed for calling "stability"
      logical,                            intent(in   ) :: thsfc_loc
      real(kind=kind_phys),               intent(in   ) :: grav
      real(kind=kind_phys), dimension(:), intent(in   ) :: prsik1, prslk1, prslki, z1
      real(kind=kind_phys), dimension(:), intent(in   ) :: ztmax_wat, ztmax_lnd, ztmax_ice

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      ! Local variables
      integer :: i, k
      real(kind=kind_phys) :: txl, txi, txo, wfrac, q0, rho
      ! For calling "stability"
      real(kind=kind_phys) :: tsurf, virtfac, tv1, thv1, tvs, z0max, ztmax

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      ! --- generate ocean/land/ice composites

      if (frac_grid) then

        do i=1, im

          ! Three-way composites (fields from sfc_diff)
          txl   = landfrac(i)            ! land fraction
          wfrac = one - txl              ! ocean fraction
          txi   = cice(i) * wfrac        ! txi = ice fraction wrt whole cell
          txo   = max(zero, wfrac-txi)   ! txo = open water fraction

         !gflx(i)   = txl*gflx_lnd(i)   + txi*gflx_ice(i)   + txo*gflx_wat(i)
          ep1d(i)   = txl*ep1d_lnd(i)   + txi*ep1d_ice(i)   + txo*ep1d_wat(i)
         !weasd(i)  = txl*weasd_lnd(i)  + txi*weasd_ice(i)  + txo*weasd_wat(i)
         !snowd(i)  = txl*snowd_lnd(i)  + txi*snowd_ice(i)  + txo*snowd_wat(i)
          weasd(i)  = txl*weasd_lnd(i)  + txi*weasd_ice(i)
          snowd(i)  = txl*snowd_lnd(i)  + txi*snowd_ice(i)
         !tprcp(i)  = txl*tprcp_lnd(i)  + txi*tprcp_ice(i)  + txo*tprcp_wat(i)

          if (.not. flag_cice(i) .and. islmsk(i) == 2) then
            evap(i) = txl*evap_lnd(i)   + wfrac*evap_ice(i)
            hflx(i) = txl*hflx_lnd(i)   + wfrac*hflx_ice(i)
            qss(i)  = txl*qss_lnd(i)    + wfrac*qss_ice(i)
            gflx(i) = txl*gflx_lnd(i)   + wfrac*gflx_ice(i)
          else
            evap(i) = txl*evap_lnd(i)   + txi*evap_ice(i)   + txo*evap_wat(i)
            hflx(i) = txl*hflx_lnd(i)   + txi*hflx_ice(i)   + txo*hflx_wat(i)
            qss(i)  = txl*qss_lnd(i)    + txi*qss_ice(i)    + txo*qss_wat(i)
            gflx(i) = txl*gflx_lnd(i)   + txi*gflx_ice(i)   + txo*gflx_wat(i)
          endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         Call stability for consistent surface properties. Currently this comes from !
!         the GFS surface layere scheme (sfc_diff), regardless of the actual surface  !
!         layer parameterization being used - to be extended in the future            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          tsfc(i)   = ( txl * cdq_lnd(i) * tsfc_lnd(i)                          &
                      + txi * cdq_ice(i) * tice(i)                              & ! DH* Ben had tsurf_ice(i), but GFS_surface_composites_post_run uses tice instead
                      + txo * cdq_wat(i) * tsfc_wat(i))                         &
                      / (txl * cdq_lnd(i) + txi * cdq_ice(i) + txo * cdq_wat(i) )
          tsurf     = ( txl * cdq_lnd(i) * tsurf_lnd(i)                         &
                      + txi * cdq_ice(i) * tsurf_ice(i)                         &
                      + txo * cdq_wat(i) * tsurf_wat(i))                        &
                      / (txl * cdq_lnd(i) + txi * cdq_ice(i) + txo * cdq_wat(i) )

          q0 = max( q1(i), qmin )
          virtfac = one + rvrdm1 * q0

          tv1 = t1(i) * virtfac ! Virtual temperature in middle of lowest layer
          if(thsfc_loc) then ! Use local potential temperature
            thv1 = t1(i) * prslki(i) * virtfac  ! Theta-v at lowest level
            tvs  = half * (tsfc(i)+tsurf) * virtfac
          else ! Use potential temperature referenced to 1000 hPa
            thv1 = t1(i) / prslk1(i) * virtfac  ! Theta-v at lowest level
            tvs  = half * (tsfc(i)+tsurf)/prsik1(i) * virtfac
          endif

          zorl(i) = exp(txl*log(zorll(i)) + txi*log(zorli(i)) + txo*log(zorlo(i)))
          z0max   = 0.01_kind_phys * zorl(i)
          ztmax   = exp(txl*log(ztmax_lnd(i)) + txi*log(ztmax_ice(i)) + txo*log(ztmax_wat(i)))

          ! Only actually need to call "stability" if multiple surface types exist...
          if(txl .eq. one) then ! 100% land
            rb(i)     = rb_lnd(i)
            ffmm(i)   = ffmm_lnd(i)
            ffhh(i)   = ffhh_lnd(i)
            fm10(i)   = fm10_lnd(i)
            fh2(i)    = fh2_lnd(i)
            cd(i)     = cd_lnd(i)
            cdq(i)    = cdq_lnd(i)
            stress(i) = stress_lnd(i)
            uustar(i) = uustar_lnd(i)
          elseif(txo .eq. one) then ! 100% open water
            rb(i)     = rb_wat(i)
            ffmm(i)   = ffmm_wat(i)
            ffhh(i)   = ffhh_wat(i)
            fm10(i)   = fm10_wat(i)
            fh2(i)    = fh2_wat(i)
            cd(i)     = cd_wat(i)
            cdq(i)    = cdq_wat(i)
            stress(i) = stress_wat(i)
            uustar(i) = uustar_wat(i)
          elseif(txi .eq. one) then ! 100% ice
            rb(i)     = rb_ice(i)
            ffmm(i)   = ffmm_ice(i)
            ffhh(i)   = ffhh_ice(i)
            fm10(i)   = fm10_ice(i)
            fh2(i)    = fh2_ice(i)
            cd(i)     = cd_ice(i)
            cdq(i)    = cdq_ice(i)
            stress(i) = stress_ice(i)
            uustar(i) = uustar_ice(i)
          else ! Mix of multiple surface types (land, water, and/or ice)
            call stability(z1(i), snowd(i), thv1, wind(i), z0max, ztmax, tvs, grav, & ! inputs
                           tv1, thsfc_loc,                                          & ! inputs
                           rb(i), ffmm(i), ffhh(i), fm10(i), fh2(i), cd(i), cdq(i), & ! outputs
                           stress(i), uustar(i))
          endif ! Checking to see if point is one or multiple surface types

          ! BWG, 2021/02/25: cmm=cd*wind, chh=cdq*wind, so use composite cd, cdq
          rho       = prsl1(i) / (rd*tv1)
          cmm(i)    =      cd(i)*wind(i)  !txl*cmm_lnd(i)    + txi*cmm_ice(i)    + txo*cmm_wat(i)
          chh(i)    = rho*cdq(i)*wind(i)  !txl*chh_lnd(i)    + txi*chh_ice(i)    + txo*chh_wat(i)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          if (dry(i)) then
            tsfcl(i) = tsfc_lnd(i)       ! over land
          elseif (wet(i)) then
            tsfcl(i) = tsfc_wat(i)       ! over water
          else
            tsfcl(i) = tice(i)           ! over ice
          endif
          if (wet(i)) then
            tsfco(i) = tsfc_wat(i)       ! over lake or ocean when uncoupled
          elseif (icy(i)) then
            tsfco(i) = tice(i)           ! over lake or ocean ice when uncoupled
          else
            tsfco(i) = tsfc_lnd(i)       ! over land
          endif
          if (icy(i)) then
            tisfc(i) = tice(i)           ! over lake or ocean ice when uncoupled
          elseif (wet(i)) then
            tisfc(i) = tsfc_wat(i)       ! over lake or ocean when uncoupled
          else
            tisfc(i) = tsfc_lnd(i)       ! over land
          endif
                                                  ! for coupled model ocean will replace this
!         if (icy(i)) tisfc(i) = tsfc_ice(i)      ! over ice when uncoupled
!         if (icy(i)) tisfc(i) = tice(i)          ! over ice when uncoupled

!         if (wet(i) .and. .not. cplflx) then
!           tsfco(i) = tsfc_wat(i)                ! over lake or ocean when uncoupled
!           tisfc(i) = tsfc_ice(i)                ! over ice when uncoupled
!         endif

!         if (.not. flag_cice(i)) then
!           if (islmsk(i) == 2) then              ! return updated lake ice thickness & concentration to global array
!             tisfc(i) = tice(i)
!           else                                  ! this would be over open ocean or land (no ice fraction)
!             hice(i)  = zero
!             cice(i)  = zero
!             tisfc(i) = tsfc(i)
!           endif
!         endif
          if (.not. icy(i)) then
            hice(i)  = zero
            cice(i)  = zero
          endif
        enddo

      else

        do i=1,im
          if (islmsk(i) == 1) then
            zorl(i)   = zorll(i)
            cd(i)     = cd_lnd(i)
            cdq(i)    = cdq_lnd(i)
            rb(i)     = rb_lnd(i)
            stress(i) = stress_lnd(i)
            ffmm(i)   = ffmm_lnd(i)
            ffhh(i)   = ffhh_lnd(i)
            uustar(i) = uustar_lnd(i)
            fm10(i)   = fm10_lnd(i)
            fh2(i)    = fh2_lnd(i)
            tsfcl(i)  = tsfc_lnd(i) ! over land
            tsfc(i)   = tsfcl(i)
            tsfco(i)  = tsfc(i)
            tisfc(i)  = tsfc(i)
            cmm(i)    = cmm_lnd(i)
            chh(i)    = chh_lnd(i)
            gflx(i)   = gflx_lnd(i)
            ep1d(i)   = ep1d_lnd(i)
            weasd(i)  = weasd_lnd(i)
            snowd(i)  = snowd_lnd(i)
           !tprcp(i)  = tprcp_lnd(i)
            evap(i)   = evap_lnd(i)
            hflx(i)   = hflx_lnd(i)
            qss(i)    = qss_lnd(i)
            hice(i)   = zero
            cice(i)   = zero
          elseif (islmsk(i) == 0) then
            zorl(i)   = zorlo(i)
            cd(i)     = cd_wat(i)
            cdq(i)    = cdq_wat(i)
            rb(i)     = rb_wat(i)
            stress(i) = stress_wat(i)
            ffmm(i)   = ffmm_wat(i)
            ffhh(i)   = ffhh_wat(i)
            uustar(i) = uustar_wat(i)
            fm10(i)   = fm10_wat(i)
            fh2(i)    = fh2_wat(i)
            tsfco(i)  = tsfc_wat(i) ! over lake (and ocean when uncoupled)
            tsfc(i)   = tsfco(i)
            tsfcl(i)  = tsfc(i)
            tisfc(i)  = tsfc(i)
            cmm(i)    = cmm_wat(i)
            chh(i)    = chh_wat(i)
            gflx(i)   = gflx_wat(i)
            ep1d(i)   = ep1d_wat(i)
            weasd(i)  = weasd_wat(i)
            snowd(i)  = snowd_wat(i)
           !tprcp(i)  = tprcp_wat(i)
            evap(i)   = evap_wat(i)
            hflx(i)   = hflx_wat(i)
            qss(i)    = qss_wat(i)
            hice(i)   = zero
            cice(i)   = zero
          else ! islmsk(i) == 2
            zorl(i)   = zorli(i)
            cd(i)     = cd_ice(i)
            cdq(i)    = cdq_ice(i)
            rb(i)     = rb_ice(i)
            ffmm(i)   = ffmm_ice(i)
            ffhh(i)   = ffhh_ice(i)
            uustar(i) = uustar_ice(i)
            fm10(i)   = fm10_ice(i)
            fh2(i)    = fh2_ice(i)
            stress(i) = stress_ice(i)
            cmm(i)    = cmm_ice(i)
            chh(i)    = chh_ice(i)
            gflx(i)   = gflx_ice(i)
            ep1d(i)   = ep1d_ice(i)
            weasd(i)  = weasd_ice(i)
            snowd(i)  = snowd_ice(i)
           !tprcp(i)  = cice(i)*tprcp_ice(i) + (one-cice(i))*tprcp_wat(i)
            qss(i)    = qss_ice(i)
            tsfc(i)   = tsfc_ice(i)
            evap(i)   = evap_ice(i)
            hflx(i)   = hflx_ice(i)
            qss(i)    = qss_ice(i)
            tisfc(i)  = tice(i)
            if (.not. flag_cice(i)) then
!             tisfc(i) = tice(i) ! over lake ice (and sea ice when uncoupled)
              zorl(i)  = cice(i) * zorli(i)   + (one - cice(i)) * zorlo(i)
              tsfc(i)  = tsfc_ice(i) ! over lake (and ocean when uncoupled)
            elseif (wet(i)) then
              if (cice(i) >= min_seaice) then ! this was already done for lake ice in sfc_sice
                txi = cice(i)
                txo = one - txi
                evap(i)   = txi * evap_ice(i)   + txo * evap_wat(i)
                hflx(i)   = txi * hflx_ice(i)   + txo * hflx_wat(i)
                tsfc(i)   = txi * tsfc_ice(i)   + txo * tsfc_wat(i)
                stress(i) = txi * stress_ice(i) + txo * stress_wat(i)
                qss(i)    = txi * qss_ice(i)    + txo * qss_wat(i)
                ep1d(i)   = txi * ep1d_ice(i)   + txo * ep1d_wat(i)
                zorl(i)   = txi * zorli(i)      + txo * zorlo(i)
              else
                evap(i)   = evap_wat(i)
                hflx(i)   = hflx_wat(i)
                tsfc(i)   = tsfc_wat(i)
                stress(i) = stress_wat(i)
                qss(i)    = qss_wat(i)
                ep1d(i)   = ep1d_wat(i)
                zorl(i)   = zorlo(i)
              endif
            endif
            if (wet(i)) then
              tsfco(i) = tsfc_wat(i)
            else
              tsfco(i) = tsfc(i)
            endif
            tsfcl(i)  = tsfc(i)
            do k=1,min(kice,km) ! store tiice in stc to reduce output in the nonfrac grid case
              stc(i,k) = tiice(i,k)
            end do
          endif

        enddo

      endif ! if (frac_grid)

      ! --- compositing done

   end subroutine GFS_surface_composites_post_run

end module GFS_surface_composites_post

!> \file GFS_surface_composites.F90
!!  Contains code related to generating composites for all GFS surface schemes.

module GFS_surface_composites_pre

   use machine, only: kind_phys

   implicit none

   private

   public GFS_surface_composites_pre_init, GFS_surface_composites_pre_finalize, GFS_surface_composites_pre_run

contains

   subroutine GFS_surface_composites_pre_init ()
   end subroutine GFS_surface_composites_pre_init

   subroutine GFS_surface_composites_pre_finalize()
   end subroutine GFS_surface_composites_pre_finalize

#if 0
!> \section arg_table_GFS_surface_composites_pre_run Argument Table
!! \htmlinclude GFS_surface_composites_pre_run.html
!!
#endif
   subroutine GFS_surface_composites_pre_run (im, frac_grid, flag_cice, cplflx, landfrac, lakefrac, oceanfrac,  &
                                 frland, dry, icy, lake, ocean, wet, cice, cimin, zorl, zorlo, zorll, zorl_ocn, &
                                 zorl_lnd, zorl_ice, snowd, snowd_ocn, snowd_lnd, snowd_ice, tprcp, tprcp_ocn,  &
                                 tprcp_lnd, tprcp_ice, uustar, uustar_lnd, uustar_ice, weasd, weasd_ocn,        &
                                 weasd_lnd, weasd_ice, ep1d_ice, tsfc, tsfco, tsfcl, tsfc_ocn, tsfc_lnd,        &
                                 tsfc_ice, tisfc, tice, tsurf, tsurf_ocn, tsurf_lnd, tsurf_ice, gflx_ice,       &
                                 errmsg, errflg)

      use machine, only: kind_phys

      implicit none

      ! Interface variables
      integer,                             intent(in   ) :: im
      logical,                             intent(in   ) :: frac_grid, cplflx
      logical, dimension(im),              intent(in   ) :: flag_cice
      logical,              dimension(im), intent(inout) :: dry, icy, lake, ocean, wet
      real(kind=kind_phys),                intent(in   ) :: cimin
      real(kind=kind_phys), dimension(im), intent(in   ) :: landfrac, lakefrac, oceanfrac, cice
      real(kind=kind_phys), dimension(im), intent(  out) :: frland
      real(kind=kind_phys), dimension(im), intent(in   ) :: zorl, snowd, tprcp, uustar, weasd

      real(kind=kind_phys), dimension(im), intent(inout) :: zorlo, zorll, tsfc, tsfco, tsfcl, tisfc, tsurf
      real(kind=kind_phys), dimension(im), intent(inout) :: snowd_ocn, snowd_lnd, snowd_ice, tprcp_ocn, &
        tprcp_lnd, tprcp_ice, zorl_ocn, zorl_lnd, zorl_ice, tsfc_ocn, tsfc_lnd, tsfc_ice, tsurf_ocn,    &
        tsurf_lnd, tsurf_ice, uustar_lnd, uustar_ice, weasd_ocn, weasd_lnd, weasd_ice, ep1d_ice, gflx_ice
      real(kind=kind_phys), dimension(im), intent(  out) :: tice

      ! CCPP error handling
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      ! Local variables
      real(kind=kind_phys), parameter :: one = 1.0d0
      integer :: i

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      do i=1,im
        frland(i) = landfrac(i)
        if (frland(i) > 0.0)                                    dry(i) = .true.
        if (cice(i) >= cimin*(1.-frland(i)) .and. frland(i)<1.) icy(i) = .true.
        if (frland(i)+cice(i) < 1.0 )                           wet(i) = .true. ! there is some open water!
      enddo  

      if (frac_grid) then
        do i=1,im
          tsfc(i) = tsfcl(i) *               frland(i)   &
                  + tisfc(i) *      cice(i)              &
                  + tsfco(i) * (one-cice(i)-frland(i))
        enddo
      elseif (cplflx) then
        do i=1,im
          if (flag_cice(i)) then
            tsfc(i) = tisfc(i) *      cice(i)            &
                    + tsfc (i) * (one-cice(i))
            icy(i) = .true.
          endif
        enddo
      endif

      if (.not. cplflx .or. .not. frac_grid) then
        do i=1,im
          zorll(i) = zorl(i)
          zorlo(i) = zorl(i)
          tsfcl(i) = tsfc(i)
          tsfco(i) = tsfc(i)
          !tisfc(i) = tsfc(i)
        enddo
      endif

      do i=1,im
        if (wet(i)) then                   ! Water
           tprcp_ocn(i) = tprcp(i)
            zorl_ocn(i) = zorlo(i)
            tsfc_ocn(i) = tsfco(i)
           tsurf_ocn(i) = tsfco(i)
!           weasd_ocn(i) = weasd(i)
!           snowd_ocn(i) = snowd(i)
           weasd_ocn(i) = 0.0
           snowd_ocn(i) = 0.0
        endif
        if (dry(i)) then                   ! Land
          uustar_lnd(i) = uustar(i)
           weasd_lnd(i) = weasd(i)
           tprcp_lnd(i) = tprcp(i)
            zorl_lnd(i) = zorll(i)
            tsfc_lnd(i) = tsfcl(i)
           tsurf_lnd(i) = tsfcl(i)
           snowd_lnd(i) = snowd(i)
        end if
        if (icy(i)) then                   ! Ice
          uustar_ice(i) = uustar(i)
           weasd_ice(i) = weasd(i)
           tprcp_ice(i) = tprcp(i)
            zorl_ice(i) = zorll(i)
!            tsfc_ice(i) = tisfc(i)
!           tsurf_ice(i) = tisfc(i)
            tsfc_ice(i) = tsfc(i)
           tsurf_ice(i) = tsfc(i)
           snowd_ice(i) = snowd(i)
            ep1d_ice(i) = 0.
            gflx_ice(i) = 0.
        end if
      enddo

     ! Assign sea ice temperature to interstitial variable
      do i = 1, im
        tice(i) = tisfc(i)
      end do

   end subroutine GFS_surface_composites_pre_run

end module GFS_surface_composites_pre


module GFS_surface_composites_post

   use machine, only: kind_phys

   implicit none

   private

   public GFS_surface_composites_post_init, GFS_surface_composites_post_finalize, GFS_surface_composites_post_run

contains

   subroutine GFS_surface_composites_post_init ()
   end subroutine GFS_surface_composites_post_init

   subroutine GFS_surface_composites_post_finalize()
   end subroutine GFS_surface_composites_post_finalize

#if 0
!> \section arg_table_GFS_surface_composites_post_run Argument Table
!! \htmlinclude GFS_surface_composites_post_run.html
!!
#endif
   subroutine GFS_surface_composites_post_run (                                                                                   &
      im, cplflx, frac_grid, flag_cice, islmsk, dry, wet, icy, landfrac, zorl, zorlo, zorll, zorl_ocn, zorl_lnd, zorl_ice,        &
      cd, cd_ocn, cd_lnd, cd_ice, cdq, cdq_ocn, cdq_lnd, cdq_ice, rb, rb_ocn, rb_lnd, rb_ice, stress, stress_ocn, stress_lnd,     &
      stress_ice, ffmm, ffmm_ocn, ffmm_lnd, ffmm_ice, ffhh, ffhh_ocn, ffhh_lnd, ffhh_ice, uustar, uustar_ocn, uustar_lnd,         &
      uustar_ice, fm10, fm10_ocn, fm10_lnd, fm10_ice, fh2, fh2_ocn, fh2_lnd, fh2_ice, tsurf, tsurf_ocn, tsurf_lnd, tsurf_ice,     &
      cmm, cmm_ocn, cmm_lnd, cmm_ice, chh, chh_ocn, chh_lnd, chh_ice, gflx, gflx_ocn, gflx_lnd, gflx_ice, ep1d, ep1d_ocn,         &
      ep1d_lnd, ep1d_ice, weasd, weasd_ocn, weasd_lnd, weasd_ice, snowd, snowd_ocn, snowd_lnd, snowd_ice, tprcp, tprcp_ocn,       &
      tprcp_lnd, tprcp_ice, evap, evap_ocn, evap_lnd, evap_ice, hflx, hflx_ocn, hflx_lnd, hflx_ice, qss, qss_ocn, qss_lnd,        &
      qss_ice, tsfc, tsfco, tsfcl, tsfc_ocn, tsfc_lnd, tsfc_ice, tisfc, tice, hice, cice, errmsg, errflg)

      use machine, only: kind_phys

      implicit none

      integer,                              intent(in) :: im
      logical,                              intent(in) :: cplflx, frac_grid
      logical, dimension(im),               intent(in) :: flag_cice, dry, wet, icy
      integer, dimension(im),               intent(in) :: islmsk
      real(kind=kind_phys), dimension(im),  intent(in) :: landfrac,                                                             &
        zorl_ocn, zorl_lnd, zorl_ice, cd_ocn, cd_lnd, cd_ice, cdq_ocn, cdq_lnd, cdq_ice, rb_ocn, rb_lnd, rb_ice, stress_ocn,    &
        stress_lnd, stress_ice, ffmm_ocn, ffmm_lnd, ffmm_ice, ffhh_ocn, ffhh_lnd, ffhh_ice, uustar_ocn, uustar_lnd, uustar_ice, &
        fm10_ocn, fm10_lnd, fm10_ice, fh2_ocn, fh2_lnd, fh2_ice, tsurf_ocn, tsurf_lnd, tsurf_ice, cmm_ocn, cmm_lnd, cmm_ice,    &
        chh_ocn, chh_lnd, chh_ice, gflx_ocn, gflx_lnd, gflx_ice, ep1d_ocn, ep1d_lnd, ep1d_ice, weasd_ocn, weasd_lnd, weasd_ice, &
        snowd_ocn, snowd_lnd, snowd_ice,tprcp_ocn, tprcp_lnd, tprcp_ice, evap_ocn, evap_lnd, evap_ice, hflx_ocn, hflx_lnd,      &
        hflx_ice, qss_ocn, qss_lnd, qss_ice, tsfc_ocn, tsfc_lnd, tsfc_ice

      real(kind=kind_phys), dimension(im),  intent(inout) :: zorl, zorlo, zorll, cd, cdq, rb, stress, ffmm, ffhh, uustar, fm10, &
        fh2, tsurf, cmm, chh, gflx, ep1d, weasd, snowd, tprcp, evap, hflx, qss, tsfc, tsfco, tsfcl, tisfc

      real(kind=kind_phys), dimension(im),  intent(in   ) :: tice ! interstitial sea ice temperature
      real(kind=kind_phys), dimension(im),  intent(inout) :: hice, cice

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      ! Local variables
      integer :: i
      real(kind=kind_phys) :: txl, txi, txo

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      ! --- generate ocean/land/ice composites

      if (frac_grid) then

        do i=1, im

          ! Three-way composites (fields from sfc_diff)
          txl = landfrac(i)
          txi = cice(i)                 ! here cice is grid fraction that is ice
          txo = 1.0 - txl - txi

          zorl(i)   = txl*zorl_lnd(i)   + txi*zorl_ice(i)   + txo*zorl_ocn(i)
          cd(i)     = txl*cd_lnd(i)     + txi*cd_ice(i)     + txo*cd_ocn(i)
          cdq(i)    = txl*cdq_lnd(i)    + txi*cdq_ice(i)    + txo*cdq_ocn(i)
          rb(i)     = txl*rb_lnd(i)     + txi*rb_ice(i)     + txo*rb_ocn(i)
          stress(i) = txl*stress_lnd(i) + txi*stress_ice(i) + txo*stress_ocn(i)
          ffmm(i)   = txl*ffmm_lnd(i)   + txi*ffmm_ice(i)   + txo*ffmm_ocn(i)
          ffhh(i)   = txl*ffhh_lnd(i)   + txi*ffhh_ice(i)   + txo*ffhh_ocn(i)
          uustar(i) = txl*uustar_lnd(i) + txi*uustar_ice(i) + txo*uustar_ocn(i)
          fm10(i)   = txl*fm10_lnd(i)   + txi*fm10_ice(i)   + txo*fm10_ocn(i)
          fh2(i)    = txl*fh2_lnd(i)    + txi*fh2_ice(i)    + txo*fh2_ocn(i)
          !tsurf(i)  = txl*tsurf_lnd(i)  + txi*tice(i)      + txo*tsurf_ocn(i)
          !tsurf(i)  = txl*tsurf_lnd(i)  + txi*tsurf_ice(i)  + txo*tsurf_ocn(i) ! not used again! Moorthi
          cmm(i)    = txl*cmm_lnd(i)    + txi*cmm_ice(i)    + txo*cmm_ocn(i)
          chh(i)    = txl*chh_lnd(i)    + txi*chh_ice(i)    + txo*chh_ocn(i)
          gflx(i)   = txl*gflx_lnd(i)   + txi*gflx_ice(i)   + txo*gflx_ocn(i)
          ep1d(i)   = txl*ep1d_lnd(i)   + txi*ep1d_ice(i)   + txo*ep1d_ocn(i)
          !weasd(i)  = txl*weasd_lnd(i)  + txi*weasd_ice(i)  + txo*weasd_ocn(i)
          !snowd(i)  = txl*snowd_lnd(i)  + txi*snowd_ice(i)  + txo*snowd_ocn(i)
          weasd(i)  = txl*weasd_lnd(i)  + txi*weasd_ice(i)
          snowd(i)  = txl*snowd_lnd(i)  + txi*snowd_ice(i)
          tprcp(i)  = txl*tprcp_lnd(i)  + txi*tprcp_ice(i)  + txo*tprcp_ocn(i)
          evap(i)   = txl*evap_lnd(i)   + txi*evap_ice(i)   + txo*evap_ocn(i)
          hflx(i)   = txl*hflx_lnd(i)   + txi*hflx_ice(i)   + txo*hflx_ocn(i)
          qss(i)    = txl*qss_lnd(i)    + txi*qss_ice(i)    + txo*qss_ocn(i)
          tsfc(i)   = txl*tsfc_lnd(i)   + txi*tice(i)       + txo*tsfc_ocn(i)
          !tsfc(i)   = txl*tsfc_lnd(i)   + txi*tsfc_ice(i)   + txo*tsfc_ocn(i)

          zorll(i) = zorl_lnd(i)
          zorlo(i) = zorl_ocn(i)
  
          if (dry(i)) tsfcl(i) = tsfc_lnd(i)      ! over land
          if (wet(i)) tsfco(i) = tsfc_ocn(i)      ! over lake or ocean when uncoupled
          tisfc(i) = tsfc(i)                      ! assume bitwise identical on non-icy points
          if (icy(i)) then
            tisfc(i) = tsfc_ice(i)                ! over ice when uncoupled
!           tisfc(i) = tice(i)                    ! over ice when uncoupled
          else
            hice(i) = 0.0
            cice(i) = 0.0
          end if

!         if (wet(i) .and. .not. cplflx) then
!           tsfco(i) = tsfc3_ocn(i)               ! over lake or ocean when uncoupled
!           tisfc(i) = tsfc3_ice(i)               ! over ice when uncoupled
!         endif

        end do

      else

        do i=1,im
          if (islmsk(i) == 1) then
            zorl(i)   = zorl_lnd(i)
            cd(i)     = cd_lnd(i)
            cdq(i)    = cdq_lnd(i)
            rb(i)     = rb_lnd(i)
            stress(i) = stress_lnd(i)
            ffmm(i)   = ffmm_lnd(i)
            ffhh(i)   = ffhh_lnd(i)
            uustar(i) = uustar_lnd(i)
            fm10(i)   = fm10_lnd(i)
            fh2(i)    = fh2_lnd(i)
            !tsurf(i) = tsurf_lnd(i)
            cmm(i)    = cmm_lnd(i)
            chh(i)    = chh_lnd(i)
            gflx(i)   = gflx_lnd(i)
            ep1d(i)   = ep1d_lnd(i)
            weasd(i)  = weasd_lnd(i)
            snowd(i)  = snowd_lnd(i)
            tprcp(i)  = tprcp_lnd(i)
            evap(i)   = evap_lnd(i)
            hflx(i)   = hflx_lnd(i)
            qss(i)    = qss_lnd(i)
            tsfc(i)   = tsfc_lnd(i)
            cmm(i)    = cmm_lnd(i)
            chh(i)    = chh_lnd(i)
          elseif (islmsk(i) == 0) then
            zorl(i)   = zorl_ocn(i)
            cd(i)     = cd_ocn(i)
            cdq(i)    = cdq_ocn(i)
            rb(i)     = rb_ocn(i)
            stress(i) = stress_ocn(i)
            ffmm(i)   = ffmm_ocn(i)
            ffhh(i)   = ffhh_ocn(i)
            uustar(i) = uustar_ocn(i)
            fm10(i)   = fm10_ocn(i)
            fh2(i)    = fh2_ocn(i)
            !tsurf(i) = tsurf_ocn(i)
            cmm(i)    = cmm_ocn(i)
            chh(i)    = chh_ocn(i)
            gflx(i)   = gflx_ocn(i)
            ep1d(i)   = ep1d_ocn(i)
            weasd(i)  = weasd_ocn(i)
            snowd(i)  = snowd_ocn(i)
            tprcp(i)  = tprcp_ocn(i)
            evap(i)   = evap_ocn(i)
            hflx(i)   = hflx_ocn(i)
            qss(i)    = qss_ocn(i)
            tsfc(i)   = tsfc_ocn(i)
            cmm(i)    = cmm_ocn(i)
            chh(i)    = chh_ocn(i)
          else
            zorl(i)   = zorl_ice(i)
            cd(i)     = cd_ice(i)
            cdq(i)    = cdq_ice(i)
            rb(i)     = rb_ice(i)
            stress(i) = stress_ice(i)
            ffmm(i)   = ffmm_ice(i)
            ffhh(i)   = ffhh_ice(i)
            uustar(i) = uustar_ice(i)
            fm10(i)   = fm10_ice(i)
            fh2(i)    = fh2_ice(i)
            !tsurf(i) = tsurf_ice(i)
            cmm(i)    = cmm_ice(i)
            chh(i)    = chh_ice(i)
            gflx(i)   = gflx_ice(i)
            ep1d(i)   = ep1d_ice(i)
            weasd(i)  = weasd_ice(i)
            snowd(i)  = snowd_ice(i)
            tprcp(i)  = tprcp_ice(i)
            evap(i)   = evap_ice(i)
            hflx(i)   = hflx_ice(i)
            qss(i)    = qss_ice(i)
            tsfc(i)   = tsfc_ice(i)
            cmm(i)    = cmm_ice(i)
            chh(i)    = chh_ice(i)
          endif

          zorll(i) = zorl_lnd(i)
          zorlo(i) = zorl_ocn(i)

          if (flag_cice(i)) then
            evap(i) = cice(i) * evap_ice(i) + (1.0-cice(i)) * evap_ocn(i)
            hflx(i) = cice(i) * hflx_ice(i) + (1.0-cice(i)) * hflx_ocn(i)
            tsfc(i) = cice(i) * tsfc_ice(i) + (1.0-cice(i)) * tsfc_ocn(i)
          endif

          if (dry(i)) tsfcl(i) = tsfc_lnd(i)      ! over land
          if (wet(i)) tsfco(i) = tsfc_ocn(i)      ! over lake or ocean when uncoupled
          tisfc(i) = tsfc(i)                      ! assume bitwise identical on non-icy points
          if (icy(i)) then
!           tisfc(i) = tsfc_ice(i)                ! over ice when uncoupled
            tisfc(i) = tice(i)                    ! over ice when uncoupled
          else
            hice(i) = 0.0
            cice(i) = 0.0
          end if

!         if (wet(i) .and. .not. cplflx) then
!           tsfco(i) = tsfc_ocn(i)                ! over lake or ocean when uncoupled
!           tisfc(i) = tsfc_ice(i)                ! over ice when uncoupled
!         endif

        end do

      end if ! if (frac_grid)

      ! --- compositing done

   end subroutine GFS_surface_composites_post_run

end module GFS_surface_composites_post

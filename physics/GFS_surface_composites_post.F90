!> \file GFS_surface_composites_post.F90
!!  Contains code related to generating composites for all GFS surface schemes.

module GFS_surface_composites_post

   use machine, only: kind_phys

   ! For consistent calculations of composite surface properties
   use sfc_diff, only: stability

   implicit none

   private

   public GFS_surface_composites_post_run

   real(kind=kind_phys), parameter :: zero = 0.0_kind_phys, one = 1.0_kind_phys, &
                                      half = 0.5_kind_phys, qmin = 1.0e-8_kind_phys

contains

!> \section arg_table_GFS_surface_composites_post_run Argument Table
!! \htmlinclude GFS_surface_composites_post_run.html
!!
   subroutine GFS_surface_composites_post_run (                                                                                   &
      im, kice, km, rd, rvrdm1, cplflx, cplwav2atm, frac_grid, flag_cice, thsfc_loc, islmsk, dry, wet, icy, wind, t1, q1, prsl1,  &
      landfrac, lakefrac, oceanfrac, zorl, zorlo, zorll, zorli, garea, frac_ice,                                                  &
      cd, cd_wat, cd_lnd, cd_ice, cdq, cdq_wat, cdq_lnd, cdq_ice, rb, rb_wat, rb_lnd, rb_ice, stress, stress_wat, stress_lnd,     &
      stress_ice, ffmm, ffmm_wat, ffmm_lnd, ffmm_ice, ffhh, ffhh_wat, ffhh_lnd, ffhh_ice, uustar, uustar_wat, uustar_lnd,         &
      uustar_ice, fm10, fm10_wat, fm10_lnd, fm10_ice, fh2, fh2_wat, fh2_lnd, fh2_ice, tsurf_wat, tsurf_lnd, tsurf_ice,            &
      cmm, cmm_wat, cmm_lnd, cmm_ice, chh, chh_wat, chh_lnd, chh_ice, gflx, gflx_wat, gflx_lnd, gflx_ice, ep1d, ep1d_wat,         &
      ep1d_lnd, ep1d_ice, weasd, weasd_lnd, weasd_ice, snowd, snowd_lnd, snowd_ice, tprcp, tprcp_wat,                             &
      tprcp_lnd, tprcp_ice, evap, evap_wat, evap_lnd, evap_ice, hflx, hflx_wat, hflx_lnd, hflx_ice, qss, qss_wat, qss_lnd,        &
      qss_ice, tsfc, tsfco, tsfcl, tsfc_wat, tisfc, hice, cice, tiice,                                                            &
      sigmaf, zvfun, lheatstrg, h0facu, h0facs, hflxq, hffac, stc, lkm, iopt_lake, iopt_lake_clm, use_lake_model,                                                               &
      grav, prsik1, prslk1, prslki, z1, ztmax_wat, ztmax_lnd, ztmax_ice, huge, errmsg, errflg)

      implicit none

      integer,                              intent(in) :: im, kice, km, lkm, iopt_lake, iopt_lake_clm
      logical,                              intent(in) :: cplflx, frac_grid, cplwav2atm, frac_ice
      logical,                              intent(in) :: lheatstrg
      logical, dimension(:),                intent(in) :: flag_cice, dry, icy
      logical, dimension(:),             intent(inout) :: wet
      integer, dimension(:),                intent(in) :: islmsk, use_lake_model
      real(kind=kind_phys), dimension(:),   intent(in) :: wind, t1, q1, prsl1, landfrac, lakefrac, oceanfrac,                   &
        cd_wat, cd_lnd, cd_ice, cdq_wat, cdq_lnd, cdq_ice, rb_wat, rb_lnd, rb_ice, stress_wat,                                  &
        stress_lnd, stress_ice, ffmm_wat, ffmm_lnd, ffmm_ice, ffhh_wat, ffhh_lnd, ffhh_ice, uustar_wat, uustar_lnd, uustar_ice, &
        fm10_wat, fm10_lnd, fm10_ice, fh2_wat, fh2_lnd, fh2_ice, tsurf_wat, tsurf_lnd, tsurf_ice, cmm_wat, cmm_lnd, cmm_ice,    &
        chh_wat, chh_lnd, chh_ice, gflx_wat, gflx_lnd, gflx_ice, ep1d_wat, ep1d_lnd, ep1d_ice, weasd_lnd, weasd_ice,            &
        snowd_lnd, snowd_ice, tprcp_wat, tprcp_lnd, tprcp_ice, evap_wat, evap_lnd, evap_ice, hflx_wat, hflx_lnd,                &
        hflx_ice, qss_wat, qss_lnd, qss_ice, tsfc_wat, zorlo, zorll, zorli, garea

      real(kind=kind_phys), dimension(:),   intent(inout) :: zorl, cd, cdq, rb, stress, ffmm, ffhh, uustar, fm10,               &
        fh2, cmm, chh, gflx, ep1d, weasd, snowd, tprcp, evap, hflx, qss, tsfc, tsfco, tsfcl, tisfc

      real(kind=kind_phys), dimension(:),   intent(inout) :: hice, cice
      real(kind=kind_phys), dimension(:),   intent(inout) :: sigmaf, zvfun, hflxq, hffac
      real(kind=kind_phys),                 intent(in   ) :: h0facu, h0facs
!     real(kind=kind_phys),                 intent(in   ) :: min_seaice
      real(kind=kind_phys),                 intent(in   ) :: rd, rvrdm1, huge

      real(kind=kind_phys), dimension(:,:), intent(in   ) :: tiice
      real(kind=kind_phys), dimension(:,:), intent(inout) :: stc

      ! Additional data needed for calling "stability"
      logical,                              intent(in   ) :: thsfc_loc
      real(kind=kind_phys),                 intent(in   ) :: grav
      real(kind=kind_phys), dimension(:),   intent(in   ) :: prsik1, prslk1, prslki, z1
      real(kind=kind_phys), dimension(:),   intent(in   ) :: ztmax_wat, ztmax_lnd, ztmax_ice

      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

      ! Local variables
      integer :: i, k
      real(kind=kind_phys) :: txl, txi, txo, wfrac, q0, rho
      ! For calling "stability"
      real(kind=kind_phys) :: tsurf, virtfac, tv1, thv1, tvs, z0max, ztmax
      real(kind=kind_phys) :: lnzorll, lnzorli, lnzorlo
!
      real(kind=kind_phys) :: tem1, tem2, gdx
      real(kind=kind_phys), parameter :: z0lo=0.1, z0up=1.0
!

      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      ! --- generate ocean/land/ice composites

       fractional_grid: if (frac_grid) then

        do i=1, im

          ! Three-way composites (fields from sfc_diff)
          txl   = landfrac(i)            ! land fraction
          wfrac = one - txl              ! ocean fraction
          txi   = cice(i) * wfrac        ! txi = ice fraction wrt whole cell
          txo   = max(zero, wfrac-txi)   ! txo = open water fraction

         !gflx(i)   = txl*gflx_lnd(i)  + txi*gflx_ice(i)  + txo*gflx_wat(i)
          ep1d(i)   = txl*ep1d_lnd(i)  + txi*ep1d_ice(i)  + txo*ep1d_wat(i)
          weasd(i)  = txl*weasd_lnd(i) + txi*weasd_ice(i)
          snowd(i)  = txl*snowd_lnd(i) + txi*snowd_ice(i)
         !tprcp(i)  = txl*tprcp_lnd(i) + txi*tprcp_ice(i) + txo*tprcp_wat(i)
!
          sigmaf(i) = txl*sigmaf(i)

          evap(i)   = txl*evap_lnd(i)  + txi*evap_ice(i)  + txo*evap_wat(i)
          hflx(i)   = txl*hflx_lnd(i)  + txi*hflx_ice(i)  + txo*hflx_wat(i)
          qss(i)    = txl*qss_lnd(i)   + txi*qss_ice(i)   + txo*qss_wat(i)
          gflx(i)   = txl*gflx_lnd(i)  + txi*gflx_ice(i)  + txo*gflx_wat(i)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         Call stability for consistent surface properties. Currently this comes from !
!         the GFS surface layere scheme (sfc_diff), regardless of the actual surface  !
!         layer parameterization being used - to be extended in the future            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          tsfc(i)   = ( txl * cdq_lnd(i) * tsfcl(i)                             &
                      + txi * cdq_ice(i) * tisfc(i)                             & ! DH* Ben had tsurf_ice(i), but GFS_surface_composites_post_run uses tice instead
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

          lnzorll = zero ; lnzorli = zero ; lnzorlo = zero
          if (zorll(i) /= huge) then
            lnzorll = log(zorll(i))
          endif
          if (zorli(i) /= huge) then
            lnzorli = log(zorli(i))
          endif
          if (zorlo(i) /= huge) then
            lnzorlo = log(zorlo(i))
          endif
          zorl(i) = exp(txl*lnzorll + txi*lnzorli + txo*lnzorlo)
 !        zorl(i) = exp(txl*log(zorll(i)) + txi*log(zorli(i)) + txo*log(zorlo(i)))
          z0max   = 0.01_kind_phys * zorl(i)
          ztmax   = exp(txl*log(ztmax_lnd(i)) + txi*log(ztmax_ice(i)) + txo*log(ztmax_wat(i)))

          ! Only actually need to call "stability" if multiple surface types exist...
          if(txl == one) then ! 100% land
            rb(i)     = rb_lnd(i)
            ffmm(i)   = ffmm_lnd(i)
            ffhh(i)   = ffhh_lnd(i)
            fm10(i)   = fm10_lnd(i)
            fh2(i)    = fh2_lnd(i)
            cd(i)     = cd_lnd(i)
            cdq(i)    = cdq_lnd(i)
            stress(i) = stress_lnd(i)
            uustar(i) = uustar_lnd(i)
          elseif(txo == one) then ! 100% open water
            rb(i)     = rb_wat(i)
            ffmm(i)   = ffmm_wat(i)
            ffhh(i)   = ffhh_wat(i)
            fm10(i)   = fm10_wat(i)
            fh2(i)    = fh2_wat(i)
            cd(i)     = cd_wat(i)
            cdq(i)    = cdq_wat(i)
            stress(i) = stress_wat(i)
            uustar(i) = uustar_wat(i)
          elseif(txi == one) then ! 100% ice
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
!
! re-compute zvfun with composite surface roughness & green vegetation fraction
!
            tem1 = (z0max - z0lo) / (z0up - z0lo)
            tem1 = min(max(tem1, zero), one)
            tem2 = max(sigmaf(i), 0.1)
            zvfun(i) = sqrt(tem1 * tem2)
            gdx = sqrt(garea(i))
!
! re-compute variables for canopy heat storage parameterization with the updated zvfun
!      in the fractional grid
!
            hflxq(i) = hflx(i)
            hffac(i) = 1.0
            if (lheatstrg) then
              if(hflx(i) > 0.) then
                hffac(i) = h0facu * zvfun(i)
              else
                hffac(i) = h0facs * zvfun(i)
              endif
              hffac(i) = 1. + hffac(i)
              hflxq(i) = hflx(i) / hffac(i)
            endif
!
            call stability(z1(i), zvfun(i), gdx, tv1, thv1, wind(i),                & ! inputs
                           z0max, ztmax, tvs, grav, thsfc_loc,                      & ! inputs
                           rb(i), ffmm(i), ffhh(i), fm10(i), fh2(i), cd(i), cdq(i), & ! outputs
                           stress(i), uustar(i))
          endif ! Checking to see if point is one or multiple surface types

          ! BWG, 2021/02/25: cmm=cd*wind, chh=cdq*wind, so use composite cd, cdq
          rho       = prsl1(i) / (rd*tv1)
          cmm(i)    =      cd(i)*wind(i)  !txl*cmm_lnd(i)    + txi*cmm_ice(i)    + txo*cmm_wat(i)
          chh(i)    = rho*cdq(i)*wind(i)  !txl*chh_lnd(i)    + txi*chh_ice(i)    + txo*chh_wat(i)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          if (dry(i)) then
          elseif (wet(i)) then
            tsfcl(i) = tsfc_wat(i)       ! over water
          else
            tsfcl(i) = tisfc(i)           ! over ice
          endif
          if (wet(i)) then
            tsfco(i) = tsfc_wat(i)       ! over lake or ocean when uncoupled
          elseif (icy(i)) then
            tsfco(i) = tisfc(i)           ! over lake or ocean ice when uncoupled
          else
            tsfco(i) = tsfcl(i)       ! over land
          endif
          if (icy(i)) then
            !tisfc(i) = tisfc(i)           ! over lake or ocean ice when uncoupled
          elseif (wet(i)) then
            tisfc(i) = tsfc_wat(i)       ! over lake or ocean when uncoupled
          else
            tisfc(i) = tsfcl(i)       ! over land
          endif
                                                  ! for coupled model ocean will replace this
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

      else ! not fractional grid

        do i=1,im

          if (use_lake_model(i)>0) then
            if(frac_ice .and. icy(i)) then
              call composite_icy(iopt_lake==iopt_lake_clm)
              call composite_wet_and_icy
            else
              call composite_wet
            endif
          else if (islmsk(i) == 1) then
          !-- land
            call composite_land
          elseif (islmsk(i) == 0) then
          !-- water
            call composite_wet
          else ! islmsk(i) == 2
          !-- ice
            call composite_icy(.false.)
            call composite_wet_and_icy
          endif
        enddo

      endif fractional_grid

      ! --- compositing done

   contains

     subroutine composite_land
            implicit none
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
            tsfc(i)   = tsfcl(i)
            tsfco(i)  = tsfc(i)
            tisfc(i)  = tsfc(i)
            cmm(i)    = cmm_lnd(i)
            chh(i)    = chh_lnd(i)
            gflx(i)   = gflx_lnd(i)
            ep1d(i)   = ep1d_lnd(i)
            weasd(i)  = weasd_lnd(i)
            snowd(i)  = snowd_lnd(i)
            evap(i)   = evap_lnd(i)
            hflx(i)   = hflx_lnd(i)
            qss(i)    = qss_lnd(i)
            hice(i)   = zero
            cice(i)   = zero
     end subroutine composite_land

     subroutine composite_wet
            implicit none
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
            weasd(i)  = zero
            snowd(i)  = zero
            evap(i)   = evap_wat(i)
            hflx(i)   = hflx_wat(i)
            qss(i)    = qss_wat(i)
            hice(i)   = zero
            cice(i)   = zero
     end subroutine composite_wet

     subroutine composite_icy(is_clm)
            implicit none
            logical, intent(in) :: is_clm
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
            if(is_clm) then
              weasd(i)  = weasd_ice(i)
              snowd(i)  = snowd_ice(i)
            else
              weasd(i)  = weasd_ice(i) * cice(i)
              snowd(i)  = snowd_ice(i) * cice(i)
            endif
            qss(i)    = qss_ice(i)
            evap(i)   = evap_ice(i)
            hflx(i)   = hflx_ice(i)
     end subroutine composite_icy

     subroutine composite_wet_and_icy
            implicit none
            txi = cice(i)
            txo = one - txi
            evap(i)   = txi * evap_ice(i)   + txo * evap_wat(i)
            hflx(i)   = txi * hflx_ice(i)   + txo * hflx_wat(i)
            tsfc(i)   = txi * tisfc(i)      + txo * tsfc_wat(i)
            stress(i) = txi * stress_ice(i) + txo * stress_wat(i)
            qss(i)    = txi * qss_ice(i)    + txo * qss_wat(i)
            ep1d(i)   = txi * ep1d_ice(i)   + txo * ep1d_wat(i)

            lnzorli = zero ; lnzorlo = zero
            if (zorli(i) /= huge) then
              lnzorli = log(zorli(i))
            endif
            if (zorlo(i) /= huge) then
              lnzorlo = log(zorlo(i))
            endif
            zorl(i) = exp(txi*lnzorli + txo*lnzorlo)
!           zorl(i)   = exp(txi*log(zorli(i)) + txo*log(zorlo(i)))
!
            if (wet(i)) then
              tsfco(i) = tsfc_wat(i)
            else
              tsfco(i) = tsfc(i)
            endif
            tsfcl(i)  = tsfc(i)
            do k=1,min(kice,km) ! store tiice in stc to reduce output in the nonfrac grid case
              stc(i,k) = tiice(i,k)
            enddo
     end subroutine composite_wet_and_icy

   end subroutine GFS_surface_composites_post_run

end module GFS_surface_composites_post

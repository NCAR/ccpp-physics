! GW SOLVERS:
!=========== SOLVER_ORODIS; SOLVER_WMSDIS, SOLVER_LSATDIS
!          + RF_DAMP if it is needed along with ugwp_tofd
!===========
! Note in contrast to dycore vertical indices: surface=1   top=levs
!
! Collection of main friction-GWD solvers
!
! subroutine ugwp_oro
! 
! subroutine gw_solver_linsatdis
! subroutine gw_solver_wmsdis
! subroutine rf_damp
!
! ===========
!            
!
     subroutine ugwp_oro(im, levs, dtp, kdt,me, lprnt, fcor, c2f2,       &
       u, v, tkin, pint, delp, pmid, pexner, gzint, gzmid, orostat,      &
       hpbl, axz, ayz, edis, kdis, dusfc, dvsfc,                         &
       dusfc_mb, dvsfc_mb, dusfc_ogw, dvsfc_ogw, dusfc_lwb, dvsfc_lwb,   &
       zmtb, zlwb, zogw, tauf_ogw, tauz_ogw, axmtb,  axlwb, axtms        )
!----------------------------------------------------------------------
! COORDE-output:  6-hour inst: U, V, T, PMSL, PS, HT (ounce)
!    3D 6-hr   aver: DYN-U, SSO-U, PBL-U, AF-U1....
!    2D 6-hr   aver: tau_SSO, tau_GWD, tau_BL; &
!              tau_sso = tau_mtb + tau_tofd + tau_lwb +tau_ogw
!    ZM 6-hr   aver: tau_RES = PS*dH/dx -zonal mean
! Experiments: Midlat 80-200km
!              LR_CTL;            ; LR_NOSSO with TOFD/TMS;
!              LR_NOGWD (MTN+TOFD); LR_GWD4 --- 4 times taub
!----------------------------------------------------------------------
      use machine ,          only : kind_phys
      use ugwp_oro_init,     only : cdmb, cleff, sigfac, hncrit, hpmin, hminmt
      use ugwp_oro_init,     only : gamm_std, sigma_std
      use ugwp_common ,      only : rgrav, grav, cpd, rd, rv, rcpd, rcpd2

     
      use ugwp_common ,      only : pi, rad_to_deg, deg_to_rad, pi2

      use cires_ugwp_module, only : kxw,  max_kdis, max_axyz

      implicit none
      logical                                     ::   lprnt
      integer                                     ::   im, levs
      integer                                     ::   me
      integer                                     ::   kdt
      real(kind_phys)                             ::   dtp
      real(kind_phys), dimension(im)              ::   hpbl                          !   pbl-height in meters
      real(kind_phys), dimension(im)              ::   fcor, c2f2
      real(kind_phys), dimension(im, 14)          ::   orostat
      real(kind_phys), dimension(im, levs)        ::   u, v, tkin, q

      real(kind_phys), dimension(im, levs)        ::   pmid, pexner, gzmid, delp
      real(kind_phys), dimension(im, levs+1)      ::   pint, gzint
 

      real(kind_phys), dimension(im, levs)        ::   axz, ayz, edis, kdis           ! total 6-hr averaged tendencies
      real(kind_phys), dimension(im, levs)        ::   krf2d
      real(kind_phys), dimension(im, levs)        ::   tauz_ogw, axmtb, axlwb, axtms  ! 3-sub components axogw = axz-(axmtb+axlwb+axtms)
      real(kind_phys), dimension(im)              ::   tauf_ogw                       ! total-source momentum flux

      real(kind_phys), dimension(im)              ::   zmtb, zlwb, zogw

      real(kind_phys), dimension(im)              ::   dusfc,      dvsfc              ! total tausfc_sso
      real(kind_phys), dimension(im)              ::   dusfc_mb,   dvsfc_mb           ! integrated  tau_mtb
      real(kind_phys), dimension(im)              ::   dusfc_ogw,  dvsfc_ogw          ! integrated  tau_ogw
      real(kind_phys), dimension(im)              ::   dusfc_lwb,  dvsfc_lwb          ! integrated  tau_lwb
      real(kind_phys), dimension(im)              ::   dusfc_tofd, dvsfc_tofd         ! integrated  tau_tofd

!
!                         mu=hprime       gamm=a/b         sigma           theta
! which stand for the standard deviation, the anisotropy, the slope and the orientation  of  the  orography.
!
      real(kind_phys) :: elvmax(im)
      real(kind_phys) :: hprime(im)

      real(kind_phys) :: theta                        !the orienatation, angle
      real(kind_phys) :: sigma                        !the slope dh/dx
      real(kind_phys) :: gamm                         !the anisotropy see ifs-oro

      real(kind_phys) :: oc, oa4(4),  clx4(4)         !kim & doyle 2005 .... attempt to do TOFD ..?
!
      integer, allocatable ::  k_elev(:), k_mtb(:), k_ogw(:), k_lee(:), k_tofd(:)

      real(kind_phys) wk(im)

      real(kind_phys) eng0, eng1
!
!
!
      real(kind_phys), dimension(levs)   :: up, vp, tp, qp, dp, zpm, pmid1, pex

      real(kind_phys), dimension(levs+1) :: taudz, rhoi, rim_z, pint1, zpi
      real(kind_phys), dimension(levs)   ::  drtau, kdis_oro
!
      real (kind_phys)                 :: elvp,  elvpd, dtaux, dtauy
      real(kind_phys)                  :: loss, mtb_fric, mbx, mby
      real(kind_phys)                  :: sigflt

      real(kind_phys)                  :: zpbl = 2000.   ! can be passed from PBL physics as in gwdps.f
!
      logical icrilv(im)
!
!----   mountain/oro gravity wave drag +TOFD
!
      real(kind=kind_phys), dimension(levs) :: utofd1, vtofd1, epstofd1, krf_tofd1
!
      real(kind=kind_phys), dimension(levs) :: drlee,  drmtb,  drlow,  drogw
      real(kind_phys)                       :: r_cpdt, acc_lim
      real(kind_phys),   dimension(im)      :: tautot, tauogw, taumtb, taulee,  taurf
      real(kind_phys)                       :: xn,     yn,     umag,   kxridge,  &
                                               tx1,    tx2
      real(kind=kind_phys),dimension(levs+1):: tau_src

      integer   ::  npt, krefj, kdswj, kotr, i, j, k
      integer   ::  ipt(im)

!
! copy 1D
!
         do i=1, im
           hprime(i) = orostat(i,  1)
           elvmax(i) = orostat(i, 14) 
!
           tautot(i) = 0.0
           tauogw(i) = 0.0
           taumtb(i) = 0.0
           taulee(i) = 0.0
           taurf(i)  = 0.0
!
           dusfc(i)      = 0.0
           dvsfc(i)      = 0.0
           dusfc_mb(i)   = 0.0
           dvsfc_mb(i)   = 0.0
           dusfc_ogw(i)  = 0.0
           dvsfc_ogw(i)  = 0.0
           dusfc_lwb(i)  = 0.0
           dvsfc_lwb(i)  = 0.0
           dusfc_tofd(i) = 0.0
           dvsfc_tofd(i) = 0.0
           tauf_ogw(i)   = 0.0
!
           zmtb(i)       = -99.
           zlwb(i)       = -99.
           zogw(i)       = -99.
           ipt(i)        = 0
         enddo
!	        print *, maxval(hprime), maxval(elvmax), ' check hprime -elevmax ugwp_oro'
!
! 3-part of oro-effects + ked_oro
!
         do k=1, levs
           do i=1, im
             axz(i,k)      = 0.0
             ayz(i,k)      = 0.0
             edis(i,k)     = 0.0
             kdis(i,k)     = 0.0
             krf2d(i,k)    = 0.0
             tauz_ogw(i,k) = 0.0
             axmtb(i:,k)   = 0.0
             axlwb(i,k)    = 0.0
             axtms(i,k)    = 0.0 
           enddo
         enddo
 
!
!   optional diag 3-parts of drag: [tx_ogw, tx_mtb, tx_lee]
!
! ----do we have orography  for mtb and gwd calculation points ?
!
        npt = 0
        do i = 1,im
          if ( (elvmax(i) > hminmt) .and. (hprime(i) > hpmin) )  then
             npt      = npt + 1
             ipt(npt) = i

          endif
        enddo
        if (npt == 0) return       ! no ororgraphy  ====> gwd/mb calculation done

!        allocate(iwklm(npt), idxzb(npt), kreflm(npt))
        allocate( k_elev(npt), k_mtb(npt), k_ogw(npt), k_lee(npt), k_tofd(npt))
        do i=1,npt
          k_ogw (i) = 2
          k_tofd(i) = 2
          k_lee (i) = 2
          k_mtb(i)  = 0
          k_elev(i) = 2
        enddo
!
!                      controls through: use ugwp_oro_init
! main ORO-loop        sigfac = n*sigma = [1.5, 2, 2.5, 4]*hprime
!


        do i = 1, npt
!
          j = ipt(i)

          elvpd =  elvmax(j)
          elvp  = min (elvpd  + sigfac * hprime(j), hncrit)

          sigma = orostat(j,13)
          gamm  = orostat(j,12)
          theta = orostat(j,11)*deg_to_rad

          if (sigma == 0.0 ) then
            sigma = sigma_std
            gamm  = gamm_std
            theta = 0.0
          endif

          oc      = orostat(j,2)
          oa4(1)  = orostat(j,3)
          oa4(2)  = orostat(j,4)
          oa4(3)  = orostat(j,5)
          oa4(4)  = orostat(j,6)
          clx4(1) = orostat(j,7) 
          clx4(2) = orostat(j,8)
          clx4(3) = orostat(j,9) 
          clx4(4) = orostat(j,10)
!
! do column-based diagnostics "more-efficient" for oro-places
!

          do k=1,levs
            up(k)    = u(j,k)
            vp(k)    = v(j,k)
            tp(k)    = tkin(j,k)
            qp(k)    = q(j,k)
            dp(k)    = delp(j,k)

            zpm(k)   = gzmid(j,k) * rgrav
            pmid1(k) = pmid(j,k)
            pex(k)   = pexner(j,k)
          enddo
          do k=1,levs+1
            zpi(k)   = gzint(j,k) * rgrav
            pint1(k) = pint(j,k)
          enddo
!
! elvp- k-index: iwklm k_elvp = index for elvmax + 4*hprime, "elevation index"
! GFS-2017
          do k=1, levs-1
            if (elvp <= zpi(k+1) .and. elvp > zpi(k)) then
              k_elev(i) =  k+1      !......simply k+1 next interface level
              exit
            endif
          enddo
!         if (elvp  .ge. 300. ) then
!           write(6,333) elvp, zpi(1), elvpd, hprime(j), sigfac, hncrit
!           pause
!         endif
!333       format(6(3x, F10.3))
!
! SSO effects: TOFD-drag/friction coefficients can be calculated
!	  
          sigflt = hprime(j)*0.01  ! turb SSo(j) ...small-scale orography < 2-5 km ....
          zpbl   = hpbl(j)
 
          call ugwp_tofd1d(levs, sigflt, elvPd, zpi(1), zpbl, up, vp, zpm,   &
                           utofd1, vtofd1, epstofd1, krf_tofd1)

          do k=1, levs
            krf2d(j,k) = krf_tofd1(k)
            axtms(j,k) = utofd1(k)
!-------
! nullify ORO-tendencies
!
            drmtb(k) = 0.0
            drlee(k) = 0.0
            drtau(k) = 0.0
            drlow(k) = 0.0
          enddo

!-------
!
! levels of          k_mtb(i)/mtb + kdswj/dwlee + krefj/ogwd inside next "subs"
!                      zmtb,         zlwb,          zogw
!                      drmtb,       drlow/drlee,    drogw
!------- 
!
! mtb : drmtb => 1-st order friction as well as TurbulentOro-Drag
!
          call ugwp_drag_mtb( k_elev(i), levs,                                        &
              elvpd, elvp, hprime(j), sigma, theta, oc,  oa4, clx4, gamm,  zpbl,      &
              up, vp, tp, qp, dp, zpm, zpi, pmid1, pint1, k_mtb(i), drmtb, taumtb(j))

          axmtb(j,1:levs) = drmtb(1:levs)*up(1:levs)
! 
!      print * ,  k_elev(i),    k_mtb(i) , taumtb(j)*1.e3, ' k_elev,    k_mtb , taumtb '
!
! tautot = taulee+tauogw + rho*drlee = -d[taulee(z)]/dz
!


          call ugwp_taub_oro(levs, k_mtb(i), kxw,  taumtb(j),  fcor(j),   &
            hprime(j) , sigma, theta, oc,  oa4, clx4, gamm, elvp,         &
            up, vp, tp, qp, dp, zpm, zpi, pmid1, pint1, xn, yn, umag,     &
            tautot(j), tauogw(j), taulee(j), drlee, tau_src,              &
            kxridge, kdswj, krefj, kotr)

!     print *, k_mtb(i), kxw, taumtb(j), fcor(j),hprime(j), ' af ugwp_taub_oro '
!     print *, kdswj, krefj, kotr, ' kdswj, krefj, kotr '


          tauf_ogw(j) = tautot(j)
          axlwb(j,1:levs) = drlee(1:levs)

          if ( k_mtb(i) > 0) zmtb(j) = zpi(k_mtb(i))- zpi(1)
          if ( krefj    > 0) zogw(j) = zpi(krefj)   - zpi(1)
          if ( kdswj    > 0) zlwb(j) = zpi(kdswj)   - zpi(1)
!         if ( k_mtb(i) > 0  .and.  zmtb(j) > zogw(j)) print *, ' zmtb > zogw ',   zmtb(j), zogw(j)
!
! tau: tauogw, kxw/kxridge  ATTENTION c2f2(j) = fcor(j)*fcor(j)/kxridge/kxridge
!
          if ( (krefj > 1) .and. ( abs(tauogw(j)) > 0.) ) then
!
            call ugwp_oro_lsatdis( krefj, levs,  tauogw(j),  tautot(j), tau_src, kxw,     &
                 fcor(j), kxridge, up, vp, tp, qp, dp, zpm, zpi, pmid1, pint1,            &
                 xn, yn, umag, drtau, kdis_oro)
!
          else
            drtau = 0.
          endif

          tauz_ogw(j,1:levs) = tau_src(1:levs)

          r_cpdt  = rcpd2/dtp
!
!
          do k = 1,levs
!
! project to x-dir & y=dir and do diagnostics
! & apply limiters and output separate oro-effects
!
            drlow(k) = drtau(k) + drlee(k)
            acc_lim  = min(abs(drlow(k)), max_axyz)
            drlow(k) = sign(acc_lim, drlow(k))

            dtaux    = drlow(k) * xn + utofd1(k)
            dtauy    = drlow(k) * yn + vtofd1(k)

            eng0     = up(k)*up(k)+vp(k)*vp(k)
            eng1     = 0.0
!
            if (k < k_mtb(i) .and. drmtb(k) /= 0 ) then
              loss     = 1.0 / (1.0+drmtb(k)*dtp)
              mtb_fric = drmtb(k)*loss
!
              mbx      =  mtb_fric * up(k)
              mby      =  mtb_fric * vp(k)
!
              ayz(j,k) = -mby        !+ ayz(j,k)
              axz(j,k) = -mbx        !+ axz(j,k)     
!
              eng1     = eng0*loss*loss +eng1
              dusfc(j) = dusfc(j) - mbx * dp(k)
              dvsfc(j) = dvsfc(j) - mby * dp(k)
            endif
!
            ayz(j,k)   = dtauy + ayz(j,k)
            axz(j,k)   = dtaux + axz(j,k)
!
            tx1  = u(j,k)  + dtaux*dtp
            tx2  = v(j,k)  + dtauy*dtp
            eng1 = tx1*tx1 + tx2*tx2 + eng1

            dusfc(j)  = dusfc(j)  + dtaux * dp(k)
            dvsfc(j)  = dvsfc(j)  + dtauy * dp(k)

            edis(j,k) = max(eng0-eng1, 0.0) * r_cpdt      !+ epstofd1(k)
            kdis(j,k) = min(kdis_oro(k), max_kdis )

          enddo
!
          dusfc(j) = -rgrav * dusfc(j)
          dvsfc(j) = -rgrav * dvsfc(j)
!
! oro-locations
!
        enddo          !          ipt - oro-loop .... "fraction of Land" in the grid box
        deallocate(k_elev, k_mtb, k_ogw, k_lee, k_tofd  )
!
      end subroutine ugwp_oro
!
!
      subroutine gw_solver_linsatdis(im, levs, dtp, kdt, me,          &
        taub, klev, if_src, nf_src, nw,  ch, naz,  spf,  xaz, yaz,    &
        fcor, c2f2, u, v, t,  q, prsi, delp, prsl, prslk, phii, phil, &
        ax, ay, eps, ked, tauz)

      use ugwp_common ,      only : rgrav, grav, cpd, rd, rv, rcpd, rcpd2
      use ugwp_common ,      only : pi, rad_to_deg, deg_to_rad, pi2

      use cires_ugwp_module, only : kxw,  max_kdis, max_axyz, max_eps
      use cires_ugwp_module, only : kvg, ktg, krad, kion

        implicit none
        integer :: im, levs
        integer :: me, kdt, nw, naz, nf_src
        real    :: dtp
        integer, dimension(im) ::  klev, if_src
        real,    dimension(im) ::  taub, fcor, c2f2

        real,   dimension(naz) ::  xaz, yaz
        real,   dimension(nw ) ::  ch, spf
!==========================
        real, dimension(im, levs)   ::  u, v, t, delp, prsl, prslk, phil, q
        real, dimension(im, levs+1) :: prsi , phii
!==========================
        real, dimension(im, levs)   :: ax, ay, eps, ked, tauz

        real, dimension(levs)   ::  u1, v1, t1, dp, pmid, zmid, pex1, &
                                    q1, rho
        real, dimension(levs+1) ::  pint , zint, ui, vi, ti,          &
                                    bn2i,  bvi, rhoi
        integer                 :: i, j, k, ksrc
        real, dimension(nw)     :: taub_spect
!       real, dimension(levs)   :: ax1, ay1, eps1
!       real, dimension(levs+1) :: ked1, tau1
        real                    :: chm, ss
        real, parameter         :: dsp = 1./20.
        logical                 :: pfirst=.true.

        save pfirst
128     Format (2x, I4, 4(2x, F10.3))

!       do i=1, nw
!         spf(i) = exp(-Ch(i)*dsp)
!       enddo
!       ss = sum(spf)
!       spf(1:nw) =  spf(1:nw)/ss

        if (pfirst ) then
          j = 1
          ksrc = klev(j)
          taub_spect(1:nw) = spf(1:nw)*taub(j)
          print *
          chm = 0.
          do i=1, nw
            write(6, 128)  i, spf(i), taub_spect(i)*1.e3, ch(i), ch(i)-chm
            chm = ch(i)
          enddo

          print *
          pause
        endif

        do j=1,im
          if (if_src(j) == 1) then
!
! compute GW-effects
! prsi, delp, prsl, prslk, phii, phil
!
            do k=1,levs
              u1(k)   = u(j,k)
              v1(k)   = v(j,k)
              t1(k)   = t(j,k)
              q1(k)   = q(j,k)       ! H2O-index -1 in tracer-array
              dp(k)   = delp(j,k)

              zmid(k) = phil(j,k) * rgrav
              pmid(k) = prsl(j,k)
!             pex1(k) = prslk(j,k)
            enddo
            do k=1,levs+1
              zint(k) = phii(j,k) * rgrav
              pint(k) = prsi(j,k)
            enddo

            call mflow_tauz(levs, u1, v1, t1, q1, dp, zmid, zint, &
                            pmid, pint, rho, ui, vi, ti,  bn2i, bvi, rhoi)
!
            ksrc = klev(j)
            taub_spect(1:nw) = spf(1:nw)*taub(j)/rhoi(ksrc)
            if (pfirst .and. j ==1 ) then

              print *, maxval(taub_spect)/kxw*bvi(ksrc)/ch(1), ' Urms '
              print *, maxval(zmid), minval(zmid) , ' zmid '
              print *, maxval(zint), minval(zint) , ' zint '
              print *, maxval(rho),  minval(rho)  , ' rho  '
              print *, maxval(rhoi), minval(rhoi) , ' rhoi  '
              print *, maxval(ti),   minval(ti)   , ' tempi  '
              print *, maxval(ui),   minval(ui)   , ' ui  '
              print *, maxval(u1),   minval(u1)   , ' ++++ u1  '
              print *, maxval(vi),   minval(vi)   , ' vi  '
              print *, maxval(v1),   minval(v1)   , ' ++++ v1  ' 
              print *, maxval(pint), minval(pint) , ' pint  '
               pause
             endif
!
             call ugwp_lsatdis_naz(levs, ksrc, nw, naz, kxw, taub_spect,      &
                                   ch, xaz, yaz, fcor(j), c2f2(j), dp,        &
                                   zmid, zint, pmid, pint, rho, ui, vi, ti,   &
                                   kvg, ktg, krad, kion, bn2i, bvi, rhoi,     &
                                   ax(j,1:levs), ay(j,1:levs), eps(j,1:levs), &
                                   ked(j,1:levs), tauz(j,1:levs))
!             kvg, ktg, krad, kion, bn2i, bvi, rhoi, ax1, ay1, eps1, ked1, tau1)

             if (pfirst .and. j ==1 ) then

               print *, maxval(taub_spect)/kxw*bvi(ksrc)/ch(1), ' Urms '
               print *, maxval(zmid), minval(zmid) , ' zmid '
               print *, maxval(zint), minval(zint) , ' zint '
               print *, maxval(rho),  minval(rho)  , ' rho  '
               print *, maxval(rhoi), minval(rhoi) , ' rhoi '
               print *, maxval(ti),   minval(ti)   , ' rhoi '
               print *, maxval(ui),   minval(ui)   , ' ui  '
               print *, maxval(vi),   minval(vi)   , ' vi  '
               print *, maxval(pint), minval(pint) , ' pint '
               pause
             endif
!
!            ax(j,:)   = ax1
!            ay(j,:)   = ay1
!            eps(j,:)  = eps1
!            ked(j,:)  = ked1(1:levs)
!            tauz(j,:) = tau1(1:levs)
          endif

        enddo
        pfirst = .false.
!
! spectral solver for discrete spectra of GWs in N-azimiths
! Linear saturation with background dissipation
!
      end subroutine gw_solver_linsatdis
!
      subroutine gw_solver_wmsdis(im, levs, dtp, kdt, me,             &
        taub, klev, if_src, nf_src, nw,  ch, naz,  spf,  xaz, yaz,    &
        fcor, c2f2, u, v, t,  q, prsi, delp, prsl, prslk, phii, phil, &
        ax, ay, eps, ked, tauz)
!     use para_taub,         only :  tau_ex
      use ugwp_common ,      only : rgrav, grav, cpd, rd, rv, rcpd, rcpd2
      use ugwp_common ,      only : pi, rad_to_deg, deg_to_rad, pi2

      use cires_ugwp_module, only : kxw,  max_kdis, max_axyz, max_eps
      use cires_ugwp_module, only : kvg, ktg, krad, kion
 
       implicit none
       integer :: im, levs, me, kdt, nw, naz, nf_src
       real    :: dtp

       integer, dimension(im) ::  klev, if_src
       real,    dimension(im) ::  taub, fcor, c2f2

       real,   dimension(naz) ::  xaz, yaz
       real,   dimension(nw ) ::  ch, spf
!==========================
       real, dimension(im, levs)   ::  u, v, t, delp, prsl, prslk, phil, q
       real, dimension(im, levs+1) :: prsi , phii
!==========================
       real, dimension(im, levs)   :: ax, ay, eps, ked, tauz

       real, dimension(levs)   :: u1, v1, t1, dp, pmid, zmid, pex1, q1, rho
       real, dimension(levs+1) :: pint , zint, ui, vi, ti, bn2i, bvi, rhoi

       integer                 :: i, j, k, ksrc
       real, dimension(nw)     :: taub_spect
!      real, dimension(levs)   :: ax1, ay1, eps1
!      real,dimension(levs+1)  :: ked1, tau1
       real                    :: tau_ex

!	print *, nf_src, 'nf_src ... gw_solver_wmsdis '
!	print *, if_src,  'if_src ... gw_solver_wmsdis '

        do j=1,im
          if (if_src(j) == 1) then
!
! compute gw-effects
! prsi, delp, prsl, prslk, phii, phil
!
            do k=1,levs
              u1(k)   = u(j,k)
              v1(k)   = v(j,k)
              t1(k)   = t(j,k)
              q1(k)   = q(j,k)       ! h2o-index -1 in tracer-array
              dp(k)   = delp(j,k)

              zmid(k) = phil(j,k)  *rgrav
              pmid(k) = prsl(j,k)
!             pex1(k) = prslk(j,k)
            enddo
            do k=1,levs+1
              zint(k) = phii(j,k)*rgrav
              pint(k) = prsi(j,k)
            enddo

            call mflow_tauz(levs, u1, v1, t1, q1, dp, zmid, zint, &
                            pmid, pint, rho, ui, vi, ti,  bn2i, bvi, rhoi)
!
! any extras bkg-arrays
!
            ksrc = klev(j)
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! more work for spectral setup for different "slopes"
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~	
            tau_ex = taub(j)
            taub_spect(1:nw) = spf(1:nw)/rhoi(ksrc) *tau_ex    ! check it ....*tau_ex(j)
 
!	 
!  call FVS93_ugwps(nw,  ch,  dch, taub_spect, spnorm,  nslope, bn2i(ksrc), bvi(ksrc), bnrho(ksrc))
!
!	print *, ' bf ugwp_wmsdis_naz ksrc', ksrc,  zmid(ksrc)

            call ugwp_wmsdis_naz(levs, ksrc, nw, naz, kxw, tau_ex, ch, xaz, yaz,    &
                                 fcor(j), c2f2(j), dp, zmid, zint, pmid, pint,      &
                                 rho, ui, vi, ti,  kvg, ktg, krad, kion, bn2i, bvi, &
                                 rhoi, ax(j,1:levs), ay(j,1:levs), eps(j,1:levs),   &
                                 ked(j,1:levs), tauz(j,1:levs))
!                kvg, ktg, krad, kion, bn2i, bvi, rhoi, ax1, ay1, eps1, ked1, tau1)
 
!	print *, ' after ugwp_wmsdis_naz ksrc', ksrc,  zint(ksrc)

!   subroutine ugwp_wmsdis_naz(levs, ksrc, nw, naz, kxw, taub_lat, ch, xaz, yaz,    &
!              fcor, c2f2, dp, zmid, zint, pmid, pint, rho, ui, vi, ti,             &
!              kvg, ktg, krad, kion, bn2i, bvi, rhoi, ax, ay, eps, ked)
 
!           ax(j,:)   = ax1
!           ay(j,:)   = ay1  
!           eps(j,:)  = eps1
!           ked(j,:)  = ked1(1:levs)
!           tauz(j,:) = tau1(1:levs)

          endif

        enddo
!
! ugwp_wmsdis_naz everything similar to linsat , except  spectral saturation
!
!
        return
      end subroutine gw_solver_wmsdis
!
!
      subroutine rf_damp(im, levs, levs_rf, dtp, rfdis, rfdist, u, v, ax, ay, eps)
        use ugwp_common, only : rcpd2

        implicit  none

        integer :: im, levs, levs_rf
        real    :: dtp
        real, dimension(levs)     :: rfdis, rfdist
        real, dimension(im, levs) :: u, v, ax, ay, eps
        real       :: ud, vd, rdtp
        integer    :: i, k
 
        rdtp  = 1.0 / dtp

        do k= levs_rf, levs
          do i=1,im
           ud = rfdis(k)*u(i,k)
           vd = rfdis(k)*u(i,k)
           ax(i,k)  = rfdist(k)*u(i,k)
           ay(i,k)  = rfdist(k)*v(i,k)
           eps(i,k) = rcpd2*(u(i,k)*u(i,k) +v(i,k)*v(i,k) -ud*ud -vd*vd)
          enddo
        enddo
      end subroutine rf_damp
!                   

!> \file moninshoc.f
!!  Contains most of the SHOC PBL/shallow convection scheme.

!> This module contains the CCPP-compliant SHOC scheme.
      module moninshoc

      contains

      subroutine moninshoc_init ()
      end subroutine moninshoc_init

      subroutine moninshoc_finalize ()
      end subroutine moninshoc_finalize

!!!!!  ==========================================================  !!!!!
! subroutine 'moninshoc' computes pbl height and applies vertical diffusion
! using the coefficient provided by the SHOC scheme (from previous step)
! 2015-05-04 - Shrinivas Moorthi - original version based on monin
! 2018-03-21 - Shrinivas Moorthi - fixed a bug related to tke vertical diffusion
!                                  and gneralized the tke location in tracer array
! 2018-03-23 - Shrinivas Moorthi - used twice the momentum diffusion coefficient
!                                  for tke as in Deardorff (1980) - added tridi1
!
!> \section arg_table_moninshoc_run Argument Table
!! \htmlinclude moninshoc_run.html
!!
      subroutine moninshoc_run (im,km,ntrac,ntcw,ncnd,dv,du,tau,rtg,
     &                          u1,v1,t1,q1,tkh,prnum,ntke,
     &                          psk,rbsoil,zorl,u10m,v10m,fm,fh,
     &                          tsea,heat,evap,stress,spd1,kpbl,
     &                          prsi,del,prsl,prslk,phii,phil,delt,
     &                          dusfc,dvsfc,dtsfc,dqsfc,dkt,hpbl,
     &                          kinver,xkzm_m,xkzm_h,xkzm_s,xkzminv,
     &                          grav,rd,cp,hvap,fv,ntoz,dt3dt_PBL,
     &                          du3dt_PBL,dv3dt_PBL,dq3dt_PBL,do3dt_PBL,
     &                          gen_tend,ldiag3d,qdiag3d,
     &                          errmsg,errflg)
!
      use machine  , only : kind_phys
      use funcphys , only : fpvs

      implicit none
!
!     arguments
!
      integer,                                  intent(in) :: im,
     &  km, ntrac, ntcw, ncnd, ntke, ntoz
      integer, dimension(im),                   intent(in) ::  kinver

      real(kind=kind_phys),                     intent(in) :: delt,
     &  xkzm_m, xkzm_h, xkzm_s, xkzminv
      real(kind=kind_phys),                     intent(in) :: grav,
     &  rd, cp, hvap, fv
      real(kind=kind_phys), dimension(im),      intent(in) :: psk,
     &  rbsoil, zorl, u10m, v10m, fm, fh, tsea, heat, evap, stress, spd1
      real(kind=kind_phys), dimension(im,km),   intent(in) :: u1, v1,
     &  t1, tkh, del, prsl, phil, prslk
      real(kind=kind_phys), dimension(im,km+1), intent(in) :: prsi, phii
      real(kind=kind_phys), dimension(im,km,ntrac), intent(in) :: q1

      real(kind=kind_phys), dimension(im,km),   intent(inout) :: du, dv,
     &  tau
      real(kind=kind_phys), dimension(im,km,ntrac), intent(inout) :: rtg

      real(kind=kind_phys), dimension(:,:),     intent(inout) :: 
     &  du3dt_PBL, dv3dt_PBL, dt3dt_PBL, dq3dt_PBL, do3dt_PBL
      logical,                                  intent(in) :: ldiag3d, 
     &  qdiag3d, gen_tend

      integer, dimension(im),                   intent(out) :: kpbl
      real(kind=kind_phys), dimension(im),      intent(out) :: dusfc,
     &  dvsfc, dtsfc, dqsfc, hpbl
      real(kind=kind_phys), dimension(im,km),   intent(out) :: prnum
      real(kind=kind_phys), dimension(im,km-1), intent(out) :: dkt

      character(len=*),                         intent(out) :: errmsg
      integer,                                  intent(out) :: errflg
!
!    locals
!
      integer, parameter :: kp = kind_phys
      integer i,is,k,kk,km1,kmpbl,kp1, ntloc
!
      logical  pblflg(im), sfcflg(im), flg(im)

      real(kind=kind_phys), dimension(im) ::  phih, phim
     &,                     rbdn, rbup, sflux, z0, crb, zol, thermal
     &,                     beta, tx1
!
      real(kind=kind_phys), dimension(im,km)  :: theta, thvx, zl, a1, ad
     &,                                          dt2odel
      real(kind=kind_phys), dimension(im,km-1):: xkzo, xkzmo, al, au
     &,                                          dku, rdzt
!
      real(kind=kind_phys) zi(im,km+1), a2(im,km*(ntrac+1))
!
      real(kind=kind_phys) dsdz2,  dsdzq,  dsdzt, dsig, dt2, rdt
     &,                    dtodsd, dtodsu, rdz,   tem,  tem1
     &,                    ttend,  utend, vtend,  qtend
     &,                    spdk2,  rbint, ri,     zol1, robn, bvf2
!
      real(kind=kind_phys), parameter ::  one=1.0_kp, zero=0.0_kp
     &,              zolcr=0.2_kp,
     &               zolcru=-0.5_kp, rimin=-100.0_kp, sfcfrac=0.1_kp,
     &               crbcon=0.25_kp, crbmin=0.15_kp,  crbmax=0.35_kp,
     &               qmin=1.0e-8_kp, zfmin=1.0d-8,    qlmin=1.0e-12_kp,
     &               aphi5=5.0_kp,   aphi16=16.0_kp,  f0=1.0e-4_kp
     &,              dkmin=zero,     dkmax=1000.0_kp
!    &,              dkmin=zero,     dkmax=1000.,     xkzminv=0.3
     &,              prmin=0.25_kp,  prmax=4.0_kp,    vk=0.4_kp,
     &               cfac=6.5_kp
      real(kind=kind_phys) :: gravi, cont, conq, gocp, go2

      gravi = one  / grav
      cont  = cp   * gravi
      conq  = hvap * gravi
      gocp  = grav / cp
      go2   = grav * 0.5_kp

! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0
!
!-----------------------------------------------------------------------
!
!     compute preliminary variables
!
      dt2   = delt
      rdt   = one / dt2
      km1   = km - 1
      kmpbl = km / 2
!
      rtg = zero
!
      do k=1,km
        do i=1,im
          zi(i,k)      = phii(i,k) * gravi
          zl(i,k)      = phil(i,k) * gravi
          dt2odel(i,k) = dt2 / del(i,k)
        enddo
      enddo
      do i=1,im
         zi(i,km+1) = phii(i,km+1) * gravi
      enddo
!
      do k = 1,km1
        do i=1,im
          rdzt(i,k)  = one / (zl(i,k+1) - zl(i,k))
          prnum(i,k) = one
        enddo
      enddo
!               Setup backgrond diffision
      do i=1,im
        prnum(i,km) = one
        tx1(i) = one / prsi(i,1)
      enddo
      do k = 1,km1
        do i=1,im
          xkzo(i,k)  = zero
          xkzmo(i,k) = zero
!         if (k < kinver(i)) then
          if (k <= kinver(i)) then
!    vertical background diffusivity for heat and momentum
            tem1       = one - prsi(i,k+1) * tx1(i)
            tem1       = min(one, exp(-tem1 * tem1 * 10.0_kp))
            xkzo(i,k)  = xkzm_h * tem1
            xkzmo(i,k) = xkzm_m * tem1
          endif
        enddo
      enddo
!
!  diffusivity in the inversion layer is set to be xkzminv (m^2/s)
!
      do k = 1,kmpbl
        do i=1,im
          if(zi(i,k+1) > 250.0_kp) then
            tem1 = (t1(i,k+1)-t1(i,k)) * rdzt(i,k)
            if(tem1 > 1.0e-5_kp) then
               xkzo(i,k)  = min(xkzo(i,k),xkzminv)
            endif
          endif
        enddo
      enddo
!
!
      do i = 1,im
         z0(i)     = 0.01_kp * zorl(i)
         kpbl(i)   = 1
         hpbl(i)   = zi(i,1)
         pblflg(i) = .true.
         sfcflg(i) = .true.
         if(rbsoil(i) > zero) sfcflg(i) = .false.
         dusfc(i)  = zero
         dvsfc(i)  = zero
         dtsfc(i)  = zero
         dqsfc(i)  = zero
      enddo
!
      do k = 1,km
        do i=1,im
          tx1(i) = zero
        enddo
        do kk=1,ncnd
          do i=1,im
            tx1(i) = tx1(i) + max(q1(i,k,ntcw+kk-1), qlmin)
          enddo
        enddo
        do i = 1,im
          theta(i,k) = t1(i,k) * psk(i) / prslk(i,k)
          thvx(i,k)  = theta(i,k)*(one+fv*max(q1(i,k,1),qmin)-tx1(i))
        enddo
      enddo
!
      do i = 1,im
         sflux(i)  = heat(i) + evap(i)*fv*theta(i,1)
         if (.not.sfcflg(i) .or. sflux(i) <= zero) pblflg(i)=.false.
         beta(i)  = dt2 / (zi(i,2)-zi(i,1))
      enddo
!
!  compute the pbl height
!
!     write(0,*)' IN moninbl u10=',u10m(1:5),' v10=',v10m(1:5)
      do i=1,im
         flg(i)  = .false.
         rbup(i) = rbsoil(i)
!
         if (pblflg(i)) then
           thermal(i) = thvx(i,1)
           crb(i) = crbcon
         else
           thermal(i) = tsea(i)*(one+fv*max(q1(i,1,1),qmin))
           tem   = max(one, sqrt(u10m(i)*u10m(i) + v10m(i)*v10m(i)))
           robn   = tem / (f0 * z0(i))
           tem1   = 1.0e-7_kp * robn
           crb(i) = max(min(0.16_kp * (tem1 ** (-0.18_kp)), crbmax),
     &                                                      crbmin)
         endif
      enddo
      do k = 1, kmpbl
        do i = 1, im
          if (.not.flg(i)) then
            rbdn(i) = rbup(i)
            spdk2   = max((u1(i,k)*u1(i,k)+v1(i,k)*v1(i,k)), one)
            rbup(i) = (thvx(i,k)-thermal(i))*phil(i,k)
     &              / (thvx(i,1)*spdk2)
            kpbl(i) = k
            flg(i)  = rbup(i) > crb(i)
          endif
        enddo
      enddo
      do i = 1,im
        if(kpbl(i) > 1) then
          k = kpbl(i)
          if (rbdn(i) >= crb(i)) then
            rbint = zero
          elseif(rbup(i) <= crb(i)) then
            rbint = one
          else
            rbint = (crb(i)-rbdn(i)) / (rbup(i)-rbdn(i))
          endif
          hpbl(i) = zl(i,k-1) + rbint*(zl(i,k)-zl(i,k-1))
          if(hpbl(i) < zi(i,kpbl(i))) kpbl(i) = kpbl(i) - 1
        else
          hpbl(i) = zl(i,1)
          kpbl(i) = 1
        endif
      enddo
!
!  compute similarity parameters
!
      do i=1,im
         zol(i) = max(rbsoil(i)*fm(i)*fm(i)/fh(i),rimin)
         if (sfcflg(i)) then
           zol(i) = min(zol(i),-zfmin)
         else
           zol(i) = max(zol(i),zfmin)
         endif
         zol1 = zol(i)*sfcfrac*hpbl(i)/zl(i,1)
         if (sfcflg(i)) then
!          phim(i) = (1.-aphi16*zol1)**(-1./4.)
!          phih(i) = (1.-aphi16*zol1)**(-1./2.)
           tem     = one / max(one - aphi16*zol1, 1.0e-8_kp)
           phih(i) = sqrt(tem)
           phim(i) = sqrt(phih(i))
         else
           phim(i) = one + aphi5*zol1
           phih(i) = phim(i)
         endif
      enddo
!
!  enhance the pbl height by considering the thermal excess
!
      do i=1,im
         flg(i)  = .true.
         if (pblflg(i)) then
           flg(i) = .false.
           rbup(i) = rbsoil(i)
         endif
      enddo
      do k = 2, kmpbl
        do i = 1, im
          if (.not.flg(i)) then
            rbdn(i) = rbup(i)
            spdk2   = max((u1(i,k)*u1(i,k)+v1(i,k)*v1(i,k)), one)
            rbup(i) = (thvx(i,k)-thermal(i)) * phil(i,k)
     &              / (thvx(i,1)*spdk2)
            kpbl(i) = k
            flg(i)  = rbup(i) > crb(i)
          endif
        enddo
      enddo
      do i = 1,im
        if (pblflg(i)) then
          k = kpbl(i)
          if(rbdn(i) >= crb(i)) then
            rbint = zero
          elseif(rbup(i) <= crb(i)) then
            rbint = one
          else
            rbint = (crb(i)-rbdn(i)) / (rbup(i)-rbdn(i))
          endif
          if (k > 1) then
            hpbl(i) = zl(i,k-1) + rbint*(zl(i,k)-zl(i,k-1))
            if(hpbl(i) < zi(i,kpbl(i))) kpbl(i) = kpbl(i) - 1
            if(kpbl(i) <= 1) then
              pblflg(i) = .false.
            endif
          else
            pblflg(i) = .false.
          endif
        endif
        if (pblflg(i)) then
          tem = phih(i)/phim(i) + cfac*vk*sfcfrac
        else
          tem = phih(i)/phim(i)
        endif
        prnum(i,1) = min(prmin,max(prmax,tem))
      enddo
!
      do i = 1, im
        if(zol(i) > zolcr) then
          kpbl(i) = 1
        endif
      enddo
!
!  compute Prandtl number above boundary layer
!
      do k = 1, km1
        kp1 = k + 1
        do i=1,im
          if(k >= kpbl(i)) then
            rdz  = rdzt(i,k)
            tem  = u1(i,k) - u1(i,kp1)
            tem1 = v1(i,k) - v1(i,kp1)
            tem  = (tem*tem + tem1*tem1) * rdz * rdz
            bvf2 = go2*(thvx(i,kp1)-thvx(i,k))*rdz / (t1(i,k)+t1(i,kp1))
            ri   = max(bvf2/tem,rimin)
            if(ri < zero) then ! unstable regime
              prnum(i,kp1) = one
            else
              prnum(i,kp1) = min(one + 2.1_kp*ri, prmax)
            endif
          elseif (k > 1) then
            prnum(i,kp1) = prnum(i,1)
          endif
!
!         prnum(i,kp1) = one
          prnum(i,kp1) = max(prmin, min(prmax, prnum(i,kp1)))
          tem      = tkh(i,kp1) * prnum(i,kp1)
          dku(i,k) = max(min(tem+xkzmo(i,k),       dkmax), xkzmo(i,k))
          dkt(i,k) = max(min(tkh(i,kp1)+xkzo(i,k), dkmax), xkzo(i,k))
        enddo
      enddo
!
!     compute tridiagonal matrix elements for heat and moisture
!
      do i=1,im
         ad(i,1) = one
         a1(i,1) = t1(i,1)   + beta(i) * heat(i)
         a2(i,1) = q1(i,1,1) + beta(i) * evap(i)
      enddo

      ntloc = 1
      if(ntrac > 1) then
        is    = 0
        do k = 2, ntrac
          if (k /= ntke) then
            ntloc = ntloc + 1
            is = is + km
            do i = 1, im
              a2(i,1+is) = q1(i,1,k)
            enddo
          endif
        enddo
      endif
!
      do k = 1,km1
        kp1 = k + 1
        do i = 1,im
          dtodsd    = dt2odel(i,k)
          dtodsu    = dt2odel(i,kp1)
          dsig      = prsl(i,k)-prsl(i,kp1)
          rdz       = rdzt(i,k)
          tem1      = dsig * dkt(i,k) * rdz
          dsdz2     = tem1 * rdz
          au(i,k)   = -dtodsd*dsdz2
          al(i,k)   = -dtodsu*dsdz2
!
          ad(i,k)   = ad(i,k)-au(i,k)
          ad(i,kp1) = one - al(i,k)
          dsdzt     = tem1 * gocp
          a1(i,k)   = a1(i,k)   + dtodsd*dsdzt
          a1(i,kp1) = t1(i,kp1) - dtodsu*dsdzt
          a2(i,kp1) = q1(i,kp1,1)
!
        enddo
      enddo
!
      if(ntrac > 1) then
        is = 0
        do kk = 2, ntrac
          if (kk /= ntke) then
            is = is + km
            do k = 1, km1
              kp1 = k + 1
              do i = 1, im
                a2(i,kp1+is) = q1(i,kp1,kk)
              enddo
            enddo
          endif
        enddo
      endif
!
!     solve tridiagonal problem for heat, moisture and tracers
!
      call tridin(im,km,ntloc,al,ad,au,a1,a2,au,a1,a2)

!
!     recover tendencies of heat and moisture
!
      do  k = 1,km
        do i = 1,im
          ttend      = (a1(i,k)-t1(i,k))   * rdt
          qtend      = (a2(i,k)-q1(i,k,1)) * rdt
          tau(i,k)   = tau(i,k)   + ttend
          rtg(i,k,1) = rtg(i,k,1) + qtend
          dtsfc(i)   = dtsfc(i)   + del(i,k)*ttend
          dqsfc(i)   = dqsfc(i)   + del(i,k)*qtend
        enddo
      enddo
      if(ldiag3d .and. .not. gen_tend) then
        do  k = 1,km
          do i = 1,im
            ttend      = (a1(i,k)-t1(i,k))
            dt3dt_PBL(i,k) = dt3dt_PBL(i,k) + ttend
          enddo
        enddo
        if(qdiag3d) then
          do  k = 1,km
            do i = 1,im
              qtend      = (a2(i,k)-q1(i,k,1))
              dq3dt_PBL(i,k) = dq3dt_PBL(i,k) + qtend
            enddo
          enddo
        endif
      endif
      do i = 1,im
        dtsfc(i)   = dtsfc(i) * cont
        dqsfc(i)   = dqsfc(i) * conq
      enddo
      if(ntrac > 1) then
        is = 0
        do kk = 2, ntrac
          if (kk /= ntke) then
            is = is + km
            do k = 1, km
              do i = 1, im
                qtend = (a2(i,k+is)-q1(i,k,kk))*rdt
                rtg(i,k,kk) = rtg(i,k,kk) + qtend
              enddo
            enddo
          endif
        enddo
        if(ldiag3d .and. ntoz>0 .and. qdiag3d .and. .not. gen_tend) then
          kk = ntoz
          is = (kk-1) * km
          do k = 1, km
            do i = 1, im
              qtend = (a2(i,k+is)-q1(i,k,kk))
              do3dt_PBL(i,k) = do3dt_PBL(i,k)+qtend
            enddo
          enddo
        endif
      endif
!
!     compute tridiagonal matrix elements for momentum
!
      do i=1,im
         ad(i,1) = one + beta(i) * stress(i) / spd1(i)
         a1(i,1) = u1(i,1)
         a2(i,1) = v1(i,1)
      enddo
!
      do k = 1,km1
        kp1 = k + 1
        do i=1,im
          dtodsd    = dt2odel(i,k)
          dtodsu    = dt2odel(i,kp1)
          dsig      = prsl(i,k)-prsl(i,kp1)
          rdz       = rdzt(i,k)
          tem1      = dsig*dku(i,k)*rdz
          dsdz2     = tem1 * rdz
          au(i,k)   = -dtodsd*dsdz2
          al(i,k)   = -dtodsu*dsdz2
!
          ad(i,k)   = ad(i,k) - au(i,k)
          ad(i,kp1) = one - al(i,k)
          a1(i,kp1) = u1(i,kp1)
          a2(i,kp1) = v1(i,kp1)
!
        enddo
      enddo

      call tridi2(im,km,al,ad,au,a1,a2,au,a1,a2)
!
!     recover tendencies of momentum
!
      do k = 1,km
        do i = 1,im
          utend    = (a1(i,k)-u1(i,k))*rdt
          vtend    = (a2(i,k)-v1(i,k))*rdt
          du(i,k)  = du(i,k)  + utend
          dv(i,k)  = dv(i,k)  + vtend
          tem      = del(i,k) * gravi
          dusfc(i) = dusfc(i) + tem * utend
          dvsfc(i) = dvsfc(i) + tem * vtend
        enddo
      enddo
      if (ldiag3d .and. .not. gen_tend) then
        do k = 1,km
          do i = 1,im
            utend    = (a1(i,k)-u1(i,k))
            vtend    = (a2(i,k)-v1(i,k))
            du3dt_PBL(i,k) = du3dt_PBL(i,k) + utend
            dv3dt_PBL(i,k) = dv3dt_PBL(i,k) + vtend
          enddo
        enddo
      endif
!
      if (ntke > 0) then    ! solve tridiagonal problem for momentum and tke
!
!     compute tridiagonal matrix elements for tke
!
        do i=1,im
           ad(i,1) = one
           a1(i,1) = q1(i,1,ntke)
        enddo
!
        do k = 1,km1
          kp1 = k + 1
          do i=1,im
            dtodsd    = dt2odel(i,k)
            dtodsu    = dt2odel(i,kp1)
            dsig      = prsl(i,k)-prsl(i,kp1)
            rdz       = rdzt(i,k)
            tem1      = dsig*dku(i,k)*(rdz+rdz)
            dsdz2     = tem1 * rdz
            au(i,k)   = -dtodsd*dsdz2
            al(i,k)   = -dtodsu*dsdz2
!
            ad(i,k)   = ad(i,k) - au(i,k)
            ad(i,kp1) = one - al(i,k)
            a1(i,kp1) = q1(i,kp1,ntke)
          enddo
        enddo

        call tridi1(im,km,al,ad,au,a1,au,a1)
!
        do k = 1, km !     recover tendencies of tke
          do i = 1, im
            qtend = (a1(i,k)-q1(i,k,ntke))*rdt
            rtg(i,k,ntke) = rtg(i,k,ntke) + qtend
          enddo
        enddo
      endif
!
      return
      end subroutine moninshoc_run

      end module moninshoc

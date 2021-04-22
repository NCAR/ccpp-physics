module cires_orowam2017


contains


      subroutine oro_wam_2017(im, levs,npt,ipt, kref,kdt,me,master, &
     &   dtp,dxres, taub, u1, v1, t1, xn, yn, bn2, rho, prsi, prsL, &
     &   grav, omega, con_rd, del, sigma, hprime, gamma, theta,  &
     &   sinlat, xlatd, taup, taud, pkdis)
! 
      USE MACHINE ,      ONLY : kind_phys
!      
      implicit none

      integer :: im, levs
      integer :: npt
      integer :: kdt, me, master
      integer :: kref(im), ipt(im)
      real(kind=kind_phys), intent(in) :: dtp, dxres
      real(kind=kind_phys), intent(in) :: taub(im)

      real(kind=kind_phys), intent(in) :: sinlat(im), xlatd(im)
      real(kind=kind_phys), intent(in), dimension(im) :: sigma,  &
     &                      hprime, gamma, theta

      real(kind=kind_phys), intent(in), dimension(im) :: xn, yn

      real(kind=kind_phys), intent(in), dimension(im, levs) ::   &
     &   u1, v1, t1,  bn2,  rho,   prsl, del
      real(kind=kind_phys), intent(in) :: grav, omega, con_rd
 
      real(kind=kind_phys), intent(in), dimension(im, levs+1) :: prsi
!
! out : taup,  taud, pkdis
!
      real(kind=kind_phys), intent(inout), dimension(im, levs+1) :: taup
      real(kind=kind_phys), intent(inout), dimension(im, levs)  :: taud
      real(kind=kind_phys), intent(inout), dimension(im, levs)  :: pkdis
      real(kind=kind_phys)     :: belps, aelps, nhills, selps
!
! multiwave oro-spectra
! locals
!
      integer :: i, j, k, isp, iw

      integer, parameter               :: nworo = 30
      real(kind=kind_phys), parameter  :: fc_flag = 0.0
      real(kind=kind_phys), parameter  :: mkzmin = 6.28e-3/50.0
      real(kind=kind_phys), parameter  :: mkz2min = mkzmin* mkzmin
      real(kind=kind_phys), parameter  :: kedmin = 1.e-3
      real(kind=kind_phys), parameter  :: kedmax = 350.,axmax=250.e-5
      real(kind=kind_phys), parameter  :: rtau   = 0.01   ! nonlin-OGW scale 1/10sec
      real(kind=kind_phys), parameter  :: Linsat2 =0.5
      real(kind=kind_phys), parameter  :: kxmin = 6.28e-3/100.
      real(kind=kind_phys), parameter  :: kxmax = 6.28e-3/5.0
      real(kind=kind_phys), parameter  :: dkx = (kxmax -kxmin)/(nworo-1)
      real(kind=kind_phys), parameter  :: kx_slope= -5./3.
      real(kind=kind_phys), parameter  :: hps =7000., rhp2 = .5/hps
      real(kind=kind_phys), parameter  :: cxmin=0.5, cxmin2=cxmin*cxmin

      real  :: akx(nworo),  cxoro(nworo), akx2(nworo)
      real  :: aspkx(nworo), c2f2(nworo) , cdf2(nworo)
      real  :: tau_sp(nworo,levs+1), wkdis(nworo, levs+1)
      real  :: tau_kx(nworo),taub_kx(nworo)
      real, dimension(nworo, levs+1)  :: wrms, akzw

      real  :: tauz(levs+1), rms_wind(levs+1)
      real  :: wave_act(nworo,levs+1)

      real  :: kxw, kzw, kzw2, kzw3, kzi, dzmet, rhoint
      real  ::  rayf, kturb
      real  ::  uz, bv, bv2,kxsp, fcor2, cf2

      real  ::  fdis
      real  ::  wfdm, wfdt, wfim, wfit
      real  ::  betadis, betam, betat, kds, cx, rhofac
      real  ::  etwk, etws, tauk, cx2sat
      real  ::  cdf1, tau_norm
!
! mean flow
!
      real, dimension(levs+1) :: uzi,rhoi,ktur, kalp, dzi

      integer :: nw, nzi, ksrc
      taud (:, :) = 0.0 ; pkdis(:,:) = 0.0 ; taup (:,:) = 0.0
      tau_sp (:,:) = 0.0 ; wrms(:,:) = 0.0
      nw  =  nworo
      nzi = levs+1

      do iw = 1, nw
!	                              !kxw =  0.25/(dxres)*iw
        kxw        = kxmin+(iw-1)*dkx
        akx(iw)    = kxw
        akx2(iw)   = kxw*kxw
        aspkx(iw)  = kxw ** (kx_slope)
        tau_kx(iw) = aspkx(iw)*dkx
      enddo

      tau_norm  = sum(tau_kx)
      tau_kx(:) =  tau_kx(:)/tau_norm

      if (kdt == 1) then
771     format( 'vay-oro19 ', 3(2x,F8.3))
        write(6,771)                         &
     &    maxval(tau_kx)*maxval(taub)*1.e3,  &
     &    minval(tau_kx), maxval(tau_kx)
      endif
!
! main loop over oro-points
!
      do i =1, npt
         j = ipt(i)

!
! estimate "nhills" => stochastic choices for OGWs
!
         if (taub(i) > 0.) then
!
! max_kxridge =min( .5*sigma(j)/hprime(j), kmax)
! ridge-dependent dkx = (max_kxridge -kxmin)/(nw-1)
! option to make grid-box variable kx-spectra kxw = kxmin+(iw-1)*dkx
!
           wave_act(1:nw, 1:levs+1) = 1.0
           ksrc = kref(i)
           tauz(1:ksrc)  = taub(i)
           taub_kx(1:nw) = tau_kx(1:nw) * taub(i)
           wkdis(:,:)    = kedmin

           call oro_meanflow(levs, nzi, u1(j,:), v1(j,:), t1(j,:),      &
     &                       prsi(j,:), prsL(j,:), grav, con_rd,        &
     &                       del(j,:), rho(i,:),                        &
     &                       bn2(i,:), uzi, rhoi,ktur, kalp,dzi,        &
     &                       xn(i), yn(i))

           fcor2 = (2*omega*sinlat(j))*(2*omega*sinlat(j))*fc_flag

           k = ksrc

           bv2    = bn2(i,k)
           uz     = uzi(k)           !u1(j,ksrc)*xn(i)+v1(j,ksrc)*yn(i)!
           kturb  = ktur(k)
           rayf   = kalp(k)
           rhoint = rhoi(k)
           dzmet  = dzi(k)
           kzw    = max(sqrt(bv2)/max(cxmin, uz), mkzmin)
!
! specify oro-kx spectra and related variables k=ksrc
!
           do iw = 1, nw
             kxw =  akx(iw)
             cxoro(iw) = 0.0 - uz
             c2f2(iw) = fcor2/akx2(iw)
             wrms(iw,k)= taub_kx(iw)/rhoint*kzw/kxw
             tau_sp(iw, k) = taub_kx(iw)
!
!
             if (cxoro(iw) > cxmin) then
               wave_act(iw,k:levs+1) = 0.                 ! crit-level
             else
               cdf2(iw) =  cxoro(iw)*cxoro(iw) -c2f2(iw)
               if ( cdf2(iw) < cxmin2)  then
                 wave_act(iw,k:levs+1) = 0.               ! coriolis cut-off
               else
                 kzw2 = max(Bv2/Cdf2(iw) - akx2(iw), mkz2min)
                 kzw = sqrt(kzw2)
                 akzw(iw,k)= kzw
                 wrms(iw,k)= taub_kx(iw)/rhoint * kzw/kxw
               endif
             endif
           enddo             !  nw-spectral loop
!
! defined abobe, k = ksrc: akx(nworo),  cxoro(nworo), tau_sp(ksrc, nworo)
! propagate upward multiwave-spectra are filtered by dissipation & instability
!
!          tau_sp(:,ksrc+1:levs+1) = tau_sp(:, ksrc)
           do k= ksrc+1, levs
             uz = uzi(k)
             bv2 =bn2(i,k)
             bv = sqrt(bv2)
             rayf = kalp(k)
             rhoint= rhoi(k)
             dzmet = dzi(k)
             rhofac = rhoi(k-1)/rhoi(k)

             do iw = 1, nworo
!
               if (wave_act(iw, k-1) <= 0.0) cycle
               cxoro(iw)= 0.0 - uz
               if ( cxoro(iw) > cxmin) then
                 wave_act(iw,k:levs+1) = 0.0     ! crit-level
               else
                 cdf2(iw) =  cxoro(iw)*cxoro(iw) -c2f2(iw)
                 if ( cdf2(iw) < cxmin2)  wave_act(iw,k:levs+1) = 0.0
               endif       
               if ( wave_act(iw,k) <= 0.0) cycle 
!
! upward propagation
!          
               kzw2 = Bv2/Cdf2(iw) - akx2(iw) 

               if (kzw2 < mkz2min) then 
                 wave_act(iw,k:levs+1) = 0.0
               else  
!
! upward propagation w/o reflection
!
                  kxw        = akx(iw)
                  kzw        = sqrt(kzw2)
                  akzw(iw,k) = kzw
                  kzw3       = kzw2*kzw

                  cx = cxoro(iw)
                  betadis = cdf2(iw) / (Cx*Cx+c2f2(iw))
                  betaM   = 1.0 / (1.0+betadis)
                  betaT   = 1.0 - BetaM
                  kds     = wkdis(iw,k-1)

                  etws  =  wrms(iw,k-1)*rhofac * kzw/akzw(iw,k-1)

                  kturb = ktur(k)+pkdis(j,k-1)
                  wfiM  = kturb*kzw2 +rayf
                  wfiT  = wfiM                   ! do updates with Pr-numbers Kv/Kt
                  cdf1  = sqrt(Cdf2(iw))
                  wfdM  = wfiM/(kxw*Cdf1)*BetaM
                  wfdT  = wfiT/(kxw*Cdf1)*BetaT
                  kzi   = 2.*kzw*(wfdM+wfdT)*dzmet
                  Fdis  = exp(-kzi)

                  etwk   = etws*Fdis
                  Cx2sat = Linsat2*Cdf2(iw)

                 if (etwk > cx2sat) then
                    Kds  =  kxw*Cdf1*rhp2/kzw3
                    etwk = cx2sat
                    wfiM = kds*kzw2
                    wfdM = wfiM/(kxw*Cdf1)
                    kzi  = 2.*kzw*(wfdm + wfdm)*dzmet
                    etwk = cx2sat*exp(-kzi)
                  endif
!	   	                         if( lat(j) eq 40.5 ) then stop
                  wkdis(iw,k) = kds
                  wrms(iw,k) =  etwk
                  tauk = etwk*kxw/kzw
                  tau_sp(iw,k) = tauk *rhoint
                  if ( tau_sp(iw,k) > tau_sp(iw,k-1))        &
     &                                tau_sp(iw,k) = tau_sp(iw,k-1)

               ENDIF  ! upward
             ENDDO    ! spectral

!......... do spectral sum of rms, wkdis, tau

             tauz(k)     = sum( tau_sp(:,k)*wave_act(:,k) )
             rms_wind(k) = sum(   wrms(:,k)*wave_act(:,k) )

             pkdis(j,k) = sum(wkdis(:,k)*wave_act(:,k))+rms_wind(k)*rtau

             if (pkdis(j,k) > kedmax) pkdis(j,k) = kedmax

          ENDDO   ! k=ksrc+1, levs

          k = ksrc
          tauz(k)      = sum(tau_sp(:,k)*wave_act(:,k))
          tauz(k)      = tauz(k+1)   ! zero momentum dep-n at k=ksrc 

          pkdis(j,k)   = sum(wkdis(:,k)*wave_act(:,k))
          rms_wind(k)  = sum(wrms(:,k)*wave_act(:,k))
          tauz(levs+1) = tauz(levs)
          taup(i, 1:levs+1) = tauz(1:levs+1)
          do  k=ksrc, levs
            taud(i,k) = ( tauz(k+1) - tauz(k))*grav/del(j,k)
!	    if (taud(i,k) .gt. 0)taud(i,k)=taud(i,k)*.01
!	    if (abs(taud(i,k)).ge.axmax)taud(i,k)=sign(taud(i,k),axmax)
          enddo
        endif                  ! taub > 0
      enddo                    ! oro-points (i, j, ipt)
!23456
      end subroutine oro_wam_2017
!-------------------------------------------------------------
!
! define mean flow  and dissipation for OGW-kx spectrum
!
!-------------------------------------------------------------      
      subroutine oro_meanflow(nz, nzi, u1, v1, t1, pint, pmid,       &
     &           grav, con_rd,                                       &
     &           delp, rho, bn2, uzi, rhoi, ktur, kalp, dzi, xn, yn)

      use ugwp_common_v1 ,  only : velmin, dw2min
      implicit none

      integer :: nz, nzi
      real, dimension(nz  ) ::  u1,  v1, t1, delp, rho, pmid
      real, dimension(nz  ) ::  bn2  ! define at the interfaces
      real, dimension(nz+1) ::  pint
      real                  ::  xn, yn
      real,intent(in) :: grav, con_rd
! output
 
      real, dimension(nz+1) ::  dzi,  uzi, rhoi, ktur, kalp

! locals
      integer :: i, j, k
      real :: ui, vi, ti, uz, vz, shr2, rdz, kamp
      real :: zgrow, zmet, rdpm, ritur, kmol, w1
      real :: rgrav, rdi
! paremeters
      real, parameter :: hps = 7000., rpspa = 1.e-5
      real, parameter :: rhps=1.0/hps
      real, parameter :: h4= 0.25/hps
      real, parameter :: rimin = 1.0/8.0, kedmin = 0.01
      real, parameter :: lturb = 30. ,    uturb = 150.0
      real, parameter :: lsc2 = lturb*lturb,usc2 = uturb*uturb
      kalp(1:nzi) = 2.e-7        ! radiative damping

      rgrav = 1.0/grav
      rdi = 1.0/con_rd

      do k=2, nz
        rdpm    = grav/(pmid(k-1)-pmid(k))
        ui      = .5*(u1(k-1)+u1(k))
        vi      = .5*(v1(k-1)+v1(k))
        uzi(k)  = Ui*xn + Vi*yn
        ti      = .5*(t1(k-1)+t1(k))
        rhoi(k) = rdi*pint(k)/ti
        rdz     = rdpm *rhoi(k)
        dzi(k)  = 1./rdz
        uz      = u1(k)-u1(k-1)
        vz      = v1(k)-v1(k-1)
        shr2    = rdz*rdz*(max(uz*uz+vz*vz, dw2min))
        zmet    = -hps*alog(pint(k)*rpspa)
        zgrow   = exp(zmet*h4)
        kmol    =  2.e-5*exp(zmet*rhps)+kedmin
        ritur   = max(bn2(k)/shr2, rimin)
        kamp    = sqrt(shr2)*lsc2 *zgrow
        w1      = 1./(1. + 5*ritur)
        ktur(k) = kamp * w1 * w1 +kmol
      enddo

      k = 1 
      uzi(k)  = uzi(k+1)
      ktur(k) = ktur(k+1)
      rhoi(k) = rdi*pint(k)/t1(k+1)
      dzi(k)  = rgrav*delp(k)/rhoi(k)
      
      k = nzi
      uzi(k)  = uzi(k-1)
      ktur(k) = ktur(k-1)
      rhoi(k) = rhoi(k-1)*.5
      dzi(k)  =  dzi(k-1)

      end subroutine oro_meanflow

end module cires_orowam2017

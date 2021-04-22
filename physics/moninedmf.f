!> \file moninedmf.f
!!  Contains most of the hybrid eddy-diffusivity mass-flux scheme except for the
!!  subroutine that calculates the mass flux and updraft properties.

!> This module contains the CCPP-compliant hybrid eddy-diffusivity mass-flux
!! scheme.
      module hedmf

      contains

!> \section arg_table_hedmf_init Argument Table
!! \htmlinclude hedmf_init.html
!!
      subroutine hedmf_init (moninq_fac,errmsg,errflg)
         use machine, only : kind_phys
         implicit none
         real(kind=kind_phys), intent(in ) :: moninq_fac
         character(len=*),     intent(out) :: errmsg
         integer,              intent(out) :: errflg
         ! Initialize CCPP error handling variables
         errmsg = ''
         errflg = 0

         if (moninq_fac == 0) then
             errflg = 1
             write(errmsg,'(*(a))') 'Logic error: moninq_fac == 0',
     &                              ' is incompatible with hedmf'
         end if
      end subroutine hedmf_init

      subroutine hedmf_finalize ()
      end subroutine hedmf_finalize


!> \defgroup HEDMF GFS Hybrid Eddy-Diffusivity Mass-Flux (HEDMF) Scheme Module
!! @{
!!  \brief  This subroutine contains all of logic for the
!! Hybrid EDMF PBL scheme except for the calculation of
!! the updraft properties and mass flux.
!!
!> \section arg_table_hedmf_run Argument Table
!! \htmlinclude hedmf_run.html
!!
!!  \section general_edmf GFS Hybrid EDMF General Algorithm
!!  -# Compute preliminary variables from input arguments.
!!  -# Calculate the first estimate of the PBL height ("Predictor step").
!!  -# Calculate Monin-Obukhov similarity parameters.
!!  -# Update thermal properties of surface parcel and recompute PBL height ("Corrector step").
!!  -# Determine whether stratocumulus layers exist and compute quantities needed for enhanced diffusion.
!!  -# Calculate the inverse Prandtl number.
!!  -# Compute diffusion coefficients below the PBL top.
!!  -# Compute diffusion coefficients above the PBL top.
!!  -# If the PBL is convective, call the mass flux scheme to replace the countergradient terms.
!!  -# Compute enhanced diffusion coefficients related to stratocumulus-topped PBLs.
!!  -# Solve for the temperature and moisture tendencies due to vertical mixing.
!!  -# Calculate heating due to TKE dissipation and add to the tendency for temperature.
!!  -# Solve for the horizontal momentum tendencies and add them to output tendency terms.
!!  \section detailed_hedmf  GFS Hybrid HEDMF Detailed Algorithm
!!  @{
      subroutine hedmf_run (im,km,ntrac,ntcw,dv,du,tau,rtg,             &
     &   u1,v1,t1,q1,swh,hlw,xmu,                                       &
     &   psk,rbsoil,zorl,u10m,v10m,fm,fh,                               &
     &   tsea,heat,evap,stress,spd1,kpbl,                               &
     &   prsi,del,prsl,prslk,phii,phil,delt,dspheat,                    &
     &   dusfc,dvsfc,dtsfc,dqsfc,hpbl,hgamt,hgamq,dkt,dku,              &
     &   kinver,xkzm_m,xkzm_h,xkzm_s,lprnt,ipr,                         &
     &   xkzminv,moninq_fac,hurr_pbl,islimsk,var_ric,                   &
     &   coef_ric_l,coef_ric_s,lssav,ldiag3d,qdiag3d,ntoz,              &
     &   du3dt_PBL,dv3dt_PBL,dt3dt_PBL,dq3dt_PBL,do3dt_PBL,             &
     &   flag_for_pbl_generic_tend,errmsg,errflg)
!
      use machine  , only : kind_phys
      use funcphys , only : fpvs
      !GJF: Note that sending these constants through the argument list
      !results in regression test failures with "PROD" mode compilation
      !flags (specifically, grav and cp)
      use physcons, grav => con_g, cp => con_cp,
     &              hvap => con_hvap, fv => con_fvirt

      implicit none
!
!     arguments
!
      logical, intent(in) :: lprnt, hurr_pbl, lssav, ldiag3d, qdiag3d
      logical, intent(in) :: flag_for_pbl_generic_tend
      integer, intent(in) :: ipr, islimsk(im)
      integer, intent(in) :: im, km, ntrac, ntcw, kinver(im), ntoz
      integer, intent(out) :: kpbl(im)

!
      real(kind=kind_phys), intent(in) :: delt, xkzm_m, xkzm_h, xkzm_s
      real(kind=kind_phys), intent(in) :: xkzminv, moninq_fac, var_ric, &
     &                     coef_ric_l, coef_ric_s
      real(kind=kind_phys), intent(inout) :: dv(im,km),     du(im,km),  &
     &                     tau(im,km),    rtg(im,km,ntrac)
      ! Only allocated if ldiag3d or qdiag3d are true
      real(kind=kind_phys), intent(inout), dimension(:,:) ::            &
     &   du3dt_PBL,dv3dt_PBL,dt3dt_PBL,dq3dt_PBL,do3dt_PBL
      real(kind=kind_phys), intent(in) ::                               &
     &                     u1(im,km),     v1(im,km),                    &
     &                     t1(im,km),     q1(im,km,ntrac),              &
     &                     swh(im,km),    hlw(im,km),                   &
     &                     xmu(im),       psk(im),                      &
     &                     rbsoil(im),    zorl(im),                     &
     &                     u10m(im),      v10m(im),                     &
     &                     fm(im),        fh(im),                       &
     &                     tsea(im),                                    &
     &                     heat(im),      evap(im),                     &
     &                     stress(im),    spd1(im)
      real(kind=kind_phys), intent(in) ::                               &
     &                     prsi(im,km+1), del(im,km),                   &
     &                     prsl(im,km),   prslk(im,km),                 &
     &                     phii(im,km+1), phil(im,km)
      real(kind=kind_phys), intent(out) ::                              &
     &                     dusfc(im),     dvsfc(im),                    &
     &                     dtsfc(im),     dqsfc(im),                    &
     &                     hpbl(im)
      real(kind=kind_phys), intent(out) ::                              &
     &                     dkt(im,km-1), dku(im,km-1)
      real(kind=kind_phys), intent(inout) ::                            &
     &                     hgamt(im),     hgamq(im)
!
      logical, intent(in) :: dspheat
!          flag for tke dissipative heating
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg

!
!    locals

      integer i,iprt,is,iun,k,kk,km1,kmpbl,latd,lond
      integer lcld(im),icld(im),kcld(im),krad(im)
      integer kx1(im), kpblx(im)
!
!     real(kind=kind_phys) betaq(im), betat(im),   betaw(im),
      real(kind=kind_phys) phih(im), phim(im),  hpblx(im),              &
     &                     rbdn(im),    rbup(im),                       &
     &                     beta(im),    sflux(im),                      &
     &                     z0(im),    crb(im),     wstar(im),           &
     &                     zol(im),   ustmin(im),  ustar(im),           &
     &                     thermal(im),wscale(im), wscaleu(im)
!
      real(kind=kind_phys) theta(im,km),thvx(im,km),  thlvx(im,km),     &
     &                     qlx(im,km),  thetae(im,km),                  &
     &                     qtx(im,km),  bf(im,km-1),  diss(im,km),      &
     &                     radx(im,km-1),                               &
     &                     govrth(im),  hrad(im),                       &
!    &                     hradm(im),   radmin(im),   vrad(im),         &
     &                     radmin(im),  vrad(im),                       &
     &                     zd(im),      zdd(im),      thlvx1(im)
!
      real(kind=kind_phys) rdzt(im,km-1),dktx(im,km-1),                 &
     &                     zi(im,km+1),  zl(im,km),                     &
     &                     xkzo(im,km-1), xkzmo(im,km-1),               &
     &                     cku(im,km-1), ckt(im,km-1),                  &
     &                     ti(im,km-1),  shr2(im,km-1),                 &
     &                     al(im,km-1),  ad(im,km),                     &
     &                     au(im,km-1),  a1(im,km),                     &
     &                     a2(im,km*ntrac)
!
      real(kind=kind_phys) tcko(im,km),  qcko(im,km,ntrac),             &
     &                     ucko(im,km),  vcko(im,km),  xmf(im,km)
!
      real(kind=kind_phys) prinv(im), rent(im)
!
      logical  pblflg(im), sfcflg(im), scuflg(im), flg(im)
      logical  ublflg(im), pcnvflg(im)
!
!  pcnvflg: true for convective(strongly unstable) pbl
!  ublflg: true for unstable but not convective(strongly unstable) pbl
!
      real(kind=kind_phys) aphi16,  aphi5,  bvf2,   wfac,
     &                     cfac,    conq,   cont, conw,
     &                     dk,      dkmax,  dkmin,
     &                     dq1,     dsdz2,  dsdzq,  dsdzt,
     &                     dsdzu,   dsdzv,
     &                     dsig,    dt2,    dthe1,  dtodsd,
     &                     dtodsu,  dw2,    dw2min,
     &                     gamcrq,  gamcrt, gocp,
     &                     gravi,   f0,
     &                     prnum,   prmax,  prmin,  pfac,  crbcon,
     &                     qmin,    tdzmin, qtend,  crbmin,crbmax,
     &                     rbint,   rdt,    rdz,    qlmin,
     &                     ri,      rimin,  rl2,    rlam,  rlamun,
     &                     rone,    rzero,  sfcfrac,
     &                     spdk2,   sri,    zol1,   zolcr, zolcru,
     &                     robn,    ttend,
     &                     utend,   vk,     vk2,
     &                     ust3,    wst3,
     &                     vtend,   zfac,   vpert,  cteit,
     &                     rentf1,  rentf2, radfac,
     &                     zfmin,   zk,     tem,    tem1,  tem2,
     &                     xkzm,    xkzmu,
     &                     ptem,    ptem1,  ptem2, tx1(im), tx2(im)
!
      real(kind=kind_phys) zstblmax,h1,     h2,     qlcr,  actei,
     &                     cldtime
      real :: ttend_fac
      
      !! for hurricane application
      real(kind=kind_phys) wspm(im,km-1)
      integer kLOC ! RGF
      real :: xDKU ! RGF

      integer, parameter :: useshape=2!0-- no change, original ALPHA adjustment,1-- shape1, 2-- shape2(adjust above sfc)
      real :: smax,ashape,sz2h, sksfc,skmax,ashape1,skminusk0, hmax
cc
      parameter(gravi=1.0/grav)
      parameter(gocp=grav/cp)
      parameter(cont=cp/grav,conq=hvap/grav,conw=1.0/grav)               ! for del in pa
!     parameter(cont=1000.*cp/grav,conq=1000.*hvap/grav,conw=1000./grav) ! for del in kpa
      parameter(rlam=30.0,vk=0.4,vk2=vk*vk)
      parameter(prmin=0.25,prmax=4.,zolcr=0.2,zolcru=-0.5)
      parameter(dw2min=0.0001,dkmin=0.0,dkmax=1000.,rimin=-100.)
      parameter(crbcon=0.25,crbmin=0.15,crbmax=0.35)
      parameter(wfac=7.0,cfac=6.5,pfac=2.0,sfcfrac=0.1)
!     parameter(qmin=1.e-8,xkzm=1.0,zfmin=1.e-8,aphi5=5.,aphi16=16.)
      parameter(qmin=1.e-8,         zfmin=1.e-8,aphi5=5.,aphi16=16.)
      parameter(tdzmin=1.e-3,qlmin=1.e-12,f0=1.e-4)
      parameter(h1=0.33333333,h2=0.66666667)
!     parameter(cldtime=500.,xkzminv=0.3)
      parameter(cldtime=500.)
!     parameter(cldtime=500.,xkzmu=3.0,xkzminv=0.3)
!     parameter(gamcrt=3.,gamcrq=2.e-3,rlamun=150.0)
      parameter(gamcrt=3.,gamcrq=0.,rlamun=150.0)
      parameter(rentf1=0.2,rentf2=1.0,radfac=0.85)
      parameter(iun=84)
!
!     parameter (zstblmax = 2500., qlcr=1.0e-5)
!     parameter (zstblmax = 2500., qlcr=3.0e-5)
!     parameter (zstblmax = 2500., qlcr=3.5e-5)
!     parameter (zstblmax = 2500., qlcr=1.0e-4)
      parameter (zstblmax = 2500., qlcr=3.5e-5)
!     parameter (actei = 0.23)
      parameter (actei = 0.7)
c
c-----------------------------------------------------------------------
c
 601  format(1x,' moninp lat lon step hour ',3i6,f6.1)
 602      format(1x,'    k','        z','        t','       th',
     1     '      tvh','        q','        u','        v',
     2     '       sp')
 603      format(1x,i5,8f9.1)
 604      format(1x,'  sfc',9x,f9.1,18x,f9.1)
 605      format(1x,'    k      zl    spd2   thekv   the1v'
     1         ,' thermal    rbup')
 606      format(1x,i5,6f8.2)
 607      format(1x,' kpbl    hpbl      fm      fh   hgamt',
     1         '   hgamq      ws   ustar      cd      ch')
 608      format(1x,i5,9f8.2)
 609      format(1x,' k pr dkt dku ',i5,3f8.2)
 610      format(1x,' k pr dkt dku ',i5,3f8.2,' l2 ri t2',
     1         ' sr2  ',2f8.2,2e10.2)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

! compute preliminary variables
!
!     iprt = 0
!     if(iprt.eq.1) then
!cc   latd = 0
!     lond = 0
!     else
!cc   latd = 0
!     lond = 0
!     endif
!
      dt2   = delt
      rdt   = 1. / dt2
      km1   = km - 1
      kmpbl = km / 2
!>  - Compute physical height of the layer centers and interfaces from the geopotential height (zi and zl)
      do k=1,km
        do i=1,im
          zi(i,k) = phii(i,k) * gravi
          zl(i,k) = phil(i,k) * gravi
        enddo
      enddo
      do i=1,im
         zi(i,km+1) = phii(i,km+1) * gravi
      enddo
!>  - Compute reciprocal of \f$ \Delta z \f$ (rdzt)
      do k = 1,km1
        do i=1,im
          rdzt(i,k) = 1.0 / (zl(i,k+1) - zl(i,k))
        enddo
      enddo
!>  - Compute reciprocal of pressure (tx1, tx2)
      do i=1,im
        kx1(i) = 1
        tx1(i) = 1.0 / prsi(i,1)
        tx2(i) = tx1(i)
      enddo
!>  - Compute background vertical diffusivities for scalars and momentum (xkzo and xkzmo)
      do k = 1,km1
        do i=1,im
          xkzo(i,k)  = 0.0
          xkzmo(i,k) = 0.0
          if (k < kinver(i)) then
!                                  vertical background diffusivity
            ptem      = prsi(i,k+1) * tx1(i)
            tem1      = 1.0 - ptem
            tem1      = tem1 * tem1 * 10.0
            xkzo(i,k) = xkzm_h * min(1.0, exp(-tem1))

!                                  vertical background diffusivity for momentum
            if (ptem >= xkzm_s) then
              xkzmo(i,k) = xkzm_m
              kx1(i)     = k + 1
            else
              if (k == kx1(i) .and. k > 1) tx2(i) = 1.0 / prsi(i,k)
              tem1 = 1.0 - prsi(i,k+1) * tx2(i)
              tem1 = tem1 * tem1 * 5.0
              xkzmo(i,k) = xkzm_m * min(1.0, exp(-tem1))
            endif
          endif
        enddo
      enddo
!     if (lprnt) then
!       print *,' xkzo=',(xkzo(ipr,k),k=1,km1)
!       print *,' xkzmo=',(xkzmo(ipr,k),k=1,km1)
!     endif
!
! diffusivity in the inversion layer is set to be xkzminv (m^2/s)
!>  - The background scalar vertical diffusivity is limited to be less than or equal to xkzminv
      do k = 1,kmpbl
        do i=1,im
!         if(zi(i,k+1) > 200..and.zi(i,k+1) < zstblmax) then
          if(zi(i,k+1) > 250.) then
            tem1 = (t1(i,k+1)-t1(i,k)) * rdzt(i,k)
            if(tem1 > 1.e-5) then
               xkzo(i,k)  = min(xkzo(i,k),xkzminv)
            endif
          endif
        enddo
      enddo
!>  - Some output variables and logical flags are initialized
      do i = 1,im
         z0(i)    = 0.01 * zorl(i)
         dusfc(i) = 0.
         dvsfc(i) = 0.
         dtsfc(i) = 0.
         dqsfc(i) = 0.
         wscale(i)= 0.
         wscaleu(i)= 0.
         kpbl(i)  = 1
         hpbl(i)  = zi(i,1)
         hpblx(i) = zi(i,1)
         pblflg(i)= .true.
         sfcflg(i)= .true.
         if(rbsoil(i) > 0.) sfcflg(i) = .false.
         ublflg(i)= .false.
         pcnvflg(i)= .false.
         scuflg(i)= .true.
         if(scuflg(i)) then
           radmin(i)= 0.
           rent(i)  = rentf1
           hrad(i)  = zi(i,1)
!          hradm(i) = zi(i,1)
           krad(i)  = 1
           icld(i)  = 0
           lcld(i)  = km1
           kcld(i)  = km1
           zd(i)    = 0.
        endif
      enddo
!>  - Compute \f$\theta\f$ (theta), \f$q_l\f$ (qlx), \f$q_t\f$ (qtx), \f$\theta_e\f$ (thetae), \f$\theta_v\f$ (thvx), \f$\theta_{l,v}\f$ (thlvx)
      do k = 1,km
        do i = 1,im
          theta(i,k) = t1(i,k) * psk(i) / prslk(i,k)
          qlx(i,k)   = max(q1(i,k,ntcw),qlmin)
          qtx(i,k)   = max(q1(i,k,1),qmin)+qlx(i,k)
          ptem       = qlx(i,k)
          ptem1      = hvap*max(q1(i,k,1),qmin)/(cp*t1(i,k))
          thetae(i,k)= theta(i,k)*(1.+ptem1)
          thvx(i,k)  = theta(i,k)*(1.+fv*max(q1(i,k,1),qmin)-ptem)
          ptem2      = theta(i,k)-(hvap/cp)*ptem
          thlvx(i,k) = ptem2*(1.+fv*qtx(i,k))
        enddo
      enddo
!>  - Initialize diffusion coefficients to 0 and calculate the total radiative heating rate (dku, dkt, radx)
      do k = 1,km1
        do i = 1,im
          dku(i,k)  = 0.
          dkt(i,k)  = 0.
          dktx(i,k) = 0.
          cku(i,k)  = 0.
          ckt(i,k)  = 0.
          tem       = zi(i,k+1)-zi(i,k)
          radx(i,k) = tem*(swh(i,k)*xmu(i)+hlw(i,k))
        enddo
      enddo
!>  - Set lcld to first index above 2.5km
      do i=1,im
         flg(i)  = scuflg(i)
      enddo
      do k = 1, km1
        do i=1,im
          if(flg(i).and.zl(i,k) >= zstblmax) then
             lcld(i)=k
             flg(i)=.false.
          endif
      enddo
      enddo
!
!  compute virtual potential temp gradient (bf) and winshear square
!>  - Compute \f$\frac{\partial \theta_v}{\partial z}\f$ (bf) and the wind shear squared (shr2)
      do k = 1, km1
      do i = 1, im
         rdz  = rdzt(i,k)
         bf(i,k) = (thvx(i,k+1)-thvx(i,k))*rdz
         ti(i,k) = 2./(t1(i,k)+t1(i,k+1))
         dw2  = (u1(i,k)-u1(i,k+1))**2
     &        + (v1(i,k)-v1(i,k+1))**2
         shr2(i,k) = max(dw2,dw2min)*rdz*rdz
      enddo
      enddo
!>  - Calculate \f$\frac{g}{\theta}\f$ (govrth), \f$\beta = \frac{\Delta t}{\Delta z}\f$ (beta), \f$u_*\f$ (ustar), total surface flux (sflux), and set pblflag to false if the total surface energy flux is into the surface
      do i = 1,im
        govrth(i) = grav/theta(i,1)
      enddo
!
      do i=1,im
         beta(i)  = dt2 / (zi(i,2)-zi(i,1))
      enddo
!
      do i=1,im
         ustar(i) = sqrt(stress(i))
      enddo
!
      do i = 1,im
         sflux(i)  = heat(i) + evap(i)*fv*theta(i,1)
         if(.not.sfcflg(i) .or. sflux(i) <= 0.) pblflg(i)=.false.
      enddo
!>  ## Calculate the first estimate of the PBL height ("Predictor step")
!!  The calculation of the boundary layer height follows Troen and Mahrt (1986) \cite troen_and_mahrt_1986 section 3. The approach is to find the level in the column where a modified bulk Richardson number exceeds a critical value.
!!
!!  The temperature of the thermal is of primary importance. For the initial estimate of the PBL height, the thermal is assumed to have one of two temperatures. If the boundary layer is stable, the thermal is assumed to have a temperature equal to the surface virtual temperature. Otherwise, the thermal is assumed to have the same virtual potential temperature as the lowest model level. For the stable case, the critical bulk Richardson number becomes a function of the wind speed and roughness length, otherwise it is set to a tunable constant.
!  compute the pbl height
!
      if (.not. (hurr_pbl .and. moninq_fac < 0.0)) then
        do i=1,im
           flg(i) = .false.
           rbup(i) = rbsoil(i)
  !
           if(pblflg(i)) then
             thermal(i) = thvx(i,1)
             crb(i) = crbcon
           else
             thermal(i) = tsea(i)*(1.+fv*max(q1(i,1,1),qmin))
             tem = sqrt(u10m(i)**2+v10m(i)**2)
             tem = max(tem, 1.)
             robn = tem / (f0 * z0(i))
             tem1 = 1.e-7 * robn
             crb(i) = 0.16 * (tem1 ** (-0.18))
             crb(i) = max(min(crb(i), crbmax), crbmin)
           endif
        enddo
      else
        do i=1,im
          flg(i) = .false.
          rbup(i) = rbsoil(i)

          ! use variable Ri for all conditions
          if(pblflg(i)) then
            thermal(i) = thvx(i,1)
          else
            thermal(i) = tsea(i)*(1.+fv*max(q1(i,1,1),qmin))
          endif
          tem = sqrt(u10m(i)**2+v10m(i)**2)
          tem = max(tem, 1.)
          robn = tem / (f0 * z0(i))
          tem1 = 1.e-7 * robn
          crb(i) = crbcon
          if (var_ric .eq. 1.) then
            if (islimsk(i) .eq. 1)  crb(I) = coef_ric_l*(tem1)**(-0.18)
            if (islimsk(i) .eq. 0)  crb(I) = coef_ric_s*(tem1)**(-0.18)
          endif
          crb(i) = max(min(crb(i), crbmax), crbmin)
        enddo
      endif

!>  Given the thermal's properties and the critical Richardson number, a loop is executed to find the first level above the surface where the modified Richardson number is greater than the critical Richardson number, using equation 10a from Troen and Mahrt (1986) \cite troen_and_mahrt_1986 (also equation 8 from Hong and Pan (1996) \cite hong_and_pan_1996):
!!  \f[
!!  h = Ri\frac{T_0\left|\vec{v}(h)\right|^2}{g\left(\theta_v(h) - \theta_s\right)}
!!  \f]
!!  where \f$h\f$ is the PBL height, \f$Ri\f$ is the Richardson number, \f$T_0\f$ is the virtual potential temperature near the surface, \f$\left|\vec{v}\right|\f$ is the wind speed, and \f$\theta_s\f$ is for the thermal. Rearranging this equation to calculate the modified Richardson number at each level, k, for comparison with the critical value yields:
!!  \f[
!!  Ri_k = gz(k)\frac{\left(\theta_v(k) - \theta_s\right)}{\theta_v(1)*\vec{v}(k)}
!!  \f]
      do k = 1, kmpbl
      do i = 1, im
        if(.not.flg(i)) then
          rbdn(i) = rbup(i)
          spdk2   = max((u1(i,k)**2+v1(i,k)**2),1.)
          rbup(i) = (thvx(i,k)-thermal(i))*
     &              (grav*zl(i,k)/thvx(i,1))/spdk2
          kpbl(i) = k
          flg(i)  = rbup(i) > crb(i)
        endif
      enddo
      enddo
!>  Once the level is found, some linear interpolation is performed to find the exact height of the boundary layer top (where \f$Ri = Ri_{cr}\f$) and the PBL height and the PBL top index are saved (hpblx and kpblx, respectively)
      do i = 1,im
        if(kpbl(i) > 1) then
          k = kpbl(i)
          if(rbdn(i) >= crb(i)) then
            rbint = 0.
          elseif(rbup(i) <= crb(i)) then
            rbint = 1.
          else
            rbint = (crb(i)-rbdn(i))/(rbup(i)-rbdn(i))
          endif
          hpbl(i) = zl(i,k-1) + rbint*(zl(i,k)-zl(i,k-1))
          if(hpbl(i) < zi(i,kpbl(i))) kpbl(i) = kpbl(i) - 1
        else
          hpbl(i) = zl(i,1)
          kpbl(i) = 1
        endif
        kpblx(i) = kpbl(i)
        hpblx(i) = hpbl(i)
      enddo
!
!  compute similarity parameters
!>  ## Calculate Monin-Obukhov similarity parameters
!!  Using the initial guess for the PBL height, Monin-Obukhov similarity parameters are calculated. They are needed to refine the PBL height calculation and for calculating diffusion coefficients.
!!
!!  First, calculate the Monin-Obukhov nondimensional stability parameter, commonly referred to as \f$\zeta\f$ using the following equation from Businger et al. (1971) \cite businger_et_al_1971 (equation 28):
!!  \f[
!!  \zeta = Ri_{sfc}\frac{F_m^2}{F_h} = \frac{z}{L}
!!  \f]
!!  where \f$F_m\f$ and \f$F_h\f$ are surface Monin-Obukhov stability functions calculated in sfc_diff.f and \f$L\f$ is the Obukhov length. Then, the nondimensional gradients of momentum and temperature (phim and phih) are calculated using equations 5 and 6 from Hong and Pan (1996) \cite hong_and_pan_1996 depending on the surface layer stability. Then, the velocity scale valid for the surface layer (\f$w_s\f$, wscale) is calculated using equation 3 from Hong and Pan (1996) \cite hong_and_pan_1996. For the neutral and unstable PBL above the surface layer, the convective velocity scale, \f$w_*\f$, is calculated according to:
!!  \f[
!!  w_* = \left(\frac{g}{\theta_0}h\overline{w'\theta_0'}\right)^{1/3}
!!  \f]
!!  and the mixed layer velocity scale is then calculated with equation 6 from Troen and Mahrt (1986) \cite troen_and_mahrt_1986
!!  \f[
!!  w_s = (u_*^3 + 7\epsilon k w_*^3)^{1/3}
!!  \f]
      do i=1,im
         zol(i) = max(rbsoil(i)*fm(i)*fm(i)/fh(i),rimin)
         if(sfcflg(i)) then
           zol(i) = min(zol(i),-zfmin)
         else
           zol(i) = max(zol(i),zfmin)
         endif
         zol1 = zol(i)*sfcfrac*hpbl(i)/zl(i,1)
         if(sfcflg(i)) then
!          phim(i) = (1.-aphi16*zol1)**(-1./4.)
!          phih(i) = (1.-aphi16*zol1)**(-1./2.)
           tem     = 1.0 / (1. - aphi16*zol1)
           phih(i) = sqrt(tem)
           phim(i) = sqrt(phih(i))
         else
           phim(i) = 1. + aphi5*zol1
           phih(i) = phim(i)
         endif
         wscale(i) = ustar(i)/phim(i)
         ustmin(i) = ustar(i)/aphi5
         wscale(i) = max(wscale(i),ustmin(i))
      enddo
      do i=1,im
        if(pblflg(i)) then
          if(zol(i) < zolcru .and. kpbl(i) > 1) then
            pcnvflg(i) = .true.
          else
            ublflg(i) = .true.
          endif
          wst3 = govrth(i)*sflux(i)*hpbl(i)
          wstar(i)= wst3**h1
          ust3 = ustar(i)**3.
          wscaleu(i) = (ust3+wfac*vk*wst3*sfcfrac)**h1
          wscaleu(i) = max(wscaleu(i),ustmin(i))
        endif
      enddo
!
! compute counter-gradient mixing term for heat and moisture
!>  ## Update thermal properties of surface parcel and recompute PBL height ("Corrector step").
!!  Next, the counter-gradient terms for temperature and humidity are calculated using equation 4 of Hong and Pan (1996) \cite hong_and_pan_1996 and are used to calculate the "scaled virtual temperature excess near the surface" (equation 9 in Hong and Pan (1996) \cite hong_and_pan_1996) so that the properties of the thermal are updated to recalculate the PBL height.
      do i = 1,im
         if(ublflg(i)) then
           hgamt(i)  = min(cfac*heat(i)/wscaleu(i),gamcrt)
           hgamq(i)  = min(cfac*evap(i)/wscaleu(i),gamcrq)
           vpert     = hgamt(i) + hgamq(i)*fv*theta(i,1)
           vpert     = min(vpert,gamcrt)
           thermal(i)= thermal(i)+max(vpert,0.)
           hgamt(i)  = max(hgamt(i),0.0)
           hgamq(i)  = max(hgamq(i),0.0)
         endif
      enddo
!
!  enhance the pbl height by considering the thermal excess
!>  The PBL height calculation follows the same procedure as the predictor step, except that it uses an updated virtual potential temperature for the thermal.
      do i=1,im
         flg(i)  = .true.
         if(ublflg(i)) then
           flg(i)  = .false.
           rbup(i) = rbsoil(i)
         endif
      enddo
      do k = 2, kmpbl
      do i = 1, im
        if(.not.flg(i)) then
          rbdn(i) = rbup(i)
          spdk2   = max((u1(i,k)**2+v1(i,k)**2),1.)
          rbup(i) = (thvx(i,k)-thermal(i))*
     &              (grav*zl(i,k)/thvx(i,1))/spdk2
          kpbl(i) = k
          flg(i)  = rbup(i) > crb(i)
        endif
      enddo
      enddo
      do i = 1,im
        if(ublflg(i)) then
           k = kpbl(i)
           if(rbdn(i) >= crb(i)) then
              rbint = 0.
           elseif(rbup(i) <= crb(i)) then
              rbint = 1.
           else
              rbint = (crb(i)-rbdn(i))/(rbup(i)-rbdn(i))
           endif
           hpbl(i) = zl(i,k-1) + rbint*(zl(i,k)-zl(i,k-1))
           if(hpbl(i) < zi(i,kpbl(i))) kpbl(i) = kpbl(i) - 1
           if(kpbl(i) <= 1) then
              ublflg(i) = .false.
              pblflg(i) = .false.
           endif
        endif
      enddo
!
!  look for stratocumulus
!>  ## Determine whether stratocumulus layers exist and compute quantities needed for enhanced diffusion
!!  - Starting at the PBL top and going downward, if the level is less than 2.5 km and \f$q_l>q_{l,cr}\f$ then set kcld = k (find the cloud top index in the PBL). If no cloud water above the threshold is found, scuflg is set to F.
      do i = 1, im
        flg(i)=scuflg(i)
      enddo
      do k = kmpbl,1,-1
      do i = 1, im
        if(flg(i) .and. k <= lcld(i)) then
          if(qlx(i,k).ge.qlcr) then
             kcld(i)=k
             flg(i)=.false.
          endif
        endif
      enddo
      enddo
      do i = 1, im
        if(scuflg(i) .and. kcld(i)==km1) scuflg(i)=.false.
      enddo
!>  - Starting at the PBL top and going downward, if the level is less than the cloud top, find the level of the minimum radiative heating rate within the cloud. If the level of the minimum is the lowest model level or the minimum radiative heating rate is positive, then set scuflg to F.
      do i = 1, im
        flg(i)=scuflg(i)
      enddo
      do k = kmpbl,1,-1
      do i = 1, im
        if(flg(i) .and. k <= kcld(i)) then
          if(qlx(i,k) >= qlcr) then
            if(radx(i,k) < radmin(i)) then
              radmin(i)=radx(i,k)
              krad(i)=k
            endif
          else
            flg(i)=.false.
          endif
        endif
      enddo
      enddo
      do i = 1, im
        if(scuflg(i) .and. krad(i) <= 1) scuflg(i)=.false.
        if(scuflg(i) .and. radmin(i)>=0.) scuflg(i)=.false.
      enddo
!>  - Starting at the PBL top and going downward, count the number of levels below the minimum radiative heating rate level that have cloud water above the threshold. If there are none, then set the scuflg to F.
      do i = 1, im
        flg(i)=scuflg(i)
      enddo
      do k = kmpbl,2,-1
      do i = 1, im
        if(flg(i) .and. k <= krad(i)) then
          if(qlx(i,k) >= qlcr) then
            icld(i)=icld(i)+1
          else
            flg(i)=.false.
          endif
        endif
      enddo
      enddo
      do i = 1, im
        if(scuflg(i) .and. icld(i) < 1) scuflg(i)=.false.
      enddo
!>  - Find the height of the interface where the minimum in radiative heating rate is located. If this height is less than the second model interface height, then set the scuflg to F.
      do i = 1, im
        if(scuflg(i)) then
           hrad(i) = zi(i,krad(i)+1)
!          hradm(i)= zl(i,krad(i))
        endif
      enddo
!
      do i = 1, im
        if(scuflg(i) .and. hrad(i)<zi(i,2)) scuflg(i)=.false.
      enddo
!>  - Calculate the hypothetical \f$\theta_v\f$ at the minimum radiative heating level that a parcel would reach due to radiative cooling after a typical cloud turnover time spent at that level.
      do i = 1, im
        if(scuflg(i)) then
          k    = krad(i)
          tem  = zi(i,k+1)-zi(i,k)
          tem1 = cldtime*radmin(i)/tem
          thlvx1(i) = thlvx(i,k)+tem1
!         if(thlvx1(i) > thlvx(i,k-1)) scuflg(i)=.false.
        endif
      enddo
!>  - Determine the distance that a parcel would sink downwards starting from the level of minimum radiative heating rate by comparing the hypothetical minimum \f$\theta_v\f$ calculated above with the environmental \f$\theta_v\f$.
      do i = 1, im
         flg(i)=scuflg(i)
      enddo
      do k = kmpbl,1,-1
      do i = 1, im
        if(flg(i) .and. k <= krad(i))then
          if(thlvx1(i) <= thlvx(i,k))then
             tem=zi(i,k+1)-zi(i,k)
             zd(i)=zd(i)+tem
          else
             flg(i)=.false.
          endif
        endif
      enddo
      enddo
!>  - Calculate the cloud thickness, where the cloud top is the in-cloud minimum radiative heating level and the bottom is determined previously.
      do i = 1, im
        if(scuflg(i))then
          kk = max(1, krad(i)+1-icld(i))
          zdd(i) = hrad(i)-zi(i,kk)
        endif
      enddo
!>  - Find the largest between the cloud thickness and the distance of a sinking parcel, then determine the smallest of that number and the height of the minimum in radiative heating rate. Set this number to \f$zd\f$. Using \f$zd\f$, calculate the characteristic velocity scale of cloud-top radiative cooling-driven turbulence.
      do i = 1, im
        if(scuflg(i))then
          zd(i) = max(zd(i),zdd(i))
          zd(i) = min(zd(i),hrad(i))
          tem   = govrth(i)*zd(i)*(-radmin(i))
          vrad(i)= tem**h1
        endif
      enddo
!
!     compute inverse prandtl number
!>  ## Calculate the inverse Prandtl number
!!  For an unstable PBL, the Prandtl number is calculated according to Hong and Pan (1996) \cite hong_and_pan_1996, equation 10, whereas for a stable boundary layer, the Prandtl number is simply \f$Pr = \frac{\phi_h}{\phi_m}\f$.
      do i = 1, im
        if(ublflg(i)) then
          tem = phih(i)/phim(i)+cfac*vk*sfcfrac
        else
          tem = phih(i)/phim(i)
        endif
        prinv(i) =  1.0 / tem
        prinv(i) = min(prinv(i),prmax)
        prinv(i) = max(prinv(i),prmin)
      enddo
      do i = 1, im
        if(zol(i) > zolcr) then
          kpbl(i) = 1
        endif
      enddo


!!! 20150915 WeiguoWang added alpha (moninq_fac) and wind-dependent modification of K by RGF
! -------------------------------------------------------------------------------------
! begin RGF modifications
! this is version MOD05

! RGF determine wspd at roughly 500 m above surface, or as close as possible,
! reuse SPDK2
!  zi(i,k) is AGL, right?  May not matter if applied only to water grid points
      if(hurr_pbl .and. moninq_fac < 0.0) then
        do i=1,im
          spdk2 = 0.
          wspm(i,1) = 0.
          do k = 1, kmpbl ! kmpbl is like a max possible pbl height
            if (zi(i,k) .le. 500. .and. zi(i,k+1) .gt. 500.) then ! find level bracketing 500 m
              spdk2 = SQRT(u1(i,k)*u1(i,k)+v1(i,k)*v1(i,k)) ! wspd near 500 m
              wspm(i,1) = spdk2/0.6  ! now the Km limit for 500 m.  just store in K=1
              wspm(i,2) = float(k)  ! height of level at gridpoint i. store in K=2
            endif
          enddo !k
        enddo ! i
      endif ! hurr_pbl and moninq_fac < 0


!     compute diffusion coefficients below pbl
!>  ## Compute diffusion coefficients below the PBL top
!!  Below the PBL top, the diffusion coefficients (\f$K_m\f$ and \f$K_h\f$) are calculated according to equation 2 in Hong and Pan (1996) \cite hong_and_pan_1996 where a different value for \f$w_s\f$ (PBL vertical velocity scale) is used depending on the PBL stability. \f$K_h\f$ is calculated from \f$K_m\f$ using the Prandtl number. The calculated diffusion coefficients are checked so that they are bounded by maximum values and the local background diffusion coefficients.
      if (.not. (hurr_pbl .and. moninq_fac < 0.0)) then
        do k = 1, kmpbl
          do i=1,im
            if(k < kpbl(i)) then
!                   zfac = max((1.-(zi(i,k+1)-zl(i,1))/
!    1               (hpbl(i)-zl(i,1))), zfmin)
              zfac = max((1.-zi(i,k+1)/hpbl(i)), zfmin)
              tem = zi(i,k+1) * (zfac**pfac) * moninq_fac ! lmh suggested by kg
              if(pblflg(i)) then
                tem1 = vk * wscaleu(i) * tem
!                       dku(i,k) = xkzmo(i,k) + tem1
!               dkt(i,k) = xkzo(i,k)  + tem1 * prinv(i)
                dku(i,k) = tem1
                dkt(i,k) = tem1 * prinv(i)
              else
                tem1 = vk * wscale(i) * tem
!                       dku(i,k) = xkzmo(i,k) + tem1
!               dkt(i,k) = xkzo(i,k)  + tem1 * prinv(i)
                dku(i,k) = tem1
                dkt(i,k) = tem1 * prinv(i)
              endif
              dku(i,k) = min(dku(i,k),dkmax)
              dku(i,k) = max(dku(i,k),xkzmo(i,k))
              dkt(i,k) = min(dkt(i,k),dkmax)
              dkt(i,k) = max(dkt(i,k),xkzo(i,k))
              dktx(i,k)= dkt(i,k)
            endif
          enddo !i
        enddo !k
      else
        !hurricane PBL case and moninq_fac < 0 (note that the i and k loop order has been switched)
        do i=1, im
          do k=1, kmpbl
            if (k < kpbl(i)) then
!             zfac = max((1.-(zi(i,k+1)-zl(i,1))/
!    1             (hpbl(i)-zl(i,1))), zfmin)
              zfac = max((1.-zi(i,k+1)/hpbl(i)), zfmin)
              tem = zi(i,k+1) * (zfac**pfac) * ABS(moninq_fac)

!!!! CHANGES FOR HEIGHT-DEPENDENT K ADJUSTMENT, WANG W
              if (useshape .ge. 1) then
                sz2h=(zi(i,k+1)-zl(i,1))/(hpbl(i)-zl(i,1))
                sz2h=max(sz2h,zfmin)
                sz2h=min(sz2h,1.0)
                zfac=(1.0-sz2h)**pfac
!                smax=0.148  !! max value of this shape function
                smax=0.148  !! max value of this shape function
                hmax=0.333  !! roughly height if max K
                skmax=hmax*(1.0-hmax)**pfac
                sksfc=min(zi(i,2)/hpbl(i),0.05)  ! surface layer top, 0.05H or ZI(2) (Zi(1)=0)
                sksfc=sksfc*(1-sksfc)**pfac

                zfac=max(zfac,zfmin)
                ashape=max(ABS(moninq_fac),0.2)  ! should not be smaller than 0.2, otherwise too much adjustment(?)
                if (useshape == 1) then 
                 ashape=(1.0 - ((sz2h*zfac/smax)**0.25) *(1.0 - ashape))
                 tem = zi(i,k+1) * (zfac) * ashape
                elseif (useshape == 2) then   !only adjus K that is > K_surface_top
                  ashape1=1.0
                  if (skmax > sksfc) then 
                    ashape1=(skmax*ashape-sksfc)/(skmax-sksfc)
                  endif
                  skminusk0 = zi(i,k+1)*zfac - hpbl(i)*sksfc
                  tem = zi(i,k+1) * (zfac) ! no adjustment
                  if (skminusk0 > 0) then   ! only adjust K which is > surface top K
                    tem = skminusk0*ashape1 + hpbl(i)*sksfc
                  endif
                endif ! useshape == 1 or 2
              endif  ! endif useshape>1
!!!! END OF CHANGES , WANG W

!!If alpha >= 0, this is the only modification of K
! if alpha = -1, the above provides the first guess for DKU, based on assumption
! alpha = +1
!               (other values of alpha < 0 can also be applied)
! if alpha > 0, the above applies the alpha suppression factor and we are
! finished

              if(pblflg(i)) then
                tem1 = vk * wscaleu(i) * tem
!                 dku(i,k) = xkzmo(i,k) + tem1
!               dkt(i,k) = xkzo(i,k)  + tem1 * prinv(i)
                dku(i,k) = tem1
                dkt(i,k) = tem1 * prinv(i)
              else
                tem1 = vk * wscale(i) * tem
!                 dku(i,k) = xkzmo(i,k) + tem1
!               dkt(i,k) = xkzo(i,k)  + tem1 * prinv(i)
                dku(i,k) = tem1
                dkt(i,k) = tem1 * prinv(i)
              endif
              dku(i,k) = min(dku(i,k),dkmax)
              dku(i,k) = max(dku(i,k),xkzmo(i,k))
              dkt(i,k) = min(dkt(i,k),dkmax)
              dkt(i,k) = max(dkt(i,k),xkzo(i,k))
              dktx(i,k)= dkt(i,k)
            endif !k < kpbl(i)
          enddo     !K loop

! possible modification of first guess DKU, under certain conditions
! (1) this applies only to columns over water
          if (islimsk(i) .eq. 0) then ! sea only
! (2) alpha test
! if alpha < 0, find alpha for each column and do the loop again
! if alpha > 0, we are finished
      
!GJF: redundant check for moninq_fac < 0?
            if (moninq_fac .lt. 0.) then      ! variable alpha test
! k-level of layer around 500 m
              kLOC = INT(wspm(i,2))
!            print *,' kLOC ',kLOC,' KPBL ',KPBL(I)

! (3) only do  this IF KPBL(I) >= kLOC.  Otherwise, we are finished, with DKU as
! if alpha = +1
              if(kpbl(i) .gt. kLOC) then
                xDKU = DKU(i,kLOC)     ! Km at k-level
! (4) DKU check.
! WSPM(i,1) is the KM cap for the 500-m level.
!  if DKU at 500-m level < WSPM(i,1), do not limit Km ANYWHERE.  Alpha =
!  abs(alpha).  No need to recalc.
!  if DKU at 500-m level > WSPM(i,1), then alpha = WSPM(i,1)/xDKU for entire
!  column
                if(xDKU .ge. wspm(i,1)) then ! ONLY if DKU at 500-m exceeds cap, otherwise already done
                  wspm(i,3) = wspm(i,1)/xDKU  ! ratio of cap to Km at k-level, store in WSPM(i,3)
                  !WSPM(i,4) = amin1(WSPM(I,3),1.0) ! this is new column alpha. cap at 1. ! should never be needed
                  wspm(i,4) = min(wspm(i,3),1.0) ! this is new column alpha. cap at 1. ! should never be needed
 !! recalculate K capped by WSPM(i,1)           
                  do k = 1, kmpbl
                    if(k < kpbl(i)) then
!                     zfac = max((1.-(zi(i,k+1)-zl(i,1))/
!    1                       (hpbl(i)-zl(i,1))), zfmin)
                      zfac = max((1.-zi(i,k+1)/hpbl(i)), zfmin)
                      tem = zi(i,k+1) * (zfac**pfac) * wspm(i,4)
!!! wang use different K shape, options!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! HANGES FOR HEIGHT-DEPENDENT K ADJUSTMENT, WANG W
                      if(useshape .ge. 1) then
                        sz2h=(zi(i,k+1)-zl(i,1))/(hpbl(i)-zl(i,1))
                        sz2h=max(sz2h,zfmin)
                        sz2h=min(sz2h,1.0)
                        zfac=(1.0-sz2h)**pfac
                        smax=0.148  !! max value of this shape function
                        hmax=0.333  !! roughly height if max K
                        skmax=hmax*(1.0-hmax)**pfac
                        sksfc=min(zi(i,2)/hpbl(i),0.05)  ! surface layer top, 0.05H or ZI(2) (Zi(1)=0)
                        sksfc=sksfc*(1-sksfc)**pfac
                                
                        zfac=max(zfac,zfmin)
                        ashape=max(wspm(i,4),0.2)  !! adjustment coef should not smaller than 0.2
                        if(useshape ==1) then 
                          ashape=(1.0 - ((sz2h*zfac/smax)**0.25)*
     &                           (1.0 - ashape))
                          tem = zi(i,k+1) * (zfac) * ashape
                        elseif (useshape == 2) then !only adjus K that is > K_surface_top 
                          ashape1=1.0
                          if (skmax > sksfc) then 
                            ashape1=(skmax*ashape-sksfc)/(skmax-sksfc)
                          endif
                          skminusk0=zi(i,k+1)*zfac - hpbl(i)*sksfc
                          tem = zi(i,k+1) * (zfac) ! no adjustment
                          if (skminusk0 > 0) then   ! only adjust K which is > surface top K
                            tem = skminusk0*ashape1 + HPBL(i)*sksfc
                          endif
                        endif  ! endif useshape=1 or 2
                      endif  ! endif useshape>1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                      if(pblflg(i)) then
                        tem1 = vk * wscaleu(i) * tem
!                        dku(i,k) = xkzmo(i,k) + tem1
!                        dkt(i,k) = xkzo(i,k)  + tem1 * prinv(i)
                        dku(i,k) = tem1
                        dkt(i,k) = tem1 * prinv(i)
                      else
                        tem1 = vk * wscale(i) * tem
!                        dku(i,k) = xkzmo(i,k) + tem1
!                        dkt(i,k) = xkzo(i,k)  + tem1 * prinv(i)
                        dku(i,k) = tem1
                        dkt(i,k) = tem1 * prinv(i)
                      endif !pblflg
                      dku(i,k) = min(dku(i,k),dkmax)
                      dku(i,k) = max(dku(i,k),xkzmo(i,k))
                      dkt(i,k) = min(dkt(i,k),dkmax)
                      dkt(i,k) = max(dkt(i,k),xkzo(i,k))
                      dktx(i,k)= dkt(i,k)
                    endif ! k < kpbl(i)
                  enddo ! K loop
                endif ! xDKU .ge. wspm(i,1)
              endif ! kpbl(i) .ge. kLOC
            endif ! moninq_fac < 0 (GJF: redundant?)
          endif ! islimsk == 0
        enddo ! I loop
      endif ! not (hurr_pbl and moninq_fac < 0)
!
! compute diffusion coefficients based on local scheme above pbl
!>  ## Compute diffusion coefficients above the PBL top
!!  Diffusion coefficients above the PBL top are computed as a function of local stability (gradient Richardson number), shear, and a length scale from Louis (1979) \cite louis_1979 :
!!  \f[
!!  K_{m,h}=l^2f_{m,h}(Ri_g)\left|\frac{\partial U}{\partial z}\right|
!!  \f]
!!  The functions used (\f$f_{m,h}\f$) depend on the local stability. First, the gradient Richardson number is calculated as
!!  \f[
!!  Ri_g=\frac{\frac{g}{T}\frac{\partial \theta_v}{\partial z}}{\frac{\partial U}{\partial z}^2}
!!  \f]
!!  where \f$U\f$ is the horizontal wind. For the unstable case (\f$Ri_g < 0\f$), the Richardson number-dependent functions are given by
!!  \f[
!!  f_h(Ri_g) = 1 + \frac{8\left|Ri_g\right|}{1 + 1.286\sqrt{\left|Ri_g\right|}}\\
!!  \f]
!!  \f[
!!  f_m(Ri_g) = 1 + \frac{8\left|Ri_g\right|}{1 + 1.746\sqrt{\left|Ri_g\right|}}\\
!!  \f]
!!  For the stable case, the following formulas are used
!!  \f[
!!  f_h(Ri_g) = \frac{1}{\left(1 + 5Ri_g\right)^2}\\
!!  \f]
!!  \f[
!!  Pr = \frac{K_h}{K_m} = 1 + 2.1Ri_g
!!  \f]
!!  The source for the formulas used for the Richardson number-dependent functions is unclear. They are different than those used in Hong and Pan (1996) \cite hong_and_pan_1996 as the previous documentation suggests. They follow equation 14 of Louis (1979) \cite louis_1979 for the unstable case, but it is unclear where the values of the coefficients \f$b\f$ and \f$c\f$ from that equation used in this scheme originate. Finally, the length scale, \f$l\f$ is calculated according to the following formula from Hong and Pan (1996) \cite hong_and_pan_1996
!!  \f[
!!  \frac{1}{l} = \frac{1}{kz} + \frac{1}{l_0}\\
!!  \f]
!!  \f[
!!  or\\
!!  \f]
!!  \f[
!!  l=\frac{l_0kz}{l_0+kz}
!!  \f]
!!  where \f$l_0\f$ is currently 30 m for stable conditions and 150 m for unstable. Finally, the diffusion coefficients are kept in a range bounded by the background diffusion and the maximum allowable values.
      do k = 1, km1
         do i=1,im
            if(k >= kpbl(i)) then
               bvf2 = grav*bf(i,k)*ti(i,k)
               ri   = max(bvf2/shr2(i,k),rimin)
               zk   = vk*zi(i,k+1)
               if(ri < 0.) then ! unstable regime
                  rl2      = zk*rlamun/(rlamun+zk)
                  dk       = rl2*rl2*sqrt(shr2(i,k))
                  sri      = sqrt(-ri)
!                 dku(i,k) = xkzmo(i,k) + dk*(1+8.*(-ri)/(1+1.746*sri))
!                 dkt(i,k) = xkzo(i,k)  + dk*(1+8.*(-ri)/(1+1.286*sri))
                  dku(i,k) = dk*(1+8.*(-ri)/(1+1.746*sri))
                  dkt(i,k) = dk*(1+8.*(-ri)/(1+1.286*sri))
               else             ! stable regime
                  rl2      = zk*rlam/(rlam+zk)
!!                tem      = rlam * sqrt(0.01*prsi(i,k))
!!                rl2      = zk*tem/(tem+zk)
                  dk       = rl2*rl2*sqrt(shr2(i,k))
                  tem1     = dk/(1+5.*ri)**2
!
                  if(k >= kpblx(i)) then
                    prnum = 1.0 + 2.1*ri
                    prnum = min(prnum,prmax)
                  else
                    prnum = 1.0
                  endif
!                 dku(i,k) = xkzmo(i,k) + tem1 * prnum
!                 dkt(i,k) = xkzo(i,k)  + tem1
                  dku(i,k) = tem1 * prnum
                  dkt(i,k) = tem1
               endif
!
               dku(i,k) = min(dku(i,k),dkmax)
               dku(i,k) = max(dku(i,k),xkzmo(i,k))
               dkt(i,k) = min(dkt(i,k),dkmax)
               dkt(i,k) = max(dkt(i,k),xkzo(i,k))
!
            endif
!
         enddo
      enddo
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  compute components for mass flux mixing by large thermals
!>  ## If the PBL is convective, call the mass flux scheme to replace the countergradient terms.
!!  If the PBL is convective, the updraft properties are initialized to be the same as the state variables and the subroutine mfpbl is called.
      do k = 1, km
        do i = 1, im
          if(pcnvflg(i)) then
            tcko(i,k) = t1(i,k)
            ucko(i,k) = u1(i,k)
            vcko(i,k) = v1(i,k)
            xmf(i,k) = 0.
          endif
        enddo
      enddo
      do kk = 1, ntrac
      do k = 1, km
        do i = 1, im
          if(pcnvflg(i)) then
            qcko(i,k,kk) = q1(i,k,kk)
          endif
        enddo
      enddo
      enddo
!>  For details of the mfpbl subroutine, step into its documentation ::mfpbl
      call mfpbl(im,im,km,ntrac,dt2,pcnvflg,
     &       zl,zi,thvx,q1,t1,u1,v1,hpbl,kpbl,
     &       sflux,ustar,wstar,xmf,tcko,qcko,ucko,vcko)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  compute diffusion coefficients for cloud-top driven diffusion
!  if the condition for cloud-top instability is met,
!    increase entrainment flux at cloud top
!
!>  ## Compute enhanced diffusion coefficients related to stratocumulus-topped PBLs
!!  If a stratocumulus layer has been identified in the PBL, the diffusion coefficients in the PBL are modified in the following way.
!!
!!  -# First, the criteria for CTEI is checked, using the threshold from equation 13 of Macvean and Mason (1990) \cite macvean_and_mason_1990. If the criteria is met, the cloud top diffusion is increased:
!!  \f[
!!  K_h^{Sc} = -c\frac{\Delta F_R}{\rho c_p}\frac{1}{\frac{\partial \theta_v}{\partial z}}
!!  \f]
!!  where the constant \f$c\f$ is set to 0.2 if the CTEI criterion is not met and 1.0 if it is.
!!
!!  -# Calculate the diffusion coefficients due to stratocumulus mixing according to equation 5 in Lock et al. (2000) \cite lock_et_al_2000 for every level below the stratocumulus top using the characteristic stratocumulus velocity scale previously calculated. The diffusion coefficient for momentum is calculated assuming a constant inverse Prandtl number of 0.75.
      do i = 1, im
        if(scuflg(i)) then
           k = krad(i)
           tem = thetae(i,k) - thetae(i,k+1)
           tem1 = qtx(i,k) - qtx(i,k+1)
           if (tem > 0. .and. tem1 > 0.) then
             cteit= cp*tem/(hvap*tem1)
             if(cteit > actei) rent(i) = rentf2
           endif
        endif
      enddo
      do i = 1, im
        if(scuflg(i)) then
           k = krad(i)
           tem1  = max(bf(i,k),tdzmin)
           ckt(i,k) = -rent(i)*radmin(i)/tem1
           cku(i,k) = ckt(i,k)
        endif
      enddo
!
      do k = 1, kmpbl
         do i=1,im
            if(scuflg(i) .and. k < krad(i)) then
               tem1=hrad(i)-zd(i)
               tem2=zi(i,k+1)-tem1
               if(tem2 > 0.) then
                  ptem= tem2/zd(i)
                  if(ptem.ge.1.) ptem= 1.
                  ptem= tem2*ptem*sqrt(1.-ptem)
                  ckt(i,k) = radfac*vk*vrad(i)*ptem
                  cku(i,k) = 0.75*ckt(i,k)
                  ckt(i,k) = max(ckt(i,k),dkmin)
                  ckt(i,k) = min(ckt(i,k),dkmax)
                  cku(i,k) = max(cku(i,k),dkmin)
                  cku(i,k) = min(cku(i,k),dkmax)
               endif
            endif
         enddo
      enddo
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!>  After \f$K_h^{Sc}\f$ has been determined from the surface to the top of the stratocumulus layer, it is added to the value for the diffusion coefficient calculated previously using surface-based mixing [see equation 6 of Lock et al. (2000) \cite lock_et_al_2000 ].
      if (.not. hurr_pbl) then
        do k = 1, kmpbl
          do i=1,im
            if(scuflg(i)) then
               dkt(i,k) = dkt(i,k)+ckt(i,k)
               dku(i,k) = dku(i,k)+cku(i,k)
               dkt(i,k) = min(dkt(i,k),dkmax)
               dku(i,k) = min(dku(i,k),dkmax)
            endif
          enddo
        enddo
      else
        do k = 1, kmpbl
          do i=1,im
            if(scuflg(i)) then
               !! if K needs to be adjusted by alpha, then no need to add this term
               if (.not. (hurr_pbl .and. moninq_fac < 0.0)) then
                 dkt(i,k) = dkt(i,k)+ckt(i,k)
                 dku(i,k) = dku(i,k)+cku(i,k)
               end if
               dkt(i,k) = min(dkt(i,k),dkmax)
               dku(i,k) = min(dku(i,k),dkmax)
            endif
          enddo
        enddo
      endif
!
!     compute tridiagonal matrix elements for heat and moisture
!
!>  ## Solve for the temperature and moisture tendencies due to vertical mixing.
!!  The tendencies of heat, moisture, and momentum due to vertical diffusion are calculated using a two-part process. First, a solution is obtained using an implicit time-stepping scheme, then the time tendency terms are "backed out". The tridiagonal matrix elements for the implicit solution for temperature and moisture are prepared in this section, with differing algorithms depending on whether the PBL was convective (substituting the mass flux term for counter-gradient term), unstable but not convective (using the computed counter-gradient terms), or stable (no counter-gradient terms).
      do i=1,im
         ad(i,1) = 1.
         a1(i,1) = t1(i,1)   + beta(i) * heat(i)
         a2(i,1) = q1(i,1,1) + beta(i) * evap(i)
      enddo

      if(ntrac >= 2) then
        do k = 2, ntrac
          is = (k-1) * km
          do i = 1, im
            a2(i,1+is) = q1(i,1,k)
          enddo
        enddo
      endif
!
      do k = 1,km1
        do i = 1,im
          dtodsd = dt2/del(i,k)
          dtodsu = dt2/del(i,k+1)
          dsig   = prsl(i,k)-prsl(i,k+1)
          rdz    = rdzt(i,k)
          tem1   = dsig * dkt(i,k) * rdz
          dsdz2     = tem1 * rdz
          au(i,k)   = -dtodsd*dsdz2
          al(i,k)   = -dtodsu*dsdz2
!
          if(pcnvflg(i) .and. k < kpbl(i)) then
             tem2      = dsig * rdz
             ptem      = 0.5 * tem2 * xmf(i,k)
             ptem1     = dtodsd * ptem
             ptem2     = dtodsu * ptem
             ad(i,k)   = ad(i,k)-au(i,k)-ptem1
             ad(i,k+1) = 1.-al(i,k)+ptem2
             au(i,k)   = au(i,k)-ptem1
             al(i,k)   = al(i,k)+ptem2
             ptem      = tcko(i,k) + tcko(i,k+1)
             dsdzt     = tem1 * gocp
             a1(i,k)   = a1(i,k)+dtodsd*dsdzt-ptem1*ptem
             a1(i,k+1) = t1(i,k+1)-dtodsu*dsdzt+ptem2*ptem
             ptem      = qcko(i,k,1) + qcko(i,k+1,1)
             a2(i,k)   = a2(i,k) - ptem1 * ptem
             a2(i,k+1) = q1(i,k+1,1) + ptem2 * ptem
          elseif(ublflg(i) .and. k < kpbl(i)) then
             ptem1 = dsig * dktx(i,k) * rdz
             tem   = 1.0 / hpbl(i)
             dsdzt = tem1 * gocp - ptem1 * hgamt(i) * tem
             dsdzq = - ptem1 * hgamq(i) * tem
             ad(i,k)   = ad(i,k)-au(i,k)
             ad(i,k+1) = 1.-al(i,k)
             a1(i,k)   = a1(i,k)+dtodsd*dsdzt
             a1(i,k+1) = t1(i,k+1)-dtodsu*dsdzt
             a2(i,k)   = a2(i,k)+dtodsd*dsdzq
             a2(i,k+1) = q1(i,k+1,1)-dtodsu*dsdzq
          else
             ad(i,k)   = ad(i,k)-au(i,k)
             ad(i,k+1) = 1.-al(i,k)
             dsdzt     = tem1 * gocp
             a1(i,k)   = a1(i,k)+dtodsd*dsdzt
             a1(i,k+1) = t1(i,k+1)-dtodsu*dsdzt
             a2(i,k+1) = q1(i,k+1,1)
          endif
!
        enddo
      enddo
!
      if(ntrac >= 2) then
        do kk = 2, ntrac
          is = (kk-1) * km
          do k = 1, km1
            do i = 1, im
              if(pcnvflg(i) .and. k < kpbl(i)) then
                dtodsd = dt2/del(i,k)
                dtodsu = dt2/del(i,k+1)
                dsig  = prsl(i,k)-prsl(i,k+1)
                tem   = dsig * rdzt(i,k)
                ptem  = 0.5 * tem * xmf(i,k)
                ptem1 = dtodsd * ptem
                ptem2 = dtodsu * ptem
                tem1  = qcko(i,k,kk) + qcko(i,k+1,kk)
                a2(i,k+is) = a2(i,k+is) - ptem1*tem1
                a2(i,k+1+is)= q1(i,k+1,kk) + ptem2*tem1
              else
                a2(i,k+1+is) = q1(i,k+1,kk)
              endif
            enddo
          enddo
        enddo
      endif
!
!     solve tridiagonal problem for heat and moisture
!
!>  The tridiagonal system is solved by calling the internal ::tridin subroutine.
      call tridin(im,km,ntrac,al,ad,au,a1,a2,au,a1,a2)

!
!     recover tendencies of heat and moisture
!
!>  After returning with the solution, the tendencies for temperature and moisture are recovered.
      do  k = 1,km
         do i = 1,im
            ttend      = (a1(i,k)-t1(i,k)) * rdt
            qtend      = (a2(i,k)-q1(i,k,1))*rdt
            tau(i,k)   = tau(i,k)+ttend
            rtg(i,k,1) = rtg(i,k,1)+qtend
            dtsfc(i)   = dtsfc(i)+cont*del(i,k)*ttend
            dqsfc(i)   = dqsfc(i)+conq*del(i,k)*qtend
            if(lssav .and. ldiag3d .and. .not.                          &
     &                flag_for_pbl_generic_tend) then
               dt3dt_PBL(i,k) = dt3dt_PBL(i,k) + ttend*delt
               if(qdiag3d) then
                  dq3dt_PBL(i,k) = dq3dt_PBL(i,k) + qtend*delt
               endif
            endif
         enddo
      enddo
      if(ntrac >= 2) then
        do kk = 2, ntrac
          is = (kk-1) * km
          do k = 1, km
            do i = 1, im
              qtend = (a2(i,k+is)-q1(i,k,kk))*rdt
              rtg(i,k,kk) = rtg(i,k,kk)+qtend
            enddo
          enddo
        enddo
        if(lssav .and. ldiag3d .and. ntoz>0 .and. qdiag3d .and.         &
     &               .not. flag_for_pbl_generic_tend) then
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
!   compute tke dissipation rate
!
!>  ## Calculate heating due to TKE dissipation and add to the tendency for temperature
!!  Following Han et al. (2016) \cite Han_2016 , turbulence dissipation contributes to the tendency of temperature in the following way. First, turbulence dissipation is calculated by equation 17 of Han et al. (2016) \cite Han_2016 for the PBL and equation 16 for the surface layer.
      if(dspheat) then
!
      do k = 1,km1
        do i = 1,im
          diss(i,k) = dku(i,k)*shr2(i,k)-grav*ti(i,k)*dkt(i,k)*bf(i,k)
!         diss(i,k) = dku(i,k)*shr2(i,k)
        enddo
      enddo
!
!     add dissipative heating at the first model layer
!
!>  Next, the temperature tendency is updated following equation 14.
      if (hurr_pbl .and. moninq_fac < 0.0) then
        ttend_fac = 0.7
      else
        ttend_fac = 0.5
      endif
      
      do i = 1,im
         tem   = govrth(i)*sflux(i)
         tem1  = tem + stress(i)*spd1(i)/zl(i,1)
         tem2  = 0.5 * (tem1+diss(i,1))
         tem2  = max(tem2, 0.)
         ttend = tem2 / cp
         tau(i,1) = tau(i,1)+ttend_fac*ttend
      enddo
!
!     add dissipative heating above the first model layer
!
      do k = 2,km1
        do i = 1,im
          tem = 0.5 * (diss(i,k-1)+diss(i,k))
          tem  = max(tem, 0.)
          ttend = tem / cp
          tau(i,k) = tau(i,k) + ttend_fac*ttend
        enddo
      enddo
!
      endif
!
!     compute tridiagonal matrix elements for momentum
!
!>  ## Solve for the horizontal momentum tendencies and add them to the output tendency terms
!!  As with the temperature and moisture tendencies, the horizontal momentum tendencies are calculated by solving tridiagonal matrices after the matrices are prepared in this section.
      do i=1,im
         ad(i,1) = 1.0 + beta(i) * stress(i) / spd1(i)
         a1(i,1) = u1(i,1)
         a2(i,1) = v1(i,1)
      enddo
!
      do k = 1,km1
        do i=1,im
          dtodsd  = dt2/del(i,k)
          dtodsu  = dt2/del(i,k+1)
          dsig    = prsl(i,k)-prsl(i,k+1)
          rdz     = rdzt(i,k)
          tem1    = dsig*dku(i,k)*rdz
          dsdz2   = tem1 * rdz
          au(i,k) = -dtodsd*dsdz2
          al(i,k) = -dtodsu*dsdz2
!
          if(pcnvflg(i) .and. k < kpbl(i)) then
             tem2      = dsig * rdz
             ptem      = 0.5 * tem2 * xmf(i,k)
             ptem1     = dtodsd * ptem
             ptem2     = dtodsu * ptem
             ad(i,k)   = ad(i,k)-au(i,k)-ptem1
             ad(i,k+1) = 1.-al(i,k)+ptem2
             au(i,k)   = au(i,k)-ptem1
             al(i,k)   = al(i,k)+ptem2
             ptem      = ucko(i,k) + ucko(i,k+1)
             a1(i,k)   = a1(i,k) - ptem1 * ptem
             a1(i,k+1) = u1(i,k+1) + ptem2 * ptem
             ptem      = vcko(i,k) + vcko(i,k+1)
             a2(i,k)   = a2(i,k) - ptem1 * ptem
             a2(i,k+1) = v1(i,k+1) + ptem2 * ptem
          else
             ad(i,k)   = ad(i,k)-au(i,k)
             ad(i,k+1) = 1.-al(i,k)
             a1(i,k+1) = u1(i,k+1)
             a2(i,k+1) = v1(i,k+1)
          endif
!
        enddo
      enddo

!
!     solve tridiagonal problem for momentum
!
      call tridi2(im,km,al,ad,au,a1,a2,au,a1,a2)
!
!     recover tendencies of momentum
!
!>  Finally, the tendencies are recovered from the tridiagonal solutions.
      do k = 1,km
         do i = 1,im
            utend = (a1(i,k)-u1(i,k))*rdt
            vtend = (a2(i,k)-v1(i,k))*rdt
            du(i,k)  = du(i,k)  + utend
            dv(i,k)  = dv(i,k)  + vtend
            dusfc(i) = dusfc(i) + conw*del(i,k)*utend
            dvsfc(i) = dvsfc(i) + conw*del(i,k)*vtend
            if(lssav .and. ldiag3d .and. .not.                          &
     &             flag_for_pbl_generic_tend) then
               du3dt_PBL(i,k) = du3dt_PBL(i,k) + utend*delt
               dv3dt_PBL(i,k) = dv3dt_PBL(i,k) + vtend*delt
            endif
!
!  for dissipative heating for ecmwf model
!
!           tem1 = 0.5*(a1(i,k)+u1(i,k))
!           tem2 = 0.5*(a2(i,k)+v1(i,k))
!           diss(i,k) = -(tem1*utend+tem2*vtend)
!           diss(i,k) = max(diss(i,k),0.)
!           ttend = diss(i,k) / cp
!           tau(i,k) = tau(i,k) + ttend
!
         enddo
      enddo
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      do i = 1, im
         hpbl(i) = hpblx(i)
         kpbl(i) = kpblx(i)
      enddo
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      return
      end subroutine hedmf_run
!> @}
!> @}

      end module hedmf

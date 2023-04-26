!>  \file samfshalcnv.f
!!  This file contains the Scale-Aware mass flux Shallow Convection scheme.

      module samfshalcnv

      use samfcnv_aerosols, only : samfshalcnv_aerosols
      use progsigma, only : progsigma_calc

      contains

      subroutine samfshalcnv_init(imfshalcnv, imfshalcnv_samf,          &
     &                            errmsg, errflg)

      integer,                   intent(in) :: imfshalcnv
      integer,                   intent(in) :: imfshalcnv_samf

      ! CCPP error handling
      character(len=*),          intent(out) :: errmsg
      integer,                   intent(out) :: errflg

      ! Consistency checks
      if (imfshalcnv/=imfshalcnv_samf) then
        write(errmsg,'(*(a))') 'Logic error: namelist choice of',       &
     &  ' shallow convection is different from SAMF'
        errflg = 1
        return
      end if
      end subroutine samfshalcnv_init

!> \defgroup SAMF_shal GFS saSAS Shallow Convection Module
!>  This subroutine contains the entirety of the SAMF shallow convection
!!  scheme.
!> @{
!!  This routine follows the \ref SAMFdeep quite closely, although it
!!  can be interpreted as only having the "static" and "feedback" control
!!  portions, since the "dynamic" control is not necessary to find the cloud
!!  base mass flux. The algorithm is simplified from SAMF deep convection by
!!  excluding convective downdrafts and being confined to operate below
!!  \f$p=0.7p_{sfc}\f$. Also, entrainment is both simpler and stronger in
!!  magnitude compared to the deep scheme.
!!
!! \section arg_table_samfshalcnv_run Argument Table
!! \htmlinclude samfshalcnv_run.html
!!
!!  \section gen_samfshalcnv GFS samfshalcnv General Algorithm
!!  -# Compute preliminary quantities needed for the static and feedback control portions of the algorithm.
!!  -# Perform calculations related to the updraft of the entraining/detraining cloud model ("static control").
!!  -# The cloud base mass flux is obtained using the cumulus updraft velocity averaged ove the whole cloud depth.
!!  -# Calculate the tendencies of the state variables (per unit cloud base mass flux) and the cloud base mass flux.
!!  -# For the "feedback control", calculate updated values of the state variables by multiplying the cloud base mass flux and the tendencies calculated per unit cloud base mass flux from the static control.
!!  \section det_samfshalcnv GFS samfshalcnv Detailed Algorithm
      subroutine samfshalcnv_run(im,km,itc,ntc,cliq,cp,cvap,            &
     &     eps,epsm1,fv,grav,hvap,rd,rv,                                &
     &     t0c,delt,ntk,ntr,delp,first_time_step,restart,               & 
     &     tmf,qmicro,progsigma,                                        &
     &     prslp,psp,phil,qtr,prevsq,q,q1,t1,u1,v1,fscav,               &
     &     rn,kbot,ktop,kcnv,islimsk,garea,                             &
     &     dot,ncloud,hpbl,ud_mf,dt_mf,cnvw,cnvc,                       &
     &     clam,c0s,c1,evef,pgcon,asolfac,hwrf_samfshal,                & 
     &     sigmain,sigmaout,errmsg,errflg)
!
      use machine , only : kind_phys
      use funcphys , only : fpvs

      implicit none
!
      integer, intent(in)  :: im, km, itc, ntc, ntk, ntr, ncloud
      integer, intent(in)  :: islimsk(:)
      real(kind=kind_phys), intent(in) :: cliq, cp, cvap,               &
     &   eps, epsm1, fv, grav, hvap, rd, rv, t0c
      real(kind=kind_phys), intent(in) ::  delt
      real(kind=kind_phys), intent(in) :: psp(:), delp(:,:),            &
     &   prslp(:,:), garea(:), hpbl(:), dot(:,:), phil(:,:),            &
     &   qmicro(:,:),tmf(:,:,:),prevsq(:,:),q(:,:)

      real(kind=kind_phys), intent(in) :: sigmain(:,:)
!
      real(kind=kind_phys), dimension(:), intent(in) :: fscav
      integer, intent(inout)  :: kcnv(:)
      ! DH* TODO - check dimensions of qtr, ntr+2 correct?  *DH
      real(kind=kind_phys), intent(inout) ::   qtr(:,:,:),              &
     &   q1(:,:), t1(:,:), u1(:,:), v1(:,:)
!
      integer, intent(out) :: kbot(:), ktop(:)
      real(kind=kind_phys), intent(out) :: rn(:),                       &
     &   cnvw(:,:), cnvc(:,:), ud_mf(:,:), dt_mf(:,:), sigmaout(:,:)
!
      real(kind=kind_phys), intent(in) :: clam,    c0s,     c1,         &
     &                     asolfac, evef, pgcon
      logical,          intent(in)  :: hwrf_samfshal,first_time_step,   &
     &     restart,progsigma
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
!
!  local variables
      integer              i,j,indx, k, kk, km1, n
      integer              kpbl(im)
!
      real(kind=kind_phys) clamd,   tkemx,   tkemn,   dtke
!
      real(kind=kind_phys) dellat,
     &                     c0l,     d0,
     &                     desdt,   dp,
     &                     dq,      dqsdp,   dqsdt,   dt,
     &                     dt2,     dtmax,   dtmin,
     &                     dxcrt,   dxcrtc0,
     &                     dv1h,    dv2h,    dv3h,
     &                     dz,      dz1,     e1,
     &                     el2orc,  elocp,   aafac,
     &                     cm,      cq,
     &                     es,      etah,    h1,      shevf,
!    &                     evfact,  evfactl,
     &                     fact1,   fact2,   factor,  dthk,
     &                     gamma,   pprime,  betaw,   tauadv,
     &                     qlk,     qrch,    qs,
     &                     rfact,   shear,   tfac,
     &                     val,     val1,    val2,
     &                     w1,      w1l,     w1s,     w2,
     &                     w2l,     w2s,     w3,      w3l,
     &                     w3s,     w4,      w4l,     w4s,
     &                     rho,     tem,     tem1,    tem2,
     &                     ptem,    ptem1
!
      integer              kb(im), kb1(im), kbcon(im), kbcon1(im),
     &                     ktcon(im), ktcon1(im), ktconn(im),
     &                     kbm(im), kmax(im)
!
      real(kind=kind_phys) aa1(im),     cina(im),
     &                     tkemean(im), clamt(im),
     &                     ps(im),      del(im,km), prsl(im,km),
     &                     umean(im),   advfac(im), gdx(im),
     &                     delhbar(im), delq(im),   delq2(im),
     &                     delqbar(im), delqev(im), deltbar(im),
!    &                     deltv(im),   dtconv(im), edt(im),
     &                     deltv(im),   dtconv(im),
     &                     pdot(im),    po(im,km),
     &                     qcond(im),   qevap(im),  hmax(im),
!    &                     rntot(im),   vshear(im),
     &                     rntot(im),
     &                     xlamud(im),  xmb(im),    xmbmax(im),
     &                     delebar(im,ntr),
     &                     delubar(im), delvbar(im)
!
      real(kind=kind_phys) c0(im), sfcpbl(im)
c
      real(kind=kind_phys) crtlame, crtlamd
!
      real(kind=kind_phys) cinpcr,  cinpcrmx,  cinpcrmn,
     &                     cinacr,  cinacrmx,  cinacrmn,
     &                     sfclfac, rhcrt
!
!  parameters for updraft velocity calculation
      real(kind=kind_phys) bet1,    cd1,     f1,      gam1,
!     &                     bb1,     bb2
     &                     bb1,     bb2,     wucb

cc

!  parameters for prognostic sigma closure
      real(kind=kind_phys) omega_u(im,km),zdqca(im,km),tmfq(im,km),
     &                     omegac(im),zeta(im,km),dbyo1(im,km),
     &                     sigmab(im),qadv(im,km)
      real(kind=kind_phys) gravinv,dxcrtas,invdelt
      logical flag_shallow
c  physical parameters
!     parameter(g=grav,asolfac=0.89)
!     parameter(g=grav)
!     parameter(elocp=hvap/cp,
!    &          el2orc=hvap*hvap/(rv*cp))
!     parameter(c0s=0.002,c1=5.e-4,d0=.01)
!     parameter(d0=.01)
      parameter(d0=.001)
!     parameter(c0l=c0s*asolfac)
!
! asolfac: aerosol-aware parameter based on Lim & Hong (2012)
!      asolfac= cx / c0s(=.002)
!      cx = min([-0.7 ln(Nccn) + 24]*1.e-4, c0s)
!      Nccn: CCN number concentration in cm^(-3)
!      Until a realistic Nccn is provided, Nccns are assumed
!      as Nccn=100 for sea and Nccn=1000 for land
!
      parameter(cm=1.0,cq=1.3)
!     parameter(fact1=(cvap-cliq)/rv,fact2=hvap/rv-fact1*t0c)
      parameter(clamd=0.1,tkemx=0.65,tkemn=0.05)
      parameter(dtke=tkemx-tkemn)
      parameter(dthk=25.,sfclfac=0.2,rhcrt=0.75)
      parameter(cinpcrmx=180.,cinpcrmn=120.)
!  shevf is an enhancing evaporation factor for shallow convection
      parameter(cinacrmx=-120.,shevf=2.0)
      parameter(dtmax=10800.,dtmin=600.)
      parameter(bet1=1.875,cd1=.506,f1=2.0,gam1=.5)
      parameter(betaw=.03,dxcrtc0=9.e3)
      parameter(h1=0.33333333)
!  progsigma
      parameter(dxcrtas=30.e3)
c  local variables and arrays
      real(kind=kind_phys) pfld(im,km),    to(im,km),     qo(im,km),
     &                     uo(im,km),      vo(im,km),     qeso(im,km),
     &                     ctr(im,km,ntr), ctro(im,km,ntr)
!  for aerosol transport
!     real(kind=kind_phys) qaero(im,km,ntc)
c  variables for tracer wet deposition,
      real(kind=kind_phys), dimension(im,km,ntc) :: chem_c, chem_pw,
     &  wet_dep
      real(kind=kind_phys), parameter :: escav   = 0.8 ! wet scavenging efficiency
!
!  for updraft velocity calculation
      real(kind=kind_phys) wu2(im,km),     buo(im,km),    drag(im,km),
     &                     wc(im)
!
!  for updraft fraction & scale-aware function
      real(kind=kind_phys) scaldfunc(im), sigmagfm(im)
!
c  cloud water
!     real(kind=kind_phys) qlko_ktcon(im), dellal(im,km), tvo(im,km),
      real(kind=kind_phys) qlko_ktcon(im), dellal(im,km),
     &                     dbyo(im,km),    zo(im,km),     xlamue(im,km),
     &                     rh(im,km),
     &                     heo(im,km),     heso(im,km),
     &                     dellah(im,km),  dellaq(im,km),
     &                     dellae(im,km,ntr),
     &                     dellau(im,km),  dellav(im,km), hcko(im,km),
     &                     ucko(im,km),    vcko(im,km),   qcko(im,km),
     &                     qrcko(im,km),   ecko(im,km,ntr),
     &                     ercko(im,km,ntr), eta(im,km),
     &                     zi(im,km),      pwo(im,km),    c0t(im,km),
     &                     sumx(im),       tx1(im),       cnvwt(im,km)
     &,                    rhbar(im)
!
!  variables for Total Variation Diminishing (TVD) flux-limiter scheme
!     on environmental subsidence and uplifting
!
      real(kind=kind_phys) q_diff(im,0:km-1), e_diff(im,0:km-1,ntr),
     &                     flxtvd(im,km-1)
      real(kind=kind_phys) rrkp, phkp
      real(kind=kind_phys) tsumn(im), tsump(im), rtnp(im)
!
      logical do_aerosols, totflg, cnvflg(im), flg(im)
!
      real(kind=kind_phys) tf, tcr, tcrf
      parameter (tf=233.16, tcr=263.16, tcrf=1.0/(tcr-tf))


c-----------------------------------------------------------------------
!
! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0

      gravinv = 1./grav
      invdelt = 1./delt

      elocp = hvap/cp
      el2orc = hvap*hvap/(rv*cp)

      fact1 = (cvap-cliq)/rv
      fact2 = hvap/rv-fact1*t0c

      if (.not.hwrf_samfshal) then
             cinacrmn=-80.
      endif

      if (progsigma) then
         dxcrt=10.e3
      else
         dxcrt=15.e3
      endif

c-----------------------------------------------------------------------
      if (.not.hwrf_samfshal) then
!>  ## Determine whether to perform aerosol transport
        do_aerosols = (itc > 0) .and. (ntc > 0) .and. (ntr > 0)
        if (do_aerosols) do_aerosols = (ntr >= itc + ntc - 3)
      endif
!
!************************************************************************
!     convert input Pa terms to Cb terms  -- Moorthi
!>  ## Compute preliminary quantities needed for the static and feedback control portions of the algorithm.
!>  - Convert input pressure terms to centibar units.
      ps   = psp   * 0.001
      prsl = prslp * 0.001
      del  = delp  * 0.001
!************************************************************************
!
      km1 = km - 1
c
c  initialize arrays
c
!>  - Initialize column-integrated and other single-value-per-column variable arrays.
!
      chem_c  = 0.
      chem_pw = 0.
      wet_dep = 0.
!
      if(hwrf_samfshal) then
       do i=1,im
        cnvflg(i) = .true.
        if(kcnv(i) == 1) cnvflg(i) = .false.
        if(cnvflg(i)) then
          kbot(i)=km+1
          ktop(i)=0
        endif
        sfcpbl(i) = sfclfac * hpbl(i)
        rn(i)=0.
        kbcon(i)=km
        ktcon(i)=1
        ktconn(i)=1
        kb(i)=km
        pdot(i) = 0.
        qlko_ktcon(i) = 0.
!       edt(i)  = 0.
        aa1(i)  = 0.
        cina(i) = 0.
!       vshear(i) = 0.
        advfac(i) = 0.
        gdx(i) = sqrt(garea(i))
        xmb(i) = 0.
          scaldfunc(i)=-1.0  ! wang initialized
          sigmagfm(i)=-1.0
       enddo

      else !gfs_samfshal
       do i=1,im
        cnvflg(i) = .true.
        if(kcnv(i) == 1) cnvflg(i) = .false.
        if(cnvflg(i)) then
          kbot(i)=km+1
          ktop(i)=0
        endif
        sfcpbl(i) = sfclfac * hpbl(i)
        rn(i)=0.
        kbcon(i)=km
        ktcon(i)=1
        ktconn(i)=1
        kb(i)=km
        pdot(i) = 0.
        qlko_ktcon(i) = 0.
!       edt(i)  = 0.0
        aa1(i)  = 0.
        cina(i) = 0.
!       vshear(i) = 0.
        gdx(i) = sqrt(garea(i))
        xmb(i) = 0.
       enddo
      endif
!!

!>  - Return to the calling routine if deep convection is present or the surface buoyancy flux is negative.
      totflg = .true.
      do i=1,im
        totflg = totflg .and. (.not. cnvflg(i))
      enddo
      if(totflg) return
!!
!>  - determine aerosol-aware rain conversion parameter over land
      do i=1,im
        if(islimsk(i) == 1) then
           c0(i) = c0s*asolfac
        else
           c0(i) = c0s
        endif
      enddo
!
!>  - determine scale-aware rain conversion parameter decreasing with decreasing grid size
      do i=1,im
        if(gdx(i) < dxcrtc0) then
          tem = gdx(i) / dxcrtc0 
          tem1 = tem**3
          c0(i) = c0(i) * tem1
        endif
      enddo
!
!>  - determine rain conversion parameter above the freezing level which exponentially decreases with decreasing temperature from Han et al.'s (2017) \cite han_et_al_2017 equation 8.
      do k = 1, km
        do i = 1, im
          if(t1(i,k) > 273.16) then
            c0t(i,k) = c0(i)
          else
            tem = d0 * (t1(i,k) - 273.16)
            tem1 = exp(tem)
            c0t(i,k) = c0(i) * tem1
          endif
        enddo
      enddo
!
!>  - Initialize convective cloud water and cloud cover to zero.
      do k = 1, km
        do i = 1, im
          cnvw(i,k) = 0.
          cnvc(i,k) = 0.
        enddo
      enddo
! hchuang code change
!>  - Initialize updraft mass fluxes to zero.
      do k = 1, km
        do i = 1, im
          ud_mf(i,k) = 0.
          dt_mf(i,k) = 0.
        enddo
      enddo
c
      dt2   = delt
!
c  model tunable parameters are all here
      aafac   = .1
!     evfact  = 0.3
!     evfactl = 0.3
!
      crtlame = 1.0e-4
      crtlamd = 3.0e-4
!
      w1l     = -8.e-3
      w2l     = -4.e-2
      w3l     = -5.e-3
      w4l     = -5.e-4
      w1s     = -2.e-4
      w2s     = -2.e-3
      w3s     = -1.e-3
      w4s     = -2.e-5
c
c  define top layer for search of the downdraft originating layer
c  and the maximum thetae for updraft
c
!>  - Determine maximum indices for the parcel starting point (kbm) and cloud top (kmax).
      do i=1,im
        kbm(i)   = km
        kmax(i)  = km
        tx1(i)   = 1.0 / ps(i)
      enddo
!
      do k = 1, km
        do i=1,im
          if (prsl(i,k)*tx1(i) > 0.70) kbm(i)   = k + 1
          if (prsl(i,k)*tx1(i) > 0.60) kmax(i)  = k + 1
        enddo
      enddo
      do i=1,im
        kbm(i)   = min(kbm(i),kmax(i))
      enddo
c
c  hydrostatic height assume zero terr and compute
c  updraft entrainment rate as an inverse function of height
c
!>  - Calculate hydrostatic height at layer centers assuming a flat surface (no terrain) from the geopotential.
      do k = 1, km
        do i=1,im
          zo(i,k) = phil(i,k) / grav
        enddo
      enddo
!>  - Calculate interface height
      if(hwrf_samfshal) then
       do k = 1, km1
        do i=1,im
          zi(i,k) = 0.5*(zo(i,k)+zo(i,k+1))
          xlamue(i,k) = clam / zi(i,k)
        enddo
       enddo
       do i=1,im
        xlamue(i,km) = xlamue(i,km1)
       enddo
      else
       do k = 1, km1
        do i=1,im
          zi(i,k) = 0.5*(zo(i,k)+zo(i,k+1))
        enddo
       enddo
      endif
c
c  pbl height
c
!>  - Find the index for the PBL top using the PBL height; enforce that it is lower than the maximum parcel starting level.
      do i=1,im
        flg(i) = cnvflg(i)
        kpbl(i)= 1
      enddo
      do k = 2, km1
        do i=1,im
          if (flg(i) .and. zo(i,k) <= hpbl(i)) then
            kpbl(i) = k
          else
            flg(i) = .false.
          endif
        enddo
      enddo
      do i=1,im
        kpbl(i)= min(kpbl(i),kbm(i))
      enddo
c
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c   convert surface pressure to mb from cb
c
!>  - Convert prsl from centibar to millibar, set normalized mass flux to 1, cloud properties to 0, and save model state variables (after advection/turbulence).
      do k = 1, km
        do i = 1, im
          if (cnvflg(i) .and. k <= kmax(i)) then
            pfld(i,k) = prsl(i,k) * 10.0
            eta(i,k)  = 1.
            rh(i,k)   = 0.
            hcko(i,k) = 0.
            qcko(i,k) = 0.
            qrcko(i,k)= 0.
            ucko(i,k) = 0.
            vcko(i,k) = 0.
            dbyo(i,k) = 0.
            pwo(i,k)  = 0.
            dellal(i,k) = 0.
            to(i,k)   = t1(i,k)
            qo(i,k)   = q1(i,k)
            uo(i,k)   = u1(i,k)
            vo(i,k)   = v1(i,k)
!           uo(i,k)   = u1(i,k) * rcs(i)
!           vo(i,k)   = v1(i,k) * rcs(i)
            wu2(i,k)  = 0.
            buo(i,k)  = 0.
            drag(i,k) = 0.
            cnvwt(i,k) = 0.
          endif
        enddo
      enddo

      do i = 1,im
         omegac(i)=0.
      enddo

      do k = 1, km
         do i = 1, im
            dbyo1(i,k)=0.
            zdqca(i,k)=0.
            omega_u(i,k)=0.
            zeta(i,k)=1.0
         enddo
      enddo
!
!  initialize tracer variables
!
      if (.not.hwrf_samfshal) then
        do n = 3, ntr+2
          kk = n-2
        do k = 1, km
          do i = 1, im
            if (cnvflg(i) .and. k <= kmax(i)) then
              ctr(i,k,kk)  = qtr(i,k,n)
              ctro(i,k,kk) = qtr(i,k,n)
              ecko(i,k,kk) = 0.
              ercko(i,k,kk) = 0.
            endif
          enddo
        enddo
        enddo
      endif
!>  - Calculate saturation specific humidity and enforce minimum moisture values.
      do k = 1, km
        do i=1,im
          if (cnvflg(i) .and. k <= kmax(i)) then
            qeso(i,k) = 0.01 * fpvs(to(i,k))      ! fpvs is in pa
            qeso(i,k) = eps * qeso(i,k) / (pfld(i,k) + epsm1*qeso(i,k))
            val1      =             1.e-8
            qeso(i,k) = max(qeso(i,k), val1)
            val2      =           1.e-10
            qo(i,k)   = max(qo(i,k), val2 )
!           qo(i,k)   = min(qo(i,k),qeso(i,k))
!           tvo(i,k)  = to(i,k) + fv * to(i,k) * qo(i,k)
          endif
        enddo
      enddo
c
c  compute moist static energy
c
!>  - Calculate moist static energy (heo) and saturation moist static energy (heso).
      do k = 1, km
        do i=1,im
          if (cnvflg(i) .and. k <= kmax(i)) then
!           tem       = grav * zo(i,k) + cp * to(i,k)
            tem       = phil(i,k) + cp * to(i,k)
            heo(i,k)  = tem  + hvap * qo(i,k)
            heso(i,k) = tem  + hvap * qeso(i,k)
c           heo(i,k)  = min(heo(i,k),heso(i,k))
          endif
        enddo
      enddo
c
c  determine level with largest moist static energy within pbl
c  this is the level where updraft starts
c
!> ## Perform calculations related to the updraft of the entraining/detraining cloud model ("static control").
!> - Find the index for a level of sfclfac*hpbl which is initial guess for the parcel starting level.
      do i=1,im
        flg(i) = cnvflg(i)
        kb1(i) = 1
      enddo
      do k = 1, km1
        do i=1,im
          if (flg(i) .and. zo(i,k) <= sfcpbl(i)) then
            kb1(i) = k
          else
            flg(i) = .false.
          endif
        enddo
      enddo
      do i=1,im
        kb1(i) = min(kb1(i),kpbl(i))
      enddo
c
!> - Search in the PBL for the level of maximum moist static energy to start the ascending parcel.
      do i=1,im
         if (cnvflg(i)) then
            hmax(i) = heo(i,kb1(i))
            kb(i) = kb1(i)
         endif
      enddo
      do k = 2, km
        do i=1,im
          if(cnvflg(i) .and. (k > kb1(i) .and. k <= kpbl(i))) then
            if(heo(i,k) > hmax(i)) then
              kb(i)   = k
              hmax(i) = heo(i,k)
            endif
          endif
        enddo
      enddo
c
!> - Calculate the temperature, water vapor mixing ratio, and pressure at interface levels.
      do k = 1, km1
        do i=1,im
          if (cnvflg(i) .and. k <= kmax(i)-1) then
            dz      = .5 * (zo(i,k+1) - zo(i,k))
            dp      = .5 * (pfld(i,k+1) - pfld(i,k))
            es      = 0.01 * fpvs(to(i,k+1))      ! fpvs is in pa
            pprime  = pfld(i,k+1) + epsm1 * es
            qs      = eps * es / pprime
            dqsdp   = - qs / pprime
            desdt   = es * (fact1 / to(i,k+1) + fact2 / (to(i,k+1)**2))
            dqsdt   = qs * pfld(i,k+1) * desdt / (es * pprime)
            gamma   = el2orc * qeso(i,k+1) / (to(i,k+1)**2)
            dt      = (grav*dz + hvap*dqsdp*dp) / (cp*(1. + gamma))
            dq      = dqsdt * dt + dqsdp * dp
            to(i,k) = to(i,k+1) + dt
            qo(i,k) = qo(i,k+1) + dq
            po(i,k) = .5 * (pfld(i,k) + pfld(i,k+1))
          endif
        enddo
      enddo
!
!> - Recalculate saturation specific humidity, moist static energy, saturation moist static energy, and horizontal momentum on interface levels. Enforce minimum specific humidity.
      do k = 1, km1
        do i=1,im
          if (cnvflg(i) .and. k <= kmax(i)-1) then
            qeso(i,k) = 0.01 * fpvs(to(i,k))      ! fpvs is in pa
            qeso(i,k) = eps * qeso(i,k) / (po(i,k) + epsm1*qeso(i,k))
            val1      =             1.e-8
            qeso(i,k) = max(qeso(i,k), val1)
            val2      =           1.e-10
            qo(i,k)   = max(qo(i,k), val2 )
!           qo(i,k)   = min(qo(i,k),qeso(i,k))
            rh(i,k)   = min(qo(i,k)/qeso(i,k), 1.)
            heo(i,k)  = .5 * grav * (zo(i,k) + zo(i,k+1)) +
     &                  cp * to(i,k) + hvap * qo(i,k)
            heso(i,k) = .5 * grav * (zo(i,k) + zo(i,k+1)) +
     &                  cp * to(i,k) + hvap * qeso(i,k)
            uo(i,k)   = .5 * (uo(i,k) + uo(i,k+1))
            vo(i,k)   = .5 * (vo(i,k) + vo(i,k+1))
          endif
        enddo
      enddo

      if (.not.hwrf_samfshal) then
       do n = 1, ntr
       do k = 1, km1
        do i=1,im
          if (cnvflg(i) .and. k <= kmax(i)-1) then
            ctro(i,k,n) = .5 * (ctro(i,k,n) + ctro(i,k+1,n))
          endif
        enddo
       enddo
       enddo
      endif
c
c  look for the level of free convection as cloud base
c
!> - Search below the index "kbm" for the level of free convection (LFC) where the condition \f$h_b > h^*\f$ is first met, where \f$h_b, h^*\f$ are the state moist static energy at the parcel's starting level and saturation moist static energy, respectively. Set "kbcon" to the index of the LFC.
      do i=1,im
        flg(i)   = cnvflg(i)
        if(flg(i)) kbcon(i) = kmax(i)
      enddo
      do k = 2, km1
        do i=1,im
          if (flg(i) .and. k < kbm(i)) then
            if(k > kb(i) .and. heo(i,kb(i)) > heso(i,k)) then
              kbcon(i) = k
              flg(i)   = .false.
            endif
          endif
        enddo
      enddo
c
      do i=1,im
        if(cnvflg(i)) then
          if(kbcon(i) == kmax(i)) cnvflg(i) = .false.
        endif
      enddo
!!
!> - If no LFC, return to the calling routine without modifying state variables.
      totflg = .true.
      do i=1,im
        totflg = totflg .and. (.not. cnvflg(i))
      enddo
      if(totflg) return
!!
!> - Determine the vertical pressure velocity at the LFC. After Han and Pan (2011) \cite han_and_pan_2011 , determine the maximum pressure thickness between a parcel's starting level and the LFC. If a parcel doesn't reach the LFC within the critical thickness, then the convective inhibition is deemed too great for convection to be triggered, and the subroutine returns to the calling routine without modifying the state variables.
      do i=1,im
        if(cnvflg(i)) then
!         pdot(i)  = 10.* dot(i,kbcon(i))
          pdot(i)  = 0.01 * dot(i,kbcon(i)) ! Now dot is in Pa/s
        endif
      enddo
c
c   turn off convection if pressure depth between parcel source level
c      and cloud base is larger than a critical value, cinpcr
c
      do i=1,im
        if(cnvflg(i)) then
          if(islimsk(i) == 1) then
            w1 = w1l
            w2 = w2l
            w3 = w3l
            w4 = w4l
          else
            w1 = w1s
            w2 = w2s
            w3 = w3s
            w4 = w4s
          endif
          if(pdot(i) <= w4) then
            tem = (pdot(i) - w4) / (w3 - w4)
          elseif(pdot(i) >= -w4) then
            tem = - (pdot(i) + w4) / (w4 - w3)
          else
            tem = 0.
          endif
          val1    =            -1.
          tem = max(tem,val1)
          val2    =             1.
          tem = min(tem,val2)
          ptem = 1. - tem
          ptem1= .5*(cinpcrmx-cinpcrmn)
          cinpcr = cinpcrmx - ptem * ptem1
          tem1 = pfld(i,kb(i)) - pfld(i,kbcon(i))

          if(tem1 > cinpcr) then
             cnvflg(i) = .false.
          endif
        endif
      enddo
!!
      totflg = .true.
      do i=1,im
        totflg = totflg .and. (.not. cnvflg(i))
      enddo
      if(totflg) return
!
! re-define kb & kbcon
!
      do i=1,im
         if (cnvflg(i)) then
            hmax(i) = heo(i,1)
            kb(i) = 1
         endif
      enddo
      do k = 2, km
        do i=1,im
          if (cnvflg(i) .and. k <= kpbl(i)) then
            if(heo(i,k) > hmax(i)) then
              kb(i)   = k
              hmax(i) = heo(i,k)
            endif
          endif
        enddo
      enddo
!
      do i=1,im
        flg(i)   = cnvflg(i)
        if(flg(i)) kbcon(i) = kmax(i)
      enddo
      do k = 2, km1
        do i=1,im
          if (flg(i) .and. k < kbm(i)) then
            if(k > kb(i) .and. heo(i,kb(i)) > heso(i,k)) then
              kbcon(i) = k
              flg(i)   = .false.
            endif
          endif
        enddo
      enddo
!
      do i=1,im
        if(cnvflg(i)) then
          if(kbcon(i) == kmax(i)) cnvflg(i) = .false.
        endif
      enddo
!!
      totflg = .true.
      do i=1,im
        totflg = totflg .and. (.not. cnvflg(i))
      enddo
      if(totflg) return
!!
      do i=1,im
        if(cnvflg(i)) then
!         pdot(i)  = 10.* dot(i,kbcon(i))
          pdot(i)  = 0.01 * dot(i,kbcon(i)) ! Now dot is in Pa/s
        endif
      enddo
!
!> - if the mean relative humidity in the subcloud layers is less than a threshold value (rhcrt), convection is not triggered.
!
      do i = 1, im
        rhbar(i) = 0.
        sumx(i) = 0.
      enddo
      do k = 1, km1
        do i = 1, im
          if (cnvflg(i)) then
            if(k >= kb(i) .and. k < kbcon(i)) then
              dz = zo(i,k+1) - zo(i,k)
              rhbar(i) = rhbar(i) + rh(i,k) * dz
              sumx(i) = sumx(i) + dz
            endif
          endif
        enddo
      enddo
      do i= 1, im
        if(cnvflg(i)) then
          rhbar(i) = rhbar(i) / sumx(i)
          if(rhbar(i) < rhcrt) then
            cnvflg(i) = .false.
          endif
        endif
      enddo
!!
      totflg = .true.
      do i=1,im
        totflg = totflg .and. (.not. cnvflg(i))
      enddo
      if(totflg) return
!!
!
! turbulent entrainment rate assumed to be proportional
!   to subcloud mean TKE
!
!c
!c  specify the detrainment rate for the updrafts
!c
      if (hwrf_samfshal) then
       do i = 1, im
        if(cnvflg(i)) then
          xlamud(i) = xlamue(i,kbcon(i))
!         xlamud(i) = crtlamd
        endif
       enddo
      else
      if(ntk > 0) then
        do i= 1, im
          if(cnvflg(i)) then
            sumx(i) = 0.
            tkemean(i) = 0.
          endif
        enddo

        do k = 1, km1
          do i = 1, im
            if(cnvflg(i)) then
              if(k >= kb(i) .and. k < kbcon(i)) then
                dz = zo(i,k+1) - zo(i,k)
                tem = 0.5 * (qtr(i,k,ntk)+qtr(i,k+1,ntk))
                tkemean(i) = tkemean(i) + tem * dz
                sumx(i) = sumx(i) + dz
              endif
            endif
          enddo
        enddo
!
        do i= 1, im
          if(cnvflg(i)) then
             tkemean(i) = tkemean(i) / sumx(i)
             if(tkemean(i) > tkemx) then
               clamt(i) = clam + clamd
             else if(tkemean(i) < tkemn) then
               clamt(i) = clam - clamd
             else
               tem = tkemx - tkemean(i)
               tem1 = 1. - 2. *  tem / dtke
               clamt(i) = clam + clamd * tem1
             endif
          endif
        enddo
!
      else
!
        do i= 1, im
          if(cnvflg(i)) then
            clamt(i)  = clam
          endif
        enddo
!
      endif
!!
!
!  assume updraft entrainment rate
!     is an inverse function of height
!
      do k = 1, km1
        do i=1,im
          if(cnvflg(i)) then
            dz = zo(i,k+1) - zo(i,k)
            xlamue(i,k) = clamt(i) / (zi(i,k) + dz)
            xlamue(i,k) = max(xlamue(i,k), crtlame)
          endif
        enddo
      enddo
      do i=1,im
        if(cnvflg(i)) then
          xlamue(i,km) = xlamue(i,km1)
        endif
      enddo
c
c  specify the detrainment rate for the updrafts
c
!! (The updraft detrainment rate is set constant and equal to the entrainment rate at cloud base.)
!!
!> - The updraft detrainment rate is vertically constant and proportional to clamt
      do i = 1, im
        if(cnvflg(i)) then
!         xlamud(i) = xlamue(i,kbcon(i))
!         xlamud(i) = crtlamd
          xlamud(i) = 0.001 * clamt(i)
        endif
      enddo
      endif    ! hwrf_samfshal
c
c  determine updraft mass flux for the subcloud layers
c
!> - Calculate the normalized mass flux for subcloud and in-cloud layers according to Pan and Wu (1995) \cite pan_and_wu_1995 equation 1:
!!  \f[
!!  \frac{1}{\eta}\frac{\partial \eta}{\partial z} = \lambda_e - \lambda_d
!!  \f]
!!  where \f$\eta\f$ is the normalized mass flux, \f$\lambda_e\f$ is the entrainment rate and \f$\lambda_d\f$ is the detrainment rate. The normalized mass flux increases upward below the cloud base and decreases upward above.
      do k = km1, 1, -1
        do i = 1, im
          if (cnvflg(i)) then
            if(k < kbcon(i) .and. k >= kb(i)) then
              dz       = zi(i,k+1) - zi(i,k)
              ptem     = 0.5*(xlamue(i,k)+xlamue(i,k+1))-xlamud(i)
              eta(i,k) = eta(i,k+1) / (1. + ptem * dz)
            endif
          endif
        enddo
      enddo
c
c  compute mass flux above cloud base
c
      do i = 1, im
        flg(i) = cnvflg(i)
      enddo
      do k = 2, km1
        do i = 1, im
         if(flg(i))then
           if(k > kbcon(i) .and. k < kmax(i)) then
              dz       = zi(i,k) - zi(i,k-1)
              ptem     = 0.5*(xlamue(i,k)+xlamue(i,k-1))-xlamud(i)
              eta(i,k) = eta(i,k-1) * (1 + ptem * dz)
              if(eta(i,k) <= 0.) then
                kmax(i) = k
                ktconn(i) = k
                kbm(i) = min(kbm(i),kmax(i))
                flg(i) = .false.
              endif
           endif
         endif
        enddo
      enddo
c
c  compute updraft cloud property
c
!> - Set cloud properties equal to the state variables at updraft starting level (kb).
      do i = 1, im
        if(cnvflg(i)) then
          indx         = kb(i)
          hcko(i,indx) = heo(i,indx)
          ucko(i,indx) = uo(i,indx)
          vcko(i,indx) = vo(i,indx)
        endif
      enddo
!  for tracers
      if (.not. hwrf_samfshal) then
      do n = 1, ntr
        do i = 1, im
          if(cnvflg(i)) then
            indx = kb(i)
            ecko(i,indx,n) = ctro(i,indx,n)
            ercko(i,indx,n) = ctro(i,indx,n)
          endif
        enddo
      enddo
      endif
c
!  cm is an enhancement factor in entrainment rates for momentum
!
!> - Calculate the cloud properties as a parcel ascends, modified by entrainment and detrainment. Discretization follows Appendix B of Grell (1993) \cite grell_1993 . Following Han and Pan (2006) \cite han_and_pan_2006, the convective momentum transport is reduced by the convection-induced pressure gradient force by the constant "pgcon", currently set to 0.55 after Zhang and Wu (2003) \cite zhang_and_wu_2003 .
      do k = 2, km1
        do i = 1, im
          if (cnvflg(i)) then
            if(k > kb(i) .and. k < kmax(i)) then
              dz   = zi(i,k) - zi(i,k-1)
              tem  = 0.5 * (xlamue(i,k)+xlamue(i,k-1)) * dz
              tem1 = 0.5 * xlamud(i) * dz
              factor = 1. + tem - tem1
              hcko(i,k) = ((1.-tem1)*hcko(i,k-1)+tem*0.5*
     &                     (heo(i,k)+heo(i,k-1)))/factor
              dbyo(i,k) = hcko(i,k) - heso(i,k)
!
              tem  = 0.5 * cm * tem
              factor = 1. + tem
              ptem = tem + pgcon
              ptem1= tem - pgcon
              ucko(i,k) = ((1.-tem)*ucko(i,k-1)+ptem*uo(i,k)
     &                     +ptem1*uo(i,k-1))/factor
              vcko(i,k) = ((1.-tem)*vcko(i,k-1)+ptem*vo(i,k)
     &                     +ptem1*vo(i,k-1))/factor
            endif
          endif
        enddo
      enddo

      if (.not.hwrf_samfshal) then
       if (do_aerosols) then
         kk = itc -3
       else
         kk = ntr
       endif
       do n = 1, kk
       do k = 2, km1
        do i = 1, im
          if (cnvflg(i)) then
            if(k > kb(i) .and. k < kmax(i)) then
              dz   = zi(i,k) - zi(i,k-1)
              tem  = 0.25 * (xlamue(i,k)+xlamue(i,k-1)) * dz
              tem  = cq * tem
              factor = 1. + tem
              ecko(i,k,n) = ((1.-tem)*ecko(i,k-1,n)+tem*
     &                     (ctro(i,k,n)+ctro(i,k-1,n)))/factor
              ercko(i,k,n) = ecko(i,k,n)
            endif
          endif
        enddo
       enddo
       enddo
       if (do_aerosols) then
         do n = 1, ntc
           kk = n + itc -3
           do k = 2, km1
             do i = 1, im
               if (cnvflg(i)) then
                 if(k > kb(i) .and. k < kmax(i)) then
                   dz = zi(i,k) - zi(i,k-1)
                   tem  = 0.25 * (xlamue(i,k)+xlamue(i,k-1)) * dz
                   tem  = cq * tem
                   factor = 1. + tem
                   ecko(i,k,kk) = ((1. - tem) * ecko(i,k-1,kk) + tem *
     &                     (ctro(i,k,kk) + ctro(i,k-1,kk))) / factor
                   ercko(i,k,kk) = ecko(i,k,kk)
                   chem_c(i,k,n) = escav * fscav(n) * ecko(i,k,kk)
                   tem = chem_c(i,k,n) / (1. + c0t(i,k) * dz)
                   chem_pw(i,k,n) = c0t(i,k) * dz * tem * eta(i,k-1)
                   ecko(i,k,kk) = tem + ecko(i,k,kk) - chem_c(i,k,n)
                 endif
               endif
             enddo
           enddo
         enddo
       endif
      endif
c
c   taking account into convection inhibition due to existence of
c    dry layers below cloud base
c
!> - With entrainment, recalculate the LFC as the first level where buoyancy is positive. The difference in pressure levels between LFCs calculated with/without entrainment must be less than a threshold (currently 25 hPa). Otherwise, convection is inhibited and the scheme returns to the calling routine without modifying the state variables. This is the subcloud dryness trigger modification discussed in Han and Pan (2011) \cite han_and_pan_2011.
      do i=1,im
        flg(i) = cnvflg(i)
        kbcon1(i) = kmax(i)
      enddo
      do k = 2, km1
      do i=1,im
        if (flg(i) .and. k < kbm(i)) then
          if(k >= kbcon(i) .and. dbyo(i,k) > 0.) then
            kbcon1(i) = k
            flg(i)    = .false.
          endif
        endif
      enddo
      enddo
      do i=1,im
        if(cnvflg(i)) then
          if(kbcon1(i) == kmax(i)) cnvflg(i) = .false.
        endif
      enddo
      do i=1,im
        if(cnvflg(i)) then
          tem = pfld(i,kbcon(i)) - pfld(i,kbcon1(i))
          if(tem > dthk) then
             cnvflg(i) = .false.
          endif
        endif
      enddo
!!
      totflg = .true.
      do i = 1, im
        totflg = totflg .and. (.not. cnvflg(i))
      enddo
      if(totflg) return
!!
c
c  calculate convective inhibition
c
!> - Calculate additional trigger condition of the convective inhibition (CIN) according to Han et al.'s (2017) \cite han_et_al_2017 equation 13.
      do k = 2, km1
        do i = 1, im
          if (cnvflg(i)) then
            if(k > kb(i) .and. k < kbcon1(i)) then
              dz1 = zo(i,k+1) - zo(i,k)
              gamma = el2orc * qeso(i,k) / (to(i,k)**2)
              rfact =  1. + fv * cp * gamma
     &                 * to(i,k) / hvap
              cina(i) = cina(i) +
!    &                 dz1 * eta(i,k) * (grav / (cp * to(i,k)))
     &                 dz1 * (grav / (cp * to(i,k)))
     &                 * dbyo(i,k) / (1. + gamma)
     &                 * rfact
              val = 0.
              cina(i) = cina(i) +
!    &                 dz1 * eta(i,k) * grav * fv *
     &                 dz1 * grav * fv *
     &                 max(val,(qeso(i,k) - qo(i,k)))
            endif
          endif
        enddo
      enddo
!> - Turn off convection if the CIN is less than a critical value (cinacr) which is inversely proportional to the large-scale vertical velocity.

      if (hwrf_samfshal) then
       do i = 1, im
        if(cnvflg(i)) then
          cinacr = cinacrmx
          if(cina(i) < cinacr) cnvflg(i) = .false.
        endif
       enddo
      else
       do i = 1, im
        if(cnvflg(i)) then
          if(islimsk(i) == 1) then
            w1 = w1l
            w2 = w2l
            w3 = w3l
            w4 = w4l
          else
            w1 = w1s
            w2 = w2s
            w3 = w3s
            w4 = w4s
          endif
          if(pdot(i) <= w4) then
            tem = (pdot(i) - w4) / (w3 - w4)
          elseif(pdot(i) >= -w4) then
            tem = - (pdot(i) + w4) / (w4 - w3)
          else
            tem = 0.
          endif

          val1    =            -1.
          tem = max(tem,val1)
          val2    =             1.
          tem = min(tem,val2)
          tem = 1. - tem
          tem1= .5*(cinacrmx-cinacrmn)
          cinacr = cinacrmx - tem * tem1
          if(cina(i) < cinacr) cnvflg(i) = .false.
         endif
       enddo
      endif
!!
      totflg = .true.
      do i=1,im
        totflg = totflg .and. (.not. cnvflg(i))
      enddo
      if(totflg) return
!!
c
c  determine first guess cloud top as the level of zero buoyancy
c    limited to the level of P/Ps=0.7
c
!> - Calculate the cloud top as the first level where parcel buoyancy becomes negative; the maximum possible value is at \f$p=0.7p_{sfc}\f$.
      do i = 1, im
        flg(i) = cnvflg(i)
        if(flg(i)) ktcon(i) = kbm(i)
      enddo
      do k = 2, km1
      do i=1,im
        if (flg(i) .and. k < kbm(i)) then
          if(k > kbcon1(i) .and. dbyo(i,k) < 0.) then
             ktcon(i) = k
             flg(i)   = .false.
          endif
        endif
      enddo
      enddo
c
c  specify upper limit of mass flux at cloud base
c
!> - Calculate the maximum value of the cloud base mass flux using the CFL-criterion-based formula of Han and Pan (2011) \cite han_and_pan_2011, equation 7.
      do i = 1, im
        if(cnvflg(i)) then
          k = kbcon(i)
          dp = 1000. * del(i,k)
          xmbmax(i) = dp / (grav * dt2)
        endif
      enddo
c
c  compute cloud moisture property and precipitation
c
!> - Set cloud moisture property equal to the enviromental moisture at updraft starting level (kb).
      do i = 1, im
        if (cnvflg(i)) then
          aa1(i) = 0.
          qcko(i,kb(i)) = qo(i,kb(i))
          qrcko(i,kb(i)) = qo(i,kb(i))
        endif
      enddo
!> - Calculate the moisture content of the entraining/detraining parcel (qcko) and the value it would have if just saturated (qrch), according to equation A.14 in Grell (1993) \cite grell_1993 . Their difference is the amount of convective cloud water (qlk = rain + condensate). Determine the portion of convective cloud water that remains suspended and the portion that is converted into convective precipitation (pwo). Calculate and save the negative cloud work function (aa1) due to water loading. Above the level of minimum moist static energy, some of the cloud water is detrained into the grid-scale cloud water from every cloud layer with a rate of 0.0005 \f$m^{-1}\f$ (dellal).
      do k = 2, km1
        do i = 1, im
          if (cnvflg(i)) then
            if(k > kb(i) .and. k < ktcon(i)) then
              dz    = zi(i,k) - zi(i,k-1)
              gamma = el2orc * qeso(i,k) / (to(i,k)**2)
              qrch = qeso(i,k)
     &             + gamma * dbyo(i,k) / (hvap * (1. + gamma))
cj
              tem  = 0.25 * (xlamue(i,k)+xlamue(i,k-1)) * dz
              tem  = cq * tem
              factor = 1. + tem
              qcko(i,k) = ((1.-tem)*qcko(i,k-1)+tem*
     &                     (qo(i,k)+qo(i,k-1)))/factor
              qrcko(i,k) = qcko(i,k)
cj
              dq = eta(i,k) * (qcko(i,k) - qrch)
c
!             rhbar(i) = rhbar(i) + qo(i,k) / qeso(i,k)
c
c  below lfc check if there is excess moisture to release latent heat
c
              if(k >= kbcon(i) .and. dq > 0.) then
                etah = .5 * (eta(i,k) + eta(i,k-1))
                dp = 1000. * del(i,k)
                if(ncloud > 0) then
                  ptem = c0t(i,k) + c1
                  qlk = dq / (eta(i,k) + etah * ptem * dz)
                  dellal(i,k) = etah * c1 * dz * qlk * grav / dp
                else
                  qlk = dq / (eta(i,k) + etah * c0t(i,k) * dz)
                endif
                buo(i,k) = buo(i,k) - grav * qlk
                qcko(i,k)= qlk + qrch
                pwo(i,k) = etah * c0t(i,k) * dz * qlk
                cnvwt(i,k) = etah * qlk * grav / dp
                zdqca(i,k)=dq/eta(i,k)
              endif
!
!  compute buoyancy and drag for updraft velocity
!
              if(k >= kbcon(i)) then
                rfact =  1. + fv * cp * gamma
     &                   * to(i,k) / hvap
                buo(i,k) = buo(i,k) + (grav / (cp * to(i,k)))
     &                   * dbyo(i,k) / (1. + gamma)
     &                   * rfact
                val = 0.
                buo(i,k) = buo(i,k) + grav * fv *
     &                     max(val,(qeso(i,k) - qo(i,k)))
                drag(i,k) = max(xlamue(i,k),xlamud(i))
              endif
!
            endif
          endif
        enddo
      enddo
c
c  calculate cloud work function
c
!     do k = 2, km1
!       do i = 1, im
!         if (cnvflg(i)) then
!           if(k >= kbcon(i) .and. k < ktcon(i)) then
!             dz1 = zo(i,k+1) - zo(i,k)
!             gamma = el2orc * qeso(i,k) / (to(i,k)**2)
!             rfact =  1. + fv * cp * gamma
!    &                 * to(i,k) / hvap
!             aa1(i) = aa1(i) +
!!   &                 dz1 * eta(i,k) * (grav / (cp * to(i,k)))
!    &                 dz1 * (grav / (cp * to(i,k)))
!    &                 * dbyo(i,k) / (1. + gamma)
!    &                 * rfact
!             val = 0.
!             aa1(i) = aa1(i) +
!!   &                 dz1 * eta(i,k) * grav * fv *
!    &                 dz1 * grav * fv *
!    &                 max(val,(qeso(i,k) - qo(i,k)))
!           endif
!         endif
!       enddo
!     enddo
!     do i = 1, im
!       if(cnvflg(i) .and. aa1(i) <= 0.) cnvflg(i) = .false.
!     enddo
!
!  calculate cloud work function
!
!> - Calculate the cloud work function according to Pan and Wu (1995) \cite pan_and_wu_1995 equation 4:
!!  \f[
!!  A_u=\int_{z_0}^{z_t}\frac{g}{c_pT(z)}\frac{\eta}{1 + \gamma}[h(z)-h^*(z)]dz
!!  \f]
!! (discretized according to Grell (1993) \cite grell_1993 equation B.10 using B.2 and B.3 of Arakawa and Schubert (1974) \cite arakawa_and_schubert_1974 and assuming \f$\eta=1\f$) where \f$A_u\f$ is the updraft cloud work function, \f$z_0\f$ and \f$z_t\f$ are cloud base and cloud top, respectively, \f$\gamma = \frac{L}{c_p}\left(\frac{\partial \overline{q_s}}{\partial T}\right)_p\f$ and other quantities are previously defined.
      do i = 1, im
        if (cnvflg(i)) then
          aa1(i) = 0.
        endif
      enddo
      do k = 2, km1
        do i = 1, im
          if (cnvflg(i)) then
            if(k >= kbcon(i) .and. k < ktcon(i)) then
              dz1 = zo(i,k+1) - zo(i,k)
              aa1(i) = aa1(i) + buo(i,k) * dz1
              dbyo1(i,k) = hcko(i,k) - heso(i,k)
            endif
          endif
        enddo
      enddo
      do i = 1, im
        if(cnvflg(i) .and. aa1(i) <= 0.) cnvflg(i) = .false.
      enddo
!!
!> - If the updraft cloud work function is negative, convection does not occur, and the scheme returns to the calling routine.
      totflg = .true.
      do i=1,im
        totflg = totflg .and. (.not. cnvflg(i))
      enddo
      if(totflg) return
!!
c
c  estimate the onvective overshooting as the level
c    where the [aafac * cloud work function] becomes zero,
c    which is the final cloud top
c    limited to the level of P/Ps=0.7
c
!> - Continue calculating the cloud work function past the point of neutral buoyancy to represent overshooting according to Han and Pan (2011) \cite han_and_pan_2011 . Convective overshooting stops when \f$ cA_u < 0\f$ where \f$c\f$ is currently 10%, or when 10% of the updraft cloud work function has been consumed by the stable buoyancy force. Overshooting is also limited to the level where \f$p=0.7p_{sfc}\f$.
      do i = 1, im
        if (cnvflg(i)) then
          aa1(i) = aafac * aa1(i)
        endif
      enddo
c
      do i = 1, im
        flg(i) = cnvflg(i)
        ktcon1(i) = kbm(i)
      enddo
      do k = 2, km1
        do i = 1, im
          if (flg(i)) then
            if(k >= ktcon(i) .and. k < kbm(i)) then
              dz1 = zo(i,k+1) - zo(i,k)
              gamma = el2orc * qeso(i,k) / (to(i,k)**2)
              rfact =  1. + fv * cp * gamma
     &                 * to(i,k) / hvap
              aa1(i) = aa1(i) +
!    &                 dz1 * eta(i,k) * (grav / (cp * to(i,k)))
     &                 dz1 * (grav / (cp * to(i,k)))
     &                 * dbyo(i,k) / (1. + gamma)
     &                 * rfact
!             val = 0.
!             aa1(i) = aa1(i) +
!!   &                 dz1 * eta(i,k) * grav * fv *
!    &                 dz1 * grav * fv *
!    &                 max(val,(qeso(i,k) - qo(i,k)))
              if(aa1(i) < 0.) then
                ktcon1(i) = k
                flg(i) = .false.
              endif
            endif
          endif
        enddo
      enddo
c
c  compute cloud moisture property, detraining cloud water
c    and precipitation in overshooting layers
c
!> - For the overshooting convection, calculate the moisture content of the entraining/detraining parcel as before. Partition convective cloud water and precipitation and detrain convective cloud water in the overshooting layers.
      do k = 2, km1
        do i = 1, im
          if (cnvflg(i)) then
            if(k >= ktcon(i) .and. k < ktcon1(i)) then
              dz    = zi(i,k) - zi(i,k-1)
              gamma = el2orc * qeso(i,k) / (to(i,k)**2)
              qrch = qeso(i,k)
     &             + gamma * dbyo(i,k) / (hvap * (1. + gamma))
cj
              tem  = 0.25 * (xlamue(i,k)+xlamue(i,k-1)) * dz
              tem  = cq * tem
              factor = 1. + tem
              qcko(i,k) = ((1.-tem)*qcko(i,k-1)+tem*
     &                     (qo(i,k)+qo(i,k-1)))/factor
              qrcko(i,k) = qcko(i,k)
cj
              dq = eta(i,k) * (qcko(i,k) - qrch)
c
c  check if there is excess moisture to release latent heat
c
              if(dq > 0.) then
                etah = .5 * (eta(i,k) + eta(i,k-1))
                dp = 1000. * del(i,k)
                if(ncloud > 0) then
                  ptem = c0t(i,k) + c1
                  qlk = dq / (eta(i,k) + etah * ptem * dz)
                  dellal(i,k) = etah * c1 * dz * qlk * grav / dp
                else
                  qlk = dq / (eta(i,k) + etah * c0t(i,k) * dz)
                endif
                qcko(i,k) = qlk + qrch
                pwo(i,k) = etah * c0t(i,k) * dz * qlk
                cnvwt(i,k) = etah * qlk * grav / dp
                zdqca(i,k)=dq/eta(i,k)
              endif
            endif
          endif
        enddo
      enddo
!
!  compute updraft velocity square(wu2)
!> - Calculate updraft velocity square(wu2) according to Han et al.'s (2017) \cite han_et_al_2017 equation 7.
!
      bb1 = 4.0
      bb2 = 0.8
!
      if (hwrf_samfshal) then
      do i = 1, im
       if (cnvflg(i)) then
         k = kbcon1(i)
         tem = po(i,k) / (rd * to(i,k))
         wucb = -0.01 * dot(i,k) / (tem * grav)
         if(wucb > 0.) then
           wu2(i,k) = wucb * wucb
         else
           wu2(i,k) = 0.
         endif
       endif
      enddo
      endif
      do k = 2, km1
        do i = 1, im
          if (cnvflg(i)) then
            if(k > kbcon1(i) .and. k < ktcon(i)) then
              dz    = zi(i,k) - zi(i,k-1)
              tem  = 0.25 * bb1 * (drag(i,k)+drag(i,k-1)) * dz
              tem1 = 0.5 * bb2 * (buo(i,k)+buo(i,k-1)) * dz
              ptem = (1. - tem) * wu2(i,k-1)
              ptem1 = 1. + tem
              wu2(i,k) = (ptem + tem1) / ptem1
              wu2(i,k) = max(wu2(i,k), 0.)
            endif
          endif
        enddo
      enddo
!
      if(progsigma)then                                                                                                                               
          do k = 2, km1
            do i = 1, im
               if (cnvflg(i)) then
                  if(k > kbcon1(i) .and. k < ktcon(i)) then
                     rho = po(i,k)*100. / (rd * to(i,k))
                     omega_u(i,k)=-1.0*sqrt(wu2(i,k))*rho*grav
                     omega_u(i,k)=MAX(omega_u(i,k),-80.)
                  endif
               endif
            enddo
         enddo
      endif    

!  compute updraft velocity averaged over the whole cumulus
!
!> - Calculate the mean updraft velocity within the cloud (wc).
      do i = 1, im
        wc(i) = 0.
        sumx(i) = 0.
      enddo
      do k = 2, km1
        do i = 1, im
          if (cnvflg(i)) then
            if(k > kbcon1(i) .and. k < ktcon(i)) then
              dz = zi(i,k) - zi(i,k-1)
              tem = 0.5 * (sqrt(wu2(i,k)) + sqrt(wu2(i,k-1)))
              wc(i) = wc(i) + tem * dz
              sumx(i) = sumx(i) + dz
            endif
          endif
        enddo
      enddo
      do i = 1, im
        if(cnvflg(i)) then
          if(sumx(i) == 0.) then
             cnvflg(i)=.false.
          else
             wc(i) = wc(i) / sumx(i)
          endif
          val = 1.e-4
          if (wc(i) < val) cnvflg(i)=.false.
        endif
      enddo
c
!> - For progsigma =T, calculate the mean updraft velocity in pressure coordinates within the cloud (wc).                                                                                        
      if(progsigma)then                                                                                                                               
         do i = 1, im
            omegac(i) = 0.
            sumx(i) = 0.
         enddo
         do k = 2, km1
            do i = 1, im
               if (cnvflg(i)) then
                  if(k > kbcon1(i) .and. k < ktcon(i)) then
                     dp = 1000. * del(i,k)
                     tem = 0.5 * (omega_u(i,k) + omega_u(i,k-1))
                     omegac(i) = omegac(i) + tem * dp
                     sumx(i) = sumx(i) + dp
                  endif
               endif
            enddo
         enddo
         do i = 1, im
            if(cnvflg(i)) then
               if(sumx(i) == 0.) then
                  cnvflg(i)=.false.
               else
                  omegac(i) = omegac(i) / sumx(i)
               endif
               val = -1.2
               if (omegac(i) > val) cnvflg(i)=.false.
            endif
         enddo

!> - For progsigma = T, calculate the xi term in Bengtsson et al. 2022 \cite Bengtsson_2022 (equation 8)
         do k = 2, km1
            do i = 1, im
               if (cnvflg(i)) then
                  if(k > kbcon1(i) .and. k < ktcon(i)) then
                     if(omega_u(i,k) .ne. 0.)then
                        zeta(i,k)=eta(i,k)*(omegac(i)/omega_u(i,k))
                     else
                        zeta(i,k)=0.
                     endif
                     zeta(i,k)=MAX(0.,zeta(i,k))
                     zeta(i,k)=MIN(1.,zeta(i,k))
                  endif
               endif
            enddo
         enddo
      endif !if progsigma

c exchange ktcon with ktcon1
c
      do i = 1, im
        if(cnvflg(i)) then
          kk = ktcon(i)
          ktcon(i) = ktcon1(i)
          ktcon1(i) = kk
        endif
      enddo
c
c  this section is ready for cloud water
c
      if(ncloud > 0) then
c
c  compute liquid and vapor separation at cloud top
c
!> - => Separate the total updraft cloud water at cloud top into vapor and condensate.
      do i = 1, im
        if(cnvflg(i)) then
          k = ktcon(i) - 1
          gamma = el2orc * qeso(i,k) / (to(i,k)**2)
          qrch = qeso(i,k)
     &         + gamma * dbyo(i,k) / (hvap * (1. + gamma))
          dq = qcko(i,k) - qrch
c
c  check if there is excess moisture to release latent heat
c
          if(dq > 0.) then
            qlko_ktcon(i) = dq
            qcko(i,k) = qrch
            zdqca(i,k) = dq
          endif
        endif
      enddo
      endif
c
     
c--- compute precipitation efficiency in terms of windshear
c
!! - Calculate the wind shear and precipitation efficiency according to equation 58 in Fritsch and Chappell (1980) \cite fritsch_and_chappell_1980 :
!! \f[
!! E = 1.591 - 0.639\frac{\Delta V}{\Delta z} + 0.0953\left(\frac{\Delta V}{\Delta z}\right)^2 - 0.00496\left(\frac{\Delta V}{\Delta z}\right)^3
!! \f]
!! where \f$\Delta V\f$ is the integrated horizontal shear over the cloud depth, \f$\Delta z\f$, (the ratio is converted to units of \f$10^{-3} s^{-1}\f$). The variable "edt" is \f$1-E\f$ and is constrained to the range \f$[0,0.9]\f$.
!     do i = 1, im
!       if(cnvflg(i)) then
!         vshear(i) = 0.
!       endif
!     enddo
!     do k = 2, km
!       do i = 1, im
!         if (cnvflg(i)) then
!           if(k > kb(i) .and. k <= ktcon(i)) then
!             shear= sqrt((uo(i,k)-uo(i,k-1)) ** 2
!    &                  + (vo(i,k)-vo(i,k-1)) ** 2)
!             vshear(i) = vshear(i) + shear
!           endif
!         endif
!       enddo
!     enddo
!     do i = 1, im
!       if(cnvflg(i)) then
!         vshear(i) = 1.e3 * vshear(i) / (zi(i,ktcon(i))-zi(i,kb(i)))
!         e1=1.591-.639*vshear(i)
!    &       +.0953*(vshear(i)**2)-.00496*(vshear(i)**3)
!         edt(i)=1.-e1
!         val =         .9
!         edt(i) = min(edt(i),val)
!         val =         .0
!         edt(i) = max(edt(i),val)
!       endif
!     enddo
c
c--- what would the change be, that a cloud with unit mass
c--- will do to the environment?
c
!> ## Calculate the tendencies of the state variables (per unit cloud base mass flux) and the cloud base mass flux.
!> - Calculate the change in moist static energy, moisture mixing ratio, and horizontal winds per unit cloud base mass flux for all layers below cloud top from equations B.14 and B.15 from Grell (1993) \cite grell_1993, and for the cloud top from B.16 and B.17.
      do k = 1, km
        do i = 1, im
          if(cnvflg(i) .and. k <= kmax(i)) then
            dellah(i,k) = 0.
            dellaq(i,k) = 0.
            dellau(i,k) = 0.
            dellav(i,k) = 0.
          endif
        enddo
      enddo
      if (.not.hwrf_samfshal) then
       do n = 1, ntr
       do k = 1, km
        do i = 1, im
          if(cnvflg(i) .and. k <= kmax(i)) then
            dellae(i,k,n) = 0.
          endif
        enddo
       enddo
       enddo
      endif
c
c--- changed due to subsidence and entrainment
c
      do k = 2, km1
        do i = 1, im
          if (cnvflg(i)) then
            if(k > kb(i) .and. k < ktcon(i)) then
              dp = 1000. * del(i,k)
              dz = zi(i,k) - zi(i,k-1)
c
              dv1h = heo(i,k)
              dv2h = .5 * (heo(i,k) + heo(i,k-1))
              dv3h = heo(i,k-1)
c
              tem  = 0.5 * (xlamue(i,k)+xlamue(i,k-1))
              tem1 = xlamud(i)

              factor = grav / dp
cj
              dellah(i,k) = dellah(i,k) +
     &     ( eta(i,k)*dv1h - eta(i,k-1)*dv3h
     &    -  tem*eta(i,k-1)*dv2h*dz
     &    +  tem1*eta(i,k-1)*.5*(hcko(i,k)+hcko(i,k-1))*dz
     &         ) * factor
cj
              tem1 = -eta(i,k) * qrcko(i,k)
              tem2 = -eta(i,k-1) * qcko(i,k-1)
              dellaq(i,k) = dellaq(i,k) + (tem1-tem2) * factor
cj
              tem1=eta(i,k)*(uo(i,k)-ucko(i,k))
              tem2=eta(i,k-1)*(uo(i,k-1)-ucko(i,k-1))
              dellau(i,k) = dellau(i,k) + (tem1-tem2) * factor
cj
              tem1=eta(i,k)*(vo(i,k)-vcko(i,k))
              tem2=eta(i,k-1)*(vo(i,k-1)-vcko(i,k-1))
              dellav(i,k) = dellav(i,k) + (tem1-tem2) * factor
cj
            endif
          endif
        enddo
      enddo
      if(.not.hwrf_samfshal) then
       do n = 1, ntr
       do k = 2, km1
        do i = 1, im
          if (cnvflg(i)) then
            if(k > kb(i) .and. k < ktcon(i)) then
              dp = 1000. * del(i,k)
cj
              tem1 = -eta(i,k) * ercko(i,k,n)
              tem2 = -eta(i,k-1) * ecko(i,k-1,n)
              dellae(i,k,n) = dellae(i,k,n) + (tem1-tem2) * grav/dp
cj
            endif
          endif
        enddo
       enddo
       enddo
      endif
c
c------- cloud top
c
      do i = 1, im
        if(cnvflg(i)) then
          indx = ktcon(i)
          dp = 1000. * del(i,indx)
          tem = eta(i,indx-1) * grav / dp
          dellah(i,indx) = tem * (hcko(i,indx-1) - heo(i,indx-1))
          dellaq(i,indx) = tem *  qcko(i,indx-1)
          dellau(i,indx) = tem * (ucko(i,indx-1) - uo(i,indx-1))
          dellav(i,indx) = tem * (vcko(i,indx-1) - vo(i,indx-1))
c
c  cloud water
c
          dellal(i,indx) = tem * qlko_ktcon(i)
        endif
      enddo
      if (.not.hwrf_samfshal) then
      do n = 1, ntr
      do i = 1, im
        if(cnvflg(i)) then
          indx = ktcon(i)
          dp = 1000. * del(i,indx)
          dellae(i,indx,n) = eta(i,indx-1) *
     &             ecko(i,indx-1,n) * grav / dp
        endif
      enddo
      enddo
      endif
!
! compute change rates due to environmental subsidence & uplift
!     using a positive definite TVD flux-limiter scheme
!
!  for moisture
!
      do k=1,km1
        do i=1,im
          if(cnvflg(i) .and. k <= ktcon(i)) then
            q_diff(i,k) = q1(i,k) - q1(i,k+1)
          endif
        enddo
      enddo
      do i=1,im
        if(cnvflg(i)) then
          if(q1(i,1) >= 0.) then
            q_diff(i,0) = max(0.,2.*q1(i,1)-q1(i,2))-
     &                        q1(i,1)
          else
            q_diff(i,0) = min(0.,2.*q1(i,1)-q1(i,2))-
     &                          q1(i,1)
          endif
        endif
      enddo
!
      flxtvd = 0.
      do k = 1, km1
        do i = 1, im
          if(cnvflg(i) .and.
     &      (k >= kb(i) .and. k < ktcon(i))) then
            if(eta(i,k) > 0.) then
              rrkp = 0.
              if(abs(q_diff(i,k)) > 1.e-22)
     &               rrkp = q_diff(i,k+1) / q_diff(i,k)
              phkp = (rrkp+abs(rrkp)) / (1.+abs(rrkp))
              tem1 = q1(i,k+1) +
     &                   phkp*(qo(i,k)-q1(i,k+1))
              flxtvd(i,k) = eta(i,k) * tem1
            endif
          endif
        enddo
      enddo
!
      do k = 2, km1
        do i = 1, im
          if(cnvflg(i) .and.
     &      (k > kb(i) .and. k <= ktcon(i))) then
             dp = 1000. * del(i,k)
             dellaq(i,k) = dellaq(i,k) +
     &           (flxtvd(i,k) - flxtvd(i,k-1)) * grav/dp
          endif
        enddo
      enddo
!
!  for tracers including TKE & ozone
!
      if (.not.hwrf_samfshal) then
!
      do n=1,ntr
        do k=1,km1
          do i=1,im
            if(cnvflg(i) .and. k <= ktcon(i)) then
              e_diff(i,k,n) = ctr(i,k,n) - ctr(i,k+1,n)
            endif
          enddo
        enddo
        do i=1,im
          if(cnvflg(i)) then
            if(ctr(i,1,n) >= 0.) then
              e_diff(i,0,n) = max(0.,2.*ctr(i,1,n)-ctr(i,2,n))-
     &                          ctr(i,1,n)
            else
              e_diff(i,0,n) = min(0.,2.*ctr(i,1,n)-ctr(i,2,n))-
     &                          ctr(i,1,n)
            endif
          endif
        enddo
      enddo
!
      do n=1,ntr
!
        flxtvd = 0.
        do k= 1, km1
          do i = 1, im
            if(cnvflg(i) .and.
     &        (k >= kb(i) .and. k < ktcon(i))) then
              if(eta(i,k) > 0.) then
                rrkp = 0.
                if(abs(e_diff(i,k,n)) > 1.e-22)
     &                 rrkp = e_diff(i,k+1,n) / e_diff(i,k,n)
                phkp = (rrkp+abs(rrkp)) / (1.+abs(rrkp))
                tem1 = ctr(i,k+1,n) +
     &                     phkp*(ctro(i,k,n)-ctr(i,k+1,n))
                flxtvd(i,k) = eta(i,k) * tem1
              endif
            endif
          enddo
        enddo
!
        do k = 2, km1
          do i = 1, im
            if(cnvflg(i) .and.
     &        (k > kb(i) .and. k <= ktcon(i))) then
               dp = 1000. * del(i,k)
               dellae(i,k,n) = dellae(i,k,n) +
     &             (flxtvd(i,k) - flxtvd(i,k-1)) * grav/dp
            endif
          enddo
        enddo
!
      enddo
!
      endif
!
!  compute convective turn-over time
!
!> - Following Bechtold et al. (2008) \cite bechtold_et_al_2008, calculate the convective turnover time using the mean updraft velocity (wc) and the cloud depth. It is also proportional to the grid size (gdx).
      do i= 1, im
        if(cnvflg(i)) then
          tem = zi(i,ktcon1(i)) - zi(i,kbcon1(i))
          dtconv(i) = tem / wc(i)
          if (.not.hwrf_samfshal) then
            tfac = 1. + gdx(i) / 75000.
            dtconv(i) = tfac * dtconv(i)
          endif
          dtconv(i) = max(dtconv(i),dtmin)
          dtconv(i) = max(dtconv(i),dt2)
          dtconv(i) = min(dtconv(i),dtmax)
        endif
      enddo
!
!> - Calculate advective time scale (tauadv) using a mean cloud layer wind speed.
      do i= 1, im
        if(cnvflg(i)) then
          sumx(i) = 0.
          umean(i) = 0.
        endif
      enddo
      do k = 2, km1
        do i = 1, im
          if(cnvflg(i)) then
            if(k >= kbcon1(i) .and. k < ktcon1(i)) then
              dz = zi(i,k) - zi(i,k-1)
              tem = sqrt(u1(i,k)*u1(i,k)+v1(i,k)*v1(i,k))
              umean(i) = umean(i) + tem * dz
              sumx(i) = sumx(i) + dz
            endif
          endif
        enddo
      enddo
      do i= 1, im
        if(cnvflg(i)) then
           umean(i) = umean(i) / sumx(i)
           umean(i) = max(umean(i), 1.)
           tauadv = gdx(i) / umean(i)
           advfac(i) = tauadv / dtconv(i)
           advfac(i) = min(advfac(i), 1.)
        endif
      enddo
c
c  compute cloud base mass flux as a function of the mean
c     updraft velcoity
c
!> - From Bengtsson et al. (2022) \cite Bengtsson_2022 prognostic closure scheme, equation 8, call progsigma_calc() to compute updraft area fraction based on a moisture budget
      if(progsigma)then
!     Initial computations, dynamic q-tendency
         if(first_time_step .and. .not.restart)then
            do k = 1,km
               do i = 1,im
                  qadv(i,k)=0.
               enddo
            enddo
         else
            do k = 1,km
               do i = 1,im
                  qadv(i,k)=(q(i,k) - prevsq(i,k))*invdelt
               enddo
            enddo
         endif

         do k = 1,km
            do i = 1,im
               tmfq(i,k)=tmf(i,k,1)
            enddo
         enddo

         flag_shallow = .true.
         call progsigma_calc(im,km,first_time_step,restart,flag_shallow,
     &        del,tmfq,qmicro,dbyo1,zdqca,omega_u,zeta,hvap,delt,
     &        qadv,kbcon1,ktcon,cnvflg,
     &        sigmain,sigmaout,sigmab)
      endif

!> - From Han et al.'s (2017) \cite han_et_al_2017 equation 6, calculate cloud base mass flux as a function of the mean updraft velcoity.
!!  As discussed in Han et al. (2017) \cite han_et_al_2017 , when dtconv is larger than tauadv, the convective mixing is not fully conducted before the cumulus cloud is advected out of the grid cell. In this case, therefore, the cloud base mass flux is further reduced in proportion to the ratio of tauadv to dtconv.

      do i= 1, im
        if(cnvflg(i)) then
          k = kbcon(i)
          rho = po(i,k)*100. / (rd*to(i,k))
          if(progsigma .and. gdx(i) < dxcrtas)then
             xmb(i) = advfac(i)*sigmab(i)*((-1.0*omegac(i))*gravinv)
          else
             xmb(i) = advfac(i)*betaw*rho*wc(i)
          endif
        endif
      enddo
!
!> - For scale-aware parameterization, the updraft fraction (sigmagfm) is first computed as a function of the lateral entrainment rate at cloud base (see Han et al.'s (2017) \cite han_et_al_2017 equation 4 and 5), following the study by Grell and Freitas (2014) \cite grell_and_freitas_2014.
      do i = 1, im
        if(cnvflg(i)) then
          tem = min(max(xlamue(i,kbcon(i)), 2.e-4), 6.e-4)
          tem = 0.2 / tem
          tem1 = 3.14 * tem * tem
          sigmagfm(i) = tem1 / garea(i)
          sigmagfm(i) = max(sigmagfm(i), 0.001)
          sigmagfm(i) = min(sigmagfm(i), 0.999)
        endif
      enddo
!
!> - Then, calculate the reduction factor (scaldfunc) of the vertical convective eddy transport of mass flux as a function of updraft fraction from the studies by Arakawa and Wu (2013) \cite arakawa_and_wu_2013 (also see Han et al.'s (2017) \cite han_et_al_2017 equation 1 and 2). The final cloud base mass flux with scale-aware parameterization is obtained from the mass flux when sigmagfm << 1, multiplied by the reduction factor (Han et al.'s (2017) \cite han_et_al_2017 equation 2).
      do i = 1, im
        if(cnvflg(i)) then
          if (gdx(i) < dxcrt) then
             if(progsigma)then
                scaldfunc(i)=(1.-sigmab(i))*(1.-sigmab(i))
             else
                scaldfunc(i) = (1.-sigmagfm(i)) * (1.-sigmagfm(i))
             endif
             scaldfunc(i) = max(min(scaldfunc(i), 1.0), 0.)
          else
            scaldfunc(i) = 1.0
          endif
          xmb(i) = xmb(i) * scaldfunc(i)
          xmb(i) = min(xmb(i),xmbmax(i))
        endif
      enddo
!
!> - Transport aerosols if present
!
!     if (.not.hwrf_samfshal) then
!      if (do_aerosols)
!    &  call samfshalcnv_aerosols(im, im, km, itc, ntc, ntr, delt,
!!   &  xlamde, xlamdd, cnvflg, jmin, kb, kmax, kbcon, ktcon, fscav,
!    &  cnvflg, kb, kmax, ktcon, fscav,
!!   &  edto, xlamd, xmb, c0t, eta, etad, zi, xlamue, xlamud, delp,
!    &  xmb, c0t, eta, zi, xlamue, xlamud, delp,
!    &  qtr, qaero)
!     endif
!
!> ## For the "feedback control", calculate updated values of the state variables by multiplying the cloud base mass flux and the tendencies calculated per unit cloud base mass flux from the static control.
!! - Recalculate saturation specific humidity.
c
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
      do k = 1, km
        do i = 1, im
          if (cnvflg(i) .and. k <= kmax(i)) then
            qeso(i,k) = 0.01 * fpvs(t1(i,k))      ! fpvs is in pa
            qeso(i,k) = eps * qeso(i,k) / (pfld(i,k) + epsm1*qeso(i,k))
            val     =             1.e-8
            qeso(i,k) = max(qeso(i,k), val )
          endif
        enddo
      enddo
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
!> - Calculate the temperature tendency from the moist static energy and specific humidity tendencies.
!> - Update the temperature, specific humidity, and horiztonal wind state variables by multiplying the cloud base mass flux-normalized tendencies by the cloud base mass flux.
!> - Accumulate column-integrated tendencies.
      do i = 1, im
        delhbar(i) = 0.
        delqbar(i) = 0.
        deltbar(i) = 0.
        delubar(i) = 0.
        delvbar(i) = 0.
        qcond(i) = 0.
      enddo
      if (.not. hwrf_samfshal) then
       do n = 1, ntr
       do i = 1, im
        delebar(i,n) = 0.
       enddo
       enddo
      endif
      do k = 1, km
        do i = 1, im
          if (cnvflg(i)) then
            if(k > kb(i) .and. k <= ktcon(i)) then
              dellat = (dellah(i,k) - hvap * dellaq(i,k)) / cp
              t1(i,k) = t1(i,k) + dellat * xmb(i) * dt2
              q1(i,k) = q1(i,k) + dellaq(i,k) * xmb(i) * dt2
!             tem = 1./rcs(i)
!             u1(i,k) = u1(i,k) + dellau(i,k) * xmb(i) * dt2 * tem
!             v1(i,k) = v1(i,k) + dellav(i,k) * xmb(i) * dt2 * tem
              u1(i,k) = u1(i,k) + dellau(i,k) * xmb(i) * dt2
              v1(i,k) = v1(i,k) + dellav(i,k) * xmb(i) * dt2
              dp = 1000. * del(i,k)
              tem = xmb(i) * dp / grav
              delhbar(i) = delhbar(i) + tem * dellah(i,k)
              delqbar(i) = delqbar(i) + tem * dellaq(i,k)
              deltbar(i) = deltbar(i) + tem * dellat
              delubar(i) = delubar(i) + tem * dellau(i,k)
              delvbar(i) = delvbar(i) + tem * dellav(i,k)
            endif
          endif
        enddo
      enddo
!
! Negative moisture is set to zero after borrowing it from
!   positive values within the mass-flux transport layers
!
      do i = 1,im
        tsumn(i) = 0.
        tsump(i) = 0.
        rtnp(i) = 1.
      enddo
      do k = 1,km1
        do i = 1,im
          if (cnvflg(i)) then
            if(k > kb(i) .and. k <= ktcon(i)) then
              tem = q1(i,k) * delp(i,k) / grav
              if(q1(i,k) < 0.) tsumn(i) = tsumn(i) + tem
              if(q1(i,k) > 0.) tsump(i) = tsump(i) + tem
            endif
          endif
        enddo
      enddo
      do i = 1,im
        if(cnvflg(i)) then
          if(tsump(i) > 0. .and. tsumn(i) < 0.) then
            if(tsump(i) > abs(tsumn(i))) then
              rtnp(i) = tsumn(i) / tsump(i)
            else
              rtnp(i) = tsump(i) / tsumn(i)
            endif
          endif
        endif
      enddo
      do k = 1,km1
        do i = 1,im
          if (cnvflg(i)) then
            if(k > kb(i) .and. k <= ktcon(i)) then
              if(rtnp(i) < 0.) then
                if(tsump(i) > abs(tsumn(i))) then
                  if(q1(i,k) < 0.) q1(i,k)= 0.
                  if(q1(i,k) > 0.) q1(i,k)=(1.+rtnp(i))*q1(i,k)
                else
                  if(q1(i,k) < 0.) q1(i,k)=(1.+rtnp(i))*q1(i,k)
                  if(q1(i,k) > 0.) q1(i,k)=0.
                endif
              endif
            endif
          endif
        enddo
      enddo
!
      if (.not.hwrf_samfshal) then
!
      indx = ntk - 2
      do n = 1, ntr
!
       do k = 1, km
         do i = 1, im
           if (cnvflg(i)) then
             if(k > kb(i) .and. k <= ktcon(i)) then
               ctr(i,k,n) = ctr(i,k,n)+dellae(i,k,n)*xmb(i)*dt2
               dp = 1000. * del(i,k)
               delebar(i,n)=delebar(i,n)+dellae(i,k,n)*xmb(i)*dp/grav
             endif
           endif
         enddo
       enddo
!
! Negative TKE, ozone, and aerosols are set to zero after borrowing them
!     from positive values within the mass-flux transport layers
!
        do i = 1,im
          tsumn(i) = 0.
          tsump(i) = 0.
          rtnp(i) = 1.
        enddo
        do k = 1,km1
          do i = 1,im
            if (cnvflg(i)) then
              if(k > kb(i) .and. k <= ktcon(i)) then
                if(n == indx) then
                  if(k > 1) then
                    dz = zi(i,k) - zi(i,k-1)
                  else
                    dz = zi(i,k)
                  endif
                  tem = ctr(i,k,n) * dz
                else
                  tem = ctr(i,k,n) * delp(i,k) / grav
                endif
                if(ctr(i,k,n) < 0.) tsumn(i) = tsumn(i) + tem
                if(ctr(i,k,n) > 0.) tsump(i) = tsump(i) + tem
              endif
            endif
          enddo
        enddo
        do i = 1,im
          if(cnvflg(i)) then
            if(tsump(i) > 0. .and. tsumn(i) < 0.) then
              if(tsump(i) > abs(tsumn(i))) then
                rtnp(i) = tsumn(i) / tsump(i)
              else
                rtnp(i) = tsump(i) / tsumn(i)
              endif
            endif
          endif
        enddo
        do k = 1,km1
        do i = 1,im
          if (cnvflg(i)) then
            if(k > kb(i) .and. k <= ktcon(i)) then
              if(rtnp(i) < 0.) then
                if(tsump(i) > abs(tsumn(i))) then
                  if(ctr(i,k,n)<0.) ctr(i,k,n)=0.
                  if(ctr(i,k,n)>0.) ctr(i,k,n)=(1.+rtnp(i))*ctr(i,k,n)
                else
                  if(ctr(i,k,n)<0.) ctr(i,k,n)=(1.+rtnp(i))*ctr(i,k,n)
                  if(ctr(i,k,n)>0.) ctr(i,k,n)=0.
                endif
              endif
            endif
          endif
        enddo
        enddo
!
        kk = n+2
        do k = 1, km
        do i = 1, im
          if (cnvflg(i)) then
            if(k > kb(i) .and. k <= ktcon(i)) then
              qtr(i,k,kk) = ctr(i,k,n)
            endif
          endif
        enddo
        enddo
!
      enddo
!
       if (do_aerosols) then
!
        do n = 1, ntc
!
!  convert wet deposition to total mass deposited over dt2 and dp
          do k = 2, km1
            do i = 1, im
              if (cnvflg(i)) then
                if(k > kb(i) .and. k < ktcon(i)) then
                  dp = 1000. * del(i,k)
                  wet_dep(i,k,n) = chem_pw(i,k,n)*grav/dp
                  wet_dep(i,k,n) = wet_dep(i,k,n)*xmb(i)*dt2*dp
                endif
              endif
            enddo
          enddo
!
          kk = n + itc - 1
          do k = 2, km1
            do i = 1, im
              if (cnvflg(i)) then
                if(k > kb(i) .and. k < ktcon(i)) then
                  dp = 1000. * del(i,k)
                  if (qtr(i,k,kk) < 0.) then
!   borrow negative mass from wet deposition
                    tem = -qtr(i,k,kk)*dp
                    if(wet_dep(i,k,n) >= tem) then
                      wet_dep(i,k,n) = wet_dep(i,k,n) - tem
                      qtr(i,k,kk) = 0.
                    else
                      wet_dep(i,k,n) = 0.
                      qtr(i,k,kk) = qtr(i,k,kk)+wet_dep(i,k,n)/dp
                    endif
                  endif
                endif
              endif
            enddo
          enddo
!
        enddo
!
       endif
!
      endif
!
!> - Recalculate saturation specific humidity using the updated temperature.
      do k = 1, km
        do i = 1, im
          if (cnvflg(i)) then
            if(k > kb(i) .and. k <= ktcon(i)) then
              qeso(i,k) = 0.01 * fpvs(t1(i,k))      ! fpvs is in pa
              qeso(i,k) = eps * qeso(i,k)/(pfld(i,k) + epsm1*qeso(i,k))
              val     =             1.e-8
              qeso(i,k) = max(qeso(i,k), val )
            endif
          endif
        enddo
      enddo
c
!> - Add up column-integrated convective precipitation by multiplying the normalized value by the cloud base mass flux.
      do i = 1, im
        rntot(i) = 0.
        delqev(i) = 0.
        delq2(i) = 0.
        flg(i) = cnvflg(i)
      enddo
      do k = km, 1, -1
        do i = 1, im
          if (cnvflg(i)) then
            if(k < ktcon(i) .and. k > kb(i)) then
              rntot(i) = rntot(i) + pwo(i,k) * xmb(i) * .001 * dt2
            endif
          endif
        enddo
      enddo
c
c evaporating rain
c
!> - Determine the evaporation of the convective precipitation and update the integrated convective precipitation.
!> - Update state temperature and moisture to account for evaporation of convective precipitation.
!> - Update column-integrated tendencies to account for evaporation of convective precipitation.
      do k = km, 1, -1
        do i = 1, im
          if (k <= kmax(i)) then
            deltv(i) = 0.
            delq(i) = 0.
            qevap(i) = 0.
            if(cnvflg(i)) then
              if(k < ktcon(i) .and. k > kb(i)) then
                rn(i) = rn(i) + pwo(i,k) * xmb(i) * .001 * dt2
              endif
            endif
            if(flg(i) .and. k < ktcon(i)) then
!             evef = edt(i) * evfact
!             if(islimsk(i) == 1) evef=edt(i) * evfactl
!             if(islimsk(i) == 1) evef=.07
              qcond(i) = shevf * evef * (q1(i,k) - qeso(i,k))
     &                 / (1. + el2orc * qeso(i,k) / t1(i,k)**2)
              dp = 1000. * del(i,k)
              factor = dp / grav
              if(rn(i) > 0. .and. qcond(i) < 0.) then
                qevap(i) = -qcond(i) * (1.-exp(-.32*sqrt(dt2*rn(i))))
                qevap(i) = min(qevap(i), rn(i)*1000.*grav/dp)
                delq2(i) = delqev(i) + .001 * qevap(i) * factor
              endif
              if(rn(i) > 0. .and. qcond(i) < 0. .and.
     &           delq2(i) > rntot(i)) then
                qevap(i) = 1000.* grav * (rntot(i) - delqev(i)) / dp
                flg(i) = .false.
              endif
              if(rn(i) > 0. .and. qevap(i) > 0.) then
                tem  = .001 * factor
                tem1 = qevap(i) * tem
                if(tem1 > rn(i)) then
                  qevap(i) = rn(i) / tem
                  rn(i) = 0.
                else
                  rn(i) = rn(i) - tem1
                endif
                q1(i,k) = q1(i,k) + qevap(i)
                t1(i,k) = t1(i,k) - elocp * qevap(i)
                deltv(i) = - elocp*qevap(i)/dt2
                delq(i) =  + qevap(i)/dt2
                delqev(i) = delqev(i) + tem * qevap(i)
              endif
              delqbar(i) = delqbar(i) + delq(i)  * factor
              deltbar(i) = deltbar(i) + deltv(i) * factor
            endif
          endif
        enddo
      enddo
cj
!     do i = 1, im
!     if(me == 31 .and. cnvflg(i)) then
!     if(cnvflg(i)) then
!       print *, ' shallow delhbar, delqbar, deltbar = ',
!    &             delhbar(i),hvap*delqbar(i),cp*deltbar(i)
!       print *, ' shallow delubar, delvbar = ',delubar(i),delvbar(i)
!       print *, ' precip =', hvap*rn(i)*1000./dt2
!       print*,'pdif= ',pfld(i,kbcon(i))-pfld(i,ktcon(i))
!     endif
!     enddo
!     do n = 1, ntr
!     do i = 1, im
!     if(me == 31 .and. cnvflg(i)) then
!     if(cnvflg(i)) then
!       print *, ' tracer delebar = ',delebar(i,n)
!     endif
!     enddo
!     enddo
cj
      do i = 1, im
        if(cnvflg(i)) then
          if(rn(i) < 0. .or. .not.flg(i)) rn(i) = 0.
          ktop(i) = ktcon(i)
          kbot(i) = kbcon(i)
          kcnv(i) = 2
        endif
      enddo
c
c  convective cloud water
c
!> - Calculate shallow convective cloud water.
      do k = 1, km
        do i = 1, im
          if (cnvflg(i)) then
            if (k >= kbcon(i) .and. k < ktcon(i)) then
              cnvw(i,k) = cnvwt(i,k) * xmb(i) * dt2
            endif
          endif
        enddo
      enddo

c
c  convective cloud cover
c
!> - Calculate convective cloud cover, which is used when pdf-based cloud fraction is used (i.e., pdfcld=.true.).
      do k = 1, km
        do i = 1, im
          if (cnvflg(i)) then
            if (k >= kbcon(i) .and. k < ktcon(i)) then
              cnvc(i,k) = 0.04 * log(1. + 675. * eta(i,k) * xmb(i))
              cnvc(i,k) = min(cnvc(i,k), 0.2)
              cnvc(i,k) = max(cnvc(i,k), 0.0)
            endif
          endif
        enddo
      enddo
c
c  cloud water
c
!> - Separate detrained cloud water into liquid and ice species as a function of temperature only.
      if (ncloud > 0) then
!
      do k = 1, km1
        do i = 1, im
          if (cnvflg(i)) then
!           if (k > kb(i) .and. k <= ktcon(i)) then
            if (k >= kbcon(i) .and. k <= ktcon(i)) then
              tem  = dellal(i,k) * xmb(i) * dt2
              tem1 = max(0.0, min(1.0, (tcr-t1(i,k))*tcrf))
              if (qtr(i,k,2) > -999.0) then
                qtr(i,k,1) = qtr(i,k,1) + tem * tem1            ! ice
                qtr(i,k,2) = qtr(i,k,2) + tem *(1.0-tem1)       ! water
              else
                qtr(i,k,1) = qtr(i,k,1) + tem
              endif
            endif
          endif
        enddo
      enddo
!
      endif
!> - Store aerosol concentrations if present
!     if (.not. hwrf_samfshal) then
!      if (do_aerosols) then
!       do n = 1, ntc
!         kk = n + itc - 1
!         do k = 1, km
!           do i = 1, im
!             if(cnvflg(i) .and. rn(i) > 0.) then
!               if (k <= kmax(i)) qtr(i,k,kk) = qaero(i,k,n)
!             endif
!           enddo
!         enddo
!       enddo
!      endif
!     endif
!
! hchuang code change
!
!> - Calculate and retain the updraft mass flux for dust transport by cumulus convection.
!
!> - Calculate the updraft convective mass flux.
      do k = 1, km
        do i = 1, im
          if(cnvflg(i)) then
            if(k >= kb(i) .and. k < ktop(i)) then
              ud_mf(i,k) = eta(i,k) * xmb(i) * dt2
            endif
          endif
        enddo
      enddo
!> - save the updraft convective mass flux at cloud top.
      do i = 1, im
        if(cnvflg(i)) then
           k = ktop(i)-1
           dt_mf(i,k) = ud_mf(i,k)
        endif
      enddo
!
!   include TKE contribution from shallow convection
!
      if (.not.hwrf_samfshal) then
      if (ntk > 0) then
!
      do k = 2, km1
        do i = 1, im
          if(cnvflg(i)) then
            if(k > kb(i) .and. k < ktop(i)) then
              tem = 0.5 * (eta(i,k-1) + eta(i,k)) * xmb(i)
              tem1 = pfld(i,k) * 100. / (rd * t1(i,k))
              if(progsigma)then
                tem2 = sigmab(i)
              else
                tem2 = max(sigmagfm(i), betaw)
              endif
              ptem = tem / (tem2 * tem1)
              qtr(i,k,ntk)=qtr(i,k,ntk)+0.5*tem2*ptem*ptem
            endif
          endif
        enddo
      enddo
!
      endif
      endif
!!
      return
      end subroutine samfshalcnv_run
!> @}
      end module samfshalcnv


      module samfcnv_aerosols

      implicit none

      private

      public :: samfdeepcnv_aerosols, samfshalcnv_aerosols

      contains

      subroutine samfdeepcnv_aerosols(im, ix, km, itc, ntc, ntr, delt,
     &  xlamde, xlamdd, cnvflg, jmin, kb, kmax, kd94, ktcon, fscav,
     &  edto, xlamd, xmb, c0t, eta, etad, zi, xlamue, xlamud, delp,
     &  qtr, qaero)

      use machine , only : kind_phys
      use physcons, only : g => con_g, qamin

      implicit none

c     -- input arguments
      integer,                                      intent(in) :: im,
     &  ix, km, itc, ntc, ntr
      real(kind=kind_phys),                         intent(in) :: delt,
     &  xlamde, xlamdd
      logical,              dimension(im),          intent(in) :: cnvflg
      integer,              dimension(im),          intent(in) :: jmin,
     &  kb, kmax, kd94, ktcon
      real(kind=kind_phys), dimension(im),          intent(in) :: edto,
     &  xlamd, xmb
      real(kind=kind_phys), dimension(ntc),         intent(in) :: fscav
      real(kind=kind_phys), dimension(im,km),       intent(in) :: c0t,
     &  eta, etad, zi, xlamue, xlamud
      real(kind=kind_phys), dimension(ix,km),       intent(in) :: delp
      real(kind=kind_phys), dimension(ix,km,ntr+2), intent(in) :: qtr
c     -- output arguments
      real(kind=kind_phys), dimension(im,km,ntc),   intent(out) :: qaero

c     -- local variables
c     -- general variables
      integer :: i, indx, it, k, kk, km1, kp1, n
      real(kind=kind_phys) :: adw, aup, dtime_max, dv1q, dv2q, dv3q,
     &  dtovdz, dz, factor, ptem, ptem1, qamax, tem, tem1
      real(kind=kind_phys), dimension(ix,km) :: xmbp
c     -- chemical transport variables
      real(kind=kind_phys), dimension(im,km,ntc) :: ctro2, ecko2, ecdo2,
     &  dellae2
c     -- additional variables for tracers for wet deposition,
      real(kind=kind_phys), dimension(im,km,ntc) :: chem_c, chem_pw,
     &  wet_dep
c     -- if reevaporation is enabled, uncomment lines below
c     real(kind=kind_phys), dimension(im,ntc) :: pwav
c     real(kind=kind_phys), dimension(im,km) :: pwdper
c     real(kind=kind_phys), dimension(im,km,ntr) :: chem_pwd
c     -- additional variables for fct
      real(kind=kind_phys), dimension(im,km) :: flx_lo, totlout, clipout

      real(kind=kind_phys), parameter :: one     = 1.0_kind_phys
      real(kind=kind_phys), parameter :: half    = 0.5_kind_phys
      real(kind=kind_phys), parameter :: quarter = 0.25_kind_phys
      real(kind=kind_phys), parameter :: zero    = 0.0_kind_phys
      real(kind=kind_phys), parameter :: epsil   = 1.e-22_kind_phys    ! prevent division by zero

c     -- begin

c     -- check if aerosols are present
      if ( ntc <= 0 .or. itc <= 0 .or. ntr <= 0 ) return
      if ( ntr < itc + ntc - 3 ) return

c     -- initialize work variables
      km1 = km - 1

      chem_c    = zero
      chem_pw   = zero
      ctro2     = zero
      dellae2   = zero
      ecdo2     = zero
      ecko2     = zero
      qaero     = zero
      wet_dep   = zero

c     -- set work arrays

      do n = 1, ntc
        it = n + itc - 1
        do k = 1, km
          do i = 1, im
            if (k <= kmax(i)) qaero(i,k,n) = max(qamin, qtr(i,k,it))
          enddo
        enddo
      enddo

      do k = 1, km
        do i = 1, im
          if (cnvflg(i)) xmbp(i,k) = g * xmb(i) / delp(i,k)
        enddo
      enddo

      do n = 1, ntc
c       -- interface level
        do k = 1, km1
          kp1 = k + 1
          do i = 1, im
            if (kp1 <= kmax(i)) ctro2(i,k,n) =
     &        half * (qaero(i,k,n) + qaero(i,kp1,n))
          enddo
        enddo
c       -- top level
        do i = 1, im
          ctro2(i,kmax(i),n) = qaero(i,kmax(i),n)
        enddo
      enddo

      do n = 1, ntc
        do k = 1, km
          do i = 1, im
            if (cnvflg(i) .and. (k <= kb(i)))
     &        ecko2(i,k,n) = ctro2(i,k,n)
          enddo
        enddo
      enddo

      do n = 1, ntc
        do i = 1, im
          if (cnvflg(i)) ecdo2(i,jmin(i),n) = ctro2(i,jmin(i),n)
        enddo
      enddo

c     do chemical tracers, first need to know how much reevaporates

c     aerosol re-evaporation is set to zero for now
c     uncomment and edit the following code to enable re-evaporation
c     chem_pwd  = zero
c     pwdper    = zero
c     pwav      = zero
c     do i = 1, im
c        do k=1,jmin(i)
c           pwdper(i,k)= -edto(i)*pwdo(i,k)/pwavo(i)
c        enddo
c     enddo
c
c     calculate include mixing ratio (ecko2), how much goes into
c     rainwater to be rained out (chem_pw), and total scavenged,
c     if not reevaporated (pwav)

      do n = 1, ntc
        do k = 2, km1
          kk = k - 1
          do i = 1, im
            if (cnvflg(i)) then
              if ((k > kb(i)) .and. (k < ktcon(i))) then
                dz   = zi(i,k) - zi(i,kk)
                tem  = half    * (xlamue(i,k)+xlamue(i,kk)) * dz
                tem1 = quarter * (xlamud(i,k)+xlamud(i,kk)) * dz
                factor = one + tem - tem1

c               if conserved (not scavenging) then
                ecko2(i,k,n) = ((one-tem1)*ecko2(i,kk,n)
     &            + half*tem*(ctro2(i,k,n)+ctro2(i,kk,n)))/factor

c               how much will be scavenged
c
c               this choice was used in GF, and is also described in a
c               successful implementation into CESM in GRL (Yu et al. 2019),
c               it uses dimesnsionless scavenging coefficients (fscav),
c               but includes henry coeffs with gas phase chemistry

c               fraction fscav is going into liquid
                chem_c(i,k,n)=fscav(n)*ecko2(i,k,n)

c               of that part is going into rain out (chem_pw)
                tem=chem_c(i,k,n)/(one+c0t(i,k)*dz)
                chem_pw(i,k,n)=c0t(i,k)*dz*tem*eta(i,kk) !etah
                ecko2(i,k,n)=tem+ecko2(i,k,n)-chem_c(i,k,n)

c               pwav needed fo reevaporation in downdraft
c               if including reevaporation, please uncomment code below
c               pwav(i,n)=pwav(i,n)+chem_pw(i,k,n)
              endif
            endif
          enddo
        enddo
        do k = 1, km1
          do i = 1, im
            if (k >= ktcon(i)) ecko2(i,k,n)=ctro2(i,k,n)
          enddo
        enddo
      enddo

c     reevaporation of some, pw and pwd terms needed later for dellae2

      do n = 1, ntc
        do k = km1, 1, -1
          kp1 = k + 1
          do i = 1, im
            if (cnvflg(i) .and. (k < jmin(i))) then
              dz = zi(i,kp1) - zi(i,k)
              if (k >= kd94(i)) then
                tem  = xlamde * dz
                tem1 = half * xlamdd * dz
              else
                tem  = xlamde * dz
                tem1 = half * (xlamd(i)+xlamdd) * dz
              endif
              factor = one + tem - tem1
              ecdo2(i,k,n) = ((one-tem1)*ecdo2(i,kp1,n)
     &          +half*tem*(ctro2(i,k,n)+ctro2(i,kp1,n)))/factor
c             if including reevaporation, please uncomment code below
c             ecdo2(i,k,n)=ecdo2(i,k,n)+pwdper(i,kp1)*pwav(i,n)
c             chem_pwd(i,k,n)=max(zero,pwdper(i,kp1)*pwav(i,n))
            endif
          enddo
        enddo
      enddo

      do n = 1, ntc
        do i = 1, im
          if (cnvflg(i)) then
c           subsidence term treated in fct routine
            dellae2(i,1,n) = edto(i)*etad(i,1)*ecdo2(i,1,n)*xmbp(i,1)
          endif
        enddo
      enddo

      do n = 1, ntc
        do i = 1, im
          if (cnvflg(i)) then
            k = ktcon(i)
            kk = k - 1
c           for the subsidence term already is considered
            dellae2(i,k,n) = eta(i,kk) * ecko2(i,kk,n) * xmbp(i,k)
          endif
        enddo
      enddo

c     --- for updraft & downdraft vertical transport
c
c     initialize maximum allowed timestep for upstream difference approach
c
      dtime_max=delt
      do k=2,km1
        kk = k - 1
        do i = 1, im
          if (kk < ktcon(i)) dtime_max = min(dtime_max,half*delp(i,kk))
        enddo
      enddo

c     now for every chemistry tracer
      do n = 1, ntc
        do k = 2, km1
          kk = k - 1
          do i = 1, im
            if (cnvflg(i) .and. (k < ktcon(i))) then
              dz = zi(i,k) - zi(i,kk)
              aup = one
              if (k <= kb(i)) aup = zero
              adw = one
              if (k > jmin(i)) adw = zero

              dv1q = half * (ecko2(i,k,n) + ecko2(i,kk,n))
              dv2q = half * (ctro2(i,k,n) + ctro2(i,kk,n))
              dv3q = half * (ecdo2(i,k,n) + ecdo2(i,kk,n))

              tem  = half * (xlamue(i,k) + xlamue(i,kk))
              tem1 = half * (xlamud(i,k) + xlamud(i,kk))

              if (k <= kd94(i)) then
                ptem  = xlamde
                ptem1 = xlamd(i) + xlamdd
              else
                ptem  = xlamde
                ptem1 = xlamdd
              endif
              dellae2(i,k,n) = dellae2(i,k,n) +
c                 detrainment from updraft
     &          ( aup*tem1*eta(i,kk)*dv1q
c                 entrainement into up and downdraft
     &          - (aup*tem*eta(i,kk)+adw*edto(i)*ptem*etad(i,k))*dv2q
c                 detrainment from downdraft
     &          + (adw*edto(i)*ptem1*etad(i,k)*dv3q) ) * dz * xmbp(i,k)

              wet_dep(i,k,n)=chem_pw(i,k,n)*g/delp(i,k)

c             sinks from where updraft and downdraft start
              if (k == jmin(i)+1) then
                dellae2(i,k,n) = dellae2(i,k,n)
     &            -edto(i)*etad(i,kk)*ctro2(i,kk,n)*xmbp(i,k)
              endif
              if (k == kb(i))then
                dellae2(i,k,n) = dellae2(i,k,n)
     &            -eta(i,k)*ctro2(i,k,n)*xmbp(i,k)
              endif
            endif
          enddo
        enddo

        do i = 1, im
          if (cnvflg(i)) then
            if (kb(i) == 1) then
              k=kb(i)
              dellae2(i,k,n) = dellae2(i,k,n)
     &          -eta(i,k)*ctro2(i,k,n)*xmbp(i,k)
            endif
          endif
        enddo

      enddo

c     for every tracer...

      do n = 1, ntc
        flx_lo  = zero
        totlout = zero
        clipout = zero
c       compute low-order mass flux, upstream
        do k = 2, km1
          kk = k - 1
          do i = 1, im
            if (cnvflg(i) .and. (kk < ktcon(i))) then
              tem = zero
              if (kk >= kb(i)  ) tem = eta(i,kk)
              if (kk <= jmin(i)) tem = tem - edto(i)*etad(i,kk)
c             low-order flux,upstream
              if (tem > zero) then
                flx_lo(i,k) = -xmb(i) * tem * qaero(i,k,n)
              elseif (tem < zero) then
                flx_lo(i,k) = -xmb(i) * tem * qaero(i,kk,n)
              endif
            endif
          enddo
        enddo

c       --- make sure low-ord fluxes don't violate positive-definiteness
        do k=1,km1
          kp1 = k + 1
          do i=1,im
            if (cnvflg(i) .and. (k <= ktcon(i))) then
c             time step / grid spacing
              dtovdz = g * dtime_max / abs(delp(i,k))
c             total flux out
              totlout(i,k)=max(zero,flx_lo(i,kp1))-min(zero,flx_lo(i,k))
              clipout(i,k)=min(one ,qaero(i,k,n)/max(epsil,totlout(i,k))
     &           / (1.0001_kind_phys*dtovdz))
            endif
          enddo
        enddo

c       recompute upstream mass fluxes
        do k = 2, km1
          kk = k - 1
          do i = 1, im
            if (cnvflg(i) .and. (kk < ktcon(i))) then
              tem = zero
              if (kk >= kb(i)  ) tem = eta(i,kk)
              if (kk <= jmin(i)) tem = tem - edto(i)*etad(i,kk)
              if (tem > zero) then
                flx_lo(i,k) = flx_lo(i,k) * clipout(i,k)
              elseif (tem < zero) then
                flx_lo(i,k) = flx_lo(i,k) * clipout(i,kk)
              endif
            endif
          enddo
        enddo

c      --- a positive-definite low-order (diffusive) solution for the subsidnce fluxes
       do k=1,km1
         kp1 = k + 1
         do i=1,im
           if (cnvflg(i) .and. (k <= ktcon(i))) then
             dtovdz = g * dtime_max / abs(delp(i,k)) ! time step /grid spacing
             dellae2(i,k,n) = dellae2(i,k,n)
     &         -(flx_lo(i,kp1)-flx_lo(i,k))*dtovdz/dtime_max
           endif
          enddo
        enddo

      enddo ! ctr

c     convert wet deposition to total mass deposited over dt and dp

      do n = 1, ntc
        do k = 1, km
          do i = 1, im
            if (cnvflg(i) .and. (k < ktcon(i)))
     &        wet_dep(i,k,n) = wet_dep(i,k,n)*xmb(i)*delt*delp(i,k)
          enddo
        enddo
      enddo

c     compute final aerosol concentrations

      do n = 1, ntc
        do k = 1, km
          do i = 1, im
            if (cnvflg(i) .and. (k <= min(kmax(i),ktcon(i)))) then
              qaero(i,k,n) = qaero(i,k,n) + dellae2(i,k,n) * delt
              if (qaero(i,k,n) < zero) then
c               add negative mass to wet deposition
                wet_dep(i,k,n) = wet_dep(i,k,n)-qaero(i,k,n)*delp(i,k)
                qaero(i,k,n) = qamin
              endif
            endif
          enddo
        enddo
      enddo

      return
      end subroutine samfdeepcnv_aerosols

      subroutine samfshalcnv_aerosols(im, ix, km, itc, ntc, ntr, delt,
     &  cnvflg, kb, kmax, kbcon, ktcon, fscav,
     &  xmb, c0t, eta, zi, xlamue, xlamud, delp,
     &  qtr, qaero)

      use machine , only : kind_phys
      use physcons, only : g => con_g, qamin

      implicit none

c     -- input arguments
      integer,                                      intent(in) :: im,
     &  ix, km, itc, ntc, ntr
      real(kind=kind_phys),                         intent(in) :: delt
!     &  xlamde, xlamdd
      logical,              dimension(im),          intent(in) :: cnvflg
!     integer,              dimension(im),          intent(in) :: jmin,
      integer,              dimension(im),          intent(in) ::
     &  kb, kmax, kbcon, ktcon
      real(kind=kind_phys), dimension(im),          intent(in) ::
     &  xmb, xlamud
      real(kind=kind_phys), dimension(ntc),         intent(in) :: fscav
      real(kind=kind_phys), dimension(im,km),       intent(in) :: c0t,
     &  eta, zi, xlamue  !, xlamud
      real(kind=kind_phys), dimension(ix,km),       intent(in) :: delp
      real(kind=kind_phys), dimension(ix,km,ntr+2), intent(in) :: qtr
c     -- output arguments
      real(kind=kind_phys), dimension(im,km,ntc),   intent(out) :: qaero

c     -- local variables
c     -- general variables
      integer :: i, indx, it, k, kk, km1, kp1, n
!     real(kind=kind_phys) :: adw, aup, dtime_max, dv1q, dv2q, dv3q,
      real(kind=kind_phys) :: aup, dtime_max, dv1q, dv2q, dv3q,
     &  dtovdz, dz, factor, ptem, ptem1, qamax, tem, tem1
      real(kind=kind_phys), dimension(ix,km) :: xmbp
c     -- chemical transport variables
      real(kind=kind_phys), dimension(im,km,ntc) :: ctro2,ecko2,dellae2
c     -- additional variables for tracers for wet deposition,
      real(kind=kind_phys), dimension(im,km,ntc) :: chem_c, chem_pw,
     &  wet_dep
c     -- if reevaporation is enabled, uncomment lines below
c     real(kind=kind_phys), dimension(im,ntc) :: pwav
c     real(kind=kind_phys), dimension(im,km) :: pwdper
c     real(kind=kind_phys), dimension(im,km,ntr) :: chem_pwd
c     -- additional variables for fct
      real(kind=kind_phys), dimension(im,km) :: flx_lo, totlout, clipout

      real(kind=kind_phys), parameter :: one     = 1.0_kind_phys
      real(kind=kind_phys), parameter :: half    = 0.5_kind_phys
      real(kind=kind_phys), parameter :: quarter = 0.25_kind_phys
      real(kind=kind_phys), parameter :: zero    = 0.0_kind_phys
      real(kind=kind_phys), parameter :: epsil   = 1.e-22_kind_phys    ! prevent division by zero
      real(kind=kind_phys), parameter :: escav   = 0.8_kind_phys       ! wet scavenging efficiency

c     -- begin

c     -- check if aerosols are present
      if ( ntc <= 0 .or. itc <= 0 .or. ntr <= 0 ) return
      if ( ntr  < itc + ntc - 3 ) return

c     -- initialize work variables
      km1 = km - 1

      chem_c    = zero
      chem_pw   = zero
      ctro2     = zero
      dellae2   = zero
      !ecdo2     = zero
      ecko2     = zero
      qaero     = zero
      wet_dep   = zero

c     -- set work arrays

      do n = 1, ntc
        it = n + itc - 1
        do k = 1, km
          do i = 1, im
            if (k <= kmax(i)) qaero(i,k,n) = max(qamin, qtr(i,k,it))
          enddo
        enddo
      enddo

      do k = 1, km
        do i = 1, im
          if (cnvflg(i)) xmbp(i,k) = g * xmb(i) / delp(i,k)
        enddo
      enddo

      do n = 1, ntc
c       -- interface level
        do k = 1, km1
          kp1 = k + 1
          do i = 1, im
            if (kp1 <= kmax(i)) ctro2(i,k,n) =
     &        half * (qaero(i,k,n) + qaero(i,kp1,n))
          enddo
        enddo
c       -- top level
        do i = 1, im
          ctro2(i,kmax(i),n) = qaero(i,kmax(i),n)
        enddo
      enddo

      do n = 1, ntc
        do k = 1, km
          do i = 1, im
            if (cnvflg(i) .and. (k <= kb(i)))
     &        ecko2(i,k,n) = ctro2(i,k,n)
          enddo
        enddo
      enddo

      !do n = 1, ntc
      !  do i = 1, im
      !    if (cnvflg(i)) ecdo2(i,jmin(i),n) = ctro2(i,jmin(i),n)
      !  enddo
      !enddo

c     do chemical tracers, first need to know how much reevaporates

c     aerosol re-evaporation is set to zero for now
c     uncomment and edit the following code to enable re-evaporation
c     chem_pwd  = zero
c     pwdper    = zero
c     pwav      = zero
c     do i = 1, im
c        do k=1,jmin(i)
c           pwdper(i,k)= -edto(i)*pwdo(i,k)/pwavo(i)
c        enddo
c     enddo
c
c     calculate include mixing ratio (ecko2), how much goes into
c     rainwater to be rained out (chem_pw), and total scavenged,
c     if not reevaporated (pwav)

      do n = 1, ntc
        do k = 2, km1
          kk = k - 1
          do i = 1, im
            if (cnvflg(i)) then
              if ((k > kb(i)) .and. (k < ktcon(i))) then
                dz   = zi(i,k) - zi(i,kk)
                tem  = half    * (xlamue(i,k)+xlamue(i,kk)) * dz
!               tem1 = quarter * (xlamud(i,k)+xlamud(i,kk)) * dz
                tem1 = quarter * (xlamud(i  )+xlamud(i   )) * dz
                factor = one + tem - tem1

c               if conserved (not scavenging) then
                ecko2(i,k,n) = ((one-tem1)*ecko2(i,kk,n)
     &            + half*tem*(ctro2(i,k,n)+ctro2(i,kk,n)))/factor

c               how much will be scavenged
c
c               this choice was used in GF, and is also described in a
c               successful implementation into CESM in GRL (Yu et al. 2019),
c               it uses dimesnsionless scavenging coefficients (fscav),
c               but includes henry coeffs with gas phase chemistry

c               fraction fscav is going into liquid
                chem_c(i,k,n)=escav*fscav(n)*ecko2(i,k,n)

c               of that part is going into rain out (chem_pw)
                tem=chem_c(i,k,n)/(one+c0t(i,k)*dz)
                chem_pw(i,k,n)=c0t(i,k)*dz*tem*eta(i,kk) !etah
                ecko2(i,k,n)=tem+ecko2(i,k,n)-chem_c(i,k,n)

c               pwav needed fo reevaporation in downdraft
c               if including reevaporation, please uncomment code below
c               pwav(i,n)=pwav(i,n)+chem_pw(i,k,n)
              endif
            endif
          enddo
        enddo
        do k = 1, km1
          do i = 1, im
            if (k >= ktcon(i)) ecko2(i,k,n)=ctro2(i,k,n)
          enddo
        enddo
      enddo

c     reevaporation of some, pw and pwd terms needed later for dellae2

!      do n = 1, ntc
!        do k = km1, 1, -1
!          kp1 = k + 1
!          do i = 1, im
!            if (cnvflg(i) .and. (k < jmin(i))) then
!              dz = zi(i,kp1) - zi(i,k)
!              if (k >= kbcon(i)) then
!                tem  = xlamde * dz
!                tem1 = half * xlamdd * dz
!              else
!                tem  = xlamde * dz
!                tem1 = half * (xlamd(i)+xlamdd) * dz
!              endif
!              factor = one + tem - tem1
!              ecdo2(i,k,n) = ((one-tem1)*ecdo2(i,kp1,n)
!     &          +half*tem*(ctro2(i,k,n)+ctro2(i,kp1,n)))/factor
c             if including reevaporation, please uncomment code below
c             ecdo2(i,k,n)=ecdo2(i,k,n)+pwdper(i,kp1)*pwav(i,n)
c             chem_pwd(i,k,n)=max(zero,pwdper(i,kp1)*pwav(i,n))
!            endif
!          enddo
!        enddo
!      enddo

!      do n = 1, ntc
!        do i = 1, im
!          if (cnvflg(i)) then
c           subsidence term treated in fct routine
!            dellae2(i,1,n) = edto(i)*etad(i,1)*ecdo2(i,1,n)*xmbp(i,1)
!          endif
!        enddo
!      enddo

      do n = 1, ntc
        do i = 1, im
          if (cnvflg(i)) then
            k = ktcon(i)
            kk = k - 1
c           for the subsidence term already is considered
            dellae2(i,k,n) = eta(i,kk) * ecko2(i,kk,n) * xmbp(i,k)
          endif
        enddo
      enddo

c     --- for updraft & downdraft vertical transport
c
c     initialize maximum allowed timestep for upstream difference approach
c
      dtime_max=delt
      do k=2,km1
        kk = k - 1
        do i = 1, im
          if (kk < ktcon(i)) dtime_max = min(dtime_max,half*delp(i,kk))
        enddo
      enddo

c     now for every chemistry tracer
      do n = 1, ntc
        do k = 2, km1
          kk = k - 1
          do i = 1, im
            if (cnvflg(i) .and. (k < ktcon(i))) then
              dz = zi(i,k) - zi(i,kk)
              aup = one
              if (k <= kb(i)) aup = zero
!              adw = one
!              if (k > jmin(i)) adw = zero

              dv1q = half * (ecko2(i,k,n) + ecko2(i,kk,n))
              dv2q = half * (ctro2(i,k,n) + ctro2(i,kk,n))
c             dv3q = half * (ecdo2(i,k,n) + ecdo2(i,kk,n))

              tem  = half * (xlamue(i,k) + xlamue(i,kk))
             !tem1 = half * (xlamud(i,k) + xlamud(i,kk))
              tem1 = half * (xlamud(i  ) + xlamud(i   ))

!              if (k <= kbcon(i)) then
!                ptem  = xlamde
!                ptem1 = xlamd(i) + xlamdd
!              else
!                ptem  = xlamde
!                ptem1 = xlamdd
!              endif
              dellae2(i,k,n) = dellae2(i,k,n) +
c                 detrainment from updraft
     &          ( aup*tem1*eta(i,kk)*dv1q
c                 entrainement into up and downdraft
!    &          - (aup*tem*eta(i,kk)+adw*edto(i)*ptem*etad(i,k))*dv2q
     &          - (aup*tem*eta(i,kk))*dv2q
c                 detrainment from downdraft
!    &          + (adw*edto(i)*ptem1*etad(i,k)*dv3q)
     &            ) * dz * xmbp(i,k)

              wet_dep(i,k,n)=chem_pw(i,k,n)*g/delp(i,k)

c             sinks from where updraft and downdraft start
!              if (k == jmin(i)+1) then
!                dellae2(i,k,n) = dellae2(i,k,n)
!     &            -edto(i)*etad(i,kk)*ctro2(i,kk,n)*xmbp(i,k)
!              endif
              if (k == kb(i))then
                dellae2(i,k,n) = dellae2(i,k,n)
     &            -eta(i,k)*ctro2(i,k,n)*xmbp(i,k)
              endif
            endif
          enddo
        enddo

        do i = 1, im
          if (cnvflg(i)) then
            if (kb(i) == 1) then
              k=kb(i)
              dellae2(i,k,n) = dellae2(i,k,n)
     &          -eta(i,k)*ctro2(i,k,n)*xmbp(i,k)
            endif
          endif
        enddo

      enddo

c     for every tracer...

      do n = 1, ntc
        flx_lo  = zero
        totlout = zero
        clipout = zero
c       compute low-order mass flux, upstream
        do k = 2, km1
          kk = k - 1
          do i = 1, im
            if (cnvflg(i) .and. (kk < ktcon(i))) then
              tem = zero
              if (kk >= kb(i)  ) tem = eta(i,kk)
!              if (kk <= jmin(i)) tem = tem - edto(i)*etad(i,kk)
c             low-order flux,upstream
              if (tem > zero) then
                flx_lo(i,k) = -xmb(i) * tem * qaero(i,k,n)
              elseif (tem < zero) then
                flx_lo(i,k) = -xmb(i) * tem * qaero(i,kk,n)
              endif
            endif
          enddo
        enddo

c       --- make sure low-ord fluxes don't violate positive-definiteness
        do k=1,km1
          kp1 = k + 1
          do i=1,im
            if (cnvflg(i) .and. (k <= ktcon(i))) then
c             time step / grid spacing
              dtovdz = g * dtime_max / abs(delp(i,k))
c             total flux out
              totlout(i,k)=max(zero,flx_lo(i,kp1))-min(zero,flx_lo(i,k))
              clipout(i,k)=min(one ,qaero(i,k,n)/max(epsil,totlout(i,k))
     &           / (1.0001_kind_phys*dtovdz))
            endif
          enddo
        enddo

c       recompute upstream mass fluxes
        do k = 2, km1
          kk = k - 1
          do i = 1, im
            if (cnvflg(i) .and. (kk < ktcon(i))) then
              tem = zero
              if (kk >= kb(i)  ) tem = eta(i,kk)
!             if (kk <= jmin(i)) tem = tem - edto(i)*etad(i,kk)
              if (tem > zero) then
                flx_lo(i,k) = flx_lo(i,k) * clipout(i,k)
              elseif (tem < zero) then
                flx_lo(i,k) = flx_lo(i,k) * clipout(i,kk)
              endif
            endif
          enddo
        enddo

c      --- a positive-definite low-order (diffusive) solution for the subsidnce fluxes
       do k=1,km1
         kp1 = k + 1
         do i=1,im
           if (cnvflg(i) .and. (k <= ktcon(i))) then
             dtovdz = g * dtime_max / abs(delp(i,k)) ! time step /grid spacing
             dellae2(i,k,n) = dellae2(i,k,n)
     &         -(flx_lo(i,kp1)-flx_lo(i,k))*dtovdz/dtime_max
           endif
          enddo
        enddo

      enddo ! ctr

c     convert wet deposition to total mass deposited over dt and dp

      do n = 1, ntc
        do k = 1, km
          do i = 1, im
            if (cnvflg(i) .and. (k < ktcon(i)))
     &        wet_dep(i,k,n) = wet_dep(i,k,n)*xmb(i)*delt*delp(i,k)
          enddo
        enddo
      enddo

c     compute final aerosol concentrations

      do n = 1, ntc
        do k = 1, km
          do i = 1, im
            if (cnvflg(i) .and. (k <= min(kmax(i),ktcon(i)))) then
              qaero(i,k,n) = qaero(i,k,n) + dellae2(i,k,n) * delt
              if (qaero(i,k,n) < zero) then
c               add negative mass to wet deposition
                wet_dep(i,k,n) = wet_dep(i,k,n)-qaero(i,k,n)*delp(i,k)
                qaero(i,k,n) = qamin
              endif
            endif
          enddo
        enddo
      enddo

      return
      end subroutine samfshalcnv_aerosols

      end module samfcnv_aerosols

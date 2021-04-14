!> \file gcm_shoc.F90
!!  Contains the Simplified Higher-Order Closure (SHOC) scheme.

!> This module contains the CCPP-compliant SHOC scheme.
module shoc
  use machine, only: kind_phys

  implicit none

  private

  public shoc_run, shoc_init, shoc_finalize
  integer, parameter :: kp = kind_phys

contains

subroutine shoc_init ()
end subroutine shoc_init

subroutine shoc_finalize ()
end subroutine shoc_finalize

!> \section arg_table_shoc_run Argument Table
!! \htmlinclude shoc_run.html
!!
subroutine shoc_run (nx, nzm, tcr, tcrf, con_cp, con_g, con_hvap, con_hfus, con_rv, con_rd,     &
                     con_pi, con_fvirt, con_eps, dtp, prsl, delp, phii, phil, u, v, omega, rhc, &
                     supice, pcrit,  cefac, cesfac, tkef1, dis_opt, hflx, evap, prnum,          &
                     gt0, gq0, ntrac, ntqv, ntcw, ntiw, ntrw, ntsw, ntgl, ntlnc, ntinc,         &
                     cld_sgs, tke, tkh, wthv_sec, errmsg, errflg)

    implicit none

    integer, intent(in) :: nx, nzm, ntrac, ntqv, ntcw, ntiw, ntrw, ntsw, ntgl, ntlnc, ntinc
    real(kind=kind_phys), intent(in) :: tcr, tcrf, con_cp, con_g, con_hvap, con_hfus, con_rv, &
                                        con_rd, con_pi, con_fvirt, con_eps,                   &
                                        dtp, supice, pcrit, cefac, cesfac, tkef1, dis_opt
  !
    real(kind=kind_phys), intent(in), dimension(nx)       :: hflx, evap
    real(kind=kind_phys), intent(in), dimension(nx,nzm)   :: prsl, delp, phil, u, v, omega, rhc, prnum
    real(kind=kind_phys), intent(in), dimension(nx,nzm+1) :: phii
   !
    real(kind=kind_phys), intent(inout), dimension(nx,nzm) :: gt0, cld_sgs, tke, tkh, wthv_sec
    real(kind=kind_phys), intent(inout), dimension(nx,nzm,ntrac) :: gq0

   character(len=*), intent(out) :: errmsg
   integer,          intent(out) :: errflg

   real(kind=kind_phys), parameter :: epsq = 1.0e-20_kp, zero=0.0_kp, one=1.0_kp

   integer :: i, k

   real(kind=kind_phys) :: tem
   real(kind=kind_phys), dimension(nx,nzm) :: qi   ! local array of suspended cloud ice
   real(kind=kind_phys), dimension(nx,nzm) :: qc   ! local array of suspended cloud water
   real(kind=kind_phys), dimension(nx,nzm) :: qsnw ! local array of suspended snowq
   real(kind=kind_phys), dimension(nx,nzm) :: qrn  ! local array of suepended rain
   real(kind=kind_phys), dimension(nx,nzm) :: qgl  ! local array of suspended graupel
   real(kind=kind_phys), dimension(nx,nzm) :: ncpl ! local array of cloud water number concentration
   real(kind=kind_phys), dimension(nx,nzm) :: ncpi ! local array of cloud ice number concentration

! Initialize CCPP error handling variables

    errmsg = ''
    errflg = 0

    if (ntiw < 0) then   ! this is valid only for Zhao-Carr scheme
      do k=1,nzm
         do i=1,nx
           qc(i,k) = gq0(i,k,ntcw)
           if (abs(qc(i,k)) < epsq) then
             qc(i,k) = zero
           endif
           tem = qc(i,k) * max(zero, MIN(one, (tcr-gt0(i,k))*tcrf))
           qi(i,k) = tem                              ! ice
           qc(i,k) = qc(i,k) - tem         ! water
           qrn(i,k)  = zero
           qsnw(i,k) = zero
           ncpl(i,k) = zero
           ncpi(i,k) = zero
         enddo
       enddo
    else
      if (ntgl > 0) then    ! graupel exists - combine graupel with snow
        do k=1,nzm
          do i=1,nx
            qc(i,k)   = gq0(i,k,ntcw)
            qi(i,k)   = gq0(i,k,ntiw)
            qrn(i,k)  = gq0(i,k,ntrw)
            qsnw(i,k) = gq0(i,k,ntsw) + gq0(i,k,ntgl)
           enddo
         enddo
      else                  ! no graupel
        do k=1,nzm
          do i=1,nx
            qc(i,k)   = gq0(i,k,ntcw)
            qi(i,k)   = gq0(i,k,ntiw)
            qrn(i,k)  = gq0(i,k,ntrw)
            qsnw(i,k) = gq0(i,k,ntsw)
           enddo
         enddo
      endif
      if (ntlnc > 0) then
        do k=1,nzm
          do i=1,nx
            ncpl(i,k) = gq0(i,k,ntlnc)
            ncpi(i,k) = gq0(i,k,ntinc)
          enddo
        enddo
      endif
    endif

    !     phy_f3d(1,1,ntot3d-2) - shoc determined sgs clouds
    !     phy_f3d(1,1,ntot3d-1) - shoc determined diffusion coefficients
    !     phy_f3d(1,1,ntot3d  ) - shoc determined  w'theta'

    call shoc_work (nx, nx, nzm, nzm+1, dtp, prsl, delp,                                &
                    phii, phil, u, v, omega, gt0, gq0(:,:,1), qi, qc, qsnw, qrn,        &
                    rhc, supice, pcrit, cefac, cesfac, tkef1, dis_opt,                  &
                    cld_sgs, tke, hflx, evap, prnum, tkh, wthv_sec,                     &
                    ntlnc, ncpl, ncpi,                                                  &
                    con_cp, con_g, con_hvap, con_hfus, con_rv, con_rd, con_pi,          &
                    con_fvirt, con_eps)

    if (ntiw < 0) then   ! this is valid only for Zhao-Carr scheme
      do k=1,nzm
         do i=1,nx
           gq0(i,k,ntcw) = qc(i,k) + qi(i,k)
         enddo
       enddo
    else
      do k=1,nzm
        do i=1,nx
          gq0(i,k,ntcw) = qc(i,k)
          gq0(i,k,ntiw) = qi(i,k)
        enddo
      enddo
      if (ntlnc > 0) then
        do k=1,nzm
          do i=1,nx
            gq0(i,k,ntlnc) = ncpl(i,k)
            gq0(i,k,ntinc) = ncpi(i,k)
          enddo
        enddo
      endif
    endif

end subroutine shoc_run

 ! Implementation of the Simplified High Order Closure (SHOC) scheme
 ! of Bogenschutz and Krueger (2013), J. Adv. Model. Earth Syst, 5, 195-211,
 ! doi: 10.1002/jame.200118. (further referred to as BK13)
 ! in a single column form suitable for use in a GCM physics package.
 ! Alex Belochitski, heavily based on the code of Peter Bogenschutz.
 ! S Moorthi - optimization, cleanup, improve and customize for gsm
 !           - improved solution for sgs-tke equation
 ! S Moorthi - 05-11-17 - modified shear production term to eliminate
 !                        spurious tke ove Antarctica.
 ! S Moorthi - 01-12-17 - added extra pressure dependent tke dissipation at
 !                        pressures below a critical value pcrit
 ! S Moorthi - 04-12-17 - fixed a bug in the definition of hl on input
 !                        replacing fac_fus by fac_sub
 ! S.Moorthi - 00-00-17 - added an alternate option for near boundary cek following
 !                        Scipion et. al., from U. Oklahoma.
 subroutine shoc_work (ix, nx, nzm, nz, dtn,                            &
                       prsl, delp, phii, phil, u, v, omega, tabs,       &
                       qwv, qi, qc, qpi, qpl, rhc, supice,              &
                       pcrit, cefac, cesfac, tkef1, dis_opt,            &
                       cld_sgs, tke, hflx, evap, prnum, tkh,            &
                       wthv_sec, ntlnc, ncpl, ncpi,                     &
                       cp, ggr, lcond, lfus, rv, rgas, pi, epsv, eps)

  use funcphys , only : fpvsl, fpvsi, fpvs    ! saturation vapor pressure for water & ice

  implicit none

  real,    intent(in) :: cp, ggr, lcond, lfus, rv, rgas, pi, epsv, eps
  integer, intent(in) :: ix      ! max number of points in the physics window in the x
  integer, intent(in) :: nx      ! Number of points in the physics window in the x

  integer, intent(in) :: nzm     ! Number of vertical layers
  integer, intent(in) :: nz      ! Number of layer interfaces  (= nzm + 1)
  integer, intent(in) :: ntlnc   ! index of liquid water number concentration
  real,    intent(in) :: dtn     ! Physics time step, s

  real,    intent(in) :: pcrit   ! pressure in Pa below which additional tke dissipation is applied
  real,    intent(in) :: cefac   ! tunable multiplier to dissipation term
  real,    intent(in) :: cesfac  ! tunable multiplier to dissipation term for bottom level
  real,    intent(in) :: tkef1   ! uncentering terms in implicit tke integration
  real,    intent(in) :: dis_opt ! when > 0 use different formula for near surface dissipation

  real,    intent(in) :: hflx(nx)
  real,    intent(in) :: evap(nx)

! The interface is talored to GFS in a sense that input variables are 2D

  real, intent(in)    :: prsl   (ix,nzm)   ! mean layer presure
  real, intent(in)    :: delp   (ix,nzm)   ! layer presure depth
  real, intent(in)    :: phii   (ix,nz )   ! interface geopotential height
  real, intent(in)    :: phil   (ix,nzm)   ! layer geopotential height
  real, intent(in)    :: u      (ix,nzm)   ! u-wind, m/s
  real, intent(in)    :: v      (ix,nzm)   ! v-wind, m/s
  real, intent(in)    :: omega  (ix,nzm)   ! omega, Pa/s
  real, intent(inout) :: tabs   (ix,nzm)   ! temperature, K
  real, intent(inout) :: qwv    (ix,nzm)   ! water vapor mixing ratio, kg/kg
  real, intent(inout) :: qc     (ix,nzm)   ! cloud water mixing ratio, kg/kg
  real, intent(inout) :: qi     (ix,nzm)   ! cloud ice   mixing ratio, kg/kg
! Anning Cheng 03/11/2016 SHOC feedback to number concentration
  real, intent(inout) :: ncpl   (nx,nzm)   ! cloud water number concentration,/m^3
  real, intent(inout) :: ncpi   (nx,nzm)   ! cloud ice   number concentration,/m^3
  real, intent(in)    :: qpl    (nx,nzm)   ! rain    mixing ratio, kg/kg
  real, intent(in)    :: qpi    (nx,nzm)   ! snow    mixing ratio, kg/kg

  real, intent(in)    :: rhc    (nx,nzm)   ! critical relative humidity
  real, intent(in)    :: supice            ! ice supersaturation parameter
  real, intent(out)   :: cld_sgs(ix,nzm)   ! sgs cloud fraction
! real, intent(inout) :: cld_sgs(nx,nzm)   ! sgs cloud fraction
  real, intent(inout) :: tke    (ix,nzm)   ! turbulent kinetic energy. m**2/s**2
! real, intent(inout) :: tk     (nx,nzm)   ! eddy viscosity
  real, intent(inout) :: tkh    (ix,nzm)   ! eddy diffusivity
  real, intent(in)    :: prnum  (nx,nzm)   ! turbulent Prandtl number
  real, intent(inout) :: wthv_sec (ix,nzm) ! Buoyancy flux, K*m/s

  real, parameter :: zero=0.0_kp,  one=1.0_kp,  half=0.5_kp, two=2.0_kp,                  &
                     three=3.0_kp, oneb3=one/three, twoby3=two/three, fourb3=twoby3+twoby3
  real, parameter :: sqrt2 = sqrt(two), twoby15 = two / 15.0_kp,                          &
                     nmin = 1.0_kp,    RI_cub = 6.4e-14_kp, RL_cub = 1.0e-15_kp,          &
                     skew_facw=1.2_kp, skew_fact=0.0_kp,                                  &
                     tkhmax=300.0_kp,  qcmin=1.0e-9_kp
  real            :: lsub, fac_cond, fac_fus, cpolv, fac_sub, ggri, kapa, gocp,           &
                     rog, sqrtpii, epsterm, onebeps, onebrvcp

! SHOC tunable parameters

  real, parameter :: lambda  = 0.04_kp
! real, parameter :: min_tke = 1.0e-6_kp  ! Minumum TKE value, m**2/s**2
  real, parameter :: min_tke = 1.0e-4_kp  ! Minumum TKE value, m**2/s**2
! real, parameter :: max_tke = 100.0_kp ! Maximum TKE value, m**2/s**2
  real, parameter :: max_tke = 40.0_kp  ! Maximum TKE value, m**2/s**2
! Maximum turbulent eddy length scale, m
! real, parameter :: max_eddy_length_scale  = 2000.0_kp
  real, parameter :: max_eddy_length_scale  = 1000.0_kp
! Maximum "return-to-isotropy" time scale, s
  real, parameter :: max_eddy_dissipation_time_scale = 2000.0_kp
  real, parameter :: Pr    = 1.0_kp           ! Prandtl number

! Constants for the TKE dissipation term based on Deardorff (1980)
  real, parameter :: pt19=0.19_kp,  pt51=0.51_kp, pt01=0.01_kp, atmin=0.01_kp, atmax=one-atmin
  real, parameter :: Cs  = 0.15_kp, epsln=1.0e-6_kp
! real, parameter :: Ck  = 0.2_kp     ! Coeff in the eddy diffusivity - TKE relationship, see Eq. 7 in BK13
  real, parameter :: Ck  = 0.1_kp     ! Coeff in the eddy diffusivity - TKE relationship, see Eq. 7 in BK13 

! real, parameter :: Ce  = Ck**3/(0.7*Cs**4)
! real, parameter :: Ce  = Ck**3/(0.7*Cs**4) * 2.2
! real, parameter :: Ce  = Ck**3/(0.7*Cs**4) * 3.0 , Ces = Ce
! real, parameter :: Ce  = Ck**3/(0.7*Cs**4) * 2.5 , Ces = Ce * 3.0 / 2.5
! real, parameter :: Ces = Ce/0.7*3.0

! real, parameter :: Ce  = Ck**3/(0.7*Cs**4), Ces = Ce*3.0/0.7 ! Commented Moor

  real, parameter :: Ce  = Ck**3/Cs**4, Ces = Ce
! real, parameter :: Ce  = Ck**3/Cs**4, Ces = Ce*3.0/0.7

! real, parameter :: vonk=0.35                ! Von Karman constant
  real, parameter :: vonk=0.4_kp              ! Von Karman constant Moorthi - as in GFS
  real, parameter :: tscale=400.0_kp          ! time scale set based off of similarity results of BK13, s
  real, parameter :: w_tol_sqd = 4.0e-04_kp   ! Min vlaue of second moment of w
! real, parameter :: w_tol_sqd = 1.0e-04_kp   ! Min vlaue of second moment of w
  real, parameter :: w_thresh  = 0.0_kp, thresh = 0.0_kp
  real, parameter :: w3_tol    = 1.0e-20_kp   ! Min vlaue of third moment of w


! These parameters are a tie-in with a microphysical scheme
! Double check their values for the Zhao-Carr scheme.
  real, parameter :: tbgmin = 233.16_kp    ! Minimum temperature for cloud water., K (ZC)
! real, parameter :: tbgmin = 258.16_kp    ! Minimum temperature for cloud water., K (ZC)
! real, parameter :: tbgmin = 253.16_kp    ! Minimum temperature for cloud water., K
  real, parameter :: tbgmax = 273.16_kp    ! Maximum temperature for cloud ice, K
  real, parameter :: a_bg   = one/(tbgmax-tbgmin)
!
! Parameters to tune the second order moments-  No tuning is performed currently

! real, parameter :: thl2tune = 2.0_kp,   qw2tune = 2.0_kp,  qwthl2tune = 2.0_kp, &
  real, parameter :: thl2tune = 1.0_kp,   qw2tune = 1.0_kp,  qwthl2tune = 1.0_kp, &
!                    thl_tol  = 1.0e-4_kp, rt_tol = 1.0e-8_kp, basetemp  = 300.0_kp
                     thl_tol  = 1.0e-2_kp, rt_tol = 1.0e-4_kp

  integer, parameter :: nitr=6

! Local variables. Note that pressure is in millibars in the SHOC code.

  real zl      (nx,nzm)  ! height of the pressure levels above surface, m
  real zi      (nx,nz)   ! height of the interface levels, m
  real adzl    (nx,nzm)  ! layer thickness i.e. zi(k+1)-zi(k) - defined at levels
  real adzi    (nx,nz)   ! level thickness i.e. zl(k)-zl(k-1) - defined at interface

  real hl      (nx,nzm)  ! liquid/ice water static energy , K
  real qv      (nx,nzm)  ! water vapor, kg/kg
  real qcl     (nx,nzm)  ! liquid water  (condensate), kg/kg
  real qci     (nx,nzm)  ! ice water  (condensate), kg/kg
  real w       (nx,nzm)  ! z-wind, m/s
  real bet     (nx,nzm)  ! ggr/tv0
  real gamaz   (nx,nzm)  ! ggr/cp*z

! Moments of the trivariate double Gaussian PDF for the SGS total water mixing ratio
! SGS liquid/ice static energy, and vertical velocity

  real qw_sec   (nx,nzm) ! Second moment total water mixing ratio, kg^2/kg^2
  real thl_sec  (nx,nzm) ! Second moment liquid/ice static energy, K^2
  real qwthl_sec(nx,nzm) ! Covariance tot. wat. mix. ratio and static energy, K*kg/kg
  real wqw_sec  (nx,nzm) ! Turbulent flux of tot. wat. mix., kg/kg*m/s
  real wthl_sec (nx,nzm) ! Turbulent flux of liquid/ice static energy, K*m/s
  real w_sec    (nx,nzm) ! Second moment of vertical velocity, m**2/s**2
  real w3       (nx,nzm) ! Third moment of vertical velocity, m**3/s**3
  real wqp_sec  (nx,nzm) ! Turbulent flux of precipitation, kg/kg*m/s

! Eddy length formulation
  real smixt    (nx,nzm) ! Turbulent length scale, m
  real isotropy (nx,nzm) ! "Return-to-isotropy" eddy dissipation time scale, s
! real isotropy_debug (nx,nzm) ! Return to isotropy scale, s without artificial limits
  real brunt    (nx,nzm) ! Moist Brunt-Vaisalla frequency, s^-1
  real conv_vel2(nx,nzm) ! Convective velocity scale cubed, m^3/s^3

  real cek(nx)

! Output of SHOC
  real diag_frac, diag_qn, diag_qi, diag_ql

! real diag_frac(nx,nzm) ! SGS cloud fraction
! real diag_qn  (nx,nzm) ! SGS cloud+ice condensate, kg/kg
! real diag_qi  (nx,nzm) ! SGS ice condensate, kg/kg
! real diag_ql  (nx,nzm) ! SGS liquid condensate, kg/kg


! Horizontally averaged variables
! real conv_vel(nzm)        ! Convective velocity scale cubed, m^3/s^3
  real wqlsb   (nzm)        ! liquid water flux, kg/kg/ m/s
  real wqisb   (nzm)        ! ice flux, kg/kg m/s
! real thlv    (nzm)        ! Grid-scale level-average virtual potential temperature
!                              (not used)


! Local variables

! real, dimension(nx,nzm) :: tkesbbuoy, tkesbshear, tkesbdiss, tkesbbuoy_debug   &
!                               tkebuoy_sgs, total_water, tscale1_debug, brunt2

  real, dimension(nx,nzm) :: total_water, brunt2, thv, tkesbdiss
  real, dimension(nx,nzm) :: def2
  real, dimension(nx)     :: denom, numer, l_inf, cldarr, thedz, thedz2

  real lstarn,    depth,    omn,         betdz,      bbb,      term,   qsatt, dqsat,        &
                  conv_var,  tkes,       skew_w,     skew_qw,  aterm,  w1_1, w1_2,  w2_1,   &
       w2_2,      w3var,     thl1_1,     thl1_2,     thl2_1,   thl2_2, qw1_1, qw1_2, qw2_1, &
       qw2_2,     ql1,       ql2,        w_ql1,      w_ql2,                                 &
       r_qwthl_1, r_wqw_1,   r_wthl_1,   testvar,    s1,    s2,    std_s1, std_s2, C1, C2,  &
       thl_first, qw_first,  w_first,    Tl1_1,      Tl1_2, betatest, pval, pkap,           &
       w2thl,     w2qw,w2ql, w2ql_1,     w2ql_2,                                            &
       thec,      thlsec,    qwsec,      qwthlsec,   wqwsec, wthlsec, thestd,dum,           &
       cqt1,      cthl1,     cqt2,       cthl2,      qn1,    qn2, qi1, qi2, omn1, omn2,     &
       basetemp2, beta1,     beta2,      qs1,        qs2,                                   &
       esval,     esval2,    om1,        om2,        epss,                                  &
       lstarn1,   lstarn2,   sqrtw2,     sqrtthl,    sqrtqt,                                &
       sqrtstd1,  sqrtstd2,  tsign,      tvar,       sqrtw2t, wqls,  wqis,                  &
       sqrtqw2_1, sqrtqw2_2, sqrtthl2_1, sqrtthl2_2, sm,      prespot,                      &
       corrtest1, corrtest2, wrk,        wrk1,       wrk2,    wrk3,    onema, pfac


  integer i,k,km1,ku,kd,ka,kb


!calculate derived constants
  lsub     = lcond+lfus
  fac_cond = lcond/cp
  fac_fus  = lfus/cp
  cpolv    = cp/lcond
  fac_sub  = lsub/cp
  ggri     = one/ggr
  kapa     = rgas/cp
  gocp     = ggr/cp
  rog      = rgas*ggri
  sqrtpii  = one/sqrt(pi+pi)
  epsterm  = rgas/rv
  onebeps  = one/epsterm
  onebrvcp = one/(rv*cp)
  epss     = eps * supice

! Map GFS variables to those of SHOC - SHOC operates on 3D fields
! Here a Y-dimension is added to the input variables, along with some unit conversions

  do k=1,nz
    do i=1,nx
      zi(i,k) = phii(i,k) * ggri
    enddo
  enddo

!
! move water from vapor to condensate if the condensate is negative
!
  do k=1,nzm
    do i=1,nx
      if (qc(i,k) < zero) then
        qwv(i,k)  = qwv(i,k) + qc(i,k)
        tabs(i,k) = tabs(i,k) - fac_cond * qc(i,k)
        qc(i,k)   = zero
      endif
      if (qi(i,k) < zero) then
        qwv(i,k)  = qwv(i,k) + qi(i,k)
        tabs(i,k) = tabs(i,k) - fac_sub  * qi(i,k)
        qi(i,k)   = zero
      endif
!
!    testing removal of ice when too warm to sustain ice 
!
!     if (qi(i,k) > zero .and. tabs(i,k) > 273.16) then
!       wrk = (tabs(i,k) - 273.16) / fac_sub
!       if (wrk < qi(i,k)) then
!         wrk       = qi(i,k) - wrk
!         qi(i,k)   = wrk
!         qwv(i,k)  = qwv(i,k) + wrk
!         tabs(i,k) = 273.16
!       else
!         tabs(i,k) = tabs(i,k) - qi(i,k) / fac_sub
!         qwv(i,k)  = qwv(i,k)  + qi(i,k)
!         qi(i,k)   = 0.0
!       endif
!     endif

    enddo
  enddo
! fill negative water vapor from below
  do k=nzm,2,-1
    km1 = k - 1
    do i=1,nx
      if (qwv(i,k) < zero) then
        qwv(i,k) = qwv(i,km1) + qwv(i,k) * delp(i,k) / delp(i,km1)
      endif
    enddo
  enddo

  do k=1,nzm
    do i=1,nx
      zl(i,k)    = phil(i,k) * ggri
      wrk        = one / prsl(i,k)
      qv(i,k)    = max(qwv(i,k), zero)
      thv(i,k)   = tabs(i,k) * (one+epsv*qv(i,k))
      w(i,k)     = - rog * omega(i,k) * thv(i,k) * wrk
      qcl(i,k)   = max(qc(i,k), zero)
      qci(i,k)   = max(qi(i,k), zero)
!
!     qpl(i,k)     = zero  ! comment or remove when using with prognostic rain/snow
!     qpi(i,k)     = zero  ! comment or remove when using with prognostic rain/snow

      wqp_sec(i,k) = zero  ! Turbulent flux of precipiation
!
      total_water(i,k) = qcl(i,k) + qci(i,k) + qv(i,k)

      prespot      = (100000.0_kp*wrk) ** kapa ! Exner function
      bet(i,k)     = ggr/(tabs(i,k)*prespot)     ! Moorthi
      thv(i,k)     = thv(i,k)*prespot            ! Moorthi
!
! Lapse rate * height = reference temperature
      gamaz(i,k) = gocp * zl(i,k)

! Liquid/ice water static energy - ! Note the the units are degrees K
      hl(i,k) = tabs(i,k) + gamaz(i,k) - fac_cond*(qcl(i,k)+qpl(i,k)) &
                                       - fac_sub *(qci(i,k)+qpi(i,k))
      w3(i,k) = zero
    enddo
  enddo

! Define vertical grid increments for later use in the vertical differentiation

  do k=2,nzm
    km1 = k - 1
    do i=1,nx
      adzi(i,k)   = zl(i,k) - zl(i,km1)
      adzl(i,km1) = zi(i,k) - zi(i,km1)
    enddo
  enddo
  do i=1,nx
    adzi(i,1)     = (zl(i,1)-zi(i,1))   !  unused in the code
    adzi(i,nz)    = adzi(i,nzm)         ! at the top - probably unused
    adzl(i,nzm)   = zi(i,nz) - zi(i,nzm)
!
    wthl_sec(i,1) = hflx(i)
    wqw_sec(i,1)  = evap(i)
  enddo


  call tke_shoc()        ! Integrate prognostic TKE equation forward in time


! diagnose second order moments of the subgrid PDF following
! Redelsperger J.L., and G. Sommeria, 1986, JAS, 43, 2619-2635 sans the use of stabilty
! weighting functions - Result is in global variables w_sec, thl_sec, qw_sec, and qwthl_sec

! call diag_moments(total_water,tke,tkh)

! Second moment of vertical velocity.
! Note that Eq 6 in BK13 gives a different expression that is dependent on
! vertical gradient of grid scale vertical velocity

  do k=1,nzm
    ku = k+1
    kd = k-1
    ka = ku
    kb = k
    if (k == 1) then
      kd = k
      kb = ka
    elseif (k == nzm) then
      ku = k
      ka = kb
    endif
    do i=1,nx
      if (tke(i,k) > zero) then
!        wrk  = half*(tkh(i,ka)+tkh(i,kb))*(w(i,ku) - w(i,kd)) &
         wrk  = half*(tkh(i,ka)*prnum(i,ka)+tkh(i,kb)*prnum(i,kb))*(w(i,ku) - w(i,kd)) &
              * sqrt(tke(i,k)) / (zl(i,ku) - zl(i,kd))
         w_sec(i,k) = max(twoby3 * tke(i,k) - twoby15 * wrk, zero)
!        w_sec(i,k) = max(twoby3 * tke(i,k), zero)
      else
         w_sec(i,k) = zero
      endif
    enddo
  enddo

  do k=2,nzm

    km1 = k - 1
    do i=1,nx

! Use backward difference in the vertical, use averaged values of "return-to-isotropy"
! time scale and diffusion coefficient

       wrk1 = one / adzi(i,k)        ! adzi(k) = (zl(k)-zl(km1))
!      wrk3 = max(tkh(i,k),pt01) * wrk1
       wrk3 = max(tkh(i,k),epsln) * wrk1

       sm   = half*(isotropy(i,k)+isotropy(i,km1))*wrk1*wrk3 ! Tau*Kh/dz^2

! SGS vertical flux liquid/ice water static energy. Eq 1 in BK13
!               No rain, snow or graupel in pdf (Annig, 08/29/2018)

       wrk1            =  hl(i,k)  - hl(i,km1)              &
                       + (qpl(i,k) - qpl(i,km1)) * fac_cond &
                       + (qpi(i,k) - qpi(i,km1)) * fac_sub
       wthl_sec(i,k) = - wrk3 * wrk1

! SGS vertical flux of total water. Eq 2 in BK13

       wrk2           = total_water(i,k) - total_water(i,km1)
       wqw_sec(i,k) = - wrk3 * wrk2

! Second moment of liquid/ice water static energy. Eq 4 in BK13

       thl_sec(i,k) = thl2tune * sm * wrk1 * wrk1

! Second moment of total water mixing ratio.  Eq 3 in BK13

       qw_sec(i,k)  = qw2tune * sm * wrk2 * wrk2

! Covariance of total water mixing ratio and liquid/ice water static energy.
! Eq 5 in BK13

       qwthl_sec(i,k) = qwthl2tune * sm * wrk1 * wrk2

    enddo   ! i  loop
  enddo     ! k  loop

!   These would be at the surface - do we need them?
  do i=1,nx
!   wthl_sec(i,1)  = wthl_sec(i,2)
!   wqw_sec(i,1)   = wqw_sec(i,2)
    thl_sec(i,1)   = thl_sec(i,2)
    qw_sec(i,1)    = qw_sec(i,2)
    qwthl_sec(i,1) = qwthl_sec(i,2)
  enddo

! Diagnose the third moment of SGS vertical velocity

  call canuto()

! Recover parameters of the subgrid PDF using diagnosed moments
! and calculate SGS cloudiness, condensation and it's effects on temeperature
! and moisture variables

  call assumed_pdf()

contains

  subroutine tke_shoc()

! This subroutine solves the TKE equation,
! Heavily based on SAM's tke_full.f90 by Marat Khairoutdinov

    real grd,betdz,Cee,lstarn, lstarp, bbb, omn, omp,qsatt,dqsat, smix,             &
         buoy_sgs,ratio,a_prod_sh,a_prod_bu,a_diss,a_prod_bu_debug, buoy_sgs_debug, &
         tscale1, wrk, wrk1, wtke, wtk2, rdtn, tkef2
    integer i,k,ku,kd,itr,k1

    rdtn = one / dtn

    call tke_shear_prod(def2)   ! Calculate shear production of TKE

! Ensure values of TKE are reasonable

    do k=1,nzm
      do i=1,nx
        tke(i,k)        = max(min_tke,tke(i,k))
        tkesbdiss(i,k)  = zero
!       tkesbshear(i,k) = zero
!       tkesbbuoy(i,k)  = zero
      enddo
    enddo

    call eddy_length()   ! Find turbulent mixing length
    call check_eddy()    ! Make sure it's reasonable

    tkef2 = one - tkef1
    do k=1,nzm
      ku = k+1
      kd = k

!     Cek = Ce * cefac

      if(k == 1) then
        ku = 2
        kd = 2
!       Cek = Ces
      elseif(k == nzm) then
        ku = k
        kd = k
!       Cek = Ces
      endif

      if (dis_opt > 0) then
        do i=1,nx
          wrk = (zl(i,k)-zi(i,1)) / adzl(i,1) + 1.5_kp
          cek(i) = (one + two / max((wrk*wrk - 3.3_kp), 0.5_kp)) * cefac
        enddo
      else
        if (k == 1) then
          cek = ces * cesfac
        else
          cek = ce  * cefac
        endif
      endif

      do i=1,nx
        grd = adzl(i,k)             !  adzl(k) = zi(k+1)-zi(k)


! TKE boyancy production term. wthv_sec (buoyancy flux) is calculated in
! assumed_pdf(). The value used here is from the previous time step

        a_prod_bu = ggr / thv(i,k) * wthv_sec(i,k)

! If wthv_sec from subgrid PDF is not available use Brunt-Vaisalla frequency from eddy_length()

!Obtain Brunt-Vaisalla frequency from diagnosed SGS buoyancy flux
!Presumably it is more precise than BV freq. calculated in  eddy_length()?

        buoy_sgs = - (a_prod_bu+a_prod_bu) / (tkh(i,ku)+tkh(i,kd) + 0.0001_kp)   ! tkh is eddy thermal diffussivity


!Compute $c_k$ (variable Cee) for the TKE dissipation term following Deardorff (1980)

        if (buoy_sgs <= zero) then
          smix = grd
        else
          smix = min(grd,max(0.1_kp*grd, 0.76_kp*sqrt(tke(i,k)/(buoy_sgs+1.0e-10_kp))))
        endif

        ratio     = smix/grd
        Cee       = Cek(i) * (pt19 + pt51*ratio) * max(one, sqrt(pcrit/prsl(i,k)))

! TKE shear production term
        a_prod_sh = half*(def2(i,ku)*tkh(i,ku)*prnum(i,ku)   &
                        + def2(i,kd)*tkh(i,kd)*prnum(i,kd))


! smixt (turb. mixing lenght) is calculated in eddy_length() 
! Explicitly integrate TKE equation forward in time
!       a_diss     = Cee/smixt(i,k)*tke(i,k)**1.5 ! TKE dissipation term
!       tke(i,k) = max(zero,tke(i,k)+dtn*(max(zero,a_prod_sh+a_prod_bu)-a_diss))

! Semi-implicitly integrate TKE equation forward in time

        wtke = tke(i,k)
        wtk2 = wtke
!       wrk  = (dtn*Cee)/smixt(i,k)
        wrk  = (dtn*Cee) / smixt(i,k)
        wrk1 = wtke + dtn*(a_prod_sh+a_prod_bu)

        do itr=1,nitr                        ! iterate for implicit solution
          wtke   = min(max(min_tke, wtke), max_tke)
          a_diss = wrk*sqrt(wtke)            ! Coefficient in the TKE dissipation term
          wtke   = wrk1 / (one+a_diss)
          wtke   = tkef1*wtke + tkef2*wtk2   ! tkef1+tkef2 = 1.0
          wtk2   = wtke
        enddo

        tke(i,k) = min(max(min_tke, wtke), max_tke)
        a_diss   = wrk*sqrt(tke(i,k))

        tscale1  = (dtn+dtn) / a_diss              ! corrected Eq 8 in BK13 -- tau = 2*tke/eps

        tkesbdiss(i,k) = rdtn*a_diss*tke(i,k) ! TKE dissipation term, epsilon


! Calculate "return-to-isotropy" eddy dissipation time scale, see Eq. 8 in BK13

        if (buoy_sgs <= zero) then
          isotropy(i,k) = min(max_eddy_dissipation_time_scale, tscale1)
        else
          isotropy(i,k) = min(max_eddy_dissipation_time_scale,          &
                           tscale1/(one+lambda*buoy_sgs*tscale1*tscale1))
        endif

! TKE budget terms

!       tkesbdiss(i,k)       = a_diss
!       tkesbshear(i,k)      = a_prod_sh
!       tkesbbuoy(i,k)       = a_prod_bu
!       tkesbbuoy_debug(i,k) = a_prod_bu_debug
!       tkebuoy_sgs(i,k)     = buoy_sgs

      enddo ! i loop
    enddo   ! k loop
    wrk = half * ck
    do k=2,nzm
      k1 = k - 1
      do i=1,nx
        tkh(i,k) = min(tkhmax, wrk * (isotropy(i,k)  * tke(i,k)    &
                                   +  isotropy(i,k1) * tke(i,k1))) ! Eddy thermal diffusivity
      enddo ! i
    enddo   ! k


  end subroutine tke_shoc


  subroutine tke_shear_prod(def2)

! Calculate TKE shear production term

    real, intent(out) :: def2(nx,nzm)

    real    rdzw, wrku, wrkv, wrkw
    integer i,k,k1

! Calculate TKE shear production term  at layer interface

    do k=2,nzm
      k1 = k - 1
      do i=1,nx
        rdzw      = one / adzi(i,k)
        wrku      = (u(i,k)-u(i,k1)) * rdzw
        wrkv      = (v(i,k)-v(i,k1)) * rdzw
!       wrkw      = (w(i,k)-w(i,k1)) * rdzw
        def2(i,k) = wrku*wrku + wrkv*wrkv !+ 2*wrkw(1) * wrkw(1)
      enddo
    enddo     ! k  loop
    do i=1,nx
!     def2(i,1) = def2(i,2)
      def2(i,1) = (u(i,1)*u(i,1) + v(i,1)*v(i,1)) / (zl(i,1)*zl(i,1))
    enddo

  end subroutine tke_shear_prod

  subroutine eddy_length()

! This subroutine computes the turbulent length scale based on a new
! formulation described in BK13

! Local variables
    real    wrk, wrk1, wrk2, wrk3
    integer i, k, kk, kl, ku, kb, kc, kli, kui

    do i=1,nx
      cldarr(i) = zero
      numer(i)  = zero
      denom(i)  = zero
    enddo

! Find the length scale outside of clouds, that includes boundary layers.

    do k=1,nzm
      do i=1,nx

! Reinitialize the mixing length related arrays to zero
!       smixt(i,k)    = one   ! shoc_mod module variable smixt
        smixt(i,k)    = epsln ! shoc_mod module variable smixt
        brunt(i,k)    = zero

!Eq. 11 in BK13 (Eq. 4.13 in Pete's dissertation)
!Outside of cloud, integrate from the surface to the cloud base
!Should the 'if' below check if the cloud liquid < a small constant instead?

        if (qcl(i,k)+qci(i,k) <= qcmin) then
          tkes       = sqrt(tke(i,k)) * adzl(i,k)
          numer(i) = numer(i) + tkes*zl(i,k) ! Numerator   in Eq. 11 in BK13
          denom(i) = denom(i) + tkes         ! Denominator in Eq. 11 in BK13
        else
          cldarr(i) = one   ! Take note of columns containing cloud.
        endif
      enddo
    enddo

! Calculate the measure of PBL depth,  Eq. 11 in BK13 (Is this really PBL depth?)
    do i=1,nx
      if (denom(i) >  zero .and. numer(i) > zero) then
        l_inf(i) = min(0.1_kp * (numer(i)/denom(i)), 100.0_kp)
      else
        l_inf(i) = 100.0_kp
      endif
    enddo

!Calculate length scale outside of cloud, Eq. 10 in BK13 (Eq. 4.12 in Pete's dissertation)
    do k=1,nzm

      kb = k-1
      kc = k+1
      if (k == 1) then
        kb = 1
        kc = 2
        thedz(:) = adzi(:,kc)
      elseif (k == nzm) then
        kb = nzm-1
        kc = nzm
        thedz(:) = adzi(:,k)
      else
        thedz(:) = adzi(:,kc) + adzi(:,k) !  = (z(k+1)-z(k-1))
      endif

      do i=1,nx

!  vars module variable bet (=ggr/tv0) ; grid module variable  adzi

        betdz = bet(i,k) / thedz(i)

        tkes = sqrt(tke(i,k))

! Compute local Brunt-Vaisalla frequency

        wrk = qcl(i,k) + qci(i,k)
        if (wrk > zero) then            ! If in the cloud

! Find the in-cloud Brunt-Vaisalla frequency

           omn = qcl(i,k) / (wrk+1.0e-20_kp) ! Ratio of liquid water to total water

! Latent heat of phase transformation based on relative water phase content
! fac_cond = lcond/cp, fac_fus = lfus/cp

           lstarn = fac_cond + (one-omn)*fac_fus

! Derivative of saturation mixing ratio over water/ice wrt temp. based on relative water phase content
           dqsat =      omn  * dtqsatw(tabs(i,k),prsl(i,k))            &
                 + (one-omn) * dtqsati(tabs(i,k),prsl(i,k))

! Saturation mixing ratio over water/ice wrt temp  based on relative water phase content

           qsatt =      omn  * qsatw(tabs(i,k),prsl(i,k))               &
                 + (one-omn) * qsati(tabs(i,k),prsl(i,k))

! liquid/ice moist static energy static energy divided by cp?

           bbb = (one + epsv*qsatt-wrk-qpl(i,k)-qpi(i,k)                &
               + 1.61_kp*tabs(i,k)*dqsat) / (one+lstarn*dqsat)

! Calculate Brunt-Vaisalla frequency using centered differences in the vertical

           brunt(i,k) = betdz*(bbb*(hl(i,kc)-hl(i,kb))                  &
                        + (bbb*lstarn - (one+lstarn*dqsat)*tabs(i,k))   &
                        * (total_water(i,kc)-total_water(i,kb))         &
                        + (bbb*fac_cond - (one+fac_cond*dqsat)*tabs(i,k))*(qpl(i,kc)-qpl(i,kb))  &
                        + (bbb*fac_sub  - (one+fac_sub*dqsat)*tabs(i,k))*(qpi(i,kc)-qpi(i,kb)) )

        else                       ! outside of cloud

! Find outside-of-cloud Brunt-Vaisalla frequency
! Only unsaturated air, rain and snow contribute to virt. pot. temp.
! liquid/ice moist static energy divided by cp?

           bbb = one + epsv*qv(i,k) - qpl(i,k) - qpi(i,k)
           brunt(i,k) = betdz*( bbb*(hl(i,kc)-hl(i,kb))                        &
                        + epsv*tabs(i,k)*(total_water(i,kc)-total_water(i,kb)) &
                        + (bbb*fac_cond-tabs(i,k))*(qpl(i,kc)-qpl(i,kb))       &
                        + (bbb*fac_sub -tabs(i,k))*(qpi(i,kc)-qpi(i,kb)) )
        endif

! Reduction of mixing length in the stable regions (where B.-V. freq. > 0) is required.
! Here we find regions of Brunt-Vaisalla freq. > 0 for later use.

        if (brunt(i,k) >= zero) then
          brunt2(i,k) = brunt(i,k)
        else
          brunt2(i,k) = zero
        endif

! Calculate turbulent length scale in the boundary layer.
! See Eq. 10 in BK13 (Eq. 4.12 in Pete's dissertation)

! Keep the length scale adequately small near the surface following Blackadar (1984)
! Note that this is not documented in BK13 and was added later for SP-CAM runs

!       if (k == 1) then
!         term = 600.*tkes
!         smixt(i,k) = term + (0.4*zl(i,k)-term)*exp(-zl(i,k)*0.01)
!       else

! tscale is the eddy turnover time scale in the boundary layer and is
! an empirically derived constant

          if (tkes > zero .and. l_inf(i) > zero) then
            wrk1 = one / (tscale*tkes*vonk*zl(i,k))
            wrk2 = one / (tscale*tkes*l_inf(i))
            wrk1 = wrk1 + wrk2 + pt01 * brunt2(i,k) / tke(i,k)
            wrk1 = sqrt(one / max(wrk1,1.0e-8_kp)) * (one/0.3_kp)
!           smixt(i,k) = min(max_eddy_length_scale, 2.8284*sqrt(wrk1)/0.3)
            smixt(i,k) = min(max_eddy_length_scale, wrk1)

!           smixt(i,k) = min(max_eddy_length_scale,(2.8284*sqrt(1./((1./(tscale*tkes*vonk*zl(i,k)))  &
!                  + (1./(tscale*tkes*l_inf(i)))+0.01*(brunt2(i,k)/tke(i,k)))))/0.3)
!         else
!           smixt(i,k) = zero
          endif

!         endif


      enddo
    enddo

! Now find the in-cloud turbulence length scale
! See Eq. 13 in BK13 (Eq. 4.18 in Pete's disseration)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Remove after coupling to subgrid PDF.
!wthv_sec = -300/ggr*brunt*tk
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! determine cubed convective velocity scale (conv_vel2) inside the cloud

!   call conv_scale()           ! inlining the relevant code

!   do i=1,nx
!     conv_vel2(i,1) = zero ! Convective velocity scale cubed
!   enddo
                                ! Integrate velocity scale in the vertical
!   do k=2,nzm
!     do i=1,nx
!       conv_vel2(i,k) = conv_vel2(i,k-1)                               &
!                      + 2.5*adzi(i,k)*bet(i,k)*wthv_sec(i,k)
!     enddo
!   enddo

    do i=1,nx

      if (cldarr(i) == 1) then ! If there's a cloud in this column 

        kl = 0
        ku = 0
        do k=2,nzm-3

! Look for the cloud base in this column  
! thresh (=0) is a  variable local to eddy_length(). Should be a module constant.
          wrk = qcl(i,k) + qci(i,k)
          if (wrk > qcmin) then
            if (kl == 0) then
               kl = k
            endif

! Look for the cloud top in this column
            if (qcl(i,k+1)+qci(i,k+1) <= qcmin) then
              ku = k
! conv_vel2 (Cubed convective velocity scale) is calculated in conv_scale()
! Use the value of conv_vel2 at the top of the cloud. 
!             conv_var = conv_vel2(i,k)** oneb3 
            endif
          endif

! Compute the mixing length scale for the cloud layer that we just found
!         if (kl > 0 .and. ku > 0 .and. ku-kl > 1) then
!         if (kl > 0 .and. ku > 0 .and. ku-kl > 0) then
          if (kl > 0 .and. ku >= kl) then
! The calculation below finds the integral in the Eq. 10 in BK13 for the current cloud
            conv_var = zero
            do kk=kl,ku
              conv_var = conv_var+ 2.5_kp*adzi(i,kk)*bet(i,kk)*wthv_sec(i,kk)
            enddo
            conv_var = conv_var ** oneb3

            if (conv_var > zero) then ! If convective vertical velocity scale > 0

              depth = (zl(i,ku)-zl(i,kl)) + adzl(i,kl)

              do kk=kl,ku
! in-cloud turbulence length scale, Eq. 13 in BK13 (Eq. 4.18)

!               wrk = conv_var/(depth*sqrt(tke(i,kk)))
!               wrk = wrk * wrk + pt01*brunt2(i,kk)/tke(i,kk)

                wrk = conv_var/(depth*depth*sqrt(tke(i,kk)))  &
                    + pt01*brunt2(i,kk)/tke(i,kk)

                smixt(i,kk) = min(max_eddy_length_scale, (one/0.3_kp)*sqrt(one/wrk))

              enddo

            endif ! If convective vertical velocity scale > 0
            kl = zero
            ku = zero
          endif ! if inside the cloud layer

        enddo   ! k=2,nzm-3
      endif     ! if in the cloudy column
    enddo       ! i=1,nx


  end subroutine eddy_length


  subroutine conv_scale()

! This subroutine calculates the cubed convective velocity scale needed
! for the definition of the length scale in clouds
! See Eq. 16 in BK13 (Eq. 4.21 in Pete's dissertation)

    integer i,  k

!!!!!!!!!
!! A bug in formulation of conv_vel
!  Obtain it by averaging conv_vel2 in the horizontal
!!!!!!!!!!

!   conv_vel(1)=zero      ! Horizontally averaged convective velocity scale cubed
    do i=1,nx
      conv_vel2(i,1) = zero ! Convective velocity scale cubed
    enddo
! Integrate velocity scale in the vertical
    do k=2,nzm
!     conv_vel(k)=conv_vel(k-1)
      do i=1,nx
!**********************************************************************
!Do not include grid-scale contribution to convective velocity scale in GCM applications
!       conv_vel(k)=conv_vel(k-1)+2.5*adzi(k)*bet(k)*(tvwle(k)+tvws(k))
!       conv_vel(k)=conv_vel(k)+2.5*adzi(i,k)*bet(i,k)*(tvws(k))
!Do not include grid-scale contribution to convective velocity scale in GCM applications
!       conv_vel2(i,k)=conv_vel2(i,k-1)+2.5*adzi(k)*bet(k)*(tvwle(k)+wthv_sec(i,k))
!**********************************************************************

        conv_vel2(i,k) = conv_vel2(i,k-1)                               &
                       + 2.5_kp*adzi(i,k)*bet(i,k)*wthv_sec(i,k)
      enddo
    enddo

  end subroutine conv_scale


  subroutine check_eddy()

! This subroutine checks eddy length values

    integer i, k, kb, ks, zend
    real    wrk
!   real zstart, zthresh, qthresh

! Temporary kludge for marine stratocumulus under very strong inversions at coarse resolution
! Placement until some explicity PBL top is put in
! Not used.
!   zthresh = 100.
!   qthresh = -6.0

    do k=1,nzm

      if (k == nzm) then
        kb = k
      else
        kb = k+1
      endif

      do i=1,nx

        wrk = 0.1_kp*adzl(i,k)
                                                            ! Minimum 0.1 of local dz
        smixt(i,k) = max(wrk, min(max_eddy_length_scale,smixt(i,k)))

! If chracteristic grid dimension in the horizontal< 1000m, set lengthscale to 
! be not larger that that.
!       if (sqrt(dx*dy) .le. 1000.) smixt(i,k)=min(sqrt(dx*dy),smixt(i,k))

        if (qcl(i,kb) == zero .and. qcl(i,k) > zero .and. brunt(i,k) > 1.0e-4_kp) then
!If just above the cloud top and atmosphere is stable, set to  0.1 of local dz
          smixt(i,k) = wrk
        endif

      enddo ! i
    enddo   ! k

  end subroutine check_eddy

  subroutine canuto()

! Subroutine impements an analytic expression for the third moment of SGS vertical velocity
! based on Canuto et at, 2001, JAS, 58, 1169-1172 (further referred to as C01)
! This allows to avoid having a prognostic equation for the third moment.
! Result is returned in a global variable w3 defined at the interface levels.

! Local variables
    integer i, k, kb, kc

    real bet2,   f0,     f1,     f2,  f3,    f4,   f5,  iso, isosqr,         &
         omega0, omega1, omega2, X0,  Y0,    X1,   Y1,  AA0, AA1, buoy_sgs2, &
                 wrk, wrk1,  wrk2, wrk3, avew
!        cond,   wrk, wrk1,  wrk2, wrk3, avew
!
! See Eq. 7 in C01 (B.7 in Pete's dissertation)
    real, parameter :: c=7.0_kp,    a0=0.52_kp/(c*c*(c-2.0_kp)), a1=0.87_kp/(c*c),         &
                       a2=0.5_kp/c, a3=0.6_kp/(c*(c-2.0_kp)), a4=2.4_kp/(3.0_kp*c+5.0_kp), &
                       a5=0.6_kp/(c*(3.0_kp*c+5.0_kp))
!Moorthi               a5=0.6_kp/(c*(3.0_kp+5.0_kp*c))

!   do k=1,nzm
    do k=2,nzm

      kb = k-1
      kc = k+1

!     if(k == 1) then
!       kb = 1
!       kc = 2
!       do i=1,nx
!         thedz(i)  = one / adzl(i,kc)
!         thedz2(i) = thedz(i)
!       enddo
!     elseif(k == nzm) then
      if(k == nzm) then
        kb = nzm-1
        kc = nzm
        do i=1,nx
          thedz(i)  = one / adzi(i,k)
          thedz2(i) = one / adzl(i,kb)
        enddo
      else
        do i=1,nx
          thedz(i)  = one / adzi(i,k)
          thedz2(i) = one / (adzl(i,k)+adzl(i,kb))
        enddo
      endif

      do i=1,nx

        iso       = half*(isotropy(i,k)+isotropy(i,kb))
        isosqr    = iso*iso ! Two-level average of "return-to-isotropy" time scale squared
        buoy_sgs2 = isosqr*half*(brunt(i,k)+brunt(i,kb))
        bet2      = half*(bet(i,k)+bet(i,kb))  !Two-level average of BV frequency squared


! Compute functions f0-f5, see Eq, 8 in C01 (B.8 in Pete's dissertation)


        avew = half*(w_sec(i,k)+w_sec(i,kb))

!<aab
! This is not a bug, but an algorithmical change.
! The line below calculates cond_w ,an estimate of the maximum allowed value of the third moment.
! It is used at the end of this subroutine to limit the value of w3.
! Here the second moment is interpolated from the layer centers to the interface, where w3 is defined.
! In the presence of strong vertical gradients of w2, the value interpolated to the interface can
! be as much as twice as as large (or as small) as the value on in layer center. When the skewness
! of W PDF is calculated in assumed_pdf(), the code there uses w2 on the layer center, and the value
! of w3 interpolated from the interfaces to the layer center. The errors introduced due to dual
! interpolation are amplified by exponentiation during the calculation of skewness
! and result in (ususally negative) values
! of skewness of W PDF that are too large ( < -10). The resulting PDF consists of two delta
! functions set far apart from each other on the W axis, and produces unrealistically
! high values of diagnosed SGS fluxes.

! The solution is to move clipping of the third moment from canuto() to right before the
! calculation of skewness in assumed_pdf()
!          cond_w = 1.2*sqrt(max(1.0e-20, 2.*avew*avew*avew))
!>aab
!

        wrk1 = bet2*iso
        wrk2 = thedz2(i)*wrk1*wrk1*iso
        wrk3 = thl_sec(i,kc) - thl_sec(i,kb)

        f0   = wrk2 * wrk1 * wthl_sec(i,k) * wrk3

        wrk  = wthl_sec(i,kc) - wthl_sec(i,kb)

        f1   = wrk2 * (wrk*wthl_sec(i,k) + half*avew*wrk3)

        wrk1 = bet2*isosqr
        f2   = thedz(i)*wrk1*wthl_sec(i,k)*(w_sec(i,k)-w_sec(i,kb))     &
             + (thedz2(i)+thedz2(i))*bet(i,k)*isosqr*wrk

        f3   = thedz2(i)*wrk1*wrk + thedz(i)*bet2*isosqr*(wthl_sec(i,k)*(tke(i,k)-tke(i,kb)))

        wrk1 = thedz(i)*iso*avew
        f4   = wrk1*(w_sec(i,k)-w_sec(i,kb) + tke(i,k)-tke(i,kb))

        f5   = wrk1*(w_sec(i,k)-w_sec(i,kb))


! Compute the "omega" terms, see Eq. 6 in C01 (B.6 in Pete's dissertation)

        omega0 = a4 / (one-a5*buoy_sgs2)
        omega1 = omega0 / (c+c)
        omega2 = omega1*f3+(5.0_kp/4.0_kp)*omega0*f4

! Compute the X0, Y0, X1, Y1 terms,  see Eq. 5 a-b in C01  (B.5 in Pete's dissertation)

        wrk1 = one / (one-(a1+a3)*buoy_sgs2)
        wrk2 = one / (one-a3*buoy_sgs2)
        X0   = wrk1 * (a2*buoy_sgs2*(one-a3*buoy_sgs2))
        Y0   = wrk2 * (two*a2*buoy_sgs2*X0)
        X1   = wrk1 * (a0*f0+a1*f1+a2*(one-a3*buoy_sgs2)*f2)
        Y1   = wrk2 * (two*a2*(buoy_sgs2*X1+(a0/a1)*f0+f1))

! Compute the A0, A1 terms,  see Eq. 5d in C01 (B.5 in Pete's dissertation)

        AA0 = omega0*X0 + omega1*Y0
        AA1 = omega0*X1 + omega1*Y1 + omega2

! Finally, we have the third moment of w, see Eq. 4c in C01 (B.4 in Pete's dissertation)
! cond_w is an estimate of third moment from second oment - If the third moment is larger
! than the estimate - limit w3.

!<aab
! Move clipping of w3 to assumed_pdf()
!       w3(i,k) = max(-cond_w, min(cond_w, (AA1-1.2*X1-1.5*f5)/(c-1.2*X0+AA0)))
        w3(i,k) = (AA1-1.2_kp*X1-1.5_kp*f5)/(c-1.2_kp*X0+AA0)
!>aab

! Implemetation of the C01 approach in this subroutine is nearly complete
! (the missing part are Eqs. 5c and 5e which are very simple)
! therefore it's easy to diagnose other third order moments obtained in C01 using this code. 

      enddo
    enddo
    do i=1,nx
      w3(i,1) = w3(i,2)
    enddo

  end subroutine canuto

  subroutine assumed_pdf()

! Compute SGS buoyancy flux, SGS cloud fraction, and SGS condensation
! using assumed analytic double-gaussian PDF for SGS vertical velocity,
! moisture, and  liquid/ice water static energy, based on the
! general approach of  Larson et al 2002, JAS, 59, 3519-3539,
! and Golaz et al 2002, JAS, 59, 3540-3551
! References in the comments in this code are given to
! the Appendix A of Pete Bogenschutz's dissertation.

! Local variables

    integer i,k,ku,kd
    real wrk, wrk1, wrk2, wrk3, wrk4, bastoeps, eps_ss1, eps_ss2, cond_w

!   bastoeps = basetemp / epsterm


! Initialize for statistics
    do k=1,nzm
      wqlsb(k) = zero
      wqisb(k) = zero
    enddo

    DO k=1,nzm

      kd = k
      ku = k + 1
!     if (k == nzm) ku = k

      DO i=1,nx

! Initialize cloud variables to zero
        diag_qn   = zero
        diag_frac = zero
        diag_ql   = zero
        diag_qi   = zero

        pval  = prsl(i,k)
        pfac  = pval * 1.0e-5_kp
        pkap  = pfac ** kapa

! Read in liquid/ice static energy, total water mixing ratio, 
! and vertical velocity to variables PDF needs
        thl_first = hl(i,k) + fac_cond*qpl(i,k) + fac_sub*qpi(i,k)
        qw_first  = total_water(i,k)
!       w_first   = half*(w(i,kd)+w(i,ku))
        w_first   = w(i,k)


! GET ALL INPUT VARIABLES ON THE SAME GRID
! Points to be computed with relation to thermo point
! Read in points that need to be averaged

        if (k < nzm) then
          w3var    = half*(w3(i,kd)+w3(i,ku))
          thlsec   = max(zero, half*(thl_sec(i,kd)+thl_sec(i,ku)) )
          qwsec    = max(zero, half*(qw_sec(i,kd)+qw_sec(i,ku)) )
          qwthlsec = half * (qwthl_sec(i,kd) + qwthl_sec(i,ku))
          wqwsec   = half * (wqw_sec(i,kd)   + wqw_sec(i,ku))
          wthlsec  = half * (wthl_sec(i,kd)  + wthl_sec(i,ku))
        else                ! at the model top assuming zeros
          w3var    = half*w3(i,k)
          thlsec   = max(zero, half*thl_sec(i,k))
          qwsec    = max(zero, half*qw_sec(i,k))
          qwthlsec = half * qwthl_sec(i,k)
          wqwsec   = half * wqw_sec(i,k)
          wthlsec  = half * wthl_sec(i,k)
        endif

!       w3var    = w3(i,k)
!       thlsec   = max(zero,thl_sec(i,k))
!       qwsec    = max(zero,qw_sec(i,k))
!       qwthlsec = qwthl_sec(i,k)
!       wqwsec   = wqw_sec(i,k)
!       wthlsec  = wthl_sec(i,k)

! Compute square roots of some variables so we don't have to do it again
        if (w_sec(i,k) > zero) then
          sqrtw2  = sqrt(w_sec(i,k))
        else
          sqrtw2  = zero
        endif
        if (thlsec > zero) then
          sqrtthl = sqrt(thlsec)
        else
          sqrtthl = zero
        endif
        if (qwsec > zero) then
          sqrtqt  = sqrt(qwsec)
        else
          sqrtqt  = zero
        endif


! Find parameters of the double Gaussian PDF of vertical velocity

! Skewness of vertical velocity
!       Skew_w = w3var / w_sec(i,k)**(3./2.)
!       Skew_w = w3var / (sqrtw2*sqrtw2*sqrtw2)     ! Moorthi

        IF (w_sec(i,k) <= w_tol_sqd) THEN ! If variance of w is too small then
                                          ! PDF is a sum of two delta functions
          Skew_w = zero
          w1_1   = w_first
          w1_2   = w_first
          w2_1   = zero
          w2_2   = zero
          aterm  = half
          onema  = half
        ELSE
!<aab
! Clip w3
          cond_w = 1.2_kp*sqrt2*max(w3_tol, sqrtw2*sqrtw2*sqrtw2)
          w3var  = max(-cond_w, min(cond_w, w3var))
!>aab

          Skew_w = w3var / (sqrtw2*sqrtw2*sqrtw2)     ! Moorthi
! Proportionality coefficients between widths of each vertical velocity 
! gaussian and the sqrt of the second moment of w
          w2_1 = 0.4_kp
          w2_2 = 0.4_kp

! Compute realtive weight of the first PDF "plume" 
! See Eq A4 in Pete's dissertaion -  Ensure 0.01 < a < 0.99

          wrk   = one - w2_1
          aterm = max(atmin,min(half*(one-Skew_w*sqrt(one/(4.0_kp*wrk*wrk*wrk+Skew_w*Skew_w))),atmax))
          onema = one - aterm

          sqrtw2t = sqrt(wrk)

! Eq. A.5-A.6
          wrk  =   sqrt(onema/aterm)
          w1_1 =   sqrtw2t * wrk
          w1_2 = - sqrtw2t / wrk

          w2_1 = w2_1 * w_sec(i,k)
          w2_2 = w2_2 * w_sec(i,k)

        ENDIF

!  Find parameters of the  PDF of liquid/ice static energy

        IF (thlsec <= thl_tol*thl_tol .or. abs(w1_2-w1_1) <= w_thresh) THEN
          thl1_1     = thl_first
          thl1_2     = thl_first
          thl2_1     = zero
          thl2_2     = zero
          sqrtthl2_1 = zero
          sqrtthl2_2 = zero
        ELSE

          corrtest1 = max(-one,min(one,wthlsec/(sqrtw2*sqrtthl)))

          thl1_1 = -corrtest1 / w1_2                 ! A.7
          thl1_2 = -corrtest1 / w1_1                 ! A.8

          wrk1   = thl1_1 * thl1_1
          wrk2   = thl1_2 * thl1_2
          wrk3   = three * (one - aterm*wrk1 - onema*wrk2)
          wrk4   = -skew_facw*Skew_w - aterm*wrk1*thl1_1 - onema*wrk2*thl1_2  ! testing - Moorthi
!         wrk4   = -skew_fact*Skew_w - aterm*wrk1*thl1_1 - onema*wrk2*thl1_2  ! testing - Moorthi
!         wrk4   =     - aterm*wrk1*thl1_1 - onema*wrk2*thl1_2
          wrk    = three * (thl1_2-thl1_1)
          if (wrk /= zero) then
            thl2_1 = thlsec * min(100.0_kp,max(zero,(thl1_2*wrk3-wrk4)/(aterm*wrk))) ! A.10
            thl2_2 = thlsec * min(100.0_kp,max(zero,(-thl1_1*wrk3+wrk4)/(onema*wrk))) ! A.11
          else
            thl2_1 = zero
            thl2_2 = zero
          endif
!
          thl1_1 = thl1_1*sqrtthl + thl_first
          thl1_2 = thl1_2*sqrtthl + thl_first

          sqrtthl2_1 = sqrt(thl2_1)
          sqrtthl2_2 = sqrt(thl2_2)

        ENDIF

!  FIND PARAMETERS FOR TOTAL WATER MIXING RATIO

        IF (qwsec <= rt_tol*rt_tol .or. abs(w1_2-w1_1) <= w_thresh) THEN
          qw1_1     = qw_first
          qw1_2     = qw_first
          qw2_1     = zero
          qw2_2     = zero
          sqrtqw2_1 = zero
          sqrtqw2_2 = zero
        ELSE

          corrtest2 = max(-one,min(one,wqwsec/(sqrtw2*sqrtqt)))

          qw1_1 = - corrtest2 / w1_2            ! A.7
          qw1_2 = - corrtest2 / w1_1            ! A.8

          tsign = abs(qw1_2-qw1_1)

!         Skew_qw = skew_facw*Skew_w

          IF (tsign > 0.4_kp) THEN
            Skew_qw = skew_facw*Skew_w
          ELSEIF (tsign <= 0.2_kp) THEN
            Skew_qw = zero
          ELSE
            Skew_qw = (skew_facw/0.2_kp) * Skew_w * (tsign-0.2_kp)
          ENDIF

          wrk1  = qw1_1 * qw1_1
          wrk2  = qw1_2 * qw1_2
          wrk3  = three * (one - aterm*wrk1 - onema*wrk2)
          wrk4  = Skew_qw - aterm*wrk1*qw1_1 - onema*wrk2*qw1_2
          wrk   = three * (qw1_2-qw1_1)

          if (wrk /= zero) then
            qw2_1 = qwsec * min(100.0_kp,max(zero,( qw1_2*wrk3-wrk4)/(aterm*wrk))) ! A.10
            qw2_2 = qwsec * min(100.0_kp,max(zero,(-qw1_1*wrk3+wrk4)/(onema*wrk))) ! A.11
          else
            qw2_1 = zero
            qw2_2 = zero
          endif
!
          qw1_1 = qw1_1*sqrtqt + qw_first
          qw1_2 = qw1_2*sqrtqt + qw_first

          sqrtqw2_1 = sqrt(qw2_1)
          sqrtqw2_2 = sqrt(qw2_2)

        ENDIF

!  CONVERT FROM TILDA VARIABLES TO "REAL" VARIABLES

        w1_1 = w1_1*sqrtw2 + w_first
        w1_2 = w1_2*sqrtw2 + w_first

!  FIND WITHIN-PLUME CORRELATIONS 

        testvar = aterm*sqrtqw2_1*sqrtthl2_1 + onema*sqrtqw2_2*sqrtthl2_2

        IF (testvar == zero) THEN
          r_qwthl_1 = zero
        ELSE
          r_qwthl_1 = max(-one,min(one,(qwthlsec-aterm*(qw1_1-qw_first)*(thl1_1-thl_first) &
                                                -onema*(qw1_2-qw_first)*(thl1_2-thl_first))/testvar)) ! A.12
        ENDIF

!  BEGIN TO COMPUTE CLOUD PROPERTY STATISTICS

!       wrk1  = gamaz(i,k) - fac_cond*qpl(i,k) - fac_sub*qpi(i,k)
!       Tl1_1 = thl1_1 - wrk1
!       Tl1_2 = thl1_2 - wrk1

        Tl1_1 = thl1_1 - gamaz(i,k)
        Tl1_2 = thl1_2 - gamaz(i,k)

! Now compute qs

! Partition based on temperature for the first plume

        IF (Tl1_1 >= tbgmax) THEN
          lstarn1  = lcond
          esval    = min(fpvsl(Tl1_1), pval)
          qs1      = eps * esval / (pval-0.378_kp*esval)
        ELSE IF (Tl1_1 <= tbgmin) THEN
          lstarn1  = lsub
          esval    = min(fpvsi(Tl1_1), pval)
          qs1      = epss * esval / (pval-0.378_kp*esval)
        ELSE
          om1      = max(zero, min(one, a_bg*(Tl1_1-tbgmin)))
          lstarn1  = lcond + (one-om1)*lfus
          esval    = min(fpvsl(Tl1_1), pval)
          esval2   = min(fpvsi(Tl1_1), pval)
          qs1      =      om1  * eps  * esval  / (pval-0.378_kp*esval)      &
                   + (one-om1) * epss * esval2 / (pval-0.378_kp*esval2)
        ENDIF

!       beta1 = (rgas/rv)*(lstarn1/(rgas*Tl1_1))*(lstarn1/(cp*Tl1_1))
!       beta1 = (lstarn1*lstarn1*onebrvcp) / (Tl1_1*Tl1_1)              ! A.18

        beta1 = lstarn1 / Tl1_1
        beta1 = beta1 * beta1 * onebrvcp


! Are the two plumes equal?  If so then set qs and beta
! in each column to each other to save computation
        IF (Tl1_1 == Tl1_2) THEN
          qs2   = qs1
          beta2 = beta1
        ELSE
          IF (Tl1_2 >= tbgmax) THEN
            lstarn2  = lcond
            esval    = min(fpvsl(Tl1_2), pval)
            qs2      = eps * esval / (pval-0.378_kp*esval)
          ELSE IF (Tl1_2 <= tbgmin) THEN
            lstarn2  = lsub
            esval    = min(fpvsi(Tl1_2), pval)
            qs2      = epss * esval / (pval-0.378_kp*esval)
          ELSE
            om2      = max(zero, min(one, a_bg*(Tl1_2-tbgmin)))
            lstarn2  = lcond + (one-om2)*lfus
            esval    = min(fpvsl(Tl1_2), pval)
            esval2   = min(fpvsi(Tl1_2), pval)
            qs2      =      om2  * eps  * esval  / (pval-0.378_kp*esval)    &
                     + (one-om2) * epss * esval2 / (pval-0.378_kp*esval2)
          ENDIF

!         beta2 = (rgas/rv)*(lstarn2/(rgas*Tl1_2))*(lstarn2/(cp*Tl1_2))   ! A.18
!         beta2 = (lstarn2*lstarn2*onebrvcp) / (Tl1_2*Tl1_2)              ! A.18

          beta2 = lstarn2 / Tl1_2
          beta2 = beta2 * beta2 * onebrvcp


        ENDIF

        qs1 = qs1 * rhc(i,k)
        qs2 = qs2 * rhc(i,k)

!  Now compute cloud stuff -  compute s term

        cqt1   = one / (one+beta1*qs1)                                    ! A.19
        wrk    = qs1 * (one+beta1*qw1_1) * cqt1
        s1     = qw1_1 - wrk                                              ! A.17
        cthl1  = cqt1*wrk*cpolv*beta1*pkap                                ! A.20

        wrk1   = cthl1 * cthl1
        wrk2   = cqt1  * cqt1
!       std_s1 = sqrt(max(zero,wrk1*thl2_1+wrk2*qw2_1-2.*cthl1*sqrtthl2_1*cqt1*sqrtqw2_1*r_qwthl_1))
        std_s1 = sqrt(max(zero, wrk1*thl2_1+wrk2*qw2_1                        &
                              - two*cthl1*sqrtthl2_1*cqt1*sqrtqw2_1*r_qwthl_1))

        qn1 = zero
        C1  = zero

        IF (std_s1 > zero) THEN
          wrk = s1 / (std_s1*sqrt2)
          C1 = max(zero, min(one, half*(one+erf(wrk))))                   ! A.15

          IF (C1 > zero) qn1 = s1*C1 + (std_s1*sqrtpii)*exp(-wrk*wrk)     ! A.16
        ELSEIF (s1 >= qcmin) THEN
          C1  = one
          qn1 = s1
        ENDIF

! now compute non-precipitating cloud condensate 

! If two plumes exactly equal, then just set many of these 
! variables to themselves to save on computation.
        IF (qw1_1 == qw1_2 .and. thl2_1 == thl2_2 .and. qs1 == qs2) THEN
          s2     = s1
          cthl2  = cthl1
          cqt2   = cqt1
          std_s2 = std_s1
          C2     = C1
          qn2    = qn1
        ELSE

          cqt2   = one / (one+beta2*qs2)
          wrk    = qs2 * (one+beta2*qw1_2) * cqt2
          s2     = qw1_2 - wrk
          cthl2  = wrk*cqt2*cpolv*beta2*pkap
          wrk1   = cthl2 * cthl2
          wrk2   = cqt2  * cqt2
!         std_s2 = sqrt(max(zero,wrk1*thl2_2+wrk2*qw2_2-2.*cthl2*sqrtthl2_2*cqt2*sqrtqw2_2*r_qwthl_1))
          std_s2 = sqrt(max(zero, wrk1*thl2_2+wrk2*qw2_2                        &
                                - two*cthl2*sqrtthl2_2*cqt2*sqrtqw2_2*r_qwthl_1))

          qn2 = zero
          C2  = zero

          IF (std_s2 > zero) THEN
            wrk = s2 / (std_s2*sqrt2)
            C2  = max(zero, min(one, half*(one+erf(wrk))))
            IF (C2 > zero) qn2 = s2*C2 + (std_s2*sqrtpii)*exp(-wrk*wrk)
          ELSEIF (s2 >= qcmin) THEN
            C2  = one
            qn2 = s2
          ENDIF

        ENDIF

! finally, compute the SGS cloud fraction
        diag_frac = aterm*C1 + onema*C2

        om1 = max(zero, min(one, (Tl1_1-tbgmin)*a_bg))
        om2 = max(zero, min(one, (Tl1_2-tbgmin)*a_bg))

        qn1 = min(qn1,qw1_1)
        qn2 = min(qn2,qw1_2)

        ql1 = qn1*om1
        ql2 = qn2*om2

        qi1 = qn1 - ql1
        qi2 = qn2 - ql2

        diag_qn = min(max(zero, aterm*qn1 + onema*qn2), total_water(i,k))
        diag_ql = min(max(zero, aterm*ql1 + onema*ql2), diag_qn)
        diag_qi = max(zero, diag_qn - diag_ql)

! Update temperature variable based on diagnosed cloud properties
        om1         = max(zero, min(one, (tabs(i,k)-tbgmin)*a_bg))
        lstarn1     = lcond + (one-om1)*lfus
        tabs(i,k) = hl(i,k) - gamaz(i,k) + fac_cond*(diag_ql+qpl(i,k)) &
                                         + fac_sub *(diag_qi+qpi(i,k)) &
                  + tkesbdiss(i,k) * (dtn/cp)      ! tke dissipative heating

! Update ncpl and ncpi Anning Cheng 03/11/2016
!       ncpl(i,k)    = diag_ql/max(qc(i,k),1.e-10)*ncpl(i,k)

! Update ncpl and ncpi Moorthi  12/12/2018
        if (ntlnc > 0) then         ! liquid and ice number concentrations predicted
          if (ncpl(i,k) > nmin) then
            ncpl(i,k) = diag_ql/max(qc(i,k),1.0e-10_kp)*ncpl(i,k)
          else
            ncpl(i,k) = max(diag_ql/(fourb3*pi*RL_cub*997.0_kp), nmin)
          endif
          if (ncpi(i,k) > nmin) then
            ncpi(i,k) = diag_qi/max(qi(i,k),1.0e-10_kp)*ncpi(i,k)
          else
            ncpi(i,k) = max(diag_qi/(fourb3*pi*RI_cub*500.0_kp), nmin)
          endif
        endif

! Update moisture fields
        qc(i,k)      = diag_ql
        qi(i,k)      = diag_qi
        qwv(i,k)     = max(zero, total_water(i,k) - diag_qn)
        cld_sgs(i,k) = diag_frac

! Compute the liquid water flux
        wqls = aterm * ((w1_1-w_first)*ql1) + onema * ((w1_2-w_first)*ql2)
        wqis = aterm * ((w1_1-w_first)*qi1) + onema * ((w1_2-w_first)*qi2)

! Compute statistics for the fluxes so we don't have to save these variables
        wqlsb(k) = wqlsb(k) + wqls
        wqisb(k) = wqisb(k) + wqis

! diagnostic buoyancy flux.  Includes effects from liquid water, ice
! condensate, liquid & ice precipitation
!       wrk = epsv * basetemp
        wrk = epsv * thv(i,k)

        bastoeps = onebeps * thv(i,k)

        if (k < nzm) then
          wthv_sec(i,k) = wthlsec + wrk*wqwsec                                    &
                        + (fac_cond-bastoeps)*wqls                                &
                        + (fac_sub-bastoeps) *wqis                                &
                        + ((lstarn1/cp)-thv(i,k))*half*(wqp_sec(i,kd)+wqp_sec(i,ku))
        else
          wthv_sec(i,k) = wthlsec + wrk*wqwsec                                    &
                        + (fac_cond-bastoeps)*wqls                                &
                        + (fac_sub-bastoeps) *wqis                                &
                        + ((lstarn1/cp)-thv(i,k))*half*wqp_sec(i,k)
        endif

!         wthv_sec(i,k) = wthlsec + wrk*wqwsec                                    &
!                       + (fac_cond-bastoeps)*wqls                                &
!                       + (fac_sub-bastoeps)*wqis                                 &
!                       + ((lstarn1/cp)-basetemp)*half*(wqp_sec(i,kd)+wqp_sec(i,ku))

      ENDDO
    ENDDO


  end subroutine assumed_pdf


! Saturation vapor pressure and mixing ratio subroutines
! Based on Flatau et al (1992), J. App. Met., 31, 1507-1513
! Code by Marat Khairoutdinov


  real function esatw(t)
    real t	! temperature (K)
    real a0,a1,a2,a3,a4,a5,a6,a7,a8
    data a0,a1,a2,a3,a4,a5,a6,a7,a8 /                       &
         6.11239921,       0.443987641,     0.142986287e-1, &
         0.264847430e-3,   0.302950461e-5,  0.206739458e-7, &
         0.640689451e-10, -0.952447341e-13,-0.976195544e-15/
    real dt
    dt    = max(-80.,t-273.16)
    esatw = a0 + dt*(a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt)))))))
  end function esatw

  real function qsatw(t,p)
!    implicit none
    real t	! temperature (K)
    real p	! pressure    (Pa)
    real esat
!   esat  = fpvs(t)
    esat  = fpvsl(t)
    qsatw = 0.622 * esat/max(esat,p-0.378*esat)
!   esat  = esatw(t)
!   qsatw = 0.622 * esat/max(esat,p-esat)
  end function qsatw


  real function esati(t)
    real t	! temperature (K)
    real a0,a1,a2,a3,a4,a5,a6,a7,a8
    data a0,a1,a2,a3,a4,a5,a6,a7,a8 /                     &
         6.11147274,     0.503160820,     0.188439774e-1, &
         0.420895665e-3, 0.615021634e-5,  0.602588177e-7, &
         0.385852041e-9, 0.146898966e-11, 0.252751365e-14/
    real dt
!    real esatw
    if(t > 273.15) then
       esati = esatw(t)
    else if(t.gt.185.) then
       dt    = t-273.16
       esati = a0 + dt*(a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt)))))))
    else   ! use some additional interpolation below 184K
       dt    = max(-100.,t-273.16)
       esati = 0.00763685 + dt*(0.000151069+dt*7.48215e-07)
    endif
  end function esati

  real function qsati(t,p)
    real t	! temperature (K)
    real p	! pressure    (Pa)
    real esat !,esati
!   esat  = fpvs(t)
    esat  = fpvsi(t)
    qsati = 0.622 * esat/max(esat,p-0.378*esat)
!   esat  = esati(t)
!   qsati = 0.622 * esat/max(esat,p-esat)
  end function qsati

  real function dtesatw(t)
    real t	! temperature (K)
    real a0,a1,a2,a3,a4,a5,a6,a7,a8
    data a0,a1,a2,a3,a4,a5,a6,a7,a8 /                        &
         0.443956472,      0.285976452e-1,   0.794747212e-3, &
         0.121167162e-4,   0.103167413e-6,   0.385208005e-9, &
        -0.604119582e-12, -0.792933209e-14, -0.599634321e-17/
    real dt
    dt      = max(-80.,t-273.16)
    dtesatw = a0 + dt* (a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt)))))))
  end function dtesatw

  real function dtqsatw(t,p)
    real t	! temperature (K)
    real p	! pressure    (Pa)
!    real dtesatw
    dtqsatw = 100.0*0.622*dtesatw(t)/p
  end function dtqsatw

  real function dtesati(t)
    real t	! temperature (K)
    real a0,a1,a2,a3,a4,a5,a6,a7,a8
    data a0,a1,a2,a3,a4,a5,a6,a7,a8 /                      &
         0.503223089,     0.377174432e-1,  0.126710138e-2, &
         0.249065913e-4,  0.312668753e-6,  0.255653718e-8, &
         0.132073448e-10, 0.390204672e-13, 0.497275778e-16/
    real dt
!    real dtesatw
    if(t > 273.15) then
       dtesati = dtesatw(t)
    else if(t > 185.) then
       dt      = t-273.16
       dtesati = a0 + dt*(a1+dt*(a2+dt*(a3+dt*(a4+dt*(a5+dt*(a6+dt*(a7+a8*dt)))))))
    else  ! use additional interpolation below 185K
       dt      = max(-100.,t-273.16)
       dtesati = 0.0013186 + dt*(2.60269e-05+dt*1.28676e-07)
    endif
  end function dtesati


  real function dtqsati(t,p)
    real t	! temperature (K)
    real p	! pressure    (Pa)
!    real dtesati
    dtqsati = 100.0*0.622*dtesati(t)/p
  end function dtqsati

end subroutine shoc_work

end module shoc

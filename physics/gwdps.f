!> \file gwdps.f
!! This file is the  parameterization of orographic gravity wave
!! drag and mountain blocking.

!> This module contains the CCPP-compliant orographic gravity wave dray scheme.
      module gwdps

      contains

!> \section arg_table_gwdps_init Argument Table
!!
      subroutine gwdps_init()
      end subroutine gwdps_init

!> \defgroup gfs_gwdps GFS Orographic Gravity Wave Drag and Mountain Blocking Scheme Module
!! \brief This subroutine includes orographic gravity wave drag and mountain
!! blocking.
!!
!> The time tendencies of zonal and meridional wind are altered to
!! include the effect of mountain induced gravity wave drag from
!! subgrid scale orography including convective breaking, shear
!! breaking and the presence of critical levels.
!!
!! \section arg_table_gwdps_run Argument Table
!! \htmlinclude gwdps_run.html
!!
!> \section gen_gwdps GFS Orographic GWD Scheme General Algorithm
!! -# Calculate subgrid mountain blocking
!! -# Calculate orographic wave drag
!!
!! The NWP model gravity wave drag (GWD) scheme in the GFS has two
!! main components: how the surface stress is computed, and then how
!! that stress is distributed over a vertical column where it may
!! interact with the models momentum.  Each of these depends on the
!! large scale environmental atmospheric state and assumptions about
!! the sub-grid scale processes.  In Alpert GWD (1987) based on linear,
!! two-dimensional non-rotating, stably stratified flow over a mountain ridge,
!! sub-grid scale gravity wave motions are assumed which propagate away
!! from the mountain.  Described in Alpert (1987), the flux measured over
!! a "low level" vertically averaged layer, in the atmosphere defines a base
!! level flux.  "Low level" was taken to be the first 1/3 of the troposphere
!! in the 1987 implementation.  This choice was meant to encompass a thick
!! low layer for vertical averages of the environmental (large scale) flow
!! quantities.  The vertical momentum flux or gravity wave stress in a
!! grid box due to a single mountain is given as in Pierrehumbert(1986) 
!! \cite pierrehumbert_1986 :
!! \f[
!! \tau =  \frac {\rho \: U^{3}\: G(F_{r})} {\Delta X \; N } 
!! \f]
!! where \f$ \Delta X \f$ is a grid increment, N is the Brunt Viasala frequency
!! \f[
!!  N(\sigma) = \frac{-g \: \sigma \:
!!  \frac{\partial\Theta}{\partial\sigma}}{\Theta \:R \:T} 
!! \f]
!! The environmental variables are calculated from a mass weighted vertical
!! average over a base layer.  \f$G(F_{r})\f$ is a monotonically increasing
!! function of Froude number :
!! \f[
!!  F_{r} = \frac{N h^{'}}{U} 
!! \f]
!! where \f$U\f$ is the wind speed calculated as a mass weighted vertical average in
!! the base layer, and  \f$h^{'}\f$, is the vertical displacement caused by the orography
!! variance.  An effective mountain length for the gravity wave processes:
!! \f[
!!  l^{*} =  \frac{\Delta X}{m} 
!! \f]
!! where \f$m\f$ is the number of mountains in a grid box, can then
!! be defined to obtain the form of the base level stress
!! \f[
!!  \tau =  \frac {\rho \: U^{3} \: G(F_{r})} {N \;l^{*}} 
!! \f]
!! giving the stress induced from the surface in a model grid box.
!! Pierrehumbert(1986) \cite pierrehumbert_1986 gives the form 
!! for the function \f$G(F_{r})\f$ as
!! \f[
!!  G(F_{r}) = \bar{G}\frac{F^{2}_{r}}{F^{2}_{r}\: + \:a^{2}} 
!! \f]
!! Where \f$ \bar{G}  \f$  is an order unity non-dimensional saturation
!! flux set to 1  and \f$a\f$ is a function of the mountain aspect ratio also
!!set to 1 in the 1987 implementation of the GFS GWD.  Typical values of
!! U=10m/s, N=0.01 1/s, l*=100km, and a=1, gives a flux of 1 Pascal and
!! if this flux is made to go to zero linearly with height then the
!! decelerations would be about 10/m/s/day which is consistent with
!! observations in PH.
!!
!! In Kim, Moorthi, Alpert's (1998, 2001) GWD currently in GFS operations,
!! the GWD scheme has the same physical basis as in Alpert (1987) with the addition
!! of enhancement factors for the amplitude, G, and mountain shape details
!! in G(Fr) to account for effects from the mountain blocking.  A factor,
!! E m’, is an enhancement factor on the stress in the Alpert '87 scheme.
!!  The E ranges from no enhancement to an upper limit of 3, E=E(OA)[1-3],
!!  and is a function of OA, the Orographic Asymmetry defined in Kim and Arakawa (1995) 
!! \cite kim_and_arakawa_1995 as
!! Orographic Asymmetry (OA): 
!! \f[
!!  OA=\frac{ \bar{x} \; - \;
!!  \sum\limits_{j=1}^{N_{b}} x_{j} \; n_{j} }{\sigma_{x}} 
!! \f]
!! where \f$N_{b}\f$ is the total number of bottom blocks in the mountain barrier,
!! \f$ \sigma_{x} \f$ is the standard deviation of the horizontal distance defined by
!!\f[
!!  \sigma_{x} = \sqrt{ \frac{\sum\limits_{j=1}^{N_{b}}
!! \; (x_{j} \; - \; \bar{x} )^2}{N_{x}} } 
!!\f]
!! where \f$N_{x}\f$ is the number of grid intervals for the large scale domain being
!! considered. So the term, E(OA)m’/  \f$ \Delta X \f$ in Kim's scheme represents
!! a multiplier on G shown in Alpert's eq (1), where m’ is the number of mountains
!! in a sub-grid scale box. Kim increased the complexity of m’ making it a
!! function of the fractional area of the sub-grid mountain and the asymmetry
!! and convexity statistics which are found from running a gravity wave
!!  model for a large number of cases:
!! \f[
!!  m^{'} = C_{m} \Delta X \left[  \frac{1 \; + \;
!!  \sum\limits_{x} L_{h} }{\Delta X}  \right]^{OA+1} 
!! \f]
!! Where, according to Kim,  \f$ \sum \frac{L_{h}}{\Delta X} \f$  is
!! the fractional area covered by the subgrid-scale orography higher than
!! a critical height  \f$ h_{c} = Fr_{c} U_{0}/N_{0} \f$ , over the
!! "low level" vertically averaged layer, for a grid box with the interval
!! \f$ \Delta X \f$.  Each \f$ L_{n}\f$  is the width of a segment of
!! orography intersection at the critical height:
!! \f[
!! Fr_{0} = \frac{N_{0} \; h^{'}}{U_{0}} 
!! \f]
!! \f[
!! G^{'}(OC,Fr_{0}) = \frac{Fr_{0}^{2}}{Fr_{0}^{2} \; + \; a^{2}} 
!! \f]
!! \f[
!! a^{2} = \frac{C_{G}}{OC} 
!! \f]
!! \f[
!!  E(OA, Fr_{0}) = (OA \; + \; 2)^{\delta} 
!! \f] 
!! and 
!! \f[  
!! \delta \; = \; \frac{C_{E} \; Fr_{0}}{Fr_{c}}  
!! \f]
!! where \f$ Fr_{c} \f$ is as in Alpert.
!!
!! This represents a closed scheme, somewhat empirical adjustments
!! to the original scheme to calculate the surface stress.
!!
!! Momentum is deposited by the sub-grid scale gravity waves break due
!! to the presence of convective mixing assumed to occur when the
!! minimum Richardson number:
!!
!! Orographic Convexity (OC) = 
!! \f[  
!! OC=\frac{ \sum\limits_{j=1}^{N_{x}}
!!  \; (h_{j} \; - \; \bar{h})^4 }{N_{x} \;\sigma_{h}^4} 
!! \f]
!!  and where 
!!  \f[
!!  \sigma_{h} = \sqrt{ \frac{\sum\limits_{j=1}^{N_{x}}
!!  \; (h_{j} \; - \; \bar{h} )^2}{N_{x}} } 
!!  \f]
!! This represents a closed scheme, somewhat empirical adjustments
!! to the original scheme to calculate the surface stress.
!!
!! Momentum is deposited by the sub-grid scale gravity waves break due
!!  to the presence of convective mixing assumed to occur when
!!  the minimum Richardson number:
!! \f[
!!  Ri_{m} = \frac{Ri(1 \; - \; Fr)}{(1 \; + \; \sqrt{Ri}Fr)^2} 
!! \f]
!! Is less than 1/4  Or if critical layers are encountered in a layer
!! the the momentum flux will vanish.  The critical layer is defined
!! when the base layer wind becomes perpendicular to the environmental
!!  wind.  Otherwise, wave breaking occurs at a level where the amplification
!!  of the wave causes the local Froude number or similarly a truncated
!!  (first term of the) Scorer parameter, to be reduced below a critical
!!  value by the saturation hypothesis (Lindzen).  This is done through
!!  eq 1 which can be written as
!! \f[
!! \tau = \rho U N k h^{'2} 
!! \f]
!! For small Froude number this is discretized in the vertical so at each
!!  level the stress is reduced by ratio of the Froude or truncated Scorer
!!  parameter, \f$ \frac{U^{2}}{N^{2}} = \frac{N \tau_{l-1}}{\rho U^{3} k} \f$ ,
!!  where the stress is from the layer below beginning with that found near
!!  the surface.   The respective change in momentum is applied in
!! that layer building up from below.
!!
!! An amplitude factor is part of the calibration of this scheme which is
!! a function of the model resolution and the vertical diffusion.  This
!! is because the vertical diffusion and the GWD account encompass
!! similar physical processes. Thus, one needs to run the model over
!! and over for various amplitude factors for GWD and vertical diffusion.
!!
!! In addition, there is also mountain blocking from lift and frictional
!!  forces.  Improved integration between how the GWD is calculated and
!! the mountain blocking of wind flow around sub-grid scale orography
!! is underway at NCEP.  The GFS already has convectively forced GWD
!!  an independent process. 
!!
!> \section det_gwdps GFS Orographic GWD Scheme Detailed Algorithm
!> @{
      subroutine gwdps_run(                                             &
     &           IM,IX,KM,A,B,C,U1,V1,T1,Q1,KPBL,                       &
     &           PRSI,DEL,PRSL,PRSLK,PHII, PHIL,DELTIM,KDT,             &
     &           HPRIME,OC,OA4,CLX4,THETA,SIGMA,GAMMA,ELVMAX,           &
     &           DUSFC,DVSFC,G, CP, RD, RV, IMX,                        &
     &           nmtvr, cdmbgwd, me, lprnt, ipr, rdxzb, errmsg, errflg)
!
!   ********************************************************************
! ----->  I M P L E M E N T A T I O N    V E R S I O N   <----------
!
!          --- Not in this code --  History of GWDP at NCEP----
!              ----------------     -----------------------
!  VERSION 3  MODIFIED FOR GRAVITY WAVES, LOCATION: .FR30(V3GWD)  *J*
!---       3.1 INCLUDES VARIABLE SATURATION FLUX PROFILE CF ISIGST
!---       3.G INCLUDES PS COMBINED W/ PH (GLAS AND GFDL)
!-----         ALSO INCLUDED IS RI  SMOOTH OVER A THICK LOWER LAYER
!-----         ALSO INCLUDED IS DECREASE IN DE-ACC AT TOP BY 1/2
!-----     THE NMC GWD INCORPORATING BOTH GLAS(P&S) AND GFDL(MIGWD)
!-----        MOUNTAIN INDUCED GRAVITY WAVE DRAG
!-----    CODE FROM .FR30(V3MONNX) FOR MONIN3
!-----        THIS VERSION (06 MAR 1987)
!-----        THIS VERSION (26 APR 1987)    3.G
!-----        THIS VERSION (01 MAY 1987)    3.9
!-----    CHANGE TO FORTRAN 77 (FEB 1989)     --- HANN-MING HENRY JUANG
!-----    20070601 ELVMAX bug fix (*j*)
!
!   VERSION 4
!                ----- This code -----
!
!-----   MODIFIED TO IMPLEMENT THE ENHANCED LOW TROPOSPHERIC GRAVITY
!-----   WAVE DRAG DEVELOPED BY KIM AND ARAKAWA(JAS, 1995).
!        Orographic Std Dev (hprime), Convexity (OC), Asymmetry (OA4)
!        and Lx (CLX4) are input topographic statistics needed.
!
!-----   PROGRAMMED AND DEBUGGED BY HONG, ALPERT AND KIM --- JAN 1996.
!-----   debugged again - moorthi and iredell --- may 1998.
!-----
!       Further Cleanup, optimization and modification
!                                       - S. Moorthi May 98, March 99.
!-----   modified for usgs orography data (ncep office note 424)
!        and with several bugs fixed  - moorthi and hong --- july 1999.
!
!-----   Modified & implemented into NRL NOGAPS
!                                       - Young-Joon Kim, July 2000
!-----
!   VERSION lm MB  (6): oz fix 8/2003
!                ----- This code -----
!
!------   Changed to include the Lott and Miller Mtn Blocking
!         with some modifications by (*j*)  4/02
!        From a Principal Coordinate calculation using the
!        Hi Res 8 minute orography, the Angle of the
!        mtn with that to the East (x) axis is THETA, the slope
!        parameter SIGMA. The anisotropy is in GAMMA - all  are input
!        topographic statistics needed.  These are calculated off-line
!        as a function of model resolution in the fortran code ml01rg2.f,
!        with script mlb2.sh.   (*j*)
!-----   gwdps_mb.f version (following lmi) elvmax < hncrit (*j*)
!        MB3a expt to enhance elvmax mtn hgt see sigfac & hncrit
!        gwdps_GWDFIX_v6.f FIXGWD GF6.0 20070608 sigfac=4.
!-----
!----------------------------------------------------------------------C
!    USE
!        ROUTINE IS CALLED FROM GBPHYS  (AFTER CALL TO MONNIN)
!
!    PURPOSE
!        USING THE GWD PARAMETERIZATIONS OF PS-GLAS AND PH-
!        GFDL TECHNIQUE.  THE TIME TENDENCIES OF U V
!        ARE ALTERED TO INCLUDE THE EFFECT OF MOUNTAIN INDUCED
!        GRAVITY WAVE DRAG FROM SUB-GRID SCALE OROGRAPHY INCLUDING
!        CONVECTIVE BREAKING, SHEAR BREAKING AND THE PRESENCE OF
!        CRITICAL LEVELS
!
!  INPUT
!        A(IX,KM)  NON-LIN TENDENCY FOR V WIND COMPONENT
!        B(IX,KM)  NON-LIN TENDENCY FOR U WIND COMPONENT
!        C(IX,KM)  NON-LIN TENDENCY FOR TEMPERATURE
!        U1(IX,KM) ZONAL WIND M/SEC  AT T0-DT
!        V1(IX,KM) MERIDIONAL WIND M/SEC AT T0-DT
!        T1(IX,KM) TEMPERATURE DEG K AT T0-DT
!        Q1(IX,KM) SPECIFIC HUMIDITY AT T0-DT
!
!        DELTIM  TIME STEP    SECS
!        SI(N)   P/PSFC AT BASE OF LAYER N
!        SL(N)   P/PSFC AT MIDDLE OF LAYER N
!        DEL(N)  POSITIVE INCREMENT OF P/PSFC ACROSS LAYER N
!        KPBL(IM) is the index of the top layer of the PBL
!        ipr & lprnt for diagnostics
!
!  OUTPUT
!        A, B    AS AUGMENTED BY TENDENCY DUE TO GWDPS
!                OTHER INPUT VARIABLES UNMODIFIED.
!  revision log:
!    May 2013  J. Wang change cleff back to opn setting
!    Jan 2014  J. Wang merge Henry and Fangin's dissipation heat in gfs to nems
!
!
!   ********************************************************************
      USE MACHINE , ONLY : kind_phys
      implicit none
!
      ! Interface variables
      integer, intent(in) :: im, ix, km, imx, kdt, ipr, me
      integer, intent(in) :: KPBL(IM) ! Index for the PBL top layer!
      real(kind=kind_phys), intent(in) ::                               &
     &                     deltim, G, CP, RD, RV, cdmbgwd(4)
      real(kind=kind_phys), intent(inout) ::                            &
     &                     A(IX,KM), B(IX,KM), C(IX,KM)
      real(kind=kind_phys), intent(in) ::                               &
     &                     U1(IX,KM),   V1(IX,KM),     T1(IX,KM),       &
     &                     Q1(IX,KM),   PRSI(IX,KM+1), DEL(IX,KM),      &
     &                     PRSL(IX,KM), PRSLK(IX,KM),  PHIL(IX,KM),     &
     &                     PHII(IX,KM+1)
      real(kind=kind_phys), intent(in) ::                               &
     &                     OC(IM), OA4(IX,4), CLX4(IX,4), HPRIME(IM)
      real(kind=kind_phys), intent(inout) :: ELVMAX(IM)
      real(kind=kind_phys), intent(in) ::                               &
     &                     THETA(IM), SIGMA(IM), GAMMA(IM)
      real(kind=kind_phys), intent(out) :: DUSFC(IM), DVSFC(IM),        &
     &                     RDXZB(IX)
      integer, intent(in) :: nmtvr
      logical, intent(in) :: lprnt
      character(len=*), intent(out) :: errmsg
      integer,          intent(out) :: errflg
!
      ! Local variables
! for lm mtn blocking
      real(kind=kind_phys) wk(IM)
      real(kind=kind_phys) bnv2lm(IM,KM),PE(IM),EK(IM),ZBK(IM),UP(IM)
      real(kind=kind_phys) DB(IM,KM),ANG(IM,KM),UDS(IM,KM)
      real(kind=kind_phys) ZLEN, DBTMP, Rtrm, PHIANG, CDmb, DBIM, ZR
      real(kind=kind_phys) ENG0, ENG1
!
!     Some constants
!
      real(kind=kind_phys) pi, dw2min, rimin, ric, bnv2min, efmin
     &,                    efmax,hpmax,hpmin, rad_to_deg, deg_to_rad
      PARAMETER (PI=3.1415926535897931)
      PARAMETER (RAD_TO_DEG=180.0/PI, DEG_TO_RAD=PI/180.0)
      PARAMETER (DW2MIN=1., RIMIN=-100., RIC=0.25, BNV2MIN=1.0E-5)
!     PARAMETER (EFMIN=0.0, EFMAX=10.0, hpmax=200.0)
      PARAMETER (EFMIN=0.0, EFMAX=10.0, hpmax=2400.0, hpmin=1.0)
!
      real(kind=kind_phys) FRC,    CE,     CEOFRC, frmax, CG, GMAX
     &,                    VELEPS, FACTOP, RLOLEV, RDI
!     &,                    CRITAC, VELEPS, FACTOP, RLOLEV, RDI
      parameter (FRC=1.0, CE=0.8, CEOFRC=CE/FRC, frmax=100., CG=0.5)
      parameter (GMAX=1.0, VELEPS=1.0, FACTOP=0.5)
!      parameter (GMAX=1.0, CRITAC=5.0E-4, VELEPS=1.0, FACTOP=0.5)
      parameter (RLOLEV=50000.0)
!     parameter (RLOLEV=500.0)
!     parameter (RLOLEV=0.5)
!
       real(kind=kind_phys) dpmin,hminmt,hncrit,minwnd,sigfac
! --- for lm mtn blocking
!     parameter (cdmb = 1.0)     !< non-dim sub grid mtn drag Amp (*j*)
      parameter (hncrit=8000.)   !< Max value in meters for ELVMAX (*j*)
!  hncrit set to 8000m and sigfac added to enhance elvmax mtn hgt
      parameter (sigfac=4.0)     !< MB3a expt test for ELVMAX factor (*j*)
      parameter (hminmt=50.)     !< min mtn height (*j*)
      parameter (minwnd=0.1)     !< min wind component (*j*)

!     parameter (dpmin=00.0)     !< Minimum thickness of the reference layer
!!    parameter (dpmin=05.0)     !< Minimum thickness of the reference layer
!     parameter (dpmin=20.0)     !< Minimum thickness of the reference layer
                                 !< in centibars
      parameter (dpmin=5000.0)   !< Minimum thickness of the reference layer
                                 !< in Pa
!
      real(kind=kind_phys) FDIR
      integer mdir
      parameter(mdir=8, FDIR=mdir/(PI+PI))
      integer nwdir(mdir)
      data nwdir/6,7,5,8,2,3,1,4/
      save nwdir
!
      LOGICAL ICRILV(IM)
!
!----   MOUNTAIN INDUCED GRAVITY WAVE DRAG
!
      real(kind=kind_phys) TAUB(IM),  XN(IM),     YN(IM),    UBAR(IM)   &
     &,                    VBAR(IM),  ULOW(IM),   OA(IM),    CLX(IM)    &
     &,                    ROLL(IM),  ULOI(IM)                          &
     &,                    DTFAC(IM), XLINV(IM),  DELKS(IM)
!    &,                    DTFAC(IM), XLINV(IM),  DELKS(IM), DELKS1(IM)
!
      real(kind=kind_phys) BNV2(IM,KM),  TAUP(IM,KM+1), ri_n(IM,KM)     &
     &,                    TAUD(IM,KM),  RO(IM,KM),     VTK(IM,KM)      &
     &,                    VTJ(IM,KM),   SCOR(IM),      VELCO(IM,KM-1)  &
     &,                    bnv2bar(im)
!
!     real(kind=kind_phys) VELKO(KM-1)
      integer   kref(IM), kint(im), iwk(im), ipt(im)
! for lm mtn blocking
      integer   iwklm(im)
!      integer   kreflm(IM), iwklm(im)
      integer   idxzb(im), ktrial, klevm1
!
      real(kind=kind_phys) gor,    gocp,  fv,    gr2,  bnv,  fr         &
     &,                    brvf,   cleff, tem,   tem1,  tem2, temc, temv&
     &,                    wdir,   ti,    rdz,   dw2,   shr2, bvf2      &
     &,                    rdelks, efact, coefm, gfobnv, onebg          &
     &,                    scork,  rscor, hd,    fro,   rim,  sira      &
     &,                    dtaux,  dtauy, pkp1log, pklog                &
     &,                    cosang, sinang, cos2a, sin2a
!
      integer kmm1, kmm2, lcap, lcapp1, kbps, kbpsp1,kbpsm1             &
     &, kmps, idir, nwd, i, j, k, klcap, kp1, kmpbl, npt, npr, kmll
!    &, kmll,kmds,ihit,jhit
!
      ! Initialize CCPP error handling variables
      errmsg = ''
      errflg = 0
!
!     parameter (cdmb = 1.0)     ! non-dim sub grid mtn drag Amp (*j*)
! non-dim sub grid mtn drag Amp (*j*)
!     cdmb = 1.0/float(IMX/192)
!     cdmb = 192.0/float(IMX)
      cdmb = 4.0 * 192.0/float(IMX)
      if (cdmbgwd(1) >= 0.0) cdmb = cdmb * cdmbgwd(1)
!
      npr = 0
      DO I = 1, IM
         DUSFC(I) = 0.
         DVSFC(I) = 0.
      ENDDO
!
      DO K = 1, KM
        DO I = 1, IM
          DB(I,K)  = 0.
          ANG(I,K) = 0.
          UDS(I,K) = 0.
        ENDDO
      ENDDO
!
      RDI    = 1.0 / RD
      onebg  = 1.0 / g
      GOR    = G/RD
      GR2    = G*GOR
      GOCP   = G/CP
      FV     = RV/RD - 1
!
!     NCNT   = 0
      KMM1   = KM - 1
      KMM2   = KM - 2
      LCAP   = KM
      LCAPP1 = LCAP + 1
!
!
      IF ( NMTVR == 14) then 
! ----  for lm and gwd calculation points
        RDXZB(:)  = 0 
        ipt = 0
        npt = 0
        DO I = 1,IM
          IF (elvmax(i) > HMINMT .and. hprime(i) > hpmin)  then
             npt      = npt + 1
             ipt(npt) = i
             if (ipr == i) npr = npt
          ENDIF
        ENDDO
        IF (npt == 0) RETURN     ! No gwd/mb calculation done!
!
!       if (lprnt) print *,' npt=',npt,' npr=',npr,' ipr=',ipr,' im=',im
!    &,' ipt(npt)=',ipt(npt)
!
! --- iwklm is the level above the height of the of the mountain.
! --- idxzb is the level of the dividing streamline.
! INITIALIZE DIVIDING STREAMLINE (DS) CONTROL VECTOR
!
        do i=1,npt
          iwklm(i)  = 2
          IDXZB(i)  = 0 
!         kreflm(i) = 0
        enddo
!       if (lprnt) 
!    &  print *,' in gwdps_lm.f npt,IM,IX,IY,km,me=',npt,IM,IX,IY,km,me
!
!
!> --- Subgrid Mountain Blocking Section
!
!..............................
!..............................
!
!  (*j*)  11/03:  test upper limit on KMLL=km - 1
!      then do not need hncrit -- test with large hncrit first.
!       KMLL  = km / 2 ! maximum mtnlm height : # of vertical levels / 2

        KMLL = kmm1

! --- No mtn should be as high as KMLL (so we do not have to start at
! --- the top of the model but could do calc for all levels).
!
          DO I = 1, npt
            j = ipt(i)
            ELVMAX(J) = min (ELVMAX(J) + sigfac * hprime(j), hncrit)
          ENDDO
!
        DO K = 1,KMLL
          DO I = 1, npt
            j = ipt(i)
! --- interpolate to max mtn height for index, iwklm(I) wk[gz]
! --- ELVMAX is limited to hncrit because to hi res topo30 orog.
            pkp1log =  phil(j,k+1) * onebg
            pklog   =  phil(j,k)   * onebg
!!!-------     ELVMAX(J) = min (ELVMAX(J) + sigfac * hprime(j), hncrit)

            if (ELVMAX(j) <=  pkp1log .and. ELVMAX(j) >=  pklog) THEN

!     print *,' in gwdps_lm.f 1  =',k,ELVMAX(j),pklog,pkp1log,me
! ---        wk for diags but can be saved and reused.  
               wk(i)    = G * ELVMAX(j) / (phil(j,k+1) - phil(j,k))
               iwklm(I) =  MAX(iwklm(I), k+1) 
!     print *,' in gwdps_lm.f 2 npt=',npt,i,j,wk(i),iwklm(i),me
            endif
!
! ---        find at prsl levels large scale environment variables
! ---        these cover all possible mtn max heights
            VTJ(I,K)  = T1(J,K)  * (1.0+FV*Q1(J,K))
            VTK(I,K)  = VTJ(I,K) / PRSLK(J,K)
            RO(I,K)   = RDI * PRSL(J,K) / VTJ(I,K) ! DENSITY Kg/M**3
          ENDDO
        ENDDO
!
! testing for highest model level of mountain top
!
!         ihit = 2
!         jhit = 0
!        do i = 1, npt
!        j=ipt(i)
!          if ( iwklm(i) .gt. ihit ) then
!            ihit = iwklm(i)
!            jhit = j
!          endif
!        enddo
!     print *, ' mb: kdt,max(iwklm),jhit,phil,me=',
!    &          kdt,ihit,jhit,phil(jhit,ihit),me
         
        klevm1 = KMLL - 1
        DO K = 1, klevm1
          kp1 = k + 1
          DO I = 1, npt
           j   = ipt(i)
           RDZ = g / (phil(j,kp1) - phil(j,k))
! ---                               Brunt-Vaisala Frequency
!> - Compute Brunt-Vaisala Frequency \f$N\f$.
            BNV2LM(I,K) = (G+G) * RDZ * (VTK(I,Kp1) - VTK(I,K))
     &                                / (VTK(I,Kp1) + VTK(I,K))
            bnv2lm(i,k) = max( bnv2lm(i,k), bnv2min )
          ENDDO
        ENDDO
!    print *,' in gwdps_lm.f 3 npt=',npt,j,RDZ,me
!
        DO I = 1, npt
          J   = ipt(i)
          DELKS(I)   = 1.0 / (PRSI(J,1) - PRSI(J,iwklm(i)))
!         DELKS1(I)  = 1.0 / (PRSI(J,1) - PRSL(J,iwklm(i)))
          UBAR (I)   = 0.0
          VBAR (I)   = 0.0
          ROLL (I)   = 0.0
          PE   (I)   = 0.0
          EK   (I)   = 0.0
          BNV2bar(I) = (PRSI(J,1)-PRSL(J,1)) * DELKS(I) * BNV2LM(I,1)
        ENDDO

! --- find the dividing stream line height
! --- starting from the level above the max mtn downward
! --- iwklm(i) is the k-index of mtn elvmax elevation
!> - Find the dividing streamline height starting from the level above
!! the maximum mountain height and processing downward.
!       DO Ktrial = KMLL, 1, -1
!         DO I = 1, npt
!            IF ( Ktrial < iwklm(I) .and. kreflm(I) == 0 ) then
!               kreflm(I) = Ktrial
!            ENDIF
!         ENDDO
!       ENDDO
!     print *,' in gwdps_lm.f 4 npt=',npt,kreflm(npt),me
!
! --- in the layer kreflm(I) to 1 find PE (which needs N, ELVMAX)
! ---  make averages, guess dividing stream (DS) line layer.
! ---  This is not used in the first cut except for testing and
! --- is the vert ave of quantities from the surface to mtn top.
!   
        DO I = 1, npt
          DO K = 1, iwklm(i)-1
            J          = ipt(i)
            RDELKS     = DEL(J,K) * DELKS(I)
            UBAR(I)    = UBAR(I)  + RDELKS * U1(J,K) ! trial Mean U below
            VBAR(I)    = VBAR(I)  + RDELKS * V1(J,K) ! trial Mean V below
            ROLL(I)    = ROLL(I)  + RDELKS * RO(I,K) ! trial Mean RO below
            if (k < iwklm(I)-1) then
              RDELKS   = (PRSL(J,K)-PRSL(J,K+1)) * DELKS(I)
            else
              RDELKS   = (PRSL(J,K)-PRSI(J,K+1)) * DELKS(I)
            endif
            BNV2bar(I) = BNV2bar(I) + BNV2lm(I,K) * RDELKS
! --- these vert ave are for diags, testing and GWD to follow (*j*).
          ENDDO
        ENDDO
!     print *,' in gwdps_lm.f 5  =',i,kreflm(npt),BNV2bar(npt),me
!
! --- integrate to get PE in the trial layer.
! --- Need the first layer where PE>EK - as soon as
! --- IDXZB is not 0 we have a hit and Zb is found.
!
        DO I = 1, npt
          J = ipt(i)
          DO K = iwklm(I), 1, -1
            PHIANG   =  atan2(V1(J,K),U1(J,K))*RAD_TO_DEG
            ANG(I,K) = ( THETA(J) - PHIANG )
            if ( ANG(I,K) >  90. ) ANG(I,K) = ANG(I,K) - 180.
            if ( ANG(I,K) < -90. ) ANG(I,K) = ANG(I,K) + 180.
            ANG(I,K) = ANG(I,K) * DEG_TO_RAD
!
!> - Compute wind speed UDS
!!\f[
!!    UDS=\max(\sqrt{U1^2+V1^2},minwnd)
!!\f]
!! where \f$ minwnd=0.1 \f$, \f$U1\f$ and \f$V1\f$ are zonal and
!! meridional wind components of model layer wind.

            UDS(I,K) =
     &          MAX(SQRT(U1(J,K)*U1(J,K) + V1(J,K)*V1(J,K)), minwnd)
! --- Test to see if we found Zb previously
            IF (IDXZB(I) == 0 ) then
              PE(I) = PE(I) + BNV2lm(I,K) * (G*ELVMAX(J) - phil(J,K))
     &                                    * (PHII(J,K+1) - PHII(J,K))
     &                                    * (onebg*onebg)
! --- KE
! --- Wind projected on the line perpendicular to mtn range, U(Zb(K)).
! --- kenetic energy is at the layer Zb
! --- THETA ranges from -+90deg |_ to the mtn "largest topo variations"
              UP(I)  =  UDS(I,K) * cos(ANG(I,K))
              EK(I)  = 0.5 *  UP(I) * UP(I)

! --- Dividing Stream lime  is found when PE =exceeds EK.
              IF (PE(I) >= EK(I)) THEN
                 IDXZB(I) = K
                 RDXZB(J) = real(K,kind=kind_phys)
              ENDIF
! --- Then mtn blocked flow is between Zb=k(IDXZB(I)) and surface
!
!> - The dividing streamline height (idxzb), of a subgrid scale
!! obstable, is found by comparing the potential (PE) and kinetic
!! energies (EK) of the upstream large scale wind and subgrid scale air
!! parcel movements. the dividing streamline is found when
!! \f$PE\geq EK\f$. Mountain-blocked flow is defined to exist between
!! the surface and the dividing streamline height (\f$h_d\f$), which
!! can be found by solving an integral equation for \f$h_d\f$:
!!\f[
!! \frac{U^{2}(h_{d})}{2}=\int_{h_{d}}^{H} N^{2}(z)(H-z)dz
!!\f]
!! where \f$H\f$ is the maximum subgrid scale elevation within the grid
!! box of actual orography, \f$h\f$, obtained from the GTOPO30 dataset
!! from the U.S. Geological Survey.
            ENDIF
          ENDDO
        ENDDO
!
!     print *,' in gwdps_lm.f 6  =',phiang,THETA(ipt(npt)),me
!     print *,' in gwdps_lm.f 7  =',IDXZB(npt),PE(npt)
!
!     if (lprnt .and. npr .gt. 0) then
!       print *,' BNV2bar,BNV2lm=',bnv2bar(npr),BNV2lm(npr,1:klevm1)
!       print *,' npr,IDXZB,UDS=',npr,IDXZB(npr),UDS(npr,:)
!       print *,' PE,UP,EK=',PE(npr),UP(npr),EK(npr)
!     endif
!
        DO I = 1, npt
          J    = ipt(i)
! --- Calc if N constant in layers (Zb guess) - a diagnostic only.
          ZBK(I) = ELVMAX(J)
     &           - SQRT(UBAR(I)*UBAR(I) + VBAR(I)*VBAR(I))/BNV2bar(I)
        ENDDO
!
!     if (lprnt .and. npr .gt. 0) then
!       print *,' iwklm,ZBK=',iwklm(npr),ZBK(npr),IDXZB(npr)
!       print *,' Zb=',PHIL(ipr),IDXZB(npr))/G
!     print *,' in gwdps_lm.f 8 npt =',npt,ZBK(npt),UP(npt),me
!     endif
!
! --- The drag for mtn blocked flow
!
        DO I = 1, npt
          J = ipt(i)
          ZLEN = 0.
!      print *,' in gwdps_lm.f 9  =',i,j,IDXZB(i),me
          IF (IDXZB(I) > 0) then
            DO K = IDXZB(I), 1, -1
              IF (PHIL(J,IDXZB(I)) >  PHIL(J,K)) then

!> - Calculate \f$ZLEN\f$, which sums up a number of contributions of
!! elliptic obstables.
!!\f[
!!    ZLEN=\sqrt{[\frac{h_{d}-z}{z+h'}]}
!!\f]
!! where \f$z\f$ is the height, \f$h'\f$ is the orographic standard
!! deviation (HPRIME).
                ZLEN = SQRT( (PHIL(J,IDXZB(I)) - PHIL(J,K))
     &                     / (PHIL(J,K )       + G*hprime(J)) )
! --- lm eq 14:
!> - Calculate the drag coefficient to vary with the aspect ratio of
!! the obstable as seen by the incident flow (see eq.14 in Lott and
!! Miller (1997) \cite lott_and_miller_1997)
!!\f[
!! R=\frac{\cos^{2}\psi+\gamma\sin^{2}\psi}{\gamma\cos^{2}\psi+\sin^{2}\psi}
!!\f]
!! where \f$\psi\f$, which is derived from THETA, is the angle between
!! the incident flow direction and the normal ridge direcion.
!! \f$\gamma\f$ is the orographic anisotropy (GAMMA).

                cosang = cos(ang(i,k))
                sinang = sin(ang(i,k))
                cos2a  = cosang * cosang
                sin2a  = sinang * sinang
                tem    = cos2a + GAMMA(J)*sin2a
                                                       ! Here Rtrm is 1.0/R
                                                       ! --------------------
                if (abs(tem) > 1.e-06) then
                  Rtrm   = (gamma(J)*cos2a + sin2a) / tem
                elseif (tem > 0.0) then
                  Rtrm   = (gamma(J)*cos2a + sin2a) * 1.0e6
                else
                  Rtrm   = - (gamma(J)*cos2a + sin2a) * 1.0e6
                endif
                ZR =  MAX( 2.0 - Rtrm, 0. )

! --- (negitive of DB -- see sign at tendency)
!> - In each model layer below the dividing streamlines, a drag from
!! the blocked flow is exerted by the obstacle on the large scale flow.
!! The drag per unit area and per unit height is written (eq.15 in
!! Lott and Miller (1997) \cite lott_and_miller_1997):
!!\f[
!! D_{b}(z)=-C_{d}\max(2-\frac{1}{R},0)\rho\frac{\sigma}{2h'}ZLEN\max(\cos\psi,\gamma\sin\psi)\frac{UDS}{2}
!!\f]
!! where \f$C_{d}\f$ is a specified constant, \f$\sigma\f$ is the
!! orographic slope.

                DBTMP = 0.25 *  CDmb * ZR * sigma(J) *
     &                  MAX(cosANG, gamma(J)*sinANG) * ZLEN / hprime(J)
                DB(I,K) =  DBTMP * UDS(I,K)
!
!               if(lprnt .and. i .eq. npr) then
!                 print *,' in gwdps_lmi.f 10 npt=',npt,i,j,idxzb(i)
!    &,           DBTMP,R' ang=',ang(i,k),' gamma=',gamma(j),' K=',K
!                 print *,' in gwdps_lmi.f 11   K=',k,ZLEN,cos(ANG(I,K))
!                 print *,' in gwdps_lmi.f 12  DB=',DB(i,k),sin(ANG(I,K))
!               endif
              endif
            ENDDO
!         if(lprnt) print *,' @K=1,ZLEN,DBTMP=',K,ZLEN,DBTMP
          endif
        ENDDO
! 
!.............................
!.............................
! end  mtn blocking section
!
      ELSEIF ( NMTVR /= 14) then
! ----  for mb not present and  gwd (nmtvr .ne .14)
        ipt = 0
        npt = 0
        DO I = 1,IM
          IF ( hprime(i) > hpmin )  then
             npt      = npt + 1
             ipt(npt) = i
             if (ipr == i) npr = npt
          ENDIF
        ENDDO
        IF (npt == 0) RETURN     ! No gwd/mb calculation done!
!
!       if (lprnt) print *,' NPR=',npr,' npt=',npt,' IPR=',IPR
!      &,' ipt(npt)=',ipt(npt)
!
        do i=1,npt
          IDXZB(i) = 0
          RDXZB(i) = 0.
        enddo
      ENDIF
!
!.............................
!.............................
!
!>-# Orographic Gravity Wave Drag Section
      KMPBL  = km / 2 ! maximum pbl height : # of vertical levels / 2
!
!  Scale cleff between IM=384*2 and 192*2 for T126/T170 and T62
!
      if (imx > 0) then
!       cleff = 1.0E-5 * SQRT(FLOAT(IMX)/384.0) !  this is inverse of CLEFF!
!       cleff = 1.0E-5 * SQRT(FLOAT(IMX)/192.0) !  this is inverse of CLEFF!
!       cleff = 0.5E-5 * SQRT(FLOAT(IMX)/192.0) !  this is inverse of CLEFF!
!       cleff = 1.0E-5 * SQRT(FLOAT(IMX)/192)/float(IMX/192)
!       cleff = 1.0E-5 / SQRT(FLOAT(IMX)/192.0) !  this is inverse of CLEFF!
        cleff = 0.5E-5 / SQRT(FLOAT(IMX)/192.0) !  this is inverse of CLEFF!
! hmhj for ndsl
! jw        cleff = 0.1E-5 / SQRT(FLOAT(IMX)/192.0) !  this is inverse of CLEFF!
!       cleff = 2.0E-5 * SQRT(FLOAT(IMX)/192.0) !  this is inverse of CLEFF!
!       cleff = 2.5E-5 * SQRT(FLOAT(IMX)/192.0) !  this is inverse of CLEFF!
      endif
      if (cdmbgwd(2) >= 0.0) cleff = cleff * cdmbgwd(2)
!
      DO K = 1,KM
        DO I =1,npt
          J         = ipt(i)
          VTJ(I,K)  = T1(J,K)  * (1.0+FV*Q1(J,K))
          VTK(I,K)  = VTJ(I,K) / PRSLK(J,K)
          RO(I,K)   = RDI * PRSL(J,K) / VTJ(I,K) ! DENSITY TONS/M**3
          TAUP(I,K) = 0.0
        ENDDO
      ENDDO
      DO K = 1,KMM1
        kp1 = k + 1
        DO I =1,npt
          J         = ipt(i)
          TI        = 2.0 / (T1(J,K)+T1(J,Kp1))
          TEM       = TI  / (PRSL(J,K)-PRSL(J,Kp1))
          RDZ       = g   / (phil(j,kp1) - phil(j,k))
          TEM1      = U1(J,K) - U1(J,Kp1)
          TEM2      = V1(J,K) - V1(J,Kp1)
          DW2       = TEM1*TEM1 + TEM2*TEM2
          SHR2      = MAX(DW2,DW2MIN) * RDZ * RDZ
          BVF2      = G*(GOCP+RDZ*(VTJ(I,Kp1)-VTJ(I,K))) * TI
          ri_n(I,K) = MAX(BVF2/SHR2,RIMIN)   ! Richardson number
!                                              Brunt-Vaisala Frequency
!         TEM       = GR2 * (PRSL(J,K)+PRSL(J,K+1)) * TEM
!         BNV2(I,K) = TEM * (VTK(I,K+1)-VTK(I,K))/(VTK(I,K+1)+VTK(I,K))
          BNV2(I,K) = (G+G) * RDZ * (VTK(I,Kp1)-VTK(I,K))
     &                            / (VTK(I,Kp1)+VTK(I,K))
          bnv2(i,k) = max( bnv2(i,k), bnv2min )
        ENDDO
      ENDDO
!      print *,' in gwdps_lm.f GWD:14  =',npt,kmm1,bnv2(npt,kmm1)
!
!     Apply 3 point smoothing on BNV2
!
!     do k=1,km
!       do i=1,im
!         vtk(i,k) = bnv2(i,k)
!       enddo
!     enddo
!     do k=2,kmm1
!       do i=1,im
!         bnv2(i,k) = 0.25*(vtk(i,k-1)+vtk(i,k+1)) + 0.5*vtk(i,k)
!       enddo
!     enddo
!
!     Finding the first interface index above 50 hPa level
!
      do i=1,npt
        iwk(i) = 2
      enddo
      DO K=3,KMPBL
        DO I=1,npt
          j   = ipt(i)
          tem = (prsi(j,1) - prsi(j,k))
          if (tem < dpmin) iwk(i) = k
        enddo
      enddo
!
!> - Calculate the reference level index: kref=max(2,KPBL+1). where
!! KPBL is the index for the PBL top layer.
      KBPS = 1
      KMPS = KM
      DO I=1,npt
        J         = ipt(i)
        kref(I)   = MAX(IWK(I), KPBL(J)+1 ) ! reference level
        DELKS(I)  = 1.0 / (PRSI(J,1) - PRSI(J,kref(I)))
!       DELKS1(I) = 1.0 / (PRSI(J,1) - PRSL(J,kref(I)))
        UBAR (I)  = 0.0
        VBAR (I)  = 0.0
        ROLL (I)  = 0.0
        KBPS      = MAX(KBPS, kref(I))
        KMPS      = MIN(KMPS, kref(I))
!
        BNV2bar(I) = (PRSI(J,1)-PRSL(J,1)) * DELKS(I) * BNV2(I,1)
      ENDDO
!      print *,' in gwdps_lm.f GWD:15  =',KBPS,KMPS
      KBPSP1 = KBPS + 1
      KBPSM1 = KBPS - 1
      DO K = 1,KBPS
        DO I = 1,npt
          IF (K < kref(I)) THEN
            J          = ipt(i)
            RDELKS     = DEL(J,K) * DELKS(I)
            UBAR(I)    = UBAR(I)  + RDELKS * U1(J,K)   ! Mean U below kref
            VBAR(I)    = VBAR(I)  + RDELKS * V1(J,K)   ! Mean V below kref
!
            ROLL(I)    = ROLL(I)  + RDELKS * RO(I,K)   ! Mean RO below kref
            if (k < kref(i)-1) then
              RDELKS     = (PRSL(J,K)-PRSL(J,K+1)) * DELKS(I)
            else
              RDELKS     = (PRSL(J,K)-PRSI(J,K+1)) * DELKS(I)
            endif
            BNV2bar(I) = BNV2bar(I) + BNV2(I,K) * RDELKS
          ENDIF
        ENDDO
      ENDDO
!      print *,' in gwdps_lm.f GWD:15B =',bnv2bar(npt)
!
!     FIGURE OUT LOW-LEVEL HORIZONTAL WIND DIRECTION AND FIND 'OA'
!
!             NWD  1   2   3   4   5   6   7   8
!              WD  W   S  SW  NW   E   N  NE  SE
!
!> - Calculate low-level horizontal wind direction, the derived
!! orographic asymmetry parameter (OA), and the derived Lx (CLX).
      DO I = 1,npt
        J      = ipt(i)
        wdir   = atan2(UBAR(I),VBAR(I)) + pi
        idir   = mod(nint(fdir*wdir),mdir) + 1
        nwd    = nwdir(idir)
        OA(I)  = (1-2*INT( (NWD-1)/4 )) * OA4(J,MOD(NWD-1,4)+1)
        CLX(I) = CLX4(J,MOD(NWD-1,4)+1)
      ENDDO
!
!-----XN,YN            "LOW-LEVEL" WIND PROJECTIONS IN ZONAL
!                                    & MERIDIONAL DIRECTIONS
!-----ULOW             "LOW-LEVEL" WIND MAGNITUDE -        (= U)
!-----BNV2             BNV2 = N**2
!-----TAUB             BASE MOMENTUM FLUX
!-----= -(RO * U**3/(N*XL)*GF(FR) FOR N**2 > 0
!-----= 0.                        FOR N**2 < 0
!-----FR               FROUDE    =   N*HPRIME / U
!-----G                GMAX*FR**2/(FR**2+CG/OC)
!
!-----INITIALIZE SOME ARRAYS
!
      DO I = 1,npt
        XN(I)     = 0.0
        YN(I)     = 0.0
        TAUB (I)  = 0.0
        ULOW (I)  = 0.0
        DTFAC(I)  = 1.0
        ICRILV(I) = .FALSE. ! INITIALIZE CRITICAL LEVEL CONTROL VECTOR
        
!
!----COMPUTE THE "LOW LEVEL" WIND MAGNITUDE (M/S)
!
        ULOW(I) = MAX(SQRT(UBAR(I)*UBAR(I) + VBAR(I)*VBAR(I)), 1.0)
        ULOI(I) = 1.0 / ULOW(I)
      ENDDO
!
      DO  K = 1,KMM1
        kp1 = k + 1
        DO  I = 1,npt
          J            = ipt(i)
          VELCO(I,K)   = 0.5 * ((U1(J,K)+U1(J,Kp1))*UBAR(I)
     &                       +  (V1(J,K)+V1(J,Kp1))*VBAR(I))
          VELCO(I,K)   = VELCO(I,K) * ULOI(I)
!         IF ((VELCO(I,K).LT.VELEPS) .AND. (VELCO(I,K).GT.0.)) THEN
!           VELCO(I,K) = VELEPS
!         ENDIF
        ENDDO
      ENDDO
!      
!
!   find the interface level of the projected wind where
!   low levels & upper levels meet above pbl
!
!     do i=1,npt
!       kint(i) = km
!     enddo
!     do k = 1,kmm1
!       do i = 1,npt
!         IF (K .GT. kref(I)) THEN
!           if(velco(i,k) .lt. veleps .and. kint(i) .eq. km) then
!             kint(i) = k+1
!           endif
!         endif
!       enddo
!     enddo
!  WARNING  KINT = KREF !!!!!!!!!
      do i=1,npt
        kint(i) = kref(i)
      enddo
!
!     if(lprnt) print *,' ubar=',ubar
!    &,' vbar=',vbar,' ulow=',ulow,' veleps=',veleps
!
      DO I = 1,npt
        J      = ipt(i)
        BNV    = SQRT( BNV2bar(I) )
        FR     = BNV     * ULOI(I) * min(HPRIME(J),hpmax)
        FR     = MIN(FR, FRMAX)
        XN(I)  = UBAR(I) * ULOI(I)
        YN(I)  = VBAR(I) * ULOI(I)
!
!     Compute the base level stress and store it in TAUB
!     CALCULATE ENHANCEMENT FACTOR, NUMBER OF MOUNTAINS & ASPECT
!     RATIO CONST. USE SIMPLIFIED RELATIONSHIP BETWEEN STANDARD
!     DEVIATION & CRITICAL HGT
!
!> - Calculate enhancement factor (E),number of mountans (m') and
!! aspect ratio constant.
!!\n As in eq.(4.9),(4.10),(4.11) in Kim and Arakawa (1995)
!! \cite kim_and_arakawa_1995, we define m' and E in such a way that they
!! depend on the geometry and location of the subgrid-scale orography
!! through OA and the nonlinearity of flow above the orography through
!! Fr. OC, which is the orographic convexity, and statistically
!! determine how protruded (sharp) the subgrid-scale orography is, is
!! included in the saturation flux G' in such a way that G' is
!! proportional to OC. The forms of E,m' and G' are:
!!\f[
!!  E(OA,F_{r_{0}})=(OA+2)^{\delta}
!!\f]
!!\f[
!!  \delta=C_{E}F_{r_{0}}/F_{r_{c}}
!!\f]
!!\f[
!!  m'(OA,CLX)=C_{m}\triangle x(1+CLX)^{OA+1}
!!\f]
!!\f[
!!  G'(OC,F_{r_{0}})=\frac{F_{r_{0}}^2}{F_{r_{0}}^2+a^{2}}
!!\f]
!!\f[
!! a^{2}=C_{G}OC^{-1}
!!\f]
!! where \f$F_{r_{c}}(=1)\f$ is the critical Froude number,
!! \f$F_{r_{0}}\f$ is the Froude number. \f$C_{E}\f$,\f$C_{m}\f$,
!! \f$C_{G}\f$ are constants.

!> - Calculate the reference-level drag \f$\tau_{0}\f$ (eq.(4.8) in
!! Kim and Arakawa (1995) \cite kim_and_arakawa_1995):
!!\f[
!! \tau_0=E\frac{m'}{\triangle x}\frac{\rho_{0}U_0^3}{N_{0}}G'
!!\f]
!! where \f$E\f$,\f$m'\f$, and \f$G'\f$ are the enhancement factor,
!! "the number of mountains", and the flux function defined above,
!! respectively.

        EFACT    = (OA(I) + 2.) ** (CEOFRC*FR)
        EFACT    = MIN( MAX(EFACT,EFMIN), EFMAX )
!
        COEFM    = (1. + CLX(I)) ** (OA(I)+1.)
!
        XLINV(I) = COEFM * CLEFF
!
        TEM      = FR    * FR * OC(J)
        GFOBNV   = GMAX  * TEM / ((TEM + CG)*BNV)  ! G/N0
!
        TAUB(I)  = XLINV(I) * ROLL(I) * ULOW(I) * ULOW(I)
     &           * ULOW(I)  * GFOBNV  * EFACT         ! BASE FLUX Tau0
!
!         tem      = min(HPRIME(I),hpmax)
!         TAUB(I)  = XLINV(I) * ROLL(I) * ULOW(I) * BNV * tem * tem
!
        K        = MAX(1, kref(I)-1)
        TEM      = MAX(VELCO(I,K)*VELCO(I,K), 0.1)
        SCOR(I)  = BNV2(I,K) / TEM  ! Scorer parameter below ref level
      ENDDO
!     if(lprnt) print *,' taub=',taub
!                                                                       
!----SET UP BOTTOM VALUES OF STRESS
!
      DO K = 1, KBPS
        DO I = 1,npt
          IF (K <= kref(I)) TAUP(I,K) = TAUB(I)
        ENDDO
      ENDDO
!
!   Now compute vertical structure of the stress.
!
      DO K = KMPS, KMM1                   ! Vertical Level K Loop!
        KP1 = K + 1
        DO I = 1, npt
!
!-----UNSTABLE LAYER IF RI < RIC
!-----UNSTABLE LAYER IF UPPER AIR VEL COMP ALONG SURF VEL <=0 (CRIT LAY)
!---- AT (U-C)=0. CRIT LAYER EXISTS AND BIT VECTOR SHOULD BE SET (.LE.)
!
          IF (K >= kref(I)) THEN
            ICRILV(I) = ICRILV(I) .OR. ( ri_n(I,K) <  RIC)
     &                            .OR. (VELCO(I,K) <= 0.0)
          ENDIF
        ENDDO
!
!> - Compute the drag above the reference level (\f$k\geq kref\f$):
!!  - Calculate the ratio of the Scorer parameter (\f$R_{scor}\f$).
!! \n From a series of experiments, Kim and Arakawa (1995)
!! \cite kim_and_arakawa_1995 found that the magnitude of drag divergence
!! tends to be underestimated by the revised scheme in low-level
!! downstream regions with wave breaking. Therefore, at low levels when
!! OA > 0 (i.e., in the "downstream" region) the saturation hypothesis
!! is replaced by the following formula based on the ratio of the
!! the Scorer parameter:
!!\f[
!! R_{scor}=\min \left[\frac{\tau_i}{\tau_{i+1}},1\right]
!!\f]
        DO I = 1,npt
          IF (K >= kref(I))   THEN
            IF (.NOT.ICRILV(I) .AND. TAUP(I,K) > 0.0 ) THEN
              TEMV = 1.0 / max(VELCO(I,K), 0.01)
!             IF (OA(I) .GT. 0. .AND.  PRSI(ipt(i),KP1).GT.RLOLEV) THEN
              IF (OA(I) > 0. .AND. kp1 < kint(i)) THEN
                SCORK   = BNV2(I,K) * TEMV * TEMV
                RSCOR   = MIN(1.0, SCORK / SCOR(I))
                SCOR(I) = SCORK
              ELSE
                RSCOR   = 1.
              ENDIF
!
!>  - The drag above the reference level is expressed as:
!!\f[
!! \tau=\frac{m'}{\triangle x}\rho NUh_d^2
!!\f]
!! where \f$h_{d}\f$ is the displacement wave amplitude. In the absence
!! of wave breaking, the displacement amplitude for the \f$i^{th}\f$
!! layer can be expressed using the drag for the layer immediately
!! below. Thus, assuming \f$\tau_i=\tau_{i+1}\f$, we can get:
!!\f[
!! h_{d_i}^2=\frac{\triangle x}{m'}\frac{\tau_{i+1}}{\rho_{i}N_{i}U_{i}}
!!\f]

              BRVF = SQRT(BNV2(I,K))        ! Brunt-Vaisala Frequency
!             TEM1 = XLINV(I)*(RO(I,KP1)+RO(I,K))*BRVF*VELCO(I,K)*0.5
              TEM1 = XLINV(I)*(RO(I,KP1)+RO(I,K))*BRVF*0.5
     &                       * max(VELCO(I,K),0.01)
              HD   = SQRT(TAUP(I,K) / TEM1)
              FRO  = BRVF * HD * TEMV
!
!    RIM is the  MINIMUM-RICHARDSON NUMBER BY SHUTTS (1985)
!
!> - The minimum Richardson number (\f$Ri_{m}\f$) or local
!! wave-modified Richardson number, which determines the onset of wave
!! breaking, is expressed in terms of \f$R_{i}\f$ and
!! \f$F_{r_{d}}=Nh_{d}/U\f$:
!!\f[
!! Ri_{m}=\frac{Ri(1-Fr_{d})}{(1+\sqrt{Ri}\cdot Fr_{d})^{2}}
!!\f]
!! see eq.(4.6) in Kim and Arakawa (1995) \cite kim_and_arakawa_1995.

              TEM2   = SQRT(ri_n(I,K))
              TEM    = 1. + TEM2 * FRO
              RIM    = ri_n(I,K) * (1.-FRO) / (TEM * TEM)
!
!    CHECK STABILITY TO EMPLOY THE 'SATURATION HYPOTHESIS'
!    OF LINDZEN (1981) EXCEPT AT TROPOSPHERIC DOWNSTREAM REGIONS
!
!>  - Check stability to employ the 'saturation hypothesis' of Lindzen
!! (1981) \cite lindzen_1981 except at tropospheric downstream regions.
!! \n Wave breaking occurs when \f$Ri_{m}<Ri_{c}=0.25\f$. Then
!! Lindzen's wave saturation hypothesis resets the displacement
!! amplitude \f$h_{d}\f$ to that corresponding to \f$Ri_{m}=0.25\f$,
!! we obtain the critical \f$h_{d}\f$(or \f$h_{c}\f$) expressed in
!! terms of the mean values of \f$U\f$, \f$N\f$, and \f$Ri\f$ (
!! eq.(4.7) in Kim and Arakawa (1995) \cite kim_and_arakawa_1995):
!!\f[
!! h_{c}=\frac{U}{N}\left\{2(2+\frac{1}{\sqrt{Ri}})^{1/2}-(2+\frac{1}{\sqrt{Ri}})\right\}
!!\f]
!! if \f$Ri_{m}\leq Ri_{c}\f$, obtain \f$\tau\f$ from the drag above
!! the reference level by using \f$h_{c}\f$ computed above; otherwise
!! \f$\tau\f$ is unchanged (note: scaled by the ratio of the Scorer
!! paramter).
!                                       ----------------------
              IF (RIM <= RIC .AND.
!    &           (OA(I) .LE. 0. .OR.  PRSI(ipt(I),KP1).LE.RLOLEV )) THEN
     &           (OA(I) <= 0. .OR.  kp1 >= kint(i) )) THEN
                 TEMC = 2.0 + 1.0 / TEM2
                 HD   = VELCO(I,K) * (2.*SQRT(TEMC)-TEMC) / BRVF
                 TAUP(I,KP1) = TEM1 * HD * HD
              ELSE
                 TAUP(I,KP1) = TAUP(I,K) * RSCOR
              ENDIF
              taup(i,kp1) = min(taup(i,kp1), taup(i,k))
            ENDIF
          ENDIF
        ENDDO
      ENDDO
!
!     DO I=1,IM
!       taup(i,km+1) = taup(i,km)
!     ENDDO
!
      IF(LCAP .LE. KM) THEN
         DO KLCAP = LCAPP1, KM+1
            DO I = 1,npt
              SIRA          = PRSI(ipt(I),KLCAP) / PRSI(ipt(I),LCAP)
              TAUP(I,KLCAP) = SIRA * TAUP(I,LCAP)
            ENDDO
         ENDDO
      ENDIF
!
!     Calculate - (g/p*)*d(tau)/d(sigma) and Decel terms DTAUX, DTAUY
!
      DO K = 1,KM
        DO I = 1,npt
          TAUD(I,K) = G * (TAUP(I,K+1) - TAUP(I,K)) / DEL(ipt(I),K)
        ENDDO
      ENDDO
!
!------LIMIT DE-ACCELERATION (MOMENTUM DEPOSITION ) AT TOP TO 1/2 VALUE
!------THE IDEA IS SOME STUFF MUST GO OUT THE 'TOP'
!
      DO KLCAP = LCAP, KM
         DO I = 1,npt
            TAUD(I,KLCAP) = TAUD(I,KLCAP) * FACTOP
         ENDDO
      ENDDO
!
!------IF THE GRAVITY WAVE DRAG WOULD FORCE A CRITICAL LINE IN THE
!------LAYERS BELOW SIGMA=RLOLEV DURING THE NEXT DELTIM TIMESTEP,
!------THEN ONLY APPLY DRAG UNTIL THAT CRITICAL LINE IS REACHED.
!
      DO K = 1,KMM1
        DO I = 1,npt
           IF (K > kref(I) .and. PRSI(ipt(i),K) >= RLOLEV) THEN
             IF(TAUD(I,K) /= 0.) THEN
               TEM = DELTIM * TAUD(I,K)
               DTFAC(I) = MIN(DTFAC(I),ABS(VELCO(I,K)/TEM))
             ENDIF
           ENDIF
        ENDDO
      ENDDO
!
!     if(lprnt .and. npr .gt. 0) then
!       print *,' before  A=',A(npr,:)
!       print *,' before  B=',B(npr,:)
!     endif

!> - Calculate outputs: A, B, DUSFC, DVSFC (see parameter description).
!!  - Below the dividing streamline height (k < idxzb), mountain
!!    blocking(\f$D_{b}\f$) is applied.
!!  - Otherwise (k>= idxzb), orographic GWD (\f$\tau\f$) is applied.
      DO K = 1,KM
        DO I = 1,npt
          J         = ipt(i)
          TAUD(I,K) = TAUD(I,K) * DTFAC(I)
          DTAUX     = TAUD(I,K) * XN(I)
          DTAUY     = TAUD(I,K) * YN(I)
          ENG0      = 0.5*(U1(j,K)*U1(j,K)+V1(J,K)*V1(J,K))
! ---  lm mb (*j*)  changes overwrite GWD
          if ( K < IDXZB(I) .AND. IDXZB(I) /= 0 ) then
            DBIM    =   DB(I,K) / (1.+DB(I,K)*DELTIM)
            A(J,K)  = - DBIM * V1(J,K) + A(J,K)
            B(J,K)  = - DBIM * U1(J,K) + B(J,K)
            ENG1    = ENG0*(1.0-DBIM*DELTIM)*(1.0-DBIM*DELTIM)
!          if ( ABS(DBIM * U1(J,K)) .gt. .01 )
!    & print *,' in gwdps_lmi.f KDT=',KDT,I,K,DB(I,K),
!    &                      dbim,idxzb(I),U1(J,K),V1(J,K),me
            DUSFC(J) = DUSFC(J) - DBIM * U1(J,K) * DEL(J,K)
            DVSFC(J) = DVSFC(J) - DBIM * V1(J,K) * DEL(J,K)
          else
!
            A(J,K)   = DTAUY     + A(J,K)
            B(J,K)   = DTAUX     + B(J,K)
            ENG1     = 0.5*(
     &                 (U1(J,K)+DTAUX*DELTIM)*(U1(J,K)+DTAUX*DELTIM)
     &               + (V1(J,K)+DTAUY*DELTIM)*(V1(J,K)+DTAUY*DELTIM))
            DUSFC(J) = DUSFC(J)  + DTAUX * DEL(J,K)
            DVSFC(J) = DVSFC(J)  + DTAUY * DEL(J,K)
          endif
          C(J,K) = C(J,K) + max(ENG0-ENG1,0.)/CP/DELTIM
        ENDDO
      ENDDO
!     if (lprnt) then
!       print *,' in gwdps_lm.f after  A=',A(ipr,:)
!       print *,' in gwdps_lm.f after  B=',B(ipr,:)
!       print *,' DB=',DB(ipr,:)
!     endif

      DO I = 1,npt
        J          = ipt(i)
!       TEM    = (-1.E3/G)
        DUSFC(J) = - onebg * DUSFC(J)
        DVSFC(J) = - onebg * DVSFC(J)
      ENDDO
!
!    MONITOR FOR EXCESSIVE GRAVITY WAVE DRAG TENDENCIES IF NCNT>0
!
!     IF(NCNT.GT.0) THEN
!        IF(LAT.GE.38.AND.LAT.LE.42) THEN
!CMIC$ GUARD 37
!           DO 92 I = 1,IM
!              IF(IKOUNT.GT.NCNT) GO TO 92
!              IF(I.LT.319.OR.I.GT.320) GO TO 92
!              DO 91 K = 1,KM
!                 IF(ABS(TAUD(I,K)) .GT. CRITAC) THEN
!                    IF(I.LE.IM) THEN
!                       IKOUNT = IKOUNT+1
!                       PRINT 123,I,LAT,KDT
!                       PRINT 124,TAUB(I),BNV(I),ULOW(I),
!    1                  GF(I),FR(I),ROLL(I),HPRIME(I),XN(I),YN(I)
!                       PRINT 124,(TAUD(I,KK),KK = 1,KM)
!                       PRINT 124,(TAUP(I,KK),KK = 1,KM+1)
!                       PRINT 124,(ri_n(I,KK),KK = 1,KM)
!                       DO 93 KK = 1,KMM1
!                          VELKO(KK) =
!    1                  0.5*((U1(I,KK)+U1(I,KK+1))*UBAR(I)+
!    2                  (V1(I,KK)+V1(I,KK+1))*VBAR(I))*ULOI(I)
!93                     CONTINUE
!                       PRINT 124,(VELKO(KK),KK = 1,KMM1)
!                       PRINT 124,(A    (I,KK),KK = 1,KM)
!                       PRINT 124,(DTAUY(I,KK),KK = 1,KM)
!                       PRINT 124,(B    (I,KK),KK = 1,KM)
!                       PRINT 124,(DTAUX(I,KK),KK = 1,KM)
!                       GO TO 92
!                    ENDIF
!                 ENDIF
!91            CONTINUE
!92         CONTINUE
!CMIC$ END GUARD 37
!123        FORMAT('  *** MIGWD PRINT *** I=',I3,' LAT=',I3,' KDT=',I3)
!124        FORMAT(2X,  10E13.6)
!        ENDIF
!     ENDIF
!
!      print *,' in gwdps_lm.f 18  =',A(ipt(1),idxzb(1))
!    &,                          B(ipt(1),idxzb(1)),me
      RETURN
      end subroutine gwdps_run
!> @}

!
!> \section arg_table_gwdps_finalize Argument Table
!!
      subroutine gwdps_finalize()
      end subroutine gwdps_finalize

      end module gwdps

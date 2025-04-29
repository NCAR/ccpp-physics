!>\file ecmwf_ngw.F90
!!

module ecmwf_ngw


contains
!------------------------------------------------------------------------------------------
!   April 2025 Adding ECMWF Non-stationary gravity wave scheme option by  Bo Yang
!------------------------------------------------------------------------------------------
!  different orientation at vertical
!  1 is the highest level for ECMWF, 1 is the lowest level for GFS
!  Vertical levels are reversed after entering the subroutine ecmwf_ngw_emc and reversed back before exiting
!  This non-stationary GWD module was obtained by Fanglin Yang from ECMWF, with permission for operational use at NCEP. We would
!  like to thank Andy Brown, Michael Sleigh, Peter Bechtold, and Nils Wedi at ECMWF for their support in porting this code to the
!  UFS.
!! Original Fortran Code by J. SCINOCCIA
! Rewritten in IFS format by A. ORR E.C.M.W.F. August 2008
! PURPOSE
! -------

! THIS ROUTINE COMPUTES NON-OROGRAPHIC GRAVITY WAVE DRAG
! AFTER SCINOCCA (2003) AND Mc LANDRESS AND SCINOCCIA (JAS 2005)
! HYDROSTATIC NON-ROTATIONAL SIMPLIFIED VERSION OF THE
! WARNER AND MCINTYRE (1996) NON-OROGRAPHIC GRAVITY WAVE PARAMETERIZATION
! CONSTANTS HAVE BEEN OPTIMIZED FOLLOWING M. ERN ET AL. (ATMOS. CHEM. PHYS. 2006)

! REFERENCE: Orr, A., P. Bechtold, J. Scinoccia, M. Ern, M. Janiskova, 2010:
! Improved middle atmosphere climate and analysis in the ECMWF forecasting system
! through a non-orographic gravity wave parametrization. J. Climate., 23, 5905-5926.

! LAUNCH SPECTRUM - GENERALIZED DESAUBIES
! INCLUDES A CRITICAL-LEVEL CORRECTION THAT PREVENTS THE
! MOMEMTUM DEPOSITION IN EACH AZIMUTH FROM DRIVING THE FLOW TO SPEEDS FASTER
! THAN THE PHASE SPEED OF THE WAVES, I.E. WHEN WAVES BREAK THEY DRAG THE MEAN
! FLOW TOWARDS THEIR PHASE SPEED - NOT PAST IT.
!!! different orientation for vertical
!!! 1 is the highest level for ECMWF, 1 is the lowest level for GFS
!---------------------------------------------------

!      subroutine cires_ugwpv1_ngw_solv2(mpi_id, master, im, levs, kdt, dtp, &
!                 tau_ngw, tm , um, vm, qm, prsl, prsi, zmet,  zmeti, prslk, &
!                 xlatd, sinlat, coslat,                                     &
!             pdudt, pdvdt, pdtdt, dked, zngw)


       subroutine ecmwf_ngw_emc(mpi_id, master, KLON, KLEV, kdt, PTSTEP, DX,        &
        tau_ngw, PTM11, PUM11, PVM11, qm1, PAPM11, PAPHM11, PGEO11,  zmeti1, prslk1, &
                 xlatd, sinlat, coslat,                                              &
             PTENU, PTENV, pdtdt, dked, zngw)


!
      use machine,          only : kind_phys


      use cires_ugwpv1_module,only :  psrc => knob_ugwp_palaunch

      use cires_ugwpv1_module,only : maxdudt, maxdtdt, max_eps, dked_min, dked_max
!!  maxdudt=250.e-5; maxdtdt=15.e-2; dked_min=0.01; dked_max=250.0
!!  max_eps=max_kdis*4.e-4=450*4.e-4

      use ugwp_common ,     only : rgrav,  grav,  cpd,    rd,  rv, rcpdl, grav2cpd,    &
                                   omega2,  rcpd,   rcpd2,  pi,    pi2, fv,            &
                                   rad_to_deg, deg_to_rad,                             &
                                   rdi,        gor,    grcp,   gocp,                   &
                                   bnv2min,  bnv2max,  dw2min, velmin, gr2,            &
                                   hpscale, rhp, rh4, grav2, rgrav2, mkzmin, mkz2min

!!! grav=con_g; rgrav=1/grav; cpd=con_cp; rd=con_rd;  rv=con_rv; rcpdl=cpd*rgrav=cpd/g
!!!grav2cpd= grav*grcp=grav*grav*rcpd= g**2/cpd omega=con_omega; omega2=2*omega1
!!  rcpd2=0.5*rcpd=1/(2*cpd)
!!   pi=con_pi; pi2=2*pi
!!   fv=con_fvirt; rad_to_deg=180.0/pi; deg_to_rad=pi2/180.0
!!   rdi= 1.0/rd; gor=grav/rd;  grcp= grav*rcpd=g/Cpd; gocp=grcp=grav*rcpd= g/cpd
!!   bnv2min=(pi2/1800.)*(pi2/1800.); bnv2max=(pi2/30.)*(pi2/30.)
!!   dw2min=1.0, velmin=sqrt(dw2min)
!!   gr2=grav*gor=grav*grav/rd=g**2/rd
!!   hpscale=7000; rhp = 1./hpscale = 1/7000; rh4=rhp2*rhp2=(0.5*rhp)**2=1/4*1/7000=1/28000
!!   rgrav2=rgrav*rgrav=1/(g**2); grav2=grav+grav=2*g; mkzmin=pi2/80.0e-3=2*pi/80.0e3
!!    mkz2min= mkzmin*mkzmin=(2*pi/80.0E3)**2


      use ugwp_wmsdis_init, only : v_kxw,  rv_kxw,   v_kxw2, tamp_mpa, tau_min, ucrit, &
                                   gw_eff,                                             &
!                                   zms,                                                &
                                   zci4, zci3, zci2,                                   &
                                   rimin, sc2, sc2u, ric
!   v_kxw=kxw=pi2/lhmet=pi2/200e3; rv_kxw=200e3/pi2
!   v_kxw2=v_kxw*v_kxw=(pi2/200e3)**2; tamp_mpa=knob_ugwp_tauamp  amplitude for GEOS-5/MERRA-2
!   tau_min=min of GW MF 0.25mPa
!   ucrit=cdmin=2e-2/mkzmax=2e-2/((2pie)/500) = 10/(2pie)
!   gw_eff=effac=1.0
!   zci4=(zms*zci(inc))**4; zci2=(zms*zci(inc))**2; zci3(inc)=(zms*zci(inc))**3
!   rimin=-10.0; sc2=lturb*lturb=30m*30m; sc2u=ulturb*ulturb=150*150


      implicit none

!in
!work


!      integer, intent(in)  :: KLAUNCH                           ! index for launch level
      integer, intent(in)  :: KLEV                            ! vertical level
      integer, intent(in)  :: KLON                              ! horiz tiles
      integer, intent(in)  :: mpi_id, master, kdt

      real(kind=kind_phys) ,intent(in)   :: PTSTEP               ! model time step
      real(kind=kind_phys) ,intent(in)   :: DX(KLON)             ! model grid size

      real(kind=kind_phys) ,intent(in)   :: tau_ngw(KLON)

      real(kind=kind_phys) ,intent(in)   :: PVM11(KLON,KLEV)       ! meridional wind
      real(kind=kind_phys) ,intent(in)   :: PUM11(KLON,KLEV)       ! zonal wind
      real(kind=kind_phys) ,intent(in)   :: qm1(KLON,KLEV)       ! spec. humidity
      real(kind=kind_phys) ,intent(in)   :: PTM11(KLON,KLEV)       ! kinetic temperature
      real(kind=kind_phys) ,intent(in)   :: PAPM11(KLON,KLEV)     ! mid-layer pressure
      real(kind=kind_phys) ,intent(in)   :: PAPHM11(KLON,KLEV+1)   !  interface pressure
      real(kind=kind_phys) ,intent(in)   :: PGEO11(KLON,KLEV)   !  full model level geopotential in meters

      real(kind=kind_phys) :: PVM1(KLON,KLEV)       ! meridional wind
      real(kind=kind_phys) :: PUM1(KLON,KLEV)       ! zonal wind
      real(kind=kind_phys) :: qm(KLON,KLEV)       ! spec. humidity
      real(kind=kind_phys) :: PTM1(KLON,KLEV)       ! kinetic temperature
      real(kind=kind_phys) :: PAPM1(KLON,KLEV)     ! mid-layer pressure
      real(kind=kind_phys) :: PAPHM1(KLON,KLEV+1)   !  interface pressure
      real(kind=kind_phys) :: PGEO1(KLON,KLEV)   !  full model level geopotential in meters

!      real(kind=kind_phys)  :: PGAW(KLON) !normalised gaussian quadrature weight/nb of longitude pts 
                                                       ! local sub-area == 4*RPI*RA**2 * PGAW
!      real(kind=kind_phys) ,intent(in)   :: PPRECIP(KLON) ! total surface precipitation




      real(kind=kind_phys) ,intent(in)   :: prslk1(KLON,KLEV)    ! mid-layer exner function
       real(kind=kind_phys) :: prslk(KLON,KLEV)    ! mid-layer exner function
!      real(kind=kind_phys) ,intent(in)   :: zmet(KLON,KLEV)     ! meters phil =philg/grav ! use PGEO1 instead

      real(kind=kind_phys) ,intent(in)   :: zmeti1(KLON,KLEV+1)  !  interface geopi/meters
      real(kind=kind_phys)  :: zmeti(KLON,KLEV+1)  !  interface geopi/meters
      real(kind=kind_phys) ,intent(in)   :: xlatd(KLON)         ! xlat_d in degrees
      real(kind=kind_phys) ,intent(in)   :: sinlat(KLON)
      real(kind=kind_phys) ,intent(in)   :: coslat(KLON)
!
! out-gw effects
!
      real(kind=kind_phys) ,intent(out) :: PTENU(KLON,KLEV)     ! zonal momentum tendency
      real(kind=kind_phys) ,intent(out) :: PTENV(KLON,KLEV)     ! meridional momentum tendency
      real(kind=kind_phys) ,intent(out) :: pdtdt(KLON,KLEV)     ! gw-heating (u*ax+v*ay)/cp and cooling


      real(kind=kind_phys) :: PFLUXU(KLON,KLEV+1)     ! zonal momentum tendency
      real(kind=kind_phys) :: PFLUXV(KLON,KLEV+1)     ! meridional momentum tendency

      real(kind=kind_phys) ,intent(out) :: dked(KLON,KLEV)      ! gw-eddy diffusion

      real(kind=kind_phys) ,intent(out) :: zngw(KLON)           ! launch height
!
!work
      INTEGER,      PARAMETER   :: IAZIDIM=4     !number of azimuths
      INTEGER,      PARAMETER   :: INCDIM=20     !number of discretized c spectral elements in launch spectrum

      REAL(kind=kind_phys),  PARAMETER    :: RA=6370.e3             !half-model level zonal velocity
      REAL(kind=kind_phys),  PARAMETER    :: GPTWO=2.0              ! 2p in equation 

      REAL(kind=kind_phys) :: ZUHM1(KLON,KLEV)             !half-model level zonal velocity
      REAL(kind=kind_phys) :: ZVHM1(KLON,KLEV)             !half-model level meridional velocity

      REAL(kind=kind_phys) :: ZBVFHM1(KLON,KLEV)           !half-model level Brunt-Vaisalla frequency
      REAL(kind=kind_phys) :: ZRHOHM1(KLON,KLEV)           !half-model level density
      REAL(kind=kind_phys) :: ZX(INCDIM)                   !coordinate transformation
      REAL(kind=kind_phys) :: ZCI(INCDIM)                  !phase speed element
      REAL(kind=kind_phys) :: ZDCI(INCDIM)
      REAL(kind=kind_phys) :: ZUI(KLON,KLEV,IAZIDIM)       !intrinsic velocity
      REAL(kind=kind_phys) :: ZUL(KLON,IAZIDIM)            !velocity in azimuthal direction at launch level
      REAL(kind=kind_phys) :: ZBVFL(KLON)                  !buoyancy at launch level
      REAL(kind=kind_phys) :: ZCOSANG(IAZIDIM)             !cos of azimuth angle
      REAL(kind=kind_phys) :: ZSINANG(IAZIDIM)             !sin of azimuth angle
      REAL(kind=kind_phys) :: ZFCT(KLON,KLEV)
      REAL(kind=kind_phys) :: ZFNORM(KLON)                 !normalisation factor (A)
      REAL(kind=kind_phys) :: ZCI_MIN(KLON,IAZIDIM)
      REAL(kind=kind_phys) :: ZTHM1(KLON,KLEV)             !temperature on half-model levels
      REAL(kind=kind_phys) :: ZFLUX(KLON,INCDIM,IAZIDIM)   !momentum flux at each vertical level and azimuth
      REAL(kind=kind_phys) :: ZPU(KLON,KLEV,IAZIDIM)       !momentum flux
      REAL(kind=kind_phys) :: ZDFL(KLON,KLEV,IAZIDIM)
      REAL(kind=kind_phys) :: ZACT(KLON,INCDIM,IAZIDIM)    !if =1 then critical level encountered
      REAL(kind=kind_phys) :: ZACC(KLON,INCDIM,IAZIDIM)
      REAL(kind=kind_phys) :: ZCRT(KLON,KLEV,IAZIDIM)

      INTEGER :: ILAUNCH                   !model level from which GW spectrum is launched
      INTEGER :: INC, JK, JL, IAZI


      REAL(KIND=kind_phys) :: ZRADTODEG, ZGELATDEG
      REAL(KIND=kind_phys) :: ZCIMIN, ZCIMAX
      REAL(KIND=kind_phys) :: ZGAM, ZPEXP, ZXMAX, ZXMIN, ZXRAN
      REAL(KIND=kind_phys) :: ZDX

      REAL(KIND=kind_phys) :: ZX1, ZX2, ZDXA, ZDXB, ZDXS
      REAL(KIND=kind_phys) :: ZANG, ZAZ_FCT, ZNORM, ZANG1, ZTX
      REAL(KIND=kind_phys) :: ZU, ZCIN, ZCPEAK
      REAL(KIND=kind_phys) :: ZCIN4, ZBVFL4, ZCIN2, ZBVFL2, ZCIN3, ZBVFL3, ZCINC
      REAL(KIND=kind_phys) :: ZATMP, ZFLUXS, ZDEP, ZFLUXSQ, ZULM, ZDFT, ZE1, ZE2
      REAL(KIND=kind_phys) :: ZMS_L,ZMS, Z0P5, Z0P0, Z50S
      REAL(KIND=kind_phys) :: ZGAUSS(KLON), ZFLUXLAUN(KLON), ZCNGL(KLON)
      REAL(KIND=kind_phys) :: ZCONS1,ZCONS2,ZDELP,ZRGPTS

!!!  try to assign values for the following 

      REAL(KIND=kind_phys) :: GSSEC
      REAL(KIND=kind_phys) :: ZGAUSSB,ZFLUXGLOB

!      REAL(KIND=kind_phys) :: PGAW(KLON)  don't need this value, set ZDX directly


!      REAL(KIND=kind_phys) ::  PGELAT(KLON) !!! use xlatd from gfs instead 
      REAL(KIND=kind_phys) ::  GCOEFF, GGAUSSA           !!! GCOEFF link to precip
!      REAL(KIND=kind_phys) ::   GGAUSSB                  !!! GGAUSSB->ZGAUSS
                                                         !!!! GGAUSSA
      REAL(KIND=kind_phys) ::  GCSTAR


      LOGICAL :: LGACALC, LGSATL, LOZPR

      integer :: NGAUSS, NSLOPE
     

      NSLOPE=1
      LGACALC=.false.
      LGSATL=.false.
      LOZPR=.true.
!      NGAUSS=4
       NGAUSS=1
!      GGAUSSA=20._kind_phys
!       GGAUSSA=10._kind_phys
        GGAUSSA=5._kind_phys
!      GGAUSSB=1.0_kind_phys
!        ZGAUSSB=0.25_kind_phys
!         ZGAUSSB=0.3_kind_phys
!          ZGAUSSB=0.35_kind_phys
          ZGAUSSB=0.38_kind_phys
!        ZGAUSSB=-0.25_kind_phys
!        ZGAUSSB=0.5_kind_phys
!        ZGAUSSB=0.3_kind_phys
      GCSTAR=1.0_kind_phys

       

!      GSSEC=(pi2/1800.)*(pi2/1800.)
       GSSEC=1.e-24
       GCOEFF=1.0_kind_phys   !!do not know the value, but never use it when NGAUSS is not equal 1
      


!!      LOGIC :: LGINDL !!LGINDL=.true. using standard atm values to calculate!! comment out
!!      REAL(KIND=kind_phys) ::  STPHI(KM),STTEM(KM), STPREH(KM)
!!        REAL(KIND=kind_phys) :: ZTHSTD,ZRHOSTD,ZBVFSTD

! Set parameters which are a function of launch height
!!! need to set the parameter ILAUCH, ZFLUXGLOB,ZGAUSSB,ZMS_L directly
!ILAUNCH=NLAUNCHL(KLAUNCH)
!ZFLUXGLOB=GFLUXLAUNL(KLAUNCH)
!ZGAUSSB=GGAUSSB(KLAUNCH)
!ZMS_L=GMSTAR_L(KLAUNCH)

!*INPUT PARAMETERS
!*       ----------------
      ZRADTODEG=57.29577951_kind_phys
      ZMS_L=2.e3_kind_phys
!       ZFLUXGLOB=3.75e-3_kind_phys
!        ZFLUXGLOB=3.7e-3_kind_phys
!         ZFLUXGLOB=3.55e-3_kind_phys
!           ZFLUXGLOB=3.65e-3_kind_phys  !standard value
!           ZFLUXGLOB=3.62e-3_kind_phys  !standard value
           ZFLUXGLOB=3.60e-3_kind_phys 

!       ZFLUXGLOB=3.2e-3_kind_phys
!        ZFLUXGLOB=3.0e-3_kind_phys
!       ZFLUXGLOB=5.0e-3_kind_phys
!        ZFLUXGLOB=4.0e-3_kind_phys
!       ZFLUXGLOB=0.0_kind_phys

      ZMS=2*pi/ZMS_L

!*       INITIALIZE FIELDS TO ZERO
!*       -------------------------

      PTENU(:,:)=0.0_kind_phys
      PTENV(:,:)=0.0_kind_phys
      pdtdt(:,:)=0.0_kind_phys
      dked(:,:)=0.0_kind_phys

       DO JK=1,KLEV+1
       DO JL=1,KLON
        PFLUXU(JL,JK)=0.0_kind_phys
        PFLUXV(JL,JK)=0.0_kind_phys
       ENDDO
       ENDDO


      DO IAZI=1,IAZIDIM
        DO JK=1,KLEV
          DO JL=1,KLON
           ZPU(JL,JK,IAZI)=0.0_kind_phys
           ZCRT(JL,JK,IAZI)=0.0_kind_phys
           ZDFL(JL,JK,IAZI)=0.0_kind_phys
          ENDDO
        ENDDO
      ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!
!       redefine ilaunch
       
      DO JK=1, KLEV
        if (PAPM11(KLON,JK) .LT. psrc) exit
      ENDDO
      ILAUNCH = max(JK-1,3)

      DO JL=1,KLON
      zngw(JL) = PGEO11(JL, ILAUNCH)
      ENDDO

      ILAUNCH = KLEV + 1 - ILAUNCH

!!!!!!!!!!!!!!!!!!!!!!!!!!!
       

!* reverse vertical coordinate to ECMWF
      DO JL=1,KLON
      PTM1(JL,:)=transfer(PTM11(JL,KLEV:1:-1),PTM11(JL,:))
      PUM1(JL,:)=transfer(PUM11(JL,KLEV:1:-1),PUM11(JL,:))
      PVM1(JL,:)=transfer(PVM11(JL,KLEV:1:-1),PVM11(JL,:))
      qm(JL,:)=transfer(qm1(JL,KLEV:1:-1),qm1(JL,:))
      PAPM1(JL,:)=transfer(PAPM11(JL,KLEV:1:-1),PAPM11(JL,:))
      PGEO1(JL,:)=transfer(PGEO11(JL,KLEV:1:-1),PGEO11(JL,:))
      prslk(JL,:)=transfer(prslk1(JL,KLEV:1:-1),prslk1(JL,:))
 
      dked(JL,:)=transfer(dked(JL,KLEV:1:-1),dked(JL,:))
      PTENU(JL,:)=transfer(PTENU(JL,KLEV:1:-1),PTENU(JL,:))
      PTENV(JL,:)=transfer(PTENV(JL,KLEV:1:-1),PTENV(JL,:))
      pdtdt(JL,:)=transfer(pdtdt(JL,KLEV:1:-1),pdtdt(JL,:))


      PAPHM1(JL,:)=transfer(PAPHM11(JL,KLEV+1:1:-1),PAPHM11(JL,:))
      zmeti(JL,:)=transfer(zmeti1(JL,KLEV+1:1:-1),zmeti1(JL,:))
      ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      


!*       INITIALIZE PARAMETERS FOR COORDINATE TRANSFORM
!*       ----------------------------------------------

! ZCIMIN,ZCIMAX - min,max intrinsic launch-level phase speed (c-U_o) (m/s)
! ZGAM - half=width of coordinate stretch

      ZCIMIN=0.50_kind_phys
      ZCIMAX=100.0_kind_phys
      ZGAM=0.25_kind_phys

      ZPEXP=GPTWO/2.0_kind_phys

! set initial min ci in each column and azimuth (used for critical levels)

      DO IAZI=1,IAZIDIM
         DO JL=1,KLON
           ZCI_MIN(JL,IAZI)=ZCIMIN
         ENDDO
      ENDDO


!*       DEFINE HALF MODEL LEVEL WINDS AND TEMPERATURE
!*       -----------------------------------

       DO JK=2,KLEV
       DO JL=1,KLON
       ZTHM1(JL,JK) =0.5_kind_phys*(PTM1(JL,JK-1)+PTM1(JL,JK))
       ZUHM1(JL,JK) =0.5_kind_phys*(PUM1(JL,JK-1)+PUM1(JL,JK))
       ZVHM1(JL,JK) =0.5_kind_phys*(PVM1(JL,JK-1)+PVM1(JL,JK))
       ENDDO
       ENDDO

       JK=1
       DO JL=1,KLON
        ZTHM1(JL,JK)=PTM1(JL,JK)
        ZUHM1(JL,JK)=PUM1(JL,JK)
        ZVHM1(JL,JK)=PVM1(JL,JK)
       ENDDO


!*       DEFINE STATIC STABILITY AND AIR DENSITY ON HALF MODEL LEVELS
!*       ------------------------------------------------------------

!       ZCONS1=1.0_kind_phys/RD
!       ZCONS2=RG**2/RCPD

        ZCONS1=1.0_kind_phys/rd
        ZCONS2=grav2cpd


       DO JK=KLEV,2,-1
         DO JL=1,KLON
!    ZDELP=PAPM1(JL,JK)-PAPM1(JL,JK-1)
            ZDELP=(PGEO1(JL,JK)-PGEO1(JL,JK-1))*grav   !! times grav since import from zmet
            ZRHOHM1(JL,JK)=PAPHM1(JL,JK)*ZCONS1/ZTHM1(JL,JK)
            ZBVFHM1(JL,JK)=ZCONS2/ZTHM1(JL,JK)*&
     &      (1.0_kind_phys+cpd*(PTM1(JL,JK)-PTM1(JL,JK-1))/ZDELP)
!     & (1.0_kind_phys-RCPD*ZRHOHM1(JL,JK)*(PTM1(JL,JK)-PTM1(JL,JK-1))/ZDELP)
           ZBVFHM1(JL,JK)=MAX(ZBVFHM1(JL,JK),GSSEC)
           ZBVFHM1(JL,JK)=SQRT(ZBVFHM1(JL,JK))
         ENDDO
       ENDDO

!*       SET UP AZIMUTH DIRECTIONS AND SOME TRIG FACTORS
!*       -----------------------------------------------

         ZANG=2*pi/IAZIDIM
         ZAZ_FCT=1.0_kind_phys


! get normalization factor to ensure that the same amount of momentum
! flux is directed (n,s,e,w) no mater how many azimuths are selected.
! note, however, the code below assumes a symmetric distribution of
! of azimuthal directions (ie 4,8,16,32,...)

      ZNORM=0.0_kind_phys
      DO IAZI=1,IAZIDIM
      ZANG1=(IAZI-1)*ZANG
      ZCOSANG(IAZI)=COS(ZANG1)
      ZSINANG(IAZI)=SIN(ZANG1)
      ZNORM=ZNORM+ABS(ZCOSANG(IAZI))
      ENDDO
      ZAZ_FCT=2._kind_phys*ZAZ_FCT/ZNORM


!*       DEFINE COORDINATE TRANSFORM
!*       -----------------------------------------------

! note that this is expresed in terms of the intrinsic phase speed
! at launch ci=c-u_o so that the transformation is identical at every
! launch site.
! See Eq. 28-30 of Scinocca 2003.

      ZXMAX=1.0_kind_phys/ZCIMIN
      ZXMIN=1.0_kind_phys/ZCIMAX

      ZXRAN=ZXMAX-ZXMIN
      ZDX=ZXRAN/REAL(INCDIM-1)
      IF(LGACALC) ZGAM=(ZXMAX-ZXMIN)/LOG(ZXMAX/ZXMIN)
!!  LGACALC=.false. ZGAM =0.25   ZX1=0.0007
!!  LGACALC=.true. ZGAM =0.37559 ZX1=0.01

      ZX1=ZXRAN/(EXP(ZXRAN/ZGAM)-1.0_kind_phys)
      ZX2=ZXMIN-ZX1




      DO INC=1,INCDIM 
      ZTX=REAL(INC-1)*ZDX+ZXMIN
      ZX(INC)=ZX1*EXP((ZTX-ZXMIN)/ZGAM)+ZX2                       !Eq. 29 of Scinocca 2003
      ZCI(INC)=1.0_kind_phys/ZX(INC)                                   !Eq. 28 of Scinocca 2003
      ZDCI(INC)=ZCI(INC)**2*(ZX1/ZGAM)*EXP((ZTX-ZXMIN)/ZGAM)*ZDX  !Eq. 30 of Scinocca 2003
      ENDDO


!*       DEFINE INTRINSIC VELOCITY (RELATIVE TO LAUNCH LEVEL VELOCITY) U(Z)-U(Zo), AND COEFFICINETS
!*       ------------------------------------------------------------------------------------------

      DO IAZI=1,IAZIDIM
       DO JL=1,KLON
       ZUL(JL,IAZI)=ZCOSANG(IAZI)*ZUHM1(JL,ILAUNCH)&
     &   +ZSINANG(IAZI)*ZVHM1(JL,ILAUNCH)
       ENDDO
      ENDDO
      DO JL=1,KLON
       ZBVFL(JL)=ZBVFHM1(JL,ILAUNCH)
      ENDDO

      DO JK=2,ILAUNCH
       DO IAZI=1,IAZIDIM
        DO JL=1,KLON
         ZU=ZCOSANG(IAZI)*ZUHM1(JL,JK)+ZSINANG(IAZI)*ZVHM1(JL,JK)
         ZUI(JL,JK,IAZI)=ZU-ZUL(JL,IAZI)
        ENDDO
       ENDDO
      ENDDO

!*       DEFINE RHO(Zo)/N(Zo)
!*       -------------------
      DO JK=2,ILAUNCH
       DO JL=1,KLON
       ZFCT(JL,JK)=ZRHOHM1(JL,JK)/ZBVFHM1(JL,JK)
       ENDDO
      ENDDO

! Optionally set ZFCT at launch level using standard atmos values, to ensure saturation is
! independent of location
!      IF (LGINDL) THEN
!        ZCONS1=1.0_kind_phys/rd
!        ZCONS2=grav2cpd
!        ZDELP=STPHI(ILAUNCH)-STPHI(ILAUNCH-1)  !!probably need to time grav depending on input
!        THSTD=0.5_kind_phys*(STTEM(ILAUNCH-1)+STTEM(ILAUNCH))
!        ZRHOSTD=STPREH(ILAUNCH-1)*ZCONS1/ZTHSTD
!        ZBVFSTD=ZCONS2/ZTHSTD*(1.0_kind_phys+rcpd*
!     & (STTEM(ILAUNCH)-STTEM(ILAUNCH-1))/ZDELP)
!
!        ZBVFSTD=MAX(ZBVFSTD,GSSEC)
!        ZBVFSTD=SQRT(ZBVFSTD)
!       DO JL=1,KLON
!        ZFCT(JL,ILAUNCH)=ZRHOSTD/ZBVFSTD
!       ENDDO
!      ENDIF

!*       SET LAUNCH MOMENTUM FLUX SPECTRAL DENSITY
!*       -----------------------------------------

! Eq. (25) of Scinocca 2003 (not including the 'A' component), and with U-Uo=0
! do this for only one azimuth since it is identical to all azimuths, and it will be renormalized
! Initial spectrum fully saturated if LGSATL

      IF(NSLOPE==1) THEN
! s=1 case
      DO INC=1,INCDIM
      ZCIN=ZCI(INC)
      ZCIN4=(ZMS*ZCIN)**4
      DO JL=1,KLON
         ZBVFL4=ZBVFL(JL)**4
         IF(LGSATL) THEN
         ZFLUX(JL,INC,1)=ZFCT(JL,ILAUNCH)*ZBVFL4*MIN(ZCIN/ZCIN4,&
     &    ZCIN/ZBVFL4)
         ELSE
         ZFLUX(JL,INC,1)=ZFCT(JL,ILAUNCH)*ZBVFL4*ZCIN/(ZBVFL4+ZCIN4)
         ENDIF
         ZACT(JL,INC,1)=1.0_kind_phys
      ENDDO
      ENDDO

      ELSEIF(NSLOPE==2) THEN
! s=2 case
      DO INC=1,INCDIM
      ZCIN=ZCI(INC)
      ZCIN4=(ZMS*ZCIN)**4
      DO JL=1,KLON
         ZBVFL4=ZBVFL(JL)**4
         ZCPEAK=ZBVFL(JL)/ZMS
         IF(LGSATL) THEN
         ZFLUX(JL,INC,1)=ZFCT(JL,ILAUNCH)*ZBVFL4*MIN&
     &   (ZCPEAK/ZCIN4,ZCIN/ZBVFL4)
         ELSE
         ZFLUX(JL,INC,1)=ZFCT(JL,ILAUNCH)*ZBVFL4*ZCIN*ZCPEAK/(ZBVFL4&
     &     *ZCPEAK+ZCIN4*ZCIN)
         ENDIF
         ZACT(JL,INC,1)=1.0_kind_phys
      ENDDO
      ENDDO

      ELSEIF(NSLOPE==-1) THEN
! s=-1 case
      DO INC=1,INCDIM
      ZCIN=ZCI(INC)
      ZCIN2=(ZMS*ZCIN)**2
      DO JL=1,KLON
        ZBVFL2=ZBVFL(JL)**2
        ZFLUX(JL,INC,1)=ZFCT(JL,ILAUNCH)*ZBVFL2*ZCIN/(ZBVFL2+ZCIN2)
        ZACT(JL,INC,1)=1.0_kind_phys
      ENDDO
      ENDDO

      ELSEIF(NSLOPE==0) THEN
! s=0 case
      DO INC=1,INCDIM
       ZCIN=ZCI(INC)
       ZCIN3=(ZMS*ZCIN)**3
      DO JL=1,KLON
        ZBVFL3=ZBVFL(JL)**3
        ZFLUX(JL,INC,1)=ZFCT(JL,ILAUNCH)*ZBVFL3*ZCIN/(ZBVFL3+ZCIN3)
        ZACT(JL,INC,1)=1.0_kind_phys
        ZACC(JL,INC,1)=1.0_kind_phys
      ENDDO
      ENDDO

      ENDIF

!*       NORMALIZE LAUNCH MOMENTUM FLUX
!*       ------------------------------

! (rho x F^H = rho_o x F_p^total)

! integrate (ZFLUX x dX)
      DO INC=1,INCDIM
        ZCINC=ZDCI(INC)
        DO JL=1,KLON
        ZPU(JL,ILAUNCH,1)=ZPU(JL,ILAUNCH,1)+ZFLUX(JL,INC,1)*ZCINC
        ENDDO
      ENDDO


!*       NORMALIZE GFLUXLAUN TO INCLUDE SENSITIVITY TO PRECIPITATION
!*       -----------------------------------------------------------

! Also other options to alter tropical values

! A=ZFNORM in Scinocca 2003.  A is independent of height.
      ZDXA=1.0_kind_phys/29.E3_kind_phys
      ZDXB=1.0_kind_phys/3.5E3_kind_phys
      DO JL=1,KLON
!!    ZDX=MAX(1.E2_kind_phys,2*RA*SQRT(pi*PGAW(JL))) !grid resolution (m)
!!      ZDX=50.E3    ! c192 for c192 usage
!!!   ZDX=25.E3    !C384
!!!      ZDX=13.E3    !C768
 

   !Scaling factor for launch flux depending on grid resolution
   ! smooth reduction below 30 km
      ZDXS=1.0_kind_phys-MIN(1.0_kind_phys,ATAN((MAX(1.0_kind_phys&
     &   /DX(JL),ZDXA)-ZDXA)/(ZDXB-ZDXA)))
       ZFLUXLAUN(JL)=ZFLUXGLOB*ZDXS
      ZFNORM(JL)=ZFLUXLAUN(JL)/ZPU(JL,ILAUNCH,1)
      ENDDO

! If LOZPR=TRUE then vary EPLAUNCH over tropics
      IF (LOZPR) THEN
       IF (NGAUSS==1) THEN
       Z50S=-50.0_kind_phys

       DO JL=1,KLON

!       ZFLUXLAUN(JL)=ZFLUXLAUN(JL)*(1.0_kind_phys+MIN&
!     &   (0.5_kind_phys,GCOEFF*PPRECIP(JL)))     !precip
!       ZFNORM(JL)=ZFLUXLAUN(JL)/ZPU(JL,ILAUNCH,1)

       ZGELATDEG=xlatd(JL)-Z50S
       ZGAUSS(JL)=ZGAUSSB*EXP((-ZGELATDEG*ZGELATDEG)&
     &  /(2*GGAUSSA*GGAUSSA))

!       ZGELATDEG=xlatd(JL)
!       ZGAUSS(JL)=-0.1_kind_phys*EXP((-ZGELATDEG*ZGELATDEG)&
!     &  /(2*GGAUSSA*GGAUSSA))+ZGAUSS(JL)


!       ZGELATDEG=xlatd(JL)
!       ZGAUSS(JL)=-0.05_kind_phys*EXP((-ZGELATDEG*ZGELATDEG)&
!     &  /(2*GGAUSSA*GGAUSSA))+ZGAUSS(JL)

       ZGELATDEG=xlatd(JL)
       ZGAUSS(JL)=0.08_kind_phys*EXP((-ZGELATDEG*ZGELATDEG)&
     &  /(2*GGAUSSA*GGAUSSA))+ZGAUSS(JL)

       ZGELATDEG=xlatd(JL)-50.0_kind_phys
       ZGAUSS(JL)= 0.0022_kind_phys*EXP((-ZGELATDEG*ZGELATDEG)&
     &  /(2*10.*10.))+ZGAUSS(JL)


       ZFLUXLAUN(JL)=(1.0_kind_phys+ZGAUSS(JL))*ZFLUXLAUN(JL)



       ZFNORM(JL)=ZFLUXLAUN(JL)/ZPU(JL,ILAUNCH,1)


       ENDDO
      ELSEIF (NGAUSS==2) THEN

        DO JL=1,KLON
!       ZGELATDEG=PGELAT(JL)*ZRADTODEG
       ZGELATDEG=xlatd(JL)
       ZGAUSS(JL)=ZGAUSSB*EXP((-ZGELATDEG*ZGELATDEG)&
     &  /(2*GGAUSSA*GGAUSSA))
       ZFLUXLAUN(JL)=(1.0_kind_phys+ZGAUSS(JL))*ZFLUXLAUN(JL)
       ZFNORM(JL)=ZFLUXLAUN(JL)/ZPU(JL,ILAUNCH,1)
       ENDDO



      ELSEIF (NGAUSS==4) THEN

! Set latitudinal dependence to optimize stratospheric winds for 36r1
      Z50S=-50.0_kind_phys
      DO JL=1,KLON
!      ZGELATDEG=PGELAT(JL)*ZRADTODEG-Z50S
       ZGELATDEG=xlatd(JL)-Z50S
      ZGAUSS(JL)=ZGAUSSB*EXP((-ZGELATDEG*ZGELATDEG)/(2*GGAUSSA*GGAUSSA))
      ZFLUXLAUN(JL)=(1.0_kind_phys+ZGAUSS(JL))*ZFLUXLAUN(JL)
      ZFNORM(JL)=ZFLUXLAUN(JL)/ZPU(JL,ILAUNCH,1)
      ENDDO

      ENDIF
      ENDIF 

      DO IAZI=1,IAZIDIM
         DO JL=1,KLON
         ZPU(JL,ILAUNCH,IAZI)=ZFLUXLAUN(JL)
         ENDDO
      ENDDO

!*       ADJUST CONSTANT ZFCT
!*       --------------------
      DO JK=2,ILAUNCH
        DO JL=1,KLON
         ZFCT(JL,JK)=ZFNORM(JL)*ZFCT(JL,JK)
        ENDDO
      ENDDO

!*       RENORMALIZE EACH SPECTRAL ELEMENT IN FIRST AZIMUTH
!*       --------------------------------------------------
      DO INC=1,INCDIM
       DO JL=1,KLON
       ZFLUX(JL,INC,1)=ZFNORM(JL)*ZFLUX(JL,INC,1)
       ENDDO
      ENDDO



!*       COPY ZFLUX INTO ALL OTHER AZIMUTHS
!*       --------------------------------

! ZACT=1 then no critical level
! ZACT=0 then critical level

      DO IAZI=2,IAZIDIM
        DO INC=1,INCDIM
         DO JL=1,KLON
         ZFLUX(JL,INC,IAZI)=ZFLUX(JL,INC,1)
         ZACT(JL,INC,IAZI)=1.0_kind_phys
         ZACC(JL,INC,IAZI)=1.0_kind_phys
         ENDDO
       ENDDO
      ENDDO

! -----------------------------------------------------------------------------

!*       BEGIN MAIN LOOP OVER LEVELS
!*       ---------------------------

!* begin IAZIDIM do-loop
!* --------------------

      DO IAZI=1,IAZIDIM

!* begin JK do-loop
!* ----------------

      DO JK=ILAUNCH-1,2,-1

!* first do critical levels
!* ------------------------

         DO JL=1,KLON
            ZCI_MIN(JL,IAZI)=MAX(ZCI_MIN(JL,IAZI),ZUI(JL,JK,IAZI))
         ENDDO

!* set ZACT to zero if critical level encountered
!* ----------------------------------------------

         Z0P5=0.5_kind_phys
         DO INC=1,INCDIM
            ZCIN=ZCI(INC)
            DO JL=1,KLON
               ZATMP=Z0P5+SIGN(Z0P5,ZCIN-ZCI_MIN(JL,IAZI))
               ZACC(JL,INC,IAZI)=ZACT(JL,INC,IAZI)-ZATMP
               ZACT(JL,INC,IAZI)=ZATMP
            ENDDO
         ENDDO

!* integrate to get critical-level contribution to mom deposition on this level, i.e. ZACC=1
!* ----------------------------------------------------------------------------------------

         DO INC=1,INCDIM
            ZCINC=ZDCI(INC)
            DO JL=1,KLON
               ZDFL(JL,JK,IAZI)=ZDFL(JL,JK,IAZI)+&
     &                ZACC(JL,INC,IAZI)*ZFLUX(JL,INC,IAZI)*ZCINC
            ENDDO
         ENDDO

!* get weighted average of phase speed in layer
                                                    
        DO JL=1,KLON
            IF(ZDFL(JL,JK,IAZI)>0.0_kind_phys) THEN
               ZATMP=ZCRT(JL,JK,IAZI)
               DO INC=1,INCDIM
                  ZATMP=ZATMP+ZCI(INC)*&
     &                   ZACC(JL,INC,IAZI)*ZFLUX(JL,INC,IAZI)*ZDCI(INC)
               ENDDO
               ZCRT(JL,JK,IAZI)=ZATMP/ZDFL(JL,JK,IAZI)
            ELSE
               ZCRT(JL,JK,IAZI)=ZCRT(JL,JK+1,IAZI)
            ENDIF
         ENDDO

!* do saturation (Eq. (26) and (27) of Scinocca 2003)
!* -------------------------------------------------

         IF(GPTWO==3.0_kind_phys) THEN
             DO INC=1,INCDIM
                ZCIN=ZCI(INC)
                ZCINC=1.0_kind_phys/ZCIN
                DO JL=1,KLON
                   ZE1=ZCIN-ZUI(JL,JK,IAZI)
                   ZE2=GCSTAR*ZFCT(JL,JK)*ZE1
                   ZFLUXSQ=ZE2*ZE2*ZE1*ZCINC
                !  ZFLUXSQ=ZE2*ZE2*ZE1/ZCIN
                   ZDEP=ZACT(JL,INC,IAZI)*(ZFLUX(JL,INC,IAZI)**2-ZFLUXSQ)
                   IF(ZDEP>0.0_kind_phys) THEN
                      ZFLUX(JL,INC,IAZI)=SQRT(ZFLUXSQ)
                   ENDIF
                ENDDO
             ENDDO
         ELSEIF(GPTWO==2.0_kind_phys) THEN
             DO INC=1,INCDIM
                ZCIN=ZCI(INC)
                ZCINC=1.0_kind_phys/ZCIN
                DO JL=1,KLON
                   ZFLUXS=GCSTAR*ZFCT(JL,JK)*&
     &                    (ZCIN-ZUI(JL,JK,IAZI))**2*ZCINC
                !  ZFLUXS=GCSTAR*ZFCT(JL,JK)*(ZCIN-ZUI(JL,JK,IAZI))**2/ZCIN
                   ZDEP=ZACT(JL,INC,IAZI)*(ZFLUX(JL,INC,IAZI)-ZFLUXS)
                   IF(ZDEP>0.0_kind_phys) THEN
                      ZFLUX(JL,INC,IAZI)=ZFLUXS
                   ENDIF
                ENDDO
             ENDDO
         ENDIF

!* integrate spectrum
!* ------------------

         DO INC=1,INCDIM
            ZCINC=ZDCI(INC)
            DO JL=1,KLON
               ZPU(JL,JK,IAZI)=ZPU(JL,JK,IAZI)+&
     &               ZACT(JL,INC,IAZI)*ZFLUX(JL,INC,IAZI)*ZCINC
            ENDDO
         ENDDO

!* end JK do-loop
!* --------------

      ENDDO
!* end IAZIDIM do-loop
!* ---------------

      ENDDO

! -----------------------------------------------------------------------------

!*       MAKE CORRECTION FOR CRITICAL-LEVEL MOMENTUM DEPOSITION
!*       ------------------------------------------------------

        Z0P0=0._kind_phys
!        ZRGPTS=1.0_kind_phys/(RG*PTSTEP)
         ZRGPTS=1.0_kind_phys/(grav*PTSTEP)
        DO IAZI=1,IAZIDIM
        DO JL=1,KLON
          ZCNGL(JL)=0.0_kind_phys
        ENDDO
        DO JK=2,ILAUNCH
        DO JL=1,KLON
         ZULM=ZCOSANG(IAZI)*PUM1(JL,JK)+ZSINANG(IAZI)*& 
     &   PVM1(JL,JK)-ZUL(JL,IAZI)
         ZDFL(JL,JK-1,IAZI)=ZDFL(JL,JK-1,IAZI)+ZCNGL(JL)
         ZDFT=MIN(ZDFL(JL,JK-1,IAZI),2.0_kind_phys*(PAPM1(JL,JK-1)-&
     & PAPM1(JL,JK))*(ZCRT(JL,JK-1,IAZI)-ZULM)*ZRGPTS)

         ZDFT=MAX(ZDFT,Z0P0)
         ZCNGL(JL)=(ZDFL(JL,JK-1,IAZI)-ZDFT)
         ZPU(JL,JK,IAZI)=ZPU(JL,JK,IAZI)-ZCNGL(JL)
        ENDDO
        ENDDO
        ENDDO


!*       SUM CONTRIBUTION FOR TOTAL ZONAL AND MERIDIONAL FLUX
!*       ---------------------------------------------------

         DO IAZI=1,IAZIDIM
         DO JK=ILAUNCH,2,-1
         DO JL=1,KLON
         PFLUXU(JL,JK)=PFLUXU(JL,JK)+ZPU(JL,JK,IAZI)*ZAZ_FCT*ZCOSANG(IAZI)
         PFLUXV(JL,JK)=PFLUXV(JL,JK)+ZPU(JL,JK,IAZI)*ZAZ_FCT*ZSINANG(IAZI)
         ENDDO
         ENDDO
         ENDDO


!*    UPDATE U AND V TENDENCIES
!*    ----------------------------

!      ZCONS1=1.0_kind_phys/RCPD
       ZCONS1=rcpd
      DO JK=1,ILAUNCH
      DO JL=1, KLON
!      ZDELP= RG/(PAPHM1(JL,JK+1)-PAPHM1(JL,JK))
       ZDELP= grav/(PAPHM1(JL,JK+1)-PAPHM1(JL,JK))
      ZE1=(PFLUXU(JL,JK+1)-PFLUXU(JL,JK))*ZDELP
      ZE2=(PFLUXV(JL,JK+1)-PFLUXV(JL,JK))*ZDELP

      if (abs(ZE1) >= maxdudt ) then
           ZE1 = sign(maxdudt, ZE1) 
      endif
      if (abs(ZE2) >= maxdudt ) then
           ZE2 = sign(maxdudt, ZE2)
      endif


      PTENU(JL,JK)=ZE1
      PTENV(JL,JK)=ZE2
!     add the tendency of dT/dt 
      ZE2=-(PUM1(JL,JK)*PTENU(JL,JK)+PVM1(JL,JK)*PTENV(JL,JK))/cpd
      if (abs(ZE2) >= max_eps) pdtdt(JL,JK) = sign(max_eps, ZE2)

!     end of the tendency of dT/dt      
      ENDDO
      ENDDO


!* reverse vertical coordinate back to GFS for outbound variables only
      DO JL=1,KLON
      dked(JL,:)=transfer(dked(JL,KLEV:1:-1),dked(JL,:))
      PTENU(JL,:)=transfer(PTENU(JL,KLEV:1:-1),PTENU(JL,:))
      PTENV(JL,:)=transfer(PTENV(JL,KLEV:1:-1),PTENV(JL,:))
      pdtdt(JL,:)=transfer(pdtdt(JL,KLEV:1:-1),pdtdt(JL,:))
      ENDDO


!      DO JL=1,KLON
!      dked(JL,:)=0.0
!      PTENU(JL,:)=0.0
!      PTENV(JL,:)=0.0
!      pdtdt(JL,:)=0.0
!      ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



       return
       end subroutine ecmwf_ngw_emc

      
end module ecmwf_ngw





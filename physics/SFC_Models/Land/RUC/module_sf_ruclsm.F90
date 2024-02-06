!>\file module_sf_ruclsm.F90
!! This file is the entity of NOAA/ESRL/GSD RUC LSM Model(WRF version 4.0).

!>\ingroup lsm_ruc_group
!! This module contains the entity of the RUC LSM model, which is a  
!! soil/veg/snowpack and ice/snowpack/land-surface model to update soil
!! moisture, soil temperature, skin temperature, snowpack water content, snowdepth,
!! and all terms of the surface energy balance and surface water balance.
MODULE module_sf_ruclsm

   use machine ,   only : kind_phys, kind_dbl_prec
   use namelist_soilveg_ruc
   use physcons,   only : rhowater, con_t0c, con_hfus, con_hvap, &
                          con_pi, con_rv, con_g, con_csol, con_tice

   implicit none

   private
   !private qsn

   public :: lsmruc, ruclsminit, rslf

!> CONSTANT PARAMETERS
!! @{
      real (kind_phys), parameter :: tfrz     = con_t0c
      real (kind_phys), parameter :: xls      = con_hvap + con_hfus 
      real (kind_phys), parameter :: piconst  = con_pi
      real (kind_phys), parameter :: r_v      = con_rv
      real (kind_phys), parameter :: grav     = con_g
      real (kind_phys), parameter :: sheatice = con_csol

      real (kind_phys), parameter :: rhoice   = 917._kind_phys ! ice density
      real (kind_phys), parameter :: sheatsn  = 2090._kind_phys ! snow heat capacity
      real (kind_phys), parameter :: P1000mb  = 100000._kind_phys

      real (kind_phys), parameter :: zero     = 0._kind_dbl_prec
      real (kind_phys), parameter :: one      = 1._kind_dbl_prec

      !-- options for snow conductivity: 1 - constant, 2 - Sturm et al.,1997
      !integer, parameter :: isncond_opt = 1
      !-- Snow fraction options
      !-- option 1: original formulation using threshold snow depth to compute snow fraction
      !integer, parameter :: isncovr_opt = 1
      !-- option 2: the tanh formulation from Niu,G.-Y.,and Yang,Z.-L., 2007,JGR,DOI:10.1029/2007JD008674.
      !integer, parameter :: isncovr_opt = 2
      !-- option 3: the tanh formulation from Niu,G.-Y.,and Yang,Z with
      !   vegetation-dependent parameters from Noah MP (personal communication with
      !   Mike Barlage)
      !integer, parameter :: isncovr_opt = 3

      !-- Mosaic_lu and mosaic_soil are defined in set_soilveg_ruc.F90 and
      !   passes to RUC LSM via namelist_soilveg_ruc.F90.

!! @}

!> VEGETATION PARAMETERS
!! @{
        INTEGER :: LUCATS
        integer, PARAMETER :: NLUS=50
        CHARACTER*8 LUTYPE
!! @}

!> SOIL PARAMETERS
!! @{
        INTEGER :: SLCATS
        INTEGER, PARAMETER :: NSLTYPE=30
        CHARACTER*8 SLTYPE
!! @}

!> LSM GENERAL PARAMETERS
!! @{
        INTEGER :: SLPCATS
        INTEGER, PARAMETER :: NSLOPE=30
        real (kind_phys) ::  SBETA_DATA,FXEXP_DATA,CSOIL_DATA,SALP_DATA,REFDK_DATA,    &
                 REFKDT_DATA,FRZK_DATA,ZBOT_DATA,  SMLOW_DATA,SMHIGH_DATA, &
                        CZIL_DATA
!! @}


CONTAINS

!-----------------------------------------------------------------
!>\ingroup lsm_ruc_group
!> The RUN LSM model is described in Smirnova et al.(1997) 
!! \cite Smirnova_1997 and Smirnova et al.(2000) \cite Smirnova_2000 
!>\section gen_lsmruc GSD RUC LSM General Algorithm
!! @{
    SUBROUTINE LSMRUC(xlat,xlon,                                 &
                   DT,init,lsm_cold_start,KTAU,iter,NSL,         &
                   graupelncv,snowncv,rainncv,raincv,            &
                   ZS,RAINBL,SNOW,SNOWH,SNOWC,FRZFRAC,frpcpn,    &
                   rhosnf,precipfr,exticeden, hgt,stdev,         &
                   Z3D,P8W,T3D,QV3D,QC3D,RHO3D,EMISBCK,          &
                   GLW,GSWdn,GSW,EMISS,CHKLOWQ, CHS,             &
                   FLQC,FLHC,rhonewsn_ex,mosaic_lu,              &
                   mosaic_soil,isncond_opt,isncovr_opt,          &
                   MAVAIL,CANWAT,VEGFRA,                         &
                   ALB,ZNT,Z0,SNOALB,ALBBCK,LAI,                 & 
                   landusef, nlcat, soilctop, nscat,             &
                   smcwlt, smcref,                               & 
                   QSFC,QSG,QVG,QCG,DEW,SOILT1,TSNAV,            &
                   TBOT,IVGTYP,ISLTYP,XLAND,                     &
                   ISWATER,ISICE,XICE,XICE_THRESHOLD,            &
                   CP,RV,RD,G0,PI,LV,STBOLT,                     &
                   SOILMOIS,SH2O,SMAVAIL,SMMAX,                  &
                   TSO,SOILT,EDIR,EC,ETT,SUBLIM,SNOH,            &
                   HFX,QFX,LH,INFILTR,                           &
                   RUNOFF1,RUNOFF2,ACRUNOFF,SFCEXC,              &
                   SFCEVP,GRDFLX,SNOWFALLAC,ACSNOW,SNOM,         &
                   SMFR3D,KEEPFR3DFLAG,                          &
                   add_fire_heat_flux,fire_heat_flux,            &
                   myj,shdmin,shdmax,rdlai2d,                    &
                   ims,ime, jms,jme, kms,kme,                    &
                   its,ite, jts,jte, kts,kte,                    &
                   errmsg, errflg)
!-----------------------------------------------------------------
   IMPLICIT NONE
!-----------------------------------------------------------------
!
! The RUC LSM model is described in:
!  Smirnova, T.G., J.M. Brown, and S.G. Benjamin, 1997: 
!     Performance of different soil model configurations in simulating 
!     ground surface temperature and surface fluxes. 
!     Mon. Wea. Rev. 125, 1870-1884.
!  Smirnova, T.G., J.M. Brown, and D. Kim, 2000: Parameterization of 
!     cold-season processes in the MAPS land-surface scheme. 
!     J. Geophys. Res. 105, 4077-4086.
!-----------------------------------------------------------------
!-- DT            time step (second)
!        init - flag for initialization
!lsm_cold_start - flag for cold start run
!        ktau - number of time step
!        NSL  - number of soil layers
!        NZS  - number of levels in soil
!        ZS   - depth of soil levels (m)
!-- RAINBL    - accumulated rain in [mm] between the PBL calls
!-- RAINNCV         one time step grid scale precipitation (mm/step)
!-- RAINCV          one time step convective precipitation (mm/step)
!        SNOW - snow water equivalent [mm]
!        FRAZFRAC - fraction of frozen precipitation
!-- PRECIPFR (mm) - time step frozen precipitation
!-- SNOWC       flag indicating snow coverage (1 for snow cover)
!-- Z3D         heights (m)
!-- P8W         3D pressure (Pa)
!-- T3D         temperature (K)
!-- QV3D        3D water vapor mixing ratio (Kg/Kg)
!        QC3D - 3D cloud water mixing ratio (Kg/Kg)
!       RHO3D - 3D air density (kg/m^3)
!-- GLW         downward long wave flux at ground surface (W/m^2)
!-- GSW         absorbed short wave flux at ground surface (W/m^2)
!-- EMISS       surface emissivity (between 0 and 1)
!        FLQC - surface exchange coefficient for moisture (kg/m^2/s)
!        FLHC - surface exchange coefficient for heat [W/m^2/s/degreeK]     
!      SFCEXC - surface exchange coefficient for heat [m/s]
!      CANWAT - CANOPY MOISTURE CONTENT (mm)
!      VEGFRA - vegetation fraction (between 0 and 100)
!         ALB - surface albedo (between 0 and 1)
!      SNOALB - maximum snow albedo (between 0 and 1)
!      ALBBCK - snow-free albedo (between 0 and 1)
!         ZNT - roughness length [m]
!-- TBOT        soil temperature at lower boundary (K)
!      IVGTYP - USGS vegetation type (24 classes)
!      ISLTYP - STASGO soil type (16 classes)
!-- XLAND       land mask (1 for land, 2 for water)
!-- CP          heat capacity at constant pressure for dry air (J/kg/K)
!-- G0          acceleration due to gravity (m/s^2)
!-- LV          latent heat of evaporation (J/kg)
!-- STBOLT      Stefan-Boltzmann constant (W/m^2/K^4)
!    SOILMOIS - soil moisture content (volumetric fraction)
!         TSO - soil temp (K)
!-- SOILT       surface temperature (K)
!-- HFX         upward heat flux at the surface (W/m^2)
!-- QFX         upward moisture flux at the surface (kg/m^2/s)
!-- LH          upward latent heat flux (W/m^2)
!   SFCRUNOFF - ground surface runoff [mm]
!   UDRUNOFF - underground runoff [mm]
!   ACRUNOFF - run-total surface runoff [mm]
!   SFCEVP - total time-step evaporation in [kg/m^2]
!   GRDFLX - soil heat flux (W/m^2: negative, if downward from surface)
!   SNOWFALLAC - run-total snowfall accumulation [mm]   
!   ACSNOW - run-toral SWE of snowfall [mm]   
!-- CHKLOWQ - is either 0 or 1 (so far set equal to 1).
!--           used only in MYJPBL. 
!-- tice - sea ice temperture (C)
!-- rhosice - sea ice density (kg m^-3)
!-- capice - sea ice volumetric heat capacity (J/m^3/K)
!-- thdifice - sea ice thermal diffusivity (m^2/s)
!--
!-- ims           start index for i in memory
!-- ime           end index for i in memory
!-- jms           start index for j in memory
!-- jme           end index for j in memory
!-- kms           start index for k in memory
!-- kme           end index for k in memory
!-------------------------------------------------------------------------
!   INTEGER,     PARAMETER            ::     nzss=5
!   INTEGER,     PARAMETER            ::     nddzs=2*(nzss-2)

   real (kind_phys),       INTENT(IN   )    ::     xlat,xlon
   real (kind_phys),       INTENT(IN   )    ::     DT
   LOGICAL,    INTENT(IN   )    ::     myj,frpcpn,init,lsm_cold_start,exticeden
   INTEGER,    INTENT(IN   )    ::     NLCAT, NSCAT 
   INTEGER,    INTENT(IN   )    ::     mosaic_lu,mosaic_soil
   INTEGER,    INTENT(IN   )    ::     isncond_opt,isncovr_opt
   INTEGER,    INTENT(IN   )    ::     ktau, iter, nsl, isice, iswater, &
                                       ims,ime, jms,jme, kms,kme,       &
                                       its,ite, jts,jte, kts,kte

!   LOGICAL,    DIMENSION( ims:ime, jms:jme ), INTENT(IN   )    ::     flag_iter, flag

   real (kind_phys),   DIMENSION( ims:ime, kms:kme, jms:jme ) , &
            INTENT(IN   )    ::                           QV3D, &
                                                          QC3D, &
                                                           p8w, &
                                                         rho3D, &
                                                           T3D, &
                                                           z3D

   real (kind_phys),   DIMENSION( ims:ime , jms:jme ),           &
               INTENT(IN   )    ::                       RAINBL, &
                                                            GLW, &
                                                          GSWdn, &
                                                            GSW, &
                                                         ALBBCK, &
                                                           FLHC, &
                                                           FLQC, &
                                                           CHS , &
                                                           XICE, &
                                                          XLAND, &
                                                         VEGFRA, &
                                                           TBOT

   real (kind_phys),       DIMENSION( ims:ime , jms:jme ),       &
               INTENT(IN   )    ::                   GRAUPELNCV, &
                                                        SNOWNCV, &
                                                         RAINCV, &
                                                        RAINNCV
   real (kind_phys),      DIMENSION( ims:ime),  INTENT(IN   )  ::   RHONEWSN_ex     !externally-calculated srf frz precip density

   real (kind_phys), DIMENSION( ims:ime , jms:jme ), INTENT(IN )::   SHDMAX
   real (kind_phys), DIMENSION( ims:ime , jms:jme ), INTENT(IN )::   SHDMIN
   real (kind_phys), DIMENSION( ims:ime , jms:jme ), INTENT(IN )::   hgt
   real (kind_phys), DIMENSION( ims:ime , jms:jme ), INTENT(IN )::   stdev
   LOGICAL, intent(in) :: add_fire_heat_flux
   real (kind_phys), DIMENSION( ims:ime , jms:jme ), INTENT(IN ):: fire_heat_flux
   LOGICAL, intent(in) :: rdlai2d

   real (kind_phys),       DIMENSION( 1:nsl), INTENT(IN   )  :: ZS

   real (kind_phys),       DIMENSION( ims:ime , jms:jme ),       &
               INTENT(INOUT)    ::                               &
                                                           SNOW, &
                                                          SNOWH, &
                                                          SNOWC, &
                                                         CANWAT, &
                                                         SNOALB, &
                                                            ALB, &
                                                            LAI, &
                                                         SMCWLT, &
                                                         SMCREF, &
                                                          EMISS, &
                                                        EMISBCK, &
                                                         MAVAIL, & 
                                                         SFCEXC, &
                                                            Z0 , &
                                                            ZNT

   real (kind_phys),       DIMENSION( ims:ime , jms:jme ),       &
               INTENT(IN   )    ::                               &
                                                        FRZFRAC

   INTEGER,    DIMENSION( ims:ime , jms:jme ),                   &
               INTENT(IN   )    ::                       IVGTYP, &
                                                         ISLTYP
   real (kind_phys),     DIMENSION( ims:ime , 1:nlcat, jms:jme ), INTENT(IN):: LANDUSEF
   real (kind_phys),     DIMENSION( ims:ime , 1:nscat, jms:jme ), INTENT(IN):: SOILCTOP

   real (kind_phys), INTENT(IN   ) ::  CP,G0,LV,STBOLT,RV,RD,PI, &
                                         XICE_threshold
 
   real (kind_phys),     DIMENSION( ims:ime , 1:nsl, jms:jme ) , &
               INTENT(INOUT)    ::                 SOILMOIS,SH2O,TSO

   real (kind_phys),       DIMENSION( ims:ime, jms:jme )       , &
               INTENT(INOUT)    ::                        SOILT, &
                                                            HFX, &
                                                            QFX, &
                                                             LH, &
                                                           EDIR, &
                                                             EC, &
                                                            ETT, &
                                                         SUBLIM, &
                                                           SNOH, &
                                                         SFCEVP, &
                                                        RUNOFF1, &
                                                        RUNOFF2, &
                                                       ACRUNOFF, &
                                                         GRDFLX, &
                                                         ACSNOW, &
                                                           SNOM, &
                                                            QVG, &
                                                            QCG, &
                                                            DEW, &
                                                           QSFC, &
                                                            QSG, &
                                                        CHKLOWQ, &
                                                         SOILT1, &
                                                          TSNAV

   real (kind_phys),       DIMENSION( ims:ime, jms:jme )       , & 
               INTENT(INOUT)    ::                      SMAVAIL, &
                                                          SMMAX

   real (kind_phys),       DIMENSION( its:ite, jts:jte )    ::   &
                                                             PC, &
                                                      SFCRUNOFF, &
                                                       UDRUNOFF, &
                                                         EMISSL, &
                                                           MSNF, &
                                                         FACSNF, &
                                                           ZNTL, &
                                                        LMAVAIL, &
                                                          SMELT, &
                                                          SNFLX, &
                                                           sflx, &
                                                            smf, &
                                                          EVAPL, &
                                                          PRCPL, &
                                                         SEAICE, &
                                                        INFILTR
! Energy and water budget variables:
   real (kind_phys),       DIMENSION( its:ite, jts:jte )    ::   &
                                                         budget, &
                                                       acbudget, &
                                                    waterbudget, &
                                                  acwaterbudget, &
                                                       smtotold, &
                                                        snowold, &
                                                      canwatold


   real (kind_phys),       DIMENSION( ims:ime, 1:nsl, jms:jme)   &
                                             ::    KEEPFR3DFLAG, &
                                                         SMFR3D

   real (kind_phys),DIMENSION( ims:ime, jms:jme ),INTENT(OUT) :: &
                                                         RHOSNF, & ! RHO of snowfall
                                                       PRECIPFR, & ! time-step frozen precip
                                                     SNOWFALLAC
!--- soil/snow properties
   real (kind_phys)                                              &
                             ::                           RHOCS, &
                                                       RHONEWSN, &
                                                          RHOSN, &
                                                      RHOSNFALL, &
                                                           BCLH, &
                                                            DQM, &
                                                           KSAT, &
                                                           PSIS, &
                                                           QMIN, &
                                                          QWRTZ, &
                                                            REF, &
                                                           WILT, &
                                                        CANWATR, &
                                                       SNOWFRAC, &
                                                          SNHEI, &
                                                           SNWE

   real (kind_phys)                                      ::  CN, &
                                                         SAT,CW, &
                                                           C1SN, &
                                                           C2SN, &
                                                         KQWRTZ, &
                                                           KICE, &
                                                            KWT


   real (kind_phys),     DIMENSION(1:NSL)             :: ZSMAIN, &
                                                         ZSHALF, &
                                                         DTDZS2

   real (kind_phys),     DIMENSION(1:2*(nsl-2))       :: DTDZS

   real (kind_phys),     DIMENSION(1:5001)            :: TBQ


   real (kind_phys),     DIMENSION( 1:nsl )          :: SOILM1D, & 
                                                          TSO1D, &
                                                        SOILICE, &
                                                        SOILIQW, &
                                                       SMFRKEEP

   real (kind_phys),     DIMENSION( 1:nsl )           :: KEEPFR
                                                
   real (kind_phys),     DIMENSION( 1:nlcat )         :: lufrac
   real (kind_phys),     DIMENSION( 1:nscat )         :: soilfrac

   real (kind_phys)                           ::            RSM, &
                                                      SNWEPRINT, &
                                                      SNHEIPRINT

   real (kind_phys)                           ::         PRCPMS, &
                                                        NEWSNMS, &
                                                      prcpncliq, &
                                                       prcpncfr, &
                                                      prcpculiq, &
                                                       prcpcufr, &
                                                           PATM, &
                                                          PATMB, &
                                                           TABS, &
                                                          QVATM, &
                                                          QCATM, &
                                                          Q2SAT, &
                                                         CONFLX, &
                                                            RHO, &
                                                           QKMS, &
                                                           TKMS, &
                                                        snowrat, &
                                                       grauprat, &
                                                         icerat, &
                                                          curat, &
                                                       INFILTRP
   real (kind_phys)      ::  cq,r61,r273,arp,brp,x,evs,eis
   real (kind_phys)      ::  cropfr, cropsm, newsm, factor

   real (kind_phys)      ::  meltfactor, ac,as, wb,rovcp
   INTEGER   ::  NROOT
   INTEGER   ::  ILAND,ISOIL,IFOREST
 
   INTEGER   ::  I,J,K,NZS,NZS1,NDDZS
   INTEGER   ::  k1,k2
   logical :: debug_print

   !-- diagnostic point
   real (kind_phys) :: testptlat, testptlon

   character(len=*),intent(out) :: errmsg
   integer,         intent(out) :: errflg

!-----------------------------------------------------------------
!   
     ! Initialize error-handling
     errflg = 0
     errmsg = ''

     debug_print = .false.
!
         rovcp = rd/cp

         NZS=NSL
         NDDZS=2*(nzs-2)

        !--
        testptlat = 35.55 !48.7074_kind_phys !39.958 !42.05 !39.0 !74.12 !29.5 
        testptlon = 278.66 !289.03_kind_phys !271.622 !286.75 !280.6 !164.0 !283.0 
        !--


!> - Table TBQ is for resolution of balance equation in vilka()
        CQ=173.15_kind_dbl_prec-.05_kind_dbl_prec
        R273=1._kind_dbl_prec/tfrz
        R61=6.1153_kind_dbl_prec*0.62198_kind_dbl_prec
        ARP=77455._kind_dbl_prec*41.9_kind_dbl_prec/461.525_kind_dbl_prec
        BRP=64._kind_dbl_prec*41.9_kind_dbl_prec/461.525_kind_dbl_prec

        DO K=1,5001
          CQ=CQ+.05_kind_dbl_prec
          EVS=EXP(17.67_kind_dbl_prec*(CQ-tfrz)/(CQ-29.65_kind_dbl_prec))
          EIS=EXP(22.514_kind_dbl_prec-6.15E3_kind_dbl_prec/CQ)
          if(CQ.ge.tfrz) then
          ! tbq is in mb
            tbq(k) = R61*evs
          else
            tbq(k) = R61*eis
          endif
        END DO

!> - Initialize soil/vegetation parameters
!--- This is temporary until SI is added to mass coordinate ---!!!!!

    if(init .and. iter == 1) then

     if( lsm_cold_start ) then
     !-- beginning of cold-start
       DO J=jts,jte
         DO i=its,ite
!
!>  - Initializing inside-snow temp if it is not defined
           IF((soilt1(i,j) .LT. 170._kind_phys) .or. (soilt1(i,j) .GT.400._kind_phys)) THEN
             IF(snowc(i,j).gt.zero) THEN
               soilt1(i,j)=min(tfrz,0.5_kind_phys*(soilt(i,j)+tso(i,1,j)) )
               IF (debug_print ) THEN
                   print *, &
                  'Temperature inside snow is initialized in RUCLSM ', soilt1(i,j),i,xlat,xlon
               ENDIF
             ELSE
               soilt1(i,j) = tso(i,1,j)
             ENDIF
           ENDIF
           tsnav(i,j) =min(zero,0.5_kind_phys*(soilt(i,j)+tso(i,1,j))-tfrz)
           !- 10feb22 - limit snow albedo at high elevations
           !- based on Roesch et al., Climate Dynamics (2001),17:933-946
           if(hgt(i,j) > 2500._kind_phys) then
             snoalb(i,j) = min(0.65_kind_phys,snoalb(i,j))
           endif

           patmb=P8w(i,kms,j)*1.e-2_kind_phys
           QSG  (i,j) = QSN(SOILT(i,j),TBQ)/PATMB
           
           if((qcg(i,j) < zero) .or. (qcg(i,j) > 0.1_kind_phys)) then
             qcg  (i,j) = qc3d(i,1,j)
             if (debug_print ) then
               print *, 'QCG is initialized in RUCLSM ', qcg(i,j),qc3d(i,1,j),i,xlat,xlon
             endif
           endif

           if((qvg(i,j) .LE. zero) .or. (qvg(i,j) .GT.0.1_kind_phys)) then
             qvg  (i,j) = qv3d(i,1,j)
             if (debug_print ) then
               print *, 'QVG is initialized in RUCLSM ', qvg(i,j),mavail(i,j),qsg(i,j),i,xlat,xlon
             endif
           endif
           qsfc(i,j) = qvg(i,j)/(1.+qvg(i,j))

           SMELT(i,j) = zero
           SNOM (i,j) = zero
           ACSNOW(i,j) = zero
           SNOWFALLAC(i,j) = zero
           PRECIPFR(i,j) = zero
           RHOSNF(i,j) = -1.e3_kind_phys ! non-zero flag
           SNFLX(i,j) = zero
           DEW  (i,j) = zero
           PC   (i,j) = zero
           zntl (i,j) = zero
           RUNOFF1(i,j) = zero
           RUNOFF2(i,j) = zero
           SFCRUNOFF(i,j) = zero
           UDRUNOFF(i,j) = zero
           ACRUNOFF(i,j) = zero
           emissl (i,j) = zero
           msnf (i,j) = zero
           facsnf (i,j) = zero
           budget(i,j) = zero
           acbudget(i,j) = zero
           waterbudget(i,j) = zero
           acwaterbudget(i,j) = zero
           smtotold(i,j)=zero
           canwatold(i,j)=zero

!>  - For RUC LSM CHKLOWQ needed for MYJPBL should 
!! 1 because is actual specific humidity at the surface, and
!! not the saturation value
           chklowq(i,j) = one
           infiltr(i,j) = zero 
           snoh  (i,j) = zero
           edir  (i,j) = zero
           ec    (i,j) = zero
           ett   (i,j) = zero
           sublim(i,j) = zero
           sflx  (i,j) = zero
           smf   (i,j) = zero
           evapl (i,j) = zero
           prcpl (i,j) = zero
         ENDDO
       ENDDO

       infiltrp = zero
       do k=1,nsl
         soilice(k)=zero
         soiliqw(k)=zero
       enddo
      endif ! cold start
     endif ! init==.true.

!-----------------------------------------------------------------

        PRCPMS = zero
        newsnms = zero
        prcpncliq = zero
        prcpculiq = zero
        prcpncfr = zero
        prcpcufr = zero

   DO J=jts,jte

      DO i=its,ite

    IF (debug_print ) THEN
       if (abs(xlat-testptlat).lt.0.2 .and.                         &
           abs(xlon-testptlon).lt.0.2)then
         print 100,'(RUC start)  i=',i,'  lat,lon=',xlat,xlon,      &
         'mavail ', mavail(i,j),' soilt',soilt(i,j),'qvg ',qvg(i,j),&
         'p8w',p8w(i,1,j),'sflay qfx',qfx(i,j),'sflay hfx',hfx(i,j),&
         'gsw ',gsw(i,j),'glw ',glw(i,j),'soilt ',soilt(i,j),       &
         'chs ',chs(i,j),'flqc ',flhc(i,j),'alb ',alb(i,j),         &
         'rainbl ',rainbl(i,j),'dt ',dt
         print *,'nzs',nzs, 'ivgtyp ',ivgtyp(i,j),'isltyp ',isltyp(i,j)
      endif
    ENDIF

         ILAND     = IVGTYP(i,j)
         ISOIL     = ISLTYP(I,J)
         TABS      = T3D(i,kms,j)
         QVATM     = QV3D(i,kms,j)
         QCATM     = QC3D(i,kms,j)
         PATM      = P8w(i,kms,j)*1.e-5_kind_phys
!> - Z3D(1) is thickness between first full sigma level and the surface, 
!! but first mass level is at the half of the first sigma level 
!! (u and v are also at the half of first sigma level)
         CONFLX    = Z3D(i,kms,j)*0.5_kind_phys
         RHO       = RHO3D(I,kms,J)
!> - Initialize snow, graupel and ice fractions in frozen precip
         snowrat = zero
         grauprat = zero
         icerat = zero
         curat = zero
       IF(FRPCPN) THEN
         prcpncliq = rainncv(i,j)*(1.-frzfrac(i,j))
         prcpncfr = rainncv(i,j)*frzfrac(i,j)
!> - Apply the same frozen precipitation fraction to convective precip
!tgs - 31 mar17 - add temperature check in case Thompson MP produces
!                 frozen precip at T > 273.
       if(frzfrac(i,j) > zero .and. tabs < tfrz) then
         prcpculiq = max(zero,raincv(i,j)*(one-frzfrac(i,j)))
         prcpcufr = max(zero,raincv(i,j)*frzfrac(i,j))
       else
          if(tabs < tfrz) then
            prcpcufr = max(zero,raincv(i,j))
            prcpculiq = zero
          else
            prcpcufr = zero
            prcpculiq = max(zero,raincv(i,j))
          endif  ! tabs < 273.
       endif  ! frzfrac > 0.
!--- 1*e-3 is to convert from mm/s to m/s
         PRCPMS   = (prcpncliq + prcpculiq)/DT*1.e-3_kind_phys
         NEWSNMS  = (prcpncfr + prcpcufr)/DT*1.e-3_kind_phys

         if((prcpncfr + prcpcufr) > zero) then
!> - Calculate snow, graupel and ice fractions in falling frozen precip
         snowrat=min(one,max(zero,snowncv(i,j)/(prcpncfr + prcpcufr)))
         grauprat=min(one,max(zero,graupelncv(i,j)/(prcpncfr + prcpcufr)))
         icerat=min(one,max(zero,(prcpncfr-snowncv(i,j)-graupelncv(i,j)) &
               /(prcpncfr + prcpcufr)))
         curat=min(one,max(zero,(prcpcufr/(prcpncfr + prcpcufr))))
         endif

       ELSE  ! .not. FRPCPN
          if (tabs.le.tfrz) then
         PRCPMS    = zero 
         NEWSNMS   = RAINBL(i,j)/DT*1.e-3_kind_phys
         !> - If here no info about constituents of frozen precipitation,
         !! suppose it is all snow
         snowrat = one 
          else
         PRCPMS    = RAINBL(i,j)/DT*1.e-3_kind_phys
         NEWSNMS   = zero 
          endif
       ENDIF

! -- save time-step water equivalent of frozen precipitation in PRECIPFR array to be used in
!    module_diagnostics
          precipfr(i,j) = NEWSNMS * DT *1.e3_kind_phys

        if   (myj)   then
         QKMS=CHS(i,j)
         TKMS=CHS(i,j)
        else
!> - Convert exchange coeff QKMS to [m/s]
         QKMS=FLQC(I,J)/RHO/MAVAIL(I,J)
!         TKMS=FLHC(I,J)/RHO/CP
         TKMS=FLHC(I,J)/RHO/(CP*(one+0.84_kind_phys*QVATM))  ! mynnsfc uses CPM
        endif
!> - Convert incoming snow and canwat from mm to m
         SNWE=SNOW(I,J)*1.E-3_kind_phys
         SNHEI=SNOWH(I,J)
         CANWATR=CANWAT(I,J)*1.E-3_kind_phys

         SNOWFRAC=SNOWC(I,J)
         RHOSNFALL=RHOSNF(I,J)

         snowold(i,j)=snwe
!-----
             zsmain(1)=zero
             zshalf(1)=zero
          do k=2,nzs
             zsmain(k)= zs(k)
             zshalf(k)=0.5_kind_phys*(zsmain(k-1) + zsmain(k))
          enddo

          do k=1,nlcat
             lufrac(k) = landusef(i,k,j)
          enddo
          do k=1,nscat
             soilfrac(k) = soilctop(i,k,j)
          enddo

!------------------------------------------------------------
!-----  DDZS and DSDZ1 are for implicit solution of soil eqns.
!-------------------------------------------------------------
        NZS1=NZS-1
!-----
    IF (debug_print ) THEN
          print *,' DT,NZS1, ZSMAIN, ZSHALF --->', dt,nzs1,zsmain,zshalf
    ENDIF

        DO  K=2,NZS1
          K1=2*K-3
          K2=K1+1
          X=DT/2./(ZSHALF(K+1)-ZSHALF(K))
          DTDZS(K1)=X/(ZSMAIN(K)-ZSMAIN(K-1))
          DTDZS2(K-1)=X
          DTDZS(K2)=X/(ZSMAIN(K+1)-ZSMAIN(K))
        END DO

        CW =4.183E6_kind_dbl_prec


!--- Constants used in Johansen soil thermal
!--- conductivity method

        KQWRTZ=7.7_kind_dbl_prec
        KICE=2.2_kind_dbl_prec
        KWT=0.57_kind_dbl_prec

!***********************************************************************
!--- Constants for snow density calculations C1SN and C2SN

        c1sn=0.026_kind_dbl_prec
        c2sn=21._kind_dbl_prec

!***********************************************************************

        NROOT= 4
!           ! rooting depth

        RHONEWSN = 200._kind_phys
       if(SNOW(i,j).gt.zero .and. SNOWH(i,j).gt.0.02_kind_phys) then
        RHOSN = SNOW(i,j)/SNOWH(i,j)
       else
        RHOSN = 300._kind_phys
       endif

    IF (debug_print ) THEN
      if(init) then
        if (abs(xlat-testptlat).lt.0.2 .and.                           &
            abs(xlon-testptlon).lt.0.2)then
           print*,'  lat,lon=',xlat,xlon
           print *,'before SOILVEGIN - z0,znt',i,z0(i,j),znt(i,j)
           print *,'ILAND, ISOIL =',i,iland,isoil
        endif
      endif
    ENDIF
 
!> - Call soilvegin() to initialize soil and surface properties
     !-- land or ice
       CALL SOILVEGIN  ( debug_print, mosaic_lu, mosaic_soil,                        &
                       soilfrac,nscat,shdmin(i,j),shdmax(i,j),                       &
                       NLCAT,ILAND,ISOIL,iswater,MYJ,IFOREST,lufrac,VEGFRA(I,J),     &
                       EMISSL(I,J),PC(I,J),MSNF(I,J),FACSNF(I,J),                    &
                       ZNT(I,J),LAI(I,J),RDLAI2D,                                    &
                       QWRTZ,RHOCS,BCLH,DQM,KSAT,PSIS,QMIN,REF,WILT,i,j,errmsg, errflg)

       !-- update background emissivity for land points, can have vegetation mosaic effect
       EMISBCK(I,J) = EMISSL(I,J)
       smcwlt(i,j)  = wilt
       smcref(i,j)  = ref

    IF (debug_print ) THEN
      if(init)then
        if (abs(xlat-testptlat).lt.0.2 .and.                           &
            abs(xlon-testptlon).lt.0.2)then
         print*,'  lat,lon=',xlat,xlon
         print *,'after SOILVEGIN - z0,znt,lai',i,z0(i,j),znt(i,j),lai(i,j)
         print *,'NLCAT,iland,EMISSL(I,J),PC(I,J),ZNT(I,J),LAI(I,J)', &
                  NLCAT,iland,EMISSL(I,J),PC(I,J),ZNT(I,J),LAI(I,J),i,j
         print *,'NSCAT,QWRTZ,RHOCS,BCLH,DQM,KSAT,PSIS,QMIN,REF,WILT',&
                 NSCAT,QWRTZ,RHOCS,BCLH,DQM,KSAT,PSIS,QMIN,REF,WILT,i,j
        endif
      endif
    ENDIF

        CN=CFACTR_DATA   ! exponent
        SAT = 5.e-4_kind_phys  ! units [m]

!-- definition of number of soil levels in the rooting zone
     IF(iforest.gt.2) THEN
!---- all vegetation types except evergreen and mixed forests
!18apr08 - define meltfactor for Egglston melting limit:
! for open areas factor is 2, and for forests - factor is 0.85
! This will make limit on snow melting smaller and let snow stay 
! longer in the forests.
         meltfactor = 2.0_kind_phys

         do k=2,nzs
         if(zsmain(k).ge.0.4_kind_phys) then
            NROOT=K
            goto  111
         endif
         enddo
     ELSE
!---- evergreen and mixed forests
!18apr08 - define meltfactor
!         meltfactor = 1.5
! 28 March 11 - Previously used value of metfactor= 1.5 needs to be further reduced 
! to compensate for low snow albedos in the forested areas. 
! Melting rate in forests will reduce.
         meltfactor = 0.85_kind_phys

         do k=2,nzs
         if(zsmain(k).ge.1.1_kind_phys) then
            NROOT=K
            goto  111
         endif
         enddo
     ENDIF
 111   continue

!-----
    IF (debug_print ) THEN
      if (abs(xlat-testptlat).lt.0.2 .and.                               &
          abs(xlon-testptlon).lt.0.2)then
         print*,'  lat,lon=',xlat,xlon
         print *,' ZNT, LAI, VEGFRA, SAT, EMIS, PC --->',                &
                   ZNT(I,J),LAI(I,J),VEGFRA(I,J),SAT,EMISSL(I,J),PC(I,J)
         print *,' ZS, ZSMAIN, ZSHALF, CONFLX, CN, SAT, --->', zs,zsmain,zshalf,conflx,cn,sat
         print *,'NROOT, meltfactor, iforest, ivgtyp, i,j ', nroot,meltfactor,iforest,ivgtyp(I,J),I,J
      endif
    ENDIF

        IF((XLAND(I,J)-1.5).GE.0._kind_phys)THEN
!-- Water 
           SMAVAIL(I,J)= one
             SMMAX(I,J)= one
             SNOW(I,J) = zero
             SNOWH(I,J)= zero
             SNOWC(I,J)= zero
           LMAVAIL(I,J)= one
! accumulated water equivalent of frozen precipitation over water [mm]
           acsnow(i,j)=acsnow(i,j)+precipfr(i,j)

           ILAND=iswater
           ISOIL=14

           patmb=P8w(i,1,j)*1.e-2_kind_phys
           qvg  (i,j) = QSN(SOILT(i,j),TBQ)/PATMB
           qsfc(i,j) = qvg(i,j)/(1.+qvg(i,j))
           CHKLOWQ(I,J)= one
           Q2SAT=QSN(TABS,TBQ)/PATMB

            DO K=1,NZS
              SOILMOIS(I,K,J)=one
              SH2O    (I,K,J)=one 
              TSO(I,K,J)= SOILT(I,J)
            ENDDO

    IF (debug_print ) THEN
      if (abs(xlat-testptlat).lt.0.2 .and.                    &
            abs(xlon-testptlon).lt.0.2)then
           PRINT*,'  water point'  
           print*,'  lat,lon=',xlat,xlon,'SOILT=', SOILT(i,j)
      endif
    ENDIF

           ELSE

! LAND POINT OR SEA ICE
       if(xice(i,j).ge.xice_threshold) then
           SEAICE(i,j)=one
       else
           SEAICE(i,j)=zero
       endif

         IF(SEAICE(I,J).GT.0.5_kind_phys)THEN
!-- Sea-ice case
    IF (debug_print ) THEN
        if (abs(xlat-testptlat).lt.0.2 .and.                    &      
            abs(xlon-testptlon).lt.0.2)then
           PRINT*,' sea-ice at water point'
           print*,'  lat,lon=',xlat,xlon
        endif
    ENDIF
            ILAND = isice
        if(nscat == 9) then
            ISOIL = 9  ! ZOBLER
        else
            ISOIL = 16 ! STATSGO
        endif
            ZNT(I,J) = 0.011_kind_phys
            ! in FV3 albedo and emiss are defined for ice
            emissl(i,j) = emisbck(i,j) ! no snow impact, old 0.98 used in WRF 
            dqm = one
            ref = one
            qmin = zero
            wilt = zero

           patmb=P8w(i,1,j)*1.e-2_kind_phys
           qvg  (i,j) = QSN(SOILT(i,j),TBQ)/PATMB
           qsg  (i,j) = qvg(i,j)
           qsfc(i,j) = qvg(i,j)/(1.+qvg(i,j))

            DO K=1,NZS
               soilmois(i,k,j) = one
               smfr3d(i,k,j)   = one
               sh2o(i,k,j)     = zero
               keepfr3dflag(i,k,j) = zero
               tso(i,k,j) = min(con_tice,tso(i,k,j))
            ENDDO
          ENDIF

!  Attention!!!!  RUC LSM uses soil moisture content minus residual (minimum
!  or dry soil moisture content for a given soil type) as a state variable.

           DO k=1,nzs
           ! soilm1d - soil moisture content minus residual [m**3/m**3]
              soilm1d (k) = min(max(zero,soilmois(i,k,j)-qmin),dqm)
              tso1d   (k) = tso(i,k,j)
              soiliqw (k) = min(max(zero,sh2o(i,k,j)-qmin),soilm1d(k))
              soilice (k) =(soilm1d (k) - soiliqw (k))/0.9_kind_phys 
           ENDDO 

           do k=1,nzs
              smfrkeep(k) = smfr3d(i,k,j)
              keepfr  (k) = keepfr3dflag(i,k,j)
           enddo

              LMAVAIL(I,J)=max(0.00001_kind_phys,min(one,soilm1d(1)/(ref-qmin)))

    IF (debug_print ) THEN
      if (abs(xlat-testptlat).lt.0.2 .and.                           &
          abs(xlon-testptlon).lt.0.2)then
        print*,'  lat,lon=',xlat,xlon
        print *,'LAND, i,j,tso1d,soilm1d,PATM,TABS,QVATM,QCATM,RHO',  &
                       i,j,tso1d,soilm1d,PATM,TABS,QVATM,QCATM,RHO
        print *,'CONFLX =',CONFLX 
        print *,'SMFRKEEP,KEEPFR   ',SMFRKEEP,KEEPFR
      endif
    ENDIF

        smtotold(i,j)=0.

      do k=1,nroot
        smtotold(i,j)=smtotold(i,j)+(qmin+soilm1d(k))*             &
                    (zshalf(k+1)-zshalf(k))
      enddo

       if (debug_print .and. abs(xlat-testptlat).lt.0.2          &
          .and. abs(xlon-testptlon).lt.0.2) then
         print *,'Old soilm1d ',i,soilm1d
       endif

        canwatold(i,j) = canwatr
!-----------------------------------------------------------------
         CALL SFCTMP (debug_print, dt,ktau,conflx,i,j,           &
                xlat, xlon, testptlat, testptlon,                &
!--- input variables
                nzs,nddzs,nroot,meltfactor,                      &   !added meltfactor
                isncond_opt,isncovr_opt,                         &
                iland,isoil,ivgtyp(i,j),isltyp(i,j),             &
                PRCPMS, NEWSNMS,SNWE,SNHEI,SNOWFRAC,             &
                exticeden,RHOSN,RHONEWSN_ex(I),RHONEWSN,         &
                RHOSNFALL,snowrat,grauprat,icerat,curat,         &
                PATM,TABS,QVATM,QCATM,RHO,                       &
                GLW(I,J),GSWdn(i,j),GSW(I,J),                    &
                EMISSL(I,J),EMISBCK(I,J),                        &
                msnf(i,j), facsnf(i,j),                          &
                QKMS,TKMS,PC(I,J),LMAVAIL(I,J),                  &
                canwatr,vegfra(I,J),alb(I,J),znt(I,J),           &
                snoalb(i,j),albbck(i,j),lai(i,j),                &
                hgt(i,j),stdev(i,j),                             &   !new
                myj,seaice(i,j),isice,                           &
                add_fire_heat_flux,fire_heat_flux(i,j),          &
!--- soil fixed fields
                QWRTZ,                                           &
                rhocs,dqm,qmin,ref,                              &
                wilt,psis,bclh,ksat,                             &
                sat,cn,zsmain,zshalf,DTDZS,DTDZS2,tbq,           &
!--- constants
                cp,rovcp,g0,lv,stbolt,cw,c1sn,c2sn,              &
                KQWRTZ,KICE,KWT,                                 &
!--- output variables
                snweprint,snheiprint,rsm,                        &
                soilm1d,tso1d,smfrkeep,keepfr,                   &
                soilt(I,J),soilt1(i,j),tsnav(i,j),dew(I,J),      &
                qvg(I,J),qsg(I,J),qcg(I,J),SMELT(I,J),           &
                SNOH(I,J),SNFLX(I,J),SNOM(I,J),SNOWFALLAC(I,J),  &
                ACSNOW(I,J),edir(I,J),ec(I,J),ett(I,J),qfx(I,J), &
                lh(I,J),hfx(I,J),sflx(I,J),sublim(I,J),          &
                evapl(I,J),prcpl(I,J),budget(i,j),runoff1(i,j),  &
                runoff2(I,J),soilice,soiliqw,infiltrp,smf(i,j))
!-----------------------------------------------------------------

! Fraction of cropland category in the grid box should not have soil moisture below 
! wilting point during the growing season.
! Let's keep soil moisture 5% above wilting point for the crop fraction of grid box.
! This change violates LSM moisture budget, but
! can be considered as a compensation for irrigation not included into LSM. 
! "Irigation" could be applied when landuse fractional information
! is available and mosaic_lu=1.
    if(mosaic_lu == 1) then
      ! greenness factor: between 0 for min greenness and 1 for max greenness.
      factor = max(zero,min(one,(vegfra(i,j)-shdmin(i,j))/max(one,(shdmax(i,j)-shdmin(i,j)))))
           if (debug_print ) then
             if (abs(xlat-testptlat).lt.0.1 .and. &
                 abs(xlon-testptlon).lt.0.1)then
                 print *,'  lat,lon=',xlat,xlon,' factor=',factor
             endif
           endif

      if((ivgtyp(i,j) == natural .or. ivgtyp(i,j) == crop) .and. factor > 0.75) then
      ! cropland or grassland, apply irrigation during the growing seaspon when fraction 
      ! of greenness is > 0.75.

        do k=1,nroot
          cropsm=1.05_kind_phys*wilt - qmin
          cropfr = min(one,lufrac(crop) + 0.4*lufrac(natural)) ! assume that 40% of natural is cropland
          newsm = cropsm*cropfr + (1.-cropfr)*soilm1d(k)
          if(soilm1d(k) < newsm) then
           IF (debug_print ) THEN
             if (abs(xlat-testptlat).lt.0.1 .and. &
                 abs(xlon-testptlon).lt.0.1)then
                 print * ,'Soil moisture is below wilting in cropland areas at time step',ktau 
                 print * ,'  lat,lon=',xlat,xlon
                 print * ,'  lufrac=',lufrac,'factor=',factor &
                         ,'lai,ivgtyp,lufrac(crop),k,soilm1d(k),cropfr,wilt,cropsm,newsm,', &
                           lai(i,j),ivgtyp(i,j),lufrac(crop),k,soilm1d(k),cropfr,wilt,cropsm,newsm
             endif
           ENDIF
            soilm1d(k) = newsm
           IF (debug_print ) THEN
             if (abs(xlat-testptlat).lt.0.1 .and. &                  
                 abs(xlon-testptlon).lt.0.1)then
               print*,'  lat,lon=',xlat,xlon
               print * ,'Added soil water to cropland areas, k,soilm1d(k)',k,soilm1d(k)
             endif
           ENDIF
          endif ! < cropsm
        enddo
      endif ! crop
    endif ! mosaic_lu

!***  DIAGNOSTICS
!--- available and maximum soil moisture content in the soil
!--- domain

        smavail(i,j) = zero
        smmax (i,j)  = zero

      !do k=1,nzs-1
      !-- root-zone soil moisture
      do k=1,nroot
        smavail(i,j)=smavail(i,j)+(qmin+soilm1d(k))*             &
                    (zshalf(k+1)-zshalf(k))
        smmax (i,j) =smmax (i,j)+(qmin+dqm)*                     &
                    (zshalf(k+1)-zshalf(k))
      enddo

     if (debug_print) then
       if (abs(xlat-testptlat).lt.0.2 .and. abs(xlon-testptlon).lt.0.2)then
         print 100,'(RUC runoff)  i=',i,'  lat,lon=',xlat,xlon,  &
         'RUNOFF1', RUNOFF1(I,J), 'RUNOFF2 ',RUNOFF2(I,J),   &
         'edir ',edir(I,J),'ec ',ec(I,J),'ett ',ett(I,J)
       endif
     endif
!--- Convert the water unit into mm
        !-- three lines below are commented because accumulation
        !   happens in sfc_drv_ruc
        ACRUNOFF(I,J)  = (RUNOFF1(I,J)+RUNOFF2(I,J))*DT*rhowater
        SMAVAIL  (I,J) = SMAVAIL(I,J) * rhowater ! mm
        SMMAX    (I,J) = SMMAX(I,J) * rhowater
        smtotold (I,J) = smtotold(I,J) * rhowater ! mm

        do k=1,nzs

             soilmois(i,k,j) = soilm1d(k) + qmin
             sh2o    (i,k,j) = min(soiliqw(k) + qmin,soilmois(i,k,j))
                  tso(i,k,j) = tso1d(k)
        enddo

        tso(i,nzs,j) = tbot(i,j)

        do k=1,nzs
             smfr3d(i,k,j) = smfrkeep(k)
           keepfr3dflag(i,k,j) = keepfr (k)
        enddo

!tgs add together dew and cloud at the ground surface
!30july13        qcg(i,j)=qcg(i,j)+dew(i,j)/qkms

        Z0       (I,J) = ZNT (I,J)
        SFCEXC   (I,J) = TKMS
        patmb=P8w(i,1,j)*1.e-2_kind_phys
        Q2SAT=QSN(TABS,TBQ)/PATMB
        QSFC(I,J) = QVG(I,J)/(one+QVG(I,J))
        ! for MYJ surface and PBL scheme
        !      if (myj) then
        ! MYJSFC expects QSFC as actual specific humidity at the surface
        IF((QVATM.GE.Q2SAT*0.95_kind_phys).AND.QVATM.LT.qvg(I,J))THEN
          CHKLOWQ(I,J)=zero
        ELSE
          CHKLOWQ(I,J)=one
        ENDIF

        if(snow(i,j)==zero) EMISSL(i,j) = EMISBCK(i,j)
        EMISS (I,J) = EMISSL(I,J)
        ! SNOW is in [mm], SNWE is in [m]; CANWAT is in mm, CANWATR is in m
        SNOW   (i,j) = SNWE*1000._kind_phys
        SNOWH  (I,J) = SNHEI 
        CANWAT (I,J) = CANWATR*1000._kind_phys

     if (debug_print) then
       if (abs(xlat-testptlat).lt.0.2 .and. abs(xlon-testptlon).lt.0.2)then
         print *,'snow(i,j),soilt(i,j),xice(i,j),tso(i,:,j)',snow(i,j),soilt(i,j),xice(i,j),tso(i,:,j)
       endif
     endif
        INFILTR(I,J) = INFILTRP

        MAVAIL (i,j) = LMAVAIL(I,J)  
    IF (debug_print ) THEN
      if (abs(xlat-testptlat).lt.0.2 .and. abs(xlon-testptlon).lt.0.2)then
        print *,' LAND, I=,J=, QFX, HFX after SFCTMP', i,j,lh(i,j),hfx(i,j)
      endif
    ENDIF
        SFCEVP (I,J) = SFCEVP (I,J) + QFX (I,J) * DT
        GRDFLX (I,J) = -one * sflx(I,J)

!tgs - SMF.NE.0. when there is phase change in the top soil layer
! The heat of soil water freezing/thawing is not computed explicitly
! and is responsible for the residual in the energy budget.

!--- SNOWC snow cover flag
       SNOWC(I,J)=SNOWFRAC

!--- RHOSNF - density of snowfall
       RHOSNF(I,J)=RHOSNFALL

! Accumulated moisture flux [kg/m^2]
       SFCEVP (I,J) = SFCEVP (I,J) + QFX (I,J) * DT

!--tgs - SMF.NE.0. when there is phase change in the top soil layer
!        The heat of freezing/thawing of soil water is not computed explicitly
!        and is responsible for the residual in the energy budget.
!        endif
!        budget(i,j)=budget(i,j)-smf(i,j)

    if (debug_print ) then
     if (abs(xlat-testptlat).lt.0.2 .and.   &
         abs(xlon-testptlon).lt.0.2)then
     !-- compute budget for a test point
       ac=zero
       as=zero
       wb=zero

       ac=canwat(i,j)-canwatold(i,j)*rhowater ! canopy water change
       as=snwe-snowold(i,j) ! SWE change
       wb = smavail(i,j)-smtotold(i,j)
       waterbudget(i,j)=rainbl(i,j)+smelt(i,j)*dt*rhowater & ! source
                      -qfx(i,j)*dt &
                      -runoff1(i,j)*dt*rhowater-runoff2(i,j)*dt*rhowater &
                      -ac-as ! - (smavail(i,j)-smtotold(i,j))

       print *,'soilm1d ',i,soilm1d
       print 100,'(RUC budgets)  i=',i,'  lat,lon=',xlat,xlon,                     &
            'budget ',budget(i,j),'waterbudget',waterbudget(i,j),                  &
            'rainbl ',rainbl(i,j),'runoff1 ',runoff1(i,j),                         &
            'smelt ',smelt(i,j)*dt*1.e3_kind_phys,'smc change ',wb,                &
            'snwe change ',as,'canw change ',ac,'runoff2 ',runoff2(i,j),           &
            'qfx*dt ',qfx(i,j)*dt,'smavail ',smavail(i,j),'smcold',smtotold(i,j)
       !--
       print *,'Smf=',smf(i,j),i,j
       print *,'SNOW,SNOWold',i,j,snwe,snowold(i,j)
       print *,'SNOW-SNOWold',i,j,max(zero,snwe-snowold(i,j))
       print *,'CANWATold, canwat ',i,j,canwatold(i,j),canwat(i,j)
       print *,'canwat(i,j)-canwatold(i,j)',max(zero,canwat(i,j)-canwatold(i,j))
     endif
    endif

 100      format (";;; ",a,i4,a,2f14.7/(4(a10,'='es14.7)))


    IF (debug_print ) THEN
     if (abs(xlat-testptlat).lt.0.2 .and.   &
         abs(xlon-testptlon).lt.0.2)then
       print *,'LAND, i,tso1d,soilm1d,soilt - end of time step',         &
                  i,tso1d,soilm1d,soilt(i,j)
       print *,'LAND, QFX, HFX after SFCTMP', i,lh(i,j),hfx(i,j)
     endif
    ENDIF

!--- end of a land or sea ice point
        ENDIF
2999  continue ! lakes
      ENDDO

   ENDDO

!-----------------------------------------------------------------
   END SUBROUTINE LSMRUC
!! @}
!-----------------------------------------------------------------

!>\ingroup lsm_ruc_group
!! This subroutine solves energy and moisture budgets.
!! - It computes density of frozen precipitation from empirical 
!! dependence on temperature at the first atmospheric level.
!! - Computes amount of liquid and frozen precipitation intercepted by 
!! the vegetation canopy.
!! - In there is snow on the ground, the snow fraction is below 0.75,
!! the snow "mosaic" approach is turned on.
!! - Updates emissivity and albedo for patch snow.
   SUBROUTINE SFCTMP (debug_print, delt,ktau,conflx,i,j,         & !--- input variables
                xlat,xlon,testptlat,testptlon,                   &
                nzs,nddzs,nroot,meltfactor,                      &
                isncond_opt,isncovr_opt,                         &
                ILAND,ISOIL,IVGTYP,ISLTYP,PRCPMS,                &
                NEWSNMS,SNWE,SNHEI,SNOWFRAC,                     &
                exticeden,RHOSN,RHONEWSN_ex,RHONEWSN,RHOSNFALL,  &
                snowrat,grauprat,icerat,curat,                   &
                PATM,TABS,QVATM,QCATM,rho,                       &
                GLW,GSWdn,GSW,EMISS,EMISBCK,msnf,facsnf,         &
                QKMS,TKMS,PC,MAVAIL,CST,VEGFRA,ALB,ZNT,          &
                ALB_SNOW,ALB_SNOW_FREE,lai,hgt,stdev,            &
                MYJ,SEAICE,ISICE,                                &
                add_fire_heat_flux,fire_heat_flux,               &
                QWRTZ,rhocs,dqm,qmin,ref,wilt,psis,bclh,ksat,    & !--- soil fixed fields
                sat,cn,zsmain,zshalf,DTDZS,DTDZS2,tbq,           &
                cp,rovcp,g0,lv,stbolt,cw,c1sn,c2sn,              & !--- constants
                KQWRTZ,KICE,KWT,                                 &
                snweprint,snheiprint,rsm,                        & !---output variables
                soilm1d,ts1d,smfrkeep,keepfr,soilt,soilt1,       &
                tsnav,dew,qvg,qsg,qcg,                           &
                SMELT,SNOH,SNFLX,SNOM,SNOWFALLAC,ACSNOW,         &
                edir1,ec1,ett1,eeta,qfx,hfx,s,sublim,            &
                evapl,prcpl,fltot,runoff1,runoff2,soilice,       &
                soiliqw,infiltr,smf)
!-----------------------------------------------------------------
       IMPLICIT NONE
!-----------------------------------------------------------------

!--- input variables

   INTEGER,  INTENT(IN   )   ::  isice,i,j,nroot,ktau,nzs ,      &
                                 nddzs                             !nddzs=2*(nzs-2)
   integer,  intent(in   )   ::  isncond_opt,isncovr_opt

   real (kind_phys),     INTENT(IN   )   ::  DELT,CONFLX,meltfactor,xlat,xlon
   real (kind_phys),     INTENT(IN   )   ::  testptlat,testptlon
   real (kind_phys),     INTENT(IN   )   ::  C1SN,C2SN,RHONEWSN_ex
   LOGICAL,  INTENT(IN   )   ::  myj, debug_print, exticeden
!--- 3-D Atmospheric variables
   real (kind_phys)                                            , &
            INTENT(IN   )    ::                            PATM, &
                                                           TABS, &
                                                          QVATM, &
                                                          QCATM
   real (kind_phys)                                            , &
            INTENT(IN   )    ::                             GLW, &
                                                            GSW, &
                                                          GSWdn, &
                                                             PC, &
                                                    msnf,facsnf, &
                                                         VEGFRA, &
                                                  ALB_SNOW_FREE, &
                                                            lai, &
                                                      hgt,stdev, &
                                                         SEAICE, &
                                                            RHO, &
                                                           QKMS, &
                                                           TKMS, &
                                                 fire_heat_flux      
   LOGICAL,   INTENT(IN   )  ::              add_fire_heat_flux      
                                                             
   INTEGER,   INTENT(IN   )  ::                          IVGTYP, ISLTYP
!--- 2-D variables
   real (kind_phys)                                            , &
            INTENT(INOUT)    ::                           EMISS, &
                                                        EMISBCK, &
                                                         MAVAIL, &
                                                       SNOWFRAC, &
                                                       ALB_SNOW, &
                                                            ALB, &
                                                            CST

!--- soil properties
   real (kind_phys)                      ::                      &
                                                          RHOCS, &
                                                           BCLH, &
                                                            DQM, &
                                                           KSAT, &
                                                           PSIS, &
                                                           QMIN, &
                                                          QWRTZ, &
                                                            REF, &
                                                            SAT, &
                                                           WILT

   real (kind_phys),     INTENT(IN   )   ::                  CN, &
                                                             CW, &
                                                             CP, &
                                                          ROVCP, &
                                                             G0, &
                                                             LV, &
                                                         STBOLT, &
                                                         KQWRTZ, &
                                                           KICE, &
                                                            KWT

   real (kind_phys),     DIMENSION(1:NZS), INTENT(IN) :: ZSMAIN, &
                                                         ZSHALF, &
                                                         DTDZS2 


   real (kind_phys),     DIMENSION(1:NDDZS), INTENT(IN) :: DTDZS

   real (kind_phys),     DIMENSION(1:5001), INTENT(IN) :: TBQ


!--- input/output variables
!-------- 3-d soil moisture and temperature
   real (kind_phys),     DIMENSION( 1:nzs )                    , &
             INTENT(INOUT)   ::                            TS1D, & 
                                                        SOILM1D, &
                                                       SMFRKEEP
   real (kind_phys),  DIMENSION( 1:nzs )                       , &
             INTENT(INOUT)   ::                          KEEPFR

   real (kind_phys),  DIMENSION(1:NZS),INTENT(INOUT) :: SOILICE, &
                                                        SOILIQW
          

   INTEGER, INTENT(INOUT)    ::                     ILAND,ISOIL
   INTEGER                   ::                     ILANDs

!-------- 2-d variables
   real (kind_phys)                                            , &
             INTENT(INOUT)   ::                             DEW, &
                                                          EDIR1, &
                                                            EC1, &
                                                           ETT1, &
                                                           EETA, &
                                                          EVAPL, &
                                                        INFILTR, &
                                                          RHOSN, &
                                                       RHONEWSN, &
                                                      rhosnfall, &
                                                        snowrat, &
                                                       grauprat, &
                                                         icerat, &
                                                          curat, &
                                                         SUBLIM, &
                                                          PRCPL, &
                                                            QVG, &
                                                            QSG, &
                                                            QCG, &
                                                            QFX, &
                                                            HFX, &
                                                          fltot, &
                                                            smf, &
                                                              S, &  
                                                        RUNOFF1, &
                                                        RUNOFF2, &
                                                         ACSNOW, &
                                                     SNOWFALLAC, &
                                                           SNWE, &
                                                          SNHEI, &
                                                          SMELT, &
                                                           SNOM, &
                                                           SNOH, &
                                                          SNFLX, &
                                                          SOILT, &
                                                         SOILT1, &
                                                          TSNAV, &
                                                            ZNT

   real (kind_phys),     DIMENSION(1:NZS)              ::        &
                                                           tice, &
                                                        rhosice, &
                                                         capice, &
                                                       thdifice, &
                                                          TS1DS, &
                                                       SOILM1DS, &
                                                      SMFRKEEPS, &
                                                       SOILIQWS, & 
                                                       SOILICES, &
                                                        KEEPFRS
!-------- 1-d variables
   real (kind_phys) :: &
                                                            DEWS, &
                                                        MAVAILS,  &
                                                          EDIR1s, &
                                                            EC1s, &
                                                            csts, &
                                                           ETT1s, &
                                                           EETAs, &
                                                          EVAPLs, &
                                                        INFILTRs, &
                                                          PRCPLS, &
                                                            QVGS, &
                                                            QSGS, &
                                                            QCGS, &
                                                            QFXS, &
                                                            HFXS, &
                                                          fltots, &
                                                        RUNOFF1S, &
                                                        RUNOFF2s, &
                                                              SS, &
                                                          SOILTs

            
                     

   real (kind_phys),  INTENT(INOUT)                     ::  RSM, &  
                                                      SNWEPRINT, &
                                                     SNHEIPRINT
!--- Local variables

   INTEGER ::  K,ILNB

   real (kind_phys)    ::  BSN, XSN                            , &
               RAINF, SNTH, NEWSN, PRCPMS, NEWSNMS             , &
               T3, UPFLUX, XINET, snowfrac2, m
   real (kind_phys)    ::  snhei_crit, snhei_crit_newsn, keep_snow_albedo, SNOWFRACnewsn
   real (kind_phys)    ::  newsnowratio, dd1

   real (kind_phys)    ::  rhonewgr,rhonewice

   real (kind_phys)    ::  RNET,GSWNEW,GSWIN,EMISSN,ZNTSN,EMISS_snowfree
   real (kind_phys)    ::  VEGFRAC, snow_mosaic, snfr, vgfr
   real    ::  cice, albice, albsn, drip, dripsn, dripliq
   real    ::  interw, intersn, infwater, intwratio

!-----------------------------------------------------------------
        integer,   parameter      ::      ilsnow=99 
        
    IF (debug_print ) THEN
        print *,' in SFCTMP',i,j,nzs,nddzs,nroot,                 &
                 SNWE,RHOSN,SNOM,SMELT,TS1D
    ENDIF

     !-- Snow fraction options
     !-- option 1: original formulation using critical snow depth to compute
     !-- snow fraction
     !-- option 2: the tanh formulation from Niu,G.-Y.,and Yang,Z.-L. 2007,JGR,DOI:10.1029/2007JD008674.
     !-- option 3: the tanh formulation from Niu,G.-Y.,and Yang,Z.-L. 2007,JGR,DOI:10.1029/2007JD008674.
     !             with vegetation dependent parameters from Noah MP (personal
     !             communication with Mike Barlage)
     !-- SNHEI_CRIT is a threshold for fractional snow in isncovr_opt=1
         snhei_crit=0.01601_kind_phys*rhowater/rhosn
         snhei_crit_newsn=0.0005_kind_phys*rhowater/rhosn
     !--
        zntsn = z0tbl(isice)
        snow_mosaic = zero
        snfr = one
        NEWSN= zero
        newsnowratio = zero
        snowfracnewsn= zero
        snowfrac2= zero
        rhonewsn = 100._kind_phys
        if(snhei == zero) snowfrac=zero
        smelt = zero
        RAINF = zero
        RSM = zero
        DD1 = zero
        INFILTR = zero
! Jul 2016 -  Avissar and Pielke (1989)
! This formulation depending on LAI defines relative contribution of the vegetation to
! the total heat fluxes between surface and atmosphere.
! With VEGFRA=100% and LAI=3, VEGFRAC=0.86 meaning that vegetation contributes
! only 86% of the total surface fluxes.
!        VGFR=0.01*VEGFRA ! % --> fraction
!        VEGFRAC=2.*lai*vgfr/(1.+2.*lai*vgfr)
        VEGFRAC=0.01_kind_phys*VEGFRA
        drip = zero
        dripsn = zero
        dripliq = zero
        smf = zero
        interw = zero
        intersn = zero
        infwater = zero

!---initialize local arrays for sea ice
          do k=1,nzs
            tice(k) = zero
            rhosice(k) = zero
            cice = zero
            capice(k) = zero
            thdifice(k) = zero
          enddo

        GSWnew=GSW
        GSWin=GSWdn !/(1.-alb)
        ALBice=ALB_SNOW_FREE
        ALBsn=alb_snow
        EMISSN = 0.99_kind_phys ! from setemis, from WRF - 0.98
        EMISS_snowfree = EMISBCK ! LEMITBL(IVGTYP)

!--- sea ice properties
!--- N.N Zubov "Arctic Ice"
!--- no salinity dependence because we consider the ice pack
!--- to be old and to have low salinity (0.0002)
       if(SEAICE.ge.0.5_kind_phys) then
          do k=1,nzs
            tice(k) = ts1d(k) - tfrz
            rhosice(k) = 917.6_kind_phys/(one-0.000165_kind_phys*tice(k))
            cice = 2115.85_kind_phys +7.7948_kind_phys*tice(k)
            capice(k) = cice*rhosice(k)
            thdifice(k) = 2.260872_kind_phys/capice(k)
           enddo
!-- SEA ICE ALB dependence on ice temperature. When ice temperature is
!-- below critical value of -10C - no change to albedo.
!-- If temperature is higher that -10C then albedo is decreasing.
!-- The minimum albedo at t=0C for ice is 0.1 less.
       ALBice = MIN(ALB_SNOW_FREE,MAX(ALB_SNOW_FREE - 0.05_kind_phys,   &
               ALB_SNOW_FREE - 0.1_kind_phys*(tice(1)+10._kind_phys)/10._kind_phys ))
       endif

    IF (debug_print ) THEN
        print *,'alb_snow_free',ALB_SNOW_FREE
        print *,'GSW,GSWnew,GLW,SOILT,EMISS,ALB,ALBice,SNWE',&
                 GSW,GSWnew,GLW,SOILT,EMISS,ALB,ALBice,SNWE
    ENDIF

	if(snhei.gt.0.0081_kind_phys*rhowater/rhosn) then
!*** Update snow density for current temperature (Koren et al 1999,doi:10.1029/1999JD900232.)
        BSN=delt/3600._kind_phys*c1sn*exp(0.08_kind_phys*min(zero,tsnav)-c2sn*rhosn*1.e-3_kind_phys)
       if(bsn*snwe*100._kind_phys.lt.1.e-4_kind_phys) goto 777
        XSN=rhosn*(exp(bsn*snwe*100._kind_phys)-one)/(bsn*snwe*100._kind_phys)
        rhosn=MIN(MAX(58.8_kind_phys,XSN),500._kind_phys)
 777   continue
      endif

      !-- snow_mosaic from the previous time step 
      if(snowfrac < 0.75_kind_phys) snow_mosaic = one

           newsn=newsnms*delt
!---- ACSNOW - run-total snowfall water [mm]
           acsnow=acsnow+newsn*rhowater

       IF(NEWSN.GT.zero) THEN

    IF (debug_print ) THEN
      print *, 'THERE IS NEW SNOW, newsn', newsn
    ENDIF

        newsnowratio = min(one,newsn/(snwe+newsn))

        !if(isncovr_opt == 2) then
        !-- update snow fraction for fresh snowfall (Swenson&Lawrence,JGR,2012)
        !   time-step snowfall [mm H2O], 0.1 - accumulation constant (unitless)
        !  snowfrac = snowfrac + tanh(0.1*newsn*1.e3)*(1.-snowfrac) ! eq. 8.1 from CLM5
        !  if(debug_print) print *,'2 - snowfrac newsn', i,j,ktau,snowfrac
        !endif

!--- 27 Feb 2014 - empirical formulations from John M. Brown
!        rhonewsn=min(250.,rhowater/max(4.179,(13.*tanh((274.15-Tabs)*0.3333))))
!--- 13 Mar 2018 - formulation from Trevor Elcott
      if (exticeden) then
        rhonewsn = rhonewsn_ex
      else
        rhonewsn=min(125._kind_phys,rhowater/max(8._kind_phys,(17._kind_phys*tanh((276.65_kind_phys-Tabs)*0.15_kind_phys))))
        rhonewgr=min(500._kind_phys,rhowater/max(2._kind_phys,(3.5_kind_phys*tanh((274.15_kind_phys-Tabs)*0.3333_kind_phys))))
        rhonewice=rhonewsn

!--- compute density of "snowfall" from weighted contribution
!                 of snow, graupel and ice fractions

         rhosnfall = min(500._kind_phys,max(58.8_kind_phys,(rhonewsn*snowrat +  &
                     rhonewgr*grauprat + rhonewice*icerat + rhonewgr*curat)))

    if (debug_print) then
    !if (abs(xlat-33.35).lt.0.2 .and. abs(xlon-272.55).lt.0.2)then
      print *,' xlat, xlon', xlat, xlon
      print *,'snow_mosaic = ',snow_mosaic
      print *,'new snow,newsnowratio,rhosnfall =',newsn,newsnowratio,rhosnfall
      print *,'snowrat,grauprat,icerat,curat,rhonewgr,rhonewsn,rhonewice',snowrat,grauprat,icerat,curat,rhonewgr,rhonewsn,rhonewice
    endif
! from now on rhonewsn is the density of falling frozen precipitation
        rhonewsn=rhosnfall
      end if

!*** Define average snow density of the snow pack considering
!*** the amount of fresh snow (eq. 9 in Koren et al.(1999) 
!*** without snow melt )
         xsn=(rhosn*snwe+rhonewsn*newsn)/                         &
             (snwe+newsn)
         rhosn=MIN(MAX(58.8_kind_phys,XSN),500._kind_phys)
       ENDIF ! end NEWSN > 0.

       IF(PRCPMS > zero) THEN

! PRCPMS is liquid precipitation rate
! RAINF is a flag used for calculation of rain water
! heat content contribution into heat budget equation. Rain's temperature
! is set equal to air temperature at the first atmospheric
! level.  

           RAINF=one
       ENDIF

        drip = zero
        intwratio= zero
     if(vegfrac > 0.01_kind_phys) then
! compute intercepted precipitation - Eq. 1 Lawrence et al.,
! J. of Hydrometeorology, 2006, CLM.
         interw=0.25_kind_phys*DELT*PRCPMS*(one-exp(-0.5_kind_phys*lai))*vegfrac
         intersn=0.25_kind_phys*NEWSN*(one-exp(-0.5_kind_phys*lai))*vegfrac
         infwater=PRCPMS - interw/delt
    if((interw+intersn) > zero) then
       intwratio=interw/(interw+intersn)
    endif

! Update water/snow intercepted by the canopy
         dd1=CST + interw + intersn
         CST=DD1
        IF(CST.GT.SAT) THEN
          CST=SAT
          DRIP=DD1-SAT
        ENDIF
     else
         CST=zero
         DRIP=zero
         interw=zero
         intersn=zero
         infwater=PRCPMS
     endif ! vegfrac > 0.01


       IF(NEWSN.GT.zero) THEN
!Update snow on the ground
         snwe=max(zero,snwe+newsn-intersn)
! Add drip to snow on the ground
      if(drip > zero) then
       if (snow_mosaic==one) then
         dripliq=drip*intwratio
         dripsn = drip - dripliq
         snwe=snwe+dripsn
         infwater=infwater+dripliq
         dripliq=zero
         dripsn = zero
       else
         snwe=snwe+drip
       endif
      endif
         snhei=snwe*rhowater/rhosn
         NEWSN=NEWSN*rhowater/rhonewsn
       ENDIF

   IF(SNHEI.GT.zero) THEN
!-- SNOW on the ground
!--- Land-use category should be changed to snow/ice for grid points with snow>0
        ILAND=ISICE
!24nov15 - based on field exp on Pleasant View soccer fields
!    if(meltfactor > 1.5) then ! all veg. types, except forests
!         SNHEI_CRIT=0.01601*1.e3/rhosn
! Petzold - 1 cm of fresh snow overwrites effects from old snow.
! Need to test SNHEI_CRIT_newsn=0.01
!         SNHEI_CRIT_newsn=0.01
!    else  ! forests
!         SNHEI_CRIT=0.02*1.e3/rhosn
!         SNHEI_CRIT_newsn=0.001*1.e3/rhosn
!    endif

      !-- update snow cover with accounting for fresh snow
        m = one ! m=1.6 in Niu&Yang, m=1 in CLM
      if(isncovr_opt == 1) then
        snowfrac=min(one,snhei/(2._kind_phys*snhei_crit))
      elseif(isncovr_opt == 2) then
        snowfrac=min(one,snhei/(2._kind_phys*snhei_crit))
        if(ivgtyp == glacier .or. ivgtyp == bare) then
        !-- sparsely vegetated or land ice
          snowfrac2 = tanh( snhei/(2.5_kind_phys * 0.2_kind_phys *(rhosn/rhonewsn)**m))
        else
          !-- Niu&Yang: znt=0.01 m for 1 degree (100km) resolution tests
          !  on 3-km scale use actual roughness, but not higher than 0.2 m.
          !  The factor is 20 for forests (~100/dx = 33.)
          snowfrac2 = tanh( snhei/(2.5_kind_phys *min(0.2_kind_phys,znt) *(rhosn/rhonewsn)**m))
        endif
        !-- snow fraction is average between method 1 and 2
        snowfrac = 0.5_kind_phys*(snowfrac+snowfrac2)
      else
      !-- isncovr_opt=3
        !m = msnf ! vegetation dependent facsnf/msnf from Noah MP
        !-- for RRFS a factor 10. was added to 'facsnf' to get reasonal values of
        !   snow cover fractions on the 3-km scale.
        !   This factor is scale dependent.
        snowfrac = tanh( snhei/(10._kind_phys * facsnf *(rhosn/rhonewsn)**m))
      endif

       if(newsn > zero ) then
         SNOWFRACnewsn=MIN(one,snowfallac*1.e-3_kind_phys/SNHEI_CRIT_newsn)
       endif

       !-- due to steep slopes and blown snow, limit snow fraction in the
       !-- mountains to 0.85 (based on Swiss weather model over the Alps)
       if(hgt > 2500._kind_phys .and. ivgtyp == glacier) snowfrac=min(0.85_kind_phys,snowfrac)

       !24nov15 - SNOWFRAC for urban category < 0.75 
       if(ivgtyp == urban) snowfrac=min(0.75_kind_phys,snowfrac)

       if(snowfrac < 0.75_kind_phys) snow_mosaic = one

       KEEP_SNOW_ALBEDO = zero
       IF (snowfracnewsn > 0.99_kind_phys .and. rhosnfall < 450._kind_phys) THEN
       ! new snow
             KEEP_SNOW_ALBEDO = one
             ! turn off separate treatment of snow covered and snow-free portions of the grid cell
             snow_mosaic=0.  ! ???
      ENDIF

    IF (debug_print ) THEN
      print *,'SNHEI_CRIT,SNOWFRAC,SNHEI_CRIT_newsn,SNOWFRACnewsn', &
               SNHEI_CRIT,SNOWFRAC,SNHEI_CRIT_newsn,SNOWFRACnewsn
    ENDIF

!-- Set znt for snow from VEGPARM table (snow/ice landuse), except for
!-- land-use types with higher roughness (forests, urban).
      IF(newsn.eq.zero .and. znt.le.0.2_kind_phys .and. IVGTYP.ne.isice) then
         if( snhei .le. 2._kind_phys*ZNT)then
         ! shallow snow
           znt=0.55_kind_phys*znt+0.45_kind_phys*z0tbl(iland)
         elseif( snhei .gt. 2._kind_phys*ZNT .and. snhei .le. 4._kind_phys*ZNT)then
           znt=0.2_kind_phys*znt+0.8_kind_phys*z0tbl(iland)
         elseif(snhei > 4._kind_phys*ZNT) then
         ! deep snow
           znt=z0tbl(iland)
         endif
       ENDIF


    IF(SEAICE .LT. 0.5_kind_phys) THEN
!----- SNOW on soil
!-- ALB dependence on snow depth
! ALB_SNOW across Canada's forested areas is very low - 0.27-0.35, this
! causes significant warm biases. Limiting ALB in these areas to be higher than 0.4
! hwlps with these biases.. 
     if( snow_mosaic == one) then
         ALBsn=alb_snow
         if(KEEP_SNOW_ALBEDO > 0.9_kind_phys .and. albsn < 0.4_kind_phys) then
         !-- Albedo correction with fresh snow and deep snow pack
         !-- will reduce warm bias in western Canada
         !-- and US West coast, where max snow albedo is low (0.3-0.5).
           !print *,'ALB increase to 0.7',alb_snow,snhei,snhei_crit,albsn,i,j
           !ALBsn = 0.7_kind_phys
         endif

         Emiss= emissn
     else
         ALBsn   = MAX(keep_snow_albedo*alb_snow,               &
                   MIN((alb_snow_free +                         &
           (alb_snow - alb_snow_free) * snowfrac), alb_snow))
         if(newsn > zero .and. KEEP_SNOW_ALBEDO > 0.9_kind_phys .and. albsn < 0.4_kind_phys) then
         !-- Albedo correction with fresh snow and deep snow pack
         !-- will reduce warm bias in western Canada
         !-- and US West coast, where max snow albedo is low (0.3-0.5).
           !print *,'ALB increase to 0.7',alb_snow,snhei,snhei_crit,albsn,i,j
           !ALBsn = 0.7_kind_phys
           !print *,'NO mosaic ALB increase to 0.7',alb_snow,snhei,snhei_crit,alb,i,j
         endif

         Emiss   = MAX(keep_snow_albedo*emissn,                 &
                   MIN((emiss_snowfree +                         &
           (emissn - emiss_snowfree) * snowfrac), emissn))
     endif ! snow_mosaic

    IF (debug_print ) THEN
    !if (abs(xlat-33.35).lt.0.2 .and. abs(xlon-272.55).lt.0.2)then
      print *,'Snow on soil ALBsn,emiss,snow_mosaic',i,j,ALBsn,emiss,snow_mosaic
    ENDIF
!28mar11  if canopy is covered with snow to 95% of its capacity and snow depth is
! higher than patchy snow treshold - then snow albedo is not less than 0.55
! (inspired by the flight from Fairbanks to Seatle)

!-- ALB dependence on snow temperature. When snow temperature is
!-- below critical value of -10C - no change to albedo.
!-- If temperature is higher that -10C then albedo is decreasing.
!-- The minimum albedo at t=0C for snow on land is 15% less than
!-- albedo of temperatures below -10C.
     if(albsn.lt.0.4_kind_phys .or. keep_snow_albedo==1) then
        ALB=ALBsn
      else
!-- change albedo when no fresh snow and snow albedo is higher than 0.5
        ALB = MIN(ALBSN,MAX(ALBSN - 0.1_kind_phys*(soilt - 263.15_kind_phys)/       &
                (tfrz-263.15_kind_phys)*ALBSN, ALBSN - 0.05_kind_phys))
      endif
    ELSE
!----- SNOW on ice
     if( snow_mosaic == one) then
         ALBsn=alb_snow
         Emiss= emissn
     else
         ALBsn   = MAX(keep_snow_albedo*alb_snow,               &
                   MIN((albice + (alb_snow - albice) * snowfrac), alb_snow))
         Emiss   = MAX(keep_snow_albedo*emissn,                 &
                   !-- emiss_snowfree=0.96 in setemis
                   MIN((emiss_snowfree +                        &
                   (emissn - emiss_snowfree) * snowfrac), emissn))
     endif

    IF (debug_print ) THEN
  print *,'Snow on ice snow_mosaic,ALBsn,emiss',i,j,ALBsn,emiss,snow_mosaic
    ENDIF
!-- ALB dependence on snow temperature. When snow temperature is
!-- below critical value of -10C - no change to albedo.
!-- If temperature is higher that -10C then albedo is decreasing.
      if(albsn.lt.alb_snow .or. keep_snow_albedo .eq.one)then
       ALB=ALBsn
      else
!-- change albedo when no fresh snow
       ALB = MIN(ALBSN,MAX(ALBSN - 0.15_kind_phys*ALBSN*(soilt - 263.15_kind_phys)/  &
                (tfrz-263.15_kind_phys), ALBSN - 0.1_kind_phys))
      endif

    ENDIF

    if (snow_mosaic==one) then 
!may 2014 - treat separately snow-free and snow-covered areas

       if(SEAICE .LT. 0.5_kind_phys) then
!  LAND
! portion not covered with snow
! compute absorbed GSW for snow-free portion

         gswnew=GSWin*(one-alb_snow_free)
!--------------
         T3      = STBOLT*SOILT*SOILT*SOILT
         UPFLUX  = T3 *SOILT
         XINET   = EMISS_snowfree*(GLW-UPFLUX)
         RNET    = GSWnew + XINET
         IF ( add_fire_heat_flux .and. fire_heat_flux >0 ) then ! JLS
          IF (debug_print ) THEN
            print *,'RNET snow-free, fire_heat_flux, xlat/xlon',RNET, fire_heat_flux,xlat,xlon
          ENDIF
            RNET = RNET + fire_heat_flux
         ENDIF

    IF (debug_print ) THEN
     print *,'Fractional snow - snowfrac=',snowfrac
     print *,'Snowfrac<1 GSWin,GSWnew -',GSWin,GSWnew,'SOILT, RNET',soilt,rnet
    ENDIF
           do k=1,nzs
          soilm1ds(k) = soilm1d(k)
          ts1ds(k) = ts1d(k)
          smfrkeeps(k) = smfrkeep(k)
          keepfrs(k) = keepfr(k)
          soilices(k) = soilice(k)
          soiliqws(k) = soiliqw(k)
            enddo
          soilts = soilt
          qvgs = qvg
          qsgs = qsg
          qcgs = qcg
          csts = cst
          mavails = mavail
          smelt=zero
          runoff1s=zero
          runoff2s=zero
       
          ilands = ivgtyp

         CALL SOIL(debug_print,xlat, xlon, testptlat, testptlon,&
!--- input variables
            i,j,iland,isoil,delt,ktau,conflx,nzs,nddzs,nroot,   &
            PRCPMS,RAINF,PATM,QVATM,QCATM,GLW,GSWnew,gswin,     &
            EMISS_snowfree,RNET,QKMS,TKMS,PC,csts,dripliq,      &
            infwater,rho,vegfrac,lai,myj,                       &
!--- soil fixed fields 
            QWRTZ,rhocs,dqm,qmin,ref,wilt,                      &
            psis,bclh,ksat,sat,cn,                              &
            zsmain,zshalf,DTDZS,DTDZS2,tbq,                     &
!--- constants
            lv,CP,rovcp,G0,cw,stbolt,tabs,                      &
            KQWRTZ,KICE,KWT,                                    &
!--- output variables for snow-free portion
            soilm1ds,ts1ds,smfrkeeps,keepfrs,                   &
            dews,soilts,qvgs,qsgs,qcgs,edir1s,ec1s,             &
            ett1s,eetas,qfxs,hfxs,ss,evapls,prcpls,fltots,runoff1s, &
            runoff2s,mavails,soilices,soiliqws,                 &
            infiltrs,smf)
        else
! SEA ICE
! portion not covered with snow
! compute absorbed GSW for snow-free portion

         gswnew=GSWin*(one-albice)
!--------------
         T3      = STBOLT*SOILT*SOILT*SOILT
         UPFLUX  = T3 *SOILT
         XINET   = EMISS_snowfree*(GLW-UPFLUX)
         RNET    = GSWnew + XINET
    IF (debug_print ) THEN
     print *,'Fractional snow - snowfrac=',snowfrac
     print *,'Snowfrac<1 GSWin,GSWnew -',GSWin,GSWnew,'SOILT, RNET',soilt,rnet
    ENDIF
            do k=1,nzs
          ts1ds(k) = ts1d(k)
            enddo
          soilts = soilt
          qvgs = qvg
          qsgs = qsg
          qcgs = qcg
          smelt=zero
          runoff1s=zero
          runoff2s=zero
 
          CALL SICE(debug_print,xlat,xlon,                      &
!--- input variables
            i,j,iland,isoil,delt,ktau,conflx,nzs,nddzs,nroot,   &
            PRCPMS,RAINF,PATM,QVATM,QCATM,GLW,GSWnew,           &
            0.98_kind_phys,RNET,QKMS,TKMS,rho,myj,              &
!--- sea ice parameters
            tice,rhosice,capice,thdifice,                       &
            zsmain,zshalf,DTDZS,DTDZS2,tbq,                     &
!--- constants
            lv,CP,rovcp,cw,stbolt,tabs,                         &
!--- output variable
            ts1ds,dews,soilts,qvgs,qsgs,qcgs,                   &
            eetas,qfxs,hfxs,ss,evapls,prcpls,fltots             &
                                                                )
           edir1 = eeta*1.e-3_kind_phys
           ec1 = zero
           ett1 = zero
           runoff1 = prcpms
           runoff2 = zero
           mavail = one
           infiltr= zero
           cst= zero
            do k=1,nzs
               soilm1d(k)=one
               soiliqw(k)=zero
               soilice(k)=one
               smfrkeep(k)=one
               keepfr(k)=zero
            enddo
        endif ! seaice < 0.5

    endif ! snow_mosaic=1.
                           
!--- recompute absorbed solar radiation and net radiation
!--- for updated value of snow albedo - ALB
         gswnew=GSWin*(one-alb)
!--------------
         T3      = STBOLT*SOILT*SOILT*SOILT
         UPFLUX  = T3 *SOILT
         XINET   = EMISS*(GLW-UPFLUX)
         RNET    = GSWnew + XINET
    IF (debug_print ) THEN
    !if (abs(xlat-testptlat).lt.0.1 .and. abs(xlon-testptlon).lt.0.1)then
        print *,'RNET=',rnet
        print *,'SNOW - I,J,newsn,snwe,snhei,GSW,GSWnew,GLW,UPFLUX,ALB',&
                 i,j,newsn,snwe,snhei,GSW,GSWnew,GLW,UPFLUX,ALB
        print *,'GSWnew',gswnew,'alb=',alb
    ENDIF

      if (SEAICE .LT. 0.5_kind_phys) then
! LAND
         IF ( add_fire_heat_flux .and. fire_heat_flux>0 ) then ! JLS
          IF (debug_print ) THEN
           print *,'RNET snow, fire_heat_flux, xlat/xlon',RNET, fire_heat_flux,xlat,xlon
          ENDIF
            RNET = RNET + fire_heat_flux
         ENDIF
           if(snow_mosaic==one)then
              snfr=one
           else
              snfr=snowfrac
           endif
         CALL SNOWSOIL (debug_print,xlat,xlon,testptlat,testptlon, & !--- input variables
            i,j,isoil,delt,ktau,conflx,nzs,nddzs,nroot,         &
            isncond_opt,isncovr_opt,                            &
            meltfactor,rhonewsn,SNHEI_CRIT,                     &  ! new
            ILAND,PRCPMS,RAINF,NEWSN,snhei,SNWE,snfr,           &
            RHOSN,PATM,QVATM,QCATM,                             &
            GLW,GSWnew,GSWin,EMISS,RNET,IVGTYP,                 &
            QKMS,TKMS,PC,CST,dripsn,infwater,                   &
            RHO,VEGFRAC,ALB,ZNT,lai,                            &
            MYJ,                                                &
!--- soil fixed fields
            QWRTZ,rhocs,dqm,qmin,ref,wilt,psis,bclh,ksat,       &
            sat,cn,zsmain,zshalf,DTDZS,DTDZS2,tbq,              & 
!--- constants
            lv,CP,rovcp,G0,cw,stbolt,tabs,                      &
            KQWRTZ,KICE,KWT,                                    &
!--- output variables
            ilnb,snweprint,snheiprint,rsm,                      &
            soilm1d,ts1d,smfrkeep,keepfr,                       &
            dew,soilt,soilt1,tsnav,qvg,qsg,qcg,                 &
            SMELT,SNOH,SNFLX,SNOM,edir1,ec1,ett1,eeta,          &
            qfx,hfx,s,sublim,prcpl,fltot,runoff1,runoff2,       &
            mavail,soilice,soiliqw,infiltr                      )
       else
! SEA ICE
           if(snow_mosaic==one)then
              snfr=one
           else
              snfr=snowfrac
           endif

         CALL SNOWSEAICE (debug_print,xlat,xlon,                &
            i,j,isoil,delt,ktau,conflx,nzs,nddzs,               &    
            isncond_opt,isncovr_opt,                            &
            meltfactor,rhonewsn,SNHEI_CRIT,                     &  ! new
            ILAND,PRCPMS,RAINF,NEWSN,snhei,SNWE,snfr,           &    
            RHOSN,PATM,QVATM,QCATM,                             &    
            GLW,GSWnew,EMISS,RNET,                              &    
            QKMS,TKMS,RHO,myj,                                  &    
!--- sea ice parameters
            ALB,ZNT,                                            &
            tice,rhosice,capice,thdifice,                       &    
            zsmain,zshalf,DTDZS,DTDZS2,tbq,                     &    
!--- constants
            lv,CP,rovcp,cw,stbolt,tabs,                         &    
!--- output variables
            ilnb,snweprint,snheiprint,rsm,ts1d,                 &    
            dew,soilt,soilt1,tsnav,qvg,qsg,qcg,                 &    
            SMELT,SNOH,SNFLX,SNOM,eeta,                         &    
            qfx,hfx,s,sublim,prcpl,fltot                        &    
                                                                )    
           edir1 = eeta*1.e-3_kind_phys
           ec1 = zero
           ett1 = zero
           runoff1 = smelt
           runoff2 = zero
           mavail = one
           infiltr = zero
           cst = zero
            do k=1,nzs
               soilm1d(k)=one
               soiliqw(k)=zero
               soilice(k)=one
               smfrkeep(k)=one
               keepfr(k)=zero
            enddo
       endif


     if (snow_mosaic==one) then
! May 2014 - now combine snow covered and snow-free land fluxes, soil temp, moist,
! etc.
        if(SEAICE .LT. 0.5_kind_phys) then
! LAND
   IF (debug_print ) THEN
    !if (abs(xlat-33.35).lt.0.2 .and. abs(xlon-272.55).lt.0.2)then
      print *,' xlat, xlon', xlat, xlon
      print *,' snowfrac = ',snowfrac
      print *,' SOILT snow on land', ktau, i,j,soilt
      print *,' SOILT on snow-free land', i,j,soilts
      print *,' ts1d,ts1ds',i,j,ts1d,ts1ds
      print *,' SNOW flux',i,j, snflx
      print *,' Ground flux on snow-covered land',i,j, s
      print *,' Ground flux on snow-free land', i,j,ss
      print *,' CSTS, CST', i,j,csts,cst
   ENDIF

            do k=1,nzs
          soilm1d(k) = soilm1ds(k)*(1.-snowfrac) + soilm1d(k)*snowfrac
          ts1d(k) = ts1ds(k)*(1.-snowfrac) + ts1d(k)*snowfrac
          smfrkeep(k) = smfrkeeps(k)*(1.-snowfrac) + smfrkeep(k)*snowfrac
       if(snowfrac > 0.5_kind_phys) then
          keepfr(k) = keepfr(k)
       else
          keepfr(k) = keepfrs(k)
       endif
          soilice(k) = soilices(k)*(1.-snowfrac) + soilice(k)*snowfrac
          soiliqw(k) = soiliqws(k)*(1.-snowfrac) + soiliqw(k)*snowfrac
            enddo
          dew = dews*(one-snowfrac) + dew*snowfrac
          soilt = soilts*(one-snowfrac) + soilt*snowfrac
          qvg = qvgs*(one-snowfrac) + qvg*snowfrac
          qsg = qsgs*(one-snowfrac) + qsg*snowfrac
          qcg = qcgs*(one-snowfrac) + qcg*snowfrac
          edir1 = edir1s*(one-snowfrac) + edir1*snowfrac
          ec1 = ec1s*(one-snowfrac) + ec1*snowfrac
          cst = csts*(one-snowfrac) + cst*snowfrac
          ett1 = ett1s*(one-snowfrac) + ett1*snowfrac
          eeta = eetas*(one-snowfrac) + eeta*snowfrac
          qfx = qfxs*(one-snowfrac) + qfx*snowfrac
          hfx = hfxs*(one-snowfrac) + hfx*snowfrac
          s = ss*(one-snowfrac) + s*snowfrac
          evapl = evapls*(one-snowfrac)
          prcpl = prcpls*(one-snowfrac) + prcpl*snowfrac
          fltot = fltots*(one-snowfrac) + fltot*snowfrac
          ALB   = MAX(keep_snow_albedo*alb,              &
                  MIN((alb_snow_free + (alb - alb_snow_free) * snowfrac), alb))

          Emiss = MAX(keep_snow_albedo*emissn,           &
                  MIN((emiss_snowfree +                  &
              (emissn - emiss_snowfree) * snowfrac), emissn))

          runoff1 = runoff1s*(one-snowfrac) + runoff1*snowfrac
          runoff2 = runoff2s*(one-snowfrac) + runoff2*snowfrac
          mavail = mavails*(one-snowfrac) + one*snowfrac
          infiltr = infiltrs*(one-snowfrac) + infiltr*snowfrac

    IF (debug_print ) THEN
    !if (abs(xlat-33.35).lt.0.2 .and. & abs(xlon-272.55).lt.0.2)then
      print *,' Ground flux combined', xlat,xlon, s
      print *,' SOILT combined on land', soilt
      print *,' TS combined on land', ts1d
    ENDIF
       else
! SEA ICE
! Now combine fluxes for snow-free sea ice and snow-covered area
    IF (debug_print ) THEN
      print *,'SOILT snow on ice', soilt
    ENDIF
            do k=1,nzs
          ts1d(k) = ts1ds(k)*(one-snowfrac) + ts1d(k)*snowfrac
            enddo
          dew = dews*(one-snowfrac) + dew*snowfrac
          soilt = soilts*(one-snowfrac) + soilt*snowfrac
          qvg = qvgs*(one-snowfrac) + qvg*snowfrac
          qsg = qsgs*(one-snowfrac) + qsg*snowfrac
          qcg = qcgs*(one-snowfrac) + qcg*snowfrac
          sublim = eeta
          eeta = eetas*(one-snowfrac) + eeta*snowfrac
          qfx = qfxs*(one-snowfrac) + qfx*snowfrac
          hfx = hfxs*(one-snowfrac) + hfx*snowfrac
          s = ss*(one-snowfrac) + s*snowfrac
          prcpl = prcpls*(one-snowfrac) + prcpl*snowfrac
          fltot = fltots*(one-snowfrac) + fltot*snowfrac
          ALB   = MAX(keep_snow_albedo*alb,              &
                  MIN((albice + (alb - alb_snow_free) * snowfrac), alb))
          Emiss = MAX(keep_snow_albedo*emissn,           &
                  MIN((emiss_snowfree +                  &
              (emissn - emiss_snowfree) * snowfrac), emissn))
          runoff1 = runoff1s*(one-snowfrac) + runoff1*snowfrac
          runoff2 = runoff2s*(one-snowfrac) + runoff2*snowfrac
    IF (debug_print ) THEN
      print *,'SOILT combined on ice', soilt
    ENDIF
       endif      
     endif ! snow_mosaic = 1.
 
      !-- 13 jan 2022
      !  update snow fraction after melting (Swenson, S.C. and Lawrence, 2012,
      !  JGR, DOI:10.1029/2012MS000165 
      !  
      !if (snwe > 0.) then
      !  if(smelt > 0.) then
        !update snow fraction after melting
          !n_melt = 200./max(10.,topo_std)
      !    snowfrac = max(0.,snowfrac - (acos(min(1.,(2.*(smelt*delt/snwe) -
      !    1.)))/piconst)**10)
          !snowfrac = 1. - (acos(min(1.,(2.*(smelt*delt/snwe) -
          !1.)))/piconst)**10.
      !    if(i==744.and.j==514 .or. i==924.and.j==568)then
           !print *,'smr,n_melt,topo_std', smr,n_melt,topo_std
      !     print *,'3 - snowfrac end', i,j,ktau,snowfrac,smelt*delt, snwe,
      !     piconst
      !    endif
      !  endif
      !else
      !  snowfrac = 0.
      !endif
      ! 
      !-- The NY07 parameterization gives more realistic snow cover fraction
      !   than SL12
      !-- 13 Jan 2022
      !-- update snow fraction after metlting (Niu, G.-Y., and Yang, Z.-L. 2007,
      !JGR,
      !   DOI:10.1029/2007JD008674)
      !   Limit on znt (<0.25) is needed to avoid very small snow fractions in the
      !   forested areas with large roughness

    IF(snhei == zero) then
      !--- all snow is melted
      iland=ivgtyp
      snowfrac = zero
      alb = alb_snow_free
      emiss = emiss_snowfree
    ELSE
    !-- update snow cover after possible melting
        m = one  ! m=1.6_kind_phys in Niu&Yang, m=1 in CLM
      if(isncovr_opt == 1) then
        snowfrac=min(one,snhei/(2._kind_phys*snhei_crit))
      elseif(isncovr_opt == 2) then
      !-- isncovr_opt=2
        snowfrac=min(one,snhei/(2._kind_phys*snhei_crit))
        if(ivgtyp == glacier .or. ivgtyp == bare) then
        !-- sparsely vegetated or land ice
          snowfrac2 = tanh( snhei/(2.5_kind_phys * 0.2_kind_phys *(rhosn/rhonewsn)**m))
        else
          !-- Niu&Yang: znt=0.01 m for 1 degree (100km) resolution tests
          !  on 3-km scale use actual roughness, but not higher than 0.2 m.
          !  The factor is 20 for forests (~100/dx = 33.)
          snowfrac2 = tanh( snhei/(2.5 *min(0.2,znt) *(rhosn/rhonewsn)**m))
        endif
        !-- snow fraction is average between method 1 and 2
        snowfrac = 0.5_kind_phys*(snowfrac+snowfrac2)
      else
      !-- isncovr_opt=3
        !m = msnf ! vegetation dependent facsnf/msnf from Noah MP
        !-- for RRFS a factor 10. was added to 'facsnf' to get reasonal values of
        !   snow cover fractions on the 3-km scale.
        !   This factor is scale dependent.
        snowfrac = tanh( snhei/(10._kind_phys * facsnf *(rhosn/rhonewsn)**m)) 
      endif

      !-- due to steep slopes and blown snow, limit snow fraction in the
      !-- mountains ( Swiss weather model)
      if(hgt > 2500._kind_phys .and. ivgtyp == glacier) snowfrac=min(0.85_kind_phys,snowfrac)

      if(ivgtyp == urban) snowfrac=min(0.75_kind_phys,snowfrac)

!  run-total accumulated snow based on snowfall and snowmelt in [mm]

       IF (debug_print ) then
       !if (abs(xlat-testptlat).lt.0.2 .and. abs(xlon-testptlon).lt.0.2)then
         print *,'Snowfallac xlat, xlon',xlat,xlon
         print *,'newsn [m],rhonewsn,newsnowratio=',newsn,rhonewsn,newsnowratio
         print *,'Time-step newsn depth [m], swe [m]',newsn,newsn*rhonewsn
         print *,'Time-step smelt: swe [m]' ,smelt*delt
         print *,'Time-step sublim: swe,[kg m-2]',sublim*delt
       endif

      snowfallac = snowfallac + newsn * 1.e3_kind_phys    ! accumulated snow depth [mm], using variable snow density

       IF (debug_print ) THEN
       !if (abs(xlat-testptlat).lt.0.2 .and. abs(xlon-testptlon).lt.0.2)then
         print *,'snowfallac,snhei,snwe',snowfallac,snhei,snwe
       endif
    ENDIF

   ELSE
!--- no snow
           snheiprint=zero
           snweprint=zero
           smelt=zero

!--------------
         T3      = STBOLT*SOILT*SOILT*SOILT
         UPFLUX  = T3 *SOILT
         XINET   = EMISS*(GLW-UPFLUX)
         RNET    = GSWnew + XINET
    IF (debug_print ) THEN
     print *,'NO snow on the ground GSWnew -',GSWnew,'RNET=',rnet
    ENDIF

       if(SEAICE .LT. 0.5_kind_phys) then
!  LAND
         IF ( add_fire_heat_flux .and. fire_heat_flux>0) then ! JLS
          IF (debug_print ) THEN
           print *,'RNET no snow, fire_heat_flux, xlat/xlon',RNET, fire_heat_flux,xlat,xlon
          endif
            RNET = RNET + fire_heat_flux
         ENDIF

         CALL SOIL(debug_print,xlat, xlon, testptlat, testptlon,&
!--- input variables
            i,j,iland,isoil,delt,ktau,conflx,nzs,nddzs,nroot,   &
            PRCPMS,RAINF,PATM,QVATM,QCATM,GLW,GSWnew,GSWin,     &
            EMISS,RNET,QKMS,TKMS,PC,cst,drip,infwater,          &
            rho,vegfrac,lai,myj,                                &
!--- soil fixed fields 
            QWRTZ,rhocs,dqm,qmin,ref,wilt,                      &
            psis,bclh,ksat,sat,cn,                              &
            zsmain,zshalf,DTDZS,DTDZS2,tbq,                     &
!--- constants
            lv,CP,rovcp,G0,cw,stbolt,tabs,                      &
            KQWRTZ,KICE,KWT,                                    &
!--- output variables
            soilm1d,ts1d,smfrkeep,keepfr,                       &
            dew,soilt,qvg,qsg,qcg,edir1,ec1,                    &
            ett1,eeta,qfx,hfx,s,evapl,prcpl,fltot,runoff1,      &
            runoff2,mavail,soilice,soiliqw,                     &
            infiltr,smf)
        else
! SEA ICE
! If current ice albedo is not the same as from the previous time step, then
! update GSW, ALB and RNET for surface energy budget
         if(ALB.ne.ALBice) GSWnew=GSW/(one-ALB)*(one-ALBice)
         alb=albice
         RNET    = GSWnew + XINET

          CALL SICE(debug_print,xlat,xlon,                       &
!--- input variables
            i,j,iland,isoil,delt,ktau,conflx,nzs,nddzs,nroot,   &
            PRCPMS,RAINF,PATM,QVATM,QCATM,GLW,GSWnew,           &
            EMISS,RNET,QKMS,TKMS,rho,myj,                       &
!--- sea ice parameters
            tice,rhosice,capice,thdifice,                       &
            zsmain,zshalf,DTDZS,DTDZS2,tbq,                     &
!--- constants
            lv,CP,rovcp,cw,stbolt,tabs,                         &
!--- output variables
            ts1d,dew,soilt,qvg,qsg,qcg,                         &
            eeta,qfx,hfx,s,evapl,prcpl,fltot                    &
                                                                )
           edir1 = eeta*1.e-3_kind_phys
           ec1 = zero
           ett1 = zero
           runoff1 = prcpms
           runoff2 = zero
           mavail = one
           infiltr = zero
           cst = zero
            do k=1,nzs
               soilm1d(k)= one
               soiliqw(k)= zero
               soilice(k)= one
               smfrkeep(k)= one
               keepfr(k)= zero
            enddo
        endif

        ENDIF

!---------------------------------------------------------------
   END SUBROUTINE SFCTMP
!---------------------------------------------------------------

!>\ingroup lsm_ruc_group
!! This function computes water vapor mixing ratio at saturation from
!! the precomputed table and a given temperature.
       FUNCTION QSN(TN,T)
!****************************************************************
   real (kind_phys),     DIMENSION(1:5001),  INTENT(IN   )   ::  T
   real (kind_phys),     INTENT(IN  )   ::  TN

      real (kind_phys)    QSN, R,R1,R2
      INTEGER I

       R=(TN-173.15_kind_dbl_prec)/.05_kind_dbl_prec+one
       I=INT(R)
       IF(I.GE.1) goto 10
       I=1
       R=1.
  10   IF(I.LE.5000) GOTO 20
       I=5000
       R=5001._kind_dbl_prec
  20   R1=T(I)
       R2=R-I
       QSN=(T(I+1)-R1)*R2 + R1
!-----------------------------------------------------------------------
  END FUNCTION QSN
!------------------------------------------------------------------------

!>\ingroup lsm_ruc_group
!> This subroutine calculates energy and moisture budget for vegetated surfaces
!! without snow, heat diffusion and Richards eqns in soil.
        SUBROUTINE SOIL (debug_print,xlat,xlon,testptlat,testptlon,&
            i,j,iland,isoil,delt,ktau,conflx,nzs,nddzs,nroot,& !--- input variables
            PRCPMS,RAINF,PATM,QVATM,QCATM,                   &
            GLW,GSW,GSWin,EMISS,RNET,                        &
            QKMS,TKMS,PC,cst,drip,infwater,rho,vegfrac,lai,  &
            myj,                                             &
            QWRTZ,rhocs,dqm,qmin,ref,wilt,psis,bclh,ksat,    & !--- soil fixed fields
            sat,cn,zsmain,zshalf,DTDZS,DTDZS2,tbq,           &
            xlv,CP,rovcp,G0_P,cw,stbolt,TABS,                & !--- constants
            KQWRTZ,KICE,KWT,                                 &
            soilmois,tso,smfrkeep,keepfr,                    & !--- output variables
            dew,soilt,qvg,qsg,qcg,                           &
            edir1,ec1,ett1,eeta,qfx,hfx,s,evapl,             &
            prcpl,fltot,runoff1,runoff2,mavail,soilice,      &
            soiliqw,infiltrp,smf)

!*************************************************************
!   Energy and moisture budget for vegetated surfaces 
!   without snow, heat diffusion and Richards eqns. in
!   soil
!
!     DELT - time step (s)
!     ktau - number of time step
!     CONFLX - depth of constant flux layer (m)
!     J,I - the location of grid point
!     IME, JME, KME, NZS - dimensions of the domain
!     NROOT - number of levels within the root zone
!     PRCPMS - precipitation rate in m/s
!     PATM - pressure [bar]
!     QVATM,QCATM - cloud and water vapor mixing ratio (kg/kg)
!                   at the first atm. level
!     GLW, GSW - incoming longwave and absorbed shortwave
!                radiation at the surface (W/m^2)
!     EMISS,RNET - emissivity of the ground surface (0-1) and net
!                  radiation at the surface (W/m^2)
!     QKMS - exchange coefficient for water vapor in the
!              surface layer (m/s)
!     TKMS - exchange coefficient for heat in the surface
!              layer (m/s)
!     PC - plant coefficient (resistance) (0-1)
!     RHO - density of atmosphere near sueface (kg/m^3)
!     VEGFRAC - greeness fraction
!     RHOCS - volumetric heat capacity of dry soil
!     DQM, QMIN - porosity minus residual soil moisture QMIN (m^3/m^3)
!     REF, WILT - field capacity soil moisture and the
!                 wilting point (m^3/m^3)
!     PSIS - matrix potential at saturation (m)
!     BCLH - exponent for Clapp-Hornberger parameterization
!     KSAT - saturated hydraulic conductivity (m/s)
!     SAT - maximum value of water intercepted by canopy (m)
!     CN - exponent for calculation of canopy water
!     ZSMAIN - main levels in soil (m)
!     ZSHALF - middle of the soil layers (m)
!     DTDZS,DTDZS2 - dt/(2.*dzshalf*dzmain) and dt/dzshalf in soil
!     TBQ - table to define saturated mixing ration
!           of water vapor for given temperature and pressure
!     SOILMOIS,TSO - soil moisture (m^3/m^3) and temperature (K)
!     DEW -  dew in kg/m^2s
!     SOILT - skin temperature (K)
!     QSG,QVG,QCG - saturated mixing ratio, mixing ratio of
!                   water vapor and cloud at the ground
!                   surface, respectively (kg/kg)
!     EDIR1, EC1, ETT1, EETA - direct evaporation, evaporation of
!            canopy water, transpiration in kg m-2 s-1 and total
!            evaporation in m s-1.
!     QFX, HFX - latent and sensible heat fluxes (W/m^2)
!     S - soil heat flux in the top layer (W/m^2)
!     RUNOFF - surface runoff (m/s)
!     RUNOFF2 - underground runoff (m)
!     MAVAIL - moisture availability in the top soil layer (0-1)
!     INFILTRP - infiltration flux from the top of soil domain (m/s)
!
!*****************************************************************
        IMPLICIT NONE
!-----------------------------------------------------------------

!--- input variables

   LOGICAL,  INTENT(IN   )   ::  debug_print
   INTEGER,  INTENT(IN   )   ::  nroot,ktau,nzs                , &
                                 nddzs                    !nddzs=2*(nzs-2)
   INTEGER,  INTENT(IN   )   ::  i,j,iland,isoil
   real (kind_phys),     INTENT(IN   )   ::  DELT,CONFLX
   real (kind_phys),     INTENT(IN   )   ::  xlat,xlon,testptlat,testptlon
   LOGICAL,  INTENT(IN   )   ::  myj
!--- 3-D Atmospheric variables
   real (kind_phys),                                             &
            INTENT(IN   )    ::                            PATM, &
                                                          QVATM, &
                                                          QCATM
!--- 2-D variables
   real (kind_phys),                                             &
            INTENT(IN   )    ::                             GLW, &
                                                            GSW, &
                                                          GSWin, &
                                                          EMISS, &
                                                            RHO, &
                                                             PC, &
                                                        VEGFRAC, &
                                                            lai, &
                                                       infwater, &
                                                           QKMS, &
                                                           TKMS

!--- soil properties
   real (kind_phys),                                             &
            INTENT(IN   )    ::                           RHOCS, &
                                                           BCLH, &
                                                            DQM, &
                                                           KSAT, &
                                                           PSIS, &
                                                           QMIN, &
                                                          QWRTZ, &
                                                            REF, &
                                                           WILT

   real (kind_phys),     INTENT(IN   )   ::                  CN, &
                                                             CW, &
                                                         KQWRTZ, &
                                                           KICE, &
                                                            KWT, &
                                                            XLV, &
                                                            g0_p


   real (kind_phys),     DIMENSION(1:NZS), INTENT(IN) :: ZSMAIN, &
                                                         ZSHALF, &
                                                         DTDZS2

   real (kind_phys),     DIMENSION(1:NDDZS), INTENT(IN)  :: DTDZS

   real (kind_phys),     DIMENSION(1:5001), INTENT(IN)  ::  TBQ


!--- input/output variables
!-------- 3-d soil moisture and temperature
   real (kind_phys),     DIMENSION( 1:nzs )                    , &
             INTENT(INOUT)   ::                             TSO, &
                                                       SOILMOIS, &
                                                       SMFRKEEP

   real (kind_phys),     DIMENSION( 1:nzs )                    , &
             INTENT(INOUT)   ::                          KEEPFR

!-------- 2-d variables
   real (kind_phys),                                             &
             INTENT(INOUT)   ::                             DEW, &
                                                            CST, &
                                                           DRIP, &
                                                          EDIR1, &
                                                            EC1, &
                                                           ETT1, &
                                                           EETA, &
                                                          EVAPL, &
                                                          PRCPL, &
                                                         MAVAIL, &
                                                            QVG, &
                                                            QSG, &
                                                            QCG, &
                                                           RNET, &
                                                            QFX, &
                                                            HFX, &
                                                              S, &
                                                            SAT, &
                                                        RUNOFF1, &
                                                        RUNOFF2, &
                                                          SOILT

!-------- 1-d variables
   real (kind_phys),  DIMENSION(1:NZS), INTENT(OUT)  :: SOILICE, &
                                                        SOILIQW

!--- Local variables

   real (kind_phys)    ::  INFILTRP, transum                   , &
               RAINF,  PRCPMS                                  , &
               TABS, T3, UPFLUX, XINET
   real (kind_phys)    ::  CP,rovcp,G0,LV,STBOLT,xlmelt,dzstop , &
               can,epot,fac,fltot,ft,fq,hft                    , &
               q1,ras,sph                                      , &
               trans,zn,ci,cvw,tln,tavln,pi                    , &
               DD1,CMC2MS,DRYCAN,WETCAN                        , &
               INFMAX,RIW, X
   real (kind_phys), DIMENSION(1:NZS) :: transp,cap,diffu,hydro, &
                                   thdif,tranf,tav,soilmoism   , &
                                   soilicem,soiliqwm,detal     , &
                                   fwsat,lwsat,told,smold

   real (kind_phys)                        ::  soiltold,smf
   real (kind_phys)    :: soilres, alfa, fex, fex_fc, fc, psit

   INTEGER ::  nzs1,nzs2,k

!-----------------------------------------------------------------

!-- define constants
        CI=RHOICE*sheatice
        XLMELT=con_hfus
        cvw=cw

        prcpl=prcpms

        smf = zero
        soiltold = soilt

        wetcan= zero
        drycan= one

!--- Initializing local arrays
        DO K=1,NZS
          TRANSP   (K)=zero
          soilmoism(k)=zero
          soilice  (k)=zero
          soiliqw  (k)=zero
          soilicem (k)=zero
          soiliqwm (k)=zero
          lwsat    (k)=zero
          fwsat    (k)=zero
          tav      (k)=zero
          cap      (k)=zero
          thdif    (k)=zero
          diffu    (k)=zero
          hydro    (k)=zero
          tranf    (k)=zero
          detal    (k)=zero
          told     (k)=zero
          smold    (k)=zero
        ENDDO

          NZS1=NZS-1
          NZS2=NZS-2
        dzstop=one/(zsmain(2)-zsmain(1))
        RAS=RHO*1.E-3_kind_phys ! rho/rhowater
        RIW=rhoice*1.e-3_kind_phys ! rhoice/rhowater

!--- Computation of volumetric content of ice in soil 

         DO K=1,NZS
!- main levels
         tln=log(tso(k)/tfrz)
         if(tln.lt.zero) then
           soiliqw(k)=(dqm+qmin)*(XLMELT*                        &
         (tso(k)-tfrz)/tso(k)/grav/psis)                         &
          **(-one/bclh)-qmin
           soiliqw(k)=max(zero,soiliqw(k))
           soiliqw(k)=min(soiliqw(k),soilmois(k))
           soilice(k)=(soilmois(k)-soiliqw(k))/RIW

!---- melting and freezing is balanced, soil ice cannot increase
       if(keepfr(k).eq.one) then
           soilice(k)=min(soilice(k),smfrkeep(k))
           soiliqw(k)=max(zero,soilmois(k)-soilice(k)*riw)
       endif

         else
           soilice(k)=zero
           soiliqw(k)=soilmois(k)
         endif

          ENDDO

          DO K=1,NZS1
!- middle of soil layers
         tav(k)=0.5_kind_phys*(tso(k)+tso(k+1))
         soilmoism(k)=0.5_kind_phys*(soilmois(k)+soilmois(k+1))
         tavln=log(tav(k)/tfrz)

         if(tavln.lt.zero) then
           soiliqwm(k)=(dqm+qmin)*(XLMELT*                       &
         (tav(k)-tfrz)/tav(k)/grav/psis)                         &
          **(-one/bclh)-qmin
           fwsat(k)=dqm-soiliqwm(k)
           lwsat(k)=soiliqwm(k)+qmin
           soiliqwm(k)=max(zero,soiliqwm(k))
           soiliqwm(k)=min(soiliqwm(k), soilmoism(k))
           soilicem(k)=(soilmoism(k)-soiliqwm(k))/riw
!---- melting and freezing is balanced, soil ice cannot increase
       if(keepfr(k).eq.one) then
           soilicem(k)=min(soilicem(k),                          &
                   0.5_kind_phys*(smfrkeep(k)+smfrkeep(k+1)))
           soiliqwm(k)=max(zero,soilmoism(k)-soilicem(k)*riw)
           fwsat(k)=dqm-soiliqwm(k)
           lwsat(k)=soiliqwm(k)+qmin
       endif

         else
           soilicem(k)=zero
           soiliqwm(k)=soilmoism(k)
           lwsat(k)=dqm+qmin
           fwsat(k)=zero
         endif

          ENDDO

          do k=1,nzs
           if(soilice(k).gt.zero) then
             smfrkeep(k)=soilice(k)
           else
             smfrkeep(k)=soilmois(k)/riw
           endif
          enddo

!******************************************************************
! SOILPROP computes thermal diffusivity, and diffusional and
!          hydraulic condeuctivities
!******************************************************************
          CALL SOILPROP( debug_print,                             &
               xlat, xlon, testptlat, testptlon,                  &
!--- input variables
               nzs,fwsat,lwsat,tav,keepfr,                        &
               soilmois,soiliqw,soilice,                          &
               soilmoism,soiliqwm,soilicem,                       &
!--- soil fixed fields
               QWRTZ,rhocs,dqm,qmin,psis,bclh,ksat,               &
!--- constants
               riw,xlmelt,CP,G0_P,cvw,ci,                         &
               kqwrtz,kice,kwt,                                   &
!--- output variables
               thdif,diffu,hydro,cap)

!********************************************************************
!--- CALCULATION OF CANOPY WATER (Smirnova et al., 1996, EQ.16) AND DEW 
 
        FQ=QKMS

        Q1=-QKMS*RAS*(QVATM - QSG)

        DEW=zero
        IF(QVATM.GE.QSG)THEN
          DEW=FQ*(QVATM-QSG)
        ENDIF

!--- WETCAN is the fraction of vegetated area covered by canopy
!--- water, and DRYCAN is the fraction of vegetated area where
!--- transpiration may take place.

          WETCAN=min(0.25_kind_phys,max(zero,(CST/SAT))**CN)
          DRYCAN=one-WETCAN

!**************************************************************
!  TRANSF computes transpiration function
!**************************************************************
           CALL TRANSF(debug_print,                           &
              xlat, xlon, testptlat, testptlon,               &
!--- input variables
              nzs,nroot,soiliqw,tabs,lai,gswin,               &
!--- soil fixed fields
              dqm,qmin,ref,wilt,zshalf,pc,iland,              &
!--- output variables
              tranf,transum)

!--- Save soil temp and moisture from the beginning of time step
          do k=1,nzs
           told(k)=tso(k)
           smold(k)=soilmois(k)
          enddo

! Sakaguchi and Zeng (2009) - dry soil resistance to evaporation
!      if (vgtype==11) then   ! MODIS wetland
        alfa=one
!      else
        fex=min(one,soilmois(1)/dqm)
        fex=max(fex,0.01_kind_phys)
        psit=psis*fex ** (-bclh)
        psit = max(-1.e5_kind_phys, psit)
        alfa=min(one,exp(G0_P*psit/r_v/SOILT))
     ! print *,'alfa=',alfa, exp(G0_P*psit/r_v/SOILT)
!      endif
        alfa=one
! field capacity
! 20jun18 - beta in Eq. (5) is called soilres in the code - it limits soil evaporation
! when soil moisture is below field capacity.  [Lee and Pielke, 1992]
! This formulation agrees with observations when top layer is < 2 cm thick.
! Soilres = 1 for snow, glaciers and wetland.
!        fc=ref  - suggested in the paper
!        fc=max(qmin,ref*0.5) ! used prior to 20jun18 change
! Switch from ref*0.5 to ref*0.25 will reduce soil resistance, increase direct
! evaporation, effects sparsely vegetated areas--> cooler during the day
!        fc=max(qmin,ref*0.25)  ! 
! For now we'll go back to ref*0.5
! 3feb21 - in RRFS testing (fv3-based), ref*0.5 gives too much direct
!          evaporation. Therefore , it is replaced with ref*0.7.
        fc=ref
        fex_fc=one
      if((soilmois(1)+qmin) > fc .or. (qvatm-qvg) > zero) then
        soilres = one
      else
        fex_fc=min(one,(soilmois(1)+qmin)/fc)
        fex_fc=max(fex_fc,0.01_kind_phys)
        soilres=0.25_kind_phys*(one-cos(piconst*fex_fc))**2._kind_phys
      endif
    IF ( debug_print ) THEN
     print *,'piconst=',piconst
     print *,'fex,psit,psis,bclh,g0_p,r_v,soilt,alfa,mavail,soilmois(1),fc,ref,soilres,fex_fc', &
              fex,psit,psis,bclh,g0_p,r_v,soilt,alfa,mavail,soilmois(1),fc,ref,soilres,fex_fc
    endif

!**************************************************************
!  SOILTEMP soilves heat budget and diffusion eqn. in soil
!**************************************************************

        CALL SOILTEMP(debug_print,xlat,xlon,testptlat,testptlon,&
!--- input variables
             i,j,iland,isoil,                                 &
             delt,ktau,conflx,nzs,nddzs,nroot,                &
             PRCPMS,RAINF,                                    &
             PATM,TABS,QVATM,QCATM,EMISS,RNET,                &
             QKMS,TKMS,PC,rho,vegfrac, lai,                   &
             thdif,cap,drycan,wetcan,                         & 
             transum,dew,mavail,soilres,alfa,                 &
!--- soil fixed fields
             dqm,qmin,bclh,zsmain,zshalf,DTDZS,tbq,           &
!--- constants
             xlv,CP,G0_P,cvw,stbolt,                          &
!--- output variables
             tso,soilt,qvg,qsg,qcg,x)

!************************************************************************

!--- CALCULATION OF DEW USING NEW VALUE OF QSG OR TRANSP IF NO DEW
        ETT1=zero
        DEW=zero

        IF(QVATM.GE.QSG)THEN
          DEW=QKMS*(QVATM-QSG)
          ETT1=zero
          DO K=1,NZS
            TRANSP(K)=zero
          ENDDO
        ELSE

          DO K=1,NROOT
            TRANSP(K)=VEGFRAC*RAS*QKMS*                       &
                    (QVATM-QSG)*                              &
                    TRANF(K)*DRYCAN/ZSHALF(NROOT+1)
               IF(TRANSP(K).GT.zero) TRANSP(K)=zero
            ETT1=ETT1-TRANSP(K)
          ENDDO
          DO k=nroot+1,nzs
            transp(k)=zero
          enddo
        ENDIF

!-- Recalculate volumetric content of frozen water in soil
         DO K=1,NZS
!- main levels
           tln=log(tso(k)/tfrz)
         if(tln.lt.zero) then
           soiliqw(k)=(dqm+qmin)*(XLMELT*                     &
          (tso(k)-tfrz)/tso(k)/grav/psis)                     & 
           **(-one/bclh)-qmin
           soiliqw(k)=max(zero,soiliqw(k))
           soiliqw(k)=min(soiliqw(k),soilmois(k))
           soilice(k)=(soilmois(k)-soiliqw(k))/riw
!---- melting and freezing is balanced, soil ice cannot increase
       if(keepfr(k).eq.one) then
           soilice(k)=min(soilice(k),smfrkeep(k))
           soiliqw(k)=max(zero,soilmois(k)-soilice(k)*riw)
       endif

         else
           soilice(k)=zero
           soiliqw(k)=soilmois(k)
         endif
         ENDDO

!*************************************************************************
! SOILMOIST solves moisture budget (Smirnova et al., 1996, EQ.22,28) 
!           and Richards eqn.
!*************************************************************************
          CALL SOILMOIST (debug_print,                         &
               xlat, xlon, testptlat, testptlon,               &
!-- input
               delt,nzs,nddzs,DTDZS,DTDZS2,RIW,                &
               zsmain,zshalf,diffu,hydro,                      &
               QSG,QVG,QCG,QCATM,QVATM,-infwater,              &
               QKMS,TRANSP,DRIP,DEW,zero,SOILICE,VEGFRAC,      &
               zero,soilres,                                   &
!-- soil properties
               DQM,QMIN,REF,KSAT,RAS,INFMAX,                   &
!-- output
               SOILMOIS,SOILIQW,MAVAIL,RUNOFF1,                &
               RUNOFF2,INFILTRP)
        
!--- KEEPFR is 1 when the temperature and moisture in soil
!--- are both increasing. In this case soil ice should not
!--- be increasing according to the freezing curve.
!--- Some part of ice is melted, but additional water is
!--- getting frozen. Thus, only structure of frozen soil is
!--- changed, and phase changes are not affecting the heat
!--- transfer. This situation may happen when it rains on the
!--- frozen soil.
 
        do k=1,nzs
       if (soilice(k).gt.zero) then
          if(tso(k).gt.told(k).and.soilmois(k).gt.smold(k)) then
              keepfr(k)=one
          else
              keepfr(k)=zero
          endif
       endif
        enddo

!--- THE DIAGNOSTICS OF SURFACE FLUXES 

          T3      = STBOLT*SOILTold*SOILTold*SOILTold
          UPFLUX  = T3 * 0.5_kind_phys*(SOILTold+SOILT)
          XINET   = EMISS*(GLW-UPFLUX)
          HFT=-TKMS*CP*RHO*(TABS-SOILT)
          HFX=-TKMS*CP*RHO*(TABS-SOILT)                        &
               *(P1000mb*0.00001_kind_phys/Patm)**ROVCP
          Q1=-QKMS*RAS*(QVATM - QSG)

          CMC2MS = zero
        IF (Q1.LE.zero) THEN
! ---  condensation
          EC1= zero
          EDIR1= zero
          ETT1= zero
     if(myj) then
!-- moisture flux for coupling with MYJ PBL
          EETA=-QKMS*RAS*(QVATM/(one+QVATM) - QSG/(one+QSG))*rhowater
          CST= CST-EETA*DELT*vegfrac
    IF (debug_print ) THEN
!!!    IF(i.eq.374.and.j.eq.310.or. EETA.gt.0.0004) then
        print *,'Cond MYJ EETA',eeta,eeta*xlv, i,j
    ENDIF
     else ! myj
!-- actual moisture flux from RUC LSM
          EETA= - RHO*DEW
          CST=CST+DELT*DEW*RAS * vegfrac
    IF (debug_print ) THEN
!    IF(i.eq.374.and.j.eq.310.or. EETA.gt.0.0004) then
       print *,'Cond RUC LSM EETA',EETA,eeta*xlv, i,j
    ENDIF
     endif ! myj
          QFX= XLV*EETA
          EETA= - RHO*DEW
        ELSE
! ---  evaporation
          EDIR1 =-soilres*(one-vegfrac)*QKMS*RAS*                      &
                  (QVATM-QVG)
          CMC2MS=CST/DELT*RAS
          EC1 = Q1 * WETCAN * vegfrac
    IF (debug_print ) THEN
     IF(i.eq.440.and.j.eq.180.or. QFX.gt.1000..or.i.eq.417.and.j.eq.540) then
       print *,'CST before update=',cst
       print *,'EC1=',EC1,'CMC2MS=',CMC2MS
     ENDIF
    ENDIF

          CST=max(zero,CST-EC1 * DELT)

     if (myj) then
!-- moisture flux for coupling with MYJ PBL
          EETA=-soilres*QKMS*RAS*(QVATM/(one+QVATM) - QVG/(one+QVG))*rhowater
     else ! myj
    IF (debug_print ) THEN
!    IF(i.eq.440.and.j.eq.180.or. QFX.gt.1000..or.i.eq.417.and.j.eq.540) then
       print *,'QKMS,RAS,QVATM/(one+QVATM),QVG/(one+QVG),QSG ', &
                QKMS,RAS,QVATM/(one+QVATM),QVG/(one+QVG),QSG
       print *,'Q1*(1.-vegfrac),EDIR1',Q1*(one-vegfrac),EDIR1
       print *,'CST,WETCAN,DRYCAN',CST,WETCAN,DRYCAN
       print *,'EC1=',EC1,'ETT1=',ETT1,'CMC2MS=',CMC2MS,'CMC2MS*ras=',CMC2MS*ras
    ENDIF
!-- actual moisture flux from RUC LSM
          EETA = (EDIR1 + EC1 + ETT1)*rhowater
    IF (debug_print ) THEN
!    IF(i.eq.374.and.j.eq.310.or. EETA.gt.0.0004) then
        print *,'RUC LSM EETA',EETA,eeta*xlv
    ENDIF
     endif ! myj
          QFX= XLV * EETA
          EETA = (EDIR1 + EC1 + ETT1)*rhowater
        ENDIF
    IF (debug_print ) THEN
     print *,'potential temp HFT ',HFT
     print *,'abs temp HFX ',HFX
    ENDIF

          EVAPL=EETA
          S=THDIF(1)*CAP(1)*DZSTOP*(TSO(1)-TSO(2))
! Energy budget
          FLTOT=RNET-HFT-XLV*EETA-S-X
    IF (debug_print ) THEN
!    IF(i.eq.440.and.j.eq.180 .or. qfx.gt.1000..or.i.eq.417.and.j.eq.540) then
       print *,'SOIL - FLTOT,RNET,HFT,QFX,S,X=',i,j,FLTOT,RNET,HFT,XLV*EETA,s,x
       print *,'edir1,ec1,ett1,mavail,qkms,qvatm,qvg,qsg,vegfrac',&
                edir1,ec1,ett1,mavail,qkms,qvatm,qvg,qsg,vegfrac
    ENDIF
    if(detal(1) .ne. zero) then
! SMF - energy of phase change in the first soil layer
         smf=fltot
    IF (debug_print ) THEN
     print *,'detal(1),xlmelt,soiliqwm(1),delt',detal(1),xlmelt,soiliqwm(1),delt
     print *,'Implicit phase change in the first layer - smf=',smf
    ENDIF
    endif


 222    CONTINUE

 1123    FORMAT(I5,8F12.3)
 1133    FORMAT(I7,8E12.4)
  123   format(i6,f6.2,7f8.1)
  122   FORMAT(1X,2I3,6F8.1,F8.3,F8.2)
!-------------------------------------------------------------------
   END SUBROUTINE SOIL
!-------------------------------------------------------------------

!>\ingroup lsm_ruc_group
!! This subroutine is called for sea ice without accumulated snow
!! on its surface. it solves heat diffusion inside ice and energy
!! budget at the surface of ice. It computes skin temperature and
!! temerature inside sea ice.
        SUBROUTINE SICE ( debug_print,xlat,xlon,                &
            i,j,iland,isoil,delt,ktau,conflx,nzs,nddzs,nroot,   & !--- input variables
            PRCPMS,RAINF,PATM,QVATM,QCATM,GLW,GSW,              &
            EMISS,RNET,QKMS,TKMS,rho,myj,                       &
            tice,rhosice,capice,thdifice,                       & !--- sea ice parameters
            zsmain,zshalf,DTDZS,DTDZS2,tbq,                     &
            xlv,CP,rovcp,cw,stbolt,tabs,                        & !--- constants
            tso,dew,soilt,qvg,qsg,qcg,                          & !--- output variables
            eeta,qfx,hfx,s,evapl,prcpl,fltot                    &
                                                                )

!*****************************************************************
!   Energy budget and  heat diffusion eqns. for
!   sea ice
!*************************************************************

        IMPLICIT NONE
!-----------------------------------------------------------------

!--- input variables

   INTEGER,  INTENT(IN   )   ::  nroot,ktau,nzs                , &
                                 nddzs                    !nddzs=2*(nzs-2)
   INTEGER,  INTENT(IN   )   ::  i,j,iland,isoil
   real (kind_phys),     INTENT(IN   )   ::  DELT,CONFLX,xlat,xlon
   LOGICAL,  INTENT(IN   )   ::  myj, debug_print
!--- 3-D Atmospheric variables
   real (kind_phys),                                             &
            INTENT(IN   )    ::                            PATM, &
                                                          QVATM, &
                                                          QCATM
!--- 2-D variables
   real (kind_phys),                                             &
            INTENT(IN   )    ::                             GLW, &
                                                            GSW, &
                                                          EMISS, &
                                                            RHO, &
                                                           QKMS, &
                                                           TKMS
!--- sea ice properties
   real (kind_phys),    DIMENSION(1:NZS)                       , &
            INTENT(IN   )    ::                                  &
                                                           tice, &
                                                        rhosice, &
                                                         capice, &
                                                       thdifice


   real (kind_phys),     INTENT(IN   )   ::                      &
                                                             CW, &
                                                            XLV


   real (kind_phys),   DIMENSION(1:NZS), INTENT(IN)  ::  ZSMAIN, &
                                                         ZSHALF, &
                                                         DTDZS2

   real (kind_phys),   DIMENSION(1:NDDZS), INTENT(IN)  :: DTDZS

   real (kind_phys),     DIMENSION(1:5001), INTENT(IN)  ::  TBQ


!--- input/output variables
!----soil temperature
   real (kind_phys),     DIMENSION( 1:nzs ),  INTENT(INOUT)   :: TSO
!-------- 2-d variables
   real (kind_phys),                                             &
             INTENT(INOUT)   ::                             DEW, &
                                                           EETA, &
                                                          EVAPL, &
                                                          PRCPL, &
                                                            QVG, &
                                                            QSG, &
                                                            QCG, &
                                                           RNET, &
                                                            QFX, &
                                                            HFX, &
                                                              S, &
                                                          SOILT

!--- Local variables
   real (kind_phys)    ::  x,x1,x2,x4,tn,denom
   real (kind_phys)    ::       RAINF,  PRCPMS                 , &
                                TABS, T3, UPFLUX, XINET

   real (kind_phys)    ::  CP,rovcp,G0,LV,STBOLT,xlmelt,dzstop , &
               epot,fltot,ft,fq,hft,ras,cvw                    

   real (kind_phys) :: FKT,D1,D2,D9,D10,DID,R211,R21,R22,R6,R7,D11, &
                  PI,H,FKQ,R210,AA,BB,PP,Q1,QS1,TS1,TQ2,TX2       , &
                  TDENOM,QGOLD,SNOH

   real (kind_phys)    ::  AA1,RHCS, icemelt


   real (kind_phys),     DIMENSION(1:NZS)  ::   cotso,rhtso

   INTEGER ::  nzs1,nzs2,k,k1,kn,kk

!-----------------------------------------------------------------

!-- define constants
        XLMELT=con_hfus
        cvw=cw

        prcpl=prcpms

          NZS1=NZS-1
          NZS2=NZS-2
        dzstop=1./(zsmain(2)-zsmain(1))
        RAS=RHO*1.E-3_kind_phys

        do k=1,nzs
           cotso(k)=zero
           rhtso(k)=zero
        enddo

        cotso(1)=zero
        rhtso(1)=TSO(NZS)

        DO 33 K=1,NZS2
          KN=NZS-K
          K1=2*KN-3
          X1=DTDZS(K1)*THDIFICE(KN-1)
          X2=DTDZS(K1+1)*THDIFICE(KN)
          FT=TSO(KN)+X1*(TSO(KN-1)-TSO(KN))                             &
             -X2*(TSO(KN)-TSO(KN+1))
          DENOM=1.+X1+X2-X2*cotso(K)
          cotso(K+1)=X1/DENOM
          rhtso(K+1)=(FT+X2*rhtso(K))/DENOM
   33  CONTINUE

!************************************************************************
!--- THE HEAT BALANCE EQUATION (Smirnova et al., 1996, EQ. 21,26)
        RHCS=CAPICE(1)
        H=one
        FKT=TKMS
        D1=cotso(NZS1)
        D2=rhtso(NZS1)
        TN=TSO(1)
        D9=THDIFICE(1)*RHCS*dzstop
        D10=TKMS*CP*RHO
        R211=.5_kind_phys*CONFLX/DELT
        R21=R211*CP*RHO
        R22=.5_kind_phys/(THDIFICE(1)*DELT*dzstop**2)
        R6=EMISS *STBOLT*.5_kind_phys*TN**4
        R7=R6/TN
        D11=RNET+R6
        TDENOM=D9*(one-D1+R22)+D10+R21+R7                             &
              +RAINF*CVW*PRCPMS
        FKQ=QKMS*RHO
        R210=R211*RHO
        AA=XLS*(FKQ+R210)/TDENOM
        BB=(D10*TABS+R21*TN+XLS*(QVATM*FKQ                            &
        +R210*QVG)+D11+D9*(D2+R22*TN)                                 &
        +RAINF*CVW*PRCPMS*max(tfrz,TABS)                              &
         )/TDENOM
        AA1=AA
        PP=PATM*rhowater
        AA1=AA1/PP
    IF (debug_print ) THEN
        PRINT *,' VILKA-SEAICE1'
        print *,'D10,TABS,R21,TN,QVATM,FKQ',                          &
                 D10,TABS,R21,TN,QVATM,FKQ
        print *,'RNET, EMISS, STBOLT, SOILT',RNET, EMISS, STBOLT, SOILT
        print *,'R210,QVG,D11,D9,D2,R22,RAINF,CVW,PRCPMS,TDENOM',     &
                 R210,QVG,D11,D9,D2,R22,RAINF,CVW,PRCPMS,TDENOM
        print *,'tn,aa1,bb,pp,fkq,r210',                              &
                 tn,aa1,bb,pp,fkq,r210
    ENDIF
        QGOLD=QSG
        CALL VILKA(TN,AA1,BB,PP,QS1,TS1,TBQ,KTAU,i,j,iland,isoil,xlat,xlon)
!--- it is saturation over sea ice
        QVG=QS1
        QSG=QS1
        TSO(1)=min(con_tice,TS1)
        QCG=zero
!--- sea ice melting is not included in this simple approach
!--- SOILT - skin temperature
          SOILT=TSO(1)
!---- Final solution for soil temperature - TSO
          DO K=2,NZS
            KK=NZS-K+1
            TSO(K)=min(con_tice,rhtso(KK)+cotso(KK)*TSO(K-1))
          END DO
!--- CALCULATION OF DEW USING NEW VALUE OF QSG OR TRANSP IF NO DEW
        DEW=zero

!--- THE DIAGNOSTICS OF SURFACE FLUXES 
          T3      = STBOLT*TN*TN*TN
          UPFLUX  = T3 *0.5_kind_phys*(TN+SOILT)
          XINET   = EMISS*(GLW-UPFLUX)
          HFT=-TKMS*CP*RHO*(TABS-SOILT)
          HFX=-TKMS*CP*RHO*(TABS-SOILT)                        &
               *(P1000mb*0.00001_kind_phys/Patm)**ROVCP
          Q1=-QKMS*RAS*(QVATM - QSG)
        IF (Q1.LE.zero) THEN
! ---  condensation
     if(myj) then
!-- moisture flux for coupling with MYJ PBL
          EETA=-QKMS*RAS*(QVATM/(1.+QVATM) - QSG/(1.+QSG))*rhowater
    IF (debug_print ) THEN
       print *,'MYJ EETA',eeta
    ENDIF
     else ! myj
!-- actual moisture flux from RUC LSM
          DEW=QKMS*(QVATM-QSG)
          EETA= - RHO*DEW
    IF (debug_print ) THEN
       print *,'RUC LSM EETA',eeta
    ENDIF
     endif ! myj
          QFX= XLS*EETA
          EETA= - RHO*DEW
        ELSE
! ---  evaporation
     if(myj) then
!-- moisture flux for coupling with MYJ PBL
          EETA=-QKMS*RAS*(QVATM/(1.+QVATM) - QVG/(1.+QVG))*rhowater
    IF (debug_print ) THEN
       print *,'MYJ EETA',eeta
    ENDIF
     else ! myj
! to convert from m s-1 to kg m-2 s-1: *rho water=1.e3************
!-- actual moisture flux from RUC LSM
          EETA = Q1*rhowater
    IF (debug_print ) THEN
       print *,'RUC LSM EETA',eeta
    ENDIF
     endif ! myj
          QFX= XLS * EETA
          EETA = Q1*rhowater
        ENDIF
          EVAPL=EETA

          S=THDIFICE(1)*CAPICE(1)*DZSTOP*(TSO(1)-TSO(2))
! heat storage in surface layer
        SNOH=zero
! There is ice melt
         X= (cp*rho*r211+rhcs*zsmain(2)*0.5_kind_phys/delt)*(SOILT-TN) +   &
            XLS*rho*r211*(QSG-QGOLD)
         X=X &
! "heat" from rain
        -RAINF*CVW*PRCPMS*(max(tfrz,TABS)-SOILT)

!-- excess energy spent on sea ice melt
        icemelt=RNET-XLS*EETA -HFT -S -X
    IF (debug_print ) THEN
        print *,'icemelt=',icemelt
    ENDIF

          FLTOT=RNET-XLS*EETA-HFT-S-X-icemelt
    IF (debug_print ) THEN
       print *,'SICE - FLTOT,RNET,HFT,QFX,S,icemelt,X=', &
                       FLTOT,RNET,HFT,XLS*EETA,s,icemelt,X
    ENDIF

!-------------------------------------------------------------------
   END SUBROUTINE SICE
!-------------------------------------------------------------------

!>\ingroup lsm_ruc_group
!! This subroutine is called for snow covered areas of land. It
!! solves energy and moisture budgets on the surface of snow, and 
!! on the interface of snow and soil. It computes skin temperature,
!! snow temperature, snow depth and snow melt.
        SUBROUTINE SNOWSOIL ( debug_print,xlat,xlon,           &
             testptlat,testptlon,                              &
             i,j,isoil,delt,ktau,conflx,nzs,nddzs,nroot,       & !--- input variables
             isncond_opt,isncovr_opt,                          &
             meltfactor,rhonewsn,SNHEI_CRIT,                   & ! new
             ILAND,PRCPMS,RAINF,NEWSNOW,snhei,SNWE,SNOWFRAC,   &
             RHOSN,                                            &
             PATM,QVATM,QCATM,                                 &
             GLW,GSW,GSWin,EMISS,RNET,IVGTYP,                  &
             QKMS,TKMS,PC,cst,drip,infwater,                   &
             rho,vegfrac,alb,znt,lai,                          &
             MYJ,                                              & !--- soil fixed fields
             QWRTZ,rhocs,dqm,qmin,ref,wilt,psis,bclh,ksat,     &
             sat,cn,zsmain,zshalf,DTDZS,DTDZS2,tbq,            &
             xlv,CP,rovcp,G0_P,cw,stbolt,TABS,                 & !--- constants
             KQWRTZ,KICE,KWT,                                  &
             ilnb,snweprint,snheiprint,rsm,                    & !--- output variables
             soilmois,tso,smfrkeep,keepfr,                     &
             dew,soilt,soilt1,tsnav,                           &
             qvg,qsg,qcg,SMELT,SNOH,SNFLX,SNOM,                &
             edir1,ec1,ett1,eeta,qfx,hfx,s,sublim,             &
             prcpl,fltot,runoff1,runoff2,mavail,soilice,       &
             soiliqw,infiltrp                                  )

!***************************************************************
!   Energy and moisture budget for snow, heat diffusion eqns.
!   in snow and soil, Richards eqn. for soil covered with snow
!
!     DELT - time step (s)
!     ktau - numver of time step
!     CONFLX - depth of constant flux layer (m)
!     J,I - the location of grid point
!     IME, JME,  NZS - dimensions of the domain
!     NROOT - number of levels within the root zone
!     PRCPMS - precipitation rate in m/s
!     NEWSNOW - pcpn in soilid form (m)
!     SNHEI, SNWE - snow height and snow water equivalent (m)
!     RHOSN - snow density (kg/m-3)
!     PATM - pressure (bar)
!     QVATM,QCATM - cloud and water vapor mixing ratio
!                   at the first atm. level (kg/kg)
!     GLW, GSW - incoming longwave and absorbed shortwave
!                radiation at the surface (W/m^2)
!     EMISS,RNET - emissivity (0-1) of the ground surface and net
!                  radiation at the surface (W/m^2)
!     QKMS - exchange coefficient for water vapor in the
!              surface layer (m/s)
!     TKMS - exchange coefficient for heat in the surface
!              layer (m/s)
!     PC - plant coefficient (resistance) (0-1)
!     RHO - density of atmosphere near surface (kg/m^3)
!     VEGFRAC - greeness fraction (0-1)
!     RHOCS - volumetric heat capacity of dry soil (J/m^3/K)
!     DQM, QMIN - porosity minus residual soil moisture QMIN (m^3/m^3)
!     REF, WILT - field capacity soil moisture and the
!                 wilting point (m^3/m^3)
!     PSIS - matrix potential at saturation (m)
!     BCLH - exponent for Clapp-Hornberger parameterization
!     KSAT - saturated hydraulic conductivity (m/s)
!     SAT - maximum value of water intercepted by canopy (m)
!     CN - exponent for calculation of canopy water
!     ZSMAIN - main levels in soil (m)
!     ZSHALF - middle of the soil layers (m)
!     DTDZS,DTDZS2 - dt/(2.*dzshalf*dzmain) and dt/dzshalf in soil
!     TBQ - table to define saturated mixing ration
!           of water vapor for given temperature and pressure
!     ilnb - number of layers in snow
!     rsm - liquid water inside snow pack (m)
!     SOILMOIS,TSO - soil moisture (m^3/m^3) and temperature (K)
!     DEW -  dew in (kg/m^2 s)
!     SOILT - skin temperature (K)
!     SOILT1 - snow temperature at 7.5 cm depth (K)
!     TSNAV - average temperature of snow pack (C)
!     QSG,QVG,QCG - saturated mixing ratio, mixing ratio of
!                   water vapor and cloud at the ground
!                   surface, respectively (kg/kg)
!     EDIR1, EC1, ETT1, EETA - direct evaporation, evaporation of
!            canopy water, transpiration (kg m-2 s-1) and total
!            evaporation in (m s-1).
!     QFX, HFX - latent and sensible heat fluxes (W/m^2)
!     S - soil heat flux in the top layer (W/m^2)
!     SUBLIM - snow sublimation (kg/m^2/s)
!     RUNOFF1 - surface runoff (m/s)
!     RUNOFF2 - underground runoff (m)
!     MAVAIL - moisture availability in the top soil layer (0-1)
!     SOILICE - content of soil ice in soil layers (m^3/m^3)
!     SOILIQW - lliquid water in soil layers (m^3/m^3)
!     INFILTRP - infiltration flux from the top of soil domain (m/s)
!     XINET - net long-wave radiation (W/m^2)
!
!*******************************************************************

        IMPLICIT NONE
!-------------------------------------------------------------------
!--- input variables
   LOGICAL,  INTENT(IN   )   ::  debug_print
   INTEGER,  INTENT(IN   )   ::  nroot,ktau,nzs     ,            &
                                 nddzs                         !nddzs=2*(nzs-2)
   INTEGER,  INTENT(IN   )   ::  i,j,isoil,isncond_opt,isncovr_opt

   real (kind_phys),     INTENT(IN   )   ::  DELT,CONFLX,PRCPMS, &
                                 RAINF,NEWSNOW,RHONEWSN,         &
                                 testptlat,testptlon,            &
                                 SNHEI_CRIT,meltfactor,xlat,xlon

   LOGICAL,    INTENT(IN   )    ::     myj

!--- 3-D Atmospheric variables
   real (kind_phys),                                             &
            INTENT(IN   )    ::                            PATM, &
                                                          QVATM, &
                                                          QCATM
!--- 2-D variables
   real (kind_phys)                                            , &
            INTENT(IN   )    ::                             GLW, &
                                                            GSW, &
                                                          GSWin, &
                                                            RHO, &
                                                             PC, &
                                                        VEGFRAC, &
                                                            lai, &
                                                       infwater, &
                                                           QKMS, &
                                                           TKMS

   INTEGER,  INTENT(IN   )   ::                          IVGTYP
!--- soil properties
   real (kind_phys)                                            , &
            INTENT(IN   )    ::                           RHOCS, &
                                                           BCLH, &
                                                            DQM, &
                                                           KSAT, &
                                                           PSIS, &
                                                           QMIN, &
                                                          QWRTZ, &
                                                            REF, &
                                                            SAT, &
                                                           WILT

   real (kind_phys),     INTENT(IN   )   ::                  CN, &
                                                             CW, &
                                                            XLV, &
                                                           G0_P, & 
                                                         KQWRTZ, &
                                                           KICE, &
                                                            KWT 


   real (kind_phys),   DIMENSION(1:NZS), INTENT(IN)  ::  ZSMAIN, &
                                                         ZSHALF, &
                                                         DTDZS2

   real (kind_phys),   DIMENSION(1:NDDZS), INTENT(IN)  :: DTDZS

   real (kind_phys),   DIMENSION(1:5001), INTENT(IN)  ::    TBQ


!--- input/output variables
!-------- 3-d soil moisture and temperature
   real (kind_phys),     DIMENSION(  1:nzs )                   , &
             INTENT(INOUT)   ::                             TSO, &
                                                       SOILMOIS, &
                                                       SMFRKEEP

   real (kind_phys),  DIMENSION( 1:nzs )                       , &
             INTENT(INOUT)   ::                          KEEPFR


   INTEGER,  INTENT(INOUT)    ::                           ILAND


!-------- 2-d variables
   real (kind_phys)                                            , &
             INTENT(INOUT)   ::                             DEW, &
                                                            CST, &
                                                           DRIP, &
                                                          EDIR1, &
                                                            EC1, &
                                                           ETT1, &
                                                           EETA, &
                                                          RHOSN, &
                                                         SUBLIM, &
                                                          PRCPL, &
                                                            ALB, &
                                                          EMISS, &
                                                            ZNT, &
                                                         MAVAIL, &
                                                            QVG, &
                                                            QSG, &
                                                            QCG, &
                                                            QFX, &
                                                            HFX, &
                                                              S, &
                                                        RUNOFF1, &
                                                        RUNOFF2, &
                                                           SNWE, &
                                                          SNHEI, &
                                                          SMELT, &
                                                           SNOM, &
                                                           SNOH, &
                                                          SNFLX, &
                                                          SOILT, &
                                                         SOILT1, &
                                                       SNOWFRAC, &
                                                          TSNAV

   INTEGER, INTENT(INOUT)    ::                            ILNB

!-------- 1-d variables
   real (kind_phys),     DIMENSION(1:NZS), INTENT(OUT)  ::  SOILICE, &
                                                            SOILIQW

   real (kind_phys),     INTENT(OUT)                         :: RSM, &
                                                          SNWEPRINT, &
                                                          SNHEIPRINT
!--- Local variables


   INTEGER ::  nzs1,nzs2,k

   real (kind_phys)    ::  INFILTRP, TRANSUM                   , &
               SNTH, NEWSN                                     , &
               TABS, T3, UPFLUX, XINET                         , &
               BETA, SNWEPR,EPDT,PP
   real (kind_phys)    ::  CP,rovcp,G0,LV,xlvm,STBOLT,xlmelt,dzstop, &
               can,epot,fac,fltot,ft,fq,hft                    , &
               q1,ras,sph                                      , &
               trans,zn,ci,cvw,tln,tavln,pi                    , &
               DD1,CMC2MS,DRYCAN,WETCAN                        , &
               INFMAX,RIW,DELTSN,H,UMVEG

   real (kind_phys), DIMENSION(1:NZS) :: transp,cap,diffu,hydro, &
                                   thdif,tranf,tav,soilmoism   , &
                                   soilicem,soiliqwm,detal     , &
                                   fwsat,lwsat,told,smold
   real (kind_phys)                        ::  soiltold, qgold

   real (kind_phys)                        ::  RNET, X

!-----------------------------------------------------------------

        cvw=cw
        XLMELT=con_hfus
!-- heat of water vapor sublimation
        XLVm=XLV+XLMELT

!--- SNOW flag -- ISICE
!--- DELTSN - is the threshold for splitting the snow layer into 2 layers.
!--- With snow density 400 kg/m^3, this threshold is equal to 7.5 cm,
!--- equivalent to 0.03 m SNWE. For other snow densities the threshold is
!--- computed using SNWE=0.03 m and current snow density.
!--- SNTH - the threshold below which the snow layer is combined with
!--- the top soil layer. SNTH is computed using snwe=0.016 m, and
!--- equals 4 cm for snow density 400 kg/m^3.

!save SOILT and QVG
       soiltold=soilt
       qgold=qvg

       x=zero

! increase thinkness of top snow layer from 3 cm SWE to 5 cm SWE
           DELTSN=0.05_kind_phys*rhowater/rhosn
           snth=0.01_kind_phys*rhowater/rhosn

! For 2-layer snow model when the snow depth is marginally higher than DELTSN,
! reset DELTSN to half of snow depth.
        IF(SNHEI.GE.DELTSN+SNTH) THEN
          if(snhei-deltsn-snth.lt.snth) deltsn=0.5_kind_phys*(snhei-snth)
    IF (debug_print ) THEN
      print *,'DELTSN is changed,deltsn,snhei,snth',i,j,deltsn,snhei,snth
    ENDIF
        ENDIF 

        CI=RHOICE*sheatice
        RAS=RHO*1.E-3_kind_dbl_prec ! rho/rhowater
        RIW=rhoice*1.e-3_kind_dbl_prec ! rhoice/rhowater
        RSM=zero

        DO K=1,NZS
          TRANSP     (K)=zero
          soilmoism  (k)=zero
          soiliqwm   (k)=zero
          soilice    (k)=zero
          soilicem   (k)=zero
          lwsat      (k)=zero
          fwsat      (k)=zero
          tav        (k)=zero
          cap        (k)=zero
          diffu      (k)=zero
          hydro      (k)=zero
          thdif      (k)=zero
          tranf      (k)=zero
          detal      (k)=zero
          told       (k)=zero
          smold      (k)=zero 
        ENDDO

        snweprint=zero
        snheiprint=zero
        prcpl=prcpms

!*** DELTSN is the depth of the top layer of snow where
!*** there is a temperature gradient, the rest of the snow layer
!*** is considered to have constant temperature


          NZS1=NZS-1
          NZS2=NZS-2
        DZSTOP=one/(zsmain(2)-zsmain(1))

!----- THE CALCULATION OF THERMAL DIFFUSIVITY, DIFFUSIONAL AND ---
!----- HYDRAULIC CONDUCTIVITY (SMIRNOVA ET AL. 1996, EQ.2,5,6) ---
!tgs - the following loop is added to define the amount of frozen
!tgs - water in soil if there is any
         DO K=1,NZS

         tln=log(tso(k)/tfrz)
         if(tln.lt.zero) then
           soiliqw(k)=(dqm+qmin)*(XLMELT*                          &
           (tso(k)-tfrz)/tso(k)/grav/psis)                         &
          **(-one/bclh)-qmin
           soiliqw(k)=max(zero,soiliqw(k))
           soiliqw(k)=min(soiliqw(k),soilmois(k))
           soilice(k)=(soilmois(k)-soiliqw(k))/riw

!---- melting and freezing is balanced, soil ice cannot increase
       if(keepfr(k).eq.1.) then
           soilice(k)=min(soilice(k),smfrkeep(k))
           soiliqw(k)=max(zero,soilmois(k)-soilice(k)*riw)
       endif

         else
           soilice(k)=zero
           soiliqw(k)=soilmois(k)
         endif

          ENDDO

          DO K=1,NZS1

         tav(k)=0.5_kind_phys*(tso(k)+tso(k+1))
         soilmoism(k)=0.5_kind_phys*(soilmois(k)+soilmois(k+1))
         tavln=log(tav(k)/tfrz)

         if(tavln.lt.zero) then
           soiliqwm(k)=(dqm+qmin)*(XLMELT*                         &
         (tav(k)-tfrz)/tav(k)/grav/psis)                           &
          **(-one/bclh)-qmin
           fwsat(k)=dqm-soiliqwm(k)
           lwsat(k)=soiliqwm(k)+qmin
           soiliqwm(k)=max(zero,soiliqwm(k))
           soiliqwm(k)=min(soiliqwm(k), soilmoism(k))
           soilicem(k)=(soilmoism(k)-soiliqwm(k))/riw
!---- melting and freezing is balanced, soil ice cannot increase
       if(keepfr(k).eq.one) then
           soilicem(k)=min(soilicem(k),                            &
                    0.5_kind_phys*(smfrkeep(k)+smfrkeep(k+1)))
           soiliqwm(k)=max(zero,soilmoism(k)-soilicem(k)*riw)
           fwsat(k)=dqm-soiliqwm(k)
           lwsat(k)=soiliqwm(k)+qmin
       endif

         else
           soilicem(k)=zero
           soiliqwm(k)=soilmoism(k)
           lwsat(k)=dqm+qmin
           fwsat(k)=zero

         endif
          ENDDO

          do k=1,nzs
           if(soilice(k).gt.zero) then
             smfrkeep(k)=soilice(k)
           else
             smfrkeep(k)=soilmois(k)/riw
           endif
          enddo

!******************************************************************
! SOILPROP computes thermal diffusivity, and diffusional and
!          hydraulic condeuctivities
!******************************************************************
          CALL SOILPROP(debug_print,                             &
               xlat, xlon, testptlat, testptlon,                 &
!--- input variables
               nzs,fwsat,lwsat,tav,keepfr,                       &
               soilmois,soiliqw,soilice,                         &
               soilmoism,soiliqwm,soilicem,                      &
!--- soil fixed fields
               QWRTZ,rhocs,dqm,qmin,psis,bclh,ksat,              & 
!--- constants
               riw,xlmelt,CP,G0_P,cvw,ci,                        &
               kqwrtz,kice,kwt,                                  &
!--- output variables
               thdif,diffu,hydro,cap)

!******************************************************************** 
!--- CALCULATION OF CANOPY WATER (Smirnova et al., 1996, EQ.16) AND DEW 
 
        SMELT=zero
        H=MAVAIL ! =1. if snowfrac=1

        FQ=QKMS


!--- If vegfrac.ne.0. then part of falling snow can be
!--- intercepted by the canopy. 

        DEW=zero
        UMVEG=one-vegfrac
        EPOT = -FQ*(QVATM-QSG) 

    IF (debug_print ) THEN
      print *,'SNWE after subtracting intercepted snow - snwe=',snwe,vegfrac,cst
    ENDIF

!-- Save SNWE from the previous time step
          SNWEPR=SNWE

!  check if all snow can evaporate during DT
         BETA=one
         EPDT = EPOT * RAS *DELT
         IF(EPDT > zero .and. SNWEPR.LE.EPDT) THEN 
            BETA=SNWEPR/EPDT
            SNWE=zero
         ENDIF

          WETCAN=min(0.25_kind_phys,max(zero,(CST/SAT))**CN)
          DRYCAN=one-WETCAN

!**************************************************************
!  TRANSF computes transpiration function
!**************************************************************
           CALL TRANSF(debug_print,                           &
              xlat, xlon, testptlat, testptlon,               &
!--- input variables
              nzs,nroot,soiliqw,tabs,lai,gswin,               &
!--- soil fixed fields
              dqm,qmin,ref,wilt,zshalf,pc,iland,              & 
!--- output variables
              tranf,transum)

!--- Save soil temp and moisture from the beginning of time step
          do k=1,nzs
           told(k)=tso(k)
           smold(k)=soilmois(k)
          enddo

!**************************************************************
! SNOWTEMP solves heat budget and diffusion eqn. in soil
!**************************************************************

    IF (debug_print ) THEN
print *, 'TSO before calling SNOWTEMP: ', tso
    ENDIF
        CALL SNOWTEMP(debug_print,xlat,xlon,testptlat,testptlon,&
!--- input variables
             i,j,iland,isoil,                                 &
             delt,ktau,conflx,nzs,nddzs,nroot,                &
             isncond_opt,isncovr_opt,                         &
             snwe,snwepr,snhei,newsnow,snowfrac,snhei_crit,   &
             beta,deltsn,snth,rhosn,rhonewsn,meltfactor,      &  ! add meltfactor
             PRCPMS,RAINF,                                    &
             PATM,TABS,QVATM,QCATM,                           &
             GLW,GSW,EMISS,RNET,                              &
             QKMS,TKMS,PC,rho,vegfrac,                        &
             thdif,cap,drycan,wetcan,cst,                     &
             tranf,transum,dew,mavail,                        &
!--- soil fixed fields
             dqm,qmin,psis,bclh,                              &
             zsmain,zshalf,DTDZS,tbq,                         &
!--- constants
             xlvm,CP,rovcp,G0_P,cvw,stbolt,                   &
!--- output variables
             snweprint,snheiprint,rsm,                        &
             tso,soilt,soilt1,tsnav,qvg,qsg,qcg,              &
             smelt,snoh,snflx,s,ilnb,x)

!************************************************************************
!--- RECALCULATION OF DEW USING NEW VALUE OF QSG OR TRANSP IF NO DEW
         DEW=zero
         ETT1=zero
         PP=PATM*rhowater
         EPOT = -FQ*(QVATM-QSG)
       IF(EPOT.GT.zero) THEN
! Evaporation
          DO K=1,NROOT
            TRANSP(K)=vegfrac*RAS*FQ*(QVATM-QSG)              &
                     *tranf(K)*DRYCAN/zshalf(NROOT+1)
            ETT1=ETT1-TRANSP(K)
          ENDDO
          DO k=nroot+1,nzs
            transp(k)=zero
          enddo

        ELSE
! Sublimation
          DEW=-EPOT
          DO K=1,NZS
            TRANSP(K)=zero
          ENDDO
        ETT1=zero
        ENDIF

!-- recalculating of frozen water in soil
         DO K=1,NZS
         tln=log(tso(k)/tfrz)
         if(tln.lt.zero) then
           soiliqw(k)=(dqm+qmin)*(XLMELT*                    &
           (tso(k)-tfrz)/tso(k)/grav/psis)                   &
          **(-one/bclh)-qmin
           soiliqw(k)=max(zero,soiliqw(k))
           soiliqw(k)=min(soiliqw(k),soilmois(k))
           soilice(k)=(soilmois(k)-soiliqw(k))/riw
!---- melting and freezing is balanced, soil ice cannot increase
       if(keepfr(k).eq.one) then
           soilice(k)=min(soilice(k),smfrkeep(k))
           soiliqw(k)=max(zero,soilmois(k)-soilice(k)*riw)
       endif

         else
           soilice(k)=zero
           soiliqw(k)=soilmois(k)
         endif
         ENDDO

!*************************************************************************
!--- TQCAN FOR SOLUTION OF MOISTURE BALANCE (Smirnova et al. 1996, EQ.22,28)
!    AND TSO,ETA PROFILES
!*************************************************************************
                CALL SOILMOIST (debug_print,xlat,xlon,testptlat,testptlon,&
!-- input
               delt,nzs,nddzs,DTDZS,DTDZS2,RIW,                    &
               zsmain,zshalf,diffu,hydro,                          &
               QSG,QVG,QCG,QCATM,QVATM,-INFWATER,                  &
               QKMS,TRANSP,zero,                                   &
               zero,SMELT,soilice,vegfrac,                         &
               snowfrac,one,                                       &
!-- soil properties
               DQM,QMIN,REF,KSAT,RAS,INFMAX,                       &
!-- output
               SOILMOIS,SOILIQW,MAVAIL,RUNOFF1,                    &
               RUNOFF2,infiltrp) 
 
!        endif

!-- Restore land-use parameters if all snow is melted
         IF(SNHEI.EQ.zero)  then
          tsnav=soilt-tfrz
         ENDIF

! 21apr2009
! SNOM [mm] goes into the passed-in ACSNOM variable in the grid derived type
        SNOM=SNOM+SMELT*DELT*rhowater
!
!--- KEEPFR is 1 when the temperature and moisture in soil
!--- are both increasing. In this case soil ice should not
!--- be increasing according to the freezing curve.
!--- Some part of ice is melted, but additional water is
!--- getting frozen. Thus, only structure of frozen soil is
!--- changed, and phase changes are not affecting the heat
!--- transfer. This situation may happen when it rains on the
!--- frozen soil.

        do k=1,nzs
       if (soilice(k).gt.zero) then
          if(tso(k).gt.told(k).and.soilmois(k).gt.smold(k)) then
              keepfr(k)=one
          else
              keepfr(k)=zero
          endif
       endif
        enddo
!--- THE DIAGNOSTICS OF SURFACE FLUXES

        T3      = STBOLT*SOILTold*SOILTold*SOILTold
        UPFLUX  = T3 *0.5_kind_phys*(SOILTold+SOILT)
        XINET   = EMISS*(GLW-UPFLUX)   
        HFX=-TKMS*CP*RHO*(TABS-SOILT)                        &
               *(P1000mb*0.00001_kind_phys/Patm)**ROVCP
    IF (debug_print ) THEN
      print *,'potential temp HFX',hfx
    ENDIF
        HFT=-TKMS*CP*RHO*(TABS-SOILT) 
    IF (debug_print ) THEN
      print *,'abs temp HFX',hft
    ENDIF
        Q1 = - FQ*RAS* (QVATM - QSG)
        CMC2MS= zero
        IF (Q1.LT.zero) THEN
! ---  condensation
        EDIR1=zero
        EC1=zero
        ETT1=zero
! ---  condensation
     if(myj) then
!-- moisture flux for coupling with MYJ PBL
          EETA=-QKMS*RAS*(QVATM/(1.+QVATM) - QSG/(1.+QSG))*rhowater
          CST= CST-EETA*DELT*vegfrac
    IF (debug_print ) THEN
      print *,'MYJ EETA cond', EETA
    ENDIF
     else ! myj
!-- actual moisture flux from RUC LSM
          DEW=QKMS*(QVATM-QSG)
          EETA= - RHO*DEW
          CST=CST+DELT*DEW*RAS * vegfrac
    IF (debug_print ) THEN
      print *,'RUC LSM EETA cond',EETA
    ENDIF
     endif ! myj
          QFX= XLVm*EETA
          EETA= - RHO*DEW
        ELSE
! ---  evaporation
        EDIR1 = Q1*UMVEG *BETA
        CMC2MS=CST/DELT*RAS
        EC1 = Q1 * WETCAN * vegfrac

        CST=max(zero,CST-EC1 * DELT)

    IF (debug_print ) THEN
     print*,'Q1,umveg,beta',Q1,umveg,beta
     print *,'wetcan,vegfrac',wetcan,vegfrac
     print *,'EC1,CMC2MS',EC1,CMC2MS
    ENDIF

     if(myj) then
!-- moisture flux for coupling with MYJ PBL
        EETA=-(QKMS*RAS*(QVATM/(one+QVATM) - QSG/(one+QSG))*rhowater)*BETA
    IF (debug_print ) THEN
      print *,'MYJ EETA', EETA*XLVm,EETA
    ENDIF
     else ! myj
! to convert from m s-1 to kg m-2 s-1: *rho water=1.e3************
!-- actual moisture flux from RUC LSM
        EETA = (EDIR1 + EC1 + ETT1)*rhowater
    IF (debug_print ) THEN
      print *,'RUC LSM EETA',EETA*XLVm,EETA
    ENDIF
     endif ! myj
        QFX= XLVm * EETA
        EETA = (EDIR1 + EC1 + ETT1)*rhowater
       ENDIF
        S=SNFLX
        sublim=Q1*rhowater !kg m-2 s-1
! Energy budget
        FLTOT=RNET-HFT-XLVm*EETA-S-SNOH-x
    IF (debug_print ) THEN
       print *,'SNOWSOIL - FLTOT,RNET,HFT,QFX,S,SNOH,X=',FLTOT,RNET,HFT,XLVm*EETA,s,SNOH,X
       print *,'edir1,ec1,ett1,mavail,qkms,qvatm,qvg,qsg,vegfrac,beta',&
                edir1,ec1,ett1,mavail,qkms,qvatm,qvg,qsg,vegfrac,beta
    ENDIF

 222     CONTINUE

 1123    FORMAT(I5,8F12.3)
 1133    FORMAT(I7,8E12.4)
  123   format(i6,f6.2,7f8.1)
 122    FORMAT(1X,2I3,6F8.1,F8.3,F8.2)

!-------------------------------------------------------------------
   END SUBROUTINE SNOWSOIL
!-------------------------------------------------------------------

!>\ingroup lsm_ruc_group
!! This subroutine is called for sea ice with accumulated snow on
!! its surface. It solves energy budget on the snow interface with 
!! atmosphere and snow interface with ice. It calculates skin 
!! temperature, snow and ice temperatures, snow depth and snow melt.
           SUBROUTINE SNOWSEAICE( debug_print,xlat,xlon,        &
            i,j,isoil,delt,ktau,conflx,nzs,nddzs,               &
            isncond_opt,isncovr_opt,                            &
            meltfactor,rhonewsn,SNHEI_CRIT,                     &  ! new
            ILAND,PRCPMS,RAINF,NEWSNOW,snhei,SNWE,snowfrac,     &
            RHOSN,PATM,QVATM,QCATM,                             &
            GLW,GSW,EMISS,RNET,                                 &
            QKMS,TKMS,RHO,myj,                                  &
            ALB,ZNT,                                            & !--- sea ice parameters
            tice,rhosice,capice,thdifice,                       &
            zsmain,zshalf,DTDZS,DTDZS2,tbq,                     &
            xlv,CP,rovcp,cw,stbolt,tabs,                        & !--- constants
            ilnb,snweprint,snheiprint,rsm,tso,                  & !--- output variables
            dew,soilt,soilt1,tsnav,qvg,qsg,qcg,                 &
            SMELT,SNOH,SNFLX,SNOM,eeta,                         &
            qfx,hfx,s,sublim,prcpl,fltot                        &
                                                                )
!***************************************************************
!   Solving energy budget for snow on sea ice and heat diffusion 
!   eqns. in snow and sea ice
!***************************************************************


        IMPLICIT NONE
!-------------------------------------------------------------------
!--- input variables

   LOGICAL,  INTENT(IN   )   ::  debug_print
   INTEGER,  INTENT(IN   )   ::  ktau,nzs     ,                  &
                                 nddzs                         !nddzs=2*(nzs-2)
   INTEGER,  INTENT(IN   )   ::  i,j,isoil,isncond_opt,isncovr_opt

   real (kind_phys),     INTENT(IN   )   ::  DELT,CONFLX,PRCPMS, &
                                 RAINF,NEWSNOW,RHONEWSN,         &
                                 meltfactor,snhei_crit,xlat,xlon
   real                      ::  rhonewcsn

   LOGICAL,  INTENT(IN   )   ::  myj
!--- 3-D Atmospheric variables
   real (kind_phys),                                             &
            INTENT(IN   )    ::                            PATM, &
                                                          QVATM, &
                                                          QCATM
!--- 2-D variables
   real (kind_phys)                                            , &
            INTENT(IN   )    ::                             GLW, &
                                                            GSW, &
                                                            RHO, &
                                                           QKMS, &
                                                           TKMS

!--- sea ice properties
   real (kind_phys),     DIMENSION(1:NZS)                      , &
            INTENT(IN   )    ::                                  &
                                                           tice, &
                                                        rhosice, &
                                                         capice, &
                                                       thdifice

   real (kind_phys),     INTENT(IN   )   ::                      &
                                                             CW, &
                                                            XLV

   real (kind_phys),     DIMENSION(1:NZS), INTENT(IN) :: ZSMAIN, &
                                                         ZSHALF, &
                                                         DTDZS2

   real (kind_phys),     DIMENSION(1:NDDZS), INTENT(IN)  :: DTDZS

   real (kind_phys),     DIMENSION(1:5001), INTENT(IN)  ::  TBQ

!--- input/output variables
!-------- 3-d soil moisture and temperature
   real (kind_phys),     DIMENSION(  1:nzs )                   , &
             INTENT(INOUT)   ::                             TSO

   INTEGER,  INTENT(INOUT)    ::                           ILAND


!-------- 2-d variables
   real (kind_phys)                                            , &
             INTENT(INOUT)   ::                             DEW, &
                                                           EETA, &
                                                          RHOSN, &
                                                         SUBLIM, &
                                                          PRCPL, &
                                                            ALB, &
                                                          EMISS, &
                                                            ZNT, &
                                                            QVG, &
                                                            QSG, &
                                                            QCG, &
                                                            QFX, &
                                                            HFX, &
                                                              S, &
                                                           SNWE, &
                                                          SNHEI, &
                                                          SMELT, &
                                                           SNOM, &
                                                           SNOH, &
                                                          SNFLX, &
                                                          SOILT, &
                                                         SOILT1, &
                                                       SNOWFRAC, &
                                                          TSNAV

   INTEGER, INTENT(INOUT)    ::                            ILNB

   real (kind_phys),     INTENT(OUT)                     :: RSM, &
                                                      SNWEPRINT, &
                                                     SNHEIPRINT
!--- Local variables


   INTEGER ::  nzs1,nzs2,k,k1,kn,kk
   real (kind_phys)    ::  x,x1,x2,dzstop,ft,tn,denom

   real (kind_phys)    ::  SNTH, NEWSN                         , &
               TABS, T3, UPFLUX, XINET                         , &
               BETA, SNWEPR,EPDT,PP
   real (kind_phys)    ::  CP,rovcp,G0,LV,xlvm,STBOLT,xlmelt   , &
               epot,fltot,fq,hft,q1,ras,ci,cvw                 , &
               RIW,DELTSN,H

   real (kind_phys)    ::  rhocsn,thdifsn,                       &
               xsn,ddzsn,x1sn,d1sn,d2sn,d9sn,r22sn

   real (kind_phys)    ::  cotsn,rhtsn,xsn1,ddzsn1,x1sn1,ftsnow,denomsn
   real (kind_phys)    ::  fso,fsn,                              &
               FKT,D1,D2,D9,D10,DID,R211,R21,R22,R6,R7,D11,      &
               FKQ,R210,AA,BB,QS1,TS1,TQ2,TX2,                   &
               TDENOM,AA1,RHCS,H1,TSOB, SNPRIM,                  &
               SNODIF,SOH,TNOLD,QGOLD,SNOHGNEW
   real (kind_phys),     DIMENSION(1:NZS)  ::  cotso,rhtso

   real (kind_phys)                   :: RNET,rsmfrac,soiltfrac,hsn,icemelt,rr
   integer                ::      nmelt

   real (kind_phys)                   :: keff, fact

!-----------------------------------------------------------------
        XLMELT=con_hfus
!-- heat of sublimation of water vapor
        XLVm=XLV+XLMELT

        !-- options for snow conductivity:
        !-- 1 - constant
        !-- opt 2 -  Sturm et al., 1997
        keff = 0.265_kind_phys

!--- SNOW flag -- ISICE
!--- DELTSN - is the threshold for splitting the snow layer into 2 layers.
!--- With snow density 400 kg/m^3, this threshold is equal to 7.5 cm,
!--- equivalent to 0.03 m SNWE. For other snow densities the threshold is
!--- computed using SNWE=0.03 m and current snow density.
!--- SNTH - the threshold below which the snow layer is combined with
!--- the top sea ice layer. SNTH is computed using snwe=0.016 m, and
!--- equals 4 cm for snow density 400 kg/m^3.

           DELTSN=0.05_kind_phys*rhowater/rhosn
           snth=0.01_kind_phys*rhowater/rhosn

! For 2-layer snow model when the snow depth is marginlly higher than DELTSN,
! reset DELTSN to half of snow depth.
        IF(SNHEI.GE.DELTSN+SNTH) THEN
          if(snhei-deltsn-snth.lt.snth) deltsn=0.5_kind_phys*(snhei-snth)
    IF (debug_print ) THEN
        print *,'DELTSN ICE is changed,deltsn,snhei,snth', &
                                  i,j, deltsn,snhei,snth
    ENDIF
        ENDIF

        CI=RHOICE*sheatice
        RAS=RHO*1.E-3_kind_dbl_prec ! rho/rhowater
        RIW=rhoice*1.e-3_kind_dbl_prec ! rhoice/rhowater
        RSM=zero

        XLMELT=con_hfus
        RHOCSN=sheatsn * RHOSN
!18apr08 - add rhonewcsn
        RHOnewCSN=sheatsn * RHOnewSN

      if(isncond_opt == 1) then
        !-- old version thdifsn = 0.265/RHOCSN
        THDIFSN = 0.265_kind_phys/RHOCSN
      else
      !-- 07Jun19 - thermal conductivity (K_eff) from Sturm et al.(1997)
      !-- keff = 10. ** (2.650 * RHOSN*1.e-3 - 1.652)
         fact = one
         if(rhosn < 156._kind_phys .or. (newsnow > zero .and. rhonewsn < 156._kind_phys)) then
           keff = 0.023_kind_phys + 0.234_kind_phys * rhosn * 1.e-3_kind_phys
         else
           keff = 0.138_kind_phys - 1.01_kind_phys * rhosn*1.e-3_kind_phys + 3.233_kind_phys * rhosn**2 * 1.e-6_kind_phys
         endif

         if(newsnow <= zero .and. snhei > one .and. rhosn > 250._kind_phys) then
         !-- some areas with large snow depth have unrealistically 
         !-- low snow density (in the Rockie's with snow depth > 1 m). 
         !-- Based on Sturm et al. keff=0.452 typical for hard snow slabs
         !-- with rhosn=488 kg/m^3. Thdifsn = 0.452/(2090*488)=4.431718e-7
         !-- In future a better compaction scheme is needed for these areas.
           thdifsn = 4.431718e-7_kind_phys 
         else
           thdifsn = keff/rhocsn * fact
         endif
      endif

        RAS=RHO*1.E-3_kind_phys

        SOILTFRAC=SOILT

        SMELT=zero
        SOH=zero
        SNODIF=zero
        SNOH=zero
        SNOHGNEW=zero
        RSM=zero
        RSMFRAC=zero
        fsn=one
        fso=zero
        cvw=cw

          NZS1=NZS-1
          NZS2=NZS-2

        QGOLD=QSG
        TNOLD=SOILT
        DZSTOP=one/(ZSMAIN(2)-ZSMAIN(1))

        snweprint=zero
        snheiprint=zero
        prcpl=prcpms

!*** DELTSN is the depth of the top layer of snow where
!*** there is a temperature gradient, the rest of the snow layer
!*** is considered to have constant temperature


        H=one
        SMELT=zero

        FQ=QKMS
        SNHEI=SNWE*rhowater/RHOSN
        SNWEPR=SNWE

!  check if all snow can evaporate during DT
         BETA=one
         EPOT = -FQ*(QVATM-QSG)
         EPDT = EPOT * RAS *DELT
         IF(EPDT.GT.zero .and. SNWEPR.LE.EPDT) THEN
            BETA=SNWEPR/max(1.e-8_kind_phys,EPDT)
            SNWE=zero
         ENDIF

!******************************************************************************
!       COEFFICIENTS FOR THOMAS ALGORITHM FOR TSO
!******************************************************************************

        cotso(1)=zero
        rhtso(1)=TSO(NZS)
        DO 33 K=1,NZS2
          KN=NZS-K
          K1=2*KN-3
          X1=DTDZS(K1)*THDIFICE(KN-1)
          X2=DTDZS(K1+1)*THDIFICE(KN)
          FT=TSO(KN)+X1*(TSO(KN-1)-TSO(KN))                           &
             -X2*(TSO(KN)-TSO(KN+1))
          DENOM=1.+X1+X2-X2*cotso(K)
          cotso(K+1)=X1/DENOM
          rhtso(K+1)=(FT+X2*rhtso(K))/DENOM
   33  CONTINUE
!--- THE NZS element in COTSO and RHTSO will be for snow
!--- There will be 2 layers in snow if it is deeper than DELTSN+SNTH
       IF(SNHEI.GE.SNTH) then
        if(snhei.le.DELTSN+SNTH) then
!-- 1-layer snow model
         ilnb=1
         snprim=max(snth,snhei)
         soilt1=tso(1)
         tsob=tso(1)
         XSN = DELT/2._kind_phys/(zshalf(2)+0.5_kind_phys*SNPRIM)
         DDZSN = XSN / SNPRIM
         X1SN = DDZSN * thdifsn
         X2 = DTDZS(1)*THDIFICE(1)
         FT = TSO(1)+X1SN*(SOILT-TSO(1))                              &
              -X2*(TSO(1)-TSO(2))
         DENOM = one + X1SN + X2 -X2*cotso(NZS1)
         cotso(NZS)=X1SN/DENOM
         rhtso(NZS)=(FT+X2*rhtso(NZS1))/DENOM
         cotsn=cotso(NZS)
         rhtsn=rhtso(NZS)
!*** Average temperature of snow pack (C)
         tsnav=0.5_kind_phys*(soilt+tso(1))                            &
                     -tfrz

        else
!-- 2 layers in snow, SOILT1 is temperasture at DELTSN depth
         ilnb=2
         snprim=deltsn
         tsob=soilt1
         XSN = DELT/2._kind_phys/(0.5_kind_phys*SNHEI)
         XSN1= DELT/2._kind_phys/(zshalf(2)+0.5_kind_phys*(SNHEI-DELTSN))
         DDZSN = XSN / DELTSN
         DDZSN1 = XSN1 / (SNHEI-DELTSN)
         X1SN = DDZSN * thdifsn
         X1SN1 = DDZSN1 * thdifsn
         X2 = DTDZS(1)*THDIFICE(1)
         FT = TSO(1)+X1SN1*(SOILT1-TSO(1))                            &
              -X2*(TSO(1)-TSO(2))
         DENOM = one + X1SN1 + X2 - X2*cotso(NZS1)
         cotso(nzs)=x1sn1/denom
         rhtso(nzs)=(ft+x2*rhtso(nzs1))/denom
         ftsnow = soilt1+x1sn*(soilt-soilt1)                          &
               -x1sn1*(soilt1-tso(1))
         denomsn = 1. + X1SN + X1SN1 - X1SN1*cotso(NZS)
         cotsn=x1sn/denomsn
         rhtsn=(ftsnow+X1SN1*rhtso(NZS))/denomsn
!*** Average temperature of snow pack (C)
         tsnav=0.5_kind_phys/snhei*((soilt+soilt1)*deltsn             &
                     +(soilt1+tso(1))*(SNHEI-DELTSN))                 &
                     -tfrz
        endif
       ENDIF

       IF(SNHEI.LT.SNTH.AND.SNHEI.GT.zero) then
!--- snow is too thin to be treated separately, therefore it
!--- is combined with the first sea ice layer.
         snprim=SNHEI+zsmain(2)
         fsn=SNHEI/snprim
         fso=one-fsn
         soilt1=tso(1)
         tsob=tso(2)
         XSN = DELT/2._kind_phys/((zshalf(3)-zsmain(2))+0.5_kind_phys*snprim)
         DDZSN = XSN /snprim
         X1SN = DDZSN * (fsn*thdifsn+fso*thdifice(1))
         X2=DTDZS(2)*THDIFICE(2)
         FT=TSO(2)+X1SN*(SOILT-TSO(2))-                              &
                       X2*(TSO(2)-TSO(3))
         denom = one + x1sn + x2 - x2*cotso(nzs-2)
         cotso(nzs1) = x1sn/denom
         rhtso(nzs1)=(FT+X2*rhtso(NZS-2))/denom
         tsnav=0.5_kind_phys*(soilt+tso(1))                          &
                     -tfrz
         cotso(nzs)=cotso(NZS1)
         rhtso(nzs)=rhtso(nzs1)
         cotsn=cotso(NZS)
         rhtsn=rhtso(NZS)
       ENDIF

!************************************************************************
!--- THE HEAT BALANCE EQUATION 
!18apr08 nmelt is the flag for melting, and SNOH is heat of snow phase changes
       nmelt=0
       SNOH=zero

        EPOT=-QKMS*(QVATM-QSG)
        RHCS=CAPICE(1)
        H=one
        FKT=TKMS
        D1=cotso(NZS1)
        D2=rhtso(NZS1)
        TN=SOILT
        D9=THDIFICE(1)*RHCS*dzstop
        D10=TKMS*CP*RHO
        R211=.5_kind_phys*CONFLX/DELT
        R21=R211*CP*RHO
        R22=.5_kind_phys/(THDIFICE(1)*DELT*dzstop**2)
        R6=EMISS *STBOLT*.5_kind_phys*TN**4
        R7=R6/TN
        D11=RNET+R6

      IF(SNHEI.GE.SNTH) THEN 
        if(snhei.le.DELTSN+SNTH) then
!--- 1-layer snow
          D1SN = cotso(NZS)
          D2SN = rhtso(NZS)
        else
!--- 2-layer snow
          D1SN = cotsn
          D2SN = rhtsn
        endif
        D9SN= THDIFSN*RHOCSN / SNPRIM
        R22SN = SNPRIM*SNPRIM*0.5_kind_phys/(THDIFSN*DELT)
      ENDIF

       IF(SNHEI.LT.SNTH.AND.SNHEI.GT.zero) then
!--- thin snow is combined with sea ice
         D1SN = D1
         D2SN = D2
         D9SN = (fsn*THDIFSN*RHOCSN+fso*THDIFICE(1)*RHCS)/           &
                 snprim
         R22SN = snprim*snprim*0.5_kind_phys                         &
                 /((fsn*THDIFSN+fso*THDIFICE(1))*delt)
      ENDIF

      IF(SNHEI.eq.zero)then
!--- all snow is sublimated
        D9SN = D9
        R22SN = R22
        D1SN = D1
        D2SN = D2
      ENDIF


!---- TDENOM for snow
        TDENOM = D9SN*(one-D1SN +R22SN)+D10+R21+R7                   &
              +RAINF*CVW*PRCPMS                                      &
              +RHOnewCSN*NEWSNOW/DELT

        FKQ=QKMS*RHO
        R210=R211*RHO
        AA=XLVM*(BETA*FKQ+R210)/TDENOM
        BB=(D10*TABS+R21*TN+XLVM*(QVATM*                             &
        (BETA*FKQ)                                                   &
        +R210*QVG)+D11+D9SN*(D2SN+R22SN*TN)                          &
        +RAINF*CVW*PRCPMS*max(tfrz,TABS)                             &
        + RHOnewCSN*NEWSNOW/DELT*min(tfrz,TABS)                      &
         )/TDENOM
        AA1=AA
        PP=PATM*1.E3_kind_phys
        AA1=AA1/PP
!18apr08  - the iteration start point
 212    continue
        BB=BB-SNOH/TDENOM
    IF (debug_print ) THEN
        print *,'VILKA-SNOW on SEAICE'
        print *,'tn,aa1,bb,pp,fkq,r210',                             &
                 tn,aa1,bb,pp,fkq,r210
        print *,'TABS,QVATM,TN,QVG=',TABS,QVATM,TN,QVG
    ENDIF

        CALL VILKA(TN,AA1,BB,PP,QS1,TS1,TBQ,KTAU,i,j,iland,isoil,xlat,xlon)
!--- it is saturation over snow
        QVG=QS1
        QSG=QS1
        QCG=zero

!--- SOILT - skin temperature of snow on ice
        SOILT=TS1
     if(nmelt==1 .and. snowfrac==one) then
       soilt = min(tfrz,soilt)
     endif

    IF (debug_print ) THEN
        print *,' AFTER VILKA-SNOW on SEAICE'
        print *,' TS1,QS1: ', ts1,qs1
    ENDIF
! Solution for temperature at 7.5 cm depth and snow-seaice interface
       IF(SNHEI.GE.SNTH) THEN
        if(snhei.gt.DELTSN+SNTH) then
!-- 2-layer snow model
          SOILT1=min(tfrz,rhtsn+cotsn*SOILT)
          TSO(1)=min(con_tice,(rhtso(NZS)+cotso(NZS)*SOILT1))
          tsob=soilt1
        else
!-- 1 layer in snow
          TSO(1)=min(con_tice,(rhtso(NZS)+cotso(NZS)*SOILT))
          SOILT1=TSO(1)
          tsob=tso(1)
        endif
       ELSEIF  (SNHEI > zero .and. SNHEI < SNTH) THEN
! blended
         TSO(2)=min(con_tice,(rhtso(NZS1)+cotso(NZS1)*SOILT))
         tso(1)=min(con_tice,(tso(2)+(soilt-tso(2))*fso))
         SOILT1=TSO(1)
         tsob=TSO(2)
       ELSE
! snow is melted
         TSO(1)=min(con_tice,SOILT)
         SOILT1=min(con_tice,SOILT)
         tsob=tso(1)
       ENDIF
!---- Final solution for TSO in sea ice
       IF (SNHEI > zero .and. SNHEI < SNTH) THEN
! blended or snow is melted
          DO K=3,NZS
            KK=NZS-K+1
            TSO(K)=min(con_tice,rhtso(KK)+cotso(KK)*TSO(K-1))
          END DO
       ELSE
          DO K=2,NZS
            KK=NZS-K+1
            TSO(K)=min(con_tice,rhtso(KK)+cotso(KK)*TSO(K-1))
          END DO
       ENDIF
!--- For thin snow layer combined with the top soil layer
!--- TSO(i,j,1) is computed by linear interpolation between SOILT
!--- and TSO(i,j,2)
!       if(SNHEI.LT.SNTH.AND.SNHEI.GT.0.)then
!          tso(1)=min(271.4,tso(2)+(soilt-tso(2))*fso)
!          soilt1=tso(1)
!          tsob = tso(2)
!       endif

      if(nmelt.eq.1) go to 220

!--- IF SOILT > tfrz F then melting of snow can happen
!    if all snow can evaporate, then there is nothing to melt
   IF(SOILT>tfrz .AND. BETA==one .AND. SNHEI>zero) THEN
!
        nmelt = 1
        soiltfrac=snowfrac*tfrz+(1.-snowfrac)*min(con_tice,SOILT)

        QSG= QSN(soiltfrac,TBQ)/PP
        T3      = STBOLT*TNold*TNold*TNold
        UPFLUX  = T3 * 0.5_kind_phys*(TNold+SOILTfrac)
        XINET   = EMISS*(GLW-UPFLUX)
         EPOT = -QKMS*(QVATM-QSG)
         Q1=EPOT*RAS

        IF (Q1.LE.zero) THEN
! ---  condensation
          DEW=-EPOT

        QFX= XLVM*RHO*DEW
        EETA=QFX/XLVM
       ELSE
! ---  evaporation
        EETA = Q1 * BETA * rhowater
! to convert from kg m-2 s-1 to m s-1: 1/rho water=1.e-3************
        QFX= - XLVM * EETA
       ENDIF

         HFX=D10*(TABS-soiltfrac)

       IF(SNHEI.GE.SNTH)then
         SOH=thdifsn*RHOCSN*(soiltfrac-TSOB)/SNPRIM
         SNFLX=SOH
       ELSE
         SOH=(fsn*thdifsn*rhocsn+fso*thdifice(1)*rhcs)*                &
              (soiltfrac-TSOB)/snprim
         SNFLX=SOH
       ENDIF
         X= (R21+D9SN*R22SN)*(soiltfrac-TNOLD) +                        &
            XLVM*R210*(QSG-QGOLD)
!-- SNOH is energy flux of snow phase change
        SNOH=RNET+QFX +HFX                                              &
                  +RHOnewCSN*NEWSNOW/DELT*(min(tfrz,TABS)-soiltfrac)  &
                  -SOH-X+RAINF*CVW*PRCPMS*                              &
                  (max(tfrz,TABS)-soiltfrac)

    IF (debug_print ) THEN
     print *,'SNOWSEAICE melt I,J,SNOH,RNET,QFX,HFX,SOH,X',i,j,SNOH,RNET,QFX,HFX,SOH,X
     print *,'RHOnewCSN*NEWSNOW/DELT*(min(tfrz,TABS)-soiltfrac)',     &
              RHOnewCSN*NEWSNOW/DELT*(min(tfrz,TABS)-soiltfrac)
     print *,'RAINF*CVW*PRCPMS*(max(tfrz,TABS)-soiltfrac)',           &
              RAINF*CVW*PRCPMS*(max(tfrz,TABS)-soiltfrac)
    ENDIF
        SNOH=AMAX1(zero,SNOH)
!-- SMELT is speed of melting in M/S
        SMELT= SNOH /XLMELT*1.E-3_kind_phys
        SMELT=AMIN1(SMELT,SNWEPR/DELT-BETA*EPOT*RAS)
        SMELT=AMAX1(zero,SMELT)

    IF (debug_print ) THEN
       print *,'1-SMELT i,j',smelt,i,j
    ENDIF
!18apr08 - Egglston limit
       SMELT= amin1 (smelt,delt/60._kind_phys* 5.6E-8_kind_phys*meltfactor*max(one,(soilt-tfrz))) ! SnowMIP
    IF (debug_print ) THEN
       print *,'2-SMELT i,j',smelt,i,j
    ENDIF

! rr - potential melting
        rr=SNWEPR/delt-BETA*EPOT*RAS
        SMELT=min(SMELT,rr)
    IF (debug_print ) THEN
      print *,'3- SMELT i,j,smelt,rr',i,j,smelt,rr
    ENDIF
        SNOHGNEW=SMELT*XLMELT*1.E3
        SNODIF=AMAX1(zero,(SNOH-SNOHGNEW))

        SNOH=SNOHGNEW

    IF (debug_print ) THEN
       print*,'soiltfrac,soilt,SNOHGNEW,SNODIF=', &
            i,j,soiltfrac,soilt,snohgnew,snodif
       print *,'SNOH,SNODIF',SNOH,SNODIF
    ENDIF

!*** From Koren et al. (1999) 13% of snow melt stays in the snow pack
        rsmfrac=min(0.18_kind_phys,(max(0.08_kind_phys,snwepr/0.10_kind_phys*0.13_kind_phys)))
       if(snhei > 0.01_kind_phys) then
        rsm=rsmfrac*smelt*delt
       else
! do not keep melted water if snow depth is less that 1 cm
        rsm=zero
       endif
!18apr08 rsm is part of melted water that stays in snow as liquid
        SMELT=AMAX1(zero,SMELT-rsm/delt)
    IF (debug_print ) THEN
       print *,'4-SMELT i,j,smelt,rsm,snwepr,rsmfrac', &
                    i,j,smelt,rsm,snwepr,rsmfrac
    ENDIF

!-- update liquid equivalent of snow depth
!-- for evaporation and snow melt
        SNWE = AMAX1(zero,(SNWEPR-                                    &
                    (SMELT+BETA*EPOT*RAS)*DELT                        &
                                         ) )
        soilt=soiltfrac
!--- If there is no snow melting then just evaporation
!--- or condensation changes SNWE
      ELSE
       if(snhei > zero.and. beta == one) then
               EPOT=-QKMS*(QVATM-QSG)
               SNWE = AMAX1(zero,(SNWEPR-                             &
                    BETA*EPOT*RAS*DELT))
       else
         snwe = zero
       endif

      ENDIF

! no iteration for snow on sea ice, because it will produce
! skin temperature higher than it is possible with snow on sea ice
!      if(nmelt.eq.1) goto 212  ! second iteration
 220  continue

       if(smelt > zero .and. rsm > zero) then
        if(snwe.le.rsm) then
    IF (debug_print ) THEN
     print *,'SEAICE SNWE<RSM snwe,rsm,smelt*delt,epot*ras*delt,beta', &
                              snwe,rsm,smelt*delt,epot*ras*delt,beta
    ENDIF
        else
!*** Update snow density on effect of snow melt, melted
!*** from the top of the snow. 13% of melted water
!*** remains in the pack and changes its density.
!*** Eq. 9 (with my correction) in Koren et al. (1999)

         xsn=(rhosn*(snwe-rsm)+rhowater*rsm)/                            &
             snwe
         rhosn=MIN(MAX(58.8_kind_phys,XSN),500._kind_phys)

        RHOCSN=sheatsn* RHOSN
        if(isncond_opt == 1) then
        !-- old version thdifsn = 0.265/RHOCSN
          THDIFSN = 0.265_kind_phys/RHOCSN
        else
      !-- 07Jun19 - thermal conductivity (K_eff) from Sturm et al.(1997)
      !-- keff = 10. ** (2.650 * RHOSN*1.e-3 - 1.652)
         fact = one
         if(rhosn < 156._kind_phys .or. (newsn > zero .and. rhonewsn < 156._kind_phys)) then
           keff = 0.023_kind_phys + 0.234_kind_phys * rhosn * 1.e-3_kind_phys
         else
           keff = 0.138_kind_phys - 1.01_kind_phys * rhosn*1.e-3_kind_phys + 3.233_kind_phys * rhosn**2 * 1.e-6_kind_phys
         endif
        
         if(newsnow <= zero .and. snhei > one .and. rhosn > 250._kind_phys) then
         !-- some areas with large snow depth have unrealistically 
         !-- low snow density (in the Rockie's with snow depth > 1 m). 
         !-- Based on Sturm et al. keff=0.452 typical for hard snow slabs
         !-- with rhosn=488 kg/m^3. Thdifsn = 0.452/(2090*488)=4.431718e-7
         !-- In future a better compaction scheme is needed for these areas.
           thdifsn = 4.431718e-7_kind_phys
         else
           thdifsn = keff/rhocsn * fact
         endif
        endif

        endif
      endif

        snweprint=snwe
!--- if VEGFRAC.ne.0. then some snow stays on the canopy
!--- and should be added to SNWE for water conservation
!                +VEGFRAC*cst
        snheiprint=snweprint*rhowater / RHOSN

    IF (debug_print ) THEN
print *, 'snweprint : ',snweprint
print *, 'D9SN,SOILT,TSOB : ', D9SN,SOILT,TSOB
    ENDIF
      IF(SNHEI.GT.zero) THEN
        if(ilnb.gt.1) then
          tsnav=0.5_kind_phys/snhei*((soilt+soilt1)*deltsn           &
                    +(soilt1+tso(1))*(SNHEI-DELTSN))                 &
                       -tfrz
        else
          tsnav=0.5_kind_phys*(soilt+tso(1)) - tfrz
        endif
      ENDIF
!--- RECALCULATION OF DEW USING NEW VALUE OF QSG
         DEW=zero
         PP=PATM*1.E3_kind_phys
         QSG= QSN(SOILT,TBQ)/PP
         EPOT = -FQ*(QVATM-QSG)
       IF(EPOT.LT.zero) THEN
! Sublimation
          DEW=-EPOT
        ENDIF

        SNOM=SNOM+SMELT*DELT*rhowater

!--- THE DIAGNOSTICS OF SURFACE FLUXES

        T3      = STBOLT*TNold*TNold*TNold
        UPFLUX  = T3 *0.5_kind_phys*(SOILT+TNold)
        XINET   = EMISS*(GLW-UPFLUX)
        HFT=-TKMS*CP*RHO*(TABS-SOILT)
        HFX=-TKMS*CP*RHO*(TABS-SOILT)                        &
               *(P1000mb*0.00001_kind_phys/Patm)**ROVCP
        Q1 = - FQ*RAS* (QVATM - QSG)
        IF (Q1.LT.zero) THEN
! ---  condensation
      if(myj) then
!-- moisture flux for coupling with MYJ PBL
          EETA=-QKMS*RAS*(QVATM/(1.+QVATM) - QSG/(1.+QSG))*rhowater
      else ! myj
!-- actual moisture flux from RUC LSM
          DEW=QKMS*(QVATM-QSG)
          EETA= - RHO*DEW
      endif ! myj
          QFX= XLVm*EETA
          EETA= - RHO*DEW
          sublim = EETA
        ELSE
! ---  evaporation
      if(myj) then
!-- moisture flux for coupling with MYJ PBL
          EETA=-QKMS*RAS*BETA*(QVATM/(1.+QVATM) - QVG/(1.+QVG))*rhowater
      else ! myj
! to convert from m s-1 to kg m-2 s-1: *rho water=1.e3************
!-- actual moisture flux from RUC LSM
          EETA = Q1*BETA*rhowater
      endif ! myj
          QFX= XLVm * EETA
          EETA = Q1*BETA*rhowater
          sublim = EETA
        ENDIF

        icemelt=zero
      IF(SNHEI.GE.SNTH)then
         S=thdifsn*RHOCSN*(soilt-TSOB)/SNPRIM
         SNFLX=S
       ELSEIF(SNHEI.lt.SNTH.and.SNHEI.GT.zero) then
         S=(fsn*thdifsn*rhocsn+fso*thdifice(1)*rhcs)*                &
              (soilt-TSOB)/snprim
         SNFLX=S
    IF (debug_print ) THEN
      print *,'SNOW is thin, snflx',i,j,snflx
    ENDIF
       ELSE 
         SNFLX=D9SN*(SOILT-TSOB)
    IF (debug_print ) THEN
      print *,'SNOW is GONE, snflx',i,j,snflx
    ENDIF
       ENDIF

        SNHEI=SNWE *rhowater / RHOSN

    IF (debug_print ) THEN
       print *,'SNHEI,SNOH',i,j,SNHEI,SNOH
    ENDIF
!
         X= (R21+D9SN*R22SN)*(soilt-TNOLD) +              &
            XLVM*R210*(QSG-QGOLD)
    IF (debug_print ) THEN
     print *,'SNOWSEAICE storage ',i,j,x
     print *,'R21,D9sn,r22sn,soiltfrac,tnold,qsg,qgold,snprim', &
              R21,D9sn,r22sn,soiltfrac,tnold,qsg,qgold,snprim
    ENDIF
         X=X &
        -RHOnewCSN*NEWSNOW/DELT*(min(tfrz,TABS)-SOILT)        &
        -RAINF*CVW*PRCPMS*(max(tfrz,TABS)-SOILT)

! -- excess energy is spent on ice melt
        icemelt = RNET-HFT-XLVm*EETA-S-SNOH-X
    IF (debug_print ) THEN
        print *,'SNOWSEAICE icemelt=',icemelt
    ENDIF

        FLTOT=RNET-HFT-XLVm*EETA-S-SNOH-x-icemelt
    IF (debug_print ) THEN
       print *,'i,j,snhei,qsg,soilt,soilt1,tso,TABS,QVATM', &
                i,j,snhei,qsg,soilt,soilt1,tso,tabs,qvatm
       print *,'SNOWSEAICE - FLTOT,RNET,HFT,QFX,S,SNOH,icemelt,snodif,X,SOILT=' &
                      ,FLTOT,RNET,HFT,XLVm*EETA,s,SNOH,icemelt,snodif,X,SOILT
    ENDIF
!-- Restore sea-ice parameters if snow is less than threshold
         IF(SNHEI.EQ.zero)  then
          tsnav=soilt-tfrz
          emiss=0.98_kind_phys
          znt=0.011_kind_phys
          alb=0.55_kind_phys
         ENDIF

!------------------------------------------------------------------------
!------------------------------------------------------------------------
   END SUBROUTINE SNOWSEAICE
!------------------------------------------------------------------------

!>\ingroup lsm_ruc_group
!> This subroutine solves energy budget equation and heat diffusion
!! equation.
           SUBROUTINE SOILTEMP( debug_print,xlat,xlon,testptlat,testptlon,&
           i,j,iland,isoil,                                 & !--- input variables
           delt,ktau,conflx,nzs,nddzs,nroot,                &
           PRCPMS,RAINF,PATM,TABS,QVATM,QCATM,              &
           EMISS,RNET,                                      &
           QKMS,TKMS,PC,RHO,VEGFRAC,lai,                    &
           THDIF,CAP,DRYCAN,WETCAN,                         &
           TRANSUM,DEW,MAVAIL,soilres,alfa,                 &
           DQM,QMIN,BCLH,                                   & !---soil fixed fields
           ZSMAIN,ZSHALF,DTDZS,TBQ,                         &
           XLV,CP,G0_P,CVW,STBOLT,                          & !--- constants
           TSO,SOILT,QVG,QSG,QCG,X)                           !---output variables

!*************************************************************
!   Energy budget equation and heat diffusion eqn are 
!   solved here and
!
!     DELT - time step (s)
!     ktau - number of time step
!     CONFLX - depth of constant flux layer (m)
!     IME, JME, KME, NZS - dimensions of the domain 
!     NROOT - number of levels within the root zone
!     PRCPMS - precipitation rate in m/s
!     COTSO, RHTSO - coefficients for implicit solution of
!                     heat diffusion equation
!     THDIF - thermal diffusivity (m^2/s)
!     QSG,QVG,QCG - saturated mixing ratio, mixing ratio of
!                   water vapor and cloud at the ground
!                   surface, respectively (kg/kg)
!     PATM -  pressure [bar]
!     QC3D,QV3D - cloud and water vapor mixing ratio
!                   at the first atm. level (kg/kg)
!     EMISS,RNET - emissivity (0-1) of the ground surface and net
!                  radiation at the surface (W/m^2)
!     QKMS - exchange coefficient for water vapor in the
!              surface layer (m/s)
!     TKMS - exchange coefficient for heat in the surface
!              layer (m/s)
!     PC - plant coefficient (resistance)
!     RHO - density of atmosphere near surface (kg/m^3)
!     VEGFRAC - greeness fraction (0-1)
!     CAP - volumetric heat capacity (J/m^3/K)
!     DRYCAN - dry fraction of vegetated area where
!              transpiration may take place (0-1)
!     WETCAN - fraction of vegetated area covered by canopy
!              water (0-1)
!     TRANSUM - transpiration function integrated over the 
!               rooting zone (m)
!     DEW -  dew in kg/m^2s
!     MAVAIL - fraction of maximum soil moisture in the top
!               layer (0-1)
!     ZSMAIN - main levels in soil (m)
!     ZSHALF - middle of the soil layers (m)
!     DTDZS - dt/(2.*dzshalf*dzmain)
!     TBQ - table to define saturated mixing ration
!           of water vapor for given temperature and pressure
!     TSO - soil temperature (K)
!     SOILT - skin temperature (K)
!
!****************************************************************

        IMPLICIT NONE
!-----------------------------------------------------------------

!--- input variables

   LOGICAL,  INTENT(IN   )   ::  debug_print
   INTEGER,  INTENT(IN   )   ::  nroot,ktau,nzs                , &
                                 nddzs                         !nddzs=2*(nzs-2)
   INTEGER,  INTENT(IN   )   ::  i,j,iland,isoil
   real (kind_phys),     INTENT(IN   )   ::  DELT,CONFLX,PRCPMS, RAINF
   real (kind_phys),     INTENT(IN   )   ::  xlat, xlon, testptlat, testptlon
   real (kind_phys),     INTENT(INOUT)   ::  DRYCAN,WETCAN,TRANSUM
!--- 3-D Atmospheric variables
   real (kind_phys),                                             &
            INTENT(IN   )    ::                            PATM, &
                                                          QVATM, &
                                                          QCATM
!--- 2-D variables
   real (kind_phys)                                            , &
            INTENT(IN   )    ::                                  &
                                                          EMISS, &
                                                            RHO, &
                                                           RNET, &  
                                                             PC, &
                                                        VEGFRAC, &
                                                            LAI, &
                                                            DEW, & 
                                                           QKMS, &
                                                           TKMS

!--- soil properties
   real (kind_phys)                                            , &
            INTENT(IN   )    ::                                  &
                                                           BCLH, &
                                                            DQM, &
                                                           QMIN
   real (kind_phys)                                            , &
            INTENT(IN   )    ::                                  &
                                                   soilres,alfa


   real (kind_phys),     INTENT(IN   )   ::                  CP, &
                                                            CVW, &
                                                            XLV, &
                                                         STBOLT, &
                                                           TABS, &
                                                           G0_P


   real (kind_phys),     DIMENSION(1:NZS), INTENT(IN) :: ZSMAIN, &
                                                         ZSHALF, &
                                                          THDIF, &
                                                            CAP

   real (kind_phys),     DIMENSION(1:NDDZS), INTENT(IN)  :: DTDZS

   real (kind_phys),     DIMENSION(1:5001), INTENT(IN)  ::  TBQ


!--- input/output variables
!-------- 3-d soil moisture and temperature
   real (kind_phys),     DIMENSION( 1:nzs )                    , &
             INTENT(INOUT)   ::                             TSO

!-------- 2-d variables
   real (kind_phys)                                            , &
             INTENT(INOUT)   ::                                  &
                                                         MAVAIL, &
                                                            QVG, &
                                                            QSG, &
                                                            QCG, &
                                                          SOILT


!--- Local variables

   real (kind_phys)    ::  x,x1,x2,x4,dzstop,can,ft,sph        , &
               tn,trans,umveg,denom,fex

   real (kind_phys)    ::  FKT,D1,D2,D9,D10,DID,R211,R21,R22,R6,R7,D11, &
                      PI,H,FKQ,R210,AA,BB,PP,Q1,QS1,TS1,TQ2,TX2       , &
                      TDENOM

   real (kind_phys)    ::  C,CC,AA1,RHCS,H1, QGOLD

   real (kind_phys),     DIMENSION(1:NZS)  ::                   cotso,rhtso

   INTEGER ::  nzs1,nzs2,k,k1,kn,kk, iter


!-----------------------------------------------------------------

        iter=0

          NZS1=NZS-1
          NZS2=NZS-2
        dzstop=1./(ZSMAIN(2)-ZSMAIN(1))

        qgold=qvg

        do k=1,nzs
           cotso(k)=zero
           rhtso(k)=zero
        enddo
!******************************************************************************
!       COEFFICIENTS FOR THOMAS ALGORITHM FOR TSO
!******************************************************************************
        cotso(1)=zero
        rhtso(1)=TSO(NZS)
        DO 33 K=1,NZS2
          KN=NZS-K
          K1=2*KN-3
          X1=DTDZS(K1)*THDIF(KN-1)
          X2=DTDZS(K1+1)*THDIF(KN)
          FT=TSO(KN)+X1*(TSO(KN-1)-TSO(KN))                             &
             -X2*(TSO(KN)-TSO(KN+1))
          DENOM=1.+X1+X2-X2*cotso(K)
          cotso(K+1)=X1/DENOM
          rhtso(K+1)=(FT+X2*rhtso(K))/DENOM
   33  CONTINUE

!************************************************************************
!--- THE HEAT BALANCE EQUATION (Smirnova et al., 1996, EQ. 21,26)

        RHCS=CAP(1)

        H=MAVAIL

        TRANS=TRANSUM*DRYCAN/ZSHALF(NROOT+1)
        CAN=WETCAN+TRANS
        UMVEG=(1.-VEGFRAC) * soilres
 2111   continue
        FKT=TKMS
        D1=cotso(NZS1)
        D2=rhtso(NZS1)
        TN=SOILT
        D9=THDIF(1)*RHCS*dzstop
        D10=TKMS*CP*RHO
        R211=.5_kind_phys*CONFLX/DELT
        R21=R211*CP*RHO
        R22=.5_kind_phys/(THDIF(1)*DELT*dzstop**2)
        R6=EMISS *STBOLT*.5_kind_phys*TN**4
        R7=R6/TN
        D11=RNET+R6
        TDENOM=D9*(one-D1+R22)+D10+R21+R7                             &
              +RAINF*CVW*PRCPMS
        FKQ=QKMS*RHO
        R210=R211*RHO
        C=VEGFRAC*FKQ*CAN
        CC=C*XLV/TDENOM
        AA=XLV*(FKQ*UMVEG+R210)/TDENOM
        BB=(D10*TABS+R21*TN+XLV*(QVATM*                               &
        (FKQ*UMVEG+C)                                                 & 
        +R210*QVG)+D11+D9*(D2+R22*TN)                                 &
        +RAINF*CVW*PRCPMS*max(tfrz,TABS)                              &
         )/TDENOM
        AA1=AA+CC
        PP=PATM*1.E3_kind_phys
        AA1=AA1/PP
        CALL VILKA(TN,AA1,BB,PP,QS1,TS1,TBQ,KTAU,i,j,iland,isoil,xlat,xlon)
        TQ2=QVATM
        TX2=TQ2*(one-H)
        Q1=TX2+H*QS1
    IF (debug_print ) THEN
        print *,'VILKA1 - TS1,QS1,TQ2,H,TX2,Q1',TS1,QS1,TQ2,H,TX2,Q1
    ENDIF
        IF(Q1.LT.QS1) GOTO 100
!--- if no saturation - goto 100
!--- if saturation - goto 90
   90   QVG=QS1
        QSG=QS1
        TSO(1)=TS1
        QCG=max(zero,Q1-QS1)
    IF (debug_print ) THEN
        print *,'90 QVG,QSG,QCG,TSO(1)',QVG,QSG,QCG,TSO(1)
    ENDIF

        GOTO 200
  100   BB=BB-AA*TX2
        AA=(AA*H+CC)/PP

        CALL VILKA(TN,AA,BB,PP,QS1,TS1,TBQ,KTAU,i,j,iland,isoil,xlat,xlon)
        Q1=TX2+H*QS1
    IF (debug_print ) THEN
!     if(i.eq.279.and.j.eq.263) then
        print *,'VILKA2 - TS1,QS1,TQ2,H,TX2,Q1',TS1,QS1,TQ2,H,TX2,Q1
    ENDIF
        IF(Q1.GE.QS1) GOTO 90
        QSG=QS1
        QVG=Q1
!   if( QS1>QVATM .and. QVATM > QVG) then
    ! very dry soil 
    ! print *,'very dry soils mavail,qvg,qs1,qvatm,ts1',i,j,mavail,qvg,qs1,qvatm,ts1
    ! QVG = QVATM
!   endif
        TSO(1)=TS1
        QCG=zero
  200   CONTINUE
    IF (debug_print ) THEN
        print *,'200 QVG,QSG,QCG,TSO(1)',QVG,QSG,QCG,TSO(1)
    ENDIF
    IF (debug_print ) THEN
     if(iter == 1) then
      print *,'QVATM,QVG,QSG,QCG,TS1',QVATM,QVG,QSG,QCG,TS1
     endif
    ENDIF

!--- SOILT - skin temperature
          SOILT=TS1

!---- Final solution for soil temperature - TSO
          DO K=2,NZS
            KK=NZS-K+1
            TSO(K)=rhtso(KK)+cotso(KK)*TSO(K-1)
          END DO

         X= (cp*rho*r211+rhcs*zsmain(2)*0.5_kind_phys/delt)*(SOILT-TN) + &
            XLV*rho*r211*(QVG-QGOLD) 

    IF (debug_print ) THEN
        print*,'SOILTEMP storage, i,j,x,soilt,tn,qvg,qvgold', &
                                  i,j,x,soilt,tn,qvg,qgold
        print *,'TEMP term (cp*rho*r211+rhcs*zsmain(2)*0.5/delt)*(SOILT-TN)',&
                 (cp*rho*r211+rhcs*zsmain(2)*0.5_kind_phys/delt)*(SOILT-TN)
        print *,'QV term XLV*rho*r211*(QVG-QGOLD)',XLV*rho*r211*(QVG-QGOLD)
    ENDIF
         X=X &
! "heat" from rain
        -RAINF*CVW*PRCPMS*(max(tfrz,TABS)-SOILT)

    IF (debug_print ) THEN
        print *,'x=',x
    ENDIF

!--------------------------------------------------------------------
   END SUBROUTINE SOILTEMP
!--------------------------------------------------------------------

!>\ingroup lsm_ruc_group
!> This subroutine solves energy bugdget equation and heat diffusion 
!! equation to obtain snow and soil temperatures.
           SUBROUTINE SNOWTEMP( debug_print,xlat,xlon,             &
           testptlat,testptlon,i,j,iland,isoil,                    & !--- input variables
           delt,ktau,conflx,nzs,nddzs,nroot,                       &
           isncond_opt,isncovr_opt,                                &
           snwe,snwepr,snhei,newsnow,snowfrac,snhei_crit,          &
           beta,deltsn,snth,rhosn,rhonewsn,meltfactor,             &  ! add meltfactor
           PRCPMS,RAINF,                                           &
           PATM,TABS,QVATM,QCATM,                                  &
           GLW,GSW,EMISS,RNET,                                     &
           QKMS,TKMS,PC,RHO,VEGFRAC,                               &
           THDIF,CAP,DRYCAN,WETCAN,CST,                            &
           TRANF,TRANSUM,DEW,MAVAIL,                               &
           DQM,QMIN,PSIS,BCLH,                                     & !--- soil fixed fields
           ZSMAIN,ZSHALF,DTDZS,TBQ,                                &
           XLVM,CP,rovcp,G0_P,CVW,STBOLT,                          & !--- constants
           SNWEPRINT,SNHEIPRINT,RSM,                               & !--- output variables
           TSO,SOILT,SOILT1,TSNAV,QVG,QSG,QCG,                     &
           SMELT,SNOH,SNFLX,S,ILNB,X)

!********************************************************************
!   Energy budget equation and heat diffusion eqn are 
!   solved here to obtain snow and soil temperatures
!
!     DELT - time step (s)
!     ktau - number of time step
!     CONFLX - depth of constant flux layer (m)
!     IME, JME, KME, NZS - dimensions of the domain 
!     NROOT - number of levels within the root zone
!     PRCPMS - precipitation rate in m/s
!     COTSO, RHTSO - coefficients for implicit solution of
!                     heat diffusion equation
!     THDIF - thermal diffusivity (W/m/K)
!     QSG,QVG,QCG - saturated mixing ratio, mixing ratio of
!                   water vapor and cloud at the ground
!                   surface, respectively (kg/kg)
!     PATM - pressure [bar]
!     QCATM,QVATM - cloud and water vapor mixing ratio
!                   at the first atm. level (kg/kg)
!     EMISS,RNET - emissivity (0-1) of the ground surface and net
!                  radiation at the surface (W/m^2)
!     QKMS - exchange coefficient for water vapor in the
!              surface layer (m/s)
!     TKMS - exchange coefficient for heat in the surface
!              layer (m/s)
!     PC - plant coefficient (resistance)
!     RHO - density of atmosphere near surface (kg/m^3)
!     VEGFRAC - greeness fraction (0-1)
!     CAP - volumetric heat capacity (J/m^3/K)
!     DRYCAN - dry fraction of vegetated area where
!              transpiration may take place (0-1) 
!     WETCAN - fraction of vegetated area covered by canopy
!              water (0-1)
!     TRANSUM - transpiration function integrated over the 
!               rooting zone (m)
!     DEW -  dew in kg/m^2/s
!     MAVAIL - fraction of maximum soil moisture in the top
!               layer (0-1)
!     ZSMAIN - main levels in soil (m)
!     ZSHALF - middle of the soil layers (m)
!     DTDZS - dt/(2.*dzshalf*dzmain)
!     TBQ - table to define saturated mixing ration
!           of water vapor for given temperature and pressure
!     TSO - soil temperature (K)
!     SOILT - skin temperature (K)
!
!*********************************************************************

        IMPLICIT NONE
!---------------------------------------------------------------------
!--- input variables

   LOGICAL,  INTENT(IN   )   ::  debug_print
   INTEGER,  INTENT(IN   )   ::  nroot,ktau,nzs                , &
                                 nddzs                             !nddzs=2*(nzs-2)

   INTEGER,  INTENT(IN   )   ::  i,j,iland,isoil,isncond_opt,isncovr_opt
   real (kind_phys),     INTENT(IN   )  ::  DELT,CONFLX,PRCPMS , &
                                 RAINF,NEWSNOW,DELTSN,SNTH     , &
                                 TABS,TRANSUM,SNWEPR           , &
                                 testptlat,testptlon           , & 
                                 rhonewsn,meltfactor,xlat,xlon,snhei_crit
   real                      ::  rhonewcsn

!--- 3-D Atmospheric variables
   real (kind_phys),                                             &
            INTENT(IN   )    ::                            PATM, &
                                                          QVATM, &
                                                          QCATM
!--- 2-D variables
   real (kind_phys)                                            , &
            INTENT(IN   )    ::                             GLW, &
                                                            GSW, &
                                                            RHO, &
                                                             PC, &
                                                        VEGFRAC, &
                                                           QKMS, &
                                                           TKMS

!--- soil properties
   real (kind_phys)                                            , &
            INTENT(IN   )    ::                                  &
                                                           BCLH, &
                                                            DQM, &
                                                           PSIS, &
                                                           QMIN

   real (kind_phys),     INTENT(IN   )   ::                  CP, &
                                                          ROVCP, &
                                                            CVW, &
                                                         STBOLT, &
                                                           XLVM, &
                                                            G0_P


   real (kind_phys),     DIMENSION(1:NZS), INTENT(IN) :: ZSMAIN, &
                                                         ZSHALF, &
                                                          THDIF, &
                                                            CAP, &
                                                          TRANF 

   real (kind_phys),     DIMENSION(1:NDDZS), INTENT(IN)  :: DTDZS

   real (kind_phys),     DIMENSION(1:5001), INTENT(IN)  ::  TBQ


!--- input/output variables
!-------- 3-d soil moisture and temperature
   real (kind_phys),     DIMENSION(  1:nzs )                   , &
             INTENT(INOUT)   ::                             TSO


!-------- 2-d variables
   real (kind_phys)                                            , &
             INTENT(INOUT)   ::                             DEW, &
                                                            CST, &
                                                          RHOSN, &
                                                          EMISS, &
                                                         MAVAIL, &
                                                            QVG, &
                                                            QSG, &
                                                            QCG, &
                                                           SNWE, &
                                                          SNHEI, &
                                                       SNOWFRAC, &
                                                          SMELT, &
                                                           SNOH, &
                                                          SNFLX, &
                                                              S, &
                                                          SOILT, &
                                                         SOILT1, &
                                                          TSNAV

   real (kind_phys),     INTENT(INOUT)                  ::   DRYCAN, WETCAN           

   real (kind_phys),     INTENT(OUT)                    ::  RSM, &
                                                      SNWEPRINT, &
                                                     SNHEIPRINT
   INTEGER,  INTENT(OUT)                    ::             ilnb
!--- Local variables


   INTEGER ::  nzs1,nzs2,k,k1,kn,kk

   real (kind_phys)    ::  x,x1,x2,x4,dzstop,can,ft,sph,         &
               tn,trans,umveg,denom

   real (kind_phys)    ::  cotsn,rhtsn,xsn1,ddzsn1,x1sn1,ftsnow,denomsn

   real (kind_phys)    ::  t3,upflux,xinet,ras,                  & 
               xlmelt,rhocsn,thdifsn,                            &
               beta,epot,xsn,ddzsn,x1sn,d1sn,d2sn,d9sn,r22sn

   real (kind_phys)    ::  fso,fsn,                              &
               FKT,D1,D2,D9,D10,DID,R211,R21,R22,R6,R7,D11,      &
               PI,H,FKQ,R210,AA,BB,PP,Q1,QS1,TS1,TQ2,TX2,        &
               TDENOM,C,CC,AA1,RHCS,H1,                          &
               tsob, snprim, sh1, sh2,                           &
               smeltg,snohg,snodif,soh,                          &
               CMC2MS,TNOLD,QGOLD,SNOHGNEW                            

   real (kind_phys),     DIMENSION(1:NZS)  ::  transp,cotso,rhtso
   real (kind_phys)                        ::             edir1, &
                                                            ec1, &
                                                           ett1, &
                                                           eeta, &
                                                            qfx, &
                                                            hfx

   real (kind_phys)       :: RNET,rsmfrac,soiltfrac,hsn,rr,keff,fact
   integer                     :: nmelt, iter

!-----------------------------------------------------------------

       iter = 0

       !-- options for snow conductivity:
       !-- 1 - constant
       !-- opt 2 -  Sturm et al., 1997
       keff = 0.265_kind_phys

       do k=1,nzs
          transp   (k)=zero
          cotso    (k)=zero
          rhtso    (k)=zero
       enddo

    IF (debug_print ) THEN
print *, 'SNOWTEMP: SNHEI,SNTH,SOILT1: ',SNHEI,SNTH,SOILT1,soilt 
    ENDIF
        XLMELT=con_hfus
        RHOCSN=sheatsn* RHOSN
        RHOnewCSN=sheatsn* RHOnewSN
        if(isncond_opt == 1) then
        !-- old version thdifsn = 0.265/RHOCSN
          THDIFSN = 0.265_kind_phys/RHOCSN
        else
        !-- 07Jun19 - thermal conductivity (K_eff) from Sturm et al.(1997)
        !-- keff = 10. ** (2.650 * RHOSN*1.e-3 - 1.652)
           fact = one
           if(rhosn < 156._kind_phys .or. (newsnow > zero .and. rhonewsn < 156._kind_phys)) then
             keff = 0.023_kind_phys + 0.234_kind_phys * rhosn * 1.e-3_kind_phys
           else
             keff = 0.138_kind_phys - 1.01_kind_phys * rhosn*1.e-3_kind_phys + 3.233_kind_phys * rhosn**2 * 1.e-6_kind_phys
             if(debug_print) then
               print *,'SnowTemp xlat,xlon,rhosn,keff', xlat,xlon,rhosn,keff,keff/rhocsn*fact
               print *,'SNOWTEMP - 0.265/rhocsn',0.265_kind_phys/rhocsn
             endif
           endif
       if ( debug_print .and. abs(xlat-testptlat).lt.0.2  .and. abs(xlon-testptlon).lt.0.2) then
           print *,'SNOWTEMP - xlat,xlon,newsnow,rhonewsn,rhosn,fact,keff',xlat,xlon,newsnow, rhonewsn,rhosn,fact,keff
       endif

         if(newsnow <= zero .and. snhei > one .and. rhosn > 250._kind_phys) then
           !-- some areas with large snow depth have unrealistically 
           !-- low snow density (in the Rockie's with snow depth > 1 m). 
           !-- Based on Sturm et al. keff=0.452 typical for hard snow slabs
           !-- with rhosn=488 kg/m^3. Thdifsn = 0.452/(2090*488)=4.431718e-7
           !-- In future a better compaction scheme is needed for these areas.
             thdifsn = 4.431718e-7_kind_phys
           else
             thdifsn = keff/rhocsn * fact
           endif
       if (debug_print .and. abs(xlat-testptlat).lt.0.2  .and. abs(xlon-testptlon).lt.0.2) then
           print *,'SNOWTEMP - thdifsn',xlat,xlon,thdifsn
           print *,'SNOWTEMP - 0.265/rhocsn',0.265_kind_phys/rhocsn
       endif

        endif

        RAS=RHO*1.E-3_kind_phys

        SOILTFRAC=SOILT

        SMELT=zero
        SOH=zero
        SMELTG=zero
        SNOHG=zero
        SNODIF=zero
        RSM = zero
        RSMFRAC = zero
        fsn=one
        fso=zero

          NZS1=NZS-1
          NZS2=NZS-2

        QGOLD=QVG
        DZSTOP=one/(ZSMAIN(2)-ZSMAIN(1))

!******************************************************************************
!       COEFFICIENTS FOR THOMAS ALGORITHM FOR TSO
!******************************************************************************
        cotso(1)=zero
        rhtso(1)=TSO(NZS)
        DO 33 K=1,NZS2
          KN=NZS-K
          K1=2*KN-3
          X1=DTDZS(K1)*THDIF(KN-1)
          X2=DTDZS(K1+1)*THDIF(KN)
          FT=TSO(KN)+X1*(TSO(KN-1)-TSO(KN))                           &
             -X2*(TSO(KN)-TSO(KN+1))
          DENOM=1.+X1+X2-X2*cotso(K)
          cotso(K+1)=X1/DENOM
          rhtso(K+1)=(FT+X2*rhtso(K))/DENOM
   33  CONTINUE
!--- THE NZS element in COTSO and RHTSO will be for snow
!--- There will be 2 layers in snow if it is deeper than DELTSN+SNTH
       IF(SNHEI.GE.SNTH) then
        if(snhei.le.DELTSN+SNTH) then
!-- 1-layer snow model
    IF (debug_print ) THEN
      print *,'1-layer - snth,snhei,deltsn',snth,snhei,deltsn
    ENDIF
         ilnb=1
         snprim=max(snth,snhei)
         tsob=tso(1)
         soilt1=tso(1)
         XSN = DELT/2._kind_phys/(zshalf(2)+0.5_kind_phys*SNPRIM)
         DDZSN = XSN / SNPRIM
         X1SN = DDZSN * thdifsn
         X2 = DTDZS(1)*THDIF(1)
         FT = TSO(1)+X1SN*(SOILT-TSO(1))                              &
              -X2*(TSO(1)-TSO(2))
         DENOM = one + X1SN + X2 -X2*cotso(NZS1)
         cotso(NZS)=X1SN/DENOM
         rhtso(NZS)=(FT+X2*rhtso(NZS1))/DENOM
         cotsn=cotso(NZS)
         rhtsn=rhtso(NZS)
!*** Average temperature of snow pack (C)
         tsnav=min(zero,0.5_kind_phys*(soilt+tso(1))-tfrz)

        else
!-- 2 layers in snow, SOILT1 is temperasture at DELTSN depth
    IF (debug_print ) THEN
      print *,'2-layer - snth,snhei,deltsn',snth,snhei,deltsn
    ENDIF
         ilnb=2
         snprim=deltsn
         tsob=soilt1
         XSN = DELT/2._kind_phys/(0.5_kind_phys*deltsn)
         XSN1= DELT/2._kind_phys/(zshalf(2)+0.5_kind_phys*(SNHEI-DELTSN))
         DDZSN = XSN / DELTSN
         DDZSN1 = XSN1 / (SNHEI-DELTSN)
         X1SN = DDZSN * thdifsn
         X1SN1 = DDZSN1 * thdifsn
         X2 = DTDZS(1)*THDIF(1)
         FT = TSO(1)+X1SN1*(SOILT1-TSO(1))                            &
              -X2*(TSO(1)-TSO(2))
         DENOM = 1. + X1SN1 + X2 - X2*cotso(NZS1)
         cotso(nzs)=x1sn1/denom
         rhtso(nzs)=(ft+x2*rhtso(nzs1))/denom
         ftsnow = soilt1+x1sn*(soilt-soilt1)                          &
               -x1sn1*(soilt1-tso(1))
         denomsn = one + X1SN + X1SN1 - X1SN1*cotso(NZS)
         cotsn=x1sn/denomsn
         rhtsn=(ftsnow+X1SN1*rhtso(NZS))/denomsn
!*** Average temperature of snow pack (C)
         tsnav=min(zero,0.5_kind_phys/snhei*((soilt+soilt1)*deltsn    &
                     +(soilt1+tso(1))*(SNHEI-DELTSN))                 &
                     -tfrz)
        endif
       ENDIF
       IF(SNHEI.LT.SNTH.AND.SNHEI.GT.zero) then
!--- snow is too thin to be treated separately, therefore it
!--- is combined with the first soil layer.
         snprim=SNHEI+zsmain(2)
         fsn=SNHEI/snprim
         fso=one-fsn
         soilt1=tso(1)
         tsob=tso(2)
         XSN = DELT/2._kind_phys/((zshalf(3)-zsmain(2))+0.5_kind_phys*snprim)
         DDZSN = XSN /snprim
         X1SN = DDZSN * (fsn*thdifsn+fso*thdif(1))
         X2=DTDZS(2)*THDIF(2)
         FT=TSO(2)+X1SN*(SOILT-TSO(2))-                              &
                       X2*(TSO(2)-TSO(3))
         denom = one + x1sn + x2 - x2*cotso(nzs-2)
         cotso(nzs1) = x1sn/denom
         rhtso(nzs1)=(FT+X2*rhtso(NZS-2))/denom
         tsnav=min(zero,0.5_kind_phys*(soilt+tso(1))                 &
                     -tfrz)
         cotso(NZS)=cotso(nzs1)
         rhtso(NZS)=rhtso(nzs1)
         cotsn=cotso(NZS)
         rhtsn=rhtso(NZS)

       ENDIF

!************************************************************************
!--- THE HEAT BALANCE EQUATION (Smirnova et al. 1996, EQ. 21,26)
!18apr08 nmelt is the flag for melting, and SNOH is heat of snow phase changes
       nmelt=0
       SNOH=zero

        ETT1=zero
        EPOT=-QKMS*(QVATM-QGOLD)
        RHCS=CAP(1)
        H=MAVAIL !1.
        TRANS=TRANSUM*DRYCAN/ZSHALF(NROOT+1)
        CAN=WETCAN+TRANS
        UMVEG=one-VEGFRAC
        FKT=TKMS
        D1=cotso(NZS1)
        D2=rhtso(NZS1)
        TN=SOILT
        D9=THDIF(1)*RHCS*dzstop
        D10=TKMS*CP*RHO
        R211=.5_kind_phys*CONFLX/DELT
        R21=R211*CP*RHO
        R22=.5_kind_phys/(THDIF(1)*DELT*dzstop**2)
        R6=EMISS *STBOLT*.5_kind_phys*TN**4
        R7=R6/TN
        D11=RNET+R6

      IF(SNHEI.GE.SNTH) THEN
        if(snhei.le.DELTSN+SNTH) then
!--- 1-layer snow
          D1SN = cotso(NZS)
          D2SN = rhtso(NZS)
    IF (debug_print ) THEN
      print *,'1 layer d1sn,d2sn',i,j,d1sn,d2sn
    ENDIF
        else
!--- 2-layer snow
          D1SN = cotsn
          D2SN = rhtsn
    IF (debug_print ) THEN
      print *,'2 layers d1sn,d2sn',i,j,d1sn,d2sn
    ENDIF
        endif
        D9SN= THDIFSN*RHOCSN / SNPRIM
        R22SN = SNPRIM*SNPRIM*0.5_kind_phys/(THDIFSN*DELT)
    IF (debug_print ) THEN
      print *,'1 or 2 layers D9sn,R22sn',d9sn,r22sn
    ENDIF
      ENDIF

       IF(SNHEI.LT.SNTH.AND.SNHEI.GT.zero) then
!--- thin snow is combined with soil
         D1SN = D1
         D2SN = D2
         D9SN = (fsn*THDIFSN*RHOCSN+fso*THDIF(1)*RHCS)/              &
                 snprim
         R22SN = snprim*snprim*0.5_kind_phys                         &
                 /((fsn*THDIFSN+fso*THDIF(1))*delt)
    IF (debug_print ) THEN
       print *,' Combined  D9SN,R22SN,D1SN,D2SN: ',D9SN,R22SN,D1SN,D2SN
    ENDIF
      ENDIF
      IF(SNHEI.eq.zero)then
!--- all snow is sublimated
        D9SN = D9
        R22SN = R22
        D1SN = D1
        D2SN = D2
    IF (debug_print ) THEN
        print *,' SNHEI = 0, D9SN,R22SN,D1SN,D2SN: ',D9SN,R22SN,D1SN,D2SN
    ENDIF
      ENDIF

 2211   continue

!18apr08  - the snow melt iteration start point
 212    continue

!---- TDENOM for snow
        TDENOM = D9SN*(one-D1SN +R22SN)+D10+R21+R7                   &
              +RAINF*CVW*PRCPMS                                      &
              +RHOnewCSN*NEWSNOW/DELT

        FKQ=QKMS*RHO
        R210=R211*RHO
        C=VEGFRAC*FKQ*CAN
        CC=C*XLVM/TDENOM
        AA=XLVM*(BETA*FKQ*UMVEG+R210)/TDENOM
        BB=(D10*TABS+R21*TN+XLVM*(QVATM*                             &
        (BETA*FKQ*UMVEG+C)                                           &
        +R210*QGOLD)+D11+D9SN*(D2SN+R22SN*TN)                        &
        +RAINF*CVW*PRCPMS*max(tfrz,TABS)                             &
        + RHOnewCSN*NEWSNOW/DELT*min(tfrz,TABS)                      &
         )/TDENOM
        AA1=AA+CC
        PP=PATM*1.E3_kind_phys
        AA1=AA1/PP
        BB=BB-SNOH/TDENOM

      IF (debug_print ) THEN
        if (abs(xlat-33.35).lt.0.2 .and. abs(xlon-272.55).lt.0.2)then
          print *,'1-', i,rnet,tabs,tn,aa1,bb,pp,ktau,newsnow,snwepr,snwe,snhei,snowfrac,soilt,soilt1,tso,rhosn
          print *,'2-', i,tdenom,fkq,vegfrac,can,R210,D10,R21,D9sn,D1sn,R22sn,R7,prcpms
        endif
      ENDIF
        CALL VILKA(TN,AA1,BB,PP,QS1,TS1,TBQ,KTAU,i,j,iland,isoil,xlat,xlon)
        TQ2=QVATM
        TX2=TQ2*(one-H)
        Q1=TX2+H*QS1
    IF (debug_print ) THEN
    !if (abs(xlat-33.35).lt.0.2 .and. abs(xlon-272.55).lt.0.2)then
     print *,'VILKA1 - TS1,QS1,TQ2,H,TX2,Q1',TS1,QS1,TQ2,H,TX2,Q1,xlat,xlon
    ENDIF
        IF(Q1.LT.QS1) GOTO 100
!--- if no saturation - goto 100
!--- if saturation - goto 90
   90   QVG=QS1
        QSG=QS1
        QCG=max(zero,Q1-QS1)
    IF (debug_print ) THEN
     print *,'90 QVG,QSG,QCG,TSO(1)',QVG,QSG,QCG,TSO(1)
    ENDIF
        GOTO 200
  100   BB=BB-AA*TX2
        AA=(AA*H+CC)/PP
        CALL VILKA(TN,AA,BB,PP,QS1,TS1,TBQ,KTAU,i,j,iland,isoil,xlat,xlon)
        Q1=TX2+H*QS1
    IF (debug_print ) THEN
    !if (abs(xlat-33.35).lt.0.2 .and. abs(xlon-272.55).lt.0.2)then
     print *,'VILKA2 - TS1,QS1,H,TX2,Q1',TS1,QS1,TQ2,H,TX2,Q1
    ENDIF
        IF(Q1.GT.QS1) GOTO 90
        QSG=QS1
        QVG=Q1
        QCG=zero
    IF (debug_print ) THEN
     print *,'No Saturation QVG,QSG,QCG,TSO(1)',QVG,QSG,QCG,TSO(1)
    ENDIF
  200   CONTINUE

!--- SOILT - skin temperature
        SOILT=TS1
     if(nmelt==1 .and. snowfrac==one .and. snwe > zero .and. SOILT > tfrz) then
     !--7feb22 on the second iteration when SNOH is known and snwe > 0. after melting,
     !-- check if the snow skin temperature is =<tfrzK
     !-- when a grid cell is fully covered with snow (snowfrac=1) 
     !-- or with partial snow cover and snow_mosaic=1 (snowfrac=1).
       if (debug_print ) then
       !if (abs(xlat-33.35).lt.0.2 .and. abs(xlon-272.55).lt.0.2)then
         print *,'soilt is too high =',soilt,xlat,xlon
         soilt = min(tfrz,soilt)
       endif
     endif

    IF (debug_print ) THEN
       !if (abs(xlat-33.35).lt.0.2 .and. abs(xlon-272.55).lt.0.2)then
       print *,'snwe,snwepr,snhei,snowfr,soilt,soilt1,tso',i,j,snwe,snwepr,snhei,snowfrac,soilt,soilt1,tso
    ENDIF
! Solution for temperature at 7.5 cm depth and snow-soil interface
       IF(SNHEI.GE.SNTH) THEN
        if(snhei.gt.DELTSN+SNTH) then
!-- 2-layer snow model
          SOILT1=rhtsn+cotsn*SOILT
          TSO(1)=rhtso(NZS)+cotso(NZS)*SOILT1
          tsob=soilt1
        else
!-- 1 layer in snow
          TSO(1)=rhtso(NZS)+cotso(NZS)*SOILT
          SOILT1=TSO(1)
          tsob=tso(1)
        endif
       ELSEIF (SNHEI > zero .and. SNHEI < SNTH) THEN
! blended 
         TSO(2)=rhtso(NZS1)+cotso(NZS1)*SOILT
         tso(1)=(tso(2)+(soilt-tso(2))*fso)
         SOILT1=TSO(1)
         tsob=TSO(2)
       ELSE
!-- very thin or zero snow. If snow is thin we suppose that
!--- tso(i,j,1)=SOILT, and later we recompute tso(i,j,1)
         TSO(1)=SOILT
         SOILT1=SOILT
         tsob=TSO(1)
       ENDIF
       if(nmelt==1.and.snowfrac==one) then
       !-- second iteration with full snow cover
         SOILT1= min(tfrz,SOILT1)
         TSO(1)= min(tfrz,TSO(1))
         tsob  = min(tfrz,tsob)
       endif

!---- Final solution for TSO
       IF (SNHEI > zero .and. SNHEI < SNTH) THEN
! blended or snow is melted
          DO K=3,NZS
            KK=NZS-K+1
            TSO(K)=rhtso(KK)+cotso(KK)*TSO(K-1)
          END DO

       ELSE
          DO K=2,NZS
            KK=NZS-K+1
            TSO(K)=rhtso(KK)+cotso(KK)*TSO(K-1)
          END DO
       ENDIF
!--- For thin snow layer combined with the top soil layer
!--- TSO(1) is recomputed by linear interpolation between SOILT
!--- and TSO(i,j,2)
!       if(SNHEI.LT.SNTH.AND.SNHEI.GT.0.)then
!          tso(1)=tso(2)+(soilt-tso(2))*fso
!          soilt1=tso(1)
!          tsob = tso(2)
!       endif


    IF (debug_print ) THEN
    !if (abs(xlat-33.35).lt.0.2 .and. abs(xlon-272.55).lt.0.2)then
      print *,'Final SOILT,SOILT1,tso,TSOB,QSG',xlat,xlon,SOILT,SOILT1,tso,TSOB,QSG,'nmelt=',nmelt
      print *,'SNWEPR-BETA*EPOT*RAS*DELT',SNWEPR-BETA*EPOT*RAS*DELT,beta,snwepr,epot
    ENDIF

     if(nmelt.eq.1) go to 220

!--- IF SOILT > tfrz F then melting of snow can happen
! if all snow can evaporate (beta<1), then there is nothing to melt
   IF(SOILT > tfrz.AND.BETA==one.AND.SNHEI>zero) THEN
     !-- snow sublimation and melting
        nmelt = 1
        soiltfrac=snowfrac*tfrz+(one-snowfrac)*SOILT
        QSG=min(QSG, QSN(soiltfrac,TBQ)/PP)
        qvg=snowfrac*qsg+(one-snowfrac)*qvg
        T3      = STBOLT*TN*TN*TN
        UPFLUX  = T3 * 0.5_kind_phys*(TN + SOILTfrac)
        XINET   = EMISS*(GLW-UPFLUX)
         EPOT = -QKMS*(QVATM-QSG)
         Q1=EPOT*RAS


        IF (Q1.LE.0..or.iter==1) THEN
! ---  condensation
          DEW=-EPOT
          DO K=1,NZS
            TRANSP(K)=0.
          ENDDO

        QFX = -XLVM*RHO*DEW
        EETA = QFX/XLVM
       ELSE
! ---  evaporation
          DO K=1,NROOT
            TRANSP(K)=-VEGFRAC*q1                                     &
                      *TRANF(K)*DRYCAN/zshalf(NROOT+1)
            ETT1=ETT1-TRANSP(K)
          ENDDO
          DO k=nroot+1,nzs
            transp(k)=0.
          enddo

        EDIR1 = Q1*UMVEG * BETA
        EC1 = Q1 * WETCAN * vegfrac
        CMC2MS=CST/DELT*RAS
        EETA = (EDIR1 + EC1 + ETT1)*rhowater
! to convert from kg m-2 s-1 to m s-1: 1/rho water=1.e-3************ 
        QFX=  XLVM * EETA
       ENDIF

         HFX=-D10*(TABS-soiltfrac)

       IF(SNHEI.GE.SNTH)then
         SOH=thdifsn*RHOCSN*(soiltfrac-TSOB)/SNPRIM
         SNFLX=SOH
       ELSE
         SOH=(fsn*thdifsn*rhocsn+fso*thdif(1)*rhcs)*                   &
              (soiltfrac-TSOB)/snprim
         SNFLX=SOH
       ENDIF

!
         X= (R21+D9SN*R22SN)*(soiltfrac-TN) +                          &
            XLVM*R210*(QVG-QGOLD)
    IF (debug_print ) THEN
      print *,'SNOWTEMP storage ',i,j,x
      print *,'R21,D9sn,r22sn,soiltfrac,tn,qsg,qvg,qgold,snprim', &
              R21,D9sn,r22sn,soiltfrac,tn,qsg,qvg,qgold,snprim
    ENDIF

!-- SNOH is energy flux of snow phase change
        SNOH=RNET-QFX -HFX - SOH - X                                    & 
                  +RHOnewCSN*NEWSNOW/DELT*(min(tfrz,TABS)-soiltfrac)  &
                  +RAINF*CVW*PRCPMS*(max(tfrz,TABS)-soiltfrac) 
        SNOH=AMAX1(0.,SNOH)
!-- SMELT is speed of melting in M/S
        SMELT= SNOH /XLMELT*1.E-3_kind_phys
    IF (debug_print ) THEN
    !if (abs(xlat-33.35).lt.0.2 .and. abs(xlon-272.55).lt.0.2)then
      print *,'1- SMELT',smelt,snoh,xlat,xlon
    ENDIF

      IF(EPOT.gt.zero .and. SNWEPR.LE.EPOT*RAS*DELT) THEN
!-- all snow can evaporate
        BETA=SNWEPR/(EPOT*RAS*DELT)
        SMELT=AMAX1(zero,AMIN1(SMELT,SNWEPR/DELT-BETA*EPOT*RAS))
        SNWE=zero
       IF (debug_print ) THEN
       !if (abs(xlat-33.35).lt.0.2 .and. abs(xlon-272.55).lt.0.2)then
         print *,'2- SMELT',xlat,xlon,snwe,smelt,rhonewsn,xlat,xlon
       ENDIF
          goto 88
      ENDIF

!18apr08 - Egglston limit
      !-- 22apr22 Do not limit snow melting for hail (rhonewsn > 450), or dense snow
      !-- (rhosn > 350.) with very warm surface temperatures (>10C)
      if( (rhosn < 350._kind_phys .or. (newsnow > zero .and. rhonewsn < 450._kind_phys)) .and. soilt < 283._kind_phys ) then
        SMELT= amin1 (smelt, delt/60._kind_phys*5.6E-8_kind_phys*meltfactor*max(one,(soilt-tfrz))) 
        IF (debug_print ) THEN
        !if (abs(xlat-33.35).lt.0.2 .and. abs(xlon-272.55).lt.0.2)then
          print *,'3- SMELT',xlat,xlon,smelt,rhosn,rhonewsn,xlat,xlon
        ENDIF
      endif

! rr - potential melting
        rr=max(zero,SNWEPR/delt-BETA*EPOT*RAS)
        if(smelt > rr) then
          SMELT = min(SMELT,rr)
          SNWE = zero
         IF (debug_print ) THEN
         !if (abs(xlat-33.35).lt.0.2 .and. abs(xlon-272.55).lt.0.2)then
           print *,'4- SMELT i,j,smelt,rr',xlat,xlon,smelt,rr
         ENDIF
        endif
 88   continue
        SNOHGNEW=SMELT*XLMELT*rhowater
        SNODIF=AMAX1(zero,(SNOH-SNOHGNEW))

        SNOH=SNOHGNEW
        IF (debug_print ) THEN
        !if (abs(xlat-33.35).lt.0.2 .and. abs(xlon-272.55).lt.0.2)then
          print *,'SNOH,SNODIF',SNOH,SNODIF
          print *,' xlat, xlon', xlat, xlon
        ENDIF

      IF( smelt > zero) then
!*** From Koren et al. (1999) 13% of snow melt stays in the snow pack
        rsmfrac=min(0.18_kind_phys,(max(0.08_kind_phys,snwepr/0.10_kind_phys*0.13_kind_phys)))
       if(snhei > 0.01_kind_phys .and. rhosn < 350._kind_phys) then
        rsm=min(snwe,rsmfrac*smelt*delt)
       else
       ! do not keep melted water if snow depth is less that 1 cm
       ! or if snow is dense
        rsm=zero
       endif
!18apr08 rsm is part of melted water that stays in snow as liquid
       if(rsm > zero) then
         SMELT=max(zero,SMELT-rsm/delt)
         IF (debug_print ) THEN
         !if (abs(xlat-33.35).lt.0.2 .and. abs(xlon-272.55).lt.0.2)then
           print *,'5- SMELT i,j,smelt,rsm,snwepr,rsmfrac', &
                        i,j,smelt,rsm,snwepr,rsmfrac
           print *,' xlat, xlon', xlat, xlon
         ENDIF
       endif ! rsm

      ENDIF ! smelt > 0

!-- update of liquid equivalent of snow depth
!-- due to evaporation and snow melt
      if(snwe > zero) then
        SNWE = AMAX1(zero,(SNWEPR-                        &
               (SMELT+BETA*EPOT*RAS)*DELT                 &
                                         ) )
        IF (debug_print ) THEN
        !if (abs(xlat-33.35).lt.0.2 .and. abs(xlon-272.55).lt.0.2)then
          print *,' Snow is melting and sublimating, snwe', xlat, xlon, SNWE
        endif
      else
       !-- all snow is sublimated or melted
         IF (debug_print ) THEN
         !if (abs(xlat-33.35).lt.0.2 .and. abs(xlon-272.55).lt.0.2)then
          print *,' all snwe is sublimated or melted', xlat, xlon, SNWE
         endif
       endif
      ELSE
      !-- NO MELTING, only sublimation
      !--- If there is no snow melting then just evaporation
      !--- or condensation changes SNWE
       if(snhei.ne.zero .and. beta == one) then
               EPOT=-QKMS*(QVATM-QSG)
               SNWE = AMAX1(zero,(SNWEPR-                  &
                    BETA*EPOT*RAS*DELT))
       else
       !-- all snow is sublibated
         snwe = zero
       endif

      ENDIF

!18apr08 - if snow melt occurred then go into iteration for energy budget 
!         solution 
     if(nmelt.eq.1) goto 212  ! second interation
 220  continue

      if(smelt > zero .and. rsm > zero) then
       if(snwe.le.rsm) then
    IF ( debug_print ) THEN
       print *,'SNWE<RSM snwe,rsm,smelt*delt,epot*ras*delt,beta', &
                     snwe,rsm,smelt*delt,epot*ras*delt,beta
    ENDIF
       else
!*** Update snow density on effect of snow melt, melted
!*** from the top of the snow. 13% of melted water
!*** remains in the pack and changes its density.
!*** Eq. 9 (with my correction) in Koren et al. (1999)
          xsn=(rhosn*(snwe-rsm)+rhowater*rsm)/                    &
              snwe
          rhosn=MIN(MAX(58.8_kind_phys,XSN),500._kind_phys)

          RHOCSN=sheatsn* RHOSN
          if(isncond_opt == 1) then
          !-- old version thdifsn = 0.265/RHOCSN
            THDIFSN = 0.265_kind_phys/RHOCSN
          else
          !-- 07Jun19 - thermal conductivity (K_eff) from Sturm et al.(1997)
          !-- keff = 10. ** (2.650 * RHOSN*1.e-3 - 1.652)
            fact = one
            if(rhosn < 156._kind_phys .or. (newsnow > zero .and. rhonewsn < 156._kind_phys)) then
              keff = 0.023_kind_phys + 0.234_kind_phys * rhosn * 1.e-3_kind_phys
            else
              keff = 0.138_kind_phys - 1.01_kind_phys * rhosn*1.e-3_kind_phys + 3.233_kind_phys * rhosn**2 * 1.e-6_kind_phys
              if(debug_print) then
                print *,'End SNOWTEMP - xlat,xlon,rhosn,keff',xlat,xlon,rhosn,keff
                print *,'End SNOWTEMP - 0.265/rhocsn',0.265/rhocsn
              endif
            endif
       if (debug_print .and. abs(xlat-testptlat).lt.0.2 .and. abs(xlon-testptlon).lt.0.2) then
           print *,'END SNOWTEMP - newsnow, rhonewsn,rhosn,fact,keff', &
                    xlat,xlon,newsnow, rhonewsn,rhosn,fact,keff,keff/rhocsn*fact
       endif

         if(newsnow <= zero .and. snhei > one .and. rhosn > 250._kind_phys) then
            !-- some areas with large snow depth have unrealistically 
            !-- low snow density (in the Rockie's with snow depth > 1 m). 
            !-- Based on Sturm et al. keff=0.452 typical for hard snow slabs
            !-- with rhosn=488 kg/m^3. Thdifsn = 0.452/(2090*488)=4.431718e-7
            !-- In future a better compaction scheme is needed for these areas.
              thdifsn = 4.431718e-7_kind_phys
            else
              thdifsn = keff/rhocsn * fact
            endif

          endif
       if (debug_print .and.  abs(xlat-testptlat).lt.0.2  .and. abs(xlon-testptlon).lt.0.2) then
           print *,'END SNOWTEMP - thdifsn',xlat,xlon,thdifsn
           print *,'END SNOWTEMP - 0.265/rhocsn',0.265/rhocsn
       endif
        endif  
       endif

!--- Compute flux in the top snow layer
       IF(SNHEI.GE.SNTH)then
         S=thdifsn*RHOCSN*(soilt-TSOB)/SNPRIM
         SNFLX=S
         S=D9*(tso(1)-tso(2))
       ELSEIF(SNHEI.lt.SNTH.and.SNHEI.GT.zero) then
         S=(fsn*thdifsn*rhocsn+fso*thdif(1)*rhcs)*                   &
              (soilt-TSOB)/snprim
         SNFLX=S
         S=D9*(tso(1)-tso(2))
       ELSE
         S=D9SN*(SOILT-TSOB)
         SNFLX=S
         S=D9*(tso(1)-tso(2))
       ENDIF

        !-- Update snow depth after melting at the interface with the atmosphere
        SNHEI=SNWE * rhowater / RHOSN

!--  If ground surface temperature
!--  is above freezing snow can melt from the bottom at the interface with soild. The following
!--  piece of code will check if bottom melting is possible.

    IF (debug_print ) THEN
    !if (abs(xlat-33.35).lt.0.2 .and. abs(xlon-272.55).lt.0.2)then
      print *,'snhei,snwe,rhosn,snowfr',snhei,snwe,rhosn,snowfrac,xlat,xlon
    endif

        IF(TSO(1).GT.tfrz .and. snhei > zero) THEN
!-- melting at the soil/snow interface
          if (snhei.GT.deltsn+snth) then
              hsn = snhei - deltsn
            IF (debug_print ) THEN
              print*,'2 layer snow - snhei,hsn',snhei,hsn
            ENDIF
          else
            IF (debug_print ) THEN
              print*,'1 layer snow or blended - snhei',snhei
            ENDIF
              hsn = snhei
          endif

         soiltfrac=snowfrac*tfrz+(one-snowfrac)*TSO(1)

         SNOHG=(TSO(1)-soiltfrac)*(cap(1)*zshalf(2)+                       &
               RHOCSN*0.5_kind_phys*hsn) / DELT
         SNOHG=AMAX1(zero,SNOHG)
         SNODIF=zero
         SMELTG=SNOHG/XLMELT*1.E-3_kind_phys
         IF (debug_print ) THEN
         !if (abs(xlat-33.35).lt.0.2 .and. abs(xlon-272.55).lt.0.2)then
           print *,' SMELTG =',smeltg,xlat,xlon
         endif
! Egglston - empirical limit on snow melt from the bottom of snow pack
      !9jun22-- the next line excludeis cases of summer hail from snowmelt limiting
      if( (rhosn < 350._kind_phys .or. (newsnow > zero .and. rhonewsn < 450._kind_phys)) .and. soilt < 283._kind_phys ) then
        SMELT=AMIN1(SMELTG, 5.8e-9_kind_phys)
      endif

! rr - potential melting
        rr=SNWE/delt
        SMELTG=AMIN1(SMELTG, rr)

        SNOHGNEW=SMELTG*XLMELT*rhowater
        SNODIF=AMAX1(zero,(SNOHG-SNOHGNEW))
    IF (debug_print ) THEN
    !if (abs(xlat-33.35).lt.0.2 .and. abs(xlon-272.55).lt.0.2)then
       print *,'TSO(1),soiltfrac,snowfrac,smeltg,SNODIF',TSO(1),soiltfrac,snowfrac,smeltg,SNODIF
       print *,' xlat, xlon', xlat, xlon
    ENDIF

        snwe=max(zero,snwe-smeltg*delt)
        SNHEI=SNWE * rhowater / RHOSN
        !-- add up all snow melt
        SMELT = SMELT + SMELTG
      
        if(snhei > zero) TSO(1) = soiltfrac

    IF (debug_print ) THEN
    !if (abs(xlat-33.35).lt.0.2 .and. abs(xlon-272.55).lt.0.2)then
       print *,'Melt from the bottom snwe,snhei',snwe,snhei
       print *,' xlat, xlon', xlat, xlon
       print *,'TSO(1),soiltfrac,snowfrac,smeltg,SNODIF',TSO(1),soiltfrac,snowfrac,smeltg,SNODIF
       print *,'Melt from the bottom snwe,snhei,snoh',snwe,snhei,snoh
       print *,' Final TSO ',tso
       if (snhei==zero) &
       print *,'Snow is all melted on the warm ground'
    ENDIF

        ENDIF ! melt on snow/soil interface

        snweprint=snwe
        snheiprint=snweprint*rhowater / RHOSN

        X= (R21+D9SN*R22SN)*(soilt-TN) +                            &
            XLVM*R210*(QSG-QGOLD)
    IF (debug_print ) THEN
    !if (abs(xlat-33.35).lt.0.2 .and. abs(xlon-272.55).lt.0.2)then
      print *,'end SNOWTEMP storage ',xlat,xlon,x
      print *,'R21,D9sn,r22sn,soiltfrac,soilt,tn,qsg,qgold,snprim', &
              R21,D9sn,r22sn,soiltfrac,soilt,tn,qsg,qgold,snprim
      print *,'snwe, snhei ',snwe,snhei
    ENDIF

         X=X &
! "heat" from snow and rain
        -RHOnewCSN*NEWSNOW/DELT*(min(tfrz,TABS)-SOILT)              &
        -RAINF*CVW*PRCPMS*(max(tfrz,TABS)-SOILT)
    IF (debug_print ) THEN
     print *,'x=',x
     print *,'SNHEI=',snhei
     print *,'SNFLX=',snflx
    ENDIF

      IF(SNHEI.GT.zero) THEN
        if(ilnb.gt.1) then
          tsnav=min(zero,0.5_kind_phys/snhei*((soilt+soilt1)*deltsn &
                    +(soilt1+tso(1))*(SNHEI-DELTSN))                &
                       -tfrz)
        else
          tsnav=min(zero,0.5_kind_phys*(soilt+tso(1)) - tfrz)
        endif
      ELSE
          tsnav= min(zero,soilt - tfrz)
      ENDIF

!------------------------------------------------------------------------
   END SUBROUTINE SNOWTEMP
!------------------------------------------------------------------------

!>\ingroup lsm_ruc_group
!! This subroutine solves moisture budget and computes soil moisture
!! and surface and sub-surface runoffs.
        SUBROUTINE SOILMOIST ( debug_print,                     &
              xlat, xlon, testptlat, testptlon,                 &
              DELT,NZS,NDDZS,DTDZS,DTDZS2,RIW,                  & !--- input parameters
              ZSMAIN,ZSHALF,DIFFU,HYDRO,                        &
              QSG,QVG,QCG,QCATM,QVATM,PRCP,                     &
              QKMS,TRANSP,DRIP,                                 &
              DEW,SMELT,SOILICE,VEGFRAC,SNOWFRAC,soilres,       &
              DQM,QMIN,REF,KSAT,RAS,INFMAX,                     & !--- soil properties
              SOILMOIS,SOILIQW,MAVAIL,RUNOFF,RUNOFF2,INFILTRP)    !--- output
!*************************************************************************
!   moisture balance equation and Richards eqn.
!   are solved here 
!   
!     DELT - time step (s)
!     IME,JME,NZS - dimensions of soil domain
!     ZSMAIN - main levels in soil (m)
!     ZSHALF - middle of the soil layers (m)
!     DTDZS -  dt/(2.*dzshalf*dzmain)
!     DTDZS2 - dt/(2.*dzshalf)
!     DIFFU - diffusional conductivity (m^2/s)
!     HYDRO - hydraulic conductivity (m/s)
!     QSG,QVG,QCG - saturated mixing ratio, mixing ratio of
!                   water vapor and cloud at the ground
!                   surface, respectively (kg/kg)
!     QCATM,QVATM - cloud and water vapor mixing ratio
!                   at the first atm. level (kg/kg)
!     PRCP - precipitation rate in m/s
!     QKMS - exchange coefficient for water vapor in the
!              surface layer (m/s)
!     TRANSP - transpiration from the soil layers (m/s)
!     DRIP - liquid water dripping from the canopy to soil (m)
!     DEW -  dew in kg/m^2s
!     SMELT - melting rate in m/s
!     SOILICE - volumetric content of ice in soil (m^3/m^3)
!     SOILIQW - volumetric content of liquid water in soil (m^3/m^3)
!     VEGFRAC - greeness fraction (0-1)
!     RAS - ration of air density to soil density
!     INFMAX - maximum infiltration rate (kg/m^2/s)
!    
!     SOILMOIS - volumetric soil moisture, 6 levels (m^3/m^3)
!     MAVAIL - fraction of maximum soil moisture in the top
!               layer (0-1)
!     RUNOFF - surface runoff (m/s)
!     RUNOFF2 - underground runoff (m)
!     INFILTRP - point infiltration flux into soil (m/s)
!            /(snow bottom runoff) (mm/s)
!
!     COSMC, RHSMC - coefficients for implicit solution of
!                     Richards equation
!******************************************************************
        IMPLICIT NONE
!------------------------------------------------------------------
!--- input variables
   LOGICAL,  INTENT(IN   )   ::  debug_print
   real (kind_phys),     INTENT(IN   )   ::  DELT
   real (kind_phys),     INTENT(IN   )   ::  xlat, xlon, testptlat, testptlon
   INTEGER,  INTENT(IN   )   ::  NZS,NDDZS

! input variables

   real (kind_phys),     DIMENSION(1:NZS), INTENT(IN) :: ZSMAIN, &
                                                         ZSHALF, &
                                                          DIFFU, &
                                                          HYDRO, &
                                                         TRANSP, &
                                                        SOILICE, &
                                                         DTDZS2

   real (kind_phys),     DIMENSION(1:NDDZS), INTENT(IN)  ::  DTDZS

   real (kind_phys),     INTENT(IN ) :: QSG,QVG,QCG,QCATM,QVATM, &
                                   QKMS,VEGFRAC,DRIP,PRCP      , &
                                   DEW,SMELT,SNOWFRAC          , &
                                   DQM,QMIN,REF,KSAT,RAS,RIW,SOILRES
                         
! output

   real (kind_phys),     DIMENSION(  1:nzs )                   , &
                         INTENT(INOUT) :: SOILMOIS,SOILIQW
                                                
   real (kind_phys),     INTENT(INOUT) :: MAVAIL,RUNOFF,RUNOFF2,INFILTRP, &
                                                        INFMAX

! local variables

   real (kind_phys),     DIMENSION( 1:nzs )  ::  COSMC,RHSMC

   real (kind_phys)    ::  DZS,R1,R2,R3,R4,R5,R6,R7,R8,R9,R10
   real (kind_phys)    ::  REFKDT,REFDK,DELT1,F1MAX,F2MAX
   real (kind_phys)    ::  F1,F2,FD,KDT,VAL,DDT,PX,FK,FKMAX
   real (kind_phys)    ::  QQ,UMVEG,INFMAX1,TRANS
   real (kind_phys)    ::  TOTLIQ,FLX,FLXSAT,QTOT
   real (kind_phys)    ::  DID,X1,X2,X4,DENOM,Q2,Q4
   real (kind_phys)    ::  dice,fcr,acrt,frzx,sum,cvfrz

   INTEGER ::  NZS1,NZS2,K,KK,K1,KN,ialp1,jj,jk

!******************************************************************************
!       COEFFICIENTS FOR THOMAS ALGORITHM FOR SOILMOIS
!******************************************************************************
          NZS1=NZS-1                                                            
          NZS2=NZS-2

 118      format(6(10Pf23.19))

           do k=1,nzs
            cosmc(k)=zero
            rhsmc(k)=zero
           enddo
 
        DID=(ZSMAIN(NZS)-ZSHALF(NZS))
        X1=ZSMAIN(NZS)-ZSMAIN(NZS1)

        DENOM=(one+DIFFU(nzs1)/X1/DID*DELT+HYDRO(NZS)/(2._kind_phys*DID)*DELT)
        COSMC(1)=DELT*(DIFFU(nzs1)/DID/X1                                 &
                    +HYDRO(NZS1)/2._kind_phys/DID)/DENOM
        RHSMC(1)=(SOILMOIS(NZS)+TRANSP(NZS)*DELT/                         &
               DID)/DENOM

!12 June 2014 - low boundary condition: 1 - zero diffusion below the lowest
! level; 2 - soil moisture at the low boundary can be lost due to the root uptake.
! So far - no interaction with the water table.

        DENOM=1.+DIFFU(nzs1)/X1/DID*DELT
        COSMC(1)=DELT*(DIFFU(nzs1)/DID/X1                                 &
                    +HYDRO(NZS1)/DID)/DENOM

        RHSMC(1)=(SOILMOIS(NZS)-HYDRO(NZS)*DELT/DID*soilmois(nzs)         &
                 +TRANSP(NZS)*DELT/DID)/DENOM

        COSMC(1)=zero
        RHSMC(1)=SOILMOIS(NZS)
!
        DO K=1,NZS2
          KN=NZS-K
          K1=2*KN-3
          X4=2.*DTDZS(K1)*DIFFU(KN-1)
          X2=2.*DTDZS(K1+1)*DIFFU(KN)
          Q4=X4+HYDRO(KN-1)*DTDZS2(KN-1)
          Q2=X2-HYDRO(KN+1)*DTDZS2(KN-1)
          DENOM=one+X2+X4-Q2*COSMC(K)
          COSMC(K+1)=Q4/DENOM
    IF (debug_print ) THEN
       if (abs(xlat-testptlat).lt.0.05 .and.                         &
           abs(xlon-testptlon).lt.0.05)then
           print *,'xlat,xlon=',xlat,xlon
           print *,'q2,soilmois(kn),DIFFU(KN),x2,HYDRO(KN+1),DTDZS2(KN-1),kn,k' &
                   ,q2,soilmois(kn),DIFFU(KN),x2,HYDRO(KN+1),DTDZS2(KN-1),kn,k
       endif
    ENDIF
          RHSMC(K+1)=(SOILMOIS(KN)+Q2*RHSMC(K)                            &
                   +TRANSP(KN)                                            &
                   /(ZSHALF(KN+1)-ZSHALF(KN))                             &
                   *DELT)/DENOM
        ENDDO

! --- MOISTURE BALANCE BEGINS HERE

          TRANS=TRANSP(1)
          UMVEG=(one-VEGFRAC)*soilres

          RUNOFF=zero
          RUNOFF2=zero
          DZS=ZSMAIN(2)
          R1=COSMC(NZS1)
          R2= RHSMC(NZS1)
          R3=DIFFU(1)/DZS
          R4=R3+HYDRO(1)*.5_kind_phys          
          R5=R3-HYDRO(2)*.5_kind_phys
          R6=QKMS*RAS
!-- Total liquid water available on the top of soil domain
!-- Without snow - 3 sources of water: precipitation,
!--         water dripping from the canopy and dew 
!-- With snow - only one source of water - snow melt

  191   format (f23.19)

        TOTLIQ=PRCP-DRIP/DELT-(one-VEGFRAC)*DEW*RAS-SMELT
    IF (debug_print ) THEN
       if (abs(xlat-testptlat).lt.0.05 .and.                         &
           abs(xlon-testptlon).lt.0.05)then
           print *,'xlat,xlon=',xlat,xlon
           print *,'UMVEG*PRCP,DRIP/DELT,UMVEG*DEW*RAS,SMELT', &
                    UMVEG*PRCP,DRIP/DELT,UMVEG*DEW*RAS,SMELT
       endif
    ENDIF

        FLX=TOTLIQ
        INFILTRP=TOTLIQ

! -----------     FROZEN GROUND VERSION    -------------------------
!   REFERENCE FROZEN GROUND PARAMETER, CVFRZ, IS A SHAPE PARAMETER OF
!   Areal (kind_phys) DISTRIBUTION FUNCTION OF SOIL ICE CONTENT WHICH EQUALS 1/CV.
!   CV IS A COEFFICIENT OF SPATIAL VARIATION OF SOIL ICE CONTENT.
!   BASED ON FIELD DATA CV DEPENDS ON Areal (kind_phys) MEAN OF FROZEN DEPTH, AND IT
!   CLOSE TO CONSTANT = 0.6 IF Areal (kind_phys) MEAN FROZEN DEPTH IS ABOVE 20 CM.
!   THAT IS WHY PARAMETER CVFRZ = 3 (INT{1/0.6*0.6})
!
!   Current logic doesn't allow CVFRZ be bigger than 3
         CVFRZ = 3._kind_phys

!-- SCHAAKE/KOREN EXPRESSION for calculation of max infiltration
         REFKDT=3._kind_phys
         REFDK=3.4341E-6_kind_phys
         DELT1=DELT/86400._kind_phys
         F1MAX=DQM*ZSHALF(2)
         F2MAX=DQM*(ZSHALF(3)-ZSHALF(2))
         F1=F1MAX*(one-SOILMOIS(1)/DQM)
         DICE=SOILICE(1)*ZSHALF(2)
         FD=F1
        do k=2,nzs1
         DICE=DICE+(ZSHALF(k+1)-ZSHALF(k))*SOILICE(K)
         FKMAX=DQM*(ZSHALF(k+1)-ZSHALF(k))
         FK=FKMAX*(one-SOILMOIS(k)/DQM)
         FD=FD+FK
        enddo
         KDT=REFKDT*KSAT/REFDK
         VAL=(1.-EXP(-KDT*DELT1))
         DDT = FD*VAL
         PX= - TOTLIQ * DELT
         IF(PX < zero) PX = zero
         IF(PX > zero) THEN
           INFMAX1 = (PX*(DDT/(PX+DDT)))/DELT
         ELSE
           INFMAX1 = zero
         ENDIF
    IF (debug_print ) THEN
      print *,'INFMAX1 before frozen part',INFMAX1
    ENDIF

! -----------     FROZEN GROUND VERSION    --------------------------
!    REDUCTION OF INFILTRATION BASED ON FROZEN GROUND PARAMETERS
!
! ------------------------------------------------------------------

         FRZX= 0.15_kind_phys*((dqm+qmin)/ref) * (0.412_kind_phys / 0.468_kind_phys)

         FCR = one
         IF ( DICE .GT. 1.E-2_kind_phys) THEN
           ACRT = CVFRZ * FRZX / DICE
           SUM = one
           IALP1 = CVFRZ - 1
           DO JK = 1,IALP1
              K = 1
              DO JJ = JK+1, IALP1
                K = K * JJ
              END DO
              SUM = SUM + (ACRT ** ( CVFRZ-JK)) / FLOAT (K)
           END DO
           FCR = one - EXP(-ACRT) * SUM
         END IF
    IF (debug_print ) THEN
          print *,'FCR--------',fcr
          print *,'DICE=',dice
    ENDIF
         INFMAX1 = INFMAX1* FCR
! -------------------------------------------------------------------

         INFMAX = MAX(INFMAX1,HYDRO(1)*SOILMOIS(1))
         INFMAX = MIN(INFMAX, -TOTLIQ)
    IF (debug_print ) THEN
      print *,'INFMAX,INFMAX1,HYDRO(1)*SOILIQW(1),-TOTLIQ', &
               INFMAX,INFMAX1,HYDRO(1)*SOILIQW(1),-TOTLIQ
    ENDIF
!----
          IF (-TOTLIQ.GT.INFMAX)THEN
            RUNOFF=-TOTLIQ-INFMAX
            FLX=-INFMAX
    IF (debug_print ) THEN
       print *,'FLX,RUNOFF1=',flx,runoff
    ENDIF
          ENDIF
! INFILTRP is total infiltration flux in M/S
          INFILTRP=FLX
! Solution of moisture budget
          R7=.5_kind_phys*DZS/DELT
          R4=R4+R7
          FLX=FLX-SOILMOIS(1)*R7
! R8 is for direct evaporation from soil, which occurs
! only from snow-free areas
          R8=UMVEG*R6*(one-snowfrac)
          QTOT=QVATM+QCATM
          R9=TRANS
          R10=QTOT-QSG

!-- evaporation regime
          IF(R10.LE.zero) THEN
            QQ=(R5*R2-FLX+R9)/(R4-R5*R1-R10*R8/(REF-QMIN))
            FLXSAT=-DQM*(R4-R5*R1-R10*R8/(REF-QMIN))                &
                   +R5*R2+R9
          ELSE
!-- dew formation regime
            QQ=(R2*R5-FLX+R8*(QTOT-QCG-QVG)+R9)/(R4-R1*R5)
            FLXSAT=-DQM*(R4-R1*R5)+R2*R5+R8*(QTOT-QVG-QCG)+R9
          END IF

          IF(QQ.LT.0.) THEN
!  print *,'negative QQ=',qq
            SOILMOIS(1)=1.e-8_kind_phys

          ELSE IF(QQ.GT.DQM) THEN
!-- saturation
            SOILMOIS(1)=DQM
    IF (debug_print ) THEN
   print *,'FLXSAT,FLX,DELT',FLXSAT,FLX,DELT,RUNOFF2
    ENDIF
            RUNOFF=RUNOFF+(FLXSAT-FLX)
          ELSE
            SOILMOIS(1)=min(dqm,max(1.e-8_kind_phys,QQ))
          END IF

    IF (debug_print ) THEN
       if (abs(xlat-testptlat).lt.0.05 .and.                         &
           abs(xlon-testptlon).lt.0.05)then
           print *,'xlat,xlon=',xlat,xlon
           print *,'SOILMOIS,SOILIQW, soilice',SOILMOIS,SOILIQW,soilice*riw
           print *,'COSMC,RHSMC',COSMC,RHSMC
       endif
    ENDIF
!--- FINAL SOLUTION FOR SOILMOIS 
!          DO K=2,NZS1
          DO K=2,NZS
            KK=NZS-K+1
            QQ=COSMC(KK)*SOILMOIS(K-1)+RHSMC(KK)

           IF (QQ.LT.zero) THEN

           ELSE IF(QQ.GT.DQM) THEN
!-- saturation
            SOILMOIS(K)=DQM
             IF(K.EQ.NZS)THEN
    IF (debug_print ) THEN
   print *,'hydro(k),QQ,DQM,k',hydro(k),QQ,DQM,k
    ENDIF
               RUNOFF2=RUNOFF2+((QQ-DQM)*(ZSMAIN(K)-ZSHALF(K)))/DELT
             ELSE
               RUNOFF2=RUNOFF2+((QQ-DQM)*(ZSHALF(K+1)-ZSHALF(K)))/DELT
             ENDIF
           ELSE
            SOILMOIS(K)=min(dqm,max(1.e-8_kind_phys,QQ))
           END IF
          END DO
    IF (debug_print ) THEN
       if (abs(xlat-testptlat).lt.0.05 .and.                         &
           abs(xlon-testptlon).lt.0.05)then
           print *,'xlat,xlon=',xlat,xlon
           print *,'END soilmois,soiliqw,soilice',soilmois,SOILIQW,soilice*riw
       endif 
    ENDIF

           MAVAIL=max(.00001_kind_phys,min(one,(SOILMOIS(1)/(REF-QMIN)*(one-snowfrac)+one*snowfrac)))
!-------------------------------------------------------------------
    END SUBROUTINE SOILMOIST
!-------------------------------------------------------------------

!>\ingroup lsm_ruc_group
!! This subroutine computes thermal diffusivity, and diffusional and 
!! hydraulic condeuctivities in soil.
            SUBROUTINE SOILPROP( debug_print,                     &
         xlat, xlon, testptlat, testptlon,                        &
         nzs,fwsat,lwsat,tav,keepfr,                              & !--- input variables
         soilmois,soiliqw,soilice,                                &
         soilmoism,soiliqwm,soilicem,                             &
         QWRTZ,rhocs,dqm,qmin,psis,bclh,ksat,                     & !--- soil fixed fields
         riw,xlmelt,CP,G0_P,cvw,ci,                               & !--- constants
         kqwrtz,kice,kwt,                                         &
         thdif,diffu,hydro,cap)                                     !--- output variables

!******************************************************************
! SOILPROP computes thermal diffusivity, and diffusional and
!          hydraulic condeuctivities
!******************************************************************
! NX,NY,NZS - dimensions of soil domain
! FWSAT, LWSAT - volumetric content of frozen and liquid water
!                for saturated condition at given temperatures (m^3/m^3)
! TAV - temperature averaged for soil layers (K)
! SOILMOIS - volumetric soil moisture at the main soil levels (m^3/m^3)
! SOILMOISM - volumetric soil moisture averaged for layers (m^3/m^3)
! SOILIQWM - volumetric liquid soil moisture averaged for layers (m^3/m^3)
! SOILICEM - volumetric content of soil ice averaged for layers (m^3/m^3)
! THDIF - thermal diffusivity for soil layers (W/m/K)
! DIFFU - diffusional conductivity (m^2/s)
! HYDRO - hydraulic conductivity (m/s)
! CAP - volumetric heat capacity (J/m^3/K)
!
!******************************************************************

        IMPLICIT NONE
!-----------------------------------------------------------------

!--- soil properties
   LOGICAL,  INTENT(IN   )   ::  debug_print
   INTEGER, INTENT(IN   )    ::                            NZS
   real (kind_phys), INTENT(IN   ) :: xlat, xlon, testptlat, testptlon

   real (kind_phys)                                            , &
            INTENT(IN   )    ::                           RHOCS, &
                                                           BCLH, &
                                                            DQM, &
                                                           KSAT, &
                                                           PSIS, &
                                                          QWRTZ, &  
                                                           QMIN

   real (kind_phys),    DIMENSION(  1:nzs )                    , &
            INTENT(IN   )    ::                        SOILMOIS, &
                                                         keepfr


   real (kind_phys),     INTENT(IN   )   ::                  CP, &
                                                            CVW, &
                                                            RIW, &  
                                                         kqwrtz, &
                                                           kice, &
                                                            kwt, &
                                                         XLMELT, &
                                                            G0_P



!--- output variables
   real (kind_phys),     DIMENSION(1:NZS)                      , &
            INTENT(INOUT)  ::      cap,diffu,hydro             , &
                                   thdif,tav                   , &
                                   soilmoism                   , &
                                   soiliqw,soilice             , &
                                   soilicem,soiliqwm           , &
                                   fwsat,lwsat

!--- local variables
   real (kind_phys),     DIMENSION(1:NZS)  ::  hk,detal,kasat,kjpl

   real (kind_phys)    ::  x,x1,x2,x4,ws,wd,fact,fach,facd,psif,ci
   real (kind_phys)    ::  tln,tavln,tn,pf,a,am,ame,h
   INTEGER ::  nzs1,k

!-- for Johansen thermal conductivity
   real (kind_phys)    ::  kzero,gamd,kdry,kas,x5,sr,ke       
               

         nzs1=nzs-1

!-- Constants for Johansen (1975) thermal conductivity
         kzero =2._kind_phys       ! if qwrtz > 0.2


         do k=1,nzs
            detal (k)=zero
            kasat (k)=zero
            kjpl  (k)=zero
            hk    (k)=zero
         enddo

           ws=dqm+qmin
           x1=xlmelt/(g0_p*psis)
           x2=x1/bclh*ws
           x4=(bclh+one)/bclh
!--- Next 3 lines are for Johansen thermal conduct.
           gamd=(one-ws)*2700._kind_phys
           kdry=(0.135_kind_phys*gamd+64.7_kind_phys)/(2700._kind_phys-0.947_kind_phys*gamd)
           !-- one more option from Christa's paper
           if(qwrtz > 0.2_kind_phys) then
             kas=kqwrtz**qwrtz*kzero**(1.-qwrtz)
           else
             kas=kqwrtz**qwrtz*3._kind_phys**(one-qwrtz)
           endif

         DO K=1,NZS1
           tn=tav(k) - tfrz
           wd=ws - riw*soilicem(k)
           psif=psis*100._kind_phys*(wd/(soiliqwm(k)+qmin))**bclh            &
                * (ws/wd)**3._kind_phys
!--- PSIF should be in [CM] to compute PF
           pf=log10(abs(psif))
           fact=one+riw*soilicem(k)
!--- HK is for McCumber thermal conductivity
         IF(PF.LE.5.2_kind_phys) THEN
           HK(K)=420._kind_phys*EXP(-(PF+2.7_kind_phys))*fact
         ELSE
           HK(K)=.1744_kind_phys*fact
         END IF

           IF(soilicem(k).NE.zero.AND.TN.LT.zero) then
!--- DETAL is taking care of energy spent on freezing or released from 
!          melting of soil water

              DETAL(K)=tfrz*X2/(TAV(K)*TAV(K))*                  &
                     (TAV(K)/(X1*TN))**X4

              if(keepfr(k).eq.one) then
                 detal(k)=zero
              endif

           ENDIF

!--- Next 10 lines calculate Johansen thermal conductivity KJPL
           kasat(k)=kas**(one-ws)*kice**fwsat(k)                  &
                    *kwt**lwsat(k)

           X5=(soilmoism(k)+qmin)/ws
         if(soilicem(k).eq.zero) then
           sr=max(0.101_kind_phys,x5)
           ke=log10(sr)+one
         else
           ke=x5
         endif

           kjpl(k)=ke*(kasat(k)-kdry)+kdry

!--- CAP -volumetric heat capacity
            CAP(K)=(one-WS)*RHOCS                                   &
                  + (soiliqwm(K)+qmin)*CVW                          &
                  + soilicem(K)*CI                                  &
                  + (dqm-soilmoism(k))*CP*1.2_kind_phys             &
            - DETAL(K)*rhowater*xlmelt

           a=RIW*soilicem(K)

        if((ws-a).lt.0.12_kind_phys)then
           diffu(K)=zero
        else
           H=max(zero,(soilmoism(K)+qmin-a)/(max(1.e-8_kind_phys,(ws-a)))) 
           facd=one
        if(a.ne.zero)facd=one-a/max(1.e-8_kind_phys,soilmoism(K))
          ame=max(1.e-8_kind_phys,ws-riw*soilicem(K))
!--- DIFFU is diffusional conductivity of soil water
          diffu(K)=-BCLH*KSAT*PSIS/ame*                            &
                  (ws/ame)**3._kind_phys                           &
                  *H**(BCLH+2._kind_phys)*facd
         endif

!--- thdif - thermal diffusivity
!           thdif(K)=HK(K)/CAP(K)
!--- Use thermal conductivity from Johansen (1975)
            thdif(K)=KJPL(K)/CAP(K)

         END DO

    IF (debug_print ) THEN
   print *,'soilice*riw,soiliqw,soilmois,ws',soilice*riw,soiliqw,soilmois,ws
    ENDIF
         DO K=1,NZS

         if((ws-riw*soilice(k)).lt.0.12_kind_phys)then
            hydro(k)=zero
         else
            fach=one
          if(soilice(k).ne.zero)                                   &
             fach=one-riw*soilice(k)/max(1.e-8_kind_phys,soilmois(k))
         am=max(1.e-8_kind_phys,ws-riw*soilice(k))
!--- HYDRO is hydraulic conductivity of soil water
          hydro(K)=min(KSAT,KSAT/am*                               & 
                  (soiliqw(K)/am)                                  &
                  **(2._kind_phys*BCLH+2._kind_phys)               &
                  * fach)
          if(hydro(k)<1.e-10_kind_phys)hydro(k)=zero
         endif

       ENDDO
    IF (debug_print ) THEN
       print *,'hydro=',hydro
    ENDIF

!-----------------------------------------------------------------------
   END SUBROUTINE SOILPROP
!-----------------------------------------------------------------------

!>\ingroup lsm_ruc_group
!> This subroutine solves the transpiration function (EQs. 18,19 in
!! Smirnova et al.(1997) \cite Smirnova_1997)
           SUBROUTINE TRANSF(  debug_print,                      &
              xlat,xlon,testptlat,testptlon,                     &
              nzs,nroot,soiliqw,tabs,lai,gswin,                  & !--- input variables
              dqm,qmin,ref,wilt,zshalf,pc,iland,                 & !--- soil fixed fields
              tranf,transum)                                       !--- output variables

!-------------------------------------------------------------------
!--- TRANF(K) - THE TRANSPIRATION FUNCTION (Smirnova et al. 1996, EQ. 18,19)
!*******************************************************************
! NX,NY,NZS - dimensions of soil domain
! SOILIQW - volumetric liquid soil moisture at the main levels (m^3/m^3)
! TRANF - the transpiration function at levels (m)
! TRANSUM - transpiration function integrated over the rooting zone (m)
!
!*******************************************************************
        IMPLICIT NONE
!-------------------------------------------------------------------

!--- input variables

   LOGICAL,  INTENT(IN   )   ::  debug_print
   INTEGER,  INTENT(IN   )   ::  nroot,nzs,iland
   real (kind_phys), INTENT(IN   ) :: xlat,xlon,testptlat,testptlon

   real (kind_phys)                                            , &
            INTENT(IN   )    ::                GSWin, TABS, lai
!--- soil properties
   real (kind_phys)                                            , &
            INTENT(IN   )    ::                             DQM, &
                                                           QMIN, &
                                                            REF, &
                                                             PC, &
                                                           WILT

   real (kind_phys),  DIMENSION(1:NZS), INTENT(IN)  :: soiliqw,  &
                                                        ZSHALF

!-- output 
   real (kind_phys),     DIMENSION(1:NZS), INTENT(OUT)  :: TRANF
   real (kind_phys),     INTENT(OUT)  ::                 TRANSUM  

!-- local variables
   real (kind_phys)    ::  totliq, did
   INTEGER ::  k

!-- for non-linear root distribution
   real (kind_phys)    ::  gx,sm1,sm2,sm3,sm4,ap0,ap1,ap2,ap3,ap4
   real (kind_phys)    ::  FTEM, PCtot, fsol, f1, cmin, cmax, totcnd
   real (kind_phys),     DIMENSION(1:NZS)   ::           PART
!--------------------------------------------------------------------

        do k=1,nzs
           part(k)=zero
           tranf(k)=zero
        enddo

        transum=zero
        totliq=soiliqw(1)+qmin
           sm1=totliq
           sm2=sm1*sm1
           sm3=sm2*sm1
           sm4=sm3*sm1
           ap0=0.299_kind_phys
           ap1=-8.152_kind_phys
           ap2=61.653_kind_phys
           ap3=-115.876_kind_phys
           ap4=59.656_kind_phys
           gx=ap0+ap1*sm1+ap2*sm2+ap3*sm3+ap4*sm4
          if(totliq.ge.ref) gx=one
          if(totliq.le.wilt) gx=zero
          if(gx.gt.one) gx=one
          if(gx.lt.zero) gx=zero
        DID=zshalf(2)
          part(1)=DID*gx
        IF(TOTLIQ.GT.REF) THEN
          TRANF(1)=DID
        ELSE IF(TOTLIQ.LE.WILT) THEN
          TRANF(1)=zero
        ELSE
          TRANF(1)=(TOTLIQ-WILT)/(REF-WILT)*DID
        ENDIF 
!-- uncomment next line for non-linear root distribution
          !TRANF(1)=part(1)

        DO K=2,NROOT
        totliq=soiliqw(k)+qmin
           sm1=totliq
           sm2=sm1*sm1
           sm3=sm2*sm1
           sm4=sm3*sm1
           gx=ap0+ap1*sm1+ap2*sm2+ap3*sm3+ap4*sm4
          if(totliq.ge.ref) gx=one
          if(totliq.le.wilt) gx=zero
          if(gx.gt.one) gx=one
          if(gx.lt.zero) gx=zero
          DID=zshalf(K+1)-zshalf(K)
          part(k)=did*gx
        IF(totliq.GE.REF) THEN
          TRANF(K)=DID
        ELSE IF(totliq.LE.WILT) THEN
          TRANF(K)=zero
        ELSE
          TRANF(K)=(totliq-WILT)                                &
                /(REF-WILT)*DID
        ENDIF
!-- uncomment next line for non-linear root distribution
          !TRANF(k)=part(k)
        END DO
    IF (debug_print ) THEN
       if (abs(xlat-testptlat).lt.0.05 .and.                         &
           abs(xlon-testptlon).lt.0.05)then
         print *,'xlat,xlon=',xlat,xlon
         print *,'soiliqw =',soiliqw,'wilt=',wilt,'qmin= ',qmin
         print *,'tranf = ',tranf
       endif
    ENDIF

! For LAI> 3 =>  transpiration at potential rate (F.Tardieu, 2013)
      if(lai > 4._kind_phys) then
        pctot=0.8_kind_phys
      else
        pctot=pc
!- 26aug16-  next 2 lines could lead to LH increase and higher 2-m Q during the day
!        pctot=min(0.8,pc*lai)
!        pctot=min(0.8,max(pc,pc*lai))
      endif
    IF ( debug_print ) THEN
       if (abs(xlat-testptlat).lt.0.05 .and.                         &
           abs(xlon-testptlon).lt.0.05)then
           print *,'xlat,xlon=',xlat,xlon
           print *,'pctot,lai,pc',pctot,lai,pc
       endif
    ENDIF
!---
!--- air temperature function
!     Avissar (1985) and AX 7/95
        IF (TABS .LE. 302.15_kind_phys) THEN
          FTEM = one / (one + EXP(-0.41_kind_phys * (TABS - 282.05_kind_phys)))
        ELSE
          FTEM = one / (one + EXP(0.5_kind_phys * (TABS - 314.0_kind_phys)))
        ENDIF
!--- incoming solar function
     cmin = one/rsmax_data
     cmax = one/rstbl(iland)
    if(lai > one) then
     cmax = lai/rstbl(iland) ! max conductance
    endif
! Noihlan & Planton (1988)
       f1=zero
!    if(lai > 0.01) then
!       f1 = 1.1/lai*gswin/rgltbl(iland)! f1=0. when GSWin=0.
!       fsol = (f1+cmin/cmax)/(1.+f1)
!       fsol=min(1.,fsol)
!    else
!       fsol=cmin/cmax
!    endif
!     totcnd = max(lai/rstbl(iland), pctot * ftem * f1) 
! Mahrer & Avissar (1982), Avissar et al. (1985)
     if (GSWin < rgltbl(iland)) then
      fsol = one / (one + exp(-0.034_kind_phys * (GSWin - 3.5_kind_phys)))
     else
      fsol = one
     endif
!--- total conductance
     totcnd =(cmin + (cmax - cmin)*pctot*ftem*fsol)/cmax

    IF ( debug_print ) THEN
       if (abs(xlat-testptlat).lt.0.05 .and.                         &
           abs(xlon-testptlon).lt.0.05)then
          print *,'xlat,xlon=',xlat,xlon
          print *,'GSWin,Tabs,lai,f1,cmax,cmin,pc,pctot,ftem,fsol',GSWin,Tabs,lai,f1,cmax,cmin,pc,pctot,ftem,fsol
          print *,'iland,RGLTBL(iland),RSTBL(iland),RSMAX_DATA,totcnd'  &
                  ,iland,RGLTBL(iland),RSTBL(iland),RSMAX_DATA,totcnd
       endif
    ENDIF

!-- TRANSUM - total for the rooting zone
          transum=zero
        DO K=1,NROOT
! linear root distribution
         TRANF(k)=max(zero,TRANF(k)*totcnd)
         transum=transum+tranf(k)
        END DO
    IF ( debug_print ) THEN
       if (abs(xlat-testptlat).lt.0.05 .and.                         &
           abs(xlon-testptlon).lt.0.05)then
         print *,'xlat,xlon=',xlat,xlon
         print *,'transum,TRANF',transum,tranf
       endif
    ENDIF

!-----------------------------------------------------------------
   END SUBROUTINE TRANSF
!-----------------------------------------------------------------

!>\ingroup lsm_ruc_group
!> This subroutine finds the solution of energy budget at the surface
!! from the pre-computed table of saturated water vapor mixing ratio 
!! and estimated surface temperature.
       SUBROUTINE VILKA(TN,D1,D2,PP,QS,TS,TT,NSTEP,ii,j,iland,isoil,xlat,xlon)
!--------------------------------------------------------------
!--- VILKA finds the solution of energy budget at the surface
!--- using table T,QS computed from Clausius-Klapeiron
!--------------------------------------------------------------
   real (kind_phys),     DIMENSION(1:5001),  INTENT(IN   )   ::  TT
   real (kind_phys),     INTENT(IN  )   ::  TN,D1,D2,PP,xlat,xlon
   INTEGER,  INTENT(IN  )   ::  NSTEP,ii,j,iland,isoil

   real (kind_phys),     INTENT(OUT  )  ::  QS, TS

   real (kind_phys)    ::  F1,T1,T2,RN
   INTEGER ::  I,I1

       I=(TN-1.7315E2_kind_dbl_prec)/.05_kind_dbl_prec+1
       T1=173.1_kind_dbl_prec+FLOAT(I)*.05_kind_dbl_prec
       F1=T1+D1*TT(I)-D2
       I1=I-F1/(.05_kind_dbl_prec+D1*(TT(I+1)-TT(I)))
       I=I1
       IF(I.GT.5000.OR.I.LT.1) GOTO 1
  10   I1=I
       T1=173.1_kind_dbl_prec+FLOAT(I)*.05_kind_dbl_prec
       F1=T1+D1*TT(I)-D2
       RN=F1/(.05_kind_dbl_prec+D1*(TT(I+1)-TT(I)))
       I=I-INT(RN)
       IF(I.GT.5000.OR.I.LT.1) GOTO 1
       IF(I1.NE.I) GOTO 10
       TS=T1-.05_kind_dbl_prec*RN
       QS=(TT(I)+(TT(I)-TT(I+1))*RN)/PP
       GOTO 20
   1   PRINT *,'     AVOST IN VILKA     Table index= ',I
       print *,'I,J=',ii,j,'LU_index = ',iland, 'Psfc[hPa] = ',pp, 'Tsfc = ',tn
       print *,'AVOST point at xlat/xlon=',xlat,xlon
   20  CONTINUE
!-----------------------------------------------------------------------
   END SUBROUTINE VILKA
!-----------------------------------------------------------------------

!>\ingroup lsm_ruc_group
!! This subroutine computes effective land and soil parameters in the
!! grid cell from the weighted contribution of soil and land categories
!! represented in the grid cell.
     SUBROUTINE SOILVEGIN  ( debug_print,mosaic_lu,mosaic_soil,        &
                             soilfrac,nscat,shdmin, shdmax,            &
                     NLCAT,IVGTYP,ISLTYP,iswater,MYJ,                  &
                     IFOREST,lufrac,vegfrac,EMISS,PC,                  &
                     MSNF,FACSNF,ZNT,LAI,RDLAI2D,                      &
                     QWRTZ,RHOCS,BCLH,DQM,KSAT,PSIS,QMIN,REF,WILT,I,J, &
                     errmsg, errflg)

!************************************************************************
!  Set-up soil and vegetation Parameters in the case when
!  snow disappears during the forecast and snow parameters
!  shold be replaced by surface parameters according to
!  soil and vegetation types in this point.
!
!        Output:
!
!
!             Soil parameters:
!               DQM: MAX soil moisture content - MIN (m^3/m^3)
!               REF:        Reference soil moisture (m^3/m^3)
!               WILT: Wilting PT soil moisture contents (m^3/m^3)
!               QMIN: Air dry soil moist content limits (m^3/m^3)
!               PSIS: SAT soil potential coefs. (m)
!               KSAT:  SAT soil diffusivity/conductivity coefs. (m/s)
!               BCLH: Soil diffusivity/conductivity exponent.
!
! ************************************************************************

   IMPLICIT NONE
!---------------------------------------------------------------------------
      integer,   parameter      ::      nsoilclas=19
      integer,   parameter      ::      nvegclas=24+3
      integer,   parameter      ::      ilsnow=99

   LOGICAL,    INTENT(IN   )    ::      debug_print
   INTEGER,    INTENT(IN   )    ::      mosaic_lu, mosaic_soil
   INTEGER,    INTENT(IN   )    ::      nlcat, nscat, iswater, i, j

!---    soiltyp classification according to STATSGO(nclasses=16)
!
!             1          SAND                  SAND
!             2          LOAMY SAND            LOAMY SAND
!             3          SANDY LOAM            SANDY LOAM
!             4          SILT LOAM             SILTY LOAM
!             5          SILT                  SILTY LOAM
!             6          LOAM                  LOAM
!             7          SANDY CLAY LOAM       SANDY CLAY LOAM
!             8          SILTY CLAY LOAM       SILTY CLAY LOAM
!             9          CLAY LOAM             CLAY LOAM
!            10          SANDY CLAY            SANDY CLAY
!            11          SILTY CLAY            SILTY CLAY
!            12          CLAY                  LIGHT CLAY
!            13          ORGANIC MATERIALS     LOAM
!            14          WATER
!            15          BEDROCK
!                        Bedrock is reclassified as class 14
!            16          OTHER (land-ice)
!            17          Playa
!            18          Lava
!            19          White Sand
!
!----------------------------------------------------------------------
         real (kind_phys)  LQMA(nsoilclas),LRHC(nsoilclas),           &
               LPSI(nsoilclas),LQMI(nsoilclas),                       &
               LBCL(nsoilclas),LKAS(nsoilclas),                       &
               LWIL(nsoilclas),LREF(nsoilclas),                       &
               DATQTZ(nsoilclas)
!-- LQMA Rawls et al.[1982]
!        DATA LQMA /0.417, 0.437, 0.453, 0.501, 0.486, 0.463, 0.398,
!     &  0.471, 0.464, 0.430, 0.479, 0.475, 0.439, 1.0, 0.20, 0.401/
!---
!-- Clapp, R. and G. Hornberger, 1978: Empirical equations for some soil
!   hydraulic properties, Water Resour. Res., 14, 601-604.

!-- Clapp et al. [1978]
     DATA LQMA /0.395, 0.410, 0.435, 0.485, 0.485, 0.451, 0.420,      &
                0.477, 0.476, 0.426, 0.492, 0.482, 0.451, 1.0,        &
                0.20,  0.435, 0.468, 0.200, 0.339/

!-- LREF Rawls et al.[1982]
!        DATA LREF /0.091, 0.125, 0.207, 0.330, 0.360, 0.270, 0.255,
!     &  0.366, 0.318, 0.339, 0.387, 0.396, 0.329, 1.0, 0.108, 0.283/

!-- Clapp et al. [1978]
        DATA LREF /0.174, 0.179, 0.249, 0.369, 0.369, 0.314, 0.299,   &
                   0.357, 0.391, 0.316, 0.409, 0.400, 0.314, 1.,      &
                   0.1,   0.249, 0.454, 0.17,  0.236/

!-- LWIL Rawls et al.[1982]
!        DATA LWIL/0.033, 0.055, 0.095, 0.133, 0.133, 0.117, 0.148,
!     &  0.208, 0.197, 0.239, 0.250, 0.272, 0.066, 0.0, 0.006, 0.029/

!-- Clapp et al. [1978]
        DATA LWIL/0.068, 0.075, 0.114, 0.179, 0.179, 0.155, 0.175,    &
                  0.218, 0.250, 0.219, 0.283, 0.286, 0.155, 0.0,      &
                  0.006, 0.114, 0.030, 0.006, 0.01/

!        DATA LQMI/0.010, 0.028, 0.047, 0.084, 0.084, 0.066, 0.067,
!     &  0.120, 0.103, 0.100, 0.126, 0.138, 0.066, 0.0, 0.006, 0.028/

!-- Carsel and Parrish [1988]
        DATA LQMI/0.045, 0.057, 0.065, 0.067, 0.034, 0.078, 0.10,     &
                  0.089, 0.095, 0.10,  0.070, 0.068, 0.078, 0.0,      &
                  0.004, 0.065, 0.020, 0.004, 0.008/

!-- LPSI Cosby et al[1984]
!        DATA LPSI/0.060, 0.036, 0.141, 0.759, 0.759, 0.355, 0.135,
!     &  0.617, 0.263, 0.098, 0.324, 0.468, 0.355, 0.0, 0.069, 0.036/
!     &  0.617, 0.263, 0.098, 0.324, 0.468, 0.355, 0.0, 0.069, 0.036/

!-- Clapp et al. [1978]
       DATA LPSI/0.121, 0.090, 0.218, 0.786, 0.786, 0.478, 0.299,     &
                 0.356, 0.630, 0.153, 0.490, 0.405, 0.478, 0.0,       &
                 0.121, 0.218, 0.468, 0.069, 0.069/

!-- LKAS Rawls et al.[1982]
!        DATA LKAS/5.83E-5, 1.70E-5, 7.19E-6, 1.89E-6, 1.89E-6,
!     &  3.67E-6, 1.19E-6, 4.17E-7, 6.39E-7, 3.33E-7, 2.50E-7,
!     &  1.67E-7, 3.38E-6, 0.0, 1.41E-4, 1.41E-5/

!-- Clapp et al. [1978]
        DATA LKAS/1.76E-4, 1.56E-4, 3.47E-5, 7.20E-6, 7.20E-6,         &
                  6.95E-6, 6.30E-6, 1.70E-6, 2.45E-6, 2.17E-6,         &
                  1.03E-6, 1.28E-6, 6.95E-6, 0.0,     1.41E-4,         &
                  3.47E-5, 1.28E-6, 1.41E-4, 1.76E-4/

!-- LBCL Cosby et al [1984]
!        DATA LBCL/2.79,  4.26,  4.74,  5.33,  5.33,  5.25,  6.66,
!     &  8.72,  8.17,  10.73, 10.39, 11.55, 5.25,  0.0, 2.79,  4.26/

!-- Clapp et al. [1978]
        DATA LBCL/4.05,  4.38,  4.90,  5.30,  5.30,  5.39,  7.12,      &
                  7.75,  8.52, 10.40, 10.40, 11.40,  5.39,  0.0,       &
                  4.05,  4.90, 11.55,  2.79,  2.79/

        DATA LRHC /1.47,1.41,1.34,1.27,1.27,1.21,1.18,1.32,1.23,       &
                   1.18,1.15,1.09,1.21,4.18,2.03,2.10,1.09,2.03,1.47/

        DATA DATQTZ/0.92,0.82,0.60,0.25,0.10,0.40,0.60,0.10,0.35,      &
                    0.52,0.10,0.25,0.00,0.,0.60,0.0,0.25,0.60,0.92/

!--------------------------------------------------------------------------
!
!     USGS Vegetation Types
!
!    1:   Urban and Built-Up Land
!    2:   Dryland Cropland and Pasture
!    3:   Irrigated Cropland and Pasture
!    4:   Mixed Dryland/Irrigated Cropland and Pasture
!    5:   Cropland/Grassland Mosaic
!    6:   Cropland/Woodland Mosaic
!    7:   Grassland
!    8:   Shrubland
!    9:   Mixed Shrubland/Grassland
!   10:   Savanna
!   11:   Deciduous Broadleaf Forest
!   12:   Deciduous Needleleaf Forest
!   13:   Evergreen Broadleaf Forest
!   14:   Evergreen Needleleaf Fores
!   15:   Mixed Forest
!   16:   Water Bodies
!   17:   Herbaceous Wetland
!   18:   Wooded Wetland
!   19:   Barren or Sparsely Vegetated
!   20:   Herbaceous Tundra
!   21:   Wooded Tundra
!   22:   Mixed Tundra
!   23:   Bare Ground Tundra
!   24:   Snow or Ice
!
!   25:   Playa
!   26:   Lava
!   27:   White Sand

! MODIS vegetation categories from VEGPARM.TBL
!    1:   Evergreen Needleleaf Forest
!    2:   Evergreen Broadleaf Forest
!    3:   Deciduous Needleleaf Forest
!    4:   Deciduous Broadleaf Forest
!    5:   Mixed Forests
!    6:   Closed Shrublands
!    7:   Open Shrublands
!    8:   Woody Savannas
!    9:   Savannas
!   10:   Grasslands
!   11:   Permanent wetlands
!   12:   Croplands
!   13:   Urban and Built-Up
!   14:   cropland/natural vegetation mosaic
!   15:   Snow and Ice
!   16:   Barren or Sparsely Vegetated
!   17:   Water
!   18:   Wooded Tundra
!   19:   Mixed Tundra
!   20:   Barren Tundra
!   21:   Lakes


!----  Below are the arrays for the vegetation parameters
         real (kind_phys) LALB(nvegclas),LMOI(nvegclas),LEMI(nvegclas), &
              LROU(nvegclas),LTHI(nvegclas),LSIG(nvegclas),             &
              LPC(nvegclas)

!************************************************************************
!----     vegetation parameters
!
!-- USGS model
!
        DATA  LALB/.18,.17,.18,.18,.18,.16,.19,.22,.20,.20,.16,.14,     &
                   .12,.12,.13,.08,.14,.14,.25,.15,.15,.15,.25,.55,     &
                   .30,.16,.60 /
        DATA LEMI/.88,4*.92,.93,.92,.88,.9,.92,.93,.94,                 &
                  .95,.95,.94,.98,.95,.95,.85,.92,.93,.92,.85,.95,      &
                  .85,.85,.90 /
!-- Roughness length is changed for forests and some others
         DATA LROU/.5,.06,.075,.065,.05,.2,.075,.1,.11,.15,.5,.5,       & 
                   .5,.5,.5,.0001,.2,.4,.05,.1,.15,.1,.065,.05,         &
                   .01,.15,.01 /

        DATA LMOI/.1,.3,.5,.25,.25,.35,.15,.1,.15,.15,.3,.3,            &
                  .5,.3,.3,1.,.6,.35,.02,.5,.5,.5,.02,.95,.40,.50,.40/
!
!---- still needs to be corrected
!
       DATA LPC /0.4,0.3,0.4,0.4,0.4,0.4,0.4,0.4,0.4,0.4,5*0.55,0.,0.55,0.55,                   &
                 0.3,0.3,0.4,0.4,0.3,0.,.3,0.,0./
!***************************************************************************


   INTEGER      ::                &
                                                         IVGTYP, &
                                                         ISLTYP

   LOGICAL,    INTENT(IN   )    ::     myj
   real (kind_phys),       INTENT(IN )      ::   SHDMAX
   real (kind_phys),       INTENT(IN )      ::   SHDMIN
   real (kind_phys),       INTENT(IN )      ::   VEGFRAC
   real (kind_phys),     DIMENSION( 1:NLCAT ),  INTENT(IN)::         LUFRAC
   real (kind_phys),     DIMENSION( 1:NSCAT ),  INTENT(IN)::         SOILFRAC

   real (kind_phys)                                            , &
            INTENT (  OUT)            ::                     pc, &
                                                           msnf, &
                                                         facsnf

   real (kind_phys)                                            , &
            INTENT (INOUT   )         ::                  emiss, &
                                                            lai, &
                                                            znt
  LOGICAL, intent(in) :: rdlai2d
!--- soil properties
   real (kind_phys)                                            , &
            INTENT(  OUT)    ::                           RHOCS, &
                                                           BCLH, &
                                                            DQM, &
                                                           KSAT, &
                                                           PSIS, &
                                                           QMIN, &
                                                          QWRTZ, &
                                                            REF, &
                                                           WILT
   INTEGER, INTENT (  OUT)   ::                         iforest
   character(len=*),intent(out) :: errmsg
   integer,         intent(out) :: errflg
   INTEGER   ::   kstart, kfin, lstart, lfin
   INTEGER   ::   k
   real (kind_phys)      ::   area,  factor, znt1, lb
   real (kind_phys),     DIMENSION( 1:NLCAT ) :: ZNTtoday, LAItoday, deltalai

!***********************************************************************
!        DATA ZS1/0.0,0.05,0.20,0.40,1.6,3.0/   ! o -  levels in soil
!        DATA ZS2/0.0,0.025,0.125,0.30,1.,2.3/   ! x - levels in soil


   ! Initialize error-handling
   errflg = 0
   errmsg = ''

   iforest = IFORTBL(IVGTYP)

    IF (debug_print ) THEN
        print *,'ifortbl(ivgtyp),ivgtyp,laitbl(ivgtyp),z0tbl(ivgtyp)', &
            ifortbl(ivgtyp),ivgtyp,laitbl(ivgtyp),z0tbl(ivgtyp)
    ENDIF

        deltalai(:) = zero

! 11oct2012 - seasonal correction on ZNT for crops and LAI for all veg. types
! factor = 1 with minimum greenness -->  vegfrac = shdmin (cold season)
! factor = 0 with maximum greenness -->  vegfrac = shdmax
! SHDMAX, SHDMIN and VEGFRAC are in % here.
      if((shdmax - shdmin) .lt. one) then
        factor = one ! min greenness
      else
        factor = one - max(zero,min(one,(vegfrac - shdmin)/max(one,(shdmax-shdmin))))
      endif

! 18sept18 - LAITBL and Z0TBL are the max values
      do k = 1,nlcat
       if(IFORTBL(k) == 1) deltalai(k)=min(0.2_kind_phys,0.8_kind_phys*LAITBL(K))
       if(IFORTBL(k) == 2 .or. IFORTBL(k) == 7) deltalai(k)=min(0.5_kind_phys,0.8_kind_phys*LAITBL(K))
       if(IFORTBL(k) == 3) deltalai(k)=min(0.45_kind_phys,0.8_kind_phys*LAITBL(K))
       if(IFORTBL(k) == 4) deltalai(k)=min(0.75_kind_phys,0.8_kind_phys*LAITBL(K))
       if(IFORTBL(k) == 5) deltalai(k)=min(0.86_kind_phys,0.8_kind_phys*LAITBL(K))

       if(k.ne.iswater) then
!-- 20aug18 - change in LAItoday based on the greenness fraction for the current day
        LAItoday(k) = LAITBL(K) - deltalai(k) * factor

         if(IFORTBL(k) == 7) then
!-- seasonal change of roughness length for crops 
           ZNTtoday(k) = Z0TBL(K) - 0.125_kind_phys * factor
         else
           ZNTtoday(k) = Z0TBL(K)
         endif
       else
        LAItoday(k) = LAITBL(K)
        ZNTtoday(k) = ZNT ! do not overwrite z0 over water with the table value
       endif
      enddo

    IF (debug_print ) THEN
        print *,'ivgtyp,factor,vegfrac,shdmin,shdmax,deltalai,laitoday(ivgtyp),znttoday(ivgtyp)', &
         i,j,ivgtyp,factor,vegfrac,shdmin,shdmax,deltalai(ivgtyp),laitoday(ivgtyp),znttoday(ivgtyp)
    ENDIF

        EMISS = zero
        ZNT   = zero
        ZNT1  = zero
        PC    = zero
        MSNF  = zero
        FACSNF= zero
        if(.not.rdlai2d) LAI = zero
        AREA  = zero
!-- mosaic approach to landuse in the grid box
! Use  Mason (1988) Eq.(15) to compute effective ZNT;
!  Lb - blending height =  L/200., where L is the length scale
! of regions with varying Z0 (Lb = 5 if L=1000 m)
        LB = 5._kind_phys
      if(mosaic_lu == 1) then
      do k = 1,nlcat
        AREA  = AREA + lufrac(k)
        EMISS = EMISS+ LEMITBL(K)*lufrac(k)
        ZNT   = ZNT  + lufrac(k)/ALOG(LB/ZNTtoday(K))**2._kind_phys
! ZNT1 - weighted average in the grid box, not used, computed for comparison
        ZNT1  = ZNT1 + lufrac(k)*ZNTtoday(K)
        if(.not.rdlai2d) LAI = LAI  + LAItoday(K)*lufrac(k)
        PC    = PC   + PCTBL(K)*lufrac(k)
        MSNF  = MSNF + MFSNO(K)*lufrac(k)
        FACSNF= FACSNF + SNCOVFAC(K)*lufrac(k)
      enddo

       if (area.gt.one) area=one
       if (area <= zero) then
          print *,'Bad area of grid box', area
          errflg = 1
          errmsg = 'ERROR(SOILVEGIN): Bad area of grid box'
          return
       endif

    IF (debug_print ) THEN
        print *,'area=',area,i,j,ivgtyp,nlcat,(lufrac(k),k=1,nlcat),EMISS,ZNT,ZNT1,LAI,PC
    ENDIF

        EMISS = EMISS/AREA
        ZNT1   = ZNT1/AREA
        ZNT = LB/EXP(SQRT(one/ZNT))
        if(.not.rdlai2d) LAI = LAI/AREA
        PC    = PC /AREA
        MSNF  = MSNF /AREA
        FACSNF= FACSNF /AREA

    IF (debug_print ) THEN 
        print *,'mosaic=',j,ivgtyp,nlcat,(lufrac(k),k=1,nlcat),EMISS,ZNT,ZNT1,LAI,PC
    ENDIF


      else
        EMISS = LEMITBL(IVGTYP)
        ZNT = ZNTtoday(IVGTYP)
        PC    = PCTBL(IVGTYP)
        MSNF  = MFSNO(IVGTYP)
        FACSNF= SNCOVFAC(IVGTYP)
        if(.not.rdlai2d) LAI = LAItoday(IVGTYP)
     endif

! parameters from SOILPARM.TBL
          RHOCS  = zero
          BCLH   = zero
          DQM    = zero
          KSAT   = zero
          PSIS   = zero
          QMIN   = zero
          REF    = zero
          WILT   = zero
          QWRTZ  = zero
          AREA   = zero
! mosaic approach
       if(mosaic_soil == 1 ) then
            do k = 1, nscat
        if(k.ne.14) then  ! STATSGO value for water
        !exclude water points from this loop
          AREA   = AREA + soilfrac(k)
          RHOCS  = RHOCS + HC(k)*1.E6_kind_phys*soilfrac(k)
          BCLH   = BCLH + BB(K)*soilfrac(k)
          DQM    = DQM + (MAXSMC(K)-                               &
                   DRYSMC(K))*soilfrac(k)
          KSAT   = KSAT + SATDK(K)*soilfrac(k)
          PSIS   = PSIS - SATPSI(K)*soilfrac(k)
          QMIN   = QMIN + DRYSMC(K)*soilfrac(k)
          REF    = REF + REFSMC(K)*soilfrac(k)
          WILT   = WILT + WLTSMC(K)*soilfrac(k)
          QWRTZ  = QWRTZ + QTZ(K)*soilfrac(k)
        endif
            enddo
       if (area.gt.one) area=one
       if (area <= zero) then
! area = 0. for water points
!          print *,'Area of a grid box', area, 'iswater = ',iswater
          RHOCS  = HC(ISLTYP)*1.E6_kind_phys
          BCLH   = BB(ISLTYP)
          DQM    = MAXSMC(ISLTYP)-                               &
                   DRYSMC(ISLTYP)
          KSAT   = SATDK(ISLTYP)
          PSIS   = - SATPSI(ISLTYP)
          QMIN   = DRYSMC(ISLTYP)
          REF    = REFSMC(ISLTYP)
          WILT   = WLTSMC(ISLTYP)
          QWRTZ  = QTZ(ISLTYP)
       else
          RHOCS  = RHOCS/AREA
          BCLH   = BCLH/AREA
          DQM    = DQM/AREA
          KSAT   = KSAT/AREA
          PSIS   = PSIS/AREA
          QMIN   = QMIN/AREA
          REF    = REF/AREA
          WILT   = WILT/AREA
          QWRTZ  = QWRTZ/AREA
       endif

! dominant category approach
        else
      if(isltyp.ne.14) then
          RHOCS  = HC(ISLTYP)*1.E6_kind_phys
          BCLH   = BB(ISLTYP)
          DQM    = MAXSMC(ISLTYP)-                               &
                   DRYSMC(ISLTYP)
          KSAT   = SATDK(ISLTYP)
          PSIS   = - SATPSI(ISLTYP)
          QMIN   = DRYSMC(ISLTYP)
          REF    = REFSMC(ISLTYP)
          WILT   = WLTSMC(ISLTYP)
          QWRTZ  = QTZ(ISLTYP)
      endif
        endif

!--------------------------------------------------------------------------
   END SUBROUTINE SOILVEGIN
!--------------------------------------------------------------------------

!>\ingroup lsm_ruc_group
!> This subroutine computes liquid and forezen soil moisture from the
!! total soil moisture, and also computes soil moisture availability in
!! the top soil layer.
  SUBROUTINE RUCLSMINIT( debug_print, landfrac, fice, min_seaice,  &
                     nzs, isltyp, ivgtyp, mavail,                  &
                     sh2o, smfr3d, tslb, smois,                    &
                     ims,ime, jms,jme, kms,kme,                    &
                     its,ite, jts,jte, kts,kte                     )

  use namelist_soilveg_ruc

#if ( WRF_CHEM == 1 )
  USE module_data_gocart_dust
#endif
   IMPLICIT NONE
   LOGICAL,  INTENT(IN   )   ::  debug_print
   real (kind_phys), DIMENSION( ims:ime),  INTENT(IN   )   :: landfrac, fice
   real (kind_phys),                       INTENT(IN   )   :: min_seaice

   INTEGER,  INTENT(IN   )   ::     &
                                    ims,ime, jms,jme, kms,kme,  &
                                    its,ite, jts,jte, kts,kte,  &
                                    nzs

   real (kind_phys), DIMENSION( ims:ime, 1:nzs, jms:jme )     , &
            INTENT(IN)    ::                              TSLB, &
                                                         SMOIS

   INTEGER, DIMENSION( ims:ime, jms:jme )                      , &
            INTENT(INOUT)    ::                   ISLTYP,IVGTYP

   real (kind_phys), DIMENSION( ims:ime, 1:nzs, jms:jme )      , &
            INTENT(OUT)    ::                            SMFR3D, &
                                                         SH2O

   real (kind_phys), DIMENSION( ims:ime, jms:jme )             , &
            INTENT(OUT)    ::                            MAVAIL

   !-- local
   real (kind_phys), DIMENSION ( 1:nzs ) :: SOILIQW

   INTEGER ::  I,J,L,itf,jtf
   real (kind_phys)    ::  RIW,XLMELT,TLN,DQM,REF,PSIS,QMIN,BCLH

   INTEGER                   :: errflag

        RIW=rhoice*1.e-3_kind_phys
        XLMELT=con_hfus

! for FIM
   itf=ite  !  min0(ite,ide-1)
   jtf=jte  !  min0(jte,jde-1)

   errflag = 0
   DO j = jts,jtf
     DO i = its,itf

       IF ( ISLTYP( i,j ) .LT. 0 ) THEN
         errflag = 1
         print *, &
         "module_sf_ruclsm.F: lsminit: out of range ISLTYP ",i,j,ISLTYP( i,j )
       ENDIF
     ENDDO
   ENDDO
   IF ( errflag .EQ. 1 ) THEN
      print *,&
      "module_sf_ruclsm.F: lsminit: out of range value "// &
                            "of ISLTYP. Is this field in the input?" 
   ENDIF

   DO J=jts,jtf
     DO I=its,itf

       ! in Zobler classification isltyp=0 for water. Statsgo classification
       ! has isltyp=14 for water
       if (isltyp(i,j) == 0) isltyp(i,j)=14
     
       if(landfrac(i) > zero ) then
       !-- land
       !-- Computate volumetric content of ice in soil
       !-- and initialize MAVAIL
         DQM    = MAXSMC   (ISLTYP(I,J)) - &
                  DRYSMC   (ISLTYP(I,J))
         REF    = REFSMC   (ISLTYP(I,J))
         PSIS   = - SATPSI (ISLTYP(I,J))
         QMIN   = DRYSMC   (ISLTYP(I,J))
         BCLH   = BB       (ISLTYP(I,J))

         mavail(i,j) = max(0.00001_kind_phys,min(one,(smois(i,1,j)-qmin)/(ref-qmin)))

         DO L=1,NZS
           !-- for land points initialize soil ice
           tln=log(TSLB(i,l,j)/tfrz)
          
           if(tln.lt.zero) then
             soiliqw(l)=(dqm+qmin)*(XLMELT*                        &
            (tslb(i,l,j)-tfrz)/tslb(i,l,j)/grav/psis)              &
            **(-one/bclh)
            soiliqw(l)=max(zero,soiliqw(l))
            soiliqw(l)=min(soiliqw(l),smois(i,l,j))
            sh2o(i,l,j)=soiliqw(l)
            smfr3d(i,l,j)=(smois(i,l,j)-soiliqw(l))/RIW
         
           else
             smfr3d(i,l,j)=zero
             sh2o(i,l,j)=smois(i,l,j)
           endif
         ENDDO

       elseif( fice(i) > min_seaice) then
       !-- ice
         mavail(i,j) = one
         DO L=1,NZS
           smfr3d(i,l,j)=one
           sh2o(i,l,j)=zero
         ENDDO
    
       else
       !-- water  ISLTYP=14
         mavail(i,j) = one
         DO L=1,NZS
           smfr3d(i,l,j)=zero
           sh2o(i,l,j)=one
         ENDDO

       endif ! land

    ENDDO
   ENDDO


  END SUBROUTINE ruclsminit
!
!-----------------------------------------------------------------
!>\ingroup lsm_ruc_group
!> This subroutine specifies vegetation related characteristics.
        SUBROUTINE RUCLSM_SOILVEGPARM( debug_print,MMINLURUC, MMINSL)
!-----------------------------------------------------------------

        IMPLICIT NONE
        LOGICAL,  INTENT(IN   )   ::  debug_print

        integer :: LUMATCH, IINDEX, LC, NUM_SLOPE
        integer :: ierr
        INTEGER , PARAMETER :: OPEN_OK = 0

        character*8 :: MMINLURUC, MMINSL
        character*128 ::  vege_parm_string
!        logical, external :: wrf_dm_on_monitor


!-----SPECIFY VEGETATION RELATED CHARACTERISTICS :
!             ALBBCK: SFC albedo (in percentage)
!                 Z0: Roughness length (m)
!               LEMI: Emissivity
!                 PC: Plant coefficient for transpiration function
! -- the rest of the parameters are read in but not used currently
!             SHDFAC: Green vegetation fraction (in percentage)
!  Note: The ALBEDO, Z0, and SHDFAC values read from the following table
!          ALBEDO, amd Z0 are specified in LAND-USE TABLE; and SHDFAC is
!          the monthly green vegetation data
!             CMXTBL: MAX CNPY Capacity (m)
!              RSMIN: Mimimum stomatal resistance (s m-1)
!              RSMAX: Max. stomatal resistance (s m-1)
!                RGL: Parameters used in radiation stress function
!                 HS: Parameter used in vapor pressure deficit functio
!               TOPT: Optimum transpiration air temperature. (K)
!             CMCMAX: Maximum canopy water capacity
!             CFACTR: Parameter used in the canopy inteception calculati
!               SNUP: Threshold snow depth (in water equivalent m) that
!                     implies 100% snow cover
!                LAI: Leaf area index (dimensionless)
!             MAXALB: Upper bound on maximum albedo over deep snow
!
!-----READ IN VEGETAION PROPERTIES FROM VEGPARM.TBL 
!                                                                       

!       IF ( wrf_dm_on_monitor() ) THEN

        OPEN(19, FILE='VEGPARM.TBL',FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
        IF(ierr .NE. OPEN_OK ) THEN
          print *,&
          'module_sf_ruclsm.F: soil_veg_gen_parm: failure opening VEGPARM.TBL'
        END IF

          print *,&
         'INPUT VEGPARM FOR ',MMINLURUC

        LUMATCH=0

 2000   FORMAT (A8)
!sms$serial begin
        READ (19,'(A)') vege_parm_string
!sms$serial end
        outer : DO 
!sms$serial begin
           READ (19,2000,END=2002)LUTYPE
           READ (19,*)LUCATS,IINDEX
!sms$serial end

            print *,&
           'VEGPARM FOR ',LUTYPE,' FOUND', LUCATS,' CATEGORIES'

           IF(LUTYPE.NE.MMINLURUC)THEN    ! Skip over the undesired table
           print *,&
              'Skipping ', LUTYPE, ' table'
              DO LC=1,LUCATS
!sms$serial begin
                 READ (19,*)
!sms$serial end
              ENDDO
              inner : DO               ! Find the next "Vegetation Parameters"
!sms$serial begin
                 READ (19,'(A)',END=2002) vege_parm_string
!sms$serial end
                 IF (TRIM(vege_parm_string) .EQ. "Vegetation Parameters") THEN
                    EXIT inner
                 END IF
               ENDDO inner
           ELSE
              LUMATCH=1
              print *,&
              'Found ', LUTYPE, ' table'
              EXIT outer                ! Found the table, read the data
           END IF

        ENDDO outer

        IF (LUMATCH == 1) then
           print *,&
           'Reading ',LUTYPE,' table'
           DO LC=1,LUCATS
!sms$serial begin
              READ (19,*)IINDEX,ALBTBL(LC),Z0TBL(LC),LEMITBL(LC),PCTBL(LC), &
                         SHDTBL(LC),IFORTBL(LC),RSTBL(LC),RGLTBL(LC),         &
                         HSTBL(LC),SNUPTBL(LC),LAITBL(LC),MAXALB(LC)
!sms$serial end
           ENDDO
!
!sms$serial begin
           READ (19,*)
           READ (19,*)TOPT_DATA
           READ (19,*)
           READ (19,*)CMCMAX_DATA
           READ (19,*)
           READ (19,*)CFACTR_DATA
           READ (19,*)
           READ (19,*)RSMAX_DATA
           READ (19,*)
           READ (19,*)BARE
           READ (19,*)
           READ (19,*)GLACIER
           READ (19,*)
           READ (19,*)NATURAL
           READ (19,*)
           READ (19,*)CROP
           READ (19,*)
           READ (19,*,iostat=ierr)URBAN
!sms$serial end
           if ( ierr /= 0 )  print *, "-------- VEGPARM.TBL READ ERROR --------"
           if ( ierr /= 0 )  print *, "Problem read URBAN from VEGPARM.TBL"
           if ( ierr /= 0 )  print *, " -- Use updated version of VEGPARM.TBL  "
           if ( ierr /= 0 )  print *,  "Problem read URBAN from VEGPARM.TBL"

        ENDIF

 2002   CONTINUE
        CLOSE (19)
!-----
    IF (debug_print ) THEN
         print *,' LEMITBL, PCTBL, Z0TBL, LAITBL --->', LEMITBL, PCTBL, Z0TBL, LAITBL
    ENDIF


        IF (LUMATCH == 0) then
!           CALL wrf_error_fatal ("Land Use Dataset '"//MMINLURUC//"' not found in VEGPARM.TBL.")
        ENDIF

!-----READ IN SOIL PROPERTIES FROM SOILPARM.TBL
!                                                                       
!      IF ( wrf_dm_on_monitor() ) THEN
        OPEN(19, FILE='SOILPARM.TBL',FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
        IF(ierr .NE. OPEN_OK ) THEN
          print *,&
          'module_sf_ruclsm.F: soil_veg_gen_parm: failure opening SOILPARM.TBL'
        END IF

        print *,'INPUT SOIL TEXTURE CLASSIFICATION = ',MMINSL

        LUMATCH=0

!sms$serial begin
        READ (19,'(A)') vege_parm_string
!sms$serial end
        outersl : DO
!sms$serial begin
           READ (19,2000,END=2003)SLTYPE
           READ (19,*)SLCATS,IINDEX
!sms$serial end

            print *,&
           'SOILPARM FOR ',SLTYPE,' FOUND', SLCATS,' CATEGORIES'

           IF(SLTYPE.NE.MMINSL)THEN    ! Skip over the undesired table
           print *,&
              'Skipping ', SLTYPE, ' table'
              DO LC=1,SLCATS
!sms$serial begin
                 READ (19,*)
!sms$serial end
              ENDDO
              innersl : DO               ! Find the next "Vegetation Parameters"
!sms$serial begin
                 READ (19,'(A)',END=2002) vege_parm_string
!sms$serial end
                 IF (TRIM(vege_parm_string) .EQ. "Soil Parameters") THEN
                    EXIT innersl
                 END IF
               ENDDO innersl
           ELSE
              LUMATCH=1
              print *,&
              'Found ', SLTYPE, ' table'
              EXIT outersl                ! Found the table, read the data
           END IF

        ENDDO outersl

        IF (LUMATCH == 1) then
     print *,'SLCATS=',SLCATS
          DO LC=1,SLCATS
!sms$serial begin
              READ (19,*) IINDEX,BB(LC),DRYSMC(LC),HC(LC),MAXSMC(LC),&
                        REFSMC(LC),SATPSI(LC),SATDK(LC), SATDW(LC),   &
                        WLTSMC(LC), QTZ(LC)
 !sms$serial end
          ENDDO
         ENDIF

 2003   CONTINUE

        CLOSE (19)

      IF(LUMATCH.EQ.0)THEN
          print *, 'SOIl TEXTURE IN INPUT FILE DOES NOT ' 
          print *, 'MATCH SOILPARM TABLE'                 
          print *, 'INCONSISTENT OR MISSING SOILPARM FILE' 
      ENDIF

!
!-----READ IN GENERAL PARAMETERS FROM GENPARM.TBL 
!                                                                       
        OPEN(19, FILE='GENPARM.TBL',FORM='FORMATTED',STATUS='OLD',IOSTAT=ierr)
        IF(ierr .NE. OPEN_OK ) THEN
          print *,&
          'module_sf_ruclsm.F: soil_veg_gen_parm: failure opening GENPARM.TBL'
        END IF

!sms$serial begin
        READ (19,*)
        READ (19,*)
        READ (19,*) NUM_SLOPE
!sms$serial end

          SLPCATS=NUM_SLOPE

          DO LC=1,SLPCATS
!sms$serial begin
              READ (19,*)SLOPE_DATA(LC)
!sms$serial end
          ENDDO

!sms$serial begin
          READ (19,*)
          READ (19,*)SBETA_DATA
          READ (19,*)
          READ (19,*)FXEXP_DATA
          READ (19,*)
          READ (19,*)CSOIL_DATA
          READ (19,*)
          READ (19,*)SALP_DATA
          READ (19,*)
          READ (19,*)REFDK_DATA
          READ (19,*)
          READ (19,*)REFKDT_DATA
          READ (19,*)
          READ (19,*)FRZK_DATA
          READ (19,*)
          READ (19,*)ZBOT_DATA
          READ (19,*)
          READ (19,*)CZIL_DATA
          READ (19,*)
          READ (19,*)SMLOW_DATA
          READ (19,*)
          READ (19,*)SMHIGH_DATA
!sms$serial end
        CLOSE (19)

!-----------------------------------------------------------------
      END SUBROUTINE RUCLSM_SOILVEGPARM
!-----------------------------------------------------------------

!>\ingroup lsm_ruc_group
!> This subroutine specifies 19 soiltyp classification according to
!! STATSGO.
  SUBROUTINE SOILIN (ISLTYP, DQM, REF, PSIS, QMIN, BCLH )

!---    soiltyp classification according to STATSGO(nclasses=16)
!
!             1          SAND                  SAND
!             2          LOAMY SAND            LOAMY SAND
!             3          SANDY LOAM            SANDY LOAM
!             4          SILT LOAM             SILTY LOAM
!             5          SILT                  SILTY LOAM
!             6          LOAM                  LOAM
!             7          SANDY CLAY LOAM       SANDY CLAY LOAM
!             8          SILTY CLAY LOAM       SILTY CLAY LOAM
!             9          CLAY LOAM             CLAY LOAM
!            10          SANDY CLAY            SANDY CLAY
!            11          SILTY CLAY            SILTY CLAY
!            12          CLAY                  LIGHT CLAY
!            13          ORGANIC MATERIALS     LOAM
!            14          WATER
!            15          BEDROCK
!                        Bedrock is reclassified as class 14
!            16          OTHER (land-ice)
! extra classes from Fei Chen
!            17          Playa
!            18          Lava
!            19          White Sand
!
!----------------------------------------------------------------------
         integer,   parameter      ::      nsoilclas=19

         integer, intent ( in)  ::                          isltyp
         real,    intent ( out) ::               dqm,ref,qmin,psis,bclh

         real (kind_phys)  LQMA(nsoilclas),LREF(nsoilclas),LBCL(nsoilclas),       &
               LPSI(nsoilclas),LQMI(nsoilclas)

!-- LQMA Rawls et al.[1982]
!        DATA LQMA /0.417, 0.437, 0.453, 0.501, 0.486, 0.463, 0.398,
!     &  0.471, 0.464, 0.430, 0.479, 0.475, 0.439, 1.0, 0.20, 0.401/
!---
!-- Clapp, R. and G. Hornberger, Empirical equations for some soil
!   hydraulic properties, Water Resour. Res., 14,601-604,1978.
!-- Clapp et al. [1978]
     DATA LQMA /0.395, 0.410, 0.435, 0.485, 0.485, 0.451, 0.420,      &
                0.477, 0.476, 0.426, 0.492, 0.482, 0.451, 1.0,        &
                0.20,  0.435, 0.468, 0.200, 0.339/

!-- Clapp et al. [1978]
        DATA LREF /0.174, 0.179, 0.249, 0.369, 0.369, 0.314, 0.299,   &
                   0.357, 0.391, 0.316, 0.409, 0.400, 0.314, 1.,      &
                   0.1,   0.249, 0.454, 0.17,  0.236/

!-- Carsel and Parrish [1988]
        DATA LQMI/0.045, 0.057, 0.065, 0.067, 0.034, 0.078, 0.10,     &
                  0.089, 0.095, 0.10,  0.070, 0.068, 0.078, 0.0,      &
                  0.004, 0.065, 0.020, 0.004, 0.008/

!-- Clapp et al. [1978]
       DATA LPSI/0.121, 0.090, 0.218, 0.786, 0.786, 0.478, 0.299,     &
                 0.356, 0.630, 0.153, 0.490, 0.405, 0.478, 0.0,       &
                 0.121, 0.218, 0.468, 0.069, 0.069/

!-- Clapp et al. [1978]
        DATA LBCL/4.05,  4.38,  4.90,  5.30,  5.30,  5.39,  7.12,      &
                  7.75,  8.52, 10.40, 10.40, 11.40,  5.39,  0.0,       &
                  4.05,  4.90, 11.55,  2.79,  2.79/


          DQM    = LQMA(ISLTYP)-                               &
                   LQMI(ISLTYP)
          REF    = LREF(ISLTYP)
          PSIS   = - LPSI(ISLTYP)
          QMIN   = LQMI(ISLTYP)
          BCLH   = LBCL(ISLTYP)

  END SUBROUTINE SOILIN

!+---+-----------------------------------------------------------------+
!>\ingroup lsm_ruc_group
!> This function calculates the liquid saturation vapor mixing ratio as 
!! a function of temperature and pressure (from Thompson scheme).
      real (kind_phys) FUNCTION RSLF(P,T)

      IMPLICIT NONE
      real (kind_phys), INTENT(IN):: P, T
      real (kind_phys):: ESL,X
      real (kind_phys), PARAMETER:: C0= .611583699E03
      real (kind_phys), PARAMETER:: C1= .444606896E02
      real (kind_phys), PARAMETER:: C2= .143177157E01
      real (kind_phys), PARAMETER:: C3= .264224321E-1
      real (kind_phys), PARAMETER:: C4= .299291081E-3
      real (kind_phys), PARAMETER:: C5= .203154182E-5
      real (kind_phys), PARAMETER:: C6= .702620698E-8
      real (kind_phys), PARAMETER:: C7= .379534310E-11
      real (kind_phys), PARAMETER:: C8=-.321582393E-13

      X=MAX(-80._kind_dbl_prec,T-273.16_kind_dbl_prec)

      ESL=C0+X*(C1+X*(C2+X*(C3+X*(C4+X*(C5+X*(C6+X*(C7+X*C8)))))))
      ESL=MIN(ESL, P*0.15_kind_dbl_prec)        ! Even with P=1050mb and T=55C, the sat. vap.  pres only contributes to ~15% of total pres.
      RSLF=.622_kind_dbl_prec*ESL/max(1.e-4_kind_dbl_prec,(P-ESL))

      END FUNCTION RSLF


END MODULE module_sf_ruclsm
